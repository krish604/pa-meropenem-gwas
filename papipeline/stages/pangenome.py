"""Stage 7 - pan-genome.

Panaroo-compatible module. The heavy lifting (graph construction, core
alignment) is delegated to Panaroo when installed; this module owns the
gene presence/absence matrix and the core/accessory partition, and can build
the partition from annotations alone for TEST mode.

Outputs: ``gene_presence_absence.tsv``, ``core_genes.tsv``,
``accessory_genes.tsv``, ``pangenome_summary.tsv``.

"accessory" means **total - core** here, not the shell bucket. On the verified
10-isolate cohort that is 5,191 against a shell of 2,850; see `partition`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError, ModeNotAllowedError, StageError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import AnnotationRecord, RunMode

LOGGER = get_logger("stages.pangenome")

GPA_COLUMNS: Tuple[str, ...] = ("gene", "sample_id", "present")
CORE_COLUMNS: Tuple[str, ...] = ("gene",)
ACCESSORY_COLUMNS: Tuple[str, ...] = ("gene", "n_samples")
SUMMARY_COLUMNS: Tuple[str, ...] = ("metric", "value")

#: A gene present in at least this fraction of samples is "core".
#:
#: Exactly 1.0, i.e. present in *every* isolate. panaroo's own `Core` bucket is
#: >= 99% of strains, so the two agree only while the cohort is small enough that
#: no achievable fraction lands between 99% and 100% - at n=10 both mean 10/10.
#: At n=967 panaroo's core includes every gene in >= 958 isolates while this
#: includes only the 967/967 genes, so the two diverge materially. Reconciling
#: them is a scientific decision, not a threshold tweak, and has not been made.
CORE_PRESENCE_THRESHOLD = 1.0


@dataclass(frozen=True)
class Pangenome:
    """The pan-genome partition and the presence/absence matrix."""

    sample_ids: Tuple[str, ...]
    genes: Tuple[str, ...]
    presence: Mapping[str, Set[str]]
    core: Tuple[str, ...]
    accessory: Tuple[str, ...]

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def long_rows(self) -> List[Dict[str, object]]:
        """One row per (gene, sample) pair - the stage's output contract.

        BREAKING CHANGE. Replaces the cohort-level `gene, n_samples, frequency`
        projection, which was lossy in a way that mattered: spec D6 requires gene
        presence/absence as a pyseer feature family, and pyseer consumes a
        per-sample matrix. Two different cohorts can share every count in an
        `n_samples`/`frequency` row while differing completely in which isolates
        carry the gene. The per-sample mapping was already held here
        (`self.presence`); only the projection discarded it.

        Long rather than wide because a wide gene x sample table grows a column
        per isolate and its header stops being readable. pyseer wants the wide
        `.Rtab` matrix (`pyseer --pres`, "as produced by roary"), so that
        transpose belongs in D6 input construction, not in this stage.
        """
        rows: List[Dict[str, object]] = []
        for gene in self.genes:
            carriers = self.presence.get(gene, set())
            for sample_id in self.sample_ids:
                rows.append(
                    {
                        "gene": gene,
                        "sample_id": sample_id,
                        "present": int(sample_id in carriers),
                    }
                )
        return rows

    def matrix_rows(self) -> List[Dict[str, object]]:
        """Cohort-level counts per gene. NOT the output contract.

        Retained because these counts are the natural way to read the table by
        eye, and because `pangenome_summary` reports the same numbers. Nothing in
        the pipeline consumes this shape - see `long_rows`.
        """
        rows: List[Dict[str, object]] = []
        for gene in self.genes:
            carriers = self.presence.get(gene, set())
            rows.append(
                {
                    "gene": gene,
                    "n_samples": len(carriers),
                    "frequency": round(len(carriers) / self.n_samples, 6)
                    if self.n_samples
                    else None,
                }
            )
        return rows

    def sample_row(self, sample_id: str) -> Dict[str, int]:
        """Binary gene x sample vector for one sample."""
        return {gene: int(sample_id in self.presence.get(gene, set())) for gene in self.genes}


def build_from_annotations(
    annotations: Mapping[str, Sequence[AnnotationRecord]],
    sample_ids: Optional[Sequence[str]] = None,
    core_threshold: float = CORE_PRESENCE_THRESHOLD,
) -> Pangenome:
    """Build a pan-genome from standardised annotations.

    A gene counts as present in a sample when the sample has a named
    annotation for it. Unnamed features are excluded from the matrix
    because they cannot be matched across samples.
    """
    ids = tuple(sample_ids) if sample_ids is not None else tuple(sorted(annotations))
    presence: Dict[str, Set[str]] = {}

    for sample_id in ids:
        for record in annotations.get(sample_id, ()):  # missing -> empty
            name = record.gene_name
            if not name:
                continue
            presence.setdefault(name, set()).add(sample_id)

    genes = tuple(sorted(presence))
    n = len(ids)
    core: List[str] = []
    accessory: List[str] = []
    for gene in genes:
        if n and len(presence[gene]) / n >= core_threshold:
            core.append(gene)
        else:
            accessory.append(gene)

    return Pangenome(
        sample_ids=ids,
        genes=genes,
        presence=presence,
        core=tuple(core),
        accessory=tuple(accessory),
    )


def read_gene_presence_absence(path: Path) -> Tuple[Tuple[str, ...], Dict[str, Set[str]]]:
    """Read the long-format gene presence/absence table.

    One row per (gene, sample) pair; `present` is 1 when the isolate carries the
    gene and 0 when it does not, and 1.0 / true are also accepted.

    Columns are located **by name**, not by position. Reading positionally is how
    a metadata column gets mistaken for a genome, which yields a core count that
    is wrong but entirely plausible.

    Returns:
        ``(sample_ids, {gene: set(carrier_sample_ids)})`` where `sample_ids` is
        ordered by first appearance in the table.

    Changed with the writer in the same commit. See
    tests/unit/test_pangenome_gpa_roundtrip.py - a reader and a writer that each
    look reasonable and jointly disagree is the failure this module has already
    had once.
    """
    path = Path(path)
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not lines:
        raise DataContractError("Gene presence/absence matrix is empty", path=str(path))

    header = [h.strip() for h in lines[0].split("\t")]
    for required in GPA_COLUMNS:
        if required not in header:
            raise DataContractError(
                f"Gene presence/absence table is missing the {required!r} column; "
                f"header was {header!r}. Declared contract is {list(GPA_COLUMNS)!r}.",
                path=str(path),
            )
    gene_i = header.index("gene")
    sample_i = header.index("sample_id")
    present_i = header.index("present")

    order: List[str] = []
    seen = set()
    presence: Dict[str, Set[str]] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) <= max(gene_i, sample_i, present_i):
            raise DataContractError(
                f"Gene presence/absence row has {len(fields)} fields, fewer than "
                f"the {len(GPA_COLUMNS)} declared columns {list(GPA_COLUMNS)!r}. "
                f"A short row means the file is truncated, not that the gene is "
                f"absent.",
                path=str(path),
            )
        gene = fields[gene_i].strip()
        sample_id = fields[sample_i].strip()
        if not gene or not sample_id:
            continue
        if sample_id not in seen:
            seen.add(sample_id)
            order.append(sample_id)
        if fields[present_i].strip().lower() in {"1", "1.0", "true"}:
            presence.setdefault(gene, set()).add(sample_id)

    # A gene present in zero samples still belongs in the mapping, with an empty
    # carrier set. Dropping it would make an absent gene indistinguishable from a
    # gene nobody has looked for.
    for line in lines[1:]:
        fields = line.split("\t")
        gene = fields[gene_i].strip() if len(fields) > gene_i else ""
        if gene:
            presence.setdefault(gene, set())

    return tuple(order), presence


def partition(
    sample_ids: Sequence[str], presence: Mapping[str, Set[str]]
) -> Pangenome:
    """Partition genes into core and accessory from a presence mapping.

    **Definition, which the emitted metrics depend on and which is not
    negotiable by accident:**

    - `core` = genes present in **every** sample (`CORE_PRESENCE_THRESHOLD`
      is exactly 1.0, not 0.99).
    - `accessory` = **every other gene**, i.e. `total - core`.

    So `n_accessory_genes` is *not* the shell bucket. panaroo splits the
    non-core remainder three ways - soft-core (>= 95% of strains), shell
    (>= 15% and < 95%, i.e. 2-9 genomes at n=10) and cloud (< 15%, i.e.
    genome-specific) - and this stage collapses all three into `accessory`.
    On the verified 10-isolate cohort: total 10,019, core 4,828, accessory
    5,191, of which shell is 2,850 and cloud 2,341. Describing those 5,191 as
    "accessory (2-9 genomes)" is wrong by 2,341; the 2-9 figure is 2,850.

    This stage emits no shell/cloud/soft-core count, so a reader of
    `pangenome_summary.tsv` cannot recover the split from it. panaroo's own
    split is available separately via
    `adapters.panaroo.parse_summary_statistics`.
    """
    genes = tuple(sorted(presence))
    n = len(sample_ids)
    core = tuple(g for g in genes if n and len(presence[g]) / n >= CORE_PRESENCE_THRESHOLD)
    accessory = tuple(g for g in genes if g not in core)
    return Pangenome(
        sample_ids=tuple(sample_ids),
        genes=genes,
        presence=dict(presence),
        core=core,
        accessory=accessory,
    )


def write_outputs(pangenome: Pangenome, out_dir: Path) -> Dict[str, Path]:
    """Write the four declared pan-genome outputs."""
    out_dir = Path(out_dir)
    paths = {
        "gene_presence_absence": out_dir / "gene_presence_absence.tsv",
        "core_genes": out_dir / "core_genes.tsv",
        "accessory_genes": out_dir / "accessory_genes.tsv",
        "pangenome_summary": out_dir / "pangenome_summary.tsv",
    }

    write_tsv(paths["gene_presence_absence"], pangenome.long_rows(), GPA_COLUMNS)
    write_tsv(paths["core_genes"], [{"gene": g} for g in pangenome.core], CORE_COLUMNS)
    write_tsv(
        paths["accessory_genes"],
        [
            {"gene": g, "n_samples": len(pangenome.presence.get(g, set()))}
            for g in pangenome.accessory
        ],
        ACCESSORY_COLUMNS,
    )

    total = len(pangenome.genes)
    write_tsv(
        paths["pangenome_summary"],
        [
            {"metric": "n_samples", "value": str(pangenome.n_samples)},
            {"metric": "n_genes_total", "value": str(total)},
            {"metric": "n_core_genes", "value": str(len(pangenome.core))},
            # total - core. NOT the shell bucket. See `partition`.
            {"metric": "n_accessory_genes", "value": str(len(pangenome.accessory))},
            {
                "metric": "core_fraction",
                "value": f"{len(pangenome.core)/total:.6f}" if total else "0",
            },
        ],
        SUMMARY_COLUMNS,
    )
    return paths


def _build_from_panaroo(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
) -> Pangenome:
    """REAL: panaroo's own presence/absence table, not the annotation table.

    The order is load-bearing and each step is a separate refusal:

    1. **The gate.** `runtime.allow_real_mode` is false in every committed
       overlay, so this is what keeps a REAL run from starting by accident. It is
       checked *first* because otherwise an operator whose gate is shut is told to
       install panaroo - which is not the problem, and fixing it would not have
       unblocked anything.
    2. **The preflight.** panaroo shells out to `cd-hit`, so `panaroo` can be
       importable while the cohort cannot be processed at all. Checked before the
       parse so that a well-formed table left on disk by an earlier run does not
       read as "this stage works".
    3. **The parse**, which refuses any sample the cohort does not contain.

    A *gate*, not a refusal of the stage: open the gate and this returns a
    partition. That distinction is load-bearing in
    `tests/integration/test_stage_taxonomy_is_self_verifying.py`, which decides
    from this function's AST whether `pangenome` refuses REAL by design, and
    which a `ModeNotAllowedError` here - rather than a mode-guarded
    `NotImplementedError` - keeps on the correct side of.
    """
    if not bool(config.runtime.get("allow_real_mode", False)):
        raise ModeNotAllowedError(
            "REAL-mode stage 7 (pangenome) is gated: runtime.allow_real_mode is "
            f"false in the machine overlay for {config.machine_name!r}. panaroo "
            "is provisioned out-of-band and its output is read from this run's "
            "intermediate directory; neither happens unless REAL mode is "
            "authorised. Set it in the overlay, or open it for one session with "
            "PIPELINE_ALLOW_REAL_MODE=1.",
            machine=config.machine_name,
        )

    from ..adapters import panaroo as panaroo_adapter

    panaroo_adapter.preflight()

    panaroo_dir = panaroo_adapter.output_dir(intermediate_root)
    samples, presence = panaroo_adapter.build_pangenome_from_panaroo(
        panaroo_dir, manifest.sample_ids
    )
    LOGGER.info("Stage 7: read panaroo's presence table from %s", panaroo_dir)
    return partition(samples, presence)


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    annotations: Optional[Mapping[str, Sequence[AnnotationRecord]]] = None,
) -> Pangenome:
    """Stage 7 entry point.

    ``TEST``/``STUB`` build the partition from the standardised annotation table,
    where a "gene" is a literal annotation gene name. ``REAL`` reads panaroo's own
    table, where a "gene" is an ortholog group id - a different thing under the
    same column name, and the reason the two paths cannot be one.

    ``write_outputs`` is not called here: `run_pipeline` owns it, so the declared
    four outputs are written in one place for both modes.

    Raises:
        StageError: if the partition has no genes. A zero-gene pan-genome is a
            missing input in every mode, not a finding - see
            `_refuse_an_empty_pangenome`.
    """
    if mode is RunMode.REAL:
        pangenome = _build_from_panaroo(config, manifest, intermediate_root)
    else:
        if annotations is None:
            from .annotation import load_from_intermediate

            annotations = load_from_intermediate(intermediate_root, manifest)
        pangenome = build_from_annotations(annotations, manifest.sample_ids)

    _refuse_an_empty_pangenome(pangenome, mode, intermediate_root)

    LOGGER.info(
        "Stage 7: %d genes | %d core | %d accessory | %d samples",
        len(pangenome.genes),
        len(pangenome.core),
        len(pangenome.accessory),
        pangenome.n_samples,
    )
    return pangenome


def _refuse_an_empty_pangenome(
    pangenome: Pangenome, mode: RunMode, source_root: Path
) -> None:
    """Refuse a zero-gene partition, in every mode.

    `docs/scientific_rules.md` section 9, on empty results:

        "An empty input file is a `DataContractError`, so a truncated or failed
        run is never mistaken for a cohort with no findings."

    A zero-gene pan-genome is that failure wearing a different hat. `write_outputs`
    renders it as a perfectly well-formed `pangenome_summary.tsv` reading
    `n_genes_total 0`, `n_core_genes 0`, `n_accessory_genes 0` - which reads as
    "this cohort has no genes", a biological claim the run has no evidence for.
    Every way of reaching zero genes is a broken input:

    - the annotation tables are not where the stage was told to read (the caller
      passed the run's own output directory instead of `tool_output_root`, which
      is how this stage reported 0/0/0 over a cohort whose committed fixture
      declares 10 genes);
    - the annotation tables are there but every feature has an empty `gene_name`,
      and unnamed features cannot be matched across isolates;
    - panaroo wrote a header-only presence table, i.e. it did not finish.

    None of those is "no findings", so all of them refuse. A `StageError` and not
    `DataContractError`: the input table itself is well-formed, it is the *stage*
    that cannot produce a partition from it, and `run.py` already names the stage
    when it wraps a failure (`_execute_stage`). The message follows the house
    template - what is missing, then what to check - as in
    `adapters/iqtree.require_alignment_upstream` and the missing-tree refusal in
    `stages/similarity`.

    Deliberately not gated on mode: the gate above (`allow_real_mode`) is a *gate*,
    not a refusal, and `test_stage_taxonomy_is_self_verifying.py` reads this
    module's AST to keep `pangenome` off `REAL_REFUSING_STAGES`. This raise is
    unguarded by `allow_real_mode`, which is what that check looks for, so the
    classification is unchanged - but the guard is on the *partition being
    empty*, not on the mode, so it applies to TEST too. A TEST run that has lost
    its fixtures should say so.
    """
    if pangenome.genes:
        return

    source_root = Path(source_root)
    if mode is RunMode.REAL:
        from ..adapters import panaroo as panaroo_adapter

        expected = source_root / "panaroo" / "gene_presence_absence.csv"
        remedy = (
            f"Check that {expected} was written and carries at least one gene "
            "row: panaroo wrote a header-only table, which means it did not "
            "finish, not that the cohort has no genes. Re-run the provisioning "
            "command and read panaroo's own log for the failure. "
            + panaroo_adapter.provision_instructions(source_root / "panaroo")
        )
    else:
        expected = source_root / "annotation" / "<sample_id>.annotation.tsv"
        remedy = (
            f"Check that {expected} exists for at least one cohort member, and "
            "that the rows carry a non-empty `gene_name` - unnamed features are "
            "excluded from the matrix because they cannot be matched across "
            "isolates, so a table of unnamed features also yields zero genes. "
            "Also check this stage was handed the tool output root and not the "
            "run's own output directory: `loader.tool_output_root` exists "
            "precisely so a run cannot read its own outputs as inputs."
        )

    raise StageError(
        f"Stage 7 (pangenome) built no pan-genome: 0 genes across "
        f"{pangenome.n_samples} samples, so there is nothing to partition into "
        f"core and accessory. {remedy} A zero-gene pan-genome is a missing "
        "input, not a cohort with no findings - scientific_rules section 9 "
        "requires an empty input to refuse rather than report zero, because "
        "0/0/0 in `pangenome_summary.tsv` is indistinguishable from a result.",
        stage="pangenome",
        mode=getattr(mode, "value", mode),
        source_root=str(source_root),
        n_samples=pangenome.n_samples,
    )
