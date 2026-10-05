"""Stage 7 REAL caller: panaroo, and the parser for what it emits.

Separated from ``papipeline.stages.pangenome`` because this is the only part that
shells out. The stage owns the contract; this owns the process and the file
format panaroo happens to produce.

Preflight before anything runs, naming the missing tool. panaroo needs two
things, and they fail differently:

  - ``panaroo`` itself must be importable. On this platform it is installed from
    a pinned git source via pip, because the bioconda recipe cannot solve on
    osx-arm64 (hard prokka dependency; see docs/environment-arm64.md section 3).
  - ``cd-hit`` must be callable. panaroo shells out to it, so it is a *runtime*
    dependency and panaroo can be importable while cd-hit is absent. That is a
    real failure mode and the reason both are checked rather than one.

Deliberately NOT wired to ``MachineConfig.tool_available``. That flag exists but
has no consumer outside config and tests - the DAG-build refusal for an
unavailable tool is tickets 04 and 11, a separate scope. Checking the real
filesystem is the honest preflight; consulting a declarative flag that nothing
enforces would be theatre.
"""

from __future__ import annotations

import csv
import io
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from ..errors import DataContractError, ToolNotAvailableError
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.panaroo")

#: panaroo's own presence/absence CSV. Column 1 is the gene id, column 2 a
#: non-unique name, column 3 free-text annotation, then one column per genome.
#: Those first three are metadata, NOT genomes - reading them as genomes is the
#: easy way to get a core count wrong, and it is why `presence_columns` locates
#: the genome columns by header rather than by position.
PRESENCE_CSV = "gene_presence_absence.csv"

#: panaroo's filtered core-gene alignment - stage 8's input. Named here because
#: stage 8 reads it and its own adapter (`adapters.gubbins`) owns the constant;
#: two owners would be two answers.
CORE_ALIGNMENT_ALN = "core_gene_alignment_filtered.aln"

#: The command that produces this directory, for the refusal messages.
#:
#: **This pipeline never runs panaroo.** There is no `subprocess` call in this
#: module and none in `stages/pangenome.py`; panaroo's output is
#: operator-provisioned, so every refusal that fires because it is absent has to
#: say how to produce it. A refusal that names a missing file and not the
#: command leaves the operator to work out which tool and which invocation, and
#: "panaroo" alone is not an answer - the flags below are the ones that matter.
#:
#: **The flags are verified, not guessed.** This is the invocation recorded in
#: `docs/environment-arm64.md` section 3, run end-to-end against real Bakta
#: GFF3 for the 10-isolate smoke cohort (exit 0, 10,019 gene families):
#:
#:     panaroo -i *.gff3 -o out -t 4 --clean-mode moderate --remove-invalid-genes
#:
#: `--remove-invalid-genes` is required, not cosmetic: panaroo *raises*
#: `ValueError: Invalid gene sequence!` on any CDS whose length is not a
#: multiple of 3, is under 34 bp, or contains an internal stop. On real Bakta
#: data that is 10 genes out of 63,359. Without the flag the command fails.
#:
#: `-o` is pointed at the run's *intermediate root*, not at the panaroo
#: subdirectory: panaroo writes into `<-o>/panaroo/`, and
#: :func:`output_dir` is defined as `<intermediate_root>/panaroo`. Pointing `-o`
#: at the subdirectory itself is the off-by-one-directory that leaves the
#: pipeline looking for files one level away from where they landed.
PROVISION_COMMAND = (
    "panaroo -i <this run's intermediate>/annotation/*.gff3 \\\n"
    "    -o <this run's intermediate> \\\n"
    "    -t <threads from the machine overlay> \\\n"
    "    --clean-mode moderate \\\n"
    "    --remove-invalid-genes"
)


def provision_instructions(output_dir: Optional[Path] = None) -> str:
    """The command and the two files, resolved against a real output directory.

    Assembled here rather than written into each message so the three refusals
    that can fire on the same missing input cannot drift apart - and so the paths
    in the message are computed from the same :func:`output_dir` the parser reads,
    not spelled out again.

    Args:
        output_dir: The directory the pipeline reads panaroo's output from. Omit
            it where the caller has no run in hand - the preflight runs before
            any path is resolved - and the directory is rendered as the operator
            would recognise it.
    """
    rendered = "<this run's intermediate>/panaroo" if output_dir is None else str(
        Path(output_dir)
    )
    return (
        "This pipeline does not run panaroo: "
        "papipeline/adapters/panaroo.py contains no subprocess call, so the "
        "output is provisioned out of band. Produce it with\n"
        f"  {PROVISION_COMMAND}\n"
        f"which writes {rendered}/ into existence. The pipeline reads two "
        f"files from it: {PRESENCE_CSV} (stage 7, the pan-genome partition) and "
        f"{CORE_ALIGNMENT_ALN} (stage 8, gubbins' input), so one command "
        "unblocks both stages. The GFF3 files come from stage 2 (Bakta) and "
        "need no adapter - Bakta already emits the embedded ##FASTA section "
        "panaroo requires."
    )

#: The wide matrix pyseer reads (`pyseer --pres`, "as produced by roary").
#: panaroo writes it itself; we keep it rather than re-deriving one.
PRESENCE_RTAB = "gene_presence_absence.Rtab"

#: Core/accessory are our partition, not panaroo's. panaroo emits a file per
#: gene family and this is the concatenation, which panaroo's own docs describe
#: as its notion of the core. It is emitted by the tool and left on disk; it is
#: NOT read by any code path here, so nothing in this pipeline cross-checks
#: anything against it. `parse_summary_statistics` is the only reader of
#: panaroo's own numbers, and it is a verification helper rather than part of
#: the run - see its docstring for why it is not wired in.
PAN_GENOME_REFERENCE = "pan_genome_reference.fa"

#: panaroo's own tally of the gene families, one per line, tab-separated:
#: ``<label>\t<range>\t<count>``. Labels are matched in full, never by prefix -
#: ``Core genes`` and ``Soft core genes`` are distinct buckets and a prefix
#: match reads whichever it hits first.
SUMMARY_STATISTICS = "summary_statistics.txt"

#: label -> attribute name. Exact labels, copied from
#: ``panaroo/generate_output.py`` ``generate_summary_stats`` (panaroo 1.8.0).
SUMMARY_LABELS: Dict[str, str] = {
    "Core genes": "core",
    "Soft core genes": "soft_core",
    "Shell genes": "shell",
    "Cloud genes": "cloud",
    "Total genes": "total",
}

#: Verified against 10 real Bakta GFF3 on 2026-10-02; see the module docstring.
CLEAN_MODE = "moderate"

#: Where panaroo's output is read from, relative to the run's intermediate root.
#:
#: Layout, not environment: the root itself comes from configuration
#: (`PipelineConfig.intermediate_root`), exactly as `amr` and `bakta` take theirs.
#: What is fixed here is the *subdirectory name*, and it is declared once so the
#: runner that eventually invokes panaroo and the reader that consumes its output
#: cannot disagree. Two independently-reasonable names would produce a stage that
#: refuses for want of a file written moments earlier.
OUTPUT_DIRNAME = "panaroo"


def output_dir(intermediate_root: Path) -> Path:
    """Where this run's panaroo output lives.

    Derived from the run's intermediate root, never from a fixed location: the
    same code has to serve TEST and REAL, two results trees and any number of
    worktrees, and a constant here would silently read a previous run's cohort.
    """
    return Path(intermediate_root) / OUTPUT_DIRNAME


def preflight() -> None:
    """Refuse by name if panaroo or cd-hit is missing.

    Raises ToolNotAvailableError naming the specific missing tool and both
    remedies, rather than letting a subprocess die later with something less
    actionable.
    """
    missing: List[str] = []

    try:
        import panaroo  # noqa: F401
    except Exception:
        missing.append(
            "panaroo (python package not importable - on this platform it is "
            "installed from source via pip; see environment/environment.yml)"
        )

    if shutil.which("cd-hit") is None:
        missing.append(
            "cd-hit (not on PATH - panaroo shells out to it; conda package "
            "`cd-hit`, pinned in environment/environment.yml)"
        )

    if missing:
        raise ToolNotAvailableError(
            "Stage 7 (pangenome) cannot verify the toolchain: "
            + "; ".join(missing)
            + ". Both are required: panaroo is importable but useless without "
            "cd-hit, so neither check is redundant."
            "\n\nNote what this refusal does and does not imply. This pipeline "
            "NEVER runs panaroo - there is no subprocess call in this module - "
            "so the preflight is a check that panaroo could be run *here*, not "
            "a step in this run. A cohort whose output was produced elsewhere is "
            f"fine. {provision_instructions()} Only fix the toolchain if you "
            "intend to produce the output on this machine."
        )


def _presence_columns(header: Sequence[str]) -> List[int]:
    """Indices of the per-genome columns in panaroo's CSV.

    Located by name, not assumed. The metadata prefix is Gene / Non-unique Gene
    name / Annotation, and panaroo is free to reorder or add to it.
    """
    metadata = {"gene", "non-unique gene name", "annotation"}
    cols = [
        i
        for i, name in enumerate(header)
        if name.strip().lower() not in metadata
    ]
    if not cols:
        raise DataContractError(
            f"{PRESENCE_CSV} has no per-genome columns; header was "
            f"{list(header)!r}. Expected the three metadata columns "
            f"(Gene, Non-unique Gene name, Annotation) followed by one column "
            f"per genome."
        )
    return cols


def parse_presence_csv(path: Path) -> Tuple[Tuple[str, ...], Dict[str, Set[str]]]:
    """Read panaroo's presence/absence CSV into sample_ids and a presence map.

    Uses `csv.reader` rather than `split(",")`. The Annotation column contains
    commas inside quotes; a naive split misaligns every field after the first
    quoted gene and silently produces wrong per-sample counts. That failure is
    invisible - the row count stays right - so it is worth naming.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    reader = csv.reader(io.StringIO(text))
    rows = [r for r in reader if r and any(c.strip() for c in r)]
    if not rows:
        raise DataContractError(f"{path} is empty; panaroo produced no table")

    header = [c.strip() for c in rows[0]]
    cols = _presence_columns(header)
    sample_ids = tuple(header[i] for i in cols)

    presence: Dict[str, Set[str]] = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        gene = row[0].strip()
        carriers = presence.setdefault(gene, set())
        for i in cols:
            if i < len(row) and row[i].strip() not in ("", "0"):
                carriers.add(header[i])

    return sample_ids, presence


@dataclass(frozen=True)
class PanarooBuckets:
    """panaroo's own gene-family tally, as written to ``summary_statistics.txt``.

    The five counts are kept separately and are NOT collapsed into
    core/non-core, because panaroo's four buckets are not a partition this
    pipeline reproduces: its ``Core`` is *>= 99%* of strains, its ``Shell`` is
    *>= 15% and < 95%* (so it excludes soft-core), and its ``Cloud`` is
    *< 15%* - not "genome-specific". Verified against 10 real Bakta GFF3 on
    2026-10-03: 4,828 core / 0 soft-core / 2,850 shell / 2,341 cloud / 10,019
    total.
    """

    core: int
    soft_core: int
    shell: int
    cloud: int
    total: int

    @property
    def non_core(self) -> int:
        """Everything that is not panaroo-core: soft-core + shell + cloud.

        4,828 core of 10,019 gives 5,191 here, of which 2,850 are shell (present
        in 2-9 of 10 genomes) and 2,341 are cloud (present in 1). Calling all
        5,191 "accessory (2-9 genomes)" conflates two buckets; `shell` alone is
        the 2-9 figure.
        """
        return self.soft_core + self.shell + self.cloud

    def as_dict(self) -> Dict[str, int]:
        return {
            "core": self.core,
            "soft_core": self.soft_core,
            "shell": self.shell,
            "cloud": self.cloud,
            "total": self.total,
        }


def parse_summary_statistics(path: Path) -> PanarooBuckets:
    """Read panaroo's ``summary_statistics.txt`` into its five separate counts.

    A verification helper. It is deliberately NOT called from the run path:
    panaroo's ``Core`` is >= 99% of strains while this pipeline's core is
    exactly 100% (`stages/pangenome.CORE_PRESENCE_THRESHOLD`), so at a cohort
    size like n=967 the two are *expected* to disagree and raising on the
    disagreement would break a correct run. Reconciling that is a decision
    about which definition of core is authoritative, and it has not been made.

    Label matching is by full equality on the first tab-separated field. The
    buckets ``Core genes`` and ``Soft core genes`` overlap as prefixes, so any
    ``startswith``/``in`` match silently reads one bucket's count as the
    other's - and at n=10 soft-core is structurally always 0 (proportion_present
    can only be a multiple of 10, and no multiple of 10 lies in [95, 99)), which
    is exactly the case where a conflation would go unnoticed.
    """
    path = Path(path)
    if not path.exists():
        raise DataContractError(
            f"panaroo did not write {SUMMARY_STATISTICS} at {path}. This file is "
            f"written on every successful run; its absence means the tool failed, "
            f"not that the cohort had no genes."
        )

    found: Dict[str, int] = {}
    unrecognised: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        label = fields[0].strip()
        attr = SUMMARY_LABELS.get(label)
        if attr is None:
            unrecognised.append(label)
            continue
        if len(fields) < 3:
            raise DataContractError(
                f"{SUMMARY_STATISTICS} line {label!r} has {len(fields)} "
                f"tab-separated fields, expected 3 "
                f"(<label>\\t<range>\\t<count>). Line was {line!r}."
            )
        found[attr] = int(fields[2].strip())

    missing = sorted(set(SUMMARY_LABELS.values()) - set(found))
    if missing:
        raise DataContractError(
            f"{SUMMARY_STATISTICS} at {path} is missing the {missing} count(s). "
            f"Recognised labels were {sorted(found)}; unrecognised lines were "
            f"{unrecognised}. A file without all five buckets is not panaroo's "
            f"summary, and reading it as one would report a partial cohort as a "
            f"complete one."
        )

    buckets = PanarooBuckets(**found)
    if buckets.non_core + buckets.core != buckets.total:
        raise DataContractError(
            f"{SUMMARY_STATISTICS} at {path} does not add up: core {buckets.core}"
            f" + soft-core {buckets.soft_core} + shell {buckets.shell} + cloud "
            f"{buckets.cloud} = {buckets.core + buckets.non_core}, but total is "
            f"{buckets.total}. A summary that does not partition its own total is "
            f"truncated or was written by something other than panaroo."
        )
    return buckets


def build_pangenome_from_panaroo(
    panaroo_dir: Path, sample_ids: Sequence[str]
) -> Tuple[Tuple[str, ...], Dict[str, Set[str]]]:
    """Parse panaroo's output, refusing any sample the cohort does not contain.

    A missing isolate is a hard failure, never a silent drop: a pangenome built
    over 9 of 10 isolates reports different core/accessory boundaries, and the
    summary numbers would look plausible while describing the wrong cohort.
    """
    panaroo_dir = Path(panaroo_dir)
    csv_path = panaroo_dir / PRESENCE_CSV
    if not csv_path.exists():
        raise DataContractError(
            f"panaroo did not write {PRESENCE_CSV} in {panaroo_dir}. Expected "
            "the file to exist after a successful run; a missing table means "
            "the tool failed, not that the cohort was empty. "
            + provision_instructions(panaroo_dir)
        )

    found, presence = parse_presence_csv(csv_path)

    expected = set(sample_ids)
    got = set(found)
    if got != expected:
        raise DataContractError(
f"panaroo's presence table does not match the cohort. "
            f"missing={sorted(expected - got)} "
            f"unexpected={sorted(got - expected)}. "
            f"Sample-ID mismatches fail loudly; they are never dropped, because "
            f"a pangenome over a different cohort yields different core "
            f"boundaries while still looking reasonable. The usual cause is a "
            f"table left on disk by an earlier run over a different manifest. "
            f"{provision_instructions(panaroo_dir)}"
        )

    return found, presence
