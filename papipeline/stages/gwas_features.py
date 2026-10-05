"""Producer for ``gwas/gwas_features.tsv`` - stage 12's binary feature matrix.

**Why this module exists.** ``papipeline.stages.gwas.run`` *reads* stage 12's
feature table (``gwas.py:1743``) and nothing in the pipeline wrote it, so the
table was operator-provisioned by default and undocumented as such. The seam
registry recorded it as ``PRODUCED`` by ``papipeline.stages.gwas``, which was
false: that module's only occurrence of the name is the read. This module is
the writer. It does not touch ``gwas.py``; that module's signature, its REAL
refusal and its engine selection are all left exactly as they were.

**The contract** (``docs/data_contract.md``, "Stage input tables"): one
``sample_id`` column, unique, then one column per feature named
``<type>__<label>``, values binary, and every ``<type>`` prefix present in
``config.gwas.feature_types``. Each prefix below is checked against that list
before anything is written, and an undeclared prefix is a refusal rather than a
warning - the contract says "must", and ``gwas.build_input`` only warns and
tests the feature anyway.

**Three families over the declared prefixes.** The prefix names the observation
being made, which is why they are not interchangeable:

``gene__<locus>``
    Whole-gene presence/absence from stage 7's pan-genome matrix. The unit is a
    gene: does this isolate carry this ortholog group?

``gene_presence_absence__<gene>_<variant>``
    One point-mutated allele from stage 4. The unit is an ALLELE, and this is
    the distinction that matters. ``oprD_V359L`` is a missense substitution in
    a gene every isolate in the cohort carries, so reading it as a gene gain or
    a gene loss is a mechanism claim that is false; folding it into a
    ``gene__<gene>``-shaped feature destroys the observation entirely. The label
    therefore carries the variant, the feature reads 1 only for the isolate
    carrying *that substitution*, and the gene's own presence is a separate
    column read from the pan-genome.

    These rows are deliberately not ``snp__``. AMRFinderPlus reports ``POINT``
    (a listed single substitution) and ``POINT_DISRUPT`` (an alignment-detected
    disruption the database does not list), and the disruption rows include
    multi-residue substitutions and indels - ``rpsJ_HKYK56del``,
    ``porB1b_GA120del`` - which are not single-nucleotide variants. Labelling
    them ``snp`` would assert a base-level resolution the evidence does not
    have, and the ``variant`` column holds protein-level allele strings, not
    coordinates.

``gene_presence_absence__<gene>``
    An intact AMR determinant's gene presence, from the same stage 4 table. Also
    a gene, so it shares the allele family's prefix; the two differ only in the
    label, and a label never carries an underscore-joined variant unless the
    table said that row was a point mutation.

Plus ``gene__oprD_absent`` and ``gene__oprD_LoF``, built by
:func:`regulators.oprd_feature_rows` and not re-derived here. That function
already exists, already owns the definition of both verdicts, and had zero
callers; calling it is what puts a producer under it. The prefix is passed in
rather than hard-coded because ``oprd_feature_columns`` documents that choice as
a schema decision belonging to the caller, and ``gene`` is the prefix the
committed stage 12 fixture uses for exactly these two columns.

**What is deliberately absent.**

* ``snp__`` features. ``config.gwas.feature_types`` declares the prefix and the
  committed fixture carries two example columns, but the contract specifies no
  variant features, ``cohort_variants``' header is PROVISIONAL with no
  consumer (``docs/design/cohort-variant-merge.md``), and ticket 16 records the
  ``--pres`` versus ``--vcf`` route as unresolved. A producer for a provisional
  table with no real instance on disk would be a guess wearing a schema.
* ``unitig__`` and ``kmer__`` features. Declared and exemplified in the fixture,
  with no source table and no producer.
* Any statistic. This module counts; it does not test. Nothing here is a
  finding, and every count it logs about a cohort carries the caution below.

**Counts on a small cohort are not findings.** Every per-cohort summary this
module emits carries ``n=<samples>; underpowered, not a finding`` explicitly.
That is a reporting requirement, not decoration: a gene carried by 2 of 10
isolates is a count, and reading it as evidence of anything is the error the
caution exists to prevent.

**Determinism.** The same inputs produce a byte-identical file. Samples are
emitted in sorted order, feature columns in a fixed ``(family rank, label)``
order, and no emitted value depends on dict or filesystem iteration order.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError, SampleIdError
from ..io.tsv import MISSING_SENTINELS, read_tsv, write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from . import regulators as reg

LOGGER = get_logger("stages.gwas_features")

# --------------------------------------------------------------------------
# Paths. One declaration each, so the reader and the writer cannot disagree.
# --------------------------------------------------------------------------

#: This stage's table, relative to the root :func:`produce_gwas_features` is
#: handed. That root is ``config.tool_output_root(mode)``, which is the same
#: directory ``gwas.run`` receives as its ``intermediate_root`` - so this path
#: is exactly the one ``gwas.py:1743`` reads.
FEATURE_TABLE_RELPATH: Tuple[str, ...] = ("gwas", "gwas_features.tsv")

#: Stage 7's gene presence/absence matrix. Stage 7 writes its four outputs into
#: the run's stage directory, so this one sits under ``stages/`` while the tool
#: input tables sit one level up.
PANGENOME_GPA_RELPATH: Tuple[str, ...] = ("stages", "gene_presence_absence.tsv")

#: Stage 4's determinant table - the *input* table the contract declares, which
#: carries the ``#`` provenance header naming the organism. Deliberately not
#: ``stages/04_amr.tsv``: that is the stage's principal output projection, and a
#: consumer that needs the point-mutation split needs the table stage 4 wrote.
AMR_TABLE_RELPATH: Tuple[str, ...] = ("amr", "amr_determinants.tsv")

#: Stage 6's variant table, spelled by the stage that reads and writes it.
REGULATOR_TABLE_RELPATH: Tuple[str, ...] = reg.REGULATOR_TABLE_RELPATH

#: The prefix passed to :func:`regulators.oprd_feature_rows`. Matches the two
#: oprD columns in the committed stage 12 fixture.
OPRD_FEATURE_TYPE: str = "gene"

#: Prefix for the whole-gene family (stage 7's matrix).
GENE_FEATURE_TYPE: str = "gene"

#: Prefix for the stage 4 family, both the allele rows and the intact
#: determinant rows.
AMR_FEATURE_TYPE: str = "gene_presence_absence"

#: Rank fixing the order families appear in, so the column order does not depend
#: on the prefixes being sorted. ``gene`` and ``gene_presence_absence`` are
#: near-synonyms as words, and an alphabetical order would let the shape of the
#: output depend on that accident.
FAMILY_RANK: Dict[str, int] = {
    GENE_FEATURE_TYPE: 0,
    AMR_FEATURE_TYPE: 1,
    OPRD_FEATURE_TYPE: 2,
}

#: Columns required of the AMR table beyond the contract's own required set.
#: ``gene`` and ``variant`` are what distinguish an allele from a whole gene, so
#: a table lacking either cannot produce this table and is refused by name.
AMR_REQUIRED: Tuple[str, ...] = ("sample_id", "antibiotic", "determinant", "gene", "variant")

#: Columns the feature table itself has. ``gwas.FEATURE_COLUMNS`` is the same
#: one-tuple; it is repeated here so this module does not import the stage whose
#: REAL refusal is not its business.
FEATURE_COLUMNS: Tuple[str, ...] = ("sample_id",)


@dataclass(frozen=True)
class InputSpec:
    """One required input, and the stage that produces it.

    Carrying the producer stage is what lets a refusal name all three things a
    reader needs: which file, where it was expected, and who was supposed to
    write it. A refusal naming only the path sends the reader to the filesystem;
    one naming the stage sends them to the code.
    """

    key: str
    relpath: Tuple[str, ...]
    producer_stage: str
    why: str

    def path(self, root: Path) -> Path:
        return Path(root).joinpath(*self.relpath)

    def missing_reason(self, root: Path) -> str:
        return (
            f"{self.key}: expected {self.path(root)} "
            f"(relative path {'/'.join(self.relpath)} under {root}), written by "
            f"stage {self.producer_stage}. {self.why}"
        )


#: Every input, in the order the emitted columns depend on them. All three are
#: required: this table has no partial form, because a feature matrix missing a
#: whole family reads as a cohort with no features of that kind rather than as a
#: run that could not read one input.
REQUIRED_INPUTS: Tuple[InputSpec, ...] = (
    InputSpec(
        key="pangenome_gene_presence_absence",
        relpath=PANGENOME_GPA_RELPATH,
        producer_stage="7 (pangenome)",
        why=(
            "Whole-gene features are read from the pan-genome matrix this "
            "stage writes; there is no other source for a gene-level feature."
        ),
    ),
    InputSpec(
        key="amr_determinants",
        relpath=AMR_TABLE_RELPATH,
        producer_stage="4 (amr)",
        why=(
            "Determinant and point-mutation allele features are read from this "
            "table. It also carries the provenance header naming the organism: "
            "without a valid organism AMRFinderPlus reports no POINT rows at "
            "all, so a table with no variants cannot be told apart from a "
            "cohort never screened for them."
        ),
    ),
    InputSpec(
        key="regulator_variants",
        relpath=REGULATOR_TABLE_RELPATH,
        producer_stage="6 (regulators)",
        why=(
            f"{GENE_FEATURE_TYPE}__{reg.OPRD_ABSENT_FEATURE} and "
            f"{GENE_FEATURE_TYPE}__{reg.OPRD_LOF_FEATURE} are read from this "
            "table through regulators.oprd_feature_rows; both verdicts are the "
            "regulator screen's, not this module's to re-derive."
        ),
    ),
)


def feature_table_path(root: Path) -> Path:
    """Where this table lives under ``root``.

    Args:
        root: The stage-input root (``config.tool_output_root(mode)``), which is
            the directory ``gwas.run`` is handed as its ``intermediate_root``.

    Returns:
        ``<root>/gwas/gwas_features.tsv``.
    """
    return Path(root).joinpath(*FEATURE_TABLE_RELPATH)


def input_paths(root: Path) -> Dict[str, Path]:
    """Every input path, keyed, for a caller that has to wire them up."""
    return {spec.key: spec.path(root) for spec in REQUIRED_INPUTS}


def _power_caution(n_samples: int) -> str:
    """The flag every per-cohort count in this module carries.

    Not decoration. A feature carried by 2 of 10 isolates is a count; reading it
    as evidence of anything is the error this exists to prevent, and a reader who
    sees the number without the flag cannot tell that it was not a test.
    """
    return (
        f"n={n_samples} isolates; underpowered, not a finding - counts only, "
        "no association was tested"
    )


def _refuse_missing_inputs(root: Path, missing: Sequence[InputSpec]) -> None:
    """Refuse naming EVERY missing input, not the first one found.

    All of them go in one message, each with its path and its producer stage.
    Refusing on the first would make an operator fix one file per run to
    discover a three-file problem, which presents as a task far larger than it
    is.
    """
    raise DataContractError(
        "Cannot build the stage 12 feature table "
        f"({'/'.join(FEATURE_TABLE_RELPATH)}): {len(missing)} of "
        f"{len(REQUIRED_INPUTS)} required inputs are absent. Each is listed "
        "with its path and the stage that produces it; all of them must exist "
        "before this table can be written.\n  - "
        + "\n  - ".join(spec.missing_reason(root) for spec in missing),
        root=str(root),
        missing=[spec.key for spec in missing],
        n_missing=len(missing),
    )


def _check_cohort(
    expected: Sequence[str],
    found: Sequence[str],
    path: Path,
    source: str,
    *,
    require_every_sample: bool,
) -> None:
    """Refuse a sample-ID mismatch, naming the samples.

    Rule 5. Which mismatches matter depends on the table, and treating them the
    same is how a correct table gets called broken:

    * ``require_every_sample=True`` for the pan-genome matrix, which is dense by
      construction - one row per (gene, sample) pair - so a missing sample is a
      truncated or foreign table. A pangenome over 9 of 10 isolates reports
      different core/accessory boundaries while looking entirely reasonable.
    * ``require_every_sample=False`` for stage 4's determinants and stage 6's
      variants, which are presence-dependent: a sample with no determinant
      produces no row, because absence is a finding rather than a missing
      record (``contracts.DENSE_PER_SAMPLE_STAGES`` says exactly this). Those
      samples are filled 0, and only a sample the manifest does not contain is a
      mismatch - that one is a table from a different cohort.
    """
    want = set(expected)
    got = set(found)
    missing = sorted(want - got)
    unexpected = sorted(got - want)
    if require_every_sample:
        problem = bool(missing) or bool(unexpected)
    else:
        problem = bool(unexpected)
    if not problem:
        return
    detail = (
        f"missing_from_table={missing} unexpected_in_table={unexpected}. "
        if require_every_sample
        else f"unexpected_in_table={unexpected}. "
        "This table is presence-dependent, so a manifest sample with no row is "
        "read as carrying none of its features, which is a finding and not a "
        "missing record. "
    )
    raise SampleIdError(
        f"{source} does not describe this run's cohort. {detail}"
        "Sample-ID mismatches fail loudly and are never dropped: a feature "
        "matrix built over a different cohort is wrong in ways no later stage "
        "can detect. The usual cause is a table left on disk by an earlier run "
        "over a different manifest.",
        path=str(path),
        source=source,
        missing=missing,
        unexpected=unexpected,
    )


def gene_features(presence: Mapping[str, Sequence[str]]) -> Dict[str, Set[str]]:
    """Whole-gene features: one ``gene__<locus>`` per gene in the matrix.

    Args:
        presence: ``{gene: carrier_sample_ids}``, as returned by
            :func:`pangenome.read_gene_presence_absence`. A gene carried by no
            sample is still a column, reading 0 throughout; dropping it would
            make an absent gene indistinguishable from one nobody looked for.

    Returns:
        ``{column_name: carriers}``. The column name is built here rather than at
        the point of use so every family hands the assembly the same shape.
    """
    return {
        f"{GENE_FEATURE_TYPE}__{gene}": set(presence[gene]) for gene in sorted(presence)
    }


def target_antibiotic(config: PipelineConfig) -> str:
    """The antibiotic under analysis, and whether it was chosen or defaulted.

    ``config.antibiotics[0]`` is a position rather than a choice, so it names
    whichever drug happens to be listed first. The project's own
    ``primary_antibiotic`` is the choice. This mirrors the resolution in
    ``gwas.build_input`` - fallback included - so the feature table and the
    analysis cannot disagree about which drug is under study, and a fallback
    that actually fires is logged, because a table filtered on the wrong drug
    comes out empty rather than looking wrong.
    """
    project = (getattr(config, "raw", None) or {}).get("project", {}) or {}
    configured = project.get("primary_antibiotic")
    if configured:
        return str(configured)
    antibiotic = str(config.antibiotics[0])
    LOGGER.warning(
        "config.project.primary_antibiotic is unset; AMR features are "
        "filtered on config.antibiotics[0] = %r, which is a position in a list "
        "rather than a choice. AMR features are empty for any antibiotic the "
        "table attributes elsewhere.",
        antibiotic,
    )
    return antibiotic


def amr_features(
    rows: Sequence[Mapping[str, object]],
    *,
    antibiotic: str,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], int]:
    """AMR features: point-mutation allele presence, and determinant presence.

    Args:
        rows: The stage 4 determinant table, already read.
        antibiotic: The antibiotic under analysis. Only determinants the table
            attributes to it become features - a ``blaIMP`` determinant is
            evidence about imipenem, and a feature matrix carrying determinants
            for other drugs would test this cohort against mechanisms the study
            never claimed.

    Returns:
        ``(allele_carriers, determinant_carriers, n_dropped_other_antibiotic)``.
        The two mappings are ``{feature_label: {sample_id}}`` - the same shape
        stage 7's matrix is read in, so the assembly treats all three families
        identically and cannot handle one of them differently by accident. They
        are kept apart so a point mutation cannot reach the table through the
        gene-presence route: the mis-modelling a prior triage found, and the
        reason this is two collections and not one.

    Raises:
        DataContractError: A row carries a ``variant`` with no ``gene``. Stage 4
            splits ``<gene>_<variant>`` as a pair, so a variant without its gene
            means the table was written by something that did not. Recovering the
            gene from the determinant would be a mechanism claim attributed to
            the wrong gene, which is worse than an absent feature.
    """
    allele_carriers: Dict[str, Set[str]] = {}
    determinant_carriers: Dict[str, Set[str]] = {}
    dropped = 0
    unsplit: List[str] = []

    for row in rows:
        if _cell(row, "antibiotic") != antibiotic:
            dropped += 1
            continue
        determinant = _cell(row, "determinant")
        gene = _cell(row, "gene")
        variant = _cell(row, "variant")
        if not determinant:
            continue
        sample_id = _cell(row, "sample_id")
        if variant:
            if not gene:
                unsplit.append(f"{sample_id}:{determinant}")
                continue
            # A point mutation: the unit is the ALLELE, so the label is the gene
            # AND the substitution and the carriers are the isolates read as
            # carrying *that* substitution.
            allele_carriers.setdefault(f"{gene}_{variant}", set()).add(sample_id)
        else:
            # An intact determinant: the unit is the GENE, so the label is the
            # gene alone and no variant ever appears in it.
            determinant_carriers.setdefault(gene or determinant, set()).add(sample_id)

    if unsplit:
        raise DataContractError(
            f"{len(unsplit)} AMR row(s) carry a `variant` with no `gene`, so "
            "the point mutation cannot be named as an allele and must not be "
            "folded into the gene-presence features: a substitution read as a "
            "gene gain is a false mechanism claim. Stage 4 splits "
            "`<gene>_<variant>` as a pair, so a table that does not came from "
            "somewhere other than this pipeline's own stage 4. Offending rows: "
            + ", ".join(sorted(unsplit)[:10]),
            n_unsplit=len(unsplit),
            examples=sorted(unsplit)[:10],
        )

    return allele_carriers, determinant_carriers, dropped


def _cell(row: Mapping[str, object], column: str) -> str:
    """One table cell as text, with the TSV missing-sentinels read as absent.

    :func:`read_tsv` already maps ``.``/``NA``/``-`` to ``None``, so this is a
    no-op for a table read from disk. It is not redundant: `build_feature_rows`
    is public and is called directly with constructed rows in tests, and a
    ``variant`` of ``"."`` read as a point mutation would put an allele column
    in the table for a gene the row says is intact - the exact mis-modelling
    this module exists to prevent, arriving through the wrong door.
    """
    value = row.get(column)
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in MISSING_SENTINELS else text


def _order_columns(columns: Mapping[str, str]) -> List[str]:
    """Total, reproducible column order: family rank first, then label.

    ``columns`` maps a full column name to its declared prefix. Sorting by
    ``(rank, name)`` is total because the names are unique, so the emitted order
    cannot depend on insertion order.
    """
    return sorted(columns, key=lambda name: (FAMILY_RANK.get(columns[name], 99), name))


def _refuse_undeclared_prefixes(
    columns: Mapping[str, str], allowed: Sequence[str]
) -> None:
    """Refuse any column whose ``<type>`` prefix is not declared in config.

    The contract says the prefix "must appear in
    ``config.gwas.feature_types``". ``gwas.build_input`` warns and tests the
    feature anyway, so a producer that emitted an undeclared prefix would ship a
    contract violation downstream in a shape the consumer had already decided to
    accept. Refusing here keeps the table conformant by construction.
    """
    undeclared = sorted(
        {prefix for prefix in columns.values() if prefix not in set(allowed)}
    )
    if not undeclared:
        return
    examples = sorted(
        name for name, prefix in columns.items() if prefix in set(undeclared)
    )[:5]
    raise DataContractError(
        f"Feature prefix(es) {undeclared} are not in "
        f"config.gwas.feature_types {sorted(allowed)}, which "
        "docs/data_contract.md requires every <type> prefix to appear in. "
        f"Examples: {examples}. Declare the prefix in config or drop the "
        "family; emitting it anyway would put a contract violation into a table "
        "whose consumer only warns.",
        undeclared=undeclared,
        allowed=sorted(allowed),
        examples=examples,
    )


def build_feature_rows(
    config: PipelineConfig,
    manifest: SampleManifest,
    *,
    gene_presence: Mapping[str, Sequence[str]],
    amr_rows: Sequence[Mapping[str, object]],
    regulator_records: Sequence[reg.RegulatorVariant],
) -> Tuple[List[Dict[str, object]], List[str], Dict[str, int]]:
    """Assemble the wide feature table from the three sources.

    Split out from :func:`produce_gwas_features` so the assembly is testable
    against synthetic inputs with no directory tree on disk, while the reading
    and the refusing stay together in one place.

    Args:
        config: The loaded pipeline configuration. Used for the declared
            feature prefixes and the target antibiotic, both of which are
            scientific choices that live in config (rule 2).
        manifest: The cohort. One output row per manifest sample.
        gene_presence: ``{gene: carrier_sample_ids}`` from stage 7's matrix.
            Must not be empty; a zero-gene pan-genome is refused. Its shape is
            the one :func:`amr_features` returns too.
        amr_rows: Stage 4's determinant table, as read.
        regulator_records: Stage 6's parsed variant records.

    Returns:
        ``(rows, columns, counts)``. ``rows`` is one dict per manifest sample in
        sorted order, ``columns`` is every feature column in the fixed order
        from :func:`_order_columns`, and ``counts`` holds the per-family feature
        totals for the log and the report.

    Raises:
        DataContractError: Stage 7's matrix contributed no genes, an AMR row
            cannot be split into gene and allele, a feature prefix is
            undeclared, or two features claim one column name.
    """
    sample_ids = sorted(manifest.sample_ids)
    if not sample_ids:
        raise DataContractError(
            "Cannot build the stage 12 feature table: the manifest names no "
            "sample, so there is no row to write.",
            n_samples=0,
        )

    if not gene_presence:
        raise DataContractError(
            "Stage 7's gene presence/absence matrix contributed no genes, so "
            "there is no whole-gene feature to write. A pan-genome with zero "
            "genes is a failed input, not a cohort with no genes: "
            "docs/scientific_rules.md section 9 requires an empty input to be "
            "refused rather than reported. The matrix is written by stage 7 "
            "(pangenome), which refuses a zero-gene partition itself; reaching "
            "here means the table on disk is not the one that stage wrote.",
            n_samples=len(sample_ids),
            input="stage 7's gene presence/absence matrix",
        )

    antibiotic = target_antibiotic(config)
    gene_carriers = gene_features(gene_presence)
    allele_carriers, determinant_carriers, dropped = amr_features(
        amr_rows, antibiotic=antibiotic
    )
    if dropped:
        LOGGER.info(
            "AMR features: %d determinant row(s) attributed to an antibiotic "
            "other than %s and excluded from the feature table",
            dropped,
            antibiotic,
        )

    oprd_rows = reg.oprd_feature_rows(
        reg.group_by_sample(regulator_records, manifest), OPRD_FEATURE_TYPE
    )
    oprd_columns = [
        f"{OPRD_FEATURE_TYPE}__{label}" for label in reg.oprd_feature_labels()
    ]

    # Every family contributes to ONE pair of mappings keyed by the full column
    # name - `presence` for carriers, `prefixes` for the declared type - so the
    # row assembly below cannot treat one family differently from another. The
    # prefixes differ, so no two families can collide on a name.
    presence: Dict[str, Set[str]] = {}
    prefixes: Dict[str, str] = {}
    for name, carriers in gene_carriers.items():
        presence[name] = set(carriers)
        prefixes[name] = GENE_FEATURE_TYPE
    for label, carriers in allele_carriers.items():
        name = f"{AMR_FEATURE_TYPE}__{label}"
        presence[name] = set(carriers)
        prefixes[name] = AMR_FEATURE_TYPE
    for label, carriers in determinant_carriers.items():
        name = f"{AMR_FEATURE_TYPE}__{label}"
        presence[name] = set(carriers)
        prefixes[name] = AMR_FEATURE_TYPE
    for name in oprd_columns:
        # The oprD pair is not a carrier set: both verdicts come from the
        # regulator screen's per-sample rows, merged in below.
        presence[name] = set()
        prefixes[name] = OPRD_FEATURE_TYPE

    _refuse_undeclared_prefixes(prefixes, list(config.gwas.feature_types))
    ordered = _order_columns(prefixes)

    oprd_by_sample = {str(row.get("sample_id")): row for row in oprd_rows}
    if len(oprd_by_sample) != len(oprd_rows):
        duplicated = sorted(
            sid for sid in oprd_by_sample
            if sum(1 for r in oprd_rows if str(r.get("sample_id")) == sid) > 1
        )
        raise DataContractError(
            "regulators.oprd_feature_rows returned two rows for one sample, so "
            "the oprD features cannot be written without one row silently "
            "winning. Offending samples: "
            + ", ".join(duplicated[:10]),
            samples=duplicated[:10],
        )

    rows: List[Dict[str, object]] = []
    for sample_id in sample_ids:
        row: Dict[str, object] = {"sample_id": sample_id}
        for name in ordered:
            row[name] = int(sample_id in presence[name])
        oprd = oprd_by_sample.get(sample_id)
        if oprd is None:
            raise DataContractError(
                "regulators.oprd_feature_rows produced no row for a manifest "
                "sample, so its oprD features cannot be written. Silence must "
                "not read as an absent gene.",
                sample_id=sample_id,
            )
        for name in oprd_columns:
            row[name] = int(oprd.get(name) or 0)
        rows.append(row)

    counts = {
        "pangenome_gene": len(gene_carriers),
        "amr_allele": len(allele_carriers),
        "amr_determinant_gene": len(determinant_carriers),
        "oprd": len(oprd_columns),
        "amr_dropped_other_antibiotic": dropped,
    }
    counts["total"] = len(ordered)

    # A column carried by every isolate, or by none, carries no information and
    # cannot be tested; the count is reported so a reader is not left wondering
    # where the constant columns went. They are kept, not dropped: choosing which
    # features to withhold is a scientific decision, and the consumer owns the
    # testability filter that already implements one. Read off the emitted rows
    # rather than off `carriers`, because the oprD pair is never a carrier set.
    # An AMR determinant gene that the pan-genome also knows becomes two
    # columns testing the same presence under two prefixes, which counts twice
    # in any multiple-testing denominator. That is a scientific decision about
    # the test set, not one to take silently, so the overlap is reported. Both
    # columns are kept: dropping either would discard a call some stage made.
    #
    # Only determinant genes can collide. An allele label is `<gene>_<variant>`
    # and so can only equal a pan-genome locus carrying that exact spelling,
    # which is not a locus name; comparing gene names to gene names is also the
    # comparison that means something, so no label is stripped to make it fit.
    redundant = sorted(set(determinant_carriers) & set(gene_presence))
    if redundant:
        LOGGER.info(
            "%d AMR determinant gene(s) are also pan-genome loci and so appear "
            "under two prefixes, e.g. %s. Both columns are written; a "
            "multiple-testing denominator that counts them twice is a decision "
            "for the analysis, not for this producer. %s",
            len(redundant),
            ", ".join(redundant[:5]),
            _power_caution(len(sample_ids)),
        )

    carried_by: Dict[str, int] = {
        name: sum(int(row[name]) for row in rows) for name in ordered
    }
    invariant = sorted(
        name
        for name in ordered
        if carried_by[name] in (0, len(sample_ids))
    )
    if invariant:
        LOGGER.info(
            "%d of %d feature column(s) are carried by every or no isolate and "
            "so cannot be tested; kept, because dropping them here would be an "
            "unstated scientific choice. %s",
            len(invariant),
            len(ordered),
            _power_caution(len(sample_ids)),
        )

    LOGGER.info(
        "Stage 12 feature table: %d samples x %d features | %s",
        len(sample_ids),
        len(ordered),
        _power_caution(len(sample_ids)),
    )
    for label, key in (
        ("pangenome gene", "pangenome_gene"),
        ("AMR point-mutation allele", "amr_allele"),
        ("AMR determinant gene", "amr_determinant_gene"),
        ("oprD", "oprd"),
    ):
        LOGGER.info("  %s features: %d", label, counts[key])
    LOGGER.info("  feature columns total: %d", counts["total"])

    return rows, ["sample_id", *ordered], counts


def _provenance(
    root: Path,
    mode: object,
    manifest: SampleManifest,
    counts: Mapping[str, int],
) -> List[str]:
    """The ``#`` banner above the header, naming every input and its producer.

    ``read_tsv`` skips ``#`` lines, so this costs the consumer nothing. It earns
    its place because a stage 12 result is otherwise unattributable: the file
    does not say which stage 4 run it came from, and a table whose variants are
    all empty has two very different explanations - a cohort screened without a
    valid organism, and a cohort screened with one and found nothing.
    """
    return [
        "Stage 12 feature matrix - binary presence/absence per sample.",
        "One column per feature, named <type>__<label>; every <type> prefix is "
        "in config.gwas.feature_types.",
        f"mode={getattr(mode, 'value', mode)} "
        f"n_samples={len(manifest.sample_ids)} "
        f"n_features={counts.get('total', 0)}",
        _power_caution(len(manifest.sample_ids)),
        f"root={root}",
        f"input: pangenome gene presence/absence <- {REQUIRED_INPUTS[0].path(root)}"
        " (stage 7)",
        f"input: AMR determinants <- {REQUIRED_INPUTS[1].path(root)} (stage 4)",
        f"input: regulator variants <- {REQUIRED_INPUTS[2].path(root)} (stage 6)",
        f"features: pangenome_gene={counts.get('pangenome_gene', 0)} "
        f"amr_point_mutation_allele={counts.get('amr_allele', 0)} "
        f"amr_determinant_gene={counts.get('amr_determinant_gene', 0)} "
        f"oprD={counts.get('oprd', 0)}",
        "amr rows excluded as another antibiotic: "
        f"{counts.get('amr_dropped_other_antibiotic', 0)}",
        "The AMR allele family holds point-mutation ALLELES, not gene gains: a "
        "label carries <gene>_<variant> only where stage 4 populated `variant`.",
        "No statistic is computed here. Nothing in this file is a finding.",
    ]


def produce_gwas_features(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode,
    intermediate_root: Path,
) -> Path:
    """Build ``<intermediate_root>/gwas/gwas_features.tsv`` and return its path.

    The producer half of stage 12's input, split out so that ``gwas.run`` stays
    a plain *reader* of the table while the wiring that feeds it lives at the
    call site. That split is the pattern
    :func:`regulators.produce_regulator_variants` established for stage 6, and
    it is what lets ``papipeline/run.py`` own every path without this module
    knowing where the run's directories are.

    Args:
        config: The loaded pipeline configuration.
        manifest: The cohort. One output row per manifest sample.
        mode: The run's mode. Recorded in the provenance banner; it does not
            branch, because the three sources are the same tables in TEST and
            REAL and a mode-gated producer would be a second, quieter way for
            the file to go missing.
        intermediate_root: The stage-input root - pass
            ``config.tool_output_root(mode)``. This is the same directory
            ``gwas.run`` reads the table from, so the two cannot disagree.

    Returns:
        The path written.

    Raises:
        DataContractError: Any required input is absent (all of them named, with
            path and producer stage), a source describes a different cohort, the
            AMR table cannot be split into gene and allele, a feature prefix is
            undeclared, or stage 7's matrix contributed no genes.
        SampleIdError: A source table's sample IDs disagree with the manifest.
    """
    from .pangenome import read_gene_presence_absence

    root = Path(intermediate_root)
    # Before the input check, so an empty cohort is reported as the empty cohort
    # it is rather than as three absent files.
    if not manifest.sample_ids:
        raise DataContractError(
            "Cannot build the stage 12 feature table: the manifest names no "
            "sample, so there is no row to write. One row per manifest sample "
            "is the contract; a table with no rows would be a header-only file "
            "reading as a cohort with no features.",
            n_samples=0,
            root=str(root),
        )
    missing = [spec for spec in REQUIRED_INPUTS if not spec.path(root).is_file()]
    if missing:
        _refuse_missing_inputs(root, missing)

    gpa_path = REQUIRED_INPUTS[0].path(root)
    amr_path = REQUIRED_INPUTS[1].path(root)
    reg_path = REQUIRED_INPUTS[2].path(root)

    pangenome_samples, presence = read_gene_presence_absence(gpa_path)
    _check_cohort(
        manifest.sample_ids,
        pangenome_samples,
        gpa_path,
        "stage 7's gene matrix",
        require_every_sample=True,
    )

    amr_rows = read_tsv(amr_path, required_columns=AMR_REQUIRED)
    _check_cohort(
        manifest.sample_ids,
        sorted({_cell(row, "sample_id") for row in amr_rows}),
        amr_path,
        "stage 4's determinant table",
        require_every_sample=False,
    )

    regulator_records = reg.load_regulator_variants(config, reg_path, root)
    _check_cohort(
        manifest.sample_ids,
        sorted({record.sample_id for record in regulator_records}),
        reg_path,
        "stage 6's regulator variant table",
        require_every_sample=False,
    )

    rows, columns, counts = build_feature_rows(
        config,
        manifest,
        gene_presence=presence,
        amr_rows=amr_rows,
        regulator_records=regulator_records,
    )
    path = feature_table_path(root)
    write_tsv(
        path,
        rows,
        columns,
        header_comment=_provenance(root, mode, manifest, counts),
    )
    LOGGER.info(
        "Wrote %s (%d features, %d of them point-mutation alleles; %s)",
        path,
        counts["total"],
        counts["amr_allele"],
        _power_caution(len(manifest.sample_ids)),
    )
    return path



__all__ = [
    "AMR_FEATURE_TYPE",
    "AMR_REQUIRED",
    "AMR_TABLE_RELPATH",
    "FEATURE_COLUMNS",
    "FEATURE_TABLE_RELPATH",
    "FAMILY_RANK",
    "GENE_FEATURE_TYPE",
    "InputSpec",
    "OPRD_FEATURE_TYPE",
    "PANGENOME_GPA_RELPATH",
    "REGULATOR_TABLE_RELPATH",
    "REQUIRED_INPUTS",
    "amr_features",
    "build_feature_rows",
    "feature_table_path",
    "gene_features",
    "input_paths",
    "produce_gwas_features",
    "target_antibiotic",
]
