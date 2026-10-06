"""Per-stage output contracts for the real pipeline.

Every stage in ``run.py`` writes one principal table into
``<intermediate>/stages``. The filename and the exact header of that table
are declared once, here, and used by both the code that writes it and the
contract that validates it.

Sharing the declaration is the whole point. A contract naming a different
file, or a column the stage no longer emits, would validate nothing while
appearing to — and that is precisely the class of bug the execution layer
exists to catch. The declarations below were taken from the headers a real
TEST run actually wrote, and
``tests/integration/test_run_observability.py::test_contracts_match_what_the_pipeline_writes``
fails if they ever drift from reality.

Column names are *not* imported from the stage modules here. Several
modules expose several ``*_COLUMNS`` constants, and a stage's principal
table is not always the record schema its module defines — the annotation
stage writes a three-column *summary*, not the ten-column record table. An
imported constant would look tidier and be wrong.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

from .validation import (
    Check,
    CheckKind,
    OutputSpec,
    covers_samples,
    has_columns,
    min_rows,
    non_empty,
    parses,
)

#: stage -> (filename as written, header as written)
STAGE_TABLES: Dict[str, Tuple[str, Tuple[str, ...]]] = {
    "validation": ("01_validation.tsv", (
        "sample_id",         "assembly_size",         "contig_count",         "n50",         "n90",         "largest_contig",         "gc_content",         "ambiguous_bases",         "basic_quality_status",         "species_confirmation",         "contamination_status",         "completeness_status",         "duplicate_status",         "flag_reasons", 
    )),
    "annotation": ("02_annotation_summary.tsv", (
        "sample_id",         "n_records",         "n_named_genes", 
    )),
    "mlst": ("03_mlst.tsv", (
        "sample_id",         "ST",         "alleles",         "MLST_status",         "mlst_scheme",         "allele_database", 
    )),
    "amr": ("04_amr.tsv", (
        "sample_id",         "antibiotic",         "determinant",         "gene",         "variant",         "determinant_type",         "mechanism",         "evidence_source",         "database",         "database_version",         "confidence",         "claim_status",         "identity_pct",         "coverage_pct", 
    )),
    "virulence": ("08_virulence.tsv", (
        "sample_id",         "virulence_factor",         "gene",         "category",         "database",         "database_version",         "confidence",         "identity_pct", 
    )),
    "pangenome": ("pangenome_summary.tsv", (
        "metric",         "value", 
    )),
    "phylogeny": ("10_phylogeny.tsv", (
        "tree",         "n_tips",         "n_internal_nodes",         "has_branch_lengths",         "matches_manifest", 
    )),
    "phenotype": ("11_phenotype.tsv", (
        "sample_id",         "antibiotic",         "phenotype",         "MIC",         "MIC_unit",         "zone_diameter",         "zone_unit",         "source", 
    )),
    "gwas": ("12_gwas.tsv", (
        "feature",         "feature_type",         "effect",         "p_value",         "adjusted_p_value",         "effect_size",         "frequency",         "lineage_distribution",         "model", 
    )),
    "convergence": ("13_convergence.tsv", (
        "determinant",         "independent_lineages",         "branch_count",         "distribution",         "convergence_category", 
    )),
    "cooccurrence": ("14_cooccurrence.tsv", (
        "feature_a",         "feature_b",         "feature_type",         "n_a",         "n_b",         "n_both",         "statistic",         "statistic_value",         "adjusted_p_value",         "interpretation_limit", 
    )),
    "reporting": ("16_figure_manifest.tsv", (
        "figure",         "kind",         "n_items", 
    )),
    # --- spec.md D1 stages with no implementation (ticket 14) -------------
    # Declared so the DAG can list spec.md's fifteen, and so a caller who asks
    # for one is told it is unbuilt rather than handed a plausible empty table.
    #
    # PROVISIONAL: these column sets are placeholders. The stages fabricate a
    # header-only table in STUB and refuse to run outside it, so no parser reads
    # them. Ticket 14 sets the real columns against the real output.
    # Derived from a real `bcftools mpileup | bcftools call -mv` run, not from
    # an expectation of bcftools. Per-isolate calls against the reference; they
    # conflate species-wide fixed differences with cohort polymorphisms, so a
    # cohort merge is still required before these are GWAS features.
    "variants": ("variants.tsv", (
        "sample_id", "chrom", "pos", "ref", "alt", "qual", "filter",
        "GT", "AC", "AN", "DP4", "MQ", "MQ0F",
    )),
    # Cohort-wide merge of the per-isolate calls. Deliberately NOT a per-sample
    # table: one row per variable position across the whole cohort, not one row
    # per isolate. The exact header is still open - docs/design/cohort-variant-
    # merge.md - so these are the columns every plausible contract needs, and
    # they are provisional. No consumer writes a row to them yet.
    # Long/tidy: one row per variable site cohort-wide, not one column per
    # isolate. pyseer is wide only via `--pres`, where samples come from the
    # header row, so the pivot is a separate concern with no consumer yet. The
    # columns now match the implemented contract in stages/cohort_variants.py.
    "cohort_variants": ("cohort_variants.tsv", (
        "chrom", "pos", "ref", "alt", "ac", "an", "an_calls", "af",
    )),
    "recombination": ("recombination.tsv", (
        "node",         "n_snps",         "mean_branch_length",         "recombination_detected", 
    )),
    # A SQUARE MATRIX, not a long pair list. Header row of sample names, one
    # row per sample carrying its distance to every other sample, `sample_id`
    # leading so the file is keyed on the one identifier every stage uses.
    #
    # This replaces a `("sample_a", "sample_b", "distance", "method")`
    # placeholder that nothing read and that duplicated every distance - 7,405
    # rows at 835 isolates against 835. The stand-in fixture written for this
    # stage was already square, so the fixture and the contract had disagreed.
    # The distance is PATRISTIC: summed branch lengths along the stage-9 tree.
    # Method is recorded in the file so a reader cannot assume a different
    # definition produced these numbers.
    "similarity": ("similarity.tsv", (
        "sample_id",
        # One column holding the sample's full distance vector. A square matrix
        # on disk as separate columns would need a header per sample, which
        # changes with the cohort; this form is keyed on the cohort instead.
        "distances",
    )),
}

#: Outputs of the steps that spec.md:351 folds into another stage.
#:
#: Still written, still consumed - the master table most of all - but no longer
#: a *stage's* principal output, so not in STAGE_TABLES. Declared here for the
#: same reason STAGE_TABLES exists: one source for the location, so the code
#: that writes a file and any contract that checks it cannot disagree about
#: where it went.
INTERNAL_TABLES: Dict[str, Tuple[str, Tuple[str, ...]]] = {
    "structural_variants": ("07_structural_variants.tsv", (
        "sample_id",         "variant_id",         "variant_type",         "position",         "affected_gene",         "size",         "evidence",         "confidence",         "call_status",         "mge", 
    )),
    "regulators": ("06_regulators.tsv", (
        "sample_id",         "gene",         "variant",         "variant_type",         "position",         "reference",         "alternate",         "effect",         "mechanism",         "confidence",         "evidence_source",         "call_status", 
    )),
    "mechanisms": ("05_mechanisms.tsv", (
        "sample_id",         "antibiotic",         "determinant",         "mechanism",         "evidence_level",         "gene",         "gene_type",         "notes", 
    )),
    "master_table": ("15_master_table.tsv", (
        "sample_id",         "antibiotic",         "phenotype",         "amr_gene",         "amr_variant",         "chromosomal_mutation",         "regulator",         "mechanism",         "structural_variant",         "mlst",         "lineage",         "virulence_profile",         "gwas_feature",         "gwas_status",         "convergence_status",         "confidence",         "evidence_notes", 
    )),
}


def internal_table_path(stage_dir: Path, name: str) -> Path:
    """Where a folded step's table lives.

    The single source of that location, as :func:`table_path` is for stages.
    Declare the name in INTERNAL_TABLES rather than inlining a path at the call
    site, or the writer and the contract drift apart - which is the class of bug
    this module exists to prevent.
    """
    try:
        filename, _columns = INTERNAL_TABLES[name]
    except KeyError:
        raise KeyError(
            f"No output is declared for folded step {name!r}. Add it to "
            "INTERNAL_TABLES rather than inventing a path here."
        ) from None
    return Path(stage_dir) / filename

#: Stages whose table is *nominally* one row per sample, and for which a
#: caller may therefore opt into a coverage check.
#:
#: This is not a claim that coverage holds. Several of these are
#: presence-dependent in practice - a sample with no antimicrobial
#: determinant, no regulator variant, no structural variant and no virulence
#: factor produces no row at all, because absence is a finding rather than a
#: missing record. Coverage is therefore never asserted by default; see the
#: note in :func:`stage_spec`.
DENSE_PER_SAMPLE_STAGES = frozenset({
    "validation", "annotation", "mlst", "mechanisms", "virulence",
    "phenotype", "integration",
})

#: Kept for callers that want the coarser "has a sample_id column" notion.
PER_SAMPLE_STAGES = DENSE_PER_SAMPLE_STAGES | frozenset({
    "amr", "regulators", "structural_variants", "variants", "similarity",
})
# `variants` was absent from this set even though spec.md D1 lists stage 6 as
# per-sample, so the taxonomy did not know the stage it had just given a
# per-isolate contract. It is here now, and `cohort_variants` deliberately is
# not: that difference is what makes the new stage cohort-wide rather than a
# rename of the old one.
# `similarity` joined on 2026-09-30 when its contract became a square matrix:
# one row per sample, so the per-sample property holds exactly.


def table_path(stage_dir: Path, stage: str) -> Path:
    """Where a stage's principal table lives.

    The single source of that location. ``run.py`` writes here and the
    contract validates here, so the two cannot disagree.
    """
    try:
        filename, _ = STAGE_TABLES[stage]
    except KeyError:
        raise KeyError(
            f"No output contract is declared for stage {stage!r}. Add it to "
            f"STAGE_TABLES rather than inventing a path here."
        ) from None
    return Path(stage_dir) / filename


def required_columns(stage: str) -> Tuple[str, ...]:
    """The header this stage's principal table is required to carry."""
    entry = STAGE_TABLES.get(stage)
    return entry[1] if entry else ()


def stage_spec(
    stage: str,
    table: Optional[Path] = None,
    *,
    stage_dir: Optional[Path] = None,
    n_samples: Optional[int] = None,
    extra_paths: Sequence[Path] = (),
    description: str = "",
    require_coverage: bool = False,
) -> OutputSpec:
    """The contract for one stage's principal table.

    The checks are about the artefact, never about biology: it exists, it
    is not empty, it parses as TSV, it carries the header the stage
    declares, and it holds at least one record. None of them can be
    satisfied by a file that merely exists.
    """
    if table is None:
        if stage_dir is None:
            raise ValueError("pass either table or stage_dir")
        table = table_path(stage_dir, stage)
    columns = required_columns(stage)
    checks = [
        non_empty(table, description=description or f"{stage}: output is non-empty"),
        parses(table, "tsv", description=f"{stage}: output is a readable TSV"),
    ]
    if columns:
        checks.append(
            has_columns(table, columns,
                        description=f"{stage}: output carries its declared header")
        )
    checks.append(min_rows(table, 1, description=f"{stage}: output has at least one record"))
    if (
        n_samples is not None
        and stage in DENSE_PER_SAMPLE_STAGES
        and require_coverage
    ):
        # Opt-in only, and off by default. Coverage is a statement about the
        # *data*: whether every isolate has a determinant call depends on
        # what was found, not on whether the stage ran. Asserting it as a
        # default execution contract would fail a correct run whose correct
        # answer is "this isolate has no AMR gene", so a caller that really
        # wants it has to ask for it.
        checks.append(
            covers_samples(table, n_samples, key="sample_id",
                           description=f"{stage}: output accounts for all {n_samples} samples")
        )
    for path in extra_paths:
        checks.append(non_empty(path, description=f"{stage}: companion artefact is non-empty"))
    return OutputSpec(stage=stage, checks=tuple(checks))
