"""Pilot-100 QC report.

The report leads with the selection method, because that is the fact a reader
most needs to be able to verify. It states the method as implemented, names
the alternative it is *not*, and reports how much the two differ.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from . import mechanisms as mech_mod

from papipeline.pilot.census import describe_reasons

LOGGER = get_logger("pilot.report")

REPORT_NAME = "PILOT100_QC_REPORT.md"

SELECTION_METHOD = (
    "FIRST 100 ASSEMBLED GENOMES PRESENT IN data/, sorted by GCA accession."
)
SELECTION_METHOD_NOT = (
    "NOT the first 100 rows of PDC_essential.tsv."
)


def _table(columns: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    if not rows:
        return "_No rows._\n"
    lines = [
        "| " + " | ".join(str(c) for c in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join("-" if v is None else str(v) for v in row) + " |"
        )
    return "\n".join(lines) + "\n"


def build_report(
    summary: Mapping[str, Any],
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    figures: Sequence[Mapping[str, Any]],
    detected: Sequence[str],
    skipped: Sequence[str],
) -> str:
    """Render the QC report as Markdown."""
    out: List[str] = []

    out.append("# Pilot-100 QC Report\n")
    out.append(
        "**Organism:** *Pseudomonas aeruginosa* | "
        "**Antibiotics:** imipenem, meropenem | "
        f"**Generated:** {summary.get('generated_utc')}\n"
    )

    out.append("## Selection method\n")
    out.append(f"**{SELECTION_METHOD}**\n")
    out.append(f"**{SELECTION_METHOD_NOT}**\n")
    out.append(
        "The cohort is defined by what is on disk in `data/`, sorted by GCA "
        "accession, and the first 100 are taken. `PDC_essential.tsv` is opened "
        "only after that selection is final, and is used solely to attach "
        "metadata to the already-selected accessions. Assemblies are reached "
        "through symlinks; `data/` is never modified.\n"
    )

    proof = summary.get("selection_proof") or {}
    out.append("### Selection is materially different from PDC row order\n")
    out.append(
        _table(
            ["Check", "Value"],
            [
                ["Assemblies discovered in data/", summary.get("assemblies_discovered")],
                ["First 100 selected (GCA order)", summary.get("selected")],
                ["PDC_essential.tsv total rows", summary.get("pdc_total_rows")],
                ["First 100 assemblies by GCA accession", proof.get("filesystem_first_accession")],
                ["100th assembly by GCA accession", proof.get("filesystem_last_accession")],
                ["Overlap with PDC row order's first 100", f"{proof.get('overlap_n')} of {summary.get('selected')}"],
                ["Identical to PDC row order's first 100?", proof.get("identical_sets")],
                ["In filesystem selection but not PDC first 100", proof.get("only_in_filesystem_selection")],
                ["In PDC first 100 but not filesystem selection", proof.get("only_in_pdc_roworder")],
            ],
        )
    )
    out.append(
        "Had the cohort been taken as the first 100 rows of "
        f"`PDC_essential.tsv` instead, {proof.get('only_in_pdc_roworder')} of "
        "those assemblies would not be in this cohort and "
        f"{proof.get('only_in_filesystem_selection')} of this cohort's "
        "assemblies would be absent. The two selections are not "
        "interchangeable.\n"
    )

    out.append("## Assembly integrity and the analysis shortfall\n")
    reasons = summary.get("exclusion_reasons") or {}
    out.append(
        "The **selection is the first 100 assemblies in `data/` in GCA order** "
        "and has not been altered. Of those 100, "
        f"**{summary.get('n_excluded')} have no usable sequence file** and are "
        f"therefore excluded from analysis, leaving "
        f"**{summary.get('n_analysable')} analysable**. Excluded assemblies are "
        "dropped, never substituted from elsewhere in the accession list.\n"
    )
    out.append(
        _table(
            ["Exclusion reason", "Assemblies"],
            [[k, v] for k, v in sorted(reasons.items())],
        )
    )
    out.append(
        f"Across all {summary.get('assemblies_discovered')} assemblies present "
        "in `data/`, only a small minority are intact. The pre-existing "
        "`run_pdc.sh` rehydrated a dehydrated NCBI download and reported "
        "success, but most files were left empty, truncated or partial. The "
        f"measured breakdown is {describe_reasons(reasons)}. "
        "**Every count in this report is out of "
        f"{summary.get('n_analysable')} analysable assemblies, not "
        f"{summary.get('n_selected')}.** "
        "Per-assembly detail is in `assembly_integrity.tsv`.\n"
    )

    out.append("## Counts\n")
    out.append(
        _table(
            ["Quantity", "Value"],
            [
                ["Total genome assemblies discovered in data/", summary.get("assemblies_discovered")],
                ["First 100 assemblies selected (SELECTION IS UNCHANGED)", summary.get("n_selected", summary.get("selected"))],
                ["Assemblies analysed (intact and within size range)", summary.get("n_analysable")],
                ["Assemblies EXCLUDED (empty / corrupt / undersized)", summary.get("n_excluded")],
                ["Number matched to PDC_essential.tsv", summary.get("pdc_matched")],
                ["Number missing from PDC_essential.tsv", summary.get("pdc_missing")],
                ["Number with imipenem metadata", summary.get("n_with_imipenem")],
                ["Number with meropenem metadata", summary.get("n_with_meropenem")],
                ["Number of genomes successfully processed", summary.get("n_processed")],
                ["Number of genomes that failed AMR detection", summary.get("n_failed")],
                ["Number of AMR gene calls", summary.get("n_gene_calls")],
                ["Number of distinct AMR genes", summary.get("n_distinct_genes")],
                ["Number of mechanism categories", summary.get("n_mechanisms")],
                ["Number of gene pairs", summary.get("n_gene_pairs")],
                ["Number of mechanism combinations", summary.get("n_mechanism_combinations")],
            ],
        )
    )

    out.append("## Tooling and provenance\n")
    out.append(
        _table(
            ["Item", "Value"],
            [
                ["AMR detection tool", summary.get("amrfinder_version")],
                ["AMR database", summary.get("amr_database")],
                ["AMR database version", summary.get("amr_database_version")],
                ["AMR identity threshold", summary.get("amr_identity_min")],
                ["AMR coverage threshold", summary.get("amr_coverage_min")],
                ["Metadata source", "PDC_essential.tsv"],
                ["Mechanism knowledge table", "config/mechanisms.tsv"],
                ["Gene-family table (ESBL/MBL)", "config/gene_families.tsv"],
            ],
        )
    )

    out.append("## Phenotype availability\n")
    for antibiotic in ("imipenem", "meropenem"):
        rows = tables.get(f"phenotype_{antibiotic}", [])
        counts: Dict[str, int] = {}
        for row in rows:
            label = str(row.get("phenotype"))
            counts[label] = counts.get(label, 0) + 1
        out.append(f"### {antibiotic.capitalize()}\n")
        out.append(
            _table(
                ["Category", "Isolates (n)"],
                [[k, counts.get(k, 0)] for k in sorted(counts)],
            )
        )

    out.append("## Mechanism category distribution\n")
    category_counts: Dict[str, int] = {}
    for row in tables.get("mechanism_summary_imipenem", []):
        category = str(row.get("Mechanism"))
        category_counts[category] = int(row.get("n_isolates") or 0)
    out.append(
        _table(
            ["Mechanism category", "Isolates (n)"],
            [[k, v] for k, v in sorted(category_counts.items(), key=lambda kv: -kv[1])],
        )
    )

    out.append("## Figures\n")
    out.append(
        _table(
            ["Figure", "Created", "Reason if skipped"],
            [[f.get("figure"), f.get("created"), f.get("reason")] for f in figures],
        )
    )

    out.append("## Scientific limitations\n")
    out.append(
        "1. **Gene detection is not phenotypic resistance.** A detected "
        "determinant is reported as DETECTED. No isolate in this pilot is "
        "reported as resistant on the basis of a gene.\n"
        "2. **OprD detection means an intact locus.** Only an explicit "
        "disruptive variant call (premature stop, frameshift, indel) is "
        "reported as reduced permeability. The `oprD` gene being present is "
        "not evidence of reduced permeability.\n"
        "3. **ESBL and MBL are nomenclature labels.** AMRFinderPlus reports a "
        "single `BETA-LACTAM` class and does not distinguish these "
        f"(database {summary.get('amr_database_version')}). The split comes "
        "from `config/gene_families.tsv`, a reviewable family-nomenclature "
        "table. It is not an enzyme activity measurement.\n"
        "4. **No MIC or zone diameter exists or was derived.** The source "
        "provides categorical R/I/S only. R/I/S is never converted into a "
        "quantitative value.\n"
        "5. **No population-structure correction.** This pilot ran no GWAS, "
        "phylogenomics or convergence analysis, so nothing here is corrected "
        "for lineage or clonal structure. Gene and mechanism frequencies are "
        "raw counts within this cohort.\n"
        "6. **Co-occurrence is association only.** Gene and mechanism pairs "
        "are counted per isolate. They are not evidence of interaction, and a "
        "pair may co-occur merely because both belong to the same clone.\n"
        f"7. **The cohort is not a random sample.** These are the "
        f"{summary.get('n_selected')} of "
        f"{summary.get('assemblies_discovered')} available assemblies in GCA "
        "accession order. Nothing here estimates prevalence in the full "
        "collection.\n"
        "8. **PDC metadata is used as given.** Isolate, biosample, collection "
        "date and location are copied from `PDC_essential.tsv` without "
        "independent verification.\n"
        "9. **No phylogenetic or convergence context.** Whether a determinant "
        "occurs in independent lineages is not assessed, so lineage-linked "
        "carriage cannot be distinguished from resistance-associated "
        "carriage.\n"
    )

    if skipped:
        out.append("## Genes with no category\n")
        out.append(
            "Genes AMRFinderPlus reported without a class or subclass, and "
            "which are therefore in no functional category. They are listed "
            "rather than dropped:\n"
        )
        out.append(_table(["Gene"], [[g] for g in sorted(skipped)]))
        out.append("")

    return "".join(out)


def write_report(path: Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    LOGGER.info("QC report written to %s", path)
    return path
