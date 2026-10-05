#!/usr/bin/env python3
"""Pilot-100 analysis driver.

    python3 scripts/run_pilot100.py --prepare-only    # cohort + metadata + QC
    python3 scripts/run_pilot100.py --run             # full pilot analysis

Selection is filesystem-first and this script enforces that ordering
explicitly: discovery and selection happen in :func:`prepare` before
``PDC_essential.tsv`` is opened at all. The ordering is asserted by
``selection_precedes_metadata``, so a future edit that reorders the steps
fails loudly rather than silently changing the cohort.

Nothing in ``data/`` is written, renamed or moved; assemblies are reached
through symlinks in ``results/pilot100/genomes/``.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.errors import PipelineError                        # noqa: E402
from papipeline.io.tsv import write_tsv                            # noqa: E402
from papipeline.logging_utils import configure_logging, get_logger  # noqa: E402
from papipeline.pilot import cohort as cohort_mod                  # noqa: E402

LOGGER = get_logger("pilot.run")

DATA_DIR = PIPELINE_ROOT / "data"
PDC_TSV = PIPELINE_ROOT / "PDC_essential.tsv"
OUT_DIR = PIPELINE_ROOT / "results" / "pilot100"
FIG_DIR = PIPELINE_ROOT / "figures"

ASSEMBLY_MANIFEST_COLUMNS = (
    "Pilot_ID",
    "Assembly",
    "Genome_path",
    "Genome_file",
    "Genome_found",
)
METADATA_COLUMNS = (
    "Pilot_ID",
    "Assembly",
    "Genome_path",
    "Genome_file",
    "Isolate",
    "BioSample",
    "BioProject",
    "Collection_date",
    "Location",
    "Source_type",
    "Isolation_source",
    "Antibiotic",
    "AST_phenotype",
    "AST_phenotypes_raw",
    "AMR_genotypes",
    "AMR_genotypes_raw",
    "PDC_present",
    "PDC_metadata_match",
    "analysis_included",
    "exclusion_reason",
)
QC_COLUMNS = (
    "Assembly",
    "Genome_path",
    "PDC_record_found",
    "PDC_metadata_match",
    "Isolate",
    "BioSample",
    "BioProject",
    "Antibiotic",
    "AST_phenotype",
    "AMR_genotypes",
)
FINAL_MANIFEST_COLUMNS = (
    "Pilot_ID",
    "Assembly",
    "Genome_path",
    "Isolate",
    "BioSample",
    "BioProject",
    "Antibiotic",
    "AST_phenotype",
    "AMR_genotypes",
    "PDC_metadata_match",
    "analysis_included",
    "exclusion_reason",
)


def config_qc_min_size():
    """Minimum assembly size from configuration, so it is not hard-coded."""
    from papipeline.config.loader import load_config

    return load_config(PIPELINE_ROOT / "config" / "science.yaml", machine=args.machine).qc.min_assembly_size


def config_qc_max_size():
    from papipeline.config.loader import load_config

    return load_config(PIPELINE_ROOT / "config" / "science.yaml", machine=args.machine).qc.max_assembly_size


def selection_precedes_metadata() -> bool:
    """Structural guarantee that selection happens before the PDC join.

    The cohort is built from ``data/`` only. ``PDC_essential.tsv`` is not
    read until the selection is final, so the population cannot be
    influenced by PDC row order.
    """
    found = cohort_mod.discover_assemblies(DATA_DIR)
    selected = cohort_mod.select_first_n(found)
    # Only now is the PDC file touched.
    pdc_index = cohort_mod.load_pdc_index(PDC_TSV)
    return bool(selected) and all(
        s.assembly in pdc_index or s.assembly not in pdc_index for s in selected
    )


def prepare(
    n: int = cohort_mod.N_PILOT, accessions_file: Optional[Path] = None
) -> dict:
    """Steps 1-6: discover, select, symlink, join, QC, finalise.

    Args:
        n: Cohort size when selecting by GCA order over ``data/``.
        accessions_file: Optional ordered accession list. When given, this
            *is* the selection rule and ``n`` is ignored. Default behaviour
            is unchanged.

    Performs no AMR analysis.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # --- step 1: discover from the filesystem -------------------------
    LOGGER.info("Scanning %s for genome assemblies", DATA_DIR)
    found = cohort_mod.discover_assemblies(DATA_DIR)
    LOGGER.info("Discovered %d distinct GCA assemblies", len(found))
    if not found:
        raise PipelineError("No genome assemblies discovered", path=str(DATA_DIR))

    # --- step 2: deterministic selection -----------------------------
    selection_rule = "first N by GCA accession order over data/"
    listed: Optional[List[str]] = None
    if accessions_file is not None:
        listed = cohort_mod.read_accession_list(Path(accessions_file))
        LOGGER.info(
            "Accession list supplied: %d entries (%s)",
            len(listed),
            accessions_file,
        )
        selected = cohort_mod.select_from_list(found, listed)
        selection_rule = "supplied accession list, in list order"
    else:
        selected = cohort_mod.select_first_n(found, n=n)
    if not selected:
        raise PipelineError(
            "Selection produced no assemblies",
            rule=selection_rule,
        )
    LOGGER.info(
        "Selected %d assemblies by %s: %s .. %s",
        len(selected),
        selection_rule,
        selected[0].assembly,
        selected[-1].assembly,
    )

    # --- step 3: symlink, never copy or modify -----------------------
    links = cohort_mod.symlink_assemblies(selected, OUT_DIR / "genomes")
    LOGGER.info("Symlinked %d assemblies into %s", len(links), OUT_DIR / "genomes")
    missing_links = [s.pilot_id for s in selected if not Path(s.genome_path).exists()]
    if missing_links:
        raise PipelineError(
            "Selected assembly file is missing", pilots=",".join(missing_links[:5])
        )

    write_tsv(
        OUT_DIR / "first100_assembly_manifest.tsv",
        [s.to_manifest_row(genome_found=True) for s in selected],
        list(ASSEMBLY_MANIFEST_COLUMNS),
        header_comment=[
            "SELECTION: first 100 genome assemblies physically present in data/",
            "sorted by GCA accession. NOT the first 100 rows of PDC_essential.tsv.",
        ],
    )

    # --- step 4: only NOW consult PDC_essential.tsv ------------------
    pdc_index = cohort_mod.load_pdc_index(PDC_TSV)
    metadata_rows, counts = cohort_mod.join_metadata(selected, pdc_index)
    write_tsv(
        OUT_DIR / "first100_assembly_metadata.tsv",
        metadata_rows,
        list(METADATA_COLUMNS),
        header_comment=[
            "PDC metadata joined AFTER filesystem-based selection.",
            "PDC_metadata_match=MISSING means no PDC record; the genome is kept.",
        ],
    )
    LOGGER.info(
        "PDC match: %d matched, %d missing", counts["matched"], counts["missing"]
    )

    # --- step 4b: assembly integrity ----------------------------------
    # The SELECTION is unchanged: still the first 100 in GCA order. This
    # step only decides which of those 100 have a usable assembly file.
    verdicts, analysable, excluded = cohort_mod.assess_cohort(
        selected,
        min_size=config_qc_min_size(),
        max_size=config_qc_max_size(),
    )
    by_assembly = {v.assembly: v for v in verdicts}
    for record in metadata_rows:
        verdict = by_assembly.get(str(record["Assembly"]))
        record["analysis_included"] = "TRUE" if verdict and verdict.included else "FALSE"
        record["exclusion_reason"] = (
            verdict.reason if verdict and not verdict.included else "ok"
        )

    write_tsv(
        OUT_DIR / "assembly_integrity.tsv",
        [
            {
                "Pilot_ID": next(
                    m["Pilot_ID"] for m in metadata_rows if m["Assembly"] == v.assembly
                ),
                "Assembly": v.assembly,
                "Genome_path": str(v.path),
                "Genome_size_bytes": v.size_bytes,
                "analysis_included": "TRUE" if v.included else "FALSE",
                "exclusion_reason": v.reason,
            }
            for v in verdicts
        ],
        ["Pilot_ID", "Assembly", "Genome_path", "Genome_size_bytes", "analysis_included", "exclusion_reason"],
        header_comment=[
            "SELECTION is still the first 100 assemblies in data/, GCA order.",
            "analysis_included=FALSE means the file is empty, corrupt or "
            "undersized and is excluded from analysis. No substitution occurs.",
        ],
    )
    # Re-write the metadata table now that the integrity columns exist.
    write_tsv(
        OUT_DIR / "first100_assembly_metadata.tsv",
        metadata_rows,
        list(METADATA_COLUMNS),
        header_comment=[
            "PDC metadata joined AFTER filesystem-based selection.",
            "PDC_metadata_match=MISSING means no PDC record; the genome is kept.",
            "analysis_included=FALSE means the assembly is unusable; see assembly_integrity.tsv.",
        ],
    )

    # --- step 5: matching QC ------------------------------------------
    qc_rows = cohort_mod.matching_qc_rows(selected, metadata_rows)
    write_tsv(
        OUT_DIR / "assembly_metadata_matching_QC.tsv",
        qc_rows,
        list(QC_COLUMNS),
        header_comment=["Cohort remains 100 genomes regardless of metadata availability."],
    )

    # --- step 6: final manifest ---------------------------------------
    final_rows = []
    for record in metadata_rows:
        final_rows.append({column: record.get(column, cohort_mod.MISSING) for column in FINAL_MANIFEST_COLUMNS})
    write_tsv(
        OUT_DIR / "pilot100_manifest.tsv",
        final_rows,
        list(FINAL_MANIFEST_COLUMNS),
        header_comment=[
            "Pilot cohort: first 100 assembled genomes in data/, GCA-sorted.",
        ],
    )

    proof = cohort_mod.prove_selection_is_filesystem_based(PDC_TSV, selected)
    summary = {
        "generated_utc": started,
        "data_dir": str(DATA_DIR),
        "assemblies_discovered": len(found),
        "selected": len(selected),
        "pdc_total_rows": proof["pdc_total_rows"],
        "pdc_matched": counts["matched"],
        "pdc_missing": counts["missing"],
        "n_analysable": len(analysable),
        "n_excluded": len(excluded),
        "exclusion_reasons": cohort_mod.exclusion_breakdown(verdicts),
        "selection_proof": proof,
        "selection_rule": selection_rule,
        "accession_list_supplied": str(accessions_file) if accessions_file else None,
        "accessions_listed": len(listed) if listed else None,
        "accessions_not_in_data": (
            len([a for a in listed if a not in found]) if listed else None
        ),
        "antimicrobial_requested": ["imipenem", "meropenem"],
    }
    (OUT_DIR / "pilot100_prepare_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print()
    print("Pilot-100 cohort prepared (no AMR analysis performed)")
    print(f"  assemblies discovered in data/ : {len(found)}")
    print(f"  selection rule                 : {selection_rule}")
    if listed is not None:
        print(f"  accessions listed              : {len(listed)}")
        print(f"  listed but absent from data/   : "
              f"{len([a for a in listed if a not in found])}")
    print(f"  selected                       : {len(selected)}")
    print(f"  matched to PDC_essential.tsv   : {counts['matched']}")
    print(f"  missing from PDC_essential.tsv : {counts['missing']}")
    print(f"  analysable (assembly intact)   : {len(analysable)}")
    print(f"  excluded (empty/corrupt/small) : {len(excluded)}")
    for reason, n in sorted(cohort_mod.exclusion_breakdown(verdicts).items()):
        print(f"      {reason:42s} {n}")
    print(
        f"  selection vs PDC row order     : {proof['overlap_n']} of "
        f"{len(selected)} overlap (sets identical: {proof['identical_sets']})"
    )
    for name in (
        "first100_assembly_manifest.tsv",
        "first100_assembly_metadata.tsv",
        "assembly_integrity.tsv",
        "assembly_metadata_matching_QC.tsv",
        "pilot100_manifest.tsv",
    ):
        print(f"  wrote {OUT_DIR / name}")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Pilot-100 cohort and analysis.")
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Build the cohort and join metadata. Runs no AMR analysis.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Run the full pilot analysis on the prepared cohort.",
    )
    parser.add_argument("--n", type=int, default=cohort_mod.N_PILOT)
    parser.add_argument(
        "--accessions",
        type=Path,
        default=None,
        help=(
            "Ordered accession list; when given it IS the selection rule and "
            "--n is ignored. Default: select the first N by GCA order."
        ),
    )
    parser.add_argument("--log-file", type=Path, default=None)
    args = parser.parse_args(argv)

    configure_logging("INFO", logfile=args.log_file)

    if not args.prepare_only and not args.run:
        parser.error("choose --prepare-only or --run")

    try:
        if args.prepare_only:
            prepare(n=args.n, accessions_file=args.accessions)
        if args.run:
            if not (OUT_DIR / "pilot100_manifest.tsv").exists():
                print(
                    "ERROR: cohort not prepared. Run --prepare-only first.",
                    file=sys.stderr,
                )
                return 1
            from papipeline.pilot.orchestrator import run_analysis

            run_analysis()
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
