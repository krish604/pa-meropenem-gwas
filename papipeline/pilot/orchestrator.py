"""Pilot-100 analysis orchestration.

``run_analysis()`` performs, in order:

1. assembly QC (stage 1 of the main pipeline, reused unchanged);
2. AMRFinderPlus detection across the cohort;
3. phenotype extraction for imipenem and meropenem from the joined PDC
   metadata;
4. gene, mechanism, pair and combination tables per antibiotic;
5. figures;
6. PowerPoint;
7. QC report.

It reuses the main pipeline's stage 1 and its typed record models rather than
reimplementing them.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig, load_config
from ..errors import PipelineError
from ..io.tsv import read_tsv, write_tsv
from ..logging_utils import get_logger
from ..manifest import Sample, SampleManifest
from ..models import RunMode
from ..stages import validation as validation_stage
from . import analysis as analysis_mod
from . import figures as figures_mod
from . import mechanisms as mech_mod
from . import pptx as pptx_mod
from . import report as report_mod
from .amr_detect import detect_cohort, database_version

LOGGER = get_logger("pilot.orchestrator")

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "pilot100"
FIG_DIR = ROOT / "figures"
AMR_DIR = OUT_DIR / "amrfinder"


def _qc_summary(members: Sequence[analysis_mod.CohortMember], config: PipelineConfig):
    """Run the main pipeline's stage 1 over the pilot assemblies."""
    manifest = SampleManifest(
        [
            Sample(m.pilot_id, m.genome_path, "pilot100")
            for m in members
        ]
    )
    records = validation_stage.run(config, manifest, RunMode.REAL)
    summary = validation_stage.summarise(records)
    rows = validation_stage.rows(records)
    return records, rows, summary


def run_analysis(
    antibiotic_labels: Sequence[str] = analysis_mod.ANTIBIOTICS,
    threads: int = 2,
    jobs: int = 4,
) -> Dict[str, Any]:
    """Run the full pilot analysis on the prepared cohort."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    AMR_DIR.mkdir(parents=True, exist_ok=True)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    config = load_config(ROOT / "config" / "science.yaml", machine=args.machine)
    families = mech_mod.gene_family_table(ROOT / "config" / "gene_families.tsv")

    # --- load the prepared cohort ------------------------------------
    selected = analysis_mod.load_cohort(
        OUT_DIR / "pilot100_manifest.tsv",
        OUT_DIR / "first100_assembly_metadata.tsv",
    )
    if not selected:
        raise PipelineError("Pilot cohort is empty; run --prepare-only first")
    LOGGER.info("Pilot cohort selected: %d assemblies", len(selected))

    # The SELECTION is the first 100 in GCA order and is never changed here.
    # Assemblies whose file is empty, corrupt or undersized cannot be
    # analysed; they are excluded with a recorded reason, never replaced.
    members = [m for m in selected if m.analysis_included]
    excluded = [m for m in selected if not m.analysis_included]
    exclusion_reasons: Dict[str, int] = {}
    for member in excluded:
        key = (member.exclusion_reason or "unknown").split("_at_offset")[0]
        for prefix in ("corrupt_binary_data", "undersized", "oversized",
                       "zero_byte_file", "file_not_found", "not_fasta_no_header"):
            if key.startswith(prefix):
                key = prefix
                break
        exclusion_reasons[key] = exclusion_reasons.get(key, 0) + 1

    LOGGER.info(
        "Analysable: %d of %d selected (excluded %s)",
        len(members),
        len(selected),
        ", ".join(f"{k}={v}" for k, v in sorted(exclusion_reasons.items())) or "none",
    )
    if excluded:
        LOGGER.error(
            "%d of the %d selected assemblies have an empty, corrupt or "
            "undersized file and are EXCLUDED from the analysis. They are not "
            "substituted. See results/pilot100/assembly_integrity.tsv.",
            len(excluded),
            len(selected),
        )
    if not members:
        raise PipelineError(
            "No analysable assemblies: every selected file is empty, corrupt "
            "or undersized. Re-download is required before any analysis.",
            n_selected=len(selected),
        )

    # --- 1. assembly QC -----------------------------------------------
    qc_records, qc_rows, qc_summary = _qc_summary(members, config)
    write_tsv(
        OUT_DIR / "assembly_QC.tsv",
        qc_rows,
        list(validation_stage.QC_COLUMNS),
        header_comment=["Stage 1 genome validation, reused from the main pipeline."],
    )
    LOGGER.info(
        "Assembly QC: %d flagged of %d", qc_summary.get("n_flagged", 0), len(members)
    )

    # --- 1b. purge reports for excluded assemblies --------------------
    # A report produced from a corrupt or partial assembly is not a result.
    # Removing it here means detect_cohort cannot silently reuse it.
    included_ids = {m.pilot_id for m in members}
    if AMR_DIR.exists():
        stale = [
            p
            for p in AMR_DIR.glob("*.amrfinder.tsv")
            if p.name.split(".")[0] not in included_ids
        ]
        for path in stale:
            path.unlink()
        if stale:
            LOGGER.warning(
                "Removed %d AMRFinderPlus report(s) belonging to excluded "
                "assemblies",
                len(stale),
            )

    # --- 2. AMR detection ---------------------------------------------
    database_dir = _locate_database()
    items = [(m.pilot_id, Path(m.genome_path)) for m in members]
    detection = detect_cohort(
        items, AMR_DIR, database_dir, "imipenem", threads=threads, jobs=jobs
    )
    analysis_mod.attach_detections(members, detection, config)

    if detection.n_failed:
        LOGGER.warning(
            "%d of %d assemblies failed AMR detection and carry no gene calls",
            detection.n_failed,
            len(items),
        )

    # --- 3-4. per-antibiotic tables -----------------------------------
    tables: Dict[str, List[Dict[str, Any]]] = {}
    gene_rows_by_antibiotic: Dict[str, List[Dict[str, Any]]] = {}
    mechanism_rows_by_antibiotic: Dict[str, List[Dict[str, Any]]] = {}
    pair_rows_by_antibiotic: Dict[str, List[Dict[str, Any]]] = {}
    n_phenotype: Dict[str, int] = {}

    for antibiotic in antibiotic_labels:
        label = antibiotic.capitalize()

        phenotype_rows, counts = analysis_mod.phenotype_rows(members, antibiotic)
        n_phenotype[antibiotic] = sum(
            n for category, n in counts.items() if category != analysis_mod.MISSING
        )
        tables[f"phenotype_{antibiotic}"] = phenotype_rows
        write_tsv(
            OUT_DIR / f"{label}_phenotype.tsv",
            phenotype_rows,
            [
                "Pilot_ID",
                "Assembly",
                "Isolate",
                "phenotype",
                "availability",
                "MIC",
                "MIC_unit",
            ],
            header_comment=[
                "Categorical susceptibility only. No MIC or zone diameter is "
                "present in the source and none is derived.",
            ],
        )

        gene_rows = analysis_mod.gene_summary_rows(members, antibiotic, config, families)
        gene_rows_by_antibiotic[antibiotic] = gene_rows
        write_tsv(
            OUT_DIR / f"{label}_gene_summary.tsv",
            gene_rows,
            [
                "Pilot_ID",
                "Assembly",
                "Isolate",
                "AST_phenotype",
                "Gene",
                "Mechanism",
                "PA_mechanism",
                "PA_evidence_level",
                "Detected",
                "Evidence",
                "AMR_class",
                "AMR_subclass",
                "n_genes_detected_in_isolate",
            ],
            header_comment=[
                f"AMRFinderPlus detection, database {database_version(database_dir)}.",
                "Detection is not phenotypic resistance.",
            ],
        )

        freq_rows = analysis_mod.gene_frequency_rows(members, antibiotic, config, families)
        write_tsv(
            OUT_DIR / f"{label}_gene_frequency.tsv",
            freq_rows,
            list(freq_rows[0].keys()) if freq_rows else ["Gene"],
        )

        mech_rows = analysis_mod.mechanism_summary_rows(members, antibiotic, families)
        mechanism_rows_by_antibiotic[antibiotic] = mech_rows
        tables[f"mechanism_summary_{antibiotic}"] = mech_rows
        write_tsv(
            OUT_DIR / f"{label}_mechanism_summary.tsv",
            mech_rows,
            list(mech_rows[0].keys()) if mech_rows else ["Mechanism"],
        )

        pair_rows = analysis_mod.gene_pair_rows(members, antibiotic)
        pair_rows_by_antibiotic[antibiotic] = pair_rows
        write_tsv(
            OUT_DIR / f"{label}_gene_pairs.tsv",
            pair_rows,
            list(pair_rows[0].keys()) if pair_rows else ["Gene_A", "Gene_B"],
            header_comment=["Association only. Not evidence of interaction or causation."],
        )

        combination_rows = analysis_mod.mechanism_combination_rows(members, antibiotic, families)
        tables[f"mechanism_combinations_{antibiotic}"] = combination_rows
        write_tsv(
            OUT_DIR / f"{label}_mechanism_combinations.tsv",
            combination_rows,
            list(combination_rows[0].keys()) if combination_rows else ["Combination"],
            header_comment=["Observed combinations only. Association only, not causation."],
        )
        LOGGER.info(
            "%s: %d gene rows, %d distinct genes, %d pairs, %d combinations",
            label,
            len(gene_rows),
            len(freq_rows),
            len(pair_rows),
            len(combination_rows),
        )

    # --- 5. figures ---------------------------------------------------
    figure_results: List[figures_mod.FigureResult] = []
    figure_paths: Dict[str, Optional[Path]] = {}
    for antibiotic in antibiotic_labels:
        results = figures_mod.build_all(
            FIG_DIR,
            antibiotic,
            gene_rows_by_antibiotic[antibiotic],
            analysis_mod.gene_frequency_rows(members, antibiotic, config, families),
            mechanism_rows_by_antibiotic[antibiotic],
            pair_rows_by_antibiotic[antibiotic],
            members,
        )
        figure_results.extend(results)
        for result in results:
            figure_paths[result.name] = result.path

    for result in figures_mod.build_shared(FIG_DIR, members, families):
        figure_results.append(result)
        figure_paths[result.name] = result.path

    write_tsv(
        OUT_DIR / "figure_manifest.tsv",
        [r.to_row() for r in figure_results],
        ["figure", "created", "path", "reason"],
    )
    for result in figure_results:
        state = "created" if result.created else f"SKIPPED ({result.reason})"
        LOGGER.info("figure %-40s %s", result.name, state)

    # --- summary ------------------------------------------------------
    all_genes: set = set()
    for member in members:
        all_genes |= member.genes
    uncategorised = sorted(
        g
        for g in all_genes
        if mech_mod.functional_category(g, None, families)
        == mech_mod.CATEGORY_UNMAPPED
    )
    mechanism_categories: set = set()
    for member in members:
        mechanism_categories |= member.categories()

    prepare_summary = json.loads(
        (OUT_DIR / "pilot100_prepare_summary.json").read_text(encoding="utf-8")
    )

    summary: Dict[str, Any] = {
        "generated_utc": generated,
        "assemblies_discovered": prepare_summary.get("assemblies_discovered"),
        "selected": prepare_summary.get("selected"),
        "n_excluded_prepare": prepare_summary.get("n_excluded"),
        "pdc_total_rows": prepare_summary.get("pdc_total_rows"),
        "pdc_matched": prepare_summary.get("pdc_matched"),
        "pdc_missing": prepare_summary.get("pdc_missing"),
        "selection_proof": prepare_summary.get("selection_proof", {}),
        "selection_overlap": prepare_summary.get("selection_proof", {}).get("overlap_n"),
        "n_selected": len(selected),
        "n_excluded": len(excluded),
        "exclusion_reasons": exclusion_reasons,
        "n_with_imipenem": n_phenotype.get("imipenem"),
        "n_with_meropenem": n_phenotype.get("meropenem"),
        "n_analysable": len(members),
        "n_processed": detection.n_successful,
        "n_failed": detection.n_failed,
        "n_gene_calls": sum(len(v) for v in detection.records.values()),
        "n_distinct_genes": len(all_genes),
        "n_mechanisms": len(mechanism_categories),
        "n_gene_pairs": len(
            pair_rows_by_antibiotic.get(antibiotic_labels[0], [])
        ),
        "n_mechanism_combinations": len(
            tables.get(f"mechanism_combinations_{antibiotic_labels[0]}", [])
        ),
        "amrfinder_version": _amrfinder_version(),
        "amr_database": f"AMRFinderPlus ({database_version(database_dir)})",
        "amr_database_version": database_version(database_dir),
        "amr_identity_min": "80 (AMRFinderPlus default)",
        "amr_coverage_min": "50 (AMRFinderPlus default)",
        "qc_summary": {
            "n_samples": qc_summary.get("n_samples"),
            "n_flagged": qc_summary.get("n_flagged"),
            "min_assembly_size": qc_summary.get("min_assembly_size"),
            "max_assembly_size": qc_summary.get("max_assembly_size"),
            "median_assembly_size": qc_summary.get("median_assembly_size"),
            "min_n50": qc_summary.get("min_n50"),
            "max_n50": qc_summary.get("max_n50"),
            "min_gc_content": qc_summary.get("min_gc_content"),
            "max_gc_content": qc_summary.get("max_gc_content"),
        },
    }
    (OUT_DIR / "pilot100_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )

    # --- 6. PowerPoint -------------------------------------------------
    deck_tables = dict(tables)
    deck_tables["qc_summary"] = [
        {"metric": k, "value": v} for k, v in summary["qc_summary"].items()
    ]
    deck_tables["mechanism_combinations"] = tables.get(
        f"mechanism_combinations_{antibiotic_labels[0]}", []
    )
    pptx_path = pptx_mod.build_deck(
        OUT_DIR / "Pilot100_Pseudomonas_AMR_Analysis.pptx",
        summary,
        deck_tables,
        figure_paths,
        antibiotic_labels=antibiotic_labels,
    )
    LOGGER.info("PowerPoint written to %s", pptx_path)

    # --- 7. QC report --------------------------------------------------
    report_text = report_mod.build_report(
        summary,
        tables,
        [r.to_row() for r in figure_results],
        detected=sorted(all_genes),
        skipped=uncategorised,
    )
    report_mod.write_report(OUT_DIR / report_mod.REPORT_NAME, report_text)

    print()
    print("Pilot-100 analysis complete")
    print(f"  selected (first 100, GCA)  : {summary['n_selected']}")
    print(f"  excluded (bad assembly)    : {summary['n_excluded']} "
          f"({', '.join(f'{k}={v}' for k, v in sorted(summary['exclusion_reasons'].items())) or 'none'})")
    print(f"  genomes processed          : {summary['n_processed']} / {summary['n_analysable']} analysable")
    print(f"  AMR gene calls             : {summary['n_gene_calls']}")
    print(f"  distinct genes             : {summary['n_distinct_genes']}")
    print(f"  mechanism categories       : {summary['n_mechanisms']}")
    print(f"  gene pairs ({antibiotic_labels[0]})       : {summary['n_gene_pairs']}")
    print(f"  mechanism combinations     : {summary['n_mechanism_combinations']}")
    print(f"  figures created            : {sum(1 for r in figure_results if r.created)}"
          f" / {len(figure_results)}")
    print(f"  outputs                    : {OUT_DIR}")
    return summary


def _locate_database() -> Path:
    """Find the provisioned AMRFinderPlus database.

    Preference order:

    1. ``$AMRFINDERPLUS_DB``, if set;
    2. the ``latest`` entry of the env's ``share/amrfinderplus``, which is
       what ``amrfinder_update`` maintains;
    3. the lexicographically highest versioned directory there;
    4. the database directory itself, if it is already versioned.

    An ``AMRFINDERPLUS_DB`` override that does not exist is a hard error
    rather than a silent fallback, so a typo cannot be masked.
    """
    from .amr_detect import is_database_dir

    override = os.environ.get("AMRFINDERPLUS_DB")
    if override:
        candidate = Path(override)
        if is_database_dir(candidate):
            return candidate
        raise PipelineError(
            "AMRFINDERPLUS_DB is set but is not a provisioned database",
            path=str(candidate),
            hint="Point it at the directory containing version.txt",
        )

    root = (
        Path.home()
        / "micromamba_pilot/envs/pilot100/share/amrfinderplus"
    )
    if not root.is_dir():
        raise PipelineError(
            "AMRFinderPlus database root not found",
            path=str(root),
            hint=(
                "Provision it with: micromamba run -n pilot100 amrfinder_update "
                "--database <env>/share/amrfinderplus"
            ),
        )

    latest = root / "latest"
    if is_database_dir(latest):
        return latest

    versioned = sorted(
        (p for p in root.iterdir() if p.is_dir() and is_database_dir(p)),
        key=lambda p: p.name,
    )
    if not versioned:
        raise PipelineError(
            "No provisioned AMRFinderPlus database found",
            path=str(root),
            hint="Run amrfinder_update to provision one",
        )
    return versioned[-1]


def _amrfinder_version() -> str:
    import shutil
    import subprocess

    executable = shutil.which("amrfinder")
    if not executable:
        return "not installed"
    try:
        out = subprocess.run(
            [executable, "--version"], capture_output=True, text=True, timeout=30
        )
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception:  # noqa: BLE001 - provenance must never break a run
        return "unknown"
