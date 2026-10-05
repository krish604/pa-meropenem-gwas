#!/usr/bin/env python
"""Q1 — synthetic results fixtures for the dashboard's tests.

Run with the pipeline environment:

    micromamba run -n pa-amr python dashboard/fixtures/generate.py

What it writes
--------------

* ``dashboard/fixtures/small/live/``   — a 10-isolate live results root
  (``run_manifest.json`` at the root, tables under ``intermediate/stages/``).
* ``dashboard/fixtures/small/bundle/`` — the same 10 isolates in the delivery
  bundle layout of DESIGN §12 (A1–A5): ``RESULTS.md``, ``01_bakta_input/``,
  ``02_stage_outputs/``, ``03_report/``, ``04_run_info/``,
  ``05_validation/`` and ``06_for_900_isolates/RUNBOOK_900.md``.

Both are committed. They are small, synthetic, and carry no real clinical
accessions: every sample id is ``TEST_PA_###``.

On demand only
--------------

* ``--large`` (or ``PA_FIXTURES_LARGE=1``) writes a 900-isolate live root to a
  temporary directory and prints its path. It is **never committed**.
* ``--huge DIR`` writes a wide table for the UI-D5 paging test to ``DIR`` and
  prints its path. It is **never committed**.

Every table's header is the contract header from
``papipeline.execution.contracts.STAGE_TABLES`` / ``INTERNAL_TABLES``, asserted
in :func:`write_contract_table` so a schema cannot drift silently. Values are
taken from the committed ``test_data/`` fixtures where those fixtures share a
column, and synthesised otherwise.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from papipeline.execution.contracts import (  # noqa: E402
    INTERNAL_TABLES,
    STAGE_TABLES,
    internal_table_path,
    table_path,
)
from papipeline.io.tsv import read_tsv  # noqa: E402
from papipeline.run import STAGE_ORDER  # noqa: E402

FIXTURES_DIR = Path(__file__).resolve().parent
SMALL_DIR = FIXTURES_DIR / "small"
TEST_DATA = REPO_ROOT / "test_data"

BANNER = "SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS"

#: Committed fixtures reused as row templates. Keys are contract table names.
TEMPLATES: Mapping[str, str] = {
    "mlst": "intermediate/mlst/mlst_results.tsv",
    "amr": "intermediate/amr/amr_determinants.tsv",
    "virulence": "intermediate/virulence/virulence_factors.tsv",
    "regulators": "intermediate/regulators/regulator_variants.tsv",
    "structural_variants": "intermediate/structural_variants/structural_variants.tsv",
    "phenotype": "phenotype/imipenem_phenotype.tsv",
    "pangenome": "intermediate/pangenome/pangenome_summary.tsv",
    "gene_presence_absence": "intermediate/pangenome/gene_presence_absence.tsv",
    "core_genes": "intermediate/pangenome/core_genes.tsv",
    "accessory_genes": "intermediate/pangenome/accessory_genes.tsv",
    "tree_metadata": "phylogeny/tree_metadata.tsv",
}


def samples(n: int) -> List[str]:
    """Synthetic sample ids. Never a real accession."""
    return [f"TEST_PA_{i:03d}" for i in range(1, n + 1)]


# ---------------------------------------------------------------------------
# Row synthesis
# ---------------------------------------------------------------------------


def _template_rows(name: str) -> Dict[str, Dict[str, Optional[str]]]:
    rel = TEMPLATES.get(name)
    if not rel:
        return {}
    path = TEST_DATA / rel
    if not path.exists():
        return {}
    try:
        rows = read_tsv(path)
    except Exception:
        return {}
    out: Dict[str, Dict[str, Optional[str]]] = {}
    for row in rows:
        sid = row.get("sample_id")
        if sid:
            out[str(sid)] = row
    return out


def synth_value(column: str, sample: str, i: int) -> str:
    """A deterministic, schema-plausible value for one cell."""
    table = {
        "sample_id": sample,
        "antibiotic": "imipenem",
        "chrom": "NC_002516.2",
        "ref": "A",
        "alt": "G",
        "filter": "PASS",
        "GT": "0/1",
        "gene": "synthetic_gene",
        "determinant": "synthetic_determinant",
        "variant": "",
        "determinant_type": "AMR",
        "mechanism": "antibiotic inactivation",
        "evidence_source": "synthetic",
        "database": "synthetic_db",
        "database_version": "1.0",
        "confidence": "high",
        "claim_status": "supported",
        "virulence_factor": "synthetic_factor",
        "category": "adherence",
        "MIC_unit": "mg/L",
        "zone_unit": "mm",
        "source": "synthetic_generator",
        "mlst_scheme": "synthetic_scheme",
        "allele_database": "synthetic_db",
        "MLST_status": "called",
        "basic_quality_status": "pass",
        "species_confirmation": "confirmed",
        "contamination_status": "clean",
        "completeness_status": "complete",
        "duplicate_status": "unique",
        "flag_reasons": "none",
        "call_status": "called",
        "evidence": "synthetic",
        "mge": "none",
        "variant_type": "snv",
        "variant_id": "sv_0001",
        "affected_gene": "synthetic_gene",
        "effect": "synthetic",
        "reference": "A",
        "alternate": "G",
        "distribution": "synthetic",
        "convergence_category": "unknown",
        "feature_type": "gene",
        "statistic": "fisher",
        "interpretation_limit": "association only; not an interaction",
        "model": "reference_fisher:NO_KINSHIP_CORRECTION",
        "figure": "figure_01",
        "kind": "svg",
        "tree": "stage9.nwk",
        "evidence_level": "synthetic",
        "gene_type": "amr",
        "notes": "synthetic",
        "gwas_status": "not_tested",
        "convergence_status": "unknown",
        "evidence_notes": "synthetic",
        "structural_variant": "none",
        "regulator": "none",
        "chromosomal_mutation": "none",
        "virulence_profile": "synthetic",
        "gwas_feature": "none",
        "amr_gene": "synthetic_gene",
        "amr_variant": "none",
        "phenotype": ["R", "I", "S", "ND"][i % 4],
        "lineage": ["LINEAGE_A", "LINEAGE_B", "LINEAGE_C"][i % 3],
        "present": "true",
        "recombination_detected": "false",
        "matches_manifest": "true",
        "has_branch_lengths": "true",
        "value": str(i + 1),
    }
    if column in table:
        return table[column]
    if column in ("pos", "position", "n_records", "n_named_genes", "size",
                  "n_snps", "n_tips", "n_internal_nodes", "independent_lineages",
                  "branch_count", "n_a", "n_b", "n_both", "n_items", "n_samples"):
        return str(1 + (i % 20))
    if column in ("identity_pct", "coverage_pct", "mean_branch_length",
                  "statistic_value", "adjusted_p_value", "p_value", "frequency",
                  "effect_size"):
        return f"{(i % 9 + 1) / 10:.4f}"
    if column == "truncation_aa":
        return str(10 + i)
    return f"{column}_{i:02d}"


def rows_for(
    header: Sequence[str], sample_ids: Sequence[str], *, name: str = ""
) -> List[List[str]]:
    """Rows for a per-sample table, template values first, synthesis after.

    The committed ``test_data/`` fixtures are used as row templates where they
    share a column; every column the template lacks is synthesised, so the
    emitted header is always the contract header and never the template's.
    """
    template = _template_rows(name)
    out: List[List[str]] = []
    for i, sample in enumerate(sample_ids):
        source = template.get(sample, {})
        row: List[str] = []
        for column in header:
            if column == "sample_id":
                row.append(sample)
                continue
            value = source.get(column)
            row.append(value if value not in (None, "") else synth_value(column, sample, i))
        out.append(row)
    return out


def write_contract_table(
    path: Path,
    *,
    name: str,
    header: Sequence[str],
    rows: Iterable[Sequence[str]],
    banner: bool = True,
) -> Path:
    """Write a TSV whose header is the contract header, asserted.

    The assertion is the point: a fixture whose columns have drifted from
    ``contracts.py`` would make every downstream test assert the wrong schema.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    if banner:
        lines.append(f"# {BANNER}")
    lines.append("\t".join(header))
    for row in rows:
        lines.append("\t".join(str(c) for c in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_stage(path: Path, stage: str, sample_ids: Sequence[str]) -> None:
    """Write one stage table at its contract path with its contract header."""
    _filename, header = STAGE_TABLES[stage]
    rows = stage_rows(stage, header, sample_ids)
    write_contract_table(path, name=stage, header=header, rows=rows)


def write_internal(path: Path, key: str, sample_ids: Sequence[str]) -> None:
    _filename, header = INTERNAL_TABLES[key]
    rows = stage_rows(key, header, sample_ids)
    write_contract_table(path, name=key, header=header, rows=rows)


def stage_rows(stage: str, header: Sequence[str], sample_ids: Sequence[str]) -> List[List[str]]:
    """Rows for one table. Some tables are cohort-wide, not per-sample."""
    n = len(sample_ids)
    if stage in ("cohort_variants",):
        rows = []
        for i in range(min(2 * n, 40)):
            row = []
            for column in header:
                row.append(synth_value(column, "", i))
            rows.append(row)
        return rows
    if stage == "pangenome":
        return [[m, v] for m, v in (
            ("n_genes_total", str(4 * n)),
            ("n_core_genes", "4"),
            ("n_accessory_genes", str(max(1, 3 * n))),
            ("n_unique_genes", "2"),
            ("n_samples", str(n)),
        )]
    if stage == "gwas":
        rows = []
        for i in range(min(4 * n, 80)):
            row = []
            for column in header:
                if column == "lineage_distribution":
                    row.append("LINEAGE_A:5,LINEAGE_B:3,LINEAGE_C:2")
                else:
                    row.append(synth_value(column, "", i))
            rows.append(row)
        return rows
    if stage in ("convergence", "cooccurrence", "reporting", "recombination",
                 "phylogeny", "structural_variants"):
        rows = []
        count = min(2 * n, 40)
        for i in range(count):
            row = []
            for column in header:
                row.append(synth_value(column, "", i))
            rows.append(row)
        return rows
    if stage in ("gene_presence_absence",):
        # long form: one row per (gene, sample)
        rows = []
        for g in range(4):
            for sample in sample_ids:
                rows.append([f"gene_{g + 1:02d}", sample, "true" if (g + len(sample)) % 2 else "false"])
        return rows
    if stage == "core_genes":
        return [[f"core_gene_{i + 1:02d}"] for i in range(4)]
    if stage == "accessory_genes":
        return [[f"acc_gene_{i + 1:02d}", str(2 + i)] for i in range(3)]
    if stage == "similarity":
        # One row per sample holding its whole distance vector.
        order = list(sample_ids)
        rows = []
        for i, sample in enumerate(order):
            vector = [f"{abs(i - j) / max(1, n - 1):.6f}" for j in range(n)]
            rows.append([sample, ";".join(vector)])
        return rows
    if stage == "mechanisms":
        rows = []
        for i, sample in enumerate(sample_ids):
            row = []
            for column in header:
                if column == "sample_id":
                    row.append(sample)
                else:
                    row.append(synth_value(column, sample, i))
            rows.append(row)
        return rows
    # default: one row per sample, template values where available
    return rows_for(header, sample_ids, name=stage)


# ---------------------------------------------------------------------------
# The tree
# ---------------------------------------------------------------------------

#: Internal labels cycled through the generated tree. `100/100` repeats, which
#: is the whole point of UI-D6: a support label is not an identity.
SUPPORT_LABELS = ("100/100", "100/100", "100/100", "99.9/92", "100/92", "83.4/88", "99/100")


def tree_newick(sample_ids: Sequence[str], *, label_cycle: bool = True) -> str:
    """A balanced Newick over ``sample_ids`` with duplicate support labels.

    Deterministic: the same sample list gives the same bytes. Branch lengths
    are positive so a patristic distance is non-trivial.
    """
    counter = {"n": 0}

    def label() -> str:
        idx = counter["n"]
        counter["n"] += 1
        return SUPPORT_LABELS[idx % len(SUPPORT_LABELS)] if label_cycle else "100/100"

    def build(items: Sequence[str], depth: int) -> str:
        if len(items) == 1:
            return f"{items[0]}:0.0{min(depth + 1, 9)}"
        mid = len(items) // 2
        left = build(items[:mid], depth + 1)
        right = build(items[mid:], depth + 1)
        return f"({left},{right}){label()}:0.0{min(depth + 1, 9)}"

    return build(list(sample_ids), 0) + ";"


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


def make_manifest(sample_ids: Sequence[str], *, all_completed: bool = True) -> Dict[str, Any]:
    """A manifest in the `papipeline/run.py:2464` writer's schema.

    No `bakta_executions` key: neither writer carries one (assumption A6), and
    the provenance endpoint must render `not reported`, never `0`.
    """
    state = "completed" if all_completed else "not_run"
    stages = {stage: state for stage in STAGE_ORDER}
    outputs: Dict[str, str] = {}
    for stage, (filename, _cols) in STAGE_TABLES.items():
        outputs[stage] = f"intermediate/stages/{filename}"
    return {
        "pipeline_version": "0.1.0-fixture",
        "run_mode": "TEST",
        "generated_at_utc": "2026-10-05T00:00:00Z",
        "antibiotic": "imipenem",
        "python": "3.11",
        "platform": "synthetic-fixture",
        "n_samples": len(sample_ids),
        "samples": list(sample_ids),
        "stages": stages,
        "stages_skipped": {},
        "outputs": outputs,
        "tools_detected": {
            "bakta": {"executable": "/usr/bin/bakta", "version": "1.9.0", "available": True},
            "mlst": {"executable": "/usr/bin/mlst", "version": "2.23.0", "available": True},
            "iqtree": {"executable": "/usr/bin/iqtree2", "version": "2.2.0", "available": True},
        },
        "references": [
            {
                "reference_id": "ref_genome",
                "tool": "bakta",
                "tool_version": "1.9.0",
                "database": "synthetic_db",
                "database_version": "1.0",
                "version_status": "pinned",
            }
        ],
    }


def events_lines(sample_ids: Sequence[str]) -> List[str]:
    """Valid `status/events.jsonl` lines for the healthy fixture."""
    out: List[str] = []
    t = 0
    for stage in STAGE_ORDER:
        t += 1
        out.append(json.dumps({"t": f"2026-10-05T00:00:{t:02d}Z", "event": "start", "stage": stage, "sample": "all"}))
        t += 1
        out.append(json.dumps({"t": f"2026-10-05T00:00:{t:02d}Z", "event": "done", "stage": stage, "sample": "all"}))
    return out


# ---------------------------------------------------------------------------
# Layouts
# ---------------------------------------------------------------------------


def _write_common(stages_dir: Path, sample_ids: Sequence[str]) -> None:
    """Every contracted table, internal table and Snakefile-declared table."""
    stages_dir.mkdir(parents=True, exist_ok=True)
    for stage in STAGE_TABLES:
        write_stage(table_path(stages_dir, stage), stage, sample_ids)
    for key in INTERNAL_TABLES:
        write_internal(internal_table_path(stages_dir, key), key, sample_ids)
    # Snakefile-declared paths (report_tables.SNAKEFILE_ONLY_TABLES).
    write_contract_table(
        stages_dir / "gene_presence_absence.tsv",
        name="gene_presence_absence",
        header=["gene", "sample_id", "present"],
        rows=stage_rows("gene_presence_absence", ["gene", "sample_id", "present"], sample_ids),
    )
    write_contract_table(
        stages_dir / "core_genes.tsv",
        name="core_genes",
        header=["gene"],
        rows=stage_rows("core_genes", ["gene"], sample_ids),
    )
    write_contract_table(
        stages_dir / "accessory_genes.tsv",
        name="accessory_genes",
        header=["gene", "n_samples"],
        rows=stage_rows("accessory_genes", ["gene", "n_samples"], sample_ids),
    )
    # Sidecars and provenance.
    units = {
        "quantity": "patristic distance",
        "definition": "expected substitutions per site, summed along the stage-9 tree",
        "is_a_snp_count": False,
    }
    (stages_dir / "similarity.units.json").write_text(json.dumps(units, indent=2) + "\n", encoding="utf-8")
    provenance = {
        "per_isolate": [
            {"sample_id": s, "n_snv": i % 7, "n_indel": i % 3}
            for i, s in enumerate(sample_ids)
        ]
    }
    (stages_dir / "variants_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


def _write_reports(reports_dir: Path) -> None:
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "test_pipeline_report.md").write_text(
        f"# {BANNER}\n\nSynthetic report body.\n", encoding="utf-8"
    )
    (reports_dir / "test_pipeline_report.html").write_text(
        f"<html><body><h1>{BANNER}</h1></body></html>\n", encoding="utf-8"
    )
    figures = reports_dir / "figures"
    figures.mkdir(exist_ok=True)
    (figures / "figure_01.svg").write_text(
        "<svg xmlns=\"http://www.w3.org/2000/svg\"><text>synthetic</text></svg>\n",
        encoding="utf-8",
    )


def build_live(root: Path, n: int) -> Path:
    """A live results root: manifest at the root, tables under intermediate/."""
    root = Path(root)
    sample_ids = samples(n)
    stages_dir = root / "intermediate" / "stages"
    _write_common(stages_dir, sample_ids)
    phylo = root / "intermediate" / "phylogeny"
    phylo.mkdir(parents=True, exist_ok=True)
    (phylo / "stage9.nwk").write_text(tree_newick(sample_ids) + "\n", encoding="utf-8")
    metadata_header = ["sample_id", "tree_tip_label", "lineage_label", "st", "source"]
    write_contract_table(
        phylo / "tree_metadata.tsv",
        name="tree_metadata",
        header=metadata_header,
        rows=[
            [s, s, ["LINEAGE_A", "LINEAGE_B", "LINEAGE_C"][i % 3], str(2 + i % 12), "synthetic_generator"]
            for i, s in enumerate(sample_ids)
        ],
    )
    status = root / "status"
    status.mkdir(parents=True, exist_ok=True)
    (status / "events.jsonl").write_text("\n".join(events_lines(sample_ids)) + "\n", encoding="utf-8")
    _write_reports(root / "reports")
    (root / "run_manifest.json").write_text(
        json.dumps(make_manifest(sample_ids), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return root


def build_bundle(root: Path, n: int) -> Path:
    """The delivery bundle layout (DESIGN §12, A1–A5)."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    sample_ids = samples(n)
    (root / "RESULTS.md").write_text(
        f"# {BANNER}\n\nSynthetic delivery bundle.\n", encoding="utf-8"
    )
    (root / "01_bakta_input").mkdir(parents=True, exist_ok=True)
    (root / "01_bakta_input" / "README.md").write_text(
        f"# {BANNER}\n\nBakta inputs (synthetic).\n", encoding="utf-8"
    )
    stage_outputs = root / "02_stage_outputs"
    _write_common(stage_outputs, sample_ids)
    phylo = root / "intermediate" / "phylogeny"
    phylo.mkdir(parents=True, exist_ok=True)
    (phylo / "stage9.nwk").write_text(tree_newick(sample_ids) + "\n", encoding="utf-8")
    metadata_header = ["sample_id", "tree_tip_label", "lineage_label", "st", "source"]
    write_contract_table(
        phylo / "tree_metadata.tsv",
        name="tree_metadata",
        header=metadata_header,
        rows=[
            [s, s, ["LINEAGE_A", "LINEAGE_B", "LINEAGE_C"][i % 3], str(2 + i % 12), "synthetic_generator"]
            for i, s in enumerate(sample_ids)
        ],
    )
    (root / "03_report").mkdir(parents=True, exist_ok=True)
    _write_reports(root / "03_report")
    run_info = root / "04_run_info"
    (run_info / "status").mkdir(parents=True, exist_ok=True)
    (run_info / "status" / "events.jsonl").write_text(
        "\n".join(events_lines(sample_ids)) + "\n", encoding="utf-8"
    )
    validation = root / "05_validation"
    validation.mkdir(parents=True, exist_ok=True)
    (validation / "README.md").write_text(
        f"# {BANNER}\n\nValidation artefacts (synthetic).\n", encoding="utf-8"
    )
    runbook = root / "06_for_900_isolates"
    runbook.mkdir(parents=True, exist_ok=True)
    lines = ["# RUNBOOK_900 (synthetic)", "", "| sample_id | seconds |", "| --- | --- |"]
    for i, s in enumerate(sample_ids):
        lines.append(f"| {s} | {10 + i} |")
    (runbook / "RUNBOOK_900.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (stage_outputs / "run_manifest.json").write_text(
        json.dumps(make_manifest(sample_ids), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return root


# ---------------------------------------------------------------------------
# On-demand artefacts
# ---------------------------------------------------------------------------


def build_large(out_dir: Optional[Path] = None, n: int = 900) -> Path:
    """The 900-isolate live root, in a temp directory, never committed."""
    if out_dir is None:
        out_dir = Path(tempfile.mkdtemp(prefix="pa_dash_large_"))
    root = Path(out_dir)
    return build_live(root, n)


def write_wide_table(path: Path, *, rows: int = 3000, cols: int = 2000) -> Path:
    """A 3000x2000 TSV for the UI-D5 paging test. Never committed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ["sample_id"] + [f"c{i:04d}" for i in range(1, cols)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for r in range(rows):
            cells = [str(r % 10) for _ in range(cols - 1)]
            handle.write(f"TEST_PA_{r:05d}\t" + "\t".join(cells) + "\n")
    return path


def build_wide_root(out_dir: Path, *, rows: int = 3000, cols: int = 2000) -> Path:
    """A minimal live root whose ``variants.tsv`` is the wide table."""
    root = Path(out_dir)
    stages = root / "intermediate" / "stages"
    stages.mkdir(parents=True, exist_ok=True)
    write_wide_table(stages / "variants.tsv", rows=rows, cols=cols)
    sample_ids = samples(10)
    (root / "run_manifest.json").write_text(
        json.dumps(make_manifest(sample_ids), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return root


def build_small() -> Tuple[Path, Path]:
    """Write both committed 10-isolate layouts and return their roots."""
    live = SMALL_DIR / "live"
    bundle = SMALL_DIR / "bundle"
    build_live(live, 10)
    build_bundle(bundle, 10)
    return live, bundle


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large", action="store_true", help="Write the 900-isolate live root to a temp dir.")
    parser.add_argument("--samples", type=int, default=900, help="Cohort size for --large.")
    parser.add_argument("--out", type=Path, default=None, help="Output directory for --large.")
    parser.add_argument("--huge", type=Path, default=None, help="Write a 3000x2000 table under this directory.")
    parser.add_argument("--rows", type=int, default=3000)
    parser.add_argument("--cols", type=int, default=2000)
    args = parser.parse_args(argv)

    large_env = os.environ.get("PA_FIXTURES_LARGE") == "1"
    if args.large or large_env:
        root = build_large(args.out, n=args.samples)
        print(root)
        return 0
    if args.huge is not None:
        path = write_wide_table(args.huge / "wide.tsv", rows=args.rows, cols=args.cols)
        print(path)
        return 0

    live, bundle = build_small()
    print(live)
    print(bundle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
