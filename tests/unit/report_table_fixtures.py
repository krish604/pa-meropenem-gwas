"""Shared builders for the report's tests.

Both `test_report_reads_run_tables.py` and `test_report_figures.py` need the
same thing: a synthetic `<intermediate>/stages` directory holding tables written
at the paths `papipeline/execution/contracts.py` declares, with the headers that
same module declares. Duplicating that in two files would let the two drift, and
a fixture that drifts from the contracts silently stops exercising the reader -
it starts exercising the refusal path instead.

Imported by name rather than as a package: `tests/unit/` has no `__init__.py`,
so pytest puts the directory on `sys.path` and a plain module import is the only
form that resolves.
"""

from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.execution.contracts import (
    INTERNAL_TABLES,
    internal_table_path,
    required_columns,
    table_path,
)
from papipeline.models import RunMode
from papipeline.stages import reporting
from papipeline.stages.report_tables import load_run_tables

SCIENCE = "config/science.yaml"

#: The ten-isolate cohort R16 names. Used wherever the literal flag matters.
COHORT_N = 10

#: Every heading `SECTION_BUILDERS` promises, so a section dropped from that
#: tuple cannot disappear from the report unnoticed.
REQUIRED_SECTIONS: Sequence[str] = (
    "Run summary and per-stage status",
    "Cohort overview",
    "AMR summary",
    "Virulence summary",
    "oprD locus verdicts",
    "Variants summary",
    "Pangenome summary",
    "Tree and similarity summary",
    "Associations: GWAS, convergence, co-occurrence",
    "Provenance of the tables read",
)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def write_table(
    path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str] = ()
) -> Path:
    """Write a TSV carrying ``columns`` as its header, in that order.

    Every table the report reads is written here with the header
    `contracts.required_columns` declares for it, because
    `load_run_tables` passes those columns to `read_tsv` as a requirement and
    will - correctly - refuse a table that carries a different one. A test
    fixture with a convenient short header would test the refusal path instead
    of the reading path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    header = list(columns) or (list(rows[0].keys()) if rows else ["sample_id"])
    lines = ["\t".join(header)]
    for row in rows:
        lines.append(
            "\t".join("" if row.get(c) is None else str(row[c]) for c in header)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_contract_table(
    root: Path, key: str, rows: Sequence[Mapping[str, object]]
) -> Path:
    """Write a stage table at its contracted path, with its contracted header."""
    return write_table(table_path(root, key), rows, required_columns(key))


def header_only(root: Path, key: str) -> Path:
    """A table carrying its contracted header and no rows at all."""
    path = table_path(root, key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(required_columns(key)) + "\n", encoding="utf-8")
    return path


def amr_row(**over: object) -> Dict[str, object]:
    row: Dict[str, object] = {
        "sample_id": "TEST_PA_001",
        "antibiotic": "imipenem",
        "determinant": "oprD_V359L",
        "gene": "oprD",
        "variant": None,
        "determinant_type": "AMR",
        "mechanism": "loss_of_function",
        "evidence_source": "amrfinderplus",
        "database": "AMRFinderPlus",
        "database_version": "2026-08-07.1",
        "confidence": None,
        "claim_status": "DETECTED",
        "identity_pct": 94.58,
        "coverage_pct": 100.0,
    }
    row.update(over)
    return row


def variant_row(**over: object) -> Dict[str, object]:
    row: Dict[str, object] = {
        "sample_id": "TEST_PA_001",
        "chrom": "NC_002516.2",
        "pos": 779,
        "ref": "C",
        "alt": "G",
        "qual": 30.4,
        "filter": ".",
        "GT": "1/1",
        "AC": 2,
        "AN": 2,
        "DP4": "0,0,1,0",
        "MQ": 60,
        "MQ0F": 0,
    }
    row.update(over)
    return row


@pytest.fixture()
def stage_dir(tmp_path: Path) -> Path:
    """An empty `<intermediate>/stages` directory."""
    root = tmp_path / "results" / "real" / "intermediate" / "stages"
    root.mkdir(parents=True)
    return root


def context(config: PipelineConfig, mode: RunMode, **over) -> reporting.ReportContext:
    ctx = reporting.ReportContext(
        mode=mode,
        config=config,
        generated_at="2026-10-04T00:00:00Z",
        antibiotic="imipenem",
        n_samples=COHORT_N,
    )
    for key, value in over.items():
        setattr(ctx, key, value)
    return ctx


def real_guard_payload(mode: RunMode) -> Dict[str, object]:
    """The declared inputs a REAL report must carry to avoid its own refusal.

    `write_report` refuses a REAL report whose three finding-bearing aggregates
    are empty, so a REAL test that reaches the HTML render has to hand over
    content for them. Only REAL needs this; the guard is not what is under test
    here (`tests/unit/test_reporting_real_refusal.py` is).
    """
    if mode is not RunMode.REAL:
        return {}
    return {
        "declared_inputs": {
            "convergence_calls": [
                reporting.ConvergenceCategory.RARE_ISOLATED
            ],
            "cooccurrence": [{"feature_a": "a", "feature_b": "b"}],
            "gwas_associations": [{"feature": "f", "adjusted_p_value": 0.01}],
        }
    }


def populated_tables(config: PipelineConfig, root: Path) -> None:
    """Write one small, synthetic table per contracted output.

    Every table carries the header `contracts.required_columns` declares for it,
    because the reader requires those columns and a convenient short header
    would make each of these assert the refusal path instead.
    """
    write_contract_table(
        root,
        "validation",
        [
            {"sample_id": f"TEST_PA_{i:03d}", "assembly_size": 6_400_000}
            for i in range(1, 4)
        ],
    )
    write_contract_table(
        root,
        "mlst",
        [
            {"sample_id": f"TEST_PA_{i:03d}", "ST": str(i), "MLST_status": "ok"}
            for i in range(1, 4)
        ],
    )
    write_contract_table(
        root,
        "amr",
        [
            amr_row(sample_id="TEST_PA_001"),
            amr_row(
                sample_id="TEST_PA_002",
                determinant="oprD_Y237TerfsTer0",
                variant="Y237TerfsTer0",
            ),
            amr_row(
                sample_id="TEST_PA_003",
                determinant="blaKPC-2",
                gene="blaKPC-2",
                mechanism="enzymatic_inactivation",
            ),
        ],
    )
    write_contract_table(
        root,
        "virulence",
        [
            {
                "sample_id": "TEST_PA_001",
                "virulence_factor": "lasA",
                "gene": "lasA",
                "category": "protease",
            }
        ],
    )
    write_contract_table(
        root,
        "variants",
        [
            variant_row(),
            variant_row(sample_id="TEST_PA_002", ref="C", alt="CT"),
            variant_row(sample_id="TEST_PA_003", ref="CTT", alt="C"),
            variant_row(sample_id="TEST_PA_003", ref="CA", alt="GT", pos=900),
        ],
    )
    write_contract_table(
        root,
        "cohort_variants",
        [
            {
                "chrom": "NC_002516.2",
                "pos": 779,
                "ref": "C",
                "alt": "G",
                "ac": 2,
                "an": 4,
                "af": 0.5,
            }
        ],
    )
    write_contract_table(
        root, "pangenome", [{"metric": "core_gene_families", "value": 4828}]
    )
    write_contract_table(
        root,
        "phylogeny",
        [
            {
                "tree": "tree.nwk",
                "n_tips": 3,
                "n_internal_nodes": 1,
                "has_branch_lengths": "true",
                "matches_manifest": "true",
            }
        ],
    )
    write_contract_table(
        root,
        "phenotype",
        [
            {"sample_id": "TEST_PA_001", "phenotype": "R", "MIC": 8},
            {"sample_id": "TEST_PA_002", "phenotype": "S", "MIC": 1},
            {"sample_id": "TEST_PA_003", "phenotype": "I", "MIC": 2},
        ],
    )
    write_contract_table(
        root,
        "gwas",
        [
            {
                "feature": "gene__oprD",
                "feature_type": "gene",
                "adjusted_p_value": 0.01,
            }
        ],
    )
    write_contract_table(
        root,
        "convergence",
        [
            {
                "determinant": "oprD",
                "independent_lineages": 3,
                "convergence_category": "rare_isolated",
            }
        ],
    )
    write_contract_table(
        root,
        "cooccurrence",
        [
            {
                "feature_a": "oprD",
                "feature_b": "blaKPC-2",
                "n_both": 2,
                "adjusted_p_value": 0.5,
            }
        ],
    )
    write_contract_table(
        root, "similarity", [{"sample_id": f"TEST_PA_{i:03d}"} for i in range(1, 4)]
    )
    from papipeline.stages.similarity import units_sidecar_path

    units_sidecar_path(table_path(root, "similarity")).write_text(
        '{"method": "patristic", "n_tips": 3, "n_pairs": 3}', encoding="utf-8"
    )
    write_table(
        internal_table_path(root, "mechanisms"),
        [{"sample_id": "TEST_PA_001", "mechanism": "loss_of_function"}],
        INTERNAL_TABLES["mechanisms"][1],
    )
    write_table(
        root / "10_alignment_summary.tsv",
        [{"metric": "core_n_sequences", "value": 3}],
        ["metric", "value"],
    )


def render(config: PipelineConfig, mode: RunMode, root: Optional[Path] = None) -> str:
    tables = load_run_tables(config, mode, stage_dir=root)
    return reporting.build_markdown(context(config, mode, run_tables=tables))
