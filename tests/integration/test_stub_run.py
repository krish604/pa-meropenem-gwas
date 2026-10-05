"""One real STUB run must produce every stage output and a well-formed event log.

This is the highest seam in the project (spec.md, Seam 1) and it is deliberately
an *execution* test, not a dry run. ``snakemake -n`` proves the DAG resolves; it
proves nothing about whether a rule can actually run, whether a stage writes
where the next rule expects, or whether the dashboard receives a single
parseable event stream. Those are the regressions that a dry run cannot see,
and the DAG was green for a long time while the workflow had never executed.

What STUB must do (spec.md D8): exercise the whole DAG with no bioinformatics
tool installed, no real parser, and no fixture on disk, emitting a tiny
plausible-shaped output per rule.

The run is hermetic. ``PIPELINE_RESULTS_ROOT`` redirects every output and
``--config status_log=`` redirects the event stream, so the test never reads or
writes the developer's real results tree or the real dashboard feed.

The results redirect is an *environment variable* rather than a ``--config`` key
on purpose. The rules shell out to a separate process that recomputes its paths
from config, so a ``--config`` override reaches the DAG and not the stage that
writes the file.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV
from papipeline.execution.contracts import STAGE_TABLES

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = PIPELINE_ROOT / "workflow" / "Snakefile"

pytestmark = pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake is not installed; this check cannot run here",
)


#: Outputs the workflow declares that are not a stage's principal table.
EXTRA_OUTPUTS = (
    "gene_presence_absence.tsv",
    "core_genes.tsv",
    "accessory_genes.tsv",
    "10_alignment_summary.tsv",
)


def _expected_stage_tables() -> set:
    return {filename for filename, _columns in STAGE_TABLES.values()}


def _stub_run(tmp_path: Path, *, cores: int = 1) -> tuple:
    """Run the whole workflow for real, in STUB mode, into a temp tree."""
    results_root = tmp_path / "results"
    status_log = tmp_path / "events.jsonl"
    environ = {**os.environ, RESULTS_ROOT_ENV: str(results_root)}
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--config",
            "mode=STUB",
            "machine=laptop",
            f"status_log={status_log}",
            "--cores", str(cores),
            "--rerun-incomplete",
        ],
        cwd=PIPELINE_ROOT,
        env=environ,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    return result, results_root, status_log


def test_a_stub_run_produces_every_stage_output(tmp_path: Path):
    """Every declared stage table must exist after one STUB run."""
    result, results_root, _ = _stub_run(tmp_path)
    assert result.returncode == 0, (
        f"the STUB run failed; nothing downstream of this is meaningful\n"
        f"--- stdout ---\n{result.stdout[-4000:]}\n"
        f"--- stderr ---\n{result.stderr[-4000:]}"
    )

    stage_dir = results_root / "stub" / "intermediate" / "stages"
    missing = [
        name
        for name in sorted(_expected_stage_tables() | set(EXTRA_OUTPUTS))
        if not (stage_dir / name).exists()
    ]
    assert not missing, (
        f"these declared outputs were not produced by a STUB run: {missing}\n"
        f"looked in {stage_dir}"
    )


def test_every_stage_table_has_its_declared_header(tmp_path: Path):
    """A stub output must be shaped like the real thing, not just present.

    A file the dashboard cannot parse is not a stub, it is a different bug.
    """
    result, results_root, _ = _stub_run(tmp_path)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]

    stage_dir = results_root / "stub" / "intermediate" / "stages"
    for stage, (_filename, columns) in sorted(STAGE_TABLES.items()):
        path = stage_dir / _filename
        header = path.read_text(encoding="utf-8").splitlines()[0]
        actual = tuple(header.split("\t"))
        assert actual == columns, (
            f"stage {stage!r} wrote a header that does not match its contract\n"
            f"  expected: {columns}\n  wrote:    {actual}"
        )


def test_a_stub_run_writes_a_well_formed_event_log(tmp_path: Path):
    """The dashboard's only input must exist and parse.

    Four fields, three event values, nothing else. A fifth field is a contract
    change the browser has to understand, and there is no second consumer to
    justify it.
    """
    result, _results_root, status_log = _stub_run(tmp_path)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]

    assert status_log.exists(), f"no event log at {status_log}"
    lines = [ln for ln in status_log.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert lines, "the event log is empty; the dashboard would show nothing"

    records = [json.loads(ln) for ln in lines]
    for record in records:
        assert set(record) == {"t", "event", "stage", "sample"}, (
            f"event record has the wrong fields: {sorted(record)}"
        )
        assert record["event"] in {"start", "done", "fail"}, (
            f"unknown event value {record['event']!r}"
        )
        assert record["stage"], "event record has no stage"

    assert {r["event"] for r in records} == {"start", "done"}, (
        "a clean run should record exactly start and done"
    )

    stages = {r["stage"] for r in records}
    assert len(stages) >= 15, f"only {len(stages)} stages reported: {sorted(stages)}"


def test_a_stub_run_needs_no_fixture_on_disk(tmp_path: Path):
    """STUB must not depend on committed fixtures or real genomes.

    The run is pointed at a results root that does not exist. If a stage still
    reached for test_data/ or data/, it would either fail or quietly succeed
    against a fixture, and STUB would mean TEST.
    """
    result, _results_root, _status_log = _stub_run(tmp_path)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    # STUB mode writes nothing under the data root.
    assert not (PIPELINE_ROOT / "data" / "stub").exists(), (
        "a STUB run wrote into data/, which means it is not fabricating"
    )


def test_run_pipeline_script_has_a_usable_machine_flag(tmp_path: Path):
    """`run_pipeline.py` reads `args.machine` but never defined the flag.

    The reporting rule shells out to it, so the whole workflow died at the last
    stage with an AttributeError from argparse - a pre-existing defect that
    only a real execution could surface, since `--dry-run` never runs a command.
    """
    script = PIPELINE_ROOT / "scripts" / "common" / "run_pipeline.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True, text=True, timeout=120,
    )
    assert "--machine" in result.stdout, (
        "run_pipeline.py reads args.machine but does not accept --machine, so "
        "main() raises AttributeError before any pipeline work happens"
    )


# --------------------------------------------------------------------------
# The directional join, wired
# --------------------------------------------------------------------------

def _test_root_with_phenotype_gap(tmp_path: Path) -> Path:
    """A throwaway pipeline root whose phenotype table is missing one sample.

    Copies config/ and test_data/ so the gap can be introduced without touching
    the committed fixtures, then drops a single phenotype row.
    """
    import shutil as _shutil

    root = tmp_path / "root"
    _shutil.copytree(PIPELINE_ROOT / "config", root / "config")
    _shutil.copytree(PIPELINE_ROOT / "test_data", root / "test_data")

    phenotype = root / "test_data" / "phenotype" / "imipenem_phenotype.tsv"
    kept = []
    dropped = None
    data_rows = 0
    for line in phenotype.read_text(encoding="utf-8").splitlines(keepends=True):
        is_comment = line.startswith("#")
        # The first non-comment line is the header; the second is the first
        # sample. Fixture identifiers are TEST_PA_*, not GCA_*.
        if not is_comment:
            data_rows += 1
            if data_rows == 2:
                dropped = line.split("\t")[0]
                continue
        kept.append(line)
    assert dropped, "no data row found to drop"
    phenotype.write_text("".join(kept), encoding="utf-8")
    return root


def test_a_missing_phenotype_row_stops_the_run(tmp_path: Path, monkeypatch):
    """The manifest is the authority on what exists.

    A genome with no phenotype row is a hard failure, not a warning: the
    pipeline cannot know which cohort was meant, and carrying the sample as
    silently-missing is how a GWAS ends up on a different set than the report
    describes. This is `papipeline.join`'s contract, and before it was wired
    the run only warned and continued.
    """
    from papipeline.config.loader import load_config
    from papipeline.errors import StageError
    from papipeline.run import run_pipeline

    root = _test_root_with_phenotype_gap(tmp_path)
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "gap-results"))
    config = load_config(root / "config" / "science.yaml")

    # run.py wraps a stage failure in StageError and preserves the detail, so
    # the operator sees why, not just which stage died.
    with pytest.raises(StageError) as excinfo:
        run_pipeline(config=config, mode="TEST")

    message = str(excinfo.value)
    assert "TEST_PA_" in message, (
        f"the failure should name the offending sample, got: {message}"
    )
    assert "phenotype row" in message, (
        f"the failure should say what is wrong, got: {message}"
    )


def test_the_join_writes_its_exclusions_for_a_real_run(tmp_path: Path, monkeypatch):
    """A clean TEST run records the join's decision, not just its verdict.

    The exclusions table is the join's audit trail: which phenotype rows were
    dropped and why. Without it, "excluded, with its reason" is a claim nobody
    can check.
    """
    from papipeline.config.loader import load_config
    from papipeline.run import run_pipeline

    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "join-results"))
    config = load_config(PIPELINE_ROOT / "config" / "science.yaml")
    result = run_pipeline(config=config, mode="TEST")

    exclusions = result.outputs.get("phenotype_exclusions")
    assert exclusions is not None, (
        "the phenotype stage ran but recorded no exclusions table, so the "
        "directional join's decisions are not auditable"
    )
    assert Path(exclusions).exists()
