"""The wiring test for the post-scan statistics package (written first).

`papipeline/downstream/` shipped complete and unit-tested
(`tests/unit/test_downstream_*.py`, 85 tests) and was called by nothing:
`.build/downstream-stats.registry-notes.md` records the verdict as "built,
unwired", in the order its own section 1 prescribes — gate, conditional,
stratified, interactions, lineage meta, evidence, report.

This file pins the wiring and the one decision the wiring had to make.

**The decision: the seven steps run in REAL, and only in REAL.** Measured, not
assumed: on the committed twenty-isolate fixtures stage 12 produces
`gene__oprD_absent` nowhere, `gene__oprD_LoF` at adjusted p = 0.901, and no
MBL feature at all, so `run_control_gate` refuses with **both** configured
controls missing (`oprD_burden` present-but-not-significant, `any_MBL` absent).
The gate is the FIRST of the seven steps and `build_report` re-runs it internally
with no flag to disable it, so wiring the package into a mode whose baseline
cannot pass would either fail every TEST run or require weakening the gate on
behalf of the fixtures — which is tuning configuration to make a test pass, the
thing `cohort_gate` already refused to do for the same run.

What that leaves testable, and this file tests all of it:

1. the runner's *body* — every one of the seven steps, in order, against the
   real `GwasInput` a TEST run builds — by calling it directly with a baseline
   that recovers the controls (the same technique `test_downstream_report.py`
   uses, because the fixtures themselves cannot);
2. that the gate is genuinely first: the fixture baseline is refused before a
   single file is written;
3. that `run.py` calls it from inside `resolved is RunMode.REAL`, read from the
   source with `ast` rather than trusting a comment;
4. that a TEST run writes no downstream artefact, so the mode decision cannot
   drift into "wired everywhere, discovered at the end of a real run".
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import List

import pytest

from papipeline.config.loader import load_config
from papipeline.downstream import controls as cg
from papipeline.downstream import report as rp
from papipeline.downstream.runner import DownstreamRun, run_downstream_analyses
from papipeline.io import read_tsv
from papipeline.models import GwasResult
from papipeline.run import run_pipeline

EVIDENCE_FILENAME = "downstream_evidence.tsv"


def _result(feature: str, adjusted_p: float, lineages=None) -> GwasResult:
    """A baseline row that recovers a control. Shape from the unit tests."""
    return GwasResult(
        feature=feature,
        feature_type="gene",
        effect=3.0,
        p_value=min(1.0, adjusted_p / 10.0),
        adjusted_p_value=adjusted_p,
        effect_size=3.0,
        frequency=0.2,
        lineage_distribution=dict(lineages or {"ST1": 4, "ST2": 4}),
        model="test",
    )


def _recovering_baseline() -> List[GwasResult]:
    return [
        _result("gene_presence_absence__blaNDM-1", 1e-6),
        _result("gene__oprD_absent", 1e-6),
    ]


@pytest.fixture(scope="module")
def test_run_and_scan(pipeline_root: Path):
    """One TEST run, with stage 12's (results, GwasInput) pair captured.

    `RunResult` does not carry them — they are run-locals of
    `run_pipeline` — so the call site is observed rather than reimplemented:
    the original `stage_gwas.run` is wrapped, called exactly as the pipeline
    calls it, and its return value kept. The pipeline itself is otherwise
    untouched, so what is captured is what a real run would hold.
    """
    import papipeline.run as run_module

    captured = {}
    original = run_module.stage_gwas.run

    def spy(*args, **kwargs):
        out = original(*args, **kwargs)
        captured["gwas"] = out
        return out

    run_module.stage_gwas.run = spy
    try:
        result = run_pipeline(
            config_path=pipeline_root / "config" / "science.yaml",
            mode="TEST",
            config=load_config(pipeline_root / "config" / "science.yaml"),
        )
    finally:
        run_module.stage_gwas.run = original

    assert "gwas" in captured, "stage 12 never ran in this TEST run"
    gwas_results, gwas_input = captured["gwas"]
    return result, gwas_results, gwas_input


class TestTheSevenStepsRun:
    def test_a_recovering_baseline_runs_every_step(self, tmp_path, test_run_and_scan):
        """The body is exercised end to end, gate through report.

        Steps, in the order section 1 of the registry notes prescribes:
        control gate, conditional scan, stratified scan, interactions (both
        tiers), lineage meta, evidence, report. All seven are asserted
        through the object the runner returns, so a step that stops being
        called shows up as a missing section rather than as a green run.
        """
        _, gwas_results, gwas_input = test_run_and_scan
        # The committed baseline cannot recover the controls (see the module
        # docstring), so the gate is given a baseline that can — the fixture
        # GwasInput is real either way, which is what the wiring is about.
        del gwas_results
        outcome = run_downstream_analyses(
            load_config(Path("config") / "science.yaml"),
            gwas_results=_recovering_baseline(),
            gwas_input=gwas_input,
            convergence=None,
            output_dir=tmp_path,
        )
        assert isinstance(outcome, DownstreamRun)
        assert [c.control for c in outcome.controls] == ["oprD_burden", "any_MBL"]
        assert all(c.recovered for c in outcome.controls)

        report = outcome.report
        assert report.tier1, "step 4 (tier 1) produced no rows"
        assert report.tier2, "step 4 (tier 2) produced no rows"
        assert report.meta, "step 5 (lineage meta) did not run"
        assert report.conditional, "step 2 (conditional scan) produced no rows"
        assert report.evidence, "step 6 (evidence) produced no rows"
        assert report.controls == tuple(outcome.controls)

    def test_the_evidence_table_is_written_and_reads_back(self, tmp_path, test_run_and_scan):
        """Step 6 and 7: the evidence table on disk is the one the report holds."""
        _, _, gwas_input = test_run_and_scan
        outcome = run_downstream_analyses(
            load_config(Path("config") / "science.yaml"),
            gwas_results=_recovering_baseline(),
            gwas_input=gwas_input,
            output_dir=tmp_path,
        )
        assert outcome.evidence_path == tmp_path / EVIDENCE_FILENAME
        rows = read_tsv(
            outcome.evidence_path,
            required_columns=list(rp.REPORT_COLUMNS),
            unique_columns=["feature"],
        )
        assert len(rows) == len(outcome.report.rows())
        assert rows, "an evidence table with no rows has nothing to review"

    def test_the_stratified_scan_holds_the_known_determinants_fixed(
        self, tmp_path, test_run_and_scan
    ):
        """Step 3 runs over the fixture cohort, not over an empty stratum."""
        _, _, gwas_input = test_run_and_scan
        outcome = run_downstream_analyses(
            load_config(Path("config") / "science.yaml"),
            gwas_results=_recovering_baseline(),
            gwas_input=gwas_input,
            output_dir=tmp_path,
        )
        stratified = outcome.report.stratified
        assert stratified is not None, "step 3 (stratified scan) did not run"
        assert set(stratified.known_columns) <= set(gwas_input.features)
        assert len(stratified.samples) < gwas_input.n_samples, (
            "the stratum is the whole cohort, so nothing was held out"
        )


class TestTheGateIsFirst:
    def test_the_fixture_baseline_is_refused_before_anything_is_written(
        self, tmp_path, test_run_and_scan
    ):
        """The measured fact the mode decision rests on, pinned so it can't rot.

        If a future fixture regenerates a significant MBL feature and an
        OprD-loss feature, this fails - which is the signal to reconsider the
        REAL-only gate in `run.py`, not a reason to weaken this assertion.
        """
        _, gwas_results, gwas_input = test_run_and_scan
        with pytest.raises(cg.ControlGateError) as exc:
            run_downstream_analyses(
                load_config(Path("config") / "science.yaml"),
                gwas_results=gwas_results,
                gwas_input=gwas_input,
                output_dir=tmp_path,
            )
        assert set(exc.value.context["missing"]) == {"any_MBL", "oprD_burden"}
        assert exc.value.context["key"] == "downstream.positive_controls"
        assert list(tmp_path.iterdir()) == [], (
            "the gate refused AFTER writing something; it is not first"
        )


class TestTheModeGate:
    def test_run_py_calls_it_from_inside_the_real_branch(
        self, pipeline_root: Path
    ):
        """Read from the source: the call exists, and it is inside REAL.

        A comment saying "REAL only" proves nothing; the `ast` walk does. The
        same technique `test_seam_run_module_references.py` uses on this file,
        for the same reason - a call site that quietly moves out of the guard
        is invisible to a TEST run, because a TEST run never reaches it.
        """
        source = (pipeline_root / "papipeline" / "run.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)

        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "run_downstream_analyses"
        ]
        assert len(calls) == 1, (
            f"expected exactly one call site, found {len(calls)}"
        )

        def names(node) -> set:
            return {
                child.id if isinstance(child, ast.Name) else child.attr
                for child in ast.walk(node)
                if isinstance(child, (ast.Name, ast.Attribute))
            }

        guarded = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.If)
            and calls[0] in list(ast.walk(node))  # somewhere under this If
            and "RunMode" in names(node.test)
            and "REAL" in names(node.test)
        ]
        assert guarded, (
            "run_downstream_analyses is called, but not from a branch whose "
            "test names RunMode.REAL"
        )

    def test_a_test_run_writes_no_downstream_artefact(
        self, test_run_and_scan, pipeline_root: Path
    ):
        """The decision on the other side of the branch, asserted on a real run."""
        result, _, _ = test_run_and_scan
        assert "downstream_evidence" not in result.outputs, (
            "a TEST run produced the downstream report; the mode gate changed"
        )
        reports_root = pipeline_root / "results" / "reports"
        assert not (reports_root / EVIDENCE_FILENAME).exists(), (
            f"{reports_root / EVIDENCE_FILENAME} exists after a TEST run"
        )
