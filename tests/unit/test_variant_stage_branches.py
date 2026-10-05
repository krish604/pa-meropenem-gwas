"""The two stage-6 branches in `run.py::_execute_stage`.

The ticket's blocker was that these stages had parsers and no dispatcher:
`run.py`'s `variants` branch wrote only the regulators internal table and its
own comment admitted it was "unreachable until ticket 14 builds the stage", and
there was no `cohort_variants` branch at all. `UNBUILT_STAGES` refused both
outside STUB, so the branches were not reachable through `run_pipeline`.

The reconciliation that followed opened that gate: both stages left
`UNBUILT_STAGES` and `standins.UNBUILT_WITH_STANDIN`, and `analysis:` in
`science.yaml` gained the two keys. These tests therefore assert the branches
*and* the opened gate — and, because opening a stage is not permission to read
a cohort, assert separately that REAL is still refused on every machine overlay.

What is tested here is that the branches are *correct now that they are
reachable*, by inspecting the dispatcher itself. That is a real seam rather
than a workaround: it is the function the Snakemake rules call and the function
the observatory wraps, so a test against it is a test against what a run will
actually do.

Two properties are load-bearing.

**The branches hand `stage_*.run()` to `run.py`, not to a private helper.**
Every other stage in the file is a thin dispatch — `record(...)` around one
`stage_X.run(...)` call. A branch that re-implemented the merge inline would
work, pass these tests, and leave two implementations of the cohort filter free
to disagree. So these assert the wiring is the same shape as its neighbours.

**`cohort_variants` consumes what `variants` produced, in the same run.**
The two branches are separate stages with a declared edge, and the merge's
denominator is the cohort. If the cohort branch recomputed calls independently
it would be reading a different thing from the one the per-isolate stage wrote.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from papipeline.execution.contracts import (
    STAGE_TABLES,
    internal_table_path,
    table_path,
)
from papipeline.run import EXECUTION_ORDER, STAGE_ORDER, UNBUILT_STAGES

REPO = Path(__file__).resolve().parents[2]
RUN_PY = REPO / "papipeline" / "run.py"


def _branch_body(stage: str) -> ast.AST:
    """The `elif stage == "<stage>":` node inside `_execute_stage`.

    Each stage in the chain is one `ast.If` whose `.test` is the comparison. Note
    the operands: in `stage == "variants"` the `left` is the *name* `stage` and
    the stage name is the comparator, so matching on `left.value` finds nothing.
    """
    tree = ast.parse(RUN_PY.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef) and node.name == "_execute_stage"):
            continue
        for inner in ast.walk(node):
            if not isinstance(inner, ast.If):
                continue
            test = inner.test
            if (
                isinstance(test, ast.Compare)
                and isinstance(test.comparators[0], ast.Constant)
                and test.comparators[0].value == stage
            ):
                return inner
    raise AssertionError(f"no branch for stage {stage!r} in _execute_stage")


def _branch_source(stage: str) -> str:
    """The branch's own statements, unparsed.

    Only `.body`, never the whole `If`: an `elif` node carries the rest of the
    chain in `.orelse`, so unparsing the node itself would fold every later
    stage into this branch and make a negative assertion ("this branch does not
    call X") pass for the wrong reason.

    `ast.unparse` normalises string quotes to single, so assertions here match
    on substrings that do not contain a quoted literal, or accept either quote.
    """
    node = _branch_body(stage)
    return "\n".join(ast.unparse(stmt) for stmt in node.body)


class TestTheBranchesExist:
    @pytest.mark.parametrize("stage", ("variants", "cohort_variants"))
    def test_the_stage_has_a_branch(self, stage):
        assert _branch_body(stage) is not None

    def test_variants_came_before_cohort_variants(self):
        """The merge consumes the calls, so its branch is downstream."""
        assert STAGE_ORDER.index("variants") < STAGE_ORDER.index("cohort_variants")
        assert EXECUTION_ORDER.index("variants") < EXECUTION_ORDER.index("cohort_variants")


class TestTheBranchesDispatchToTheStages:
    """Not to a private helper. A second implementation is the failure mode."""

    def test_variants_calls_the_stage_entry_point(self):
        assert "stage_variants.run(" in _branch_source("variants")

    def test_cohort_variants_calls_the_stage_entry_point(self):
        assert "stage_cohort_variants.run(" in _branch_source("cohort_variants")

    def test_neither_branch_reimplements_the_merge(self):
        """A branch that filtered rows itself would own a second cohort filter."""
        source = _branch_source("cohort_variants")
        assert "merge_calls" not in source, (
            "the cohort branch must delegate to stage_cohort_variants.run(); "
            "calling merge_calls here would leave two implementations of the "
            "minor-allele-frequency filter free to disagree"
        )

    def test_variants_still_folds_in_the_regulators_screen(self):
        """spec.md:351 — `regulators` is an internal step of `variants`.

        Built and tested long before the calling existed, and preserved. The
        branch owns it now, so the branch must still write it or the screen
        silently stops running.
        """
        source = _branch_source("variants")
        assert "stage_regulators.run(" in source
        assert "regulators" in source


class TestTheBranchesWriteTheDeclaredOutputs:
    def test_variants_writes_the_contract_table(self):
        assert "table_path(stage_dir, 'variants')" in _branch_source("variants")

    def test_cohort_variants_writes_the_contract_table(self):
        assert (
            "table_path(stage_dir, 'cohort_variants')"
            in _branch_source("cohort_variants")
        )

    def test_variants_writes_the_internal_regulators_table(self):
        assert (
            "internal_table_path(stage_dir, 'regulators')"
            in _branch_source("variants")
        )

    def test_the_columns_written_are_the_declared_contract(self):
        """Not a hand-copied column list.

        The contract in `STAGE_TABLES` was derived from a real bcftools run; a
        literal list in the branch would be a second statement of it, free to
        drift from the validator that checks the file it wrote.
        """
        assert "stage_variants.PER_ISOLATE_COLUMNS" in _branch_source("variants")
        assert (
            "stage_cohort_variants.MERGE_COLUMNS" in _branch_source("cohort_variants")
        )

    @pytest.mark.parametrize("stage", ("variants", "cohort_variants"))
    def test_each_branch_marks_its_own_stage(self, stage):
        assert f"mark('{stage}')" in _branch_source(stage)


class TestTheCohortBranchConsumesThePerIsolateStage:
    def test_it_is_given_the_calls_from_this_run(self):
        """One run's calls, not a re-read that could disagree."""
        assert "variant_calls" in _branch_source("cohort_variants")

    def test_the_per_isolate_branch_publishes_the_calls(self):
        """The per-isolate branch must leave its calls in scope for the merge."""
        assert "variant_calls" in _branch_source("variants")

    def test_the_shared_state_is_declared_nonlocal(self):
        """A closure variable written in one branch must be declared to persist.

        `_execute_stage` is a closure over the run's mutable state and declares
        every one of those names `nonlocal`. Without `variant_calls` in that
        list, assigning it in the `variants` branch raises `UnboundLocalError`
        or silently shadows, and the merge would see an empty cohort.
        """
        tree = ast.parse(RUN_PY.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == "_execute_stage":
                declared = set()
                for inner in ast.walk(node):
                    if isinstance(inner, ast.Nonlocal):
                        declared.update(inner.names)
                assert "variant_calls" in declared
                return
        raise AssertionError("_execute_stage not found")


class TestWhatEachBranchDoesInStub:
    """The decision the ticket asks for, pinned.

    STUB fabricates rather than executes, and `_execute_stage` returns before any
    stage branch on that path. So both branches do nothing in STUB, and the
    tables come from `stub.fabricate` instead. That is the correct division: the
    alternative - running the real aligner in STUB - would defeat the mode,
    which exists to exercise the DAG with no tool installed and no fixture on
    disk.

    The consequence is asserted rather than assumed: `variants` and
    `cohort_variants` are in `STAGE_TABLES`, so `stub.fabricate` writes both
    header-only tables and the DAG resolves without these branches.
    """

    @pytest.mark.parametrize("stage", ("variants", "cohort_variants"))
    def test_stub_fabricates_the_stage_table(self, stage):
        """Declared, so STUB writes a header-only table for it."""
        from papipeline.execution.contracts import required_columns

        assert stage in STAGE_TABLES
        assert required_columns(stage), f"{stage} declares no header to fabricate"

    def test_stub_writes_both_tables_without_the_branches(self, tmp_path):
        """The mechanism, exercised: no aligner, no fixture, both tables exist."""
        from papipeline.stub import fabricate

        written = fabricate("variants", tmp_path)
        assert written["variants"].is_file()
        header = written["variants"].read_text(encoding="utf-8").splitlines()
        assert header and header[0].split("\t") == list(
            STAGE_TABLES["variants"][1]
        )
        # Zero rows: a fabricated row would erase the difference between "no
        # data" and "data that happens to be empty".
        assert len(header) == 1

        cohort = fabricate("cohort_variants", tmp_path)
        assert cohort["cohort_variants"].is_file()

    def test_the_stub_path_returns_before_the_branches(self):
        """Structural: the STUB guard is an unconditional early `return`.

        Located by source order rather than by nesting level, because the stage
        chain lives inside a `try:` and the guard directly in the function body
        - so the two are siblings at different depths, and only their order in
        the file says which runs first.
        """
        source = RUN_PY.read_text(encoding="utf-8")
        guard = source.index("stub.fabricate(stage, stage_dir)")
        first_branch = source.index('elif stage == "variants":')
        assert guard < first_branch, (
            "the STUB guard must precede the stage chain, or STUB would run the "
            "real caller"
        )
        # And the guard returns before the chain, so nothing below it executes.
        between = source[guard:first_branch]
        assert "return" in between, (
            "the STUB block does not return before the stage chain; without the "
            "return, STUB would fall through and align real genomes"
        )


class TestVariantStagesReachableButRealGated:
    """The reconciliation's invariant: enabled in STUB/TEST, still shut for REAL.

    This class previously asserted the opposite - that both stages stayed in
    `UNBUILT_STAGES` - which was true while the branches had no caller. It
    asserted the shape of that intermediate state, and opening the gate
    necessarily falsifies it.

    The distinction it now pins is the one that actually matters, and it is a
    distinction between two *different* gates that are easy to conflate:

    * **the stage gate**, `analysis.<stage>` in ``science.yaml`` plus
      membership of `UNBUILT_STAGES`. Open, or the Snakemake rule executes
      nothing and exits non-zero, failing the DAG. Opening it says "this stage
      exists and is wanted" - nothing about what data it may read;
    * **the mode gate**, `runtime.allow_real_mode` in the machine overlay.
      Still false on all three overlays, so a REAL run of these stages is
      refused before any genome is touched.

    A stage that is enabled is not thereby authorised, and the second gate is
    the one that stops a REAL cohort being analysed on a laptop. So both are
    asserted here, and the second is asserted per-overlay rather than from the
    loaded laptop config alone: a single read would pass while another overlay
    had been flipped.
    """

    @pytest.mark.parametrize("stage", ("variants", "cohort_variants"))
    def test_the_stage_gate_is_open(self, stage):
        """No longer refused as unbuilt."""
        assert stage not in UNBUILT_STAGES, (
            f"{stage} is still refused as unbuilt; the branches have a caller "
            "and the stage has an implementation, so the refusal is now false"
        )

    @pytest.mark.parametrize("stage", ("variants", "cohort_variants"))
    def test_the_stage_is_enabled_in_the_science_config(self, config, stage):
        """A stage absent from `analysis` is disabled, whatever UNBUILT_STAGES says.

        This is not hypothetical: the first attempt at the reconciliation left
        these keys out, the stage was reported "disabled in configuration", and
        the STUB rule failed because `run_stage.py` executed zero stages and
        exited non-zero. The gate has two halves and this is the one that fails
        quietly.
        """
        assert config.stage_enabled(stage) is True, (
            f"analysis.{stage} is not true in science.yaml, so the stage is "
            "disabled and its Snakemake rule will execute nothing and fail"
        )

    def test_stub_fabricates_both_without_reaching_the_caller(self, tmp_path):
        """Reachable in STUB, and still doing no real work.

        Reachability is the new state; fabricating rather than executing is the
        property that must survive it, or STUB would align real genomes.
        """
        from papipeline.stub import fabricate

        for stage in ("variants", "cohort_variants"):
            written = fabricate(stage, tmp_path)
            assert written[stage].is_file()
            assert len(written[stage].read_text(encoding="utf-8").splitlines()) == 1, (
                "STUB writes a header and no rows; a fabricated row would erase "
                "the difference between no data and empty data"
            )

    @pytest.mark.parametrize(
        "overlay", ["laptop.yaml", "bigmachine.yaml", "smoke.yaml"]
    )
    def test_the_mode_gate_is_still_shut_on_every_overlay(self, overlay):
        """REAL stays refused, and this reads all three rather than one.

        Reading the overlay the tests happen to load would pass while another
        had been flipped - and `allow_real_mode` is the single flag standing
        between this and a REAL cohort being analysed on a machine that cannot
        hold it.
        """
        import yaml

        path = REPO / "config" / "machines" / overlay
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert data["runtime"]["allow_real_mode"] is False, (
            f"{overlay} has allow_real_mode true. REAL mode is authorised "
            "separately and in writing; enabling a stage in science.yaml is not "
            "that authorisation."
        )

    def test_real_mode_is_still_refused_through_the_public_interface(self, config):
        """End to end, not just the flag: `run_pipeline(mode='REAL')` still raises.

        The flag being false is the precondition; this is the behaviour, and it
        is the one a caller actually meets.
        """
        from papipeline.errors import ModeNotAllowedError
        from papipeline.run import run_pipeline

        with pytest.raises(ModeNotAllowedError) as excinfo:
            run_pipeline(config=config, mode="REAL")
        message = str(excinfo.value)
        assert "allow_real_mode" in message

    def test_nothing_reads_real_genomes_for_these_stages_in_test_mode(
        self, config, tmp_path, monkeypatch
    ):
        """TEST reaches the branches without the caller or a reference.

        The stages are now enabled, so "enabled" must not have quietly become
        "invokes minimap2". A stage that reached the aligner in TEST would need
        a reference TEST does not have, and would fail confusingly rather than
        obviously.
        """
        from papipeline.models import RunMode, Sample
        from papipeline.manifest import SampleManifest
        from papipeline.stages import variants as stage_variants

        def explode(*args, **kwargs):
            raise AssertionError("TEST mode must not invoke the aligner")

        monkeypatch.setattr(stage_variants, "call_isolate", explode, raising=False)
        calls = stage_variants.run(
            config,
            SampleManifest(samples=[Sample(sample_id="TEST_PA_001")]),
            RunMode.TEST,
            data_root=REPO / "test_data",
            workdir=tmp_path,
        )
        assert calls["TEST_PA_001"]
