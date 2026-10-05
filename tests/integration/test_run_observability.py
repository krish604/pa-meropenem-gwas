"""The production path, observed.

These tests exist to prove one specific claim: that a **normal** pipeline
run — the same ``run_pipeline`` the CLI calls — now feeds the observatory.
Not ``fixture.py``, which drives the runner directly. The production
orchestrator.

So every test here calls ``run_pipeline`` and then reads the execution
store, the event bus and the observatory API. If the wiring regresses, these
fail.

What is also checked, because the integration is only safe if it is
inert when off:

- with the observatory disabled, the pipeline behaves exactly as before and
  writes no execution state;
- no stage runs twice, with or without observation;
- a stage that fails still stops the run, and is never reported SUCCEEDED;
- a stage whose outputs are wrong is INVALID even though the stage itself
  completed without error.
"""

from __future__ import annotations

import json
import logging
import shutil
import sys
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import StageError
from papipeline.execution import StageState
from papipeline.execution.contracts import (
    STAGE_TABLES,
    internal_table_path,
    required_columns,
    table_path,
)
from papipeline.models import RunMode
from papipeline.observatory import read_snapshot
from papipeline.observatory.api import payload
from papipeline.observatory.events import EventBus
from papipeline.observatory.wiring import ObservatorySettings, config_hash, create
from papipeline.run import EXECUTION_ORDER, PREREQUISITES, run_pipeline

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "config" / "science.yaml"


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


@pytest.fixture()
def observed(tmp_path):
    """An enabled observatory handle on a temporary store."""
    handle = create(load_config(CONFIG), db_path=tmp_path / "obs.db",
                    run_key="observed-run", enable=True)
    yield handle
    handle.close()


@pytest.fixture(autouse=True)
def _own_results_tree(tmp_path, monkeypatch):
    """Redirect this file's runs into ``tmp_path`` and nothing else.

    **Why this is here: a neighbour's output file was being read as this run's
    own.** These tests wrote into the repository's real ``results/test`` tree,
    which every other TEST-mode run in the suite also writes. When stage 4
    raises, ``execution/runner.py``'s ``_partial_output`` decides FAILED vs
    INCOMPLETE by asking whether a declared output exists with an mtime inside a
    **one-second** window around the attempt (``began_at = time.time() - 1.0``,
    a deliberate margin for filesystem timestamp resolution). A table another
    test wrote less than a second earlier therefore reads as "this attempt died
    part-way through", and the stage is reported INCOMPLETE.

    That is not hypothetical and it is not new: it was latent at
    ``e4690cb`` and became reachable when FIX6-PROBE removed the eager
    ``--version`` sweep, which had been costing 3.6 s per ``run_pipeline`` call
    and slowing the suite enough to keep the tests more than a second apart.
    The fix for the ordering bug made a pre-existing race reachable. Verified at
    baseline by exercising ``_partial_output`` directly: no table -> FAILED, a
    table 2 s old -> FAILED, a table written moments earlier -> INCOMPLETE.

    So this fixture is the honest fix rather than a widened assertion: a test
    that reads a shared mutable tree cannot make a claim about its own run. It
    is R1/R10-neutral - **no assertion in this file is changed**; each one gets
    a private tree, which is what it was already assuming it had.
    """
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "results"))


def stage_dir() -> Path:
    """This run's stage-table directory, from the CONFIGURED results root.

    Derived from the loaded config rather than spelled as
    ``REPO / "results" / "test" / ...``: the four call sites that need this path
    used to name the repository's real tree by hand, which is precisely what
    `_own_results_tree` redirects away from. A literal here would keep reading
    the shared tree while the runs wrote somewhere else, and the tests that
    corrupt or compare a table would act on a file no run had written.
    """
    config = load_config(CONFIG)
    return Path(config.intermediate_root(RunMode.TEST)) / "stages"


def run(mode="TEST", observatory=None, only=None, skip=None, write_html=False,
        resume=False):
    return run_pipeline(
        config_path=CONFIG, mode=mode, only=only, skip=skip,
        write_html=write_html, config=load_config(CONFIG),
        observatory=observatory, resume=resume,
    )


# ── the pipeline still works with the observatory off ─────────────

def _runnable() -> int:
    """Stages that actually execute in a TEST run.

    The spec's fifteen, minus the three with no implementation (ticket 14),
    which are disabled and unrequested and so are skipped before they are ever
    marked. Derived rather than hardcoded, so this file cannot drift from the
    taxonomy the way a literal 16 already did.
    """
    from papipeline.run import STAGE_ORDER, UNBUILT_STAGES

    return len(STAGE_ORDER) - len(UNBUILT_STAGES)


class TestObservatoryIsOptional:
    def test_a_normal_run_completes_with_the_observatory_off(self, tmp_path):
        """The default path. No observatory, no store, no behaviour change."""
        result = run()
        assert len(result.stage_status) == _runnable()
        assert all(v == "completed" for v in result.stage_status.values())

    def test_no_execution_store_is_written_when_disabled(self, tmp_path):
        db = tmp_path / "should-not-exist.db"
        handle = create(load_config(CONFIG), db_path=db, run_key="x", enable=False)
        assert handle.enabled is False
        assert handle.store is None
        run(observatory=handle)
        assert not db.exists(), "a disabled observatory must not create a store"
        handle.close()

    def test_the_observatory_is_not_a_pipeline_dependency(self, tmp_path):
        """Observability is not required for the science."""
        handle = create(load_config(CONFIG), db_path=tmp_path / "off.db", enable=False)
        result = run(observatory=handle)
        assert result.n_samples == 20
        handle.close()

    def test_configuration_defaults_to_disabled(self):
        settings = ObservatorySettings.from_config(load_config(CONFIG))
        assert settings.enabled is False

    def test_an_explicit_argument_overrides_configuration(self, tmp_path):
        settings = ObservatorySettings.from_config(
            load_config(CONFIG), db_path=tmp_path / "x.db", run_key="r", enable=True)
        assert settings.enabled is True
        assert settings.run_key == "r"


# ── the production path feeds the observatory ──────────────────────

class TestProductionRunIsObserved:
    def test_a_real_run_persists_every_stage(self, observed):
        result = run(observatory=observed)
        assert len(result.stage_status) == _runnable()

        snap = read_snapshot(observed.store, "observed-run")
        assert snap.exists is True
        assert len(snap.tasks) == _runnable()
        assert snap.counts == {StageState.SUCCEEDED.value: _runnable()}

    def test_every_stage_is_present_in_the_snapshot(self, observed):
        run(observatory=observed)
        snap = read_snapshot(observed.store, "observed-run")
        recorded = {t["stage"] for t in snap.tasks}
        # The three unbuilt stages never run in TEST, so they are never
        # observed. Every stage that *does* run must be in the snapshot, or the
        # observatory is silently missing work.
        from papipeline.run import UNBUILT_STAGES

        expected = set(EXECUTION_ORDER) - set(UNBUILT_STAGES)
        assert recorded == expected, (
            f"missing from the snapshot: {sorted(expected - recorded)}; "
            f"unexpected: {sorted(recorded - expected)}"
        )

    def test_each_task_carries_real_provenance(self, observed):
        run(observatory=observed)
        snap = read_snapshot(observed.store, "observed-run")
        for task in snap.tasks:
            assert task["command"], "a task must record what executed it"
            assert task["config_hash"], "a task must record the configuration hash"
            assert task["started_at"] and task["ended_at"]
            assert task["elapsed_seconds"] is not None
            assert task["validation_state"] == StageState.SUCCEEDED.value

    def test_the_events_are_real_transitions(self, observed):
        run(observatory=observed)
        names = [e.name for e in observed.bus.recent()]
        assert names.count("TASK_STARTED") == _runnable()
        assert names.count("TASK_COMPLETED") == _runnable()
        # At least a start and a completion per stage that ran. A literal 32
        # encoded the sixteen-stage taxonomy and drifted the moment three stages
        # became unbuilt.
        assert observed.bus.sequence >= 2 * _runnable()

    def test_the_api_serves_the_real_run(self, observed):
        run(observatory=observed)
        body = payload(observed.store, "observed-run")
        assert len(body["tasks"]) == _runnable()
        assert body["exists"] is True
        assert body["counts"].get("SUCCEEDED") == _runnable()
        # The idle overlay is a statement about *activity*, so a finished run
        # is idle while still being fully reported.
        assert body["idle"] is True

    def test_a_second_run_is_a_rerun_and_never_emits_task_resumed(self, observed):
        """A deliberate second execution is not a resume.

        This is the rule that was previously wrong: every invocation that
        found a prior SUCCEEDED row announced TASK_RESUMED, so an ordinary
        re-run claimed to have resumed work it had in fact redone.
        """
        run(observatory=observed)
        first = observed.bus.sequence
        run(observatory=observed)

        assert observed.classification.kind.value == "RE_RUN"
        names = [e.name for e in observed.bus.recent(first)]
        # No resume event anywhere in the second run.
        assert "TASK_RESUMED" not in names
        assert "RUN_RESUMED" not in names
        assert "RUN_STARTED" in names
        assert observed.bus.sequence > first

    def test_attempt_numbering_continues_across_a_rerun(self, observed):
        run(observatory=observed)
        run(observatory=observed)
        snap = read_snapshot(observed.store, "observed-run")
        assert all(t["attempt"] == 2 for t in snap.tasks)
        history = observed.store.attempts("observed-run", "amr", "")
        assert [a.attempt for a in history] == [1, 2]
        assert [a.state.value for a in history] == ["SUCCEEDED", "SUCCEEDED"]

    def test_a_resume_skips_validated_stages_and_emits_task_resumed(self, observed):
        """The genuine resume case."""
        run(observatory=observed)
        first = observed.bus.sequence
        result = run(observatory=observed, resume=True)

        assert observed.classification.kind.value == "RESUME"
        names = [e.name for e in observed.bus.recent()][first:]
        assert "RUN_RESUMED" in names
        assert names.count("TASK_RESUMED") == _runnable()
        # Nothing was re-executed, because everything still re-validated.
        assert all(v == "skipped_validated"
                   for v in result.stage_status.values())

    def test_downstream_availability_is_announced_from_the_store(self, observed):
        run(observatory=observed)
        unblocked = [e for e in observed.bus.recent() if e.name == "STAGE_UNBLOCKED"]
        stages = {e.stage for e in unblocked}
        # Every stage that declares prerequisites and had them satisfied.
        # Stages that declare prerequisites and had them satisfied. `integration`
        # folded into `reporting`; `cooccurrence` now declares them.
        assert "cooccurrence" in stages
        assert "gwas" in stages
        for event in unblocked:
            assert event.payload["prerequisites"], "a prerequisite list is required"


# ── no stage runs twice ────────────────────────────────────────────

class TestNoDoubleExecution:
    def test_each_stage_body_runs_exactly_once_observed(self, observed, monkeypatch):
        calls = {}
        import papipeline.run as run_module

        modules = {
                "validation": "stage_validation", "annotation": "stage_annotation",
                "mlst": "stage_mlst", "amr": "stage_amr",
                "virulence": "stage_virulence",
                "pangenome": "stage_pangenome", "phylogeny": "stage_phylo",
                "gwas": "stage_gwas", "convergence": "stage_convergence",
                "cooccurrence": "stage_cooccurrence",
                # `integration` folded into `reporting` (spec.md:351), and the
                # `reporting` branch has no run() to wrap - see below.
            "phenotype": "stage_phenotype",
        }
        # reporting has no run(); it is a write_report() stage.
        modules.pop("reporting", None)

        for stage, module_name in modules.items():
            target = getattr(run_module, module_name, None)
            if target is None or not hasattr(target, "run"):
                continue
            original = target.run
            key = stage

            def counted(*args, _o=original, _k=key, **kwargs):
                calls[_k] = calls.get(_k, 0) + 1
                return _o(*args, **kwargs)

            monkeypatch.setattr(target, "run", counted)

        run(observatory=observed)
        assert calls, "no stage was intercepted"
        for stage, count in calls.items():
            assert count == 1, f"{stage} ran {count} times"

    def test_the_store_records_one_task_per_stage(self, observed):
        run(observatory=observed)
        rows = observed.store.list_for_run("observed-run")
        assert len(rows) == _runnable()
        assert len({(r["stage"], r["subject"]) for r in rows}) == _runnable()

    def test_attempt_history_has_one_entry_per_run(self, observed):
        run(observatory=observed)
        for stage in ("amr", "cooccurrence", "gwas"):
            history = observed.store.attempts("observed-run", stage, "")
            assert len(history) == 1, f"{stage} logged {len(history)} attempts"


# ── failure behaviour is unchanged and truthful ────────────────────

class TestFailureBehaviour:
    def test_a_failing_stage_stops_the_run(self, observed, monkeypatch):
        import papipeline.stages.amr as stage_amr

        def explode(*args, **kwargs):
            raise RuntimeError("deliberate test failure")

        monkeypatch.setattr(stage_amr, "run", explode)
        with pytest.raises(StageError) as caught:
            run(observatory=observed)
        assert "amr" in str(caught.value)

    def test_a_failing_stage_is_never_reported_succeeded(self, observed, monkeypatch):
        import papipeline.stages.amr as stage_amr

        def explode(*args, **kwargs):
            raise RuntimeError("deliberate test failure")

        monkeypatch.setattr(stage_amr, "run", explode)
        with pytest.raises(StageError):
            run(observatory=observed)
        snap = read_snapshot(observed.store, "observed-run")
        amr = next(n for n in snap.nodes if n["stage"] == "amr")
        assert amr["state"] != StageState.SUCCEEDED.value
        assert amr["state"] == StageState.FAILED.value
        row = observed.store.get("observed-run", "amr", "")
        assert row["state"] == "FAILED"
        assert "deliberate test failure" in (row["error"] or "")

    def test_stages_before_the_failure_still_succeeded(self, observed, monkeypatch):
        import papipeline.stages.amr as stage_amr

        def explode(*args, **kwargs):
            raise RuntimeError("deliberate test failure")

        monkeypatch.setattr(stage_amr, "run", explode)
        with pytest.raises(StageError):
            run(observatory=observed)
        snap = read_snapshot(observed.store, "observed-run")
        done = {n["stage"]: n["state"] for n in snap.nodes if n["total"]}
        assert done["validation"] == StageState.SUCCEEDED.value
        assert done["mlst"] == StageState.SUCCEEDED.value
        assert done["amr"] == StageState.FAILED.value
        assert "gwas" not in done, "a stage after the failure must not be recorded"

    def test_a_stage_that_writes_nothing_is_incomplete_not_succeeded(self, observed, tmp_path):
        """The stage completes without error but produces no output."""
        import papipeline.stages.amr as stage_amr

        def write_nothing(*args, **kwargs):
            # Empties the table the contract will look for, then returns.
            table_path(stage_dir(), "amr").write_text("")
            return {}

        original = stage_amr.run
        try:
            stage_amr.run = write_nothing
            with pytest.raises(StageError):
                run(observatory=observed)
        finally:
            stage_amr.run = original
        row = observed.store.get("observed-run", "amr", "")
        assert row["state"] in (StageState.INCOMPLETE.value, StageState.INVALID.value)
        assert row["state"] != StageState.SUCCEEDED.value

    def test_a_stage_with_a_wrong_header_is_invalid(self, observed, monkeypatch):
        """A plausible table with the wrong columns is INVALID, not success.

        The corruption is applied to ``write_tsv`` so it lands *after* the
        stage wrote its own table; replacing the stage function would be
        overwritten by the stage's own write.
        """
        import papipeline.run as run_module

        original = run_module.write_tsv
        target = table_path(stage_dir(), "amr")

        def write_then_corrupt(path, rows, columns, **kwargs):
            written = original(path, rows, columns, **kwargs)
            if Path(path) == target:
                Path(path).write_text("wrong\tcolumns\n1\t2\n")
            return written

        monkeypatch.setattr(run_module, "write_tsv", write_then_corrupt)
        with pytest.raises(StageError):
            run(observatory=observed)
        row = observed.store.get("observed-run", "amr", "")
        assert row["state"] == StageState.INVALID.value
        assert row["validation_state"] == StageState.INVALID.value
        assert "columns" in (row["validation_detail"] or "")

    def test_no_retry_is_attempted_for_a_scientific_stage(self, observed, monkeypatch):
        """The pipeline always stopped at the first failure. It still does."""
        import papipeline.stages.amr as stage_amr
        calls = {"n": 0}

        def explode(*args, **kwargs):
            calls["n"] += 1
            raise RuntimeError("no second attempt please")

        monkeypatch.setattr(stage_amr, "run", explode)
        with pytest.raises(StageError):
            run(observatory=observed)
        assert calls["n"] == 1


# ── multi-stage: A completes, downstream becomes available ─────────

class TestMultiStageFlow:
    def test_upstream_completion_precedes_downstream_and_is_recorded(self, observed):
        # `mechanisms` folded into `cooccurrence` (spec.md:351), so the owner
        # is what a `--only` selection names now.
        run(only=["cooccurrence"], observatory=observed)

        snap = read_snapshot(observed.store, "observed-run")
        recorded = {t["stage"]: t for t in snap.tasks}
        # --only expands prerequisites, so the declared dependencies ran.
        from papipeline.run import UNBUILT_STAGES

        for dependency in PREREQUISITES["cooccurrence"]:
            if dependency in UNBUILT_STAGES:
                # Declared, unbuilt, and in TEST not runnable. The edge is real
                # in the DAG (Snakemake satisfies it from a stand-in fixture);
                # at the run level there is nothing to observe.
                assert dependency not in recorded
                continue
            assert dependency in recorded, f"{dependency} did not run"
            assert recorded[dependency]["state"] == StageState.SUCCEEDED.value
        assert recorded["cooccurrence"]["state"] == StageState.SUCCEEDED.value

    def test_the_scientific_result_is_identical_observed_or_not(self, observed, tmp_path):
        """The decisive regression check: observation changes nothing."""
        plain = run()
        observed_run = run(observatory=observed)
        assert plain.stage_status == observed_run.stage_status
        assert plain.n_samples == observed_run.n_samples
        for key, path in plain.outputs.items():
            other = observed_run.outputs[key]
            assert Path(path).read_bytes() == Path(other).read_bytes(), (
                f"{key} differs between an observed and an unobserved run"
            )

    def test_the_master_table_is_byte_identical(self, observed):
        """A stronger check on the stage that consumes everything."""
        run(observatory=observed)
        # `integration` folded into `reporting` (spec.md:351), so the master
        # table is an internal table rather than a stage's principal output.
        observed_table = internal_table_path(stage_dir(), "master_table")
        first = observed_table.read_bytes()
        run()
        assert observed_table.read_bytes() == first


# ── contracts cannot drift from what the pipeline writes ───────────

class TestContractsMatchThePipeline:
    def test_every_declared_file_exists_after_a_run(self, observed):
        run(observatory=observed)
        from papipeline.run import UNBUILT_STAGES

        stage_tables = stage_dir()
        for stage, (filename, columns) in STAGE_TABLES.items():
            if stage in UNBUILT_STAGES:
                # No implementation, so no output. Asserted in STUB instead -
                # see test_stage_taxonomy.py. A 0-row file here would be a
                # claim the stage never made.
                assert not (stage_tables / filename).exists(), (
                    f"{stage} has no implementation but wrote {filename}; an "
                    "empty table is a finding and this stage looked for nothing"
                )
                continue
            path = stage_tables / filename
            assert path.exists(), f"{stage}: {filename} was not written"
            header = path.read_text().splitlines()[0].split("\t")
            assert tuple(header) == columns, (
                f"{stage} header drifted from its contract\n"
                f"  declared {columns}\n  actual   {tuple(header)}"
            )

    def test_table_path_points_at_the_declared_filename(self, tmp_path):
        for stage, (filename, _) in STAGE_TABLES.items():
            assert table_path(tmp_path, stage).name == filename

    def test_every_stage_has_a_contract(self):
        assert set(STAGE_TABLES) == set(EXECUTION_ORDER)

    def test_a_stage_without_a_contract_is_refused(self, tmp_path):
        with pytest.raises(KeyError):
            table_path(tmp_path, "not_a_stage")
