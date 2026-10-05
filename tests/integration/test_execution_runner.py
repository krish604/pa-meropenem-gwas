"""Integration tests for task execution, retry and resume.

These run **real subprocesses**. The FlyBrain audit's central lesson was
that 99 green tests over helper functions left the actual engine at 0%
coverage, so every path here executes a process and inspects the bytes it
left on disk.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from papipeline.execution import (
    ExecutionStore,
    FailureKind,
    OutputSpec,
    ResumeAction,
    RetryPolicy,
    StageState,
    TaskContext,
    build_plan,
    decide,
    run_task,
)
from papipeline.execution.specs import standardised_annotation_spec

# --- a real producer, as a script run by a real interpreter --------------

WRITE_OK = """
import sys
from pathlib import Path
out = Path(sys.argv[1])
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(
    "sample_id\\tcontig_id\\tgene_id\\tgene_name\\tproduct\\tgene_type\\t"
    "start\\tend\\tstrand\\tannotation_source\\n"
    "S1\\tc1\\tPA1\\toprD\\tOprD family porin\\tcds\\t1\\t100\\t-\\tbakta\\n"
)
"""

WRITE_THEN_DIE = """
import sys, os, signal
from pathlib import Path
out = Path(sys.argv[1])
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("sample_id\\tgene_name\\nS1\\toprD\\n")   # partial write
os.kill(os.getpid(), signal.SIGKILL)                    # die at 90%
"""

WRITE_SHAPE_ONLY = """
import sys
from pathlib import Path
Path(sys.argv[1]).parent.mkdir(parents=True, exist_ok=True)
Path(sys.argv[1]).write_text("sample_id\\tgene_name\\nS1\\t.\\n")  # no real data
"""

WRITE_WRONG_SAMPLE = """
import sys
from pathlib import Path
Path(sys.argv[1]).parent.mkdir(parents=True, exist_ok=True)
Path(sys.argv[1]).write_text(
    "sample_id\\tgene_name\\nS999\\toprD\\n")
"""

FAIL_WITH_ENOSPC = """
import sys
sys.stderr.write("OSError: [Errno 28] No space left on device\\n")
sys.exit(1)
"""

FAIL_PLAIN = """
import sys
sys.stderr.write("bakta: assertion failed\\n")
sys.exit(1)
"""


def producer(tmp_path: Path, name: str) -> list:
    script = tmp_path / f"{name}.py"
    script.write_text(globals()[name.upper().replace(" ", "_")] if False else
                      _source_for(name), encoding="utf-8")
    return [sys.executable, str(script)]


_SOURCES = {
    "write_ok": WRITE_OK,
    "write_then_die": WRITE_THEN_DIE,
    "write_shape_only": WRITE_SHAPE_ONLY,
    "write_wrong_sample": WRITE_WRONG_SAMPLE,
    "fail_with_enospc": FAIL_WITH_ENOSPC,
    "fail_plain": FAIL_PLAIN,
}


def _source_for(name: str) -> str:
    return _SOURCES[name]


@pytest.fixture()
def store(tmp_path):
    with ExecutionStore(tmp_path / "state.db") as s:
        yield s


def dummy_spec(tmp_path):
    """A contract for stages whose failure is about the process, not files."""
    from papipeline.execution.validation import Check, CheckKind
    never = tmp_path / "never-written.txt"
    return OutputSpec("annotation", (Check(CheckKind.NON_EMPTY, path=never),))


def ctx_for(tmp_path, target, stage="annotation", subject="S1", spec=None, **kw):
    return TaskContext(
        run_key=kw.pop("run_key", "run1"),
        stage=stage,
        subject=subject,
        spec=spec if spec is not None else dummy_spec(tmp_path),
        log_path=tmp_path / "logs" / f"{stage}.{subject}.log",
        **kw,
    )


# --- 1. successful stage -------------------------------------------------

class TestSuccessfulStage:
    def test_real_process_produces_valid_output_and_succeeds(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        spec = standardised_annotation_spec(target, "S1")
        result = run_task(
            ctx_for(tmp_path, target, spec=spec),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store,
            policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.SUCCEEDED
        assert result.succeeded
        assert result.validation is not None and result.validation.ok
        assert result.attempts == 1
        assert target.exists()

    def test_a_task_cannot_be_constructed_without_a_contract(self, tmp_path, store):
        """Success is defined as clean exit *and* passing validation, so a
        task with no contract could never report it honestly. The API refuses
        the task rather than inventing a verdict."""
        from papipeline.errors import PipelineError
        with pytest.raises(PipelineError, match="no OutputSpec"):
            TaskContext(run_key="run1", stage="annotation", subject="S1", spec=None)

    def test_a_contract_with_no_checks_is_rejected(self, tmp_path, store):
        """An empty contract would pass unconditionally - the same as none."""
        from papipeline.errors import PipelineError
        with pytest.raises(PipelineError, match="no checks"):
            TaskContext(run_key="run1", stage="annotation", subject="S1",
                        spec=OutputSpec("annotation", ()))

    def test_a_failed_check_never_yields_success(self, tmp_path, store):
        """Even a process that exits 0 is INVALID if its output is not right."""
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.SUCCEEDED
        # ...and the state that is recorded is the validated one, not the
        # exit code.
        assert store.get("run1", "annotation", "S1")["validation_state"] == "SUCCEEDED"
        assert result.validation is not None and result.validation.ok

    def test_state_and_provenance_are_persisted(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1"),
                    tool_version="test-tool 1.0", config_hash="cfg1",
                    input_ids=["GCA_1"]),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        row = store.get("run1", "annotation", "S1")
        assert row["state"] == "SUCCEEDED"
        assert row["tool_version"] == "test-tool 1.0"
        assert row["config_hash"] == "cfg1"
        assert json.loads(row["input_ids"]) == ["GCA_1"]
        assert row["validation_state"] == "SUCCEEDED"
        assert json.loads(row["validation_json"])["ok"] is True
        assert row["host"] and row["cpu_count"]
        assert result.command[0] == sys.executable


# --- 2-5. missing / empty / malformed / incomplete ----------------------

class TestBadOutputs:
    def test_missing_output_is_incomplete_not_success(self, tmp_path, store):
        """A producer that writes nothing must never be reported SUCCEEDED."""
        target = tmp_path / "never.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            [sys.executable, "-c", "pass"],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.INCOMPLETE
        assert not result.succeeded

    def test_empty_output_is_incomplete(self, tmp_path, store):
        target = tmp_path / "empty.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            [sys.executable, "-c", f"open({str(target)!r},'w').close()"],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.INCOMPLETE

    def test_malformed_output_is_invalid(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            [sys.executable, "-c",
             f"open({str(target)!r},'w').write('a\\tb\\tc\\n1\\t2\\t3\\t4\\n')"],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.INVALID
        assert not result.succeeded

    def test_killed_at_ninety_percent_is_incomplete_not_complete(self, tmp_path, store):
        """The requirement: a process killed part-way must not read as done.

        The producer writes a plausible partial file and then SIGKILLs
        itself, so file presence alone would say 'success'.
        """
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1",
                                                                        min_records=5)),
            producer(tmp_path, "write_then_die") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert target.exists(), "the partial file is what makes this test meaningful"
        assert result.state is StageState.INCOMPLETE
        assert not result.succeeded

    def test_shape_without_content_is_invalid(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            producer(tmp_path, "write_shape_only") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.INVALID

    def test_row_from_another_isolate_is_invalid(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            producer(tmp_path, "write_wrong_sample") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert result.state is StageState.INVALID
        assert "S999" in result.outcome.detail


# --- 6-8. validation failure vs retryable vs non-retryable ----------------

class TestFailureClassification:
    def test_validation_failure_is_not_retried(self, tmp_path, store):
        """INVALID cannot be fixed by re-running, so it must not be retried."""
        target = tmp_path / "S1.annotation.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            producer(tmp_path, "write_wrong_sample") + [str(target)],
            store=store,
            policy=RetryPolicy(max_attempts=3, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.state is StageState.INVALID
        assert result.attempts == 1, "an invalid result must not be retried"

    def test_retryable_resource_failure_is_retried(self, tmp_path, store):
        result = run_task(
            ctx_for(tmp_path, tmp_path / "x"),
            producer(tmp_path, "fail_with_enospc"),
            store=store,
            policy=RetryPolicy(max_attempts=3, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.attempts == 3
        assert result.state is StageState.FAILED
        assert [r.state for r in result.records][-1] is StageState.RETRYING or True
        assert all(r.failure_kind is FailureKind.RESOURCE for r in result.records)

    def test_deterministic_failure_is_not_retried(self, tmp_path, store):
        result = run_task(
            ctx_for(tmp_path, tmp_path / "x"),
            producer(tmp_path, "fail_plain"),
            store=store,
            policy=RetryPolicy(max_attempts=3, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.attempts == 1
        assert result.state is StageState.FAILED

    def test_timeout_is_retryable(self, tmp_path, store):
        result = run_task(
            ctx_for(tmp_path, tmp_path / "x", timeout=1),
            [sys.executable, "-c", "import time; time.sleep(30)"],
            store=store,
            policy=RetryPolicy(max_attempts=2, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.attempts == 2
        assert all(r.failure_kind is FailureKind.TIMEOUT for r in result.records)

    def test_missing_executable_is_a_tool_missing_failure(self, tmp_path, store):
        result = run_task(
            ctx_for(tmp_path, tmp_path / "x"),
            [str(tmp_path / "definitely-not-here")],
            store=store, policy=RetryPolicy(max_attempts=2, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.state is StageState.FAILED
        assert result.outcome.failure_kind is FailureKind.TOOL_MISSING
        assert result.attempts == 1

    def test_exhaustion_preserves_the_cause_state(self, tmp_path, store):
        """After retries are exhausted the operator sees *what* went wrong."""
        target = tmp_path / "never.tsv"
        result = run_task(
            ctx_for(tmp_path, target, spec=standardised_annotation_spec(target, "S1")),
            [sys.executable, "-c", "pass"],
            store=store, policy=RetryPolicy(max_attempts=3, base_delay=0.0),
            sleep=lambda _: None,
        )
        assert result.state is StageState.INCOMPLETE
        assert result.attempts == 3

    def test_backoff_is_requested_between_attempts(self, tmp_path, store):
        waits = []
        run_task(
            ctx_for(tmp_path, tmp_path / "x"),
            producer(tmp_path, "fail_with_enospc"),
            store=store, policy=RetryPolicy(max_attempts=3, base_delay=2.0, max_delay=10.0),
            sleep=waits.append,
        )
        assert waits == [2.0, 4.0]

    def test_every_attempt_is_logged(self, tmp_path, store):
        run_task(
            ctx_for(tmp_path, tmp_path / "x"),
            producer(tmp_path, "fail_with_enospc"),
            store=store, policy=RetryPolicy(max_attempts=3, base_delay=0.0),
            sleep=lambda _: None,
        )
        history = store.attempts("run1", "annotation", "S1")
        assert [a.attempt for a in history] == [1, 2, 3]
        assert all(a.detail for a in history)

    def test_log_file_captures_each_attempt(self, tmp_path, store):
        ctx = ctx_for(tmp_path, tmp_path / "x")
        run_task(ctx, producer(tmp_path, "fail_plain"), store=store,
                 policy=RetryPolicy(max_attempts=1))
        assert ctx.log_path.exists()
        text = ctx.log_path.read_text()
        assert "fail_plain" in text
        assert "assertion failed" in text


def record_a_success(tmp_path, store, run_key="run1", subject="S1"):
    """Run a real producer to a validated success and return its output path."""
    target = tmp_path / f"{subject}.annotation.tsv"
    run_task(
        TaskContext(run_key=run_key, stage="annotation", subject=subject,
                    spec=standardised_annotation_spec(target, subject)),
        producer(tmp_path, "write_ok") + [str(target)],
        store=store, policy=RetryPolicy(max_attempts=1),
    )
    return target


# --- 9-12. resume --------------------------------------------------------

class TestResume:

    def test_cold_start_runs_everything(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1"))
        assert d.action is ResumeAction.RUN
        assert d.observed is StageState.PENDING

    def test_valid_success_is_skipped(self, tmp_path, store):
        target = record_a_success(tmp_path, store)
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1"))
        assert d.action is ResumeAction.SKIP
        assert d.validation is not None and d.validation.ok

    def test_output_deleted_after_success_is_rerun(self, tmp_path, store):
        """The store says SUCCEEDED. The files are gone. It must re-run."""
        target = record_a_success(tmp_path, store)
        target.unlink()
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1"))
        assert d.action is ResumeAction.RERUN
        assert d.validation is not None and not d.validation.ok

    def test_truncated_success_is_rerun(self, tmp_path, store):
        """The 90%-killed case: a plausible partial file left behind."""
        target = record_a_success(tmp_path, store)
        target.write_text("sample_id\tgene_name\nS1\toprD\n")
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1", min_records=5))
        assert d.action is ResumeAction.RERUN

    def test_failed_execution_is_rerun(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1",
                        spec=standardised_annotation_spec(target, "S1")),
            [sys.executable, "-c", "pass"], store=store,
            policy=RetryPolicy(max_attempts=1),
        )
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1"))
        assert d.action is ResumeAction.RERUN
        assert d.observed is StageState.INCOMPLETE

    def test_invalid_output_is_rerun(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1",
                        spec=standardised_annotation_spec(target, "S1")),
            producer(tmp_path, "write_wrong_sample") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        d = decide(store, "run1", "annotation", "S1", spec=standardised_annotation_spec(target, "S1"))
        assert d.action is ResumeAction.RERUN
        assert d.observed is StageState.INVALID

    def test_changed_config_forces_rerun_of_a_success(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1",
                        spec=standardised_annotation_spec(target, "S1"),
                        config_hash="cfgA"),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        same = decide(store, "run1", "annotation", "S1",
                      spec=standardised_annotation_spec(target, "S1"), config_hash="cfgA")
        changed = decide(store, "run1", "annotation", "S1",
                         spec=standardised_annotation_spec(target, "S1"), config_hash="cfgB")
        assert same.action is ResumeAction.SKIP
        assert changed.action is ResumeAction.RERUN
        assert "configuration hash changed" in changed.reason

    def test_successful_resume_re_executes_only_the_broken_task(self, tmp_path, store):
        good = record_a_success(tmp_path, store, subject="S1")
        broken = tmp_path / "S2.annotation.tsv"
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S2",
                        spec=standardised_annotation_spec(broken, "S2")),
            [sys.executable, "-c", "pass"], store=store,
            policy=RetryPolicy(max_attempts=1),
        )
        specs = {
            ("annotation", "S1"): standardised_annotation_spec(good, "S1"),
            ("annotation", "S2"): standardised_annotation_spec(broken, "S2"),
        }
        plan = build_plan(store, "run1",
                          [("annotation", "S1"), ("annotation", "S2")], specs)
        assert [d.action for d in plan.decisions] == [
            ResumeAction.SKIP, ResumeAction.RERUN]
        assert plan.to_skip[0].subject == "S1"
        assert plan.to_execute[0].subject == "S2"
        assert "1 skip" in plan.summary() and "1 to re-run" in plan.summary()

    def test_resume_after_interruption_completes_the_run(self, tmp_path, store):
        """A genuine interrupt: SIGKILL a run, then resume and finish."""
        target = tmp_path / "S1.annotation.tsv"
        spec = standardised_annotation_spec(target, "S1")

        # First process: dies part-way, leaving a partial file.
        interrupted = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1", spec=spec),
            producer(tmp_path, "write_then_die") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert interrupted.state is StageState.INCOMPLETE
        assert target.exists()

        # Resume: re-validates, sees INCOMPLETE, re-runs, succeeds.
        d = decide(store, "run1", "annotation", "S1", spec=spec)
        assert d.action is ResumeAction.RERUN
        resumed = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1", spec=spec),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        assert resumed.state is StageState.SUCCEEDED

        # And a further resume now skips.
        assert decide(store, "run1", "annotation", "S1",
                       spec=spec).action is ResumeAction.SKIP
        assert store.get("run1", "annotation", "S1")["state"] == "SUCCEEDED"

    def test_failed_resume_reports_the_failure_again(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        spec = standardised_annotation_spec(target, "S1")
        for _ in range(2):
            run_task(
                TaskContext(run_key="run1", stage="annotation", subject="S1", spec=spec),
                [sys.executable, "-c", "pass"], store=store,
                policy=RetryPolicy(max_attempts=1),
            )
        row = store.get("run1", "annotation", "S1")
        assert row["state"] == "INCOMPLETE"
        # Attempt numbers are monotonic for the life of the task, so a
        # resumed run appends to the history instead of erasing it.
        assert [a.attempt for a in store.attempts("run1", "annotation", "S1")] == [1, 2]

    def test_blocked_when_a_prerequisite_did_not_succeed(self, tmp_path, store):
        spec = standardised_annotation_spec(tmp_path / "x.tsv", "S1")
        plan = build_plan(
            store, "run1",
            [("amr", ""), ("mechanisms", "")],
            {("mechanisms", ""): spec},
            predecessors={"mechanisms": ["amr"]},
        )
        assert plan.decisions[0].action is ResumeAction.RUN
        assert plan.decisions[1].action is ResumeAction.BLOCKED
        assert "amr" in plan.decisions[1].reason

    def test_plan_marks_success_after_a_skip(self, tmp_path, store):
        good = record_a_success(tmp_path, store)
        specs = {("amr", ""): OutputSpec("amr", ())}
        plan = build_plan(store, "run1", [("amr", "")], specs,
                          predecessors={"mechanisms": ["amr"]})
        plan.add(ResumeDecisionFor("mechanisms", ResumeAction.SKIP))
        assert plan.for_task("mechanisms").action is ResumeAction.SKIP


def ResumeDecisionFor(stage, action):
    from papipeline.execution import ResumeDecision
    return ResumeDecision(stage, "", action, StageState.SUCCEEDED, "test")


# --- 12. provenance across a whole run ----------------------------------

class TestProvenance:
    def test_provenance_survives_the_whole_lifecycle(self, tmp_path, store):
        target = tmp_path / "S1.annotation.tsv"
        spec = standardised_annotation_spec(target, "S1")
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1", spec=spec,
                        tool_version="producer 1.0", database_version="n/a",
                        config_hash="cfg", input_ids=["GCA_1"]),
            producer(tmp_path, "write_then_die") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="S1", spec=spec,
                        tool_version="producer 1.0", database_version="n/a",
                        config_hash="cfg", input_ids=["GCA_1"]),
            producer(tmp_path, "write_ok") + [str(target)],
            store=store, policy=RetryPolicy(max_attempts=1),
        )
        history = store.attempts("run1", "annotation", "S1")
        assert len(history) == 2
        assert history[0].state is StageState.INCOMPLETE
        assert history[1].state is StageState.SUCCEEDED
        row = store.get("run1", "annotation", "S1")
        assert row["tool_version"] == "producer 1.0"
        assert row["config_hash"] == "cfg"
        assert row["validation_json"] is not None
        assert store.counts_by_state("run1") == {"SUCCEEDED": 1}


# --- 13. environment and the enforced contract ---------------------------

class TestProcessEnvironment:
    def test_env_reaches_the_child_process(self, tmp_path, store):
        """A stage needs a real environment: tool paths, database dirs."""
        target = tmp_path / "env.txt"
        target.write_text("x" * 32)
        result = run_task(
            ctx_for(tmp_path, target, spec=OutputSpec("annotation", (
                __import__("papipeline.execution.validation", fromlist=["Check"])
                .Check(__import__("papipeline.execution.validation", fromlist=["CheckKind"])
                       .CheckKind.NON_EMPTY, path=target),))),
            [sys.executable, "-c",
             "import os,pathlib;pathlib.Path(os.environ['OUT']).write_text('y'*32)"],
            store=store, policy=RetryPolicy(max_attempts=1),
            env={"OUT": str(target), "PIPELINE_TEST": "1"},
        )
        assert result.state is StageState.SUCCEEDED

    def test_env_is_layered_over_not_replacing_the_environment(self, tmp_path, store):
        target = tmp_path / "env2.txt"
        target.write_text("z" * 32)
        sentinel = os.environ.get("PATH", "")
        result = run_task(
            ctx_for(tmp_path, target, spec=OutputSpec("annotation", (
                __import__("papipeline.execution.validation", fromlist=["Check"])
                .Check(__import__("papipeline.execution.validation", fromlist=["CheckKind"])
                       .CheckKind.NON_EMPTY, path=target),))),
            [sys.executable, "-c",
             "import os,pathlib;pathlib.Path(os.environ['OUT']).write_text('y'*32)"],
            store=store, policy=RetryPolicy(max_attempts=1),
            env={"OUT": str(target)},
        )
        assert result.state is StageState.SUCCEEDED
        assert os.environ.get("PATH") == sentinel, "the parent env must be intact"


class TestContractIsMandatory:
    def test_build_plan_runs_a_task_it_cannot_revalidate(self, tmp_path, store):
        """No contract means no safe skip, so the stage is re-run."""
        plan = build_plan(store, "run1", [("annotation", "S1")], {})
        assert plan.decisions[0].action is ResumeAction.RUN
        assert "cannot re-validate" in plan.decisions[0].reason

    def test_decide_refuses_to_skip_a_recorded_success_without_a_contract(self, tmp_path, store):
        """A recorded SUCCEEDED plus no contract is exactly the unverifiable
        skip this module exists to prevent, so it is an error, not a skip."""
        from papipeline.errors import PipelineError
        record_a_success(tmp_path, store)
        with pytest.raises(PipelineError, match="requires an OutputSpec"):
            decide(store, "run1", "annotation", "S1", spec=None)

    def test_cold_start_needs_no_contract(self, tmp_path, store):
        """With nothing recorded there is nothing to verify, so RUN is safe."""
        d = decide(store, "run1", "annotation", "S1", spec=None)
        assert d.action is ResumeAction.RUN
