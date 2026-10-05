"""Unit tests for the execution state model, retry policy and state store."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.errors import PipelineError
from papipeline.execution import (
    ExecutionOutcome,
    ExecutionStore,
    FailureKind,
    RetryPolicy,
    StageState,
    classify_process_failure,
)
from papipeline.execution.retry import AttemptRecord


# --- state model ---------------------------------------------------------


class TestStageState:
    def test_all_seven_states_exist(self):
        assert {s.value for s in StageState} == {
            "PENDING", "RUNNING", "SUCCEEDED", "FAILED",
            "INCOMPLETE", "INVALID", "RETRYING",
        }

    def test_only_succeeded_counts_as_success(self):
        """'Files exist' must never mean success."""
        assert StageState.SUCCEEDED.is_success
        for state in StageState:
            if state is not StageState.SUCCEEDED:
                assert not state.is_success, f"{state} must not report success"

    def test_failed_incomplete_invalid_are_distinct(self):
        assert len({StageState.FAILED, StageState.INCOMPLETE,
                    StageState.INVALID}) == 3

    def test_three_failure_states_need_rerun(self):
        for state in (StageState.FAILED, StageState.INCOMPLETE, StageState.INVALID):
            assert state.needs_rerun
        assert not StageState.SUCCEEDED.needs_rerun

    def test_running_and_retrying_are_not_terminal(self):
        assert not StageState.RUNNING.is_terminal
        assert not StageState.RETRYING.is_terminal

    def test_coerce_tolerates_case_and_whitespace(self):
        assert StageState.coerce("succeeded") is StageState.SUCCEEDED
        assert StageState.coerce("  INVALID  ") is StageState.INVALID

    def test_coerce_rejects_unknown_state(self):
        with pytest.raises(PipelineError):
            StageState.coerce("BANANA")


class TestFailureKind:
    def test_invalid_is_never_retryable(self):
        """Re-running an identical command reproduces identical wrong bytes."""
        assert not FailureKind.INVALID.retryable

    @pytest.mark.parametrize("kind", [FailureKind.INCOMPLETE, FailureKind.RESOURCE,
                                      FailureKind.TIMEOUT])
    def test_transient_kinds_are_retryable(self, kind):
        assert kind.retryable

    @pytest.mark.parametrize("kind", [FailureKind.EXECUTION, FailureKind.TOOL_MISSING])
    def test_deterministic_kinds_are_not_retryable(self, kind):
        assert not kind.retryable

    def test_kind_maps_to_the_right_state(self):
        assert FailureKind.INVALID.state is StageState.INVALID
        assert FailureKind.INCOMPLETE.state is StageState.INCOMPLETE
        assert FailureKind.EXECUTION.state is StageState.FAILED
        assert FailureKind.TIMEOUT.state is StageState.FAILED


class TestClassifyProcessFailure:
    def test_plain_nonzero_exit_is_a_deterministic_execution_failure(self):
        outcome = classify_process_failure(1, "bakta: some assertion failed")
        assert outcome.state is StageState.FAILED
        assert outcome.failure_kind is FailureKind.EXECUTION
        assert not outcome.retryable

    def test_disk_full_is_retryable(self):
        outcome = classify_process_failure(1, "OSError: [Errno 28] No space left on device")
        assert outcome.failure_kind is FailureKind.RESOURCE
        assert outcome.retryable

    def test_too_many_open_files_is_retryable(self):
        outcome = classify_process_failure(1, "OSError: [Errno 24] Too many open files EMFILE")
        assert outcome.failure_kind is FailureKind.RESOURCE

    def test_timeout_is_retryable(self):
        outcome = classify_process_failure(None, "subprocess timed out after 3600s")
        assert outcome.failure_kind is FailureKind.TIMEOUT
        assert outcome.retryable

    @pytest.mark.parametrize("code,detail", [
        (127, "bakta: command not found"),
        (127, ""),
        (1, "bakta: no such file or directory"),
    ])
    def test_missing_tool_is_deterministic(self, code, detail):
        outcome = classify_process_failure(code, detail)
        assert outcome.failure_kind is FailureKind.TOOL_MISSING
        assert not outcome.retryable

    def test_refuses_to_classify_a_zero_exit(self):
        """A zero exit is not a failure, and this helper must not pretend.

        Silently returning FAILED would let a future caller mark a
        successful run broken; silently returning SUCCEEDED would let one
        skip a run whose output was never validated.
        """
        with pytest.raises(PipelineError) as excinfo:
            classify_process_failure(0, "all good")
        assert "exit code 0" in str(excinfo.value)


# --- retry policy --------------------------------------------------------


class TestRetryPolicy:
    def test_max_attempts_of_one_disables_retry(self):
        policy = RetryPolicy(max_attempts=1)
        failed = ExecutionOutcome(state=StageState.FAILED,
                                   failure_kind=FailureKind.TIMEOUT)
        assert not policy.should_retry(failed, attempt=1)

    def test_transient_failure_is_retried_up_to_the_limit(self):
        policy = RetryPolicy(max_attempts=3)
        failed = ExecutionOutcome(state=StageState.FAILED,
                                   failure_kind=FailureKind.TIMEOUT)
        assert policy.should_retry(failed, attempt=1)
        assert policy.should_retry(failed, attempt=2)
        assert not policy.should_retry(failed, attempt=3)

    def test_invalid_is_never_retried(self):
        policy = RetryPolicy(max_attempts=5)
        invalid = ExecutionOutcome(state=StageState.INVALID,
                                   failure_kind=FailureKind.INVALID)
        assert not policy.should_retry(invalid, attempt=1)

    def test_there_is_no_way_to_opt_into_retrying_invalid(self):
        """The brief forbids retrying a scientifically invalid result.

        An escape hatch would let exactly that happen, and it would be
        delivered by a flag whose default nobody reads. The fix for an
        over-strict validator is to correct the spec.
        """
        with pytest.raises(TypeError):
            RetryPolicy(retry_invalid=True)  # type: ignore[call-arg]

    def test_backoff_grows_exponentially_and_is_capped(self):
        policy = RetryPolicy(base_delay=2.0, max_delay=10.0)
        assert policy.delay_for(1) == 2.0
        assert policy.delay_for(2) == 4.0
        assert policy.delay_for(3) == 8.0
        assert policy.delay_for(4) == 10.0  # capped
        assert policy.delay_for(9) == 10.0

    def test_zero_base_delay_disables_waiting(self):
        assert RetryPolicy(base_delay=0.0).delay_for(1) == 0.0

    def test_next_state_is_retrying_while_attempts_remain(self):
        policy = RetryPolicy(max_attempts=3)
        transient = ExecutionOutcome(state=StageState.FAILED,
                                    failure_kind=FailureKind.TIMEOUT)
        assert policy.next_state(transient, attempt=1) is StageState.RETRYING

    def test_next_state_preserves_the_cause_once_exhausted(self):
        policy = RetryPolicy(max_attempts=2)
        incomplete = ExecutionOutcome(state=StageState.INCOMPLETE,
                                      failure_kind=FailureKind.INCOMPLETE)
        assert policy.next_state(incomplete, attempt=2) is StageState.INCOMPLETE

    def test_success_never_retries(self):
        ok = ExecutionOutcome(state=StageState.SUCCEEDED)
        assert not RetryPolicy(max_attempts=3).should_retry(ok, attempt=1)

    @pytest.mark.parametrize("bad", [{"max_attempts": 0}, {"base_delay": -1}])
    def test_invalid_policy_is_rejected(self, bad):
        with pytest.raises(ValueError):
            RetryPolicy(**bad)


# --- state store ---------------------------------------------------------


@pytest.fixture()
def store(tmp_path):
    with ExecutionStore(tmp_path / "state.db") as s:
        yield s


class TestExecutionStore:
    def test_records_every_required_field(self, store):
        """The phase brief lists these fields explicitly."""
        out = Path("/tmp/out.tsv")
        store.record(
            "run1", "annotation", "S1",
            state=StageState.SUCCEEDED, attempt=2,
            max_attempts=3,
            started_at="2026-09-27T10:00:00Z", ended_at="2026-09-27T10:06:00Z",
            elapsed_seconds=360.0,
            command=["bakta", "g.fna", "--out", "/tmp/out"],
            tool_version="1.12.1", database_version="v6.0-light",
            config_hash="abc123", input_ids=["GCA_1"],
            output_paths=[out],
            validation_state="SUCCEEDED", validation_detail="all checks passed",
            validation_json=json.dumps({"ok": True}),
        )
        row = store.get("run1", "annotation", "S1")
        assert row is not None
        assert row["stage"] == "annotation"
        assert row["subject"] == "S1"
        assert row["state"] == "SUCCEEDED"
        assert row["attempt"] == 2
        assert row["started_at"] == "2026-09-27T10:00:00Z"
        assert row["ended_at"] == "2026-09-27T10:06:00Z"
        assert row["elapsed_seconds"] == 360.0
        assert json.loads(row["command"])[0] == "bakta"
        assert row["tool_version"] == "1.12.1"
        assert row["database_version"] == "v6.0-light"
        assert row["config_hash"] == "abc123"
        assert json.loads(row["input_ids"]) == ["GCA_1"]
        assert json.loads(row["output_paths"]) == [str(out)]
        assert row["validation_state"] == "SUCCEEDED"
        assert row["validation_detail"] == "all checks passed"
        assert json.loads(row["validation_json"])["ok"] is True

    def test_records_host_and_resource_information(self, store):
        store.record("run1", "validation", "S1", state=StageState.SUCCEEDED)
        row = store.get("run1", "validation", "S1")
        assert row["host"]
        assert row["cpu_count"] is not None
        assert row["updated_at"]

    def test_records_error_information(self, store):
        store.record("run1", "amr", "S2", state=StageState.FAILED,
                     failure_kind=FailureKind.EXECUTION, error="exit 1: boom")
        row = store.get("run1", "amr", "S2")
        assert row["failure_kind"] == "execution"
        assert row["error"] == "exit 1: boom"

    def test_second_record_updates_in_place(self, store):
        store.record("run1", "mlst", "S1", state=StageState.RUNNING, attempt=1)
        store.record("run1", "mlst", "S1", state=StageState.SUCCEEDED, attempt=1)
        rows = store.list_for_run("run1")
        assert len(rows) == 1
        assert rows[0]["state"] == "SUCCEEDED"

    def test_subjects_are_independent_tasks(self, store):
        store.record("run1", "annotation", "S1", state=StageState.SUCCEEDED)
        store.record("run1", "annotation", "S2", state=StageState.FAILED)
        assert store.counts_by_state("run1") == {"SUCCEEDED": 1, "FAILED": 1}

    def test_cohort_level_task_uses_empty_subject(self, store):
        store.record("run1", "pangenome", "", state=StageState.SUCCEEDED)
        assert store.get("run1", "pangenome") is not None

    def test_run_keys_are_isolated(self, store):
        store.record("runA", "amr", "S1", state=StageState.SUCCEEDED)
        store.record("runB", "amr", "S1", state=StageState.FAILED)
        assert store.get("runA", "amr", "S1")["state"] == "SUCCEEDED"
        assert store.get("runB", "amr", "S1")["state"] == "FAILED"

    def test_attempt_history_accumulates_across_retries(self, store):
        """Provenance of a late success must include the earlier failures."""
        store.record("run1", "annotation", "S1", state=StageState.RUNNING)
        for attempt, state, kind in (
            (1, StageState.FAILED, FailureKind.TIMEOUT),
            (2, StageState.INCOMPLETE, FailureKind.INCOMPLETE),
            (3, StageState.SUCCEEDED, None),
        ):
            store.record_attempt("run1", "annotation", "S1", AttemptRecord(
                attempt=attempt, state=state, failure_kind=kind,
                detail=f"attempt {attempt}", exit_code=None, elapsed_seconds=1.0,
            ))
        history = store.attempts("run1", "annotation", "S1")
        assert [a.attempt for a in history] == [1, 2, 3]
        assert [a.state for a in history] == [
            StageState.FAILED, StageState.INCOMPLETE, StageState.SUCCEEDED,
        ]

    def test_attempt_requires_a_recorded_task(self, store):
        with pytest.raises(PipelineError):
            store.record_attempt("run1", "ghost", "S1", AttemptRecord(
                attempt=1, state=StageState.FAILED, failure_kind=None,
                detail="", exit_code=1, elapsed_seconds=0.0,
            ))

    def test_set_state_moves_state_without_clobbering_provenance(self, store):
        store.record("run1", "gwas", "", state=StageState.RUNNING,
                     tool_version="pyseer 1.1.2", config_hash="h1")
        store.set_state("run1", "gwas", "", StageState.RETRYING,
                        failure_kind=FailureKind.TIMEOUT, error="slow")
        row = store.get("run1", "gwas", "")
        assert row["state"] == "RETRYING"
        assert row["tool_version"] == "pyseer 1.1.2"  # preserved
        assert row["config_hash"] == "h1"  # preserved

    def test_state_survives_reopening_the_database(self, tmp_path):
        path = tmp_path / "state.db"
        with ExecutionStore(path) as s:
            s.record("run1", "amr", "S1", state=StageState.SUCCEEDED,
                     tool_version="4.2.7")
        with ExecutionStore(path) as s:
            row = s.get("run1", "amr", "S1")
            assert row is not None
            assert row["state"] == "SUCCEEDED"
            assert row["tool_version"] == "4.2.7"

    def test_missing_task_reads_as_none(self, store):
        assert store.get("nope", "nothing", "") is None
