"""Task execution: run, validate, retry, record.

This is the one place that decides a task's state. The order is fixed and
matters:

1. run the command,
2. **validate the declared outputs**,
3. classify the outcome,
4. apply retry policy,
5. persist.

Validation comes before the success decision on purpose. A process that
exits 0 having written the wrong bytes is not a success, and the only way
to know is to look at what it wrote.
"""

from __future__ import annotations

import os
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence

from ..errors import PipelineError
from .retry import AttemptRecord, RetryPolicy
from .state import (
    ExecutionOutcome,
    FailureKind,
    StageState,
    TaskEvent,
    classify_process_failure,
)
from .store import ExecutionStore
from .validation import OutputSpec, ValidationResult, validate

DEFAULT_TIMEOUT_SECONDS = 3600


@dataclass
class TaskContext:
    """Everything about a task except the command to run.

    Attributes:
        run_key: Identifies the cohort/run. State is scoped to it.
        stage: Stage name, matching ``papipeline.run.STAGE_ORDER``.
        subject: Sample id for a per-sample task, ``""`` for cohort-level.
        spec: The declared output contract.
        tool_version: Recorded for provenance; never used for control flow.
        database_version: Recorded for provenance.
        config_hash: Hash of the effective configuration, so a changed config
            invalidates a cached success on resume.
        input_ids: Identifiers of the inputs, recorded for provenance.
        timeout: Wall-clock ceiling for one attempt.
        log_path: Where this task's output is captured.
    """

    run_key: str
    stage: str
    spec: OutputSpec
    subject: str = ""
    tool_version: Optional[str] = None
    database_version: Optional[str] = None
    config_hash: Optional[str] = None
    input_ids: Sequence[str] = ()
    timeout: int = DEFAULT_TIMEOUT_SECONDS
    log_path: Optional[Path] = None

    def __post_init__(self) -> None:
        if self.spec is None:
            raise PipelineError(
                f"Task {self.key} declares no OutputSpec. SUCCEEDED means a "
                "clean exit *and* a passing output contract, so a task without "
                "a contract has no way to report success honestly.",
                stage=self.stage, subject=self.subject, run_key=self.run_key,
            )
        if not self.spec.checks:
            raise PipelineError(
                f"Task {self.key} declares an OutputSpec with no checks. An "
                "empty contract would pass unconditionally, which is the same "
                "as no contract at all.",
                stage=self.stage, subject=self.subject, run_key=self.run_key,
            )

    @property
    def key(self) -> tuple:
        return (self.run_key, self.stage, self.subject)


@dataclass
class TaskResult:
    """The outcome of :func:`run_task`, across all attempts."""

    stage: str
    subject: str
    state: StageState
    attempts: int
    outcome: ExecutionOutcome
    validation: Optional[ValidationResult] = None
    records: List[AttemptRecord] = field(default_factory=list)
    command: tuple = ()
    outputs: tuple = ()

    @property
    def succeeded(self) -> bool:
        return self.state.is_success

    def to_row(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "subject": self.subject,
            "state": self.state.value,
            "attempts": self.attempts,
            "exit_code": self.outcome.exit_code,
            "elapsed_seconds": round(self.outcome.elapsed_seconds, 3),
            "validation_ok": None if self.validation is None else self.validation.ok,
            "detail": self.outcome.detail,
        }


def _partial_output(spec: OutputSpec, since: Optional[float] = None) -> str:
    """Describe declared outputs *this attempt* left behind.

    Returns a short description, or ``""`` when the attempt wrote nothing.
    This is what tells "died part-way through" from "never got started": a
    process that failed *and* left output behind has produced bytes that
    downstream stages would otherwise read as authoritative.

    ``since`` is a wall-clock time and it matters. Without it, a table left
    by a **previous** run of the same task is indistinguishable from one
    this attempt half-wrote, so a stage that raised immediately would be
    reported INCOMPLETE instead of FAILED. Only files this attempt touched
    count.
    """
    present: List[str] = []
    for check in spec.checks:
        if check.path is None or check.kind.name not in ("EXISTS", "NON_EMPTY"):
            continue
        if not check.path.exists():
            continue
        if since is not None:
            try:
                if check.path.stat().st_mtime < since:
                    continue  # stale: left behind by an earlier run
            except OSError:
                continue
        if check.path.stat().st_size == 0:
            present.append(f"{check.path.name} is empty")
        else:
            present.append(f"{check.path.name} ({check.path.stat().st_size} B)")
    return ", ".join(present[:3])


def _call_once(
    fn: Callable[[], Any],
    timeout: int,
    log_path: Optional[Path],
) -> tuple:
    """Run an in-process callable, with the same contract as a subprocess.

    Fifteen of the sixteen pipeline stages are Python functions, not
    external tools, so ``run_task`` has to be able to execute a callable
    and still travel the identical state / validation / retry / event path.
    That is what this does: it is an *invoker*, not a second engine. The
    caller still passes one OutputSpec, still gets one verdict, and the
    store still records the same rows.

    A callable cannot be killed the way a process can, so ``timeout`` is
    recorded but not enforced; a runaway stage is the interpreter's
    problem, as it always was. This is stated rather than pretended.
    """
    started = time.monotonic()
    began_at = time.time() - 1.0  # a small margin for filesystem timestamp resolution
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        result = fn()
    except Exception as exc:  # noqa: BLE001 - classified, never swallowed
        elapsed = time.monotonic() - started
        detail = f"{type(exc).__name__}: {exc}"
        if log_path is not None:
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write(f"\n! {detail}\n")
        return 1, elapsed, detail, began_at
    elapsed = time.monotonic() - started
    detail = ""
    if log_path is not None:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(f"\n+ completed in {elapsed:.3f}s\n")
    return 0, elapsed, detail, began_at


def _run_once(command: Sequence[str], timeout: int, log_path: Optional[Path],
              env: Optional[Mapping[str, str]] = None,
              pid_sink: Optional[Callable[[int], None]] = None):
    """Run one command, capturing stdout+stderr to ``log_path``.

    Returns ``(returncode, elapsed, tail)``. Never raises for a non-zero
    exit; only an inability to start the process is signalled, by raising
    ``OSError`` which the caller classifies.

    ``env`` is layered over the current environment rather than replacing
    it, so a stage can be given a database directory or a tool path without
    silently losing ``PATH``.
    """
    child_env = None
    if env is not None:
        child_env = {**os.environ, **env}
    started = time.monotonic()
    began_at = time.time() - 1.0
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handle = log_path.open("a", encoding="utf-8")
        # The command is recorded, so a log is self-describing: a reader can
        # tell what produced it without cross-referencing the store.
        handle.write(f"\n$ {' '.join(command)}\n")
        handle.flush()
    else:
        handle = None

    # Popen rather than subprocess.run, for exactly one reason: a caller
    # measuring the child needs its pid, and run() does not expose one. The
    # process is still waited on, still timed out, and still captured.
    process = subprocess.Popen(
        list(command),
        stdout=handle if handle is not None else subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=child_env,
    )
    if pid_sink is not None:
        try:
            pid_sink(process.pid)
        except Exception:  # noqa: BLE001 - telemetry cannot break a run
            pass
    timed_out = False
    try:
        stdout, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        process.kill()
        stdout, _ = process.communicate()
    finally:
        if handle is not None:
            handle.close()

    class _Done:
        pass

    proc = _Done()
    proc.returncode = process.returncode
    proc.stdout = stdout.decode("utf-8", "replace") if isinstance(stdout, bytes) else stdout
    if timed_out:
        raise subprocess.TimeoutExpired(list(command), timeout)
    elapsed = time.monotonic() - started
    tail = ""
    if log_path is not None and log_path.exists():
        text = log_path.read_text(encoding="utf-8", errors="replace")
        tail = text[-2000:]
    elif proc.stdout:
        tail = str(proc.stdout)[-2000:]
    return proc.returncode, elapsed, tail, began_at


def run_task(
    ctx: TaskContext,
    command: "Sequence[str] | Callable[[], Any]",
    *,
    store: Optional[ExecutionStore] = None,
    policy: Optional[RetryPolicy] = None,
    sleep: Optional[Callable[[float], None]] = None,
    attempt_hook: Optional[Callable[[int, StageState], None]] = None,
    env: Optional[Mapping[str, str]] = None,
    event_sink: Optional[Callable[["TaskEvent"], None]] = None,
    pid_sink: Optional[Callable[[int], None]] = None,
) -> TaskResult:
    """Execute one task to a terminal state, retrying per policy.

    ``command`` is either an argument vector for an external tool or a
    zero-argument callable for an in-process stage. Both are executed here,
    through this one function, so both get identical state persistence,
    output validation, retry policy and event emission. There is no second
    task engine, and a stage is never executed twice.

    The success decision is made *after* validation, never from the exit
    code alone. Every attempt is written to the store, including the ones
    that failed, so the provenance of a late success is complete.
    """
    def announce(state: StageState, attempt_no: int, detail: str = "") -> None:
        """Publish a transition, if anyone is listening.

        Emission must never affect the run, so a broken subscriber is
        swallowed here rather than propagating into the execution path.
        """
        if event_sink is None:
            return
        event = TaskEvent(
            stage=ctx.stage, subject=ctx.subject, run_key=ctx.run_key,
            state=state, attempt=attempt_no, detail=detail,
            exit_code=outcome.exit_code if outcome else None,
            elapsed_seconds=outcome.elapsed_seconds if outcome else 0.0,
            validation_state=(validation.state.value
                              if validation is not None and state is StageState.SUCCEEDED
                              else None),
            validation_detail=(validation.detail
                               if validation is not None and state is StageState.SUCCEEDED
                               else ""),
        )
        try:
            event_sink(event)
        except Exception:  # noqa: BLE001 - an observer cannot break a run
            pass

    policy = policy or RetryPolicy()
    # A callable is recorded as a single, honest token rather than being
    # stringified into something that looks like a shell command.
    if callable(command):
        command_label: tuple = (f"<in-process> {getattr(command, '__qualname__', 'callable')}",)
    else:
        command_label = tuple(str(c) for c in command)
    records: List[AttemptRecord] = []
    attempt = 0
    outcome: Optional[ExecutionOutcome] = None
    validation: Optional[ValidationResult] = None
    started_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    if store is not None:
        store.record(
            ctx.run_key, ctx.stage, ctx.subject,
            state=StageState.RUNNING, attempt=0,
            max_attempts=policy.max_attempts, started_at=started_at,
            command=list(command_label), tool_version=ctx.tool_version,
            database_version=ctx.database_version, config_hash=ctx.config_hash,
            input_ids=list(ctx.input_ids),
        )

    while True:
        attempt += 1
        if attempt > 1 and store is not None:
            store.set_state(ctx.run_key, ctx.stage, ctx.subject, StageState.RETRYING)
        if attempt_hook is not None:
            attempt_hook(attempt, StageState.RUNNING)
        announce(StageState.RUNNING, attempt)

        try:
            if callable(command):
                returncode, elapsed, tail, began_at = _call_once(
                    command, ctx.timeout, ctx.log_path)
            else:
                returncode, elapsed, tail, began_at = _run_once(
                    command, ctx.timeout, ctx.log_path, env, pid_sink)
        except subprocess.TimeoutExpired:
            outcome = ExecutionOutcome(
                state=StageState.FAILED, failure_kind=FailureKind.TIMEOUT,
                detail=f"exceeded {ctx.timeout}s", command=command_label, elapsed_seconds=ctx.timeout,
            )
        except OSError as exc:
            outcome = ExecutionOutcome(
                state=StageState.FAILED, failure_kind=FailureKind.TOOL_MISSING,
                detail=f"{type(exc).__name__}: {exc}", command=command_label,
            )
        except Exception as exc:  # noqa: BLE001
            outcome = ExecutionOutcome(
                state=StageState.FAILED, failure_kind=FailureKind.EXECUTION,
                detail=f"{type(exc).__name__}: {exc}", command=command_label,
            )
        else:
            if returncode == 0:
                # Exit 0 is necessary, not sufficient. Validate.
                if ctx.spec is not None:
                    validation = validate(ctx.spec)
                    if validation.ok:
                        outcome = ExecutionOutcome(
                            state=StageState.SUCCEEDED, command=command_label,
                            elapsed_seconds=elapsed, exit_code=0,
                            detail=validation.detail,
                            output_paths=tuple(str(p) for p in ctx.spec.paths()),
                        )
                    else:
                        kind = (FailureKind.INCOMPLETE
                                if validation.state is StageState.INCOMPLETE
                                else FailureKind.INVALID)
                        outcome = ExecutionOutcome(
                            state=validation.state, failure_kind=kind,
                            detail=validation.detail, command=command_label,
                            elapsed_seconds=elapsed, exit_code=0,
                        )
                else:
                    outcome = ExecutionOutcome(
                        state=StageState.SUCCEEDED, command=command_label,
                        elapsed_seconds=elapsed, exit_code=0,
                        detail="completed; outputs not yet validated",
                    )
            else:
                outcome = classify_process_failure(returncode, tail, command)
                outcome = ExecutionOutcome(
                    state=outcome.state, failure_kind=outcome.failure_kind,
                    detail=tail or f"exit {returncode}", command=command_label,
                    elapsed_seconds=elapsed, exit_code=returncode,
                )
                # A process that died *with partial output on disk* got far
                # enough to be worth retrying; one that produced nothing at
                # all failed before doing any work, which is usually a
                # deterministic problem. The distinction decides whether a
                # retry is offered, so it is made explicitly rather than
                # inferred from the return code alone.
                partial = _partial_output(ctx.spec, began_at)
                if partial and outcome.failure_kind is not FailureKind.TOOL_MISSING:
                    outcome = ExecutionOutcome(
                        state=StageState.INCOMPLETE,
                        failure_kind=FailureKind.INCOMPLETE,
                        detail=(f"process exited {returncode} with partial output: "
                                f"{partial}; {outcome.detail}")[-600:],
                        command=command_label, elapsed_seconds=elapsed,
                        exit_code=returncode,
                    )

        record = AttemptRecord(
            attempt=attempt, state=outcome.state, failure_kind=outcome.failure_kind,
            detail=outcome.detail, exit_code=outcome.exit_code,
            elapsed_seconds=outcome.elapsed_seconds, command=command_label,
        )
        records.append(record)
        if store is not None:
            record = store.record_attempt(
                ctx.run_key, ctx.stage, ctx.subject, record, ctx.log_path)
            # The store owns attempt numbering, so a resumed task continues
            # its history rather than restarting at 1.
            attempt = record.attempt

        if outcome.succeeded:
            final = StageState.SUCCEEDED
            break

        if policy.should_retry(outcome, attempt):
            delay = policy.sleep_for(attempt)
            if store is not None:
                store.record(
                    ctx.run_key, ctx.stage, ctx.subject,
                    state=StageState.RETRYING, attempt=attempt,
                    max_attempts=policy.max_attempts,
                    command=list(command_label), tool_version=ctx.tool_version,
                    database_version=ctx.database_version, config_hash=ctx.config_hash,
                    input_ids=list(ctx.input_ids), failure_kind=outcome.failure_kind,
                    error=outcome.detail,
                    validation_state=(validation.state.value if validation else None),
                    validation_detail=(validation.detail if validation else None),
                    validation_json=(validation.to_json() if validation else None),
                )
            announce(StageState.RETRYING, attempt,
                     f"retrying after {outcome.failure_kind or outcome.state.value}")
            if sleep is not None:
                sleep(delay)
            else:
                policy.wait(attempt)
            continue

        final = policy.next_state(outcome, attempt)
        break

    ended_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    if store is not None:
        store.record(
            ctx.run_key, ctx.stage, ctx.subject,
            state=final, attempt=attempt, max_attempts=policy.max_attempts,
            started_at=started_at, ended_at=ended_at,
            elapsed_seconds=outcome.elapsed_seconds,
            command=list(command_label), tool_version=ctx.tool_version,
            database_version=ctx.database_version, config_hash=ctx.config_hash,
            input_ids=list(ctx.input_ids),
            output_paths=list(outcome.output_paths) if outcome.output_paths else
                          ([str(p) for p in ctx.spec.paths()] if ctx.spec else []),
            validation_state=(validation.state.value if validation else None),
            validation_detail=(validation.detail if validation else None),
            validation_json=(validation.to_json() if validation else None),
            failure_kind=outcome.failure_kind,
            error=(None if outcome.succeeded else outcome.detail),
        )

    announce(final, attempt, outcome.detail)

    return TaskResult(
        stage=ctx.stage, subject=ctx.subject, state=final, attempts=attempt,
        outcome=outcome, validation=validation, records=records,
        command=command_label,
        outputs=tuple(outcome.output_paths) if outcome.output_paths else
                (tuple(str(p) for p in ctx.spec.paths()) if ctx.spec else ()),
    )
