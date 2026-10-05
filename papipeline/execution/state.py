"""Execution state model.

The rule this module exists to enforce:

    "Files exist" must NEVER, by itself, mean "the stage succeeded".

A stage is ``SUCCEEDED`` only when its process completed *and* its declared
outputs passed declarative validation. Everything else is a distinct,
named state, so that a resume decision can tell the difference between
"never ran", "ran and failed", "ran and was killed", and "ran and produced
something that is not what we asked for".

No scientific algorithm, threshold or interpretation lives here. This module
describes *execution*, not biology.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Optional

from ..errors import PipelineError


class StageState(str, enum.Enum):
    """Every state a task can be in.

    ``INCOMPLETE`` and ``INVALID`` are deliberately distinct from ``FAILED``:

    * ``FAILED``       - the process itself did not succeed (non-zero exit,
                         timeout, missing binary, exception).
    * ``INCOMPLETE``   - the process may have succeeded, but required output
                         is missing or truncated. A run killed at 90% lands
                         here, never in ``SUCCEEDED``.
    * ``INVALID``      - output is present and non-empty, but failed
                         validation: wrong schema, unparseable, degenerate,
                         or scientifically impossible for the request.

    Conflating these is what makes a run untrustworthy: a resume that treats
    any of them as success silently propagates a broken result.
    """

    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCEEDED = "SUCCEEDED"
    FAILED = "FAILED"
    INCOMPLETE = "INCOMPLETE"
    INVALID = "INVALID"
    RETRYING = "RETRYING"

    @property
    def is_terminal(self) -> bool:
        """True when no further automatic transition is expected."""
        return self in _TERMINAL

    @property
    def is_success(self) -> bool:
        """True only for a validated success. File presence is not consulted."""
        return self is StageState.SUCCEEDED

    @property
    def needs_rerun(self) -> bool:
        """True when the task must be executed again to make progress."""
        return self in _NEEDS_RERUN

    @classmethod
    def coerce(cls, value: object) -> "StageState":
        """Parse a persisted state, tolerating case and unknown values."""
        if isinstance(value, cls):
            return value
        text = str(value or "").strip().upper()
        try:
            return cls(text)
        except ValueError:
            raise PipelineError(
                "Unknown execution state", state=str(value), known=[s.value for s in cls]
            )


_TERMINAL: FrozenSet[StageState] = frozenset(
    {StageState.SUCCEEDED, StageState.FAILED, StageState.INCOMPLETE, StageState.INVALID}
)
_NEEDS_RERUN: FrozenSet[StageState] = frozenset(
    {StageState.PENDING, StageState.FAILED, StageState.INCOMPLETE, StageState.INVALID}
)


class FailureKind(str, enum.Enum):
    """Why a task attempt did not succeed.

    The distinction that matters is *transient* versus *deterministic*. A
    timeout or a full disk will succeed on a retry; a validation failure
    will not, and retrying it forever burns hours and hides a real defect.
    """

    #: The process did not succeed: non-zero exit, signal, exception.
    EXECUTION = "execution"
    #: Required output missing, empty, or truncated.
    INCOMPLETE = "incomplete"
    #: Output present but failed declarative validation.
    INVALID = "invalid"
    #: Environment-level problem likely to clear: ENOSPC, EMFILE, EAGAIN.
    RESOURCE = "resource"
    #: Wall-clock timeout.
    TIMEOUT = "timeout"
    #: A required executable is missing or unusable.
    TOOL_MISSING = "tool_missing"

    @property
    def retryable(self) -> bool:
        """Whether re-running the same inputs could plausibly succeed.

        ``INVALID`` is never retryable: the tool ran, produced something,
        and that something is wrong. Re-running the identical command with
        the identical inputs reproduces the identical wrong output.
        """
        return self in _RETRYABLE

    @property
    def state(self) -> StageState:
        """The stage state this failure kind implies."""
        return _KIND_TO_STATE[self]


_RETRYABLE: FrozenSet[FailureKind] = frozenset(
    {FailureKind.INCOMPLETE, FailureKind.RESOURCE, FailureKind.TIMEOUT}
)
_KIND_TO_STATE: Dict[FailureKind, StageState] = {
    FailureKind.EXECUTION: StageState.FAILED,
    FailureKind.INCOMPLETE: StageState.INCOMPLETE,
    FailureKind.INVALID: StageState.INVALID,
    FailureKind.RESOURCE: StageState.FAILED,
    FailureKind.TIMEOUT: StageState.FAILED,
    FailureKind.TOOL_MISSING: StageState.FAILED,
}

#: Exit codes and OS error names that indicate a transient environment
#: problem rather than a deterministic failure. Kept small and explicit so
#: the classification is auditable rather than clever.
TRANSIENT_ERRNO: FrozenSet[int] = frozenset({28, 75})  # ENOSPC, EFBIG
TRANSIENT_OSERROR: FrozenSet[str] = frozenset(
    {"EMFILE", "ENFILE", "EAGAIN", "ETXTBSY", "EINTR"}
)


@dataclass(frozen=True)
class ExecutionOutcome:
    """The result of one attempt, before retry policy is applied.

    Attributes:
        state: The state this attempt produced.
        failure_kind: Why, when the attempt did not succeed.
        detail: Human-readable explanation. Never used for control flow.
        command: The exact argv, for provenance.
        elapsed_seconds: Wall time of the attempt.
        exit_code: Process exit code, or ``None`` if it never ran.
        output_paths: Paths the task declared or produced.
    """

    state: StageState
    failure_kind: Optional[FailureKind] = None
    detail: str = ""
    command: tuple = ()
    elapsed_seconds: float = 0.0
    exit_code: Optional[int] = None
    output_paths: tuple = ()

    @property
    def succeeded(self) -> bool:
        return self.state.is_success

    @property
    def retryable(self) -> bool:
        return bool(self.failure_kind and self.failure_kind.retryable)

    def as_row(self) -> Dict[str, object]:
        return {
            "state": self.state.value,
            "failure_kind": self.failure_kind.value if self.failure_kind else "",
            "detail": self.detail,
            "elapsed_seconds": round(self.elapsed_seconds, 3),
            "exit_code": self.exit_code,
        }


def classify_process_failure(
    exit_code: Optional[int], detail: str = "", command: tuple = ()
) -> ExecutionOutcome:
    """Turn a non-successful process into a classified outcome.

    Raises:
        PipelineError: if ``exit_code`` is 0. Classifying a successful
            process as a failure is a footgun: a caller that used this
            helper to decide success would mark a run broken, and a caller
            that used it to decide failure would skip a run that produced
            perfectly good output. A zero exit is not a failure, and
            success still has to be earned by validating the outputs - so
            this function refuses the input rather than guessing.
    """
    if exit_code == 0:
        raise PipelineError(
            "classify_process_failure called with exit code 0",
            hint="a zero exit is not a failure; validate the outputs instead",
            detail=detail[:200],
        )
    text = (detail or "").lower()
    kind = FailureKind.EXECUTION
    if exit_code == 127 or ("no such file" in text and
                            ("bakta" in text or "amrfinder" in text or
                             "not found" in text)):
        kind = FailureKind.TOOL_MISSING
    else:
        for errno in TRANSIENT_ERRNO:
            if f"errno {errno}" in text or f"error {errno}" in text:
                kind = FailureKind.RESOURCE
                break
        else:
            for name in TRANSIENT_OSERROR:
                if name.lower() in text:
                    kind = FailureKind.RESOURCE
                    break
            else:
                if "timed out" in text or "timeout" in text:
                    kind = FailureKind.TIMEOUT
    return ExecutionOutcome(
        state=kind.state, failure_kind=kind, detail=detail, command=command,
        exit_code=exit_code,
    )


@dataclass(frozen=True)
class TaskEvent:
    """A state transition, announced as it happens.

    The execution layer publishes these so a live view can follow a run
    without polling. They are notifications, not state: the store remains
    authoritative, and a consumer that misses one loses freshness rather
    than correctness.

    This lives here, in the execution layer, so that the runner has no
    dependency on any particular consumer. An observer that wants richer
    event names or a transport adapts this; it does not change it.
    """

    stage: str
    subject: str
    run_key: str
    state: StageState
    attempt: int = 0
    detail: str = ""
    exit_code: Optional[int] = None
    elapsed_seconds: float = 0.0
    #: Set only on the transition that follows a passing contract check.
    validation_state: Optional[str] = None
    validation_detail: str = ""

    @property
    def validated(self) -> bool:
        """True on the transition where outputs were checked and passed."""
        return self.state is StageState.SUCCEEDED and bool(self.validation_state)

    def to_row(self) -> Dict[str, object]:
        return {
            "stage": self.stage, "subject": self.subject,
            "run_key": self.run_key, "state": self.state.value,
            "attempt": self.attempt, "detail": self.detail,
            "exit_code": self.exit_code, "elapsed_seconds": self.elapsed_seconds,
            "validation_state": self.validation_state,
            "validation_detail": self.validation_detail,
        }


@dataclass
class StateTransition:
    """One recorded movement through the state machine, for audit."""

    frm: StageState
    to: StageState
    attempt: int
    at: str
    reason: str = ""
    detail: str = ""
    fields: Dict[str, object] = field(default_factory=dict)
