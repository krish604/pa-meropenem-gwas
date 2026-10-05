"""Controlled retry.

The policy is deliberately boring, because the failure mode it must avoid
is a retry loop that burns a day re-running a command that cannot succeed.

Rules:

* only *retryable* failure kinds are retried;
* ``INVALID`` is never retried, because the tool ran, produced something,
  and that something is wrong - the identical command on identical inputs
  reproduces the identical wrong output;
* backoff grows exponentially and is capped;
* every attempt is recorded, so the provenance of a failure survives the
  retry that eventually succeeded;
* the state after exhaustion is the *cause*, not a generic failure, so an
  operator can tell a timeout from a validation rejection.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional

from .state import ExecutionOutcome, FailureKind, StageState

#: Default ceiling. Three attempts is enough to ride out a transient
#: resource or timeout condition and not enough to hide a real defect.
DEFAULT_MAX_ATTEMPTS = 3
DEFAULT_BASE_DELAY_SECONDS = 2.0
DEFAULT_MAX_DELAY_SECONDS = 60.0


@dataclass(frozen=True)
class RetryPolicy:
    """Bounded retry with exponential backoff.

    Attributes:
        max_attempts: Total attempts including the first. ``1`` disables retry.
        base_delay: First backoff, in seconds. Doubled each retry.
        max_delay: Backoff ceiling.

    There is deliberately no option to retry ``INVALID`` output. Re-running
    an identical command on identical inputs reproduces identical bytes, so
    a retry of an invalid result cannot succeed and only delays the report
    of a real defect. If a validator turns out to be over-strict, the fix is
    to correct the spec, not to re-run the tool.
    """

    max_attempts: int = DEFAULT_MAX_ATTEMPTS
    base_delay: float = DEFAULT_BASE_DELAY_SECONDS
    max_delay: float = DEFAULT_MAX_DELAY_SECONDS

    def __post_init__(self) -> None:
        if self.max_attempts < 1:
            raise ValueError("max_attempts must be >= 1")
        if self.base_delay < 0 or self.max_delay < 0:
            raise ValueError("delays must be non-negative")

    def should_retry(self, outcome: ExecutionOutcome, attempt: int) -> bool:
        """Whether a further attempt is permitted after ``attempt``."""
        if attempt >= self.max_attempts:
            return False
        kind = outcome.failure_kind
        if kind is None:
            return False
        return kind.retryable

    def delay_for(self, attempt: int) -> float:
        """Backoff before attempt ``attempt + 1``. Capped at ``max_delay``."""
        if self.base_delay <= 0:
            return 0.0
        return min(self.base_delay * (2 ** max(0, attempt - 1)), self.max_delay)

    def next_state(self, outcome: ExecutionOutcome, attempt: int) -> StageState:
        """The state to record now.

        While retries remain the task is ``RETRYING``; once exhausted the
        cause state is preserved so the operator sees *what* went wrong.
        """
        if outcome.succeeded:
            return StageState.SUCCEEDED
        if self.should_retry(outcome, attempt):
            return StageState.RETRYING
        return outcome.state

    def sleep_for(self, attempt: int) -> float:
        """Seconds to wait before the next attempt."""
        return self.delay_for(attempt)

    def wait(self, attempt: int) -> None:
        delay = self.sleep_for(attempt)
        if delay > 0:
            time.sleep(delay)


@dataclass(frozen=True)
class AttemptRecord:
    """One attempt, for the audit trail."""

    attempt: int
    state: StageState
    failure_kind: Optional[FailureKind]
    detail: str
    exit_code: Optional[int]
    elapsed_seconds: float
    command: tuple = ()

    def to_row(self) -> dict:
        return {
            "attempt": self.attempt,
            "state": self.state.value,
            "failure_kind": self.failure_kind.value if self.failure_kind else "",
            "detail": self.detail,
            "exit_code": self.exit_code,
            "elapsed_seconds": round(self.elapsed_seconds, 3),
        }
