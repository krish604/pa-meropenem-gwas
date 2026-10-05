"""What kind of execution is this?

The pipeline has four distinct situations that were previously conflated
under the word "re-run", and conflating them is how a resume ends up
silently redoing work, or a fresh run claiming to have resumed something.

Run identity is the discriminator, and it is the *run key*:

============  ==========================================================
Situation     Definition
============  ==========================================================
NEW RUN       The run key has no recorded state. Nothing is reused.
RE-RUN        The run key has state, but the caller is **not** resuming:
              everything executes again from scratch. A deliberate new
              execution, not a continuation.
RESUME        The run key has state *and* the caller asked to resume it.
              Tasks whose outputs re-validate are SKIPped; everything else
              executes. The run continues.
RETRY         One task, second attempt, within a single invocation. Not a
              run-level concept at all - it is an attempt-level one.
============  ==========================================================

The rule that matters most, because it is easy to get backwards:

    ``TASK_RESUMED`` is emitted only when a run is genuinely RESUMED.

A re-run of the same key executes every stage and emits ``RUN_STARTED``,
never ``TASK_RESUMED``. A first execution of a key emits ``RUN_STARTED``
too. Only an explicit resume emits ``RUN_RESUMED``, and only the tasks it
actually skips emit ``TASK_RESUMED``.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from .state import StageState
from .store import ExecutionStore


class RunKind(str, enum.Enum):
    """How this invocation relates to what came before."""

    NEW_RUN = "NEW_RUN"
    RE_RUN = "RE_RUN"
    RESUME = "RESUME"


#: What a caller is asking for. Deliberately explicit: the default is to
#: execute everything, because skipping work on a hunch is far worse than
#: redoing it.
class Intent(str, enum.Enum):
    EXECUTE = "EXECUTE"
    RESUME = "RESUME"


@dataclass
class RunClassification:
    """The verdict, plus the evidence for it.

    The evidence is carried rather than discarded so the observatory and the
    logs can say *why* this was a resume, and so a test can assert the
    reason rather than only the label.
    """

    kind: RunKind
    run_key: str
    prior_rows: int = 0
    prior_succeeded: int = 0
    prior_incomplete: int = 0
    prior_failed: int = 0
    prior_invalid: int = 0
    reason: str = ""

    def to_row(self) -> Dict[str, Any]:
        return {
            "kind": self.kind.value,
            "run_key": self.run_key,
            "prior_rows": self.prior_rows,
            "prior_succeeded": self.prior_succeeded,
            "prior_incomplete": self.prior_incomplete,
            "prior_failed": self.prior_failed,
            "prior_invalid": self.prior_invalid,
            "reason": self.reason,
        }


def classify(
    store: Optional[ExecutionStore],
    run_key: str,
    *,
    intent: Intent = Intent.EXECUTE,
) -> RunClassification:
    """Decide what kind of execution this is, from the store alone.

    ``store is None`` means no state is being kept, so there is nothing to
    reuse and nothing to resume: that is a NEW RUN by definition, even if
    the caller asked to resume.
    """
    if store is None:
        return RunClassification(
            kind=RunKind.NEW_RUN, run_key=run_key, prior_rows=0,
            reason="no execution store is open, so no state can be reused",
        )

    rows = store.list_for_run(run_key)
    if not rows:
        return RunClassification(
            kind=RunKind.NEW_RUN, run_key=run_key, prior_rows=0,
            reason=f"run key {run_key!r} has no recorded state",
        )

    tally = {s.value: 0 for s in StageState}
    for row in rows:
        state = str(row.get("state") or "")
        if state in tally:
            tally[state] += 1

    verdict = RunClassification(
        kind=RunKind.NEW_RUN, run_key=run_key, prior_rows=len(rows),
        prior_succeeded=tally[StageState.SUCCEEDED.value],
        prior_incomplete=tally[StageState.INCOMPLETE.value],
        prior_failed=tally[StageState.FAILED.value],
        prior_invalid=tally[StageState.INVALID.value],
    )

    if intent is Intent.RESUME:
        verdict.kind = RunKind.RESUME
        unfinished = (
            verdict.prior_incomplete + verdict.prior_failed + verdict.prior_invalid
        )
        verdict.reason = (
            f"resume requested; {verdict.prior_succeeded} task(s) recorded "
            f"SUCCEEDED and {unfinished} need re-execution"
        )
    else:
        verdict.kind = RunKind.RE_RUN
        verdict.reason = (
            f"run key {run_key!r} already has {verdict.prior_rows} task(s); "
            f"executing everything again without reusing state"
        )
    return verdict


@dataclass
class SkipDecision:
    """Whether one task can be skipped, and on what evidence."""

    stage: str
    subject: str
    skip: bool
    reason: str
    state: str = ""

    def to_row(self) -> Dict[str, Any]:
        return {
            "stage": self.stage, "subject": self.subject, "skip": self.skip,
            "reason": self.reason, "state": self.state,
        }


def resumable(
    store: Optional[ExecutionStore],
    run_key: str,
    stage: str,
    subject: str,
    spec: Any,
    *,
    config_hash: Optional[str] = None,
) -> SkipDecision:
    """Should this task be skipped during a resume?

    The rule, in order:

    1. no prior row -> execute;
    2. prior row is not ``SUCCEEDED`` -> execute (it never finished);
    3. the configuration hash changed -> execute (the recorded success
       describes a different run);
    4. the outputs do not re-validate **now** -> execute (files were deleted,
       truncated or corrupted since);
    5. otherwise -> skip.

    Step 4 is the point of the whole exercise. A recorded ``SUCCEEDED`` is a
    claim about the past; the filesystem is the only authority on the
    present.
    """
    from .resume import decide

    if store is None:
        return SkipDecision(stage, subject, False, "no execution store")
    row = store.get(run_key, stage, subject)
    if row is None:
        return SkipDecision(stage, subject, False, "no prior state for this task")
    state = str(row.get("state") or "")
    if state != StageState.SUCCEEDED.value:
        return SkipDecision(
            stage, subject, False, f"prior state was {state}, not SUCCEEDED", state)
    if config_hash is not None:
        recorded = row.get("config_hash")
        if recorded and recorded != config_hash:
            return SkipDecision(
                stage, subject, False,
                f"configuration changed since the recorded success "
                f"({recorded} -> {config_hash})", state)

    decision = decide(store, run_key, stage, subject, spec=spec)
    if decision.action.value == "SKIP":
        return SkipDecision(
            stage, subject, True,
            f"prior SUCCEEDED and outputs re-validate now: {decision.reason}",
            state)
    return SkipDecision(
        stage, subject, False,
        f"prior SUCCEEDED but outputs re-validate as "
        f"{getattr(decision.observed, 'value', decision.observed)}: "
        f"{decision.reason}", state)
