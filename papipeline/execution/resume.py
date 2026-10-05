"""Resume decisions.

A resume must answer one question per task: *can I trust that this already
ran correctly?* The answer is never "the row says SUCCEEDED" and never
"the files are there". It is re-derived:

* the recorded state is consulted for **why** it last ran, and
* the declared outputs are **re-validated now**.

That is what makes a process killed at 90% safe. Such a run leaves files
behind. Without re-validation the resume sees files, believes success, and
carries a truncated annotation into every downstream stage.

The five outcomes the brief asks for map onto :class:`ResumeAction`:

===========================  ==========================================
recorded / observed          action
===========================  ==========================================
valid success                ``SKIP``    - re-validated now, keep
invalid output               ``RERUN``   - INVALID
incomplete output            ``RERUN``   - INCOMPLETE
failed execution             ``RERUN``   - FAILED
not started                  ``RUN``     - PENDING
prerequisite not satisfied   ``BLOCKED`` - cannot run yet
===========================  ==========================================
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from ..errors import PipelineError
from .state import StageState
from .store import ExecutionRow, ExecutionStore
from .validation import OutputSpec, ValidationResult, validate


class ResumeAction(str, enum.Enum):
    SKIP = "SKIP"
    RUN = "RUN"
    RERUN = "RERUN"
    BLOCKED = "BLOCKED"

    @property
    def executes(self) -> bool:
        return self in (ResumeAction.RUN, ResumeAction.RERUN)


@dataclass(frozen=True)
class ResumeDecision:
    """Why a task will or will not run."""

    stage: str
    subject: str
    action: ResumeAction
    observed: StageState
    reason: str
    validation: Optional[ValidationResult] = None
    detail: str = ""

    def to_row(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "subject": self.subject,
            "action": self.action.value,
            "observed": self.observed.value,
            "reason": self.reason,
            "validation_ok": None if self.validation is None else self.validation.ok,
            "detail": self.detail,
        }


@dataclass
class ResumePlan:
    """The decisions for one run, plus the counts an operator needs."""

    decisions: List[ResumeDecision] = field(default_factory=list)

    def add(self, decision: ResumeDecision) -> None:
        self.decisions.append(decision)

    def for_task(self, stage: str, subject: str = "") -> Optional[ResumeDecision]:
        for d in self.decisions:
            if d.stage == stage and d.subject == subject:
                return d
        return None

    @property
    def to_execute(self) -> List[ResumeDecision]:
        return [d for d in self.decisions if d.action.executes]

    @property
    def to_skip(self) -> List[ResumeDecision]:
        return [d for d in self.decisions if d.action is ResumeAction.SKIP]

    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {a.value: 0 for a in ResumeAction}
        for d in self.decisions:
            out[d.action.value] += 1
        return out

    def summary(self) -> str:
        c = self.counts()
        return (f"{c['SKIP']} skip, {c['RUN']} to run, {c['RERUN']} to re-run, "
                f"{c['BLOCKED']} blocked")


def decide(
    store: Optional[ExecutionStore],
    run_key: str,
    stage: str,
    subject: str = "",
    *,
    spec: OutputSpec,
    config_hash: Optional[str] = None,
    revalidate: bool = True,
) -> ResumeDecision:
    """Decide what to do with one task on resume.

    Args:
        store: The execution store, or ``None`` for a cold start.
        run_key: Run/cohort key.
        stage: Stage name.
        subject: Sample id, or ``""`` for a cohort-level task.
        spec: The declared output contract. Without one the recorded state is
            trusted, because there is nothing to re-validate against; the
            decision says so in its reason rather than implying a check ran.
        config_hash: When it differs from the recorded one, a previously
            successful task must re-run: its outputs were produced under a
            different configuration.
        revalidate: Set ``False`` only for a cheap pre-flight listing.
    """
    row: Optional[ExecutionRow] = (
        store.get(run_key, stage, subject) if store is not None else None
    )

    if row is None:
        return ResumeDecision(stage, subject, ResumeAction.RUN, StageState.PENDING,
                              "no recorded state")

    recorded = StageState.coerce(row.get("state"))
    recorded_config = row.get("config_hash")

    if recorded is StageState.SUCCEEDED and config_hash and recorded_config \
            and recorded_config != config_hash:
        return ResumeDecision(
            stage, subject, ResumeAction.RERUN, recorded,
            "configuration hash changed since the recorded success",
            detail=f"{recorded_config} -> {config_hash}",
        )

    if recorded is not StageState.SUCCEEDED:
        return ResumeDecision(
            stage, subject, ResumeAction.RERUN, recorded,
            f"last recorded state was {recorded.value}",
            detail=row.get("error") or row.get("validation_detail") or "",
        )

    # Recorded SUCCEEDED. Do not believe it without looking.
    if spec is None:
        raise PipelineError(
            "decide() requires an OutputSpec. Skipping on a recorded "
            "SUCCEEDED without re-validating the files is exactly the failure "
            "this module exists to prevent.",
            run_key=run_key, stage=stage, subject=subject,
        )
    if not revalidate:
        return ResumeDecision(
            stage, subject, ResumeAction.SKIP, recorded,
            "recorded SUCCEEDED; re-validation was explicitly disabled by the "
            "caller, so this skip is unverified",
        )

    result = validate(spec)
    if result.ok:
        return ResumeDecision(
            stage, subject, ResumeAction.SKIP, recorded,
            "recorded SUCCEEDED and outputs re-validated now", validation=result,
        )
    return ResumeDecision(
        stage, subject, ResumeAction.RERUN, recorded,
        f"recorded SUCCEEDED but outputs re-validate as {result.state.value}",
        validation=result, detail=result.detail,
    )


def build_plan(
    store: Optional[ExecutionStore],
    run_key: str,
    tasks: Sequence[tuple],
    specs: Optional[Dict[tuple, OutputSpec]] = None,
    *,
    config_hash: Optional[str] = None,
    predecessors: Optional[Dict[str, Sequence[str]]] = None,
) -> ResumePlan:
    """Decide for a whole run, honouring stage dependencies.

    Args:
        tasks: ``(stage, subject)`` pairs, in execution order.
        specs: Optional ``{(stage, subject): OutputSpec}`` for re-validation.
        config_hash: Compared against each row's recorded hash.
        predecessors: ``{stage: [dependency stages]}``, taken from
            ``papipeline.run.PREREQUISITES``. A task whose prerequisite did
            not succeed is ``BLOCKED`` rather than run, so a resume never
            executes a stage on top of a failed one.
    """
    specs = specs or {}
    plan = ResumePlan()
    satisfied: Dict[str, bool] = {}
    for stage, subject in tasks:
        deps = list((predecessors or {}).get(stage, ()))
        unmet = [d for d in deps if not satisfied.get(d, False)]
        if unmet:
            plan.add(ResumeDecision(
                stage, subject, ResumeAction.BLOCKED, StageState.PENDING,
                f"prerequisite(s) not succeeded: {', '.join(unmet)}",
            ))
            satisfied[stage] = False
            continue
        spec = specs.get((stage, subject))
        if spec is None:
            # A task with no contract cannot be re-validated, so it cannot be
            # resumed safely. Run it rather than guess.
            plan.add(ResumeDecision(
                stage, subject, ResumeAction.RUN, StageState.PENDING,
                "no output contract declared; cannot re-validate a recorded "
                "success, so this stage will be run",
            ))
            satisfied[stage] = False
            continue
        decision = decide(store, run_key, stage, subject, spec=spec,
                          config_hash=config_hash)
        plan.add(decision)
        satisfied[stage] = decision.action is ResumeAction.SKIP
    return plan
