"""Execution state, output validation, retry and resume.

Infrastructure only. Nothing in this package interprets biology, sets a
scientific threshold, or alters an analysis. It exists so that a stage's
execution state is *truthful*:

* a task is ``SUCCEEDED`` only when its outputs passed a declared contract;
* ``FAILED``, ``INCOMPLETE`` and ``INVALID`` are distinct, because they
  need different responses and only one of them is worth retrying;
* a resume re-validates rather than trusting a row or a file listing.

See ``docs/execution_state_and_validation.md``.
"""

from .retry import AttemptRecord, RetryPolicy
from .resume import ResumeAction, ResumeDecision, ResumePlan, build_plan, decide
from .runner import TaskContext, TaskResult, run_task
from .state import (
    ExecutionOutcome,
    FailureKind,
    StageState,
    TaskEvent,
    classify_process_failure,
)
from .store import ExecutionStore, host_info
from .validation import (
    Check,
    CheckKind,
    CheckResult,
    OutputSpec,
    SiblingSpec,
    ValidationResult,
    covers_samples,
    exists,
    has_columns,
    has_real_values,
    identifies_file,
    identifiers_known,
    min_rows,
    non_empty,
    not_constant_matrix,
    parses,
    validate,
)

__all__ = [
    "AttemptRecord",
    "Check",
    "CheckKind",
    "CheckResult",
    "ExecutionOutcome",
    "ExecutionStore",
    "FailureKind",
    "OutputSpec",
    "ResumeAction",
    "ResumeDecision",
    "ResumePlan",
    "RetryPolicy",
    "SiblingSpec",
    "StageState",
    "TaskEvent",
    "TaskContext",
    "TaskResult",
    "ValidationResult",
    "build_plan",
    "classify_process_failure",
    "covers_samples",
    "decide",
    "exists",
    "has_columns",
    "has_real_values",
    "host_info",
    "identifies_file",
    "identifiers_known",
    "min_rows",
    "non_empty",
    "not_constant_matrix",
    "parses",
    "run_task",
    "validate",
]
