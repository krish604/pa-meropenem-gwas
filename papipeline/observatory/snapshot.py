"""A read-only projection of the real execution store.

Every value here comes from a row that :mod:`papipeline.execution` actually
wrote. Where the engine does not record something, this module returns
``None`` and the UI shows that it is unavailable — it does not substitute a
plausible number. That rule is the whole point: a dashboard that invents
its own telemetry is worse than no dashboard, because it looks like
evidence.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from ..execution.state import StageState
from ..execution.store import ExecutionRow, ExecutionStore
from .graph import build_graph, label_for

#: Aggregating many subjects into one stage node needs a rule, and the rule
#: has to favour the state an operator must act on. Ordered most to least
#: urgent: something is happening now, then something is being retried,
#: then something is wrong, then something is merely unfinished. A stage
#: with 64 successes and 1 failure is a failure, not a success.
STATE_PRECEDENCE: Sequence[StageState] = (
    StageState.RUNNING,
    StageState.RETRYING,
    StageState.INVALID,
    StageState.FAILED,
    StageState.INCOMPLETE,
    StageState.PENDING,
    StageState.SUCCEEDED,
)

_JSON_COLUMNS = ("command", "input_ids", "output_paths", "validation_json")


def _decode(row: ExecutionRow) -> Dict[str, Any]:
    """One row, with the JSON columns parsed and unknown states tolerated."""
    out: Dict[str, Any] = dict(row)
    for column in _JSON_COLUMNS:
        raw = out.get(column)
        if isinstance(raw, str) and raw:
            try:
                out[column] = json.loads(raw)
            except (ValueError, TypeError):
                out[column] = None
        elif raw == "":
            out[column] = None
    state = out.get("state")
    try:
        out["state"] = StageState.coerce(state).value
    except Exception:  # a state this build does not know about
        out["state"] = str(state) if state else StageState.PENDING.value
    return out


def aggregate_state(counts: Dict[str, int], total: int) -> str:
    """Roll subject states up into one stage state.

    ``total == 0`` means the stage has never been attempted, which is
    ``PENDING`` and is genuinely different from "attempted and unfinished".
    """
    if total <= 0:
        return StageState.PENDING.value
    for state in STATE_PRECEDENCE:
        if counts.get(state.value, 0) > 0:
            return state.value
    return StageState.PENDING.value


def task_detail(row: ExecutionRow) -> Dict[str, Any]:
    """Everything the left panel shows about one task."""
    decoded = _decode(row)
    attempts = decoded.get("attempt") or 0
    max_attempts = decoded.get("max_attempts") or 1
    elapsed = decoded.get("elapsed_seconds")
    return {
        "run_key": decoded.get("run_key"),
        "stage": decoded.get("stage"),
        "stage_label": label_for(str(decoded.get("stage") or "")),
        "subject": decoded.get("subject") or "",
        "state": decoded.get("state"),
        "attempt": attempts,
        "max_attempts": max_attempts,
        "failure_kind": decoded.get("failure_kind"),
        "started_at": decoded.get("started_at"),
        "ended_at": decoded.get("ended_at"),
        "elapsed_seconds": elapsed if elapsed is not None else None,
        "command": decoded.get("command"),
        "tool_version": decoded.get("tool_version"),
        "database_version": decoded.get("database_version"),
        "config_hash": decoded.get("config_hash"),
        "input_ids": decoded.get("input_ids"),
        "output_paths": decoded.get("output_paths"),
        "validation_state": decoded.get("validation_state"),
        "validation_detail": decoded.get("validation_detail"),
        "validation": decoded.get("validation_json"),
        "error": decoded.get("error"),
        # Deliberately not claimed here: the execution table has no log
        # column, so the store's log_path() lookup is the only honest source.
        # The API enriches task detail with it.
        "log_path": None,
        "host": decoded.get("host"),
        "cpu_count": decoded.get("cpu_count"),
        "memory_mb": decoded.get("memory_mb"),
        "updated_at": decoded.get("updated_at"),
    }


@dataclass
class Snapshot:
    """A consistent read of one run at one moment."""

    run_key: str
    exists: bool
    nodes: List[Dict[str, Any]] = field(default_factory=list)
    edges: List[Dict[str, Any]] = field(default_factory=list)
    tasks: List[Dict[str, Any]] = field(default_factory=list)
    counts: Dict[str, int] = field(default_factory=dict)
    subjects: List[str] = field(default_factory=list)
    total_expected: int = 0
    elapsed_seconds: float = 0.0

    @property
    def idle(self) -> bool:
        """True when there is genuinely nothing to show.

        No store, no rows, or nothing running and nothing in flight. The
        UI shows PIPELINE IDLE rather than an empty graph pretending to be
        a live run.
        """
        if not self.exists or not self.tasks:
            return True
        live = (StageState.RUNNING.value, StageState.RETRYING.value)
        return not any(t["state"] in live for t in self.tasks)

    def to_row(self) -> Dict[str, Any]:
        return {
            "run_key": self.run_key,
            "exists": self.exists,
            "idle": self.idle,
            "nodes": self.nodes,
            "edges": self.edges,
            "tasks": self.tasks,
            "counts": self.counts,
            "subjects": self.subjects,
            "total_expected": self.total_expected,
            "elapsed_seconds": self.elapsed_seconds,
        }


def read_snapshot(
    store: Optional[ExecutionStore],
    run_key: str,
    order: Optional[Sequence[str]] = None,
) -> Snapshot:
    """Project the store into a snapshot.

    Stages with no rows are still present as nodes, in ``PENDING`` with
    ``total == 0``. A stage that has never run is part of the pipeline's
    structure and hiding it would misdescribe the system.
    """
    from .graph import default_order
    order = order if order is not None else default_order()
    nodes, edges = build_graph(order)
    snap = Snapshot(run_key=run_key, exists=store is not None, edges=[e.to_row() for e in edges])

    rows: List[ExecutionRow] = []
    if store is not None:
        try:
            rows = store.list_for_run(run_key)
        except Exception:
            # A missing or unreadable database is an idle observatory, not a
            # crash. The UI is the wrong place to surface a storage error.
            rows = []
            snap.exists = False

    by_stage: Dict[str, List[Dict[str, Any]]] = {}
    subjects = set()
    for raw in rows:
        detail = task_detail(raw)
        by_stage.setdefault(str(detail["stage"]), []).append(detail)
        if detail["subject"]:
            subjects.add(str(detail["subject"]))
        elapsed = detail.get("elapsed_seconds")
        if elapsed:
            snap.elapsed_seconds += float(elapsed)

    snap.tasks = [t for group in by_stage.values() for t in group]
    snap.tasks.sort(key=lambda t: (t["stage"], t["subject"]))
    snap.subjects = sorted(subjects)

    overall: Dict[str, int] = {}
    for node in nodes:
        group = by_stage.get(node.stage, [])
        counts: Dict[str, int] = {}
        for task in group:
            counts[task["state"]] = counts.get(task["state"], 0) + 1
            overall[task["state"]] = overall.get(task["state"], 0) + 1
        row = node.to_row()
        row["counts"] = counts
        row["total"] = len(group)
        row["state"] = aggregate_state(counts, len(group))
        row["subjects"] = sorted(str(t["subject"]) for t in group)
        snap.nodes.append(row)

    snap.counts = overall
    # "Expected" means the cohort the run is working on: every subject that
    # appears anywhere in the store. Stages that legitimately have one
    # cohort-level row carry subject "" and are not counted as genomes.
    snap.total_expected = len([s for s in subjects if s])
    return snap
