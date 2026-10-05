"""Throughput and ETA, derived only from what the store actually recorded.

The engine does not report progress, so there is no honest way to answer
"what percent of this stage is done". What it *does* record is
``elapsed_seconds`` per finished task and the wall-clock span of the run.
From those, and nothing else, we can say how fast the cohort is finishing
work and roughly when the remainder will be done.

Both functions return ``None`` rather than a number when the evidence is
insufficient. A dashboard that shows a confident ETA from two data points
is lying with a decimal point.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

#: Below this many finished tasks, a rate is noise rather than a rate.
MIN_SAMPLES_FOR_RATE = 3


@dataclass
class Throughput:
    """A measured rate and what it implies, or an explicit unknown."""

    completed: int
    failed: int
    remaining: int
    #: Finished tasks per minute, or ``None`` when it cannot be measured.
    tasks_per_minute: Optional[float] = None
    #: Minutes until the remaining work is done, or ``None``.
    eta_minutes: Optional[float] = None
    #: Why a rate or ETA is missing, for display rather than invention.
    basis: str = ""
    span_seconds: float = 0.0

    def to_row(self) -> Dict[str, Any]:
        return {
            "completed": self.completed,
            "failed": self.failed,
            "remaining": self.remaining,
            "tasks_per_minute": self.tasks_per_minute,
            "eta_minutes": self.eta_minutes,
            "basis": self.basis,
            "span_seconds": self.span_seconds,
        }


def _parse_stamp(value: Optional[str]) -> Optional[float]:
    """Epoch seconds from the store's ``%Y-%m-%dT%H:%M:%SZ`` stamps."""
    if not value:
        return None
    text = str(value).strip()
    for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%dT%H:%M:%S"):
        try:
            return time.mktime(time.strptime(text, fmt)) - time.timezone
        except ValueError:
            continue
    return None


def run_span_seconds(tasks: Sequence[Dict[str, Any]]) -> float:
    """Observed wall-clock span of the run, from real start/end stamps.

    Falls back to the summed durations when stamps are missing, and says so
    in the basis rather than pretending it is the same measurement.
    """
    starts = [t for t in (_parse_stamp(x.get("started_at")) for x in tasks) if t]
    ends = [t for t in (_parse_stamp(x.get("ended_at")) for x in tasks) if t]
    if starts and ends:
        span = max(ends) - min(starts)
        if span > 0:
            return float(span)
    return float(sum(float(x.get("elapsed_seconds") or 0.0) for x in tasks))


def compute(
    tasks: Sequence[Dict[str, Any]],
    total_expected: int,
    *,
    now: Optional[float] = None,
    live_states: Sequence[str] = ("RUNNING", "RETRYING"),
) -> Throughput:
    """Rate and ETA for the cohort.

    ``total_expected`` is the number of subjects the run is working on;
    ``remaining`` is what is not yet finished. Running tasks are excluded
    from the completed count, because a task that has not finished has not
    demonstrated a duration.
    """
    now = time.time() if now is None else now
    finished = [
        t for t in tasks
        if t.get("state") == "SUCCEEDED" and (t.get("elapsed_seconds") or 0) > 0
    ]
    failed = [t for t in tasks if t.get("state") in ("FAILED", "INVALID", "INCOMPLETE")]
    live = [t for t in tasks if t.get("state") in live_states]

    remaining = max(0, int(total_expected) - len(finished) - len(failed))
    out = Throughput(
        completed=len(finished),
        failed=len(failed),
        remaining=remaining,
        span_seconds=run_span_seconds(tasks),
    )

    if not finished:
        out.basis = (
            "no task has completed yet, so there is no measured rate"
        )
        return out
    if len(finished) < MIN_SAMPLES_FOR_RATE:
        out.basis = (
            f"only {len(finished)} task(s) finished; a rate from "
            f"{len(finished)} sample(s) is not meaningful"
        )
        return out

    # Prefer real wall-clock span; otherwise use the summed durations, and
    # record which of the two the number came from.
    span = out.span_seconds
    if span <= 0:
        span = float(sum(float(t["elapsed_seconds"]) for t in finished))
        out.basis = "derived from summed task durations, not wall-clock span"
    else:
        out.basis = (
            f"wall-clock span of the run over {len(finished)} completed tasks"
        )

    minutes = span / 60.0
    if minutes <= 0:
        out.basis = "tasks completed too quickly to derive a rate"
        return out

    rate = len(finished) / minutes
    out.tasks_per_minute = round(rate, 4)
    if remaining > 0:
        out.eta_minutes = round(remaining / rate, 1)
    else:
        out.eta_minutes = 0.0
        out.basis += "; nothing remaining, so the ETA is zero"
    if live:
        out.basis += f"; {len(live)} task(s) in flight"
    return out


def median_duration(tasks: Sequence[Dict[str, Any]]) -> Optional[float]:
    """Median seconds a finished task took, or ``None`` if there are none."""
    values = sorted(
        float(t["elapsed_seconds"]) for t in tasks
        if t.get("elapsed_seconds") and t.get("state") == "SUCCEEDED"
    )
    if not values:
        return None
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return (values[mid - 1] + values[mid]) / 2.0
