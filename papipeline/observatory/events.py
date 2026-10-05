"""Real execution events, emitted by the code that actually executes.

There is no second state system here. An event is a *notification* that
:mod:`papipeline.execution.runner` derived from a state transition it had
already made and written to the store; the event carries no state of its
own and the UI never treats it as the source of truth. The store remains
authoritative and the UI re-reads it, so a dropped event costs a frame of
freshness and nothing else.

Events reach this module through one additive, optional parameter on
``run_task`` (``event_sink``). No existing behaviour changes when it is
``None``, which is why it is safe to add before the orchestrator is wired.
"""

from __future__ import annotations

import json
import threading
import time
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, Optional

from ..execution.state import StageState

#: Event names. These are the transitions the execution layer can actually
#: make; there is no TASK_PROGRESS because the engine reports no progress.
TASK_CREATED = "TASK_CREATED"
TASK_STARTED = "TASK_STARTED"
TASK_COMPLETED = "TASK_COMPLETED"
TASK_VALIDATED = "TASK_VALIDATED"
TASK_FAILED = "TASK_FAILED"
TASK_INVALID = "TASK_INVALID"
TASK_INCOMPLETE = "TASK_INCOMPLETE"
TASK_RETRY = "TASK_RETRY"
TASK_RESUMED = "TASK_RESUMED"

#: Run-level events. These bracket a whole execution. They are what makes a
#: second run distinguishable from a resume in the event stream: a re-run
#: emits RUN_STARTED, and only a resume emits RUN_RESUMED.
RUN_STARTED = "RUN_STARTED"
RUN_RESUMED = "RUN_RESUMED"
RUN_COMPLETED = "RUN_COMPLETED"
RUN_FAILED = "RUN_FAILED"

#: Emitted for a task that a resume deliberately did not re-execute.
TASK_SKIPPED = "TASK_SKIPPED"

ALL_EVENTS = (
    TASK_CREATED, TASK_STARTED, TASK_COMPLETED, TASK_VALIDATED,
    TASK_FAILED, TASK_INVALID, TASK_INCOMPLETE, TASK_RETRY, TASK_RESUMED,
    TASK_SKIPPED, RUN_STARTED, RUN_RESUMED, RUN_COMPLETED, RUN_FAILED,
)


@dataclass(frozen=True)
class ExecutionEvent:
    """One real transition, as it happened."""

    name: str
    stage: str
    subject: str
    run_key: str
    state: str
    attempt: int = 0
    at: float = field(default_factory=time.time)
    detail: str = ""
    payload: Dict[str, Any] = field(default_factory=dict)
    #: Position in the emitting bus. Assigned by :meth:`EventBus.emit`, so a
    #: client can ask for "everything after N" without relying on wall-clock
    #: time, which is not unique within a burst.
    seq: int = 0

    def to_row(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "stage": self.stage,
            "subject": self.subject,
            "run_key": self.run_key,
            "state": self.state,
            "attempt": self.attempt,
            "seq": self.seq,
            "at": self.at,
            "detail": self.detail,
            "payload": dict(self.payload),
        }

    def to_json(self) -> str:
        return json.dumps(self.to_row(), sort_keys=True)


def event_for_state(
    state: StageState,
    *,
    run_key: str,
    stage: str,
    subject: str,
    attempt: int,
    detail: str = "",
    payload: Optional[Dict[str, Any]] = None,
) -> ExecutionEvent:
    """Name a transition. The mapping is total and explicit.

    Mapping ``SUCCEEDED`` to ``TASK_COMPLETED`` and then, separately,
    ``TASK_VALIDATED`` is deliberate: they are two different facts. A task
    that exited cleanly and a task whose outputs passed their contract are
    not the same event, and the UI shows them as different rows.
    """
    name = {
        StageState.RUNNING: TASK_STARTED,
        StageState.RETRYING: TASK_RETRY,
        StageState.SUCCEEDED: TASK_COMPLETED,
        StageState.FAILED: TASK_FAILED,
        StageState.INVALID: TASK_INVALID,
        StageState.INCOMPLETE: TASK_INCOMPLETE,
        StageState.PENDING: TASK_CREATED,
    }.get(state, TASK_STARTED)
    return ExecutionEvent(
        name=name,
        stage=stage,
        subject=subject,
        run_key=run_key,
        state=state.value,
        attempt=attempt,
        detail=detail,
        payload=dict(payload or {}),
    )


class EventBus:
    """Fan-out to subscribers, with a bounded replay buffer.

    Bounded on purpose: a long run emits a great many events and the UI
    only needs enough to fill its history panel on connect. A subscriber
    that falls behind loses the oldest events rather than growing without
    limit, and the client learns its buffer position so it can re-read the
    store instead of assuming continuity.
    """

    def __init__(self, capacity: int = 500) -> None:
        self._lock = threading.Lock()
        self._subscribers: List[Callable[[ExecutionEvent], None]] = []
        self._buffer: List[ExecutionEvent] = []
        self._sequence = 0
        self._capacity = max(1, capacity)
        #: Diagnostics for the API: how many events this bus has carried.
        self.emitted = 0

    @property
    def sequence(self) -> int:
        with self._lock:
            return self._sequence

    def subscribe(self, callback: Callable[[ExecutionEvent], None]) -> Callable[[], None]:
        with self._lock:
            self._subscribers.append(callback)

        def unsubscribe() -> None:
            with self._lock:
                if callback in self._subscribers:
                    self._subscribers.remove(callback)

        return unsubscribe

    def emit(self, event: ExecutionEvent) -> None:
        with self._lock:
            self._sequence += 1
            event = replace(event, seq=self._sequence)
            self._buffer.append(event)
            if len(self._buffer) > self._capacity:
                del self._buffer[: len(self._buffer) - self._capacity]
            subscribers = list(self._subscribers)
            self.emitted += 1
        for callback in subscribers:
            try:
                callback(event)
            except Exception:  # noqa: BLE001
                # One broken subscriber must not stop the others, and must
                # never propagate back into the execution engine.
                continue

    def recent(self, since: int = 0, limit: int = 200) -> List[ExecutionEvent]:
        """Buffered events after a sequence number.

        A client reconnecting passes the last sequence it saw. If that
        number has fallen out of the buffer the gap is reported explicitly
        so the client re-reads the store rather than showing a partial
        history as though it were complete.
        """
        with self._lock:
            if not since:
                return self._buffer[-limit:]
            oldest = self._sequence - len(self._buffer)
            if since < oldest:
                # The requested point has already fallen out of the buffer.
                return []
            return [e for e in self._buffer if e.seq > since][-limit:]

    def gap_since(self, since: int) -> int:
        """How many events were dropped between ``since`` and the buffer."""
        with self._lock:
            if not self._buffer:
                return 0
            oldest = self._sequence - len(self._buffer)
            return max(0, oldest - since)

    def sink(self) -> Callable[[ExecutionEvent], None]:
        """A ``run_task``-compatible event sink for this bus."""

        def emit(event: ExecutionEvent) -> None:
            self.emit(event)

        return emit


#: Process-wide bus. A module-level singleton is appropriate here: the
#: observatory is a view onto one pipeline, and a second bus would mean a
#: second, competing story about what is running.
BUS = EventBus()


def set_bus(bus: EventBus) -> EventBus:
    """Install a bus, for tests and for embedding the app in a host process."""
    global BUS
    BUS = bus
    return bus


def adapt(event: Any) -> ExecutionEvent:
    """Turn an execution-layer :class:`TaskEvent` into an observatory event.

    A ``SUCCEEDED`` transition produces *two* events when validation
    actually ran: ``TASK_COMPLETED`` for the process and ``TASK_VALIDATED``
    for the contract. That is a list because the mapping is not always
    one-to-one, and collapsing them would hide the distinction the whole
    execution layer exists to make.
    """
    from ..execution.state import TaskEvent

    if not isinstance(event, TaskEvent):
        raise TypeError(f"expected a TaskEvent, got {type(event).__name__}")

    common = dict(
        run_key=event.run_key, stage=event.stage, subject=event.subject,
        attempt=event.attempt,
    )
    out = [ExecutionEvent(
        name=event_for_state(event.state, detail=event.detail, **common).name,
        state=event.state.value, detail=event.detail,
        payload={"exit_code": event.exit_code,
                 "elapsed_seconds": event.elapsed_seconds},
        **common,
    )]
    if event.validated:
        out.append(ExecutionEvent(
            name=TASK_VALIDATED,
            state=event.state.value,
            detail=event.validation_detail or "outputs satisfy their contract",
            payload={"validation_state": event.validation_state},
            **common,
        ))
    return out[0] if len(out) == 1 else _joined(out)


def _joined(events: List[ExecutionEvent]) -> ExecutionEvent:
    """A completed-then-validated pair, as one observable transition."""
    first = events[0]
    return ExecutionEvent(
        name=TASK_COMPLETED,
        stage=first.stage, subject=first.subject, run_key=first.run_key,
        state=first.state, attempt=first.attempt, at=first.at,
        detail=first.detail,
        payload={**first.payload, "validated": True,
                 "validation": events[1].detail},
    )


def bus_sink(bus: Optional[EventBus] = None) -> Callable[[Any], None]:
    """A ``run_task(event_sink=...)`` callable that publishes to ``bus``.

    Pass this straight to ``run_task`` and real transitions reach every SSE
    subscriber::

        run_task(ctx, cmd, store=store, event_sink=bus_sink())
    """

    def sink(event: Any) -> None:
        target = bus if bus is not None else BUS
        target.emit(adapt(event))

    return sink
