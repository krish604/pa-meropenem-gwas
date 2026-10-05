"""The observatory's HTTP surface: read-only views of real state, plus SSE.

Every route reads the execution store or the operating system. There is no
route that invents a value, and no route that mutates pipeline state — the
observatory observes. The one write-shaped thing it can do is start a
sandboxed *demo* execution, which is clearly labelled and writes only to
its own database and directory.

The SSE stream pushes two kinds of message and nothing else:

``snapshot``
    the full current picture, so a reconnecting client is immediately
    correct without replaying history;
``event``
    a real transition that just happened.

A client is told its buffer position so it can detect a gap and re-read
rather than silently presenting a partial history.
"""

from __future__ import annotations

import asyncio
import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

from ..execution.store import SCHEMA_VERSION, ExecutionStore
from . import events as ev
from . import graph as g
from . import metrics as mx
from . import snapshot as sn
from . import throughput as tp

STATIC_ROOT = Path(__file__).parent / "static"
DEFAULT_RUN_KEY = os.environ.get("PAPIPELINE_OBSERVATORY_RUN", "observatory")
DEFAULT_DB = os.environ.get(
    "PAPIPELINE_OBSERVATORY_DB", str(Path.home() / ".local/share/papipeline/observatory.db")
)
#: How often the stream re-reads the store. Events are pushed as they
#: happen; this interval is the backstop for anything an event cannot
#: describe, and the heartbeat.
POLL_SECONDS = float(os.environ.get("PAPIPELINE_OBSERVATORY_POLL", "1.0"))


def _open_store(path: Optional[str]) -> Optional[ExecutionStore]:
    """Open the store read-only-ish, or ``None`` when there is not one yet.

    A missing database is a normal state for a fresh machine, not an error.
    """
    target = Path(path or DEFAULT_DB)
    if not target.exists():
        return None
    try:
        return ExecutionStore(target)
    except Exception:  # noqa: BLE001
        return None


def build_snapshot(store: Optional[ExecutionStore], run_key: str) -> sn.Snapshot:
    return sn.read_snapshot(store, run_key)


def payload(
    store: Optional[ExecutionStore],
    run_key: str,
    *,
    include_tasks: bool = True,
) -> Dict[str, Any]:
    """The whole picture in one object, as the client receives it."""
    snap = build_snapshot(store, run_key)
    host = mx.sample_host()
    rates = tp.compute(snap.tasks, snap.total_expected)
    body: Dict[str, Any] = snap.to_row()
    body["structure"] = g.graph_payload()
    body["host"] = host.to_row()
    body["throughput"] = rates.to_row()
    body["median_task_seconds"] = tp.median_duration(snap.tasks)
    body["store"] = {
        "path": str(store.path) if store is not None else None,
        "schema_version": SCHEMA_VERSION,
        "present": store is not None,
    }
    # A store holding rows for other run keys is a mismatch, and the sixteen
    # PENDING nodes below are otherwise indistinguishable from a run that has
    # not started. `mismatch` lets a client say so instead of drawing an empty
    # pipeline that looks healthy.
    available = store.run_keys() if store is not None else []
    body["run_key_mismatch"] = {
        "requested": run_key,
        "available": available,
        "mismatch": bool(available) and run_key not in available,
    }
    body["capabilities"] = {
        "psutil": mx.psutil_available(),
        "per_task_cpu": False,
        "per_task_progress": False,
        "note": (
            "This build records host metrics and per-task elapsed time only. "
            "The execution engine stores no per-task CPU, memory or progress, "
            "so those fields are reported as unavailable rather than "
            "estimated."
        ),
    }
    if not include_tasks:
        body.pop("tasks", None)
    return body


def create_app(
    db_path: Optional[str] = None,
    run_key: str = DEFAULT_RUN_KEY,
    *,
    bus: Optional[ev.EventBus] = None,
) -> FastAPI:
    """Build the app. Injected store path and bus keep it testable."""
    app = FastAPI(
        title="PAPipeline Computational Observatory",
        version="0.1.0",
        description=(
            "Read-only observability over PAPipeline execution state. "
            "Every value is read from the execution store or the host; "
            "nothing here is simulated."
        ),
    )
    app.state.db_path = db_path
    app.state.run_key = run_key
    app.state.bus = bus if bus is not None else ev.BUS
    app.state.started_at = time.time()

    def store() -> Optional[ExecutionStore]:
        opened = _open_store(app.state.db_path)
        if opened is None:
            return None
        return opened

    # -- meta ------------------------------------------------------------

    @app.get("/api/health")
    def health() -> Dict[str, Any]:
        opened = store()
        # Which keys this database actually holds. Reported so that pointing the
        # server at the wrong run key is visible here rather than showing up
        # only as an all-PENDING snapshot.
        available = opened.run_keys() if opened is not None else []
        return {
            "status": "ok",
            "run_key": app.state.run_key,
            "database": app.state.db_path,
            "database_present": Path(app.state.db_path).exists() if app.state.db_path else False,
            "event_sequence": app.state.bus.sequence,
            "psutil": mx.psutil_available(),
            "available_run_keys": available,
            # False when the database holds runs but not the requested one.
            # True when the key matches, and also when the database is empty or
            # absent - "no runs yet" is not a mismatch.
            "run_key_found": bool(available) and app.state.run_key in available,
        }

    @app.get("/api/structure")
    def structure() -> Dict[str, Any]:
        """The pipeline's own declared graph, independent of any run."""
        return g.graph_payload()

    @app.get("/api/snapshot")
    def snapshot_route(
        run: Optional[str] = Query(default=None),
        tasks: bool = Query(default=True),
    ) -> Dict[str, Any]:
        opened = store()
        return payload(
            opened, run or app.state.run_key, include_tasks=tasks
        )

    @app.get("/api/host")
    def host_route() -> Dict[str, Any]:
        return mx.sample_host().to_row()

    @app.get("/api/tasks")
    def tasks_route(
        run: Optional[str] = Query(default=None),
        stage: Optional[str] = Query(default=None),
        state: Optional[str] = Query(default=None),
    ) -> List[Dict[str, Any]]:
        """One flat list, for filtering in the client."""
        opened = store()
        key = run or app.state.run_key
        rows: List[Dict[str, Any]] = []
        if opened is not None:
            for row in opened.list_for_run(key):
                detail = sn.task_detail(row)
                if stage and detail["stage"] != stage:
                    continue
                if state and detail["state"] != state:
                    continue
                rows.append(detail)
        rows.sort(key=lambda r: (r["stage"], r["subject"]))
        return rows

    @app.get("/api/tasks/{stage}")
    @app.get("/api/tasks/{stage}/{subject}")
    def task_route(stage: str, subject: str = "") -> Dict[str, Any]:
        """One task, with its real attempt history.

        Both spellings are served: a cohort-level stage has no subject, and
        the client should not have to know that to ask for it.
        """
        """One task, with its real attempt history."""
        opened = store()
        if opened is None:
            raise HTTPException(404, "no execution store")
        row = opened.get(app.state.run_key, stage, subject)
        if row is None:
            raise HTTPException(404, f"no recorded task {stage}/{subject}")
        detail = sn.task_detail(row)
        history: List[Dict[str, Any]] = []
        if opened is not None:
            # The execution table has no log column: the runner records the
            # path per attempt. Enrich the detail so the left panel can offer
            # the log that actually exists.
            detail["log_path"] = opened.log_path(
                app.state.run_key, stage, subject)
            for record in opened.attempts(app.state.run_key, stage, subject):
                history.append({
                    "attempt": record.attempt,
                    "state": record.state.value,
                    "failure_kind": (record.failure_kind.value
                                     if record.failure_kind else None),
                    "detail": record.detail,
                    "exit_code": record.exit_code,
                    "elapsed_seconds": record.elapsed_seconds,
                })
        detail["attempts_history"] = history
        return detail

    @app.get("/api/logs/{stage}")
    @app.get("/api/logs/{stage}/{subject}")
    def logs_route(
        stage: str,
        subject: str = "",
        tail: int = Query(default=2000, ge=1, le=200000),
    ) -> Dict[str, Any]:
        """A task's real log, read from the file the runner wrote.

        The path comes from the store, not from the request, so this route
        cannot be used to read arbitrary files.
        """
        opened = store()
        if opened is None:
            raise HTTPException(404, "no execution store")
        row = opened.get(app.state.run_key, stage, subject)
        if row is None:
            raise HTTPException(404, "no recorded task")
        # The runner records the log path per attempt, not per task row.
        path = opened.log_path(app.state.run_key, stage, subject) or row.get("log_path")
        if not path:
            return {"available": False,
                    "reason": "this task recorded no log path", "text": ""}
        target = Path(str(path))
        if not target.exists():
            return {"available": False,
                    "reason": f"log not present at {target}", "text": ""}
        try:
            text = target.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            return {"available": False, "reason": str(exc), "text": ""}
        return {
            "available": True,
            "path": str(target),
            "bytes": target.stat().st_size,
            "text": text[-tail:],
            "truncated": len(text) > tail,
        }

    @app.get("/api/events")
    def events_route(
        since: int = Query(default=0, ge=0),
        run: Optional[str] = Query(default=None),
    ) -> Dict[str, Any]:
        """Buffered events, for a client that would rather poll."""
        bus: ev.EventBus = app.state.bus
        key = run or app.state.run_key
        return {
            "sequence": bus.sequence,
            "gap": bus.gap_since(since),
            "events": [e.to_row() for e in bus.recent(since)
                       if e.run_key == key],
        }

    # -- the live stream -------------------------------------------------

    @app.get("/api/stream")
    async def stream(request: Request) -> StreamingResponse:
        """Server-sent events: snapshots and real transitions.

        SSE rather than a websocket because the traffic is strictly
        server-to-client. That drops the dependency on a socket library
        and makes a reconnect a plain HTTP retry, which is what a
        long-running pipeline on a laptop actually needs.
        """
        bus: ev.EventBus = app.state.bus
        loop = asyncio.get_event_loop()
        queue: "asyncio.Queue[Optional[Dict[str, Any]]]" = asyncio.Queue(maxsize=256)

        def on_event(event: ev.ExecutionEvent) -> None:
            # A bus may be shared by several observed runs. A client is only
            # ever shown its own run's transitions.
            if event.run_key != app.state.run_key:
                return

            def put() -> None:
                try:
                    queue.put_nowait({"type": "event", "data": event.to_row()})
                except asyncio.QueueFull:
                    # The client is too slow to keep up. Dropping the event
                    # is safe: the next snapshot re-states the truth.
                    pass

            try:
                loop.call_soon_threadsafe(put)
            except RuntimeError:
                pass

        unsubscribe = bus.subscribe(on_event)

        async def body():
            try:
                yield _sse({
                    "type": "hello",
                    "run_key": app.state.run_key,
                    "sequence": bus.sequence,
                    "capabilities": payload(store(), app.state.run_key,
                                            include_tasks=False)["capabilities"],
                })
                last = 0.0
                while True:
                    if await request.is_disconnected():
                        break
                    try:
                        message = await asyncio.wait_for(queue.get(), timeout=POLL_SECONDS)
                    except asyncio.TimeoutError:
                        message = None
                    if message is None:
                        # Backstop snapshot + heartbeat. This is what makes
                        # a dropped or coalesced event harmless.
                        now = time.monotonic()
                        if now - last >= POLL_SECONDS:
                            last = now
                            yield _sse({
                                "type": "snapshot",
                                "data": payload(store(), app.state.run_key),
                            })
                        else:
                            yield ": keep-alive\n\n"
                        continue
                    yield _sse(message)
            finally:
                unsubscribe()

        return StreamingResponse(
            body(),
            media_type="text/event-stream",
            headers={
                "Cache-Control": "no-cache, no-transform",
                "Connection": "keep-alive",
                "X-Accel-Buffering": "no",
            },
        )

    # -- the app itself --------------------------------------------------

    if STATIC_ROOT.is_dir():
        app.mount("/static", StaticFiles(directory=str(STATIC_ROOT)), name="static")

        @app.get("/", response_class=HTMLResponse)
        def index() -> HTMLResponse:
            return HTMLResponse((STATIC_ROOT / "index.html").read_text(encoding="utf-8"))

        @app.get("/healthz", response_class=PlainTextResponse)
        def healthz() -> str:
            return "ok"

    return app


def _sse(message: Dict[str, Any]) -> str:
    """One SSE frame.

    The whole payload goes in ``data:`` as a single JSON line rather than
    splitting named fields, so a client parses one object and cannot
    desynchronise from a field it did not know about.
    """
    return f"event: {message['type']}\ndata: {json.dumps(message.get('data', {}), default=str)}\n\n"


app = create_app()
