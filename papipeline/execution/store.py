"""Durable execution state.

Every stage that runs through :mod:`papipeline.execution` leaves a row here:
its state, the command that produced it, the versions and configuration in
force, the inputs it consumed, the outputs it declared, and the verdict of
the output contract.

Two properties are load-bearing:

1. **Every attempt is kept**, in ``attempt_log``. A task that succeeded on
   its third try must still be able to show what happened on the first two,
   so a late success carries complete provenance.
2. **Access is serialised.** The benchmark and any future concurrent runner
   use one store from several threads. ``check_same_thread=False`` permits
   that without making it safe - interleaved cursor use on a single
   connection has crashed the interpreter - so every statement runs under a
   lock.
"""

from __future__ import annotations

import functools
import json
import os
import platform
import socket
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from ..errors import PipelineError
from .retry import AttemptRecord
from .state import FailureKind, StageState

SCHEMA_VERSION = 1

SCHEMA = """
CREATE TABLE IF NOT EXISTS execution (
    id                INTEGER PRIMARY KEY,
    run_key           TEXT NOT NULL,
    stage             TEXT NOT NULL,
    subject           TEXT NOT NULL DEFAULT '',
    state             TEXT NOT NULL,
    failure_kind      TEXT,
    attempt           INTEGER NOT NULL DEFAULT 0,
    max_attempts      INTEGER NOT NULL DEFAULT 1,
    started_at        TEXT,
    ended_at          TEXT,
    elapsed_seconds   REAL,
    command           TEXT,
    tool_version      TEXT,
    database_version  TEXT,
    config_hash       TEXT,
    input_ids         TEXT,
    output_paths      TEXT,
    validation_state  TEXT,
    validation_detail TEXT,
    validation_json   TEXT,
    error             TEXT,
    host              TEXT,
    cpu_count         INTEGER,
    memory_mb         INTEGER,
    updated_at        TEXT NOT NULL,
    UNIQUE (run_key, stage, subject)
);

CREATE TABLE IF NOT EXISTS attempt_log (
    id               INTEGER PRIMARY KEY,
    execution_id     INTEGER NOT NULL REFERENCES execution(id) ON DELETE CASCADE,
    attempt          INTEGER NOT NULL,
    state            TEXT NOT NULL,
    failure_kind     TEXT,
    detail           TEXT,
    exit_code        INTEGER,
    elapsed_seconds  REAL,
    command          TEXT,
    log_path         TEXT,
    recorded_at      TEXT NOT NULL,
    UNIQUE (execution_id, attempt)
);

CREATE INDEX IF NOT EXISTS idx_execution_run_state
    ON execution (run_key, state);
CREATE INDEX IF NOT EXISTS idx_attempt_execution
    ON attempt_log (execution_id, attempt);
"""

#: One row's worth of persisted state, as a plain mapping.
ExecutionRow = Dict[str, Any]


def host_info() -> Dict[str, Any]:
    """Static capacity snapshot for the machine this run is on."""
    memory_mb: Optional[int] = None
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        memory_mb = int(pages * page_size / (1024 * 1024))
    except (ValueError, OSError, AttributeError):
        memory_mb = None
    return {
        "host": socket.gethostname(),
        "cpu_count": os.cpu_count(),
        "memory_mb": memory_mb,
        "platform": platform.platform(),
    }


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _synchronized(method):
    """Serialise a method against every other statement on this store.

    A decorator rather than a ``with`` block inside each method, so the
    bodies stay exactly as written and no body has to be re-indented.
    """

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        with self._lock:
            return method(self, *args, **kwargs)

    return wrapper


class ExecutionStore:
    """SQLite-backed execution state.

    One row per ``(run_key, stage, subject)``; the full attempt history
    alongside it. Opens the database on construction, so the schema is
    present whether or not anything is written.
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        if str(self.path) != ":memory:":
            self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._conn = sqlite3.connect(
            str(self.path), check_same_thread=False, timeout=30.0,
        )
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA foreign_keys=ON")
        self._conn.execute("PRAGMA busy_timeout=30000")
        self._conn.executescript(SCHEMA)
        self._conn.commit()

    @property
    def schema_version(self) -> int:
        return SCHEMA_VERSION

    def close(self) -> None:
        with self._lock:
            try:
                self._conn.close()
            except sqlite3.Error:
                pass

    def __enter__(self) -> "ExecutionStore":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # -- reads ------------------------------------------------------------

    @_synchronized
    def get(self, run_key: str, stage: str, subject: str = "") -> Optional[ExecutionRow]:
        cur = self._conn.execute(
            "SELECT * FROM execution WHERE run_key=? AND stage=? AND subject=?",
            (run_key, stage, subject),
        )
        row = cur.fetchone()
        return dict(row) if row else None

    @_synchronized
    def list_for_run(self, run_key: str) -> List[ExecutionRow]:
        cur = self._conn.execute(
            "SELECT * FROM execution WHERE run_key=? ORDER BY stage, subject", (run_key,)
        )
        return [dict(r) for r in cur.fetchall()]

    @_synchronized
    def run_keys(self) -> List[str]:
        """Run keys this store actually holds rows for, sorted.

        Exists so a reader can tell "this run has not started" apart from "you
        are looking at the wrong key". Without it a `--run` typo, or a default
        that disagrees with the overlay, renders as sixteen PENDING stages -
        indistinguishable from a run that has not started, and so a quiet wrong
        answer rather than an error.
        """
        cur = self._conn.execute("SELECT DISTINCT run_key FROM execution")
        return sorted(r["run_key"] for r in cur.fetchall())

    @_synchronized
    def attempts(self, run_key: str, stage: str, subject: str = "") -> List[AttemptRecord]:
        cur = self._conn.execute(
            """SELECT a.* FROM attempt_log a
               JOIN execution e ON e.id = a.execution_id
               WHERE e.run_key=? AND e.stage=? AND e.subject=?
               ORDER BY a.attempt""",
            (run_key, stage, subject),
        )
        return [
            AttemptRecord(
                attempt=int(r["attempt"]),
                state=StageState.coerce(r["state"]),
                failure_kind=(
                    FailureKind(r["failure_kind"]) if r["failure_kind"] else None
                ),
                detail=r["detail"] or "",
                exit_code=r["exit_code"],
                elapsed_seconds=float(r["elapsed_seconds"] or 0.0),
                command=tuple(json.loads(r["command"])) if r["command"] else (),
            )
            for r in cur.fetchall()
        ]

    @_synchronized
    def counts_by_state(self, run_key: str) -> Dict[str, int]:
        cur = self._conn.execute(
            "SELECT state, COUNT(*) AS n FROM execution WHERE run_key=? GROUP BY state",
            (run_key,),
        )
        return {r["state"]: int(r["n"]) for r in cur.fetchall()}

    @_synchronized
    def log_path(self, run_key: str, stage: str, subject: str = "") -> Optional[str]:
        """The log file recorded for a task, or ``None``.

        The ``execution`` table has no log column: the runner passes the
        path to :meth:`record_attempt`, which is where it is actually
        recorded. So this reads the most recent attempt that carried one,
        rather than looking for a column that was never written.
        """
        found = self._conn.execute(
            """SELECT a.log_path FROM attempt_log a
               JOIN execution e ON e.id = a.execution_id
               WHERE e.run_key = ? AND e.stage = ? AND e.subject = ?
                 AND a.log_path IS NOT NULL AND a.log_path != ''
               ORDER BY a.attempt DESC LIMIT 1""",
            (run_key, stage, subject),
        ).fetchone()
        return str(found["log_path"]) if found is not None else None

    # -- writes -----------------------------------------------------------

    @_synchronized
    def record(
        self,
        run_key: str,
        stage: str,
        subject: str = "",
        *,
        state: StageState,
        attempt: int = 0,
        max_attempts: int = 1,
        started_at: Optional[str] = None,
        ended_at: Optional[str] = None,
        elapsed_seconds: Optional[float] = None,
        command: Sequence[str] = (),
        tool_version: Optional[str] = None,
        database_version: Optional[str] = None,
        config_hash: Optional[str] = None,
        input_ids: Sequence[str] = (),
        output_paths: Sequence[Path] = (),
        validation_state: Optional[str] = None,
        validation_detail: Optional[str] = None,
        validation_json: Optional[str] = None,
        failure_kind: Optional[FailureKind] = None,
        error: Optional[str] = None,
        include_host: bool = True,
    ) -> ExecutionRow:
        """Insert or update the row for ``(run_key, stage, subject)``.

        Upsert rather than insert, so a re-run of the same task updates in
        place and the attempt history accumulates instead of forking.
        """
        host = host_info() if include_host else {}
        payload = {
            "run_key": run_key,
            "stage": stage,
            "subject": subject,
            "state": state.value,
            "failure_kind": failure_kind.value if failure_kind else None,
            "attempt": int(attempt),
            "max_attempts": int(max_attempts),
            "started_at": started_at,
            "ended_at": ended_at,
            "elapsed_seconds": elapsed_seconds,
            "command": json.dumps(list(command)) if command else None,
            "tool_version": tool_version,
            "database_version": database_version,
            "config_hash": config_hash,
            "input_ids": json.dumps(list(input_ids)) if input_ids else None,
            "output_paths": (
                json.dumps([str(p) for p in output_paths]) if output_paths else None
            ),
            "validation_state": validation_state,
            "validation_detail": validation_detail,
            "validation_json": validation_json,
            "error": error,
            "host": host.get("host"),
            "cpu_count": host.get("cpu_count"),
            "memory_mb": host.get("memory_mb"),
            "updated_at": _now(),
        }
        columns = ", ".join(payload)
        placeholders = ", ".join(f":{k}" for k in payload)
        updates = ", ".join(
            f"{k}=excluded.{k}" for k in payload
            if k not in ("run_key", "stage", "subject")
        )
        sql = (
            f"INSERT INTO execution ({columns}) VALUES ({placeholders}) "
            f"ON CONFLICT(run_key, stage, subject) DO UPDATE SET {updates}"
        )
        self._conn.execute(sql, payload)
        self._conn.commit()
        row = self.get(run_key, stage, subject)
        assert row is not None
        return row

    @_synchronized
    def next_attempt(self, run_key: str, stage: str, subject: str = "") -> int:
        """The next attempt number for a task, across every invocation.

        Attempt numbers are monotonic for the lifetime of the ``(run, stage,
        subject)`` triple, not per call. A task that is interrupted and later
        resumed is a *continuation*, so its history must read 1, 2, 3, 4
        rather than restarting and overwriting earlier evidence.
        """
        row = self.get(run_key, stage, subject)
        if row is None:
            raise PipelineError(
                "Cannot allocate an attempt for an unrecorded task",
                run_key=run_key, stage=stage, subject=subject,
            )
        found = self._conn.execute(
            "SELECT COALESCE(MAX(attempt), 0) + 1 AS n FROM attempt_log WHERE execution_id = ?",
            (row["id"],),
        ).fetchone()
        return int(found["n"])

    @_synchronized
    def record_attempt(
        self,
        run_key: str,
        stage: str,
        subject: str,
        record: AttemptRecord,
        log_path: Optional[Path] = None,
    ) -> AttemptRecord:
        """Append one attempt to the history. Never overwrites.

        The attempt number is allocated here rather than trusted from the
        caller, so that re-running a task appends to its history instead of
        replacing the previous run's evidence.
        """
        row = self.get(run_key, stage, subject)
        if row is None:
            raise PipelineError(
                "Cannot log an attempt for an unrecorded task",
                run_key=run_key, stage=stage, subject=subject,
            )
        attempt = self.next_attempt(run_key, stage, subject)
        from dataclasses import replace

        stored = replace(record, attempt=attempt)
        self._conn.execute(
            """INSERT INTO attempt_log
               (execution_id, attempt, state, failure_kind, detail, exit_code,
                elapsed_seconds, command, log_path, recorded_at)
               VALUES (?,?,?,?,?,?,?,?,?,?)""",
            (
                row["id"], attempt, record.state.value,
                record.failure_kind.value if record.failure_kind else None,
                record.detail, record.exit_code, record.elapsed_seconds,
                json.dumps(list(record.command)) if record.command else None,
                str(log_path) if log_path else None, _now(),
            ),
        )
        self._conn.commit()
        return stored

    @_synchronized
    def set_state(
        self,
        run_key: str,
        stage: str,
        subject: str,
        state: StageState,
        *,
        failure_kind: Optional[FailureKind] = None,
        error: Optional[str] = None,
    ) -> None:
        """Transition a task's state without touching its history."""
        self._conn.execute(
            """UPDATE execution SET state=?, failure_kind=?, error=?, updated_at=?
               WHERE run_key=? AND stage=? AND subject=?""",
            (
                state.value,
                failure_kind.value if failure_kind else None,
                error,
                _now(), run_key, stage, subject,
            ),
        )
        self._conn.commit()
