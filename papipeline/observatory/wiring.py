"""Turning the observatory on inside a real pipeline run.

The observatory is strictly optional. When it is off, this module hands back
a disabled handle whose methods do nothing, and the pipeline runs exactly as
it did before - which is the property that makes it safe to leave in.

When it is on, the handle owns three things the pipeline would otherwise not
have: the execution store it writes to, the event bus it publishes
transitions to, and the per-stage :class:`~papipeline.execution.TaskContext`
that carries a stage's real inputs, outputs and contract.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from ..execution.contracts import stage_spec, table_path
from ..execution.state import TaskEvent
from ..execution.store import ExecutionStore
from ..execution.validation import OutputSpec
from .events import BUS, ExecutionEvent, EventBus, TASK_RESUMED, bus_sink

DEFAULT_DB = Path.home() / ".local" / "share" / "papipeline" / "observatory.db"
DEFAULT_RUN_KEY = "pipeline"


@dataclass
class ObservatorySettings:
    """Whether to observe, and against what."""

    enabled: bool = False
    db_path: Path = DEFAULT_DB
    run_key: str = DEFAULT_RUN_KEY

    @classmethod
    def from_config(
        cls,
        config: Any,
        *,
        db_path: Optional[Path] = None,
        run_key: Optional[str] = None,
        enable: Optional[bool] = None,
    ) -> "ObservatorySettings":
        """Read the setting from configuration, with an explicit override.

        Configuration wins over nothing and is overridden by an explicit
        argument, so ``--observatory`` on the command line beats
        ``runtime.observatory.enabled: false`` in the file.
        """
        runtime = dict(getattr(config, "runtime", {}) or {})
        section = runtime.get("observatory") or {}
        if not isinstance(section, dict):
            section = {}
        on = enable
        if on is None:
            on = bool(section.get("enabled", False))
        # `expanduser()` because the overlays write `~/.local/share/...`, and a
        # bare `Path()` treats that as a *relative* directory named "~". The run
        # then created `./~/.local/share/papipeline/observatory.db` inside the
        # repository while `/api/health` looked in the real home and reported
        # `database_present: false` - the store existed, in the wrong place,
        # because the writer and the reader resolved the same string differently.
        return cls(
            enabled=bool(on),
            db_path=Path(db_path or section.get("db") or DEFAULT_DB).expanduser(),
            run_key=str(run_key or section.get("run_key") or DEFAULT_RUN_KEY),
        )


def config_hash(config: Any) -> str:
    """A stable hash of the effective configuration.

    Recorded on every task so that a changed configuration invalidates a
    recorded success, which is what makes resume honest. Serialised
    defensively because config objects hold paths and enums.
    """

    def default(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, (set, frozenset, tuple)):
            return sorted(str(v) for v in value)
        if hasattr(value, "__dataclass_fields__"):
            return {k: default(v) for k, v in vars(value).items()}
        if hasattr(value, "value"):
            return value.value
        return str(value)

    try:
        payload = json.dumps(default(config), sort_keys=True, default=str)
    except (TypeError, ValueError):
        payload = repr(config)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


class Observatory:
    """The enabled handle. Use :meth:`create` rather than constructing it."""

    def __init__(
        self,
        settings: ObservatorySettings,
        *,
        store: Optional[ExecutionStore] = None,
        bus: Optional[EventBus] = None,
    ) -> None:
        self.settings = settings
        self._store = store
        # Own bus by default. Sharing the process-wide singleton would leak
        # one run's events into another run's stream whenever two observed
        # pipelines coexist, and would make a test suite order-dependent.
        self.bus = bus if bus is not None else EventBus()
        self._config_hash: Optional[str] = None
        self.resumed: list = []
        self.skipped: List[tuple] = []
        self.classification = None

    @property
    def enabled(self) -> bool:
        return self.settings.enabled

    @property
    def store(self) -> Optional[ExecutionStore]:
        return self._store

    @property
    def event_sink(self):
        return bus_sink(self.bus)

    def open(self, config: Any) -> "Observatory":
        """Open the store and freeze the configuration hash for this run."""
        if not self.enabled:
            return self
        self.settings.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._store = ExecutionStore(self.settings.db_path)
        self._config_hash = config_hash(config)
        return self

    def set_config_hash(self, value: str) -> None:
        self._config_hash = value

    @property
    def hash(self) -> str:
        return self._config_hash or ""

    def close(self) -> None:
        if self._store is not None:
            self._store.close()
            self._store = None

    def __enter__(self) -> "Observatory":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # ── per-task construction ──────────────────────────────────────

    def spec_for(
        self,
        stage: str,
        *,
        stage_dir: Path,
        n_samples: Optional[int] = None,
        extra_paths: Sequence[Path] = (),
    ) -> OutputSpec:
        """The real contract for a stage's principal table."""
        return stage_spec(
            stage, stage_dir=stage_dir, n_samples=n_samples,
            extra_paths=extra_paths,
        )

    def context(
        self,
        stage: str,
        *,
        subject: str = "",
        spec: Optional[OutputSpec] = None,
        stage_dir: Optional[Path] = None,
        n_samples: Optional[int] = None,
        input_ids: Sequence[str] = (),
        tool_version: Optional[str] = None,
        database_version: Optional[str] = None,
        log_dir: Optional[Path] = None,
        timeout: int = 3600,
    ):
        """A :class:`TaskContext` carrying this stage's real identity.

        The task identity is ``(run_key, stage, subject)`` - the same triple
        the store is keyed on - so it is stable across runs and resume
        remains possible. No random component is introduced, deliberately.
        """
        from ..execution import TaskContext

        if spec is None:
            if stage_dir is None:
                raise ValueError("pass spec or stage_dir")
            spec = self.spec_for(stage, stage_dir=stage_dir, n_samples=n_samples)
        log_path = None
        if log_dir is not None:
            suffix = f".{subject}" if subject else ""
            log_path = Path(log_dir) / f"{stage}{suffix}.log"
        return TaskContext(
            run_key=self.settings.run_key,
            stage=stage,
            subject=subject,
            spec=spec,
            tool_version=tool_version,
            database_version=database_version,
            config_hash=self._config_hash,
            input_ids=list(input_ids),
            timeout=timeout,
            log_path=log_path,
        )

    # ── resume ─────────────────────────────────────────────────────

    # ── run level ──────────────────────────────────────────────────

    def announce_run(
        self, classification, *, mode: str = "", stages: Sequence[str] = ()
    ) -> None:
        """Bracket a whole execution with the correct run-level event.

        A first execution and a deliberate re-run both announce
        ``RUN_STARTED`` - they execute the same work. Only a genuine resume
        announces ``RUN_RESUMED``, because only a resume reuses state.
        """
        if not self.enabled:
            return
        from .events import RUN_RESUMED, RUN_STARTED
        from ..execution.semantics import RunKind

        is_resume = classification.kind is RunKind.RESUME
        self.bus.emit(ExecutionEvent(
            name=RUN_RESUMED if is_resume else RUN_STARTED,
            stage="", subject="", run_key=self.settings.run_key,
            state="RUNNING", attempt=0,
            detail=classification.reason,
            payload={
                "kind": classification.kind.value,
                "prior_rows": classification.prior_rows,
                "prior_succeeded": classification.prior_succeeded,
                "prior_incomplete": classification.prior_incomplete,
                "mode": mode,
                "stages": list(stages),
            },
        ))

    def announce_run_finished(
        self, classification, *, executed: Sequence[str] = (),
        skipped: Sequence[str] = (), failed: Optional[str] = None,
        detail: str = "",
    ) -> None:
        """Close a run with ``RUN_COMPLETED`` or ``RUN_FAILED``."""
        if not self.enabled:
            return
        from .events import RUN_COMPLETED, RUN_FAILED
        from ..execution.semantics import RunKind

        name = RUN_FAILED if failed else RUN_COMPLETED
        self.bus.emit(ExecutionEvent(
            name=name, stage="", subject="", run_key=self.settings.run_key,
            state="FAILED" if failed else "SUCCEEDED", attempt=0,
            detail=failed or detail or "run finished",
            payload={
                "kind": classification.kind.value,
                "executed": list(executed),
                "skipped": list(skipped),
                "failed_stage": failed,
            },
        ))

    def announce_skipped(self, stage: str, subject: str = "", reason: str = "") -> None:
        """A task a resume validated and therefore did not re-execute.

        This is the only place ``TASK_RESUMED`` comes from. It is emitted
        once per task actually skipped, so the count of these events is the
        count of work a resume genuinely saved.
        """
        if not self.enabled or self._store is None:
            return
        from .events import TASK_RESUMED

        previous = self._store.get(self.settings.run_key, stage, subject)
        if previous is None or previous.get("state") != "SUCCEEDED":
            # Nothing to have resumed; do not claim otherwise.
            return
        self.bus.emit(ExecutionEvent(
            name=TASK_RESUMED,
            stage=stage, subject=subject, run_key=self.settings.run_key,
            state="SUCCEEDED",
            attempt=int(previous.get("attempt") or 0),
            detail=reason or "validated by a previous attempt in this run; not re-executed",
            payload={"previous_state": previous.get("state"),
                     "previous_ended_at": previous.get("ended_at"),
                     "previous_attempt": previous.get("attempt")},
        ))
        self.skipped.append((stage, subject))

    def mark_unblocked(self, stage: str, prerequisites: Sequence[str]) -> None:
        """Announce that a stage's declared prerequisites are satisfied.

        Emitted only when the named stages are genuinely recorded
        ``SUCCEEDED`` in this run, so it reflects the store rather than the
        order of a Python loop.
        """
        if not self.enabled or self._store is None:
            return
        satisfied = []
        for name in prerequisites:
            rows = [r for r in self._store.list_for_run(self.settings.run_key)
                    if r.get("stage") == name]
            if rows and all(r.get("state") == "SUCCEEDED" for r in rows):
                satisfied.append(name)
        if not satisfied:
            return
        self.bus.emit(ExecutionEvent(
            name="STAGE_UNBLOCKED",
            stage=stage,
            subject="",
            run_key=self.settings.run_key,
            state="PENDING",
            attempt=0,
            detail=f"prerequisites satisfied: {', '.join(satisfied)}",
            payload={"prerequisites": list(prerequisites)},
        ))


def create(
    config: Any,
    *,
    db_path: Optional[Path] = None,
    run_key: Optional[str] = None,
    enable: Optional[bool] = None,
) -> Observatory:
    """Build a handle, opened if it is enabled."""
    settings = ObservatorySettings.from_config(
        config, db_path=db_path, run_key=run_key, enable=enable)
    return Observatory(settings).open(config)
