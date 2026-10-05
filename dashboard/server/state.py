"""UI-D1 — the state directory is the only write location.

The results root is read-only. Everything the dashboard needs to remember goes
here, outside the results root by construction, and everything here is derived:
deleting the whole directory costs one index rebuild and no information.

    ~/.pa_dashboard/
      index/          <stage>.<key16>.idx.json    UI-D5, disposable
      logs/           dashboard.log
      session.json    opened root, opened-at, transport

`PA_DASH_STATE_DIR` overrides the location, so a test never writes into the
developer's home directory.

There is deliberately no API here that writes anywhere else. The containment
check in `security.py` (UI-D4) would reject a write path under the opened
results root even if a caller invented one.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Optional

#: Overrides the state directory. Read once per process.
STATE_DIR_ENV = "PA_DASH_STATE_DIR"

#: The default location when the environment says nothing.
DEFAULT_STATE_DIR = Path("~/.pa_dashboard")

#: Subdirectory names. Declared here so the layout in DESIGN §10 is written
#: down once rather than spelled at four call sites.
INDEX_DIR = "index"
LOGS_DIR = "logs"
SESSION_FILE = "session.json"


def default_state_dir() -> Path:
    """The state directory, `PA_DASH_STATE_DIR` honoured.

    `~` is expanded here rather than left to the filesystem so the logged path
    is the real one.
    """
    override = os.environ.get(STATE_DIR_ENV)
    if override:
        return Path(override).expanduser()
    return DEFAULT_STATE_DIR.expanduser()


@dataclass(frozen=True)
class StateDir:
    """One opened state directory.

    Resolved once, at construction, and `resolve()`d so a symlinked state
    directory is followed rather than written through (DESIGN §11, last row).
    """

    root: Path

    @classmethod
    def open(cls, root: Optional[Path] = None) -> "StateDir":
        chosen = Path(root).expanduser() if root is not None else default_state_dir()
        # `resolve()` on a not-yet-existing path is fine on Python 3.6+ and
        # resolves the existing prefix, which is what we want: the parent
        # directory may not exist yet.
        resolved = chosen.resolve()
        (resolved / INDEX_DIR).mkdir(parents=True, exist_ok=True)
        (resolved / LOGS_DIR).mkdir(parents=True, exist_ok=True)
        return cls(root=resolved)

    # -- layout ---------------------------------------------------------
    @property
    def index_dir(self) -> Path:
        return self.root / INDEX_DIR

    @property
    def logs_dir(self) -> Path:
        return self.root / LOGS_DIR

    @property
    def session_path(self) -> Path:
        return self.root / SESSION_FILE

    @property
    def log_path(self) -> Path:
        return self.logs_dir / "dashboard.log"

    # -- the only writes -----------------------------------------------
    def write_json(self, path: Path, payload: Mapping[str, Any]) -> Path:
        """Write JSON atomically, into the state directory only.

        Temp file in the same directory then `os.replace`, so a kill mid-write
        leaves the previous file rather than a truncated one. This is the
        "temporary file in the same directory, then `os.replace`" rule from
        UI-D5, applied to `session.json` too.
        """
        resolved = Path(path)
        if not self._contains(resolved):
            # A defence in depth, not the primary control: the primary control
            # is that no caller builds a state path from anything the client
            # sent. This turns a future mistake into a refusal.
            raise ValueError(
                f"refusing to write {resolved}: outside the state directory "
                f"{self.root}. UI-D1 permits writes to the state directory and "
                f"nowhere else."
            )
        resolved.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(resolved.parent),
            prefix=f".{resolved.name}.",
            suffix=".tmp",
            delete=False,
        )
        tmp = Path(handle.name)
        try:
            with handle:
                json.dump(payload, handle, indent=2, sort_keys=True, default=str)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(str(tmp), str(resolved))
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        return resolved

    def _contains(self, candidate: Path) -> bool:
        try:
            candidate = Path(candidate).resolve()
        except OSError:
            return False
        return candidate == self.root or self.root in candidate.parents

    # -- session.json ---------------------------------------------------
    def write_session(
        self,
        *,
        root: Optional[Path],
        kind: str,
        manifest_writer: str,
        transport: str = "sse",
    ) -> Path:
        """Record what is open. Rewritten when the root changes.

        Derived, like everything else here: deleting it loses no information,
        because the next start writes it again.
        """
        payload = {
            "opened_root": str(root) if root is not None else None,
            "opened_at_utc": datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat()
            .replace("+00:00", "Z"),
            "source_kind": kind,
            "manifest_writer": manifest_writer,
            "transport": transport,
        }
        return self.write_json(self.session_path, payload)

    def read_session(self) -> Optional[Mapping[str, Any]]:
        """The last session record, or None if there is not one yet."""
        if not self.session_path.exists():
            return None
        try:
            payload = json.loads(self.session_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None
        return payload if isinstance(payload, dict) else None


__all__ = [
    "DEFAULT_STATE_DIR",
    "INDEX_DIR",
    "LOGS_DIR",
    "SESSION_FILE",
    "STATE_DIR_ENV",
    "StateDir",
    "default_state_dir",
]