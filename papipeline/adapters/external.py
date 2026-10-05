"""External tool detection and adapters.

No bioinformatics tool is a hard import dependency of this package. Each
adapter declares what it needs; :func:`detect_tools` reports what is actually
present; and a stage refuses to run a tool it needs rather than silently
degrading.

This is also where the "no implicit database updates" rule lives: an adapter
is constructed with an explicit database version, and a request to update a
database is only honoured when the caller passes ``allow_update=True``, which
in turn is gated on configuration.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from ..errors import ToolExecutionError, ToolNotAvailableError
from ..logging_utils import get_logger

LOGGER = get_logger("adapters")


@dataclass(frozen=True)
class ToolStatus:
    """Availability of one external tool."""

    name: str
    executable: Optional[str]
    version: Optional[str]
    available: bool

    def to_row(self) -> Dict[str, Optional[str]]:
        return {
            "tool": self.name,
            "executable": self.executable,
            "version": self.version,
            "available": "true" if self.available else "false",
        }


#: The version recorded for a tool whose version was NOT read by executing it.
#:
#: A presence sweep may not run a binary just to fill this column in: the sweep
#: covers every tool the project has ever heard of, and executing all of them
#: would run tools this project is forbidden to run (Bakta and Panaroo). So
#: presence is reported and the version is reported as unknown, which is the
#: honest value - a tool that was found on PATH and never asked.
#:
#: Chosen over `None` because `None` in a provenance row reads as "nothing to
#: report" and is indistinguishable from a tool that was absent. UNKNOWN says
#: the tool is there and its version was deliberately not measured.
UNPROBED_VERSION = "UNKNOWN"


#: Tool name -> the version flag to try.
VERSION_FLAGS: Mapping[str, Sequence[str]] = {
    "bakta": ("--version",),
    "amrfinder": ("--version",),
    "mlst": ("--version",),
    "rgi": ("--version",),
    "panaroo": ("--version",),
    "snp-sites": ("--version",),
    "iqtree2": ("--version",),
    "iqtree": ("--version",),
    "pyseer": ("--version",),
    "makeblastdb": ("-version",),
    "blastn": ("-version",),
    "minimap2": ("--version",),
    "seqkit": ("--version",),
    "snakemake": ("--version",),
    "datasets": ("--version",),
}


def probe_version(name: str, executable: str, timeout: int = 20) -> Optional[str]:
    """Return the first line of ``<tool> <flag>``, or ``None`` on any failure."""
    for flag in VERSION_FLAGS.get(name, ("--version",)):
        try:
            completed = subprocess.run(
                [executable, flag],
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except (subprocess.TimeoutExpired, OSError):
            continue
        output = (completed.stdout or completed.stderr or "").strip()
        if output:
            return output.splitlines()[0].strip()
    return None


def detect_tools(
    names: Sequence[str], *, execute: bool = False
) -> Dict[str, ToolStatus]:
    """Locate each named tool. Presence only, unless ``execute`` says otherwise.

    **The default runs NOTHING.** ``shutil.which`` resolves a name without
    executing the binary; the only subprocess this function can start is the one
    behind ``execute=True``, and the default is the safe one on purpose - a
    sweep is a plausible mistake and an accidental ``--version`` on Bakta is the
    thing this signature exists to make impossible by omission.

    Every status still reports presence either way. Without ``execute`` the
    version is :data:`UNPROBED_VERSION` rather than invented.

    ``execute=True`` is reserved for a caller that is about to INVOKE that tool
    - a stage that has decided it needs it, and will run it, in this call.
    :func:`require_tool` callers are in that position.
    """
    statuses: Dict[str, ToolStatus] = {}
    for name in names:
        executable = shutil.which(name)
        if executable is None:
            statuses[name] = ToolStatus(name, None, None, False)
            LOGGER.debug("Tool not found on PATH: %s", name)
            continue
        version = probe_version(name, executable) if execute else UNPROBED_VERSION
        statuses[name] = ToolStatus(name, executable, version, True)
        LOGGER.debug("Tool found: %s -> %s (%s)", name, executable, version)
    return statuses


def probe_on_demand(
    name: str, statuses: Dict[str, ToolStatus]
) -> ToolStatus:
    """Read ONE tool's version by executing it, updating ``statuses`` in place.

    **The sanctioned way to run a tool's version flag.** The caller must be a
    stage that has already decided to INVOKE that tool. Every other reader takes
    the presence sweep's :data:`UNPROBED_VERSION` and runs nothing.

    In place, deliberately: the status map is what ``run_manifest.json`` records
    under ``tools_detected``, and after a stage has genuinely run a tool its
    measured version is the truthful thing to record there. Nothing else in the
    map changes, so a reader can still tell measured from unmeasured by
    comparing against :data:`UNPROBED_VERSION`.

    Returns the updated status. A tool that is not present is left exactly as it
    was: there is nothing to execute, and executing anything to find that out
    would be the defect this function exists to prevent.
    """
    status = statuses.get(name)
    if status is None or not status.available:
        return ToolStatus(name, None, None, False)
    version = probe_version(name, status.executable or name)
    statuses[name] = ToolStatus(name, status.executable, version, True)
    LOGGER.debug("Probed on demand: %s -> %s (%s)", name, status.executable, version)
    return statuses[name]


def detect_all() -> Dict[str, ToolStatus]:
    """Report EVERY tool the pipeline knows about, WITHOUT executing any of them.

    **Presence only, by ``shutil.which``, for every tool in
    :data:`VERSION_FLAGS`.** No subprocess is started.

    Why this is not a version sweep: this function is reached from
    ``run_pipeline`` before stage 1, over a list of every tool the project has
    ever heard of, and asking each one for its version means running each one.
    That included ``bakta --version`` and ``panaroo --version`` on every run,
    including runs configured with ``annotation.reuse_tool_output: require``
    precisely so Bakta would never be executed - a REAL 10-isolate run's tripwire
    logged both (``GUARDS.md``: pids 93671, 93679). Reuse mode decides whether a
    tool is INVOKED; a version probe does not consult it and so cannot honour it.

    Nothing is dropped from the report: every tool in :data:`VERSION_FLAGS` still
    appears, with its executable and availability, and its version is
    :data:`UNPROBED_VERSION` rather than a guess. A stage that is about to
    invoke a tool reads that tool's version with :func:`probe_on_demand`, or
    detects a named set with ``detect_tools(..., execute=True)``.
    """
    return detect_tools(sorted(VERSION_FLAGS), execute=False)


def require_tool(statuses: Mapping[str, ToolStatus], name: str, stage: str) -> ToolStatus:
    """Return a tool status or raise a clear, actionable error.

    Raises:
        ToolNotAvailableError: The tool is not installed. The message names
            the environment file that provides it.
    """
    status = statuses.get(name)
    if status is None or not status.available:
        raise ToolNotAvailableError(
            "Required tool is not installed",
            tool=name,
            stage=stage,
            hint=(
                "Install the pinned environment first: "
                "micromamba env create -f environment/environment.yml"
            ),
        )
    return status


def require_optional_tool(
    statuses: Mapping[str, ToolStatus], name: str, stage: str
) -> Optional[ToolStatus]:
    """Like :func:`require_tool` but returns ``None`` when absent."""
    status = statuses.get(name)
    if status is None or not status.available:
        LOGGER.info("Optional tool not available, stage %s will be skipped: %s", stage, name)
        return None
    return status


@dataclass
class CommandResult:
    """Captured result of an external command."""

    command: Sequence[str]
    returncode: int
    stdout: str
    stderr: str
    duration_s: float

    def to_row(self) -> Dict[str, object]:
        return {
            "command": " ".join(self.command),
            "returncode": self.returncode,
            "duration_s": round(self.duration_s, 3),
        }


@dataclass
class ToolAdapter:
    """Base class for external tool adapters.

    Subclasses set :attr:`tool_name` and implement :meth:`run`. The adapter
    never mutates a database: :meth:`ensure_database` is opt-in and refuses
    unless ``allow_update`` was set by configuration.
    """

    tool_name: str = "unknown"
    database: str = "none"
    database_version: str = "UNPINNED"

    def __init__(
        self,
        status: ToolStatus,
        database_version: Optional[str] = None,
        allow_database_update: bool = False,
        workdir: Optional[Path] = None,
    ) -> None:
        self.status = status
        if database_version:
            self.database_version = database_version
        self.allow_database_update = allow_database_update
        self.workdir = Path(workdir) if workdir else None
        self.command_log: List[CommandResult] = []

    # -- database handling ------------------------------------------------
    def ensure_database(self, requested_version: Optional[str] = None) -> str:
        """Return the pinned database version.

        Raises:
            ToolExecutionError: An update was requested but is not permitted.
                Databases are provisioned out of band, never mid-run.
        """
        if requested_version and requested_version != self.database_version:
            if not self.allow_database_update:
                raise ToolExecutionError(
                    "Database version change requested during an analysis run",
                    tool=self.tool_name,
                    pinned=self.database_version,
                    requested=requested_version,
                    hint=(
                        "Provision the new database out of band and update "
                        "config/references.tsv"
                    ),
                )
            raise ToolExecutionError(
                "In-run database updates are not implemented",
                tool=self.tool_name,
                hint="Update config/references.tsv and re-run",
            )
        return self.database_version

    # -- execution --------------------------------------------------------
    def _execute(
        self,
        command: Sequence[str],
        cwd: Optional[Path] = None,
        check: bool = True,
    ) -> CommandResult:
        """Run a command, log it, and return the captured result."""
        import time

        LOGGER.info("exec: %s", " ".join(command))
        start = time.monotonic()
        try:
            completed = subprocess.run(
                list(command),
                cwd=str(cwd) if cwd else None,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise ToolExecutionError(
                "Failed to launch external tool",
                tool=self.tool_name,
                command=" ".join(command),
                error=str(exc),
            ) from exc
        duration = time.monotonic() - start

        result = CommandResult(
            command=list(command),
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            duration_s=duration,
        )
        self.command_log.append(result)

        if check and result.returncode != 0:
            raise ToolExecutionError(
                "External tool exited non-zero",
                tool=self.tool_name,
                command=" ".join(command),
                returncode=result.returncode,
                stderr=result.stderr[-2000:],
            )
        return result

    def provenance(self) -> Dict[str, Optional[str]]:
        """Provenance fields recorded for every run."""
        return {
            "tool": self.tool_name,
            "tool_version": self.status.version,
            "database": self.database,
            "database_version": self.database_version,
            "commands": "; ".join(" ".join(r.command) for r in self.command_log) or None,
        }
