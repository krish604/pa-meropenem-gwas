"""Derive per-stage state from a run's own logs when no manifest exists.

**Why this module exists.** `run_manifest.json` is written only on success
(`papipeline/run.py:_write_run_manifest`). A REAL run that fails at stage 7
therefore leaves **no** manifest, and the delivery bundle that records it
(`artifacts/logs/full_run.log`, `PER_STAGE_OUTCOMES.md`) carries the per-stage
record instead. A dashboard that reads only the manifest would render every
stage `not_run` with the pre-run-manifest sentence, which is false: seven
stages completed.

This module parses the artefacts the failed run *did* leave:

* `artifacts/logs/full_run.log` — the `--- stage <name> ---` start markers and
  the `Stage <name> ended FAILED: <verbatim>` line;
* `artifacts/logs/stage_<name>.log` — one file per stage the run reached;
* the stage table's presence under `artifacts/stage_tables/`.

It returns the same `{states, reasons}` shape `badges.classify` already reads
from `run_manifest.stages` / `stages_skipped`, so no page learns a second
vocabulary. A stage that ran but whose table is missing is left `completed`
here and downgraded to `not_assessed` by `classify`, which supplies the
`TableRead` reason.

**Nothing is invented.** A stage with no start marker and no stage log is
`not_run`; when a stage failed, the stages after it name that stage as the
blocker rather than receiving a generic absence.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from papipeline.execution.contracts import STAGE_TABLES
from papipeline.run import STAGE_ORDER

#: The three lines this parser reads. Anchored enough to ignore the many other
#: log lines that mention a stage name (`Stage 4 provenance: ...`,
#: `Stage 7: screened ...`).
_STAGE_MARKER = re.compile(r"--- stage ([A-Za-z0-9_]+) ---")
_STAGE_FAILED = re.compile(r"Stage ([A-Za-z0-9_]+) ended FAILED: (.*)")
_RUN_HEADER = re.compile(r"=== (\w+) mode \| (\d+) samples \| antibiotic=(\w+)")
_OPRD_LOCUS = re.compile(
    r"oprD locus (\S+): (resolved|refused:[A-Za-z0-9_]+) \(coverage ([\d.]+)%\)"
)

#: The delivery layout's log directory, relative to the bundle root.
DELIVERY_LOG_DIR = Path("artifacts") / "logs"
FULL_RUN_LOG = "full_run.log"
STAGE_LOG_PREFIX = "stage_"
STAGE_LOG_SUFFIX = ".log"


@dataclass(frozen=True)
class DerivedStageState:
    """Per-stage state derived from a run's logs, not from a manifest.

    ``records`` is True when there was evidence to derive from at all, so a
    caller can tell "derived and this stage did not run" from "nothing to
    derive from". ``states`` and ``reasons`` carry the same shape the manifest
    reader supplies to `badges.classify`.
    """

    records: bool
    states: Mapping[str, str] = field(default_factory=dict)
    reasons: Mapping[str, str] = field(default_factory=dict)
    ran: Tuple[str, ...] = ()
    source: str = ""
    n_stage_logs: int = 0
    run_mode: Optional[str] = None
    antibiotic: Optional[str] = None
    n_samples: Optional[int] = None

    def as_dict(self) -> Dict[str, Any]:
        return {
            "records": self.records,
            "states": dict(self.states),
            "reasons": dict(self.reasons),
            "ran": list(self.ran),
            "source": self.source,
            "n_stage_logs": self.n_stage_logs,
            "run_mode": self.run_mode,
            "antibiotic": self.antibiotic,
            "n_samples": self.n_samples,
        }


def derive(root: Path, stage_dir: Path, *, layout: str = "") -> Optional[DerivedStageState]:
    """Parse a delivery bundle's logs into a per-stage record.

    Returns None unless the delivery layout's `full_run.log` is present, so a
    live tree or an assumed-layout bundle is never given a derived record it
    does not have.
    """
    root = Path(root)
    stage_dir = Path(stage_dir)
    if layout != "delivery":
        return None
    full = root / DELIVERY_LOG_DIR / FULL_RUN_LOG
    if not full.is_file():
        return None

    try:
        text = full.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None

    started: set = set()
    failed: Dict[str, str] = {}
    run_mode: Optional[str] = None
    antibiotic: Optional[str] = None
    n_samples: Optional[int] = None
    for line in text.splitlines():
        marker = _STAGE_MARKER.search(line)
        if marker:
            started.add(marker.group(1))
        failure = _STAGE_FAILED.search(line)
        if failure:
            failed[failure.group(1)] = failure.group(2).strip()
        header = _RUN_HEADER.search(line)
        if header:
            run_mode = header.group(1)
            n_samples = int(header.group(2))
            antibiotic = header.group(3)

    logs_dir = full.parent
    reached: set = set()
    n_stage_logs = 0
    if logs_dir.is_dir():
        for path in sorted(logs_dir.glob(f"{STAGE_LOG_PREFIX}*{STAGE_LOG_SUFFIX}")):
            name = path.name[len(STAGE_LOG_PREFIX): -len(STAGE_LOG_SUFFIX)]
            if name in STAGE_TABLES:
                reached.add(name)
                n_stage_logs += 1

    states: Dict[str, str] = {}
    reasons: Dict[str, str] = {}
    blocker_stage: Optional[str] = None
    blocker_reason: Optional[str] = None
    for stage in STAGE_ORDER:
        if stage in failed:
            blocker_stage = stage
            blocker_reason = failed[stage]

    for stage in STAGE_ORDER:
        if stage in failed:
            states[stage] = "failed"
            # Verbatim, from the log's own `ended FAILED:` line.
            reasons[stage] = failed[stage]
        elif stage in started or stage in reached:
            # Ran. `classify` downgrades to `not_assessed` when the principal
            # table is missing, and supplies the TableRead reason itself.
            states[stage] = "completed"
        else:
            states[stage] = "not_run"
            if blocker_stage is not None and blocker_reason:
                # The failure line already names the stage number
                # (`ToolNotAvailableError: Stage 7 (pangenome) ...`), so it is
                # carried verbatim rather than re-derived from STAGE_ORDER's
                # 0-based index, which numbers pangenome 8.
                reasons[stage] = (
                    f"not run: blocked by the failure of stage {blocker_stage}: "
                    f"{blocker_reason}"
                )
            else:
                reasons[stage] = (
                    "not run: artifacts/logs/ holds no `--- stage "
                    f"{stage} ---` start marker and no stage_{stage}.log"
                )

    return DerivedStageState(
        records=True,
        states=states,
        reasons=reasons,
        ran=tuple(s for s in STAGE_ORDER if s in started or s in reached),
        source=str(full),
        n_stage_logs=n_stage_logs,
        run_mode=run_mode,
        antibiotic=antibiotic,
        n_samples=n_samples,
    )


def locus_coverage(root: Path) -> Tuple[List[Dict[str, Any]], Optional[str]]:
    """The separate oprD locus-coverage gate, parsed from `full_run.log`.

    This is a **different instrument** from the tblastn structural call: the
    coverage gate asks "did an alignment span this locus", the structural call
    asks "is the ORF intact in the assembly". The run reports both; this
    function surfaces the first so the two are never silently merged.

    Returns `(items, reason)`; `items` is empty with a non-empty `reason` when
    the log is absent or holds no gate lines.
    """
    root = Path(root)
    full = root / DELIVERY_LOG_DIR / FULL_RUN_LOG
    if not full.is_file():
        return [], (
            "no locus-coverage gate lines were read: "
            f"{full} is not present in this source"
        )
    try:
        text = full.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        return [], f"{full} is unreadable: {exc}"
    items: List[Dict[str, Any]] = []
    for line in text.splitlines():
        match = _OPRD_LOCUS.search(line)
        if not match:
            continue
        items.append(
            {
                "sample_id": match.group(1),
                "verdict": match.group(2),
                "coverage_pct": float(match.group(3)),
                "display_state": (
                    "not_assessed"
                    if match.group(2).startswith("refused:")
                    else "intact"
                ),
            }
        )
    if not items:
        return [], f"{full} holds no `oprD locus <id>: <verdict> (coverage N%)` line"
    return items, None


__all__ = [
    "DELIVERY_LOG_DIR",
    "FULL_RUN_LOG",
    "DerivedStageState",
    "derive",
    "locus_coverage",
]
