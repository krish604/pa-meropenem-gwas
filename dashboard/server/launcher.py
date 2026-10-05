"""K4 — the launcher. Disabled unless `--allow-launch`, and allowlisted besides.

The default action is a **dry run**: `snakemake -n full_run`, which shows the
plan and starts nothing. Expensive tools are flagged by name — bakta, panaroo,
gubbins, iqtree, pyseer — with per-isolate timings from `RUNBOOK_900.md` **when
available**. That file does not exist anywhere yet (assumption A5), so the
timings degrade honestly to `timings not recorded` (UI-D2). They are never
invented from the cohort size.

**A REAL run requires ALL of:**

1. the typed phrase, exactly `run real samples`;
2. the machine overlay `config/bigmachine.yaml`;
3. the reuse mode shown, so a reader can see what was reused before what ran;
4. `PIPELINE_ALLOW_REAL_MODE` set **only in the child process environment**,
   never in this process and never in the dashboard's own environment.

Nothing else is executable. The command is an **allowlist**, not a template:
there is no shell interpolation of anything a client sent, and the run mode
cannot be chosen by a request body. `Runner` is the seam — QA (Q6) substitutes a
fake runner; this module never invokes bakta, panaroo, gubbins, iqtree, pyseer
or snakemake itself.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .security import TOKEN_ENV

#: The exact phrase. Compared with equality after stripping surrounding
#: whitespace and nothing else — no case folding, no prefix match, no
#: "did you mean".
REAL_PHRASE = "run real samples"

#: The environment variable the pipeline's own gate reads
#: (`config/loader.py:53`). Set on the child only.
ALLOW_REAL_MODE_ENV = "PIPELINE_ALLOW_REAL_MODE"

#: The machine overlay a REAL run must be launched under. The laptop overlay
#: caps the cohort at 20 samples and says why.
REQUIRED_OVERLAY = "bigmachine"

#: The only command the dry run may name. Allowlisted, never built from a
#: string.
SNAKEMAKE = "snakemake"
FULL_RULE = "full_run"

#: Tools whose per-isolate cost dominates a run. Flagged so the reader sees
#: where the time goes before starting, not after.
EXPENSIVE_TOOLS: Tuple[str, ...] = (
    "bakta",
    "panaroo",
    "gubbins",
    "iqtree",
    "pyseer",
)

#: When no runbook supplied them.
TIMINGS_NOT_RECORDED = "timings not recorded"

#: Where a runbook may be, probed in this order (assumption A1).
RUNBOOK_CANDIDATES: Tuple[str, ...] = (
    "06_for_900_isolates/RUNBOOK_900.md",
    "RUNBOOK_900.md",
    "04_run_info/RUNBOOK_900.md",
)

#: The pipeline CLI's own entry point, used so the launcher does not have to
#: guess how a run is invoked.
RUNNER_MODULE = "papipeline.cli"

#: The launcher is off unless the process was started with `--allow-launch`.
ALLOW_LAUNCH_ENV = "PA_DASH_ALLOW_LAUNCH"


class LauncherRefused(Exception):
    """A launch was refused. `message` is the sentence the client sees."""

    def __init__(self, message: str, detail: Optional[Mapping[str, Any]] = None) -> None:
        super().__init__(message)
        self.message = message
        self.detail = dict(detail or {})


@dataclass
class RunnerResult:
    """What a runner produced."""

    argv: Tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    dry_run: bool = True

    def as_dict(self) -> Dict[str, Any]:
        return {
            "argv": list(self.argv),
            "returncode": self.returncode,
            "stdout": self.stdout[-20000:],
            "stderr": self.stderr[-20000:],
            "dry_run": self.dry_run,
        }


class Runner:
    """The seam. QA (Q6) substitutes a fake; the default refuses to execute.

    Split out because the instruction for this phase is explicit: do not run
    bakta, panaroo, gubbins, iqtree, pyseer or snakemake while developing.
    """

    def __init__(self, *, executable: Optional[Sequence[str]] = None) -> None:
        self.executable = tuple(executable or (SNAKEMAKE,))

    def dry_run(self, repo_root: Path, extra: Sequence[str] = ()) -> RunnerResult:
        """The plan, without starting anything.

        `-n` is snakemake's dry-run flag and `--dry-run` its long form. Both are
        stated in the argv the caller can read, so a reader of the log can see
        that nothing was executed.
        """
        argv = (*self.executable, "--dry-run", FULL_RULE, *extra)
        return RunnerResult(argv=argv, returncode=0, stdout="", stderr="", dry_run=True)

    def execute(self, argv: Sequence[str], env: Mapping[str, str], cwd: Path) -> RunnerResult:
        """Refused in this build.

        Deliberate: Phase 1's launcher is a *preflight* surface. The command a
        real launch would run is chosen here, allowlisted, and shown — and the
        endpoint that would start it does not exist. `would_start` is therefore
        always false, which is what the openapi says.
        """
        raise LauncherRefused(
            "This build does not start a run. `POST /api/launcher/preflight` "
            "shows the plan and starts nothing; a launch is a human command "
            "run in the repository, not a button in a dashboard."
        )


def _fake_runner_factory(argv_builder: Any) -> Runner:
    return argv_builder


def capabilities(*, allow_launch: bool, loopback: bool) -> Dict[str, Any]:
    """What the launcher may legally do here.

    `can_launch` is false in every case in this build, and `reason` says why.
    The allowed/forbidden lists are the honest inventory, so a client can show
    the boundary rather than discovering it.
    """
    reason: str
    if not allow_launch:
        reason = (
            "The launcher is disabled: the server was started without "
            "--allow-launch. With the flag, the only action available is a "
            "dry run (`snakemake --dry-run full_run`), which shows the plan and "
            "starts nothing."
        )
    else:
        reason = (
            "This build starts nothing. With --allow-launch the launcher shows "
            "a dry-run plan and the gates a real run would have to pass; which "
            "command a real launch would run is not inferred here from the "
            "pipeline's CLI."
        )
    return {
        "can_launch": False,
        "reason": reason,
        "allow_launch_enabled": bool(allow_launch),
        "allowed_actions": [
            "POST /api/launcher/preflight: report the checks and the dry-run plan",
            f"show the allowlisted dry-run command {SNAKEMAKE} --dry-run {FULL_RULE}",
            "show per-isolate timings when a runbook supplied them",
        ],
        "forbidden_actions": [
            "starting a REAL run",
            "starting a STUB or TEST run",
            f"setting {ALLOW_REAL_MODE_ENV} in this process",
            f"accepting a token as a query parameter (it is the {TOKEN_ENV} header only)",
            "running bakta, panaroo, gubbins, iqtree or pyseer directly",
            "interpolating a client-supplied string into a shell command",
        ],
        "expensive_tools": list(EXPENSIVE_TOOLS),
        "real_phrase": REAL_PHRASE,
        "required_overlay": REQUIRED_OVERLAY,
    }


def locate_runbook(root: Optional[Path]) -> Tuple[Optional[Path], List[Dict[str, Any]]]:
    """Find `RUNBOOK_900.md`, and name every path probed.

    It does not exist in this repository yet (assumption A5), so the probe list
    is the deliverable as much as the path is.
    """
    probes: List[Dict[str, Any]] = []
    if root is None:
        probes.append({"n": 1, "path": None, "ok": False, "reason": "no results root is open"})
        return None, probes
    root = Path(root)
    for index, relative in enumerate(RUNBOOK_CANDIDATES, start=1):
        path = root / relative
        ok = path.exists()
        probes.append(
            {
                "n": index,
                "path": str(path),
                "ok": ok,
                "reason": "matched" if ok else "no file at this path",
            }
        )
        if ok:
            return path, probes
    return None, probes


def per_isolate_timings(root: Optional[Path]) -> Tuple[Optional[Dict[str, float]], str, List[Dict[str, Any]]]:
    """Per-isolate timings, or `timings not recorded`.

    Returns `(timings, timings_source, probes)`. `timings` is an OBJECT keyed on
    `sample_id` when a runbook supplied one and `None` otherwise — never an
    estimate derived from the cohort size, which would be an invented
    measurement (UI-D2, and the openapi's own note on this field).
    """
    path, probes = locate_runbook(root)
    if path is None:
        return None, "none", probes
    import re

    timings: Dict[str, float] = {}
    # A markdown table row: `| SAMPLE | 123.4 |` or `SAMPLE: 123.4`.
    row = re.compile(r"^\s*\|?\s*([A-Za-z0-9][A-Za-z0-9_.\-]*)\s*\|?\s*[:=]?\s*"
                     r"([0-9]+(?:\.[0-9]+)?)\s*(?:s|sec|seconds|min|m)?\s*\|?\s*$")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return None, "none", probes + [
            {"n": len(probes) + 1, "path": str(path), "ok": False,
             "reason": f"unreadable: {exc}"}
        ]
    for line in text.splitlines():
        match = row.match(line)
        if match:
            timings[match.group(1)] = float(match.group(2))
    if not timings:
        return None, "none", probes + [
            {"n": len(probes) + 1, "path": str(path), "ok": False,
             "reason": (
                 f"{path} was found but carries no per-isolate timing rows this "
                 f"reader recognises, so the timings are not reported rather "
                 f"than guessed"
             )}
        ]
    return timings, "runbook", probes


def status(*, log_path: Optional[Path], timings_source: str) -> Dict[str, Any]:
    """Is a run in progress, and what were the observed timings.

    `running` is derived from the event log: a `start` with no later terminal
    event for that stage. With no log at all it is `None` — "unknown", not
    "not running". Nothing is running because nothing was ever started.
    """
    from .events import aggregate, read_from
    from .papipeline_refs import STAGE_ORDER

    running: Optional[bool] = None
    last_event: Optional[Dict[str, Any]] = None
    if log_path is not None and Path(log_path).exists():
        replay = read_from(log_path, 0)
        state, _totals, _unmapped = aggregate(replay.events, stage_order=STAGE_ORDER)
        running = any(entry.running_samples for entry in state.values())
        if replay.events:
            last_event = replay.events[-1].as_dict()
    return {
        "running": running,
        "last_event": last_event,
        "per_isolate_timings": None,
        "timings_source": timings_source,
        "per_isolate_timings_label": (
            TIMINGS_NOT_RECORDED if timings_source == "none" else "recorded by a runbook"
        ),
    }


def _real_gates(
    *,
    phrase: Optional[str],
    overlay: Optional[str],
    repo_root: Optional[Path],
) -> Tuple[List[Dict[str, Any]], bool]:
    """The four gates a REAL run must pass, each reported."""
    checks: List[Dict[str, Any]] = []

    checks.append(
        {
            "name": "typed_phrase",
            "passed": phrase is not None and phrase == REAL_PHRASE,
            "detail": (
                f"the phrase {REAL_PHRASE!r} was typed exactly"
                if phrase is not None and phrase == REAL_PHRASE
                else (
                    f"a REAL run requires the phrase {REAL_PHRASE!r} typed "
                    f"exactly. Nothing shorter opens this gate: 'run real', "
                    f"'yes' or a capitalised variant are not the phrase, and a "
                    f"near-match is not consent to run on real genomes."
                )
            ),
        }
    )

    overlay_path = (
        Path(repo_root) / "config" / "machines" / f"{REQUIRED_OVERLAY}.yaml"
        if repo_root is not None
        else None
    )
    overlay_ok = bool(overlay) and overlay == REQUIRED_OVERLAY and (
        overlay_path is not None and overlay_path.exists()
    )
    checks.append(
        {
            "name": "machine_overlay",
            "passed": overlay_ok,
            "detail": (
                f"the {REQUIRED_OVERLAY!r} overlay is selected and present at "
                f"{overlay_path}"
                if overlay_ok
                else (
                    f"a REAL run requires the machine overlay "
                    f"config/machines/{REQUIRED_OVERLAY}.yaml. The laptop "
                    f"overlay caps a cohort at 20 samples and says why; it is "
                    f"not the overlay a real-genome run is launched under."
                )
            ),
        }
    )

    reuse_ok = reuse_mode_known(repo_root)
    checks.append(
        {
            "name": "reuse_mode_shown",
            "passed": reuse_ok,
            "detail": (
                "the annotation reuse mode is read from configuration and shown "
                "before anything runs, so a reader sees what will be reused "
                "rather than what will be executed"
                if reuse_ok
                else (
                    "the annotation reuse setting could not be read from "
                    "configuration, so it cannot be shown. A real run is not "
                    "started without it."
                )
            ),
        }
    )

    env_ok = os.environ.get(ALLOW_REAL_MODE_ENV) == "1"
    checks.append(
        {
            "name": "allow_real_mode_env",
            "passed": False,  # deliberately: never satisfied in this process
            "detail": (
                f"{ALLOW_REAL_MODE_ENV} is set on the CHILD process only, never "
                f"in the dashboard's own environment. This check reports what "
                f"the child would receive; it is not a gate this process can "
                f"satisfy for itself, and the dashboard does not set it."
            ),
        }
    )
    return checks, all(c["passed"] for c in checks)


def reuse_mode_known(repo_root: Optional[Path]) -> bool:
    """Whether the annotation reuse setting is readable and therefore showable."""
    if repo_root is None:
        return False
    path = Path(repo_root) / "config" / "science.yaml"
    if not path.exists():
        return False
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return False
    return "reuse_tool_output" in text


def reuse_mode(repo_root: Optional[Path]) -> Optional[str]:
    """The configured reuse mode, or None when it cannot be read."""
    if repo_root is None:
        return None
    path = Path(repo_root) / "config" / "science.yaml"
    try:
        import yaml

        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return None
    annotation = payload.get("annotation") or {}
    value = annotation.get("reuse_tool_output")
    return str(value) if value is not None else None


def preflight(
    *,
    allow_launch: bool,
    repo_root: Optional[Path],
    results_root: Optional[Path],
    log_path: Optional[Path],
    runner: Optional[Runner] = None,
    mode: str = "dry_run",
    phrase: Optional[str] = None,
    overlay: Optional[str] = None,
) -> Dict[str, Any]:
    """The checks and the plan. Starts nothing.

    `would_start` is always false — the openapi declares it `const: false` and
    this implementation honours that rather than shipping a schema the code
    contradicts.
    """
    timings, timings_source, runbook_probes = per_isolate_timings(results_root)
    checks: List[Dict[str, Any]] = [
        {
            "name": "launch_enabled",
            "passed": bool(allow_launch),
            "detail": (
                "--allow-launch was passed, so the dry-run plan below is shown"
                if allow_launch
                else (
                    "the launcher is disabled because --allow-launch was not "
                    "passed. Nothing is shown and nothing would run."
                )
            ),
        },
        {
            "name": "runbook_timings",
            "passed": timings_source == "runbook",
            "detail": (
                f"per-isolate timings were read from a runbook ({len(timings or {})} isolates)"
                if timings_source == "runbook"
                else (
                    f"{TIMINGS_NOT_RECORDED}. Probed: "
                    + "; ".join(
                        f"{p['path'] or '(unset)'} ({p['reason']})" for p in runbook_probes
                    )
                    + ". No per-isolate estimate is derived from the cohort "
                    "size; an invented timing would be a measurement nobody made."
                )
            ),
        },
        {
            "name": "expensive_tools_flagged",
            "passed": True,
            "detail": (
                "these tools dominate a run's wall clock and are flagged before "
                "anything starts: " + ", ".join(EXPENSIVE_TOOLS)
            ),
        },
    ]

    plan: Optional[Dict[str, Any]] = None
    if allow_launch:
        effective = runner or Runner()
        result = effective.dry_run(Path(repo_root) if repo_root else Path("."))
        plan = {
            "kind": "dry_run",
            "argv": list(result.argv),
            "note": (
                "`--dry-run` is snakemake's plan-only flag. This command was "
                "NOT executed by the dashboard; it is what a human would run "
                "to see the plan."
            ),
            "expensive_tools": list(EXPENSIVE_TOOLS),
            "per_isolate_timings": timings,
            "per_isolate_timings_label": (
                TIMINGS_NOT_RECORDED if timings_source == "none" else "recorded by a runbook"
            ),
        }

    real_checks: List[Dict[str, Any]] = []
    real_ok = False
    if str(mode).strip() == "real":
        real_checks, real_ok = _real_gates(
            phrase=phrase, overlay=overlay, repo_root=repo_root
        )

    payload: Dict[str, Any] = {
        "would_start": False,
        "checks": checks,
        "plan": plan,
        "runbook_probes": runbook_probes,
        "reuse_mode": reuse_mode(repo_root),
    }
    if str(mode).strip() == "real":
        payload["real_gates"] = real_checks
        payload["real_gates_all_passed"] = real_ok
        payload["reason"] = (
            (
                f"a REAL run would need the phrase {REAL_PHRASE!r} typed "
                f"exactly, the {REQUIRED_OVERLAY!r} overlay, and "
                f"{ALLOW_REAL_MODE_ENV} set on the child process only. Even "
                f"with all of them, this build has no endpoint that starts a "
                f"run: a launch is a human command run in the repository."
            )
        )
    return payload


def log_launch(state_dir: Any, record: Mapping[str, Any]) -> Path:
    """Append one launch attempt to the state directory's log.

    The state directory is the only write location (UI-D1). Every launch is
    recorded — including the ones refused — because a refused real-genome launch
    is exactly the event a reader of the log needs to see.
    """
    import json

    path = Path(state_dir.logs_dir / "launcher.log")
    entry = {
        "t": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z"),
        **dict(record),
    }
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry, sort_keys=True, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    return path


__all__ = [
    "ALLOW_REAL_MODE_ENV",
    "ALLOW_LAUNCH_ENV",
    "EXPENSIVE_TOOLS",
    "REAL_PHRASE",
    "REQUIRED_OVERLAY",
    "RUNBOOK_CANDIDATES",
    "TIMINGS_NOT_RECORDED",
    "LauncherRefused",
    "Runner",
    "RunnerResult",
    "capabilities",
    "log_launch",
    "locate_runbook",
    "per_isolate_timings",
    "preflight",
    "reuse_mode",
    "status",
]