#!/usr/bin/env python3
"""Run one stage, announcing it to the status log on the way in and out.

Every rule in the workflow goes through here, so the event stream is a
by-product of running the pipeline rather than a parallel thing that can drift
out of step with it. Three events per rule invocation, or two plus a failure:

    start   before the command runs
    done    after it exits 0
    fail    after it exits non-zero

The stage's own exit code is propagated, and a `fail` event is emitted for any
non-zero exit, not just the ones Snakemake would consider errors. A stage that
dies is exactly the thing a dashboard exists to show.

If the event log cannot be written, the stage still runs and still returns its
own exit code. The log is observability, not the science, and a full disk must
not be allowed to destroy a run. That is reported on stderr, loudly.

The sample argument is ``all`` here because these rules are cohort-level: one
invocation covers every isolate. Per-isolate events arrive when a stage fans out
per sample, and this script already takes the flag for exactly that.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Sequence

EMITTER = Path(__file__).resolve().parent / "emit.py"


def emit(log: Path, event: str, stage: str, sample: str) -> None:
    """Append one event, never raising."""
    result = subprocess.run(
        [sys.executable, str(EMITTER), "--log", str(log),
         "--event", event, "--stage", stage, "--sample", sample],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        sys.stderr.write(
            f"run_rule: could not record the {event!r} event for stage "
            f"{stage!r}: {result.stderr.strip() or result.returncode}\n"
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, help="Stage name to record.")
    parser.add_argument(
        "--sample", default="all",
        help="Sample identifier, or 'all' for a cohort-level rule.",
    )
    parser.add_argument(
        "--status-log", type=Path, default=Path("status/events.jsonl"),
        help="The JSON Lines file the dashboard reads.",
    )
    parser.add_argument(
        "command", nargs=argparse.REMAINDER,
        help="The command to run, after a literal '--'.",
    )
    args = parser.parse_args(argv)

    command = [c for c in args.command if c != "--"]
    if not command:
        print("run_rule: no command given", file=sys.stderr)
        return 2

    emit(args.status_log, "start", args.stage, args.sample)
    result = subprocess.run(command, check=False)
    code = result.returncode

    emit(
        args.status_log,
        "done" if code == 0 else "fail",
        args.stage,
        args.sample,
    )
    return code


if __name__ == "__main__":
    raise SystemExit(main())
