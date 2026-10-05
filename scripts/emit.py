#!/usr/bin/env python3
"""Append one status event per rule transition to ``status/events.jsonl``.

NOT currently read by anything. An earlier version of this docstring called
the file "the dashboard's only input ... it replays this file on start and
then tails it". That is not true of the code as it stands: the observatory
serves /api/events from an in-process EventBus fed by the runner, and nothing
in papipeline/ opens this file. There are two observability mechanisms and
they do not meet, and this one has no consumer yet. The format is kept
deliberately minimal so that bridging them, when it happens, is a small change:

    {"t": "<iso8601>", "event": "start|done|fail", "stage": "...", "sample": "..."}

``event`` has exactly three values and the record has exactly four fields. A
fifth field, or a fourth value, is a contract the dashboard has to understand,
and there is no second consumer to justify that complexity.

Aggregate rules record the literal sample ``all``, which keeps the number of
lines proportional to the work done rather than to the cohort size.

Written to be called once per transition from a rule's ``shell:`` directive, so
it must be cheap, must not import the pipeline, and must fail loudly but
without taking the run down with it.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
from pathlib import Path
from typing import Sequence

#: The only accepted transitions. A fourth value is a contract change, not a
#: convenience, so it is refused rather than tolerated.
EVENTS = ("start", "done", "fail")

#: The sentinel an aggregate rule uses for its cohort-wide transitions.
AGGREGATE_SAMPLE = "all"

FIELDS = ("t", "event", "stage", "sample")


def _now() -> str:
    """An aware, UTC, second-resolution timestamp.

    A naive timestamp is ambiguous between machines, and this log is read by a
    browser on a different machine from the one that wrote it.
    """
    return (
        dt.datetime.now(dt.timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def build_record(event: str, stage: str, sample: str) -> dict:
    """Validate the three caller-supplied fields and stamp the record.

    Raises:
        ValueError: A field is missing, or ``event`` is not one of the three.
    """
    if event not in EVENTS:
        raise ValueError(
            f"event must be one of {', '.join(EVENTS)}; got {event!r}"
        )
    if not stage or not stage.strip():
        raise ValueError("stage is required and may not be empty")
    if not sample or not sample.strip():
        raise ValueError(
            "sample is required and may not be empty; use "
            f"{AGGREGATE_SAMPLE!r} for an aggregate rule"
        )
    return {
        "t": _now(),
        "event": event,
        "stage": stage.strip(),
        "sample": sample.strip(),
    }


def append_event(log_path: Path, record: dict) -> None:
    """Append one JSON object as one line.

    The file is opened in append mode and never truncated, because the log is
    the history of a run: a re-run adds to it rather than erasing what the
    previous run did. A trailing newline is written even when the file did not
    end with one, so a run killed mid-write cannot leave the next record glued
    onto a partial one.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(record, separators=(",", ":"), sort_keys=False) + "\n"

    # A run killed mid-write can leave a partial final line. Probing in a
    # separate handle keeps the append handle write-only, which is what
    # "append" means.
    needs_newline = False
    if log_path.exists() and log_path.stat().st_size > 0:
        with log_path.open("rb") as probe:
            probe.seek(-1, os.SEEK_END)
            needs_newline = probe.read(1) != b"\n"

    with log_path.open("a", encoding="utf-8") as handle:
        if needs_newline:
            handle.write("\n")
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log",
        type=Path,
        default=Path("status/events.jsonl"),
        help="JSON Lines file to append to.",
    )
    parser.add_argument("--event", required=True, help=f"One of: {', '.join(EVENTS)}")
    parser.add_argument("--stage", required=True, help="Stage name.")
    parser.add_argument(
        "--sample",
        required=True,
        help=f"Sample identifier, or {AGGREGATE_SAMPLE!r} for an aggregate rule.",
    )
    args = parser.parse_args(argv)

    try:
        record = build_record(args.event, args.stage, args.sample)
    except ValueError as exc:
        print(f"emit.py: {exc}", file=sys.stderr)
        return 2

    try:
        append_event(args.log, record)
    except OSError as exc:
        # The log is observability, not the science. A run must not be lost
        # because the status directory is unwritable, but it must be said.
        print(f"emit.py: could not append to {args.log}: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
