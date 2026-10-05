"""K3 — the event stream, and the polling fallback.

**The dashboard reads `status/events.jsonl`, and it is the only mechanism it
uses** (DESIGN §6.1). The observatory's in-process `EventBus` is not read:
its store upserts on `(run_key, stage, subject)` while `run.py:370` passes
`subject=""`, so a 900-isolate run holds 16 rows rather than the ~14,400 cells
the grid needs.

Four rules from the contract (`scripts/emit.py`):

- exactly four fields — `t`, `event`, `stage`, `sample`;
- exactly three `event` values — `start`, `done`, `fail`;
- `sample` is a sample id or the literal `all` for an aggregate rule;
- a line that is not valid JSON, or that is missing a field, or that carries a
  fifth field or a fourth value, is **skipped with a warning naming its line
  number**. It is never coerced into the shape.

The file is appended with `flush()` + `fsync()` on every record, so a run killed
mid-write leaves a replayable log, which is what `emit.py:85-93` was written to
guarantee.

Read-only throughout: opened `'r'`, tailed by byte offset, never truncated,
never written. The SSE generator bounds its queue and closes its handle on
disconnect, so a client that goes away leaks nothing.
"""

from __future__ import annotations

import asyncio
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, AsyncIterator, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

from .source import EVENT_LOG_RELATIVE

#: The three accepted transitions. A fourth is a contract change, not a
#: convenience, so it is skipped rather than tolerated.
EVENTS: Tuple[str, ...] = ("start", "done", "fail")

#: The four fields, and no others.
FIELDS: Tuple[str, ...] = ("t", "event", "stage", "sample")

#: The sentinel an aggregate rule uses for its cohort-wide transitions.
AGGREGATE_SAMPLE = "all"

#: Snakemake rule name -> code stage name (DESIGN §6.2).
#:
#: The log carries RULE names, which are not always the code stage names
#: (`rule validate` vs `validation`). An unmapped name is reported in the
#: stream's own diagnostics rather than dropped.
RULE_TO_STAGE: Mapping[str, str] = {
    # rule names that match the code stage names exactly
    "validation": "validation",
    "annotation": "annotation",
    "mlst": "mlst",
    "amr": "amr",
    "variants": "variants",
    "cohort_variants": "cohort_variants",
    "recombination": "recombination",
    "similarity": "similarity",
    "virulence": "virulence",
    "pangenome": "pangenome",
    "phylogeny": "phylogeny",
    "phenotype": "phenotype",
    "gwas": "gwas",
    "convergence": "convergence",
    "cooccurrence": "cooccurrence",
    "reporting": "reporting",
    # spec.md D1 names, which differ for four stages (§1)
    "validate": "validation",
    "annotate": "annotation",
    "combination": "cooccurrence",
    "report": "reporting",
    # aggregate rules that are not stages. Mapped to None below.
}

#: Rule names that exist in `workflow/Snakefile` and are not stages. Recorded
#: rather than reported as unmapped, because `all` and `clean_results` are
#: ordinary and their absence from the stage axis is not a defect.
NON_STAGE_RULES: frozenset = frozenset({
    "all", "provenance", "synthetic_fixtures", "full_run",
    "clean_results", "list_stages",
})


def map_stage(stage: str) -> Optional[str]:
    """The code stage a log line names, or None when it names no stage."""
    return RULE_TO_STAGE.get(str(stage).strip())


@dataclass
class Event:
    """One line, validated."""

    t: str
    event: str
    stage: str
    sample: str
    line: int
    mapped_stage: Optional[str] = None
    offset: int = 0

    def as_dict(self) -> Dict[str, Any]:
        return {
            "t": self.t,
            "event": self.event,
            "stage": self.stage,
            "sample": self.sample,
            "line": self.line,
            "mapped_stage": self.mapped_stage,
        }


@dataclass
class Replay:
    """The result of reading the log from an offset."""

    cursor: int = 0
    events: List[Event] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    truncated: bool = False

    @property
    def n_skipped(self) -> int:
        return len(self.warnings)


def _parse_line(raw: str, line_number: int, offset: int) -> Tuple[Optional[Event], str]:
    """Validate one line, or return the warning that says why it was skipped.

    The strictness is the contract: a fifth field or a fourth `event` value is
    not a tolerance to absorb, and neither is a missing field.
    """
    try:
        payload = json.loads(raw)
    except ValueError as exc:
        return None, (
            f"line {line_number}: not valid JSON ({exc}); skipped rather than "
            f"coerced"
        )
    if not isinstance(payload, dict):
        return None, (
            f"line {line_number}: not a JSON object; the contract is four "
            f"fields on one line"
        )
    keys = set(payload)
    if keys != set(FIELDS):
        missing = sorted(set(FIELDS) - keys)
        extra = sorted(keys - set(FIELDS))
        return None, (
            f"line {line_number}: field set is {sorted(keys)}, which is not the "
            f"contract's {list(FIELDS)}"
            + (f"; missing {missing}" if missing else "")
            + (f"; carries an extra field {extra}" if extra else "")
            + ". A fifth field is a contract change, not a tolerance."
        )
    event = payload["event"]
    if event not in EVENTS:
        return None, (
            f"line {line_number}: event {event!r} is not one of "
            f"{', '.join(EVENTS)}; a fourth value is a contract change, not a "
            f"tolerance"
        )
    for name in ("t", "stage", "sample"):
        if not isinstance(payload[name], str) or not payload[name].strip():
            return None, (
                f"line {line_number}: {name} is absent or empty; the contract "
                f"requires it"
            )
    return (
        Event(
            t=payload["t"],
            event=event,
            stage=payload["stage"],
            sample=payload["sample"],
            line=line_number,
            mapped_stage=map_stage(payload["stage"]),
            offset=offset,
        ),
        "",
    )


def read_from(path: Optional[Path], cursor: int = 0, *, limit: int = 5000) -> Replay:
    """Read the log from ``cursor`` to its current end.

    Opened read-only. A log that shrank — a run deleted it, or it was rotated
    — resets the cursor to 0 rather than seeking into the middle of a file that
    no longer holds those bytes.
    """
    result = Replay(cursor=cursor)
    if path is None:
        result.warnings.append(
            "no event log was found; probes are listed in the response. An "
            "absent log is not an empty log: an empty log reads as 'nothing "
            "has happened'."
        )
        return result
    path = Path(path)
    if not path.exists():
        result.warnings.append(f"no event log at {path}")
        return result
    try:
        size = path.stat().st_size
    except OSError as exc:
        result.warnings.append(f"{path} cannot be stat'd: {exc}")
        return result
    start = cursor if 0 <= cursor <= size else 0
    if start != cursor:
        result.warnings.append(
            f"cursor {cursor} is past the end of {path} ({size} bytes); the log "
            f"was replaced or truncated, so the replay restarts from 0"
        )
    try:
        with open(path, "r", encoding="utf-8", newline="", errors="replace") as handle:
            handle.seek(start)
            offset = start
            line_number = 0
            # Count absolute line numbers from the byte offset so a warning
            # names the line in the FILE, not the line in this read.
            line_number = _lines_before(path, start)
            for raw in handle:
                byte_length = len(raw.encode("utf-8", errors="replace"))
                line_number += 1
                start_of_line = offset
                offset += byte_length
                if not raw.endswith("\n"):
                    # A partial final line: a run killed mid-write. Left for the
                    # next read rather than parsed as a truncated record.
                    result.truncated = True
                    offset = start_of_line
                    break
                stripped = raw.strip()
                if not stripped:
                    continue
                event, warning = _parse_line(stripped, line_number, start_of_line)
                if warning:
                    result.warnings.append(warning)
                    continue
                assert event is not None
                result.events.append(event)
                if len(result.events) >= limit:
                    result.truncated = True
                    break
            result.cursor = offset
    except OSError as exc:
        result.warnings.append(f"{path} could not be read: {exc}")
    return result


def _lines_before(path: Path, offset: int) -> int:
    """How many complete lines precede ``offset``. One cheap pass."""
    if offset <= 0:
        return 0
    count = 0
    try:
        with open(path, "rb") as handle:
            remaining = offset
            while remaining > 0:
                chunk = handle.read(min(65536, remaining))
                if not chunk:
                    break
                count += chunk.count(b"\n")
                remaining -= len(chunk)
    except OSError:
        return 0
    return count


def find_log(root: Optional[Path], repo_root: Optional[Path] = None) -> Tuple[Optional[Path], List[Dict[str, Any]]]:
    """Where the event log is, and every path probed.

    `emit.py`'s default is the repository-relative `status/events.jsonl`; a
    results-root-relative one is also probed because a bundle may have
    collected it.
    """
    probes: List[Dict[str, Any]] = []
    candidates: List[Optional[Path]] = [
        (Path(root) / EVENT_LOG_RELATIVE) if root is not None else None,
        (Path(repo_root) / EVENT_LOG_RELATIVE) if repo_root is not None else None,
    ]
    for candidate in candidates:
        if candidate is None:
            probes.append(
                {"n": len(probes) + 1, "path": None, "ok": False,
                 "reason": "not configured"}
            )
            continue
        ok = candidate.exists()
        probes.append(
            {
                "n": len(probes) + 1,
                "path": str(candidate),
                "ok": ok,
                "reason": "matched" if ok else "no file at this path",
            }
        )
        if ok:
            return candidate, probes
    return None, probes


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------


@dataclass
class StageEventState:
    """What the log says about one stage."""

    started: int = 0
    done: int = 0
    failed: int = 0
    last_event: Optional[Event] = None
    running_samples: Dict[str, str] = field(default_factory=dict)

    @property
    def badge_hint(self) -> Optional[str]:
        """The badge the log alone implies, or None when it implies nothing.

        Advisory: the authoritative badge comes from
        `badges.classify`, which also reads the manifest and the table. This is
        what `/api/events/state` carries for the metric cards.
        """
        if self.failed:
            return "failed"
        if self.running_samples:
            return "running"
        if self.done:
            return "completed"
        if self.started:
            return "running"
        return None


def aggregate(
    events: Sequence[Event],
    *,
    stage_order: Sequence[str],
    prior: Optional[Mapping[str, StageEventState]] = None,
) -> Tuple[Dict[str, StageEventState], Dict[str, int], List[str]]:
    """Fold events into per-stage state, plus totals and unmapped warnings."""
    state: Dict[str, StageEventState] = {
        name: (prior[name] if prior and name in prior else StageEventState())
        for name in stage_order
    }
    unmapped: List[str] = []
    totals = {"start": 0, "done": 0, "fail": 0}
    for event in events:
        totals[event.event] = totals.get(event.event, 0) + 1
        stage = event.mapped_stage
        if stage is None:
            if event.stage not in NON_STAGE_RULES:
                unmapped.append(
                    f"line {event.line}: stage {event.stage!r} maps to no code "
                    f"stage name. The log carries Snakemake rule names, so a new "
                    f"rule needs an entry in the rule->stage table."
                )
            continue
        entry = state.setdefault(stage, StageEventState())
        entry.last_event = event
        if event.event == "start":
            entry.started += 1
            entry.running_samples[event.sample] = event.t
        elif event.event == "done":
            entry.done += 1
            entry.running_samples.pop(event.sample, None)
        else:
            entry.failed += 1
            entry.running_samples.pop(event.sample, None)
    return state, totals, unmapped


# ---------------------------------------------------------------------------
# SSE
# ---------------------------------------------------------------------------

#: Frame order: retry once, then a full snapshot, then deltas on a 100 ms
#: timer. A comment heartbeat every 15 s keeps an idle run from looking dead.
RETRY_MS = 2000
BATCH_INTERVAL_S = 0.1
HEARTBEAT_S = 15.0

#: The queue is bounded. An unbounded one turns a fast producer and a slow
#: client into a memory leak, and a run that emits 14,400 events would do it
#: in seconds.
QUEUE_MAX = 4096


def _frame(*, event: Optional[str], data: str, event_id: Optional[int] = None) -> bytes:
    """One SSE frame, encoded."""
    parts: List[str] = []
    if event:
        parts.append(f"event: {event}")
    if event_id is not None:
        parts.append(f"id: {event_id}")
    parts.append(f"data: {data}")
    return ("\n".join(parts) + "\n\n").encode("utf-8")


def _comment(text: str) -> bytes:
    """An SSE comment: a heartbeat the client ignores and a proxy keeps alive."""
    return f": {text}\n\n".encode("utf-8")


async def sse_frames(
    log_path: Optional[Path],
    *,
    initial: Mapping[str, Any],
    last_event_id: Optional[int] = None,
    poll_interval: float = BATCH_INTERVAL_S,
) -> AsyncIterator[bytes]:
    """Replay the log, then tail it.

    Frame order, per §6.3: `retry: 2000` once at connect; one `snapshot`
    carrying the full replayed state **before** any `delta`, so a client that
    connects mid-run is never blank; `delta` per new line batched on a 100 ms
    timer; `id: <byte offset>` on every event so `Last-Event-ID` resumes
    exactly; a comment heartbeat every 15 s.

    Read-only: `read_from` opens the file `'r'`, seeks, and closes. Nothing
    here truncates or writes.
    """
    yield f"retry: {RETRY_MS}\n\n".encode("utf-8")

    replay = read_from(log_path, cursor=last_event_id or 0)
    snapshot = dict(initial)
    snapshot["cursor"] = replay.cursor
    snapshot["n_events"] = len(replay.events)
    snapshot["n_lines_skipped"] = replay.n_skipped
    snapshot["skipped_warnings"] = replay.warnings[:200]
    yield _frame(event="snapshot", data=json.dumps(snapshot, default=str), event_id=replay.cursor)

    cursor = replay.cursor
    queue: asyncio.Queue = asyncio.Queue(maxsize=QUEUE_MAX)
    dropped = 0

    async def pump() -> None:
        """Read new lines into the bounded queue, batching on the timer."""
        nonlocal cursor, dropped
        while True:
            await asyncio.sleep(poll_interval)
            try:
                step = await asyncio.to_thread(read_from, log_path, cursor)
            except Exception as exc:  # a read failure must not kill the stream
                await _offer(queue, {"warning": f"event log read failed: {exc}"})
                continue
            for warning in step.warnings:
                await _offer(queue, {"warning": warning})
            if step.events:
                payload = [e.as_dict() for e in step.events]
                try:
                    queue.put_nowait({"delta": payload})
                except asyncio.QueueFull:
                    # Bounded: the oldest batch goes, and the drop is announced
                    # rather than hidden. A client that reconnects with
                    # Last-Event-ID resumes from a byte offset, so nothing is
                    # permanently lost.
                    try:
                        queue.get_nowait()
                        dropped += 1
                        queue.put_nowait({"dropped": dropped})
                        queue.put_nowait({"delta": payload})
                    except (asyncio.QueueEmpty, asyncio.QueueFull):
                        pass
            if step.cursor != cursor:
                cursor = step.cursor
            if step.truncated:
                await _offer(queue, {"truncated": True, "cursor": cursor})

    async def _offer(q: asyncio.Queue, item: Any) -> None:
        try:
            q.put_nowait(item)
        except asyncio.QueueFull:
            pass

    pump_task = asyncio.create_task(pump())
    heartbeat = 0.0
    try:
        while True:
            try:
                item = await asyncio.wait_for(queue.get(), timeout=HEARTBEAT_S)
            except asyncio.TimeoutError:
                heartbeat += HEARTBEAT_S
                yield _comment("heartbeat")
                continue
            if "delta" in item:
                yield _frame(
                    event="delta",
                    data=json.dumps(item["delta"], default=str),
                    event_id=cursor,
                )
            elif "warning" in item:
                yield _comment(f"warning: {item['warning']}")
            elif "dropped" in item:
                yield _comment(
                    f"dropped {item['dropped']} batch(es) to keep the queue "
                    f"bounded; reconnect with Last-Event-ID to resume exactly"
                )
            elif item.get("truncated"):
                yield _frame(
                    event="delta",
                    data=json.dumps(
                        {"truncated": True, "cursor": item.get("cursor", cursor)},
                        default=str,
                    ),
                    event_id=cursor,
                )
    except asyncio.CancelledError:
        # The client disconnected. Cancelling the pump closes the file handle
        # it holds; nothing is left open behind a dropped connection.
        raise
    finally:
        pump_task.cancel()
        try:
            await pump_task
        except (asyncio.CancelledError, Exception):
            pass


__all__ = [
    "AGGREGATE_SAMPLE",
    "EVENTS",
    "FIELDS",
    "NON_STAGE_RULES",
    "QUEUE_MAX",
    "RULE_TO_STAGE",
    "Event",
    "Replay",
    "StageEventState",
    "aggregate",
    "find_log",
    "map_stage",
    "read_from",
    "sse_frames",
]