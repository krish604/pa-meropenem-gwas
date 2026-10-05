"""The status event log is the dashboard's only input.

Ticket 09. Every rule announces that it started, succeeded or failed by
appending one JSON object per line to ``status/events.jsonl``. The format is
deliberately tiny: a fourth field, or a fourth value of ``event``, would be a
contract the dashboard has to know about, and there is no second consumer to
justify the complexity.

The properties that matter:

* exactly the four fields, and exactly three values of ``event``;
* aggregate rules record the literal sample ``all``, so event volume stays
  proportionate to the work rather than to the cohort size;
* append-only, so a re-run does not erase the history of the run before it;
* a malformed or truncated line is survivable, because a run killed mid-write
  still has to leave a dashboard that can be rebuilt.

Written before the implementation.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
EMITTER = PIPELINE_ROOT / "scripts" / "emit.py"

FIELDS = {"t", "event", "stage", "sample"}
EVENTS = {"start", "done", "fail"}


def _emit(log: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(EMITTER), "--log", str(log), *args],
        capture_output=True, text=True, timeout=60,
    )


def _lines(log: Path) -> list[dict]:
    return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# the shape of one event
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("event", sorted(EVENTS))
def test_the_three_events_are_accepted(tmp_path: Path, event: str):
    log = tmp_path / "events.jsonl"
    assert _emit(log, "--event", event, "--stage", "annotate",
                 "--sample", "GCA_000000001.1").returncode == 0
    row = _lines(log)[0]
    assert set(row) == FIELDS, f"exactly four fields, got {sorted(row)}"
    assert row["event"] == event
    assert row["stage"] == "annotate"
    assert row["sample"] == "GCA_000000001.1"


def test_the_timestamp_is_iso8601_and_utc(tmp_path: Path):
    import datetime as dt

    log = tmp_path / "events.jsonl"
    _emit(log, "--event", "start", "--stage", "validate", "--sample", "all")
    stamp = _lines(log)[0]["t"]
    parsed = dt.datetime.fromisoformat(stamp)
    assert parsed.tzinfo is not None, "a naive timestamp is ambiguous across machines"


def test_an_unknown_event_value_is_refused(tmp_path: Path):
    """Four values of `event` would be a contract change, not a convenience."""
    log = tmp_path / "events.jsonl"
    result = _emit(log, "--event", "retry", "--stage", "annotate", "--sample", "x")
    assert result.returncode != 0
    assert "retry" in (result.stderr + result.stdout)
    assert not log.exists() or _lines(log) == []


def test_a_missing_stage_is_refused(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    result = _emit(log, "--event", "start", "--sample", "x")
    assert result.returncode != 0
    assert "stage" in (result.stderr + result.stdout).lower()


def test_a_missing_sample_is_refused(tmp_path: Path):
    """`all` is a value, not a licence to omit the field."""
    log = tmp_path / "events.jsonl"
    result = _emit(log, "--event", "start", "--stage", "annotate")
    assert result.returncode != 0
    assert "sample" in (result.stderr + result.stdout).lower()


# ---------------------------------------------------------------------------
# aggregate rules
# ---------------------------------------------------------------------------


def test_an_aggregate_rule_records_the_literal_all(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    _emit(log, "--event", "start", "--stage", "pangenome", "--sample", "all")
    assert _lines(log)[0]["sample"] == "all"


# ---------------------------------------------------------------------------
# append-only
# ---------------------------------------------------------------------------


def test_events_accumulate_rather_than_replace(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    _emit(log, "--event", "start", "--stage", "annotate", "--sample", "GCA_000000001.1")
    _emit(log, "--event", "done", "--stage", "annotate", "--sample", "GCA_000000001.1")
    rows = _lines(log)
    assert [r["event"] for r in rows] == ["start", "done"]


def test_a_lifecycle_is_start_done_in_order(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    for event in ("start", "done"):
        _emit(log, "--event", event, "--stage", "annotate", "--sample", "GCA_000000001.1")
    assert [r["event"] for r in _lines(log)] == ["start", "done"]


def test_the_log_directory_is_created_on_demand(tmp_path: Path):
    log = tmp_path / "nested" / "status" / "events.jsonl"
    assert _emit(log, "--event", "start", "--stage", "validate",
                 "--sample", "all").returncode == 0
    assert log.is_file()


# ---------------------------------------------------------------------------
# survivability: a run killed mid-write
# ---------------------------------------------------------------------------


def test_a_truncated_final_line_does_not_stop_the_next_write(tmp_path: Path):
    """A crashed run leaves a partial line; the log must stay usable."""
    log = tmp_path / "events.jsonl"
    _emit(log, "--event", "start", "--stage", "annotate", "--sample", "GCA_000000001.1")
    with log.open("a", encoding="utf-8") as handle:
        handle.write('{"t": "2026-01-01T00:00:00", "event": "do')  # no newline
    assert _emit(log, "--event", "fail", "--stage", "annotate",
                 "--sample", "GCA_000000001.1").returncode == 0
    raw = log.read_text(encoding="utf-8")
    assert raw.count("\n") >= 2
    # the first and the last complete lines are intact
    complete = [json.loads(l) for l in raw.splitlines() if l.strip().startswith("{")
                and l.strip().endswith("}")]
    assert [r["event"] for r in complete] == ["start", "fail"]


def test_one_line_is_always_one_json_object(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    for i in range(5):
        _emit(log, "--event", "done", "--stage", "variants", "--sample", f"GCA_{i:09d}.1")
    for line in log.read_text(encoding="utf-8").splitlines():
        json.loads(line)  # raises if a line is not exactly one object


def test_a_stage_or_sample_containing_a_quote_cannot_break_the_line(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    _emit(log, "--event", "fail", "--stage", 'we"ird', "--sample", "x")
    row = _lines(log)[0]
    assert row["stage"] == 'we"ird'
    assert row["sample"] == "x"


# ---------------------------------------------------------------------------
# the file must not be committed
# ---------------------------------------------------------------------------


def test_the_status_directory_is_git_ignored():
    import subprocess as sp

    out = sp.run(["git", "check-ignore", "-q", "status/events.jsonl"],
                 cwd=PIPELINE_ROOT, capture_output=True, text=True)
    assert out.returncode == 0, "status/events.jsonl must never be committed"
