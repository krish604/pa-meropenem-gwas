"""`scripts/run_rule.py` is the event emitter's only production caller.

Ticket 09, promoted from `built` to `wired`. The emitter was finished and
well-tested while nothing invoked it, so the dashboard's input had no producer.
Wrapping the command is the wiring: one place that announces a stage on the way
in and out, so the event stream cannot drift out of step with what actually ran.

The properties that matter, beyond "it appends three lines":

* the stage's own exit code is propagated unchanged, so wrapping a rule cannot
  turn a failure into a success;
* a failing command still produces a `fail` event, because a stage that dies is
  precisely what a dashboard exists to show;
* an unwritable status log does not stop the run. The log is observability, not
  the science.

Written before the implementation.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
RUN_RULE = PIPELINE_ROOT / "scripts" / "run_rule.py"


def _run(log: Path, *command: str, stage: str = "annotate", sample: str = "all"):
    return subprocess.run(
        [sys.executable, str(RUN_RULE), "--stage", stage, "--sample", sample,
         "--status-log", str(log), "--", *command],
        capture_output=True, text=True, timeout=300,
    )


def _events(log: Path):
    return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# the happy path
# ---------------------------------------------------------------------------


def test_a_succeeding_command_emits_start_then_done(tmp_path: Path):
    log = tmp_path / "status" / "events.jsonl"
    result = _run(log, sys.executable, "-c", "pass")
    assert result.returncode == 0
    events = _events(log)
    assert [e["event"] for e in events] == ["start", "done"]
    assert all(e["stage"] == "annotate" for e in events)
    assert all(e["sample"] == "all" for e in events)


def test_the_status_directory_is_created(tmp_path: Path):
    log = tmp_path / "deeply" / "nested" / "status" / "events.jsonl"
    assert _run(log, sys.executable, "-c", "pass").returncode == 0
    assert log.is_file()


# ---------------------------------------------------------------------------
# failure
# ---------------------------------------------------------------------------


def test_a_failing_command_emits_start_then_fail_and_propagates_the_code(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    result = _run(log, sys.executable, "-c", "raise SystemExit(7)")
    assert result.returncode == 7, "the stage's own exit code must survive wrapping"
    assert [e["event"] for e in _events(log)] == ["start", "fail"]


def test_a_signalled_command_is_still_reported_as_a_failure(tmp_path: Path):
    """A crash is not an orderly exit, and the dashboard must still see it."""
    log = tmp_path / "events.jsonl"
    result = _run(log, sys.executable, "-c", "import os, signal; os.kill(os.getpid(), signal.SIGKILL)")
    events = _events(log)
    assert events[0]["event"] == "start"
    assert events[-1]["event"] == "fail", f"expected a fail event, got {events}"
    assert result.returncode != 0


# ---------------------------------------------------------------------------
# the log is observability, not the science
# ---------------------------------------------------------------------------


def test_an_unwritable_status_log_does_not_stop_the_run(tmp_path: Path):
    """A full disk must not be able to destroy a run."""
    unwritable = tmp_path / "afile"
    unwritable.write_text("not a directory", encoding="utf-8")
    log = unwritable / "events.jsonl"        # cannot be created under a file
    result = _run(log, sys.executable, "-c", "pass")
    assert result.returncode == 0, "the command's own result must win"
    assert "could not record" in result.stderr


def test_an_unwritable_log_does_not_mask_a_command_failure(tmp_path: Path):
    unwritable = tmp_path / "afile"
    unwritable.write_text("not a directory", encoding="utf-8")
    log = unwritable / "events.jsonl"
    result = _run(log, sys.executable, "-c", "raise SystemExit(3)")
    assert result.returncode == 3


# ---------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------


def test_no_command_is_refused(tmp_path: Path):
    log = tmp_path / "events.jsonl"
    result = subprocess.run(
        [sys.executable, str(RUN_RULE), "--stage", "x", "--status-log", str(log)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode != 0
    assert not log.exists() or _events(log) == []


def test_a_per_sample_stage_records_that_sample(tmp_path: Path):
    """The cohort default is 'all'; a per-isolate rule must be able to say so."""
    log = tmp_path / "events.jsonl"
    _run(log, sys.executable, "-c", "pass", sample="GCA_000000001.1")
    assert all(e["sample"] == "GCA_000000001.1" for e in _events(log))


def test_two_rules_do_not_interleave_into_one_record(tmp_path: Path):
    """Concurrent workers share the log, so a record must be written whole."""
    log = tmp_path / "events.jsonl"
    for stage in ("annotate", "mlst", "amr", "phenotype"):
        _run(log, sys.executable, "-c", "pass", stage=stage)
    rows = _events(log)
    assert len(rows) == 8
    assert sorted({r["stage"] for r in rows}) == ["amr", "annotate", "mlst", "phenotype"]


# ---------------------------------------------------------------------------
# the workflow actually uses it
# ---------------------------------------------------------------------------


def test_the_workflow_actually_routes_its_rules_through_this_wrapper():
    """The production-caller proof.

    Unit tests on `run_rule` say the wrapper works. They do not say the workflow
    uses it, and an emitter nothing invokes is precisely the defect this ticket
    was about. So assert the Snakefile routes every stage through it, and that
    no stage shells out to a script directly any more.
    """
    snakefile = (PIPELINE_ROOT / "workflow" / "Snakefile").read_text(encoding="utf-8")

    assert "run_rule.py" in snakefile, "the workflow does not reference run_rule.py"

    import re

    rules = re.split(r"\nrule ", snakefile)[1:]
    stage_rules = [
        r for r in rules
        if re.search(r"^    shell:", r, re.M) and "run_rule.py" in r
    ]
    assert len(stage_rules) >= 15, (
        f"expected every stage rule to be wrapped; found {len(stage_rules)}"
    )

    # no rule may declare two execution keywords: Snakemake refuses to load it
    for chunk in rules:
        name = chunk.split(":")[0].strip()
        keywords = [
            kw for kw in ("run:", "shell:", "script:", "notebook:", "wrapper:")
            if re.search(rf"^    {re.escape(kw)}", chunk, re.M)
        ]
        assert len(keywords) <= 1, f"rule {name} declares {keywords}"


def test_the_status_log_path_comes_from_the_machine_overlay():
    """Nothing environmental is hard-coded, including where the log goes."""
    from papipeline.config.loader import load_config

    config = load_config(
        PIPELINE_ROOT / "config" / "science.yaml", machine="laptop"
    )
    assert config.machine.status_dir().name == "status"
    assert config.machine.status_dir() == PIPELINE_ROOT / "status"
