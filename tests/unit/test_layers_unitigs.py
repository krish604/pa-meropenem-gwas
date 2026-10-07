"""The unitig adapter: stub in TEST, refusal in REAL, flags never invented.

``unitig-caller`` is the tool this package is allowed to wrap, and AGENTS.md
rule 1 says its flags are verified with ``--help`` rather than assumed. The
tool is not installed on this machine, so :data:`FLAGS_VERIFIED` is ``False``
and these tests pin the three consequences:

* presence is *reported* by the adapter without executing anything;
* TEST and STUB runs write a stub table carrying the isolate set and no unitig
  features - the stub proves the interface, not the biology;
* a REAL run refuses, and the refusal says which command was never run, so the
  next reader knows exactly what has to be verified before this adapter can
  invoke anything.

Synthetic only: sample IDs are hand-built and nothing under ``data/``, ``db/``
or ``PDC_essential.tsv`` is touched.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import pytest

from papipeline.adapters.external import UNPROBED_VERSION
from papipeline.errors import PipelineError
from papipeline.io.tsv import read_tsv
from papipeline.layers.unitigs import (
    FLAGS_VERIFIED,
    UNITIG_TABLE_FILENAME,
    UNITIG_TOOL,
    UNVERIFIED_REASON,
    VERIFIED_FLAGS,
    UnverifiedFlagsError,
    run_unitig_screen,
    unitig_tool_status,
)
from papipeline.models import RunMode

SAMPLES: Tuple[str, ...] = tuple(f"TEST_PA_{i:03d}" for i in range(1, 6))


def _banner(path: Path) -> List[str]:
    return [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("#")
    ]


def _header_line(path: Path) -> List[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        return line.split("\t")
    raise AssertionError(f"{path} has no header row")


class TestDetectionNeverExecutesAnything:
    def test_absence_is_reported_without_running_a_command(self, monkeypatch):
        monkeypatch.setattr(
            "papipeline.adapters.external.shutil.which", lambda name: None
        )
        status = unitig_tool_status()
        assert status.name == UNITIG_TOOL
        assert status.available is False
        assert status.executable is None
        assert status.version is None

    def test_presence_is_reported_and_the_version_left_unprobed(self, monkeypatch):
        monkeypatch.setattr(
            "papipeline.adapters.external.shutil.which",
            lambda name: "/opt/bin/unitig-caller",
        )
        status = unitig_tool_status()
        assert status.available is True
        # A presence sweep never runs a binary it has not been told to run.
        assert status.version == UNPROBED_VERSION


class TestTheFlagStateIsRecordedHonestly:
    def test_no_flag_is_claimed_as_verified(self):
        assert FLAGS_VERIFIED is False
        assert VERIFIED_FLAGS == ()

    def test_the_reason_names_the_command_that_could_not_be_run(self):
        assert UNITIG_TOOL in UNVERIFIED_REASON
        assert "--help" in UNVERIFIED_REASON

    def test_the_reason_says_the_tool_is_absent_rather_than_inventing_flags(self):
        assert "not installed" in UNVERIFIED_REASON or "absent" in UNVERIFIED_REASON


class TestStubModes:
    @pytest.mark.parametrize("mode", [RunMode.TEST, RunMode.STUB])
    def test_a_stub_table_carrying_the_isolate_set_is_written(
        self, mode, tmp_path: Path
    ):
        path = run_unitig_screen(mode=mode, sample_ids=SAMPLES, output_dir=tmp_path)
        assert path == tmp_path / UNITIG_TABLE_FILENAME
        rows = read_tsv(path)
        assert [row["sample_id"] for row in rows] == list(SAMPLES)
        assert _header_line(path) == ["sample_id"]

    @pytest.mark.parametrize("mode", [RunMode.TEST, RunMode.STUB])
    def test_the_banner_says_it_is_a_stub_and_what_was_not_run(
        self, mode, tmp_path: Path
    ):
        path = run_unitig_screen(mode=mode, sample_ids=SAMPLES, output_dir=tmp_path)
        banner = "\n".join(_banner(path))
        assert "stub" in banner.lower()
        assert UNITIG_TOOL in banner
        assert "--help" in banner

    def test_a_stub_table_carries_no_unitig_features(self, tmp_path: Path):
        path = run_unitig_screen(
            mode=RunMode.TEST, sample_ids=SAMPLES, output_dir=tmp_path
        )
        rows = read_tsv(path)
        assert all(list(row) == ["sample_id"] for row in rows)


class TestRealModeIsRefused:
    def test_real_mode_refuses_to_invoke_an_unverified_command(
        self, tmp_path: Path
    ):
        with pytest.raises(UnverifiedFlagsError) as excinfo:
            run_unitig_screen(mode=RunMode.REAL, sample_ids=SAMPLES, output_dir=tmp_path)
        message = str(excinfo.value)
        assert UNITIG_TOOL in message
        assert "--help" in message
        # Nothing is written: a refusal that still leaves an artefact behind
        # would read as a screen that ran.
        assert not (tmp_path / UNITIG_TABLE_FILENAME).exists()

    def test_the_refusal_is_a_pipeline_error_so_a_run_stops_cleanly(self):
        assert issubclass(UnverifiedFlagsError, PipelineError)

    def test_the_refusal_quotes_the_recorded_reason(self):
        with pytest.raises(UnverifiedFlagsError) as excinfo:
            run_unitig_screen(mode=RunMode.REAL, sample_ids=SAMPLES, output_dir=Path("."))
        assert UNVERIFIED_REASON.split(".")[0] in str(excinfo.value)
