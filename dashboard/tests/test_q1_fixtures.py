"""Q1 — the fixtures generator.

The generated small fixtures are committed; the large one is not. These tests
assert the properties the rest of the suite depends on: contract headers, a
duplicate-support tree, synthetic ids only, and the bundle's A1–A5 directories.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

from conftest import FIXTURES_DIR, SMALL_BUNDLE, SMALL_LIVE
from dashboard.fixtures import generate as gen

from papipeline.execution.contracts import INTERNAL_TABLES, STAGE_TABLES, table_path


def _header(path: Path) -> list[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        if line.strip():
            return line.split("\t")
    return []


def test_small_live_and_bundle_are_committed():
    assert SMALL_LIVE.is_dir()
    assert SMALL_BUNDLE.is_dir()
    assert (SMALL_LIVE / "run_manifest.json").exists()
    assert (SMALL_BUNDLE / "02_stage_outputs" / "run_manifest.json").exists()


@pytest.mark.parametrize("stage", sorted(STAGE_TABLES))
def test_live_stage_tables_match_the_contract(stage, small_live):
    _filename, header = STAGE_TABLES[stage]
    path = table_path(small_live / "intermediate" / "stages", stage)
    assert path.exists(), path
    assert _header(path) == list(header)


@pytest.mark.parametrize("key", sorted(INTERNAL_TABLES))
def test_live_internal_tables_match_the_contract(key, small_live):
    from papipeline.execution.contracts import internal_table_path

    _filename, header = INTERNAL_TABLES[key]
    path = internal_table_path(small_live / "intermediate" / "stages", key)
    assert path.exists(), path
    assert _header(path) == list(header)


def test_small_fixtures_use_synthetic_ids_only():
    pattern = re.compile(r"\bPDT\d|\bGCF_\d|SAMN\d|SAMEA\d")
    for path in FIXTURES_DIR.rglob("*"):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        assert not pattern.search(text), f"real accession-like id in {path}"


def test_tree_has_duplicate_support_labels(small_live):
    tree = (small_live / "intermediate" / "phylogeny" / "stage9.nwk").read_text(encoding="utf-8")
    assert tree.count("100/100") >= 3
    # A tip label must survive byte-for-byte; underscores are part of the id.
    assert "TEST_PA_001" in tree
    assert "TEST PA 001" not in tree


def test_bundle_has_the_assumed_layout(small_bundle):
    for name in ("01_bakta_input", "02_stage_outputs", "03_report", "04_run_info",
                 "05_validation", "06_for_900_isolates"):
        assert (small_bundle / name).is_dir(), name
    assert (small_bundle / "06_for_900_isolates" / "RUNBOOK_900.md").exists()
    assert (small_bundle / "RESULTS.md").exists()


@pytest.mark.skipif(
    os.environ.get("PA_FIXTURES_LARGE") == "1",
    reason=(
        "this test guards the NORMAL run: its first assertion is that "
        "PA_FIXTURES_LARGE is unset. The full suite is also run with "
        "PA_FIXTURES_LARGE=1 to exercise the perf gates, and that run "
        "deliberately sets the variable this test guards against."
    ),
)
def test_large_build_is_gated_and_never_committed(tmp_path):
    assert os.environ.get("PA_FIXTURES_LARGE") != "1", (
        "PA_FIXTURES_LARGE must not be set for a normal run"
    )
    root = gen.build_large(tmp_path / "large", n=900)
    assert (root / "intermediate" / "stages" / "03_mlst.tsv").exists()
    assert not str(root).startswith(str(FIXTURES_DIR))
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_samples"] == 900
