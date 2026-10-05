"""The pilot report must present the census it measured, not one it remembers.

`pilot/report.py` and `pilot/pptx.py` both interpolated a literal figure of 835
assemblies, plus a breakdown (233 zero-byte, 51 truncated, 426 corrupt, 93
intact) that walks through how 835 became 65. Those numbers are a property of
one particular download, not of the run being reported. Re-run the pilot against
a different `data/` and the prose is a lie that still looks authoritative.

The census is already available: `run_pilot100.py` writes
`pilot100_prepare_summary.json` with `exclusion_reasons` from
`exclusion_breakdown`, and both builders receive that `summary` mapping. So the
fix is to read the mapping, not to route data anywhere.

Two properties this pins:

1. Every key in `exclusion_reasons` appears in the rendered prose, with its
   count. A reason that silently vanishes is worse than a wrong number, because
   the report then looks complete.
2. Keys are rendered through an explicit phrase mapping, never by reformatting
   the snake_case identifier. An unmapped key is still rendered - visibly, with
   its count - rather than dropped or folded into an anonymous total.

The fixture below is deliberately *not* the 835 run: different counts, a subset
of the real categories, and two keys outside the mapping. A swapped-constant
implementation cannot pass it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.pilot import pptx as pptx_mod
from papipeline.pilot import report as report_mod

SUMMARY = {
    "generated_utc": "2026-09-29 00:00:00",
    "assemblies_discovered": 812,
    "n_selected": 100,
    "n_excluded": 44,
    "n_analysable": 56,
    "pdc_matched": 100,
    "pdc_missing": 0,
    "amr_database_version": "test",
    "exclusion_reasons": {
        # real categories, counts that match nothing in the 835 run
        "zero_byte_file": 3,
        "corrupt_binary_data": 21,
        # NOT in the mapping: must still appear, with its count
        "assembly_unreadable_by_bakta": 4,
        # neither in the mapping nor a known category
        "other": 16,
    },
}

def _has_pptx() -> bool:
    import importlib.util

    return importlib.util.find_spec("pptx") is not None


TABLES: dict = {}
FIGURES: list = []
DETECTED: list = []
SKIPPED: list = []


@pytest.fixture()
def summary_file(tmp_path: Path) -> Path:
    path = tmp_path / "pilot100_prepare_summary.json"
    path.write_text(json.dumps(SUMMARY), encoding="utf-8")
    return path


def _report(summary_file: Path) -> str:
    summary = json.loads(summary_file.read_text(encoding="utf-8"))
    return report_mod.build_report(summary, TABLES, FIGURES, DETECTED, SKIPPED)


class TestExclusionRendering:
    def test_every_excluded_assembly_is_accounted_for_in_prose(self, summary_file):
        """The measured total must be the sum of the rendered reasons."""
        reasons = SUMMARY["exclusion_reasons"]
        assert sum(reasons.values()) == SUMMARY["n_excluded"]

    def test_each_reason_and_count_appears_in_the_prose(self, summary_file):
        """No key may vanish, and the sentence must carry the count next to it.

        Asserting the key appears *anywhere* is not enough: the exclusion table
        already printed raw keys, so such a test passes even if the sentence
        drops or rewrites them. The count has to be adjacent in the prose.
        """
        from papipeline.pilot.census import describe_reasons

        text = _report(summary_file)
        for key, count in SUMMARY["exclusion_reasons"].items():
            sentence = describe_reasons({key: count})
            assert sentence in text, f"{sentence!r} absent from the report"

    def test_unmapped_key_is_named_in_the_prose(self, summary_file):
        """An unmapped key stays visible; it is not dropped or aggregated."""
        from papipeline.pilot.census import describe_reasons

        assert "assembly_unreadable_by_bakta" in describe_reasons(
            {"assembly_unreadable_by_bakta": 4}
        )
        assert "assembly_unreadable_by_bakta" in _report(summary_file)

    def test_mapped_keys_use_the_phrase_not_the_identifier(self):
        """Rendering comes from the explicit table, never from the identifier."""
        from papipeline.pilot.census import describe_reasons

        rendered = describe_reasons({"corrupt_binary_data": 2, "zero_byte_file": 1})
        assert rendered == "corrupt: 2, zero-byte: 1"

    def test_prose_accounts_for_every_excluded_assembly(self):
        """The rendered phrases must sum to the reported excluded total."""
        from papipeline.pilot.census import describe_reasons

        assert sum(SUMMARY["exclusion_reasons"].values()) == SUMMARY["n_excluded"]
        rendered = describe_reasons(SUMMARY["exclusion_reasons"])
        for key, count in SUMMARY["exclusion_reasons"].items():
            assert f": {count}" in rendered or f": {count:,}" in rendered

    def test_no_literal_assembly_total_from_the_835_run(self, summary_file):
        """The remembered figure must not survive anywhere in the prose."""
        text = _report(summary_file)
        assert "835" not in text

    def test_discovered_total_is_read_from_the_summary(self, summary_file):
        """The denominator comes from the run, not from a constant."""
        text = _report(summary_file)
        assert str(SUMMARY["assemblies_discovered"]) in text

    def test_selection_phrase_reads_the_total_from_the_summary(self):
        """The 'N of M' phrase takes both figures as arguments, not constants."""
        from papipeline.pilot.census import describe_selection

        assert describe_selection(100, 812) == "100 of 812 available assemblies"

    @pytest.mark.skipif(
        pytest.importorskip.__module__ is None or not _has_pptx(),
        reason="python-pptx is not installed in this environment",
    )
    def test_deck_agrees_with_the_report(self, summary_file, tmp_path):
        """Both renderers must read the same mapping; no duplicate prose."""
        from pptx import Presentation

        out = tmp_path / "deck.pptx"
        pptx_mod.build_deck(out, SUMMARY, TABLES, FIGURES)
        text = "\n".join(
            shape.text_frame.text
            for slide in Presentation(str(out)).slides
            for shape in slide.shapes
            if shape.has_text_frame
        )
        assert "835" not in text
        assert str(SUMMARY["assemblies_discovered"]) in text
