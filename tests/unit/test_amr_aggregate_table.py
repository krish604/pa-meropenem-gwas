"""A REAL `amr` run must write the aggregate table, not only reports/.

**The gap this closes.** `_run_by_calling_tool` returned per-isolate records and
wrote AMRFinderPlus reports into `reports/`, but never wrote
`intermediate/amr/amr_determinants.tsv`. Two consumers depended on that file:

* `load_amr_table` reads it - the TEST path, satisfied by the fixture generator;
* the Snakemake rule declares it as `amr`'s input, via `AMR_IN`.

So in REAL it existed for neither, and the DAG could not be built at all:

    MissingInputException in rule amr
      results/real/intermediate/amr/amr_determinants.tsv

A completed REAL run had produced ten genuine AMRFinderPlus reports into a
directory nothing read. This is the same shape as the stage-2 annotation defect:
a tool-output table with a producer in one mode only.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.models import AmrDeterminant, ClaimStatus
from papipeline.stages import amr as stage_amr


def _determinant(sample_id="PDT_A", determinant="blaOXA-1"):
    return AmrDeterminant(
        sample_id=sample_id,
        antibiotic="imipenem",
        determinant=determinant,
        gene=determinant,
        variant=None,
        determinant_type="AMR",
        mechanism="beta-lactam",
        evidence_source="amrfinderplus",
        database="AMRFinderPlus",
        database_version="2026-08-07.1",
        confidence=99.5,
        claim_status=ClaimStatus.DETECTED,
        identity_pct=99.9,
        coverage_pct=100.0,
    )


class TestTheAggregateTableIsWritten:
    def test_it_carries_the_columns_the_loader_requires(self, tmp_path):
        """`REQUIRED` is what `load_amr_table` validates against.

        A table that exists but lacks these columns is the same failure as no
        table: the DAG resolves and the read fails part-way through a run.
        """
        out = tmp_path / "amr_determinants.tsv"
        stage_amr._write_determinant_table(
            out, {"PDT_A": [_determinant()], "PDT_B": [_determinant("PDT_B")]}
        )
        header = out.read_text().splitlines()[0].split("\t")
        missing = [c for c in stage_amr.REQUIRED if c not in header]
        assert not missing, f"the aggregate table omits required columns: {missing}"

    def test_every_detected_determinant_reaches_the_table(self, tmp_path):
        out = tmp_path / "amr_determinants.tsv"
        grouped = {
            "PDT_A": [_determinant("PDT_A", "blaOXA-1"), _determinant("PDT_A", "ampC")],
            "PDT_B": [_determinant("PDT_B", "OprD")],
        }
        stage_amr._write_determinant_table(out, grouped)
        assert len(out.read_text().strip().splitlines()) == 4  # header + 3

    def test_an_empty_screen_still_writes_the_header(self, tmp_path):
        """A screen that found nothing must differ from a screen never run.

        Silently skipping the write makes "no determinants" and "stage 4 did not
        execute" the same file state, and the second is a failure the DAG should
        surface rather than absorb.
        """
        out = tmp_path / "amr_determinants.tsv"
        stage_amr._write_determinant_table(out, {"PDT_A": []})
        text = out.read_text()
        assert out.exists()
        assert text.strip(), "header must be written even with no rows"
        for column in stage_amr.REQUIRED:
            assert column in text.splitlines()[0].split("\t")

    def test_it_creates_missing_parent_directories(self, tmp_path):
        out = tmp_path / "nested" / "deeper" / "amr_determinants.tsv"
        stage_amr._write_determinant_table(out, {"PDT_A": []})
        assert out.exists()

    def test_the_table_is_readable_by_the_real_loader(self, tmp_path):
        """Round trip, not just a header check.

        The point of writing it with the loader's exact columns is that
        `load_amr_table` can read it. Asserting the file's shape alone would pass
        with a column order or naming the loader rejects.
        """
        out = tmp_path / "amr_determinants.tsv"
        stage_amr._write_determinant_table(
            out, {"PDT_A": [_determinant("PDT_A", "ampC")]}
        )
        loaded = stage_amr.load_amr_table(
            out,
            antibiotic="imipenem",
            strict_mapping=False,
            mechanism_map={},
        )
        assert len(loaded) == 1
        assert loaded[0].determinant == "ampC"

class TestTheCallSiteIsWired:
    """The tests above exercise the writer directly, so they pass whether or not
    the REAL path calls it. That is not a hypothetical gap: removing the call
    from `_run_by_calling_tool` left all five green, and the call is the entire
    reason the file exists in REAL.

    Asserted on the source, which is the pattern the stage-taxonomy and
    stand-in tests already use for exactly this problem.
    """

    def test_the_real_path_writes_the_aggregate(self):
        from pathlib import Path as _Path

        source = (
            _Path(__file__).resolve().parents[2] / "papipeline" / "stages" / "amr.py"
        ).read_text(encoding="utf-8")
        body = source.split("def _run_by_calling_tool(")[1].split("\ndef ")[0]
        assert "_write_determinant_table(" in body, (
            "_run_by_calling_tool no longer writes amr_determinants.tsv, so a REAL "
            "run produces reports/ that nothing in the DAG reads and the declared "
            "AMR_IN input exists only in TEST"
        )

    def test_it_is_written_before_the_return(self):
        from pathlib import Path as _Path

        source = (
            _Path(__file__).resolve().parents[2] / "papipeline" / "stages" / "amr.py"
        ).read_text(encoding="utf-8")
        body = source.split("def _run_by_calling_tool(")[1].split("\ndef ")[0]
        assert body.index("_write_determinant_table(") < body.rindex("return grouped"), (
            "the aggregate must be written before the records are returned, or an "
            "early return leaves the file unwritten"
        )
