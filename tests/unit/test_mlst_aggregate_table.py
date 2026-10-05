"""Stage 3 must write the aggregate the DAG declares as its input.

**The gap.** `load_mlst` reads `intermediate/mlst/mlst_results.tsv`, and the
Snakemake rule declares it as `MLST_IN`. The REAL branch of `_run_by_calling_tool`
returned per-isolate records and never wrote it, so the file had producers in
TEST only - the synthetic fixture generator and a legacy script - and none in
REAL. Structurally identical to the stage-4 `amr_determinants.tsv` gap fixed in
69eb707, and the sixth instance of that class found today.

It stayed hidden because the TEST path is satisfied by the fixture generator, so
no test noticed, and no REAL run had reached stage 3 since.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.models import MlstCall
from papipeline.stages import mlst as stage_mlst


def _call(sample_id="PDT_A", st="235", alleles=None):
    return MlstCall(
        sample_id=sample_id,
        sequence_type=st,
        alleles=alleles if alleles is not None else {k: str(i) for i, k in enumerate(["acsA","aroE","guaA","mutL","nuoD","ppsA","trpE"], 1)},
        mlst_status="typed",
        mlst_scheme="pseudomonas_aeruginosa",
        allele_database="pubmlst:pseudomonas_aeruginosa",
    )


class TestTheAggregateIsWritten:
    def test_it_carries_the_columns_load_mlst_validates(self, tmp_path):
        """`REQUIRED` is what `load_mlst` checks. A file that exists but lacks
        these is the same failure as no file: the DAG resolves and the read
        fails part-way through a run."""
        out = tmp_path / "mlst_results.tsv"
        stage_mlst._write_mlst_results(out, {"PDT_A": _call()})
        header = out.read_text().splitlines()[0].split("\t")
        missing = [c for c in stage_mlst.REQUIRED if c not in header]
        assert not missing, f"aggregate omits required columns: {missing}"

    def test_it_is_readable_by_the_real_loader(self, tmp_path, config):
        """Round trip, not just a header check.

        The point of matching the loader's columns is that `load_mlst` can read
        it; asserting the shape alone would pass with a column order it rejects.
        """
        out = tmp_path / "mlst_results.tsv"
        stage_mlst._write_mlst_results(out, {"PDT_A": _call(st="235")})
        loaded = stage_mlst.load_mlst(config, out)
        assert isinstance(loaded, dict) and len(loaded) == 1
        assert next(iter(loaded.values())).sequence_type == "235"

    def test_an_empty_screen_still_writes_the_header(self, tmp_path):
        """"Typed nothing" must differ from "stage 3 did not run".

        Skipping the write makes those the same file state, and the second is a
        failure the DAG should surface rather than absorb.
        """
        out = tmp_path / "mlst_results.tsv"
        stage_mlst._write_mlst_results(out, {"PDT_A": None})
        assert out.exists() and out.read_text().strip(), "header must be written"

    def test_an_untyped_sample_is_written_with_its_status(self, tmp_path):
        """A cohort member with no sequence type is still a cohort member."""
        untyped = MlstCall(
            sample_id="PDT_B", sequence_type="", alleles={},
            mlst_status="no_call", mlst_scheme="pseudomonas_aeruginosa",
            allele_database="pubmlst:pseudomonas_aeruginosa",
        )
        out = tmp_path / "mlst_results.tsv"
        stage_mlst._write_mlst_results(out, {"PDT_B": untyped})
        assert "PDT_B" in out.read_text()

    def test_it_creates_missing_parent_directories(self, tmp_path):
        out = tmp_path / "a" / "b" / "mlst_results.tsv"
        stage_mlst._write_mlst_results(out, {"PDT_A": _call()})
        assert out.exists()


class TestTheCallSiteIsWired:
    """Exercising the writer proves nothing if REAL never calls it."""

    def test_the_real_path_writes_the_aggregate(self):
        from pathlib import Path as _Path

        source = (
            _Path(__file__).resolve().parents[2] / "papipeline" / "stages" / "mlst.py"
        ).read_text(encoding="utf-8")
        body = source.split("def _run_by_calling_tool(")[1].split("\ndef ")[0]
        assert "_write_mlst_results(" in body, (
            "_run_by_calling_tool no longer writes mlst_results.tsv, so a REAL "
            "run produces records the DAG cannot read and the declared MLST_IN "
            "input exists only in TEST"
        )

    def test_it_is_written_before_the_return(self):
        from pathlib import Path as _Path

        source = (
            _Path(__file__).resolve().parents[2] / "papipeline" / "stages" / "mlst.py"
        ).read_text(encoding="utf-8")
        body = source.split("def _run_by_calling_tool(")[1].split("\ndef ")[0]
        assert body.index("_write_mlst_results(") < body.rindex("return calls")