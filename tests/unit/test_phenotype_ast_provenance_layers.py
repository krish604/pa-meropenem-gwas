"""Which layer holds `ast_method` / `ast_standard` / `ast_edition`? Both.

There was a live disagreement about this, and both claims were repeated as
settled fact. Each names a **different layer**, and each is true of the layer
it names, which is why neither could be dismissed:

* "silently dropped at the write by `extrasaction='ignore'`" - true of the
  **file**;
* "stored on the model and exported" - true of the **model**.

So this file probes both layers and pins each, rather than picking a winner.
The pair of tests is the answer: the values are in memory and gone from disk,
and which of those a consumer sees depends entirely on whether it is handed
the `PhenotypeCall` objects or the table.

This matters beyond tidiness. ``provenance_report`` reads those three fields to
decide whether the cohort's calls may be pooled, and it can only do so because
stage 11 holds them in memory. Once written to ``11_phenotype.tsv`` they are
gone, so nothing downstream of the write can audit the gap - including any
consumer that reads the table back. A future change that adds the columns to
the contract will not break anything; it will make a second, on-disk copy of
the provenance exist. That is worth noticing rather than discovering.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.execution.contracts import STAGE_TABLES
from papipeline.io.tsv import read_tsv, write_tsv
from papipeline.models import Phenotype, PhenotypeCall

#: The stage-11 principal table: filename and the header, declared once.
PHENOTYPE_TABLE, PHENOTYPE_COLUMNS = STAGE_TABLES["phenotype"]

AST_PROVENANCE_FIELDS = ("ast_method", "ast_standard", "ast_edition")


def _call(**kw) -> PhenotypeCall:
    """A call with every AST provenance field POPULATED.

    Populated on purpose. A probe built from the real smoke data would find
    every value empty and could not tell "the model drops them" apart from
    "the source never supplied them" - which is the whole question.
    """
    base = dict(
        sample_id="GCA_000000001.1",
        antibiotic="imipenem",
        phenotype=Phenotype.R,
        mic=8.0,
        mic_unit="mg/L",
        source="probe",
    )
    base.update(kw)
    return PhenotypeCall(
        ast_method=base.pop("ast_method", "broth microdilution"),
        ast_standard=base.pop("ast_standard", "CLSI M100"),
        ast_edition=base.pop("ast_edition", "35th ed."),
        **base,
    )


class TestTheModelHoldsThem:
    """Layer 1 - `PhenotypeCall`. Both claims about the model are true."""

    def test_the_fields_exist_and_are_populated(self):
        call = _call()
        assert call.ast_method == "broth microdilution"
        assert call.ast_standard == "CLSI M100"
        assert call.ast_edition == "35th ed."

    @pytest.mark.parametrize("field", AST_PROVENANCE_FIELDS)
    def test_to_row_exports_each_field(self, field):
        row = _call().to_row()
        assert field in row, (
            f"{field} is a dataclass field but `to_row()` omits it, so the "
            "value would be lost even to an in-memory consumer"
        )
        assert row[field]

    def test_to_row_carries_them_alongside_the_contract_columns(self):
        """The export is wider than the table. Stated, because the gap between
        the two is the whole subject of this file."""
        exported = set(_call().to_row())
        declared = set(PHENOTYPE_COLUMNS)
        assert set(AST_PROVENANCE_FIELDS) <= exported
        assert not set(AST_PROVENANCE_FIELDS) & declared


class TestTheFileDropsThem:
    """Layer 2 - `11_phenotype.tsv`. The other claim is true here too."""

    def test_the_contract_header_omits_the_ast_columns(self):
        """The declaration, not the write, is where they are lost.

        `io.tsv.write_tsv` does not choose columns; its caller passes them.
        `run.py:1116-1125` passes exactly this list, so the omission is here in
        the contract and the two cannot drift apart unnoticed.
        """
        assert not [
            c for c in PHENOTYPE_COLUMNS if c.startswith("ast_")
        ], (
            "the stage-11 contract has grown ast_* columns. If that is "
            "deliberate, `provenance_report`'s on-disk consumers now exist and "
            "this file's premise is out of date"
        )

    def test_writing_a_row_drops_them(self, tmp_path: Path):
        out = write_tsv(
            tmp_path / PHENOTYPE_TABLE, [_call().to_row()], list(PHENOTYPE_COLUMNS)
        )
        header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header == list(PHENOTYPE_COLUMNS)
        assert not [c for c in header if c.startswith("ast_")]

    def test_reading_the_table_back_cannot_recover_them(self, tmp_path: Path):
        """The consequence, stated as a test rather than a comment.

        `provenance_report` counts provenance gaps from `PhenotypeCall` fields.
        Anything that reads the table instead sees a row with no provenance at
        all and cannot tell that from a source that recorded none - the two
        are the same on disk.
        """
        out = write_tsv(
            tmp_path / PHENOTYPE_TABLE, [_call().to_row()], list(PHENOTYPE_COLUMNS)
        )
        row = read_tsv(out, required_columns=("sample_id",))[0]
        assert not [k for k in row if k.startswith("ast_")]
        assert row["MIC"] == "8", "the contract columns survive, so this is a loss"

    def test_the_write_is_the_drop_not_the_reader(self, tmp_path: Path):
        """`extrasaction="ignore"` is the mechanism, named.

        Without it `csv.DictWriter` would raise on an unlisted key, and the
        stage would fail loudly instead of writing a narrower table. Silent is
        the hazard: the run succeeds and the provenance is gone.
        """
        import csv

        with pytest.raises(ValueError):
            csv.DictWriter(
                open(tmp_path / "x.tsv", "w", newline=""),
                fieldnames=list(PHENOTYPE_COLUMNS),
                delimiter="\t",
            ).writerow(_call().to_row())


class TestTheRealSmokeTableCarriesNone:
    """Layer 0 - the source. So the in-memory values are empty regardless."""

    smoke = Path(__file__).resolve().parents[2] / "db" / "smoke_phenotypes" / (
        "imipenem_phenotype.tsv"
    )

    @pytest.mark.skipif(
        not smoke.is_file(),
        reason="db/smoke_phenotypes/ is uncommitted by AGENTS.md rule 4 - "
               "absent in a fresh clone",
    )
    def test_the_source_has_no_ast_columns_at_all(self):
        rows = read_tsv(self.smoke, required_columns=("sample_id",))
        assert rows
        assert not [k for k in rows[0] if k.startswith("ast_")]