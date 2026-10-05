"""Tests for the panaroo adapter: the preflight and the presence/absence parser.

Two failure modes here are invisible if you do not look for them, and both are
tested deliberately rather than incidentally.

1. **Quoted commas in the Annotation column.** panaroo writes

       group_1375,C4-dicarboxylate ABC transporter,"protein binding, ATP binding",A,B,C

   A naive ``line.split(",")`` misaligns every field after the first quoted gene,
   so the per-genome cells get read from the wrong columns. The row count is
   unchanged and the header is unchanged, so nothing looks wrong - the counts are
   just false. This exact failure happened during the investigation that proved
   panaroo works on arm64: a shell one-liner using ``awk -F','`` reported
   ``core=0, specific=9324`` for a cohort whose true figures were 4,828 core and
   0 genome-specific. It looked exactly like a catastrophic clustering failure
   and would have caused a working tool to be written off. Hence
   ``csv.reader``, hence ``test_quoted_commas_do_not_shift_the_presence_columns``.

2. **panaroo importable but cd-hit absent.** panaroo shells out to cd-hit, so
   checking only the import would pass a preflight for a run that cannot succeed.
   Hence both checks, and ``test_the_preflight_names_cd_hit_separately``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import panaroo as adapter
from papipeline.errors import DataContractError, ToolNotAvailableError

THREE = ["S1", "S2", "S3"]

HEADER = "Gene,Non-unique Gene name,Annotation," + ",".join(THREE)


def _csv(rows: list[str]) -> str:
    return "\n".join([HEADER] + rows) + "\n"


class TestPreflight:
    def test_it_refuses_when_panaroo_is_not_importable(self, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def blocked(name, *args, **kwargs):
            if name == "panaroo":
                raise ImportError("no module named panaroo")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", blocked)
        monkeypatch.setattr(adapter.shutil, "which", lambda name: "/usr/bin/cd-hit")

        with pytest.raises(ToolNotAvailableError) as excinfo:
            adapter.preflight()
        assert "panaroo" in str(excinfo.value)

    def test_it_refuses_when_cd_hit_is_not_on_path(self, monkeypatch):
        """The case that a single check would miss: panaroo imports fine."""
        monkeypatch.setattr(adapter.shutil, "which", lambda name: None)
        # Ensure the panaroo import itself succeeds so only cd-hit is missing.
        monkeypatch.setitem(__import__("sys").modules, "panaroo", object())

        with pytest.raises(ToolNotAvailableError) as excinfo:
            adapter.preflight()
        assert "cd-hit" in str(excinfo.value)

    def test_it_names_both_when_both_are_missing(self, monkeypatch):
        """Not an either/or. A message naming one tool sends the reader off to
        fix one tool, fixes it, and hits the same wall."""
        import builtins

        real_import = builtins.__import__

        def blocked(name, *args, **kwargs):
            if name == "panaroo":
                raise ImportError("no module named panaroo")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", blocked)
        monkeypatch.setattr(adapter.shutil, "which", lambda name: None)

        with pytest.raises(ToolNotAvailableError) as excinfo:
            adapter.preflight()
        message = str(excinfo.value)
        assert "panaroo" in message and "cd-hit" in message

    def test_it_says_panaroo_alone_is_not_enough(self, monkeypatch):
        """cd-hit missing while panaroo imports - the message must not imply
        installing panaroo would fix it."""
        monkeypatch.setattr(adapter.shutil, "which", lambda name: None)
        monkeypatch.setitem(__import__("sys").modules, "panaroo", object())

        with pytest.raises(ToolNotAvailableError) as excinfo:
            adapter.preflight()
        assert "cd-hit" in str(excinfo.value)


class TestPresenceColumnDetection:
    def test_it_locates_genome_columns_past_the_three_metadata_columns(self):
        cols = adapter._presence_columns(["Gene", "Non-unique Gene name", "Annotation"] + THREE)
        assert cols == [3, 4, 5]

    def test_it_is_case_insensitive_on_metadata_names(self):
        cols = adapter._presence_columns(["gene", "NON-UNIQUE GENE NAME", "annotation"] + THREE)
        assert cols == [3, 4, 5]

    def test_it_refuses_when_there_are_no_genome_columns(self):
        with pytest.raises(DataContractError):
            adapter._presence_columns(["Gene", "Non-unique Gene name", "Annotation"])


class TestParsePresenceCsv:
    def test_it_returns_sample_ids_and_a_presence_map(self, tmp_path):
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(
            _csv([
                "g1,,,1,1,0",
                "g2,,,1,0,0",
            ])
        )
        samples, presence = adapter.parse_presence_csv(path)
        assert samples == tuple(THREE)
        assert presence["g1"] == {"S1", "S2"}
        assert presence["g2"] == {"S1"}

    def test_quoted_commas_do_not_shift_the_presence_columns(self, tmp_path):
        """The core=0/specific=9324 trap, as a regression test.

        A gene whose annotation contains a comma is the trigger. If anything
        parses this with split(",") instead of csv.reader, g2 loses its real
        carriers and gains phantom ones - silently, with the row count intact.
        """
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(
            _csv([
                'g1,"ABC transporter","protein binding, ATP binding",1,1,0',
                'g2,,,1,0,0',
            ])
        )
        samples, presence = adapter.parse_presence_csv(path)
        assert samples == tuple(THREE), "header must still yield three genomes"
        assert presence["g1"] == {"S1", "S2"}, (
            "a comma inside the Annotation column shifted the genome columns"
        )
        assert presence["g2"] == {"S1"}

    def test_an_absent_cell_is_not_presence(self, tmp_path):
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(_csv(["g1,,,0,,0"]))
        _, presence = adapter.parse_presence_csv(path)
        assert presence["g1"] == set()

    def test_an_empty_file_is_refused(self, tmp_path):
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text("")
        with pytest.raises(DataContractError):
            adapter.parse_presence_csv(path)

    def test_a_header_only_file_yields_no_genes_rather_than_failing(self, tmp_path):
        """Not a crash: panaroo can legitimately emit a header with no rows, and
        an empty cohort is a fact, not a malformed input."""
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(HEADER + "\n")
        samples, presence = adapter.parse_presence_csv(path)
        assert samples == tuple(THREE)
        assert presence == {}


class TestBuildPangenomeFromPanaroo:
    def test_it_refuses_a_cohort_mismatch(self, tmp_path):
        """A missing isolate must fail loudly.

        A pangenome over 9 of 10 isolates reports different core/accessory
        boundaries while every number still looks plausible, so a silent drop
        here would quietly change the science.
        """
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(_csv(["g1,,,1,1,0"]))  # only S1..S3 present

        with pytest.raises(DataContractError) as excinfo:
            adapter.build_pangenome_from_panaroo(tmp_path, ["S1", "S2", "S4"])
        message = str(excinfo.value)
        assert "S4" in message, "the isolate that went missing must be named"

    def test_it_names_an_unexpected_sample_too(self, tmp_path):
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(_csv(["g1,,,1,1,0"]))
        with pytest.raises(DataContractError) as excinfo:
            adapter.build_pangenome_from_panaroo(tmp_path, ["S1", "S2"])
        assert "S3" in str(excinfo.value)

    def test_it_refuses_when_panaroo_wrote_no_table(self, tmp_path):
        """A missing table means the tool failed, not that the cohort was empty."""
        with pytest.raises(DataContractError) as excinfo:
            adapter.build_pangenome_from_panaroo(tmp_path, THREE)
        assert adapter.PRESENCE_CSV in str(excinfo.value)

    def test_it_accepts_a_matching_cohort(self, tmp_path):
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(_csv(["g1,,,1,1,0", "g2,,,1,1,1"]))
        samples, presence = adapter.build_pangenome_from_panaroo(tmp_path, THREE)
        assert samples == tuple(THREE)
        assert presence["g2"] == set(THREE)
