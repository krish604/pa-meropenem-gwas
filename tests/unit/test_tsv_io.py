"""Tests for the strict TSV reader and writer.

Covers the data-contract rules the whole pipeline depends on: malformed rows,
missing values, duplicate identifiers, composite keys and comment handling.
"""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError, DuplicateSampleError
from papipeline.io.tsv import read_tsv, write_fasta, write_lines, write_tsv


class TestReadTsvBasics:
    def test_reads_simple_table(self, write_tsv_file):
        path = write_tsv_file("a.tsv", "sample_id\tvalue\nS1\t1\nS2\t2\n")
        rows = read_tsv(path)
        assert rows == [{"sample_id": "S1", "value": "1"}, {"sample_id": "S2", "value": "2"}]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(DataContractError, match="not found"):
            read_tsv(tmp_path / "nope.tsv")

    def test_empty_file_raises(self, write_tsv_file):
        path = write_tsv_file("empty.tsv", "")
        with pytest.raises(DataContractError, match="empty"):
            read_tsv(path)

    def test_header_only_raises(self, write_tsv_file):
        path = write_tsv_file("header_only.tsv", "a\tb\n")
        with pytest.raises(DataContractError, match="no data rows"):
            read_tsv(path)

    def test_comments_only_raises(self, write_tsv_file):
        path = write_tsv_file("comments.tsv", "# just a comment\n# another\n")
        with pytest.raises(DataContractError, match="only comments"):
            read_tsv(path)

    def test_blank_lines_are_skipped(self, write_tsv_file):
        path = write_tsv_file("blanks.tsv", "a\tb\nS1\t1\n\n\nS2\t2\n")
        assert len(read_tsv(path)) == 2


class TestMissingValues:
    @pytest.mark.parametrize("sentinel", [".", "NA", "na", "nan", "None", "null", "-", ""])
    def test_sentinels_become_none(self, write_tsv_file, sentinel):
        path = write_tsv_file("m.tsv", f"sample_id\tvalue\nS1\t{sentinel}\n")
        rows = read_tsv(path)
        assert rows[0]["value"] is None, f"{sentinel!r} should be treated as missing"

    def test_real_values_are_preserved(self, write_tsv_file):
        path = write_tsv_file("v.tsv", "sample_id\tvalue\nS1\tNDM\n")
        assert read_tsv(path)[0]["value"] == "NDM"

    def test_short_row_pads_with_none(self, write_tsv_file):
        """A short row pads to None rather than dropping the column."""
        path = write_tsv_file("short.tsv", "a\tb\tc\n1\t2\n")
        rows = read_tsv(path)
        assert rows[0] == {"a": "1", "b": "2", "c": None}


class TestMalformedTsv:
    def test_too_many_fields_raises(self, write_tsv_file):
        path = write_tsv_file("long.tsv", "a\tb\n1\t2\t3\n")
        with pytest.raises(DataContractError, match="more fields than the header"):
            read_tsv(path)

    def test_too_many_fields_reports_line(self, write_tsv_file):
        path = write_tsv_file("long2.tsv", "a\tb\n1\t2\n3\t4\t5\n")
        with pytest.raises(DataContractError) as exc:
            read_tsv(path)
        assert exc.value.context["line"] == 3

    def test_duplicate_header_raises(self, write_tsv_file):
        path = write_tsv_file("dupheader.tsv", "a\ta\n1\t2\n")
        with pytest.raises(DataContractError, match="duplicate column names"):
            read_tsv(path)

    def test_missing_required_column_raises(self, write_tsv_file):
        path = write_tsv_file("req.tsv", "a\tb\n1\t2\n")
        with pytest.raises(DataContractError, match="missing required columns"):
            read_tsv(path, required_columns=("sample_id",))

    def test_min_columns_enforced(self, write_tsv_file):
        path = write_tsv_file("min.tsv", "a\tb\tc\n1\n")
        with pytest.raises(DataContractError, match="too few fields"):
            read_tsv(path, min_columns=2)


class TestUniqueness:
    def test_duplicate_single_column_raises(self, write_tsv_file):
        path = write_tsv_file("dup.tsv", "sample_id\tx\nS1\t1\nS1\t2\n")
        with pytest.raises(DuplicateSampleError, match="single-column"):
            read_tsv(path, unique_columns=("sample_id",))

    def test_duplicate_reports_both_lines(self, write_tsv_file):
        path = write_tsv_file("dup2.tsv", "sample_id\tx\nS1\t1\nS2\t2\nS1\t3\n")
        with pytest.raises(DuplicateSampleError) as exc:
            read_tsv(path, unique_columns=("sample_id",))
        assert exc.value.context["first_line"] == 2
        assert exc.value.context["second_line"] == 4

    def test_unique_column_on_absent_column_raises(self, write_tsv_file):
        path = write_tsv_file("absent.tsv", "a\tb\n1\t2\n")
        with pytest.raises(DataContractError, match="absent column"):
            read_tsv(path, unique_columns=("nope",))

    def test_composite_key_allows_repeated_single_column(self, write_tsv_file):
        """The same sample may appear once per antibiotic."""
        path = write_tsv_file(
            "composite.tsv",
            "sample_id\tantibiotic\tphenotype\nS1\timipenem\tR\nS1\tmeropenem\tS\n",
        )
        rows = read_tsv(path, unique_together=[("sample_id", "antibiotic")])
        assert len(rows) == 2

    def test_composite_key_duplicate_raises(self, write_tsv_file):
        path = write_tsv_file(
            "composite_dup.tsv",
            "sample_id\tantibiotic\tphenotype\nS1\timipenem\tR\nS1\timipenem\tS\n",
        )
        with pytest.raises(DuplicateSampleError, match="composite"):
            read_tsv(path, unique_together=[("sample_id", "antibiotic")])

    def test_partially_populated_composite_key_raises(self, write_tsv_file):
        """A key missing one component cannot be deduplicated safely."""
        path = write_tsv_file("partial.tsv", "sample_id\tantibiotic\nS1\t.\nS2\timipenem\n")
        with pytest.raises(DataContractError, match="Incompletely populated"):
            read_tsv(path, unique_together=[("sample_id", "antibiotic")])


class TestComments:
    def test_comment_lines_are_skipped(self, write_tsv_file):
        path = write_tsv_file(
            "c.tsv",
            "# provenance line\n# another\na\tb\n1\t2\n# trailing\n",
        )
        rows = read_tsv(path)
        assert rows == [{"a": "1", "b": "2"}]

    def test_comment_prefix_can_be_disabled(self, write_tsv_file):
        path = write_tsv_file("c2.tsv", "#a\tb\n1\t2\n")
        # With comments disabled the '#a' line becomes the header.
        rows = read_tsv(path, comment_prefix=None)
        assert rows[0]["#a"] == "1"


class TestWriteTsv:
    def test_round_trip_preserves_values(self, tmp_path):
        path = tmp_path / "out.tsv"
        write_tsv(path, [{"a": "1", "b": None}], ["a", "b"])
        rows = read_tsv(path)
        assert rows[0]["a"] == "1"
        assert rows[0]["b"] is None

    def test_column_order_is_stable(self, tmp_path):
        path = tmp_path / "order.tsv"
        write_tsv(path, [{"b": "2", "a": "1"}], ["a", "b"])
        header = path.read_text().splitlines()[0]
        assert header == "a\tb"

    def test_provenance_header_round_trips(self, tmp_path):
        """A file written with comments can be read back."""
        path = tmp_path / "prov.tsv"
        write_tsv(path, [{"a": "1"}], ["a"], header_comment=["generated", "seed=1"])
        assert path.read_text().startswith("# generated")
        assert read_tsv(path) == [{"a": "1"}]

    def test_missing_key_becomes_sentinel(self, tmp_path):
        path = tmp_path / "miss.tsv"
        write_tsv(path, [{"a": "1"}], ["a", "b"])
        assert "\t." in path.read_text().splitlines()[1]

    def test_lists_are_joined(self, tmp_path):
        path = tmp_path / "list.tsv"
        write_tsv(path, [{"a": ["x", "y"]}], ["a"])
        assert read_tsv(path)[0]["a"] == "x;y"

    def test_creates_parent_directories(self, tmp_path):
        path = tmp_path / "deep" / "nested" / "out.tsv"
        write_tsv(path, [{"a": "1"}], ["a"])
        assert path.exists()


class TestOtherWriters:
    def test_write_lines(self, tmp_path):
        path = write_lines(tmp_path / "l.txt", ["a", "b"])
        assert path.read_text() == "a\nb\n"

    def test_write_fasta_wraps(self, tmp_path):
        path = write_fasta(tmp_path / "f.fasta", [("h1", "ACGT" * 30)])
        lines = path.read_text().splitlines()
        assert lines[0] == ">h1"
        assert all(len(line) <= 80 for line in lines[1:])
