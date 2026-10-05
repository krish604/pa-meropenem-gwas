"""Tests for the FASTA reader and assembly summary statistics."""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError
from papipeline.io.fasta import (
    AMBIGUITY_CODES,
    assembly_stats,
    file_integrity,
    read_fasta,
    summarise_alignment,
)
from papipeline.io.tsv import write_fasta


@pytest.fixture
def simple_fasta(tmp_path):
    path = tmp_path / "g.fna"
    write_fasta(
        path,
        [
            ("contig_1 length=500", "ACGT" * 125),
            ("contig_2 length=300", "GGCC" * 75),
            ("contig_3 length=100", "ATAT" * 25),
        ],
    )
    return path


class TestReadFasta:
    def test_reads_records(self, simple_fasta):
        records = list(read_fasta(simple_fasta))
        assert len(records) == 3
        assert records[0].header == "contig_1 length=500"

    def test_record_id_is_first_token(self, simple_fasta):
        assert list(read_fasta(simple_fasta))[0].record_id == "contig_1"

    def test_sequence_length(self, simple_fasta):
        assert len(list(read_fasta(simple_fasta))[0]) == 500

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(DataContractError, match="not found"):
            list(read_fasta(tmp_path / "absent.fna"))

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "e.fna"
        path.write_text("")
        with pytest.raises(DataContractError, match="empty"):
            list(read_fasta(path))

    def test_sequence_before_header_raises(self, tmp_path):
        path = tmp_path / "bad.fna"
        path.write_text("ACGTACGT\n>h1\nACGT\n")
        with pytest.raises(DataContractError, match="does not start with a header"):
            list(read_fasta(path))

    def test_empty_header_raises(self, tmp_path):
        path = tmp_path / "empty_header.fna"
        path.write_text(">\nACGT\n")
        with pytest.raises(DataContractError, match="empty header"):
            list(read_fasta(path))

    def test_lowercase_is_normalised(self, tmp_path):
        path = tmp_path / "lower.fna"
        path.write_text(">h1\nacgtACGT\n")
        assert list(read_fasta(path))[0].sequence == "ACGTACGT"

    def test_blank_lines_are_ignored(self, tmp_path):
        path = tmp_path / "blank.fna"
        path.write_text(">h1\n\nACGT\n\nACGT\n")
        assert list(read_fasta(path))[0].sequence == "ACGTACGT"


class TestAssemblyStats:
    def test_total_size_and_counts(self, simple_fasta):
        summary = assembly_stats(simple_fasta)
        assert summary.assembly_size == 900
        assert summary.contig_count == 3
        assert summary.largest_contig == 500

    def test_n50_and_n90(self, simple_fasta):
        """Lengths 500/300/100, total 900.

        N50 needs a cumulative 450: the 500 contig alone reaches it.
        N90 needs a cumulative 810: 500+300=800 falls short, so all three
        contigs are required and N90 is the shortest, 100.
        """
        summary = assembly_stats(simple_fasta)
        assert summary.n50 == 500
        assert summary.n90 == 100

    def test_n50_single_contig(self, tmp_path):
        path = tmp_path / "one.fna"
        write_fasta(path, [("c1", "ACGT" * 100)])
        summary = assembly_stats(path)
        assert summary.n50 == 400
        assert summary.n90 == 400

    def test_gc_content(self, tmp_path):
        path = tmp_path / "gc.fna"
        write_fasta(path, [("c1", "GGGGCCCC")])
        assert assembly_stats(path).gc_content == 100.0

    def test_zero_gc_content(self, tmp_path):
        path = tmp_path / "at.fna"
        write_fasta(path, [("c1", "ATATATAT")])
        assert assembly_stats(path).gc_content == 0.0

    def test_ambiguous_bases_counted(self, tmp_path):
        path = tmp_path / "amb.fna"
        write_fasta(path, [("c1", "ACGTNNNNR")])
        summary = assembly_stats(path)
        assert summary.ambiguous_bases == 5

    def test_ambiguity_codes_are_the_iupac_set(self):
        assert set("NRYKMSWBDHV") <= AMBIGUITY_CODES

    def test_invalid_characters_reported(self, tmp_path):
        path = tmp_path / "bad_chars.fna"
        path.write_text(">c1\nACGTZZZ\n")
        assert assembly_stats(path).invalid_characters == 3

    def test_metrics_are_none_for_no_records(self, tmp_path):
        """Unmeasurable metrics are None, never 0."""
        path = tmp_path / "hdr_only.fna"
        path.write_text(">c1\n")
        summary = assembly_stats(path)
        assert summary.assembly_size == 0
        assert summary.n50 is None
        assert summary.gc_content is None

    def test_to_row_keys(self, simple_fasta):
        row = assembly_stats(simple_fasta).to_row()
        assert set(row) >= {"assembly_size", "contig_count", "n50", "n90", "gc_content"}


class TestSummariseAlignment:
    def test_counts_sequences_and_sites(self, tmp_path):
        path = tmp_path / "aln.fasta"
        write_fasta(path, [("s1", "ACGTAC"), ("s2", "ACGTAC"), ("s3", "ACGTAC")])
        summary = summarise_alignment(path)
        assert summary["n_sequences"] == 3
        assert summary["n_sites"] == 6
        assert summary["missing_fraction"] == 0.0

    def test_missing_fraction(self, tmp_path):
        """'AC-GT' is 5 sites with 1 gap; 2 sequences = 10 sites, 2 gaps."""
        path = tmp_path / "aln2.fasta"
        write_fasta(path, [("s1", "AC-GT"), ("s2", "AC-GT")])
        assert summarise_alignment(path)["missing_fraction"] == pytest.approx(0.2)

    def test_missing_fraction_all_gaps(self, tmp_path):
        path = tmp_path / "aln3.fasta"
        write_fasta(path, [("s1", "----"), ("s2", "----")])
        assert summarise_alignment(path)["missing_fraction"] == 1.0

    def test_ragged_alignment_reports_longest(self, tmp_path):
        path = tmp_path / "ragged.fasta"
        path.write_text(">s1\nACGT\n>s2\nACGTAC\n")
        assert summarise_alignment(path)["n_sites"] == 6

    def test_empty_alignment(self, tmp_path):
        path = tmp_path / "empty_aln.fasta"
        path.write_text("")
        summary = summarise_alignment(path)
        assert summary["n_sequences"] == 0
        assert summary["missing_fraction"] is None


class TestTruncatedAndCorruptFiles:
    """Real bulk-downloaded genome sets contain truncated and garbled files.

    These are the failure modes discovered when the pilot ran against the
    actual data/, where 46 of the first 100 assemblies proved unreadable.
    """

    def test_zero_byte_file_raises(self, tmp_path):
        path = tmp_path / "empty.fna"
        path.write_bytes(b"")
        with pytest.raises(DataContractError, match="empty"):
            list(read_fasta(path))

    def test_binary_garbage_raises_a_typed_error(self, tmp_path):
        """A non-UTF-8 byte must not escape as UnicodeDecodeError."""
        path = tmp_path / "corrupt.fna"
        path.write_bytes(b">contig_1\n" + b"ACGT" * 4 + b"\xbf\x0e-\xbf\x8e\x8d")
        with pytest.raises(DataContractError) as exc:
            list(read_fasta(path))
        assert "corrupt" in str(exc.value) or "non-text" in str(exc.value)

    def test_corrupt_error_names_the_offset(self, tmp_path):
        path = tmp_path / "corrupt2.fna"
        path.write_bytes(b">c\nACGT\n" + b"\xbf" * 3)
        with pytest.raises(DataContractError) as exc:
            list(read_fasta(path))
        assert "byte_offset" in exc.value.context

    def test_truncation_at_4mib_boundary(self, tmp_path):
        """The observed corruption signature: valid FASTA, then garbage."""
        path = tmp_path / "half.fna"
        path.write_bytes(b">c\n" + b"ACGT" * 64 + b"\n" + b"\xbf" * 32)
        with pytest.raises(DataContractError, match="corrupt"):
            list(read_fasta(path))

    def test_header_only_file_yields_one_empty_record(self, tmp_path):
        """A header with no sequence is degenerate, not corrupt.

        It yields a single zero-length record, so downstream QC sees size 0
        and fails the size threshold rather than raising.
        """
        path = tmp_path / "hdr.fna"
        path.write_bytes(b">c1\n")
        records = list(read_fasta(path))
        assert len(records) == 1
        assert records[0].sequence == ""

    def test_no_records_raises(self, tmp_path):
        path = tmp_path / "blank.fna"
        path.write_bytes(b"\n\n\n")
        with pytest.raises(DataContractError):
            list(read_fasta(path))

    def test_gzipped_is_detected_not_parsed(self, tmp_path):
        import gzip

        path = tmp_path / "gz.fna.gz"
        with gzip.open(path, "wb") as handle:
            handle.write(b">c\nACGT\n")
        info = file_integrity(path)
        assert info["is_gzipped"] is True
        assert info["size_bytes"] > 0


class TestFileIntegrity:
    def test_missing_file(self, tmp_path):
        info = file_integrity(tmp_path / "absent.fna")
        assert info["exists"] is False
        assert info["size_bytes"] is None

    def test_empty_file(self, tmp_path):
        path = tmp_path / "e.fna"
        path.write_bytes(b"")
        info = file_integrity(path)
        assert info["is_empty"] is True
        assert info["looks_like_fasta"] is False

    def test_valid_fasta(self, tmp_path):
        path = tmp_path / "g.fna"
        path.write_text(">c1 desc\nACGT\n")
        info = file_integrity(path)
        assert info["exists"] is True
        assert info["is_empty"] is False
        assert info["looks_like_fasta"] is True
        assert info["first_line"].startswith(">c1")

    def test_binary_is_not_fasta(self, tmp_path):
        path = tmp_path / "b.fna"
        path.write_bytes(b"\x00\x01\x02")
        assert file_integrity(path)["looks_like_fasta"] is False


class TestValidationQuarantinesBadFiles:
    """One bad file must not abort a cohort, and must not count as a pass."""

    def _qc(self, config, path):
        from papipeline.stages.validation import validate_sample

        return validate_sample(config, "S1", path)

    def test_unreadable_file_is_quarantined(self, config, tmp_path):
        path = tmp_path / "corrupt.fna"
        path.write_bytes(b">c\nACGT\n" + b"\xbf" * 8)
        record = self._qc(config, path)
        assert record.basic_quality_status == "unreadable_input"
        assert "unreadable" in record.flag_reasons

    def test_empty_file_is_quarantined(self, config, tmp_path):
        path = tmp_path / "empty.fna"
        path.write_bytes(b"")
        record = self._qc(config, path)
        assert record.basic_quality_status == "unreadable_input"

    def test_missing_file_is_missing_input(self, config, tmp_path):
        record = self._qc(config, tmp_path / "absent.fna")
        assert record.basic_quality_status == "missing_input"

    def test_good_file_still_validates(self, config, tmp_path):
        from papipeline.io.tsv import write_fasta

        path = tmp_path / "good.fna"
        write_fasta(path, [("c1", "ACGT" * 500)])
        record = self._qc(config, path)
        assert record.basic_quality_status in ("pass", "flagged")
        assert record.assembly_size == 2000

    def test_cohort_run_continues_past_a_bad_file(self, config, tmp_path):
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample
        from papipeline.io.tsv import write_fasta
        from papipeline.stages import validation

        good = tmp_path / "good.fna"
        write_fasta(good, [("c1", "ACGT" * 500)])
        bad = tmp_path / "bad.fna"
        bad.write_bytes(b">c\nACGT\n\xbf\xbf")
        empty = tmp_path / "empty.fna"
        empty.write_bytes(b"")
        manifest = SampleManifest(
            [
                Sample("S_ok", str(good)),
                Sample("S_bad", str(bad)),
                Sample("S_empty", str(empty)),
            ]
        )
        records = validation.run(config, manifest, RunMode.TEST)
        assert len(records) == 3
        summary = validation.summarise(records)
        assert summary["n_unreadable"] == 2
        assert summary["n_usable"] == 1
