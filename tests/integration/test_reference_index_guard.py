"""A `.fai` index that lies about the reference, caught before anything uses it.

Found the hard way. `bgzip` + `samtools faidx` on the pinned PAO1 reference
produced an index reporting **64,400 bases for a 6,264,404-base file**. Every
tool downstream then read a reference two orders of magnitude too small, and
produced confident, entirely wrong output: minimap2 emitted 64 alignments for a
6.4 Mb genome, and bcftools called 5,485 "variants" in a 1 Mb window, all with
`QUAL=30.4183`.

Nothing errored. The numbers came out. That is the whole danger — D5 already
warns that a wrong reference is invisible downstream, and this is a way of
producing one that passes every other check.

The asymmetry this pins: indexing the **uncompressed** FASTA gives the correct
6,264,404, while indexing the **bgzipped** copy gives 64,400. So the correct
path is proven correct, the truncating path is proven caught, and the two are not
treated as interchangeable. Code that reaches for whichever is convenient must
call the check.

The check is on the *index*, not the FASTA: the source file is byte-identical in
both cases and verified by D5. It is the index that is wrong, so verifying the
file again would prove nothing.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.reference import verify_fai_index

REAL_REFERENCE = Path(
    "db/reference/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna"
)
EXPECTED_BASES = 6_264_404

pytestmark = pytest.mark.skipif(
    not REAL_REFERENCE.is_file(),
    reason="the pinned reference is not downloaded on this machine",
)


def _bases_in_fai(fai: Path) -> int:
    return sum(int(line.split("\t")[1]) for line in fai.read_text().splitlines() if line)


def _index_with(fasta: Path, out_dir: Path) -> Path:
    """Index `fasta` with samtools, bgzipped or not as the caller arranged it."""
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / fasta.name
    shutil.copy(fasta, target)
    subprocess.run(["samtools", "faidx", str(target)], check=True, capture_output=True)
    return Path(f"{target}.fai")


@pytest.fixture()
def uncompressed_index(tmp_path):
    return _index_with(REAL_REFERENCE, tmp_path / "plain")


@pytest.fixture()
def bgzipped_index(tmp_path):
    """Index through bgzip, which is the path that truncates."""
    out_dir = tmp_path / "gz"
    out_dir.mkdir(parents=True)
    target = out_dir / REAL_REFERENCE.name
    with REAL_REFERENCE.open("rb") as src, target.open("wb") as dst:
        subprocess.run(["bgzip", "-c"], stdin=src, stdout=dst, check=True)
    subprocess.run(["samtools", "faidx", str(target)], check=True, capture_output=True)
    return Path(f"{target}.fai")


class TestTheCorrectPathIsCorrect:
    def test_uncompressed_index_reports_every_base(self, uncompressed_index):
        """Proves the good path, so the guard is not rejecting a working setup."""
        assert _bases_in_fai(uncompressed_index) == EXPECTED_BASES

    def test_it_passes_the_guard(self, uncompressed_index):
        assert verify_fai_index(uncompressed_index, expected_bases=EXPECTED_BASES) == EXPECTED_BASES


class TestTheTruncatingPathIsCaught:
    def test_bgzipped_index_really_does_truncate(self, bgzipped_index):
        """Documents the hazard itself.

        If a future samtools fixes this, this fails and the guard below can be
        revisited. Until then it is the reason the check exists, and a test that
        merely asserted "the guard works" would pass even if the hazard vanished
        and the guard had become theatre.
        """
        assert _bases_in_fai(bgzipped_index) < EXPECTED_BASES

    def test_the_guard_refuses_a_short_index(self, bgzipped_index):
        with pytest.raises(DataContractError) as excinfo:
            verify_fai_index(bgzipped_index, expected_bases=EXPECTED_BASES)
        message = str(excinfo.value)
        assert "64" in message or str(_bases_in_fai(bgzipped_index)) in message

    def test_the_refusal_names_both_numbers(self, bgzipped_index):
        """A refusal the user cannot act on is just an obstacle."""
        with pytest.raises(DataContractError) as excinfo:
            verify_fai_index(bgzipped_index, expected_bases=EXPECTED_BASES)
        message = str(excinfo.value)
        assert str(EXPECTED_BASES) in message
        assert str(_bases_in_fai(bgzipped_index)) in message

    def test_it_explains_that_bgzip_indexing_is_the_cause(self, bgzipped_index):
        with pytest.raises(DataContractError) as excinfo:
            verify_fai_index(bgzipped_index, expected_bases=EXPECTED_BASES)
        assert "bgzip" in str(excinfo.value).lower() or "compress" in str(excinfo.value).lower()


    def test_a_multi_contig_index_is_summed_not_maxed(self, tmp_path):
        """Several contigs must be totalled, not reduced to the largest.

        PAO1 is a single contig, so a `max` over contig lengths would give the
        right answer here and silently under-report every multi-contig reference
        that follows. The fixture has two contigs whose sum differs from the
        larger of them.
        """
        fai = tmp_path / "multi.fai"
        fai.write_text(
            "contig_a\t100\t0\t0\t0\ncontig_b\t300\t0\t0\t0\n"
        )
        assert verify_fai_index(fai, expected_bases=400) == 400
        with pytest.raises(DataContractError):
            verify_fai_index(fai, expected_bases=300)

class TestRefusals:
    def test_a_missing_index_is_refused(self, tmp_path):
        with pytest.raises(DataContractError):
            verify_fai_index(tmp_path / "absent.fai", expected_bases=EXPECTED_BASES)

    def test_an_empty_index_is_refused(self, tmp_path):
        empty = tmp_path / "empty.fai"
        empty.write_text("")
        with pytest.raises(DataContractError):
            verify_fai_index(empty, expected_bases=EXPECTED_BASES)

    def test_a_zero_length_contig_is_refused(self, tmp_path):
        """A zero-length entry is the shape a truncated index takes."""
        fai = tmp_path / "zero.fai"
        fai.write_text("NC_002516.2\t0\t0\t0\t0\n")
        with pytest.raises(DataContractError):
            verify_fai_index(fai, expected_bases=EXPECTED_BASES)

    def test_an_index_of_only_zeros_is_refused(self, tmp_path):
        """Every contig reporting zero is refused, not treated as 'no contigs'."""
        fai = tmp_path / "allzero.fai"
        fai.write_text("a\t0\t0\t0\t0\nb\t0\t0\t0\t0\n")
        with pytest.raises(DataContractError):
            verify_fai_index(fai, expected_bases=0)
