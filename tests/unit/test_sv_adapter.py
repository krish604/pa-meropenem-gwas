"""The assembly-based SV caller.

Every test here is for behaviour that would otherwise be a plausible number
rather than a failure: calling assembly fragmentation as variation, reporting a
lopsided-vs-two-sided gap correctly, or emitting `confirmed` from evidence that
cannot support it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import sv as caller


def _block(
    ref_start, ref_end, query_start, query_end, identity=99.0,
    query_name="ctg1", ref_name="NC_002516.2",
):
    return caller.AlignmentBlock(
        ref_start=ref_start, ref_end=ref_end,
        query_start=query_start, query_end=query_end,
        identity_pct=identity, ref_name=ref_name, query_name=query_name,
    )


def _row(block, ref_len=6264404, query_len=1228470):
    """A `show-coords -T -l -H` row, 11 columns, as nucmer 4.0.1 emits it."""
    return "\t".join(
        str(x)
        for x in (
            block.ref_start, block.ref_end,
            block.query_start, block.query_end,
            block.query_length, block.query_length,
            block.identity_pct, ref_len, query_len,
            block.ref_name, block.query_name,
        )
    )


class TestParsing:
    def test_it_reads_the_verified_column_layout(self):
        """Columns verified against real output, not assumed.

        1 S1 2 E1 3 S2 4 E2 5 LEN1 6 LEN2 7 %IDY 8 LEN R 9 LEN Q 10 ref 11 query
        """
        blocks = caller.parse_show_coords(_row(_block(100, 200, 10, 110)))
        assert len(blocks) == 1
        b = blocks[0]
        assert (b.ref_start, b.ref_end) == (100, 200)
        assert (b.query_start, b.query_end) == (10, 110)
        assert b.identity_pct == 99.0
        assert b.query_name == "ctg1"

    def test_a_short_row_is_dropped_not_fatal(self):
        assert caller.parse_show_coords("1\t2\t3\n") == []

    def test_an_unparsable_number_is_dropped(self):
        assert caller.parse_show_coords(
            _row(_block(1, 2, 3, 4)).replace("\t99.0\t", "\tNA\t")
        ) == []


class TestCommandConstruction:
    def test_every_nucmer_flag_is_explicit(self, tmp_path):
        """A default that changes between mummer releases would silently change
        which variants are found, and a pin in the overlay must mean something."""
        cmd = caller.nucmer_command(
            program="nucmer", reference=tmp_path / "ref.fna",
            query=tmp_path / "q.fna", prefix=tmp_path / "p",
            minmatch=20, mincluster=65, maxgap=90, maxmatch=True, threads=4,
        )
        for flag, value in (
            ("--minmatch", "20"), ("--mincluster", "65"),
            ("--maxgap", "90"), ("--threads", "4"),
        ):
            assert cmd[cmd.index(flag) + 1] == value
        assert "--maxmatch" in cmd
        assert cmd[-2:] == [str(tmp_path / "ref.fna"), str(tmp_path / "q.fna")]

    def test_maxmatch_is_omitted_when_false(self, tmp_path):
        cmd = caller.nucmer_command(
            program="nucmer", reference=tmp_path / "r", query=tmp_path / "q",
            prefix=tmp_path / "p", minmatch=20, mincluster=65, maxgap=90,
            maxmatch=False, threads=1,
        )
        assert "--maxmatch" not in cmd

    def test_show_coords_sorts_by_query(self, tmp_path):
        """`-q` is load-bearing: gaps are found within one contig, so blocks
        must arrive grouped and ordered by query coordinate. The `-r` default
        interleaves contigs and the gap logic then compares unrelated blocks."""
        cmd = caller.show_coords_command(
            program="show-coords", delta=tmp_path / "p.delta"
        )
        assert "-q" in cmd and "-T" in cmd and "-l" in cmd and "-H" in cmd


class TestGapCalling:
    def test_an_insertion_in_the_assembly_is_found(self):
        """Query coordinate advances, reference does not."""
        blocks = [_block(1000, 2000, 1000, 2000), _block(2001, 3000, 5000, 5999)]
        out = caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        )
        assert len(out) == 1
        assert out[0].variant_type == "insertion"
        assert out[0].size == 5000 - 2000 - 1

    def test_a_deletion_in_the_assembly_is_found(self):
        """Reference coordinate advances, query does not."""
        blocks = [_block(1000, 2000, 1000, 2000), _block(9000, 10000, 2001, 3000)]
        out = caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        )
        assert len(out) == 1
        assert out[0].variant_type == "deletion"
        assert out[0].size == 9000 - 2000 - 1

    def test_a_gap_large_in_both_coordinates_is_not_a_variant(self):
        """A repeat or a contig join advances both.

        Reporting it would be reporting the aligner's uncertainty as biology.
        """
        blocks = [_block(1000, 2000, 1000, 2000), _block(9000, 10000, 9000, 10000)]
        assert caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        ) == []

    def test_ordinary_divergence_is_not_a_variant(self):
        blocks = [_block(1000, 2000, 1000, 2000), _block(2005, 3000, 2004, 2999)]
        assert caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        ) == []

    def test_a_gap_below_the_threshold_is_not_reported(self):
        blocks = [_block(1000, 2000, 1000, 2000), _block(2001, 2500, 2500, 2999)]
        assert caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        ) == []

    def test_every_call_is_a_candidate(self):
        """Assembly-only evidence cannot meet rule 7's bar for `confirmed`.

        `confirmed_only()` is empty for every REAL cohort by construction, and
        that is a property of the input rather than something to tune away.
        """
        blocks = [_block(1000, 2000, 1000, 2000), _block(2001, 3000, 5000, 5999)]
        out = caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        )
        assert out and all(c.call_status == "candidate" for c in out)
        assert not any(c.call_status == "confirmed" for c in out)

    def test_a_noisy_flanking_block_suppresses_the_call(self):
        """A breakpoint is only credible between two confident blocks."""
        blocks = [
            _block(1000, 2000, 1000, 2000, identity=99.0),
            _block(2001, 3000, 5000, 5999, identity=70.0),
        ]
        assert caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        ) == []

    def test_gaps_are_not_compared_across_contigs(self):
        """The last block of one contig and the first of the next are not a
        deletion; they are two contigs."""
        blocks = [
            _block(1000, 2000, 1000, 2000, query_name="ctg1"),
            _block(9000, 10000, 2001, 3000, query_name="ctg2"),
        ]
        assert caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        ) == []

    def test_variant_ids_are_unique_within_a_sample(self):
        blocks = [
            _block(1000, 2000, 1000, 2000, query_name="c1"),
            _block(2001, 3000, 5000, 5999, query_name="c1"),
            _block(1000, 2000, 1000, 2000, query_name="c2"),
            _block(2001, 3000, 5000, 5999, query_name="c2"),
        ]
        out = caller.call_from_gaps(
            blocks, sample_id="S1", min_sv_size=1000,
            min_flanking_identity_pct=90.0,
        )
        assert len({c.variant_id for c in out}) == len(out) == 2


class TestUnassessableRegions:
    def test_contig_ends_are_unassessable_not_insertions(self, tmp_path):
        """The fragmentation trap.

        A draft assembly breaks where a read spans a repeat, so unaligned
        sequence at a contig end is routine and carries no information about the
        reference. Rule 7 has a state for exactly this - a contig gap - and
        conflating it with "no variant present" is what it forbids.
        """
        blocks = [_block(1000, 2000, 5000, 6000)]
        regions = caller.unassessable_regions(
            blocks,
            query_contig_lengths={"ctg1": 20000},
            min_sv_size=1000,
        )
        # leading 1..4999 and trailing 6001..20000
        assert len(regions) == 2
        assert regions[0] == ("ctg1", 1, 4999)

    def test_a_fully_unaligned_contig_is_one_unassessable_region(self):
        regions = caller.unassessable_regions(
            [], query_contig_lengths={"ctg1": 5000}, min_sv_size=1000
        )
        assert regions == [("ctg1", 0, 5000)]

    def test_a_fully_aligned_contig_has_no_unassessable_region(self):
        blocks = [_block(1, 20000, 1, 20000)]
        assert caller.unassessable_regions(
            blocks, query_contig_lengths={"ctg1": 20000}, min_sv_size=1000
        ) == []


class TestContigLengths:
    def test_it_sums_wrapped_sequence_lines(self, tmp_path):
        fasta = tmp_path / "g.fna"
        fasta.write_text(">c1 first\nACGTAC\nGTACGT\n>c2\nAAAA\n", encoding="utf-8")
        assert caller.contig_lengths(fasta) == {"c1": 12, "c2": 4}


class TestScreenAssembly:
    def _runner(self, results):
        calls = []

        def runner(command):
            calls.append(command)
            return results[len(calls) - 1]

        return runner, calls

    def _result(self, stdout="", returncode=0):
        from papipeline.adapters.external import CommandResult

        return CommandResult(
            command=["x"], returncode=returncode, stdout=stdout,
            stderr="", duration_s=0.1,
        )

    def test_a_missing_delta_is_reported_not_ignored(self, tmp_path):
        """nucmer exiting 0 with no delta means the screen found nothing to
        call; treating that as "no SVs" would be a fabricated negative."""
        runner, _ = self._runner([self._result()])
        with pytest.raises(caller.SvPreflightError):
            caller.screen_assembly(
                sample_id="S1", genome=tmp_path / "g.fna",
                reference=tmp_path / "ref.fna", workdir=tmp_path / "w",
                nucmer="nucmer", show_coords="show-coords",
                minmatch=20, mincluster=65, maxgap=90, maxmatch=True,
                threads=1, min_sv_size=1000, min_flanking_identity_pct=90.0,
                runner=runner,
            )

    def test_a_nonzero_show_coords_is_raised(self, tmp_path):
        delta = tmp_path / "w" / "S1.delta"
        delta.parent.mkdir(parents=True)
        delta.write_text("x", encoding="utf-8")
        runner, _ = self._runner([self._result(), self._result(returncode=1)])
        from papipeline.errors import ToolExecutionError

        with pytest.raises(ToolExecutionError):
            caller.screen_assembly(
                sample_id="S1", genome=tmp_path / "g.fna",
                reference=tmp_path / "ref.fna", workdir=tmp_path / "w",
                nucmer="nucmer", show_coords="show-coords",
                minmatch=20, mincluster=65, maxgap=90, maxmatch=True,
                threads=1, min_sv_size=1000, min_flanking_identity_pct=90.0,
                runner=runner,
            )