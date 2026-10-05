"""The oprD locus resolver: what it resolves, and what it refuses.

The cases here are the ones that decided the design. Two of them are refusals
that a naive best-hit rule would get wrong, and both are asserted against the
real measurements taken from the ten smoke isolates rather than invented:

* a two-fragment tie is **not** a tie. Five isolates have PA0958 split across two
  abutting CDS on one contig; their per-CDS coverages are 38-53% but the unions
  are 67-100%. Refusing them would report five isolates that carry the locus as
  unlocatable, so the fragments must merge.
* a refusal is **not** ``absent``. This is the property the whole module exists
  to protect, and it is asserted directly.

Real-data regression values are recorded per test with the command that produced
them, per the repo's convention of pinning a measurement rather than a belief.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import oprd_locus as oprd
from papipeline.adapters.oprd_locus import Verdict
from papipeline.errors import PipelineError, ToolExecutionError

QUERY_LENGTH = 443  # PAO1 PA0958, 443 residues + stop. See the module docstring.


def _hit(
    subject_id,
    *,
    identity=97.0,
    qstart=1,
    qend=None,
    contig="contig_1",
    strand="+",
    qlen=QUERY_LENGTH,
    bitscore=400.0,
):
    """One blastp HSP against a 443-residue reference."""
    qend = qend if qend is not None else qstart
    return oprd.Hit(
        subject_id=subject_id,
        contig=contig,
        strand=strand,
        identity_pct=identity,
        alignment_length=qend - qstart + 1,
        query_length=qlen,
        subject_length=qend - qstart + 1,
        query_start=qstart,
        query_end=qend,
        bitscore=bitscore,
    )


def _paralog(
    name, identity=36.0, contig="contig_9", strand="+",
    qstart=1, qend=QUERY_LENGTH,
):
    """The family member that dominates the search: ~36% identity, full length."""
    return _hit(
        name, identity=identity, qstart=qstart, qend=qend,
        contig=contig, strand=strand, bitscore=230.0,
    )


class TestAClearBestHitResolves:
    def test_a_single_full_length_hit_resolves(self):
        result = oprd.resolve("S1", [_hit("GENE_1", qstart=1, qend=443)])
        assert result.verdict == Verdict.RESOLVED
        assert result.is_resolved and not result.is_refused
        assert result.gene_ids == ("GENE_1",)
        assert result.contig == "contig_1"
        assert result.coverage_pct == 100.0
        assert result.identity_pct == 97.0

    def test_a_paralog_alongside_the_true_locus_does_not_contend(self):
        """PDT000167136.1 measured: DCAFIL_02540 at 97.97% plus four paralogs.

        The paralogs sit at 33.7-36.6% over 76-108% of the query, so a naive
        "highest identity wins" is fine here but a naive "more hits is better"
        would pick the paralog family. Neither is allowed to unseat the orthologue.
        """
        result = oprd.resolve(
            "PDT000167136.1",
            [
                _hit("DCAFIL_02540", identity=97.97, qstart=1, qend=443),
                _paralog("DCAFIL_00733", identity=34.55),
                _paralog("DCAFIL_02412", identity=36.58),
                _paralog("DCAFIL_03519", identity=36.28),
                _paralog("DCAFIL_04775", identity=36.00),
            ],
        )
        assert result.verdict == Verdict.RESOLVED
        assert result.gene_ids == ("DCAFIL_02540",)
        assert result.n_candidate_hits == 1
        assert result.n_paralog_hits == 4
        assert result.n_total_hits == 5

    def test_near_but_not_quite_full_coverage_still_resolves(self):
        """PDT000292998.1 measured: 95.91% identity over 416/443 = 93.91%."""
        result = oprd.resolve(
            "PDT000292998.1", [_hit("LAGMFP_03568", identity=95.91, qstart=1, qend=416)]
        )
        assert result.verdict == Verdict.RESOLVED
        assert result.coverage_pct == pytest.approx(93.91, abs=0.01)


class TestASplitLocusIsNotATie:
    """Five of the ten isolates. The rule the brief proposed would fail all five."""

    def test_two_abutting_fragments_merge_into_one_locus(self):
        """PDT000167133.1 measured: MIKKFE_01919 q1-236 at 97.46%,
        MIKKFE_01920 q238-443 at 98.06%; union 442/443 = 99.77%.

        Per-CDS each fragment is 53% and 46% - both under the 80% floor. Only the
        union clears it. Per-CDS coverage calls this isolate absent.
        """
        result = oprd.resolve(
            "PDT000167133.1",
            [
                _hit("MIKKFE_01919", identity=97.46, qstart=1, qend=236),
                _hit("MIKKFE_01920", identity=98.06, qstart=238, qend=443),
                _paralog("MIKKFE_00212", identity=34.55),
            ],
        )
        assert result.verdict == Verdict.RESOLVED
        assert result.gene_ids == ("MIKKFE_01919", "MIKKFE_01920")
        assert result.coverage_pct == pytest.approx(99.77, abs=0.01)

    def test_fragments_with_no_gap_at_all_merge(self):
        """PDT000167135.1 measured: q1-212 and q210-443, union 443/443 = 100%."""
        result = oprd.resolve(
            "PDT000167135.1",
            [
                _hit("JHMFCP_04448", identity=98.29, qstart=210, qend=443),
                _hit("JHMFCP_04449", identity=97.17, qstart=1, qend=212),
            ],
        )
        assert result.verdict == Verdict.RESOLVED
        assert result.coverage_pct == 100.0

    def test_the_split_is_not_how_two_genes_are_told_apart(self):
        """Same query spans, different contigs: that IS two loci, and it refuses.

        The fragments merge on contig+strand, never on query overlap alone, so a
        pair of separate genes cannot be silently concatenated into one call.
        """
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=236, contig="contig_1"),
                _hit("GENE_B", qstart=238, qend=443, contig="contig_7"),
            ],
        )
        # Neither cluster clears the floor on its own (53.2% and 46.4%), so the
        # verdict is INSUFFICIENT rather than AMBIGUOUS. It is still a refusal,
        # and crucially it is not a resolution built by silently joining the two.
        assert result.verdict == Verdict.INSUFFICIENT
        assert result.is_refused
        assert result.gene_ids == ("GENE_A",)

    def test_opposite_strands_do_not_merge(self):
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=236, contig="contig_1", strand="+"),
                _hit("GENE_B", qstart=238, qend=443, contig="contig_1", strand="-"),
            ],
        )
        assert result.is_refused
        assert result.verdict != Verdict.RESOLVED
        assert len(result.clusters) == 2


class TestNoHitRefusesWithANamedReason:
    def test_paralogs_only_is_a_refusal_not_an_absence(self):
        """The case the whole module exists for.

        PDT000424983.1's proteome is dominated by the porin family at 33.7-42.6%.
        Nothing here supports "the locus is present" and nothing supports "it is
        absent", so the verdict is a named refusal and never ``absent``.
        """
        result = oprd.resolve(
            "PDT000424983.1",
            [
                _paralog("NCHPHM_01114", identity=35.78, contig="contig_26"),
                _paralog("NCHPHM_05423", identity=36.75, contig="contig_108"),
                _paralog("NCHPHM_05564", identity=33.68, contig="contig_114"),
            ],
        )
        assert result.verdict == Verdict.NO_HIT
        assert result.is_refused
        assert "NOT evidence that oprD is absent" in result.reason
        assert "42.64" not in result.reason  # the max seen was 36.75, not invented
        assert "36.75" in result.reason

    def test_an_empty_search_refuses_rather_than_resolving(self):
        result = oprd.resolve("S1", [])
        assert result.verdict == Verdict.NO_HIT
        assert result.n_candidate_hits == 0
        assert result.coverage_pct == 0.0

    def test_the_floor_sits_inside_the_measured_gap(self):
        """45.16% is the highest non-orthologue hit and 88.61% the lowest
        orthologue fragment, measured over all ten isolates. The floor must stay
        strictly between them or one of the two tiers changes meaning.
        """
        assert 45.16 < oprd.MIN_IDENTITY_PCT < 88.61

    def test_an_identical_paralog_does_not_rescue_a_refusal(self):
        result = oprd.resolve(
            "S1", [_paralog("P", identity=69.99, qstart=1, qend=443)]
        )
        assert result.verdict == Verdict.NO_HIT


class TestATieRefuses:
    def test_two_full_length_loci_on_different_contigs_refuse(self):
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=443, contig="contig_1"),
                _hit("GENE_B", qstart=1, qend=443, contig="contig_2"),
            ],
        )
        assert result.verdict == Verdict.AMBIGUOUS
        assert "cannot be told apart" in result.reason
        assert result.gene_ids == ()

    def test_a_tie_never_reports_either_gene_as_oprd(self):
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=443, contig="contig_1"),
                _hit("GENE_B", qstart=1, qend=443, contig="contig_2"),
            ],
        )
        assert "GENE_A" not in result.gene_ids
        assert "GENE_B" not in result.gene_ids

    def test_a_second_partial_copy_near_half_the_floor_refuses(self):
        """PDT000294804.1's shape, escalated: a 45% second copy at 88% identity
        is a paralog problem, and calling oprD from it would be a guess."""
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=443, contig="contig_1"),
                _hit(
                    "GENE_B", identity=88.0, qstart=1, qend=200, contig="contig_5"
                ),
            ],
        )
        assert result.verdict == Verdict.AMBIGUOUS_SECOND
        assert result.is_refused
        assert "paralog assignment problem" in result.reason

    def test_a_second_copy_below_half_the_floor_does_not_refuse(self):
        result = oprd.resolve(
            "S1",
            [
                _hit("GENE_A", qstart=1, qend=443, contig="contig_1"),
                _hit("GENE_B", identity=88.0, qstart=1, qend=100, contig="contig_5"),
            ],
        )
        assert result.verdict == Verdict.RESOLVED
        assert result.gene_ids == ("GENE_A",)


class TestInsufficientCoverageIsItsOwnRefusal:
    """PDT000294804.1 measured: union 72.91%. Present signal, under-determined."""

    def test_a_fragment_below_the_floor_is_not_present_and_not_absent(self):
        result = oprd.resolve(
            "PDT000294804.1",
            [
                _hit("CFOAGI_04608", identity=98.25, qstart=1, qend=171),
                _hit("CFOAGI_04606", identity=97.37, qstart=292, qend=443),
            ],
        )
        assert result.verdict == Verdict.INSUFFICIENT
        assert result.is_refused
        assert result.coverage_pct == pytest.approx(72.91, abs=0.01)
        assert "neither shown present nor shown absent" in result.reason

    def test_a_single_short_fragment_is_the_weakest_call_in_the_cohort(self):
        """PDT000294805.1 measured: NGELIL_00073, q207-443, 88.61% over 53.5%.

        This is the call S flagged as least certain, and it is right to refuse it.
        88.61% is 11 points below the rest of the cohort; at 53.5% coverage there
        is not enough of the protein to say whether this is a divergent OprD or a
        drifted paralog.
        """
        result = oprd.resolve(
            "PDT000294805.1", [_hit("NGELIL_00073", identity=88.61, qstart=207, qend=443)]
        )
        assert result.verdict == Verdict.INSUFFICIENT
        assert result.coverage_pct == pytest.approx(53.50, abs=0.01)

    def test_the_coverage_floor_sits_inside_the_measured_gap(self):
        """Measured unions: 72.91% then 93.91%, nothing between. See module."""
        assert 72.91 < oprd.MIN_COVERAGE_PCT < 93.91


class TestClustering:
    def test_overlapping_fragments_are_counted_once(self):
        cluster = oprd.cluster_hits(
            [_hit("A", qstart=1, qend=300), _hit("B", qstart=200, qend=443)]
        )[0]
        assert cluster.coverage_pct == 100.0
        assert len(cluster.covered) == QUERY_LENGTH

    def test_clusters_are_ordered_by_coverage(self):
        clusters = oprd.cluster_hits(
            [
                _hit("BIG", qstart=1, qend=443, contig="contig_1"),
                _hit("SMALL", qstart=1, qend=50, contig="contig_2"),
            ]
        )
        assert [c.contig for c in clusters] == ["contig_1", "contig_2"]

    def test_a_hit_below_the_floor_is_not_a_candidate_at_all(self):
        result = oprd.resolve(
            "S1", [_hit("P", identity=36.26, qstart=1, qend=443)]
        )
        assert result.n_candidate_hits == 0


class TestParsingAndLocations:
    def test_the_outfmt_requests_query_coordinates(self):
        """Without qstart/qend the fragments cannot be merged, so this is part of
        the contract rather than an optimisation."""
        assert "qstart" in oprd.OUTFMT_FIELDS and "qend" in oprd.OUTFMT_FIELDS

    def test_a_row_of_the_wrong_width_is_refused_not_parsed(self):
        with pytest.raises(PipelineError) as excinfo:
            oprd.parse_blast_table("GENE_1\t97.0\t443\n", QUERY_LENGTH)
        assert "outfmt" in str(excinfo.value)

    def test_a_non_numeric_field_is_refused(self):
        with pytest.raises(PipelineError) as excinfo:
            oprd.parse_blast_table(
                "\t".join(
                    ["GENE_1", "not-a-number", "443", "443", "441", "1", "443",
                     "1", "441", "887"]
                ) + "\n",
                QUERY_LENGTH,
            )
        assert "not numeric" in str(excinfo.value)

    def test_a_well_formed_row_parses(self):
        hits = oprd.parse_blast_table(
            "\t".join(
                ["DCAFIL_02540", "97.97", "443", "443", "443", "1", "443", "1",
                 "443", "887"]
            ) + "\n",
            QUERY_LENGTH,
        )
        assert len(hits) == 1
        assert hits[0].subject_id == "DCAFIL_02540"
        assert hits[0].identity_pct == 97.97
        assert hits[0].query_coverage_pct == 100.0

    def test_cds_locations_come_from_the_isolate_gff(self, tmp_path):
        gff = tmp_path / "iso.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "contig_1\tBakta\tCDS\t100\t500\t.\t+\t0\tID=cds-A;locus_tag=GENE_A\n"
            "contig_7\tBakta\tCDS\t200\t600\t.\t-\t0\tID=cds-B;locus_tag=GENE_B\n",
            encoding="utf-8",
        )
        found = oprd.load_cds_locations(gff)
        assert found == {"GENE_A": ("contig_1", "+"), "GENE_B": ("contig_7", "-")}

    def test_a_missing_gff_is_named(self, tmp_path):
        with pytest.raises(PipelineError) as excinfo:
            oprd.load_cds_locations(tmp_path / "absent.gff3")
        assert "Isolate GFF3 not found" in str(excinfo.value)

    def test_locations_turn_two_fragments_into_one_cluster(self):
        """The end-to-end mechanism the tie test depends on."""
        hits = [
            _hit("GENE_A", qstart=1, qend=236, contig=""),
            _hit("GENE_B", qstart=238, qend=443, contig=""),
        ]
        located = oprd.attach_locations(
            hits, {"GENE_A": ("contig_1", "+"), "GENE_B": ("contig_1", "+")}
        )
        result = oprd.resolve("S1", located)
        assert result.verdict == Verdict.RESOLVED
        assert len(oprd.cluster_hits(located)) == 1


class TestBlastpInvocation:
    def test_the_command_is_blastp_with_the_declared_outfmt(self):
        command = oprd.blastp_command(
            program="blastp", query=Path("q.faa"), subject=Path("s.faa"), threads=2
        )
        assert command[0] == "blastp"
        assert "-query" in command and "-subject" in command
        assert command[command.index("-outfmt") + 1] == oprd.OUTFMT
        assert command[command.index("-num_threads") + 1] == "2"

    def test_max_target_seqs_is_a_positive_integer(self):
        """blastp 2.17.0 rejects 0 with `expected >=1`. Verified, not assumed."""
        assert oprd.MAX_TARGET_SEQS >= 1
        command = oprd.blastp_command(
            program="blastp", query=Path("q.faa"), subject=Path("s.faa"), threads=1
        )
        assert int(command[command.index("-max_target_seqs") + 1]) >= 1

    def test_a_non_zero_exit_raises_with_the_stderr(self, tmp_path):
        (tmp_path / "q.faa").write_text(">x\nMKV\n", encoding="utf-8")

        def runner(command):
            return 1, "", "Error: something specific"

        with pytest.raises(ToolExecutionError) as excinfo:
            oprd.run_blastp(
                tmp_path / "q.faa", tmp_path / "s.faa", runner=runner
            )
        assert "something specific" in str(excinfo.value)

    def test_a_clean_run_parses_its_own_output(self, tmp_path):
        (tmp_path / "q.faa").write_text(">x\nMKV\n", encoding="utf-8")
        row = "\t".join(
            ["GENE_A", "97.0", "3", "3", "3", "1", "3", "1", "3", "40"]
        )

        def runner(command):
            return 0, row + "\n", ""

        hits = oprd.run_blastp(
            tmp_path / "q.faa", tmp_path / "s.faa", runner=runner
        )
        assert len(hits) == 1
        assert hits[0].query_length == 3


class TestReferenceProteinDerivation:
    """The query is derived from the two pinned files, not downloaded separately."""

    GFF = (
        "##gff-version 3\n"
        "NC_002516.2\tRefSeq\tgene\t1043983\t1045314\t.\t-\t.\t"
        "ID=gene-PA0958;Name=oprD;gene=oprD;locus_tag=PA0958\n"
        "NC_002516.2\tRefSeq\tCDS\t1043983\t1045314\t.\t-\t0\t"
        "ID=cds-NP_249649.1;Parent=gene-PA0958;Name=NP_249649.1;gbkey=CDS;"
        "gene=oprD;locus_tag=PA0958;product=porin D;protein_id=NP_249649.1;"
        "transl_table=11\n"
    )

    def _fasta(self, tmp_path, contig="NC_002516.2"):
        """A FASTA whose NC_002516.2 record spans the PA0958 CDS.

        The bases are placeholders - this asserts the pin's *structure* is read
        correctly, not PAO1's sequence, which the real-data blast runs cover.
        The record is written whole rather than in slices so no read offset can
        be got wrong.

        PA0958 is on the **minus** strand, so the plus-strand bases are the
        reverse complement of the coding sequence. The coding sequence here is
        `"ATG" + "GCT"*442 + "TAA"` - 444 codons ending in a stop, exactly PA0958's
        shape - whose reverse complement is `"TTA" + "AGC"*442 + "CAT"`. Writing
        the coding sequence in directly would yield a protein starting with S, and
        the module would - correctly - refuse it.
        """
        path = tmp_path / "ref.fna"
        record = list("N" * (1043982))              # plus-strand 1..1043982
        record.append("TTA" + "AGC" * 442 + "CAT")  # 1332 bases, minus-strand CDS
        wrapped = "\n".join(
            "".join(record[i:i + 60]) for i in range(0, len(record), 60)
        )
        path.write_text(f">{contig} ASM676v1\n{wrapped}\n", encoding="utf-8")
        return path

    def test_a_minus_strand_locus_is_translated_in_its_reading_frame(self, tmp_path):
        """If the strand were ignored this protein would start with S, not M."""
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        protein = oprd.reference_protein(gff, self._fasta(tmp_path), "PA0958")
        assert protein.startswith("MA")

    def test_the_terminal_stop_is_dropped_so_a_perfect_match_is_100_percent(self, tmp_path):
        """1332 bp is 444 codons: 443 residues and a stop.

        Keeping the stop would cap coverage at 443/444 = 99.77% and quietly move
        the coverage floor, so the query is 443 residues.
        """
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        protein = oprd.reference_protein(gff, self._fasta(tmp_path), "PA0958")
        assert len(protein) == 443
        assert "*" not in protein

    def test_an_internal_stop_is_refused(self, tmp_path):
        """A pseudogene or a frameshifted pin makes every identity number
        meaningless, so it is refused rather than translated."""
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        path = self._fasta(tmp_path)
        # The gene is on the minus strand, so its codon 1 is the LAST plus-strand
        # codon: seq[1043982:1043985] == "TTA" is the reverse complement of the
        # terminal TAA. Gene codon 100 therefore starts 297 bases further along
        # the plus strand, at 0-based 1045014; making that "TTA" puts a TAA stop
        # in the middle of the reading frame rather than at its end.
        lines = path.read_text(encoding="utf-8").splitlines()
        seq = "".join(l for l in lines if not l.startswith(">"))
        assert seq[1043982:1043985] == "TTA", "terminal stop codon, revcomp"
        assert seq[1045014:1045017] == "AGC", "gene codon 100"
        seq = seq[:1045014] + "TTA" + seq[1045017:]
        path.write_text(
            ">NC_002516.2 ASM676v1\n"
            + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
            + "\n",
            encoding="utf-8",
        )
        with pytest.raises(PipelineError) as excinfo:
            oprd.reference_protein(gff, path, "PA0958")
        assert "internal stop" in str(excinfo.value)

    def test_the_structure_is_derived_from_the_pin(self, tmp_path):
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        protein = oprd.reference_protein(gff, self._fasta(tmp_path), "PA0958")
        assert len(protein) == 443  # 1332 bp / 3 = 444 codons, stop removed
        assert protein.startswith("M")

    def test_a_missing_locus_is_named_not_guessed(self, tmp_path):
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        with pytest.raises(PipelineError) as excinfo:
            oprd.reference_protein(gff, self._fasta(tmp_path), "PA9999")
        assert "PA9999" in str(excinfo.value)

    def test_a_multi_segment_cds_is_refused(self, tmp_path):
        gff = tmp_path / "ref.gff"
        gff.write_text(
            self.GFF
            + "NC_002516.2\tRefSeq\tCDS\t1043983\t1044783\t.\t-\t0\t"
            "ID=cds-2;Parent=gene-PA0958;locus_tag=PA0958;transl_table=11\n",
            encoding="utf-8",
        )
        with pytest.raises(PipelineError) as excinfo:
            oprd.reference_protein(gff, self._fasta(tmp_path), "PA0958")
        assert "CDS features" in str(excinfo.value)

    def test_an_unimplemented_transl_table_is_refused(self, tmp_path):
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF.replace("transl_table=11", "transl_table=4"), "utf-8")
        with pytest.raises(PipelineError) as excinfo:
            oprd.reference_protein(gff, self._fasta(tmp_path), "PA0958")
        assert "transl_table" in str(excinfo.value)

    def test_a_contig_absent_from_the_fasta_is_refused(self, tmp_path):
        gff = tmp_path / "ref.gff"
        gff.write_text(self.GFF, encoding="utf-8")
        with pytest.raises(PipelineError) as excinfo:
            oprd.reference_protein(
                gff, self._fasta(tmp_path, contig="SOMETHING_ELSE"), "PA0958"
            )
        assert "not found" in str(excinfo.value)