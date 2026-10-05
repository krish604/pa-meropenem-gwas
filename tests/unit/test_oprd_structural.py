"""Structural classification of the oprD locus: verdicts, and what they rest on.

The module this tests answers two different questions. ``resolve`` asks *is this
the orthologue locus*; :func:`classify_structure` asks *is the gene WHOLE*.
Only the second one can say ``disrupted``, and only a recorded negative search
can say ``absent`` - which is the whole reason the classifier takes a
:class:`LocusSearch` rather than a bare hit count.

Real-data regression values come from the ten smoke isolates, measured with
``tblastn`` of the PAO1 PA0958 protein (444 aa) against each assembly:

    tblastn -query PA0958_protein -db <isolate_index> \
      -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart \
              send evalue bitscore qseq sseq" \
      -evalue 1e-3 -max_target_seqs 5000 -db_gencode 11 -num_threads 4

    -db_gencode 11 is bacterial code; blast's default is 1 (standard), which is
    wrong for P. aeruginosa. -max_target_seqs must be >= 1 in blast 2.17.

Every isolate returned a hit at 94.59-98.68% identity and 99.32-100% query
coverage, so **none of the ten is ``absent``**, and the classification is driven
entirely by the reading frame. Per-isolate ORF lengths below are the assembly's
own translation from the OprD start codon to its first in-frame stop.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from papipeline import viz
from papipeline.adapters import oprd_locus as oprd
from papipeline.adapters.oprd_locus import (
    LocusEvidence,
    LocusSearch,
    StructuralCall,
    StructuralVerdict,
    classify_structure,
    sensitivity_table,
)

#: PAO1 PA0958's **codon** count, i.e. 443 residues plus the terminal stop.
#:
#: The measurement this module records puts the reference at 444 including its
#: own '*', so 444 is what every caller passes as ``reference_aa_length``. Note
#: the previous revision's comment here said "444 residues" and that was wrong
#: by one - ``reference_protein()`` in the adapter already documented 443.
REFERENCE_AA = 444
REFERENCE_RESIDUES = REFERENCE_AA - 1


def _search(n_hits: int = 2, *, recorded: bool = True) -> LocusSearch:
    return LocusSearch(
        program="tblastn",
        query="PAO1 PA0958 oprD protein (444 aa)",
        database="db/smoke_genomes/PDT000167136.1.fna (BLAST nucleotide index)",
        parameters=(("evalue", "1e-3"), ("max_target_seqs", "5000"),
                    ("db_gencode", "11")),
        command=("tblastn -query PA0958.faa -db index "
                 "-evalue 1e-3 -max_target_seqs 5000 -db_gencode 11"),
        n_hits=n_hits,
    )


@dataclass(frozen=True)
class _Measured:
    """One isolate as tblastn and the assembly's own translation measured it."""

    sample_id: str
    orf_aa: int
    identity_pct: float
    coverage_pct: float
    stop_codon: str | None
    stop_position_aa: int | None
    split: bool
    contig: str
    strand: str
    footprint: tuple
    #: True where the lesion is a FRAMESHIFT rather than a plain nonsense codon.
    #: Measured, not assumed. Re-measured from the assemblies in this round by
    #: rebuilding each isolate's oprD CDS from its own contig and aligning it to
    #: the PAO1 PA0958 CDS nucleotide by nucleotide; the previous revision had
    #: this set for only four of the eight frameshifted isolates and called the
    #: other four "nonsense". All eight carry an upstream indel whose net length
    #: is not a multiple of three, and in each of them the aligned reference
    #: codon at the stop is ITSELF the stop codon - zero nucleotides changed -
    #: so the stop is a consequence of the frame, not a substitution:
    #:
    #:     sample          net frame at the stop   codon   aligned ref codon
    #:     PDT000167133.1  -4                      TGA     TGA
    #:     PDT000167135.1  +1                      TAA     TAA
    #:     PDT000292995.1  -2                      TGA     TGA
    #:     PDT000292998.1  +5                      TGA     TGA
    #:     PDT000294804.1  +1361 (a 1361 nt insertion inside the locus)
    #:                                               TGA     -
    #:     PDT000294805.1  -2                      TGA     TGA
    #:     PDT000311294.1  -2                      TGA     TGA
    #:     PDT000424983.1  -11                     TGA     -
    #:
    #: Correcting the flag changes lesion_type from premature_stop to
    #: frameshift_with_premature_stop for those four and changes NO verdict.
    frameshift: bool = False
    #: Whether the recorded stop IS PAO1's own terminator, reached in frame.
    #: Measured by aligning the isolate's CDS to the reference CDS: True when
    #: the isolate's reading frame runs through the reference's terminal codon.
    #: False for every premature stop, by definition. This is the field RULING 2
    #: turns on, and `None` is the honest "not asked" - see
    #: `TestRuling2AnInFrameDeletionIsNotDisruption`.
    stop_is_reference_terminator: bool = False


# Measured on the ten smoke isolates. `orf_aa` counts residues from the OprD ATG
# to the first in-frame stop; PAO1 is 444 aa with its own stop at codon 444.
#
# `stop_is_reference_terminator` is True for exactly the two isolates whose
# reading frame reaches PAO1's own terminator: PDT000167136.1 with no indel at
# all, and PDT000034122.1 with a net -6 nt IN-FRAME deletion (DEL 5 at reference
# nt 1132-1136 plus DEL 1 at 1149, 17 nt apart, inside the 12-nt tandem repeat
# `CTACGG`x2 at reference nt 1143-1154) that still ends on the native TAA at
# isolate codon 442. Under RULING 2 the second of those is `intact`.
MEASURED = (
    _Measured("PDT000034122.1", 441, 94.59, 100.00, "TAA", 442, False,
              "JFJU01000003.1", "+", (510136, 511461),
              stop_is_reference_terminator=True),
    _Measured("PDT000167133.1", 236, 98.55, 99.77, "TGA", 237, True,
              "MPBS01000001.1", "+", (2127559, 2128886), frameshift=True),
    _Measured("PDT000167135.1", 218, 98.30, 100.00, "TAA", 219, True,
              "MPBQ01000002.1", "-", (388669, 390001), frameshift=True),
    _Measured("PDT000167136.1", 443, 97.97, 100.00, "TAA", 444, False,
              "MPBO01000003.1", "+", (1528885, 1530216),
              stop_is_reference_terminator=True),
    _Measured("PDT000292995.1", 372, 98.68, 99.77, "TGA", 373, False,
              "CP027172.1", "-", (2093476, 2094805), frameshift=True),
    _Measured("PDT000292998.1", 435, 95.91, 100.00, "TGA", 436, True,
              "CP027166.1", "-", (3778827, 3780163), frameshift=True),
    _Measured("PDT000294804.1", 197, 98.25, 100.00, "TGA", 198, False,
              "CP027538.1", "-", (4881206, 4883907), frameshift=True),
    _Measured("PDT000294805.1", 262, 96.60, 100.00, "TGA", 263, True,
              "PSQQ01000001.1", "+", (75456, 76767), frameshift=True),
    _Measured("PDT000311294.1", 262, 96.60, 100.00, "TGA", 263, True,
              "CP029148.1", "+", (3930548, 3931859), frameshift=True),
    _Measured("PDT000424983.1", 189, 98.00, 99.32, "TGA", 190, False,
              "RXUU01000038.1", "-", (83794, 85108), frameshift=True),
)


def _evidence(m: _Measured) -> LocusEvidence:
    return LocusEvidence(
        identity_pct=m.identity_pct,
        coverage_pct=m.coverage_pct,
        orf_aa_length=m.orf_aa,
        reference_aa_length=REFERENCE_AA,
        internal_stop_codon=m.stop_codon,
        internal_stop_position_aa=m.stop_position_aa,
        stop_is_reference_terminator=m.stop_is_reference_terminator,
        frameshift=m.frameshift,
        split_abutting_fragments=m.split,
        at_contig_edge=False,
        contig=m.contig,
        strand=m.strand,
        footprint_1based=m.footprint,
        n_hsp=2,
    )


def _classify(m: _Measured, **kwargs) -> StructuralCall:
    return classify_structure(m.sample_id, _evidence(m), search=_search(),
                              **kwargs)


class TestRealIsolateMeasurements:
    """The ten isolates, through the classifier, at the module's defaults."""

    def test_no_isolate_is_absent(self):
        """Every isolate returned a tblastn hit, so `absent` is unsupportable.

        This is the load-bearing negative: the previous round's symbol test
        reported `absent` for all ten, and all ten carry the locus at 94.59%
        identity or better.
        """
        verdicts = {_classify(m).verdict for m in MEASURED}
        assert StructuralVerdict.ABSENT not in verdicts

    @pytest.mark.parametrize(
        "sample_id,expected,lesion",
        [
            # lesion type decides; ORF length does not appear in this table.
            ("PDT000167136.1", StructuralVerdict.INTACT, oprd.LesionType.NONE),
            # 441 aa, stop at codon 442 - but that stop IS the native TAA,
            # reached in frame behind a net -6 nt in-frame deletion, so under
            # RULING 2 there is no lesion at all. See
            # TestRuling2AnInFrameDeletionIsNotDisruption.
            ("PDT000034122.1", StructuralVerdict.INTACT, oprd.LesionType.NONE),
            ("PDT000292998.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),                # 435 aa, frameshift + stop
            ("PDT000292995.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),                # 372 aa, frameshift + stop
            ("PDT000167133.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
            ("PDT000167135.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
            ("PDT000294804.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
            ("PDT000294805.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
            ("PDT000311294.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
            ("PDT000424983.1", StructuralVerdict.DISRUPTED,
             oprd.LesionType.BOTH),
        ],
    )
    def test_verdict_per_isolate(self, sample_id, expected, lesion):
        m = next(x for x in MEASURED if x.sample_id == sample_id)
        call = _classify(m)
        assert call.verdict == expected
        assert call.evidence.lesion_type == lesion

    def test_the_orf_length_banding_and_ruling_2_disagree_only_on_frameshifts(self):
        """Where the retired 0.90 band and RULING 2 part company, and why.

        The band compared `orf_length_fraction` against 0.90. That number cannot
        see a frameshift: a frameshifted locus can still reconstruct most of the
        reference protein's residues before its shifted reading frame dies, so
        the band called it `intact`. PDT000292998.1 is the case - 435 of 444
        codons, 0.9797, comfortably over the line, with a +5 nt insertion at
        reference codon 402 and a stop at codon 436 that exists only because of
        it.

        Where the band and RULING 2 agree is an in-frame deletion: the fraction
        stays high and RULING 2 calls it `intact` too. PDT000034122.1 is that
        case, at 0.9932. The band reached the right answer there, but by
        comparing a length to a chosen number rather than by observing that the
        reading frame never leaves the native frame - which is why the band also
        had to be 0.90 and not 0.98, and why it moved three verdicts when it was
        swept.
        """
        for sample_id, frac, expected in (
            ("PDT000292998.1", 0.9797, StructuralVerdict.DISRUPTED),
            ("PDT000034122.1", 0.9932, StructuralVerdict.INTACT),
        ):
            m = next(x for x in MEASURED if x.sample_id == sample_id)
            assert round(_evidence(m).orf_length_fraction, 4) == frac
            assert frac > 0.90            # the retired band said intact for both
            assert _classify(m).verdict == expected
        # The band cannot see a frameshift; the lesion can.
        assert next(m for m in MEASURED
                    if m.sample_id == "PDT000292998.1").frameshift is True
        assert next(m for m in MEASURED
                    if m.sample_id == "PDT000034122.1").frameshift is False

    @pytest.mark.parametrize(
        "sample_id,truncation",
        [
            ("PDT000167136.1", 0),     # complete: stop AT codon 444
            ("PDT000034122.1", 2),
            ("PDT000292998.1", 8),     # the reference C-terminal octamer
            ("PDT000292995.1", 71),
            ("PDT000167133.1", 207),
            ("PDT000167135.1", 225),
            ("PDT000294804.1", 246),
            ("PDT000294805.1", 181),
            ("PDT000311294.1", 181),
            ("PDT000424983.1", 254),
        ],
    )
    def test_truncation_is_recorded_as_evidence(self, sample_id, truncation):
        """Reported, never decided on. See `test_the_orf_length_cannot_move_a
        _verdict` for the second half of that claim."""
        m = next(x for x in MEASURED if x.sample_id == sample_id)
        call = _classify(m)
        assert call.evidence.truncation_aa == truncation
        assert call.as_row()["truncation_aa"] == truncation
        assert call.as_row()["lesion_type"] == call.evidence.lesion_type

    def test_a_complete_orf_has_zero_truncation(self):
        """The count comes from the lesion, so a terminator stop yields 0.

        `reference_aa_length - orf_aa_length` would say 1 for this isolate,
        because the reference figure of 444 counts the reference's own stop
        codon. Counting from the stop codon instead is exact.
        """
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        ev = _evidence(m)
        assert ev.orf_aa_length == REFERENCE_RESIDUES
        assert ev.reference_aa_length - ev.orf_aa_length == 1   # the naive count
        assert ev.truncation_aa == 0                           # the exact one

    def test_truncated_orf_names_its_stop_codon_and_position(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000167133.1")
        call = _classify(m)
        assert call.verdict == StructuralVerdict.DISRUPTED
        assert "TGA" in call.reason
        assert "237" in call.reason
        assert call.evidence.internal_stop_position_aa == 237
        assert call.evidence.orf_aa_length == 236

    def test_complete_orf_does_not_treat_the_terminator_as_a_lesion(self):
        """PAO1's own stop sits at codon 444; that is the terminator.

        Getting this wrong calls every complete oprD gene disrupted, because
        every ORF ends on a stop codon.
        """
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        call = _classify(m)
        assert call.verdict == StructuralVerdict.INTACT
        assert "internal_stop_codon" not in call.reason

    def test_every_verdict_carries_its_evidence(self):
        for m in MEASURED:
            call = _classify(m)
            assert call.evidence.identity_pct == m.identity_pct
            assert call.evidence.coverage_pct == m.coverage_pct
            assert call.evidence.orf_aa_length == m.orf_aa
            assert call.evidence.contig == m.contig
            assert call.evidence.internal_stop_codon == m.stop_codon
            assert call.search is not None and call.search.is_recorded
            row = call.as_row()
            assert row["structural_verdict"] == call.verdict
            assert row["search_recorded"] is True
            assert "db_gencode=11" in row["search_parameters"]

    def test_no_isolate_sits_at_a_contig_edge(self):
        """The nearest locus edge measured was 75455 bp.

        If any split locus were an assembly break rather than biology, this is
        where it would show.
        """
        for m in MEASURED:
            assert _evidence(m).at_contig_edge is False


class TestTheFourVerdicts:
    def test_intact_requires_a_long_orf_and_no_shortening_lesion(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        assert _classify(m).verdict == StructuralVerdict.INTACT

    def test_disrupted_from_a_short_orf(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000424983.1")
        call = _classify(m)
        assert call.verdict == StructuralVerdict.DISRUPTED
        assert call.is_loss_of_function

    def test_disrupted_from_a_frameshift_with_no_stop(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        ev = LocusEvidence(
            identity_pct=97.9, coverage_pct=100.0, orf_aa_length=120,
            reference_aa_length=REFERENCE_AA, frameshift=True,
            split_abutting_fragments=True, contig="c1", strand="+",
            footprint_1based=(1000, 1400),
        )
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.DISRUPTED
        assert "frameshift" in call.reason
        assert "abutting" in call.reason

    def test_absent_only_from_a_recorded_zero_hit_search(self):
        ev = LocusEvidence(
            identity_pct=None, coverage_pct=0.0, orf_aa_length=None,
            reference_aa_length=REFERENCE_AA, contig="", strand="",
        )
        call = classify_structure("s1", ev, search=_search(n_hits=0))
        assert call.verdict == StructuralVerdict.ABSENT
        assert call.is_absent
        assert "0 hits" in call.reason

    def test_not_assessed_when_no_search_was_recorded(self):
        ev = LocusEvidence(coverage_pct=0.0, orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)
        call = classify_structure("s1", ev, search=None)
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("no_recorded_search")
        assert not call.is_absent

    def test_not_assessed_when_the_locus_was_not_located(self):
        ev = LocusEvidence(identity_pct=38.0, coverage_pct=41.0,
                           orf_aa_length=180, reference_aa_length=REFERENCE_AA,
                           contig="c1", strand="+", footprint_1based=(1, 600))
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("locus_not_located")

    def test_not_assessed_when_no_orf_was_measured(self):
        ev = LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=None, reference_aa_length=REFERENCE_AA,
                           contig="c1", strand="+")
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("orf_not_measured")


class TestAbsenceRequiresARecordedSearch:
    """`absent` is the one claim that cannot be re-derived from the outputs."""

    def test_absent_is_unreachable_without_a_search_object(self):
        ev = LocusEvidence(coverage_pct=0.0, identity_pct=None,
                           orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)
        call = classify_structure("s1", ev, search=None)
        assert call.verdict != StructuralVerdict.ABSENT

    def test_absent_is_unreachable_when_the_search_still_found_hits(self):
        ev = LocusEvidence(coverage_pct=0.0, identity_pct=None,
                           orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)
        call = classify_structure("s1", ev, search=_search(n_hits=4))
        assert call.verdict != StructuralVerdict.ABSENT

    def test_the_recorded_search_is_carried_into_the_row(self):
        ev = LocusEvidence(coverage_pct=0.0, identity_pct=None,
                           orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)
        search = _search(n_hits=0)
        row = classify_structure("s1", ev, search=search).as_row()
        assert row["search_recorded"] is True
        assert row["search_command"] == search.command
        assert row["search_database"] == search.database
        assert row["search_n_hits"] == 0

    def test_an_unrecorded_search_object_is_not_treated_as_recorded(self):
        search = LocusSearch(program="", query="", database="",
                             parameters=(), command="", n_hits=0)
        assert search.is_recorded is False


class TestNotAssessedIsReachableAndNeverBecomesAbsent:
    def test_not_assessed_is_a_real_verdict_not_a_fallback(self):
        assert "not_assessed" in StructuralVerdict.ALL

    def test_there_is_no_undecided_length_band_to_fall_into(self):
        """The `orf_length_in_undecided_band` outcome is gone, with the banding.

        It existed to refuse to round a length to the nearer side. With no
        length threshold there is nothing to round: every located, translated
        locus with no lesion is `intact` and every one with a lesion is
        `disrupted`. The refusal survives in the four named reasons that remain.
        """
        assert not hasattr(oprd, "INTACT_MIN_ORF_LENGTH_FRACTION")
        assert not hasattr(oprd, "DISRUPTED_MAX_ORF_LENGTH_FRACTION")
        import inspect

        params = inspect.signature(classify_structure).parameters
        assert "intact_min_orf_fraction" not in params
        assert "disrupted_max_orf_fraction" not in params
        for m in MEASURED:
            call = _classify(m)
            assert "undecided_band" not in call.reason
            assert call.is_decided

    @pytest.mark.parametrize("m", MEASURED, ids=lambda m: m.sample_id)
    def test_no_measured_isolate_is_absent_or_undecided_at_the_defaults(self, m):
        """Every isolate gets a decisive verdict; none needs the undecided band."""
        call = _classify(m)
        assert call.verdict in (StructuralVerdict.INTACT,
                                StructuralVerdict.DISRUPTED)
        assert call.is_decided

    def test_every_not_assessed_reason_is_named(self):
        """A NAMED reason, so an undecided sample is never silently absent."""
        seen = set()
        recorded = _search()
        for kwargs, ev in (
            (dict(search=None), LocusEvidence(coverage_pct=0.0)),
            (dict(search=recorded),
             LocusEvidence(identity_pct=38.0, coverage_pct=41.0,
                           orf_aa_length=100,
                           reference_aa_length=REFERENCE_AA)),
            (dict(search=recorded),
             LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)),
            (dict(search=recorded),
             LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=300, reference_aa_length=REFERENCE_AA,
                           internal_stop_codon="TGA",
                           internal_stop_position_aa=None)),
            (dict(search=recorded),
             LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=300, reference_aa_length=REFERENCE_AA,
                           at_contig_edge=True)),
        ):
            call = classify_structure("s1", ev, **kwargs)
            assert call.verdict == StructuralVerdict.NOT_ASSESSED
            assert not call.is_absent
            seen.add(call.reason.split(":")[0])
        assert seen == {"no_recorded_search", "locus_not_located",
                        "orf_not_measured", "stop_position_unknown",
                        "locus_truncated_by_contig_edge"}


class TestNamedThresholdsAndSensitivity:
    def test_the_only_surviving_thresholds_are_presence_thresholds(self):
        """Aligned span and identity answer "is it there". Nothing else is cut."""
        assert oprd.STRUCTURAL_MIN_COVERAGE_PCT == oprd.MIN_COVERAGE_PCT
        assert oprd.STRUCTURAL_MIN_IDENTITY_PCT == oprd.MIN_IDENTITY_PCT
        assert oprd.MIN_COVERAGE_PCT == 80.0
        assert oprd.MIN_IDENTITY_PCT == 70.0

    def test_the_defaults_are_not_hidden_in_an_expression(self):
        import inspect

        sig = inspect.signature(classify_structure)
        assert sig.parameters["min_identity_pct"].default == \
            oprd.STRUCTURAL_MIN_IDENTITY_PCT
        assert sig.parameters["min_coverage_pct"].default == \
            oprd.STRUCTURAL_MIN_COVERAGE_PCT

    def test_the_orf_length_cannot_move_a_verdict(self):
        """FLATNESS, proved directly rather than asserted about a sweep.

        `orf_aa_length` is set to fourteen values spanning 1 aa to 1000 aa -
        including the real measurement for every isolate - and the verdict and
        the lesion type must not change. This is the property the old banding
        could not have: under 0.90, PDT000292995.1 at 0.8378 was `intact` and at
        0.85 `disrupted`.
        """
        probes = (1, 50, 100, 200, 300, 372, 400, 435, 441, 442, 443, 444,
                  500, 1000)
        for m in MEASURED:
            verdicts, lesions, truncations = set(), set(), set()
            for length in probes:
                ev = LocusEvidence(
                    identity_pct=m.identity_pct, coverage_pct=m.coverage_pct,
                    orf_aa_length=length, reference_aa_length=REFERENCE_AA,
                    internal_stop_codon=m.stop_codon,
                    internal_stop_position_aa=m.stop_position_aa,
                    frameshift=m.frameshift,
                    split_abutting_fragments=m.split,
                    contig=m.contig, strand=m.strand,
                    footprint_1based=m.footprint, n_hsp=2)
                call = classify_structure(m.sample_id, ev, search=_search())
                verdicts.add(call.verdict)
                lesions.add(ev.lesion_type)
                truncations.add(ev.truncation_aa)
            assert len(verdicts) == 1, (m.sample_id, verdicts)
            assert len(lesions) == 1, (m.sample_id, lesions)
            if next(iter(lesions)) != oprd.LesionType.NONE:
                # Truncation is counted FROM the lesion, so where there is a
                # lesion it is invariant too. Where there is none it falls back
                # to the encoded length and legitimately tracks it - which is
                # why it is an attribute and not a verdict.
                assert len(truncations) == 1, (m.sample_id, truncations)

    def test_the_lesion_alone_decides(self):
        """One evidence object, every lesion type, verdict follows the lesion."""
        def only(**over) -> StructuralCall:
            base = dict(identity_pct=98.0, coverage_pct=100.0,
                        orf_aa_length=400, reference_aa_length=REFERENCE_AA,
                        contig="c1", strand="+",
                        footprint_1based=(1000, 2300))
            base.update(over)
            return classify_structure("s1", LocusEvidence(**base),
                                      search=_search())

        # Same length in all four. Only the lesion differs.
        assert only().verdict == StructuralVerdict.INTACT
        assert only(internal_stop_codon="TGA",
                    internal_stop_position_aa=401).verdict == \
            StructuralVerdict.DISRUPTED
        assert only(frameshift=True).verdict == StructuralVerdict.DISRUPTED
        assert only(frameshift=True, internal_stop_codon="TAG",
                    internal_stop_position_aa=401).verdict == \
            StructuralVerdict.DISRUPTED

    def test_a_stop_at_the_reference_terminator_is_not_a_lesion(self):
        ev = LocusEvidence(identity_pct=98.0, coverage_pct=100.0,
                           orf_aa_length=REFERENCE_RESIDUES,
                           reference_aa_length=REFERENCE_AA,
                           internal_stop_codon="TAA",
                           internal_stop_position_aa=REFERENCE_AA,
                           stop_is_reference_terminator=True)
        assert ev.is_terminal_stop is True
        assert ev.has_premature_stop is False
        assert ev.lesion_type == oprd.LesionType.NONE
        assert classify_structure("s1", ev,
                                  search=_search()).verdict == \
            StructuralVerdict.INTACT


class TestRuling2AnInFrameDeletionIsNotDisruption:
    """RULING 2: `disrupted` means a FRAMESHIFT or a PREMATURE STOP.

    An in-frame deletion that preserves the NATIVE stop is not a frameshift and
    not a premature stop, so it is not `disrupted` - however many residues it
    removes.

    The concrete case is PDT000034122.1. Measured on its own assembly, its oprD
    is 1326 nt where PAO1's is 1332: a net -6 nt **in-frame** deletion, delivered
    as two indels 17 nt apart (DEL 5 at reference nt 1132-1136 and DEL 1 at 1149)
    inside a perfect 12-nt tandem repeat `CTACGG`x2 at reference nt 1143-1154.
    The isolate's ORF is 441 residues and its stop codon is TAA, identical to
    PAO1's, in frame, and terminal.

    The previous revision called this `disrupted` with
    `lesion_type=premature_stop`, because it decided "is this the terminator?"
    by comparing the stop's codon INDEX (442) against the reference's codon
    COUNT (444). That comparison is only valid when the isolate's reading frame
    is exactly as long as the reference's - which is precisely what an in-frame
    deletion breaks. The reason string even read "truncating 2 reference
    residues" while the stated rule was that truncation cannot decide anything.
    """

    def _ev(self, **over) -> LocusEvidence:
        base = dict(identity_pct=94.59, coverage_pct=100.00,
                    orf_aa_length=441, reference_aa_length=REFERENCE_AA,
                    internal_stop_codon="TAA",
                    internal_stop_position_aa=442,
                    frameshift=False, split_abutting_fragments=False,
                    at_contig_edge=False, contig="JFJU01000003.1", strand="+",
                    footprint_1based=(510136, 511461), n_hsp=1)
        base.update(over)
        return LocusEvidence(**base)

    def test_the_native_stop_preserved_behind_an_in_frame_deletion_is_intact(self):
        ev = self._ev(stop_is_reference_terminator=True)
        assert ev.is_terminal_stop is True
        assert ev.has_premature_stop is False
        assert ev.lesion_type == oprd.LesionType.NONE
        call = classify_structure("PDT000034122.1", ev, search=_search())
        assert call.verdict == StructuralVerdict.INTACT

    def test_the_measured_isolate_row_moves_to_intact(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000034122.1")
        assert _classify(m).verdict == StructuralVerdict.INTACT
        assert _classify(m).evidence.lesion_type == oprd.LesionType.NONE

    def test_the_truncation_is_still_reported_as_evidence(self):
        """Intact does not mean unmeasured: the two missing residues are still
        visible in the row, and still decide nothing."""
        call = classify_structure("PDT000034122.1",
                                  self._ev(stop_is_reference_terminator=True),
                                  search=_search())
        assert call.evidence.truncation_aa == 2
        assert call.as_row()["truncation_aa"] == 2
        assert call.as_row()["orf_aa_length"] == 441

    def test_a_nonsense_substitution_at_a_non_native_codon_is_still_disrupted(self):
        """The invariant this fix must not weaken.

        Same unbroken reading frame, same locus, but the stop codon was created
        by a substitution at a codon that is a normal residue in PAO1. That is a
        premature stop and it is `disrupted`.
        """
        ev = self._ev(internal_stop_codon="TAG",
                      internal_stop_position_aa=120,
                      orf_aa_length=119,
                      stop_is_reference_terminator=False)
        assert ev.is_terminal_stop is False
        assert ev.has_premature_stop is True
        assert ev.lesion_type == oprd.LesionType.PREMATURE_STOP
        assert classify_structure("s1", ev, search=_search()).verdict == \
            StructuralVerdict.DISRUPTED

    def test_a_frameshift_is_still_disrupted_behind_a_terminator_reading(self):
        """An in-frame verdict is not bought by claiming the stop is native."""
        ev = self._ev(frameshift=True, stop_is_reference_terminator=True)
        assert ev.lesion_type == oprd.LesionType.FRAMESHIFT
        assert classify_structure("s1", ev, search=_search()).verdict == \
            StructuralVerdict.DISRUPTED

    def test_an_unmeasured_terminator_falls_back_to_the_index_test(self):
        """Absent the measurement the module stays conservative.

        `stop_is_reference_terminator=None` means the question was not asked, so
        the old index comparison stands. Defaulting it to True would manufacture
        intact calls out of unmeasured evidence, which is the one thing this
        module exists to prevent.
        """
        ev = self._ev(stop_is_reference_terminator=None)
        assert ev.is_terminal_stop is False
        assert classify_structure("s1", ev, search=_search()).verdict == \
            StructuralVerdict.DISRUPTED

    def test_the_measurement_cannot_rescue_a_stop_the_frame_never_reaches(self):
        """The index test and the measurement agree in the normal case.

        A complete ORF satisfies both, so requiring both cannot regress the
        one isolate that was already intact.
        """
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        ev = _evidence(m)
        assert ev.stop_is_reference_terminator is True
        assert ev.internal_stop_position_aa >= REFERENCE_AA
        assert ev.is_terminal_stop is True

    def test_a_split_locus_is_not_a_lesion_type(self):
        """The reading frame runs through the assembly gap; the protein is one.

        Split is recorded on the evidence and named in the reason so a reader
        sees it, but it is an annotation property, not a reading-frame lesion,
        so on its own it cannot make a locus `disrupted`.
        """
        m = next(x for x in MEASURED if x.sample_id == "PDT000167136.1")
        ev = LocusEvidence(
            identity_pct=m.identity_pct, coverage_pct=m.coverage_pct,
            orf_aa_length=m.orf_aa, reference_aa_length=REFERENCE_AA,
            internal_stop_codon=m.stop_codon,
            internal_stop_position_aa=m.stop_position_aa,
            split_abutting_fragments=True, contig=m.contig, strand=m.strand,
            footprint_1based=m.footprint, n_hsp=2)
        assert ev.lesion_type == oprd.LesionType.NONE
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.INTACT
        assert "not a reading-frame lesion" in call.reason

    def test_a_frameshift_is_disrupted_with_no_orf_length_measured(self):
        """The lesion is sufficient. The old code asked for a length first and
        returned `orf_not_measured` on this input."""
        ev = LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=None, reference_aa_length=REFERENCE_AA,
                           frameshift=True, contig="c1", strand="+",
                           footprint_1based=(1000, 1400))
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.DISRUPTED
        assert call.evidence.truncation_aa is None
        assert "frameshift" in call.reason

    def test_a_contig_edge_cannot_be_called_intact(self):
        """Unknown extent is not intact, and it is not `disrupted` either."""
        ev = LocusEvidence(identity_pct=97.0, coverage_pct=99.0,
                           orf_aa_length=400, reference_aa_length=REFERENCE_AA,
                           at_contig_edge=True, contig="c1", strand="+")
        call = classify_structure("s1", ev, search=_search())
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("locus_truncated_by_contig_edge")
        assert not call.is_absent

    def test_sensitivity_table_reports_every_isolate_at_every_threshold(self):
        calls = [_classify(m) for m in MEASURED]
        rows = sensitivity_table(calls)
        assert len(rows) == len(MEASURED)
        verdicts = {r[k] for r in rows for k in r
                    if k.startswith(("id_floor_", "cov_floor_"))}
        assert verdicts <= StructuralVerdict.ALL
        assert StructuralVerdict.ABSENT not in verdicts
        for row in rows:
            for key in ("lesion_type",):
                assert row[key] in oprd.LesionType.ALL

    def test_the_table_has_no_orf_length_sweep_left_to_move(self):
        """The swept columns are presence columns. That is the whole claim."""
        calls = [_classify(m) for m in MEASURED]
        rows = sensitivity_table(calls)
        swept = {k for r in rows for k in r
                 if k.startswith(("id_floor_", "cov_floor_"))}
        assert swept, "the presence sweep disappeared entirely"
        assert not any(k.startswith("orf") and "floor" in k
                       for r in rows for k in r)
        # The length is still reported, as evidence, on every row.
        for row in rows:
            assert "orf_aa_length" in row
            assert "orf_length_fraction" in row
            assert "truncation_aa" in row
            assert "lesion_type" in row

    def test_the_only_isolate_that_moves_is_the_weakest_on_identity(self):
        """One movement, and it is the PRESENCE floor, not the frame.

        PDT000034122.1 sits at 94.59% identity; raising the identity floor to 95
        stops the locus being located at all, so it leaves the verdict space
        rather than being re-graded. No isolate moves on its reading frame.
        """
        calls = [_classify(m) for m in MEASURED]
        rows = sensitivity_table(calls)
        moved = {r["sample_id"] for r in rows
                 if len({r[k] for k in r
                         if k.startswith(("id_floor_", "cov_floor_"))}) > 1}
        assert moved == {"PDT000034122.1"}
        row = next(r for r in rows if r["sample_id"] == "PDT000034122.1")
        assert row["id_floor_95"] == StructuralVerdict.NOT_ASSESSED
        assert row["id_floor_95"] != StructuralVerdict.ABSENT
        for r in rows:
            assert r["verdict"] == _classify(
                next(m for m in MEASURED if m.sample_id == r["sample_id"])
            ).verdict

    def test_identity_floor_95_pushes_the_weakest_isolate_out_of_a_verdict(self):
        """The presence floor, and it leaves the verdict space rather than
        flipping within it. Never `absent`: a low identity is a paralog, not a
        deletion."""
        m = next(x for x in MEASURED if x.sample_id == "PDT000034122.1")
        assert m.identity_pct == 94.59
        assert _classify(m, min_identity_pct=90.0).verdict == \
            StructuralVerdict.INTACT
        assert _classify(m, min_identity_pct=95.0).verdict == \
            StructuralVerdict.NOT_ASSESSED
        assert _classify(m, min_identity_pct=95.0).reason.startswith(
            "locus_not_located")


class TestStatusFunctionTakesTheStructuralVerdict:
    """`oprd_status_per_sample` consumes the classifier's vocabulary directly."""

    def test_absent_state_comes_only_from_an_absent_verdict(self):
        ev = LocusEvidence(coverage_pct=0.0, orf_aa_length=None,
                           reference_aa_length=REFERENCE_AA)
        call = classify_structure("s1", ev, search=_search(n_hits=0))
        status = viz.oprd_status_per_sample({"s1": []}, None, {"s1": call})
        assert status["s1"] == "absent"

    def test_loss_of_function_state_comes_only_from_a_disrupted_verdict(self):
        m = next(x for x in MEASURED if x.sample_id == "PDT000167133.1")
        call = _classify(m)
        assert call.verdict == "disrupted"
        status = viz.oprd_status_per_sample({"s1": []}, None, {"s1": call})
        assert status["s1"] == "disrupted"

    def test_a_refused_resolution_is_not_assessed_not_absent(self):
        @dataclass
        class _Refused:
            sample_id: str = "s1"
            verdict: str = "refused:insufficient_coverage"
            is_resolved: bool = False

        status = viz.oprd_status_per_sample({"s1": []}, None,
                                            {"s1": _Refused()})
        assert status["s1"] == "not_assessed"

    def test_the_symbol_test_can_no_longer_report_absent(self):
        """The regression this whole change exists to prevent.

        `absent` used to be the else-branch of `"oprD" in genes`, which reported
        the gene missing for eight of ten isolates that all carry it.
        """
        status = viz.oprd_status_per_sample({"s1": []}, {"s1": ["mexR"]}, None)
        assert status["s1"] == "not_assessed"
        assert status["s1"] != "absent"

    def test_the_resolver_stays_optional(self):
        """Kept optional on purpose: it changes what a caller reports, not
        whether it runs. Both call shapes work."""
        assert viz.oprd_status_per_sample({"s1": []}, {"s1": ["oprD"]}) == \
            {"s1": "intact"}
        assert viz.oprd_status_per_sample({"s1": []}) == {}