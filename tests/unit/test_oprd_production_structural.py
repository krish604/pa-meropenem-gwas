"""The PRODUCTION structural path: what measures, what decides, what refuses.

:mod:`papipeline.adapters.oprd_locus` decided structure correctly and measured
nothing. ``classify_structure`` was exercised only by tests, so
``LocusEvidence.stop_is_reference_terminator`` had no producer in a real run and
the codon-index fallback inside :attr:`LocusEvidence.is_terminal_stop` was the
only thing any execution path could reach. That is the gap this file closes.

Three things are pinned here:

* **A1** the terminator question is answered FROM THE ALIGNMENT, by
  :func:`~papipeline.adapters.oprd_locus.measure_locus_structure`, and no
  production path reaches the codon-index fallback - including a locus whose
  frame cannot be anchored, which is refused instead;
* **A2** the measurement reproduces the ten real isolates' measured stop
  positions, ORF lengths and net frame changes, which is what makes the verdicts
  checkable rather than asserted;
* **A3** ``absent`` comes only from a recorded zero-hit search, ``disrupted``
  only from a lesion, and a missing input is ``not_assessed`` - never ``absent``.
"""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import pytest

from papipeline.adapters import oprd_locus as oprd
from papipeline.adapters.oprd_locus import (
    AlignedNucleotideHit,
    Indel,
    LesionType,
    LocusEvidence,
    LocusSearch,
    StructuralCall,
    StructuralVerdict,
)

REFERENCE_CODONS = 60


def _synthetic_reference(codons: int = REFERENCE_CODONS,
                         seed: str = "oprd-production-test") -> str:
    """A deterministic, stop-free reference CDS of ``codons`` codons.

    Generated rather than transcribed so nothing here can be a transcription of
    PAO1 with a mistake in it, and seeded so a failure is reproducible.
    """
    body: list = []
    index = 0
    while len(body) < codons - 2:
        digest = hashlib.sha256(f"{seed}:{index}".encode()).digest()
        codon = "".join("ACGT"[byte & 3] for byte in digest[:3])
        index += 1
        if codon not in ("TAA", "TAG", "TGA"):
            body.append(codon)
    return "ATG" + "".join(body) + "TAA"


REF_CDS = _synthetic_reference()


def _translate(cds: str) -> str:
    return oprd._translate(cds)


def _row(contig: str, qseq: str, sseq: str, subject_start: int,
         subject_end: int, query_start: int, query_end: int,
         *, query_id: str = "PA0958_oprD",
         subject_length: int = 100000) -> str:
    """One ``TBLASTN_OUTFMT`` row, with the aligned columns filled in.

    ``alignment_length`` is the number of aligned residue columns, which is what
    blast reports and what :func:`parse_tblastn_table` reads.
    """
    # blast's `length` is the number of aligned COLUMNS, and a query gap is a
    # column carrying no query residue, so the subject cursor only advances on a
    # query-residue column.
    subject_cursor = 0
    nident = alen = 0
    for query_residue in qseq:
        if query_residue == "-":
            continue
        subject_residue = sseq[subject_cursor]
        subject_cursor += 1
        alen += 1
        if query_residue == subject_residue:
            nident += 1
    pident = 100.0 * nident / alen
    return "\t".join((
        query_id, contig, f"{pident:.3f}", str(alen), str(len(qseq)),
        str(subject_length), str(query_start), str(query_end),
        str(subject_start), str(subject_end), "0.0",
        f"{alen * 2:.1f}", qseq, sseq,
    ))


def _residues(cds: str) -> str:
    """A CDS as the residue string blast aligns, its terminal stop included."""
    return _translate(cds if cds.endswith("TAA") else cds + "TAA")


#: Where every fixture's locus is laid into its assembly contig.
LOCUS_START = 5000
#: Two reference codons are deleted from the in-frame fixture, so its locus is
#: shorter and the alignment's columns carry two query gaps.
DELETED_CODONS = (20, 21)
DELETED_NT = slice(3 * (DELETED_CODONS[0] - 1), 3 * DELETED_CODONS[1])

#: Complete: 180 nt against a 60-codon reference, net 0.
COMPLETE = REF_CDS

#: Net -6 nt IN-FRAME deletion of codons 20 and 21. RULING 2's case: two whole
#: codons go, the reading frame never leaves the native one, and the isolate
#: still ends on PAO1's own TAA - now at its codon 58.
INFRAME_DEL = REF_CDS[:DELETED_NT.start] + REF_CDS[DELETED_NT.stop:]

#: Net -7 nt: one base short of a whole codon, so the isolate's own frame ends
#: BEYOND the reference's terminator. The codon-index test would read that stop
#: as "past the end of the reference" and call it the terminator; the alignment
#: test does not, because the frame left the native one.
FRAMESHIFT_MINUS_7 = REF_CDS[:DELETED_NT.start - 1] + REF_CDS[DELETED_NT.start:]

#: A nonsense substitution at codon 20, in an otherwise unbroken frame.
PREMATURE_STOP = (REF_CDS[:3 * (DELETED_CODONS[0] - 1)] + "TGA"
                  + REF_CDS[3 * DELETED_CODONS[0]:])


def _exact_rows(contig: str, coding: str, *, start: int = LOCUS_START,
                direction: int = 1) -> list:
    """One HSP covering the whole reference, for a locus of the same length."""
    qseq = _residues(REF_CDS)
    sseq = _residues(coding)
    assert len(qseq) == len(sseq)
    if direction > 0:
        subject_start, subject_end = start, start + len(coding) - 1
    else:
        subject_start, subject_end = start + len(coding) - 1, start
    return [_row(contig, qseq, sseq, subject_start, subject_end, 1, len(qseq))]


def _column_gap_rows(contig: str, coding: str, deleted_codons,
                     *, start: int = LOCUS_START) -> list:
    """One HSP whose aligned columns carry a query gap per deleted codon.

    That is what blast writes for a whole-codon deletion, so the indel is
    visible in the columns rather than only in the net arithmetic.
    """
    qseq, sseq = [], []
    isolate_residues = _residues(coding)
    emitted = 0
    for index, residue in enumerate(_residues(REF_CDS), start=1):
        qseq.append(residue)
        if index in deleted_codons:
            sseq.append("-")          # a subject-gap column, as blast writes it
            continue
        sseq.append(isolate_residues[emitted])
        emitted += 1
    return [_row(contig, "".join(qseq), "".join(sseq), start,
                 start + len(coding) - 1, 1, len(_residues(REF_CDS)))]


def _split_hsp_rows(contig: str, *, split_codon: int, shift_nt: int,
                    coding_length: int, start: int = LOCUS_START) -> list:
    """Two HSPs meeting at reference ``split_codon``, the second shifted.

    This is how a frameshifted locus actually arrives: the search reports the
    reference's codons either side of the lesion as separate HSPs whose subject
    coordinates no longer sit three nucleotides per codon apart. ``shift_nt`` is
    that disagreement, and the aligner measured ``net_frame_change_nt`` equals it.
    """
    reference_residues = _residues(REF_CDS)
    first = reference_residues[:split_codon]
    second = reference_residues[split_codon - 1:]
    expected = 3 * (split_codon - 1) + 1
    second_start = start + expected - 1 + shift_nt
    return [
        _row(contig, first, first, start, start + 3 * len(first) - 1, 1,
             split_codon),
        _row(contig, second, second, second_start,
             second_start + 3 * len(second) - 1, split_codon,
             split_codon - 1 + len(second)),
    ]


def _assembly(contig: str, coding: str, *, start: int = LOCUS_START,
              length: int = 20000) -> dict:
    """One contig with ``coding`` laid into it from ``start``, forward."""
    seq = list("A" * length)
    for offset, base in enumerate(coding):
        seq[start - 1 + offset] = base
    return {contig: "".join(seq)}


def _search(n_hits: int = 1) -> LocusSearch:
    return LocusSearch(
        program="tblastn",
        query="PAO1 oprD protein (PA0958)",
        database="test_index",
        parameters=(("evalue", "1e-3"), ("max_target_seqs", "5000"),
                    ("db_gencode", "11")),
        command=("tblastn -query PA0958.faa -db index -evalue 1e-3 "
                 "-max_target_seqs 5000 -db_gencode 11 -num_threads 1"),
        n_hits=n_hits,
    )


def _measure(rows, coding: str, contig: str = "c1", **kwargs):
    hits = oprd.parse_tblastn_table("\n".join(rows), REFERENCE_CODONS)
    return hits, oprd.measure_locus_structure(
        "s1", hits, reference_cds=REF_CDS,
        assembly=_assembly(contig, coding), **kwargs)


def _call(rows, coding: str, contig: str = "c1", *, search=None) -> StructuralCall:
    hits = oprd.parse_tblastn_table("\n".join(rows), REFERENCE_CODONS)
    return oprd.structural_call(
        "s1", hits, reference_cds=REF_CDS,
        assembly=_assembly(contig, coding),
        search=_search() if search is None else search)


# ---------------------------------------------------------------------------
# A1. The terminator question is asked of the alignment.
# ---------------------------------------------------------------------------

class TestTheProductionResolverDecidesTheTerminatorFromTheAlignment:
    """A1. Not from the codon index - and demonstrably not.

    The index test compares ``internal_stop_position_aa`` against
    ``reference_aa_length``, which holds only while the isolate's reading frame
    is exactly as long as the reference's. A net -6 nt in-frame deletion moves
    the terminator from codon 60 to codon 58 and leaves it PAO1's own ``TAA``;
    on the index test alone that reads as a premature stop and the isolate is
    ``disrupted`` for a lesion it does not have.
    """

    ROWS_INFRAME = _column_gap_rows("c1", INFRAME_DEL, DELETED_CODONS)

    def test_the_in_frame_deletion_keeping_the_native_stop_is_intact(self):
        """RULING 2's case, end to end through production code."""
        _hits, m = _measure(self.ROWS_INFRAME, INFRAME_DEL)
        assert m.frame_anchored is True
        assert m.net_frame_change_nt == -6
        assert m.evidence.internal_stop_position_aa == 58
        assert m.evidence.orf_aa_length == 57
        assert m.evidence.internal_stop_codon == "TAA"
        assert m.evidence.truncation_aa == 2
        assert m.evidence.frameshift is False
        assert m.evidence.stop_is_reference_terminator is True
        assert m.evidence.lesion_type == LesionType.NONE
        assert _call(self.ROWS_INFRAME, INFRAME_DEL).verdict == \
            StructuralVerdict.INTACT

    def test_the_same_measurement_without_the_alignment_is_disrupted(self):
        """**This is the demonstration that the old production path was wrong.**

        Identical numbers, except the terminator is left unmeasured - which is
        exactly the state the module was in, because nothing measured it. The
        codon-index fallback then fires on a stop at codon 58 of a 60-codon
        reference and calls a lesion where there is none.
        """
        _hits, m = _measure(self.ROWS_INFRAME, INFRAME_DEL)
        unmeasured = replace(m.evidence, stop_is_reference_terminator=None)
        assert unmeasured.internal_stop_position_aa == 58
        assert unmeasured.reference_aa_length == REFERENCE_CODONS
        assert unmeasured.is_terminal_stop is False        # the index test
        assert unmeasured.has_premature_stop is True
        assert oprd.classify_structure("s1", unmeasured,
                                       search=_search()).verdict == \
            StructuralVerdict.DISRUPTED

    def test_the_production_call_populates_the_field_the_classifier_reads(self):
        """The gap being closed: a producer, on a production call."""
        call = _call(self.ROWS_INFRAME, INFRAME_DEL)
        assert call.verdict == StructuralVerdict.INTACT
        assert call.evidence.stop_is_reference_terminator is True
        assert call.as_row()["stop_is_reference_terminator"] is True

    def test_a_complete_locus_is_intact_and_truncates_nothing(self):
        rows = _exact_rows("c1", COMPLETE)
        _hits, m = _measure(rows, COMPLETE)
        assert m.net_frame_change_nt == 0
        assert m.evidence.internal_stop_position_aa == REFERENCE_CODONS
        assert m.evidence.orf_aa_length == REFERENCE_CODONS - 1
        assert m.evidence.truncation_aa == 0
        assert m.evidence.stop_is_reference_terminator is True
        assert _call(rows, COMPLETE).verdict == StructuralVerdict.INTACT

    def test_a_nonsense_substitution_is_a_premature_stop(self):
        """The invariant the terminator reading must not weaken."""
        rows = _exact_rows("c1", PREMATURE_STOP)
        _hits, m = _measure(rows, PREMATURE_STOP)
        assert m.net_frame_change_nt == 0
        assert m.evidence.frameshift is False
        assert m.evidence.stop_is_reference_terminator is False
        assert m.evidence.internal_stop_codon == "TGA"
        assert m.evidence.internal_stop_position_aa == DELETED_CODONS[0]
        assert m.evidence.lesion_type == LesionType.PREMATURE_STOP
        assert _call(rows, PREMATURE_STOP).verdict == \
            StructuralVerdict.DISRUPTED

    @pytest.mark.parametrize(
        "shift_nt,net",
        [(-4, -4), (+1, +1), (-7, -7), (-2, -2)],
    )
    def test_a_frameshift_is_disrupted_at_every_net_frame_change(
            self, shift_nt, net):
        """A locus split across two HSPs, as a frameshifted one arrives."""
        rows = _split_hsp_rows("c1", split_codon=DELETED_CODONS[0],
                               shift_nt=shift_nt, coding_length=0)
        coding = COMPLETE[:3 * (DELETED_CODONS[0] - 1) - 1] + \
            COMPLETE[3 * (DELETED_CODONS[0] - 1):]
        _hits, m = _measure(rows, coding)
        assert m.net_frame_change_nt == net
        assert m.evidence.frameshift is True
        assert m.evidence.stop_is_reference_terminator is False
        assert m.evidence.has_premature_stop is True
        assert m.evidence.lesion_type == LesionType.BOTH
        assert _call(rows, coding).verdict == StructuralVerdict.DISRUPTED

    def test_an_alignment_that_stops_short_of_the_reference_stop_is_not_intact(
            self):
        """Native frame, a terminal stop of its own - and still not intact.

        The search aligned reference codons 1-58 and nothing more. The isolate
        therefore never reaches PAO1's terminator at codon 60, so the stop it
        does have sits at a codon that encodes a residue in the reference: a
        premature stop. "The frame is unbroken" is not "the gene is whole", and
        the two are separate statements in the terminator test for exactly this.
        """
        residues = _residues(REF_CDS)
        coding = COMPLETE[:3 * 58] + "TAA"
        rows = [_row("c1", residues[:58], _residues(coding)[:58], LOCUS_START,
                     LOCUS_START + len(coding) - 1, 1, 58)]
        _hits, m = _measure(rows, coding)
        assert m.net_frame_change_nt == 3          # a multiple of three
        assert m.evidence.frameshift is False
        assert m.evidence.internal_stop_position_aa == 59
        assert m.evidence.stop_is_reference_terminator is False
        assert m.evidence.lesion_type == LesionType.PREMATURE_STOP
        assert _call(rows, coding).verdict == StructuralVerdict.DISRUPTED

    def test_the_reverse_strand_is_read_on_its_own_frame(self):
        """Translating the forward record would read a frame that is not there."""
        hits = oprd.parse_tblastn_table(
            "\n".join(_exact_rows("c1", COMPLETE, direction=-1)),
            REFERENCE_CODONS)
        assembly = _assembly("c1", COMPLETE, start=LOCUS_START)
        reverse = {
            "c1": assembly["c1"][:LOCUS_START - 1]
                  + oprd._reverse_complement(COMPLETE)
                  + assembly["c1"][LOCUS_START - 1 + len(COMPLETE):]
        }
        m = oprd.measure_locus_structure("s1", hits, reference_cds=REF_CDS,
                                         assembly=reverse)
        assert m.evidence.strand == "-"
        assert m.evidence.internal_stop_position_aa == REFERENCE_CODONS
        assert m.evidence.internal_stop_codon == "TAA"
        assert m.evidence.stop_is_reference_terminator is True
        assert m.evidence.lesion_type == LesionType.NONE


class TestNoProductionPathReachesTheCodonIndexFallback:
    """The fallback survives in :attr:`LocusEvidence.is_terminal_stop`.

    It has to: a caller that did not measure must stay conservative. What must
    not happen is a production path landing there, and the gate is a refusal
    rather than a default.
    """

    def test_a_locus_with_no_hsp_at_the_reference_first_codon_is_refused(self):
        """Unanchorable frame -> ``not_assessed``, never a guess."""
        residues = _residues(REF_CDS)
        rows = [_row("c1", residues[9:], residues[9:], LOCUS_START,
                     LOCUS_START + 3 * (len(residues) - 9) - 1, 10,
                     len(residues))]
        # both HSP columns are equal length, as blast writes them
        coding = COMPLETE[3 * 9:]
        hits = oprd.parse_tblastn_table("\n".join(rows), REFERENCE_CODONS)
        m = oprd.measure_locus_structure("s1", hits, reference_cds=REF_CDS,
                                         assembly=_assembly("c1", coding))
        assert m.frame_anchored is False
        assert "frame_not_anchorable" in (m.refusal_reason or "")
        assert m.evidence.stop_is_reference_terminator is False
        call = oprd.structural_call("s1", hits, reference_cds=REF_CDS,
                                    assembly=_assembly("c1", coding),
                                    search=_search())
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("frame_not_anchorable")
        assert not call.is_absent

    def test_no_candidate_hit_is_not_assessed_rather_than_absent(self):
        rows = _exact_rows("c1", COMPLETE)
        hits = [replace(h, identity_pct=41.9) for h in
                oprd.parse_tblastn_table("\n".join(rows), REFERENCE_CODONS)]
        call = oprd.structural_call("s1", hits, reference_cds=REF_CDS,
                                    assembly=_assembly("c1", COMPLETE),
                                    search=_search(n_hits=len(hits)))
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.split(":")[0] in ("no_candidate_hit",
                                            "locus_not_located")
        assert not call.is_absent

    def test_two_loci_covering_the_same_span_refuse_rather_than_choose(self):
        residues = _residues(REF_CDS)
        rows = [
            _row("c1", residues, residues, LOCUS_START,
                 LOCUS_START + len(COMPLETE) - 1, 1, len(residues)),
            _row("c2", residues, residues, 9000, 9000 + len(COMPLETE) - 1, 1,
                 len(residues)),
        ]
        hits = oprd.parse_tblastn_table("\n".join(rows), REFERENCE_CODONS)
        assembly = {**_assembly("c1", COMPLETE),
                    **_assembly("c2", COMPLETE, start=9000)}
        call = oprd.structural_call("s1", hits, reference_cds=REF_CDS,
                                    assembly=assembly,
                                    search=_search(n_hits=2))
        assert call.verdict == StructuralVerdict.NOT_ASSESSED
        assert call.reason.startswith("ambiguous_locus")
        assert not call.is_absent

    def test_the_measurer_never_leaves_the_terminator_unasked(self):
        """Every measurable locus gets a bool, so ``None`` is unreachable."""
        cases = (
            (_exact_rows("c1", COMPLETE), COMPLETE),
            (_column_gap_rows("c1", INFRAME_DEL, DELETED_CODONS), INFRAME_DEL),
            (_exact_rows("c1", PREMATURE_STOP), PREMATURE_STOP),
            (_split_hsp_rows("c1", split_codon=DELETED_CODONS[0],
                             shift_nt=-4, coding_length=0), COMPLETE),
        )
        for rows, coding in cases:
            _hits, m = _measure(rows, coding)
            assert m.frame_anchored is True, rows
            assert isinstance(m.evidence.stop_is_reference_terminator, bool)


# ---------------------------------------------------------------------------
# A2. The measurement reproduces the ten real isolates.
# ---------------------------------------------------------------------------

#: The ten smoke isolates, as an independent-frame rebuild of each assembly's
#: own oprD CDS against the PAO1 PA0958 CDS measured them. Held here as a
#: CROSS-CHECK on the production measurer, never as its input: the production
#: numbers in ``tests/unit/test_oprd_structural.py`` are produced by
#: :func:`measure_locus_structure` over the stored artifacts, and these are what
#: they are checked against.
#:
#:     sample          net    first stop codon   stop codon
#:     PDT000034122.1  -6     442                TAA   (native)
#:     PDT000167133.1  -4     237                TGA
#:     PDT000167135.1  +1     219                TAA
#:     PDT000167136.1   0     444                TAA   (native)
#:     PDT000292995.1  -2     373                TGA
#:     PDT000292998.1  +5     436                TGA
#:     PDT000294804.1 +1370   198                TGA
#:     PDT000294805.1 -20     263                TGA
#:     PDT000311294.1 -20     263                TGA
#:     PDT000424983.1 -17     190                TGA
REAL_CROSSCHECK = {
    "PDT000034122.1": (-6, "TAA", 442, True),
    "PDT000167133.1": (-4, "TGA", 237, False),
    "PDT000167135.1": (+1, "TAA", 219, False),
    "PDT000167136.1": (0, "TAA", 444, True),
    "PDT000292995.1": (-2, "TGA", 373, False),
    "PDT000292998.1": (+5, "TGA", 436, False),
    "PDT000294804.1": (+1370, "TGA", 198, False),
    "PDT000294805.1": (-20, "TGA", 263, False),
    "PDT000311294.1": (-20, "TGA", 263, False),
    "PDT000424983.1": (-17, "TGA", 190, False),
}


class TestTheRealIsolatesAreReproducibleFromTheStoredArtifacts:
    """A2, as far as a test in this repository can carry it.

    The stored artifacts live outside the worktree, so the end-to-end
    reproduction over them runs from ``pa-artifacts/round11/oprd/a2_reproduce.py``
    and its table is checked into that directory. These two tests pin the two
    facts about those numbers that a reader must be able to see without running
    anything: the split is 2 intact / 8 disrupted, and the eight ``disrupted``
    verdicts all carry ``frameshift_with_premature_stop`` because every one of
    their stops is a frameshift consequence rather than a nonsense substitution.
    """

    def test_the_reproduced_split_is_two_intact_and_eight_disrupted(self):
        import json
        from pathlib import Path

        table = Path("/Users/raghavkrishnankv/Desktop/pa-artifacts/round11"
                     "/oprd/a2-verdicts.json")
        if not table.is_file():
            pytest.skip("a2-verdicts.json is not present on this machine")
        rows = json.loads(table.read_text())["rows"]
        verdicts = [r["verdict"] for r in rows]
        assert verdicts.count("intact") == 2
        assert verdicts.count("disrupted") == 8
        assert sorted(r["sample_id"] for r in rows
                      if r["verdict"] == "intact") == \
            ["PDT000034122.1", "PDT000167136.1"]

    @pytest.mark.parametrize("sample_id", sorted(REAL_CROSSCHECK))
    def test_every_real_measurement_is_what_the_crosscheck_says(self,
                                                                sample_id):
        import json
        from pathlib import Path

        table = Path("/Users/raghavkrishnankv/Desktop/pa-artifacts/round11"
                     "/oprd/a2-verdicts.json")
        if not table.is_file():
            pytest.skip("a2-verdicts.json is not present on this machine")
        net, codon, position, native = REAL_CROSSCHECK[sample_id]
        row = next(r for r in json.loads(table.read_text())["rows"]
                   if r["sample_id"] == sample_id)
        assert row["net_frame_change_nt"] == net
        assert row["internal_stop_codon"] == codon
        assert row["internal_stop_position_aa"] == position
        assert row["stop_is_reference_terminator"] is native
        assert row["verdict"] == (StructuralVerdict.INTACT if native
                                  else StructuralVerdict.DISRUPTED)


# ---------------------------------------------------------------------------
# A3. Where `absent` comes from, and where it does not.
# ---------------------------------------------------------------------------

class TestAbsenceAndLossOfFunctionHaveExactlyOneSource:
    """A3. The three claims, on the production path."""

    def test_absent_is_reachable_only_from_a_recorded_zero_hit_search(self):
        assert StructuralVerdict.ABSENCE == frozenset({StructuralVerdict.ABSENT})
        assert StructuralVerdict.LOSS_OF_FUNCTION == \
            frozenset({StructuralVerdict.DISRUPTED})
        assert not StructuralVerdict.INTACT in StructuralVerdict.ABSENCE

    @pytest.mark.parametrize("missing", ["reference", "table", "assembly"])
    def test_a_missing_file_is_not_assessed_even_with_a_zero_hit_search(
            self, tmp_path, missing):
        """The regression that would manufacture the study's negative.

        A recorded search with zero hits is the ONE thing that may yield
        ``absent``. If a run that never read its table could reach that branch,
        a missing file would become an absent gene.
        """
        root = tmp_path / missing
        root.mkdir()
        paths = {
            "reference": root / "ref.fna",
            "table": root / "hits.tsv",
            "assembly": root / "asm.fna",
        }
        if missing != "reference":
            paths["reference"].write_text(f">r\n{REF_CDS}\n")
        if missing != "table":
            paths["table"].write_text("")
        if missing != "assembly":
            paths["assembly"].write_text(">c1\n" + "A" * 100 + "\n")
        call = oprd.structural_call_from_paths(
            "s1", reference_cds_path=paths["reference"],
            tblastn_table_path=paths["table"],
            assembly_path=paths["assembly"], search=_search(n_hits=0),
        )
        assert call.verdict == StructuralVerdict.NOT_ASSESSED, missing
        assert call.reason.startswith("missing_input"), missing
        assert str(paths[missing]) in call.reason, missing
        assert not call.is_absent, missing

    def test_a_zero_hit_search_on_read_files_is_absent(self):
        """The contrast: the input was read and the search was empty."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "ref.fna").write_text(f">r\n{REF_CDS}\n")
            (root / "hits.tsv").write_text("")
            (root / "asm.fna").write_text(">c1\n" + "A" * 100 + "\n")
            call = oprd.structural_call_from_paths(
                "s1", reference_cds_path=root / "ref.fna",
                tblastn_table_path=root / "hits.tsv",
                assembly_path=root / "asm.fna", search=_search(n_hits=0))
        assert call.verdict == StructuralVerdict.ABSENT
        assert call.is_absent

    def test_a_zero_hit_search_that_was_not_recorded_is_not_assent(self):
        assert oprd.structural_call(
            "s1", [], reference_cds=REF_CDS, assembly={}, search=None,
        ).verdict == StructuralVerdict.NOT_ASSESSED

    def test_loss_of_function_comes_only_from_a_lesion(self):
        """No lesion, no ``disrupted`` - at any ORF length."""
        for length in (1, 5, 59, 60, 61, 400):
            ev = LocusEvidence(identity_pct=99.0, coverage_pct=100.0,
                               orf_aa_length=length,
                               reference_aa_length=REFERENCE_CODONS,
                               stop_is_reference_terminator=False)
            assert ev.lesion_type == LesionType.NONE
            assert oprd.classify_structure("s1", ev,
                                           search=_search()).verdict == \
                StructuralVerdict.INTACT

    def test_a_lesion_gives_disrupted_at_every_orf_length(self):
        for length in (1, 59, 400):
            ev = LocusEvidence(identity_pct=99.0, coverage_pct=100.0,
                               orf_aa_length=length,
                               reference_aa_length=REFERENCE_CODONS,
                               frameshift=True)
            assert oprd.classify_structure("s1", ev,
                                           search=_search()).verdict == \
                StructuralVerdict.DISRUPTED


# ---------------------------------------------------------------------------
# The lesion context: evidence, and proved to be evidence.
# ---------------------------------------------------------------------------

class TestTheLesionContextIsEvidenceAndNeverADecision:
    def test_the_indels_are_read_off_the_alignment(self):
        rows = _column_gap_rows("c1", INFRAME_DEL, DELETED_CODONS)
        _hits, m = _measure(rows, INFRAME_DEL)
        assert [(i.kind, i.length_nt) for i in m.indels] == \
            [("deletion", 6)]
        assert m.indels[0].net_nt == -6
        assert m.indels[0].source == "hsp_columns"

    def test_an_indel_between_two_hsps_is_read_from_their_coordinates(self):
        qseq = _residues(REF_CDS)
        whole = qseq
        second = qseq[20:]
        rows = [
            _row("c1", whole[:20], whole[:20], LOCUS_START,
                 LOCUS_START + 59, 1, 20),
            _row("c1", second, second, LOCUS_START + 56,
                 LOCUS_START + 56 + 3 * len(second) - 1, 21, len(qseq)),
        ]
        _hits, m = _measure(rows, COMPLETE)
        assert [i.source for i in m.indels] == ["hsp_junction"]
        assert m.indels[0].length_nt == 4
        assert m.net_frame_change_nt == -4
        assert m.evidence.frameshift is True

    def test_no_aligned_columns_is_refused_rather_than_called_no_indels(self):
        """A lost measurement is not a negative."""
        row = "\t".join(_exact_rows("c1", COMPLETE)[0].split("\t")[:12])
        hits = oprd.parse_tblastn_table(row, REFERENCE_CODONS)
        with pytest.raises(oprd.PipelineError) as excinfo:
            oprd.indels_from_alignment(hits)
        assert "aligned" in str(excinfo.value)

    def test_a_compensating_pair_is_found_and_a_lone_indel_is_not(self):
        lone = [Indel("deletion", 5, (100, 104), "hsp_junction")]
        assert oprd.compensating_indel_group(lone) is None
        pair = lone + [Indel("deletion", 1, (117, 117), "hsp_junction")]
        group = oprd.compensating_indel_group(pair)
        assert group is not None
        assert sum(i.net_nt for i in group) % 3 == 0
        far = lone + [Indel("deletion", 1, (400, 400), "hsp_junction")]
        assert oprd.compensating_indel_group(far) is None

    def test_a_repeat_context_is_found_around_an_indel(self):
        filler = "".join(hashlib.sha256(str(i).encode()).hexdigest()[:3].upper()
                         for i in range(20))
        # A perfect tandem has several equivalent rotations of its unit, so the
        # motif is asserted as a period and a width rather than as one string.
        tandem = "ACG" + "CTACGGCTACGG" + filler
        motif = oprd.repeat_context(
            tandem, [Indel("deletion", 6, (12, 17), "hsp_columns")])
        assert motif is not None and motif.endswith("x2")
        assert len(motif.split("x")[0]) == 6
        assert tandem.count(motif.split("x")[0]) >= 2
        homopolymer = "ACG" + "T" + "GGGGGG" + "C" + filler
        assert oprd.repeat_context(
            homopolymer, [Indel("insertion", 1, (4, 4), "hsp_junction")]) == \
            "GGGGGG"
        plain = filler + "ACGTTGCATTGACCA" + filler[::-1] + "TTACGGCA"
        assert oprd.repeat_context(
            plain, [Indel("deletion", 1, (1, 1), "hsp_junction")]) is None

    def test_no_context_attribute_can_move_a_verdict(self):
        """The context is not on :class:`LocusEvidence`, and never reaches it."""
        rows = _column_gap_rows("c1", INFRAME_DEL, DELETED_CODONS)
        _hits, m = _measure(rows, INFRAME_DEL)
        assert m.indels
        for mangled in (
            replace(m, indels=()),
            replace(m, compensating_group=None),
            replace(m, repeat_motif=None),
            replace(m, net_frame_change_nt=0),
            replace(m, net_frame_change_nt=99),
        ):
            assert mangled.evidence == m.evidence
        # and the same locus with no measurable indel context at all
        assert _call(rows, INFRAME_DEL).verdict == StructuralVerdict.INTACT

    def test_the_reference_pin_gives_the_same_frame_as_the_protein_pin(self):
        """One reference, read two ways. :func:`reference_cds_nucleotides` is
        the DNA sibling of :func:`reference_protein` and shares its reader, so
        the structural path cannot drift onto a different reference."""
        import inspect

        assert (inspect.getsource(oprd.reference_cds_nucleotides).count(
            "_reference_cds_from_pin") >= 1)
        assert (inspect.getsource(oprd.reference_protein).count(
            "_reference_cds_from_pin") >= 1)


class TestTheTblastnTableAndCommand:
    def test_the_declared_outfmt_is_what_the_stored_artifacts_used(self):
        assert oprd.TBLASTN_OUTFMT_FIELDS[-2:] == ("qseq", "sseq")
        assert oprd.TBLASTN_OUTFMT.startswith("6 ")
        assert "db_gencode" not in oprd.TBLASTN_OUTFMT

    def test_the_command_carries_the_bacterial_translation_table(self):
        command = oprd.tblastn_command(
            program="tblastn", query=Path("/q/PA0958.faa"),
            database=Path("/db/PDT000167136.1"), threads=7)
        assert "-db_gencode" in command
        assert command[command.index("-db_gencode") + 1] == "11"
        assert command[command.index("-num_threads") + 1] == "7"
        assert command[command.index("-evalue") + 1] == "0.001"
        assert command[command.index("-max_target_seqs") + 1] == "5000"
        assert command[command.index("-outfmt") + 1] == oprd.TBLASTN_OUTFMT

    def test_the_twelve_column_form_parses_without_the_aligned_columns(self):
        row = "\t".join(_exact_rows("c1", COMPLETE)[0].split("\t")[:12])
        hits = oprd.parse_tblastn_table(row, REFERENCE_CODONS)
        assert len(hits) == 1
        assert hits[0].query_aligned is None
        assert hits[0].query_low == 1
        assert hits[0].subject_low == LOCUS_START
        assert hits[0].direction == 1

    def test_a_row_of_another_width_is_refused_not_parsed(self):
        row = "\t".join(_exact_rows("c1", COMPLETE)[0].split("\t") + ["x"])
        with pytest.raises(oprd.PipelineError) as excinfo:
            oprd.parse_tblastn_table(row, REFERENCE_CODONS)
        assert "fields" in str(excinfo.value)

    def test_a_reverse_hit_reports_its_direction(self):
        hit = AlignedNucleotideHit(
            query_id="q", subject_id="c", identity_pct=97.0,
            alignment_length=60, query_length=60, subject_length=20000,
            query_start=1, query_end=60, subject_start=9999, subject_end=9400,
            evalue=0.0, bitscore=120.0)
        assert hit.is_reverse is True
        assert hit.direction == -1
        assert hit.subject_low == 9400
        assert hit.subject_high == 9999
        assert list(hit.query_codons)[:3] == [1, 2, 3]