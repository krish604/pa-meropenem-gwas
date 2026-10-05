"""Per-locus coverage: `not_assessed` must never be reported as `no variant`.

The producer screens a locus by intersecting calls with it. An intersection is
silent about loci the aligner never reached, so an unaligned locus and a clean
locus produce byte-identical output. This file pins the mechanism that
distinguishes them, and pins it against the alternatives that look equivalent
and are not.

Three things are pinned here:

* the **threshold boundary** - at what covered fraction a locus stops being
  assessed, and that the boundary itself is inclusive;
* **supplementary-included counting** - 8 of the 10 real smoke isolates reach
  oprD only on a supplementary alignment, so excluding those records would turn
  the check into a constant false alarm;
* **the distinction itself** - a below-threshold locus is `not_assessed` with a
  NAMED reason, its records are suppressed, and no code path can turn that into
  "there is none".

The measured coverage figures are real: they come from stage 6's own alignments
in `results/real/intermediate/variants/*.sorted.bam` against the pinned PAO1
GFF, and are reproduced here as the constants the sensitivity claim rests on.
"""

from pathlib import Path

import pytest

from papipeline.adapters import gff
from papipeline.manifest import Sample, SampleManifest
from papipeline.models import RegulatorVariant, VariantType
from papipeline.stages import regulators as R
from papipeline.stages.regulators import (
    COVERAGE_NOT_ASSESSED_REASON,
    MIN_LOCUS_COVERAGE_FRACTION,
    NOT_ASSESSED,
    covered_reference_positions,
    locus_assessment,
    locus_coverage_fraction,
    locus_coverage_rows,
    not_assessed_reason,
    oprd_features,
    oprd_status,
    screen_calls,
)

# ---------------------------------------------------------------------------
# Measured facts, pinned.
#
# oprD is NC_002516.2:1043983-1045314 on the minus strand, 1,332 nt, read from
# the pinned GFF. ampR is NC_002516.2:4592990-4593880, 891 nt. Both verified by
# grepping the GFF, not copied from any prior report.
# ---------------------------------------------------------------------------
OPRD_START, OPRD_END, OPRD_LEN = 1043983, 1045314, 1332
AMPR_START, AMPR_END, AMPR_LEN = 4592990, 4593880, 891

#: PDT000034122.1's single alignment over ampR starts at 4593022, i.e. 32 nt
#: into the CDS, so 859 of 891 positions are covered. This is the ONLY pair out
#: of 100 (isolate, locus) pairs on the 10-isolate smoke set that is not at
#: exactly 1.0000, and therefore it is the single observation the whole
#: threshold's behaviour rests on. Named here so tightening the constant is a
#: decision about this measurement and not an unexplained change in output.
MEASURED_AMPR_COVERAGE_PDT000034122 = round(859 / AMPR_LEN, 4)
MEASURED_FULL_COVERAGE = 1.0


class _Record:
    """Minimal stand-in exposing what the coverage helper reads.

    Deliberately not a ``pysam.AlignedSegment``: the helper's contract is the
    four attributes below, and a test that built a real BAM would be testing
    pysam as much as the helper.
    """

    def __init__(self, start0, end, *, supplementary=False, secondary=False):
        self.reference_start = start0
        self.reference_end = end
        self.is_supplementary = supplementary
        self.is_secondary = secondary


# ===========================================================================
# The named constant and its default
# ===========================================================================


def test_threshold_is_a_named_float_with_an_explicit_default():
    """The threshold is a constant, not a literal typed at a call site."""
    assert isinstance(MIN_LOCUS_COVERAGE_FRACTION, float)
    assert MIN_LOCUS_COVERAGE_FRACTION == 0.90
    # 0 < threshold <= 1: a threshold above 1 would make every locus
    # unassessable, and 0 or below would assert nothing.
    assert 0.0 < MIN_LOCUS_COVERAGE_FRACTION <= 1.0


def test_reason_is_a_named_string_and_is_not_empty():
    assert COVERAGE_NOT_ASSESSED_REASON == "locus_coverage_below_minimum"
    assert isinstance(COVERAGE_NOT_ASSESSED_REASON, str)
    assert COVERAGE_NOT_ASSESSED_REASON.strip() == COVERAGE_NOT_ASSESSED_REASON


# ===========================================================================
# Threshold boundary
# ===========================================================================


@pytest.mark.parametrize(
    "coverage,expected_not_assessed",
    [
        (1.0, False),
        # The measured ampR pair. 0.9641 is ABOVE the 0.90 default, so the
        # shipped threshold assesses it; 0.8999 is below and is not assessed.
        (0.9641, False),
        (MIN_LOCUS_COVERAGE_FRACTION, False),   # exactly at: assessed
        (0.9001, False),
        (0.8999, True),
        (0.0, True),
    ],
)
def test_boundary_is_inclusive_at_the_threshold(coverage, expected_not_assessed):
    """`coverage >= threshold` passes; one ulp below fails.

    Inclusive on purpose: the documented default has to be able to include a
    fully covered locus without being written as 0.9000001.
    """
    got = not_assessed_reason(coverage)
    assert (got is not None) is expected_not_assessed
    if expected_not_assessed:
        assert got == COVERAGE_NOT_ASSESSED_REASON


def test_explicit_threshold_overrides_the_default():
    """The threshold is a parameter, so its sensitivity is checkable directly."""
    # The shipped default assesses the measured pair.
    assert not_assessed_reason(0.9641, 0.90) is None
    # One ulp above the measured value flips it.
    assert not_assessed_reason(0.9641, 0.9642) == COVERAGE_NOT_ASSESSED_REASON
    # At 0.9641 exactly, the measured value meets the floor.
    assert not_assessed_reason(0.9641, 0.9641) is None
    # And an unambiguously partial locus is caught by the default.
    assert not_assessed_reason(0.5, 0.90) == COVERAGE_NOT_ASSESSED_REASON


def test_unknown_coverage_is_not_a_pass_and_not_a_failure():
    """`coverage=None` yields no reason - and that is the safe direction.

    It does not fabricate a passing value, because asserting a measurement
    nobody made is the error this whole mechanism exists to prevent. The
    absence is visible in `locus_coverage_rows` as `coverage=None`.
    """
    assert not_assessed_reason(None) is None
    rows = locus_coverage_rows(["S1"], ["oprD"], None)
    assert rows == [{
        "sample_id": "S1", "gene": "oprD", "coverage": None,
        "min_locus_coverage": MIN_LOCUS_COVERAGE_FRACTION,
        "not_assessed": False, "reason": None,
    }]


def test_measured_sensitivity_of_the_default_threshold():
    """The sensitivity of the default, stated from the measured distribution.

    100 (isolate, locus) pairs on the 10-isolate smoke set: 99 at exactly
    1.0000 and one at 0.9641. So with `coverage >= threshold`:
      any threshold in (0, 0.9641]   changes 0 classifications
      any threshold in (0.9641, 1.00] changes exactly 1 - PDT000034122.1 x ampR
    A threshold of 1.00 changes one pair and no more; it does not cascade.
    """
    observed = [MEASURED_FULL_COVERAGE] * 99 + [MEASURED_AMPR_COVERAGE_PDT000034122]
    assert len(observed) == 100
    lowest = min(observed)

    for threshold in (0.10, 0.50, 0.90, lowest):
        flipped = [c for c in observed if not_assessed_reason(c, threshold)]
        assert flipped == [], f"threshold {threshold} flipped {flipped}"
    for threshold in (0.9642, 0.98, 1.00):
        flipped = [c for c in observed if not_assessed_reason(c, threshold)]
        assert flipped == [MEASURED_AMPR_COVERAGE_PDT000034122], threshold


# ===========================================================================
# Supplementary-included counting
# ===========================================================================


def test_supplementary_alignments_are_counted():
    """A locus covered only by a supplementary alignment is covered.

    Measured: 8 of the 10 smoke isolates reach oprD on a record with SAM flag
    2048. Excluding those would report most of the cohort's loci as unseen.
    """
    only_supplementary = [_Record(OPRD_START - 1, OPRD_END, supplementary=True)]
    assert covered_reference_positions(
        only_supplementary, OPRD_START, OPRD_END) == set(
            range(OPRD_START, OPRD_END + 1))
    assert locus_coverage_fraction(
        only_supplementary, OPRD_START, OPRD_END) == MEASURED_FULL_COVERAGE


def test_excluding_supplementary_is_what_would_break_it():
    """The exclusion is demonstrated, not asserted, and it is not the default."""
    only_supplementary = [_Record(OPRD_START - 1, OPRD_END, supplementary=True)]
    assert locus_coverage_fraction(
        only_supplementary, OPRD_START, OPRD_END,
        include_supplementary=False) == 0.0
    assert not_assessed_reason(locus_coverage_fraction(
        only_supplementary, OPRD_START, OPRD_END,
        include_supplementary=False)) == COVERAGE_NOT_ASSESSED_REASON


def test_secondary_alignments_are_also_counted_by_default():
    secondary = [_Record(OPRD_START - 1, OPRD_END, secondary=True)]
    assert locus_coverage_fraction(secondary, OPRD_START, OPRD_END) == 1.0
    assert locus_coverage_fraction(
        secondary, OPRD_START, OPRD_END, include_secondary=False) == 0.0


def test_union_of_primary_and_supplementary_is_the_union():
    """Coverage is the union of records, so a split alignment does not lose
    the part the primary record stops short of."""
    primary = _Record(OPRD_START - 1, OPRD_START + 99)
    supplementary = _Record(OPRD_START + 99, OPRD_END, supplementary=True)
    assert locus_coverage_fraction(
        [primary, supplementary], OPRD_START, OPRD_END) == MEASURED_FULL_COVERAGE
    # And each on its own is short.
    assert locus_coverage_fraction(
        [primary], OPRD_START, OPRD_END) == 100 / OPRD_LEN


def test_measured_alignment_boundary_inside_ampR():
    """The one measured pair below 1.0000, reproduced from its CIGAR geometry.

    PDT000034122.1's record over ampR begins at 4593022 with a hard clip, so
    positions 4592990-4593021 - the first 32 nt of the CDS, about 10 codons -
    are covered by nothing. 859/891 = 0.9641.
    """
    record = _Record(4593022 - 1, AMPR_END)          # 1-based 4593022..4593880
    covered = covered_reference_positions([record], AMPR_START, AMPR_END)
    assert len(covered) == 859
    assert min(covered) == 4593022 and max(covered) == AMPR_END
    assert locus_coverage_fraction(
        [record], AMPR_START, AMPR_END) == pytest.approx(
        MEASURED_AMPR_COVERAGE_PDT000034122, abs=5e-5)
    assert round(locus_coverage_fraction(
        [record], AMPR_START, AMPR_END), 4) == MEASURED_AMPR_COVERAGE_PDT000034122
    # ABOVE the shipped default, so this pair IS assessed as shipped. The test
    # says so explicitly, so nobody tightens the constant to 0.98 and silently
    # changes this isolate's ampR verdict without reading why it would.
    assert not_assessed_reason(
        locus_coverage_fraction([record], AMPR_START, AMPR_END)) is None
    # Raising the threshold above the measured value is what excludes it.
    assert not_assessed_reason(
        locus_coverage_fraction([record], AMPR_START, AMPR_END), 0.9642
    ) == COVERAGE_NOT_ASSESSED_REASON


def test_inverted_interval_reads_as_zero_coverage_not_an_exception():
    assert locus_coverage_fraction([], 100, 50) == 0.0
    assert covered_reference_positions([], 100, 50) == set()


# ===========================================================================
# The distinction: `not_assessed` with a named reason, never "no variant"
# ===========================================================================


def _record(sample="S1", gene="oprD", vtype=VariantType.SNV, position=1044100):
    return RegulatorVariant(
        sample_id=sample, gene=gene, variant=f"{gene}_v", variant_type=vtype,
        position=position, reference="A", alternate="T", effect="x",
        mechanism="m", confidence=None, evidence_source="test",
        call_status="candidate",
    )


#: Every label this module can produce for a locus. `not_assessed` is the only
#: one permitted for a below-threshold locus, and none of these may stand in
#: for it.
STATUS_LABELS = ("absent", "disrupted", "variant", NOT_ASSESSED)


def test_below_threshold_is_not_assessed_with_the_named_reason():
    assessment = locus_assessment("oprD", [], coverage=0.5)
    assert assessment["status"] == NOT_ASSESSED
    assert assessment["reason"] == COVERAGE_NOT_ASSESSED_REASON
    assert assessment["coverage"] == 0.5
    # The named reason, compared as a constant, so a reworded message cannot
    # pass as the same reason.
    assert assessment["reason"] is R.COVERAGE_NOT_ASSESSED_REASON
    assert assessment["min_locus_coverage"] == MIN_LOCUS_COVERAGE_FRACTION


def test_below_threshold_overrides_a_disruptive_record():
    """A call at an under-covered locus is not evidence about that locus.

    The record says `disrupted`; the coverage says the locus was not read. The
    coverage wins, because a truncating call in a locus the aligner only partly
    read still asserts the locus's sequence, which is the unsupported claim.
    """
    records = [_record(vtype=VariantType.FRAMESHIFT)]
    assert R.gene_status(records)["oprD"] == "disrupted"
    assessment = locus_assessment("oprD", records, coverage=0.5)
    assert assessment["status"] == NOT_ASSESSED
    assert assessment["reason"] == COVERAGE_NOT_ASSESSED_REASON
    # The record is still counted, so the suppression is visible rather than
    # looking like the call never happened.
    assert assessment["records"] == 1


def test_below_threshold_is_never_reported_as_no_variant():
    """The core assertion, spelled out over every status this module emits."""
    for records in ([], [_record()], [_record(vtype=VariantType.FRAMESHIFT)],
                    [_record(vtype=VariantType.GENE_ABSENCE)]):
        status = R.gene_status(records, coverage={"oprD": 0.5}).get("oprD")
        assert status == NOT_ASSESSED, status
        assert status != "variant"
        assert status != "absent"
        assert status != "disrupted"
        # `oprd_status` collapses to the same three declared values plus
        # not_assessed, and must land on not_assessed.
        assert oprd_status(records, coverage=0.5) == NOT_ASSESSED
        # Both features read 0 - never 1 on under-covered evidence.
        assert oprd_features(records, coverage=0.5) == {
            R.OPRD_ABSENT_FEATURE: 0, R.OPRD_LOF_FEATURE: 0}


def test_absent_is_not_reachable_through_the_coverage_path():
    """A GENE_ABSENCE call under-covered does not become `absent`.

    `absent` is a claim that the locus is gone from the genome. Coverage
    evidence cannot support that claim, so the coverage path must not be able
    to produce it, from any record type.
    """
    for records in ([], [_record(vtype=VariantType.GENE_ABSENCE)]):
        got = R.gene_status(records, coverage={"oprD": 0.1}).get("oprD")
        assert got == NOT_ASSESSED
        assert got != "absent"


def test_at_or_above_threshold_behaves_exactly_as_before():
    """No behaviour change when coverage is sufficient or not supplied."""
    records = [_record(vtype=VariantType.FRAMESHIFT)]
    baseline = R.gene_status(records)
    assert R.gene_status(records, coverage={"oprD": 1.0}) == baseline
    assert R.gene_status(records, coverage={"oprD": 0.90}) == baseline
    assert R.gene_status(records) == baseline
    assert oprd_features(records) == oprd_features(records, coverage=1.0)


def test_a_gene_absent_from_the_coverage_mapping_is_still_assessed():
    """Coverage evidence for one locus must not silence another."""
    records = [_record(vtype=VariantType.FRAMESHIFT)]
    status = R.gene_status(records, coverage={"ampR": 0.5})
    assert status["oprD"] == "disrupted"       # not mentioned -> assessed
    assert status["ampR"] == NOT_ASSESSED


# ===========================================================================
# The producer: records suppressed, accounting carries the reason
# ===========================================================================


class _Spec:
    """Just enough of RegulatorSpec for `gff.classify` and the record build."""

    def __init__(self):
        self.evidence_level = "candidate"
        self.variant_classes = ("SNV", "INDEL", "FRAMESHIFT", "PREMATURE_STOP",
                                "GENE_ABSENCE", "UNKNOWN")
        self.screenable = True
        self.locus_tag = "PA0958"
        self.mechanism = "test_mechanism"
        self.gene = "oprD"


def _manifest():
    return SampleManifest([Sample(sample_id=s) for s in ("S1", "S2")])


def _screen(calls, coverage, threshold=MIN_LOCUS_COVERAGE_FRACTION):
    interval = gff.GeneInterval(
        locus_tag="PA0958", contig="NC_002516.2", start=OPRD_START,
        end=OPRD_END, strand="-", name="oprD")
    reference = "".join(
        [("ACGT"[i % 4]) for i in range(OPRD_END + 10)])
    contigs = {"NC_002516.2": reference}
    return screen_calls(
        calls, _manifest(),
        intervals=[interval], gene_by_tag={"PA0958": "oprD"},
        specs={"oprD": _Spec()}, contigs=contigs,
        locus_coverage=coverage, min_locus_coverage=threshold,
    )


def _calls(sample, position, ref="A", alt="T"):
    return [{"chrom": "NC_002516.2", "pos": str(position), "ref": ref, "alt": alt}]


def test_records_are_suppressed_below_threshold_and_accounted():
    calls = {"S1": _calls("S1", 1044100), "S2": _calls("S2", 1044100)}
    covered = {("S1", "oprD"): 1.0, ("S2", "oprD"): 0.5}
    records, accounting = _screen(calls, covered)

    assert [r.sample_id for r in records] == ["S1"]
    assert accounting["records_suppressed_below_coverage"] == 1
    entry = accounting["loci_not_assessed"]["S2:oprD"]
    assert entry["reason"] == COVERAGE_NOT_ASSESSED_REASON
    assert entry["coverage"] == 0.5
    assert entry["min_locus_coverage"] == MIN_LOCUS_COVERAGE_FRACTION
    assert entry["suppressed_records"] == 1
    assert accounting["coverage_evidence"] == "supplied"
    assert accounting["min_locus_coverage"] == MIN_LOCUS_COVERAGE_FRACTION


def test_no_coverage_evidence_leaves_the_old_behaviour_and_says_so():
    """Absent input: unchanged behaviour, and the key says `absent`."""
    calls = {"S1": _calls("S1", 1044100), "S2": _calls("S2", 1044100)}
    records, accounting = _screen(calls, None)
    assert len(records) == 2
    assert "loci_not_assessed" not in accounting
    assert "records_suppressed_below_coverage" not in accounting
    assert "coverage_evidence" not in accounting


def test_full_coverage_emits_an_empty_not_assessed_map_not_a_missing_one():
    calls = {"S1": _calls("S1", 1044100), "S2": _calls("S2", 1044100)}
    covered = {("S1", "oprD"): 1.0, ("S2", "oprD"): 1.0}
    records, accounting = _screen(calls, covered)
    assert len(records) == 2
    assert accounting["loci_not_assessed"] == {}
    assert accounting["records_suppressed_below_coverage"] == 0


def test_threshold_parameter_changes_the_verdict_on_the_same_evidence():
    """Sensitivity is a parameter change, observable without editing code."""
    calls = {"S1": _calls("S1", 1044100), "S2": _calls("S2", 1044100)}
    covered = {("S1", "oprD"): 1.0, ("S2", "oprD"): 0.9641}
    for threshold, expect_s2_emitted in ((0.90, True), (0.9641, True),
                                         (0.9642, False), (1.00, False)):
        records, accounting = _screen(calls, covered, threshold=threshold)
        got_s2 = any(r.sample_id == "S2" for r in records)
        assert got_s2 is expect_s2_emitted, threshold
        assert ("S2:oprD" in accounting["loci_not_assessed"]) is (
            not expect_s2_emitted), threshold


def test_coverage_rows_expose_every_pair_including_the_unassessed_ones():
    rows = locus_coverage_rows(
        ["S1", "S2"], ["oprD", "ampR"],
        {("S1", "oprD"): 1.0, ("S2", "oprD"): 0.5, ("S1", "ampR"): 0.96},
    )
    assert len(rows) == 4
    by_key = {(r["sample_id"], r["gene"]): r for r in rows}
    assert by_key[("S2", "oprD")]["not_assessed"] is True
    assert by_key[("S2", "oprD")]["reason"] == COVERAGE_NOT_ASSESSED_REASON
    assert by_key[("S1", "oprD")]["not_assessed"] is False
    # A pair with no evidence at all is still a row, with coverage None. It is
    # not silently dropped, because a dropped row reads as a passing one.
    assert by_key[("S2", "ampR")]["coverage"] is None
    assert by_key[("S2", "ampR")]["not_assessed"] is False


def test_every_status_label_is_from_the_declared_vocabulary():
    """Guards against a future label sneaking in beside `not_assessed`."""
    observed = set()
    for coverage in (None, 1.0, 0.9641, 0.5, 0.0):
        for records in ([], [_record()],
                        [_record(vtype=VariantType.FRAMESHIFT)],
                        [_record(vtype=VariantType.GENE_ABSENCE)]):
            got = R.gene_status(records, coverage={"oprD": coverage}).get("oprD")
            if got is not None:
                observed.add(got)
            observed.add(oprd_status(records, coverage=coverage))
    assert observed <= set(STATUS_LABELS), observed
    assert NOT_ASSESSED in observed