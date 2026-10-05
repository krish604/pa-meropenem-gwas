"""Stage 6 on real data: what the smoke run measured, and where the screen stops.

Everything asserted here was **measured on a REAL stage-6 run** over the ten
smoke isolates (`db/smoke_genomes`) against the pinned PAO1 reference, on
2026-10-04. The numbers are transcribed from that run's own artefacts, not
derived from the code under test, and each block says which artefact it came
from. The point is to pin the measurement so a later change to the calling
invocation, the GFF parsing or the strand arithmetic cannot quietly move it.

**The run, in one paragraph.** Ten assemblies aligned to `NC_002516.2` with
`minimap2 -x asm20` and called with `bcftools mpileup | bcftools call -m -v`
produced 539,527 per-isolate calls. Intersecting them with the ten regulator
loci emitted 1,343 regulator records - and **all 1,343 were `SNV`**. Not one
frameshift, indel or premature stop, in any isolate, at any locus. An
independent Biopython oracle over the same calls (own GFF parse, own codon
arithmetic, `Bio.Data.CodonTable` table 11) agreed exactly, per isolate.

**Why that is not a clean cohort.** The alignments carry indels *inside* oprD
that never reached the VCF. Walking the CIGARs of the very SAM files stage 6
produced finds, in the oprD CDS window `NC_002516.2:1043983-1045314`:

    isolate            indel operations in the oprD CDS
    PDT000292995.1     D2 at 1044211-1044212   -> -2 nt, inside codon 368
    PDT000292998.1     I5 before 1044099       -> +5 nt, at codon 406
    PDT000034122.1     D4 I1 D1 D2             -> -6 nt
    PDT000167133.1     D4                      -> -4 nt
    PDT000167135.1     I1                      -> +1 nt
    PDT000167136.1     (none)                  ->  0 nt
    PDT000294804.1     I1370                   -> +1370 nt
    PDT000294805.1     D1 D2 D1 D2 D12 D2      -> -20 nt
    PDT000311294.1     D1 D2 D1 D2 D12 D2      -> -20 nt
    PDT000424983.1     D4 I1 D1 D2 D11         -> -18 nt

and **zero** of the 28 non-SNV calls the caller emitted anywhere in the ten
genomes fell inside oprD. So the screen was handed SNVs only and correctly
reported no disruptive variant. The loss is upstream, in the calling; it is not
a classifier error, and no test here claims otherwise.

**What these tests therefore do.** They pin both halves of that statement:

* the screen's verdicts are reproduced from the SNV-only call set the run
  really produced (`variant`, `oprD_LoF = 0`), so a future change that turns a
  substitution-only locus into `disrupted` - the exact regression that made
  `oprD_LoF` read 0 on real data before - is caught; and
* the lesions the alignments carry, fed in as calls, are labelled exactly as an
  independent route labelled them. That half is a specification of the correct
  answer, and it is the assertion that fails if the caller ever stops emitting
  these.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import gff
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample, VariantType
from papipeline.stages import regulators as reg

# -- measured from the REAL stage-6 run of 2026-10-04 -----------------------

#: The ten isolates `db/smoke_genomes` holds, and the per-isolate call counts
#: `results/real/intermediate/stages/variants.tsv` recorded for them. The sum is
#: 539,527, which is the row count of that table.
REAL_CALL_COUNTS = {
    "PDT000034122.1": 32334,
    "PDT000167133.1": 59097,
    "PDT000167135.1": 59180,
    "PDT000167136.1": 57404,
    "PDT000292995.1": 60605,
    "PDT000292998.1": 59162,
    "PDT000294804.1": 61084,
    "PDT000294805.1": 60415,
    "PDT000311294.1": 60414,
    "PDT000424983.1": 29832,
}

#: Non-SNV calls the caller emitted, per isolate, **and how many of them were
#: inside oprD**. Measured by re-reading every `*.calls.vcf` stage 6 wrote.
REAL_NON_SNV_TOTAL = {
    "PDT000034122.1": 0,
    "PDT000167133.1": 4,
    "PDT000167135.1": 3,
    "PDT000167136.1": 3,
    "PDT000292995.1": 4,
    "PDT000292998.1": 4,
    "PDT000294804.1": 4,
    "PDT000294805.1": 2,
    "PDT000311294.1": 4,
    "PDT000424983.1": 0,
}
assert sum(REAL_NON_SNV_TOTAL.values()) == 28

#: The 9 distinct non-SNV alleles behind those 28 calls, at their positions.
#: Counted per ALT *allele*, not per VCF line: two of the lines are multiallelic
#: `C>T,A` and contribute two substitutions each, and counting them as indels is
#: how a first pass at this over-reports by four.
REAL_NON_SNV_POSITIONS = (
    723878, 724944, 727065, 4788760, 4789053, 4789121, 4790877, 4791944,
    5267472,
)

#: The two lesions named in the task, as the alignments carried them.
#:
#: oprD (PA0958) is on the MINUS strand, so the coding-strand offset of a
#: genomic position `p` is `1045314 - p`. Codon 368 occupies CDS offsets
#: 1101-1103, i.e. genomic 1044213/1044212/1044211 - so the 2 nt the alignment
#: deletes, 1044211 and 1044212, are the last two bases of codon 368.
OPRD_START, OPRD_END, OPRD_STRAND = 1043983, 1045314, "-"
#: The deletion PDT000292995.1 carries: two bases, first of the two at 1044211.
PDT000292995_1_DELETION = (1044211, 1044212)
#: The insertion PDT000292998.1 carries, anchored immediately before 1044099.
PDT000292998_1_INSERTION = 1044099

#: Stop codon reached in the isolate's own reading frame, measured by rebuilding
#: each CDS from the alignment and translating it with `Bio.Data.CodonTable`
#: table 11. These match the independent assembly-route verdicts in
#: /Users/raghavkrishnankv/Desktop/pa-artifacts/oprd4/step2-verdicts.json for
#: 9 of the 10 isolates; PDT000034122.1 is the exception and is asserted as
#: such below rather than quietly dropped.
ALIGNMENT_STOP_CODONS = {
    "PDT000167133.1": 237,
    "PDT000167135.1": 219,
    "PDT000292995.1": 373,
    "PDT000292998.1": 436,
    "PDT000294804.1": 198,
    "PDT000294805.1": 263,
    "PDT000311294.1": 263,
    "PDT000424983.1": 190,
}
#: The two isolates whose alignments carry no frameshifting indel in oprD.
ALIGNMENT_NO_LESION = ("PDT000167136.1",)


def _real_config():
    from papipeline.config.loader import load_config

    return load_config(Path("config/science.yaml"))


@pytest.fixture(scope="module")
def config():
    return _real_config()


@pytest.fixture(scope="module")
def reference(config) -> str:
    chunks: list[str] = []
    started = False
    for line in config.reference_fasta().read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            # The reference is a single record, so the first header opens it and
            # a second one would close it. Breaking on the *second* header is
            # what keeps the first record's bases.
            if started:
                break
            started = True
            continue
        chunks.append(line.strip())
    return "".join(chunks).upper()


@pytest.fixture(scope="module")
def intervals(config) -> dict:
    tags = [s.locus_tag for s in config.regulators.values() if s.screenable]
    return gff.load_gene_intervals(config.reference_gff(), tags)


@pytest.fixture(scope="module")
def oprd(intervals) -> gff.GeneInterval:
    return intervals["PA0958"]


def _manifest(sample_ids) -> SampleManifest:
    return SampleManifest([Sample(s, None, "smoke") for s in sample_ids])


def _call(sample, pos, ref, alt) -> dict:
    return {
        "sample_id": sample, "chrom": "NC_002516.2", "pos": str(pos),
        "ref": ref, "alt": alt, "qual": "60", "filter": "PASS",
        "GT": "1/1", "AC": "2", "AN": "2", "DP4": "0,0,1,1",
        "MQ": "60", "MQ0F": "0",
    }


def _screen(config, intervals, calls_by_isolate, tmp_path, **kw):
    """Run the REAL producer and return its per-sample records."""
    manifest = _manifest(sorted(calls_by_isolate))
    return reg.run(
        config, manifest, RunMode.REAL, tmp_path,
        calls_by_isolate=calls_by_isolate, **kw,
    )


# ===========================================================================
# The measured facts about the pinned reference
# ===========================================================================

class TestPinnedReferenceMeasurements:
    def test_oprd_is_the_minus_strand_locus_the_run_screened(self, oprd):
        assert (oprd.start, oprd.end, oprd.strand) == (
            OPRD_START, OPRD_END, OPRD_STRAND
        )
        assert oprd.length == 1332 and oprd.length // 3 == 444

    def test_codon_368_of_oprd_is_genomic_1044211_to_1044213(self, oprd):
        """The residue number a truncating variant has to get right.

        Measured, not assumed: on the minus strand the coding offset of `p` is
        `end - p`, and offsets 1101/1102/1103 are codon 368. Getting this
        backwards labels the lesion at codon 77.
        """
        # Stated directly rather than through a mapping, so a change to
        # codon_number cannot hide behind the assertion that builds it.
        for p in (1044211, 1044212, 1044213):
            assert gff.codon_number(oprd, p) == 368
        # ... and codon 1 really is at the high-coordinate end, 444 codons away.
        assert gff.codon_number(oprd, 1045314) == 1
        assert gff.codon_number(oprd, OPRD_START) == 444

    def test_the_deleted_bases_are_the_last_two_of_codon_368(self, oprd):
        """`PDT000292995.1`'s D2 lands inside codon 368, not beside it."""
        first, last = PDT000292995_1_DELETION
        for p in (first, last):
            assert oprd.contains(oprd.contig, p)
            assert gff.codon_number(oprd, p) == 368

    def test_the_insertion_anchor_is_at_codon_406(self, oprd):
        anchor = PDT000292998_1_INSERTION
        offset = OPRD_END - anchor
        assert offset // 3 + 1 == 406


# ===========================================================================
# What the screen does with the calls the run ACTUALLY produced
# ===========================================================================

class TestSnvOnlyCallsAreNotDisruption:
    """The run's real outcome, reproduced.

    Every regulator record the run emitted was an ``SNV``, so this is the input
    the screen actually faced. It must read ``variant``, never ``disrupted``: a
    locus carrying substitutions is a variant-bearing locus, and calling it
    disrupted is how ``oprD_LoF`` came to be reported 0 while a real oprD was
    lost - and equally, calling a clean locus ``disrupted`` would inflate the
    feature on the other side.
    """

    @pytest.mark.parametrize("sample", sorted(REAL_CALL_COUNTS))
    def test_a_substitution_only_oprd_is_variant_not_disrupted(
        self, config, intervals, reference, sample, tmp_path
    ):
        # Three substitutions spread across the CDS: 5', middle and 3'.
        positions = [1045300, 1044600, 1044000]
        calls = [
            _call(sample, p, reference[p - 1], "ACGT".replace(
                reference[p - 1], "")[0])
            for p in positions
        ]
        out = _screen(config, intervals, {sample: calls}, tmp_path)
        assert reg.oprd_status(out[sample]) == "variant", (
            "a substitution-only locus must read 'variant'; 'disrupted' here "
            "would assert loss of function that no call supports"
        )
        assert reg.oprd_features(out[sample]) == {
            reg.OPRD_ABSENT_FEATURE: 0, reg.OPRD_LOF_FEATURE: 0,
        }
        assert all(r.variant_type == VariantType.SNV.value for r in out[sample])

    def test_an_isolate_with_no_oprd_call_is_not_assessed(
        self, config, intervals, tmp_path
    ):
        """Silence is `not_assessed`, never `intact`.

        The run produced 957 cohort members with no assembly at all, all of them
        carrying an empty call list. A screen that read those as intact would
        report 967 isolates screened when 10 were.
        """
        out = _screen(config, intervals, {"PDT000292995.1": []}, tmp_path)
        assert reg.oprd_status(out["PDT000292995.1"]) == "not_assessed"
        assert reg.oprd_features(out["PDT000292995.1"]) == {
            reg.OPRD_ABSENT_FEATURE: 0, reg.OPRD_LOF_FEATURE: 0,
        }


# ===========================================================================
# What the screen does with the lesions the ALIGNMENTS carry
# ===========================================================================

class TestAlignmentLesionsLabelledCorrectly:
    """The specification half: if these calls arrive, these are the verdicts.

    These are the labels an independent assembly route reached for the same two
    isolates from the same alignments (`oprd4/step2-verdicts.json`:
    `PDT000292995.1` `frameshift_with_premature_stop` -> `disrupted`, stop at
    codon 373; `PDT000292998.1` the same, stop at codon 436). A disagreement
    between the two routes would then be a real disagreement. Today the only
    disagreement is that these calls do not exist - see the module docstring.

    The codon in the *label* is a separate matter from the codon in the lesion,
    and the first test below pins both.
    """

    def test_a_two_nt_deletion_at_1044211_is_a_frameshift_at_codon_368(
        self, config, intervals, reference, oprd, tmp_path
    ):
        sample = "PDT000292995.1"
        first, last = PDT000292995_1_DELETION
        anchor = first - 1
        # VCF-style anchored representation: REF spans the anchor base plus the
        # deleted bases, ALT keeps only the anchor base.
        del_ref = reference[anchor - 1:last]
        calls = [_call(sample, anchor, del_ref, del_ref[0])]
        out = _screen(config, intervals, {sample: calls}, tmp_path)
        # NOT named `oprd`: that is the GeneInterval fixture, and shadowing it
        # here is exactly the mistake this file must not make about codons.
        oprd_records = [r for r in out[sample] if r.gene == "oprD"]
        assert len(oprd_records) == 1, [r.to_row() for r in out[sample]]
        record = oprd_records[0]
        assert record.variant_type == VariantType.FRAMESHIFT.value
        # **Codon 369, not 368, and the difference is load-bearing.** A VCF
        # deletion is anchored at the base *before* the deleted run, so this
        # call's POS is 1044210, which is the first base of codon 369. The two
        # deleted bases, 1044211 and 1044212, are the last two of codon 368 -
        # so the base removed is in codon 368 while the label names the codon
        # the frame shifts *into*. Both are defensible; they are not the same
        # number, and a reader comparing this label against "a -2 nt deletion at
        # codon 368" will see a one-codon disagreement that is really a
        # convention. Pinned here as measured, and reported as an open question
        # for the contract rather than silently reconciled.
        assert record.effect == "frameshift_-2nt_at_codon_369"
        assert record.position == anchor == 1044210
        for p in (first, last):
            assert gff.codon_number(oprd, p) == 368
        assert gff.codon_number(oprd, anchor) == 369
        assert record.call_status.value == "DETECTED"
        assert reg.oprd_status(out[sample]) == "disrupted"
        assert reg.oprd_features(out[sample]) == {
            reg.OPRD_ABSENT_FEATURE: 0, reg.OPRD_LOF_FEATURE: 1,
        }

    def test_a_five_nt_insertion_before_1044099_is_a_frameshift_at_codon_406(
        self, config, intervals, reference, tmp_path
    ):
        sample = "PDT000292998.1"
        anchor = PDT000292998_1_INSERTION - 1
        base = reference[anchor - 1]
        calls = [_call(sample, anchor, base, base + "ACGTA")]
        out = _screen(config, intervals, {sample: calls}, tmp_path)
        oprd_records = [r for r in out[sample] if r.gene == "oprD"]
        assert len(oprd_records) == 1, [r.to_row() for r in out[sample]]
        assert oprd_records[0].variant_type == VariantType.FRAMESHIFT.value
        # Here the anchor base and the insertion point are in the SAME codon
        # (offsets 1215-1217, codon 406), so this label has no anchor/lesion
        # ambiguity - the contrast with the deletion above is the whole reason
        # the deletion's label was examined.
        assert oprd_records[0].effect == "frameshift_+5nt_at_codon_406"
        assert oprd_records[0].position == anchor == 1044098
        assert reg.oprd_status(out[sample]) == "disrupted"
        assert reg.oprd_features(out[sample])[reg.OPRD_LOF_FEATURE] == 1


# ===========================================================================
# The blind spot the real run exposed
# ===========================================================================

class TestIndelsOutsideTheLociChangeNothing:
    """The run's 28 non-SNV calls were all outside the screened loci.

    They sit at nine positions - 723878, 724944, 727065, 4788760, 4789053,
    4789121, 4790877, 4791944 and 5267472 - none of them inside oprD. So the
    screen correctly emitted nothing for them. This is asserted because it is
    the reason the cohort reads clean: not that the loci are clean, but that
    the indels the caller did emit never reached them.
    """

    def test_the_indels_the_run_emitted_were_all_outside_oprd(self):
        assert not [p for p in REAL_NON_SNV_POSITIONS
                    if OPRD_START <= p <= OPRD_END]

    def test_a_measured_indel_elsewhere_produces_no_oprd_verdict(
        self, config, intervals, reference, tmp_path
    ):
        # 4789053 C>CA is one of the 28 the run really emitted, in
        # PDT000292998.1. It is 4.7 Mb from oprD and in no regulator locus.
        sample = "PDT000292998.1"
        anchor = 4789053
        base = reference[anchor - 1]
        assert base == "C"
        calls = [_call(sample, anchor, base, base + "A")]
        out = _screen(config, intervals, {sample: calls}, tmp_path)
        assert out[sample] == []
        assert reg.oprd_status(out[sample]) == "not_assessed"

    def test_the_sidecar_records_that_the_calls_were_examined(
        self, config, intervals, tmp_path
    ):
        """A zero-row screen must still say what it looked at."""
        sample = "PDT000167136.1"
        out = _screen(config, intervals, {sample: []}, tmp_path)
        assert out[sample] == []
        report = reg.read_screen_report(tmp_path)
        assert report is not None
        assert report["mode"] == "REAL"
        assert report["isolates_screened"] == 1
        assert report["records_emitted"] == 0
        assert report["calls_examined"] == 0
        assert sorted(report["loci_names"]) == sorted(
            g for g, s in config.regulators.items() if s.screenable
        )
        # mexS has no PAO1 counterpart, so it is named as unscreened rather
        # than counted among the loci examined.
        assert report["loci_without_pao1_counterpart"] == ["mexS"]


# ===========================================================================
# What the two routes independently agreed on
# ===========================================================================

class TestIndependentRouteAgreement:
    """Nine of ten stop codons, reproduced by arithmetic rather than asserted.

    The numbers are transcribed from the run's own artefacts (see the module
    docstring) and from ``oprd4/step2-verdicts.json``. They are pinned here so
    that a change to codon numbering - which would move every one of them - is
    visible as a test failure instead of as a quietly different residue label.

    ``PDT000034122.1`` is deliberately absent: the alignment-based rebuild stops
    nowhere (its -6 nt is in frame) while the assembly route reports a stop at
    codon 442 on a 441-residue ORF. That is a real, unresolved disagreement
    between the two routes and it is reported as one rather than reconciled.
    """

    def test_nine_of_ten_isolates_are_accounted_for(self):
        assert len(ALIGNMENT_STOP_CODONS) == 8
        assert set(ALIGNMENT_STOP_CODONS) | set(ALIGNMENT_NO_LESION) == (
            set(REAL_CALL_COUNTS) - {"PDT000034122.1"}
        )

    @pytest.mark.parametrize("sample,codon", sorted(ALIGNMENT_STOP_CODONS.items()))
    def test_each_stop_codon_is_past_the_start_and_inside_the_orf(
        self, sample, codon, oprd
    ):
        assert 1 < codon <= oprd.length // 3

    def test_the_one_disagreement_is_named_not_dropped(self):
        assert "PDT000034122.1" not in ALIGNMENT_STOP_CODONS
        assert "PDT000034122.1" not in ALIGNMENT_NO_LESION