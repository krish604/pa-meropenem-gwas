"""The stage-6 regulator screen's REAL producer.

Three things are under test, and they are tested for different reasons.

**The intersection and its refusals.** The screen is the only thing standing
between a cohort of called variants and a claim about oprD, mexR and the rest,
so "which calls does it emit, and what does it refuse rather than emit" is the
scientific content. The refusals matter more than the happy path: every one of
them exists because the alternative was a table that looked well formed and said
something false.

**Empty is not not-run.** A screen that runs and finds nothing, and a screen that
never ran, both leave a table with no rows. Downstream that difference is the
difference between "this cohort is clean" and "nothing was ever looked at", so
the two must not be confusable. The mechanism is the ``regulator_screen.json``
sidecar, and both branches are asserted.

**Against the real pinned reference, not a fixture.** ``PA0958`` (oprD) is on the
**minus** strand, and both bugs this file guards against were strand bugs that a
plus-strand fixture cannot see: an unreachable ``PREMATURE_STOP`` past a gene's
first codon, and codon numbering running backwards. So the locus coordinates,
strand and sequence come from ``db/reference/GCF_000006765.1``, and the expected
numbers were derived from that file rather than from the code under test.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.adapters import gff
from papipeline.errors import DataContractError
from papipeline.manifest import SampleManifest
from papipeline.models import ClaimStatus, RunMode, Sample, VariantType
from papipeline.stages import regulators as reg

# -- coordinates measured from the pinned reference, not from the code ---------
# PA0958 (oprD): gene span == CDS span, 1043983..1045314, MINUS strand, 1332 nt,
# 444 codons. Verified by dumping both feature types out of the real GFF.
OPRD_START, OPRD_END, OPRD_STRAND = 1043983, 1045314, "-"
OPRD_NT = OPRD_END - OPRD_START + 1
assert OPRD_NT == 1332 and OPRD_NT // 3 == 444

# PA4600 (nfxB): PLUS strand, so the same assertions run on the other strand.
NFXB_START, NFXB_END, NFXB_STRAND = 5155561, 5156124, "+"

#: Coding-strand position (1-based) <-> genome coordinate, per strand.
def _coding_to_genome(interval: gff.GeneInterval, coding_pos: int) -> int:
    """Genome coordinate of 1-based base ``coding_pos`` on the coding strand."""
    if interval.strand == "+":
        return interval.start + coding_pos - 1
    return interval.end - coding_pos + 1


def _genome_to_coding(interval: gff.GeneInterval, genome_pos: int) -> int:
    """Inverse of :func:`_coding_to_genome`."""
    if interval.strand == "+":
        return genome_pos - interval.start + 1
    return interval.end - genome_pos + 1


def _genome_bases(reference: str, start: int, length: int) -> str:
    """Plus-strand reference bases, 1-based inclusive, uppercased."""
    return reference[start - 1: start - 1 + length].upper()


@pytest.fixture(scope="module")
def reference() -> str:
    """The pinned PAO1 sequence, read once for this module."""
    config = _real_config()
    fasta = config.reference_fasta()
    name: str | None = None
    chunks: list[str] = []
    for line in fasta.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    return "".join(chunks).upper()


def _real_config():
    from papipeline.config.loader import load_config

    return load_config(Path("config/science.yaml"))


@pytest.fixture(scope="module")
def oprd(reference) -> gff.GeneInterval:
    intervals = gff.load_gene_intervals(_real_config().reference_gff(), ["PA0958"])
    return intervals["PA0958"]


def _two_sample_manifest() -> SampleManifest:
    return SampleManifest(
        [Sample("TEST_A_01", None, "test"), Sample("TEST_A_02", None, "test")]
    )


def _call(sample, contig, pos, ref, alt) -> dict:
    """One row in PER_ISOLATE_COLUMNS shape."""
    return {
        "sample_id": sample, "chrom": contig, "pos": str(pos),
        "ref": ref, "alt": alt, "qual": "60", "filter": "PASS",
        "GT": "1/1", "AC": "1", "AN": "1", "DP4": "0,0,30,29",
        "MQ": "60", "MQ0F": "0",
    }


# ===========================================================================
# The intersection
# ===========================================================================

class TestIntersection:
    def test_a_frameshift_in_oprd_is_emitted_with_the_right_codon(
        self, config, reference, oprd
    ):
        """A -2 nt deletion at coding codon 368 of oprD.

        The residue number is the point. oprD is on the minus strand, so a
        counter anchored at the gene's low coordinate calls codon 368 "codon 77"
        and codon 1 "codon 444" - and the truncation site is exactly the residue
        a reader has to be able to look up.
        """
        # Coding codon 368 occupies coding bases 1102-1104, which on the minus
        # strand are genome 1044213, 1044212, 1044211. Delete two of them.
        first = _coding_to_genome(oprd, 1102)
        ref = _genome_bases(reference, first, 3)
        alt = ref[0]
        assert (len(alt) - len(ref)) == -2

        records, accounting = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, first, ref, alt)]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )

        assert len(records) == 1
        record = records[0]
        assert record.gene == "oprD"
        assert record.variant_type == VariantType.FRAMESHIFT.value
        assert record.effect == "frameshift_-2nt_at_codon_368"
        assert record.position == first
        assert record.reference == ref and record.alternate == alt
        assert record.mechanism == "reduced_permeability"
        assert record.evidence_source == reg.SCREEN_EVIDENCE_SOURCE
        assert record.call_status is config.regulator("oprD").evidence_level
        assert accounting["records_emitted"] == 1
        assert accounting["calls_outside_screened_loci"] == 0

    def test_a_plus_5nt_insertion_is_a_frameshift_too(
        self, config, reference, oprd
    ):
        """The brief's second expectation: +5 nt in PDT000292998.1.

        5 is not a multiple of 3, so an insertion of it shifts the frame. The
        assembly route called that isolate intact, which is a real disagreement
        and is recorded as such in the final report - it is not something this
        test can settle, because it has no variant calls to feed on.

        VCF puts POS on the base *before* the insertion point and carries that
        base in REF, so the codon named in `effect` is POS's codon. The codon is
        therefore asserted against an independently computed
        :func:`gff.codon_number` rather than a hand-counted constant; the one
        hand-counted codon in this file is the headline 368 above.
        """
        pos = _coding_to_genome(oprd, 600)  # last base of coding codon 200
        ref = _genome_bases(reference, pos, 1)
        alt = ref + "ACGTA"
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, pos, ref, alt)]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )
        assert len(records) == 1
        assert records[0].variant_type == VariantType.FRAMESHIFT.value
        assert records[0].effect == (
            f"frameshift_+5nt_at_codon_{gff.codon_number(oprd, pos)}"
        )
        assert gff.codon_number(oprd, pos) == 200

    def test_an_inframe_indel_is_not_a_frameshift(
        self, config, reference, oprd
    ):
        pos = _coding_to_genome(oprd, 301)  # first base of coding codon 101
        ref = _genome_bases(reference, pos, 1)
        alt = ref + "ACG"  # +3 nt: the frame is preserved
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, pos, ref, alt)]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )
        assert len(records) == 1
        assert records[0].variant_type == VariantType.INDEL.value
        assert records[0].effect == "inframe_indel_+3nt_at_codon_101"

    def test_an_inframe_indel_does_not_count_as_loss_of_function(
        self, config, reference, oprd
    ):
        """`gene_status` routes through DISRUPTIVE_TYPES, so an INDEL is 'variant'.

        This is the distinction the two-stage model rests on, and it is the
        reason the class has to be right: an in-frame indel reaching
        `oprD_LoF` would manufacture resistance for a gene that is still
        translated.
        """
        pos = _coding_to_genome(oprd, 301)
        ref = _genome_bases(reference, pos, 1)
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, pos, ref, ref + "ACG")]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )
        assert reg.oprd_status(records) == "variant"
        assert reg.oprd_features(records) == {"oprD_absent": 0, "oprD_LoF": 0}

    def test_a_frameshift_does_count_as_loss_of_function(
        self, config, reference, oprd
    ):
        pos = _coding_to_genome(oprd, 1102)
        ref = _genome_bases(reference, pos, 3)
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, pos, ref, ref[0])]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )
        assert reg.oprd_status(records) == "disrupted"
        assert reg.oprd_features(records) == {"oprD_absent": 0, "oprD_LoF": 1}

    def test_a_stop_gaining_substitution_is_found_on_the_minus_strand(
        self, config, reference, oprd
    ):
        """PREMATURE_STOP past a gene's first codon, on the minus strand.

        Both bugs guarded here were invisible to a plus-strand, codon-1 fixture:
        the within-codon index was the offset from the gene *start*, so the stop
        check was skipped for every variant past the first three bases; and the
        VCF allele was substituted into an already-reverse-complemented codon
        without being complemented first. Together they made PREMATURE_STOP
        unreachable at oprD, mexR, ampD and ampR - four of the ten loci - while
        every existing test passed, because the one case they exercised was
        codon 1 on the plus strand, where the wrong expression lands on the
        right answer.
        """
        assert oprd.strand == "-", "this test is only meaningful on the minus strand"
        found = []
        for number in range(2, 445):
            first_coding = (number - 1) * 3 + 1
            codon = "".join(
                reference[_coding_to_genome(oprd, first_coding + k) - 1].translate(
                    str.maketrans("ACGT", "TGCA")
                )
                for k in range(3)
            )
            if codon in gff.STOP_CODONS:
                # oprD's codon 444 is its own stop. Substituting it is a no-op,
                # and classify correctly answers UNKNOWN, so there is nothing to
                # find here.
                continue
            target = {"AA": "TAA", "AG": "TAG", "GA": "TGA"}.get(codon[1:])
            if target is None:
                continue
            genome = _coding_to_genome(oprd, first_coding)
            vcf_ref = reference[genome - 1]
            vcf_alt = gff.reverse_complement(target[0])
            records, _ = reg.screen_calls(
                {"TEST_A_01": [_call("TEST_A_01", oprd.contig, genome, vcf_ref, vcf_alt)]},
                _two_sample_manifest(),
                intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
                specs={"oprD": config.regulator("oprD")},
                contigs={oprd.contig: reference},
            )
            assert len(records) == 1, f"codon {number} produced no record"
            if records[0].variant_type == VariantType.PREMATURE_STOP.value:
                found.append((number, records[0].effect))

        assert found, "no stop-gaining substitution detected anywhere in oprD"
        # Deep codons must be found: the old gate `0 <= offset < 3` allowed
        # only the first codon, so codon 296 among others was reported as SNV.
        assert any(number > 100 for number, _ in found), (
            f"only shallow codons detected: {sorted(n for n, _ in found)[:5]}"
        )
        for number, effect in found:
            assert f"at_codon_{number}_" in effect, (
                f"codon {number} labelled {effect!r}: the minus-strand codon "
                f"number is running backwards"
            )

    def test_a_call_outside_every_screened_locus_is_not_emitted(
        self, config, reference, oprd
    ):
        records, accounting = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", oprd.contig, 10, "A", "T")]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")},
            contigs={oprd.contig: reference},
        )
        assert records == []
        assert accounting["calls_outside_screened_loci"] == 1

    def test_a_position_in_an_overlap_is_reported_for_both_loci(
        self, reference
    ):
        """Dropping one would under-report the study's central mechanisms."""
        a = gff.GeneInterval("PA_A", "C", 100, 500, "+", name="oprD")
        b = gff.GeneInterval("PA_B", "C", 400, 900, "+", name="ampR")
        specs = {
            "oprD": _spec("oprD", "reduced_permeability"),
            "ampR": _spec("ampR", "efflux"),
        }
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", "C", 450, "A", "T")]},
            _two_sample_manifest(),
            intervals=[a, b], gene_by_tag={"PA_A": "oprD", "PA_B": "ampR"},
            specs=specs, contigs={"C": "A" * 1000},
        )
        assert {r.gene for r in records} == {"oprD", "ampR"}

    def test_the_gene_symbol_comes_from_config_not_the_gff_name(
        self, reference
    ):
        """mexZ and dacB carry their locus tag in the GFF's Name=, not a symbol.

        ``gff.assign`` returns ``interval.name or interval.locus_tag``, so
        trusting it puts "PA2020" in the gene column - which then fails to join
        ``config.regulator("mexZ")`` and gets logged as "not on the regulator
        screen list" for a locus that is on it. The mapping is taken from
        ``config/regulators.tsv`` instead, which is the authority for the symbol.
        """
        interval = gff.GeneInterval("PA2020", "C", 100, 400, "+", name="PA2020")
        records, _ = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", "C", 200, "A", "T")]},
            _two_sample_manifest(),
            intervals=[interval], gene_by_tag={"PA2020": "mexZ"},
            specs={"mexZ": _spec("mexZ", "efflux")},
            contigs={"C": "A" * 1000},
        )
        assert [r.gene for r in records] == ["mexZ"]

    def test_a_class_outside_screen_for_is_counted_not_silently_dropped(
        self, reference
    ):
        """Every locus lists GENE_ABSENCE; nothing can emit it from per-site calls.

        A silently dropped class reads as a screened-and-clean locus, so it is
        counted per gene and carried into the sidecar.
        """
        interval = gff.GeneInterval("PA_A", "C", 100, 400, "+", name="oprD")
        records, accounting = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", "C", 200, "A", "T")]},
            _two_sample_manifest(),
            intervals=[interval], gene_by_tag={"PA_A": "oprD"},
            specs={"oprD": _spec("oprD", "reduced_permeability", classes=("FRAMESHIFT",))},
            contigs={"C": "A" * 1000},
        )
        assert records == []
        assert accounting["calls_out_of_declared_scope"] == {"oprD:SNV": 1}

    def test_a_draft_scaffold_contig_is_counted_not_assigned(
        self, reference, oprd
    ):
        """The coordinate-frame bug this whole design exists to prevent.

        Stage 6 calls are PAO1 coordinates. A call on an assembly contig must
        never resolve against a PAO1 gene, and must not vanish either.
        """
        records, accounting = reg.screen_calls(
            {"TEST_A_01": [_call("TEST_A_01", "JFJU01000001.1", OPRD_START + 5, "A", "T")]},
            _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": _spec("oprD", "reduced_permeability")},
            contigs={oprd.contig: reference},
        )
        assert records == []
        assert accounting["calls_on_unknown_contig"] == {"JFJU01000001.1": 1}

    def test_an_unparseable_position_is_counted_not_crashed_on(
        self, reference, oprd
    ):
        bad = _call("TEST_A_01", oprd.contig, OPRD_START + 5, "A", "T")
        bad["pos"] = "not-a-number"
        records, accounting = reg.screen_calls(
            {"TEST_A_01": [bad]}, _two_sample_manifest(),
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": _spec("oprD", "reduced_permeability")},
            contigs={oprd.contig: reference},
        )
        assert records == []
        assert accounting["calls_unparseable_position"] == 1

    def test_output_is_deterministically_ordered(
        self, config, reference, oprd
    ):
        manifest = SampleManifest(
            [Sample("TEST_A_02", None, "test"), Sample("TEST_A_01", None, "test")]
        )
        a = _call("TEST_A_01", oprd.contig, _coding_to_genome(oprd, 300), "A", "T")
        b = _call("TEST_A_01", oprd.contig, _coding_to_genome(oprd, 100), "C", "G")
        calls = {
            "TEST_A_01": [a, b],
            "TEST_A_02": [_call("TEST_A_02", oprd.contig, _coding_to_genome(oprd, 200), "G", "A")],
        }
        first, _ = reg.screen_calls(
            calls, manifest, intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")}, contigs={oprd.contig: reference},
        )
        second, _ = reg.screen_calls(
            {k: list(v) for k, v in reversed(list(calls.items()))}, manifest,
            intervals=[oprd], gene_by_tag={"PA0958": "oprD"},
            specs={"oprD": config.regulator("oprD")}, contigs={oprd.contig: reference},
        )
        assert [(r.sample_id, r.position) for r in first] == [
            (r.sample_id, r.position) for r in second
        ]
        assert [r.position for r in first if r.sample_id == "TEST_A_01"] == sorted(
            r.position for r in first if r.sample_id == "TEST_A_01"
        )


def _spec(gene, mechanism, classes=("SNV", "INDEL", "FRAMESHIFT", "PREMATURE_STOP")):
    from papipeline.config.loader import RegulatorSpec

    return RegulatorSpec(
        gene=gene, mechanism=mechanism, gene_role="permeability",
        screen_for=tuple(classes), promoter_screen=False,
        reference_length=None, evidence_level=ClaimStatus.DETECTED,
        locus_tag="PA_A",
    )


# ===========================================================================
# Refusals
# ===========================================================================

class TestRefusals:
    """Each of these used to be a silently empty table.

    A refusal names what is missing, so the operator knows what to do. An empty
    table names nothing and is read as a cohort carrying no regulator variant.
    """

    def _manifest(self):
        return SampleManifest([Sample("TEST_A_01", None, "test")])

    def test_no_calls_is_refused_by_name(self, config, tmp_path):
        with pytest.raises(DataContractError) as excinfo:
            reg.run(config, self._manifest(), RunMode.REAL, tmp_path)
        message = str(excinfo.value)
        assert "calls_by_isolate" in message
        assert "variant_calls" in message, "must say what the caller should pass"
        assert not (tmp_path / "regulators" / "regulator_variants.tsv").exists()

    def test_a_manifest_sample_with_no_calls_is_refused(
        self, config, tmp_path
    ):
        """Rule 5: a sample-ID mismatch is a hard failure, never a silent drop."""
        manifest = SampleManifest(
            [Sample("TEST_A_01", None, "test"), Sample("TEST_A_02", None, "test")]
        )
        with pytest.raises(DataContractError) as excinfo:
            reg.run(
                config, manifest, RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": []},
            )
        message = str(excinfo.value)
        assert "TEST_A_02" in message
        assert "never screened" in message
        assert not (tmp_path / "regulators" / "regulator_variants.tsv").exists()

    def test_a_call_from_a_sample_the_manifest_never_names_is_refused(
        self, config, tmp_path
    ):
        with pytest.raises(DataContractError) as excinfo:
            reg.run(
                config, self._manifest(), RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": [], "GHOST_1": []},
            )
        assert "GHOST_1" in str(excinfo.value)
        assert not (tmp_path / "regulators" / "regulator_variants.tsv").exists()

    def test_a_missing_reference_is_refused_by_name(self, config, tmp_path):
        with pytest.raises(DataContractError) as excinfo:
            reg.run(
                config, self._manifest(), RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": []},
                reference_fasta=tmp_path / "absent.fna",
            )
        message = str(excinfo.value)
        assert "absent.fna" in message
        assert "coordinate" in message.lower()

    def test_a_missing_gff_is_refused_by_name(self, config, tmp_path, monkeypatch):
        monkeypatch.setattr(
            type(config), "reference_gff", lambda self: tmp_path / "absent.gff"
        )
        # The adapter raises plain PipelineError here, not DataContractError;
        # asserted against the type it actually raises rather than the one that
        # would read better.
        with pytest.raises(gff.PipelineError) as excinfo:
            reg.run(
                config, self._manifest(), RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": []},
            )
        assert "absent.gff" in str(excinfo.value)

    def test_no_screenable_locus_is_refused_rather_than_writing_nothing(
        self, config, tmp_path
    ):
        """Would otherwise write a zero-row table: a fabricated cohort-wide negative."""
        import dataclasses

        specs = {
            gene: dataclasses.replace(spec, locus_tag=None)
            for gene, spec in config.regulators.items()
        }
        config = dataclasses.replace(config, regulators=specs)
        with pytest.raises(DataContractError) as excinfo:
            reg.run(
                config, self._manifest(), RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": []},
            )
        assert "locus_tag" in str(excinfo.value)
        assert not (tmp_path / "regulators" / "regulator_variants.tsv").exists()

    def test_a_locus_tag_absent_from_the_gff_is_refused(
        self, config, tmp_path
    ):
        """The adapter's own refusal, reached through the producer."""
        import dataclasses

        broken = dict(config.regulators)
        for gene, spec in broken.items():
            if spec.screenable:
                broken[gene] = dataclasses.replace(spec, locus_tag="PA9999999")
                break
        config = dataclasses.replace(config, regulators=broken)
        # Plain PipelineError from the adapter, not DataContractError.
        with pytest.raises(gff.PipelineError) as excinfo:
            reg.run(
                config, self._manifest(), RunMode.REAL, tmp_path,
                calls_by_isolate={"TEST_A_01": []},
            )
        assert "PA9999999" in str(excinfo.value)


# ===========================================================================
# Empty vs not-run
# ===========================================================================

class TestEmptyIsNotNotRun:
    """The distinction the whole sidecar mechanism exists for.

    Two states leave a header-only table: a screen that ran and found nothing,
    and a screen that never ran. Only one of them is a result about the cohort.
    """

    def _write_header_only(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(reg.REGULATOR_COLUMNS) + "\n", encoding="utf-8")
        return path

    def test_a_header_only_table_with_no_sidecar_is_refused(self, config, tmp_path):
        """Its emptiness is evidence of nothing, so it must not read as a negative."""
        table = self._write_header_only(
            tmp_path / "regulators" / "regulator_variants.tsv"
        )
        with pytest.raises(DataContractError) as excinfo:
            reg.load_regulator_variants(config, table, intermediate_root=tmp_path)
        message = str(excinfo.value)
        assert "empty_table_without_provenance" in message
        assert "never ran" in message

    def test_a_header_only_table_with_a_sidecar_is_reported_as_a_result(
        self, config, tmp_path
    ):
        table = self._write_header_only(
            tmp_path / "regulators" / "regulator_variants.tsv"
        )
        reg._write_screen_report(
            tmp_path / "regulators" / reg.SCREEN_REPORT_NAME,
            {
                "mode": "REAL", "isolates_in_manifest": 10, "isolates_screened": 10,
                "loci_screened": 10, "records_emitted": 0,
                "calls_outside_screened_loci": 41234,
            },
        )
        with pytest.raises(DataContractError) as excinfo:
            reg.load_regulator_variants(config, table, intermediate_root=tmp_path)
        message = str(excinfo.value)
        assert "screen_ran_and_found_nothing" in message
        assert "10 of 10 isolates" in message
        assert "41234" in message

    def test_the_real_producer_writes_a_table_and_a_sidecar_when_it_finds_nothing(
        self, config, tmp_path
    ):
        """End to end: a clean cohort is a recorded positive, not a silent file."""
        manifest = SampleManifest([Sample("TEST_A_01", None, "test")])
        grouped = reg.run(
            config, manifest, RunMode.REAL, tmp_path,
            calls_by_isolate={"TEST_A_01": []},
        )
        assert grouped == {"TEST_A_01": []}

        table = tmp_path / "regulators" / "regulator_variants.tsv"
        assert table.is_file()
        body = [
            line for line in table.read_text(encoding="utf-8").splitlines()
            if not line.startswith("#")
        ]
        assert len(body) == 1, "header only: the screen ran and found nothing"
        assert body[0].split("\t") == list(reg.REGULATOR_COLUMNS)

        report = json.loads(
            (tmp_path / "regulators" / reg.SCREEN_REPORT_NAME).read_text("utf-8")
        )
        assert report["records_emitted"] == 0
        assert report["isolates_screened"] == 1
        assert report["isolates_in_manifest"] == 1
        assert report["loci_screened"] == 10
        # mexS has no PAO1 counterpart; recorded as unscreened, not dropped.
        assert report["loci_without_pao1_counterpart"] == ["mexS"]
        assert "mexZ" in report["loci_names"]

    def test_a_producer_run_is_readable_back_by_the_loader(self, config, tmp_path):
        """What the producer writes must satisfy the contract the reader enforces."""
        manifest = SampleManifest([Sample("TEST_A_01", None, "test")])
        oprd = gff.load_gene_intervals(config.reference_gff(), ["PA0958"])["PA0958"]
        reg.run(
            config, manifest, RunMode.REAL, tmp_path,
            calls_by_isolate={
                "TEST_A_01": [
                    _call("TEST_A_01", oprd.contig, _coding_to_genome(oprd, 1102),
                          "A", "T")
                ]
            },
        )
        records = reg.load_regulator_variants(
            config,
            tmp_path / "regulators" / "regulator_variants.tsv",
            intermediate_root=tmp_path,
        )
        assert len(records) == 1
        assert records[0].gene == "oprD"
        assert reg.oprd_status(records) in ("disrupted", "variant")

    def test_the_table_carries_provenance_comments_the_reader_skips(
        self, config, tmp_path
    ):
        manifest = SampleManifest([Sample("TEST_A_01", None, "test")])
        reg.run(
            config, manifest, RunMode.REAL, tmp_path,
            calls_by_isolate={"TEST_A_01": []},
        )
        text = (tmp_path / "regulators" / "regulator_variants.tsv").read_text("utf-8")
        comments = [l for l in text.splitlines() if l.startswith("#")]
        assert any("mode=REAL" in c for c in comments)
        assert any(reg.SCREEN_REPORT_NAME in c for c in comments)


# ===========================================================================
# The strand bugs, pinned at their source
# ===========================================================================

class TestCodonArithmeticOnBothStrands:
    """`adapters/gff` is imported by the producer now, so its arithmetic is ours.

    Both of these were latent while the module had no caller.
    """

    def test_codon_number_is_one_based_at_the_five_prime_end(self, oprd):
        assert gff.codon_number(oprd, _coding_to_genome(oprd, 1)) == 1
        assert gff.codon_number(oprd, _coding_to_genome(oprd, 3)) == 1
        assert gff.codon_number(oprd, _coding_to_genome(oprd, 4)) == 2

    def test_codon_number_368_is_368_and_not_77_on_the_minus_strand(self, oprd):
        """Counting from interval.start inverts every codon on a minus-strand gene."""
        genome = _coding_to_genome(oprd, (368 - 1) * 3 + 1)
        assert gff.codon_number(oprd, genome) == 368
        assert (genome - oprd.start) // 3 + 1 == 77, "the naive expression, kept honest"

    def test_codon_number_agrees_on_the_plus_strand_too(self):
        interval = gff.GeneInterval("PA_N", "C", 1000, 2000, "+")
        assert gff.codon_number(interval, 1000) == 1
        assert gff.codon_number(interval, 1002) == 1
        assert gff.codon_number(interval, 1003) == 2

    def test_codon_index_is_the_base_within_the_codon_not_within_the_gene(self):
        plus = gff.GeneInterval("A", "C", 100, 200, "+")
        minus = gff.GeneInterval("A", "C", 100, 200, "-")
        assert [gff.codon_index(plus, p) for p in (100, 101, 102, 103)] == [0, 1, 2, 0]
        # Mirrored, because codon_at reverse-complements on the minus strand.
        assert [gff.codon_index(minus, p) for p in (100, 101, 102, 103)] == [2, 1, 0, 2]

    def test_stop_gaining_is_reachable_past_the_first_codon_on_both_strands(
        self, reference
    ):
        """The exhaustive check, small enough to keep in the suite.

        Every single-base substitution at every codon of one plus-strand and one
        minus-strand locus. The invariant: `classify` reports PREMATURE_STOP if
        and only if the substituted codon is a stop.
        """
        stops = gff.STOP_CODONS
        for tag, gene in (("PA0958", "oprD"), ("PA4600", "nfxB")):
            interval = gff.load_gene_intervals(
                _real_config().reference_gff(), [tag]
            )[tag]
            coding = reference[interval.start - 1: interval.end]
            if interval.strand == "-":
                coding = gff.reverse_complement(coding)
            checked = 0
            for number in range(1, len(coding) // 3 + 1):
                codon = coding[(number - 1) * 3:(number - 1) * 3 + 3]
                for offset in range(3):
                    for base in "ACGT":
                        if base == codon[offset]:
                            continue
                        mutated = codon[:offset] + base + codon[offset + 1:]
                        coding_pos = (number - 1) * 3 + offset + 1
                        genome = _coding_to_genome(interval, coding_pos)
                        # The reference base is already the plus-strand base, so
                        # it is the VCF REF unchanged. Only the ALT is built from
                        # the coding strand and needs complementing.
                        observed = reference[genome - 1]
                        allele = gff.reverse_complement(base) if interval.strand == "-" else base
                        call = gff.VariantCall("S", interval.contig, genome, observed, allele)
                        got = gff.classify(call, reference, interval)
                        checked += 1
                        assert got == (
                            VariantType.PREMATURE_STOP.value if mutated in stops
                            else VariantType.SNV.value
                        ), (
                            f"{gene} codon {number} base {offset + 1} "
                            f"{codon}->{mutated} at {genome}: got {got}"
                        )
            # 3 bases x 3 alternatives x every codon: 1,692 for nfxB (188
            # codons), 3,996 for oprD (444). Asserted per locus so a silently
            # shortened sweep cannot pass as a clean run.
            assert checked == (len(coding) // 3) * 9, (
                f"{gene}: checked {checked}, expected "
                f"{(len(coding) // 3) * 9}"
            )