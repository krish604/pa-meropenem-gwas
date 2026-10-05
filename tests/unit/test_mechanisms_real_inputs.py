"""Stage 5 mechanisms: the variant-allele defect, and a REAL path that reads.

Two separate round-12 obligations, in one file because they are the same
argument: **this stage may only ever lower a claim, and it may only do so on
evidence it actually read.**

Part 1 - `determinant.variant` was never read (TRIAGE class (b) rank 1).
`interpret_determinant` keyed on `determinant.gene or determinant.split("~")[0]`
and looked the gene up. `oprD_V359L` therefore produced byte-identical output to
intact `oprD`: mechanism `reduced_permeability`, evidence `DETECTED`, and the
knowledge table's own note saying *"Presence of oprD indicates an intact
locus"*. One row, two contradictory statements, one of them the thing
`config/mechanisms.tsv:31` explicitly forbids. Round 11 made it live by adding
`--organism`, which is where the 55 POINT rows came from.

Part 2 - there was **no REAL path at all**. No loader, no mode branch; a
standalone or Snakemake invocation had nothing to read. `run` now re-reads its
four inputs from their contracted paths in REAL and refuses, naming the input,
rather than emitting a table whose emptiness reads as "no mechanism found".

The tests here exercise the real code paths with synthetic tables written into
`tmp_path` (ruling R12). No bioinformatics tool is invoked and no fixture is
added to `test_data/`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.errors import StageError
from papipeline.io.tsv import write_tsv
from papipeline.manifest import SampleManifest
from papipeline.models import (
    AmrDeterminant,
    ClaimStatus,
    RegulatorVariant,
    RunMode,
    Sample,
    StructuralCallStatus,
    StructuralVariant,
)
from papipeline.stages import mechanisms as stage_mech

AMR_COLUMNS = [
    "sample_id", "antibiotic", "determinant", "gene", "variant",
    "determinant_type", "mechanism", "evidence_source", "database",
    "database_version", "confidence", "claim_status", "identity_pct",
    "coverage_pct",
]

REGULATOR_COLUMNS = [
    "sample_id", "gene", "variant", "variant_type", "position", "reference",
    "alternate", "effect", "mechanism", "confidence", "evidence_source",
    "call_status",
]

SV_COLUMNS = [
    "sample_id", "variant_id", "variant_type", "position", "affected_gene",
    "size", "evidence", "confidence", "call_status", "mge",
]

ANNOTATION_COLUMNS = [
    "sample_id", "contig_id", "gene_id", "gene_name", "product", "gene_type",
    "start", "end", "strand", "annotation_source",
]


def _determinant(gene="oprD", sample="S1", variant=None, **kwargs):
    determinant = kwargs.pop("determinant", None) or (
        f"{gene}_{variant}" if variant else gene
    )
    defaults = dict(
        sample_id=sample,
        antibiotic="imipenem",
        determinant=determinant,
        gene=gene,
        variant=variant,
        determinant_type="AMR",
        mechanism=None,
        evidence_source="amrfinderplus",
        database="AMRFinderPlus",
        database_version="2026-08-07.1",
        confidence="99.0",
        claim_status=ClaimStatus.DETECTED,
    )
    defaults.update(kwargs)
    return AmrDeterminant(**defaults)


def _manifest(*sample_ids):
    return SampleManifest(
        samples=[Sample(sample_id=s) for s in (sample_ids or ("S1",))]
    )


# ---------------------------------------------------------------------------
# Part 1 - a variant allele is not a gene detection
# ---------------------------------------------------------------------------


class TestAVariantAlleleIsNotAnIntactLocus:
    """The defect was the *presence claim*, not the strength of the claim.

    `ClaimStatus` is ordered SUPPORTED > ASSOCIATED > PREDICTED > DETECTED >
    UNKNOWN (`models.CLAIM_STATUS_PRECEDENCE`), so `PREDICTED` is stronger than
    `DETECTED` and cannot be used to hold a claim down - `ceiling_for(config,
    "oprD", PREDICTED)` returns DETECTED, because DETECTED is oprD's ceiling.
    An allele *was* detected, so DETECTED is true of the row. What was false was
    the row also asserting, in the knowledge table's own words, that presence of
    oprD indicates an intact locus - for a locus the same row calls altered.
    """

    def test_the_variant_row_does_not_carry_the_presence_claim(self, config):
        """The contradiction, asserted absent.

        `config/mechanisms.tsv:31` reads *"Presence of oprD indicates an intact
        locus, not susceptibility. Loss-of-function must come from stage 6 or
        stage 7."* Emitting that beside `oprD_V359L` is the defect.
        """
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD", variant="V359L"), "imipenem"
        )
        assert "Presence of oprD indicates an intact locus" not in call.notes

    def test_the_variant_row_says_what_it_actually_observed(self, config):
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD", variant="V359L"), "imipenem"
        )
        assert "variant allele" in call.notes
        assert "stage 6 or stage 7" in call.notes

    def test_intact_oprd_still_carries_the_presence_claim(self, config):
        """Unchanged behaviour for a real detection - this is not a blanket edit."""
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD"), "imipenem"
        )
        assert "Presence of oprD indicates an intact locus" in call.notes
        assert call.evidence_level is ClaimStatus.DETECTED

    def test_the_two_rows_are_no_longer_identical(self, config):
        """The defect in one assertion."""
        intact = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD"), "imipenem"
        )
        variant = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD", variant="V359L"), "imipenem"
        )
        assert intact.notes != variant.notes
        assert intact.mechanism == variant.mechanism

    def test_the_variant_claim_is_the_tools_own_not_a_substitute(self, config):
        """An allele was detected, so DETECTED stands. Nothing is invented."""
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD", variant="V359L"), "imipenem"
        )
        assert call.evidence_level is ClaimStatus.DETECTED

    def test_a_variant_at_another_gene_keeps_its_own_note(self, config):
        """`nalC` variants are what its knowledge-table note says it records.

        Substituting the allele note for every gene would be a different wrong
        answer: it would discard a note the table asks to be shown.
        """
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="nalC", variant="G71E"), "imipenem"
        )
        assert "nalC loss-of-function variants are reported" in call.notes
        assert call.evidence_level is ClaimStatus.DETECTED
        assert call.mechanism == "efflux"

    def test_the_gene_ceiling_still_applies_to_an_allele(self, config):
        """The central guarantee is untouched by the new branch."""
        call = stage_mech.interpret_determinant(
            config,
            _determinant(
                gene="oprD", variant="V359L", claim_status=ClaimStatus.SUPPORTED
            ),
            "imipenem",
        )
        assert call.evidence_level is ClaimStatus.DETECTED


class TestAStageSixVariantCarriesNoPresenceClaimEither:
    """The same defect one layer down, and it is the LARGER half.

    `interpret_variant` is the stage-6 route - the one the module docstring
    calls the route by which loss-of-function *should* arrive. It appended
    `spec.notes` unconditionally, so every oprD row from the stage 6 screen
    carried `config/mechanisms.tsv`'s *"Presence of oprD indicates an intact
    locus"* beside a description of a variant in oprD. Measured on the real
    screen output for the 10 smoke isolates: **769** oprD rows, 769 of them
    carrying the claim (`round12/cooc/o5_report.json`). Fixing only
    `interpret_determinant` would have left the larger half in place.
    """

    @staticmethod
    def _variant(gene="oprD", sample="S1", variant_type="SNV"):
        return RegulatorVariant(
            sample_id=sample,
            gene=gene,
            variant=f"{gene}:p.V359L",
            variant_type=variant_type,
            position=1076,
            reference="G",
            alternate="T",
            effect="missense",
            mechanism=None,
            confidence=None,
            evidence_source="screen",
            call_status=ClaimStatus.DETECTED,
        )

    def test_an_oprd_variant_does_not_carry_the_presence_claim(self, config):
        call = stage_mech.interpret_variant(config, self._variant(), "imipenem")
        assert call.mechanism == "reduced_permeability"
        assert "Presence of oprD indicates an intact locus" not in call.notes
        assert "stage 6 or stage 7" in call.notes

    def test_a_gene_disruption_still_says_what_it_saw(self, config):
        """The mechanism class is unchanged; only the presence claim is dropped."""
        # RegulatorVariant is frozen, so the type is a parameter rather than an
        # assignment - a frozen model should not be poked with setattr in a test
        # that is not testing setattr.
        call = stage_mech.interpret_variant(
            config,
            self._variant(variant_type="GENE_DISRUPTION"),
            "imipenem",
        )
        assert call.mechanism == "reduced_permeability"
        assert call.evidence_level is ClaimStatus.DETECTED

    def test_a_variant_at_another_gene_keeps_its_own_note(self, config):
        call = stage_mech.interpret_variant(
            config, self._variant(gene="mexZ"), "imipenem"
        )
        assert "variant SNV detected in mexZ" in call.notes
        assert "Promoter and coding variants both relevant" in call.notes

    def test_an_intact_oprd_determinant_still_carries_the_claim(self, config):
        """The counterpart, so the test above cannot pass by a blanket removal."""
        call = stage_mech.interpret_determinant(
            config, _determinant(gene="oprD"), "imipenem"
        )
        assert "Presence of oprD indicates an intact locus" in call.notes


class TestAVariantSuppressesTheIntactLocusCertification:
    def test_annotate_declines_when_told_about_the_allele(self, config):
        """`locus_intact` asserts an intact locus; an altered one is not that.

        Declined rather than emitted at a lower status: there is no
        "intact_but-altered" status, and either of the two real ones would state
        something the variant contradicts.
        """
        assert (
            stage_mech.annotate_intact_locus(
                config, "S1", "oprD", "imipenem", "annotation",
                variant_allele="V359L",
            )
            is None
        )

    def test_and_still_emits_when_there_is_no_allele(self, config):
        call = stage_mech.annotate_intact_locus(
            config, "S1", "oprD", "imipenem", "annotation", variant_allele=None
        )
        assert call is not None
        assert call.mechanism == "locus_intact"

    def test_run_does_not_certify_a_locus_that_carries_a_variant(self, config):
        """The end-to-end symptom: one sample, one table, no contradiction."""
        calls = stage_mech.run(
            config,
            _manifest("S1"),
            RunMode.TEST,
            "imipenem",
            amr_calls={"S1": [_determinant(gene="oprD", variant="V359L")]},
            annotated_genes={"S1": ["oprD", "mexR"]},
        )["S1"]
        assert [c.mechanism for c in calls] == ["reduced_permeability"]
        assert "Presence of oprD indicates an intact locus" not in calls[0].notes

    def test_run_does_certify_an_unmutated_locus(self, config):
        calls = stage_mech.run(
            config,
            _manifest("S1"),
            RunMode.TEST,
            "imipenem",
            amr_calls={"S1": [_determinant(gene="oprD")]},
            annotated_genes={"S1": ["oprD"]},
        )["S1"]
        assert sorted(c.mechanism for c in calls) == [
            "locus_intact",
            "reduced_permeability",
        ]


# ---------------------------------------------------------------------------
# Part 2 - the REAL disk path
# ---------------------------------------------------------------------------


def _write_real_tree(root: Path, *, amr=True, regulators=True, sv=True,
                     annotation=True) -> Path:
    """A complete REAL intermediate tree, or a hole in one of the four inputs."""
    stage_dir = root / "stages"
    stage_dir.mkdir(parents=True, exist_ok=True)

    if amr:
        write_tsv(
            stage_dir / "04_amr.tsv",
            [
                {
                    "sample_id": "S1", "antibiotic": "imipenem",
                    "determinant": "oprD_V359L", "gene": "oprD",
                    "variant": "V359L", "determinant_type": "AMR",
                    "mechanism": "", "evidence_source": "amrfinderplus",
                    "database": "AMRFinderPlus",
                    "database_version": "2026-08-07.1", "confidence": "99.0",
                    "claim_status": "DETECTED", "identity_pct": "99.0",
                    "coverage_pct": "100.0",
                }
            ],
            AMR_COLUMNS,
        )
    if regulators:
        regulators_dir = root / "regulators"
        regulators_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(
            regulators_dir / "regulator_variants.tsv",
            [
                {
                    "sample_id": "S1", "gene": "mexZ", "variant": "c.-12T>C",
                    "variant_type": "SNV", "position": "12",
                    "reference": "T", "alternate": "C", "effect": "",
                    "mechanism": "", "confidence": "",
                    "evidence_source": "screen", "call_status": "DETECTED",
                }
            ],
            REGULATOR_COLUMNS,
        )
        (regulators_dir / "regulator_screen.json").write_text(
            '{"isolates_screened": 1, "isolates_in_manifest": 1}',
            encoding="utf-8",
        )
    if sv:
        write_tsv(
            stage_dir / "07_structural_variants.tsv",
            [
                {
                    "sample_id": "S1", "variant_id": "SV1",
                    "variant_type": "INSERTION", "position": "100",
                    "affected_gene": "mexT", "size": "1200",
                    "evidence": "nucmer", "confidence": "",
                    "call_status": "CONFIRMED", "mge": "",
                }
            ],
            SV_COLUMNS,
        )
    if annotation:
        annotation_dir = root / "annotation"
        annotation_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(
            annotation_dir / "S1.annotation.tsv",
            [
                {
                    "sample_id": "S1", "contig_id": "c1", "gene_id": "b1",
                    "gene_name": "oprD", "product": "porin OprD",
                    "gene_type": "CDS", "start": "1", "end": "900",
                    "strand": "+", "annotation_source": "bakta",
                },
                {
                    "sample_id": "S1", "contig_id": "c1", "gene_id": "b2",
                    "gene_name": "mexR", "product": "repressor MexR",
                    "gene_type": "CDS", "start": "1000", "end": "1500",
                    "strand": "-", "annotation_source": "bakta",
                },
            ],
            ANNOTATION_COLUMNS,
        )
    return root


class TestRealPathsAreTheContractedOnes:
    def test_the_paths_are_read_from_contracts_not_inlined(self, config, tmp_path):
        from papipeline.execution.contracts import internal_table_path, table_path

        paths = stage_mech.real_input_paths(config, tmp_path, RunMode.REAL)
        stage_dir = tmp_path / "stages"
        assert paths["amr_calls"] == str(table_path(stage_dir, "amr"))
        assert paths["structural_variants"] == str(
            internal_table_path(stage_dir, "structural_variants")
        )
        assert paths["regulator_variants"].endswith(
            "regulators/regulator_variants.tsv"
        )
        assert paths["annotated_genes"].endswith("/annotation")

    def test_an_unknown_antibiotic_still_refuses_before_anything_is_read(
        self, config, tmp_path
    ):
        with pytest.raises(Exception):
            stage_mech.run(
                config, _manifest("S1"), RunMode.REAL, "nonsense",
                amr_calls={}, intermediate_root=tmp_path,
            )


class TestRealRunReadsFromDisk:
    def test_a_complete_tree_produces_calls_with_no_argument_supplied(
        self, config, tmp_path
    ):
        """The point of the loader: `amr_calls=None` and the table still fills.

        Before the round-12 change there was no loader at all, so a standalone
        invocation had nothing to read and this call could only have returned
        an empty dict - which reads as "no mechanism found".
        """
        root = _write_real_tree(tmp_path)
        calls = stage_mech.run(
            config,
            _manifest("S1"),
            RunMode.REAL,
            "imipenem",
            amr_calls={},
            intermediate_root=root,
        )["S1"]
        assert calls, "a complete tree produced no mechanism calls"
        assert {c.gene for c in calls} == {"oprD", "mexZ", "mexT"}

    def test_the_oprd_variant_reaches_the_table_through_the_loader(self, config, tmp_path):
        """The round-11 POINT rows are the case this fix exists for.

        55 POINT rows across the 10 smoke isolates, oprD among them, so a REAL
        run now emits `reduced_permeability` at `PREDICTED` for `oprD_V359L`
        where it emitted `DETECTED` before - and no `locus_intact` beside it.
        """
        root = _write_real_tree(tmp_path)
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.REAL, "imipenem",
            amr_calls={}, intermediate_root=root,
        )["S1"]
        oprd = [c for c in calls if c.gene == "oprD"]
        assert [c.determinant for c in oprd] == ["oprD_V359L"]
        assert oprd[0].mechanism == "reduced_permeability"
        # The sample carries oprD_V359L and an oprD annotation; the table must
        # not also certify the locus intact.
        assert not [c for c in oprd if c.mechanism == "locus_intact"]

    def test_the_structural_and_regulator_inputs_are_both_consumed(self, config, tmp_path):
        root = _write_real_tree(tmp_path)
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.REAL, "imipenem",
            amr_calls={}, intermediate_root=root,
        )["S1"]
        determinants = {c.determinant for c in calls}
        assert "SV1" in determinants
        assert "c.-12T>C" in determinants  # MechanismCall.determinant is the variant

    def test_the_annotation_input_is_what_supplies_gene_names(self, config, tmp_path):
        """`02_annotation_summary.tsv` carries counts, so this is the only source.

        Worth a test because the alternative - a summary table - looks like it
        should work and yields nothing.
        """
        root = _write_real_tree(tmp_path)
        loaded = stage_mech.load_annotated_gene_names(
            root / "annotation", _manifest("S1")
        )
        assert loaded == {"S1": ["mexR", "oprD"]}

    def test_amr_calls_is_required_rather_than_omissible(self, config):
        """The seam guard, restated at the stage.

        `convergence` and `cooccurrence` make their equivalent optional. This one
        does not: omitting it would reach the stage body with no stage-4 result
        and produce an empty table that reads as "no mechanism found". Pass `{}`
        to mean "there is nothing".
        """
        import inspect

        parameters = inspect.signature(stage_mech.run).parameters
        assert parameters["amr_calls"].default is inspect.Parameter.empty

    def test_the_caller_arguments_are_ignored_in_real(self, config, tmp_path):
        """A stale in-memory argument must not win over what is on disk."""
        root = _write_real_tree(tmp_path)
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.REAL, "imipenem",
            amr_calls={"S1": [_determinant(gene="mexR", sample="S1")]},
            intermediate_root=root,
        )["S1"]
        assert "mexR" not in {c.determinant for c in calls}
        assert "oprD_V359L" in {c.determinant for c in calls}


class TestRealRunRefusesRatherThanEmittingAnEmptyTable:
    @pytest.mark.parametrize(
        "hole", ["amr", "regulators", "sv", "annotation"],
    )
    def test_an_absent_input_is_named(self, config, tmp_path, hole):
        root = _write_real_tree(tmp_path, **{hole: False})
        with pytest.raises(StageError) as excinfo:
            stage_mech.run(
                config, _manifest("S1"), RunMode.REAL, "imipenem",
                amr_calls={}, intermediate_root=root,
            )
        message = str(excinfo.value)
        assert "mechanisms" in message.lower()
        # The reader is told which stage input and which stage produces it.
        expected = {
            "amr": "amr_calls",
            "regulators": "regulator_variants",
            "sv": "structural_variants",
            "annotation": "annotated_genes",
        }[hole]
        assert expected in message

    def test_a_present_but_empty_table_is_a_result_not_a_refusal(
        self, config, tmp_path
    ):
        """The distinction the check is drawn on.

        An absent table means the producing stage did not run; a table that is
        present with no rows means it ran and found nothing. Refusing on the
        second throws away a real answer - "this cohort carries no structural
        variant" is a finding about the cohort. `regulators.
        load_regulator_variants` draws the same line with the provenance
        sidecar; there is none here, so presence is what is available.
        """
        root = _write_real_tree(tmp_path)
        # Header only, no rows: the shape a producer writes when it finds
        # nothing.
        for name in ("04_amr.tsv", "07_structural_variants.tsv"):
            path = root / "stages" / name
            header = path.read_text(encoding="utf-8").splitlines()[0]
            path.write_text(header + "\n", encoding="utf-8")
        # The regulator and annotation inputs are still populated, so the point
        # is that the run PROCEEDS rather than refusing on the two empty ones.
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.REAL, "imipenem",
            amr_calls={}, intermediate_root=root,
        )["S1"]
        determinants = {c.determinant for c in calls}
        assert "oprD_V359L" not in determinants   # its table is empty
        assert "SV1" not in determinants          # its table is empty too
        assert "c.-12T>C" in determinants          # the regulator table survived

    def test_a_present_but_empty_annotation_directory_is_also_a_result(
        self, config, tmp_path
    ):
        root = _write_real_tree(tmp_path)
        for path in (root / "annotation").glob("*.tsv"):
            path.unlink()
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.REAL, "imipenem",
            amr_calls={}, intermediate_root=root,
        )
        assert {c.determinant for c in calls["S1"]} >= {"oprD_V359L"}
        assert not [c for c in calls["S1"] if c.mechanism == "locus_intact"]

    def test_no_intermediate_root_is_a_refusal_not_a_fallback(self, config):
        with pytest.raises(StageError) as excinfo:
            stage_mech.run(
                config, _manifest("S1"), RunMode.REAL, "imipenem",
                amr_calls={"S1": [_determinant()]}, intermediate_root=None,
            )
        assert "intermediate_root" in str(excinfo.value)

    def test_stub_still_raises_rather_than_computing(self, config):
        """The fabrication guard is untouched."""
        with pytest.raises(NotImplementedError):
            stage_mech.run(
                config, _manifest("S1"), RunMode.STUB, "imipenem",
                amr_calls={"S1": [_determinant()]},
            )

    def test_test_mode_still_computes_from_arguments(self, config):
        """The mode branch is a branch, not a replacement."""
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.TEST, "imipenem",
            amr_calls={"S1": [_determinant(gene="oprD")]},
        )
        assert calls["S1"]


class TestUnmappedInputIsCountedNotGuessed:
    def test_a_gene_outside_the_knowledge_table_is_not_invented(self, config):
        calls = stage_mech.run(
            config, _manifest("S1"), RunMode.TEST, "imipenem",
            amr_calls={"S1": [_determinant(gene="notAKnownGene")]},
        )
        assert calls["S1"] == []
