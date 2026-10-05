"""Tests for the AMR parser and mechanism mapping (stages 4 and 5).

These encode the pipeline's central scientific guarantee: a detected
determinant is reported as DETECTED and nothing more, and the presence of a
gene is never reported as reduced permeability.
See docs/scientific_rules.md rules 1, 2 and 3.
"""

from __future__ import annotations

import pytest

from papipeline.errors import UnknownGeneError
from papipeline.knowledge import (
    ceiling_for,
    mechanism_map_for_antibiotic,
    unmapped_genes,
)
from papipeline.models import (
    AmrDeterminant,
    ClaimStatus,
    RegulatorVariant,
    StructuralCallStatus,
    StructuralVariant,
    VariantType,
)
from papipeline.stages import amr as stage_amr
from papipeline.stages import mechanisms as stage_mech

AMR_HEADER = (
    "sample_id\tantibiotic\tdeterminant\tgene\tvariant\tdeterminant_type\t"
    "mechanism\tevidence_source\tdatabase\tdatabase_version\tconfidence\n"
)


def _amr_table(rows: str) -> str:
    return AMR_HEADER + rows


def _determinant(gene="blaNDM-1", sample="S1", **kwargs):
    defaults = dict(
        sample_id=sample,
        antibiotic="imipenem",
        determinant=gene,
        gene=gene,
        variant=None,
        determinant_type="acquired",
        mechanism=None,
        evidence_source="amrfinderplus",
        database="AMRFinderPlus",
        database_version="2024-01",
        confidence="99.0",
    )
    defaults.update(kwargs)
    return AmrDeterminant(**defaults)


class TestBaseGene:
    def test_strips_allele_suffix(self):
        assert stage_amr.base_gene("blaOXA-48~SYNTHETIC_ALLELE") == "blaOXA-48"

    def test_no_suffix(self):
        assert stage_amr.base_gene("oprD") == "oprD"

    def test_empty(self):
        assert stage_amr.base_gene("") == ""


class TestParseAmrfinderStdout:
    STDOUT = (
        "Protein identifier\tGene symbol\tSequence Name\tScope\tElement Type\t"
        "Element Subtype\tClass\tSubclass\tMethod\tTarget\tIdentity\t"
        "Resistance mechanism\tCoverage\n"
        "contig_1_00001\tblaNDM-1\tblaNDM-1\tcore\tAMR\tbeta-lactamase\t"
        "B\tMetallo-beta-lactamase\tAMRFinderPlus\tBETA-LACTAM\t99.5\t.\t100.0\n"
    )

    def test_parses_one_determinant(self):
        results = stage_amr.parse_amrfinder_stdout(
            self.STDOUT, "S1", "imipenem", "AMRFinderPlus", "2024-01"
        )
        assert len(results) == 1
        assert results[0].gene == "blaNDM-1"

    def test_records_database_provenance(self):
        results = stage_amr.parse_amrfinder_stdout(
            self.STDOUT, "S1", "imipenem", "AMRFinderPlus", "2024-01"
        )
        assert results[0].database == "AMRFinderPlus"
        assert results[0].database_version == "2024-01"

    def test_status_is_detected_never_more(self):
        results = stage_amr.parse_amrfinder_stdout(
            self.STDOUT, "S1", "imipenem", "AMRFinderPlus", "2024-01"
        )
        assert results[0].claim_status is ClaimStatus.DETECTED

    def test_identity_parsed(self):
        results = stage_amr.parse_amrfinder_stdout(
            self.STDOUT, "S1", "imipenem", "AMRFinderPlus", "2024-01"
        )
        assert results[0].identity_pct == pytest.approx(99.5)

    def test_empty_output(self):
        assert (
            stage_amr.parse_amrfinder_stdout("", "S1", "imipenem", "DB", "v1") == []
        )


class TestLoadAmrTable:
    def test_known_gene_gets_mechanism_and_antibiotic(self, config, write_tsv_file):
        path = write_tsv_file(
            "amr.tsv",
            _amr_table("S1\timipenem\tblaNDM-1\tblaNDM-1\t.\tAMR\t.\tamrfinderplus\tDB\tv1\t99\n"),
        )
        results = stage_amr.load_amr_table(
            path, "imipenem", mechanism_map=mechanism_map_for_antibiotic(config, "imipenem")
        )
        assert results[0].mechanism == "acquired_determinant"
        assert results[0].antibiotic == "imipenem"

    def test_unknown_gene_becomes_unmapped(self, config, write_tsv_file):
        """An unmapped determinant is not attributed to imipenem."""
        path = write_tsv_file(
            "amr.tsv",
            _amr_table("S1\timipenem\tblaFAKE-1\tblaFAKE-1\t.\tAMR\t.\tamrfinderplus\tDB\tv1\t99\n"),
        )
        results = stage_amr.load_amr_table(
            path, "imipenem", mechanism_map=mechanism_map_for_antibiotic(config, "imipenem")
        )
        assert results[0].antibiotic == stage_amr.UNMAPPED_ANTIBIOTIC
        assert results[0].mechanism is None

    def test_unmapped_is_not_fabricated(self, config, write_tsv_file):
        path = write_tsv_file(
            "amr.tsv",
            _amr_table("S1\timipenem\tblaFAKE-1\tblaFAKE-1\t.\tAMR\t.\tamrfinderplus\tDB\tv1\t99\n"),
        )
        results = stage_amr.load_amr_table(
            path, "imipenem", mechanism_map=mechanism_map_for_antibiotic(config, "imipenem")
        )
        assert results[0].claim_status is ClaimStatus.DETECTED

    def test_required_columns_enforced(self, config, write_tsv_file):
        path = write_tsv_file("amr.tsv", "sample_id\tdeterminant\nS1\tblaNDM-1\n")
        with pytest.raises(Exception):
            stage_amr.load_amr_table(path, "imipenem")

    def test_loads_synthetic_fixtures(self, config, intermediate_root):
        results = stage_amr.load_amr_table(
            intermediate_root / "amr" / "amr_determinants.tsv",
            "imipenem",
            mechanism_map=mechanism_map_for_antibiotic(config, "imipenem"),
        )
        assert results
        assert all(r.database_version for r in results)


class TestKnowledgeLookups:
    def test_known_gene_resolves(self, config):
        assert config.mechanism_for_gene("oprD").mechanism == "reduced_permeability"

    def test_allele_suffix_resolves(self, config):
        assert config.mechanism_for_gene("blaOXA-48~ALLELE").gene == "blaOXA-48"

    def test_unknown_gene_raises(self, config):
        """An unknown gene is surfaced, not silently dropped."""
        with pytest.raises(UnknownGeneError, match="not present"):
            config.mechanism_for_gene("geneThatDoesNotExist")

    def test_unknown_genes_helper(self, config):
        assert unmapped_genes(config, ["oprD", "mexR", "nope"]) == ["nope"]

    def test_efflux_genes_are_classified(self, config):
        from papipeline.knowledge import efflux_regulator_genes

        efflux = efflux_regulator_genes(config)
        assert "mexR" in efflux
        assert "nfxB" in efflux
        assert "oprD" not in efflux

    @pytest.mark.parametrize(
        "gene", ["oprD", "mexR", "nalC", "nalD", "mexZ", "nfxB", "mexT", "mexS", "ampD", "ampR", "dacB"]
    )
    def test_every_required_locus_is_configured(self, config, gene):
        assert gene in config.mechanisms
        assert gene in config.regulators

    @pytest.mark.parametrize(
        "gene", ["blaOXA-48", "blaNDM-1", "blaVIM-1", "blaKPC-2", "blaTEM-1"]
    )
    def test_acquired_determinants_are_configured(self, config, gene):
        assert config.mechanisms[gene].mechanism_class == "acquired_determinant"


class TestCeilingFor:
    def test_cannot_exceed_ceiling(self, config):
        """A detection may not be reported as an association."""
        assert (
            ceiling_for(config, "oprD", ClaimStatus.SUPPORTED) is ClaimStatus.DETECTED
        )
        assert (
            ceiling_for(config, "oprD", ClaimStatus.ASSOCIATED) is ClaimStatus.DETECTED
        )

    def test_lower_claims_pass_through(self, config):
        assert ceiling_for(config, "oprD", ClaimStatus.UNKNOWN) is ClaimStatus.UNKNOWN

    def test_all_ceilings_are_detected(self, config):
        """No gene in the table may be configured above DETECTED."""
        for gene, spec in config.mechanisms.items():
            assert spec.claim_ceiling is ClaimStatus.DETECTED, gene


class TestInterpretDeterminant:
    def test_known_gene_becomes_mechanism(self, config):
        call = stage_mech.interpret_determinant(
            config, _determinant("blaNDM-1"), "imipenem"
        )
        assert call.mechanism == "acquired_determinant"
        assert call.evidence_level is ClaimStatus.DETECTED

    def test_unknown_gene_yields_none(self, config):
        assert (
            stage_mech.interpret_determinant(
                config, _determinant("blaFAKE-1"), "imipenem"
            )
            is None
        )

    def test_notes_state_detection_is_not_resistance(self, config):
        call = stage_mech.interpret_determinant(
            config, _determinant("blaNDM-1"), "imipenem"
        )
        assert "does not establish phenotypic resistance" in call.notes

    def test_input_claim_is_clamped(self, config):
        """Even if a caller passes SUPPORTED, the ceiling applies."""
        call = stage_mech.interpret_determinant(
            config,
            _determinant("blaNDM-1", claim_status=ClaimStatus.SUPPORTED),
            "imipenem",
        )
        assert call.evidence_level is ClaimStatus.DETECTED


class TestInterpretVariant:
    def _variant(self, gene="oprD", variant_type=VariantType.GENE_DISRUPTION, sample="S1"):
        return RegulatorVariant(
            sample_id=sample,
            gene=gene,
            variant=f"{gene}_v1",
            variant_type=variant_type.value,
            position=100,
            reference="A",
            alternate="T",
            effect="test",
            mechanism=None,
            confidence="high",
        )

    def test_oprd_disruption_maps_to_permeability(self, config):
        call = stage_mech.interpret_variant(
            config, self._variant("oprD", VariantType.GENE_DISRUPTION), "imipenem"
        )
        assert call.mechanism == "reduced_permeability"

    def test_oprd_absence_maps_to_permeability(self, config):
        call = stage_mech.interpret_variant(
            config, self._variant("oprD", VariantType.GENE_ABSENCE), "imipenem"
        )
        assert call.mechanism == "reduced_permeability"

    def test_regulator_variant_maps_to_efflux(self, config):
        call = stage_mech.interpret_variant(
            config, self._variant("nfxB", VariantType.SNV), "imipenem"
        )
        assert call.mechanism == "efflux"

    def test_unknown_gene_yields_none(self, config):
        assert (
            stage_mech.interpret_variant(
                config, self._variant("unknownGene", VariantType.SNV), "imipenem"
            )
            is None
        )

    def test_status_is_detected(self, config):
        call = stage_mech.interpret_variant(
            config, self._variant("oprD", VariantType.GENE_DISRUPTION), "imipenem"
        )
        assert call.evidence_level is ClaimStatus.DETECTED


class TestInterpretStructural:
    def _sv(self, gene="oprD", status=StructuralCallStatus.CONFIRMED, sample="S1"):
        return StructuralVariant(
            sample_id=sample,
            variant_id="SV1",
            variant_type="insertion_sequence",
            position=100,
            affected_gene=gene,
            size=1200,
            evidence="contig_breakpoint",
            confidence="high",
            call_status=status,
        )

    def test_confirmed_maps_to_detected(self, config):
        call = stage_mech.interpret_structural(config, self._sv(), "imipenem")
        assert call.evidence_level is ClaimStatus.DETECTED

    def test_candidate_is_predicted_not_detected(self, config):
        """Rule 7: a candidate SV is never reported as confirmed."""
        call = stage_mech.interpret_structural(
            config, self._sv(status=StructuralCallStatus.CANDIDATE), "imipenem"
        )
        assert call.evidence_level is ClaimStatus.PREDICTED

    def test_not_assessable_yields_no_call(self, config):
        """'Cannot assess' must not become 'no mechanism'."""
        assert (
            stage_mech.interpret_structural(
                config, self._sv(status=StructuralCallStatus.NOT_ASSESSABLE), "imipenem"
            )
            is None
        )

    def test_no_affected_gene_yields_none(self, config):
        sv = self._sv()
        without_gene = StructuralVariant(
            sample_id=sv.sample_id,
            variant_id=sv.variant_id,
            variant_type=sv.variant_type,
            position=sv.position,
            affected_gene=None,
            size=sv.size,
            evidence=sv.evidence,
            confidence=sv.confidence,
            call_status=sv.call_status,
        )
        assert stage_mech.interpret_structural(config, without_gene, "imipenem") is None

    def test_unknown_gene_yields_none(self, config):
        assert (
            stage_mech.interpret_structural(
                config, self._sv(gene="unknownGene"), "imipenem"
            )
            is None
        )


class TestIntactLocus:
    def test_oprd_presence_is_locus_intact_not_permeability(self, config):
        """Rule 1: detecting oprD is evidence of an intact locus."""
        call = stage_mech.annotate_intact_locus(
            config, "S1", "oprD", "imipenem", "annotation"
        )
        assert call.mechanism == "locus_intact"
        assert "not evidence of reduced permeability" in call.notes

    def test_non_oprd_gene_yields_none(self, config):
        assert (
            stage_mech.annotate_intact_locus(
                config, "S1", "mexR", "imipenem", "annotation"
            )
            is None
        )

    def test_intact_locus_is_still_only_detected(self, config):
        call = stage_mech.annotate_intact_locus(
            config, "S1", "oprD", "imipenem", "annotation"
        )
        assert call.evidence_level is ClaimStatus.DETECTED
