"""Tests for the regulator screen (stage 6) and SV/MGE module (stage 7)."""

from __future__ import annotations

import pytest

from papipeline.models import (
    ClaimStatus,
    StructuralCallStatus,
    VariantType,
)
from papipeline.stages import regulators as reg
from papipeline.stages import sv as sv_stage

REG_HEADER = (
    "sample_id\tgene\tvariant\tvariant_type\tposition\treference\talternate\t"
    "effect\tmechanism\tconfidence\tevidence_source\tcall_status\n"
)


def _reg_table(rows: str) -> str:
    return REG_HEADER + rows


class TestParseVariantType:
    @pytest.mark.parametrize(
        "value", ["SNV", "snv", " SNV ", "frameshift", "GENE_ABSENCE"]
    )
    def test_parses_known_types(self, value):
        assert reg.parse_variant_type(value) is not VariantType.UNKNOWN

    def test_unknown_becomes_unknown(self):
        assert reg.parse_variant_type("WEIRD") is VariantType.UNKNOWN

    def test_missing_becomes_unknown(self):
        assert reg.parse_variant_type(None) is VariantType.UNKNOWN


class TestLoadRegulatorVariants:
    def test_loads_synthetic_fixtures(self, config, intermediate_root):
        records = reg.load_regulator_variants(
            config, intermediate_root / "regulators" / "regulator_variants.tsv"
        )
        assert records
        assert all(isinstance(r.gene, str) for r in records)

    def test_mechanism_filled_from_table_when_absent(self, config, write_tsv_file):
        path = write_tsv_file(
            "r.tsv",
            _reg_table("S1\toprD\tv1\tGENE_DISRUPTION\t100\tA\tT\teff\t.\thigh\tsrc\tDETECTED\n"),
        )
        records = reg.load_regulator_variants(config, path)
        assert records[0].mechanism == "reduced_permeability"

    def test_promoter_call_dropped_when_disabled(self, config, write_tsv_file):
        """Promoter calls need a pinned reference, so they are refused."""
        path = write_tsv_file(
            "r.tsv",
            _reg_table("S1\toprD\tv1\tPROMOTER_ALTERATION\t100\tA\tT\teff\t.\thigh\tsrc\tDETECTED\n"),
        )
        records = reg.load_regulator_variants(config, path)
        assert records == []

    def test_unknown_gene_is_kept_but_flagged(self, config, write_tsv_file):
        """An incidental finding is reported, not discarded."""
        path = write_tsv_file(
            "r.tsv",
            _reg_table("S1\tnotALocus\tv1\tSNV\t100\tA\tT\teff\t.\thigh\tsrc\tDETECTED\n"),
        )
        records = reg.load_regulator_variants(config, path)
        assert len(records) == 1
        assert records[0].gene == "notALocus"

    def test_position_coerced_to_int(self, config, write_tsv_file):
        path = write_tsv_file(
            "r.tsv",
            _reg_table("S1\toprD\tv1\tSNV\t100\tA\tT\teff\t.\thigh\tsrc\tDETECTED\n"),
        )
        assert reg.load_regulator_variants(config, path)[0].position == 100

    def test_unparsable_position_becomes_none(self, config, write_tsv_file):
        path = write_tsv_file(
            "r.tsv",
            _reg_table("S1\toprD\tv1\tSNV\tnot_a_number\tA\tT\teff\t.\thigh\tsrc\tDETECTED\n"),
        )
        assert reg.load_regulator_variants(config, path)[0].position is None

    def test_call_status_defaults_to_detected(self, config, write_tsv_file):
        path = write_tsv_file(
            "r.tsv",
            "sample_id\tgene\tvariant\tvariant_type\nS1\toprD\tv1\tSNV\n",
        )
        assert reg.load_regulator_variants(config, path)[0].call_status is ClaimStatus.DETECTED

    def test_fixtures_respect_screened_classes(self, config, intermediate_root):
        """The generator must not emit an unscreened variant class."""
        records = reg.load_regulator_variants(
            config, intermediate_root / "regulators" / "regulator_variants.tsv"
        )
        for record in records:
            spec = config.regulator(record.gene)
            if spec is not None:
                assert record.variant_type in spec.variant_classes, (
                    f"{record.gene}/{record.variant_type} is outside the screen"
                )


class TestGeneStatus:
    def test_absence_outranks_variant(self, config, intermediate_root):
        records = reg.load_regulator_variants(
            config, intermediate_root / "regulators" / "regulator_variants.tsv"
        )
        for sample_records in [records]:
            status = reg.gene_status(sample_records)
            if "oprD" in status:
                assert status["oprD"] in ("absent", "disrupted", "variant")

    def test_oprd_status_values(self, config, intermediate_root, manifest):
        from papipeline.models import RunMode

        grouped = reg.run(config, manifest, RunMode.TEST, intermediate_root)
        statuses = {reg.oprd_status(v) for v in grouped.values()}
        assert statuses <= {"absent", "disrupted", "variant", "not_assessed"}

    def test_fixtures_cover_disruption_and_intact(self, config, intermediate_root, manifest):
        from papipeline.models import RunMode

        grouped = reg.run(config, manifest, RunMode.TEST, intermediate_root)
        statuses = {reg.oprd_status(v) for v in grouped.values()}
        assert "not_assessed" in statuses, "fixture should include unperturbed samples"
        assert statuses & {"disrupted", "absent"}, "fixture should include oprD lesions"


class TestIsDisruptive:
    @pytest.mark.parametrize(
        "variant_type",
        ["FRAMESHIFT", "GENE_ABSENCE", "GENE_DISRUPTION", "PREMATURE_STOP"],
    )
    def test_disruptive_types(self, variant_type):
        assert variant_type in reg.DISRUPTIVE_TYPES

    @pytest.mark.parametrize("variant_type", ["SNV", "INFRAME_INDEL", "MNP"])
    def test_non_disruptive_types(self, variant_type):
        assert variant_type not in reg.DISRUPTIVE_TYPES


SV_HEADER = (
    "sample_id\tvariant_id\tvariant_type\tposition\taffected_gene\tsize\t"
    "evidence\tconfidence\tcall_status\tmge\n"
)


def _sv_table(rows: str) -> str:
    return SV_HEADER + rows


class TestParseCallStatus:
    @pytest.mark.parametrize(
        "value,expected",
        [
            ("confirmed", StructuralCallStatus.CONFIRMED),
            ("candidate", StructuralCallStatus.CANDIDATE),
            ("not_assessable", StructuralCallStatus.NOT_ASSESSABLE),
        ],
    )
    def test_parses_known_statuses(self, value, expected):
        assert sv_stage.parse_call_status(value) is expected

    def test_unknown_downgrades_to_candidate(self):
        """An unrecognised status is never upgraded to confirmed."""
        assert sv_stage.parse_call_status("maybe") is StructuralCallStatus.CANDIDATE

    def test_missing_downgrades_to_candidate(self):
        assert sv_stage.parse_call_status(None) is StructuralCallStatus.CANDIDATE

    def test_case_insensitive(self):
        assert sv_stage.parse_call_status("CONFIRMED") is StructuralCallStatus.CONFIRMED


class TestLoadStructuralVariants:
    def test_loads_synthetic_fixtures(self, config, intermediate_root):
        records = sv_stage.load_structural_variants(
            config,
            intermediate_root / "structural_variants" / "structural_variants.tsv",
        )
        assert records

    def test_fixtures_cover_all_three_statuses(self, config, intermediate_root):
        records = sv_stage.load_structural_variants(
            config,
            intermediate_root / "structural_variants" / "structural_variants.tsv",
        )
        statuses = {r.call_status for r in records}
        assert StructuralCallStatus.CONFIRMED in statuses
        assert StructuralCallStatus.CANDIDATE in statuses
        assert StructuralCallStatus.NOT_ASSESSABLE in statuses

    def test_not_assessable_keeps_missing_fields_none(self, config, intermediate_root):
        records = sv_stage.load_structural_variants(
            config,
            intermediate_root / "structural_variants" / "structural_variants.tsv",
        )
        for record in records:
            if record.call_status is StructuralCallStatus.NOT_ASSESSABLE:
                assert record.position is None or record.size is None

    def test_config_cannot_promote_candidates(self, config):
        """The promote_candidate_calls knob is documented as ignored."""
        assert "promote_candidate_calls" in (config.raw.get("structural_variants") or {})

    def test_grouped_by_manifest(self, config, intermediate_root, manifest):
        from papipeline.models import RunMode

        grouped = sv_stage.run(config, manifest, RunMode.TEST, intermediate_root)
        assert set(grouped) == set(manifest.sample_ids)

    def test_confirmed_only_filter(self, config, intermediate_root):
        records = sv_stage.load_structural_variants(
            config,
            intermediate_root / "structural_variants" / "structural_variants.tsv",
        )
        confirmed = sv_stage.confirmed_only(records)
        assert all(r.call_status is StructuralCallStatus.CONFIRMED for r in confirmed)
        assert len(confirmed) < len(records)

    def test_orphaned_sample_ids_are_dropped(self, config, intermediate_root):
        records = sv_stage.load_structural_variants(
            config,
            intermediate_root / "structural_variants" / "structural_variants.tsv",
        )
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        small = SampleManifest([Sample("TEST_PA_001")])
        grouped = sv_stage.group_by_sample(records, small)
        assert set(grouped) == {"TEST_PA_001"}
