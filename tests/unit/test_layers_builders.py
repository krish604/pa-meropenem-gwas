"""Layer feature builders: ticket 15's L1-L5 encoding over synthetic sources.

Every test here hands the builders a hand-built :class:`Sources` of synthetic
rows and asserts the exact feature columns and values that come out. Nothing
reads ``data/``, ``db/`` or ``PDC_essential.tsv``; nothing runs a real-mode
pipeline.

Decisions pinned here:

* **Layer prefixes are load-bearing.** ``L1_``-``L5_`` is how a reader, a
  provenance banner and the dictionary locate a feature's evidence, so every
  builder must emit it and ``build_layers`` must not lose it.
* **The oprD loss tiers are distinct features, and a missense is neither.**
  ``oprD_absent``, ``oprD_LoF_tier1`` (stop, frameshift, IS insertion, large
  deletion), ``oprD_LoF_tier2`` (probable loss) and the ``oprD_off_any``
  composite are built separately; any other oprD call stays an individual
  feature and never sets the composite.
* **Composites are loss-based.** ``L2_PDC_high`` needs a PDC allele *and* an
  ampD/ampR/dacB loss; the pump-high composites need a loss among their group.
  A missense never counts towards either.
* **The groupings that are not spelled out by the ticket come from the
  knowledge tables**, not from this module: the efflux pump groups are read
  from ``config/mechanisms.tsv`` ``biological_role``.
* **The partial-call vocabulary is the pipeline's own** - the three states are
  asserted to be a subset of ``pilot.pdc_fields.KNOWN_CALL_STATES`` rather than
  restated from memory.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Set, Tuple

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.layers.builders import (
    build_layer1,
    build_layer2,
    build_layer3,
    build_layer4,
    build_layer5,
    build_layers,
)
from papipeline.layers.features import Feature, LayerBundle
from papipeline.layers.inputs import (
    PARTIAL_CALL_STATES,
    Sources,
    partial_call_flags,
)
from papipeline.models import ClaimStatus, RegulatorVariant

SAMPLES: Tuple[str, ...] = tuple(f"TEST_PA_{i:03d}" for i in range(1, 13))


def _expect(on: Sequence[str] | Set[str]) -> Dict[str, int]:
    """The binary column a builder should produce for ``on``."""
    on = set(on)
    return {sample_id: int(sample_id in on) for sample_id in SAMPLES}


def _sources(**kwargs) -> Sources:
    kwargs.setdefault("antibiotic", "imipenem")
    kwargs.setdefault("sample_ids", SAMPLES)
    return Sources(**kwargs)


def _amr(
    sample_id: str,
    gene: str,
    *,
    variant: str = "SYNTHETIC_ALLELE",
    antibiotic: str = "imipenem",
    determinant: str | None = None,
    **extra,
) -> Dict[str, str]:
    row = {
        "sample_id": sample_id,
        "antibiotic": antibiotic,
        "determinant": determinant or f"{gene}~SYNTHETIC_ALLELE",
        "gene": gene,
        "variant": variant,
    }
    row.update(extra)
    return row


def _reg(
    sample_id: str,
    gene: str,
    variant_type: str,
    *,
    variant: str | None = None,
    effect: str = "synthetic_effect",
) -> RegulatorVariant:
    return RegulatorVariant(
        sample_id=sample_id,
        gene=gene,
        variant=variant or f"SYNTHETIC_{gene}_{variant_type}",
        variant_type=variant_type,
        position=None,
        reference=None,
        alternate=None,
        effect=effect,
        mechanism=None,
        confidence="high",
        evidence_source="synthetic",
        call_status=ClaimStatus.DETECTED,
    )


def _sv(
    sample_id: str,
    gene: str,
    variant_type: str,
    call_status: str = "confirmed",
) -> Dict[str, str]:
    return {
        "sample_id": sample_id,
        "variant_id": f"{sample_id}_SV_1",
        "variant_type": variant_type,
        "call_status": call_status,
        "affected_gene": gene,
    }


def _by_name(features: List[Feature]) -> Dict[str, Feature]:
    return {feature.name: feature for feature in features}


class TestLayer1Acquired:
    def test_any_mbl_unions_the_three_families_the_ticket_names(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaNDM-1"),
                _amr("TEST_PA_002", "blaVIM-1"),
                _amr("TEST_PA_003", "blaIMP-18"),
                _amr("TEST_PA_004", "blaNDM-5"),
            )
        )
        features, values = build_layer1(sources, SAMPLES)
        assert values["L1_any_MBL"] == _expect(
            {"TEST_PA_001", "TEST_PA_002", "TEST_PA_003", "TEST_PA_004"}
        )
        assert {feature.layer for feature in features} == {1}
        assert all(feature.name.startswith("L1_") for feature in features)

    def test_each_carbapenemase_group_is_its_own_column(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaKPC-2"),
                _amr("TEST_PA_002", "blaGES-1"),
                _amr("TEST_PA_003", "blaOXA-48"),
                _amr("TEST_PA_004", "blaOXA-23"),
            )
        )
        _, values = build_layer1(sources, SAMPLES)
        assert values["L1_KPC"] == _expect({"TEST_PA_001"})
        assert values["L1_GES"] == _expect({"TEST_PA_002"})
        assert values["L1_OXA"] == _expect({"TEST_PA_003", "TEST_PA_004"})
        assert values["L1_any_MBL"] == _expect(set())

    def test_genes_outside_the_ticketed_groups_never_set_a_group_column(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaTEM-1"),
                _amr("TEST_PA_002", "blaCTX-M-15"),
                _amr("TEST_PA_003", "blaSPM-1"),
            )
        )
        features, values = build_layer1(sources, SAMPLES)
        assert {feature.name for feature in features} >= {
            "L1_any_MBL",
            "L1_KPC",
            "L1_GES",
            "L1_OXA",
        }
        for column in ("L1_any_MBL", "L1_KPC", "L1_GES", "L1_OXA"):
            assert values[column] == _expect(set())
        assert not any("TEM" in name or "CTX" in name or "SPM" in name for name in values)

    def test_the_gene_symbol_falls_back_to_the_determinant(self):
        sources = _sources(
            amr_rows=(_amr("TEST_PA_001", "", determinant="blaVIM-1~SYNTHETIC"),)
        )
        _, values = build_layer1(sources, SAMPLES)
        assert values["L1_any_MBL"] == _expect({"TEST_PA_001"})

    def test_rows_another_antibiotic_claims_never_count(self):
        sources = _sources(
            amr_rows=(_amr("TEST_PA_001", "blaNDM-1", antibiotic="meropenem"),)
        )
        _, values = build_layer1(sources, SAMPLES)
        assert values["L1_any_MBL"] == _expect(set())

    def test_the_provenance_token_names_the_source_table(self):
        sources = _sources(amr_rows=(_amr("TEST_PA_001", "blaKPC-2"),))
        features, _ = build_layer1(sources, SAMPLES)
        assert any(
            "amr_determinants.tsv" in feature.source_tokens for feature in features
        )


class TestLayer2Pdc:
    def test_each_observed_pdc_allele_gets_its_own_column(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaPDC-3", variant="V216I"),
                _amr("TEST_PA_002", "blaPDC-17", variant="D146del"),
            )
        )
        features, values = build_layer2(sources, SAMPLES)
        assert values["L2_blaPDC-3"] == _expect({"TEST_PA_001"})
        assert values["L2_blaPDC-17"] == _expect({"TEST_PA_002"})
        assert all(feature.name.startswith("L2_") for feature in features)

    def test_residue_positions_become_flags_and_unparsable_variants_do_not(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaPDC-3", variant="V216I"),
                _amr("TEST_PA_002", "blaPDC-3", variant="-"),
                _amr("TEST_PA_003", "blaPDC-3", variant="D146del"),
            )
        )
        _, values = build_layer2(sources, SAMPLES)
        assert values["L2_PDC_pos_216"] == _expect({"TEST_PA_001"})
        assert values["L2_PDC_pos_146"] == _expect({"TEST_PA_003"})
        position_columns = {name for name in values if name.startswith("L2_PDC_pos_")}
        assert position_columns == {"L2_PDC_pos_216", "L2_PDC_pos_146"}

    def test_a_non_pdc_gene_never_produces_an_allele_or_position(self):
        sources = _sources(amr_rows=(_amr("TEST_PA_001", "blaKPC-2", variant="V216I"),))
        _, values = build_layer2(sources, SAMPLES)
        assert not any(name.startswith("L2_PDC_pos_") for name in values)
        assert not any(name.startswith("L2_bla") for name in values)

    def test_regulator_loss_and_missense_are_separate_features(self):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "ampD", "FRAMESHIFT"),
                _reg("TEST_PA_002", "ampR", "SNV"),
                _reg("TEST_PA_003", "dacB", "GENE_ABSENCE"),
            )
        )
        _, values = build_layer2(sources, SAMPLES)
        assert values["L2_ampD_loss"] == _expect({"TEST_PA_001"})
        assert values["L2_ampD_missense"] == _expect(set())
        assert values["L2_ampR_loss"] == _expect(set())
        assert values["L2_ampR_missense"] == _expect({"TEST_PA_002"})
        assert values["L2_dacB_loss"] == _expect({"TEST_PA_003"})
        assert values["L2_dacB_missense"] == _expect(set())

    def test_a_missense_is_never_counted_as_a_loss(self):
        sources = _sources(regulator_records=(_reg("TEST_PA_001", "ampD", "SNV"),))
        _, values = build_layer2(sources, SAMPLES)
        assert values["L2_ampD_loss"] == _expect(set())
        assert values["L2_ampD_missense"] == _expect({"TEST_PA_001"})

    def test_pdc_high_needs_an_allele_and_a_regulator_loss(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaPDC-3", variant="V216I"),
                _amr("TEST_PA_002", "blaPDC-3", variant="V216I"),
            ),
            regulator_records=(
                _reg("TEST_PA_001", "ampD", "FRAMESHIFT"),
                _reg("TEST_PA_003", "ampD", "FRAMESHIFT"),
            ),
        )
        _, values = build_layer2(sources, SAMPLES)
        assert values["L2_PDC_high"] == _expect({"TEST_PA_001"})

    def test_only_the_ticketed_pdc_regulators_feed_the_composite(self):
        sources = _sources(
            amr_rows=(_amr("TEST_PA_001", "blaPDC-3", variant="V216I"),),
            regulator_records=(_reg("TEST_PA_001", "mexR", "FRAMESHIFT"),),
        )
        _, values = build_layer2(sources, SAMPLES)
        assert values["L2_PDC_high"] == _expect(set())

    def test_the_composite_is_always_defined_and_the_dictionary_is_populated(self):
        sources = _sources()
        features, values = build_layer2(sources, SAMPLES)
        assert "L2_PDC_high" in values
        assert all(all(v == 0 for v in values["L2_PDC_high"].values()) for _ in (0,))
        catalogue = _by_name(features)
        assert "ampD" in catalogue["L2_ampD_loss"].rule
        assert "regulator_variants.tsv" in catalogue["L2_ampD_loss"].source_tokens
        assert "amr_determinants.tsv" in catalogue["L2_PDC_high"].source_tokens


class TestLayer3OprdTiers:
    def test_the_four_canonical_features_exist_even_with_no_call(self):
        _, values = build_layer3(_sources(), SAMPLES)
        assert set(values) >= {
            "L3_oprD_absent",
            "L3_oprD_LoF_tier1",
            "L3_oprD_LoF_tier2",
            "L3_oprD_off_any",
        }
        for column in (
            "L3_oprD_absent",
            "L3_oprD_LoF_tier1",
            "L3_oprD_LoF_tier2",
            "L3_oprD_off_any",
        ):
            assert values[column] == _expect(set())

    def test_stop_frameshift_and_disruption_are_tier_one_absence_is_its_own(
        self,
    ):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "oprD", "PREMATURE_STOP"),
                _reg("TEST_PA_002", "oprD", "FRAMESHIFT"),
                _reg(
                    "TEST_PA_003",
                    "oprD",
                    "GENE_DISRUPTION",
                    effect="disrupted_by_insertion",
                ),
                _reg("TEST_PA_004", "oprD", "GENE_ABSENCE"),
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_LoF_tier1"] == _expect(
            {"TEST_PA_001", "TEST_PA_002", "TEST_PA_003"}
        )
        assert values["L3_oprD_absent"] == _expect({"TEST_PA_004"})
        assert values["L3_oprD_off_any"] == _expect(
            {"TEST_PA_001", "TEST_PA_002", "TEST_PA_003", "TEST_PA_004"}
        )

    def test_a_confirmed_structural_insertion_or_deletion_is_tier_one(self):
        sources = _sources(
            sv_rows=(
                _sv("TEST_PA_001", "oprD", "insertion_sequence", "confirmed"),
                _sv("TEST_PA_002", "oprD", "deletion", "confirmed"),
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_LoF_tier1"] == _expect(
            {"TEST_PA_001", "TEST_PA_002"}
        )
        assert values["L3_oprD_absent"] == _expect(set())

    def test_a_candidate_structural_variant_is_probable_loss_only(self):
        sources = _sources(
            sv_rows=(
                _sv("TEST_PA_001", "oprD", "deletion", "candidate"),
                _sv("TEST_PA_002", "oprD", "insertion_sequence", "candidate"),
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_LoF_tier2"] == _expect({"TEST_PA_001", "TEST_PA_002"})
        assert values["L3_oprD_LoF_tier1"] == _expect(set())

    def test_not_assessable_is_neither_absence_nor_loss(self):
        sources = _sources(
            sv_rows=(_sv("TEST_PA_001", "oprD", "deletion", "not_assessable"),)
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_absent"] == _expect(set())
        assert values["L3_oprD_LoF_tier1"] == _expect(set())
        assert values["L3_oprD_LoF_tier2"] == _expect(set())
        assert values["L3_oprD_off_any"] == _expect(set())

    def test_a_structural_call_at_another_gene_never_reaches_oprd(self):
        sources = _sources(
            sv_rows=(_sv("TEST_PA_001", "mexXY", "deletion", "confirmed"),)
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_off_any"] == _expect(set())

    def test_an_inframe_indel_is_probable_loss_not_tier_one(self):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "oprD", "INFRAME_INDEL"),
                _reg("TEST_PA_002", "oprD", "INDEL"),
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_LoF_tier2"] == _expect({"TEST_PA_001", "TEST_PA_002"})
        assert values["L3_oprD_LoF_tier1"] == _expect(set())

    def test_a_missense_stays_individual_and_never_sets_the_composite(self):
        sources = _sources(
            regulator_records=(_reg("TEST_PA_001", "oprD", "SNV", variant="V359L"),)
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_V359L"] == _expect({"TEST_PA_001"})
        assert values["L3_oprD_off_any"] == _expect(set())
        assert values["L3_oprD_LoF_tier1"] == _expect(set())
        assert values["L3_oprD_absent"] == _expect(set())

    def test_every_other_oprd_call_becomes_its_own_feature(self):
        sources = _sources(
            regulator_records=(
                _reg(
                    "TEST_PA_001",
                    "oprD",
                    "PROMOTER_ALTERATION",
                    variant="PROM_1",
                ),
                _reg("TEST_PA_002", "oprD", "MNP", variant="A100G"),
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_PROM_1"] == _expect({"TEST_PA_001"})
        assert values["L3_oprD_A100G"] == _expect({"TEST_PA_002"})
        assert values["L3_oprD_off_any"] == _expect(set())

    def test_the_optional_structural_calls_table_feeds_absence_and_tier_one(self):
        sources = _sources(
            oprd_structural_rows=(
                {
                    "sample_id": "TEST_PA_001",
                    "structural_verdict": "absent",
                    "lesion_type": None,
                },
                {
                    "sample_id": "TEST_PA_002",
                    "structural_verdict": "disrupted",
                    "lesion_type": "frameshift",
                },
                {
                    "sample_id": "TEST_PA_003",
                    "structural_verdict": "disrupted",
                    "lesion_type": "premature_stop",
                },
                {
                    "sample_id": "TEST_PA_004",
                    "structural_verdict": "intact",
                    "lesion_type": None,
                },
                {
                    "sample_id": "TEST_PA_005",
                    "structural_verdict": "not_assessed",
                    "lesion_type": None,
                },
            )
        )
        _, values = build_layer3(sources, SAMPLES)
        assert values["L3_oprD_absent"] == _expect({"TEST_PA_001"})
        assert values["L3_oprD_LoF_tier1"] == _expect(
            {"TEST_PA_002", "TEST_PA_003"}
        )
        assert values["L3_oprD_off_any"] == _expect(
            {"TEST_PA_001", "TEST_PA_002", "TEST_PA_003"}
        )

    def test_the_sources_are_named_in_the_dictionary(self):
        features, _ = build_layer3(
            _sources(
                regulator_records=(_reg("TEST_PA_001", "oprD", "FRAMESHIFT"),),
                sv_rows=(_sv("TEST_PA_002", "oprD", "deletion", "confirmed"),),
            ),
            SAMPLES,
        )
        catalogue = _by_name(features)
        tier1 = catalogue["L3_oprD_LoF_tier1"]
        assert "regulator_variants.tsv" in tier1.source_tokens
        assert "structural_variants.tsv" in tier1.source_tokens
        assert catalogue["L3_oprD_absent"].rule
        assert catalogue["L3_oprD_off_any"].rule


class TestLayer4Efflux:
    def test_loss_and_missense_are_separate_for_each_named_regulator(
        self, config: PipelineConfig
    ):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "mexR", "FRAMESHIFT"),
                _reg("TEST_PA_002", "mexZ", "SNV"),
                _reg("TEST_PA_003", "nfxB", "PREMATURE_STOP"),
                _reg("TEST_PA_004", "nalC", "GENE_ABSENCE"),
            )
        )
        _, values = build_layer4(config, sources, SAMPLES)
        assert values["L4_mexR_loss"] == _expect({"TEST_PA_001"})
        assert values["L4_mexR_missense"] == _expect(set())
        assert values["L4_mexZ_loss"] == _expect(set())
        assert values["L4_mexZ_missense"] == _expect({"TEST_PA_002"})
        assert values["L4_nfxB_loss"] == _expect({"TEST_PA_003"})
        assert values["L4_nalC_loss"] == _expect({"TEST_PA_004"})
        assert values["L4_nalD_loss"] == _expect(set())

    def test_the_five_named_regulators_each_get_both_columns(
        self, config: PipelineConfig
    ):
        features, values = build_layer4(config, _sources(), SAMPLES)
        for gene in ("mexR", "nalC", "nalD", "nfxB", "mexZ"):
            assert f"L4_{gene}_loss" in values
            assert f"L4_{gene}_missense" in values
        assert all(feature.name.startswith("L4_") for feature in features)

    def test_regulators_outside_the_ticketed_five_are_not_encoded(
        self, config: PipelineConfig
    ):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "mexT", "GENE_ABSENCE"),
                _reg("TEST_PA_002", "mexS", "PREMATURE_STOP"),
                _reg("TEST_PA_003", "oprD", "FRAMESHIFT"),
            )
        )
        _, values = build_layer4(config, sources, SAMPLES)
        assert not any("mexT" in name or "mexS" in name for name in values)
        assert not any("oprD" in name for name in values)

    def test_the_pump_groups_are_read_from_the_mechanism_table(
        self, config: PipelineConfig
    ):
        features, values = build_layer4(config, _sources(), SAMPLES)
        pump_columns = {name for name in values if name.startswith("L4_pump_high_")}
        # config/mechanisms.tsv biological_role: mexR represses the mexAB
        # operon; nalC/nalD/mexZ/nfxB are MexXY machinery or its repressors.
        assert pump_columns == {"L4_pump_high_mexAB", "L4_pump_high_MexXY"}
        catalogue = _by_name(features)
        assert "mechanisms.tsv" in catalogue["L4_pump_high_mexAB"].source_tokens

    def test_a_loss_sets_its_pump_group_and_a_missense_does_not(
        self, config: PipelineConfig
    ):
        sources = _sources(
            regulator_records=(
                _reg("TEST_PA_001", "mexR", "FRAMESHIFT"),
                _reg("TEST_PA_002", "mexZ", "SNV"),
                _reg("TEST_PA_003", "nfxB", "PREMATURE_STOP"),
            )
        )
        _, values = build_layer4(config, sources, SAMPLES)
        assert values["L4_pump_high_mexAB"] == _expect({"TEST_PA_001"})
        # mexZ's SNV alone must not read as pump-high; nfxB's stop must.
        assert values["L4_pump_high_MexXY"] == _expect({"TEST_PA_003"})


class TestLayer5Targets:
    def test_each_target_gene_gets_an_any_variant_column(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "gyrA", variant="S83L"),
                _amr("TEST_PA_002", "parC", variant="S80L"),
                _amr("TEST_PA_003", "ftsI", variant="P533L"),
            )
        )
        _, values = build_layer5(sources, SAMPLES, "imipenem")
        assert values["L5_gyrA"] == _expect({"TEST_PA_001"})
        assert values["L5_parC"] == _expect({"TEST_PA_002"})
        assert values["L5_ftsI"] == _expect({"TEST_PA_003"})
        assert values["L5_any_QRDR"] == _expect({"TEST_PA_001", "TEST_PA_002"})

    def test_the_qrdr_count_counts_distinct_positions_on_the_two_targets(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "gyrA", variant="S83L"),
                _amr("TEST_PA_001", "gyrA", variant="D87N"),
                _amr("TEST_PA_002", "parC", variant="S80L"),
                _amr("TEST_PA_003", "ftsI", variant="P533L"),
            )
        )
        features, values = build_layer5(sources, SAMPLES, "imipenem")
        assert values["L5_QRDR_count"]["TEST_PA_001"] == 2
        assert values["L5_QRDR_count"]["TEST_PA_002"] == 1
        # ftsI is a target but not a quinolone target: no QRDR position.
        assert values["L5_QRDR_count"]["TEST_PA_003"] == 0
        catalogue = _by_name(features)
        assert catalogue["L5_QRDR_count"].value_kind == "count"
        assert catalogue["L5_any_QRDR"].value_kind == "binary"

    def test_a_variant_that_parses_to_no_residue_still_sets_the_gene_column(self):
        sources = _sources(amr_rows=(_amr("TEST_PA_001", "gyrA", variant="unknown"),))
        _, values = build_layer5(sources, SAMPLES, "imipenem")
        assert values["L5_gyrA"] == _expect({"TEST_PA_001"})
        assert values["L5_any_QRDR"] == _expect(set())
        assert values["L5_QRDR_count"]["TEST_PA_001"] == 0

    def test_the_same_position_on_two_isolates_counts_once_each(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "gyrA", variant="S83L"),
                _amr("TEST_PA_002", "gyrA", variant="S83L"),
            )
        )
        _, values = build_layer5(sources, SAMPLES, "imipenem")
        assert values["L5_QRDR_count"] == {
            sample_id: int(sample_id in {"TEST_PA_001", "TEST_PA_002"})
            for sample_id in SAMPLES
        }

    def test_rows_for_another_antibiotic_never_count(self):
        sources = _sources(
            amr_rows=(_amr("TEST_PA_001", "gyrA", variant="S83L", antibiotic="ceftazidime"),)
        )
        _, values = build_layer5(sources, SAMPLES, "imipenem")
        assert values["L5_gyrA"] == _expect(set())

    def test_every_feature_carries_the_layer_five_prefix_and_its_table(self):
        sources = _sources(amr_rows=(_amr("TEST_PA_001", "parC", variant="S80L"),))
        features, _ = build_layer5(sources, SAMPLES, "imipenem")
        assert all(feature.name.startswith("L5_") and feature.layer == 5 for feature in features)
        assert all(
            "amr_determinants.tsv" in feature.source_tokens for feature in features
        )


class TestPartialCallFlag:
    def test_the_three_recorded_states_set_the_flag(self):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaKPC-2", call_state="PARTIAL"),
                _amr("TEST_PA_002", "blaKPC-2", call_state="MISTRANSLATION"),
                _amr("TEST_PA_003", "blaKPC-2", call_state="HMM"),
                _amr("TEST_PA_004", "blaKPC-2", call_state="COMPLETE"),
            )
        )
        flags = partial_call_flags(sources)
        assert flags == {
            **{sample_id: 0 for sample_id in SAMPLES},
            "TEST_PA_001": 1,
            "TEST_PA_002": 1,
            "TEST_PA_003": 1,
        }

    def test_the_flag_looks_at_every_source_table(self):
        sources = _sources(
            regulator_records=(_reg("TEST_PA_001", "oprD", "FRAMESHIFT"),),
            regulator_rows=(
                {
                    "sample_id": "TEST_PA_001",
                    "gene": "oprD",
                    "variant": "V1",
                    "variant_type": "FRAMESHIFT",
                    "call_state": "PARTIAL",
                },
            ),
            sv_rows=(_sv("TEST_PA_002", "oprD", "deletion", "candidate"),),
        )
        flags = partial_call_flags(sources)
        assert flags["TEST_PA_001"] == 1
        assert flags["TEST_PA_002"] == 0

    def test_a_cohort_with_no_call_state_column_reports_missing_not_zero(self):
        sources = _sources(amr_rows=(_amr("TEST_PA_001", "blaKPC-2"),))
        flags = partial_call_flags(sources)
        assert set(flags.values()) == {None}

    def test_the_states_are_the_pipeline_vocabulary(self):
        from papipeline.pilot.pdc_fields import KNOWN_CALL_STATES

        assert PARTIAL_CALL_STATES == frozenset(
            {"PARTIAL", "MISTRANSLATION", "HMM"}
        )
        # KNOWN_CALL_STATES is a tuple in pilot.pdc_fields; the assertion is
        # the subset relation, not the container type.
        assert set(PARTIAL_CALL_STATES) <= set(KNOWN_CALL_STATES)


class TestTheComposedBundle:
    def test_all_five_layers_reach_the_bundle_with_prefixes(
        self, config: PipelineConfig
    ):
        bundle = build_layers(config, _sources(), SAMPLES)
        assert isinstance(bundle, LayerBundle)
        assert bundle.sample_ids == SAMPLES
        assert {feature.layer for feature in bundle.features} == {1, 2, 3, 4, 5}
        assert all(
            feature.name.startswith(f"L{feature.layer}_")
            for feature in bundle.features
        )
        assert all(
            len(values) == len(SAMPLES) for values in bundle.values.values()
        )

    def test_values_are_integers_zero_or_the_count(
        self, config: PipelineConfig
    ):
        sources = _sources(
            amr_rows=(
                _amr("TEST_PA_001", "blaNDM-1"),
                _amr("TEST_PA_002", "gyrA", variant="S83L"),
            ),
            regulator_records=(_reg("TEST_PA_003", "oprD", "FRAMESHIFT"),),
        )
        bundle = build_layers(config, sources, SAMPLES)
        for name, values in bundle.values.items():
            for sample_id, value in values.items():
                assert isinstance(value, int), name
                assert value >= 0, name
        assert bundle.values["L1_any_MBL"]["TEST_PA_001"] == 1
        assert bundle.values["L3_oprD_LoF_tier1"]["TEST_PA_003"] == 1
        assert bundle.values["L5_QRDR_count"]["TEST_PA_002"] == 1

    def test_no_source_row_is_silently_dropped_from_the_cohort(
        self, config: PipelineConfig
    ):
        sources = _sources(
            amr_rows=(_amr("TEST_PA_012", "blaKPC-2"),),
            regulator_records=(_reg("TEST_PA_012", "oprD", "GENE_ABSENCE"),),
        )
        bundle = build_layers(config, sources, SAMPLES)
        assert bundle.values["L1_KPC"]["TEST_PA_012"] == 1
        assert bundle.values["L3_oprD_absent"]["TEST_PA_012"] == 1
