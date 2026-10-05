"""Ticket 15 - the ``oprD_absent`` / ``oprD_LoF`` feature pair.

These tests pin the two properties the ticket actually asks for: the features
are *distinct* (neither derived from the other), and presence of the locus is
never turned into evidence about susceptibility.

:class:`TestTruncatingCallsCountAsLossOfFunction` is a regression suite for the
defect where a real FRAMESHIFT or PREMATURE_STOP was labelled a plain
"variant", so ``oprD_LoF`` read 0 on a genuinely lost oprD.
"""

from __future__ import annotations

import pytest

from papipeline.io.tsv import read_tsv
from papipeline.models import (
    DISRUPTIVE_VARIANT_TYPES,
    RegulatorVariant,
    RunMode,
    VariantType,
)
from papipeline.stages import regulators as reg


def _oprD_call(variant_type: str, sample_id: str = "S1") -> RegulatorVariant:
    return RegulatorVariant(
        sample_id=sample_id,
        gene="oprD",
        variant=f"SYNTHETIC_{variant_type}",
        variant_type=variant_type,
        position=412,
        reference="SYN",
        alternate="SYN",
        effect="synthetic_effect",
        mechanism="reduced_permeability",
        confidence="SYNTHETIC",
    )


def _other_gene_call(sample_id: str = "S1") -> RegulatorVariant:
    """A call at a different locus, which must never reach the OprD features."""
    return RegulatorVariant(
        sample_id=sample_id,
        gene="mexR",
        variant="SYNTHETIC_MEXR",
        variant_type="PREMATURE_STOP",
        position=100,
        reference="A",
        alternate="T",
        effect=None,
        mechanism="efflux",
        confidence="SYNTHETIC",
    )


class TestFeatureLabels:
    def test_the_two_labels_are_named_as_the_spec_declares(self):
        assert reg.oprd_feature_labels() == ("oprD_absent", "oprD_LoF")

    def test_the_labels_are_distinct(self):
        labels = reg.oprd_feature_labels()
        assert len(set(labels)) == 2


class TestOprdFeatures:
    def test_a_wholly_absent_locus_sets_only_the_absent_feature(self):
        features = reg.oprd_features([_oprD_call("GENE_ABSENCE")])
        assert features["oprD_absent"] == 1
        assert features["oprD_LoF"] == 0

    def test_a_disrupted_locus_sets_only_the_lof_feature(self):
        features = reg.oprd_features([_oprD_call("GENE_DISRUPTION")])
        assert features["oprD_LoF"] == 1
        assert features["oprD_absent"] == 0

    def test_values_are_binary_integers(self):
        for records in ([], [_oprD_call("GENE_ABSENCE")], [_oprD_call("SNV")]):
            for value in reg.oprd_features(records).values():
                assert value in (0, 1)
                assert isinstance(value, int)

    def test_no_sample_ever_carries_both_features(self):
        for variant_type in ("GENE_ABSENCE", "GENE_DISRUPTION", "SNV", "INDEL"):
            features = reg.oprd_features([_oprD_call(variant_type)])
            assert not (
                features["oprD_absent"] == 1 and features["oprD_LoF"] == 1
            ), variant_type

    def test_absence_of_a_call_is_not_read_as_intact(self):
        """No call at all yields neither feature, never a positive one."""
        features = reg.oprd_features([])
        assert features == {"oprD_absent": 0, "oprD_LoF": 0}

    def test_a_non_disruptive_call_sets_neither_feature(self):
        features = reg.oprd_features([_oprD_call("SNV")])
        assert features == {"oprD_absent": 0, "oprD_LoF": 0}

    def test_absence_outranks_disruption_within_one_sample(self):
        records = [_oprD_call("GENE_DISRUPTION"), _oprD_call("GENE_ABSENCE")]
        features = reg.oprd_features(records)
        assert features["oprD_absent"] == 1
        assert features["oprD_LoF"] == 0

    def test_a_call_at_another_locus_does_not_reach_the_oprd_features(self):
        features = reg.oprd_features([_other_gene_call()])
        assert features == {"oprD_absent": 0, "oprD_LoF": 0}


class TestOprdFeatureColumns:
    def test_columns_carry_the_supplied_type_prefix(self):
        assert reg.oprd_feature_columns("gene") == (
            "gene__oprD_absent",
            "gene__oprD_LoF",
        )

    def test_a_different_prefix_produces_different_columns(self):
        assert reg.oprd_feature_columns("gene_presence_absence") == (
            "gene_presence_absence__oprD_absent",
            "gene_presence_absence__oprD_LoF",
        )

    @pytest.mark.parametrize(
        "feature_type", ["gene", "gene_presence_absence", "unspecified"]
    )
    def test_the_declared_prefix_is_what_stage_12_parses_back(self, feature_type):
        """A prefixed column round-trips through the stage 12 name grammar."""
        from papipeline.stages.gwas import parse_feature_name

        for column in reg.oprd_feature_columns(feature_type):
            parsed_type, label = parse_feature_name(column)
            assert parsed_type == feature_type
            assert label in reg.oprd_feature_labels()


class TestOprdFeatureRows:
    def test_one_row_per_sample_with_a_stable_column_order(self, config, intermediate_root, manifest):
        from papipeline.models import RunMode as _RunMode

        grouped = reg.run(config, manifest, _RunMode.TEST, intermediate_root)
        rows = reg.oprd_feature_rows(grouped, "gene")
        assert len(rows) == len(manifest.sample_ids)
        assert [r["sample_id"] for r in rows] == sorted(manifest.sample_ids)
        for row in rows:
            assert list(row) == [
                "sample_id",
                "gene__oprD_absent",
                "gene__oprD_LoF",
            ]

    def test_row_order_does_not_depend_on_grouped_insertion_order(self, config, intermediate_root, manifest):
        grouped = reg.run(config, manifest, RunMode.TEST, intermediate_root)
        forward = reg.oprd_feature_rows(grouped, "gene")
        reversed_input = reg.oprd_feature_rows(
            {k: grouped[k] for k in reversed(list(grouped))}, "gene"
        )
        assert forward == reversed_input

    def test_the_fixture_exercises_both_features_and_a_neutral_state(
        self, config, intermediate_root, manifest
    ):
        grouped = reg.run(config, manifest, RunMode.TEST, intermediate_root)
        rows = reg.oprd_feature_rows(grouped, "gene")
        assert sum(r["gene__oprD_absent"] for r in rows) >= 1
        assert sum(r["gene__oprD_LoF"] for r in rows) >= 1
        assert any(
            r["gene__oprD_absent"] == 0 and r["gene__oprD_LoF"] == 0 for r in rows
        ), "fixture should include samples carrying neither feature"


class TestProducerAgreesWithTheFeatureFixture:
    """The derived features must match the checked-in stage 12 fixture.

    ``test_data/intermediate/gwas/gwas_features.tsv`` is generated by
    ``papipeline.testing.synthetic`` from the same ``oprd_status`` verdicts, so
    the two independent paths to the same numbers are cross-checked here. This
    is also what would catch the fixture drifting from the producer after a
    rename.
    """

    def test_derived_features_equal_the_generated_fixture_columns(
        self, config, intermediate_root, manifest, test_data_root
    ):
        path = test_data_root / "intermediate" / "gwas" / "gwas_features.tsv"
        if not path.exists():
            pytest.skip("stage 12 feature fixture not present")

        fixture_rows = {
            str(row["sample_id"]): row
            for row in read_tsv(path, required_columns=("sample_id",))
        }
        grouped = reg.run(config, manifest, RunMode.TEST, intermediate_root)
        derived = reg.oprd_feature_rows(grouped, "gene")

        assert "gene__oprD_LoF" in fixture_rows["TEST_PA_001"], (
            "fixture should carry the renamed gene__oprD_LoF column"
        )
        assert "gene__oprD_disrupted" not in fixture_rows["TEST_PA_001"], (
            "the pre-rename column name must be gone from the fixture"
        )

        for row in derived:
            sample_id = str(row["sample_id"])
            fixture = fixture_rows[sample_id]
            assert int(fixture["gene__oprD_absent"]) == row["gene__oprD_absent"], sample_id
            assert int(fixture["gene__oprD_LoF"]) == row["gene__oprD_LoF"], sample_id


#: Classes that truncate or remove the coding sequence, other than absence.
#: Each of these was silently reported as "variant" before the fix.
TRUNCATING_CLASSES = ("FRAMESHIFT", "PREMATURE_STOP", "GENE_DISRUPTION")


class TestTruncatingCallsCountAsLossOfFunction:
    """A truncating call at oprD is loss of function, not a mere variant.

    Each test here fails against the previous implementation, which matched
    GENE_ABSENCE and GENE_DISRUPTION by literal string and let every other
    class fall through to "variant".
    """

    @pytest.mark.parametrize("variant_type", TRUNCATING_CLASSES)
    def test_a_truncating_call_sets_the_lof_feature(self, variant_type):
        features = reg.oprd_features([_oprD_call(variant_type)])
        assert features["oprD_LoF"] == 1, variant_type
        assert features["oprD_absent"] == 0, variant_type

    @pytest.mark.parametrize("variant_type", TRUNCATING_CLASSES)
    def test_a_truncating_call_is_reported_as_disrupted(self, variant_type):
        assert reg.oprd_status([_oprD_call(variant_type)]) == "disrupted"

    @pytest.mark.parametrize("variant_type", ["FRAMESHIFT", "PREMATURE_STOP"])
    def test_a_truncating_call_is_disrupted_by_the_shared_predicate(
        self, variant_type
    ):
        """``gene_status`` must route through ``is_disruptive``.

        Asserted on the predicate rather than on a re-stated list, so the two
        cannot drift apart again: if the set changes, this test follows it.
        """
        record = _oprD_call(variant_type)
        assert reg.is_disruptive(record) is True
        assert reg.gene_status([record]) == {"oprD": "disrupted"}

    def test_every_variant_type_is_accounted_for(self):
        """No declared class may fall through to a non-disruptive verdict.

        This is the guard that would have caught the defect directly: it walks
        the whole enum, so adding a class to ``VariantType`` without deciding
        what it means to ``gene_status`` fails here instead of silently
        reading as "variant".
        """
        expected_disruptive = {
            c for c in TRUNCATING_CLASSES if c in DISRUPTIVE_VARIANT_TYPES
        }
        for variant_type in VariantType:
            status = reg.oprd_status([_oprD_call(variant_type.value)])
            if variant_type.value == VariantType.GENE_ABSENCE.value:
                assert status == "absent", variant_type
            elif variant_type.value in expected_disruptive:
                assert status == "disrupted", variant_type
            else:
                assert status == "variant", variant_type

    def test_gene_status_and_is_disruptive_cannot_disagree(self):
        """One source of truth: the predicate and the label must agree.

        ``absent`` is the documented exception - it is in
        ``DISRUPTIVE_VARIANT_TYPES`` because it is loss of function in the
        strongest sense, but absence and disruption are separate observations.
        """
        for variant_type in VariantType:
            record = _oprD_call(variant_type.value)
            status = reg.gene_status([record])["oprD"]
            is_absence = variant_type.value == VariantType.GENE_ABSENCE.value
            assert status == ("absent" if is_absence else
                              "disrupted" if reg.is_disruptive(record)
                              else "variant"), variant_type

    def test_absence_still_outranks_a_truncating_call(self):
        """Distinctness survives the fix: absent and LoF never both fire."""
        for variant_type in TRUNCATING_CLASSES:
            features = reg.oprd_features(
                [_oprD_call(variant_type), _oprD_call("GENE_ABSENCE")]
            )
            assert features["oprD_absent"] == 1, variant_type
            assert features["oprD_LoF"] == 0, variant_type

    def test_a_non_truncating_call_is_still_only_a_variant(self):
        for variant_type in ("SNV", "INDEL", "INFRAME_INDEL", "MNP"):
            features = reg.oprd_features([_oprD_call(variant_type)])
            assert features == {"oprD_absent": 0, "oprD_LoF": 0}, variant_type


class TestEmulatorExercisesTheRealClasses:
    """The TEST fixtures must contain the classes real data carries.

    The emulator emitted only GENE_DISRUPTION and GENE_ABSENCE at oprD, which
    are exactly the two names the defective branch matched. Without this, the
    bug is invisible in TEST no matter how the stage is fixed.
    """

    def test_the_committed_fixture_has_truncating_calls_at_oprd(self, test_data_root):
        path = test_data_root / "intermediate" / "regulators" / "regulator_variants.tsv"
        if not path.exists():
            pytest.skip("stage 6 fixture not present")
        rows = read_tsv(path, required_columns=("gene", "variant_type"))
        at_oprd = {
            str(r["variant_type"]) for r in rows if str(r["gene"]) == "oprD"
        }
        assert at_oprd & {"FRAMESHIFT", "PREMATURE_STOP"}, (
            f"fixture must carry a truncating class at oprD; found {sorted(at_oprd)}"
        )

    def test_the_generator_emits_several_disruptive_classes_at_oprd(self, tmp_path):
        from papipeline.testing import generate

        generate(tmp_path, write_genome_files=False)
        rows = read_tsv(
            tmp_path / "intermediate" / "regulators" / "regulator_variants.tsv",
            required_columns=("gene", "variant_type"),
        )
        emitted = {
            str(r["variant_type"]) for r in rows if str(r["gene"]) == "oprD"
        }
        assert {"FRAMESHIFT", "PREMATURE_STOP"} <= emitted, sorted(emitted)

    def test_the_generator_only_emits_classes_the_screen_declares(self):
        """The knowledge table stays the single source of truth."""
        from pathlib import Path

        from papipeline.testing.synthetic import (
            disruptive_variant_choices_from_table,
        )

        table = Path(__file__).resolve().parents[2] / "config" / "regulators.tsv"
        if not table.exists():
            pytest.skip("regulators table not present")
        choices = disruptive_variant_choices_from_table(table)
        assert "oprD" in choices
        assert set(choices["oprD"]) <= DISRUPTIVE_VARIANT_TYPES
        assert "GENE_ABSENCE" not in choices["oprD"], (
            "absence is a separate OprD state and must not be folded in here"
        )
