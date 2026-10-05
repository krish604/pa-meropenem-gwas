"""Tests for the phenotype interface (stage 11).

The most important tests here are the negative ones: they assert that the
pipeline refuses to fabricate quantitative data and refuses to guess a
category. See docs/scientific_rules.md rules 1 and 6.
"""

from __future__ import annotations

import pytest

from papipeline.errors import PhenotypeError
from papipeline.models import Phenotype
from papipeline.stages import phenotype as stage

HEADER = "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\n"


def _table(rows: str) -> str:
    return HEADER + rows


class TestParsePhenotypeRow:
    def test_accepts_every_allowed_value(self, config):
        for value in ("R", "I", "S", "SDD", "ND"):
            row = {
                "sample_id": "S1",
                "antibiotic": "imipenem",
                "phenotype": value,
            }
            call = stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )
            assert call.phenotype.value == value

    def test_lowercase_is_normalised(self, config):
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": "r"}
        call = stage.parse_phenotype_row(
            row, config.allowed_phenotypes, config.antibiotics
        )
        assert call.phenotype is Phenotype.R

    def test_invalid_phenotype_raises(self, config):
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": "RESISTANT"}
        with pytest.raises(PhenotypeError, match="not in the allowed set"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    @pytest.mark.parametrize("value", ["", ".", "NA", None])
    def test_missing_phenotype_raises(self, config, value):
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": value}
        with pytest.raises(PhenotypeError):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    def test_missing_sample_id_raises(self, config):
        row = {"antibiotic": "imipenem", "phenotype": "R"}
        with pytest.raises(PhenotypeError, match="no sample_id"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    def test_invalid_antibiotic_raises(self, config):
        row = {"sample_id": "S1", "antibiotic": "ceftazidime", "phenotype": "R"}
        with pytest.raises(PhenotypeError, match="not configured"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    def test_missing_antibiotic_raises(self, config):
        row = {"sample_id": "S1", "phenotype": "R"}
        with pytest.raises(PhenotypeError, match="no antibiotic"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )


class TestNoFabricatedQuantitation:
    """R/I/S must never acquire an MIC the laboratory did not measure."""

    def test_ric_s_carry_no_mic(self, config):
        for value in ("R", "I", "S"):
            row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": value}
            call = stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )
            assert call.mic is None, f"{value} must not acquire an MIC"
            assert call.zone_diameter is None

    def test_measured_mic_is_preserved(self, config):
        row = {
            "sample_id": "S1",
            "antibiotic": "imipenem",
            "phenotype": "R",
            "MIC": "16",
            "MIC_unit": "mg/L",
        }
        call = stage.parse_phenotype_row(
            row, config.allowed_phenotypes, config.antibiotics
        )
        assert call.mic == 16.0
        assert call.mic_unit == "mg/L"

    def test_unparsable_mic_raises(self, config):
        row = {
            "sample_id": "S1",
            "antibiotic": "imipenem",
            "phenotype": "R",
            "MIC": "high",
        }
        with pytest.raises(PhenotypeError, match="not a number"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    def test_mic_with_nd_raises(self, config):
        """ND means not determined, so it cannot carry an MIC."""
        row = {
            "sample_id": "S1",
            "antibiotic": "imipenem",
            "phenotype": "ND",
            "MIC": "8",
        }
        with pytest.raises(PhenotypeError, match="no determinate category"):
            stage.parse_phenotype_row(
                row, config.allowed_phenotypes, config.antibiotics
            )

    def test_zone_diameter_only_from_measurement(self, config):
        row = {"sample_id": "S1", "antibiotic": "imipenem", "phenotype": "S"}
        call = stage.parse_phenotype_row(
            row, config.allowed_phenotypes, config.antibiotics
        )
        assert call.zone_diameter is None


class TestLoadPhenotype:
    def test_loads_synthetic_fixtures(self, config, phenotype_dir, manifest):
        calls = stage.load_phenotype(
            config, phenotype_dir, "imipenem", sample_ids=manifest.sample_ids
        )
        assert len(calls) == 20
        assert all(c.antibiotic == "imipenem" for c in calls)

    def test_fixtures_exercise_every_allowed_value(self, config, phenotype_dir):
        calls = stage.load_phenotype(config, phenotype_dir, "imipenem")
        seen = {c.phenotype.value for c in calls}
        assert seen == {"R", "I", "S", "SDD", "ND"}

    def test_unconfigured_antibiotic_raises(self, config, phenotype_dir):
        from papipeline.errors import AntibioticError

        with pytest.raises(AntibioticError):
            stage.load_phenotype(config, phenotype_dir, "ceftazidime")

    def test_duplicate_sample_raises(self, config, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "imipenem_phenotype.tsv",
            _table("S1\timipenem\tR\t.\t.\nS1\timipenem\tS\t.\t.\n"),
        )
        with pytest.raises(Exception):
            stage.load_phenotype(config, tmp_path, "imipenem")

    def test_invalid_phenotype_in_file_raises(self, config, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "imipenem_phenotype.tsv", _table("S1\timipenem\tVERY-RESISTANT\t.\t.\n")
        )
        with pytest.raises(PhenotypeError):
            stage.load_phenotype(config, tmp_path, "imipenem")

    def test_unconfigured_antibiotic_in_file_raises(
        self, config, write_tsv_file, tmp_path
    ):
        path = write_tsv_file(
            "imipenem_phenotype.tsv", _table("S1\tceftazidime\tR\t.\t.\n")
        )
        with pytest.raises(PhenotypeError, match="unconfigured antibiotic"):
            stage.load_phenotype(config, tmp_path, "imipenem")

    def test_missing_samples_are_absent_not_defaulted(
        self, config, write_tsv_file, tmp_path
    ):
        """A sample with no record is absent, never defaulted to S."""
        path = write_tsv_file(
            "imipenem_phenotype.tsv", _table("S1\timipenem\tR\t.\t.\n")
        )
        calls = stage.load_phenotype(
            config, tmp_path, "imipenem", sample_ids=["S1", "S2", "S3"]
        )
        assert {c.sample_id for c in calls} == {"S1"}

    def test_missing_file_raises(self, config, tmp_path):
        with pytest.raises(Exception):
            stage.load_phenotype(config, tmp_path, "imipenem")


class TestBinarise:
    """Mapping a category to 1/0 for GWAS must never guess."""

    @staticmethod
    def _call(value: str):
        from papipeline.models import PhenotypeCall

        return PhenotypeCall(
            sample_id="S1", antibiotic="imipenem", phenotype=Phenotype(value)
        )

    def test_resistant_maps_to_one(self):
        assert stage.binarise(self._call("R"), ["R"], ["S"]) == 1

    def test_susceptible_maps_to_zero(self):
        assert stage.binarise(self._call("S"), ["R"], ["S"]) == 0

    @pytest.mark.parametrize("value", ["I", "SDD", "ND"])
    def test_intermediate_categories_are_excluded_not_guessed(self, value):
        """I must never be silently folded into R or S."""
        assert stage.binarise(self._call(value), ["R"], ["S"]) is None

    def test_unlisted_category_is_excluded(self):
        assert stage.binarise(self._call("SDD"), ["R"], []) is None


class TestPhenotypeMatrix:
    def test_long_format_matrix(self, config, phenotype_dir):
        """The reported view.

        `phenotype_matrix` is now the CONTINUOUS trait (log2 MIC), which the
        association model consumes. The categorical view still has to exist for
        the report, so it lives under its own name rather than being lost when
        the trait changed. See spec decision D4a.
        """
        calls = stage.load_phenotype(config, phenotype_dir, "imipenem")
        matrix = stage.phenotype_category_matrix(calls)
        assert len(matrix) == 20
        assert all("imipenem" in v for v in matrix.values())

    def test_matrix_keys_are_sample_ids(self, config, phenotype_dir):
        calls = stage.load_phenotype(config, phenotype_dir, "imipenem")
        matrix = stage.phenotype_category_matrix(calls)
        assert all(k.startswith("TEST_PA_") for k in matrix)
