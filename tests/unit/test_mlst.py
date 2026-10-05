"""Tests for the MLST parser (stage 3).

The recurring theme: an incomplete profile must be reported as incomplete,
never as a guessed sequence type.
"""

from __future__ import annotations

import pytest

from papipeline.stages import mlst as stage
from papipeline.models import MlstStatus

HEADER = "sample_id\tST\talleles\tMLST_status\tmlst_scheme\tallele_database\n"
FULL_ALLELES = "abcZ:1;adk:2;aroE:3;guaA:4;gyrB:5;ircD:6;toxA:7"
PARTIAL_ALLELES = "abcZ:1;adk:2;aroE:3"


class TestParseAlleles:
    def test_parses_semicolon_pairs(self):
        assert stage.parse_alleles(FULL_ALLELES)["abcZ"] == "1"

    def test_accepts_comma_separator(self):
        assert stage.parse_alleles("abcZ:1,adk:2") == {"abcZ": "1", "adk": "2"}

    def test_empty_gives_empty_dict(self):
        assert stage.parse_alleles("") == {}
        assert stage.parse_alleles(None) == {}

    def test_skips_malformed_token(self):
        assert stage.parse_alleles("abcZ:1;garbage;adk:2") == {"abcZ": "1", "adk": "2"}

    def test_skips_incomplete_pair(self):
        assert stage.parse_alleles("abcZ:1;:2;adk:3") == {"abcZ": "1", "adk": "3"}

    @pytest.mark.parametrize("allele", ["0", "NA", "-", "None", "none"])
    def test_null_sentinel_allele_dropped(self, allele):
        """An unmatched locus is dropped, not recorded as allele 0."""
        assert stage.parse_alleles(f"abcZ:{allele}") == {}

    def test_allele_value_zero_dropped(self):
        assert stage.parse_alleles("abcZ:0") == {}


class TestNormaliseScheme:
    def test_alias_resolved(self, config):
        assert stage.normalise_scheme(config, "paeruginosa") == "pseudomonas_aeruginosa"

    def test_canonical_name_unchanged(self, config):
        assert (
            stage.normalise_scheme(config, "pseudomonas_aeruginosa")
            == "pseudomonas_aeruginosa"
        )

    def test_missing_scheme_uses_fallback(self, config):
        assert stage.normalise_scheme(config, None) == "unknown"


class TestParseRow:
    def test_complete_profile_is_typed(self, config):
        call = stage.parse_row(
            {
                "sample_id": "S1",
                "ST": "12",
                "alleles": FULL_ALLELES,
                "MLST_status": "typed",
            },
            config,
        )
        assert call.sequence_type == "12"
        assert call.mlst_status == MlstStatus.TYPED.value
        assert len(call.alleles) == 7

    def test_partial_profile_keeps_alleles(self, config):
        call = stage.parse_row(
            {
                "sample_id": "S1",
                "ST": "NA",
                "alleles": PARTIAL_ALLELES,
                "MLST_status": "partial",
            },
            config,
        )
        assert call.sequence_type is None
        assert call.mlst_status == MlstStatus.PARTIAL.value
        assert len(call.alleles) == 3

    def test_typed_but_too_few_loci_is_demoted(self, config):
        """Fewer loci than the minimum is partial, whatever the row claims."""
        call = stage.parse_row(
            {
                "sample_id": "S1",
                "ST": "12",
                "alleles": PARTIAL_ALLELES,
                "MLST_status": "typed",
            },
            config,
        )
        assert call.mlst_status == MlstStatus.PARTIAL.value
        assert call.sequence_type is None

    def test_typed_without_st_is_demoted(self, config):
        call = stage.parse_row(
            {
                "sample_id": "S1",
                "ST": ".",
                "alleles": FULL_ALLELES,
                "MLST_status": "typed",
            },
            config,
        )
        assert call.mlst_status == MlstStatus.PARTIAL.value
        assert call.sequence_type is None

    def test_no_alleles_becomes_no_call(self, config):
        call = stage.parse_row(
            {"sample_id": "S1", "ST": ".", "alleles": ".", "MLST_status": "typed"},
            config,
        )
        assert call.mlst_status == MlstStatus.NO_CALL.value

    def test_unknown_status_becomes_no_call(self, config):
        call = stage.parse_row(
            {
                "sample_id": "S1",
                "ST": "12",
                "alleles": FULL_ALLELES,
                "MLST_status": "weird",
            },
            config,
        )
        assert call.mlst_status == MlstStatus.NO_CALL.value

    def test_allele_string_round_trip(self, config):
        call = stage.parse_row(
            {"sample_id": "S1", "ST": "12", "alleles": FULL_ALLELES, "MLST_status": "typed"},
            config,
        )
        assert stage.parse_alleles(call.allele_string()) == call.alleles

    def test_missing_sample_id_raises(self, config):
        with pytest.raises(ValueError, match="no sample_id"):
            stage.parse_row({"ST": "1", "alleles": FULL_ALLELES}, config)


class TestLoadMlst:
    def test_loads_synthetic_fixtures(self, config, intermediate_root):
        calls = stage.load_mlst(
            config, intermediate_root / "mlst" / "mlst_results.tsv"
        )
        assert len(calls) == 20

    def test_fixtures_include_typed_and_partial(self, config, intermediate_root):
        calls = stage.load_mlst(
            config, intermediate_root / "mlst" / "mlst_results.tsv"
        )
        statuses = {c.mlst_status for c in calls.values()}
        assert "typed" in statuses
        assert "partial" in statuses

    def test_partial_calls_have_no_sequence_type(self, config, intermediate_root):
        calls = stage.load_mlst(
            config, intermediate_root / "mlst" / "mlst_results.tsv"
        )
        for call in calls.values():
            if call.mlst_status == MlstStatus.PARTIAL.value:
                assert call.sequence_type is None

    def test_duplicate_sample_raises(self, config, write_tsv_file):
        path = write_tsv_file(
            "dup.tsv",
            HEADER + f"S1\t1\t{FULL_ALLELES}\ttyped\tx\ty\nS1\t2\t{FULL_ALLELES}\ttyped\tx\ty\n",
        )
        with pytest.raises(Exception):
            stage.load_mlst(config, path)

    def test_missing_required_column_raises(self, config, write_tsv_file):
        path = write_tsv_file("bad.tsv", "sample_id\tST\nS1\t1\n")
        with pytest.raises(Exception):
            stage.load_mlst(config, path)


class TestAlignToManifest:
    def test_untyped_sample_maps_to_none(self, config, intermediate_root, manifest):
        calls = stage.load_mlst(
            config, intermediate_root / "mlst" / "mlst_results.tsv"
        )
        aligned = stage.align_to_manifest(calls, manifest)
        assert set(aligned) == set(manifest.sample_ids)
        assert all(v is not None for v in aligned.values())

    def test_absent_sample_maps_to_none(self, config, manifest, sample_manifest):
        aligned = stage.align_to_manifest({}, sample_manifest)
        assert list(aligned.values()) == [None, None, None]
