"""Tests for sample ID validation and the master manifest."""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError, SampleIdError
from papipeline.manifest import SampleManifest, discover_manifest, validate_sample_id
from papipeline.models import Sample


class TestValidateSampleId:
    @pytest.mark.parametrize(
        "sample_id",
        ["TEST_PA_001", "S1", "sample-1", "sample_1", "GCA_000710625.1", "a" * 128],
    )
    def test_accepts_valid_identifiers(self, sample_id):
        assert validate_sample_id(sample_id) == sample_id

    @pytest.mark.parametrize("sample_id", [None, "", "   "])
    def test_rejects_empty(self, sample_id):
        with pytest.raises(SampleIdError):
            validate_sample_id(sample_id)

    def test_rejects_surrounding_whitespace(self):
        """Whitespace is not silently stripped; it signals a real problem."""
        with pytest.raises(SampleIdError, match="whitespace"):
            validate_sample_id(" S1 ")

    @pytest.mark.parametrize("sample_id", ["a", "a" * 129])
    def test_rejects_bad_length(self, sample_id):
        with pytest.raises(SampleIdError, match="length"):
            validate_sample_id(sample_id)

    @pytest.mark.parametrize(
        "sample_id", ["a/b", "a\\b", "a b", "a\tb", "a;b", "a,b", "a|b", "a#b", "a'b"]
    )
    def test_rejects_forbidden_characters(self, sample_id):
        with pytest.raises(SampleIdError):
            validate_sample_id(sample_id)

    @pytest.mark.parametrize("sample_id", ["sample!", "sample@1", "a$b", "a%b"])
    def test_rejects_non_alphanumeric(self, sample_id):
        with pytest.raises(SampleIdError):
            validate_sample_id(sample_id)


class TestSampleManifest:
    def test_accepts_distinct_samples(self, sample_manifest):
        assert len(sample_manifest) == 3
        assert sample_manifest.sample_ids == ["TEST_A_01", "TEST_A_02", "TEST_A_03"]

    def test_rejects_duplicate_sample_ids(self):
        with pytest.raises(DataContractError, match="Duplicate sample_id"):
            SampleManifest([Sample("S1"), Sample("S1")])

    def test_rejects_invalid_sample_id(self):
        with pytest.raises(SampleIdError):
            SampleManifest([Sample("bad id")])

    def test_get_and_require(self, sample_manifest):
        assert sample_manifest.get("TEST_A_01").sample_id == "TEST_A_01"
        assert sample_manifest.get("NOPE") is None
        with pytest.raises(DataContractError, match="not present"):
            sample_manifest.require("NOPE")

    def test_index(self, sample_manifest):
        assert set(sample_manifest.index()) == set(sample_manifest.sample_ids)

    def test_iteration(self, sample_manifest):
        assert [s.sample_id for s in sample_manifest] == sample_manifest.sample_ids


class TestDiscoverManifest:
    def test_loads_synthetic_fixtures(self, manifest):
        assert len(manifest) == 20
        assert manifest.sample_ids[0] == "TEST_PA_001"
        assert manifest.sample_ids[-1] == "TEST_PA_020"

    def test_every_sample_has_an_assembly_path(self, manifest):
        assert all(s.assembly_path for s in manifest)

    def test_assembly_files_exist(self, manifest):
        from pathlib import Path

        missing = [
            s.sample_id
            for s in manifest
            if not Path(s.assembly_path).exists()
        ]
        assert missing == []

    def test_missing_metadata_raises(self, tmp_path):
        with pytest.raises(DataContractError, match="metadata file not found"):
            discover_manifest(tmp_path)

    def test_duplicate_metadata_rows_raise(self, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "sample_metadata.tsv",
            "sample_id\tassembly_path\nS1\ta.fna\nS1\tb.fna\n",
        )
        with pytest.raises(DataContractError):
            discover_manifest(tmp_path)

    def test_invalid_sample_id_raises(self, write_tsv_file, tmp_path):
        path = write_tsv_file("sample_metadata.tsv", "sample_id\nbad id\n")
        with pytest.raises(SampleIdError):
            discover_manifest(tmp_path)

    def test_require_files_rejects_absent_assembly(self, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "sample_metadata.tsv", "sample_id\tassembly_path\nS1\tmissing.fna\n"
        )
        with pytest.raises(DataContractError, match="does not exist"):
            discover_manifest(tmp_path, require_files=True)

    def test_tolerates_absent_assembly_by_default(self, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "sample_metadata.tsv", "sample_id\tassembly_path\nS1\tmissing.fna\n"
        )
        manifest = discover_manifest(tmp_path)
        assert len(manifest) == 1
