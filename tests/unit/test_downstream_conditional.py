"""Deliverable 2: the conditional and stratified scans.

The conditional scan fits every candidate feature with the known layer features
as covariates; a feature that IS one of the covariates is reported as
``skipped_covariate`` rather than silently dropped. The stratified scan is the
same association machinery restricted to the isolates that lack every known
determinant, and it refuses to run when that stratum is smaller than the
configured group minimum.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import conditional as cd
from papipeline.downstream import known_determinants as kd
from papipeline.errors import DataContractError
from papipeline.manifest import SampleManifest
from papipeline.models import Phenotype, PhenotypeCall, Sample
from papipeline.stages import gwas as stage

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture(scope="module")
def known(config) -> kd.KnownDeterminants:
    return kd.load_known_determinants(config)


def _make_input(config, phenotypes: dict, columns: dict, lineages=None):
    """Build a ``GwasInput`` from plain dicts, reusing stage 9's loader."""
    sample_ids = sorted(phenotypes)
    calls = [
        PhenotypeCall(
            sample_id=s, antibiotic="imipenem", phenotype=Phenotype(phenotypes[s])
        )
        for s in sample_ids
    ]
    manifest = SampleManifest([Sample(s, None, "test") for s in sample_ids])
    rows = [
        {"sample_id": s, **{c: columns[c].get(s, 0) for c in sorted(columns)}}
        for s in sample_ids
    ]
    return stage.build_input(config, manifest, calls, rows, lineages=lineages)


def _balanced_phenotypes(n: int) -> dict:
    half = n // 2
    return {f"S{i:04d}": ("R" if i < half else "S") for i in range(n)}


class TestCovariateSelection:
    def test_only_known_layer_features_become_covariates(self, config, known):
        phenotypes = _balanced_phenotypes(60)
        columns = {
            "gene__oprD_absent": {s: int(i % 2) for i, s in enumerate(sorted(phenotypes))},
            "gene_presence_absence__blaKPC-2": {
                s: int(i % 3 == 0) for i, s in enumerate(sorted(phenotypes))
            },
            "gene__femA_novel": {s: int(i % 5 == 0) for i, s in enumerate(sorted(phenotypes))},
        }
        gi = _make_input(config, phenotypes, columns)
        covariates = cd.covariate_columns(gi, known, config=config)
        assert "gene__oprD_absent" in covariates
        assert "gene_presence_absence__blaKPC-2" in covariates
        assert "gene__femA_novel" not in covariates

    def test_extra_covariates_come_from_the_config_key(self, config, known):
        phenotypes = _balanced_phenotypes(60)
        columns = {
            "gene__femA_novel": {s: int(i % 2) for i, s in enumerate(sorted(phenotypes))},
            "gene__oprD_absent": {s: int(i % 3) for i, s in enumerate(sorted(phenotypes))},
        }
        gi = _make_input(config, phenotypes, columns)
        tuned = replace(
            config, raw={**config.raw, "downstream": {"covariates": ["gene__femA_novel"]}}
        )
        covariates = cd.covariate_columns(gi, known, config=tuned)
        assert "gene__femA_novel" in covariates
        assert "gene__oprD_absent" in covariates

    def test_an_extra_covariate_missing_from_the_matrix_is_named(self, config, known):
        phenotypes = _balanced_phenotypes(60)
        columns = {"gene__oprD_absent": {s: 0 for s in phenotypes}}
        gi = _make_input(config, phenotypes, columns)
        tuned = replace(
            config, raw={**config.raw, "downstream": {"covariates": ["gene__ghost"]}}
        )
        with pytest.raises(DataContractError) as exc:
            cd.covariate_columns(gi, known, config=tuned)
        assert "downstream.covariates" in str(exc.value)
        assert "gene__ghost" in str(exc.value)


class TestConditionalScan:
    @pytest.fixture(scope="class")
    def fitted(self, config, known):
        phenotypes = _balanced_phenotypes(60)
        sample_ids = sorted(phenotypes)
        strong = {
            s: int(phenotypes[s] == "R" if i > 1 else 0) for i, s in enumerate(sample_ids)
        }
        noise = {s: i % 2 for i, s in enumerate(sample_ids)}
        rare = {s: int(i == 0) for i, s in enumerate(sample_ids)}
        columns = {
            "gene__strong": strong,
            "gene__noise": noise,
            "gene__rare": rare,
            "gene__oprD_absent": noise,
        }
        gi = _make_input(config, phenotypes, columns)
        covariates = cd.covariate_columns(gi, known, config=config)
        results = cd.conditional_scan(gi, covariates)
        return gi, covariates, {r.feature: r for r in results}

    def test_a_feature_that_is_a_covariate_is_skipped_by_name(self, fitted):
        _, _, by_feature = fitted
        assert by_feature["gene__oprD_absent"].status == "skipped_covariate"
        assert by_feature["gene__oprD_absent"].p_value is None

    def test_a_still_associated_feature_reports_adjusted_p(self, fitted):
        _, covariates, by_feature = fitted
        result = by_feature["gene__strong"]
        assert result.status == "tested"
        assert covariates == result.covariates
        assert result.p_value is not None
        assert result.adjusted_p_value is not None
        assert result.adjusted_p_value >= result.p_value
        assert result.adjusted_p_value < 0.05

    def test_a_single_carrier_feature_is_reported_not_testable(self, fitted):
        _, _, by_feature = fitted
        assert by_feature["gene__rare"].status == "not_testable"

    def test_the_marginal_p_of_a_nulled_feature_is_not_significant(self, fitted):
        _, _, by_feature = fitted
        result = by_feature["gene__noise"]
        assert result.status == "tested"
        assert result.adjusted_p_value > 0.05


class TestStratifiedScan:
    def test_only_isolates_lacking_every_known_column_survive(self, config, known):
        phenotypes = _balanced_phenotypes(60)
        sample_ids = sorted(phenotypes)
        columns = {
            "gene_presence_absence__blaKPC-2": {s: int(i < 6) for i, s in enumerate(sample_ids)},
            "gene__femA_novel": {
                s: int(10 <= i < 16) for i, s in enumerate(sample_ids)
            },
        }
        gi = _make_input(config, phenotypes, columns)
        scan = cd.stratified_scan(
            gi, ["gene_presence_absence__blaKPC-2"], config
        )
        carriers = [s for s in sample_ids if columns["gene_presence_absence__blaKPC-2"][s]]
        assert set(scan.samples).isdisjoint(carriers)
        assert set(scan.samples) == set(sample_ids) - set(carriers)
        assert scan.results
        assert all(r.model.startswith("reference_fisher") for r in scan.results)

    def test_the_stratified_subset_keeps_enough_of_both_groups(self, config):
        phenotypes = _balanced_phenotypes(60)
        sample_ids = sorted(phenotypes)
        columns = {"gene_presence_absence__blaKPC-2": {s: int(i < 6) for i, s in enumerate(sample_ids)}}
        gi = _make_input(config, phenotypes, columns)
        scan = cd.stratified_scan(gi, ["gene_presence_absence__blaKPC-2"], config)
        # 60 - 6 = 54, split 24 R / 30 S.
        assert len(scan.samples) == 54
        assert sum(1 for s in scan.samples if phenotypes[s] == "R") == 24

    def test_a_stratum_below_the_group_minimum_is_refused_with_the_key(self, config):
        phenotypes = _balanced_phenotypes(60)
        sample_ids = sorted(phenotypes)
        columns = {"gene_presence_absence__blaKPC-2": {s: int(i < 45) for i, s in enumerate(sample_ids)}}
        gi = _make_input(config, phenotypes, columns)
        with pytest.raises(DataContractError) as exc:
            cd.stratified_scan(gi, ["gene_presence_absence__blaKPC-2"], config)
        assert "gwas.min_samples_per_group" in str(exc.value)

    def test_lineage_unknown_is_excluded_from_the_stratum_definition(self, config):
        """The stratum is defined by known determinants only; lineages do not enter."""
        phenotypes = _balanced_phenotypes(60)
        sample_ids = sorted(phenotypes)
        columns = {"gene__femA_novel": {s: 0 for s in sample_ids}}
        lineages = {s: ("unknown" if i % 2 else "ST1") for i, s in enumerate(sample_ids)}
        gi = _make_input(config, phenotypes, columns, lineages=lineages)
        scan = cd.stratified_scan(gi, [], config)
        assert set(scan.samples) == set(sample_ids)
