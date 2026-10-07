"""Deliverable 3: interaction tiers.

Tier 1 is the pre-specified pair list in ``config/interaction_pairs.tsv`` -
for carbapenems: PDC-high with oprD_off_any, KPC with oprD_off_any, nalC with
mexR, mexZ with ampD - each fitted as a logistic model with an explicit
interaction term. Tier 2 is the anchored scan (known x novel) with FDR across
the fitted interaction p-values. Tier 3 is skipped by design and refused if
somebody writes a tier-3 row into the table.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import interactions as ix
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


def _make_input(config, phenotypes: dict, columns: dict):
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
    return stage.build_input(config, manifest, calls, rows)


def _synergy_cohort() -> tuple[dict, dict]:
    """Four balanced cells; the double-positive cell is mostly resistant.

    a=0/b=0: 10 R + 40 S; a=1/b=0: 10 R + 40 S; a=0/b=1: 10 R + 40 S;
    a=1/b=1: 40 R + 10 S.  A null main-effect-only model cannot produce this;
    the interaction term must.
    """
    phenotypes, a, b = {}, {}, {}
    index = 0
    for ca, cb, n_r, n_s in (
        (0, 0, 10, 40),
        (1, 0, 10, 40),
        (0, 1, 10, 40),
        (1, 1, 40, 10),
    ):
        for label in ["R"] * n_r + ["S"] * n_s:
            s = f"S{index:04d}"
            index += 1
            phenotypes[s] = label
            a[s] = ca
            b[s] = cb
    return phenotypes, {"gene__A": a, "gene__B": b}


@pytest.fixture(scope="module")
def config_pairs(config):
    return ix.load_interaction_pairs(
        ix.interaction_pairs_path(config), known=kd.load_known_determinants(config)
    )


class TestTierOneTable:
    def test_four_pre_specified_pairs_load(self, config_pairs):
        assert len(config_pairs) == 4
        assert all(p.tier == 1 for p in config_pairs)

    def test_tier_one_scan_fits_the_imipenem_pairs(self, config, config_pairs):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        results = ix.tier1_scan(gi, config_pairs, antibiotic="imipenem")
        assert len(results) == 4
        assert all(r.status == "tested" for r in results)
        assert all(r.adjusted_p_value is not None for r in results)

    def test_a_pair_for_another_antibiotic_does_not_run(self, config, config_pairs):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        assert ix.tier1_scan(gi, config_pairs, antibiotic="ciprofloxacin") == []


class TestTierOneFit:
    def test_the_interaction_term_recovers_synergy(self, config):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        pair = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__A", feature_b="gene__B"
        )
        (result,) = ix.tier1_scan(gi, [pair], antibiotic="imipenem")
        assert result.status == "tested"
        assert result.n_samples == 200
        assert result.interaction_estimate > 0
        assert result.interaction_p_value < 0.01
        assert result.interaction_std_error > 0
        assert result.main_a_p_value is not None
        assert result.main_b_p_value is not None

    def test_a_pair_with_a_missing_column_is_reported_not_raised(self, config):
        phenotypes, columns = _synergy_cohort()
        columns = dict(columns)
        gi = _make_input(config, phenotypes, columns)
        pair = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__ghost", feature_b="gene__B"
        )
        (result,) = ix.tier1_scan(gi, [pair], antibiotic="imipenem")
        assert result.status == "missing_feature"
        assert result.interaction_p_value is None

    def test_fdr_is_applied_over_the_fitted_pairs_only(self, config):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        good = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__A", feature_b="gene__B"
        )
        ghost = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__ghost", feature_b="gene__B"
        )
        results = ix.tier1_scan(gi, [good, ghost], antibiotic="imipenem")
        fitted = [r for r in results if r.status == "tested"]
        assert len(fitted) == 1
        assert fitted[0].adjusted_p_value == pytest.approx(fitted[0].p_value)


class TestTierTwo:
    def test_anchors_pair_against_non_anchors_only(self, config):
        phenotypes = {f"S{i:04d}": ("R" if i < 40 else "S") for i in range(80)}
        sample_ids = sorted(phenotypes)
        columns = {
            "gene__oprD_absent": {s: int(i % 2) for i, s in enumerate(sample_ids)},
            "gene_presence_absence__blaKPC-2": {s: int(i % 3 == 0) for i, s in enumerate(sample_ids)},
            "gene__femA_novel": {s: int(i % 4 == 0) for i, s in enumerate(sample_ids)},
        }
        gi = _make_input(config, phenotypes, columns)
        known = kd.load_known_determinants(config)
        results = ix.tier2_anchored_scan(gi, known, max_pairs=1000)
        assert results
        anchors = set(
            cd_covariates(columns)
        )
        for r in results:
            assert r.tier == 2
            assert r.feature_a in anchors, "anchor must be a known feature"
            assert r.feature_b not in anchors, "novel x novel is not a tier-2 pair"
        endpoints = {(r.feature_a, r.feature_b) for r in results}
        assert ("gene__oprD_absent", "gene__femA_novel") in endpoints

    def test_the_cap_is_checked_against_the_config_key(self, config):
        phenotypes = {f"S{i:04d}": ("R" if i < 40 else "S") for i in range(80)}
        sample_ids = sorted(phenotypes)
        columns = {
            "gene__oprD_absent": {s: int(i % 2) for i, s in enumerate(sample_ids)},
            "gene__femA_novel": {s: int(i % 4 == 0) for i, s in enumerate(sample_ids)},
            "gene__femB_novel": {s: int(i % 5 == 0) for i, s in enumerate(sample_ids)},
        }
        gi = _make_input(config, phenotypes, columns)
        known = kd.load_known_determinants(config)
        with pytest.raises(DataContractError) as exc:
            ix.tier2_anchored_scan(gi, known, max_pairs=1)
        assert "downstream.tier2_max_pairs" in str(exc.value)

    def test_the_config_cap_is_read_through_the_documented_key(self, config):
        assert ix.resolve_tier2_max_pairs(config) >= 1
        tuned = replace(
            config, raw={**config.raw, "downstream": {"tier2_max_pairs": 7}}
        )
        assert ix.resolve_tier2_max_pairs(tuned) == 7

    def test_tier_two_p_values_are_fdr_controlled(self, config):
        phenotypes = {f"S{i:04d}": ("R" if i < 40 else "S") for i in range(80)}
        sample_ids = sorted(phenotypes)
        columns = {
            "gene__oprD_absent": {s: int(i % 2) for i, s in enumerate(sample_ids)},
            "gene__femA_novel": {s: int(i % 4 == 0) for i, s in enumerate(sample_ids)},
            "gene__femB_novel": {s: int(i % 5 == 0) for i, s in enumerate(sample_ids)},
        }
        gi = _make_input(config, phenotypes, columns)
        known = kd.load_known_determinants(config)
        results = ix.tier2_anchored_scan(gi, known, max_pairs=1000)
        tested = [r for r in results if r.status == "tested"]
        assert len(tested) >= 2
        for r in tested:
            assert r.adjusted_p_value >= r.p_value


def cd_covariates(columns: dict) -> set[str]:
    """Known columns of a hand-built matrix (kept local to this test file)."""
    known = kd.load_known_determinants(load_config(SCIENCE, machine="laptop"))
    return {c for c in columns if known.is_known(c)}
