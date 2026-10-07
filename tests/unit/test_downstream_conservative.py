"""The conservative path in the tier-1 fit, pinned so it cannot drift.

Two situations make a logistic model with an interaction term uninformative,
and both must produce NO evidence rather than manufactured evidence:

1. **A pre-specified pair whose columns are not in this matrix.** The tier-1
   table is a commitment: every pair in ``config/interaction_pairs.tsv`` has to
   appear in the output of a tier-1 scan, even when the layer encoder did not
   produce one of its columns in this run. Such a row is reported
   ``tested`` with p = 1.0 - a null result, never a significant one - and the
   note says why. A pair that is NOT pre-specified (a locally constructed one)
   with a missing column is reported ``missing_feature`` instead, so the
   distinction is the table, not the column.

2. **A product column that is identically zero or collinear.** The interaction
   term carries no information; the likelihood-ratio answer is p = 1.0. The
   main effects are still fitted, because they are estimable.

Both rules are one-sided: they can only ever weaken a claim.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import interactions as ix
from papipeline.downstream import known_determinants as kd
from papipeline.manifest import SampleManifest
from papipeline.models import Phenotype, PhenotypeCall, Sample
from papipeline.stages import gwas as stage

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture(scope="module")
def config_pairs(config):
    return ix.load_interaction_pairs(
        ix.interaction_pairs_path(config), known=kd.load_known_determinants(config)
    )


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
    """The cohort the tier-1 spec test uses: 200 isolates, synergy in cell 4."""
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


class TestDeclaredPairsWithoutColumns:
    def test_a_pre_specified_pair_is_reported_tested_with_no_evidence(
        self, config, config_pairs
    ):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        results = ix.tier1_scan(gi, config_pairs, antibiotic="imipenem")
        assert len(results) == 4
        for result in results:
            assert result.status == "tested"
            assert result.p_value == 1.0
            assert result.interaction_p_value == 1.0
            assert result.interaction_estimate == 0.0
            assert result.adjusted_p_value is not None
            assert result.adjusted_p_value >= result.p_value
            assert result.n_samples == gi.n_samples
            assert result.note

    def test_the_conservative_path_can_only_report_a_null_result(self, config, config_pairs):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        results = ix.tier1_scan(gi, config_pairs, antibiotic="imipenem")
        threshold = config.gwas.significance_threshold
        assert all(r.p_value > threshold for r in results)
        assert all(r.adjusted_p_value > threshold for r in results)

    def test_an_undeclared_pair_with_a_missing_column_is_still_missing_feature(
        self, config
    ):
        phenotypes, columns = _synergy_cohort()
        gi = _make_input(config, phenotypes, columns)
        pair = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__ghost", feature_b="gene__B"
        )
        (result,) = ix.tier1_scan(gi, [pair], antibiotic="imipenem")
        assert result.status == "missing_feature"
        assert result.p_value is None
        assert result.interaction_p_value is None

    def test_a_declared_flag_comes_from_the_table_not_the_matrix(self, config):
        """Loading with the knowledge table is what marks a pair pre-specified."""
        known = kd.load_known_determinants(config)
        pairs = ix.load_interaction_pairs(ix.interaction_pairs_path(config), known=known)
        assert all(p.declared is True for p in pairs)
        plain = ix.load_interaction_pairs(ix.interaction_pairs_path(config))
        assert all(p.declared is None for p in plain)


class TestDegenerateInteractionTerm:
    @pytest.fixture(scope="class")
    def disjoint_cohort(self, config):
        """Two mutually exclusive features: the product column is all zeros."""
        phenotypes = {f"S{i:04d}": ("R" if i < 40 else "S") for i in range(80)}
        sample_ids = sorted(phenotypes)
        a = {s: int(i % 4 in (0, 1)) for i, s in enumerate(sample_ids)}
        b = {s: int(i % 4 in (2, 3)) for i, s in enumerate(sample_ids)}
        assert all(a[s] * b[s] == 0 for s in sample_ids)
        gi = _make_input(config, phenotypes, {"gene__A": a, "gene__B": b})
        pair = ix.InteractionPair(
            tier=1, antibiotic="imipenem", feature_a="gene__A", feature_b="gene__B"
        )
        (result,) = ix.tier1_scan(gi, [pair], antibiotic="imipenem")
        return gi, result

    def test_the_pair_is_tested_but_the_interaction_is_not_estimable(self, disjoint_cohort):
        gi, result = disjoint_cohort
        assert result.status == "tested"
        assert result.interaction_p_value == 1.0
        assert result.interaction_estimate == 0.0
        assert result.interaction_std_error is None
        assert result.n_samples == gi.n_samples

    def test_the_main_effects_are_still_estimated(self, disjoint_cohort):
        _gi, result = disjoint_cohort
        assert result.main_a_p_value is not None
        assert result.main_b_p_value is not None
        assert 0.0 <= result.main_a_p_value <= 1.0
        assert 0.0 <= result.main_b_p_value <= 1.0

    def test_a_degenerate_pair_is_still_in_the_fdr_family(self, disjoint_cohort):
        _gi, result = disjoint_cohort
        assert result.adjusted_p_value is not None
        assert result.adjusted_p_value >= result.interaction_p_value
