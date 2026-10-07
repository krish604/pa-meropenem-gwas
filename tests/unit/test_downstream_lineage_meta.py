"""Deliverable 4: the lineage-stratified meta-analysis.

The cohort association is pooled ACROSS lineages, not within one of them, so
the feature's effect is estimated separately inside each lineage and the
per-lineage estimates are combined by fixed-effect inverse variance. Cochran's
Q says how much of the variation between lineages is more than sampling error
would produce, and I squared says what share of it that is.

Every expected value below is recomputed in this file from the textbook
formulas with ``math`` only - deliberately not by calling the code under test,
so a change in the implementation cannot quietly redefine the answer.

Models carry the repository's ``:NO_KINSHIP_CORRECTION`` stamp: a per-lineage
stratum is a subset of the same structured cohort, and the reference test does
not correct for kinship.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import lineage_meta as lm
from papipeline.errors import DataContractError
from papipeline.manifest import SampleManifest
from papipeline.models import Phenotype, PhenotypeCall, Sample
from papipeline.stages import gwas as stage

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


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


def _two_lineage_tables() -> dict:
    """`(lineage, carries, phenotype)` for 40 isolates, two lineages.

    ST1: carriers 8 R / 2 S, non-carriers 2 R / 8 S.
    ST2: carriers 4 R / 6 S, non-carriers 6 R / 4 S.
    """
    rows = []
    rows += [("ST1", 1, p) for p in ["R"] * 8 + ["S"] * 2]
    rows += [("ST1", 0, p) for p in ["R"] * 2 + ["S"] * 8]
    rows += [("ST2", 1, p) for p in ["R"] * 4 + ["S"] * 6]
    rows += [("ST2", 0, p) for p in ["R"] * 6 + ["S"] * 4]
    return {"rows": rows, "feature": "gene__femA_novel"}


def _input_from_rows(config, rows, feature: str, lineages: dict | None = None):
    phenotypes = {f"S{i:04d}": p for i, (_lineage, _carries, p) in enumerate(rows)}
    columns = {feature: {f"S{i:04d}": carries for i, (_l, carries, _p) in enumerate(rows)}}
    labels = (
        lineages
        if lineages is not None
        else {f"S{i:04d}": lineage for i, (lineage, _c, _p) in enumerate(rows)}
    )
    return _make_input(config, phenotypes, columns, lineages=labels)


# --------------------------------------------------------------------------
# The formulas this file checks against, written out rather than imported.
# --------------------------------------------------------------------------


def _lnor(carrier_r: float, carrier_s: float, non_r: float, non_s: float) -> float:
    """Log odds ratio of the 2x2 table rows=outcome, cols=carriage."""
    return math.log((carrier_r * non_s) / (non_r * carrier_s))


def _lnor_se(carrier_r: float, carrier_s: float, non_r: float, non_s: float) -> float:
    return math.sqrt(
        1 / carrier_r + 1 / carrier_s + 1 / non_r + 1 / non_s
    )


def _two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def _chi2_sf_1(q: float) -> float:
    """Survival function of chi-squared with one degree of freedom."""
    return math.erfc(math.sqrt(q / 2.0))


@pytest.fixture(scope="module")
def two_lineages(config):
    table = _two_lineage_tables()
    gi = _input_from_rows(config, table["rows"], table["feature"])
    return gi, table["feature"]


class TestStratumEstimates:
    def test_each_stratum_reports_the_textbook_log_odds_ratio(self, config, two_lineages):
        gi, feature = two_lineages
        (st1, st2) = lm.stratum_estimates(gi, feature)
        assert [e.lineage for e in (st1, st2)] == ["ST1", "ST2"]

        # ST1: 8 R carriers, 2 S carriers, 2 R non-carriers, 8 S non-carriers.
        expected_st1 = _lnor(8, 2, 2, 8)
        assert st1.effect == pytest.approx(expected_st1)
        assert st1.std_error == pytest.approx(_lnor_se(8, 2, 2, 8))
        assert st1.p_value == pytest.approx(_two_sided_p(expected_st1 / _lnor_se(8, 2, 2, 8)))
        assert st1.corrected is False
        assert st1.n_samples == 20

        # ST2: 4 R carriers, 6 S carriers, 6 R non-carriers, 4 S non-carriers.
        expected_st2 = _lnor(4, 6, 6, 4)
        assert st2.effect == pytest.approx(expected_st2)
        assert st2.std_error == pytest.approx(_lnor_se(4, 6, 6, 4))
        assert st2.p_value == pytest.approx(_two_sided_p(expected_st2 / _lnor_se(4, 6, 6, 4)))
        assert st2.corrected is False

    def test_a_zero_cell_gets_the_haldane_anscombe_half(self, config):
        """A lineage whose carriers are all resistant still yields an estimate."""
        rows = []
        rows += [("ST1", 1, p) for p in ["R"] * 8 + ["S"] * 2]
        rows += [("ST1", 0, p) for p in ["R"] * 2 + ["S"] * 8]
        rows += [("ST2", 1, p) for p in ["R"] * 5]
        rows += [("ST2", 0, p) for p in ["R"] * 3 + ["S"] * 7]
        gi = _input_from_rows(config, rows, "gene__femA_novel")
        estimates = {e.lineage: e for e in lm.stratum_estimates(gi, "gene__femA_novel")}

        st2 = estimates["ST2"]
        assert st2.corrected is True
        corrected = (5.5, 0.5, 3.5, 7.5)  # carrier_r, carrier_s, non_r, non_s
        assert st2.effect == pytest.approx(_lnor(*corrected))
        assert st2.std_error == pytest.approx(_lnor_se(*corrected))

    def test_a_stratum_with_one_outcome_class_is_dropped(self, config):
        """All-resistant lineage: the odds ratio is undefined, not huge."""
        rows = []
        rows += [("ST1", 1, p) for p in ["R"] * 8 + ["S"] * 2]
        rows += [("ST1", 0, p) for p in ["R"] * 2 + ["S"] * 8]
        rows += [("ST2", 1, p) for p in ["R"] * 10]
        rows += [("ST2", 0, p) for p in ["R"] * 10]
        gi = _input_from_rows(config, rows, "gene__femA_novel")
        estimates = lm.stratum_estimates(gi, "gene__femA_novel")
        assert [e.lineage for e in estimates] == ["ST1"]

    def test_strata_carry_the_no_kinship_stamp(self, config, two_lineages):
        gi, feature = two_lineages
        for estimate in lm.stratum_estimates(gi, feature):
            assert estimate.model.endswith(":NO_KINSHIP_CORRECTION")


class TestPooling:
    @pytest.fixture(scope="class")
    def pooled(self, config, two_lineages):
        gi, feature = two_lineages
        (result,) = lm.lineage_meta_analysis(gi, config)
        return result

    def test_the_pool_is_the_inverse_variance_mean(self, two_lineages, pooled):
        gi, _feature = two_lineages
        e1, e2 = lm.stratum_estimates(gi, "gene__femA_novel")
        w1, w2 = 1 / e1.std_error**2, 1 / e2.std_error**2
        expected = (w1 * e1.effect + w2 * e2.effect) / (w1 + w2)
        expected_se = math.sqrt(1 / (w1 + w2))
        assert pooled.status == "tested"
        assert pooled.pooled_effect == pytest.approx(expected)
        assert pooled.pooled_std_error == pytest.approx(expected_se)
        assert pooled.p_value == pytest.approx(_two_sided_p(expected / expected_se))

    def test_cochran_q_and_i_squared_match_the_formula(self, two_lineages, pooled):
        gi, _feature = two_lineages
        e1, e2 = lm.stratum_estimates(gi, "gene__femA_novel")
        w1, w2 = 1 / e1.std_error**2, 1 / e2.std_error**2
        mean = (w1 * e1.effect + w2 * e2.effect) / (w1 + w2)
        q = w1 * (e1.effect - mean) ** 2 + w2 * (e2.effect - mean) ** 2
        expected_i2 = max(0.0, (q - 1.0) / q) * 100.0
        assert pooled.q_statistic == pytest.approx(q)
        assert pooled.q_p_value == pytest.approx(_chi2_sf_1(q))
        assert pooled.i_squared == pytest.approx(expected_i2)
        assert pooled.n_lineages == 2

    def test_opposite_lineage_effects_are_reported_as_heterogeneous(
        self, two_lineages, pooled
    ):
        """ST1 says the feature raises resistance, ST2 says it lowers it."""
        assert pooled.q_statistic > 1.0
        assert pooled.q_p_value < 0.05
        assert pooled.i_squared > 50.0

    def test_meta_models_carry_the_no_kinship_stamp(self, pooled):
        assert pooled.model.endswith(":NO_KINSHIP_CORRECTION")

    def test_a_feature_absent_from_the_matrix_is_refused_by_name(self, two_lineages):
        """A typo in an explicit feature list must not come back as an empty row."""
        gi, _feature = two_lineages
        with pytest.raises(DataContractError) as exc:
            lm.lineage_meta_analysis(
                gi, load_config(SCIENCE, machine="laptop"), features=["gene__nope"]
            )
        assert "gene__nope" in str(exc.value)


class TestInsufficiencyAndRefusal:
    def test_a_feature_carried_in_one_lineage_only_is_not_pooled(self, config):
        rows = []
        rows += [("ST1", 1, p) for p in ["R"] * 8 + ["S"] * 2]
        rows += [("ST1", 0, p) for p in ["R"] * 2 + ["S"] * 8]
        # ST2 carries nothing: a stratum with no carriers has no odds ratio.
        rows += [("ST2", 0, p) for p in ["R"] * 10 + ["S"] * 10]
        gi = _input_from_rows(config, rows, "gene__femA_novel")
        (result,) = lm.lineage_meta_analysis(gi, config)
        assert result.status == "insufficient_lineages"
        assert result.pooled_effect is None
        assert result.p_value is None
        assert result.q_statistic is None
        assert result.i_squared is None
        assert result.n_lineages == 1

    def test_a_cohort_without_lineage_labels_is_refused(self, config):
        phenotypes = {f"S{i:04d}": ("R" if i < 20 else "S") for i in range(40)}
        columns = {"gene__femA_novel": {s: int(i % 2) for i, s in enumerate(sorted(phenotypes))}}
        gi = _make_input(config, phenotypes, columns)
        with pytest.raises(DataContractError) as exc:
            lm.lineage_meta_analysis(gi, config)
        assert "lineage.method" in str(exc.value)

    def test_a_single_lineage_cohort_is_refused(self, config):
        phenotypes = {f"S{i:04d}": ("R" if i < 20 else "S") for i in range(40)}
        sample_ids = sorted(phenotypes)
        columns = {"gene__femA_novel": {s: int(i % 2) for i, s in enumerate(sample_ids)}}
        lineages = {s: "ST1" for s in sample_ids}
        gi = _make_input(config, phenotypes, columns, lineages=lineages)
        with pytest.raises(DataContractError) as exc:
            lm.lineage_meta_analysis(gi, config)
        assert "lineage.method" in str(exc.value)
