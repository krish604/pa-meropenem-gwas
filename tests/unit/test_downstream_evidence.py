"""Deliverable 5: the novelty filter and the evidence tiers.

``config/known_determinants.tsv`` is what "already known" means, so a feature
it resolves to is not a novelty candidate. Everything else is graded A to D:

* **A** - survives conditioning AND is seen in enough independent lineages;
* **B** - survives conditioning but is lineage-bound (or below the lineage
  floor);
* **C** - significant in the baseline scan only, i.e. the association did not
  survive conditioning on the known determinants;
* **D** - not significant.

The claim attached to a tier is a repository claim level and nothing else:
``SUPPORTED``, ``ASSOCIATED`` or ``UNKNOWN``. There is no CAUSAL level to
reach for - ``papipeline.models.ClaimStatus`` has no such member, and this
module never invents one.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import evidence as ev
from papipeline.downstream import known_determinants as kd
from papipeline.models import ClaimStatus, ConvergenceCall, ConvergenceCategory, GwasResult

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture(scope="module")
def known(config) -> kd.KnownDeterminants:
    return kd.load_known_determinants(config)


def _result(feature: str, adjusted_p, lineages=None) -> GwasResult:
    return GwasResult(
        feature=feature,
        feature_type="gene",
        effect=2.0,
        p_value=1e-4 if adjusted_p is None else min(1.0, adjusted_p / 10.0),
        adjusted_p_value=adjusted_p,
        effect_size=2.0,
        frequency=0.2,
        # `is None` rather than `or`: an explicit {} must stay empty, because
        # one test below is exactly about a result with no distribution.
        lineage_distribution=dict({"ST1": 4, "ST2": 4} if lineages is None else lineages),
        model="test",
    )


class TestNoveltyFilter:
    def test_known_and_novel_are_partitioned(self, known):
        features = [
            "gene__oprD_absent",
            "gene_presence_absence__blaNDM-1",
            "L2_PDC_high",
            "gene__femA_novel",
            "gene__strong",
        ]
        known_hits, novel = ev.novelty_partition(features, known)
        assert known_hits == [
            "gene__oprD_absent",
            "gene_presence_absence__blaNDM-1",
            "L2_PDC_high",
        ]
        assert novel == ["gene__femA_novel", "gene__strong"]

    def test_a_novel_feature_records_no_known_match(self, known):
        record = ev.classify_feature(
            "gene__femA_novel",
            conditional_p=1e-4,
            baseline_p=1e-4,
            independent_lineages=3,
            known=known,
            config=load_config(SCIENCE, machine="laptop"),
        )
        assert record.novelty == "novel"
        assert record.known_matches == ()

    def test_a_known_feature_names_the_determinant_it_matched(self, known):
        record = ev.classify_feature(
            "gene_presence_absence__blaNDM-1",
            conditional_p=1e-4,
            baseline_p=1e-4,
            independent_lineages=3,
            known=known,
            config=load_config(SCIENCE, machine="laptop"),
        )
        assert record.novelty == "known"
        assert record.known_matches == ("blaNDM-1",)


class TestEvidenceTiers:
    def _classify(self, known, *, conditional_p, baseline_p, lineages, config=None):
        return ev.classify_feature(
            "gene__femA_novel",
            conditional_p=conditional_p,
            baseline_p=baseline_p,
            independent_lineages=lineages,
            known=known,
            config=config or load_config(SCIENCE, machine="laptop"),
        )

    def test_tier_a_needs_conditioning_and_lineages(self, known):
        record = self._classify(
            known, conditional_p=1e-4, baseline_p=1e-5, lineages=3
        )
        assert record.evidence_tier == "A"
        assert record.claim_status == ClaimStatus.SUPPORTED.value
        assert record.basis == "conditional_and_convergent"

    def test_tier_b_when_the_signal_is_single_lineage(self, known):
        record = self._classify(
            known, conditional_p=1e-4, baseline_p=1e-5, lineages=1
        )
        assert record.evidence_tier == "B"
        assert record.claim_status == ClaimStatus.ASSOCIATED.value

    def test_tier_c_when_only_the_baseline_scan_was_significant(self, known):
        record = self._classify(
            known, conditional_p=0.4, baseline_p=1e-4, lineages=3
        )
        assert record.evidence_tier == "C"
        assert record.claim_status == ClaimStatus.ASSOCIATED.value
        assert record.basis == "marginal_only"

    def test_tier_d_when_neither_scan_is_significant(self, known):
        record = self._classify(
            known, conditional_p=0.4, baseline_p=0.6, lineages=3
        )
        assert record.evidence_tier == "D"
        assert record.claim_status == ClaimStatus.UNKNOWN.value
        assert record.basis == "not_significant"

    def test_a_missing_adjusted_p_is_not_significance(self, known):
        record = self._classify(
            known, conditional_p=None, baseline_p=None, lineages=3
        )
        assert record.evidence_tier == "D"
        assert record.claim_status == ClaimStatus.UNKNOWN.value

    def test_the_threshold_comes_from_the_config(self, known, config):
        tuned = replace(
            config, raw={**config.raw, "gwas": {**config.raw["gwas"], "significance_threshold": 0.01}}
        )
        # The reader must see the OVERRIDDEN value. Note what it is read from:
        # `dataclasses.replace(config, raw=...)` does not rebuild the derived
        # `gwas` field, so `tuned.gwas.significance_threshold` is still the
        # committed 0.05 while `tuned.raw[...]` is 0.01 - which is why the
        # threshold is read from the raw mapping.
        assert ev.significance_threshold(tuned) == 0.01
        borderline = self._classify(known, conditional_p=0.02, baseline_p=0.02, lineages=3, config=tuned)
        assert borderline.evidence_tier == "D"

    def test_the_lineage_floor_comes_from_the_config(self, known, config):
        tuned = replace(
            config,
            raw={
                **config.raw,
                "convergence": {**config.raw["convergence"], "min_independent_lineages": 3},
            },
        )
        assert ev.lineage_floor(tuned) == 3
        record = self._classify(
            known, conditional_p=1e-4, baseline_p=1e-4, lineages=2, config=tuned
        )
        assert record.evidence_tier == "B"


class TestLineageCounting:
    def test_lineages_are_counted_from_the_carrier_distribution(self, config, known):
        record = ev.classify_feature(
            "gene__femA_novel",
            conditional_p=1e-4,
            baseline_p=1e-4,
            independent_lineages=ev.independent_lineages(
                _result("gene__femA_novel", 1e-4, {"ST1": 3, "ST2": 2, "ST3": 1})
            ),
            known=known,
            config=config,
        )
        assert record.independent_lineages == 3

    def test_the_unknown_sentinel_is_not_a_lineage(self, config, known):
        count = ev.independent_lineages(
            _result("gene__femA_novel", 1e-4, {"unknown": 8, "ST1": 2})
        )
        assert count == 1

    def test_a_convergence_call_overrides_the_derived_count(self, config, known):
        call = ConvergenceCall(
            determinant="gene__femA_novel",
            independent_lineages=1,
            branch_count=4,
            distribution={"ST1": 4},
            convergence_category=ConvergenceCategory.LINEAGE_ASSOCIATED,
        )
        count = ev.independent_lineages(
            _result("gene__femA_novel", 1e-4, {"ST1": 3, "ST2": 2}),
            convergence={"gene__femA_novel": call},
        )
        assert count == 1

    def test_a_result_without_a_distribution_counts_zero(self):
        assert ev.independent_lineages(_result("gene__strong", 1e-4, {})) == 0


class TestBuildEvidenceTable:
    def test_the_table_covers_both_scans(self, config, known):
        baseline = [
            _result("gene__femA_novel", 1e-4),
            _result("gene__oprD_absent", 1e-4),
            _result("gene__noise", 0.8),
        ]
        conditional = [
            _result("gene__femA_novel", 1e-4),
            _result("gene__oprD_absent", 1e-4),
            _result("gene__noise", 0.8),
        ]
        table = ev.build_evidence_table(
            conditional_results=conditional,
            baseline_results=baseline,
            known=known,
            config=config,
        )
        by_feature = {r.feature: r for r in table}
        assert set(by_feature) == {
            "gene__femA_novel",
            "gene__oprD_absent",
            "gene__noise",
        }
        # Two lineages in the distribution, floor is 2 -> tier A.
        assert by_feature["gene__femA_novel"].evidence_tier == "A"
        assert by_feature["gene__femA_novel"].novelty == "novel"
        assert by_feature["gene__oprD_absent"].evidence_tier == "A"
        assert by_feature["gene__oprD_absent"].novelty == "known"
        assert by_feature["gene__noise"].evidence_tier == "D"
        assert by_feature["gene__noise"].claim_status == ClaimStatus.UNKNOWN.value

    def test_a_feature_only_the_baseline_scan_saw_is_tier_c(self, config, known):
        table = ev.build_evidence_table(
            conditional_results=[],
            baseline_results=[_result("gene__femA_novel", 1e-4)],
            known=known,
            config=config,
        )
        (record,) = table
        assert record.evidence_tier == "C"
        assert record.conditional_adjusted_p is None
        assert record.baseline_adjusted_p == pytest.approx(1e-4)


class TestClaimLevels:
    def test_every_emitted_claim_is_a_repository_level(self, config, known):
        records = [
            ev.classify_feature(
                feature,
                conditional_p=conditional_p,
                baseline_p=baseline_p,
                independent_lineages=lineages,
                known=known,
                config=config,
            )
            for feature in ("gene__femA_novel", "gene__oprD_absent")
            for conditional_p, baseline_p, lineages in (
                (1e-4, 1e-4, 3),
                (1e-4, 1e-4, 1),
                (0.4, 1e-4, 3),
                (0.4, 0.6, 3),
            )
        ]
        allowed = {status.value for status in ClaimStatus}
        assert {r.claim_status for r in records} <= allowed
        assert "CAUSAL" not in {r.claim_status for r in records}
        assert "CAUSAL" not in allowed

    def test_only_the_four_tiers_exist(self):
        assert ev.EVIDENCE_TIERS == ("A", "B", "C", "D")
        assert set(ev.TIER_CLAIMS) == set(ev.EVIDENCE_TIERS)
