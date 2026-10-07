"""Deliverable 6: the downstream report.

The report is the only place the downstream package puts a claim in front of a
reader, so this file pins three things about it:

1. the positive-control gate runs FIRST - a report that cannot show it
   recovered the controls is not produced at all, and the refusal names the
   control it failed to recover;
2. every claim on a row is one of the repository's levels - ``DETECTED``,
   ``PREDICTED``, ``ASSOCIATED``, ``SUPPORTED``, ``UNKNOWN`` - and no other
   level, in particular no ``CAUSAL``, exists anywhere in the output;
3. what is written to disk is what can be read back.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import controls as cg
from papipeline.downstream import known_determinants as kd
from papipeline.downstream import report as rp
from papipeline.io.tsv import read_tsv
from papipeline.models import ClaimStatus, GwasResult

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


def _result(feature: str, adjusted_p: float, lineages=None) -> GwasResult:
    return GwasResult(
        feature=feature,
        feature_type="gene",
        effect=3.0,
        p_value=min(1.0, adjusted_p / 10.0),
        adjusted_p_value=adjusted_p,
        effect_size=3.0,
        frequency=0.2,
        lineage_distribution=dict(lineages or {"ST1": 4, "ST2": 4}),
        model="test",
    )


def _recovering_baseline() -> list:
    return [
        _result("gene_presence_absence__blaNDM-1", 1e-6),
        _result("gene__oprD_absent", 1e-6),
    ]


class TestTheGateRunsFirst:
    def test_a_control_the_baseline_missed_stops_the_report(self, config):
        """One control recovered is not enough: the report never starts."""
        with pytest.raises(cg.ControlGateError) as exc:
            rp.build_report(config, baseline_results=[_result("gene__oprD_absent", 1e-6)])
        assert "any_MBL" in str(exc.value)
        assert set(exc.value.context["missing"]) == {"any_MBL"}

    def test_an_empty_baseline_stops_the_report(self, config):
        with pytest.raises(cg.ControlGateError) as exc:
            rp.build_report(config, baseline_results=[])
        assert set(exc.value.context["missing"]) == {"oprD_burden", "any_MBL"}

    def test_a_recovering_baseline_produces_a_report(self, config):
        report = rp.build_report(config, baseline_results=_recovering_baseline())
        assert {o.control for o in report.controls} == {"oprD_burden", "any_MBL"}
        assert all(o.recovered for o in report.controls)


class TestClaimLevels:
    @pytest.fixture(scope="class")
    def report(self, config):
        return rp.build_report(
            config,
            baseline_results=_recovering_baseline(),
            conditional_results=[
                _result("gene__femA_novel", 1e-4),
                _result("gene__noise", 0.8),
            ],
        )

    def test_every_row_carries_a_repository_claim_level(self, report):
        allowed = {status.value for status in ClaimStatus}
        rows = report.rows()
        assert rows, "a report with no rows has nothing to review"
        assert {row["claim_status"] for row in rows} <= allowed
        assert report.claim_levels() <= allowed

    def test_no_causal_level_reaches_the_output(self, report):
        assert "CAUSAL" not in report.claim_levels()
        rendered = " ".join(str(row) for row in report.rows())
        assert "CAUSAL" not in rendered
        assert "causal" not in rendered.lower()

    def test_a_novel_tier_a_row_is_the_strongest_things_get(self, report):
        by_feature = {row["feature"]: row for row in report.rows()}
        assert by_feature["gene__femA_novel"]["evidence_tier"] == "A"
        assert by_feature["gene__femA_novel"]["claim_status"] == "SUPPORTED"
        assert by_feature["gene__femA_novel"]["novelty"] == "novel"
        assert by_feature["gene__noise"]["claim_status"] == "UNKNOWN"


class TestSectionsTravelIntoTheReport:
    def test_tier_one_tier_two_and_meta_are_carried_through(self, config):
        from papipeline.downstream import interactions as ix
        from papipeline.downstream import lineage_meta as lm

        tier1 = ix.InteractionResult(
            tier=1,
            antibiotic="imipenem",
            feature_a="PDC_high",
            feature_b="oprD_off_any",
            status="tested",
            p_value=1e-3,
            adjusted_p_value=2e-3,
            interaction_p_value=1e-3,
            interaction_estimate=1.5,
            interaction_std_error=0.5,
            n_samples=200,
        )
        meta = lm.MetaResult(
            feature="gene__femA_novel",
            status="tested",
            n_lineages=2,
            pooled_effect=0.6,
            pooled_std_error=0.3,
            p_value=0.045,
            q_statistic=0.4,
            q_p_value=0.53,
            i_squared=0.0,
        )
        report = rp.build_report(
            config,
            baseline_results=_recovering_baseline(),
            tier1=[tier1],
            tier2=[tier1],
            meta=[meta],
        )
        assert [r.feature_a for r in report.tier1] == ["PDC_high"]
        assert [r.feature_a for r in report.tier2] == ["PDC_high"]
        assert [m.feature for m in report.meta] == ["gene__femA_novel"]


class TestRoundTrip:
    def test_what_is_written_reads_back(self, config, tmp_path):
        report = rp.build_report(
            config,
            baseline_results=_recovering_baseline(),
            conditional_results=[_result("gene__femA_novel", 1e-4)],
        )
        path = report.write(tmp_path / "downstream_evidence.tsv")
        rows = read_tsv(path, required_columns=rp.REPORT_COLUMNS, unique_columns=["feature"])
        assert len(rows) == len(report.rows())
        by_feature = {row["feature"]: row for row in rows}
        record = by_feature["gene__femA_novel"]
        assert record["novelty"] == "novel"
        assert record["claim_status"] == "SUPPORTED"
        assert record["evidence_tier"] == "A"
        assert float(record["conditional_adjusted_p"]) == pytest.approx(1e-4)
