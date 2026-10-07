"""The positive-control gate (deliverable 1).

The baseline scan must recover the configured controls - by default
``oprD_burden`` and ``any_MBL`` - or the downstream analyses that follow it
must not run. A failure names EVERY missing control, says whether it was
absent from the association table or present but not significant, and names
the configuration key that lists the controls.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import controls as cg
from papipeline.models import GwasResult

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


def _result(feature: str, adjusted_p: float, effect: float = 3.0) -> GwasResult:
    return GwasResult(
        feature=feature,
        feature_type="gene",
        effect=effect,
        p_value=min(1.0, adjusted_p / 10.0),
        adjusted_p_value=adjusted_p,
        effect_size=effect,
        frequency=0.2,
        lineage_distribution={"ST1": 4, "ST2": 4},
        model="test",
    )


class TestGatePasses:
    def test_both_constituent_level_controls_recover(self, config):
        results = [
            _result("gene_presence_absence__blaNDM-1", 1e-6),
            _result("gene__oprD_absent", 1e-6),
        ]
        outcomes = cg.run_control_gate(config, results)
        assert {o.control for o in outcomes} == {"oprD_burden", "any_MBL"}
        assert all(o.recovered for o in outcomes)
        assert all(o.reason == "recovered" for o in outcomes)

    def test_a_layer_suffixed_control_matches_by_suffix(self, config):
        results = [
            _result("L1_any_MBL", 1e-4),
            _result("L3_oprd_burden", 1e-4),
        ]
        outcomes = cg.run_control_gate(config, results)
        assert all(o.recovered for o in outcomes)


class TestGateFailsLoudly:
    def test_a_missing_control_is_named_in_the_message_and_the_context(self, config):
        results = [_result("gene__oprD_absent", 1e-6)]
        with pytest.raises(cg.ControlGateError) as exc:
            cg.run_control_gate(config, results)
        message = str(exc.value)
        assert "any_MBL" in message
        missing = set(exc.value.context["missing"])
        assert missing == {"any_MBL"}, missing
        assert exc.value.context["key"] == "downstream.positive_controls"

    def test_present_but_not_significant_is_not_recovery(self, config):
        results = [
            _result("gene_presence_absence__blaNDM-1", 0.5),
            _result("gene__oprD_absent", 0.5),
        ]
        with pytest.raises(cg.ControlGateError) as exc:
            cg.run_control_gate(config, results)
        message = str(exc.value)
        assert "not significant" in message
        assert set(exc.value.context["missing"]) == {"any_MBL", "oprD_burden"}

    def test_the_gate_runs_before_anything_else_can_be_reported(self, config):
        """An empty association table recovers nothing and must fail."""
        with pytest.raises(cg.ControlGateError) as exc:
            cg.run_control_gate(config, [])
        assert set(exc.value.context["missing"]) == {"any_MBL", "oprD_burden"}


class TestControlConfiguration:
    def test_the_default_control_list_is_the_two_documented_controls(self, config):
        assert cg.positive_controls(config) == ["oprD_burden", "any_MBL"]

    def test_the_control_list_comes_from_the_config_key(self, config):
        tuned = replace(config, raw={**config.raw, "downstream": {"positive_controls": ["any_MBL"]}})
        assert cg.positive_controls(tuned) == ["any_MBL"]
        specs = cg.resolve_controls(tuned)
        assert [s.name for s in specs] == ["any_MBL"]

    def test_an_empty_control_list_is_refused_not_obeyed(self, config):
        tuned = replace(config, raw={**config.raw, "downstream": {"positive_controls": []}})
        with pytest.raises(Exception) as exc:
            cg.positive_controls(tuned)
        assert "downstream.positive_controls" in str(exc.value)

    def test_outcomes_report_the_matched_features(self, config):
        results = [
            _result("gene_presence_absence__blaNDM-1", 1e-6),
            _result("gene__oprD_absent", 1e-6),
        ]
        outcomes = {o.control: o for o in cg.evaluate_controls(
            results, cg.resolve_controls(config), config.gwas.significance_threshold
        )}
        assert "gene_presence_absence__blaNDM-1" in outcomes["any_MBL"].matched_features
        assert "gene__oprD_absent" in outcomes["oprD_burden"].matched_features
