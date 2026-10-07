"""The ``layers`` configuration block: drop thresholds are configuration.

Ticket 15 plus the layer-encoding package drop a feature present in more than
98 percent of isolates or in fewer than 5 isolates, and log what was dropped.
Those two numbers are scientific constants, so they live in
``config/science.yaml`` under ``layers`` - the same pattern as
``cohort_gate`` - and never in code:

* the shipped configuration states both keys, so a reader can see and edit
  them without opening a module;
* a stripped-down configuration falls back to documented code defaults;
* an invalid value is a :class:`ConfigError` whose message names the key a
  human must edit, so the remedy travels with the failure.

Nothing here reads a path, a thread count or a threshold from anywhere else.
Synthetic only: every configuration under test is a copy of the repository's
own ``config/`` written into ``tmp_path``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import ConfigError
from papipeline.layers.settings import (
    DEFAULT_MAX_PREVALENCE,
    DEFAULT_MIN_CARRIERS,
    LayerSettings,
    layer_settings,
)

PIPELINE_ROOT = Path(__file__).resolve().parents[2]

#: The keys, spelled as literals: a refusal must name a string a human can
#: paste into ``config/science.yaml``, so the test pins the spelling itself
#: rather than importing the constant the implementation also defines.
MIN_CARRIERS_KEY = "layers.min_carriers"
MAX_PREVALENCE_KEY = "layers.max_prevalence"

#: The shipped values this package adds to ``config/science.yaml``.
SHIPPED_MIN_CARRIERS = "min_carriers: 5"
SHIPPED_MAX_PREVALENCE = "max_prevalence: 0.98"


def _config_under(root: Path, mutate=None) -> PipelineConfig:
    """A self-contained copy of the repository configuration under ``root``.

    The ``layers`` keys are edited through the same YAML the run reads, so a
    change to the shipped spelling fails this helper loudly rather than
    silently testing against a stale literal.
    """
    root = Path(root)
    shutil.copytree(PIPELINE_ROOT / "config", root / "config")
    science = root / "config" / "science.yaml"
    text = science.read_text(encoding="utf-8")
    assert SHIPPED_MIN_CARRIERS in text, (
        "config/science.yaml no longer carries the `layers` "
        f"`{SHIPPED_MIN_CARRIERS}` line this test rewrites"
    )
    assert SHIPPED_MAX_PREVALENCE in text, (
        "config/science.yaml no longer carries the `layers` "
        f"`{SHIPPED_MAX_PREVALENCE}` line this test rewrites"
    )
    if mutate is not None:
        text = mutate(text)
    science.write_text(text, encoding="utf-8")
    return load_config(science, machine=None)


class TestTheShippedConfiguration:
    def test_both_thresholds_are_stated_in_science_yaml(self):
        text = (PIPELINE_ROOT / "config" / "science.yaml").read_text(
            encoding="utf-8"
        )
        assert SHIPPED_MIN_CARRIERS in text
        assert SHIPPED_MAX_PREVALENCE in text
        assert "layers:" in text

    def test_the_shipped_values_are_read_back(self, tmp_path: Path):
        config = _config_under(tmp_path)
        assert layer_settings(config) == LayerSettings(
            min_carriers=5, max_prevalence=0.98
        )

    def test_a_changed_threshold_is_read_from_the_file_not_from_code(
        self, tmp_path: Path
    ):
        config = _config_under(
            tmp_path, lambda text: text.replace("min_carriers: 5", "min_carriers: 7")
        )
        assert layer_settings(config).min_carriers == 7


class TestTheDefaults:
    def test_the_code_defaults_match_the_shipped_values(self):
        assert DEFAULT_MIN_CARRIERS == 5
        assert DEFAULT_MAX_PREVALENCE == 0.98

    def test_a_section_without_the_keys_uses_the_documented_defaults(
        self, tmp_path: Path
    ):
        def _strip(text: str) -> str:
            assert f"  {SHIPPED_MIN_CARRIERS}\n" in text
            assert f"  {SHIPPED_MAX_PREVALENCE}\n" in text
            return text.replace(f"  {SHIPPED_MIN_CARRIERS}\n", "").replace(
                f"  {SHIPPED_MAX_PREVALENCE}\n", ""
            )

        config = _config_under(tmp_path, _strip)
        assert layer_settings(config) == LayerSettings(
            min_carriers=DEFAULT_MIN_CARRIERS,
            max_prevalence=DEFAULT_MAX_PREVALENCE,
        )


class TestInvalidValuesAreRefused:
    def test_a_non_numeric_carrier_count_names_the_key(self, tmp_path: Path):
        config = _config_under(
            tmp_path, lambda text: text.replace("min_carriers: 5", "min_carriers: five")
        )
        with pytest.raises(ConfigError) as excinfo:
            layer_settings(config)
        assert MIN_CARRIERS_KEY in str(excinfo.value)

    def test_a_negative_carrier_count_names_the_key(self, tmp_path: Path):
        config = _config_under(
            tmp_path, lambda text: text.replace("min_carriers: 5", "min_carriers: -1")
        )
        with pytest.raises(ConfigError) as excinfo:
            layer_settings(config)
        assert MIN_CARRIERS_KEY in str(excinfo.value)

    def test_a_prevalence_outside_the_unit_interval_names_the_key(
        self, tmp_path: Path
    ):
        config = _config_under(
            tmp_path,
            lambda text: text.replace("max_prevalence: 0.98", "max_prevalence: 1.5"),
        )
        with pytest.raises(ConfigError) as excinfo:
            layer_settings(config)
        assert MAX_PREVALENCE_KEY in str(excinfo.value)

    def test_a_non_numeric_prevalence_names_the_key(self, tmp_path: Path):
        config = _config_under(
            tmp_path,
            lambda text: text.replace("max_prevalence: 0.98", "max_prevalence: high"),
        )
        with pytest.raises(ConfigError) as excinfo:
            layer_settings(config)
        assert MAX_PREVALENCE_KEY in str(excinfo.value)

    def test_a_section_that_is_not_a_mapping_is_refused_naming_the_key(
        self, tmp_path: Path
    ):
        def _scalar(text: str) -> str:
            assert "layers:\n" in text
            # A scalar cannot have an indented body, so the two keys come out
            # too - the helper asserted their spelling before this ran.
            text = text.replace(f"  {SHIPPED_MIN_CARRIERS}\n", "")
            text = text.replace(f"  {SHIPPED_MAX_PREVALENCE}\n", "")
            return text.replace("layers:\n", "layers: enabled\n", 1)

        config = _config_under(tmp_path, _scalar)
        with pytest.raises(ConfigError) as excinfo:
            layer_settings(config)
        assert "layers" in str(excinfo.value)


class TestTheBoundaryRuleItself:
    def test_the_carrier_floor_is_strict(self):
        settings = LayerSettings(min_carriers=5, max_prevalence=0.98)
        assert settings.drop_reason(4, 20) is not None
        assert settings.drop_reason(5, 20) is None

    def test_the_prevalence_ceiling_is_strict(self):
        settings = LayerSettings(min_carriers=5, max_prevalence=0.98)
        # 19/20 = 0.95, at or below the ceiling: kept.
        assert settings.drop_reason(19, 20) is None
        # 20/20 = 1.0, more than 0.98: dropped as background.
        assert settings.drop_reason(20, 20) is not None

    def test_a_drop_reason_names_the_config_key(self):
        settings = LayerSettings(min_carriers=5, max_prevalence=0.98)
        assert MIN_CARRIERS_KEY in settings.drop_reason(1, 20)
        assert MAX_PREVALENCE_KEY in settings.drop_reason(20, 20)

    def test_a_drop_reason_reports_the_counts(self):
        settings = LayerSettings(min_carriers=5, max_prevalence=0.98)
        assert "1" in settings.drop_reason(1, 20)
        assert "20" in settings.drop_reason(20, 20)
