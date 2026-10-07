"""Drop thresholds for the layer feature matrices, read from configuration.

Ticket 15 says a feature present in more than 98 percent of isolates or in
fewer than 5 isolates is dropped, and the drop list is logged. Those two
numbers are scientific constants, so they live in ``config/science.yaml`` under
the ``layers`` section - the same pattern ``papipeline.cohort_gate`` uses for
its own thresholds - and never in this module beyond the documented code
defaults that apply when a key is stripped out.

Reading goes through ``config.raw`` rather than a typed field because
``science.yaml`` is free to carry top-level sections the loader does not model;
a machine overlay that redeclares one is rejected by the loader, but a science
file is the place new scientific keys are added. An invalid value is a
:class:`~papipeline.errors.ConfigError` whose message names the key, so the
remedy travels with the failure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional

from ..config.loader import PipelineConfig
from ..errors import ConfigError

#: The section of the science configuration this module reads.
CONFIG_SECTION = "layers"

#: The two keys, in full. A refusal message names these strings verbatim so a
#: reader can paste them into ``config/science.yaml``.
MIN_CARRIERS_KEY = "layers.min_carriers"
MAX_PREVALENCE_KEY = "layers.max_prevalence"

#: Code defaults, matching the shipped values. They exist so a stripped-down
#: configuration still has a rule rather than no rule.
DEFAULT_MIN_CARRIERS = 5
DEFAULT_MAX_PREVALENCE = 0.98


@dataclass(frozen=True)
class LayerSettings:
    """The drop rule, resolved once per run."""

    min_carriers: int = DEFAULT_MIN_CARRIERS
    max_prevalence: float = DEFAULT_MAX_PREVALENCE

    def drop_reason(self, carriers: int, n_isolates: int) -> Optional[str]:
        """Why a feature with this prevalence is dropped, or ``None`` to keep it.

        Both boundaries are strict, matching the ticket's wording: *fewer
        than* 5 isolates, *more than* 98 percent. The reason string names the
        config key rather than only the number, so a log line says which line
        of ``science.yaml`` to change.
        """
        if n_isolates <= 0:
            raise ValueError("drop_reason needs at least one isolate")
        if carriers < self.min_carriers:
            return (
                f"carriers={carriers} < {MIN_CARRIERS_KEY}={self.min_carriers}"
            )
        if carriers / n_isolates > self.max_prevalence:
            return (
                f"carriers={carriers}/{n_isolates} > "
                f"{MAX_PREVALENCE_KEY}={self.max_prevalence}"
            )
        return None


def _section(config: PipelineConfig) -> Mapping[str, Any]:
    raw = getattr(config, "raw", None) or {}
    section = raw.get(CONFIG_SECTION)
    if section is None:
        return {}
    if not isinstance(section, Mapping):
        raise ConfigError(
            f"`{CONFIG_SECTION}` must be a mapping of thresholds, e.g. "
            f"{MIN_CARRIERS_KEY} and {MAX_PREVALENCE_KEY}",
            key=CONFIG_SECTION,
        )
    return section


def _carriers(section: Mapping[str, Any]) -> int:
    raw = section.get("min_carriers")
    if raw is None:
        return DEFAULT_MIN_CARRIERS
    if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
        raise ConfigError(
            f"{MIN_CARRIERS_KEY} must be a non-negative integer: the number of "
            "isolates a feature must be present in to be kept",
            key=MIN_CARRIERS_KEY,
            value=raw,
        )
    return raw


def _prevalence(section: Mapping[str, Any]) -> float:
    raw = section.get("max_prevalence")
    if raw is None:
        return DEFAULT_MAX_PREVALENCE
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise ConfigError(
            f"{MAX_PREVALENCE_KEY} must be a number greater than 0 and at most 1: "
            "the fraction of isolates above which a feature is background",
            key=MAX_PREVALENCE_KEY,
            value=raw,
        )
    value = float(raw)
    if not 0 < value <= 1:
        raise ConfigError(
            f"{MAX_PREVALENCE_KEY} must be greater than 0 and at most 1",
            key=MAX_PREVALENCE_KEY,
            value=value,
        )
    return value


def layer_settings(config: PipelineConfig) -> LayerSettings:
    """Resolve ``layers.*`` from the science configuration."""
    section = _section(config)
    return LayerSettings(
        min_carriers=_carriers(section), max_prevalence=_prevalence(section)
    )
