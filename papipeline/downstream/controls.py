"""Deliverable 1: the positive-control gate.

Before anything downstream reports a finding, the baseline association table
has to show that the scan can still recover what this project already knows:
``oprD_burden`` (the OprD-loss features) and ``any_MBL``, by default, listed in
``config/science.yaml`` under ``downstream.positive_controls``.

A control is *recovered* when a feature that matches it is significant at
``gwas.significance_threshold``. Everything else is a failure, and the failure
names EVERY missing control and says which of two things went wrong: the
control was **absent from the association table** (no feature matched it at
all) or it was **present but not significant**. An empty table recovers
nothing, so the gate refuses an empty table rather than passing vacuously.

What counts as a match, and why it is not simply "any feature named after the
gene": a control that NAMES a determinant (``any_MBL``, a gene) is recovered by
detecting that determinant, including its bare ``gene__blaNDM-1`` column. A
control named after a STATE of a gene - ``oprD_burden``, which resolves to
``oprD`` through the *stem* rule - is recovered only by features of that state
(``gene__oprD_absent``, ``L3_oprD_off_any``), never by the bare presence
column. ``config/mechanisms.tsv`` is explicit that detecting ``oprD``
sequence is not detecting OprD loss, so the bare presence column must not be
able to pass this gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

from ..errors import ConfigError, PipelineError
from .known_determinants import (
    KIND_COMPOSITE,
    KIND_GENE,
    KnownDeterminants,
    load_known_determinants,
    normalise_label,
)

#: The key the control list is read from, and the key a refusal names.
CONFIG_KEY = "downstream.positive_controls"

#: Controls used when the science file does not list any. Stated here rather
#: than invented per call so the default is one place to read.
DEFAULT_POSITIVE_CONTROLS: Tuple[str, ...] = ("oprD_burden", "any_MBL")

#: Resolution paths that mean "this control IS this determinant", so detecting
#: the determinant bare (``gene__blaNDM-1``) counts as recovery.
_DIRECT_KINDS = frozenset({KIND_GENE, KIND_COMPOSITE, "family"})


class ControlGateError(PipelineError):
    """The baseline scan did not recover every configured control."""


@dataclass(frozen=True)
class ControlSpec:
    """One configured control: its name, its constituents, and how it matches.

    ``bare_member_match`` is False for a stem-resolved control (``oprD_burden``
    -> ``oprD``), which is what keeps ``gene__oprD`` - an intact porin - from
    satisfying an OprD-*loss* control.
    """

    name: str
    members: Tuple[str, ...] = ()
    bare_member_match: bool = False

    def matches(self, feature: str) -> bool:
        label = normalise_label(feature).casefold()
        if label == self.name.casefold():
            return True
        for member in self.members:
            key = member.casefold()
            if label == key and self.bare_member_match:
                return True
            if label.startswith(key + "_"):
                return True
        return False


@dataclass(frozen=True)
class ControlOutcome:
    """The verdict for one control."""

    control: str
    recovered: bool
    reason: str
    matched_features: Tuple[str, ...] = ()
    best_adjusted_p: Optional[float] = None


def positive_controls(config) -> List[str]:
    """The configured control names, in the order the science file lists them.

    Raises:
        ConfigError: the key is present but empty, or is not a list. An empty
            list is a refusal rather than an instruction to check nothing.
    """
    section = config.raw.get("downstream") or {}
    configured = section.get("positive_controls", list(DEFAULT_POSITIVE_CONTROLS))
    if isinstance(configured, str) or not isinstance(configured, (list, tuple)):
        raise ConfigError(
            f"{CONFIG_KEY} must be a list of control names",
            key=CONFIG_KEY,
            configured=repr(configured),
        )
    names = [str(name).strip() for name in configured]
    if not names or any(not name for name in names):
        raise ConfigError(
            f"{CONFIG_KEY} is empty; the positive-control gate needs at least "
            "one control, because a gate that checks nothing passes vacuously",
            key=CONFIG_KEY,
        )
    return names


def resolve_controls(
    config, known: Optional[KnownDeterminants] = None
) -> List[ControlSpec]:
    """Turn the configured names into specs, resolving their constituents."""
    knowledge = known if known is not None else load_known_determinants(config)
    specs: List[ControlSpec] = []
    for name in positive_controls(config):
        kind = knowledge.match_kind(name)
        members = tuple(knowledge.members_for(name) or ())
        specs.append(
            ControlSpec(
                name=name,
                members=members,
                bare_member_match=kind in _DIRECT_KINDS,
            )
        )
    return specs


def evaluate_controls(
    results: Iterable,
    specs: Sequence[ControlSpec],
    threshold: float,
) -> List[ControlOutcome]:
    """Judge each control against an association table.

    ``results`` is anything with ``feature`` and ``adjusted_p_value`` - the
    baseline stage-12 rows.
    """
    rows = list(results)
    outcomes: List[ControlOutcome] = []
    for spec in specs:
        matched = [row for row in rows if spec.matches(row.feature)]
        # A matched row the engine did not test reports no adjusted p; None is
        # not comparable with a float, so the "best" is taken over the rows
        # that actually have one.
        values = [
            row.adjusted_p_value
            for row in matched
            if row.adjusted_p_value is not None
        ]
        significant = [
            row for row in matched if row.adjusted_p_value is not None and row.adjusted_p_value <= threshold
        ]
        if significant:
            reason = "recovered"
        elif matched:
            reason = "not significant"
        else:
            reason = "absent"
        outcomes.append(
            ControlOutcome(
                control=spec.name,
                recovered=bool(significant),
                reason=reason,
                matched_features=tuple(row.feature for row in matched),
                best_adjusted_p=min(values) if values else None,
            )
        )
    return outcomes


def _explanation(outcome: ControlOutcome, threshold: float) -> str:
    if outcome.reason == "absent":
        return "absent from the association table (no feature matched it)"
    best = outcome.best_adjusted_p
    reported = "no adjusted p" if best is None else f"best adjusted p={best:.3g}"
    return f"present but not significant ({reported} > {threshold:g})"


def run_control_gate(config, results: Iterable) -> List[ControlOutcome]:
    """Recover every configured control or raise, naming each one that failed.

    Raises:
        ControlGateError: at least one control was not recovered. The message
            lists every missing control and why; ``context["missing"]`` holds
            their names and ``context["key"]`` the configuration key.
    """
    specs = resolve_controls(config)
    threshold = config.gwas.significance_threshold
    outcomes = evaluate_controls(results, specs, threshold)
    missing = [outcome for outcome in outcomes if not outcome.recovered]
    if not missing:
        return outcomes

    lines = [
        f"  {outcome.control}: {_explanation(outcome, threshold)}"
        for outcome in missing
    ]
    message = (
        "positive-control gate failed: the baseline association table did not "
        f"recover {len(missing)} of {len(outcomes)} configured control(s); "
        "the analyses that depend on it must not run\n" + "\n".join(lines)
    )
    raise ControlGateError(
        message,
        missing=tuple(sorted(outcome.control for outcome in missing)),
        key=CONFIG_KEY,
        threshold=threshold,
    )
