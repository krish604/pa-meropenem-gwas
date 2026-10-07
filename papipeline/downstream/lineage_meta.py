"""Deliverable 4: the lineage-stratified meta-analysis.

The cohort association pools across lineages, so a feature that is a hit only
because it swept one lineage looks like a hit anywhere else. Here the effect
is estimated separately INSIDE each lineage - the log odds ratio of the
stratum's own 2x2 table - and the per-lineage estimates are combined by
fixed-effect inverse variance. Nothing is pooled within a lineage first and
nothing is dropped for being lineage-bound: the strata are the point.

Cochran's Q says how much of the between-lineage spread exceeds what sampling
error alone would produce, and I squared is that share as a percentage. Both
travel with the pooled estimate, because a pooled number without a statement
about its heterogeneity is a number a reader cannot check.

Two kinds of honesty are enforced here:

* a stratum with no carriers, no non-carriers, or one outcome class only has
  no odds ratio at all - the estimate is undefined, not enormous - so that
  stratum is left out and the count of what was used is reported
  (``n_lineages``). Fewer than two usable strata is ``insufficient_lineages``,
  not a one-lineage "meta-analysis";
* a cohort with fewer than two lineage labels is refused outright, naming
  ``lineage.method``, because "pooled across lineages" with one lineage is a
  claim about something that was not measured.

A zero cell is corrected Haldane-Anscombe (+0.5 to all four cells) and the
stratum says ``corrected is True``, so a corrected estimate can never pass for
an uncorrected one.

Everything carries the repository's ``:NO_KINSHIP_CORRECTION`` stamp: each
stratum is a subset of the same structured cohort, and no model here corrects
for kinship.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from ..errors import DataContractError
from ..stages.convergence import UNKNOWN_LINEAGE

#: Per-stratum model stamp.
STRATUM_MODEL = "stratum_log_odds_ratio:NO_KINSHIP_CORRECTION"

#: Pooled model stamp.
META_MODEL = "fixed_effect_inverse_variance:NO_KINSHIP_CORRECTION"

#: Statuses a :class:`MetaResult` can carry.
STATUS_TESTED = "tested"
STATUS_INSUFFICIENT = "insufficient_lineages"

#: The config key a refusal about missing or too-thin lineage labels names.
LINEAGE_METHOD_KEY = "lineage.method"

#: Added to every cell of a 2x2 table that contains a zero.
HALDANE_ANSCOMBE = 0.5


@dataclass(frozen=True)
class StratumEstimate:
    """One lineage's own estimate of one feature's effect."""

    lineage: str
    n_samples: int
    n_positive: int
    n_negative: int
    n_carriers: int
    effect: float
    std_error: float
    p_value: float
    corrected: bool
    model: str = STRATUM_MODEL


@dataclass(frozen=True)
class MetaResult:
    """One feature's pooled result, with the heterogeneity that qualifies it."""

    feature: str
    status: str
    n_lineages: int = 0
    estimates: Tuple[StratumEstimate, ...] = ()
    pooled_effect: Optional[float] = None
    pooled_std_error: Optional[float] = None
    p_value: Optional[float] = None
    q_statistic: Optional[float] = None
    q_p_value: Optional[float] = None
    i_squared: Optional[float] = None
    model: str = META_MODEL


def _two_sided_p(z: float) -> float:
    """Two-sided normal p from ``z``, without importing a distribution."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def _chi2_sf(value: float, df: int) -> float:
    """Survival function of chi-squared on ``df`` degrees of freedom.

    Cochran's Q on k strata has k-1 df, so this is never a fixed shape: for
    two strata it is the 1-df case (``erfc(sqrt(Q/2))``), and for three or
    more a 1-df answer would be badly wrong. scipy is a declared dependency of
    this pipeline, so the distribution comes from there rather than from a
    formula written out here.
    """
    if df < 1:
        return 1.0
    if value <= 0.0:
        return 1.0
    from scipy.stats import chi2

    return float(chi2.sf(value, df))


def _lineage_labels(gi) -> List[str]:
    """The lineage labels worth stratifying over, sorted for determinism.

    ``unknown`` - the sentinel :data:`papipeline.stages.convergence.UNKNOWN_LINEAGE`
    - is a missing label, not a lineage, and a cohort of nothing but missing
    labels has no strata to pool.
    """
    return sorted(
        {
            label
            for label in gi.lineages.values()
            if label and label != UNKNOWN_LINEAGE
        }
    )


def _refuse_thin_lineages(gi, config) -> None:
    """Refuse a cohort that has fewer than two lineages to pool across."""
    labels = _lineage_labels(gi)
    if len(labels) >= 2:
        return
    section = config.raw.get("lineage") or {}
    configured = section.get("method")
    raise DataContractError(
        "The lineage-stratified meta-analysis pools ACROSS lineages and needs "
        f"at least two; this cohort has {len(labels)}. Set {LINEAGE_METHOD_KEY} "
        "in config/science.yaml so stage 9 labels the isolates, or run the "
        "feature-level scans, which do not need lineages.",
        method=LINEAGE_METHOD_KEY,
        configured_method=str(configured) if configured else "unset",
        n_lineages=len(labels),
        n_labelled=sum(
            1
            for sample_id in gi.sample_ids
            if gi.lineages.get(sample_id)
            and gi.lineages[sample_id] != UNKNOWN_LINEAGE
        ),
        n_samples=gi.n_samples,
    )


def stratum_estimates(gi, feature: str) -> Tuple[StratumEstimate, ...]:
    """Estimate ``feature``'s effect inside each lineage, usable strata only.

    A stratum is usable when it has at least one carrier, at least one
    non-carrier and both outcome classes present - the four cells of an odds
    ratio that can actually be formed. Unusable strata are skipped silently
    here (the caller reports how many were usable in ``n_lineages``); a
    feature that is simply not in the matrix is a refusal, because a typo must
    not come back as an empty tuple.

    Returns:
        Estimates in sorted lineage order.
    """
    if feature not in gi.features:
        raise DataContractError(
            "Feature is absent from the feature matrix, so there is nothing to "
            "estimate per lineage",
            feature=feature,
            n_features=len(gi.features),
        )

    vector = gi.features[feature]
    groups: Dict[str, List[str]] = {}
    for sample_id in gi.sample_ids:
        lineage = gi.lineages.get(sample_id) or UNKNOWN_LINEAGE
        if lineage == UNKNOWN_LINEAGE:
            continue
        groups.setdefault(lineage, []).append(sample_id)

    estimates: List[StratumEstimate] = []
    for lineage in sorted(groups):
        samples = groups[lineage]
        carrier_r = sum(
            1 for s in samples if vector.get(s, 0) and gi.binary_outcome[s] == 1
        )
        carrier_s = sum(
            1 for s in samples if vector.get(s, 0) and gi.binary_outcome[s] == 0
        )
        non_r = sum(
            1 for s in samples if not vector.get(s, 0) and gi.binary_outcome[s] == 1
        )
        non_s = sum(
            1 for s in samples if not vector.get(s, 0) and gi.binary_outcome[s] == 0
        )
        n_positive = carrier_r + non_r
        n_negative = carrier_s + non_s
        n_carriers = carrier_r + carrier_s
        n_non = non_r + non_s

        if n_carriers < 1 or n_non < 1 or n_positive < 1 or n_negative < 1:
            # No odds ratio is definable: no carriers, no non-carriers, or one
            # outcome class only. Dropping the stratum is the honest answer;
            # a zero-corrected estimate from a degenerate table is not.
            continue

        cells = [carrier_r, carrier_s, non_r, non_s]
        corrected = min(cells) == 0
        if corrected:
            cells = [cell + HALDANE_ANSCOMBE for cell in cells]
        carrier_r, carrier_s, non_r, non_s = cells

        effect = math.log((carrier_r * non_s) / (non_r * carrier_s))
        std_error = math.sqrt(
            1.0 / carrier_r + 1.0 / carrier_s + 1.0 / non_r + 1.0 / non_s
        )
        estimates.append(
            StratumEstimate(
                lineage=lineage,
                n_samples=len(samples),
                n_positive=n_positive,
                n_negative=n_negative,
                n_carriers=n_carriers,
                effect=effect,
                std_error=std_error,
                p_value=_two_sided_p(effect / std_error),
                corrected=corrected,
            )
        )
    return tuple(estimates)


def _pool(estimates: Sequence[StratumEstimate]) -> Dict[str, Optional[float]]:
    """Fixed-effect inverse-variance pooling, with Q and I squared.

    Returns every pooled quantity in one dict so the statuses can be assigned
    by the caller: ``i_squared`` is undefined when Q is zero (no spread at
    all), which is why it is computed here rather than divided by zero later.
    """
    weights = [1.0 / (estimate.std_error**2) for estimate in estimates]
    total_weight = sum(weights)
    pooled_effect = sum(
        w * estimate.effect for w, estimate in zip(weights, estimates)
    ) / total_weight
    pooled_std_error = math.sqrt(1.0 / total_weight)
    mean = pooled_effect
    q_statistic = sum(
        w * (estimate.effect - mean) ** 2 for w, estimate in zip(weights, estimates)
    )
    df = len(estimates) - 1
    if q_statistic > 0.0:
        i_squared = max(0.0, (q_statistic - df) / q_statistic) * 100.0
    else:
        i_squared = 0.0
    return {
        "pooled_effect": pooled_effect,
        "pooled_std_error": pooled_std_error,
        "p_value": _two_sided_p(pooled_effect / pooled_std_error),
        "q_statistic": q_statistic,
        "q_p_value": _chi2_sf(q_statistic, df),
        "i_squared": i_squared,
    }


def lineage_meta_analysis(
    gi, config, features: Optional[Sequence[str]] = None
) -> List[MetaResult]:
    """Pool each feature across lineages.

    Args:
        gi: the constructed GWAS input, with stage-9 lineages on it.
        config: the pipeline configuration; only its ``lineage.method``
            setting is read, to say in the refusal how this cohort was meant
            to be labelled.
        features: explicit feature list. Every name must be in the matrix - a
            typo must be refused by name, not returned as an empty row.

    Raises:
        DataContractError: fewer than two lineage labels in the cohort, or an
            explicitly named feature that is not in the matrix.

    Returns:
        One :class:`MetaResult` per feature, sorted by feature name.
    """
    _refuse_thin_lineages(gi, config)

    if features is None:
        selected = list(gi.testable_features())
    else:
        selected = list(features)
        absent = [name for name in selected if name not in gi.features]
        if absent:
            raise DataContractError(
                "Feature is absent from the feature matrix, so it has no "
                "lineages to pool",
                feature=",".join(absent),
                n_features=len(gi.features),
            )

    results: List[MetaResult] = []
    for feature in sorted(set(selected)):
        estimates = stratum_estimates(gi, feature)
        if len(estimates) < 2:
            results.append(
                MetaResult(
                    feature=feature,
                    status=STATUS_INSUFFICIENT,
                    n_lineages=len(estimates),
                    estimates=estimates,
                )
            )
            continue
        pooled = _pool(estimates)
        results.append(
            MetaResult(
                feature=feature,
                status=STATUS_TESTED,
                n_lineages=len(estimates),
                estimates=estimates,
                pooled_effect=pooled["pooled_effect"],
                pooled_std_error=pooled["pooled_std_error"],
                p_value=pooled["p_value"],
                q_statistic=pooled["q_statistic"],
                q_p_value=pooled["q_p_value"],
                i_squared=pooled["i_squared"],
            )
        )
    return results
