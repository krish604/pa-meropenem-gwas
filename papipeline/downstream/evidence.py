"""Deliverable 5: the novelty filter and the evidence tiers.

What is "already known" is ``config/known_determinants.tsv``, derived from
``config/mechanisms.tsv`` - so a feature that table resolves to is not a
novelty candidate, and one it does not resolve to is. The partition is a
filter, not a claim: novelty says nothing about whether the association is
real, only about whether this project had already written it down.

Grading runs on two scans at once, because they answer different questions:

* **A** - significant AFTER conditioning on the known layer features, and seen
  in at least the configured number of independent lineages. The claim that
  travels with it is ``SUPPORTED``, and nothing stronger: causation is out of
  scope for this pipeline (``ClaimStatus`` has no CAUSAL member).
* **B** - significant after conditioning but below the lineage floor, i.e.
  the association is there in this cohort but is carried by too few lineages
  to separate from lineage. ``ASSOCIATED``.
* **C** - significant in the baseline scan ONLY: the association did not
  survive conditioning on what is already known. ``ASSOCIATED``, and the
  ``marginal_only`` basis says why it is the weaker claim.
* **D** - not significant in either scan. ``UNKNOWN``.

A missing adjusted p-value is not a significant one: a feature the engine
skipped (too rare, or a covariate) grades D rather than borrowing a p-value
from somewhere.

Two thresholds are read from the raw science mapping rather than from the
derived config fields, because ``dataclasses.replace(config, raw=...)`` is how
this repository overrides a science key and it does NOT rebuild the derived
sub-configs - a reader that took ``config.gwas`` would silently keep the
committed value for an overridden run.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple

from ..models import ClaimStatus
from ..stages.convergence import UNKNOWN_LINEAGE
from .known_determinants import KnownDeterminants

#: The evidence tiers, in strength order. Four, and only four.
EVIDENCE_TIERS: Tuple[str, ...] = ("A", "B", "C", "D")

#: The repository claim level each tier may report. No tier reports a level
#: ``ClaimStatus`` does not have.
TIER_CLAIMS: Mapping[str, str] = {
    "A": ClaimStatus.SUPPORTED.value,
    "B": ClaimStatus.ASSOCIATED.value,
    "C": ClaimStatus.ASSOCIATED.value,
    "D": ClaimStatus.UNKNOWN.value,
}

#: Why a tier got the claim it got - recorded on the row so the reasoning
#: travels with the number instead of living in a module docstring.
TIER_BASIS: Mapping[str, str] = {
    "A": "conditional_and_convergent",
    "B": "conditional_only",
    "C": "marginal_only",
    "D": "not_significant",
}

#: Novelty vocabulary: a feature is one or the other, never "sort of novel".
NOVEL = "novel"
KNOWN = "known"


@dataclass(frozen=True)
class EvidenceRecord:
    """One feature's novelty, tier, claim and the two p-values behind them."""

    feature: str
    novelty: str
    evidence_tier: str
    claim_status: str
    basis: str
    known_matches: Tuple[str, ...] = ()
    independent_lineages: int = 0
    conditional_adjusted_p: Optional[float] = None
    baseline_adjusted_p: Optional[float] = None


def significance_threshold(config) -> float:
    """``gwas.significance_threshold``, read from the raw science mapping.

    Falls back to the derived field only when the key is absent, which the
    committed science file never is - the fallback exists so a hand-built
    config with no raw section still gets a threshold rather than an error.
    """
    section = config.raw.get("gwas") or {}
    raw_value = section.get("significance_threshold")
    if raw_value is None:
        return float(config.gwas.significance_threshold)
    return float(raw_value)


def lineage_floor(config) -> int:
    """``convergence.min_independent_lineages``: the tier-A lineage floor.

    Same raw-first reading as :func:`significance_threshold`, for the same
    reason.
    """
    section = config.raw.get("convergence") or {}
    raw_value = section.get("min_independent_lineages")
    if raw_value is None:
        return int(config.convergence.min_independent_lineages)
    return int(raw_value)


def novelty_partition(
    features: Sequence[str], known: KnownDeterminants
) -> Tuple[List[str], List[str]]:
    """Split features into (known, novel), preserving the caller's order.

    Order is preserved on both sides because the report and the tests both
    read these as lists; sorting would silently reorder a caller's table.
    """
    known_hits: List[str] = []
    novel: List[str] = []
    for feature in features:
        (known_hits if known.is_known(feature) else novel).append(feature)
    return known_hits, novel


def independent_lineages(result, convergence: Optional[Mapping] = None) -> int:
    """How many lineages carry this feature.

    A convergence call, when there is one, is the authority: it is the stage
    that actually looked at branches, and a distribution can disagree with it
    when the two were computed on different inputs. Otherwise the count is the
    number of lineages with at least one carrier in
    ``result.lineage_distribution``, excluding ``unknown`` - a missing label
    is not a lineage - and excluding lineages with zero carriers, which are
    present in the distribution only as a zero.

    Returns:
        0 for a result with no distribution, which is "not observed in any
        lineage", not "all of them".
    """
    if convergence:
        call = convergence.get(getattr(result, "feature", None))
        if call is not None:
            return int(call.independent_lineages)
    distribution = getattr(result, "lineage_distribution", None) or {}
    return sum(
        1
        for lineage, count in distribution.items()
        if lineage and lineage != UNKNOWN_LINEAGE and count and count > 0
    )


def classify_feature(
    feature: str,
    *,
    conditional_p: Optional[float],
    baseline_p: Optional[float],
    independent_lineages: int,
    known: KnownDeterminants,
    config,
) -> EvidenceRecord:
    """Grade one feature. Keyword-only, because every input is a statement.

    Args:
        conditional_p: adjusted p from the conditional scan, or ``None`` when
            that scan did not test the feature.
        baseline_p: adjusted p from the baseline scan, likewise.
        independent_lineages: from :func:`independent_lineages`.
        known: the novelty filter.
        config: supplies the significance threshold and the lineage floor.
    """
    threshold = significance_threshold(config)
    floor = lineage_floor(config)

    def significant(value: Optional[float]) -> bool:
        return value is not None and value <= threshold

    if significant(conditional_p):
        tier = "A" if independent_lineages >= floor else "B"
    elif significant(baseline_p):
        tier = "C"
    else:
        tier = "D"

    is_known = known.is_known(feature)
    return EvidenceRecord(
        feature=feature,
        novelty=KNOWN if is_known else NOVEL,
        evidence_tier=tier,
        claim_status=TIER_CLAIMS[tier],
        basis=TIER_BASIS[tier],
        known_matches=tuple(known.members_for(feature) or ()) if is_known else (),
        independent_lineages=independent_lineages,
        conditional_adjusted_p=conditional_p,
        baseline_adjusted_p=baseline_p,
    )


def build_evidence_table(
    *,
    conditional_results: Iterable,
    baseline_results: Iterable,
    known: KnownDeterminants,
    config,
    convergence: Optional[Mapping] = None,
) -> List[EvidenceRecord]:
    """Grade the union of both scans, one record per feature.

    The lineage count comes from the BASELINE result when the feature is in
    both scans, because the baseline scan is the full-cohort one and its
    distribution covers every carrier; the conditional record is used only
    when the baseline never saw the feature.

    Features are ordered baseline-first then conditional, which is the order a
    reader checks: what the plain scan found, then what conditioning did to it.
    """
    baseline = list(baseline_results)
    conditional = list(conditional_results)
    baseline_by_feature = {row.feature: row for row in baseline}
    conditional_by_feature = {row.feature: row for row in conditional}

    order: List[str] = list(dict.fromkeys(
        [row.feature for row in baseline]
        + [row.feature for row in conditional]
    ))

    records: List[EvidenceRecord] = []
    for feature in order:
        baseline_row = baseline_by_feature.get(feature)
        conditional_row = conditional_by_feature.get(feature)
        lineage_source = baseline_row if baseline_row is not None else conditional_row
        records.append(
            classify_feature(
                feature,
                conditional_p=(
                    conditional_row.adjusted_p_value
                    if conditional_row is not None
                    else None
                ),
                baseline_p=(
                    baseline_row.adjusted_p_value if baseline_row is not None else None
                ),
                independent_lineages=(
                    independent_lineages(lineage_source, convergence)
                    if lineage_source is not None
                    else 0
                ),
                known=known,
                config=config,
            )
        )
    return records
