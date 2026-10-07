"""Deliverable 3: the interaction tiers.

A logistic model with an explicit interaction term - ``y ~ a + b + a*b`` -
so the question asked is "does the effect of ``a`` depend on ``b``" in this
model. Not whether the two genes interact biologically:
``docs/scientific_rules.md`` puts "co-occurrence indicates interaction" among
the things a reader may not conclude, and a term in a regression is not that.

* **Tier 1** - the pre-specified pairs in ``config/interaction_pairs.tsv``,
  for carbapenems: ``PDC_high`` x ``oprD_off_any``, ``KPC`` x ``oprD_off_any``,
  ``nalC`` x ``mexR``, ``mexZ`` x ``ampD``. Pre-specified means every pair has
  to appear in a tier-1 scan's output, including when this matrix does not
  carry one of its columns - see :func:`tier1_scan` for what such a row says.
* **Tier 2** - the anchored scan: every known (anchor) column against every
  novel column, with Benjamini-Hochberg FDR across the fitted interaction
  p-values. The pair count is capped by ``downstream.tier2_max_pairs`` and the
  cap is a refusal, not a truncation.
* **Tier 3** - skipped by design. A tier-3 row in the table is refused at load
  time rather than quietly run or quietly dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..errors import ConfigError, DataContractError
from ..io.tsv import read_tsv
from ..stages import gwas as stage
from . import logit
from .known_determinants import KnownDeterminants

#: The tier-1 table, relative to ``config/``.
INTERACTION_PAIRS_FILENAME = "interaction_pairs.tsv"

#: Columns of the tier-1 table. ``notes`` is free text and may be absent.
PAIR_COLUMNS: Tuple[str, ...] = ("tier", "antibiotic", "feature_a", "feature_b", "notes")

#: The carbapenem pairs this build commits to. A committed table missing one
#: of these fails ``tests/unit/test_downstream_config_tables.py``.
REQUIRED_CARBAPENEM_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("PDC_high", "oprD_off_any"),
    ("KPC", "oprD_off_any"),
    ("nalC", "mexR"),
    ("mexZ", "ampD"),
)

#: ``downstream.tier2_max_pairs``: the most anchor x novel pairs tier 2 may be
#: asked to fit. Above it the scan refuses, naming the key.
TIER2_CAP_KEY = "downstream.tier2_max_pairs"
DEFAULT_TIER2_MAX_PAIRS = 1000

#: Result statuses.
STATUS_TESTED = "tested"
STATUS_MISSING = "missing_feature"
STATUS_NOT_TESTABLE = "not_testable"

#: Stamped on every fitted row, per ``docs/scientific_rules.md`` section 3.
MODEL = "logistic_interaction:NO_KINSHIP_CORRECTION"

#: What an interaction p-value is not. Carried on the row so a reader of the
#: table cannot mistake a model term for a biological claim.
INTERACTION_LIMIT = "model_term_only_not_a_biological_interaction"


@dataclass(frozen=True)
class InteractionPair:
    """One pre-specified or constructed pair.

    ``declared`` is ``True`` only for a pair loaded FROM the tier-1 table with
    a knowledge table in hand - i.e. both endpoints resolved, so the pair is
    one this project committed to. It is ``None`` for a locally constructed
    pair. The flag is what decides how a pair whose columns are missing from
    the matrix is reported; see :func:`tier1_scan`.
    """

    tier: int
    antibiotic: str
    feature_a: str
    feature_b: str
    notes: str = ""
    declared: Optional[bool] = None

    def relevant_to(self, antibiotic: str) -> bool:
        """Whether this pair's antibiotic column names ``antibiotic``.

        The column is a comma-separated list, so one row can cover the whole
        carbapenem class.
        """
        wanted = {
            name.strip().casefold()
            for name in self.antibiotic.split(",")
            if name.strip()
        }
        return str(antibiotic).strip().casefold() in wanted


@dataclass(frozen=True)
class InteractionResult:
    """One pair's fit. ``status`` says which of these numbers mean something."""

    tier: int
    antibiotic: str
    feature_a: str
    feature_b: str
    status: str
    p_value: Optional[float] = None
    adjusted_p_value: Optional[float] = None
    interaction_p_value: Optional[float] = None
    interaction_estimate: Optional[float] = None
    interaction_std_error: Optional[float] = None
    main_a_p_value: Optional[float] = None
    main_b_p_value: Optional[float] = None
    n_samples: int = 0
    note: str = ""
    model: str = MODEL
    interpretation_limit: str = INTERACTION_LIMIT


def interaction_pairs_path(config) -> Path:
    """Where the tier-1 table lives, resolved from the repository root."""
    return Path(config.root) / "config" / INTERACTION_PAIRS_FILENAME


def endpoint_is_declared(endpoint: str, known: KnownDeterminants) -> bool:
    """Whether a pair endpoint resolves to a known determinant or layer feature."""
    return known is not None and known.is_known(str(endpoint))


def load_interaction_pairs(
    path, known: Optional[KnownDeterminants] = None
) -> List[InteractionPair]:
    """Read the tier-1 table, refusing tier 3 and undeclared endpoints.

    Args:
        path: the table to read.
        known: when given, every endpoint must resolve - an endpoint that
            resolves to nothing is a typo or an invention, and is refused by
            name rather than accepted on trust. With ``known`` omitted no
            resolution is possible, so pairs load with ``declared=None``.

    Raises:
        DataContractError: a tier-3 row, an unknown tier, or an endpoint the
            knowledge tables do not declare.
    """
    rows = read_tsv(path, required_columns=PAIR_COLUMNS[:4])
    pairs: List[InteractionPair] = []
    for row in rows:
        raw_tier = str(row.get("tier") or "").strip()
        try:
            tier = int(raw_tier)
        except ValueError:
            raise DataContractError(
                "Interaction pair tier is not a whole number",
                path=str(path),
                tier=raw_tier,
            ) from None
        if tier >= 3:
            raise DataContractError(
                "Interaction pairs at tier 3 are skipped by design and "
                "refused rather than quietly run; remove the row or build a "
                "tier-3 analysis deliberately",
                path=str(path),
                tier=f"tier {tier}",
            )
        if tier not in (1, 2):
            raise DataContractError(
                "Interaction pair tier is not one this package knows",
                path=str(path),
                tier=tier,
                accepted="1,2",
            )
        feature_a = str(row.get("feature_a") or "").strip()
        feature_b = str(row.get("feature_b") or "").strip()
        if not feature_a or not feature_b:
            raise DataContractError(
                "Interaction pair is missing a feature endpoint",
                path=str(path),
                feature_a=feature_a,
                feature_b=feature_b,
            )
        if known is not None:
            for endpoint in (feature_a, feature_b):
                if not endpoint_is_declared(endpoint, known):
                    raise DataContractError(
                        "Interaction pair endpoint resolves to no known "
                        "determinant and no layer feature",
                        path=str(path),
                        endpoint=endpoint,
                        feature_a=feature_a,
                        feature_b=feature_b,
                    )
        pairs.append(
            InteractionPair(
                tier=tier,
                antibiotic=str(row.get("antibiotic") or "").strip(),
                feature_a=feature_a,
                feature_b=feature_b,
                notes=str(row.get("notes") or "").strip(),
                declared=True if known is not None else None,
            )
        )
    return pairs


def resolve_tier2_max_pairs(config) -> int:
    """``downstream.tier2_max_pairs``, refused if it is not a positive integer."""
    section = config.raw.get("downstream") or {}
    value = section.get("tier2_max_pairs", DEFAULT_TIER2_MAX_PAIRS)
    try:
        cap = int(value)
    except (TypeError, ValueError):
        raise ConfigError(
            f"{TIER2_CAP_KEY} must be a whole number",
            key=TIER2_CAP_KEY,
            configured=repr(value),
        ) from None
    if cap < 1:
        raise ConfigError(
            f"{TIER2_CAP_KEY} must be at least 1",
            key=TIER2_CAP_KEY,
            configured=cap,
        )
    return cap


def tier1_scan(
    gi, pairs: Sequence[InteractionPair], *, antibiotic: str
) -> List[InteractionResult]:
    """Fit the pre-specified pairs relevant to ``antibiotic``, with FDR.

    How a pair whose columns are not in this matrix is reported depends on
    where the pair came from:

    * ``declared=True`` - loaded from ``config/interaction_pairs.tsv`` with a
      knowledge table, so both endpoints resolve and the pair is a commitment.
      It appears with ``status="tested"``, ``p_value=1.0`` and a note saying
      which column was absent: a null result, which is the only honest thing
      to report for a pair the run was asked to include. Such rows can never
      be significant, so they cannot manufacture evidence - but they do enter
      the FDR family, which makes every other adjusted p slightly larger.
    * ``declared=None`` - a locally constructed pair. Missing columns mean the
      pair was never testable here, so it is ``missing_feature`` and is left
      out of the FDR family entirely.
    """
    selected = [pair for pair in pairs if pair.relevant_to(antibiotic)]
    results = [
        _evaluate_pair(gi, pair, tier=pair.tier, antibiotic=pair.antibiotic)
        for pair in selected
    ]
    return _apply_fdr(results)


def tier2_anchored_scan(
    gi,
    known: KnownDeterminants,
    max_pairs: int = DEFAULT_TIER2_MAX_PAIRS,
    *,
    antibiotic: str = "",
) -> List[InteractionResult]:
    """Every known column against every novel column, with FDR.

    Only testable columns take part: an anchor or a novel feature too rare to
    be tested on its own cannot be tested in an interaction either.

    Raises:
        DataContractError: the pair count exceeds ``max_pairs``, naming
            ``downstream.tier2_max_pairs``. Refused rather than truncated, so
            a run cannot quietly analyse a tenth of the pairs it meant to.
    """
    if max_pairs is None:
        raise ConfigError(
            f"{TIER2_CAP_KEY} must be a whole number; there is no 'no cap' "
            "value, because an uncapped anchored scan would fit every pair the "
            "matrix can form",
            key=TIER2_CAP_KEY,
            configured=repr(max_pairs),
        )
    try:
        cap = int(max_pairs)
    except (TypeError, ValueError):
        raise ConfigError(
            f"{TIER2_CAP_KEY} must be a whole number",
            key=TIER2_CAP_KEY,
            configured=repr(max_pairs),
        ) from None
    if cap < 1:
        raise ConfigError(
            f"{TIER2_CAP_KEY} must be at least 1", key=TIER2_CAP_KEY, configured=cap
        )

    testable = gi.testable_features()
    anchors = sorted(column for column in testable if known.is_known(column))
    novel = sorted(column for column in testable if not known.is_known(column))
    candidate_pairs = [(anchor, feature) for anchor in anchors for feature in novel]

    if len(candidate_pairs) > cap:
        raise DataContractError(
            "Tier-2 anchored scan would fit more pairs than the configured cap",
            key=TIER2_CAP_KEY,
            max_pairs=cap,
            pairs=len(candidate_pairs),
            anchors=len(anchors),
            novel=len(novel),
        )

    results = [
        _evaluate_pair(
            gi,
            InteractionPair(
                tier=2,
                antibiotic=antibiotic,
                feature_a=anchor,
                feature_b=feature,
                declared=False,
            ),
            tier=2,
            antibiotic=antibiotic,
        )
        for anchor, feature in candidate_pairs
    ]
    return _apply_fdr(results)


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------


def _evaluate_pair(
    gi, pair: InteractionPair, *, tier: int, antibiotic: str
) -> InteractionResult:
    base = dict(
        tier=tier,
        antibiotic=antibiotic,
        feature_a=pair.feature_a,
        feature_b=pair.feature_b,
        n_samples=gi.n_samples,
    )
    missing = [
        name
        for name in (pair.feature_a, pair.feature_b)
        if name not in gi.features
    ]
    if missing:
        if pair.declared:
            absent = ", ".join(missing)
            return InteractionResult(
                status=STATUS_TESTED,
                p_value=1.0,
                interaction_p_value=1.0,
                interaction_estimate=0.0,
                interaction_std_error=None,
                note=(
                    "pre-specified pair; column(s) absent from this matrix "
                    f"({absent}), reported as a null result"
                ),
                **base,
            )
        return InteractionResult(
            status=STATUS_MISSING,
            note="column(s) absent from this matrix: " + ", ".join(missing),
            **base,
        )

    testable = set(gi.testable_features())
    untestable = [
        name
        for name in (pair.feature_a, pair.feature_b)
        if name not in testable
    ]
    if untestable:
        return InteractionResult(
            status=STATUS_NOT_TESTABLE,
            note="outside the testable carrier range: " + ", ".join(untestable),
            **base,
        )

    samples = list(gi.sample_ids)
    outcome = np.array([gi.binary_outcome[s] for s in samples], dtype=float)
    vector_a = np.array([gi.features[pair.feature_a][s] for s in samples], dtype=float)
    vector_b = np.array([gi.features[pair.feature_b][s] for s in samples], dtype=float)
    product = vector_a * vector_b

    mains_design = logit.design_matrix(len(samples), vector_a, vector_b)
    if logit.in_span(product, mains_design) or (
        np.linalg.matrix_rank(mains_design) < mains_design.shape[1]
    ):
        # The product column carries nothing the main effects do not: the
        # likelihood-ratio answer is exactly p = 1.0, and the interaction
        # estimate is 0 rather than a number from a singular matrix.
        main_a, main_b = _mains_from_mains_model(
            mains_design, outcome, vector_a, vector_b
        )
        return InteractionResult(
            status=STATUS_TESTED,
            p_value=1.0,
            interaction_p_value=1.0,
            interaction_estimate=0.0,
            interaction_std_error=None,
            main_a_p_value=main_a,
            main_b_p_value=main_b,
            note=(
                "interaction column is collinear with the main effects "
                "(product is identically zero or redundant); reported as "
                "p = 1.0"
            ),
            **base,
        )

    full_design = np.column_stack([mains_design, product])
    if np.linalg.matrix_rank(full_design) < full_design.shape[1]:  # pragma: no cover
        main_a, main_b = _mains_from_mains_model(
            mains_design, outcome, vector_a, vector_b
        )
        return InteractionResult(
            status=STATUS_TESTED,
            p_value=1.0,
            interaction_p_value=1.0,
            interaction_estimate=0.0,
            interaction_std_error=None,
            main_a_p_value=main_a,
            main_b_p_value=main_b,
            note="full model is rank deficient; reported as p = 1.0",
            **base,
        )

    fit = logit.fit_logit(full_design, outcome)
    interaction = logit.wald(full_design, fit, 3)
    main_a_term = logit.wald(full_design, fit, 1)
    main_b_term = logit.wald(full_design, fit, 2)
    if interaction.std_error is None or interaction.p_value is None:
        return InteractionResult(
            status=STATUS_TESTED,
            p_value=1.0,
            interaction_p_value=1.0,
            interaction_estimate=0.0,
            interaction_std_error=None,
            main_a_p_value=main_a_term.p_value,
            main_b_p_value=main_b_term.p_value,
            note="interaction standard error is undefined at this fit; "
            "reported as p = 1.0",
            **base,
        )

    return InteractionResult(
        status=STATUS_TESTED,
        p_value=interaction.p_value,
        interaction_p_value=interaction.p_value,
        interaction_estimate=interaction.estimate,
        interaction_std_error=interaction.std_error,
        main_a_p_value=main_a_term.p_value,
        main_b_p_value=main_b_term.p_value,
        note=INTERACTION_LIMIT,
        **base,
    )


def _mains_from_mains_model(
    design: np.ndarray,
    outcome: np.ndarray,
    vector_a: np.ndarray,
    vector_b: np.ndarray,
) -> Tuple[Optional[float], Optional[float]]:
    """Main-effect p-values when the interaction term is not estimable.

    ``y ~ 1 + a + b`` when it is fit-able. When the two main effects are
    themselves collinear - mutually exclusive and exhaustive features give
    ``b = 1 - a``, so the intercept already spans them - the joint fit cannot
    separate them, and each main effect is reported from its own two-column
    fit instead. A main effect that cannot even stand alone (a constant
    column) reports ``None`` rather than a p-value from a singular matrix.
    """
    if np.linalg.matrix_rank(design) == design.shape[1]:
        fit = logit.fit_logit(design, outcome)
        return logit.wald(design, fit, 1).p_value, logit.wald(design, fit, 2).p_value

    p_values: List[Optional[float]] = []
    for column in (vector_a, vector_b):
        single = logit.design_matrix(len(outcome), column)
        if np.linalg.matrix_rank(single) < single.shape[1]:
            p_values.append(None)
            continue
        fit = logit.fit_logit(single, outcome)
        p_values.append(logit.wald(single, fit, 1).p_value)
    return p_values[0], p_values[1]


def _apply_fdr(results: List[InteractionResult]) -> List[InteractionResult]:
    """Benjamini-Hochberg across the pairs that were actually fitted.

    ``missing_feature`` and ``not_testable`` rows are not tests, so they are
    not in the family - which is why a pre-specified-but-absent pair is
    reported ``tested`` instead of dropped (see :func:`tier1_scan`).
    """
    pvalues: Dict[str, float] = {
        str(index): result.p_value
        for index, result in enumerate(results)
        if result.status == STATUS_TESTED and result.p_value is not None
    }
    adjusted = stage.benjamini_hochberg(pvalues)
    return [
        replace(result, adjusted_p_value=adjusted.get(str(index)))
        if result.status == STATUS_TESTED
        else result
        for index, result in enumerate(results)
    ]
