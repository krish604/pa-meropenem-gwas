"""The binary-outcome fits the conditional and interaction scans share.

Two models, both plain logistic regression:

* **conditional scan** - a likelihood-ratio test of the model with one extra
  feature column against the model with only the known covariates;
* **interaction tiers** - one model carrying both main effects and their
  product, so the interaction term's Wald standard error comes from the
  observed information at that point.

``numpy`` and ``scipy`` are declared dependencies of this repository
(``environment/environment.yml``); nothing here adds one. ``statsmodels`` is
NOT a declared dependency and is deliberately not imported.

What makes this module worth having rather than calling a fitter directly:
the uninformative designs this package meets in practice - a feature that lies
exactly in the span of the covariates, a product column that is identically
zero, a feature carried by one isolate - must produce a NULL RESULT, not an
exception and not a spuriously small p-value. :func:`in_span`,
:func:`lrt_pvalue`'s clamp and :func:`wald`'s undefined-standard-error guard
are the three places that guarantee it.

A note on ``fit.converged``: under quasi-separation the maximum likelihood is
approached at a coefficient of unbounded magnitude, so BFGS regularly reports
``Desired error not necessarily achieved due to precision loss`` while sitting
at the likelihood supremum - which is exactly the point the likelihood-ratio
test needs. Callers therefore read ``log_likelihood`` and are told
``converged`` is informational, never a reason to drop a fit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import chi2, norm

#: Optimiser iteration cap. Large enough that a 60-sample design reaches its
#: supremum, small enough that a pathological matrix cannot hang a run.
MAX_ITERATIONS = 500

#: Below this residual, a column is treated as lying in the span of the basis.
#: The conditional test's null case is an EXACT duplicate column, so this only
#: has to be tighter than float noise on a 0/1 matrix.
SPAN_TOLERANCE = 1e-9

#: ``|coefficient|`` at or above which a fit is recorded as separated. Nothing
#: raises on this: a separated fit still has a likelihood, which is what the
#: likelihood-ratio test compares.
SEPARATION_BOUND = 30.0

#: How the model is described in a result's ``model`` field. The suffix is the
#: repository's own convention for "this p-value is not corrected for
#: population structure" - see ``docs/scientific_rules.md`` section 3.
MODEL = "logistic:NO_KINSHIP_CORRECTION"


@dataclass(frozen=True)
class LogitFit:
    """One fitted logistic model.

    ``converged`` is the optimiser's own verdict and is informational.
    ``separated`` says the coefficients ran to the bound, i.e. the maximum
    likelihood is not finite - the likelihood is still the right one to
    compare.
    """

    coefficients: Tuple[float, ...]
    log_likelihood: float
    converged: bool
    separated: bool

    @property
    def n_parameters(self) -> int:
        return len(self.coefficients)


@dataclass(frozen=True)
class WaldTerm:
    """One coefficient's Wald estimate, standard error and two-sided p-value.

    ``std_error`` is ``None`` when the observed information is singular for
    that coefficient - the estimate is then not identifiable, and a p-value
    computed from a pseudo-inverse would be a number without a meaning.
    """

    estimate: float
    std_error: Optional[float]
    p_value: Optional[float]


def design_matrix(n_samples: int, *columns: Sequence[float]) -> np.ndarray:
    """Intercept column followed by ``columns``, as a float matrix.

    The sample count is an argument rather than inferred, so an intercept-only
    design (no covariates at all) is a valid call instead of an empty one.
    """
    matrix = np.column_stack(
        [np.ones(int(n_samples))]
        + [np.asarray(c, dtype=float) for c in columns]
    )
    return matrix


def log_likelihood(design: np.ndarray, y: np.ndarray, coefficients: np.ndarray) -> float:
    """Log-likelihood, computed through ``logaddexp`` so large eta cannot NaN."""
    eta = design @ coefficients
    return float(np.sum(y * eta - np.logaddexp(0.0, eta)))


def fit_logit(design: np.ndarray, y: np.ndarray) -> LogitFit:
    """Fit ``P(y=1) = sigmoid(design @ beta)`` by BFGS on the negative
    log-likelihood, with an analytic gradient."""
    design = np.asarray(design, dtype=float)
    y = np.asarray(y, dtype=float)

    def objective(beta: np.ndarray) -> float:
        eta = design @ beta
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    def gradient(beta: np.ndarray) -> np.ndarray:
        eta = design @ beta
        return design.T @ (expit(eta) - y)

    result = minimize(
        objective,
        np.zeros(design.shape[1]),
        jac=gradient,
        method="BFGS",
        options={"maxiter": MAX_ITERATIONS, "gtol": 1e-10},
    )
    coefficients = np.asarray(result.x, dtype=float)
    separated = bool(np.max(np.abs(coefficients)) >= SEPARATION_BOUND) if coefficients.size else False
    return LogitFit(
        coefficients=tuple(float(v) for v in coefficients),
        log_likelihood=-objective(coefficients),
        converged=bool(result.success),
        separated=separated,
    )


def in_span(column: Sequence[float], basis: np.ndarray) -> bool:
    """Whether ``column`` is a linear combination of ``basis``'s columns.

    The basis normally already carries an intercept, so "the feature adds
    nothing the covariates do not already have" is a span question, and
    answering it before fitting is what turns a rank-deficient optimisation
    into an exact ``p = 1.0``.
    """
    values = np.asarray(column, dtype=float)
    basis = np.asarray(basis, dtype=float)
    if basis.size == 0 or basis.shape[1] == 0:
        return bool(np.max(np.abs(values)) <= SPAN_TOLERANCE)
    coefficients, *_ = np.linalg.lstsq(basis, values, rcond=None)
    residual = values - basis @ coefficients
    return bool(np.max(np.abs(residual)) <= SPAN_TOLERANCE)


def lrt_pvalue(full: LogitFit, reduced: LogitFit, degrees_of_freedom: int = 1) -> float:
    """Likelihood-ratio p-value for ``full`` against ``reduced``.

    The statistic is clamped at zero: under an exactly nested pair the two
    likelihoods are equal up to float noise, and a tiny negative statistic
    would become a p-value above 1.
    """
    statistic = 2.0 * (full.log_likelihood - reduced.log_likelihood)
    if statistic < 0.0:
        statistic = 0.0
    return float(chi2.sf(statistic, degrees_of_freedom))


def wald(design: np.ndarray, fit: LogitFit, index: int) -> WaldTerm:
    """Wald estimate, standard error and two-sided p-value for one coefficient.

    The standard error comes from the observed information at the fitted
    point. When that matrix is singular the pseudo-inverse is used, which
    keeps the well-identified coefficients correct; a coefficient the
    information matrix cannot separate comes back with ``std_error=None``
    rather than a fabricated p-value.
    """
    coefficients = np.asarray(fit.coefficients, dtype=float)
    eta = design @ coefficients
    probabilities = expit(eta)
    weights = probabilities * (1.0 - probabilities)
    information = (design.T * weights) @ design
    try:
        covariance = np.linalg.inv(information)
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(information)
    variance = float(covariance[index, index])
    estimate = float(coefficients[index])
    if not np.isfinite(variance) or variance <= 0.0:
        return WaldTerm(estimate=estimate, std_error=None, p_value=None)
    std_error = float(np.sqrt(variance))
    z_score = estimate / std_error
    return WaldTerm(
        estimate=estimate,
        std_error=std_error,
        p_value=float(2.0 * norm.sf(abs(z_score))),
    )
