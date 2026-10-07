"""Deliverable 2: the conditional and the stratified scan.

**Conditional** - refit every candidate feature WITH the known layer features
as covariates, and compare against the covariate-only model by a
likelihood-ratio test. A feature that IS one of the covariates is reported
``skipped_covariate`` rather than silently dropped, so the output says which
columns were held fixed. A feature with too few carriers for
:func:`papipeline.stages.gwas.GwasInput.testable_features` is reported
``not_testable``. Both are rows, not omissions.

A feature lying exactly in the span of the covariates - the null case, and the
one a naive fit turns into a rank-deficient optimisation - gets the exact
likelihood-ratio answer of ``p = 1.0`` rather than a number from a singular
matrix. Multiple-testing correction is the repository's own
:func:`papipeline.stages.gwas.benjamini_hochberg`, applied across the features
that were actually tested.

**Stratified** - the same association machinery restricted to the isolates
that carry NONE of the known determinants: the "is anything left once the
known biology is accounted for" question. The stratum is rebuilt through
:func:`papipeline.stages.gwas.build_input`, so the group floor
(``gwas.min_samples_per_group``) and every input refusal apply to the subset
exactly as they apply to the cohort - a stratum too small to test is refused,
not analysed. Lineages do not enter the stratum definition; they are carried
through for the lineage distribution on each result.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..errors import ConfigError, DataContractError
from ..manifest import SampleManifest
from ..models import GwasResult, Phenotype, PhenotypeCall, Sample
from ..stages import gwas as stage
from . import logit
from .known_determinants import KnownDeterminants

#: Result statuses a conditional scan can emit.
STATUS_TESTED = "tested"
STATUS_SKIPPED = "skipped_covariate"
STATUS_NOT_TESTABLE = "not_testable"

#: The key extra covariates are read from, and the key a refusal names.
COVARIATES_KEY = "downstream.covariates"


@dataclass(frozen=True)
class ConditionalResult:
    """One feature's conditional verdict."""

    feature: str
    status: str
    p_value: Optional[float] = None
    adjusted_p_value: Optional[float] = None
    covariates: List[str] = field(default_factory=list)
    effect: Optional[float] = None
    n_samples: int = 0
    model: str = logit.MODEL


@dataclass(frozen=True)
class StratifiedScan:
    """The subset analysed, and what came out of it."""

    samples: Tuple[str, ...]
    known_columns: Tuple[str, ...]
    results: Tuple[GwasResult, ...] = ()


def _extra_covariates(config) -> List[str]:
    """``downstream.covariates`` from the science file, validated as a list."""
    if config is None:
        return []
    section = config.raw.get("downstream") or {}
    configured = section.get("covariates") or []
    if isinstance(configured, str) or not isinstance(configured, (list, tuple)):
        raise ConfigError(
            f"{COVARIATES_KEY} must be a list of column names",
            key=COVARIATES_KEY,
            configured=repr(configured),
        )
    return [str(name).strip() for name in configured if str(name).strip()]


def covariate_columns(
    gi, known: KnownDeterminants, config=None
) -> List[str]:
    """The columns the conditional scan holds fixed, sorted for determinism.

    Known layer features present in the matrix, plus whatever
    ``downstream.covariates`` adds. A configured covariate that is NOT in the
    matrix is a refusal: silently ignoring it would refit the model the
    operator asked for and report it as the model they asked for.
    """
    extra = _extra_covariates(config)
    missing = [name for name in extra if name not in gi.features]
    if missing:
        raise DataContractError(
            "A configured covariate is absent from the feature matrix",
            key=COVARIATES_KEY,
            columns=",".join(missing),
            n_features=len(gi.features),
        )
    known_columns = [column for column in gi.features if known.is_known(column)]
    return sorted(set(known_columns) | set(extra))


def conditional_scan(gi, covariates: Sequence[str]) -> List[ConditionalResult]:
    """Fit every testable feature against the covariate-only model.

    Args:
        gi: the constructed GWAS input.
        covariates: columns to hold fixed, normally from
            :func:`covariate_columns`.

    Returns:
        One :class:`ConditionalResult` per covariate and per testable feature,
        with Benjamini-Hochberg adjusted p-values across the tested ones.
    """
    covariate_list = list(covariates)
    absent = [name for name in covariate_list if name not in gi.features]
    if absent:
        raise DataContractError(
            "A covariate is absent from the feature matrix",
            key=COVARIATES_KEY,
            columns=",".join(absent),
            n_features=len(gi.features),
        )

    samples = list(gi.sample_ids)
    outcome = np.array([gi.binary_outcome[s] for s in samples], dtype=float)
    vectors = {
        name: np.array([gi.features[name][s] for s in samples], dtype=float)
        for name in set(covariate_list) | set(gi.testable_features())
    }
    reduced_design = logit.design_matrix(
        len(samples), *[vectors[name] for name in covariate_list]
    )
    reduced_fit = logit.fit_logit(reduced_design, outcome)

    testable = set(gi.testable_features())
    pvalues: Dict[str, float] = {}
    effects: Dict[str, Optional[float]] = {}
    statuses: Dict[str, str] = {}

    # Every column in the matrix gets a row, plus every covariate: a feature
    # too rare to test is reported `not_testable` rather than left out, so the
    # reader can tell "not tested" from "not there".
    for name in sorted(set(gi.features) | set(covariate_list)):
        if name in covariate_list:
            statuses[name] = STATUS_SKIPPED
            continue
        if name not in testable:
            statuses[name] = STATUS_NOT_TESTABLE
            continue
        vector = vectors[name]
        if logit.in_span(vector, reduced_design):
            # The covariates already carry everything the feature carries:
            # the two models are the same model, so the likelihood ratio is 0.
            pvalues[name] = 1.0
            effects[name] = None
        else:
            full_design = np.column_stack([reduced_design, vector])
            full_fit = logit.fit_logit(full_design, outcome)
            pvalues[name] = logit.lrt_pvalue(full_fit, reduced_fit)
            effects[name] = _odds_ratio(full_fit.coefficients[-1])
        statuses[name] = STATUS_TESTED

    adjusted = stage.benjamini_hochberg(pvalues)
    return [
        ConditionalResult(
            feature=name,
            status=status,
            p_value=pvalues.get(name),
            adjusted_p_value=adjusted.get(name),
            covariates=list(covariate_list),
            effect=effects.get(name),
            n_samples=gi.n_samples,
        )
        for name, status in sorted(statuses.items())
    ]


def _odds_ratio(coefficient: float) -> Optional[float]:
    """``exp(coefficient)``, or ``None`` when it overflows to infinity."""
    try:
        value = float(np.exp(coefficient))
    except (OverflowError, FloatingPointError):
        return None
    return value if np.isfinite(value) else None


def stratified_scan(gi, known_columns: Sequence[str], config) -> StratifiedScan:
    """Run the reference association inside the isolates lacking every known
    determinant.

    Args:
        gi: the constructed GWAS input for the whole cohort.
        known_columns: columns whose carriers are held out. A column named
            here but absent from the matrix is refused rather than ignored -
            ignoring it would leave carriers in a stratum that was supposed to
            exclude them.
        config: the pipeline configuration; the subset goes back through
            ``build_input`` so the group floor applies to it.

    Raises:
        DataContractError: a named column is missing from the matrix, or the
            stratum cannot clear ``gwas.min_samples_per_group``.
    """
    columns = list(known_columns)
    missing = [name for name in columns if name not in gi.features]
    if missing:
        raise DataContractError(
            "A known-determinant column is absent from the feature matrix, so "
            "the stratum cannot be defined",
            columns=",".join(missing),
            n_features=len(gi.features),
        )

    stratum = [
        sample_id
        for sample_id in gi.sample_ids
        if all(int(gi.features[name][sample_id]) == 0 for name in columns)
    ]

    if not config.gwas.positive or not config.gwas.negative:
        raise ConfigError(
            "gwas.outcome.positive and gwas.outcome.negative are both needed "
            "to rebuild the stratum's phenotype calls",
            key="gwas.outcome",
        )
    labels = {1: config.gwas.positive[0], 0: config.gwas.negative[0]}
    antibiotic = str(
        (config.raw.get("project") or {}).get("primary_antibiotic")
        or config.antibiotics[0]
    )

    calls = [
        PhenotypeCall(
            sample_id=sample_id,
            antibiotic=antibiotic,
            phenotype=Phenotype(labels[gi.binary_outcome[sample_id]]),
        )
        for sample_id in stratum
    ]
    manifest = SampleManifest(
        [Sample(sample_id=sample_id, assembly_path=None, source="downstream-stratified")
         for sample_id in stratum]
    )
    rows = [
        {
            "sample_id": sample_id,
            **{
                column: gi.features[column][sample_id]
                for column in gi.features
            },
        }
        for sample_id in stratum
    ]
    lineages = {
        sample_id: gi.lineages[sample_id]
        for sample_id in stratum
        if sample_id in gi.lineages
    }

    subset = stage.build_input(config, manifest, calls, rows, lineages=lineages)
    results = stage.ReferenceEngine().run(subset, config)
    return StratifiedScan(
        samples=tuple(stratum),
        known_columns=tuple(columns),
        results=tuple(results),
    )
