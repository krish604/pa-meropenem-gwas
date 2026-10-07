"""The one entry point that runs all seven downstream steps, in order.

The steps themselves are owned by the submodules in this package; this module
is the *order*, the arguments each step is handed, and the single file the run
writes. It exists so the order the registry notes prescribe
(``.build/downstream-stats.registry-notes.md`` section 1) is stated in code
once, and so ``papipeline/run.py`` can call one function instead of seven.

Order, and why it is not free to choose:

1. **control gate** - :func:`papipeline.downstream.controls.run_control_gate`.
   Before anything downstream reports a finding, the baseline association table
   has to show the scan still recovers what this project already knows. It
   raises :class:`ControlGateError` naming every control it failed to recover,
   so nothing below it runs.
2. **conditional scan** - known determinants held fixed.
3. **stratified scan** - the reference association inside the isolates lacking
   every known determinant. It rebuilds its own subset through
   ``stages.gwas.build_input``, so ``gwas.min_samples_per_group`` applies to the
   stratum rather than to the cohort it was carved out of.
4. **interactions** - tier 1 (pre-specified pairs from
   ``config/interaction_pairs.tsv``) then tier 2 (every known column against
   every novel column, capped by ``downstream.tier2_max_pairs``).
5. **lineage meta** - fixed-effect pooling per feature across lineages; needs
   stage-9 lineages on the ``GwasInput`` and refuses with ``lineage.method``
   when there are fewer than two.
6. **evidence** - novelty and tiers A-D over both scans, with stage 13's
   convergence call as the authority for the independent-lineage count when
   the determinant is one it classified.
7. **report** - :func:`papipeline.downstream.report.build_report`, which runs
   the gate *again* on its own baseline and only then renders. The gate runs
   twice by construction: the first refusal stops the work, the second stops a
   report built from an argument someone else supplied.

What is written: the evidence table, and nothing else. ``DownstreamReport``
carries the tier-1, tier-2 and meta sections as objects; they are their own
tables, no column order has been fixed for them, and inventing one here would
be a format decision made in a wiring commit. See
``docs/data_contract.md`` for the evidence table.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..logging_utils import get_logger
from ..stages.gwas_features import target_antibiotic
from .conditional import conditional_scan, covariate_columns, stratified_scan
from .controls import ControlOutcome, run_control_gate
from .evidence import build_evidence_table
from .interactions import (
    interaction_pairs_path,
    load_interaction_pairs,
    resolve_tier2_max_pairs,
    tier1_scan,
    tier2_anchored_scan,
)
from .lineage_meta import lineage_meta_analysis
from .report import DownstreamReport, build_report
from .known_determinants import load_known_determinants

LOGGER = get_logger("downstream")

#: The one file a run writes. Named for its content rather than for the step
#: that produced it, because steps 1-7 produce one table between them.
EVIDENCE_FILENAME = "downstream_evidence.tsv"


@dataclass(frozen=True)
class DownstreamRun:
    """What the seven steps produced, for the caller to record."""

    controls: Tuple[ControlOutcome, ...]
    report: DownstreamReport
    evidence_path: Path


def run_downstream_analyses(
    config: PipelineConfig,
    *,
    gwas_results: Sequence[Any],
    gwas_input: Any,
    convergence: Optional[Mapping[str, Any]] = None,
    output_dir: Path,
) -> DownstreamRun:
    """Run the seven steps over one cohort's stage-12 results and write the table.

    Args:
        config: Loaded pipeline configuration. Every science key this reads
            comes from the raw mapping (``downstream.*``,
            ``gwas.significance_threshold``, ``convergence.min_independent_lineages``)
            the way the rest of the package reads them, so a caller that
            ``dataclasses.replace``\ s ``raw`` sees the override.
        gwas_results: stage 12's association rows, in memory. Anything with
            ``feature`` and ``adjusted_p_value`` - the gate and the evidence
            table are the only steps that read it directly.
        gwas_input: the stage-12 :class:`~papipeline.stages.gwas.GwasInput`,
            carrying the feature matrix, the outcome and stage-9's lineages.
            Steps 2-5 are defined on it.
        convergence: stage 13's calls keyed by determinant, optional. Used as
            the authority for a feature's independent-lineage count when there
            is a call for it; a key that never matches simply falls back to
            the lineage distribution, which is the documented behaviour of
            :func:`papipeline.downstream.evidence.independent_lineages`.
        output_dir: where the evidence table is written. The caller owns the
            choice - the package does not pick a results directory.

    Returns:
        The control verdicts, the report object (which carries every section),
        and the path written.

    Raises:
        ControlGateError: the baseline did not recover every configured
            positive control. Raised before any other step, so no output file
            is written when it fails.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1, first. Everything below reports a finding, and a finding the
    # scan cannot be trusted to produce is worse than no finding.
    controls = tuple(run_control_gate(config, gwas_results))
    LOGGER.info(
        "downstream: positive-control gate passed (%s)",
        ", ".join(f"{c.control}=recovered" for c in controls),
    )

    known = load_known_determinants(config)
    antibiotic = target_antibiotic(config)

    # Step 2 - conditional, known determinants held fixed.
    covariates = covariate_columns(gwas_input, known, config)
    conditional = conditional_scan(gwas_input, covariates)

    # Step 3 - stratified. Only the known columns that are actually in THIS
    # matrix: naming one that is not would refuse (a column absent from the
    # matrix cannot define a stratum), and the fixture cohorts carry a subset
    # of the knowledge table, as a real cohort does.
    known_columns = [f for f in gwas_input.features if known.is_known(f)]
    stratified = stratified_scan(gwas_input, known_columns, config)

    # Step 4 - both interaction tiers, over the project's own antibiotic so
    # the tier-1 table's relevance filter and the run agree on the drug.
    pairs = load_interaction_pairs(interaction_pairs_path(config), known=known)
    tier1 = tier1_scan(gwas_input, pairs, antibiotic=antibiotic)
    tier2 = tier2_anchored_scan(
        gwas_input,
        known,
        max_pairs=resolve_tier2_max_pairs(config),
        antibiotic=antibiotic,
    )

    # Step 5 - lineage meta, pooled across stage-9's labels.
    meta = lineage_meta_analysis(gwas_input, config)

    # Step 6 - grade the union of both scans.
    evidence = build_evidence_table(
        conditional_results=conditional,
        baseline_results=gwas_results,
        known=known,
        config=config,
        convergence=convergence,
    )

    # Step 7 - the report, which gates again on the baseline it is given.
    report = build_report(
        config,
        baseline_results=gwas_results,
        conditional_results=conditional,
        stratified=stratified,
        tier1=tier1,
        tier2=tier2,
        meta=meta,
        known=known,
        evidence=evidence,
    )
    evidence_path = report.write(output_dir / EVIDENCE_FILENAME)

    LOGGER.info(
        "downstream: evidence=%d features | conditional=%d | stratified=%d "
        "samples in %d known-column(s) | tier1=%d pairs | tier2=%d pairs | "
        "meta=%d features -> %s",
        len(report.evidence),
        len(report.conditional),
        len(stratified.samples),
        len(stratified.known_columns),
        len(report.tier1),
        len(report.tier2),
        len(report.meta),
        evidence_path,
    )
    return DownstreamRun(
        controls=controls, report=report, evidence_path=evidence_path
    )
