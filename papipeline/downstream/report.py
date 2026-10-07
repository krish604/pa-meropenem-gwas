"""Deliverable 6: the downstream report.

The only place this package puts a claim in front of a reader. Three rules
govern it:

1. **the positive-control gate runs first.** ``build_report`` calls
   :func:`papipeline.downstream.controls.run_control_gate` before it looks at
   anything else, so a baseline that did not recover every configured control
   produces no report at all - not an empty one, not a report with a warning
   in it. The refusal names the control and why.
2. **every claim on a row is a repository claim level.** ``DETECTED``,
   ``PREDICTED``, ``ASSOCIATED``, ``SUPPORTED``, ``UNKNOWN`` - the members of
   :class:`papipeline.models.ClaimStatus` and no others. There is no
   ``CAUSAL`` level to reach for, and none is invented here: the tiers can
   only ever attach ``SUPPORTED``, ``ASSOCIATED`` or ``UNKNOWN``.
3. **what is written is what reads back.** :meth:`DownstreamReport.write`
   goes through the repository's own :func:`papipeline.io.tsv.write_tsv`, so
   the file on disk is readable by :func:`papipeline.io.tsv.read_tsv` with
   :data:`REPORT_COLUMNS` and nothing is lost in translation.

The file this writer produces is the evidence table: one row per feature with
its novelty, tier, claim and the two adjusted p-values behind them. The
interaction and meta-analysis sections travel on the report object rather than
in that file, because they have their own tables - see the registry notes for
why they are not written by a rule yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ..errors import DataContractError
from ..io.tsv import write_tsv
from ..models import ClaimStatus
from .conditional import StratifiedScan
from .controls import ControlOutcome, run_control_gate
from .evidence import EvidenceRecord, build_evidence_table
from .interactions import InteractionResult
from .known_determinants import load_known_determinants
from .lineage_meta import MetaResult

#: The evidence file's columns, in write order. The test that round-trips a
#: written report reads the file back with exactly these.
REPORT_COLUMNS: Tuple[str, ...] = (
    "feature",
    "novelty",
    "evidence_tier",
    "claim_status",
    "basis",
    "conditional_adjusted_p",
    "baseline_adjusted_p",
    "independent_lineages",
    "known_matches",
)


@dataclass(frozen=True)
class DownstreamReport:
    """A gated, graded downstream result - what the run would tell a reader."""

    controls: Tuple[ControlOutcome, ...]
    evidence: Tuple[EvidenceRecord, ...]
    conditional: Tuple[object, ...] = ()
    stratified: Optional[StratifiedScan] = None
    tier1: Tuple[InteractionResult, ...] = ()
    tier2: Tuple[InteractionResult, ...] = ()
    meta: Tuple[MetaResult, ...] = ()

    def rows(self) -> List[dict]:
        """The evidence table as plain rows, one per feature.

        Values are what :func:`papipeline.io.tsv.write_tsv` can write: floats,
        ints, strings, ``None`` for a p-value no scan produced, and the known
        matches as a sequence (``write_tsv`` joins sequences with ``;``).
        """
        return [
            {
                "feature": record.feature,
                "novelty": record.novelty,
                "evidence_tier": record.evidence_tier,
                "claim_status": record.claim_status,
                "basis": record.basis,
                "conditional_adjusted_p": record.conditional_adjusted_p,
                "baseline_adjusted_p": record.baseline_adjusted_p,
                "independent_lineages": record.independent_lineages,
                "known_matches": record.known_matches,
            }
            for record in self.evidence
        ]

    def claim_levels(self) -> set:
        """Every claim level that appears on a row.

        Restricted to the rows, because those are the only place this report
        asserts anything about a feature; the sections carry statistics, not
        claims.
        """
        return {row["claim_status"] for row in self.rows()}

    def write(self, path) -> Path:
        """Write the evidence table to ``path`` and return the path written."""
        return write_tsv(Path(path), self.rows(), REPORT_COLUMNS)


def build_report(
    config,
    *,
    baseline_results: Sequence,
    conditional_results: Sequence = (),
    stratified: Optional[StratifiedScan] = None,
    tier1: Sequence[InteractionResult] = (),
    tier2: Sequence[InteractionResult] = (),
    meta: Sequence[MetaResult] = (),
    known=None,
    evidence: Optional[Sequence[EvidenceRecord]] = None,
) -> DownstreamReport:
    """Gate the baseline scan, grade it, and carry every section into one object.

    Args:
        config: the pipeline configuration.
        baseline_results: stage-12 association rows. Checked against the
            positive controls first, and the check is what decides whether a
            report exists.
        conditional_results: conditional-scan rows, used for the conditional
            column of the evidence table.
        stratified: the known-determinant-free scan, if one ran.
        tier1: pre-specified interaction results.
        tier2: anchored interaction results.
        meta: lineage meta-analysis results.
        known: the novelty filter; loaded from ``config`` when omitted.
        evidence: precomputed evidence records. Supplied records are used as
            they are - the caller has already graded them - but the gate still
            runs first, because a report is gated either way.

    Raises:
        ControlGateError: the baseline did not recover every configured
            control. Raised before any other work.
    """
    baseline = list(baseline_results)
    # First thing, before anything else can be measured or rendered.
    controls = tuple(run_control_gate(config, baseline))

    knowledge = known if known is not None else load_known_determinants(config)
    records = (
        list(evidence)
        if evidence is not None
        else build_evidence_table(
            conditional_results=list(conditional_results),
            baseline_results=baseline,
            known=knowledge,
            config=config,
        )
    )

    # Belt and braces for rule 2 above: the tiers only ever attach a level
    # ClaimStatus defines, and this refuses a record that says otherwise
    # rather than writing a level the repository does not have.
    allowed = {status.value for status in ClaimStatus}
    stray = sorted({record.claim_status for record in records} - allowed)
    if stray:
        raise DataContractError(
            "A downstream claim level is not a member of ClaimStatus, so it "
            "cannot be reported",
            claim=",".join(stray),
            allowed=",".join(sorted(allowed)),
        )

    return DownstreamReport(
        controls=controls,
        evidence=tuple(records),
        conditional=tuple(conditional_results),
        stratified=stratified,
        tier1=tuple(tier1),
        tier2=tuple(tier2),
        meta=tuple(meta),
    )
