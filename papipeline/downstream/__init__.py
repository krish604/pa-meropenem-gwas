"""Post-scan statistics: the control gate, the conditional and stratified
scans, the interaction tiers, the lineage meta-analysis, the novelty filter and
the evidence tiers that feed the report.

What each submodule is for:

===============================  =========================================
:mod:`~papipeline.downstream.known_determinants`
                                 what "already known" means, derived from
                                 ``config/mechanisms.tsv`` into
                                 ``config/known_determinants.tsv``
:mod:`~papipeline.downstream.controls`
                                 deliverable 1: the positive-control gate
:mod:`~papipeline.downstream.conditional`
                                 deliverable 2: conditional and stratified
                                 scans
:mod:`~papipeline.downstream.interactions`
                                 deliverable 3: interaction tiers 1 and 2
                                 (tier 3 is skipped by design)
:mod:`~papipeline.downstream.lineage_meta`
                                 deliverable 4: fixed-effect meta-analysis of
                                 per-lineage estimates
:mod:`~papipeline.downstream.evidence`
                                 deliverable 5: novelty filter and tiers A-D
:mod:`~papipeline.downstream.report`
                                 deliverable 6: the table a reader sees
:mod:`~papipeline.downstream.logit`
                                 the shared binary-outcome fits
===============================  =========================================

Two rules the whole package obeys. Claim levels come from
:class:`papipeline.models.ClaimStatus` and there is no ``CAUSAL`` member to
reach for. And every fitted model is stamped ``:NO_KINSHIP_CORRECTION``, the
same stamp :class:`papipeline.stages.gwas.ReferenceEngine` uses, because
nothing here corrects for population structure either.

None of this is wired into a stage yet; see
``.build/downstream-stats.registry-notes.md``.
"""

from __future__ import annotations

__all__ = [
    "conditional",
    "controls",
    "evidence",
    "interactions",
    "known_determinants",
    "lineage_meta",
    "logit",
    "report",
]
