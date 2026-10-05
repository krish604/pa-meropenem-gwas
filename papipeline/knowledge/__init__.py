"""Knowledge-table queries.

The mechanism and regulator tables are editable data, so every lookup here is
explicit about what happens when a gene is absent: :class:`UnknownGeneError` is
raised (or a caller-requested ``None`` is returned), never a silent pass-through.
That is what keeps the knowledge base auditable.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from ..config.loader import MechanismSpec, PipelineConfig, RegulatorSpec
from ..errors import UnknownGeneError
from ..logging_utils import get_logger
from ..models import ClaimStatus

LOGGER = get_logger("knowledge")


def genes_for_antibiotic(
    config: PipelineConfig, antibiotic: str
) -> List[MechanismSpec]:
    """All mechanism rows relevant to ``antibiotic``."""
    config.require_antibiotic(antibiotic)
    return [
        spec for spec in config.mechanisms.values() if spec.relevant_to(antibiotic)
    ]


def mechanism_map_for_antibiotic(
    config: PipelineConfig, antibiotic: str
) -> Dict[str, MechanismSpec]:
    """Gene -> MechanismSpec for the given antibiotic."""
    return {spec.gene: spec for spec in genes_for_antibiotic(config, antibiotic)}


def regulator_genes(config: PipelineConfig) -> List[str]:
    """Loci screened by stage 6, in table order."""
    return list(config.regulators.keys())


def classify_determinant(
    config: PipelineConfig, gene: str, antibiotic: str
) -> MechanismSpec:
    """Map a detected gene to its mechanism.

    Raises:
        UnknownGeneError: The gene is not in the knowledge table.
    """
    spec = config.mechanism_for_gene(gene)
    if not spec.relevant_to(antibiotic):
        LOGGER.debug(
            "Gene %s is not mapped to %s (mapped to %s)",
            gene,
            antibiotic,
            spec.antibiotic,
        )
    return spec


def ceiling_for(
    config: PipelineConfig, gene: str, requested: ClaimStatus
) -> ClaimStatus:
    """Clamp a requested claim status to the gene's configured ceiling.

    A raw detection can never be reported as an association or as supported
    purely because the gene is known to be important. Only stage 12/13 may
    raise a claim above the ceiling.
    """
    spec = config.mechanism_for_gene(gene)
    order = {
        ClaimStatus.UNKNOWN: 0,
        ClaimStatus.DETECTED: 1,
        ClaimStatus.PREDICTED: 2,
        ClaimStatus.ASSOCIATED: 3,
        ClaimStatus.SUPPORTED: 4,
    }
    if order[requested] > order[spec.claim_ceiling]:
        LOGGER.debug(
            "Clamping claim for %s from %s to ceiling %s",
            gene,
            requested.value,
            spec.claim_ceiling.value,
        )
        return spec.claim_ceiling
    return requested


def unmapped_genes(
    config: PipelineConfig, genes: Sequence[str]
) -> List[str]:
    """Return the subset of ``genes`` absent from the knowledge table."""
    return [g for g in genes if not config.has_gene(g)]


def require_all_known(config: PipelineConfig, genes: Sequence[str]) -> None:
    """Raise if any gene is unknown. Used at stage boundaries."""
    missing = unmapped_genes(config, genes)
    if missing:
        raise UnknownGeneError(
            "One or more determinants are absent from the knowledge table",
            genes=",".join(missing),
            hint="Add rows to config/mechanisms.tsv",
        )


def mechanism_summary(
    config: PipelineConfig, antibiotic: str
) -> List[Tuple[str, str, str]]:
    """(mechanism, mechanism_class, n_genes) for reporting."""
    counts: Dict[Tuple[str, str], int] = {}
    for spec in genes_for_antibiotic(config, antibiotic):
        key = (spec.mechanism, spec.mechanism_class)
        counts[key] = counts.get(key, 0) + 1
    return sorted((m, mc, str(n)) for (m, mc), n in counts.items())


def regulator_with_promoter(config: PipelineConfig) -> List[RegulatorSpec]:
    """Regulator loci whose promoter region is in scope per the table."""
    return [r for r in config.regulators.values() if r.promoter_screen]


def genes_by_mechanism_class(
    config: PipelineConfig, mechanism_class: str
) -> List[str]:
    """Knowledge-table genes whose ``mechanism_class`` matches.

    The mechanism class lives in ``config/mechanisms.tsv``, so this is the
    single place that joins a regulator locus to its class rather than
    duplicating the mapping.
    """
    return sorted(
        gene
        for gene, spec in config.mechanisms.items()
        if spec.mechanism_class == mechanism_class and gene in config.regulators
    )


def efflux_regulator_genes(config: PipelineConfig) -> List[str]:
    """Regulator loci that act on efflux.

    Used by figure 7 so the figure lists every efflux locus that was
    screened, including those with no variants, rather than only the ones
    that happened to fire.
    """
    return genes_by_mechanism_class(config, "efflux")
