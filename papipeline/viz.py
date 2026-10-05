"""Stage 16 (part 1) - figure data preparation.

Every figure is produced from a *prepared data table* built here. The
separation matters: the preparation functions are pure, testable, and
contain no plotting, and the plotting code contains no biology. Nothing in
this module hard-codes a result, a gene name, or a category - all of that
arrives from the pipeline outputs.

The thirteen figures required by the specification map to the builders
below:

===  ==========================================================
  #  builder
===  ==========================================================
  1  :func:`phenotype_distribution_table`
  2  :func:`determinant_distribution_table`
  3  :func:`mechanism_distribution_table`
  4  :func:`gene_by_phenotype_table`
  5  :func:`mechanism_by_phenotype_table`
  6  :func:`oprd_by_phenotype_table`
  7  :func:`efflux_regulator_by_phenotype_table`
  8  :func:`amr_heatmap_matrix`
  9  :func:`mechanism_heatmap_matrix`
 10  :func:`tree_annotation_table`
 11  :func:`gwas_association_table`
 12  :func:`convergence_table`
 13  :func:`cooccurrence_network_table`
===  ==========================================================
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from .logging_utils import get_logger
from .models import (
    AmrDeterminant,
    ClaimStatus,
    ConvergenceCall,
    GwasResult,
    MasterRecord,
    MechanismCall,
    MlstCall,
    PhenotypeCall,
    Phenotype,
    RegulatorVariant,
    StructuralCallStatus,
    StructuralVariant,
)

LOGGER = get_logger("viz")

MULTI_VALUE_SEPARATOR = ","

#: OprD states, in display order. Derived statuses not in this list are
#: appended after the known ones so a new status is visible, not dropped.
OPRD_STATE_ORDER: Tuple[str, ...] = ("intact", "disrupted", "absent", "not_assessed")

#: States a structural locus verdict is allowed to speak in. A verdict outside
#: this set is not trusted to set a status, because `oprd_status_per_sample`
#: feeds the study's central negative and only the four states above are
#: displayable. Mirrors `adapters.oprd_locus.StructuralVerdict` without
#: importing it, so this module keeps no dependency on the adapter.
_STRUCTURAL_STATES: frozenset = frozenset(OPRD_STATE_ORDER)

PHENOTYPE_ORDER: Tuple[str, ...] = ("R", "I", "SDD", "S", "ND")


def _split(value: Optional[str]) -> List[str]:
    if not value:
        return []
    return [v for v in value.split(MULTI_VALUE_SEPARATOR) if v]


def ordered_phenotypes(present: Iterable[str]) -> List[str]:
    """Phenotype categories in a stable, clinically meaningful order."""
    seen = set(present)
    ordered = [p for p in PHENOTYPE_ORDER if p in seen]
    ordered.extend(sorted(seen - set(PHENOTYPE_ORDER)))
    return ordered


# --------------------------------------------------------------------------
# 1-3: distributions
# --------------------------------------------------------------------------


def phenotype_distribution_table(
    phenotype: Mapping[str, PhenotypeCall]
) -> Tuple[List[str], Counter]:
    """Figure 1: imipenem phenotype distribution."""
    counter: Counter = Counter()
    for call in phenotype.values():
        counter[call.phenotype.value] += 1
    return ordered_phenotypes(counter.keys()), counter


def determinant_distribution_table(
    amr: Mapping[str, Sequence[AmrDeterminant]]
) -> Tuple[List[str], Counter]:
    """Figure 2: imipenem AMR determinant distribution."""
    counter: Counter = Counter()
    for determinants in amr.values():
        seen: Set[str] = set()
        for determinant in determinants:
            name = determinant.gene or determinant.determinant
            if name and name not in seen:
                seen.add(name)
                counter[name] += 1
    return sorted(counter), counter


def mechanism_distribution_table(
    mechanisms: Mapping[str, Sequence[MechanismCall]]
) -> Tuple[List[str], Counter]:
    """Figure 3: imipenem resistance mechanism distribution."""
    counter: Counter = Counter()
    for calls in mechanisms.values():
        for mechanism in {c.mechanism for c in calls if c.mechanism}:
            counter[mechanism] += 1
    return sorted(counter), counter


# --------------------------------------------------------------------------
# 4-7: cross-tabulations
# --------------------------------------------------------------------------


def crosstab(
    feature_sets: Mapping[str, Set[str]],
    feature_name: str,
    phenotypes: Mapping[str, Optional[str]],
) -> Tuple[List[str], List[str], Dict[Tuple[str, str], int]]:
    """Build a feature x phenotype contingency table.

    Args:
        feature_sets: sample -> set of features carrying that sample.
        feature_name: The feature to tabulate.
        phenotypes: sample -> phenotype category (``None`` when missing).

    Returns:
        ``(phenotype_order, [feature_name], {(phenotype, feature): n})``.

    The category axis spans every phenotype observed in the cohort, not
    only those where this feature was seen. A category in which the feature
    is absent must still appear with a zero count, otherwise the figure
    would silently omit a phenotype. Samples with no phenotype appear under
    no category at all; they are not folded into a neighbouring bin.
    """
    counter: Counter = Counter()
    for sample_id, category in phenotypes.items():
        if category is None:
            continue
        if feature_name in feature_sets.get(sample_id, set()):
            counter[(category, feature_name)] += 1
    order = ordered_phenotypes({c for c in phenotypes.values() if c is not None})
    return order, [feature_name], dict(counter)


def _category_totals(phenotypes: Mapping[str, Optional[str]]) -> Counter:
    """How many samples sit in each phenotype category, cohort-wide."""
    return Counter(c for c in phenotypes.values() if c is not None)


def _invert(
    per_item: Mapping[str, Sequence[str]]
) -> Dict[str, Set[str]]:
    """Transpose a ``feature -> samples`` mapping into ``sample -> features``.

    :func:`crosstab` indexes by sample, so the per-gene carrier sets have to
    be transposed before being passed to it.
    """
    out: Dict[str, Set[str]] = {}
    for feature, samples in per_item.items():
        for sample_id in samples:
            out.setdefault(sample_id, set()).add(feature)
    return out


def gene_by_phenotype_table(
    amr: Mapping[str, Sequence[AmrDeterminant]],
    phenotype: Mapping[str, PhenotypeCall],
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Figure 4: AMR gene x imipenem phenotype, per gene."""
    gene_samples: Dict[str, Set[str]] = {}
    for sample_id, determinants in amr.items():
        for determinant in determinants:
            name = determinant.gene or determinant.determinant
            if name:
                gene_samples.setdefault(name, set()).add(sample_id)

    sample_to_genes = _invert(gene_samples)
    phenotypes = {sid: call.phenotype.value for sid, call in phenotype.items()}
    totals = _category_totals(phenotypes)
    out: Dict[str, Dict[str, Dict[str, int]]] = {}
    for gene in sorted(gene_samples):
        order, _features, counts = crosstab(sample_to_genes, gene, phenotypes)
        out[gene] = {
            category: {
                "carriers": counts.get((category, gene), 0),
                "total": totals.get(category, 0),
            }
            for category in order
        }
    return out


def mechanism_by_phenotype_table(
    mechanisms: Mapping[str, Sequence[MechanismCall]],
    phenotype: Mapping[str, PhenotypeCall],
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Figure 5: mechanism x imipenem phenotype, per mechanism."""
    mech_samples: Dict[str, Set[str]] = {}
    for sample_id, calls in mechanisms.items():
        for call in calls:
            if call.mechanism:
                mech_samples.setdefault(call.mechanism, set()).add(sample_id)

    sample_to_mechs = _invert(mech_samples)
    phenotypes = {sid: call.phenotype.value for sid, call in phenotype.items()}
    totals = _category_totals(phenotypes)
    out: Dict[str, Dict[str, Dict[str, int]]] = {}
    for mechanism in sorted(mech_samples):
        order, _features, counts = crosstab(sample_to_mechs, mechanism, phenotypes)
        out[mechanism] = {
            category: {
                "carriers": counts.get((category, mechanism), 0),
                "total": totals.get(category, 0),
            }
            for category in order
        }
    return out


def oprd_status_per_sample(
    variants: Mapping[str, Sequence[RegulatorVariant]],
    annotated_genes: Optional[Mapping[str, Sequence[str]]] = None,
    locus_resolution: Optional[Mapping[str, Any]] = None,
) -> Dict[str, str]:
    """Derive OprD status per sample.

    ``absent`` and ``disrupted`` are the two states this study's central negative
    rests on, so each has exactly one legitimate source and no other:

    * ``disrupted`` comes only from a ``disrupted`` structural verdict - which
      the classifier issues for a frameshift or a premature stop, decided by
      lesion type and never by ORF length - or from the variant screen's own
      ``GENE_DISRUPTION`` call, which is an independent evidence source with its
      own contract.
    * ``absent`` comes only from an ``absent`` structural verdict, and
      :class:`~papipeline.adapters.oprd_locus.StructuralVerdict` only permits
      that from a *recorded* search that returned no hit. Not locating a locus
      is not evidence of deleting it.

    ``intact`` is only assigned when the locus is confirmed present, because the
    absence of a variant call is not by itself evidence of an intact gene.

    **The annotation-symbol test can no longer produce ``absent``.** It used to,
    and on real data that was a fabricated negative: ``"oprD" in genes`` asks
    whether Bakta wrote the symbol, and on real data Bakta writes it on the
    OprD/OprP/OprQ paralog family - across ten isolates, 34 of the 36
    ``gene=oprD`` CDS features are paralogs at 33.7-36.8% identity to PAO1
    PA0958, and 8 of the 10 isolates carry the true locus as an unlabelled CDS.
    It reported ``absent`` for eight isolates that all have the gene. It now
    yields ``intact`` when the symbol is present and ``not_assessed`` otherwise,
    which is the most that a string test on an annotation symbol can support.

    A resolution from :mod:`papipeline.adapters.oprd_locus` is keyed on protein
    or DNA identity against the pinned PAO1 oprD sequence instead. It may be a
    :class:`~papipeline.adapters.oprd_locus.LocusResolution` (blastp against the
    isolate proteome) or a
    :class:`~papipeline.adapters.oprd_locus.StructuralCall` (tblastn against the
    assembly, which also judges the reading frame). A structural verdict is
    authoritative where it speaks; a refusal is ``not_assessed``, never
    ``absent``.

    ``locus_resolution`` still defaults to ``None`` so that callers which do not
    run a locus search keep working unchanged. That default is retained
    deliberately: it changes what those callers *report* (see below), not
    whether they run. **It does change figures.** With no resolution supplied,
    every sample that previously came out ``absent`` from the symbol test now
    comes out ``not_assessed``. That is the intended correction - the old value
    was the fabricated negative - but it is a visible change and callers that
    depended on the old output must be re-checked, not assumed compatible.
    """
    status: Dict[str, str] = {}
    for sample_id, calls in variants.items():
        oprd = [c for c in calls if c.gene == "oprD"]
        if not oprd:
            continue
        types = {c.variant_type for c in oprd}
        if "GENE_ABSENCE" in types:
            status[sample_id] = "absent"
        elif "GENE_DISRUPTION" in types:
            status[sample_id] = "disrupted"
        else:
            status[sample_id] = "variant"
    for sample_id, resolution in (locus_resolution or {}).items():
        if sample_id in status:
            continue
        verdict = getattr(resolution, "verdict", None)
        # A structural verdict speaks in the same vocabulary as this function
        # and is authoritative where it is decisive.
        if verdict in _STRUCTURAL_STATES:
            status[sample_id] = verdict
            continue
        status[sample_id] = (
            "intact" if getattr(resolution, "is_resolved", False)
            else "not_assessed"
        )
    for sample_id, genes in (annotated_genes or {}).items():
        if sample_id in status:
            continue
        if locus_resolution is not None and sample_id in locus_resolution:
            continue
        # A symbol on an annotation row is evidence the gene was called, and
        # nothing more. It is not evidence the gene is whole, and its absence is
        # not evidence the gene is gone.
        status[sample_id] = "intact" if "oprD" in set(genes) else "not_assessed"
    return status


def oprd_by_phenotype_table(
    oprd_status: Mapping[str, str],
    phenotype: Mapping[str, PhenotypeCall],
) -> Tuple[List[str], List[str], Dict[Tuple[str, str], int]]:
    """Figure 6: OprD status x imipenem phenotype."""
    counts: Counter = Counter()
    for sample_id, call in phenotype.items():
        state = oprd_status.get(sample_id, "not_assessed")
        counts[(call.phenotype.value, state)] += 1

    states = [s for s in OPRD_STATE_ORDER if s in {k[1] for k in counts}]
    states.extend(sorted({k[1] for k in counts} - set(OPRD_STATE_ORDER)))
    phenotypes = ordered_phenotypes({k[0] for k in counts})
    return phenotypes, states, dict(counts)


def efflux_regulator_by_phenotype_table(
    variants: Mapping[str, Sequence[RegulatorVariant]],
    efflux_regulator_genes: Sequence[str],
    phenotype: Mapping[str, PhenotypeCall],
) -> Tuple[List[str], List[str], Dict[Tuple[str, str], int]]:
    """Figure 7: efflux-regulator alterations x imipenem phenotype.

    A sample is counted as carrying an alteration if it has a variant in any
    of the configured efflux regulator loci. Genes with no variants at all
    are still listed, with zero counts, so the figure shows the loci that
    were screened rather than only the ones that fired.
    """
    wanted = set(efflux_regulator_genes)
    counts: Counter = Counter()
    for sample_id, call in phenotype.items():
        altered = {
            c.gene for c in variants.get(sample_id, ()) if c.gene in wanted
        }
        for gene in sorted(altered):
            counts[(call.phenotype.value, gene)] += 1

    genes = list(efflux_regulator_genes)
    phenotypes = ordered_phenotypes({k[0] for k in counts})
    return phenotypes, genes, dict(counts)


# --------------------------------------------------------------------------
# 8-9: heatmaps
# --------------------------------------------------------------------------


def presence_matrix(
    feature_sets: Mapping[str, Set[str]],
    sample_ids: Sequence[str],
    features: Optional[Sequence[str]] = None,
) -> Tuple[List[str], List[str], List[List[int]]]:
    """Build a sample x feature binary matrix.

    Rows are features (so the matrix can be fed to a heatmap directly with
    one row per feature) and columns are samples.
    """
    names = list(features) if features is not None else sorted(
        {name for members in feature_sets.values() for name in members}
    )
    matrix = [
        [1 if name in feature_sets.get(sample_id, set()) else 0 for sample_id in sample_ids]
        for name in names
    ]
    return names, list(sample_ids), matrix


def amr_heatmap_matrix(
    amr: Mapping[str, Sequence[AmrDeterminant]], sample_ids: Sequence[str]
) -> Tuple[List[str], List[str], List[List[int]]]:
    """Figure 8: AMR gene heatmap matrix."""
    feature_sets: Dict[str, Set[str]] = {}
    for sample_id, determinants in amr.items():
        feature_sets[sample_id] = {
            d.gene or d.determinant for d in determinants if (d.gene or d.determinant)
        }
    return presence_matrix(feature_sets, sample_ids)


def mechanism_heatmap_matrix(
    mechanisms: Mapping[str, Sequence[MechanismCall]], sample_ids: Sequence[str]
) -> Tuple[List[str], List[str], List[List[int]]]:
    """Figure 9: mechanism heatmap matrix."""
    feature_sets: Dict[str, Set[str]] = {}
    for sample_id, calls in mechanisms.items():
        feature_sets[sample_id] = {c.mechanism for c in calls if c.mechanism}
    return presence_matrix(feature_sets, sample_ids)


# --------------------------------------------------------------------------
# 10-13: tree, GWAS, convergence, network
# --------------------------------------------------------------------------


def tree_annotation_table(
    sample_ids: Sequence[str],
    amr: Mapping[str, Sequence[AmrDeterminant]],
    mechanisms: Mapping[str, Sequence[MechanismCall]],
    phenotype: Mapping[str, PhenotypeCall],
    lineages: Mapping[str, str],
    mlst: Mapping[str, Optional[MlstCall]],
) -> List[Dict[str, object]]:
    """Figure 10: per-tip annotation for the phylogeny figure.

    The tip order is the manifest order, which stage 10 has already
    validated against the tree.
    """
    rows: List[Dict[str, object]] = []
    for sample_id in sample_ids:
        call = phenotype.get(sample_id)
        mlst_call = mlst.get(sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "lineage": lineages.get(sample_id),
                "phenotype": call.phenotype.value if call else None,
                "mlst": mlst_call.sequence_type if mlst_call else None,
                "amr_genes": sorted(
                    {d.gene or d.determinant for d in amr.get(sample_id, ())}
                    - {None}
                ),
                "mechanisms": sorted({m.mechanism for m in mechanisms.get(sample_id, ())}),
            }
        )
    return rows


def gwas_association_table(
    results: Sequence[GwasResult], config=None, top: int = 20
) -> List[Dict[str, object]]:
    """Figure 11: GWAS association plot data, strongest association first.

    Rows are ordered by adjusted p-value. Lineage-confounded features are
    annotated so a reader can see why a feature is not a resistance
    association.
    """
    threshold = getattr(getattr(config, "gwas", None), "significance_threshold", 0.05)
    confound_threshold = getattr(
        getattr(config, "gwas", None), "lineage_confound_threshold", 0.9
    )

    ordered = sorted(
        results,
        key=lambda r: (r.adjusted_p_value if r.adjusted_p_value is not None else 1.0),
    )[:top]

    rows: List[Dict[str, object]] = []
    for result in ordered:
        total = sum(result.lineage_distribution.values())
        top_share = (
            max(result.lineage_distribution.values()) / total if total else 0.0
        )
        rows.append(
            {
                "feature": result.feature,
                "feature_type": result.feature_type,
                "p_value": result.p_value,
                "adjusted_p_value": result.adjusted_p_value,
                "effect": result.effect,
                "frequency": result.frequency,
                "dominant_lineage_share": round(top_share, 4),
                "lineage_linked": top_share >= confound_threshold,
                "passes_threshold": (
                    result.adjusted_p_value is not None
                    and result.adjusted_p_value <= threshold
                ),
                "model": result.model,
            }
        )
    return rows


def convergence_table(
    calls: Sequence[ConvergenceCall],
) -> Tuple[List[str], Counter]:
    """Figure 12: convergence category distribution."""
    counter: Counter = Counter(call.convergence_category.value for call in calls)
    return sorted(counter), counter


def cooccurrence_network_table(
    pairs: Sequence, threshold: float = 0.05
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Figure 13: co-occurrence network as nodes and edges.

    Only pairs passing the adjusted-p threshold become edges. A gene-gene
    and a gene-mechanism pair live in the same edge list but carry their
    ``feature_type`` so the figure can style them differently; nothing about
    the network implies direction or causation.
    """
    edges = [
        {
            "source": p.feature_a,
            "target": p.feature_b,
            "feature_type": p.feature_type,
            "weight": p.statistic_value,
            "adjusted_p_value": p.adjusted_p_value,
            "n_both": p.n_both,
            "interpretation_limit": p.interpretation_limit,
        }
        for p in pairs
        if p.adjusted_p_value is not None and p.adjusted_p_value <= threshold
    ]
    degree: Counter = Counter()
    for edge in edges:
        degree[edge["source"]] += 1
        degree[edge["target"]] += 1
    nodes = [
        {"name": name, "degree": count, "feature_type": _node_type(name, edges)}
        for name, count in sorted(degree.items())
    ]
    return nodes, edges


def _node_type(name: str, edges: Sequence[Mapping[str, object]]) -> str:
    types = {
        str(edge["feature_type"])
        for edge in edges
        if edge["source"] == name or edge["target"] == name
    }
    return "/".join(sorted(types)) if types else "unknown"
