"""Stage 15 - integration.

Builds the master table: one row per sample, joining genotype, mechanism and
phenotype into a single record.

The central design decision is what goes in the many-to-one columns. A sample
can carry several genes, several variants and several mechanisms, so those
columns hold a delimited list rather than a single "winner". Choosing one
determinant per sample and calling it *the* determinant would be an
interpretive act, and an indefensible one.

The ``confidence`` column is derived, not inherited:

* ``SUPPORTED`` - the feature is GWAS-significant *and* convergent across
  independent lineages;
* ``ASSOCIATED`` - GWAS-significant, or convergent, but not both;
* ``PREDICTED`` - inferred without direct observation, e.g. a candidate SV;
* ``DETECTED`` - directly observed, with no statistical support claimed;
* ``UNKNOWN`` - nothing to report for this sample.

``confidence`` never exceeds the evidence available, and it is never
"causal" because that value does not exist in the vocabulary.

**A REAL run refuses a manifest sample with no lineage label.** Not because a
lineage is needed to fill a column - ``lineage`` is one column of seventeen -
but because every lineage-aware claim above this table is computed from it:
convergence counts *independent lineages*, co-occurrence stratifies by them, and
the GWAS lineage-confound check needs more than one. A master table that
publishes ``lineage=None`` for a sample and ``lineage="unknown"`` for another
presents those two as the same fact, and the three downstream readers cannot
tell which is which - one is "we do not know", the other is "we know and it is
nothing". Under the old producer EVERY sample read as ``"unknown"``, so this
was not an edge case; it was the whole table.

The check names the offending samples rather than their count, because the count
tells a reader nothing they can act on and the list is what they can.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..errors import StageError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import (
    AmrDeterminant,
    ClaimStatus,
    ConvergenceCall,
    GwasResult,
    MASTER_TABLE_COLUMNS,
    MasterRecord,
    MechanismCall,
    MlstCall,
    PhenotypeCall,
    RegulatorVariant,
    RunMode,
    StructuralCallStatus,
    StructuralVariant,
    VirulenceFactor,
    strongest_claim,
)

LOGGER = get_logger("stages.integration")

MULTI_VALUE_SEPARATOR = ","

#: The sentinel `stages.phylogeny` writes when stage 3 had no sequence type.
#: Imported rather than restated: the refusing stages and the producer must
#: match on one value, and two copies of the literal `"unknown"` is how they
#: would stop matching on each other.
from .phylogeny import UNKNOWN_LINEAGE  # noqa: E402  (avoids an import cycle)


def require_complete_lineages(
    manifest: SampleManifest, lineages: Mapping[str, str]
) -> None:
    """Refuse a manifest sample whose `lineage_label` is absent or the sentinel.

    "Absent" and `"unknown"` are refused together and reported together, because
    to a reader of the master table they are the same defect: one says the
    crosswalk never reached this sample, the other says it reached it and found
    no sequence type. Neither is a lineage, and publishing either while the
    other exists makes the column unreadable.

    Raises:
        StageError: naming every offending sample.
    """
    missing = [
        sample_id
        for sample_id in manifest.sample_ids
        if not str(lineages.get(sample_id) or "").strip()
        or str(lineages.get(sample_id)).strip() == UNKNOWN_LINEAGE
    ]
    if not missing:
        return
    raise StageError(
        f"{len(missing)} of {len(manifest)} manifest sample(s) have no "
        f"lineage_label, so this master table cannot be built: "
        f"{', '.join(missing[:15])}"
        + (" ..." if len(missing) > 15 else "")
        + ". `lineage` is not a decorative column - convergence counts "
        "INDEPENDENT LINEAGES from it, co-occurrence stratifies by it, and the "
        "GWAS lineage-confound check needs more than one distinct value, so a "
        "table with a blank or sentinel entry reports those findings against a "
        "grouping that was never established. Under ruling R4 a lineage_label "
        "IS the MLST sequence type from stage 3, so an unlabelled sample means "
        "stage 3 produced no sequence type for it: check "
        "`intermediate/mlst/mlst_results.tsv` for that sample's row and its "
        "`MLST_status`. These samples are named, not counted, so they can be "
        "looked up.",
        missing=",".join(missing[:15]),
        n_missing=len(missing),
        n_samples=len(manifest),
        sentinel=UNKNOWN_LINEAGE,
    )


def _join(values: Iterable[str]) -> Optional[str]:
    """Join a set of values, or return ``None`` when empty.

    ``None`` is used rather than an empty string so that "nothing detected"
    is distinguishable from "detected as blank" downstream.
    """
    ordered = sorted({v for v in values if v})
    return MULTI_VALUE_SEPARATOR.join(ordered) if ordered else None


def sample_confidence(
    config: PipelineConfig,
    sample_id: str,
    mechanisms: Sequence[MechanismCall],
    gwas_results: Sequence[GwasResult],
    convergence: Sequence[ConvergenceCall],
    structural: Sequence[StructuralVariant],
) -> ClaimStatus:
    """Derive the per-sample confidence from the available evidence."""
    threshold = config.gwas.significance_threshold
    significant_features = {
        r.feature for r in gwas_results if r.adjusted_p_value is not None and r.adjusted_p_value <= threshold
    }
    convergent = {
        c.determinant
        for c in convergence
        if c.convergence_category.value == "recurrent_convergent"
    }

    sample_mechanism_determinants = {m.determinant for m in mechanisms}

    has_supported = bool(sample_mechanism_determinants & significant_features & convergent)
    if has_supported:
        return ClaimStatus.SUPPORTED

    has_associated = bool(sample_mechanism_determinants & (significant_features | convergent))
    if has_associated:
        return ClaimStatus.ASSOCIATED

    statuses: List[ClaimStatus] = [m.evidence_level for m in mechanisms]
    for sv in structural:
        if sv.call_status is StructuralCallStatus.CANDIDATE:
            statuses.append(ClaimStatus.PREDICTED)
        elif sv.call_status is StructuralCallStatus.CONFIRMED:
            statuses.append(ClaimStatus.DETECTED)
    return strongest_claim(statuses)


def build_master_table(
    config: PipelineConfig,
    manifest: SampleManifest,
    antibiotic: str,
    phenotype: Mapping[str, PhenotypeCall],
    amr: Mapping[str, Sequence[AmrDeterminant]],
    mechanisms: Mapping[str, Sequence[MechanismCall]],
    regulator_variants: Mapping[str, Sequence[RegulatorVariant]],
    structural: Mapping[str, Sequence[StructuralVariant]],
    mlst: Mapping[str, Optional[MlstCall]],
    lineages: Mapping[str, str],
    virulence: Mapping[str, Sequence[VirulenceFactor]],
    gwas_results: Sequence[GwasResult],
    convergence: Sequence[ConvergenceCall],
    mode: Optional[RunMode] = None,
) -> List[MasterRecord]:
    """Assemble the master table, one record per manifest sample.

    `mode` gates the lineage completeness check. It is optional and defaults to
    no check, because every existing caller reaches this function without one
    and a caller that has not said which mode it is in has not claimed to be a
    production run. `run` passes it.
    """
    if mode is RunMode.REAL:
        require_complete_lineages(manifest, lineages)

    records: List[MasterRecord] = []

    for sample_id in manifest.sample_ids:
        pheno = phenotype.get(sample_id)
        sample_amr = amr.get(sample_id, ())
        sample_mech = mechanisms.get(sample_id, ())
        sample_variants = regulator_variants.get(sample_id, ())
        sample_sv = structural.get(sample_id, ())
        sample_vf = virulence.get(sample_id, ())

        mlst_call = mlst.get(sample_id)
        lineage = lineages.get(sample_id)

        regulator_names = {
            v.gene for v in sample_variants if config.regulator(v.gene) is not None
        }

        confirmed_sv = [
            sv.variant_id
            for sv in sample_sv
            if sv.call_status is StructuralCallStatus.CONFIRMED
        ]
        candidate_sv = [
            sv.variant_id
            for sv in sample_sv
            if sv.call_status is StructuralCallStatus.CANDIDATE
        ]
        if candidate_sv and not confirmed_sv:
            # A candidate is reported, but labelled so it cannot be mistaken
            # for a confirmed call.
            confirmed_sv = [f"candidate:{vid}" for vid in candidate_sv]

        evidence_notes = _build_notes(
            sample_amr, sample_variants, sample_sv, pheno is not None
        )

        records.append(
            MasterRecord(
                sample_id=sample_id,
                antibiotic=antibiotic,
                phenotype=pheno.phenotype.value if pheno else None,
                amr_gene=_join(d.gene or d.determinant for d in sample_amr),
                amr_variant=_join(
                    d.determinant for d in sample_amr if d.variant
                ),
                chromosomal_mutation=_join(v.variant for v in sample_variants),
                regulator=_join(regulator_names),
                mechanism=_join(m.mechanism for m in sample_mech),
                structural_variant=_join(confirmed_sv),
                mlst=mlst_call.sequence_type if mlst_call else None,
                lineage=lineage,
                virulence_profile=_join(v.virulence_factor for v in sample_vf),
                gwas_feature=None,
                gwas_status=None,
                convergence_status=None,
                confidence=sample_confidence(
                    config,
                    sample_id,
                    sample_mech,
                    gwas_results,
                    convergence,
                    sample_sv,
                ),
                evidence_notes=evidence_notes,
            )
        )

    return records


def _build_notes(
    amr: Sequence[AmrDeterminant],
    variants: Sequence[RegulatorVariant],
    structural: Sequence[StructuralVariant],
    has_phenotype: bool,
) -> Optional[str]:
    """Short, factual provenance string per sample.

    States what was observed and, equally importantly, what was *not*
    concluded.
    """
    parts: List[str] = []
    if amr:
        parts.append(f"determinant_detected:{len(amr)}")
    if variants:
        parts.append(f"variant_detected:{len(variants)}")
    confirmed = sum(
        1 for sv in structural if sv.call_status is StructuralCallStatus.CONFIRMED
    )
    candidates = sum(
        1 for sv in structural if sv.call_status is StructuralCallStatus.CANDIDATE
    )
    if confirmed:
        parts.append(f"sv_confirmed:{confirmed}")
    if candidates:
        parts.append(f"sv_candidate_not_promoted:{candidates}")
    if not has_phenotype:
        parts.append("phenotype_missing")
    parts.append("detection_is_not_resistance")
    return ";".join(parts)


def annotate_gwas(
    records: Sequence[MasterRecord],
    gwas_results: Sequence[GwasResult],
    config: PipelineConfig,
) -> List[MasterRecord]:
    """Attach the strongest GWAS feature and its status per sample.

    Rewrites the records (they are frozen dataclasses, so new instances are
    created via ``dataclasses.replace``).
    """
    from dataclasses import replace

    threshold = config.gwas.significance_threshold
    significant = {
        r.feature: r
        for r in gwas_results
        if r.adjusted_p_value is not None and r.adjusted_p_value <= threshold
    }
    if not significant:
        return list(records)

    out: List[MasterRecord] = []
    for record in records:
        genes = set((record.amr_gene or "").split(MULTI_VALUE_SEPARATOR)) - {""}
        matches = sorted(genes & set(significant))
        if not matches:
            out.append(
                replace(record, gwas_status="no_significant_associated_feature")
            )
            continue
        result = significant[matches[0]]
        out.append(
            replace(
                record,
                gwas_feature=matches[0],
                gwas_status=(
                    f"associated:adj.p={result.adjusted_p_value:.4g}"
                    if result.adjusted_p_value is not None
                    else "associated"
                ),
            )
        )
    return out


def annotate_convergence(
    records: Sequence[MasterRecord], convergence: Sequence[ConvergenceCall]
) -> List[MasterRecord]:
    """Attach the convergence category of a sample's convergent determinant."""
    from dataclasses import replace

    by_determinant = {c.determinant: c for c in convergence}
    out: List[MasterRecord] = []
    for record in records:
        genes = set((record.amr_gene or "").split(MULTI_VALUE_SEPARATOR)) - {""}
        categories = {
            by_determinant[g].convergence_category.value
            for g in genes
            if g in by_determinant
        }
        if not categories:
            out.append(replace(record, convergence_status="not_assessed"))
        elif len(categories) == 1:
            out.append(replace(record, convergence_status=next(iter(categories))))
        else:
            out.append(replace(record, convergence_status="mixed:" + ",".join(sorted(categories))))
    return out


def write_master_table(records: Sequence[MasterRecord], path) -> None:
    """Write the master table in the declared column order."""
    write_tsv(path, [r.to_row() for r in records], list(MASTER_TABLE_COLUMNS))


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    antibiotic: str,
    phenotype: Mapping[str, PhenotypeCall],
    amr: Mapping[str, Sequence[AmrDeterminant]],
    mechanisms: Mapping[str, Sequence[MechanismCall]],
    regulator_variants: Mapping[str, Sequence[RegulatorVariant]],
    structural: Mapping[str, Sequence[StructuralVariant]],
    mlst: Mapping[str, Optional[MlstCall]],
    lineages: Mapping[str, str],
    virulence: Mapping[str, Sequence[VirulenceFactor]],
    gwas_results: Sequence[GwasResult],
    convergence: Sequence[ConvergenceCall],
) -> List[MasterRecord]:
    """Stage 15 entry point."""
    records = build_master_table(
        config,
        manifest,
        antibiotic,
        phenotype,
        amr,
        mechanisms,
        regulator_variants,
        structural,
        mlst,
        lineages,
        virulence,
        gwas_results,
        convergence,
        mode=mode,
    )
    records = annotate_gwas(records, gwas_results, config)
    records = annotate_convergence(records, convergence)

    by_confidence: Dict[str, int] = {}
    for record in records:
        by_confidence[record.confidence.value] = (
            by_confidence.get(record.confidence.value, 0) + 1
        )
    LOGGER.info(
        "Stage 15: %d master records | confidence %s",
        len(records),
        ", ".join(f"{k}={v}" for k, v in sorted(by_confidence.items())),
    )
    return records
