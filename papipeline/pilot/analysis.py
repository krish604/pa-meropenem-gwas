"""Pilot-100 analysis: phenotype join, gene and mechanism tables, pairs.

Consumes the AMRFinderPlus sweep (:mod:`papipeline.pilot.amr_detect`) and the
PDC metadata join (:mod:`papipeline.pilot.cohort`) and produces every
gene-level and mechanism-level deliverable for imipenem and meropenem.

Scientific constraints, all inherited from the main pipeline:

* a detected gene is ``DETECTED``; it is never a resistance call
  (docs/scientific_rules.md #1);
* R/I/S is never converted into an MIC or a zone diameter (#6);
* co-occurrence is association only, never causation (#2);
* a gene absent from ``config/mechanisms.tsv`` gets no *P. aeruginosa*
  mechanism and is reported as ``UNMAPPED`` rather than dropped;
* missing phenotype stays missing; it is never imputed.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..errors import PipelineError
from ..io.tsv import read_tsv, write_tsv
from ..logging_utils import get_logger
from ..models import Phenotype
from . import mechanisms as mech_mod
from . import pdc_fields

LOGGER = get_logger("pilot.analysis")

#: Antibiotics analysed, in report order.
ANTIBIOTICS: Tuple[str, ...] = ("imipenem", "meropenem")

MISSING = "MISSING"
NOT_DETECTED = "not_detected"


# --------------------------------------------------------------------------
# Cohort loading
# --------------------------------------------------------------------------


@dataclass
class CohortMember:
    """One pilot sample with its metadata and detection results."""

    pilot_id: str
    assembly: str
    genome_path: str
    isolate: str
    biosample: str
    bioproject: str
    pdc_match: str
    pdc_present: str
    ast: pdc_fields.AstField
    genotypes: pdc_fields.GenotypeField
    #: False when the assembly file is empty, corrupt or undersized. Set by
    #: --prepare-only from the integrity check; such samples are excluded
    #: from analysis and never substituted.
    analysis_included: bool = True
    exclusion_reason: str = "ok"
    genes: Set[str] = field(default_factory=set)
    gene_annotations: Dict[str, Dict[str, str]] = field(default_factory=dict)
    #: OprD locus recorded by the PDC genotype field. AMRFinderPlus does not
    #: report oprD (a porin is not an AMR gene), so OprD status comes from
    #: the PDC genotype calls and is labelled as such.
    oprd_locus_reported: bool = False
    oprd_variants: List[str] = field(default_factory=list)
    oprd_disrupted: bool = False

    def phenotype(self, antibiotic: str) -> Optional[Phenotype]:
        return self.ast.get(antibiotic)

    def phenotype_label(self, antibiotic: str) -> str:
        value = self.phenotype(antibiotic)
        return value.value if value else MISSING

    def categories(
        self, families: Optional[Mapping[str, str]] = None
    ) -> Set[str]:
        """Functional mechanism categories present for this isolate."""
        out: Set[str] = set()
        for gene in self.genes:
            out.add(
                mech_mod.functional_category(
                    gene, self.gene_annotations.get(gene), families
                )
            )
        if self.oprd_locus_reported:
            out.add(mech_mod.CATEGORY_OPRD)
        return out


def load_cohort(
    manifest_path: Path, metadata_path: Path
) -> List[CohortMember]:
    """Load the prepared cohort and parse both PDC free-text fields."""
    manifest_rows = read_tsv(manifest_path, required_columns=("Pilot_ID", "Assembly"))
    metadata_rows = {
        str(row["Assembly"]): row
        for row in read_tsv(
            metadata_path, required_columns=("Pilot_ID", "Assembly")
        )
    }

    members: List[CohortMember] = []
    for row in manifest_rows:
        assembly = str(row["Assembly"])
        record = metadata_rows.get(assembly, {})
        members.append(
            CohortMember(
                pilot_id=str(row["Pilot_ID"]),
                assembly=assembly,
                genome_path=str(row.get("Genome_path") or ""),
                isolate=str(record.get("Isolate") or MISSING),
                biosample=str(record.get("BioSample") or MISSING),
                bioproject=str(record.get("BioProject") or MISSING),
                pdc_match=str(record.get("PDC_metadata_match") or MISSING),
                pdc_present=str(record.get("PDC_present") or MISSING),
                analysis_included=str(record.get("analysis_included") or "TRUE").upper()
                == "TRUE",
                exclusion_reason=str(record.get("exclusion_reason") or "ok"),
                ast=pdc_fields.parse_ast_field(record.get("AST_phenotypes_raw")),
                genotypes=pdc_fields.parse_genotype_field(
                    record.get("AMR_genotypes_raw")
                ),
            )
        )

    LOGGER.info("Loaded %d pilot cohort members", len(members))
    return members


def attach_detections(
    members: Sequence[CohortMember], detection, config
) -> None:
    """Attach AMRFinderPlus detections to each cohort member, in place."""
    for member in members:
        records = detection.records.get(member.pilot_id, [])
        member.genes = {r.gene for r in records if r.gene}
        member.gene_annotations = detection.class_map.get(member.pilot_id, {})
        names = [c.gene for c in member.genotypes.calls]
        member.oprd_variants = mech_mod.oprd_calls(names)
        member.oprd_locus_reported = bool(member.oprd_variants)
        member.oprd_disrupted = bool(mech_mod.oprd_disrupting_genes(names))


# --------------------------------------------------------------------------
# Step 9: phenotype
# --------------------------------------------------------------------------


def phenotype_rows(
    members: Sequence[CohortMember], antibiotic: str
) -> Tuple[List[Dict[str, object]], Counter]:
    """Per-assembly phenotype rows and the category distribution.

    A category outside the allowed set, or a duplicated antibiotic key, is
    ``MISSING`` rather than coerced.
    """
    rows: List[Dict[str, object]] = []
    counts: Counter = Counter()
    for member in members:
        label = member.phenotype_label(antibiotic)
        counts[label] += 1
        rows.append(
            {
                "Pilot_ID": member.pilot_id,
                "Assembly": member.assembly,
                "Isolate": member.isolate,
                "phenotype": label,
                "availability": member.ast.availability(antibiotic),
                # Explicitly empty: no MIC or zone value exists in this source
                # and none is derived. See scientific rule 6.
                "MIC": ".",
                "MIC_unit": ".",
            }
        )
    return rows, counts


# --------------------------------------------------------------------------
# Step 10: gene tables
# --------------------------------------------------------------------------


def gene_summary_rows(
    members: Sequence[CohortMember],
    antibiotic: str,
    config,
    families: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, object]]:
    """One row per (isolate, gene) for the antibiotic.

    Only detected genes are listed. Absence of a row means the gene was not
    detected; the per-isolate gene count column makes that explicit.
    """
    rows: List[Dict[str, object]] = []
    for member in members:
        for gene in sorted(member.genes):
            assignment = mech_mod.assign(
                gene,
                member.gene_annotations.get(gene),
                config,
                oprd_disrupted=member.oprd_disrupted,
                families=families,
            )
            rows.append(
                {
                    "Pilot_ID": member.pilot_id,
                    "Assembly": member.assembly,
                    "Isolate": member.isolate,
                    "AST_phenotype": member.phenotype_label(antibiotic),
                    "Gene": gene,
                    "Mechanism": assignment.functional_category,
                    "PA_mechanism": assignment.pa_mechanism or MISSING,
                    "PA_evidence_level": assignment.pa_evidence_level.value,
                    "Detected": "TRUE",
                    "Evidence": assignment.evidence,
                    "AMR_class": assignment.amr_class or MISSING,
                    "AMR_subclass": assignment.amr_subclass or MISSING,
                    "n_genes_detected_in_isolate": len(member.genes),
                }
            )
    return rows


def gene_frequency_rows(
    members: Sequence[CohortMember],
    antibiotic: str,
    config,
    families: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, object]]:
    """Per-gene frequency, split by phenotype category.

    Carriers are tracked by isolate ID, not by phenotype label: counting
    labels would cap every gene at the number of distinct categories rather
    than the number of isolates carrying it.

    An isolate with no usable category for this antibiotic is bucketed as
    ``MISSING`` rather than dropped, so overall frequency stays visible and
    the missing split is explicit in the table.
    """
    carriers_by_gene: Dict[str, Set[str]] = defaultdict(set)
    category_of: Dict[str, str] = {}
    category_totals: Counter = Counter()

    for member in members:
        label = member.phenotype_label(antibiotic)
        category_of[member.pilot_id] = label
        category_totals[label] += 1
        for gene in member.genes:
            carriers_by_gene[gene].add(member.pilot_id)

    total_isolates = sum(category_totals.values())
    rows: List[Dict[str, object]] = []
    for gene in sorted(carriers_by_gene):
        assignment = mech_mod.assign(
            gene, _first_annotation(members, gene), config, families=families
        )
        carriers = carriers_by_gene[gene]
        per_category: Dict[str, int] = {}
        for category in category_totals:
            per_category[category] = sum(
                1 for pilot_id in carriers if category_of[pilot_id] == category
            )
        row: Dict[str, object] = {
            "Gene": gene,
            "Mechanism": assignment.functional_category,
            "PA_mechanism": assignment.pa_mechanism or MISSING,
            "n_isolates_detected": len(carriers),
            "n_isolates_in_cohort": total_isolates,
            "frequency": round(len(carriers) / total_isolates, 6)
            if total_isolates
            else None,
            "n_categories_present": sum(1 for v in per_category.values() if v > 0),
        }
        for category in sorted(per_category):
            row[f"n_{category}"] = per_category[category]
            row[f"freq_{category}"] = (
                round(per_category[category] / category_totals[category], 6)
                if category_totals[category]
                else None
            )
        rows.append(row)
    return rows


def _first_annotation(
    members: Sequence[CohortMember], gene: str
) -> Optional[Dict[str, str]]:
    """First available class/subclass annotation for a gene across the cohort."""
    for member in members:
        annotation = member.gene_annotations.get(gene)
        if annotation:
            return annotation
    return None


# --------------------------------------------------------------------------
# Step 11: mechanism tables
# --------------------------------------------------------------------------


def mechanism_summary_rows(
    members: Sequence[CohortMember],
    antibiotic: str,
    families: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, object]]:
    """Per-mechanism counts, split by phenotype category.

    A mechanism is counted once per isolate.
    """
    category_totals: Counter = Counter()
    category_of: Dict[str, str] = {}
    carriers_by_mech: Dict[str, Set[str]] = defaultdict(set)
    genes_per_mech: Dict[str, Set[str]] = defaultdict(set)
    carriers_per_mech: Dict[str, Set[str]] = defaultdict(set)

    for member in members:
        label = member.phenotype_label(antibiotic)
        category_of[member.pilot_id] = label
        category_totals[label] += 1
        for category in member.categories(families):
            # Carriers tracked by isolate ID, not by phenotype label.
            carriers_by_mech[category].add(member.pilot_id)
        for gene in member.genes:
            category = mech_mod.functional_category(
                gene, member.gene_annotations.get(gene), families
            )
            genes_per_mech[category].add(gene)
            carriers_per_mech[category].add(member.pilot_id)
        if member.oprd_locus_reported:
            carriers_per_mech[mech_mod.CATEGORY_OPRD].add(member.pilot_id)

    rows: List[Dict[str, object]] = []
    for mechanism in sorted(carriers_by_mech):
        carriers = carriers_by_mech[mechanism]
        row: Dict[str, object] = {
            "Mechanism": mechanism,
            "n_isolates": len(carriers),
            "n_isolates_in_cohort": sum(category_totals.values()),
            "genes_in_category": ",".join(sorted(genes_per_mech[mechanism])),
            "n_genes_in_category": len(genes_per_mech[mechanism]),
        }
        for category in sorted(category_totals):
            n = sum(1 for pilot_id in carriers if category_of[pilot_id] == category)
            row[f"n_{category}"] = n
            row[f"freq_{category}"] = (
                round(n / category_totals[category], 6)
                if category_totals[category]
                else None
            )
        rows.append(row)

    n_oprd_reported = sum(1 for m in members if m.oprd_locus_reported)
    n_oprd_disrupted = sum(1 for m in members if m.oprd_disrupted)

    def _oprd_row(label, pilot_ids):
        """An OprD row with the same R/S/I breakdown as every other row."""
        carriers = set(pilot_ids)
        row: Dict[str, object] = {
            "Mechanism": label,
            "n_isolates": len(carriers),
            "n_isolates_in_cohort": sum(category_totals.values()),
            "genes_in_category": "oprD",
            "n_genes_in_category": 1,
        }
        for category in sorted(category_totals):
            n = sum(1 for p in carriers if category_of[p] == category)
            row[f"n_{category}"] = n
            row[f"freq_{category}"] = (
                round(n / category_totals[category], 6)
                if category_totals[category]
                else None
            )
        return row

    if n_oprd_reported:
        rows.append(
            _oprd_row(
                "OprD locus reported (PDC genotype call)",
                [m.pilot_id for m in members if m.oprd_locus_reported],
            )
        )
    if n_oprd_disrupted:
        rows.append(
            _oprd_row(
                "PA:reduced_permeability (OprD disruption)",
                [m.pilot_id for m in members if m.oprd_disrupted],
            )
        )
    return rows


# --------------------------------------------------------------------------
# Step 12: gene pairs
# --------------------------------------------------------------------------


def gene_pair_rows(
    members: Sequence[CohortMember], antibiotic: str
) -> List[Dict[str, object]]:
    """Unordered gene pairs co-occurring in the same isolate.

    Each pair is counted **once per isolate** regardless of how many copies
    or contigs carry either gene. Co-occurrence is association only.
    """
    pair_isolates: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for member in members:
        genes = sorted(member.genes)
        for i, gene_a in enumerate(genes):
            for gene_b in genes[i + 1 :]:
                pair_isolates[(gene_a, gene_b)].add(member.pilot_id)

    rows: List[Dict[str, object]] = []
    for (gene_a, gene_b), isolates in sorted(
        pair_isolates.items(), key=lambda kv: (-len(kv[1]), kv[0])
    ):
        rows.append(
            {
                "Gene_A": gene_a,
                "Gene_B": gene_b,
                "n_isolates_with_both": len(isolates),
                "n_isolates_cohort": len(members),
                "frequency": round(len(isolates) / len(members), 6) if members else None,
                "interpretation_limit": "co_occurrence_association_only_not_causal",
            }
        )
    return rows


# --------------------------------------------------------------------------
# Step 13: mechanism combinations
# --------------------------------------------------------------------------


def mechanism_combination_rows(
    members: Sequence[CohortMember],
    antibiotic: str,
    families: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, object]]:
    """Observed mechanism combinations, e.g. ``PDC + OprD + Efflux``.

    Only combinations actually observed are emitted. Ordered by frequency
    then alphabetically, so the table is stable.
    """
    combination_isolates: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)
    combination_genes: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)

    for member in members:
        categories = member.categories(families)
        if not categories:
            continue
        ordered = tuple(c for c in mech_mod.CATEGORY_ORDER if c in categories)
        for size in range(1, len(ordered) + 1):
            for subset in _combinations(ordered, size):
                combination_isolates[subset].add(member.pilot_id)
                for gene in member.genes:
                    if mech_mod.functional_category(
                        gene, member.gene_annotations.get(gene), families
                    ) in subset:
                        combination_genes[subset].add(gene)

    rows: List[Dict[str, object]] = []
    for combination, isolates in sorted(
        combination_isolates.items(), key=lambda kv: (-len(kv[1]), kv[0])
    ):
        rows.append(
            {
                "Combination": " + ".join(combination),
                "n_mechanisms": len(combination),
                "n_isolates": len(isolates),
                "n_isolates_cohort": len(members),
                "frequency": round(len(isolates) / len(members), 6) if members else None,
                "genes": ",".join(sorted(combination_genes[combination])),
                "interpretation_limit": "co_occurrence_association_only_not_causal",
            }
        )
    return rows


def _combinations(items: Sequence[str], size: int) -> Iterable[Tuple[str, ...]]:
    from itertools import combinations

    return combinations(items, size)
