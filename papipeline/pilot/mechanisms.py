"""Mechanism categorisation for the pilot cohort.

The critical constraint: **no gene-to-mechanism relationship is invented.**

Two separate, differently-sourced layers are used, and they are never mixed:

1. **Functional category** (ESBL, MBL, other beta-lactamase, other AMR).
   Taken verbatim from AMRFinderPlus' own ``Class`` and ``Subclass`` fields.
   These are the reference database's statements about a gene, recorded
   with the database name and version. This pipeline adds no biology.

2. ***P. aeruginosa* mechanism** (reduced permeability, efflux, AmpC
   regulation, acquired determinant). Taken **only** from
   ``config/mechanisms.tsv``. A gene absent from that table is reported as
   ``UNMAPPED`` and receives no *P. aeruginosa* mechanism, however important
   it may be.

OprD is handled by a third, explicit rule that follows the project's
scientific rules (docs/scientific_rules.md #1):

* ``oprD`` detected as a gene means an **intact locus** -> ``locus_intact``;
* only an oprD-disrupting variant or an oprD-absent call means
  ``reduced_permeability``.

An oprD *gene* hit from AMRFinderPlus is therefore never by itself reported
as reduced permeability. OprD-disrupting calls, when present, come from the
PDC genotype field, which records them explicitly (e.g. ``oprD_W339Ter``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..logging_utils import get_logger
from ..models import ClaimStatus

LOGGER = get_logger("pilot.mechanisms")

#: Category labels required by the pilot specification.
CATEGORY_PDC = "PDC"
CATEGORY_OPRD = "OprD/Porin"
CATEGORY_EFFLUX = "Efflux"
CATEGORY_ESBL = "ESBL"
CATEGORY_MBL = "MBL"
CATEGORY_OTHER_BETALACTAMASE = "Other beta-lactamase"
CATEGORY_OTHER_AMR = "Other AMR"

#: Canonical category order for reports and figures.
CATEGORY_ORDER: Tuple[str, ...] = (
    CATEGORY_PDC,
    CATEGORY_OPRD,
    CATEGORY_EFFLUX,
    CATEGORY_ESBL,
    CATEGORY_MBL,
    CATEGORY_OTHER_BETALACTAMASE,
    CATEGORY_OTHER_AMR,
)

#: UNMAPPED is reported, never silently dropped.
CATEGORY_UNMAPPED = "UNMAPPED"

#: AMRFinderPlus ``Class`` values that denote an extended-spectrum beta-
#: lactamase. These are the database's own class labels, matched verbatim.
ESBL_CLASSES: frozenset = frozenset({"ESBL"})

#: AMRFinderPlus ``Class`` values that denote a metallo-beta-lactamase.
MBL_CLASSES: frozenset = frozenset({"METALLOBETA-LACTAM", "METALLO-BETA-LACTAM"})

#: Subclass values that denote a beta-lactamase of some kind.
BETALACTAM_SUBCLASSES: frozenset = frozenset(
    {
        "BETA-LACTAM",
        "CARBAPENEM",
        "CEPHALOSPORIN",
        "PENICILLIN",
        "MONOBACTAM",
    }
)

#: Subclass values that denote an efflux component.
EFFLUX_SUBCLASSES: frozenset = frozenset({"EFFLUX"})

#: Gene-name families that identify the *P. aeruginosa* PDC (Pseudomonas
#: peptidoglycan class C beta-lactamase) locus. This is a name-prefix match
#: on the database's own gene symbol, not a functional inference.
PDC_GENE_PREFIX = "blaPDC"

#: Loci whose *loss* is the reduced-permeability mechanism.
OPRD_LOCUS = "oprD"

#: Matches a trailing allele number, e.g. ``blaVIM-6`` -> ``blaVIM``.
ALLELE_SUFFIX_RE = re.compile(r"-\d+$")


def gene_family_table(path: Optional[Path] = None) -> Dict[str, str]:
    """Load the gene-family table used to split ESBL / MBL.

    AMRFinderPlus does not label ESBL or MBL, so the split comes from this
    editable, reviewable table. Loading failures are not silently ignored:
    an unreadable table would silently degrade every beta-lactam to
    "Other beta-lactamase", so the error is raised.
    """
    if path is None:
        path = Path(__file__).resolve().parents[2] / "config" / "gene_families.tsv"
    path = Path(path)
    if not path.exists():
        LOGGER.warning("Gene family table not found at %s", path)
        return {}
    from ..io.tsv import read_tsv

    rows = read_tsv(
        path,
        required_columns=("family", "gene_prefix"),
        unique_columns=("gene_prefix",),
    )
    table = {str(row["gene_prefix"]).strip(): str(row["family"]).strip() for row in rows}
    LOGGER.info("Loaded %d gene-family entries from %s", len(table), path.name)
    return table


def family_for_gene(gene: str, families: Mapping[str, str]) -> Optional[str]:
    """Family for a gene symbol, by prefix match after allele-number stripping.

    Strips a trailing ``-<digits>`` allele number, then matches the longest
    table prefix. This is what lets one row cover a whole numbered family:
    ``blaCTX`` covers ``blaCTX-M-15``, and ``blaVIM`` covers ``blaVIM-6``.

    Longest-prefix wins so a more specific row can override a shorter one.
    """
    if not gene or not families:
        return None
    name = ALLELE_SUFFIX_RE.sub("", gene.strip())
    best: Optional[str] = None
    best_length = -1
    for prefix, family in families.items():
        if name.startswith(prefix) and len(prefix) > best_length:
            best, best_length = family, len(prefix)
    return best


def _norm(value: Optional[str]) -> str:
    return (value or "").strip().upper()


def functional_category(
    gene: str,
    annotation: Optional[Mapping[str, str]] = None,
    families: Optional[Mapping[str, str]] = None,
) -> str:
    """Assign a functional category.

    Precedence, highest first:

    1. ``blaPDC*`` -> PDC, by gene-name prefix;
    2. ``oprD`` -> OprD/Porin, by gene name;
    3. ESBL / MBL, from the editable gene-family table, because
       AMRFinderPlus does not distinguish them;
    4. ESBL / MBL, from the database's ``Class`` when it does distinguish
       them (defensive: unused on the database version used here);
    5. efflux, from the database's ``Subclass``;
    6. other beta-lactamase, from the database's ``Class``/``Subclass``;
    7. other AMR, if the database supplied any annotation;
    8. ``UNMAPPED``, if it supplied none.

    Args:
        gene: The AMRFinderPlus element symbol.
        annotation: ``{class, subclass, type}`` from the same report row.
        families: gene prefix -> family, from :func:`gene_family_table`.
    """
    annotation = annotation or {}
    gene_class = _norm(annotation.get("class"))
    subclass = _norm(annotation.get("subclass"))
    name = (gene or "").strip()

    if name.startswith(PDC_GENE_PREFIX):
        return CATEGORY_PDC
    if name == OPRD_LOCUS:
        return CATEGORY_OPRD

    family = family_for_gene(name, families or {})
    if family == CATEGORY_ESBL:
        return CATEGORY_ESBL
    if family == CATEGORY_MBL:
        return CATEGORY_MBL

    if gene_class in ESBL_CLASSES:
        return CATEGORY_ESBL
    if gene_class in MBL_CLASSES:
        return CATEGORY_MBL
    if subclass in EFFLUX_SUBCLASSES:
        return CATEGORY_EFFLUX
    if subclass in BETALACTAM_SUBCLASSES or gene_class in BETALACTAM_SUBCLASSES:
        return CATEGORY_OTHER_BETALACTAMASE
    if subclass or gene_class:
        return CATEGORY_OTHER_AMR
    return CATEGORY_UNMAPPED


def efflux_genes(config) -> List[str]:
    """Efflux loci from the knowledge table (regulators and components)."""
    return sorted(
        gene
        for gene, spec in config.mechanisms.items()
        if spec.mechanism_class == "efflux"
    )


@dataclass(frozen=True)
class MechanismAssignment:
    """One gene's categorisation, with its provenance."""

    gene: str
    functional_category: str
    pa_mechanism: Optional[str]
    pa_evidence_level: ClaimStatus
    evidence: str
    database: Optional[str]
    database_version: Optional[str]
    amr_class: Optional[str]
    amr_subclass: Optional[str]

    def to_row(self) -> Dict[str, object]:
        return {
            "Gene": self.gene,
            "Mechanism": self.functional_category,
            "PA_mechanism": self.pa_mechanism or CATEGORY_UNMAPPED,
            "PA_evidence_level": self.pa_evidence_level.value,
            "Evidence": self.evidence,
            "Database": self.database,
            "Database_version": self.database_version,
            "AMR_class": self.amr_class,
            "AMR_subclass": self.amr_subclass,
        }


def assign(
    gene: str,
    annotation: Optional[Mapping[str, str]],
    config,
    oprd_disrupted: bool = False,
    families: Optional[Mapping[str, str]] = None,
) -> MechanismAssignment:
    """Categorise one detected gene.

    Args:
        gene: AMRFinderPlus element symbol.
        annotation: ``{class, subclass, type}`` from the report.
        config: Pipeline configuration, for the knowledge table.
        oprd_disrupted: Whether an OprD-disrupting call exists for this
            isolate. Derived from an explicit variant call, never inferred
            from the presence of the oprD gene.
        families: gene prefix -> family, for the ESBL/MBL split.
    """
    annotation = annotation or {}
    category = functional_category(gene, annotation, families)
    evidence = f"AMRFinderPlus class={_norm(annotation.get('class')) or '-'}|subclass={_norm(annotation.get('subclass')) or '-'}"

    pa_mechanism: Optional[str] = None
    level = ClaimStatus.DETECTED

    # OprD: presence is an intact locus, not reduced permeability.
    if gene == OPRD_LOCUS:
        if oprd_disrupted:
            pa_mechanism = "reduced_permeability"
            evidence += "|oprd_disrupting_variant_present"
        else:
            pa_mechanism = "locus_intact"
            evidence += "|oprd_locus_detected_no_disruption_reported"
    elif config.has_gene(gene):
        spec = config.mechanism_for_gene(gene)
        pa_mechanism = spec.mechanism
        level = spec.claim_ceiling
    else:
        # No entry in config/mechanisms.tsv. No PA mechanism is assigned.
        pa_mechanism = None

    return MechanismAssignment(
        gene=gene,
        functional_category=category,
        pa_mechanism=pa_mechanism,
        pa_evidence_level=level,
        evidence=evidence,
        database="AMRFinderPlus",
        database_version=None,
        amr_class=annotation.get("class"),
        amr_subclass=annotation.get("subclass"),
    )


def oprd_calls(genotype_gene_names: Iterable[str]) -> List[str]:
    """Every OprD-related gene token in a set of PDC genotype calls.

    Matching is case-insensitive: the locus is written ``oprD`` in the
    knowledge table but appears in mixed case in the data.
    """
    return [
        token.strip()
        for token in genotype_gene_names
        if token.strip().lower().startswith(OPRD_LOCUS.lower())
    ]


def oprd_disrupting_genes(genotype_gene_names: Iterable[str]) -> Set[str]:
    """OprD variant calls that imply a disrupted locus.

    Only explicit disruptive variant syntax counts: a premature stop
    (``Ter``), a frameshift (``fs``), or an insertion/deletion. A
    missense-only change is *not* treated as disruption, because the
    project's rule is to report the variant rather than assume its effect.
    """
    disruptive: Set[str] = set()
    for name in oprd_calls(genotype_gene_names):
        tail = name[len(OPRD_LOCUS) :]
        if any(
            marker in tail
            for marker in ("fsTer", "fs", "Ter", "_del", "del", "ins", "Dup")
        ):
            disruptive.add(name)
    return disruptive


def oprd_locus_present(genotype_gene_names: Iterable[str]) -> bool:
    """Whether the PDC genotype field records any oprD feature.

    ``True`` means the field reports something at the locus, which for this
    data source is a variant call. It does not by itself mean the locus is
    functional.
    """
    return bool(oprd_calls(genotype_gene_names))
