"""Typed domain model and controlled vocabularies.

This module is the single source of truth for the vocabulary the rest of the
pipeline is allowed to emit. The scientific rules in ``docs/scientific_rules.md``
are enforced here in code, not just in prose:

* gene presence is never promoted to phenotypic resistance;
* association is never promoted to causation;
* R/I/S is never converted to a quantitative MIC or zone diameter;
* unknown values stay unknown (``None``), they are never defaulted.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Sequence

# --------------------------------------------------------------------------
# Controlled vocabularies
# --------------------------------------------------------------------------


class ClaimStatus(str, Enum):
    """Status of a claim linking a genotype to a phenotype.

    The pipeline must always state which of these applies. There is no
    "causal" member: causation is out of scope for this pipeline.
    """

    DETECTED = "DETECTED"
    PREDICTED = "PREDICTED"
    ASSOCIATED = "ASSOCIATED"
    SUPPORTED = "SUPPORTED"
    UNKNOWN = "UNKNOWN"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


#: Ordering used when collapsing several claim statuses into one per sample.
CLAIM_STATUS_PRECEDENCE: Sequence[ClaimStatus] = (
    ClaimStatus.SUPPORTED,
    ClaimStatus.ASSOCIATED,
    ClaimStatus.PREDICTED,
    ClaimStatus.DETECTED,
    ClaimStatus.UNKNOWN,
)


def strongest_claim(statuses: Sequence[ClaimStatus]) -> ClaimStatus:
    """Return the strongest claim status in ``statuses``.

    Falls back to ``UNKNOWN`` for an empty sequence so that "no evidence"
    is represented explicitly instead of being dropped.
    """
    if not statuses:
        return ClaimStatus.UNKNOWN
    present = set(statuses)
    for candidate in CLAIM_STATUS_PRECEDENCE:
        if candidate in present:
            return candidate
    return ClaimStatus.UNKNOWN


class Phenotype(str, Enum):
    """Qualitative susceptibility categories. Never converted to MIC."""

    R = "R"
    I = "I"
    S = "S"
    SDD = "SDD"
    ND = "ND"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value

    @property
    def is_quantifiable(self) -> bool:
        """Only ``ND`` and ``SDD`` carry no resistance signal at all."""
        return self in (Phenotype.R,)


class RunMode(str, Enum):
    """Execution mode.

    ``STUB`` runs the whole DAG with every rule emitting a tiny, plausible
    fake output and no bioinformatics tool being invoked. It is the default mode
    for development and CI, and it is the mode the workflow's own dry run uses.

    ``TEST`` runs the real code paths over the 20 committed synthetic fixtures.

    ``REAL`` is gated: both the configuration flag and an explicit instruction
    from the user must open before it may read real genomes.
    """

    STUB = "STUB"
    TEST = "TEST"
    REAL = "REAL"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class VariantType(str, Enum):
    """Variant classes emitted by the regulator/chromosomal screen."""

    SNV = "SNV"
    MNP = "MNP"
    INDEL = "INDEL"
    FRAMESHIFT = "FRAMESHIFT"
    PREMATURE_STOP = "PREMATURE_STOP"
    GENE_DISRUPTION = "GENE_DISRUPTION"
    GENE_ABSENCE = "GENE_ABSENCE"
    PROMOTER_ALTERATION = "PROMOTER_ALTERATION"
    INFRAME_INDEL = "INFRAME_INDEL"
    UNKNOWN = "UNKNOWN"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


#: Variant classes expected to impair a gene's function.
#:
#: This lives here rather than in ``stages.regulators`` because two modules
#: need it and must agree: the stage that decides whether a locus is
#: disrupted, and the synthetic emulator that decides which classes to emit.
#: A second, separate copy in either place is how the two drifted apart once
#: already - ``stages.regulators`` held this set while ``gene_status()``
#: ignored it and hard-coded two type names instead, so a real FRAMESHIFT or
#: PREMATURE_STOP was reported as a plain "variant" and loss-of-function read
#: as 0 on real data.
#:
#: Membership is a statement about mechanism, not about any one gene: every
#: class here truncates or removes the coding sequence. GENE_ABSENCE is
#: included because it is loss of function in the strongest sense, but callers
#: that report absence and disruption as *separate* observations must test for
#: GENE_ABSENCE first and only then fall back to this set.
DISRUPTIVE_VARIANT_TYPES = frozenset(
    {
        VariantType.FRAMESHIFT.value,
        VariantType.GENE_ABSENCE.value,
        VariantType.GENE_DISRUPTION.value,
        VariantType.PREMATURE_STOP.value,
    }
)


class StructuralCallStatus(str, Enum):
    """Explicitly separates confirmed structural variants from candidates."""

    CONFIRMED = "confirmed"
    CANDIDATE = "candidate"
    NOT_ASSESSABLE = "not_assessable"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class SvType(str, Enum):
    DELETION = "deletion"
    INSERTION = "insertion"
    REARRANGEMENT = "rearrangement"
    INSERTION_SEQUENCE = "insertion_sequence"
    GENE_DISRUPTION = "gene_disruption"
    PROMOTER_REGION_CHANGE = "promoter_region_change"
    MGE_ASSOCIATION = "mge_association"
    INVERSION = "inversion"
    DUPLICATION = "duplication"
    UNKNOWN = "unknown"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class MechanismClass(str, Enum):
    """Mechanism categories used by ``config/mechanisms.tsv``."""

    REDUCED_PERMEABILITY = "reduced_permeability"
    EFFLUX = "efflux"
    AMPC_REGULATION = "ampc_regulation"
    ACQUIRED_DETERMINANT = "acquired_determinant"
    STRUCTURAL = "structural"
    ENZYME_MODIFICATION = "enzyme_modification"
    UNKNOWN = "unknown"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class GeneRole(str, Enum):
    """Role a gene plays with respect to susceptibility."""

    EFFLUX_REGULATOR = "efflux_regulator"
    PERMEABILITY = "permeability"
    BETA_LACTAMASE_REGULATOR = "beta_lactamase_regulator"
    TARGET_GENE = "target_gene"
    ACQUIRED_DETERMINANT = "acquired_determinant"
    UNKNOWN = "unknown"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class ConvergenceCategory(str, Enum):
    """Categories from stage 13."""

    LINEAGE_ASSOCIATED = "lineage_associated"
    RECURRENT_CONVERGENT = "recurrent_convergent"
    WIDESPREAD_BACKGROUND = "widespread_background"
    RARE_ISOLATED = "rare_isolated"
    UNKNOWN = "unknown"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


class MlstStatus(str, Enum):
    TYPED = "typed"
    PARTIAL = "partial"
    ALLELE_INCOMPLETE = "allele_incomplete"
    NO_CALL = "no_call"

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.value


# --------------------------------------------------------------------------
# Record dataclasses
# --------------------------------------------------------------------------


def _clean(value: Any) -> Any:
    """Normalise empty strings and pandas NA to ``None`` (missing stays missing)."""
    if value is None:
        return None
    if isinstance(value, float) and value != value:  # NaN
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"na", "nan", "none", "null", "."}:
        return None
    return text


@dataclass(frozen=True)
class Sample:
    """A single analysis unit. ``sample_id`` is the join key everywhere.

    The identifier is the versioned assembly accession, which is already the
    genome directory name and the FASTA stem, so the assembly join needs no
    lookup table. The remaining fields are *provenance*, not join keys: they are
    carried so the chain from assembly back to isolate is inspectable, and they
    are required to be unambiguous (at most one assembly per value).
    """

    sample_id: str
    assembly_path: Optional[str] = None
    source: str = "unknown"
    assembly: Optional[str] = None
    isolate_id: Optional[str] = None
    biosample: Optional[str] = None
    isolate_name: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "assembly_path": self.assembly_path,
            "source": self.source,
            "assembly": self.assembly,
            "isolate_id": self.isolate_id,
            "biosample": self.biosample,
            "isolate_name": self.isolate_name,
        }


@dataclass(frozen=True)
class AssemblyQC:
    """Stage 1 output. Missing metrics stay ``None`` rather than becoming 0."""

    sample_id: str
    assembly_size: Optional[int] = None
    contig_count: Optional[int] = None
    n50: Optional[int] = None
    n90: Optional[int] = None
    largest_contig: Optional[int] = None
    gc_content: Optional[float] = None
    ambiguous_bases: Optional[int] = None
    basic_quality_status: str = "unknown"
    species_confirmation: Optional[str] = None
    contamination_status: Optional[str] = None
    completeness_status: Optional[str] = None
    duplicate_status: Optional[str] = None
    flag_reasons: str = ""

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class AnnotationRecord:
    """Stage 2 standardised annotation row."""

    sample_id: str
    contig_id: str
    gene_id: str
    gene_name: Optional[str]
    product: Optional[str]
    gene_type: str
    start: Optional[int]
    end: Optional[int]
    strand: Optional[str]
    annotation_source: str

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MlstCall:
    """Stage 3 output. ``alleles`` keeps per-locus detail for traceability."""

    sample_id: str
    sequence_type: Optional[str]
    alleles: Dict[str, str] = field(default_factory=dict)
    mlst_status: str = MlstStatus.NO_CALL.value
    mlst_scheme: Optional[str] = None
    allele_database: Optional[str] = None

    def allele_string(self) -> str:
        return ";".join(f"{k}:{v}" for k, v in sorted(self.alleles.items()))

    def to_row(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "ST": self.sequence_type,
            "alleles": self.allele_string(),
            "MLST_status": self.mlst_status,
            "mlst_scheme": self.mlst_scheme,
            "allele_database": self.allele_database,
        }


@dataclass(frozen=True)
class AmrDeterminant:
    """Stage 4 output.

    ``evidence_source`` records the tool, ``database``/``database_version`` the
    reference. ``confidence`` is the tool-reported confidence, never an
    inferred probability of resistance.
    """

    sample_id: str
    antibiotic: str
    determinant: str
    gene: Optional[str]
    variant: Optional[str]
    determinant_type: str
    mechanism: Optional[str]
    evidence_source: str
    database: str
    database_version: str
    confidence: Optional[str]
    claim_status: ClaimStatus = ClaimStatus.DETECTED
    identity_pct: Optional[float] = None
    coverage_pct: Optional[float] = None

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["claim_status"] = self.claim_status.value
        return row


@dataclass(frozen=True)
class MechanismCall:
    """Stage 5 output: a determinant/variant interpreted as a mechanism."""

    sample_id: str
    antibiotic: str
    determinant: str
    mechanism: str
    evidence_level: ClaimStatus
    gene: Optional[str] = None
    gene_type: Optional[str] = None
    notes: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["evidence_level"] = self.evidence_level.value
        return row


@dataclass(frozen=True)
class RegulatorVariant:
    """Stage 6 output for resistance-associated chromosomal loci."""

    sample_id: str
    gene: str
    variant: str
    variant_type: str
    position: Optional[int]
    reference: Optional[str]
    alternate: Optional[str]
    effect: Optional[str]
    mechanism: Optional[str]
    confidence: Optional[str]
    evidence_source: str = "synthetic"
    call_status: ClaimStatus = ClaimStatus.DETECTED

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["call_status"] = self.call_status.value
        return row


@dataclass(frozen=True)
class StructuralVariant:
    """Stage 7 output. ``call_status`` must never be upgraded silently."""

    sample_id: str
    variant_id: str
    variant_type: str
    position: Optional[int]
    affected_gene: Optional[str]
    size: Optional[int]
    evidence: str
    confidence: Optional[str]
    call_status: StructuralCallStatus
    mge: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["call_status"] = self.call_status.value
        return row


@dataclass(frozen=True)
class VirulenceFactor:
    """Stage 8 output. Independent of AMR analysis by construction."""

    sample_id: str
    virulence_factor: str
    gene: str
    category: Optional[str]
    database: str
    database_version: str
    confidence: Optional[str]
    identity_pct: Optional[float] = None

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PhenotypeCall:
    """Stage 11 input/output.

    ``mic`` and ``mic_unit`` are only populated from measured data. They are
    never synthesised from R/I/S.
    """

    sample_id: str
    antibiotic: str
    phenotype: Phenotype
    mic: Optional[float] = None
    mic_unit: Optional[str] = None
    zone_diameter: Optional[float] = None
    zone_unit: Optional[str] = None
    source: str = "unknown"
    # AST provenance, carried per row where the source provides it. An MIC from
    # one laboratory's method is not automatically comparable with another's, so
    # the run reports how many rows lack each field rather than assuming they
    # are all alike.
    ast_method: Optional[str] = None
    ast_standard: Optional[str] = None
    ast_edition: Optional[str] = None
    # The continuous trait the association model consumes. Computed from a
    # measured MIC only; never from a category.
    log2_mic: Optional[float] = None
    # True only when this row's category was *derived* from a measured MIC
    # against a configured standard. False means the source call stands.
    reinterpreted_from_mic: bool = False

    def to_row(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "antibiotic": self.antibiotic,
            "phenotype": self.phenotype.value,
            "MIC": self.mic,
            "MIC_unit": self.mic_unit,
            "zone_diameter": self.zone_diameter,
            "zone_unit": self.zone_unit,
            "source": self.source,
            "ast_method": self.ast_method,
            "ast_standard": self.ast_standard,
            "ast_edition": self.ast_edition,
            "log2_MIC": self.log2_mic,
            "reinterpreted_from_MIC": self.reinterpreted_from_mic,
        }


@dataclass(frozen=True)
class GwasResult:
    """Stage 12 output. ``adjusted_p_value`` is always present alongside p."""

    feature: str
    feature_type: str
    effect: Optional[float]
    p_value: Optional[float]
    adjusted_p_value: Optional[float]
    effect_size: Optional[float]
    frequency: Optional[float]
    lineage_distribution: Dict[str, int] = field(default_factory=dict)
    model: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["lineage_distribution"] = ";".join(
            f"{k}={v}" for k, v in sorted(self.lineage_distribution.items())
        )
        return row


@dataclass(frozen=True)
class ConvergenceCall:
    """Stage 13 output. Convergence is always reported with lineage context."""

    determinant: str
    independent_lineages: int
    branch_count: int
    distribution: Dict[str, int]
    convergence_category: ConvergenceCategory

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["distribution"] = ";".join(
            f"{k}={v}" for k, v in sorted(self.distribution.items())
        )
        row["convergence_category"] = self.convergence_category.value
        return row


@dataclass(frozen=True)
class Cooccurrence:
    """Stage 14 output. Association statistics only, never causal claims."""

    feature_a: str
    feature_b: str
    feature_type: str
    n_a: int
    n_b: int
    n_both: int
    statistic: str
    statistic_value: Optional[float]
    adjusted_p_value: Optional[float]
    interpretation_limit: str = "association_only_not_causal"

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MasterRecord:
    """Stage 15 integration row: the genotype-mechanism-phenotype model."""

    sample_id: str
    antibiotic: str
    phenotype: Optional[str] = None
    amr_gene: Optional[str] = None
    amr_variant: Optional[str] = None
    chromosomal_mutation: Optional[str] = None
    regulator: Optional[str] = None
    mechanism: Optional[str] = None
    structural_variant: Optional[str] = None
    mlst: Optional[str] = None
    lineage: Optional[str] = None
    virulence_profile: Optional[str] = None
    gwas_feature: Optional[str] = None
    gwas_status: Optional[str] = None
    convergence_status: Optional[str] = None
    confidence: ClaimStatus = ClaimStatus.UNKNOWN
    evidence_notes: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        row = asdict(self)
        row["confidence"] = self.confidence.value
        return row


#: Column order for the stage 15 master table.
MASTER_TABLE_COLUMNS: Sequence[str] = tuple(MasterRecord.__dataclass_fields__.keys())


def records_to_rows(records: Sequence[Any]) -> List[Dict[str, Any]]:
    """Convert a sequence of record dataclasses to a list of dict rows."""
    return [r.to_row() for r in records]


def coerce_optional_str(value: Any) -> Optional[str]:
    """Public wrapper around the missing-value normaliser."""
    return _clean(value)


def normalise_mapping(values: Mapping[str, Any]) -> Dict[str, Optional[str]]:
    """Apply missing-value normalisation across a mapping."""
    return {k: _clean(v) for k, v in values.items()}
