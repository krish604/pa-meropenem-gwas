"""The five layer feature builders (ticket 15 plus the L1-L5 encoding).

Each builder turns the source tables into named features over one cohort, and
each feature carries the four things ``feature_dictionary.tsv`` will publish:
its name, its layer, the source tokens it was read from, and the rule that
produces its values.

What the ticket fixes as *science* lives here as documented constants: the
carbapenemase groups, the PDC regulators, the five efflux regulators, the
target genes, and the oprD loss tiers. What is *environment* - thresholds,
paths, threads - does not: those come from ``config/science.yaml`` through
:mod:`papipeline.layers.settings`. Two groupings the ticket does not spell out
come from the knowledge tables rather than from this module:

* the efflux pump groups are read from ``config/mechanisms.tsv``
  ``biological_role`` (the first ``mex*`` operon token), so adding a regulator
  is a table edit;
* the "disruptive" predicate is
  :func:`papipeline.stages.regulators.is_disruptive` over
  ``models.DISRUPTIVE_VARIANT_TYPES`` - the pipeline's single definition, not
  a second copy.

Layer 5's QRDR features use the residue positions AMRFinderPlus emits in PAO1
numbering. The repository defines no QRDR *range* anywhere, so
``any_QRDR``/``QRDR_count`` count parsable residue positions on the two
quinolone targets rather than inventing a codon list; see
``.build/layer-encoder.registry-notes.md``.
"""

from __future__ import annotations

import re
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..errors import SampleIdError
from ..logging_utils import get_logger
from ..models import SvType, VariantType
from ..stages import regulators as reg_stage
from .features import COUNT, Feature, LayerBundle
from .inputs import Sources

LOGGER = get_logger("layers")

#: Layer 1: the acquired carbapenemase groups the ticket names. ``any_MBL``
#: unions exactly these three families - SPM and other metallo-enzymes are not
#: folded in, because the ticket says NDM/VIM/IMP.
MBL_FAMILIES: Tuple[str, ...] = ("NDM", "VIM", "IMP")
GROUP_FAMILIES: Tuple[str, ...] = ("KPC", "GES", "OXA")

#: Layer 2: the chromosomal ampC regulators that raise PDC expression.
PDC_REGULATORS: Tuple[str, ...] = ("ampD", "ampR", "dacB")

#: Layer 3.
OPRD_GENE = "oprD"

#: Layer 4: the five regulators the ticket names. ``mexT`` and ``mexS`` are in
#: ``config/regulators.tsv`` but not in the ticket, so they are not encoded.
L4_REGULATOR_GENES: Tuple[str, ...] = ("mexR", "nalC", "nalD", "nfxB", "mexZ")

#: Layer 5: the targets, and the two quinolone targets whose residue positions
#: the QRDR features count.
L5_TARGET_GENES: Tuple[str, ...] = ("gyrA", "parC", "ftsI")
QRDR_TARGET_GENES: Tuple[str, ...] = ("gyrA", "parC")

#: oprD tier 1 - a documented lesion: stop, frameshift, IS insertion or large
#: deletion. ``GENE_DISRUPTION`` is the screen's own knockout class and sits
#: here with them; the fixture's disruption rows are insertion-driven.
TIER1_REGULATOR_TYPES = frozenset(
    {
        VariantType.FRAMESHIFT.value,
        VariantType.PREMATURE_STOP.value,
        VariantType.GENE_DISRUPTION.value,
    }
)

#: oprD tier 2 - probable loss: a coding indel the screen did not resolve to a
#: truncation. Probable rather than documented, which is why it is its own
#: tier (``docs/scientific_rules.md``: an in-frame deletion that preserves the
#: native stop is not loss of function).
TIER2_REGULATOR_TYPES = frozenset(
    {VariantType.INDEL.value, VariantType.INFRAME_INDEL.value}
)

#: Structural variants that remove or break the locus, versus every other SV
#: type - a duplication or a promoter-region change is not loss of the gene.
STRUCTURAL_LOSS_TYPES = frozenset(
    {
        SvType.INSERTION.value,
        SvType.INSERTION_SEQUENCE.value,
        SvType.DELETION.value,
        SvType.GENE_DISRUPTION.value,
    }
)
SV_CONFIRMED = "confirmed"
SV_CANDIDATE = "candidate"
SV_NOT_ASSESSABLE = "not_assessable"

#: The canonical layer-3 features. An individual oprD variant whose name would
#: collide with one of these is dropped with a warning rather than silently
#: overwriting a tier.
OPRD_CANONICAL = (
    f"L3_{OPRD_GENE}_absent",
    f"L3_{OPRD_GENE}_LoF_tier1",
    f"L3_{OPRD_GENE}_LoF_tier2",
    f"L3_{OPRD_GENE}_off_any",
)

#: Sequence variants that count as a missense rather than a lesion.
MISSENSE_TYPES = frozenset({VariantType.SNV.value, VariantType.MNP.value})

#: Residue positions AMRFinderPlus writes as amino-acid changes: ``S83L``,
#: ``V359Lfs*45``, ``D146del``, ``146ins``. Anything else - an allele name, a
#: contig coordinate, ``-`` - parses to no position rather than to a wrong one.
_VARIANT_POSITION_PATTERNS: Tuple["re.Pattern[str]", ...] = (
    re.compile(r"^[A-Za-z](\d+)[A-Za-z*]"),
    re.compile(r"^[A-Za-z]?(\d+)(?:del|ins|dup)", re.IGNORECASE),
)

#: The operon token read out of ``config/mechanisms.tsv`` ``biological_role``.
_PUMP_TOKEN = re.compile(r"mex[A-Za-z]+", re.IGNORECASE)

_AMR_SOURCE = ("amr_determinants.tsv",)
_REGULATOR_SOURCE = ("regulator_variants.tsv",)
_SV_SOURCE = ("structural_variants.tsv",)
_STRUCTURAL_SOURCE = ("oprd_locus/structural_calls.tsv",)
_MECHANISM_SOURCE = ("mechanisms.tsv",)


def _cell(row: Mapping[str, object], column: str) -> str:
    value = row.get(column)
    if value is None:
        return ""
    return str(value).strip()


def _column(sample_ids: Sequence[str]) -> Dict[str, int]:
    return {sample_id: 0 for sample_id in sample_ids}


def _check_sample(sample_id: str, sample_ids: Sequence[str], source: str) -> None:
    if sample_id in sample_ids:
        return
    raise SampleIdError(
        f"{source} names sample {sample_id!r}, which is not in this bundle's "
        "cohort. Sample-ID mismatches fail loudly and are never dropped.",
        source=source,
        sample_id=sample_id,
    )


def _amr_rows(
    sources: Sources, antibiotic: Optional[str] = None
) -> List[Mapping[str, object]]:
    """AMR rows for the antibiotic under analysis."""
    target = sources.antibiotic if antibiotic is None else antibiotic
    return [row for row in sources.amr_rows if _cell(row, "antibiotic") == target]


def _gene_token(row: Mapping[str, object]) -> str:
    """The gene symbol, upper-cased, with a leading ``bla`` removed.

    AMRFinderPlus writes ``blaNDM-1``; the group prefixes the ticket names are
    the family part. When the ``gene`` cell is empty the ``determinant`` cell
    carries the same symbol before its ``~`` qualifier.
    """
    gene = _cell(row, "gene") or _cell(row, "determinant").split("~")[0]
    gene = gene.strip()
    if gene[:3].lower() == "bla":
        gene = gene[3:]
    return gene.upper()


def _residue_position(variant: str) -> Optional[int]:
    """The residue number in an amino-acid variant, or ``None``."""
    text = (variant or "").strip()
    if not text or text in {"-", "."}:
        return None
    for pattern in _VARIANT_POSITION_PATTERNS:
        match = pattern.match(text)
        if match:
            return int(match.group(1))
    return None


def _sanitised_feature_name(prefix: str, token: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z._-]", "_", token).strip("_")
    return f"{prefix}{cleaned}" if cleaned else ""


# ---------------------------------------------------------------------------
# Layer 1 - acquired determinants
# ---------------------------------------------------------------------------


def build_layer1(
    sources: Sources, sample_ids: Sequence[str]
) -> Tuple[List[Feature], Dict[str, Dict[str, int]]]:
    """``L1_any_MBL`` plus one column per carbapenemase group the ticket names."""
    sample_ids = tuple(sample_ids)
    values: Dict[str, Dict[str, int]] = {
        "L1_any_MBL": _column(sample_ids),
        "L1_KPC": _column(sample_ids),
        "L1_GES": _column(sample_ids),
        "L1_OXA": _column(sample_ids),
    }
    for row in _amr_rows(sources):
        token = _gene_token(row)
        if not token:
            continue
        sample_id = _cell(row, "sample_id")
        _check_sample(sample_id, sample_ids, "amr_determinants.tsv")
        if token.startswith(MBL_FAMILIES):
            values["L1_any_MBL"][sample_id] = 1
        for family in GROUP_FAMILIES:
            if token.startswith(family):
                values[f"L1_{family}"][sample_id] = 1

    antibiotic = sources.antibiotic
    features = [
        Feature(
            name="L1_any_MBL",
            layer=1,
            source_tokens=_AMR_SOURCE,
            rule=(
                f"1 when any {antibiotic} determinant row's gene symbol (a "
                f"leading 'bla' removed) begins {' or '.join(MBL_FAMILIES)}; "
                "the ticket's metallo-beta-lactamase union, and nothing else"
            ),
        ),
        *(
            Feature(
                name=f"L1_{family}",
                layer=1,
                source_tokens=_AMR_SOURCE,
                rule=(
                    f"1 when any {antibiotic} determinant row's gene symbol "
                    f"(a leading 'bla' removed) begins {family}"
                ),
            )
            for family in GROUP_FAMILIES
        ),
    ]
    return features, values


# ---------------------------------------------------------------------------
# Layer 2 - PDC allele, positions and regulators
# ---------------------------------------------------------------------------


def build_layer2(
    sources: Sources, sample_ids: Sequence[str]
) -> Tuple[List[Feature], Dict[str, Dict[str, int]]]:
    """PDC allele columns, residue-position flags, ampD/ampR/dacB, PDC-high."""
    sample_ids = tuple(sample_ids)
    values: Dict[str, Dict[str, int]] = {}
    order: List[str] = []

    def add(name: str) -> Dict[str, int]:
        column = _column(sample_ids)
        values[name] = column
        order.append(name)
        return column

    alleles: Dict[str, Dict[str, int]] = {}
    positions: Dict[int, Dict[str, int]] = {}
    for row in _amr_rows(sources):
        gene = _cell(row, "gene") or _cell(row, "determinant").split("~")[0]
        if "PDC" not in gene.upper():
            continue
        sample_id = _cell(row, "sample_id")
        _check_sample(sample_id, sample_ids, "amr_determinants.tsv")
        if gene not in alleles:
            alleles[gene] = add(f"L2_{gene}")
        alleles[gene][sample_id] = 1
        position = _residue_position(_cell(row, "variant"))
        if position is not None:
            if position not in positions:
                positions[position] = add(f"L2_PDC_pos_{position}")
            positions[position][sample_id] = 1

    regulator_values: Dict[str, Dict[str, Dict[str, int]]] = {}
    for gene in PDC_REGULATORS:
        loss = add(f"L2_{gene}_loss")
        missense = add(f"L2_{gene}_missense")
        regulator_values[gene] = {"loss": loss, "missense": missense}
        for record in sources.regulator_records:
            if record.gene != gene:
                continue
            _check_sample(record.sample_id, sample_ids, "regulator_variants.tsv")
            if reg_stage.is_disruptive(record):
                loss[record.sample_id] = 1
            elif record.variant_type in MISSENSE_TYPES:
                missense[record.sample_id] = 1

    high = add("L2_PDC_high")
    for sample_id in sample_ids:
        has_allele = any(
            column[sample_id] for column in alleles.values()
        )
        has_loss = any(
            regulator_values[gene]["loss"][sample_id] for gene in PDC_REGULATORS
        )
        high[sample_id] = int(has_allele and has_loss)

    features: List[Feature] = []
    for gene in sorted(alleles):
        features.append(
            Feature(
                name=f"L2_{gene}",
                layer=2,
                source_tokens=_AMR_SOURCE,
                rule=f"1 when a {sources.antibiotic} determinant row reports {gene}",
            )
        )
    for position in sorted(positions):
        features.append(
            Feature(
                name=f"L2_PDC_pos_{position}",
                layer=2,
                source_tokens=_AMR_SOURCE,
                rule=(
                    "1 when a PDC row's amino-acid variant parses to residue "
                    f"{position} (the numbering AMRFinderPlus emits)"
                ),
            )
        )
    disruptive = ", ".join(sorted(reg_stage.DISRUPTIVE_TYPES))
    for gene in PDC_REGULATORS:
        features.append(
            Feature(
                name=f"L2_{gene}_loss",
                layer=2,
                source_tokens=_REGULATOR_SOURCE,
                rule=(
                    f"1 when {gene} carries a disruptive variant class "
                    f"({disruptive}), i.e. loss of function"
                ),
            )
        )
        features.append(
            Feature(
                name=f"L2_{gene}_missense",
                layer=2,
                source_tokens=_REGULATOR_SOURCE,
                rule=(
                    f"1 when {gene} carries an SNV or MNP; never counted as a "
                    "loss and never part of L2_PDC_high"
                ),
            )
        )
    features.append(
        Feature(
            name="L2_PDC_high",
            layer=2,
            source_tokens=_AMR_SOURCE + _REGULATOR_SOURCE,
            rule=(
                "1 when a PDC allele is detected and at least one of "
                f"{', '.join(PDC_REGULATORS)} carries a loss-of-function call"
            ),
        )
    )
    return features, values


# ---------------------------------------------------------------------------
# Layer 3 - oprD absence, loss tiers, individual calls
# ---------------------------------------------------------------------------


def build_layer3(
    sources: Sources, sample_ids: Sequence[str]
) -> Tuple[List[Feature], Dict[str, Dict[str, int]]]:
    """``oprD_absent``, the two LoF tiers, the ``off_any`` composite, and every
    other oprD call as its own feature."""
    sample_ids = tuple(sample_ids)
    values: Dict[str, Dict[str, int]] = {}
    order: List[str] = []

    def add(name: str) -> Dict[str, int]:
        column = _column(sample_ids)
        values[name] = column
        order.append(name)
        return column

    absent = add(f"L3_{OPRD_GENE}_absent")
    tier1 = add(f"L3_{OPRD_GENE}_LoF_tier1")
    tier2 = add(f"L3_{OPRD_GENE}_LoF_tier2")

    individual: Dict[str, Dict[str, int]] = {}
    for record in sources.regulator_records:
        if record.gene != OPRD_GENE:
            continue
        sample_id = record.sample_id
        _check_sample(sample_id, sample_ids, "regulator_variants.tsv")
        if record.variant_type == VariantType.GENE_ABSENCE.value:
            absent[sample_id] = 1
        elif record.variant_type in TIER1_REGULATOR_TYPES:
            tier1[sample_id] = 1
        elif record.variant_type in TIER2_REGULATOR_TYPES:
            tier2[sample_id] = 1
        else:
            name = _sanitised_feature_name(f"L3_{OPRD_GENE}_", record.variant)
            if not name:
                continue
            if name in OPRD_CANONICAL:
                LOGGER.warning(
                    "layers: oprD variant %r sanitises to the canonical "
                    "feature %s; keeping the canonical feature and dropping "
                    "this one",
                    record.variant,
                    name,
                )
                continue
            if name not in individual:
                individual[name] = add(name)
            individual[name][sample_id] = 1

    for row in sources.sv_rows:
        if _cell(row, "affected_gene") != OPRD_GENE:
            continue
        status = _cell(row, "call_status").lower()
        if status == SV_NOT_ASSESSABLE:
            continue
        variant_type = _cell(row, "variant_type")
        if variant_type not in STRUCTURAL_LOSS_TYPES:
            LOGGER.debug(
                "layers: %s %s at %s is not a loss type; not encoded as an "
                "oprD loss",
                status,
                variant_type,
                OPRD_GENE,
            )
            continue
        sample_id = _cell(row, "sample_id")
        _check_sample(sample_id, sample_ids, "structural_variants.tsv")
        if status == SV_CONFIRMED:
            tier1[sample_id] = 1
        elif status == SV_CANDIDATE:
            tier2[sample_id] = 1

    if sources.oprd_structural_rows is not None:
        for row in sources.oprd_structural_rows:
            sample_id = _cell(row, "sample_id")
            _check_sample(sample_id, sample_ids, "oprd_locus/structural_calls.tsv")
            # 'disrupted' means a frameshift or a premature stop and nothing
            # else (papipeline.adapters.oprd_locus), which is tier 1 exactly.
            verdict = _cell(row, "structural_verdict").lower()
            if verdict == "absent":
                absent[sample_id] = 1
            elif verdict == "disrupted":
                tier1[sample_id] = 1

    off_any = add(f"L3_{OPRD_GENE}_off_any")
    for sample_id in sample_ids:
        off_any[sample_id] = int(
            bool(absent[sample_id] or tier1[sample_id] or tier2[sample_id])
        )

    structural_tokens = (
        (_STRUCTURAL_SOURCE,) if sources.oprd_structural_rows is not None else ()
    )
    features = [
        Feature(
            name=f"L3_{OPRD_GENE}_absent",
            layer=3,
            source_tokens=_REGULATOR_SOURCE + structural_tokens,
            rule=(
                "1 when oprD has no intact copy: a GENE_ABSENCE call from the "
                "regulator screen, or an 'absent' verdict in the optional oprD "
                "structural table"
            ),
        ),
        Feature(
            name=f"L3_{OPRD_GENE}_LoF_tier1",
            layer=3,
            source_tokens=_REGULATOR_SOURCE + _SV_SOURCE + structural_tokens,
            rule=(
                "1 when oprD carries a tier-1 lesion: a stop, frameshift or "
                "gene-disruption call, a confirmed structural insertion or "
                "deletion, or a structural verdict of 'disrupted' "
                "(frameshift or premature stop)"
            ),
        ),
        Feature(
            name=f"L3_{OPRD_GENE}_LoF_tier2",
            layer=3,
            source_tokens=_REGULATOR_SOURCE + _SV_SOURCE,
            rule=(
                "1 when oprD carries a tier-2 probable loss: an in-frame or "
                "unresolved coding indel, or a candidate structural variant"
            ),
        ),
        Feature(
            name=f"L3_{OPRD_GENE}_off_any",
            layer=3,
            source_tokens=_REGULATOR_SOURCE + _SV_SOURCE + structural_tokens,
            rule=(
                "1 when L3_oprD_absent, L3_oprD_LoF_tier1 or "
                "L3_oprD_LoF_tier2 is 1; a missense or any other call outside "
                "the tiers never sets it"
            ),
        ),
    ]
    for name in sorted(individual):
        record_variant = name[len(f"L3_{OPRD_GENE}_") :]
        features.append(
            Feature(
                name=name,
                layer=3,
                source_tokens=_REGULATOR_SOURCE,
                rule=(
                    f"1 when the oprD call keyed {record_variant!r} was "
                    "observed; it lies outside the absence and loss tiers, so "
                    "it is an individual feature and never part of "
                    "L3_oprD_off_any"
                ),
            )
        )
    return features, values


# ---------------------------------------------------------------------------
# Layer 4 - efflux regulators and pump-high composites
# ---------------------------------------------------------------------------


def pump_groups(config: PipelineConfig) -> Dict[str, Tuple[str, Tuple[str, ...]]]:
    """``display token -> (token key, member genes)`` for the ticketed genes.

    The grouping is read from ``config/mechanisms.tsv`` ``biological_role``:
    the first ``mex*`` token in the role text is the operon the gene acts on.
    Keys are case-folded so ``MexXY`` and ``mexXY`` join into one group; the
    display token keeps the casing the table first wrote.
    """
    groups: Dict[str, Tuple[str, List[str]]] = {}
    for gene in L4_REGULATOR_GENES:
        spec = config.mechanisms.get(gene)
        role = spec.biological_role if spec is not None else ""
        match = _PUMP_TOKEN.search(role or "")
        token = match.group(0) if match is not None else gene
        key = token.lower()
        if key in groups:
            display, members = groups[key]
            members.append(gene)
        else:
            groups[key] = (token, [gene])
    return {key: (display, tuple(members)) for key, (display, members) in groups.items()}


def build_layer4(
    config: PipelineConfig,
    sources: Sources,
    sample_ids: Sequence[str],
) -> Tuple[List[Feature], Dict[str, Dict[str, int]]]:
    """Loss and missense columns for the five regulators, plus pump-high."""
    sample_ids = tuple(sample_ids)
    values: Dict[str, Dict[str, int]] = {}
    order: List[str] = []

    def add(name: str) -> Dict[str, int]:
        column = _column(sample_ids)
        values[name] = column
        order.append(name)
        return column

    losses: Dict[str, Dict[str, int]] = {}
    for gene in L4_REGULATOR_GENES:
        loss = add(f"L4_{gene}_loss")
        missense = add(f"L4_{gene}_missense")
        losses[gene] = loss
        for record in sources.regulator_records:
            if record.gene != gene:
                continue
            _check_sample(record.sample_id, sample_ids, "regulator_variants.tsv")
            if reg_stage.is_disruptive(record):
                loss[record.sample_id] = 1
            elif record.variant_type in MISSENSE_TYPES:
                missense[record.sample_id] = 1

    disruptive = ", ".join(sorted(reg_stage.DISRUPTIVE_TYPES))
    features: List[Feature] = []
    for gene in L4_REGULATOR_GENES:
        features.append(
            Feature(
                name=f"L4_{gene}_loss",
                layer=4,
                source_tokens=_REGULATOR_SOURCE,
                rule=(
                    f"1 when {gene} carries a disruptive variant class "
                    f"({disruptive}), i.e. loss of function"
                ),
            )
        )
        features.append(
            Feature(
                name=f"L4_{gene}_missense",
                layer=4,
                source_tokens=_REGULATOR_SOURCE,
                rule=(
                    f"1 when {gene} carries an SNV or MNP; a missense is kept "
                    "as its own feature and never counts towards a pump-high "
                    "composite"
                ),
            )
        )

    for key, (display, members) in pump_groups(config).items():
        column = add(f"L4_pump_high_{display}")
        for sample_id in sample_ids:
            column[sample_id] = int(
                any(losses[gene][sample_id] for gene in members)
            )
        features.append(
            Feature(
                name=f"L4_pump_high_{display}",
                layer=4,
                source_tokens=_REGULATOR_SOURCE + _MECHANISM_SOURCE,
                rule=(
                    f"1 when any of {', '.join(members)} carries a loss-of-"
                    "function call; the group is read from "
                    "config/mechanisms.tsv biological_role (first mex* "
                    "token), and a missense never counts"
                ),
            )
        )
    return features, values


# ---------------------------------------------------------------------------
# Layer 5 - target genes and QRDR positions
# ---------------------------------------------------------------------------


def build_layer5(
    sources: Sources, sample_ids: Sequence[str], antibiotic: str
) -> Tuple[List[Feature], Dict[str, Dict[str, int]]]:
    """Any-variant columns for gyrA/parC/ftsI plus the two QRDR features."""
    sample_ids = tuple(sample_ids)
    values: Dict[str, Dict[str, int]] = {}

    def add(name: str) -> Dict[str, int]:
        column = _column(sample_ids)
        values[name] = column
        return column

    gene_columns = {gene: add(f"L5_{gene}") for gene in L5_TARGET_GENES}
    qrdr_pairs: Dict[str, Set[Tuple[str, int]]] = {
        sample_id: set() for sample_id in sample_ids
    }

    for row in _amr_rows(sources, antibiotic):
        gene = _cell(row, "gene")
        if gene not in L5_TARGET_GENES:
            continue
        sample_id = _cell(row, "sample_id")
        _check_sample(sample_id, sample_ids, "amr_determinants.tsv")
        gene_columns[gene][sample_id] = 1
        position = _residue_position(_cell(row, "variant"))
        if gene in QRDR_TARGET_GENES and position is not None:
            qrdr_pairs[sample_id].add((gene, position))

    any_qrdr = add("L5_any_QRDR")
    count = add("L5_QRDR_count")
    for sample_id in sample_ids:
        count[sample_id] = len(qrdr_pairs[sample_id])
        any_qrdr[sample_id] = int(bool(qrdr_pairs[sample_id]))

    features = [
        Feature(
            name=f"L5_{gene}",
            layer=5,
            source_tokens=_AMR_SOURCE,
            rule=(
                f"1 when any {antibiotic} determinant row reports a variant "
                f"in {gene} (PAO1 numbering as emitted by AMRFinderPlus)"
            ),
        )
        for gene in L5_TARGET_GENES
    ]
    features.append(
        Feature(
            name="L5_any_QRDR",
            layer=5,
            source_tokens=_AMR_SOURCE,
            rule=(
                "1 when any gyrA or parC variant parses to a residue position "
                "(PAO1 numbering as emitted by AMRFinderPlus). The repository "
                "defines no QRDR codon range, so no range is invented here: "
                "any quinolone-target substitution counts"
            ),
        )
    )
    features.append(
        Feature(
            name="L5_QRDR_count",
            layer=5,
            value_kind=COUNT,
            source_tokens=_AMR_SOURCE,
            rule=(
                "the number of distinct (gene, residue position) pairs among "
                "parsable gyrA/parC variants; 0 when there are none"
            ),
        )
    )
    return features, values


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def build_layers(
    config: PipelineConfig, sources: Sources, sample_ids: Sequence[str]
) -> LayerBundle:
    """Run all five builders and bundle their features and values."""
    sample_ids = tuple(sample_ids)
    features: List[Feature] = []
    values: Dict[str, Dict[str, int]] = {}

    for builder_features, builder_values in (
        build_layer1(sources, sample_ids),
        build_layer2(sources, sample_ids),
        build_layer3(sources, sample_ids),
        build_layer4(config, sources, sample_ids),
        build_layer5(sources, sample_ids, sources.antibiotic),
    ):
        for feature in builder_features:
            if feature.name in values:
                raise ValueError(
                    f"duplicate feature name {feature.name!r} across layers"
                )
        features.extend(builder_features)
        values.update(builder_values)

    for feature in features:
        column = values[feature.name]
        missing = [sample_id for sample_id in sample_ids if sample_id not in column]
        if missing:
            raise SampleIdError(
                f"feature {feature.name} has no value for "
                f"{len(missing)} isolate(s): {missing}",
                feature=feature.name,
                missing=missing,
            )

    LOGGER.info(
        "layers: built %d features over %d isolates (%s)",
        len(features),
        len(sample_ids),
        ", ".join(
            f"L{layer}={len([f for f in features if f.layer == layer])}"
            for layer in (1, 2, 3, 4, 5)
        ),
    )
    return LayerBundle(sample_ids=sample_ids, features=features, values=values)
