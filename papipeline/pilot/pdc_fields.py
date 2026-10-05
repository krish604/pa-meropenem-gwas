"""Parsers for the two PDC_essential.tsv free-text fields.

Both fields are comma-separated ``key=value`` lists:

* ``AST phenotypes``  -> ``imipenem=R,meropenem=S,ciprofloxacin=I``
* ``AMR genotypes``    -> ``blaPDC-3=COMPLETE,oprD_V359L=POINT``

Parsing rules that protect the analysis from silently wrong numbers:

* a duplicated key inside one record is **ambiguous** and is recorded as
  such, never resolved by taking the first or last value;
* a category outside the configured allowed set is reported as invalid
  rather than coerced;
* R/I/S is never turned into an MIC or a zone diameter, and no MIC is
  present in this source to begin with;
* an unparsable token is skipped and counted, never guessed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from ..logging_utils import get_logger
from ..models import Phenotype

LOGGER = get_logger("pilot.pdc")

#: Categories accepted for a susceptibility call.
ALLOWED_CATEGORIES: Tuple[str, ...] = ("R", "I", "S", "SDD", "ND")

#: Genotype call states seen in PDC_essential.tsv. Recorded, not judged.
KNOWN_CALL_STATES: Tuple[str, ...] = (
    "COMPLETE",
    "PARTIAL",
    "PARTIAL_END_OF_CONTIG",
    "POINT",
    "HMM",
    "MISTRANSLATION",
)

#: Suffixes in PDC gene identifiers that denote a protein-level change
#: rather than a distinct gene. Retained for reporting, never invented.
VARIANT_MARKERS: Tuple[str, ...] = ("fsTer", "Ter", "fs", "del", "ins", "dup")


@dataclass
class AstField:
    """Parsed ``AST phenotypes`` for one isolate."""

    values: Dict[str, str] = field(default_factory=dict)
    ambiguous: Dict[str, List[str]] = field(default_factory=dict)
    invalid: Dict[str, List[str]] = field(default_factory=dict)
    unparsed: int = 0

    def get(self, antibiotic: str) -> Optional[Phenotype]:
        """Return the category for ``antibiotic``, or ``None`` if unavailable.

        An ambiguous or invalid entry returns ``None``: the isolate is then
        counted as lacking a usable value rather than being assigned one.
        """
        key = antibiotic.strip().lower()
        if key in self.ambiguous or key in self.invalid:
            return None
        raw = self.values.get(key)
        if raw is None:
            return None
        return Phenotype(raw)

    def availability(self, antibiotic: str) -> str:
        """Why a value is or is not available, for the QC report."""
        key = antibiotic.strip().lower()
        if key in self.ambiguous:
            return "AMBIGUOUS_DUPLICATE_KEY"
        if key in self.invalid:
            return "INVALID_CATEGORY"
        if key in self.values:
            return "PRESENT"
        return "ABSENT"


def parse_ast_field(text: Optional[str]) -> AstField:
    """Parse a PDC ``AST phenotypes`` value."""
    result = AstField()
    if not text:
        return result
    seen: Dict[str, List[str]] = {}
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            result.unparsed += 1
            continue
        antibiotic, _, category = token.partition("=")
        antibiotic = antibiotic.strip().lower()
        category = category.strip().upper()
        if not antibiotic or not category:
            result.unparsed += 1
            continue
        seen.setdefault(antibiotic, []).append(category)

    for antibiotic, categories in seen.items():
        distinct = sorted(set(categories))
        if len(distinct) > 1:
            # Conflicting values for one antibiotic: ambiguous, not resolved.
            result.ambiguous[antibiotic] = categories
            LOGGER.warning(
                "Ambiguous AST value for %s: %s", antibiotic, ",".join(categories)
            )
        elif distinct[0] not in ALLOWED_CATEGORIES:
            result.invalid[antibiotic] = categories
            LOGGER.warning(
                "AST category %r for %s is outside the allowed set %s",
                distinct[0],
                antibiotic,
                ",".join(ALLOWED_CATEGORIES),
            )
        else:
            result.values[antibiotic] = distinct[0]
    return result


@dataclass(frozen=True)
class GenotypeCall:
    """One ``gene=STATE`` entry from PDC."""

    gene: str
    state: str
    is_variant: bool

    @property
    def base_gene(self) -> str:
        """Gene symbol with any variant suffix removed.

        ``oprD_V359L`` -> ``oprD``; ``blaOXA-486`` -> ``blaOXA-486``
        (the numeric suffix is part of the allele name, not a variant).
        """
        if "_" not in self.gene:
            return self.gene
        head, _, tail = self.gene.partition("_")
        if any(marker in tail for marker in VARIANT_MARKERS) or any(
            ch.isdigit() for ch in tail
        ):
            return head
        return self.gene


@dataclass
class GenotypeField:
    """Parsed ``AMR genotypes`` for one isolate."""

    calls: List[GenotypeCall] = field(default_factory=list)
    unparsed: int = 0
    unknown_states: Set[str] = field(default_factory=set)

    @property
    def genes(self) -> List[str]:
        return [c.gene for c in self.calls]

    def base_genes(self) -> Set[str]:
        return {c.base_gene for c in self.calls}

    def variant_calls(self) -> List[GenotypeCall]:
        return [c for c in self.calls if c.is_variant]

    def state_of(self, gene: str) -> Optional[str]:
        for call in self.calls:
            if call.gene == gene:
                return call.state
        return None


def _looks_like_variant(gene: str) -> bool:
    """Whether a PDC gene token denotes a specific variant.

    ``oprD_V359L`` and ``mexR_I24AfsTer94`` are variants; ``blaOXA-486``
    and ``mexA`` are gene identifiers.
    """
    if "_" not in gene:
        return False
    head, _, tail = gene.partition("_")
    if not head or not tail:
        return False
    if any(marker in tail for marker in VARIANT_MARKERS):
        return True
    # Amino-acid style change codes, e.g. V359L, G307D, S278P.
    return bool(len(tail) >= 3 and tail[0].isalpha() and tail[-1].isalpha())


def parse_genotype_field(text: Optional[str]) -> GenotypeField:
    """Parse a PDC ``AMR genotypes`` value."""
    result = GenotypeField()
    if not text:
        return result
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            result.unparsed += 1
            continue
        gene, _, state = token.partition("=")
        gene = gene.strip()
        state = state.strip().upper()
        if not gene or not state:
            result.unparsed += 1
            continue
        if state not in KNOWN_CALL_STATES:
            result.unknown_states.add(state)
        result.calls.append(
            GenotypeCall(gene=gene, state=state, is_variant=_looks_like_variant(gene))
        )
    return result


def antibiotic_label(antibiotic: str) -> str:
    """Display label used in output filenames and slide titles."""
    return antibiotic.strip().capitalize()
