"""Resolving PAO1-coordinate variants to regulator genes.

**Why reference coordinates.** Stage 6 calls variants against PAO1, so its
positions are reference coordinates. The cohort's own assemblies carry
*Bakta's* annotations, and those are assembly-local contigs - `JFJU01000001.1`
and thousands more - on draft scaffolds whose coordinate system is unrelated to
PAO1's. Intersecting the two directly would silently attribute variants to
unrelated genes. The mapping therefore runs through the pinned PAO1 GFF, which is
in the same frame as the variants.

**What a hit means.** A variant *in a screened locus*, reported with the
mechanism configured for that locus. Never phenotypic resistance: establishing
that needs stage 12, robustness needs stage 13.

**`GENE_ABSENCE` is not implemented here.** Absence is a property of the
assembly, not of a variant call, and it comes from the gene presence/absence
stage. Every locus lists it in `screen_for`, so the filter will pass such a class
through when that stage exists; until then it is simply never emitted, and
`screen_for` is not a promise this module fulfils alone.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..errors import PipelineError

#: Codons that terminate translation.
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class GeneInterval:
    """A gene's genomic span on one reference contig."""

    locus_tag: str
    contig: str
    start: int
    end: int
    strand: str
    name: str = ""
    product: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def contains(self, contig: str, position: int) -> bool:
        return contig == self.contig and self.start <= position <= self.end


@dataclass(frozen=True)
class VariantCall:
    """One row of stage 6's per-isolate table, in reference coordinates."""

    sample_id: str
    contig: str
    position: int
    reference: str
    alternate: str


def load_gene_intervals(
    gff_path: Path, tags: Iterable[str]
) -> Dict[str, GeneInterval]:
    """Index the named locus tags from the pinned GFF.

    Only `gene` features are read. CDS rows share the same `locus_tag`, so
    counting both makes a locus look like it has two features - which is how the
    locus resolver caught itself once already.

    Raises if the GFF is absent, or yields nothing for the tags asked for: a
    silent empty index reports every isolate as carrying no regulator variant,
    which is a fabricated negative on the study's central mechanisms.
    """
    path = Path(gff_path)
    if not path.is_file():
        raise PipelineError(
            f"PAO1 GFF not found at {path}. Stage 6's variants are in PAO1 "
            "coordinates, so resolving them to genes needs the pinned reference "
            "annotation. There is no fallback: the cohort's own assemblies are "
            "in a different coordinate system, and intersecting them would "
            "attribute variants to unrelated genes.",
            path=str(path),
        )
    wanted = {t.strip() for t in tags if t and t.strip()}
    if not wanted:
        return {}
    found: Dict[str, GeneInterval] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        attrs = dict(
            part.split("=", 1) for part in fields[8].split(";") if "=" in part
        )
        tag = attrs.get("locus_tag", "")
        if tag not in wanted or tag in found:
            continue
        found[tag] = GeneInterval(
            locus_tag=tag,
            contig=fields[0],
            start=int(fields[3]),
            end=int(fields[4]),
            strand=fields[6],
            name=attrs.get("Name", ""),
            product=attrs.get("product", ""),
        )
    missing = sorted(wanted - set(found))
    if missing:
        raise PipelineError(
            f"{len(missing)} configured locus tag(s) are absent from the pinned "
            f"GFF: {missing}. config/regulators.tsv and the reference annotation "
            "disagree; regenerate the table with scripts/resolve_locus_tags.py.",
            missing=",".join(missing),
        )
    return found


def gene_at(
    intervals: Sequence[GeneInterval], contig: str, position: int
) -> Optional[GeneInterval]:
    """The gene containing a position, or None."""
    for interval in intervals:
        if interval.contains(contig, position):
            return interval
    return None


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def codon_at(
    reference_sequence: str, interval: GeneInterval, position: int
) -> Optional[str]:
    """The reference codon covering `position`, or None if it runs off the end.

    Position is 1-based inclusive; the sequence is 0-based. On the reverse strand
    the codon is the reverse complement, so a stop is a stop regardless of which
    strand the gene is transcribed from - which is the whole reason this reads
    the strand at all.
    """
    offset = position - interval.start
    if offset < 0 or offset > interval.length - 1:
        return None
    if interval.strand == "-":
        # On the minus strand the codon reads right-to-left in plus-strand
        # coordinates, so walk back from `position` toward higher coordinates.
        codon_start = offset - (offset % 3)
        start = interval.start + codon_start
        if start + 3 > interval.end + 1:
            return None
        plus = reference_sequence[start - 1: start + 2]
        return reverse_complement(plus) if len(plus) == 3 else None
    start = position - (offset % 3)
    if start + 3 > interval.end + 1:
        return None
    chunk = reference_sequence[start - 1: start + 2]
    return chunk if len(chunk) == 3 else None


def classify(
    call: VariantCall,
    reference_sequence: str = "",
    interval: Optional[GeneInterval] = None,
) -> str:
    """The variant class, derived from the allele strings.

    `PREMATURE_STOP` needs reference context to be honest, so it is reported
    only when a reference sequence and the gene's span are both available.
    Without them a substitution is reported as `SNV`, which **understates** it -
    the deliberate direction, because a fabricated `PREMATURE_STOP` would be a
    mechanism claim the data does not support.
    """
    ref = (call.reference or "").upper()
    alt = (call.alternate or "").upper()
    if not ref or not alt or ref == alt:
        return "UNKNOWN"
    if len(ref) == 1 and len(alt) == 1:
        if reference_sequence and interval is not None:
            codon = codon_at(reference_sequence, interval, call.position)
            if codon:
                index = codon_index(interval, call.position)
                # `codon` is written on the coding strand; `alt` is a VCF
                # allele, which on the minus strand is the complement. Mirroring
                # the index without complementing the base compares a base from
                # one strand against a codon from the other, so a genuine
                # stop-gaining substitution on any minus-strand locus - oprD,
                # mexR, ampD, ampR, four of the ten screenable loci - could
                # never be recognised. Verified by brute force over all 444
                # oprD codons before and after this line.
                coding_alt = alt if interval.strand == "+" else alt.translate(
                    _COMPLEMENT
                )
                mutated = codon[:index] + coding_alt + codon[index + 1:]
                if mutated in STOP_CODONS:
                    return "PREMATURE_STOP"
        return "SNV"
    delta = len(alt) - len(ref)
    return "INDEL" if delta % 3 == 0 else "FRAMESHIFT"


def codon_index(interval: GeneInterval, position: int) -> int:
    """Which base of its codon ``position`` is, 0-based from the codon start.

    Strand-aware, and this is the whole function. :func:`codon_at` reports the
    codon already written 5'->3' *on the coding strand* - on the minus strand
    that is the reverse complement of the plus-strand block - so the base to
    substitute has to be indexed within that codon, not within the gene.

    It used to be `call.position - interval.start`, the offset from the gene
    start, gated on `0 <= index < 3`. That gate is only true for the first three
    bases of a gene, so for every variant past them the stop check was skipped
    and the call degraded to a bare ``SNV``. Proven on the pinned reference:
    ``PA0958`` (oprD, minus strand), coding codon 296 is ``AAG``; a
    ``NC_002516.2:1044429 T>A`` makes it ``TAG``, a stop, and ``classify``
    returned ``'SNV'`` because the offset was 446.

    The existing unit test only exercised codon 1, where the wrong expression
    happens to land on the right answer, so it agreed with itself and hid this.

    On the minus strand ``codon_at`` reads the plus-strand block
    ``offset - offset % 3`` and reverse-complements it, so plus-strand index
    ``offset % 3`` becomes coding index ``2 - (offset % 3)``.
    """
    offset = position - interval.start
    within = offset % 3
    return within if interval.strand == "+" else 2 - within


def codon_number(interval: GeneInterval, position: int) -> int:
    """Which codon ``position`` falls in, 1-based from the gene's 5' end.

    Strand-aware, for the same reason :func:`codon_index` is: on the minus
    strand the gene is transcribed from its *high* coordinate end, so counting
    codons from ``interval.start`` counts them backwards.

    ``offset // 3 + 1`` is right on the plus strand and exactly inverted on the
    minus one. Measured on the pinned reference, ``PA0958`` (oprD, minus strand,
    1,332 nt): its first codon ``ATG`` sits at 1045314 with offset 1331 and
    ``offset // 3 + 1`` calls it codon 444; its 368th codon ``GAT`` sits at
    1044213 with offset 230 and the same expression calls it codon 77.

    :func:`codon_at` does not need this - it reports the codon's bases, and the
    ``offset - offset % 3`` block it reads is the right block on either strand -
    but anything that *names* a codon has to agree with the transcript's
    direction, or it labels the truncation site with the wrong residue number.
    """
    offset = position - interval.start
    if interval.strand == "-":
        offset = interval.length - 1 - offset
    return offset // 3 + 1


def variant_label(call: VariantCall) -> str:
    return f"{call.reference}>{call.alternate}"


def assign(
    call: VariantCall,
    intervals: Sequence[GeneInterval],
) -> List[Tuple[str, GeneInterval]]:
    """Every screened locus containing this variant.

    A position inside an overlap is reported for each gene covering it. Dropping
    one would under-report, and these are the study's central mechanisms.
    """
    hits: List[Tuple[str, GeneInterval]] = []
    for interval in intervals:
        if interval.contains(call.contig, call.position):
            hits.append((interval.name or interval.locus_tag, interval))
    return hits