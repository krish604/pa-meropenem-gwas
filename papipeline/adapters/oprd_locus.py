"""Resolving an isolate's oprD locus to the PAO1 PA0958 orthologue, or refusing.

**Why protein identity and not the annotation's symbol.** Stage 2 reports
Bakta's ``gene=`` symbol, and ``papipeline/viz.py`` decides ``oprD`` carriage with
a string test on it. That test is wrong on real data. Across the ten smoke
isolates, Bakta writes ``gene=oprD`` on 36 CDS features and 34 of them are the
OprD/OprP/OprQ porin family at 33.7-36.8% identity to PAO1 PA0958; only 2 are
the orthologue. Worse, 8 of the 10 isolates carry the true locus as a CDS Bakta
leaves **unlabelled**, so a symbol test misses the locus where it is present and
finds it where it is a paralog. Coordinates cannot substitute either: the isolate
GFFs are assembly-local ``contig_N`` and carry no PAO1 locus tags.

**Why the fragments must be merged before anything is decided.** Five of the ten
isolates have the locus split into two CDS by a draft-assembly gap, the two pieces
tiling the reference protein in one orientation. Measured on the ten isolates,
the per-CDS coverages are 38.5/53.2, 52.7/47.8, 82.7/16.2, 38.5/34.2 and 33.8/33.1
percent, but the *unions* are 99.8, 100.0, 99.1, 72.9 and 67.0. A rule that looks
at one best hit calls those first three absent. A rule that calls a tie ambiguous
refuses them. Both are wrong, so this module merges collinear hits within one
``(contig, strand)`` and decides on the union.

**What a refusal means, which is the point of the module.** Every outcome other
than ``resolved`` is a *refusal with a named reason*, and none of them is
``absent``. An isolate whose paralogs are the only thing found has not been shown
to lack oprD; it has been shown not to have been assessed. Collapsing the two
would manufacture the study's central negative. ``papipeline/viz.py`` already
carries ``not_assessed`` as a distinct state for exactly this reason.

**The thresholds are read off the measured distribution, not chosen.** Over all
ten isolates at ``E <= 10`` the identity distribution is bimodal with nothing
between the modes: the highest non-orthologue hit is 45.16% and the lowest
orthologue fragment is 88.61%, a 43.45-point gap. ``MIN_IDENTITY_PCT = 70.0``
sits inside it, so the exact value is not load-bearing here - but a floor that
high would discard a genuinely divergent OprD, and that limitation is recorded
rather than hidden. The union-coverage distribution has the same shape, with
nothing between 72.91% and 93.91%; ``MIN_COVERAGE_PCT = 80.0`` sits inside that
gap and is the load-bearing threshold, because it is what separates "the locus is
recognisably the locus" from "a fragment of it was recovered".

``resolve`` is a pure function over parsed blastp tabular output, so the decision
logic is testable without running blast. :func:`run_blastp` is the only impure
part and takes an injectable ``runner``, as the virulence adapter does.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

from ..errors import PipelineError, ToolExecutionError

#: Lowest identity at which a hit is considered a candidate orthologue fragment.
#: Chosen inside the measured 45.16-88.61% gap; see the module docstring.
MIN_IDENTITY_PCT = 70.0

#: Fraction of the reference protein a single cluster must reconstruct before the
#: locus is called resolved. Chosen inside the measured 72.91-93.91% gap.
MIN_COVERAGE_PCT = 80.0

#: blastp e-value for the search. Loose on purpose: a paralog sitting at E=1e-30
#: is evidence about which genes the family has, and this module is asked to
#: *refuse* rather than to filter, so the search must not hide anything.
EVALUE = 10.0

#: blastp caps reported targets per query; 200 is far above the 24-27 HSPs any
#: single isolate produced. `0` is ILLEGAL in blastp 2.17.0 (`expected >=1`).
MAX_TARGET_SEQS = 200

#: blastp tabular columns. Query start/end are mandatory, not optional: without
#: them the fragments of one split locus cannot be told from two real genes.
OUTFMT_FIELDS = (
    "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "bitscore",
)
OUTFMT = "6 " + " ".join(OUTFMT_FIELDS)


class Verdict:
    """The outcomes. Only ``RESOLVED`` asserts the locus is present."""

    RESOLVED = "resolved"
    NO_HIT = "refused:no_orthologous_hit"
    INSUFFICIENT = "refused:insufficient_coverage"
    AMBIGUOUS = "refused:ambiguous_locus"
    AMBIGUOUS_SECOND = "refused:ambiguous_second_locus"

    #: Every refusal. None of these means the locus is absent.
    REFUSALS = frozenset({NO_HIT, INSUFFICIENT, AMBIGUOUS, AMBIGUOUS_SECOND})


@dataclass(frozen=True)
class Hit:
    """One blastp HSP: reference protein against one isolate CDS."""

    subject_id: str
    contig: str
    strand: str
    identity_pct: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    bitscore: float

    @property
    def query_coverage_pct(self) -> float:
        return 100.0 * self.alignment_length / self.query_length

    @property
    def query_low(self) -> int:
        return min(self.query_start, self.query_end)

    @property
    def query_high(self) -> int:
        return max(self.query_start, self.query_end)


@dataclass(frozen=True)
class Cluster:
    """Hits on one contig in one orientation: candidate fragments of one locus."""

    contig: str
    strand: str
    covered: Tuple[int, ...]
    hits: Tuple[Hit, ...]

    @property
    def coverage_pct(self) -> float:
        return 100.0 * len(self.covered) / self.hits[0].query_length

    @property
    def min_identity_pct(self) -> float:
        return min(h.identity_pct for h in self.hits)

    @property
    def gene_ids(self) -> Tuple[str, ...]:
        return tuple(h.subject_id for h in self.hits)


@dataclass(frozen=True)
class LocusResolution:
    """What the resolver concluded about one isolate's oprD locus."""

    sample_id: str
    verdict: str
    reason: str
    gene_ids: Tuple[str, ...] = ()
    contig: str = ""
    strand: str = ""
    coverage_pct: float = 0.0
    identity_pct: Optional[float] = None
    n_total_hits: int = 0
    n_candidate_hits: int = 0
    n_paralog_hits: int = 0
    clusters: Tuple[Cluster, ...] = field(default=(), repr=False)

    @property
    def is_resolved(self) -> bool:
        return self.verdict == Verdict.RESOLVED

    @property
    def is_refused(self) -> bool:
        return self.verdict in Verdict.REFUSALS


def parse_blast_table(text: str, query_length: int) -> List[Hit]:
    """Parse ``OUTFMT`` tabular output into hits.

    The subject's contig and strand are not in the tabular output - blastp is
    handed a FASTA, not a GFF - so they arrive as keyword arguments populated
    from the isolate's own GFF3 by :func:`load_cds_locations`. They default to
    ``""``, which makes every hit one cluster; a caller that has a GFF should pass
    it, because the cluster split is what separates a split locus from two loci.
    """
    hits: List[Hit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(OUTFMT_FIELDS):
            raise PipelineError(
                f"blastp output row has {len(fields)} fields, expected "
                f"{len(OUTFMT_FIELDS)} ({' '.join(OUTFMT_FIELDS)}). A different "
                f"-outfmt would make coverage and the fragment merge silently "
                f"wrong, so this is refused rather than parsed. Row: {line[:200]}",
                fields=len(fields),
            )
        try:
            hits.append(
                Hit(
                    subject_id=fields[0],
                    contig="",
                    strand="",
                    identity_pct=float(fields[1]),
                    alignment_length=int(fields[2]),
                    query_length=query_length,
                    subject_length=int(fields[4]),
                    query_start=int(fields[5]),
                    query_end=int(fields[6]),
                    bitscore=float(fields[9]),
                )
            )
        except ValueError as exc:
            raise PipelineError(
                f"blastp output row is not numeric where it must be: {exc}. "
                f"Row: {line[:200]}",
                row=line[:200],
            ) from exc
    return hits


def load_cds_locations(gff_path: Path) -> Dict[str, Tuple[str, str]]:
    """Map each CDS ``locus_tag`` to its ``(contig, strand)`` in the isolate GFF3.

    Only ``CDS`` features are read, and a repeated ``locus_tag`` keeps the first,
    matching :func:`papipeline.adapters.gff.load_gene_intervals`'s handling of the
    shared-tag problem.
    """
    path = Path(gff_path)
    if not path.is_file():
        raise PipelineError(
            f"Isolate GFF3 not found at {path}. Without contig and strand the "
            f"hits cannot be grouped into clusters, and the two fragments of one "
            f"split locus become indistinguishable from two separate loci.",
            path=str(path),
        )
    out: Dict[str, Tuple[str, str]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        attrs = dict(
            part.split("=", 1) for part in fields[8].split(";") if "=" in part
        )
        tag = attrs.get("locus_tag") or attrs.get("ID")
        if tag and tag not in out:
            out[tag] = (fields[0], fields[6])
    return out


def attach_locations(
    hits: Sequence[Hit], locations: Mapping[str, Tuple[str, str]]
) -> List[Hit]:
    """Re-issue hits carrying their own contig and strand."""
    return [
        Hit(
            subject_id=h.subject_id,
            contig=locations.get(h.subject_id, ("", ""))[0],
            strand=locations.get(h.subject_id, ("", ""))[1],
            identity_pct=h.identity_pct,
            alignment_length=h.alignment_length,
            query_length=h.query_length,
            subject_length=h.subject_length,
            query_start=h.query_start,
            query_end=h.query_end,
            bitscore=h.bitscore,
        )
        for h in hits
    ]


def cluster_hits(hits: Sequence[Hit]) -> List[Cluster]:
    """Group candidate hits into one cluster per ``(contig, strand)``.

    Coverage within a cluster is the **union** of the query intervals its hits
    span, not the best single hit and not their sum. The union is what makes a
    split locus read as one locus: two fragments covering 38.5% and 34.2% cover
    72.9% between them and neither covers it alone, while two genuinely separate
    genes sit on different contigs and never merge however well they match.
    """
    grouped: Dict[Tuple[str, str], List[Hit]] = {}
    for hit in hits:
        grouped.setdefault((hit.contig, hit.strand), []).append(hit)

    clusters: List[Cluster] = []
    for (contig, strand), members in grouped.items():
        covered: set = set()
        for hit in members:
            covered.update(range(hit.query_low, hit.query_high + 1))
        clusters.append(
            Cluster(
                contig=contig,
                strand=strand,
                covered=tuple(sorted(covered)),
                hits=tuple(sorted(members, key=lambda h: h.query_low)),
            )
        )
    clusters.sort(key=lambda c: -c.coverage_pct)
    return clusters


def resolve(
    sample_id: str,
    hits: Sequence[Hit],
    *,
    min_identity_pct: float = MIN_IDENTITY_PCT,
    min_coverage_pct: float = MIN_COVERAGE_PCT,
) -> LocusResolution:
    """Decide one isolate's oprD locus from its blastp hits, or refuse.

    The order of the checks is the order of what can be asserted. A cluster that
    reconstructs the reference protein is the only thing that supports "the
    locus is present", so the ambiguity checks come first: two such clusters is a
    tie with no honest tiebreak, and one such cluster plus a substantial second
    one is a call about a paralog, not about oprD. Only then is coverage asked.

    Returns a :class:`LocusResolution` whose ``verdict`` is
    :data:`Verdict.RESOLVED` or one of the named refusals. It never raises for a
    scientific outcome - only for malformed input, which is a different thing.
    """
    candidates = [h for h in hits if h.identity_pct >= min_identity_pct]
    clusters = cluster_hits(candidates)
    complete = [c for c in clusters if c.coverage_pct >= min_coverage_pct]

    def outcome(verdict: str, reason: str, **kwargs) -> LocusResolution:
        return LocusResolution(
            sample_id=sample_id,
            verdict=verdict,
            reason=reason,
            n_total_hits=len(hits),
            n_candidate_hits=len(candidates),
            n_paralog_hits=len(hits) - len(candidates),
            clusters=tuple(clusters),
            **kwargs,
        )

    if not candidates:
        return outcome(
            Verdict.NO_HIT,
            f"No CDS in {sample_id} reaches {min_identity_pct}% identity to the "
            f"PAO1 oprD (PA0958) protein across {len(hits)} reported HSPs. The "
            f"highest identity seen was "
            f"{max((h.identity_pct for h in hits), default=0.0):.2f}%, which is "
            f"the OprD/OprP/OprQ porin family. This is NOT evidence that oprD is "
            f"absent: it is evidence the locus was not located.",
        )

    if len(complete) > 1:
        return outcome(
            Verdict.AMBIGUOUS,
            f"{len(complete)} separate loci in {sample_id} each reconstruct at "
            f"least {min_coverage_pct}% of the PAO1 oprD protein "
            f"({', '.join(f'{c.contig}({c.strand}) {c.coverage_pct:.1f}%' for c in complete)}). "
            f"Two independent loci matching the reference equally well cannot be "
            f"told apart from this evidence, so neither is reported as oprD.",
        )

    if complete:
        best = complete[0]
        others = [c for c in clusters if c is not best]
        substantial = [
            c for c in others if c.coverage_pct >= min_coverage_pct / 2.0
        ]
        if substantial:
            return outcome(
                Verdict.AMBIGUOUS_SECOND,
                f"{sample_id} has one near-complete oprD locus on "
                f"{best.contig}({best.strand}) at {best.coverage_pct:.1f}% "
                f"coverage, and a second locus on "
                f"{substantial[0].contig}({substantial[0].strand}) at "
                f"{substantial[0].coverage_pct:.1f}% coverage. A second copy this "
                f"close to a full match makes this a paralog assignment problem, "
                f"not an oprD carriage call.",
            )
        return outcome(
            Verdict.RESOLVED,
            f"oprD resolved to {', '.join(best.gene_ids)} on "
            f"{best.contig}({best.strand}): {best.coverage_pct:.1f}% of the PAO1 "
            f"PA0958 protein at {best.min_identity_pct:.2f}% identity or better.",
            gene_ids=best.gene_ids,
            contig=best.contig,
            strand=best.strand,
            coverage_pct=best.coverage_pct,
            identity_pct=best.min_identity_pct,
        )

    top = clusters[0]
    return outcome(
        Verdict.INSUFFICIENT,
        f"{sample_id} carries {len(candidates)} hit(s) at or above "
        f"{min_identity_pct}% identity to PAO1 oprD, but the best collinear "
        f"cluster ({top.contig}({top.strand}), {', '.join(top.gene_ids)}) "
        f"reconstructs only {top.coverage_pct:.1f}% of the reference protein, "
        f"below the {min_coverage_pct}% floor. The locus is neither shown present "
        f"nor shown absent; a draft assembly that recovers part of a gene cannot "
        f"settle it.",
        gene_ids=top.gene_ids,
        contig=top.contig,
        strand=top.strand,
        coverage_pct=top.coverage_pct,
        identity_pct=top.min_identity_pct,
    )


def reference_cds_nucleotides(
    gff_path: Path, fasta_path: Path, locus_tag: str
) -> str:
    """The PAO1 reference **CDS nucleotides** for ``locus_tag``, from the pin.

    The sibling of :func:`reference_protein` and derived from the same two
    pinned files, for one reason: the structural screen aligns the isolate's
    nucleotides against the reference's, so a reference read twice - once as
    protein for the blastp resolver and once as DNA for the tblastn structural
    path - could drift, and nothing would catch it. Both now come from
    :func:`_reference_cds_from_pin`.

    The returned string INCLUDES the terminal stop codon: PAO1 PA0958 is 1332 nt
    = 444 codons = 443 residues plus its own ``TAA``. The structural path needs
    the stop, because "is this stop the reference's own terminator" is the
    question RULING 2 turns on.
    """
    return _reference_cds_from_pin(gff_path, fasta_path, locus_tag)


def _reference_cds_from_pin(
    gff_path: Path, fasta_path: Path, locus_tag: str
) -> str:
    """Read the reference CDS for ``locus_tag`` off the pinned GFF and FASTA.

    Every refusal :func:`reference_protein` documents lives here, because the
    structural path inherits all of them: a partial, multi-segment, frameshifted
    or pseudogene pin would silently change what "100% coverage" means and what
    the reference's terminator is.
    """
    gff = Path(gff_path)
    fasta = Path(fasta_path)
    for path, what in ((gff, "PAO1 GFF"), (fasta, "PAO1 FASTA")):
        if not path.is_file():
            raise PipelineError(
                f"{what} not found at {path}. The oprD query protein is derived "
                f"from the pinned reference, and there is no substitute: an "
                f"unpinned query makes every identity number unreproducible.",
                path=str(path),
            )

    contig: Optional[str] = None
    region: Optional[Tuple[int, int, str]] = None
    cds: List[Tuple[int, int, str, int]] = []
    for line in gff.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        attrs = dict(
            part.split("=", 1) for part in fields[8].split(";") if "=" in part
        )
        if attrs.get("locus_tag") != locus_tag:
            continue
        if fields[2] == "gene" and contig is None:
            contig = fields[0]
            region = (int(fields[3]), int(fields[4]), fields[6])
        elif fields[2] == "CDS":
            cds.append((int(fields[3]), int(fields[4]), fields[6],
                        int(attrs.get("transl_table", "11"))))

    if contig is None or region is None:
        raise PipelineError(
            f"No gene feature for {locus_tag} in {gff}. The pinned reference and "
            f"the locus screen disagree; regenerate with "
            f"scripts/resolve_locus_tags.py rather than guessing the tag.",
            locus_tag=locus_tag,
        )
    if len(cds) != 1:
        raise PipelineError(
            f"{locus_tag} has {len(cds)} CDS features in {gff}, expected 1. A "
            f"multi-segment CDS would make the translated protein a chimera and "
            f"every coverage figure against it meaningless.",
            locus_tag=locus_tag,
            n_cds=len(cds),
        )
    start, end, strand, table = cds[0]
    if (end - start + 1) % 3 != 0:
        raise PipelineError(
            f"{locus_tag} CDS is {end - start + 1} bp, not a multiple of three. "
            f"Translating it would produce a frameshifted protein.",
            locus_tag=locus_tag,
            length=end - start + 1,
        )
    if table != 11:
        raise PipelineError(
            f"{locus_tag} CDS declares transl_table={table}; only table 11 is "
            f"implemented. Bacterial start codons differ between tables, so "
            f"guessing would change the protein.",
            locus_tag=locus_tag,
            transl_table=table,
        )

    sequence = _read_contig(fasta, contig)
    cds_seq = sequence[start - 1:end]
    if strand == "-":
        cds_seq = _reverse_complement(cds_seq)
    return str(cds_seq).upper()


def reference_protein(gff_path: Path, fasta_path: Path, locus_tag: str) -> str:
    """The PAO1 reference protein for ``locus_tag``, translated from the pin.

    Derived from the two files already pinned rather than from a separately
    downloaded proteome. A third artifact could drift from the first two and
    there would be nothing to catch it; these two are the reference the whole
    pipeline calls variants against, so deriving from them keeps one frame.

    The returned string is the protein **without** its terminal stop, so it is
    443 residues for PA0958 and a perfect match scores 100% coverage.

    Raises if the locus is absent, is not a single CDS, has an internal stop, or
    its length is not a multiple of three: a partial, frameshifted or
    pseudogene pin would silently change what "100% coverage" means.
    """
    cds_seq = _reference_cds_from_pin(gff_path, fasta_path, locus_tag)
    protein = _translate(cds_seq)
    locus_tag = str(locus_tag)
    if protein.startswith("M") is False:
        raise PipelineError(
            f"{locus_tag} does not start with methionine (got "
            f"{protein[:1]!r}). Either the strand or the table is wrong.",
            locus_tag=locus_tag,
        )
    # The internal-stop check comes FIRST. A pseudogene typically ends in a stop
    # as well, so testing `endswith("*")` first would strip that stop and never
    # look at the one in the middle - which is the one that matters.
    if "*" in protein[:-1]:
        raise PipelineError(
            f"{locus_tag} translates with an internal stop at position "
            f"{protein.index('*') + 1}, so it is not one protein. A pseudogene or "
            f"a frameshifted pin would make every identity number meaningless.",
            locus_tag=locus_tag,
        )
    # The terminal stop is dropped HERE ONLY, because this is the blastp
    # resolver's query and its coverage is measured against the reference
    # PROTEIN: 443 residues is the whole protein, and keeping the stop would cap
    # a perfect match at 99.77% and quietly move the presence coverage floor.
    #
    # It is NOT the tblastn query. The structural screen measures against the
    # reference CDS, whose codon count INCLUDES the stop, so its query is
    # :func:`tblastn_query_protein` and is written by :func:`write_tblastn_query`.
    # Handing this function's 443-residue string to tblastn is what made
    # `stop_is_reference_terminator` unsatisfiable in production: `qend <= qlen`
    # can never exceed the query's own length, so an alignment of a 443-residue
    # query cannot reach reference codon 444. Both counts now come from
    # :func:`reference_codon_count`, so they cannot drift apart again.
    if protein.endswith("*"):
        protein = protein[:-1]
    return protein


def reference_codon_count(reference_cds: str) -> int:
    """The reference's codon count, **terminal stop included**. The definition.

    PAO1 PA0958 is 1332 nt = 444 codons = 443 residues plus its own ``TAA``, so
    this is 444. Two callers need that number and they must not compute it
    separately, because they are compared against each other:

    * :func:`measure_locus_structure` compares the alignment's high query
      coordinate against it to decide whether the search reached the
      reference's own terminator, and
    * :func:`tblastn_query_protein` builds the query whose length must be
      exactly that, or the first comparison is unsatisfiable for every isolate
      and an intact oprD becomes unreachable.

    One definition, so the query cannot be shorter than the reference it is
    measured against. Callers pass the reference CDS nucleotides **with** the
    stop - :func:`reference_cds_nucleotides` - never a protein.
    """
    return len(reference_cds) // 3


def tblastn_query_protein(reference_cds: str) -> str:
    """The tblastn query: the reference protein **with** its terminal stop.

    444 characters for PA0958, the last of them ``*``, which is what the stored
    round-9 tables were produced from and what
    :func:`measure_locus_structure`'s ``ref_high == reference_codons`` test
    requires. :func:`reference_protein` is the *other* query - the blastp
    resolver's - and is one character shorter by design; see its comment.

    Refuses if the reference is not a whole number of codons, or if its
    translation does not end in a stop. A query missing its terminator is not a
    smaller query, it is one that can never satisfy the terminator test, and the
    failure would otherwise show up only as every isolate reading
    ``disrupted``.
    """
    if len(reference_cds) % 3:
        raise PipelineError(
            f"The reference CDS is {len(reference_cds)} nt, which is not a whole "
            f"number of codons. A reference whose reading frame does not close "
            f"cannot be searched against, and "
            f"{len(reference_cds) // 3} codons of it would be measured as though "
            f"they were the whole locus.",
            length=len(reference_cds),
        )
    protein = _translate(reference_cds)
    if not protein.endswith("*"):
        raise PipelineError(
            f"The reference CDS translates to {len(protein)} residues and does "
            f"not end in a stop codon. The tblastn query must carry the "
            f"reference's own terminator, because the structural screen asks "
            f"whether the search reached it; a query without one can never "
            f"satisfy that question.",
            residues=len(protein),
        )
    return protein


def write_tblastn_query(path: Path, reference_cds: str) -> Path:
    """Write the tblastn query FASTA to ``path`` and return it.

    The production query for the structural screen, and the only place it is
    built: :func:`run_tblastn` searches the file it is handed and deliberately
    does not pad or trim it, so a query written anywhere else is a query whose
    length nothing holds against :func:`reference_codon_count`.

    Written to disk rather than assembled in memory because tblastn takes a path,
    and writing it there keeps the exact bytes searched auditable after the run -
    the same reason :func:`resolve_isolate` writes its blastp query.
    """
    query = Path(path)
    query.parent.mkdir(parents=True, exist_ok=True)
    query.write_text(
        f">{TBLASTN_QUERY_ID}\n{tblastn_query_protein(reference_cds)}\n",
        encoding="utf-8",
    )
    return query


_COMPLEMENT = str.maketrans("ACGTNacgtnRYKMBVDHrykmbvdh", "TGCANtgcanYRMKVBHDyrmkvbhd")


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate(cds: str) -> str:
    return "".join(_CODONS[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))


def _read_contig(fasta: Path, contig: str) -> str:
    """One record's sequence from a FASTA, read without a third-party parser."""
    wanted = contig.split()[0]
    chunks: List[str] = []
    seen = False
    with Path(fasta).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if seen:
                    break
                seen = line[1:].split()[0] == wanted
                continue
            if seen:
                chunks.append(line.strip())
    if not seen:
        raise PipelineError(
            f"Contig {wanted} not found in {fasta}. The GFF's contig and the "
            f"FASTA's must be the same assembly release, or every coordinate "
            f"shifts silently.",
            contig=contig,
        )
    return "".join(chunks)


def blastp_command(
    *,
    program: str,
    query: Path,
    subject: Path,
    threads: int,
    evalue: float = EVALUE,
) -> List[str]:
    """Tabular blastp of the reference protein against one isolate's proteome.

    Query start/end are requested because :func:`cluster_hits` merges fragments
    by their span on the reference; without them a split locus is two unmergeable
    hits. blastp exits 0 and prints nothing when it finds nothing, which reads as
    "no orthologue" - the same silent-empty-output shape the virulence adapter
    documents for blastn.
    """
    return [
        program,
        "-query", str(query),
        "-subject", str(subject),
        "-outfmt", OUTFMT,
        "-evalue", str(evalue),
        "-max_target_seqs", str(MAX_TARGET_SEQS),
        "-num_threads", str(threads),
    ]


def run_blastp(
    query: Path,
    subject: Path,
    *,
    program: str = "blastp",
    threads: int = 1,
    runner=None,
    query_length: Optional[int] = None,
) -> List[Hit]:
    """Run blastp and return its hits.

    ``runner`` is injectable so the command construction and the failure path are
    testable without blast installed, matching
    :func:`papipeline.adapters.virulence.screen_isolate`.
    """
    if runner is None:
        def runner(command: Sequence[str]) -> Tuple[int, str, str]:
            completed = subprocess.run(
                list(command), capture_output=True, text=True, check=False
            )
            return completed.returncode, completed.stdout, completed.stderr

    length = query_length if query_length is not None else len(
        "".join(
            line.strip()
            for line in Path(query).read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
            if line and not line.startswith(">")
        )
    )
    code, out, err = runner(
        blastp_command(
            program=program, query=Path(query), subject=Path(subject),
            threads=threads,
        )
    )
    if code != 0:
        raise ToolExecutionError(
            f"blastp failed for {Path(subject).name} (exit {code}): "
            f"{err.strip()[:400]}",
            command=" ".join(
                blastp_command(
                    program=program, query=Path(query), subject=Path(subject),
                    threads=threads,
                )
            ),
        )
    return parse_blast_table(out, query_length=length)


def resolve_isolate(
    sample_id: str,
    proteome: Path,
    gff3: Path,
    reference: str,
    *,
    workdir: Path,
    program: str = "blastp",
    threads: int = 1,
    runner=None,
    min_identity_pct: float = MIN_IDENTITY_PCT,
    min_coverage_pct: float = MIN_COVERAGE_PCT,
) -> LocusResolution:
    """End-to-end for one isolate: translate the query, blastp it, decide.

    The query is written to ``workdir`` rather than assembled in memory because
    blastp takes a path, and writing it there keeps the exact bytes searched
    auditable after the run.
    """
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    query = work / "reference_protein.faa"
    query.write_text(f">PA0958 oprD\n{reference}\n", encoding="utf-8")

    hits = run_blastp(
        query, Path(proteome), program=program, threads=threads, runner=runner,
        query_length=len(reference),
    )
    located = attach_locations(hits, load_cds_locations(Path(gff3)))
    return resolve(
        sample_id,
        located,
        min_identity_pct=min_identity_pct,
        min_coverage_pct=min_coverage_pct,
    )

# ---------------------------------------------------------------------------
# Structural classification
#
# Everything above answers "is this the orthologue locus?". Nothing above
# answers "is the gene WHOLE?", which is the question the study's central
# negative actually turns on: an oprD locus that is present but whose reading
# frame dies at codon 237 is not an intact oprD, and calling it `intact` because
# blast found it would be the same category error as calling a paralog `intact`.
#
# The classification is by LESION TYPE and it is structural: the verdict is
# driven by what broke the reading frame - a frameshift, a premature stop, or
# neither - and never by how long the resulting ORF is. Truncation length is
# recorded as evidence on the verdict (LocusEvidence.truncation_aa) so a reader
# can see how much of the protein is missing without that number being able to
# decide anything. `sensitivity_table` sweeps the two PRESENCE thresholds and
# nothing else; there is no structural threshold left to sweep.
#
# RULING 2 is what "disrupted" means here, stated once so it cannot be re-read:
# a frameshift OR a premature stop, and nothing else. An in-frame deletion that
# preserves the NATIVE stop is neither, so it is `intact` however many residues
# it removes. That case needed a measurement rather than a comparison - see
# LocusEvidence.stop_is_reference_terminator and is_terminal_stop, and the
# reasoning in report-step1.md in this round's artefact directory.
# ---------------------------------------------------------------------------

#: A cluster must reach this identity to be considered the orthologue locus at
#: all. Defaults to the module's existing identity floor; structural
#: classification inherits it rather than inventing a second one.
#:
#: THIS IS A PRESENCE THRESHOLD, not a structure threshold. It answers "is the
#: locus there at all", which is the only question the aligned span is for.
STRUCTURAL_MIN_IDENTITY_PCT = MIN_IDENTITY_PCT

#: Fraction of the reference protein the locus must reconstruct before the
#: question "is it whole?" is even meaningful. Below this the locus is not
#: located well enough to say anything structural about it.
#:
#: THIS IS A PRESENCE THRESHOLD. Measured 99.32-100% across the ten smoke
#: isolates, so it is not sensitive on real data; it exists so a partial
#: assembly cannot be mistaken for a whole locus.
STRUCTURAL_MIN_COVERAGE_PCT = MIN_COVERAGE_PCT

#: tblastn search parameters. `db_gencode` 11 is bacterial code; blast's default
#: is 1 (standard), which is wrong for P. aeruginosa. `max_target_seqs` must be
#: >= 1 in blast 2.17; `0` is rejected as legacy syntax.
TBLASTN_EVALUE = 1.0e-3
TBLASTN_DB_GENCODE = 11
TBLASTN_MAX_TARGET_SEQS = 5000
#: `qseq`/`sseq` are requested and they are LOAD-BEARING, not decoration. The
#: aligned columns are what place an indel at a reference position:
#: :func:`indels_from_alignment` reads the gap columns of each HSP to find the
#: deletions and insertions inside it, and the coordinate geometry to find the
#: one between two HSPs. Without them the lesion positions are not derivable,
#: and a repeat-context-aware aligner asked to rediscover them from the raw
#: sequences scatters inside exactly the homopolymers and tandem duplications
#: this screen is looking for. These fourteen fields are the outfmt the stored
#: round artifacts were produced with; see ``tblastn/*.cmd``.
TBLASTN_OUTFMT_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qseq", "sseq",
)
TBLASTN_OUTFMT = "6 " + " ".join(TBLASTN_OUTFMT_FIELDS)

#: The aligned-column field names, so a row's last two columns are named rather
#: than counted.
TBLASTN_ALIGNED_QUERY_FIELD = "qseq"
TBLASTN_ALIGNED_SUBJECT_FIELD = "sseq"

#: The FASTA id the tblastn query carries, and therefore the ``qseqid`` every row
#: of its table reports. It is the PA0958 locus tag joined by an underscore
#: because blast truncates a FASTA header at the first whitespace, so a header
#: like ``>PA0958 oprD`` would reach the table as ``qseqid=PA0958`` and silently
#: rename the query the stored artifacts were produced with.
TBLASTN_QUERY_ID = "PA0958_oprD"


class StructuralVerdict:
    """The four structural outcomes. Every one of them carries its evidence."""

    INTACT = "intact"
    DISRUPTED = "disrupted"
    ABSENT = "absent"
    NOT_ASSESSED = "not_assessed"

    ALL = frozenset({INTACT, DISRUPTED, ABSENT, NOT_ASSESSED})

    #: The only verdict that licenses an `absent` in the figures.
    ABSENCE = frozenset({ABSENT})
    LOSS_OF_FUNCTION = frozenset({DISRUPTED})


class LesionType:
    """What KIND of lesion the reading frame carries, if any.

    This vocabulary replaced the ORF-length bands, and it is the whole of the
    classification. The rule the module implements is RULING 2:

    * :data:`FRAMESHIFT`, :data:`PREMATURE_STOP`, :data:`BOTH` -> ``disrupted``
    * :data:`NONE` -> ``intact``

    **`disrupted` means a frameshift or a premature stop, and nothing else.**
    In particular an in-frame deletion that preserves the NATIVE stop is
    neither, so it is not ``disrupted`` - not however many residues it removes,
    and not at any ORF-length fraction. The distinction is carried by
    :attr:`LocusEvidence.stop_is_reference_terminator`: a stop that the
    measurement shows is the reference's own terminator is not a lesion, while
    a stop created at a codon that encodes a normal residue in the reference is.

    There is deliberately no fifth value for "short but unbroken". A locus that
    is two residues shorter than PAO1 with the reading frame intact is not in a
    band; it is `intact` with a truncation attribute, and truncating it is a
    statement about the measurement, not a threshold this module is entitled to
    set.

    A locus split across two abutting CDS by a draft-assembly gap is NOT a lesion
    type. The reading frame runs through the gap and the protein is one protein;
    the split is a property of the annotation, and it is carried as evidence so
    a reader can see it, not as a verdict input.
    """

    NONE = "no_lesion"
    FRAMESHIFT = "frameshift"
    PREMATURE_STOP = "premature_stop"
    BOTH = "frameshift_with_premature_stop"

    ALL = frozenset({NONE, FRAMESHIFT, PREMATURE_STOP, BOTH})

    #: The lesion types that make the reading frame unable to encode the
    #: reference protein. `NONE` is not in it.
    DISRUPTING = frozenset({FRAMESHIFT, PREMATURE_STOP, BOTH})


@dataclass(frozen=True)
class LocusSearch:
    """The search that was run, recorded so a negative result is auditable.

    ``absent`` is only reachable through this object. A negative claim about a
    locus is the one claim in this pipeline that cannot be re-derived from the
    outputs alone, so the command, the database and the parameters travel with
    the verdict rather than being left in a log.
    """

    program: str
    query: str
    database: str
    parameters: Tuple[Tuple[str, str], ...]
    command: str
    n_hits: int

    @property
    def is_recorded(self) -> bool:
        return bool(self.command and self.database and self.program)

    def as_row(self) -> Dict[str, Any]:
        return {
            "program": self.program,
            "query": self.query,
            "database": self.database,
            "parameters": " ".join(f"{k}={v}" for k, v in self.parameters),
            "command": self.command,
            "n_hits": self.n_hits,
        }


@dataclass(frozen=True)
class LocusEvidence:
    """Everything the verdict was decided on. No verdict cites anything else."""

    identity_pct: Optional[float] = None
    coverage_pct: float = 0.0
    orf_aa_length: Optional[int] = None
    reference_aa_length: Optional[int] = None
    internal_stop_codon: Optional[str] = None
    internal_stop_position_aa: Optional[int] = None
    stop_is_reference_terminator: Optional[bool] = None
    frameshift: bool = False
    split_abutting_fragments: bool = False
    at_contig_edge: bool = False
    contig: str = ""
    strand: str = ""
    footprint_1based: Tuple[int, int] = (0, 0)
    n_hsp: int = 0

    @property
    def orf_length_fraction(self) -> Optional[float]:
        """Reported, never decided on.

        Kept because it is a true and useful description of the measurement, and
        because removing it would hide that the number still exists. It is NOT
        read by :func:`classify_structure` and no threshold is derived from it:
        :func:`test_the_orf_length_cannot_move_a_verdict` proves the verdict is
        invariant across every value of ``orf_aa_length``.

        Note the denominator convention: callers pass the reference's codon
        count *including* its terminal stop (444 for PAO1 PA0958, which is 443
        residues), so this fraction is understated by 1/444 = 0.225 points
        against a residue count. That is a fixture convention, not a decision.
        """
        if not self.orf_aa_length or not self.reference_aa_length:
            return None
        return self.orf_aa_length / self.reference_aa_length

    @property
    def is_terminal_stop(self) -> bool:
        """True when the recorded stop IS the reference's own terminator.

        PAO1 PA0958's stop sits at codon 444 of a 444-codon reference. A stop
        there is how a gene ends, not a lesion: every ORF ends on a stop codon,
        so treating the terminator as damage calls every complete oprD gene
        disrupted.

        **Two tests, and the codon INDEX alone is not one of them.** Comparing
        ``internal_stop_position_aa`` against ``reference_aa_length`` only holds
        while the isolate's reading frame is exactly as long as the reference's,
        which is the one assumption an in-frame deletion breaks. A net -6 nt
        in-frame deletion moves the terminator from codon 444 to codon 442 while
        leaving it the reference's own TAA; on the index test alone that reads
        as a premature stop at 442, and the isolate is called `disrupted` for a
        lesion it does not have. That is RULING 2's case, and
        ``TestRuling2AnInFrameDeletionIsNotDisruption`` pins it.

        So the index test is kept **only** as the fallback for an unmeasured
        terminator, and the load-bearing test is
        :attr:`stop_is_reference_terminator`, which the measurer sets by
        aligning the isolate's CDS to the reference CDS and asking whether the
        ORF runs through the reference's own terminator in frame.

        ``stop_is_reference_terminator=None`` means the question was not asked.
        That case resolves to the index test, which is the pre-existing
        conservative behaviour: an unasked question must never manufacture an
        ``intact``. The measurements do not have to agree - a terminus recorded
        as native wins, because it is the direct observation - but they never
        have to be reconciled by the caller either.
        """
        if self.internal_stop_codon is None or self.internal_stop_position_aa is None:
            return False
        if self.stop_is_reference_terminator is not None:
            return bool(self.stop_is_reference_terminator)
        return (
            bool(self.reference_aa_length)
            and self.internal_stop_position_aa >= self.reference_aa_length
        )

    @property
    def has_premature_stop(self) -> bool:
        """A stop codon before the reference's own terminator."""
        return (
            self.internal_stop_codon is not None
            and self.internal_stop_position_aa is not None
            and bool(self.reference_aa_length)
            and not self.is_terminal_stop
        )

    @property
    def stop_position_unknown(self) -> bool:
        """A stop codon was seen but not placed.

        An unplaced stop cannot be told from the terminator, so it cannot be
        classified. This routes to ``not_assessed`` rather than being resolved
        in either direction by guesswork.
        """
        return (
            self.internal_stop_codon is not None
            and self.internal_stop_position_aa is None
        )

    @property
    def lesion_type(self) -> str:
        """The lesion vocabulary value for this evidence. Derived, not stored.

        Derived rather than stored so a caller cannot assert a lesion type that
        contradicts the booleans it also supplied.
        """
        premature = self.has_premature_stop
        if self.frameshift and premature:
            return LesionType.BOTH
        if self.frameshift:
            return LesionType.FRAMESHIFT
        if premature:
            return LesionType.PREMATURE_STOP
        return LesionType.NONE

    @property
    def truncation_aa(self) -> Optional[int]:
        """Reference residues this ORF does not encode. EVIDENCE, not a decision.

        Counted from the lesion rather than from a length fraction, because the
        lesion is what truncates: a premature stop at codon ``stop_pos`` removes
        the ``reference_aa_length - stop_pos`` residues that lie C-terminal to
        it. With the 444-codon PAO1 reference this gives 0 for a complete ORF
        (stop at 444) and 8 for PDT000292998.1 (stop at 436) - which is the
        reference's C-terminal octamer, independently counted.

        :func:`classify_structure` reads this only to print it. Changing it
        cannot change a verdict.
        """
        if not self.reference_aa_length:
            return None
        if self.has_premature_stop:
            return max(
                0, self.reference_aa_length - int(self.internal_stop_position_aa)
            )
        if self.orf_aa_length is None:
            return None
        # No stop to count from, so fall back to the encoded length. The `- 1`
        # is the reference's terminal stop codon, which encodes no residue.
        return max(0, self.reference_aa_length - 1 - self.orf_aa_length)

    def as_row(self) -> Dict[str, Any]:
        return {
            "identity_pct": self.identity_pct,
            "coverage_pct": self.coverage_pct,
            "orf_aa_length": self.orf_aa_length,
            "reference_aa_length": self.reference_aa_length,
            "orf_length_fraction": self.orf_length_fraction,
            "lesion_type": self.lesion_type,
            "truncation_aa": self.truncation_aa,
            "internal_stop_codon": self.internal_stop_codon,
            "internal_stop_position_aa": self.internal_stop_position_aa,
            "stop_is_reference_terminator": self.stop_is_reference_terminator,
            "frameshift": self.frameshift,
            "split_abutting_fragments": self.split_abutting_fragments,
            "at_contig_edge": self.at_contig_edge,
            "contig": self.contig,
            "strand": self.strand,
            "footprint_start": self.footprint_1based[0],
            "footprint_end": self.footprint_1based[1],
            "n_hsp": self.n_hsp,
        }


@dataclass(frozen=True)
class StructuralCall:
    """A structural verdict plus the evidence and the search behind it."""

    sample_id: str
    verdict: str
    reason: str
    evidence: LocusEvidence = field(default_factory=LocusEvidence)
    search: Optional[LocusSearch] = None

    @property
    def is_absent(self) -> bool:
        return self.verdict in StructuralVerdict.ABSENCE

    @property
    def is_loss_of_function(self) -> bool:
        return self.verdict in StructuralVerdict.LOSS_OF_FUNCTION

    @property
    def is_decided(self) -> bool:
        return self.verdict in (StructuralVerdict.INTACT,
                                StructuralVerdict.DISRUPTED,
                                StructuralVerdict.ABSENT)

    def as_row(self) -> Dict[str, Any]:
        row = {
            "sample_id": self.sample_id,
            "structural_verdict": self.verdict,
            "structural_reason": self.reason,
        }
        row.update(self.evidence.as_row())
        row["search_recorded"] = bool(self.search and self.search.is_recorded)
        if self.search:
            row["search_command"] = self.search.command
            row["search_database"] = self.search.database
            row["search_parameters"] = " ".join(
                f"{k}={v}" for k, v in self.search.parameters)
            row["search_n_hits"] = self.search.n_hits
        return row


def classify_structure(
    sample_id: str,
    evidence: LocusEvidence,
    *,
    search: Optional[LocusSearch] = None,
    min_identity_pct: float = STRUCTURAL_MIN_IDENTITY_PCT,
    min_coverage_pct: float = STRUCTURAL_MIN_COVERAGE_PCT,
) -> StructuralCall:
    """Classify one locus structurally, by LESION TYPE.

    **The verdict is driven by what broke the reading frame, not by how long the
    ORF is.** A frameshift or a premature stop is ``disrupted``. No lesion is
    ``intact``. Truncation length is carried as
    :attr:`LocusEvidence.truncation_aa` and read only to print it.

    There is no ORF-length parameter on this signature and there is no band to
    be in. An earlier revision banded ``intact`` against ``disrupted`` at
    ``INTACT_MIN_ORF_LENGTH_FRACTION`` / ``DISRUPTED_MAX_ORF_LENGTH_FRACTION``
    (both 0.90) and carried a third ``orf_length_in_undecided_band`` outcome for
    lengths between them. That put a chosen number between an observation and a
    verdict, and it swept: three of the ten smoke isolates changed verdict as
    the band moved. A length fraction also cannot see a frameshift, which is the
    one lesion that matters most here - PDT000292998.1 reconstructs 0.9797 of the
    reference and is frameshifted at codon 402.

    **What the band agreed on, and why that was luck.** It agreed with RULING 2
    on PDT000034122.1 (0.9932, `intact`), because an in-frame deletion keeps the
    fraction high. It agreed by comparing a length to a chosen number rather than
    by observing that the reading frame never left the native frame, which is why
    it was 0.90 and not 0.98 and why moving it moved verdicts.

    The two surviving thresholds are **presence** thresholds. The aligned span
    answers "is the locus there at all", and only that:

    The order of the checks is the order of what may be asserted.

    1. ``absent`` is reachable ONLY from a recorded search that returned no hit.
       Without :class:`LocusSearch` there is no search to audit, so the answer is
       ``not_assessed`` with the reason ``no_recorded_search`` - never
       ``absent``. This is the whole point of the function.
    2. A hit below the identity or coverage floor has not located the orthologue
       locus, so nothing structural can be said about it.
    3. A stop codon that was seen but not placed cannot be classified at all.
    4. A frameshift, a premature stop, or both is ``disrupted``.
    5. A locus that runs off the end of its contig has an unknown extent, so its
       completeness cannot be asserted: ``not_assessed``.
    6. A locus with no lesion but no measured ORF has not been translated, so
       "no lesion observed" is not "no lesion exists": ``not_assessed``.
    7. Anything else is ``intact``.

    Every ``not_assessed`` carries a NAMED reason and is never ``absent``.
    """
    def call(verdict: str, reason: str) -> StructuralCall:
        return StructuralCall(sample_id=sample_id, verdict=verdict,
                              reason=reason, evidence=evidence, search=search)

    located = evidence.coverage_pct >= min_coverage_pct and (
        evidence.identity_pct is None
        or evidence.identity_pct >= min_identity_pct)

    # 1. absence, and only ever from a recorded negative search
    if search is None:
        return call(
            StructuralVerdict.NOT_ASSESSED,
            "no_recorded_search: no tblastn search was recorded for "
            f"{sample_id}, so a negative result cannot be audited and the "
            "locus is not called absent.",
        )
    if not located and evidence.coverage_pct == 0.0 and search.n_hits == 0:
        return call(
            StructuralVerdict.ABSENT,
            f"absent: {search.program} of the PAO1 oprD protein against "
            f"{search.database} returned {search.n_hits} hits "
            f"({search.parameters}). No orthologous locus is present in this "
            "assembly.",
        )

    # 2. the locus was not located well enough to judge structurally. This is
    #    the aligned span doing its one job.
    if not located:
        id_text = ("identity unmeasured"
                   if evidence.identity_pct is None
                   else f"{evidence.identity_pct:.2f}% identity")
        return call(
            StructuralVerdict.NOT_ASSESSED,
            f"locus_not_located: best cluster at "
            f"{evidence.coverage_pct:.1f}% coverage and {id_text}, below the "
            f"{min_coverage_pct}%/{min_identity_pct}% presence floors. Structure "
            "cannot be assessed on a locus that was not located.",
        )

    # 3. a stop with no position cannot be told from the terminator
    if evidence.stop_position_unknown:
        return call(
            StructuralVerdict.NOT_ASSESSED,
            f"stop_position_unknown: a {evidence.internal_stop_codon} stop was "
            f"recorded for {sample_id} but not its codon, and a stop at codon "
            f"{evidence.reference_aa_length} is the reference's own terminator "
            "rather than a lesion. Placing it is required before this can be "
            "called intact or disrupted.",
        )

    # 4. the lesion list, for the reason string. EVERY lesion observed is
    #    carried, not just the first: an isolate can be both frameshifted and
    #    premature-stopped, and the verdict has to be auditable against all of
    #    the evidence rather than the top of it.
    ref_len = evidence.reference_aa_length or 0
    stop_pos = evidence.internal_stop_position_aa
    truncation = evidence.truncation_aa
    lesions: List[str] = []
    if evidence.has_premature_stop:
        lesions.append(
            f"premature stop {evidence.internal_stop_codon} at codon {stop_pos} "
            f"of {ref_len}, truncating {truncation} reference residues")
    if evidence.frameshift:
        lesions.append(
            f"frameshift across "
            f"{evidence.footprint_1based[1] - evidence.footprint_1based[0] + 1} bp "
            f"on {evidence.contig}({evidence.strand})")
    if evidence.split_abutting_fragments:
        # NOT a lesion: the reading frame runs through the assembly gap and the
        # protein is one protein. Recorded because it changes how the locus must
        # be read, not because it changes the verdict.
        lesions.append(
            "split into abutting same-strand same-contig fragments on "
            f"{evidence.contig}({evidence.strand}) spanning "
            f"{evidence.footprint_1based[0]}-{evidence.footprint_1based[1]} "
            "(annotation split, not a reading-frame lesion)")
    if evidence.at_contig_edge:
        lesions.append(
            f"locus reaches a contig edge on {evidence.contig}"
            f"({evidence.strand})")
    lesion = "; ".join(lesions)

    # 5. LESION TYPE decides. Nothing below this line reads a length.
    lesion_type = evidence.lesion_type
    if lesion_type in LesionType.DISRUPTING:
        return call(
            StructuralVerdict.DISRUPTED,
            f"lesion_type={lesion_type}: the reading frame cannot encode the "
            f"reference protein, so the locus is disrupted"
            + (f". Observed: {lesion}." if lesion else ".")
            + f" Truncation is reported as evidence, not as the basis of this "
            f"call"
            + (f": {truncation} of {ref_len} reference residues."
               if truncation is not None else "."),
        )

    # 6. an unknown extent cannot be called intact
    if evidence.at_contig_edge:
        return call(
            StructuralVerdict.NOT_ASSESSED,
            f"locus_truncated_by_contig_edge: no lesion was found, but the "
            f"locus on {evidence.contig}({evidence.strand}) reaches a contig "
            "end, so how much of it exists is unknown. Calling it intact would "
            "assert a completeness that was not measured."
            + (f" Observed: {lesion}." if lesion else ""),
        )

    # 7. no ORF measured means no lesion was looked for
    if evidence.orf_aa_length is None:
        return call(
            StructuralVerdict.NOT_ASSESSED,
            "orf_not_measured: the locus was located and no lesion was "
            "recorded, but no ORF length was measured, so the reading frame was "
            "never translated and 'no lesion observed' is not 'no lesion "
            "exists'."
            + (f" Observed: {lesion}." if lesion else ""),
        )

    # 8. no lesion, extent known, frame translated: intact.
    return call(
        StructuralVerdict.INTACT,
        f"intact: lesion_type={lesion_type}, so the reading frame encodes the "
        f"reference protein. ORF of {evidence.orf_aa_length} aa against a "
        f"{ref_len} codon reference"
        + (f", truncating {truncation} reference residue(s)"
           if truncation else ", with no truncation")
        + f", at {evidence.coverage_pct:.1f}% aligned span and "
        + (f"{evidence.identity_pct:.2f}% identity"
           if evidence.identity_pct is not None else "identity unmeasured")
        + f" on {evidence.contig}({evidence.strand})"
        + (f". Note: {lesion}." if lesion else "."),
    )


def sensitivity_table(
    calls: Sequence[StructuralCall],
    *,
    identity_floors: Sequence[float] = (70.0, 80.0, 90.0, 95.0),
    coverage_floors: Sequence[float] = (80.0, 90.0, 95.0, 99.0),
) -> List[Dict[str, Any]]:
    """How each isolate's verdict moves across the two SURVIVING thresholds.

    Only the **presence** thresholds are swept, because only they are
    thresholds: ``min_identity_pct`` and ``min_coverage_pct`` decide whether the
    locus was located at all, which is the one question the aligned span answers.

    **The ORF-length columns are gone, and that is the result, not an omission.**
    The previous revision swept ``INTACT_MIN_ORF_LENGTH_FRACTION`` across
    0.80-0.99 and three of the ten isolates changed verdict as it moved, which
    is what a threshold in the classification path looks like from the outside.
    There is no ORF-length threshold left to sweep. The length is still reported
    per row - ``orf_aa_length``, ``orf_length_fraction``, ``truncation_aa`` and
    ``lesion_type`` - as evidence, and
    ``test_the_orf_length_cannot_move_a_verdict`` shows the verdict is invariant
    when ``orf_aa_length`` is set to every value in a wide range.

    A verdict may still move along these rows, and when it does it is because
    the locus stopped being *located*, not because its reading frame was
    re-graded.
    """
    rows: List[Dict[str, Any]] = []
    for call in calls:
        row: Dict[str, Any] = {
            "sample_id": call.sample_id,
            "lesion_type": call.evidence.lesion_type,
            "orf_aa_length": call.evidence.orf_aa_length,
            "orf_length_fraction": call.evidence.orf_length_fraction,
            "truncation_aa": call.evidence.truncation_aa,
            "verdict": call.verdict,
        }
        for floor in identity_floors:
            redone = classify_structure(
                call.sample_id, call.evidence, search=call.search,
                min_identity_pct=floor,
            )
            row[f"id_floor_{floor:.0f}"] = redone.verdict
        for floor in coverage_floors:
            redone = classify_structure(
                call.sample_id, call.evidence, search=call.search,
                min_coverage_pct=floor,
            )
            row[f"cov_floor_{floor:.0f}"] = redone.verdict
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# THE PRODUCTION STRUCTURAL PATH
#
# Everything above `sensitivity_table` decides; nothing above it MEASURES. This
# section is the production populator: it turns a tblastn table plus an assembly
# into the :class:`LocusEvidence` that :func:`classify_structure` reads, so
# `stop_is_reference_terminator` and `frameshift` have a producer in a real run
# and not only in a test table.
#
# **The terminator question is asked of the ALIGNMENT.** It used to be asked of
# the codon index - `internal_stop_position_aa >= reference_aa_length` - which is
# only true while the isolate's reading frame is exactly as long as the
# reference's, i.e. precisely the assumption an in-frame deletion breaks. Three
# statements, all read off the alignment, replace it:
#
#   1. the reading frame is NATIVE across the locus, i.e. the aligned span
#      consumed `3 * reference_codons` nucleotides, so `net_frame_change_nt` is a
#      multiple of three;
#   2. the alignment REACHES the reference's terminal codon, i.e. the highest
#      aligned reference codon is the reference's last;
#   3. the isolate's own reading frame runs to the end of the aligned span
#      before it stops, i.e. its first in-frame stop is the final codon.
#
# Together those say the isolate's stop IS the reference's own terminator,
# reached in frame, without ever comparing a codon number to a codon count. See
# `TestTheProductionResolverDecidesTheTerminatorFromTheAlignment`.
#
# **A production path must never leave the question unasked.**
# `stop_is_reference_terminator=None` still falls back to the codon index inside
# :attr:`LocusEvidence.is_terminal_stop`, which is the right conservative
# behaviour for a caller that did not measure. :func:`structural_call` therefore
# REFUSES with `frame_not_anchorable` when no HSP reaches the reference's first
# codon, instead of handing `classify_structure` an unmeasured terminator. There
# is no route from this section into the fallback.
# ---------------------------------------------------------------------------

#: Reference codons either side of an indel that are inspected for local frame
#: compensation, and for the repeat context the indel sits in.
#:
#: 30 nt is the neighbourhood the measured lesion table used, and it is kept so
#: the two are comparable: compensating pairs on these ten isolates are 17-26 nt
#: apart, and nothing between 27 and 30 nt is claimed either way.
INDEL_NEIGHBOURHOOD_NT = 30

#: A run of one base this long, overlapping an indel, is a homopolymer context.
#:
#: Measured, not chosen to fit. **The motif is searched on the REFERENCE CDS**,
#: which is what the locus's coordinates are stated in, and the shortest run any
#: of the ten isolates presents at an indel is a 5-mer ``G`` run at reference nt
#: 631 - PDT000167135.1's single inserted ``G`` lands inside it and makes it a
#: 6-mer in the isolate. So the isolate's own tract is at least as long as the
#: one reported here, never shorter.
HOMOPOLYMER_MIN_RUN_NT = 5

#: Tandem duplication unit lengths searched around an indel.
#:
#: 5 nt is the smallest unit measured at an indel on these isolates
#: (`GTCCG` x3 in PDT000292998.1); 12 nt is the largest (`CTACGG` x2 in
#: PDT000034122.1). Both ends of the range are the measured ones.
TANDEM_MIN_UNIT_NT = 5
TANDEM_MAX_UNIT_NT = 12


@dataclass(frozen=True)
class AlignedNucleotideHit:
    """One tblastn HSP: the reference protein against one assembly region.

    ``subject_start``/``subject_end`` are blast's own 1-based subject
    coordinates, so ``subject_start > subject_end`` on the minus strand and
    :attr:`direction` is what orients the locus. ``query_aligned`` and
    ``subject_aligned`` are the HSP's aligned columns, which is where an indel
    inside one HSP is read off; they are ``None`` only if the search was run
    without them, and :func:`indels_from_alignment` says so rather than guessing.
    """

    query_id: str
    subject_id: str
    identity_pct: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bitscore: float
    query_aligned: Optional[str] = None
    subject_aligned: Optional[str] = None

    @property
    def query_low(self) -> int:
        return min(self.query_start, self.query_end)

    @property
    def query_high(self) -> int:
        return max(self.query_start, self.query_end)

    @property
    def subject_low(self) -> int:
        return min(self.subject_start, self.subject_end)

    @property
    def subject_high(self) -> int:
        return max(self.subject_start, self.subject_end)

    @property
    def is_reverse(self) -> bool:
        return self.subject_start > self.subject_end

    @property
    def direction(self) -> int:
        """+1 if the gene reads along increasing subject coordinate, else -1."""
        return -1 if self.is_reverse else 1

    @property
    def query_codons(self) -> range:
        return range(self.query_low, self.query_high + 1)


def parse_tblastn_table(text: str, query_codon_length: int) -> List[
    AlignedNucleotideHit
]:
    """Parse ``TBLASTN_OUTFMT`` tabular output into HSPs.

    Two widths are accepted and only two: the fourteen declared columns, and the
    twelve columns without the aligned sequences. Anything else is refused rather
    than parsed, for the reason :func:`parse_blast_table` gives - a different
    ``-outfmt`` would make the reading frame silently wrong.

    The twelve-column form parses, and its hits carry no aligned columns, so
    :func:`indels_from_alignment` cannot place an indel within them. That is a
    lost measurement rather than a wrong one, and the production caller is
    expected to run the declared outfmt.
    """
    n_declared = len(TBLASTN_OUTFMT_FIELDS)
    allowed = (n_declared, n_declared - 2)
    hits: List[AlignedNucleotideHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) not in allowed:
            raise PipelineError(
                f"tblastn output row has {len(fields)} fields, expected "
                f"{n_declared} ({' '.join(TBLASTN_OUTFMT_FIELDS)}) or "
                f"{n_declared - 2} without the aligned sequences. A different "
                f"-outfmt would place the indels and therefore the reading "
                f"frame silently wrong, so this is refused rather than parsed. "
                f"Row: {line[:200]}",
                fields=len(fields),
            )
        has_columns = len(fields) == n_declared
        try:
            hits.append(
                AlignedNucleotideHit(
                    query_id=fields[0],
                    subject_id=fields[1],
                    identity_pct=float(fields[2]),
                    alignment_length=int(fields[3]),
                    query_length=int(fields[4]),
                    subject_length=int(fields[5]),
                    query_start=int(fields[6]),
                    query_end=int(fields[7]),
                    subject_start=int(fields[8]),
                    subject_end=int(fields[9]),
                    evalue=float(fields[10]),
                    bitscore=float(fields[11]),
                    query_aligned=fields[12] if has_columns else None,
                    subject_aligned=fields[13] if has_columns else None,
                )
            )
        except ValueError as exc:
            raise PipelineError(
                f"tblastn output row is not numeric where it must be: {exc}. "
                f"Row: {line[:200]}",
                row=line[:200],
            ) from exc
    return hits


def tblastn_command(
    *,
    program: str,
    query: Path,
    database: Path,
    threads: int,
    evalue: float = TBLASTN_EVALUE,
    max_target_seqs: int = TBLASTN_MAX_TARGET_SEQS,
) -> List[str]:
    """The tblastn command the structural path runs.

    Exactly the command that produced the stored round artifacts, including
    ``-db_gencode 11``: blast's default translation table is 1, the standard
    code, which is wrong for *P. aeruginosa*, and a table-1 translation of a
    bacterial locus invents stops. ``-num_threads`` is passed in and never
    hard-coded.
    """
    return [
        program,
        "-query", str(query),
        "-db", str(database),
        "-outfmt", TBLASTN_OUTFMT,
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-db_gencode", str(TBLASTN_DB_GENCODE),
        "-num_threads", str(threads),
    ]


@dataclass(frozen=True)
class Indel:
    """One insertion or deletion inside the locus, in reference coordinates.

    ``length_nt`` is positive for both kinds; :attr:`kind` says which, and
    :attr:`net_nt` is the signed effect on the reading frame (negative removes
    nucleotides, positive adds them).

    ``reference_nt_1based`` is the reference CDS nucleotides the indel covers.
    For a deletion that is the deleted span. For an insertion it is the single
    reference nucleotide the inserted bases sit in front of, because an
    insertion covers no reference nucleotide and saying otherwise would put a
    fabricated coordinate on it.
    """

    kind: str
    length_nt: int
    reference_nt_1based: Tuple[int, int]
    source: str

    @property
    def is_deletion(self) -> bool:
        return self.kind == "deletion"

    @property
    def net_nt(self) -> int:
        return -self.length_nt if self.is_deletion else self.length_nt

    @property
    def position_nt(self) -> int:
        return self.reference_nt_1based[0]


def indels_from_alignment(
    hits: Sequence[AlignedNucleotideHit],
    *,
    reference_codon_low: int = 1,
) -> Tuple[Indel, ...]:
    """Every indel the alignment itself shows, in reference coordinates.

    Two sources, and together they are exhaustive for a locus whose HSPs the
    search reported:

    * **inside one HSP**, the aligned columns. A column with a query residue and
      a subject gap is a deletion of that residue's nucleotides at that
      reference codon; a column with a query gap and a subject residue is an
      insertion in front of that reference codon. Consecutive columns of the same
      kind are merged, so a 12 nt deletion reads as one indel and not as four.
    * **between two HSPs**, the coordinates. HSP ``B``'s first aligned
      nucleotide sits where the reference says it should if the frame is native;
      the difference between where it does sit and where it should is an indel at
      the reference codon where the two HSPs meet. This is how a 5 nt insertion
      is found at reference codon 402 and a 1361 nt one at codon 172.

    **Why the alignment and not a fresh pairwise alignment.** The indels this
    screen cares about sit in homopolymers and tandem duplications, which are
    precisely where a re-aligner scatters one indel into a dozen. Anchoring on
    blast's own columns keeps a single inserted ``G`` in a ``G`` run a single
    inserted ``G``, and :func:`indels_from_alignment` is exact wherever blast
    aligned, which on this locus is everywhere - the search reports 99.32-100%
    query coverage.

    Raises if a hit carries no aligned columns, because the alternative is
    reporting "no indels" for a hit whose indels were never looked at.
    """
    ordered = sorted(hits, key=lambda h: (h.query_low, h.subject_low))
    missing = [h for h in ordered if h.query_aligned is None
               or h.subject_aligned is None]
    if missing:
        raise PipelineError(
            f"{len(missing)} of {len(ordered)} tblastn HSPs carry no aligned "
            f"columns, so the indels inside them cannot be read. Re-run the "
            f"search with the declared outfmt "
            f"({' '.join(TBLASTN_OUTFMT_FIELDS)}). Reporting 'no indels' here "
            f"would be a fabricated negative.",
            n_hits=len(missing),
        )

    indels: List[Indel] = []
    previous: Optional[AlignedNucleotideHit] = None
    for hit in ordered:
        direction = hit.direction
        if previous is not None:
            # Where this HSP's first aligned nucleotide SITS, in the coding-frame
            # coordinates the previous HSP established, minus where it is EXPECTED
            # from the reference alone. The difference is an indel in between.
            observed = _coding_position(previous, hit)
            expected = 3 * (hit.query_low - reference_codon_low) + 1
            delta = observed - expected
            if delta:
                # Located to the reference codon at which the two HSPs meet.
                # The exact nucleotides are not knowable from the alignment -
                # nothing was aligned there - so the codon is the finest honest
                # position, and it is the codon the compensation test uses.
                junction_nt = 3 * (hit.query_low - reference_codon_low) + 1
                indels.append(Indel(
                    kind="deletion" if delta < 0 else "insertion",
                    length_nt=abs(delta),
                    reference_nt_1based=(junction_nt, junction_nt),
                    source="hsp_junction",
                ))
        indels.extend(_indels_inside(hit))
        previous = hit
    indels.sort(key=lambda i: (i.reference_nt_1based[0], i.kind))
    return tuple(indels)


def _reference_nt_high(hit: AlignedNucleotideHit) -> int:
    """The highest reference CDS nucleotide this hit's alignment reaches."""
    return 3 * hit.query_high


def _coding_position(previous: AlignedNucleotideHit,
                     hit: AlignedNucleotideHit) -> int:
    """``hit``'s first aligned nucleotide, in ``previous``'s coding coordinates.

    Coding coordinates start at 1 on the nucleotide aligned to the reference's
    first aligned codon and increase along the gene, whichever subject strand
    that is.
    """
    anchor = previous.subject_start
    # `previous.subject_start` is the coordinate of its FIRST aligned nucleotide
    # on the gene, which blast reports as subject_start on the plus strand and as
    # subject_end on the minus strand. Normalise both to a signed offset.
    first = anchor
    offset = hit.subject_start - first
    if previous.direction != hit.direction:
        raise PipelineError(
            f"HSPs on {hit.subject_id} are reported on both strands "
            f"({previous.subject_start}-{previous.subject_end} and "
            f"{hit.subject_start}-{hit.subject_end}); the locus orientation is "
            f"ambiguous and the reading frame cannot be anchored.",
            subject_id=hit.subject_id,
        )
    return offset * previous.direction + 1


def _indels_inside(hit: AlignedNucleotideHit) -> List[Indel]:
    """Indels read off one HSP's aligned columns."""
    query, subject = hit.query_aligned or "", hit.subject_aligned or ""
    if len(query) != len(subject):
        raise PipelineError(
            f"tblastn HSP on {hit.subject_id} reports {len(query)} aligned query "
            f"residues against {len(subject)} aligned subject residues. The gap "
            f"columns cannot be read and are not guessed at.",
            subject_id=hit.subject_id,
        )
    out: List[Indel] = []
    # First reference CDS nucleotide of each aligned query codon.
    codon_first_nt = 3 * (hit.query_low - 1) + 1
    for column, (qc, sc) in enumerate(zip(query, subject)):
        if qc != "-" and sc != "-":
            continue
        kind = "deletion" if sc == "-" else "insertion"
        low = codon_first_nt + 3 * column
        span = (min(low, _reference_nt_high(hit)),
                min(low + 2, _reference_nt_high(hit)))
        if out and out[-1].kind == kind and out[-1].source == "hsp_columns" \
                and out[-1].reference_nt_1based[1] >= span[0] - 3:
            prior = out[-1]
            out[-1] = Indel(
                kind=kind,
                length_nt=prior.length_nt + 3,
                reference_nt_1based=(prior.reference_nt_1based[0], span[1]),
                source="hsp_columns",
            )
            continue
        out.append(Indel(kind=kind, length_nt=3,
                         reference_nt_1based=span, source="hsp_columns"))
    return out


def compensating_indel_group(
    indels: Sequence[Indel], *, window_nt: int = INDEL_NEIGHBOURHOOD_NT,
) -> Optional[Tuple[Indel, ...]]:
    """A pair of nearby indels whose combined length restores the frame.

    **The definition is local frame compensation, and it is evidence, never a
    decision.** Two or more indels within ``window_nt`` of each other whose net
    length is a multiple of three is a mechanism that can hide a frameshift -
    which is exactly why it is worth naming per isolate - but it does not make
    a locus ``intact``. PDT000034122.1 has such a pair and is ``intact`` for a
    different reason entirely: its reading frame never leaves the native one and
    it ends on PAO1's own ``TAA``. Nothing in
    :func:`classify_structure` reads this.

    Returns the first such group in reference order, or ``None``.
    """
    ordered = sorted(indels, key=lambda i: i.position_nt)
    for first in range(len(ordered)):
        group = [ordered[first]]
        for candidate in ordered[first + 1:]:
            if candidate.position_nt - ordered[first].position_nt > window_nt:
                break
            group.append(candidate)
        if len(group) < 2:
            continue
        if sum(i.net_nt for i in group) % 3 == 0:
            return tuple(group)
    return None


def repeat_context(
    reference_cds: str,
    indels: Sequence[Indel],
    *,
    window_nt: int = INDEL_NEIGHBOURHOOD_NT,
    homopolymer_min_run_nt: int = HOMOPOLYMER_MIN_RUN_NT,
    tandem_min_unit_nt: int = TANDEM_MIN_UNIT_NT,
    tandem_max_unit_nt: int = TANDEM_MAX_UNIT_NT,
) -> Optional[str]:
    """A homopolymer or tandem-repeat motif overlapping an indel, or ``None``.

    Reported per isolate because indels are not evenly distributed: three of the
    ten isolates present an indel inside a repeat, and a reviewer is entitled to
    know which. ``None`` means "no such motif was found near any indel", which is
    a statement about the reference sequence and not about the isolate's
    mechanism. Nothing reads this to decide anything.
    """
    best: Optional[Tuple[int, str]] = None
    for indel in indels:
        low = max(1, indel.position_nt - window_nt)
        high = min(len(reference_cds), indel.position_nt + window_nt)
        motif = _longest_homopolymer(
            reference_cds, low, high, homopolymer_min_run_nt)
        if motif:
            best = _keep_longest(best, motif)
        motif = _longest_tandem(
            reference_cds, low, high, tandem_min_unit_nt, tandem_max_unit_nt)
        if motif:
            best = _keep_longest(best, motif)
    return best[1] if best else None


def _keep_longest(current: Optional[Tuple[int, str]],
                  candidate: Tuple[int, str]) -> Tuple[int, str]:
    return candidate if current is None or candidate[0] > current[0] else current


def _longest_homopolymer(sequence: str, low: int, high: int,
                         min_run_nt: int) -> Optional[Tuple[int, str]]:
    """The longest run of one base in ``sequence[low:high]``, if long enough."""
    best: Optional[Tuple[int, str]] = None
    run = 0
    previous = ""
    for index in range(low, high + 1):
        base = sequence[index - 1]
        run = run + 1 if base == previous else 1
        previous = base
        if run >= min_run_nt and (best is None or run > best[0]):
            best = (run, base * run)
    return best


def _longest_tandem(sequence: str, low: int, high: int,
                    min_unit_nt: int, max_unit_nt: int
                    ) -> Optional[Tuple[int, str]]:
    """The longest tandem duplication in ``sequence[low:high]``, if any.

    Reported as ``unit x copies`` and scored on ``unit * copies``, the width of
    the repeated tract, so a 5-mer x3 outranks a 6-mer x2. A homopolymer is also
    a tandem duplication of a 1-mer, so this cannot return nothing where the
    homopolymer search did; both are tried and the wider tract wins.
    """
    best: Optional[Tuple[int, str]] = None
    for start in range(low, high + 1):
        for unit in range(min_unit_nt, max_unit_nt + 1):
            if start + 2 * unit - 1 > high:
                break
            motif = sequence[start - 1:start - 1 + unit]
            copies = 1
            while (start + (copies + 1) * unit - 1 <= high
                   and sequence[start - 1 + copies * unit:
                                start - 1 + (copies + 1) * unit] == motif):
                copies += 1
            if copies >= 2 and (best is None or copies * unit > best[0]):
                best = (copies * unit, f"{motif}x{copies}")
    return best


def read_assembly(path: Path) -> Dict[str, str]:
    """Every record of an isolate assembly, keyed by the first token of its id.

    ``_read_contig`` reads one record and is what the reference pin uses; the
    structural path needs the locus contig chosen from the alignment first and
    the assembly read afterwards, and reading the file twice to find out which
    record matters is a way to read a different release than the coordinates
    were computed against.
    """
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0] if line[1:].split() else ""
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


@dataclass(frozen=True)
class LocusMeasurement:
    """What production measured about one isolate's oprD reading frame.

    :attr:`evidence` is what :func:`classify_structure` reads and nothing else
    reads it. :attr:`indels`, :attr:`net_frame_change_nt`,
    :attr:`compensating_group` and :attr:`repeat_motif` are context for a reader
    and are never consulted by the classifier - asserted by
    `test_the_lesion_context_cannot_move_a_verdict`.
    """

    sample_id: str
    evidence: LocusEvidence
    indels: Tuple[Indel, ...] = ()
    net_frame_change_nt: int = 0
    compensating_group: Optional[Tuple[Indel, ...]] = None
    repeat_motif: Optional[str] = None
    frame_anchored: bool = False
    refusal_reason: Optional[str] = None

    @property
    def has_compensating_indel(self) -> bool:
        return self.compensating_group is not None


def measure_locus_structure(
    sample_id: str,
    hits: Sequence[AlignedNucleotideHit],
    *,
    reference_cds: str,
    assembly: Mapping[str, str],
    min_identity_pct: float = STRUCTURAL_MIN_IDENTITY_PCT,
) -> LocusMeasurement:
    """Measure one isolate's oprD reading frame from its tblastn alignment.

    The steps, and why each is here:

    1. keep the HSPs at or above the **presence** identity floor, which is the
       same floor :func:`resolve` uses for the same reason - below it the hit is
       a paralog;
    2. group them by contig and take the contig whose HSPs cover the most
       reference codons, which is what separates one locus from a paralog on a
       neighbouring scaffold. A tie is a refusal, never a coin flip;
    3. anchor the reading frame on the HSP that reaches the reference's **first**
       codon. Without it the frame cannot be established and the measurement
       refuses rather than guessing, because an unanchored frame is exactly how
       a stop gets called the terminator when it is not;
    4. read the isolate's own nucleotides along the aligned span, in the gene's
       frame, and translate;
    5. the net frame change is ``aligned nucleotides - 3 * aligned reference
       codons``. That is the frameshift measurement, and it comes from the
       alignment's arithmetic rather than from a free re-alignment;
    6. the terminator test is the three alignment statements above, and the
       codon index is never consulted.
    """
    reference_codons = reference_codon_count(reference_cds)
    candidates = [h for h in hits if h.identity_pct >= min_identity_pct]

    if not candidates:
        return LocusMeasurement(
            sample_id=sample_id,
            evidence=LocusEvidence(
                identity_pct=None, coverage_pct=0.0,
                reference_aa_length=reference_codons,
            ),
            refusal_reason=(
                f"no_candidate_hit: no HSP reaches the {min_identity_pct}% "
                f"identity presence floor in {sample_id}, so the orthologue "
                f"locus was not located."
            ),
        )

    by_contig: Dict[str, List[AlignedNucleotideHit]] = {}
    for hit in candidates:
        by_contig.setdefault(hit.subject_id, []).append(hit)

    def codon_union(hs: Sequence[AlignedNucleotideHit]) -> set:
        covered: set = set()
        for hit in hs:
            covered.update(hit.query_codons)
        return covered

    scored = sorted(
        ((len(codon_union(hs)), contig, hs) for contig, hs in by_contig.items()),
        key=lambda item: (-item[0], item[1]),
    )
    if len(scored) > 1 and scored[0][0] == scored[1][0]:
        return LocusMeasurement(
            sample_id=sample_id,
            evidence=LocusEvidence(
                identity_pct=None, coverage_pct=0.0,
                reference_aa_length=reference_codons,
            ),
            refusal_reason=(
                f"ambiguous_locus: {scored[0][1]} and {scored[1][1]} each cover "
                f"{scored[0][0]} of the {reference_codons} reference codons in "
                f"{sample_id}. A tie between two loci is not resolvable from "
                f"this evidence."
            ),
        )

    _n_covered, contig, locus = scored[0]
    covered = codon_union(locus)
    coverage_pct = 100.0 * len(covered) / reference_codons
    identity_pct = min(h.identity_pct for h in locus)
    ref_low = min(h.query_low for h in locus)
    ref_high = max(h.query_high for h in locus)
    foot_low = min(h.subject_low for h in locus)
    foot_high = max(h.subject_high for h in locus)
    strands = {h.direction for h in locus}
    if len(strands) != 1:
        return LocusMeasurement(
            sample_id=sample_id,
            evidence=LocusEvidence(
                identity_pct=identity_pct, coverage_pct=coverage_pct,
                reference_aa_length=reference_codons, contig=contig,
                footprint_1based=(foot_low, foot_high), n_hsp=len(locus),
            ),
            refusal_reason=(
                f"mixed_strands: the {len(locus)} HSPs on {contig} are reported "
                f"on both strands, so the locus has no single reading frame."
            ),
        )

    anchors = [h for h in locus if h.query_low == ref_low]
    if ref_low != 1 or not anchors:
        return LocusMeasurement(
            sample_id=sample_id,
            evidence=LocusEvidence(
                identity_pct=identity_pct, coverage_pct=coverage_pct,
                orf_aa_length=None,
                reference_aa_length=reference_codons,
                # False, not None: `None` means "the question was not asked" and
                # this measurer always asks it. Asserting False here keeps the
                # unasked case unreachable from production evidence, so the
                # codon-index fallback cannot be reached even by a caller that
                # skips the refusal below.
                stop_is_reference_terminator=False,
                contig=contig,
                strand="+" if locus[0].direction > 0 else "-",
                footprint_1based=(foot_low, foot_high), n_hsp=len(locus),
            ),
            refusal_reason=(
                f"frame_not_anchorable: the best locus on {contig} covers "
                f"reference codons {ref_low}-{ref_high} but no HSP reaches the "
                f"reference's first codon, so the isolate's reading frame cannot "
                f"be anchored and 'is this stop the reference's terminator' "
                f"cannot be asked. Refusing rather than answering it from a "
                f"codon index."
            ),
        )

    sequence = assembly.get(contig)
    if sequence is None:
        raise PipelineError(
            f"Contig {contig} is not in the assembly passed to "
            f"measure_locus_structure. The tblastn subject ids and the assembly "
            f"records must be the same release, or every coordinate read here "
            f"is meaningless.",
            contig=contig,
            sample_id=sample_id,
        )
    if foot_high > len(sequence) or foot_low < 1:
        raise PipelineError(
            f"The aligned locus on {contig} spans {foot_low}-{foot_high} but the "
            f"record is {len(sequence)} nt. Coordinates and sequence disagree.",
            contig=contig,
            sample_id=sample_id,
        )

    direction = locus[0].direction
    start = anchors[0].subject_start
    stop_coordinate = (max(h.subject_high for h in locus) if direction > 0
                       else min(h.subject_low for h in locus))
    span_low = min(start, stop_coordinate)
    span_high = max(start, stop_coordinate)
    coding = sequence[span_low - 1:span_high]
    if direction < 0:
        # The gene is on the minus strand, so its own reading frame runs along
        # the reverse complement of the aligned span, not along the forward
        # record. Translating the forward slice would read a frame that does not
        # exist in the isolate.
        coding = _reverse_complement(coding)
    protein = _translate(coding)
    stop_index = protein.find("*")
    net_frame_change_nt = len(coding) - 3 * (ref_high - ref_low + 1)
    frameshift = net_frame_change_nt % 3 != 0
    frame_is_native = not frameshift

    # THE TERMINATOR TEST. Three statements about the alignment, and no codon
    # index. See the section header above.
    #
    # The first statement is satisfiable only because the tblastn query is
    # `reference_codon_count(reference_cds)` residues long - written by
    # `write_tblastn_query`. `qend <= qlen` always, so a query one residue
    # shorter than `reference_codons` can never reach the reference's last codon
    # and this test reads False for every isolate, including an intact one.
    reaches_reference_terminator = ref_high == reference_codons
    runs_to_the_aligned_end = stop_index == len(protein) - 1
    stop_is_reference_terminator = bool(
        stop_index is not None
        and reaches_reference_terminator
        and frame_is_native
        and runs_to_the_aligned_end
    )

    at_contig_edge = foot_low == 1 or foot_high == len(sequence)
    evidence = LocusEvidence(
        identity_pct=identity_pct,
        coverage_pct=coverage_pct,
        orf_aa_length=stop_index,
        reference_aa_length=reference_codons,
        internal_stop_codon=(coding[3 * stop_index:3 * stop_index + 3]
                             if stop_index is not None else None),
        internal_stop_position_aa=(stop_index + 1
                                  if stop_index is not None else None),
        stop_is_reference_terminator=stop_is_reference_terminator,
        frameshift=frameshift,
        # The aligned span is a search artefact, not an annotation split. A
        # locus split across two abutting CDS by a draft-assembly gap is
        # recorded from the isolate GFF, and either way it is not a lesion.
        split_abutting_fragments=False,
        at_contig_edge=at_contig_edge,
        contig=contig,
        strand="+" if direction > 0 else "-",
        footprint_1based=(foot_low, foot_high),
        n_hsp=len(locus),
    )

    indels = indels_from_alignment(locus)
    group = compensating_indel_group(indels)
    return LocusMeasurement(
        sample_id=sample_id,
        evidence=evidence,
        indels=indels,
        net_frame_change_nt=net_frame_change_nt,
        compensating_group=group,
        repeat_motif=repeat_context(reference_cds, indels),
        frame_anchored=True,
    )


def structural_call(
    sample_id: str,
    hits: Sequence[AlignedNucleotideHit],
    *,
    reference_cds: str,
    assembly: Mapping[str, str],
    search: Optional[LocusSearch] = None,
    min_identity_pct: float = STRUCTURAL_MIN_IDENTITY_PCT,
    min_coverage_pct: float = STRUCTURAL_MIN_COVERAGE_PCT,
) -> StructuralCall:
    """Measure, then classify. The production entry point for one isolate.

    Two refusals happen **here**, before :func:`classify_structure` is reached,
    and both exist so no production path can reach
    :attr:`LocusEvidence.is_terminal_stop`'s codon-index fallback:

    * the measurement refused (no candidate hit, an ambiguous locus, mixed
      strands) - carried through verbatim;
    * the reading frame could not be anchored - ``frame_not_anchorable``.

    Neither can be ``absent``. Only a recorded search with zero hits can be.
    """
    measurement = measure_locus_structure(
        sample_id, hits, reference_cds=reference_cds, assembly=assembly,
        min_identity_pct=min_identity_pct,
    )
    # `no_candidate_hit` is the ONE refusal handed to the classifier rather
    # than answered here, because it is the case a recorded empty search has to
    # reach: only `classify_structure` may emit `absent`, and only from a
    # recorded search with no hits. Every other refusal - an ambiguous locus,
    # mixed strands, an unanchorable frame - is answered here, with its own
    # name, so none of them can be re-read as `absent` by a caller that goes on
    # to read `verdict` loosely.
    refused = measurement.refusal_reason is not None \
        and not measurement.refusal_reason.startswith("no_candidate_hit")
    if refused:
        return StructuralCall(
            sample_id=sample_id, verdict=StructuralVerdict.NOT_ASSESSED,
            reason=measurement.refusal_reason,
            evidence=measurement.evidence, search=search,
        )
    return classify_structure(
        sample_id, measurement.evidence, search=search,
        min_identity_pct=min_identity_pct, min_coverage_pct=min_coverage_pct,
    )


def structural_call_from_paths(
    sample_id: str,
    *,
    reference_cds_path: Path,
    tblastn_table_path: Path,
    assembly_path: Path,
    search: Optional[LocusSearch] = None,
    min_identity_pct: float = STRUCTURAL_MIN_IDENTITY_PCT,
) -> StructuralCall:
    """``structural_call`` over files, with every missing input named.

    **A missing input is ``not_assessed`` and can never be ``absent``.** If the
    reference, the tblastn table or the assembly is not on disk this returns
    ``not_assessed`` with a reason naming the file, and it does so *before* the
    search object is consulted. That ordering is the point: a recorded search
    with zero hits is the one thing that may yield ``absent``, so a run that
    never read its table must not be able to reach that branch - which is
    precisely how a missing file would otherwise become a fabricated negative.
    """
    for path, what in (
        (Path(reference_cds_path), "reference oprD CDS"),
        (Path(tblastn_table_path), "tblastn table"),
        (Path(assembly_path), "isolate assembly"),
    ):
        if not Path(path).is_file():
            return StructuralCall(
                sample_id=sample_id, verdict=StructuralVerdict.NOT_ASSESSED,
                reason=f"missing_input: the {what} was not found at {path}. An "
                       f"isolate that was never searched has not been shown to "
                       f"lack oprD; it has been shown not to have been "
                       f"assessed.",
                evidence=LocusEvidence(), search=search,
            )
    reference_cds = "".join(
        line.strip()
        for line in Path(reference_cds_path).read_text(
            encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()
    table = Path(tblastn_table_path).read_text(encoding="utf-8",
                                               errors="replace")
    hits = parse_tblastn_table(
        table, query_codon_length=reference_codon_count(reference_cds))
    return structural_call(
        sample_id, hits, reference_cds=reference_cds,
        assembly=read_assembly(Path(assembly_path)), search=search,
        min_identity_pct=min_identity_pct,
    )
