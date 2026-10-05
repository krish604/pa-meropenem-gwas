"""FASTA reading and assembly summary statistics.

Pure Python and streaming: an 8 Mb *P. aeruginosa* genome is read in one
pass with O(1) memory beyond the sequence being accumulated, and no third
party bioinformatics dependency is required for stage 1.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Set

from ..errors import DataContractError
from ..logging_utils import get_logger

LOGGER = get_logger("io.fasta")

#: IUPAC ambiguity codes counted as "ambiguous bases".
AMBIGUITY_CODES: Set[str] = set("NRYKMSWBDHVX*-")
VALID_NUCLEOTIDES: Set[str] = set("ACGTUN")


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA entry."""

    header: str
    sequence: str

    @property
    def record_id(self) -> str:
        """First whitespace-delimited token of the header."""
        return self.header.split()[0] if self.header.split() else ""

    def __len__(self) -> int:
        return len(self.sequence)


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    """Stream a FASTA file record by record.

    The file is read as bytes and decoded with explicit error reporting, so a
    truncated or garbled file produces a :class:`DataContractError` naming
    the file and byte offset, rather than an opaque ``UnicodeDecodeError``.
    Real assembly files in the wild are frequently truncated mid-transfer,
    and one bad file must not abort a cohort-wide run.

    Raises:
        DataContractError: File missing, empty, not decodable as FASTA, or
            the first record has no header line.
    """
    path = Path(path)
    if not path.exists():
        raise DataContractError("FASTA file not found", path=str(path))
    if path.stat().st_size == 0:
        raise DataContractError("FASTA file is empty (0 bytes)", path=str(path))

    header: Optional[str] = None
    chunks: List[str] = []
    offset = 0

    with path.open("rb") as handle:
        for raw in handle:
            line_bytes = len(raw)
            try:
                line = raw.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise DataContractError(
                    "FASTA file contains non-text bytes; it is truncated or "
                    "corrupt",
                    path=str(path),
                    byte_offset=offset + exc.start,
                    bad_byte=hex(raw[exc.start]),
                ) from exc
            offset += line_bytes
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header, "".join(chunks))
                header = line[1:].strip()
                if not header:
                    raise DataContractError(
                        "FASTA record has an empty header",
                        path=str(path),
                    )
                chunks = []
            else:
                if header is None:
                    raise DataContractError(
                        "FASTA file does not start with a header line",
                        path=str(path),
                    )
                chunks.append(line.upper())

    if header is not None:
        yield FastaRecord(header, "".join(chunks))
    if header is None:
        raise DataContractError(
            "FASTA file contains no sequence records", path=str(path)
        )


def file_integrity(path: Path) -> Dict[str, object]:
    """Cheap structural check of an assembly file.

    Does not parse the whole genome. Returns a dict with:

    ``exists``      - the path is present
    ``size_bytes``  - file size, or ``None`` when absent
    ``is_empty``    - zero bytes
    ``first_line``  - the first line, truncated for reporting
    ``looks_like_fasta`` - starts with ``>``
    ``is_gzipped``  - gzip magic number
    """
    path = Path(path)
    if not path.exists():
        return {
            "exists": False,
            "size_bytes": None,
            "is_empty": None,
            "first_line": None,
            "looks_like_fasta": None,
            "is_gzipped": None,
        }
    size = path.stat().st_size
    with path.open("rb") as handle:
        head = handle.read(512)
    return {
        "exists": True,
        "size_bytes": size,
        "is_empty": size == 0,
        "first_line": head[:200].decode("utf-8", errors="replace").splitlines()[0]
        if head
        else None,
        "looks_like_fasta": head.startswith(b">"),
        "is_gzipped": head.startswith(b"\x1f\x8b"),
    }


def count_bases(sequence: str) -> "Counter":
    """Count base composition. Returns a ``collections.Counter``."""
    from collections import Counter

    return Counter(sequence)


def assembly_stats(path: Path) -> "AssemblySummary":
    """Compute stage 1 metrics in a single streaming pass.

    Returns an :class:`AssemblySummary` with ``None`` for any metric that
    could not be computed (empty file, no records) rather than 0, so that
    "not measured" is never confused with "measured as zero".
    """
    lengths: List[int] = []
    total = 0
    gc = 0
    ambiguous = 0
    n_records = 0
    non_acgtn = 0

    for record in read_fasta(path):
        seq = record.sequence
        n_records += 1
        lengths.append(len(seq))
        total += len(seq)
        gc += seq.count("G") + seq.count("C")
        for base in seq:
            if base in AMBIGUITY_CODES:
                ambiguous += 1
            elif base not in VALID_NUCLEOTIDES:
                non_acgtn += 1

    return AssemblySummary(
        path=str(path),
        assembly_size=total if n_records else None,
        contig_count=n_records or None,
        n50=_nxx(lengths, 0.50),
        n90=_nxx(lengths, 0.90),
        largest_contig=max(lengths) if lengths else None,
        gc_content=round(100.0 * gc / total, 4) if total else None,
        ambiguous_bases=ambiguous if n_records else None,
        invalid_characters=non_acgtn if n_records else None,
    )


@dataclass(frozen=True)
class AssemblySummary:
    """Raw assembly metrics. Interpretation happens in stage 1."""

    path: str
    assembly_size: Optional[int]
    contig_count: Optional[int]
    n50: Optional[int]
    n90: Optional[int]
    largest_contig: Optional[int]
    gc_content: Optional[float]
    ambiguous_bases: Optional[int]
    invalid_characters: Optional[int] = None

    def to_row(self) -> dict:
        return {
            "path": self.path,
            "assembly_size": self.assembly_size,
            "contig_count": self.contig_count,
            "n50": self.n50,
            "n90": self.n90,
            "largest_contig": self.largest_contig,
            "gc_content": self.gc_content,
            "ambiguous_bases": self.ambiguous_bases,
            "invalid_characters": self.invalid_characters,
        }


def _nxx(lengths: Sequence[int], fraction: float) -> Optional[int]:
    """Return the NXX contig length. ``None`` for an empty assembly."""
    if not lengths:
        return None
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    if total == 0:
        return None
    threshold = total * fraction
    cumulative = 0
    for length in ordered:
        cumulative += length
        if cumulative >= threshold:
            return length
    return ordered[-1]


def summarise_alignment(path: Path) -> dict:
    """Summarise an alignment: number of sequences, sites, and missing data.

    An empty or absent alignment is reported as an empty summary rather than
    raising: a stage may legitimately find no alignment, and that is a result
    to report, not a crash.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {"n_sequences": 0, "n_sites": 0, "missing_fraction": None}

    records = list(read_fasta(path))
    if not records:
        return {"n_sequences": 0, "n_sites": 0, "missing_fraction": None}
    lengths = {len(r.sequence) for r in records}
    n_sites = max(lengths)
    if len(lengths) > 1:
        LOGGER.warning(
            "Alignment has ragged sequence lengths: %s", sorted(lengths)
        )
    total_sites = sum(len(r.sequence) for r in records)
    gap_chars = sum(
        r.sequence.count("-") + r.sequence.count("?") + r.sequence.count("N")
        for r in records
    )
    return {
        "n_sequences": len(records),
        "n_sites": n_sites,
        "missing_fraction": round(gap_chars / total_sites, 6) if total_sites else None,
    }
