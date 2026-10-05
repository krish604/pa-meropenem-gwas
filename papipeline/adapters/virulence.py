"""Virulence screening by raw blastn against VFDB.

**Why raw blast rather than abricate.** abricate is the conventional wrapper, and
it has no installable osx-arm64 build: its dependency chain reaches
`perl-lwp-protocol-https -> perl-test-requiresinternet -> perl-socket`, and
`perl-socket` has no osx-arm64 build in any channel checked. That is not a
wrinkle to work around - but abricate's own dependency list is
`blast >=2.7`, so blast *is* the mechanism abricate uses. This calls blastn
directly and parses VFDB's headers, which carry the factor, gene and category
that abricate would otherwise have supplied.

There is no `abricate-vfdb` conda package either (checked across bioconda
osx-arm64 and osx-64), so the sequences are provisioned under the overlay's
`db_root`. Nothing here downloads anything: `allow_database_update` is false and
a missing database is a preflight failure, not a cue to fetch.

**Carriage only.** A hit says a factor is present, nothing about whether it is
expressed, functional, or reached a phenotype. The stage is deliberately
independent of the AMR stage and consults no mechanism tables.
"""

from __future__ import annotations

import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from ..errors import PipelineError, ToolExecutionError
from ..logging_utils import get_logger

#: `[Type IV pili (VF0082) - Adherence (VFC0001)]`
#:
#: Anchored on the `(VFnnnn) - ` token rather than on a bare hyphen. A first
#: attempt split on `\s*-\s*`, which also matches the hyphen *inside* a factor
#: name: `HSI-2 (VF0943) - Effector delivery system (VFC0086)` became
#: gene `HSI`, category `2 (VF0943) - Effector delivery system (VFC0086)`. That
#: is silently wrong rather than a parse failure, and it put a fabricated
#: category in a scientific output.
_BRACKETED = re.compile(r"\[([^\[\]]+)\]")
_FACTOR_BLOCK = re.compile(
    r"^\s*(?P<factor>.+?)\s+\((?P<vf_id>VF\d+)\)\s+-\s+"
    r"(?P<category>.+?)\s+\((?P<category_id>VFC\d+)\)\s*$"
)


LOGGER = get_logger("adapters.virulence")


class VirulencePreflightError(PipelineError):
    """The database or the tool is missing. Named so the cause is unambiguous."""


@dataclass(frozen=True)
class VfdbEntry:
    """One VFDB sequence and the annotation abricate would have supplied."""

    accession: str
    gene: str
    factor: str
    category: str
    vf_id: str
    category_id: str
    description: str


@dataclass(frozen=True)
class BlastHit:
    accession: str
    identity_pct: float
    alignment_length: int
    subject_length: int
    evalue: float

    @property
    def coverage_pct(self) -> float:
        """Percent of the *factor* covered, not of the genome."""
        if self.subject_length <= 0:
            return 0.0
        return round(100.0 * self.alignment_length / self.subject_length, 2)


def parse_vfdb_headers(path: Path) -> Dict[str, VfdbEntry]:
    """Map accession -> annotation, reading headers only.

    Headers are `vfdb~~~<locus>~~~<refseq>~~~ (<locus>) <product>
    [<factor> (VFnnnn) - <category> (VFCnnnn)] [<organism>]`.

    The sequences themselves are never loaded: `makeblastdb` reads the file, and
    holding 4,592 sequences in memory to read their titles would be waste. A
    header that cannot be parsed is skipped rather than guessed at - a factor
    with no gene would otherwise be reported as an anonymous carriage.
    """
    entries: Dict[str, VfdbEntry] = {}
    skipped: List[str] = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            fields = header.split()
            if not fields:
                continue
            accession = fields[0]
            # Match the annotation block specifically. Anchoring the whole
            # header against `$` fails, because a second bracketed field
            # ([organism]) follows the annotation.
            block = None
            for candidate in _BRACKETED.finditer(header):
                block = _FACTOR_BLOCK.match(candidate.group(1))
                if block is not None:
                    break
            if block is None:
                # No annotation means no factor name and no category, and this
                # stage reports carriage of named factors. An unnamed hit would
                # be a row a reader cannot interpret, so it is skipped here and
                # counted, rather than invented.
                skipped.append(accession)
                continue
            # VFDB's second header field is the locus/gene symbol (`fimV`,
            # `PA1663`) and is the only per-sequence gene identifier present.
            # The locus symbol is the second `~~~`-separated component of the
            # accession (`vfdb~~~PA1663~~~NP_250354~~~` -> `PA1663`), which is
            # reliable across every header shape. Taking it from the prose
            # fields instead yields the product description - "unknown",
            # "transcriptional" - because the second field is `(locus)` for some
            # entries and a product word for others.
            parts = accession.split("~~~")
            gene = parts[1].strip("*") if len(parts) > 1 else ""
            entries[accession] = VfdbEntry(
                accession=accession,
                gene=gene,
                factor=block.group("factor").strip(),
                category=block.group("category").strip(),
                vf_id=block.group("vf_id"),
                category_id=block.group("category_id"),
                description=header,
            )
    if not entries:
        raise VirulencePreflightError(
            "VFDB sequences file contains no parsable factor annotations, so "
            f"nothing could be screened ({len(skipped)} headers unparsed): {path}",
            path=str(path),
        )
    LOGGER.info(
        "VFDB: %d annotated factors%s", len(entries),
        f", {len(skipped)} headers skipped as unannotated" if skipped else "",
    )
    return entries


def blast_db_path(sequences: Path) -> Path:
    """`makeblastdb` writes `<name>.nin` beside the input.

    Built by string concatenation rather than `with_suffix`, which would
    strip a dotted stem as though the tail were a file extension - the same
    trap that made stage 7 look for `PDT000034122.delta`.
    """
    return Path(f"{sequences}.nin")


def ensure_blast_db(
    sequences: Path,
    *,
    makeblastdb: str = "makeblastdb",
    runner=None,
) -> Path:
    """Build the BLAST index if absent; return the database name to search.

    Returns the *stem*, not the `.nin` file. `blastn -db` takes a database name
    and locates the index itself; handing it `sequences.nin` directly fails with
    "No alias or index file found ... in search path", which reads like a
    provisioning fault rather than a naming one.

    Never rebuilds an existing index: rebuilding unconditionally would rewrite a
    database the operator pinned, and the pin is the point of putting it under
    `db_root`.
    """
    index = blast_db_path(sequences)
    if index.exists():
        return Path(sequences)
    command = [makeblastdb, "-in", str(sequences), "-dbtype", "nucl"]
    result = _run(command, runner=runner, what="makeblastdb")
    if result.returncode != 0:
        raise VirulencePreflightError(
            "makeblastdb failed, so VFDB cannot be searched: "
            f"{result.stderr.strip()[:400]}",
            command=" ".join(command),
        )
    if not index.exists():
        raise VirulencePreflightError(
            f"makeblastdb reported success but wrote no index at {index}",
            path=str(index),
        )
    # Return the database name as given, not via with_suffix: blastn locates the
    # index itself and a dotted stem would lose its tail.
    return Path(sequences)


def blastn_command(
    *,
    program: str,
    genome: Path,
    database: Path,
    threads: int,
    evalue: float,
) -> List[str]:
    """Tabular blastn against the VFDB index.

    `-outfmt 6 qseqid sacc stitle pident length slen evalue`. The subject
    length is `slen`, **not** `qlen`. In blastn's tabular output `qlen` is
    the length of the *query* - here a whole 1.2 Mb genome - so computing
    coverage from it yields ~0.1% for every hit and a 90%% coverage
    threshold silently rejects all 4,592 factors. That failure is invisible:
    blastn exits 0 and reports zero findings, which reads as "this isolate
    carries no virulence factors" rather than as a broken screen. `slen` is
    the factor's own length, so coverage means what it says.
    """
    return [
        program,
        "-query", str(genome),
        "-db", str(database),
        "-outfmt", "6 qseqid sacc stitle pident length slen evalue",
        "-evalue", str(evalue),
        "-num_threads", str(threads),
    ]


def parse_blast_table(text: str) -> List[BlastHit]:
    hits: List[BlastHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        try:
            hits.append(
                BlastHit(
                    accession=parts[1],
                    identity_pct=float(parts[3]),
                    alignment_length=int(parts[4]),
                    subject_length=int(parts[5]),
                    evalue=float(parts[6]),
                )
            )
        except ValueError:
            # A malformed row is dropped rather than aborting the screen: blast
            # emitting an unexpected column would otherwise lose the whole
            # isolate's carriage.
            continue
    return hits


def select_hits(
    hits: Iterable[BlastHit],
    entries: Dict[str, VfdbEntry],
    *,
    min_identity_pct: float,
    min_coverage_pct: float,
) -> List[BlastHit]:
    """One best hit per VFDB factor, passing both thresholds.

    One hit per factor, because VFDB holds several sequences under a single
    factor - `Type IV pili (VF0082)` covers fimV, xcpQ, algG and others - so
    without the collapse a pilus locus contributes several carriage rows and a
    count of factors is really a count of sequences. Keying `best` by factor is
    what enforces it; an earlier version truncated the whole list instead, which
    reported exactly one factor for every genome.

    Best by identity then coverage, so a full-length exact match is preferred
    over a fragment of the same factor.
    """
    best: Dict[str, Tuple[str, BlastHit]] = {}
    for hit in hits:
        if hit.identity_pct < min_identity_pct:
            continue
        if hit.coverage_pct < min_coverage_pct:
            continue
        entry = entries.get(hit.accession)
        if entry is None:
            continue
        key = entry.factor
        current = best.get(key)
        if current is None or (hit.identity_pct, hit.coverage_pct) > (
            current[1].identity_pct,
            current[1].coverage_pct,
        ):
            best[key] = (hit.accession, hit)

    chosen: List[BlastHit] = [hit for _key, (_accession, hit) in best.items()]
    chosen.sort(
        key=lambda h: (
            entries[h.accession].category if h.accession in entries else "",
            entries[h.accession].factor if h.accession in entries else "",
            h.accession,
        )
    )
    return chosen


def screen_isolate(
    *,
    genome: Path,
    entries: Dict[str, VfdbEntry],
    database: Path,
    program: str,
    threads: int,
    evalue: float,
    min_identity_pct: float,
    min_coverage_pct: float,
    database_name: str,
    database_version: str,
    runner=None,
) -> List[Dict[str, object]]:
    """Every passing VFDB factor in one genome, as carriage rows."""
    result = _run(
        blastn_command(
            program=program,
            genome=genome,
            database=database,
            threads=threads,
            evalue=evalue,
        ),
        runner=runner,
        what="blastn",
    )
    if result.returncode != 0:
        raise ToolExecutionError(
            f"blastn failed for {genome.name}: {result.stderr.strip()[:400]}",
            command=" ".join(result.command),
        )
    hits = parse_blast_table(result.stdout)

    # Checked before selection, because `select_hits` groups by factor and has
    # to skip an accession it cannot resolve. Checking afterwards would mean the
    # skip had already swallowed the evidence, and an index out of step with the
    # headers would under-count carriage with nothing raised.
    unknown = sorted({h.accession for h in hits if h.accession not in entries})
    if unknown:
        raise VirulencePreflightError(
            f"blastn reported {len(unknown)} subject(s) absent from the VFDB "
            f"header set (first: {unknown[0]!r}); the BLAST index and the "
            "sequences file are out of step. Rebuild with ensure_blast_db().",
            accession=unknown[0],
        )

    rows: List[Dict[str, object]] = []
    for hit in select_hits(
        hits,
        entries,
        min_identity_pct=min_identity_pct,
        min_coverage_pct=min_coverage_pct,
    ):
        entry = entries[hit.accession]
        rows.append(
            {
                "sample_id": "",
                "virulence_factor": entry.factor,
                "gene": entry.gene,
                "category": entry.category,
                "vf_id": entry.vf_id,
                "category_id": entry.category_id,
                "database": database_name,
                "database_version": database_version,
                "confidence": None,
                "identity_pct": hit.identity_pct,
                "coverage_pct": hit.coverage_pct,
                "vfdb_accession": entry.accession,
            }
        )
    return rows


def _run(command: Sequence[str], *, runner, what: str):
    if runner is not None:
        return runner(command=command)
    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command), capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise ToolExecutionError(
            f"could not launch {what}",
            command=" ".join(command),
            error=str(exc),
        ) from exc
    from .external import CommandResult

    return CommandResult(
        command=list(command),
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        duration_s=round(time.monotonic() - started, 3),
    )


def database_version_on_disk(sequences: Path) -> str:
    """Version the screen actually used.

    From the file's own digest rather than a hard-coded string: VFDB is
    re-released, and a pinned value that no longer describes the bytes on disk
    is worse than no pin.
    """
    import hashlib

    digest = hashlib.sha256()
    with Path(sequences).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()[:16]}"