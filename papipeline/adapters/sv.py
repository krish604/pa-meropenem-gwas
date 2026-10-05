"""Assembly-based structural-variant calling from `show-coords` output.

The caller is deliberately narrow. This pipeline's only input is assembled
genomes, so there is no read depth, no coverage and no split-read or
discordant-pair evidence - the evidence class `sniffles`, `delly` and `gridss`
all require, and which none of them can substitute for. What an
assembly-to-reference alignment does support is a *breakpoint with aligned
flanks*: two confidently-matched blocks on the same contig with a gap between
them. Everything here is that and nothing more.

Two failure modes drove the design, and both produce a plausible-looking table
rather than an error:

* **Assembly fragmentation read as variation.** A draft assembly breaks wherever
  a read spans a repeat, so unaligned query sequence is routine. Calling it an
  insertion would report contig breaks as structural variation - in this cohort,
  hundreds of them. Unaligned contig ends are therefore `not_assessable`, which
  is what scientific rule 7 already names a contig gap, and emit no call.
* **`confirmed` reached without evidence.** Nothing here can satisfy rule 7's
  bar for `confirmed`, so every call is a candidate. `confirmed_only()` is empty
  for every REAL cohort by construction, and that is a property of the input,
  not a bug to be tuned away.
"""

from __future__ import annotations

import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from ..errors import PipelineError, ToolExecutionError
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.sv")


class SvPreflightError(PipelineError):
    """A tool or the reference is missing. Named so the cause is unambiguous."""


@dataclass(frozen=True)
class AlignmentBlock:
    """One aligned block from `show-coords -T -l -H -q`.

    Column positions verified against real nucmer 4.0.1 output, not assumed:

    1 S1 2 E1 3 S2 4 E2 5 LEN1 6 LEN2 7 %IDY 8 LEN R 9 LEN Q 10 ref 11 query
    """

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    identity_pct: float
    ref_name: str
    query_name: str

    @property
    def query_length(self) -> int:
        return self.query_end - self.query_start + 1

    @property
    def ref_length(self) -> int:
        return self.ref_end - self.ref_start + 1


@dataclass(frozen=True)
class SvCandidate:
    sample_id: str
    variant_id: str
    variant_type: str
    position: Optional[int]
    size: int
    evidence: str
    call_status: str
    affected_gene: Optional[str] = None
    confidence: Optional[str] = None
    mge: Optional[str] = None


def nucmer_command(
    *,
    program: str,
    reference: Path,
    query: Path,
    prefix: Path,
    minmatch: int,
    mincluster: int,
    maxgap: int,
    maxmatch: bool,
    threads: int,
) -> List[str]:
    """Align one assembly against the reference.

    Every flag is passed explicitly rather than left at its default. A default
    that changes between mummer releases would silently change which variants
    are found, and a pin in the overlay is supposed to mean something.
    """
    command = [
        program,
        "--minmatch", str(minmatch),
        "--mincluster", str(mincluster),
        "--maxgap", str(maxgap),
        "--threads", str(threads),
        "--prefix", str(prefix),
    ]
    if maxmatch:
        command.append("--maxmatch")
    command.extend([str(reference), str(query)])
    return command


def show_coords_command(*, program: str, delta: Path) -> List[str]:
    """Tab-delimited, with sequence lengths, no header, sorted by query.

    `-q` is load-bearing: gaps are found *within* one query contig, so the
    blocks must arrive grouped and ordered by query coordinate. `-r`, the
    default, sorts by reference and interleaves contigs.
    """
    return [program, "-T", "-l", "-H", "-q", str(delta)]


def parse_show_coords(text: str) -> List[AlignmentBlock]:
    """Parse `show-coords -T -l -H -q` output.

    Malformed lines are dropped rather than fatal: show-coords pads short rows,
    and one truncated line should not cost an isolate's whole screen.
    """
    blocks: List[AlignmentBlock] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        try:
            blocks.append(
                AlignmentBlock(
                    ref_start=int(parts[0]),
                    ref_end=int(parts[1]),
                    query_start=int(parts[2]),
                    query_end=int(parts[3]),
                    identity_pct=float(parts[6]),
                    ref_name=parts[9].strip(),
                    query_name=parts[10].strip(),
                )
            )
        except (ValueError, IndexError):
            continue
    return blocks


def call_from_gaps(
    blocks: Sequence[AlignmentBlock],
    *,
    sample_id: str,
    min_sv_size: int,
    min_flanking_identity_pct: float,
) -> List[SvCandidate]:
    """Candidates from gaps between consecutive blocks on one query contig.

    A gap is only reported when the two flanking alignments are both confident
    and the gap is lopsided - one coordinate advances and the other does not.
    A gap large in *both* is a repeat or a contig join, not a variant, and a gap
    with both advances small is ordinary divergence.

    With `query_gap >> ref_gap` the assembly has sequence the reference lacks
    (insertion). With `ref_gap >> query_gap` the assembly lacks sequence the
    reference has (deletion).
    """
    by_contig: Dict[str, List[AlignmentBlock]] = {}
    for block in blocks:
        by_contig.setdefault(block.query_name, []).append(block)

    candidates: List[SvCandidate] = []
    for contig, contig_blocks in sorted(by_contig.items()):
        ordered = sorted(contig_blocks, key=lambda b: b.query_start)
        for left, right in zip(ordered, ordered[1:]):
            if left.identity_pct < min_flanking_identity_pct:
                continue
            if right.identity_pct < min_flanking_identity_pct:
                continue
            query_gap = right.query_start - left.query_end - 1
            ref_gap = right.ref_start - left.ref_end - 1
            if query_gap < 0 or ref_gap < 0:
                # The blocks overlap in at least one frame; nucmer has already
                # called this one alignment.
                continue

            # A gap may legitimately be zero in one frame - that *is* what a
            # lopsided gap looks like, and requiring both to be positive would
            # discard every true insertion and deletion.
            if query_gap >= min_sv_size and query_gap > 4 * ref_gap:
                variant_type, size, side = "insertion", query_gap, "assembly"
            elif ref_gap >= min_sv_size and ref_gap > 4 * query_gap:
                variant_type, size, side = "deletion", ref_gap, "reference"
            else:
                # Lopsidedness is what separates a variant from a repeat. A gap
                # large in both coordinates is neither, and reporting it would
                # be reporting the aligner's uncertainty.
                continue

            candidates.append(
                SvCandidate(
                    sample_id=sample_id,
                    variant_id=(
                        f"{sample_id}:{contig}:{left.query_end + 1}:{variant_type}"
                    ),
                    variant_type=variant_type,
                    position=left.query_end + 1,
                    size=size,
                    evidence=(
                        f"nucmer:{variant_type}_in_{side};flanking_identity="
                        f"{min(left.identity_pct, right.identity_pct):.1f};"
                        f"ref_gap={ref_gap};query_gap={query_gap}"
                    ),
                    call_status="candidate",
                )
            )
    return candidates


def unassessable_regions(
    blocks: Sequence[AlignmentBlock],
    *,
    query_contig_lengths: Dict[str, int],
    min_sv_size: int,
) -> List[tuple]:
    """Unaligned contig ends, reported as `not_assessable` and never as calls.

    A draft assembly breaks at repeats, so unaligned sequence at a contig end
    carries no information about what the reference does or does not have. Rule
    7 already has a state for this - "the region could not be evaluated (contig
    gap, no coverage)" - and conflating it with "no variant present" is what it
    forbids.
    """
    by_contig: Dict[str, List[AlignmentBlock]] = {}
    for block in blocks:
        by_contig.setdefault(block.query_name, []).append(block)

    regions: List[tuple] = []
    for contig, length in sorted(query_contig_lengths.items()):
        contig_blocks = by_contig.get(contig) or []
        if not contig_blocks:
            regions.append((contig, 0, length))
            continue
        start = min(b.query_start for b in contig_blocks)
        end = max(b.query_end for b in contig_blocks)
        if start - 1 >= min_sv_size:
            regions.append((contig, 1, start - 1))
        if length - end >= min_sv_size:
            regions.append((contig, end + 1, length))
    return regions


def contig_lengths(fasta: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    with Path(fasta).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths.setdefault(name, 0)
            elif name:
                lengths[name] += len(line.strip())
    return lengths


def screen_assembly(
    *,
    sample_id: str,
    genome: Path,
    reference: Path,
    workdir: Path,
    nucmer: str,
    show_coords: str,
    minmatch: int,
    mincluster: int,
    maxgap: int,
    maxmatch: bool,
    threads: int,
    min_sv_size: int,
    min_flanking_identity_pct: float,
    runner=None,
) -> tuple:
    """Candidates and unassessable-region count for one assembly."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    # Not `workdir / sample_id` with `.with_suffix()`: isolate IDs contain
    # dots (`PDT000034122.1`), so with_suffix(".delta") would strip `.1` as
    # though it were a file extension and look for `PDT000034122.delta`.
    # It reported "nucmer wrote no delta" on a successful alignment.
    prefix = workdir / sample_id

    _run(
        nucmer_command(
            program=nucmer,
            reference=reference,
            query=genome,
            prefix=prefix,
            minmatch=minmatch,
            mincluster=mincluster,
            maxgap=maxgap,
            maxmatch=maxmatch,
            threads=threads,
        ),
        runner=runner,
        what="nucmer",
    )
    delta = Path(f"{prefix}.delta")
    if not delta.exists():
        raise SvPreflightError(
            f"nucmer wrote no delta file for {sample_id} at {delta}; "
            "alignment produced nothing to call",
            path=str(delta),
        )

    result = _run(
        show_coords_command(program=show_coords, delta=delta),
        runner=runner,
        what="show-coords",
    )
    if result.returncode != 0:
        raise ToolExecutionError(
            f"show-coords failed for {sample_id}: {result.stderr.strip()[:400]}",
            command=" ".join(result.command),
        )

    blocks = parse_show_coords(result.stdout)
    candidates = call_from_gaps(
        blocks,
        sample_id=sample_id,
        min_sv_size=min_sv_size,
        min_flanking_identity_pct=min_flanking_identity_pct,
    )
    regions = unassessable_regions(
        blocks,
        query_contig_lengths=contig_lengths(genome),
        min_sv_size=min_sv_size,
    )
    return candidates, regions


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