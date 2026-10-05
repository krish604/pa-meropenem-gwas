"""Running the tblastn search the oprD structural screen needs.

**Why this module exists.** :mod:`papipeline.adapters.oprd_locus` can decide
whether an isolate's oprD locus is WHOLE - ``intact`` or ``disrupted`` - but it
can only do that from an aligned-nucleotide HSP table, and no run has ever
produced one. Round 12's WIRE.md section A1 established that with a runtime
demonstration: ``structural_call_from_paths`` wants
``tblastn_table_path``, every isolate came back ``not_assessed`` with a reason
naming a file that was never written, and merging that into
``oprd_locus_resolution`` took all ten isolates from ``intact`` to
``not_assessed`` while adding no information. The decision logic was correct and
unreachable. This module is the missing producer.

**tblastn is the search, not blastn.** The query is a protein - the pinned PAO1
PA0958 OprD sequence, derived by ``reference_protein`` from the pinned GFF and
FASTA - and the subject is an isolate's assembled nucleotides. A translated
search is what places an indel at a reference position, because the aligned
columns it returns (``qseq``/``sseq``) are the reading frame's columns.
``virulence.screen_isolate`` uses blastn because ITS query is DNA; the two are
not interchangeable and the choice here is not a preference.

**The database is built here, because ``tblastn -db`` cannot read a FASTA.**
Measured, not assumed - pointing ``-db`` at the assembly FASTA:

    $ tblastn -query .../PA0958_oprD.faa -db .../asm/PDT000034122.1.fna \
        -outfmt "6 qseqid sseqid pident length" -evalue 1e-3 \
        -max_target_seqs 5000 -db_gencode 11 -num_threads 4
    BLAST Database error: No alias or index file found for nucleotide database
    [.../asm/PDT000034122.1.fna] in search path

    exit status 2, zero rows on stdout.

So the signature takes ``assembly_fasta`` and this module formats it with
``makeblastdb -dbtype nucl`` before searching. That is a fact about BLAST+, and
it is the reason this file exists rather than a three-line wrapper around
``oprd_locus.tblastn_command``.

**The search is reproducible to the byte against the stored artifacts.** The
round-9 shell script (``pa-artifacts/oprd2-locus/run_tblastn.sh``) built the
database the same way and recorded every parameter; re-running it here over the
same assemblies reproduces all ten stored tables byte for byte. Two details of
those recorded parameters are load-bearing and are re-verified rather than
copied:

* ``-db_gencode 11`` - bacterial code. blast's default is 1, the standard code,
  which is wrong for *P. aeruginosa* and invents stops when translating the
  isolate's locus.
* ``-max_target_seqs 5000`` - ``0`` is rejected by blast 2.17 as legacy syntax
  (``expected >=1``). Not passed, not defaulted to zero.

The recorded ``.cmd`` files render the outfmt with escaped spaces
(``6\\ qseqid\\ sseqid\\ ...``), which is an artefact of ``printf '%q '`` on an
argv element that legitimately contains spaces, not a different outfmt. This
module passes ``oprd_locus.TBLASTN_OUTFMT`` as a single argv element, which is
the field list those files mean.

**``out_tsv`` is exactly what ``structural_call_from_paths`` consumes**, because
that function is the only consumer and it re-parses the file with
``parse_tblastn_table``. So the file is validated with that parser *before* it
is put in place: a table this module produced is either parseable or was never
written. An empty table is valid and means zero hits - a recorded search that
found nothing, which is the one condition under which ``oprd_locus`` will answer
``absent``. Writing a malformed table and discovering it three stages later would
turn a tool problem into a silent ``not_assessed``.

**The search is recorded on the result, not left in a log.**
``oprd_locus.LocusSearch`` is the only route by which an ``absent`` is reachable,
so :attr:`TblastnRunInfo.search` builds one from what actually ran: the argv, the
database prefix, the parameters, and the row count. A caller that passes
``result.search`` into ``structural_call_from_paths`` therefore cannot
manufacture an ``absent`` it did not measure.

Nothing here computes a statistic across isolates. Every number on this record
belongs to one isolate's own search.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from ..errors import PipelineError, ToolExecutionError, ToolNotAvailableError
from ..logging_utils import get_logger
from . import oprd_locus as oprd
from .external import CommandResult

LOGGER = get_logger("adapters.tblastn")

#: The two tools this adapter runs. Named together because the refusal a caller
#: sees has to be able to say which one is absent - a search that cannot be
#: formatted and a search that cannot be run are different fixes.
TBLASTN = "tblastn"
MAKEBLASTDB = "makeblastdb"

#: ``-version`` is the flag both accept; ``--version`` is blastn's alias and is
#: not relied on. Read from ``blastp -version`` on this machine, not guessed.
VERSION_FLAG = "-version"


# ---------------------------------------------------------------------------
# Tool location
# ---------------------------------------------------------------------------


def resolve_tool(name: str, search_dirs: Sequence[Path]) -> Path:
    """Locate ``name``, or refuse naming it.

    ``PATH`` first, then ``search_dirs`` in the order given, which is the order
    :meth:`papipeline.config.loader.MachineConfig.tool_candidates` documents:
    ``PATH`` is authoritative and the machine overlay's directories are the
    documented escape hatch for a tool that is not on it. Nothing else is
    consulted and no path is hard-coded - the same rule
    ``gubbins.resolve_binary`` states, for the same reason.

    Raises:
        ToolNotAvailableError: ``name`` is in neither ``PATH`` nor ``search_dirs``.
            The message names the tool and lists where it looked, because a
            refusal that names no tool leaves the operator guessing which of the
            two to install.
    """
    found = shutil.which(name)
    if found:
        return Path(found)
    for directory in search_dirs:
        candidate = Path(directory) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    raise ToolNotAvailableError(
        f"{name} is required to measure an isolate's oprD reading frame and was "
        f"not found. It looked on PATH and in "
        f"{[str(d) for d in search_dirs]}. "
        f"Install the NCBI BLAST+ {name} that the environment pins (verify the "
        f"version with `{name} {VERSION_FLAG}` rather than assuming one), or add "
        f"the directory holding it to paths.tool_search_dirs in this machine's "
        f"overlay. This is not the same as an isolate that lacks oprD: it means "
        f"the locus was never searched.",
        tool=name,
        searched=[str(d) for d in search_dirs],
    )


def tool_version(executable: Path, name: str = TBLASTN) -> str:
    """The version string the tool itself reports, first line.

    ``unknown`` rather than a refusal when it cannot be read. The version is
    provenance recorded on every result; failing a completed search because a
    provenance string was unavailable would discard real evidence over a
    cosmetic loss. ``gubbins.binary_version`` documents the same choice for the
    same reason.
    """
    try:
        completed = subprocess.run(
            [str(executable), VERSION_FLAG],
            capture_output=True, text=True, check=False, timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    text = (completed.stdout or "") + (completed.stderr or "")
    for line in text.splitlines():
        if line.strip():
            return line.strip()
    return "unknown"


# ---------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TblastnRunInfo:
    """What one tblastn search did, recorded so a negative result is auditable.

    ``oprd_locus.LocusSearch`` is the object that licenses an ``absent``, and it
    is reachable only from a recorded search. This record carries the same facts
    plus the things that make *this* run re-checkable: the exact argv, both
    tools' version strings, the sha256 of the table as written, and the
    ``makeblastdb`` argv that produced the database.

    :attr:`search` builds the ``LocusSearch``. Callers should pass it to
    ``structural_call_from_paths`` rather than constructing one by hand - a
    hand-built search can claim a database and a hit count that were never
    observed.
    """

    out_tsv: Path
    query: Path
    assembly: Path
    database: str
    command: List[str]
    parameters: Tuple[Tuple[str, str], ...]
    tool_version: str
    makeblastdb_command: List[str]
    makeblastdb_version: str
    n_hits: int
    table_sha256: str
    executable: Path = Path(TBLASTN)

    @property
    def tool(self) -> str:
        return TBLASTN

    @property
    def command_line(self) -> str:
        """The argv as one string, for the record and for re-running."""
        return " ".join(self.command)

    @property
    def search(self) -> oprd.LocusSearch:
        """The recorded search, in the shape ``oprd_locus`` requires."""
        return oprd.LocusSearch(
            program=TBLASTN,
            query=str(self.query),
            database=self.database,
            parameters=tuple(self.parameters),
            command=self.command_line,
            n_hits=self.n_hits,
        )

    def as_row(self) -> Dict[str, Any]:
        """One flat row, for a stage table that wants the provenance inline."""
        return {
            "tool": self.tool,
            "tool_version": self.tool_version,
            "executable": str(self.executable),
            "query": str(self.query),
            "assembly": str(self.assembly),
            "database": self.database,
            "outfmt": oprd.TBLASTN_OUTFMT,
            "parameters": " ".join(f"{k}={v}" for k, v in self.parameters),
            "command": self.command_line,
            "makeblastdb_version": self.makeblastdb_version,
            "n_hits": self.n_hits,
            "table_sha256": self.table_sha256,
        }


#: Suffix of the directory the derived BLAST database is written into, beside
#: the table it serves. It is named so that a reader of the work root can tell a
#: derived index from a pipeline output at a glance.
DATABASE_DIR_SUFFIX = ".blastdb"


def database_prefix(out_tsv: Path) -> Path:
    """Where the BLAST index for ``out_tsv``'s isolate is built.

    Beside the table, not beside the assembly: the assembly is an input under
    ``data/`` and writing an index next to it would mutate the input tree. The
    database is derived from the assembly and is safe to delete; the table is
    not, which is why they are named differently.

    ``blast_db_path`` in :mod:`papipeline.adapters.virulence` builds its name by
    string concatenation for the same reason this does not use ``with_suffix``:
    ``with_suffix`` replaces the last dotted segment, so a name like
    ``PDT000034122.1`` - which is what an isolate id gives a table called after
    the isolate - would lose its tail.
    """
    table = Path(out_tsv)
    return table.parent / f"{table.name}{DATABASE_DIR_SUFFIX}" / table.name


# ---------------------------------------------------------------------------
# The runner protocol
# ---------------------------------------------------------------------------


def _run(
    command: Sequence[str], *, runner: Optional[Callable[..., Any]], what: str
) -> CommandResult:
    """Execute ``command``, or hand it to an injected ``runner``.

    The injected protocol is :mod:`papipeline.adapters.virulence`'s: the runner is
    called as ``runner(command=command)`` and returns a
    :class:`papipeline.adapters.external.CommandResult`. That is how the whole
    suite tests a tool it does not have installed, and reusing it means a caller
    that already fakes blastn needs no second fake here.
    """
    if runner is not None:
        result = runner(command=command)
        if not isinstance(result, CommandResult):
            raise PipelineError(
                f"The injected runner returned {type(result).__name__} for "
                f"{what}, not a CommandResult. The house runner protocol is "
                f"`runner(command=...) -> CommandResult`; see "
                f"papipeline.adapters.virulence._run. A runner that returned "
                f"something else would be read as a tool failure with nothing "
                f"raised.",
                what=what,
            )
        return result
    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command), capture_output=True, text=True, check=False,
        )
    except OSError as exc:
        raise ToolExecutionError(
            f"could not launch {what}",
            command=" ".join(command), error=str(exc),
        ) from exc
    return CommandResult(
        command=list(command),
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        duration_s=round(time.monotonic() - started, 3),
    )


# ---------------------------------------------------------------------------
# The search
# ---------------------------------------------------------------------------


def run_tblastn(
    query_faa: Path,
    assembly_fasta: Path,
    out_tsv: Path,
    *,
    tool_search_dirs: Sequence[Path],
    threads: int = 1,
    runner: Optional[Callable[..., Any]] = None,
) -> TblastnRunInfo:
    """Search one isolate's assembly for the reference oprD protein. Write the table.

    ``out_tsv`` is the table ``oprd_locus.structural_call_from_paths`` reads:
    ``oprd_locus.TBLASTN_OUTFMT`` columns, tab-separated, written only after
    ``oprd_locus.parse_tblastn_table`` has accepted it. An empty file is a valid
    result meaning zero hits, and is the only shape that may later be read as
    ``absent``.

    Args:
        query_faa: The reference protein, FASTA. Written by
            ``oprd_locus.resolve_isolate`` as ``reference_protein.faa``; derived
            from the pinned PAO1 GFF and FASTA by ``reference_protein``.
        assembly_fasta: The isolate's assembled nucleotides. Formatted into a
            BLAST database here, because ``tblastn -db`` refuses a FASTA.
        out_tsv: Where the table is written. Parent directories are created.
        tool_search_dirs: Extra directories to search, in order, after ``PATH``.
            Normally ``MachineConfig.resolved_tool_search_dirs()``.
        threads: Passed to ``-num_threads``. Must be >= 1; blast rejects 0 and a
            silent 0 would mean a search that did not use what it was given.
        runner: Injected command runner, ``runner(command=...) ->
            CommandResult``. ``None`` runs the tools for real.

    Returns:
        TblastnRunInfo: the argv, both tools' versions, the database and
            parameters, the row count and the sha256 of the table as written.

    Raises:
        ToolNotAvailableError: ``tblastn`` or ``makeblastdb`` is missing. The
            message names which one.
        PipelineError: An input is not a file, ``threads`` is below 1, the table
            the tool produced does not parse, or an injected runner returned
            the wrong type.
        ToolExecutionError: ``makeblastdb`` or ``tblastn`` exited non-zero.
    """
    query = Path(query_faa)
    assembly = Path(assembly_fasta)
    table_path = Path(out_tsv)

    for path, what in (
        (query, "reference oprD protein (query)"),
        (assembly, "isolate assembly"),
    ):
        if not path.is_file():
            raise PipelineError(
                f"tblastn cannot run: the {what} was not found at {path}. The "
                f"structural oprD screen measures the isolate's reading frame "
                f"against the reference's, so both files are the search itself - "
                f"there is nothing to fall back on.",
                path=str(path),
            )

    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise PipelineError(
            f"tblastn was given threads={threads!r}; it needs an integer >= 1. "
            f"blast rejects 0, and passing it on would either fail the tool or "
            f"run a search that did not use the thread count the run allocated.",
            threads=threads,
        )

    tblastn = resolve_tool(TBLASTN, tool_search_dirs)
    makeblastdb = resolve_tool(MAKEBLASTDB, tool_search_dirs)

    prefix = database_prefix(table_path)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    # Always rebuilt, never reused. `virulence.ensure_blast_db` skips an existing
    # index because there the pinned sequences file is the thing being protected;
    # here the index is derived from a per-isolate assembly, and a stale index
    # beside a changed assembly would search the wrong sequence while reporting
    # success. makeblastdb overwrites in place and costs ~0.03 s per isolate.
    index_command = [
        str(makeblastdb),
        "-in", str(assembly),
        "-dbtype", "nucl",
        "-out", str(prefix),
    ]
    indexed = _run(index_command, runner=runner, what=MAKEBLASTDB)
    if indexed.returncode != 0:
        raise ToolExecutionError(
            f"makeblastdb failed for {assembly.name} (exit "
            f"{indexed.returncode}), so tblastn had no nucleotide database to "
            f"search: {indexed.stderr.strip()[:400]}",
            command=" ".join(index_command),
        )
    if not Path(f"{prefix}.nin").is_file():
        raise ToolExecutionError(
            f"makeblastdb exited 0 but wrote no nucleotide index at "
            f"{prefix}.nin, so tblastn would have searched nothing.",
            command=" ".join(index_command), path=f"{prefix}.nin",
        )

    command = oprd.tblastn_command(
        program=str(tblastn),
        query=query,
        database=Path(prefix),
        threads=threads,
    )
    result = _run(command, runner=runner, what=TBLASTN)
    if result.returncode != 0:
        raise ToolExecutionError(
            f"tblastn failed for {assembly.name} (exit {result.returncode}): "
            f"{result.stderr.strip()[:400]}",
            command=" ".join(command),
        )

    table = result.stdout or ""
    # Validated before it is put in place. `structural_call_from_paths` re-parses
    # this file with the same parser, so a row it would refuse must not be
    # written: a malformed table discovered there surfaces as a bare
    # PipelineError from a decision function, naming no tool and no stage.
    #
    # `query_codon_length` is the query's own residue count, and this function
    # is not given the reference CDS - only the query file - so it counts what
    # was handed rather than deriving it. `parse_tblastn_table` does not read
    # the argument as of this commit, and if a future parser uses it, the
    # terminator off-by-one that made the production query unsatisfiable becomes
    # visible in that commit rather than silently wrong here. The query itself is
    # written by `oprd_locus.write_tblastn_query`, which counts the reference's
    # codons through `oprd_locus.reference_codon_count` - the same definition
    # `measure_locus_structure` uses - so the two cannot drift apart.
    residues = sum(
        len(line.strip())
        for line in query.read_text(encoding="utf-8", errors="replace")
        .splitlines()
        if line.strip() and not line.startswith(">")
    )
    hits = oprd.parse_tblastn_table(table, query_codon_length=residues)

    table_path.parent.mkdir(parents=True, exist_ok=True)
    staging = table_path.parent / f".{table_path.name}.partial"
    staging.write_text(table, encoding="utf-8")
    os.replace(staging, table_path)

    LOGGER.info(
        "tblastn %s: %d HSP rows -> %s",
        assembly.stem, len(hits), table_path,
    )
    return TblastnRunInfo(
        out_tsv=table_path,
        query=query,
        assembly=assembly,
        database=str(prefix),
        command=list(command),
        parameters=(
            ("outfmt", oprd.TBLASTN_OUTFMT),
            ("evalue", str(oprd.TBLASTN_EVALUE)),
            ("max_target_seqs", str(oprd.TBLASTN_MAX_TARGET_SEQS)),
            ("db_gencode", str(oprd.TBLASTN_DB_GENCODE)),
            ("num_threads", str(threads)),
        ),
        tool_version=tool_version(tblastn),
        makeblastdb_command=list(index_command),
        makeblastdb_version=tool_version(makeblastdb, MAKEBLASTDB),
        n_hits=len(hits),
        table_sha256=hashlib.sha256(table.encode("utf-8")).hexdigest(),
        executable=tblastn,
    )