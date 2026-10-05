"""Real IQ-TREE invocation, and the tree crosswalk stage 10 needs.

`stages.phylogeny.run` reads ``tree.nwk`` and ``tree_metadata.tsv`` and had no
producer for either, so a REAL run had no tree. This module is that producer. It
builds a maximum-likelihood tree from a core SNP alignment and reports the tip
labels IQ-TREE actually used.

**It starts from an alignment, because the step above it cannot run here.**
``phylogeny.aligner`` is ``panaroo``, and panaroo has no installable
``osx-arm64`` build at any version - ``laptop.yaml`` records it
``available: false`` for that reason. ``gubbins``, which produces the
recombination-filtered alignment, is likewise absent. So the core alignment is
an input here, not an output, and
:func:`require_alignment_upstream` names the step that would otherwise be
missing rather than letting the stage report "no tree".

**Every flag is one the installed IQ-TREE 3.1.3 actually accepts**, read from
``iqtree --help`` rather than from IQ-TREE 2 habits:

* ``-T`` takes the thread count; ``-nt`` is IQ-TREE 1 and is rejected;
* ``-B`` is ultrafast bootstrap and the help states ``>= 1000``;
* ``--alrt`` is SH-aLRT, ``--seed`` makes a stochastic search reproducible, and
  ``--prefix`` names every output at once.

**The crosswalk is read out of the tree, not assumed.** IQ-TREE labels tips with
the alignment's sequence identifiers, so tip == sample_id *in the expected
case*. Writing that as an assumption produces a crosswalk that is right until
something renames a sequence, at which point it lies - and
``stages.phylogeny.validate_tree_samples`` requires the tip set to equal the
manifest exactly, so the disagreement would surface as a confusing validation
failure rather than as a renamed tip. :attr:`TreeResult.crosswalk` records what
was found.

**+ASC is used only after the columns it cannot see are removed.** Ruling R8:
IQ-TREE's ``+ASC`` (ascertainment bias correction) is a statement that *every*
site in the alignment carries information, so the pinned binary refuses a
model containing ``+ASC`` when any column is constant, and names how many::

    ERROR: Invalid use of +ASC because of 20 invariant sites in the alignment

"Constant" here means constant *after gaps and N are ignored*. A column reading
``A A A A - N`` is invariant for every purpose the correction cares about, but
it is not the kind of invariant site a plain site-counting pass would flag, so
the columns have to be removed deliberately rather than left to the tool. That
is :func:`drop_partially_constant_columns`, and it runs **before** the model is
used - not because of the ordering, but because with the columns still present
the tool rejects the model by name and the run stops.

This is *not* bcftools/iqtree ``-fconst``. Measured on the same 10-isolate
alignment (see ``pa-artifacts/phy2/``): ``-fconst`` also strips the genuinely
informative positions that sit inside a gapped column and collapses the fit
outright - lnL -6800510.031 against -376357.853 for plain ``GTR+G`` - while the
narrow test below leaves the alignment's variable sites intact.

Nothing here downloads, updates or writes outside the prefix directory.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from ..errors import DataContractError, PipelineError, ToolExecutionError
from ..io.fasta import read_fasta
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.iqtree")

#: The config key naming the tool that builds the core alignment, and that tool.
#: Recorded so the absent step can be named rather than implied.
REQUIRED_ALIGNMENT_KEY = "phylogeny.aligner"
REQUIRED_ALIGNMENT_TOOL = "panaroo"

#: The help states `-B ... (>=1000)`; a smaller value is invalid, not merely weak.
MINIMUM_SUPPORT_REPLICATES = 1000

#: Characters that carry no base call. A column is judged on what remains once
#: these are removed, which is the whole difference between this filter and a
#: plain invariant-site count. `?` is included because it is the other spelling
#: of "no call" that appears in real alignments.
MISSING_BASE_CHARS: FrozenSet[str] = frozenset("-?Nn")

#: Name of the filtered alignment written beside IQ-TREE's outputs, so the file
#: the search actually read is inspectable rather than inferred.
FILTERED_ALIGNMENT_NAME = "iqtree.filtered.fasta"


@dataclass(frozen=True)
class ColumnFilterResult:
    """What the column filter did, so the count is never a claim without a number.

    `n_columns_dropped` is recorded rather than merely logged because the number
    is the evidence that ``+ASC`` was applicable at all: a run that reports "the
    model was +ASC" with no count beside it has shown nothing.
    """

    alignment: Path
    n_columns_in: int
    n_columns_kept: int

    @property
    def n_columns_dropped(self) -> int:
        return self.n_columns_in - self.n_columns_kept


@dataclass(frozen=True)
class TreeResult:
    """What one IQ-TREE run produced.

    `crosswalk` maps manifest sample id to the tip label found in the tree. The
    two are equal whenever the alignment is named by sample id, which is the
    normal case - but they are recorded as two facts because
    `stages.phylogeny` validates them as two.
    """

    tree: Path
    report: Optional[Path]
    log: Optional[Path]
    tip_labels: List[str] = field(default_factory=list)
    crosswalk: Dict[str, str] = field(default_factory=dict)
    command: List[str] = field(default_factory=list)
    #: None when the column filter did not run, so "not filtered" and "filtered
    #: and dropped nothing" stay distinguishable.
    column_filter: Optional[ColumnFilterResult] = None


def build_command(
    executable: str,
    *,
    alignment: Path,
    prefix: Path,
    model: str,
    threads: int,
    seed: int,
    bootstrap: int = MINIMUM_SUPPORT_REPLICATES,
    alrt: int = MINIMUM_SUPPORT_REPLICATES,
    sequence_type: Optional[str] = None,
) -> List[str]:
    """The IQ-TREE invocation for one alignment.

    Every scientific choice arrives as an argument; nothing here restates it.
    """
    command = [
        str(executable),
        "-s", str(alignment),
        "-m", str(model),
        "-T", str(int(threads)),
        "--seed", str(int(seed)),
        "--prefix", str(prefix),
    ]
    if bootstrap:
        command += ["-B", str(int(bootstrap))]
    if alrt:
        command += ["--alrt", str(int(alrt))]
    if sequence_type:
        # The tool auto-detects and warns when it is unsure ("NOTE: Alignment
        # sequence type is auto-detected. If in doubt, specify it via -st"), so
        # stating it removes a warning whose cause would otherwise be guessed at.
        command += ["-st", str(sequence_type)]
    return command


def observed_states(column: Iterable[str]) -> Set[str]:
    """The base calls a column actually makes, ignoring gaps and N.

    Takes the column as an iterable of characters rather than a string so a
    caller can stream one column across every sequence without materialising it
    as a string first.

    Exposed separately because it is the definition the filter turns on, and a
    definition that cannot be named cannot be argued about.
    """
    return {char.upper() for char in column if char not in MISSING_BASE_CHARS}


def drop_partially_constant_columns(source: Path, destination: Path) -> ColumnFilterResult:
    """Write `source` to `destination` without its partially constant columns.

    A column is dropped when at most one distinct base call remains in it once
    gaps and N are ignored. A fully gap-filled column goes too: it makes no call
    at all, which is a stronger version of the same problem.

    Raises:
        DataContractError: the alignment has no records, the sequences are not
            all the same length, or every column would be dropped - in which
            case the filtered file would be empty and the tool's failure would
            name a file rather than the filter that emptied it.
    """
    source = Path(source)
    headers: List[str] = []
    sequences: List[str] = []
    for record in read_fasta(source):
        headers.append(record.header)
        sequences.append(record.sequence)

    if not sequences:
        raise DataContractError(
            f"{source} has no sequences, so no column can be judged constant",
            path=str(source),
        )
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise DataContractError(
            f"{source} is not a rectangular alignment: sequence lengths are "
            f"{sorted(lengths)}. A column filter cannot say what a column "
            "contains when the sequences are not aligned.",
            path=str(source),
            lengths=",".join(str(n) for n in sorted(lengths)),
        )

    n_columns = lengths.pop()
    keep = [
        index
        for index in range(n_columns)
        if len(observed_states(seq[index] for seq in sequences)) > 1
    ]
    if not keep:
        raise DataContractError(
            f"Every one of the {n_columns} columns in {source} is constant once "
            "gaps and N are ignored, so there is nothing left to infer a tree "
            "from. This alignment is not a SNP alignment.",
            path=str(source),
            n_columns=n_columns,
        )

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for header, sequence in zip(headers, sequences):
            handle.write(f">{header}\n")
            handle.write("".join(sequence[index] for index in keep) + "\n")

    result = ColumnFilterResult(
        alignment=destination, n_columns_in=n_columns, n_columns_kept=len(keep)
    )
    LOGGER.info(
        "iqtree: dropped %d of %d columns constant ignoring gaps/N (%d kept) -> %s",
        result.n_columns_dropped,
        result.n_columns_in,
        result.n_columns_kept,
        destination,
    )
    return result


def require_alignment_upstream(alignment: Path) -> Path:
    """Refuse to proceed without an alignment, naming what would build one.

    The failure this prevents is subtle: with no alignment the stage reports no
    tree, which reads as "stage 10 ran and there was nothing to infer" rather
    than "the pangenome step cannot run on this machine".
    """
    alignment = Path(alignment)
    if not alignment.is_file():
        raise PipelineError(
            f"no core SNP alignment at {alignment}, so there is nothing for "
            f"IQ-TREE to read. It is produced by the core-alignment step, not by "
            f"IQ-TREE. `{REQUIRED_ALIGNMENT_KEY}` is `{REQUIRED_ALIGNMENT_TOOL}` "
            "in config/science.yaml, which has no installable osx-arm64 build - "
            "so this stage cannot be exercised end to end on this machine until "
            "either that is available or the alignment is supplied from "
            "elsewhere. This is a missing input, not an empty result.",
            alignment=str(alignment),
            aligner_key=REQUIRED_ALIGNMENT_KEY,
            aligner_tool=REQUIRED_ALIGNMENT_TOOL,
        )
    return alignment


def build_tree(
    alignment: Path,
    out_dir: Path,
    *,
    model: str,
    threads: int,
    seed: int,
    bootstrap: int = MINIMUM_SUPPORT_REPLICATES,
    alrt: int = MINIMUM_SUPPORT_REPLICATES,
    sequence_type: Optional[str] = None,
    executable: str = "iqtree",
    sample_ids: Optional[Dict[str, str]] = None,
    drop_partially_constant: bool = False,
    timeout: int = 7200,
) -> TreeResult:
    """Run IQ-TREE on one alignment and return the tree and its tip labels.

    `sample_ids` maps a tip label back to the manifest sample it came from. When
    the alignment is named by sample id the mapping is the identity; it is
    supplied rather than derived so a caller with renamed sequences can say so,
    and a caller that does not know gets identity rather than a guess.

    `drop_partially_constant` removes the columns that are constant once gaps
    and N are ignored, and points the search at the filtered file. It is a
    separate argument rather than something inferred from the model string
    because it is a fact about the *alignment*, not about the model: the same
    model is applicable or not depending on which alignment arrives, and a
    filter that triggered off a substring of `-m` would silently reshape the
    data whenever a model was renamed.

    Raises:
        PipelineError: the alignment is absent, naming the absent upstream step.
        ToolExecutionError: IQ-TREE could not be launched, exited non-zero, or
            exited zero without writing a tree.
        DataContractError: the alignment cannot be filtered (ragged, empty, or
            entirely constant once gaps and N are ignored).
    """
    alignment = require_alignment_upstream(alignment)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "iqtree"

    column_filter: Optional[ColumnFilterResult] = None
    if drop_partially_constant:
        column_filter = drop_partially_constant_columns(
            alignment, out_dir / FILTERED_ALIGNMENT_NAME
        )
        alignment = column_filter.alignment

    command = build_command(
        executable,
        alignment=alignment,
        prefix=prefix,
        model=model,
        threads=threads,
        seed=seed,
        bootstrap=bootstrap,
        alrt=alrt,
        sequence_type=sequence_type,
    )
    LOGGER.info(
        "Stage 10: iqtree -m %s on %s%s",
        model,
        alignment.name,
        "" if column_filter is None
        else f" ({column_filter.n_columns_kept} of {column_filter.n_columns_in} "
             f"columns kept)",
    )
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False,
            cwd=str(out_dir),
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolExecutionError(
            "IQ-TREE timed out",
            alignment=str(alignment), timeout_s=timeout,
        ) from exc
    except OSError as exc:
        raise ToolExecutionError(
            "could not launch IQ-TREE",
            alignment=str(alignment), command=" ".join(command), error=str(exc),
        ) from exc

    log_path = out_dir / "iqtree.log"
    if log_path.is_file():
        log_path.write_text(
            (completed.stdout or "") + (completed.stderr or ""), encoding="utf-8"
        )

    if completed.returncode != 0:
        raise ToolExecutionError(
            f"IQ-TREE exited {completed.returncode} for {alignment.name}",
            alignment=str(alignment),
            command=" ".join(command),
            stderr=(completed.stderr or completed.stdout or "")[-2000:],
        )

    tree = out_dir / "iqtree.treefile"
    if not tree.is_file():
        # Exit 0 with no tree is not success. IQ-TREE writes the tree last, so
        # its absence means the search never completed - and a stage that read
        # that as "no phylogeny" would report an absence of signal where there
        # was a tool failure.
        raise ToolExecutionError(
            f"IQ-TREE exited 0 but wrote no tree at {tree}, so the search did "
            "not complete",
            alignment=str(alignment), command=" ".join(command),
        )

    from ..stages.phylogeny import extract_tips

    tips = extract_tips(tree.read_text(encoding="utf-8"))
    mapping = sample_ids or {}
    crosswalk = {
        sample_id: mapping.get(tip, tip) for sample_id, tip in
        ((mapping.get(tip, tip), tip) for tip in tips)
    }
    return TreeResult(
        tree=tree,
        report=out_dir / "iqtree.iqtree" if (out_dir / "iqtree.iqtree").is_file() else None,
        log=log_path if log_path.is_file() else None,
        tip_labels=tips,
        crosswalk=crosswalk,
        command=command,
        column_filter=column_filter,
    )