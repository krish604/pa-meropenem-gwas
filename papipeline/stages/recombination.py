"""Stage 8 — recombination masking over the core gene alignment.

Between the pangenome (7) and the phylogeny (9), per `spec.md` D1. It asks
Gubbins where recombination happened and produces two things the pipeline can
use: a **recombination-masked polymorphic-sites alignment**, which stages 9 and
12 read as their SNP input, and a **per-node table** naming the tree branches
that carry recombination.

**It does not produce the tree.** Gubbins builds one internally because it needs
one to assign SNPs to branches. That tree is kept as an extra output and is not
the contract output: stage 9 runs its own inference on the masked alignment, and
two files must not claim the same path.

**`node` is not `sample_id`.** The table is per *branch*: gubbins labels tips
with the alignment's sequence identifiers and internal nodes `Node_<n>`, and
both kinds appear as rows. This is why the stage is deliberately absent from
`DENSE_PER_SAMPLE_STAGES` — the property "one row per sample" does not hold, and
asserting it would be asserting something false about a table that reads like a
sample table.

**`mean_branch_length` is computed here, not read from Gubbins.**
`per_branch_statistics.csv` has thirteen columns — SNP counts, block counts,
base counts, `r/m`, `rho/theta`, genome length, clonal frame — and not one of
them is a branch length. So this stage derives it from the node-labelled tree:

    the arithmetic mean of the lengths of all edges in the clade subtended by
    the node, including the node's own edge

Defined for every row: a tip's clade is its single edge, an internal node's is
its whole subtree. See `mean_branch_lengths` for why not the narrower "mean of
the immediate children".

**A row with `n_snps = 0` is not a branch where recombination was looked for
and not found.** On the verified 10-isolate cohort two isolates are identical at
every column where both carry a called base — they differ at 77 columns, all of
them gap placement — so Gubbins has no SNP to assign to either terminal branch
and reports zero for both. Counting those as "branches without recombination"
reports an absent measurement as an absent finding.

**`run()` gates REAL on `allow_real_mode`; the REAL work happens one level up.**
`run.py` dispatches this stage through `run.derive_recombination_tables`, which
drives `papipeline.adapters.gubbins.run_gubbins` (the source-built binary, with
a self-test) and then calls :func:`block_rows` and :func:`mean_branch_lengths`
here. So the REAL path runs gubbins rather than parsing a directory, and this
entry point parses caller-supplied gubbins-shaped output, which is what the
committed TEST fixture under `test_data/intermediate/gubbins/` holds. `TEST` runs
the real parser and the real writer over that directory, with no external tool -
the same shape as `stages.similarity.run(tree_path=…)`. The gate on `run()` is
the second of two: `run_gubbins` checks `allow_real_mode` before it looks for
anything.

The design record, including the reasoning behind every choice here, is
``docs/design/recombination-contract.md``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError, ModeNotAllowedError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode

LOGGER = get_logger("stages.recombination")

#: The stage's declared contract table, `execution/contracts.py`
#: `STAGE_TABLES["recombination"]`. Re-declared here so the writer cannot emit
#: a header the contract does not check; the test that both agree is what keeps
#: the two from drifting.
BLOCK_COLUMNS: Tuple[str, ...] = (
    "node", "n_snps", "mean_branch_length", "recombination_detected",
)

#: Where the masked alignment is written, and the basename it must have.
#:
#: Not a choice available to this stage: `stages/phylogeny.py` reads
#: `phylogeny_dir / "core_snp_alignment.fasta"` and hands it straight to
#: IQ-TREE, and `stages/gwas.py` records it as the source of its SNP alignment.
#: An alignment written anywhere else is one no consumer reads.
SNP_ALIGNMENT_NAME = "core_snp_alignment.fasta"


def _clade_edge_lengths(newick: str) -> Dict[str, List[float]]:
    """Every edge in the clade subtended by each labelled node.

    Returned as a mapping from node label to that node's own edge length plus
    the edge lengths of everything beneath it.

    **Internal nodes must carry their gubbins label.** gubbins' tree labels an
    internal node *after* the closing parenthesis —
    `(a:1,b:1)Node_1:5` — which is not standard Newick, so a generic parser
    either drops the label or reads `Node_1` as the label of the preceding
    child. Both failure modes are silent and both put the wrong SNP count
    against the wrong branch length, which is the one number a reader will
    check. Parsed here for that specific shape rather than through a library,
    because the library is right about Newick and this is not Newick.

    Node labels are matched by full equality against the keys of
    `per_branch_statistics.csv`. No pattern matching: `Node_1` must not be
    allowed to match `Node_11`, and a prefix or regex match is how that
    happens.
    """
    source = newick.strip().rstrip(";").strip()
    position = 0
    # (label, own_edge_length_or_None, [child tuples])
    root: List[Any] = ["", None, []]

    def read_label() -> str:
        nonlocal position
        start = position
        while position < len(source) and source[position] not in "(),:;":
            position += 1
        return source[start:position].strip()

    def read_length() -> Optional[float]:
        nonlocal position
        if position >= len(source) or source[position] != ":":
            return None
        position += 1
        start = position
        while position < len(source) and source[position] not in "(),:;":
            position += 1
        token = source[start:position]
        try:
            return float(token)
        except ValueError as exc:
            raise DataContractError(
                "Newick branch length is not a number",
                token=token, tree=newick[:200],
            ) from exc

    def parse_node():
        nonlocal position
        if position < len(source) and source[position] == "(":
            position += 1
            children = []
            while True:
                children.append(parse_node())
                if position < len(source) and source[position] == ",":
                    position += 1
                    continue
                if position < len(source) and source[position] == ")":
                    position += 1
                    break
                raise DataContractError(
                    "Unbalanced parentheses in the gubbins tree",
                    position=position,
                )
            # gubbins puts the internal node's label AFTER the ')', and its own
            # edge length after that. Reading the label here is what pairs each
            # `Node_n` in the tree with its row in per_branch_statistics.csv.
            label = read_label()
            return [label, read_length(), children]
        label = read_label()
        length = read_length()
        return [label, length if length is not None else 0.0, []]

    root = parse_node()

    lengths: Dict[str, List[float]] = {}

    def collect(node) -> List[float]:
        label, own, children = node
        collected: List[float] = [own if own is not None else 0.0]
        for child in children:
            collected.extend(collect(child))
        if label:
            lengths.setdefault(label, []).extend(collected)
        return collected

    collect(root)
    return lengths


def mean_branch_lengths(newick: str) -> Dict[str, float]:
    """Mean edge length over each node's whole clade, keyed by node label.

    **Why the clade and not the immediate children.** "Mean branch length" for a
    node most often means the mean of the edges descending directly from it.
    That definition is undefined for a tip — a tip has no children — so it would
    leave `mean_branch_length` empty on all ten tip rows of the verified cohort.
    An empty cell in a branch-length column is exactly the kind of hole a later
    reader fills with a zero, and a zero there reads as "no evolution on this
    branch", which is a claim about biology.

    Taking the whole clade is defined for tips and internal nodes alike and is
    the quantity worth having: it is monotone in how much change the subtree
    carries, which is what you compare against `n_snps` to see whether long
    branches are the recombinant ones.

    The alternative reading — "the length of the edge subtending this node" —
    would be a rename of `mean_branch_length` rather than a mean of anything.
    """
    return {
        label: (sum(values) / len(values))
        for label, values in _clade_edge_lengths(newick).items()
        if values
    }


def block_rows(
    statistics: Sequence[Mapping[str, Any]],
    branch_lengths: Mapping[str, float],
) -> List[Dict[str, Any]]:
    """Join gubbins' per-branch statistics to per-node branch lengths.

    Every node in the statistics file must be present in the tree, and vice
    versa. A node in one and not the other means the two files describe
    different trees, and the pairing is then wrong rather than merely missing —
    so this refuses instead of writing a partial join.

    The root is included, and its branch length is legitimately 0.0: gubbins
    midpoint-roots, and a midpoint root on an already-rooted tree leaves the
    root edge with no length. That is a real zero, not a missing value.
    """
    if not branch_lengths:
        raise DataContractError(
            "the gubbins tree yielded no labelled nodes, so no branch length "
            "can be attached to any row of the per-branch statistics"
        )

    missing = [str(row["node"]) for row in statistics if str(row["node"]) not in branch_lengths]
    if missing:
        raise DataContractError(
            f"{len(missing)} node(s) in per_branch_statistics.csv are not in the "
            f"tree gubbins wrote: {', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else "")
            + ". The two files describe different trees, so every branch length "
            "after this point would be attached to the wrong node. This is not a "
            "missing value to fill in - it is a tree that does not match its own "
            "statistics.",
            missing=",".join(missing[:10]),
            n_missing=len(missing),
        )

    extra = sorted(set(branch_lengths) - {str(row["node"]) for row in statistics})
    if extra:
        raise DataContractError(
            f"the tree has {len(extra)} labelled node(s) the per-branch "
            f"statistics do not list: {', '.join(extra[:10])}"
            + (" ..." if len(extra) > 10 else "")
            + ". Refusing rather than dropping them, because a dropped node's "
            "SNP count would look like a node with no recombination.",
            extra=",".join(extra[:10]),
        )

    return [
        {
            "node": str(row["node"]),
            "n_snps": int(row["n_snps"]),
            # Not rounded here. `io.tsv.write_tsv` formats floats at `%.6g`, and
            # two places deciding precision is how a value ends up rounded twice
            # — or not at all, depending on which caller wrote the file.
            "mean_branch_length": float(branch_lengths[str(row["node"])]),
            "recombination_detected": int(row["recombination_detected"]),
        }
        for row in statistics
    ]


def _stub_rows() -> List[Dict[str, Any]]:
    """STUB mode's table: the declared header and nothing else.

    Not shared with TEST. A stub table that carried plausible rows would be the
    exact failure this repository's stub mode exists to avoid — a
    header-only-plus-one-fabricated-row table reads as a result.
    """
    return []


def _write_alignment(source: Path, destination: Path) -> Path:
    """Copy gubbins' masked SNP alignment to the name the consumers read.

    Copied rather than symlinked: the phylogeny and gwas stages read it in a
    different results tree, and a symlink across trees is a file that resolves
    on the machine that made the run and not on the one that reads it.
    """
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        raise DataContractError(
            f"gubbins reported success but its masked alignment is absent at "
            f"{source}. Without it there is no SNP alignment for the phylogeny "
            "or GWAS stage, and this stage would otherwise look like a success "
            "that produced nothing.",
            path=str(source),
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    intermediate_root: Optional[Path] = None,
    phylogeny_dir: Optional[Path] = None,
    gubbins_dir: Optional[Path] = None,
    out_path: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    """Stage 8 entry point.

    Args:
        config: Loaded configuration. Supplies the REAL gate and the tool
            search directories.
        manifest: The cohort. Unused on the TEST path, which reads whatever
            cohort the supplied gubbins output describes; taken so the eventual
            ``run_pipeline`` call has the signature the other stages share.
        mode: ``TEST`` parses caller-supplied Gubbins output. ``STUB``
            fabricates the header. ``REAL`` is gated on
            ``runtime.allow_real_mode`` and is not the dispatch path - see the
            module docstring.
        intermediate_root: The run's intermediate root. Unused on the TEST path.
        phylogeny_dir: Where the masked alignment is written.
        gubbins_dir: In ``TEST``, a directory holding Gubbins-shaped output. In
            ``REAL``, where gubbins will be run.
        out_path: Where to write ``recombination.tsv``. Omit to skip writing.

    Returns:
        One row per tree node: ``node``, ``n_snps``, ``mean_branch_length``,
        ``recombination_detected``.

    Raises:
        NotImplementedError: In ``REAL`` mode while ``runtime.allow_real_mode``
            is false. That is a *gate on the run*, not a refusal of the stage,
            and it is not what `run_pipeline` calls in REAL: the dispatch path
            is :func:`papipeline.run.derive_recombination_tables`, which drives
            :func:`papipeline.adapters.gubbins.run_gubbins` and
            :func:`stages.recombination.build_outputs`.
        DataContractError: In ``TEST``, ``gubbins_dir`` was not supplied, or the
            files it holds do not agree with each other.
    """
    if mode is RunMode.REAL:
        # A *gate on the run*, not a refusal of the stage: open it and the
        # stage proceeds. Same shape as `stages/similarity.py:477`, and the
        # taxonomy derives "refuses REAL" from the AST precisely so that these
        # two are not confused - a gate guarded by `allow_real_mode` must not
        # be listed in `run.REAL_REFUSING_STAGES`.
        #
        # Previously this raised unconditionally, on the grounds that
        # `run.py` had no dispatch branch for stage 8. There is one now:
        # `run.derive_recombination_tables`, which runs gubbins through the
        # adapter and calls `block_rows` + `mean_branch_lengths` directly rather
        # than this entry point. `run_gubbins` checks the same flag before it
        # looks for anything, so the check here is the second of two, not the
        # only one.
        if not bool(config.runtime.get("allow_real_mode", False)):
            raise NotImplementedError(
                "REAL-mode stage 8 (recombination) is gated: "
                "runtime.allow_real_mode is false. Set it in the machine overlay, "
                "or open it for one session with PIPELINE_ALLOW_REAL_MODE=1. The "
                "committed overlays keep it shut so a REAL run cannot begin by "
                "accident. The stage itself is built: `run.derive_recombination_"
                "tables` drives `adapters.gubbins.run_gubbins`, and a binary that "
                "fails its self-test is refused by name rather than run. "
                f"(mode={getattr(mode, 'value', mode)})"
            )

    if mode is RunMode.STUB:
        return _stub_rows()

    if gubbins_dir is None:
        # Same reasoning as `stages.similarity.run(tree_path=...)`: a default
        # resolved from configuration would let a TEST run parse a directory
        # nobody chose, and the resulting table would be attributed to this
        # stage's contract without anyone having selected the input.
        raise DataContractError(
            "TEST-mode recombination needs `gubbins_dir`: it reads gubbins-shaped "
            "output from a caller-supplied directory. There is no default, "
            "because resolving one from configuration would let a TEST run "
            "report a table it measured from an input nobody chose. "
            f"(mode={getattr(mode, 'value', mode)})"
        )

    from ..adapters import gubbins as gubbins_adapter

    gubbins_dir = Path(gubbins_dir)
    statistics = gubbins_adapter.parse_per_branch_statistics(
        gubbins_dir / f"{gubbins_adapter.PREFIX}{gubbins_adapter.PER_BRANCH_STATISTICS_SUFFIX}"
    )
    node_labelled = gubbins_dir / f"{gubbins_adapter.PREFIX}{gubbins_adapter.NODE_LABELLED_TREE_SUFFIX}"
    if not node_labelled.is_file():
        raise DataContractError(
            f"gubbins' node-labelled tree is absent at {node_labelled}. It is "
            "written on every successful run and is the only place a branch "
            "length for each node exists - per_branch_statistics.csv has no "
            "branch-length column - so without it the declared table cannot be "
            "filled. This is a missing input, not an absence of recombination.",
            path=str(node_labelled),
        )

    lengths = mean_branch_lengths(node_labelled.read_text(encoding="utf-8"))
    rows = block_rows(statistics, lengths)

    LOGGER.info(
        "Stage 8: %d node rows, %d with recombination detected",
        len(rows), sum(r["recombination_detected"] for r in rows),
    )

    if phylogeny_dir is not None:
        build_outputs(
            gubbins_dir=gubbins_dir,
            phylogeny_dir=Path(phylogeny_dir),
        )

    if out_path is not None:
        write_tsv(Path(out_path), rows, BLOCK_COLUMNS)
    return rows


def build_outputs(
    *,
    gubbins_dir: Path,
    phylogeny_dir: Path,
) -> Dict[str, Path]:
    """Copy gubbins' masked alignment and tree to where the pipeline reads them.

    The alignment takes the name `stages.phylogeny` and `stages.gwas` already
    look for. **The tree does not get renamed to `tree.nwk`.** That path belongs
    to stage 9's own IQ-TREE inference; two files claiming it would leave a
    downstream reader unable to tell which inference its distances were measured
    on. Gubbins' tree is kept under its own name as an extra output.

    Not re-rooted. Gubbins midpoint-roots its own tree; see
    ``docs/design/recombination-contract.md`` section 8 for why re-rooting a
    zero-length root edge is what trips dendropy's bare ``AssertionError``.
    """
    from ..adapters import gubbins as gubbins_adapter

    gubbins_dir, phylogeny_dir = Path(gubbins_dir), Path(phylogeny_dir)
    alignment = _write_alignment(
        gubbins_dir / f"{gubbins_adapter.PREFIX}{gubbins_adapter.SNP_ALIGNMENT_SUFFIX}",
        phylogeny_dir / SNP_ALIGNMENT_NAME,
    )
    written = {"snp_alignment": alignment}

    for key, suffix in (
        ("final_tree", gubbins_adapter.FINAL_TREE_SUFFIX),
        ("node_labelled_tree", gubbins_adapter.NODE_LABELLED_TREE_SUFFIX),
    ):
        source = gubbins_dir / f"{gubbins_adapter.PREFIX}{suffix}"
        if source.is_file():
            destination = phylogeny_dir / f"gubbins_{suffix.lstrip('.')}"
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
            written[key] = destination
    return written


__all__ = [
    "BLOCK_COLUMNS",
    "SNP_ALIGNMENT_NAME",
    "block_rows",
    "build_outputs",
    "mean_branch_lengths",
    "run",
]