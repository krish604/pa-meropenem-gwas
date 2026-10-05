"""Stage 10 - a patristic distance matrix from the stage-9 phylogeny.

The stage answers one question: *how far apart are two isolates, according to
the tree the pipeline already built?* The answer is a **patristic distance** -
the sum of the branch lengths along the path between two tips.

**Not a distance matrix read back from IQ-TREE.** IQ-TREE builds the tree from
a distance matrix, so re-reading one and calling the result "phylogenetic
similarity" would be circular: the two would agree by construction while
measuring nothing the tree does not already say. Reading the tree also means
one definition of distance flows through stages 9 and 10 rather than two that
happen to agree.

**Output is a square matrix**, `sample_id` leading: one row per sample, one
distance column per sample, zero diagonal, symmetric. The long pair list it
replaced duplicated every distance - 7,405 rows at 835 isolates against 835 -
and was read by nothing. The `method` recorded alongside is `patristic` so a
reader cannot assume a different definition produced these numbers.

**A node's label is never its identity.** This is the invariant the whole module
is built around, and it is structural rather than a convention: the parser
assigns every node a unique pre-order ``node_id``, ancestry is a walk over
``parent`` pointers between node objects, and depths are looked up by that id.
No dictionary on the path from Newick text to a number is keyed by a label.

It has to be that way. IQ-TREE writes its branch-support values as internal node
labels, and support values saturate: on a real 10-tip *Pseudomonas* tree
``100/100`` occurs **three times**, at three different depths. An earlier
version keyed its parent map, its depth map and its ancestor chains on the
label, so those three nodes shared one key, last writer won, the parent map
acquired a cycle, and a ``seen`` guard quietly truncated the ancestor chains. On
that tree it raised; on a variant where the collision sat away from the root it
raised nothing at all and returned 15 of 45 pairs wrong, some of them negative.
The committed fixtures cannot show this - ``test_data/phylogeny/tree.nwk`` names
its internal nodes ``n1``..``n18`` - so the tests in
``test_similarity_duplicate_labels.py`` run against the real tree instead.

A tip label is used for exactly one thing: naming a row of the returned matrix,
which the data contract requires. That boundary is checked -
:func:`_require_unique_tip_labels` refuses two tips sharing a label, because two
tips with one name cannot be two rows. Tip labels come from the manifest, which
is required to hold no duplicates, so this cannot fire on a well-formed cohort.

**The matrix's units are recorded, in a sidecar.** A patristic distance is a sum
of branch lengths, and IQ-TREE's branch lengths are expected substitutions per
site. A reader who does not know that will read `0.0093` as nine SNPs, which is
wrong by the length of the alignment - and the error is invisible, because both
quantities are small decimals and the matrix is symmetric either way. So
:func:`write_units_sidecar` writes the unit, the definition and the explicit
"this is not a SNP count" beside the matrix, following the convention
`stages/regulators.py` set with `regulator_screen.json`: a table with no
sidecar is the ambiguous case.

Not a `#` comment above the header, which `write_tsv` supports and which
`read_tsv` skips: the contract fixes this file's shape as one header row plus
one row per sample, and a reader parsing the matrix by index must not have to
know that a comment may be there. The sidecar is a separate file, so the matrix
stays exactly the shape the contract describes.

The tip *parser* is reused from :mod:`papipeline.stages.phylogeny` so the two
stages cannot disagree about which strings are tips - see `extract_tips`. The
tree *structure* parser is this module's own, because that one needs branch
lengths and parent links, which `extract_tips` discards.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode
from .phylogeny import extract_tips

LOGGER = get_logger("stages.similarity")

#: How the distances were obtained. Recorded so a reader is never left guessing
#: whether these came from the tree, from IQ-TREE, or from something else.
METHOD = "patristic"

#: Leading column of the matrix.
ID_COLUMN = "sample_id"

#: `distances` key in each returned row.
DISTANCES = "distances"

#: Sidecar name, derived from the matrix's own name so the pair cannot drift.
UNITS_SIDECAR_SUFFIX = ".units.json"

#: The unit, stated so it can be quoted rather than paraphrased.
DISTANCE_UNIT = "expected substitutions per site"

#: What a patristic distance is NOT, in the words a reader is most likely to
#: reach for. This is the whole reason the sidecar exists: the numbers look like
#: SNP counts and are not, and nothing in a bare square matrix says so.
NOT_A_SNP_COUNT = (
    "These are NOT SNP counts. A SNP count is an integer number of differing "
    "sites; a patristic distance is a sum of branch lengths, and IQ-TREE's "
    "branch lengths are expected substitutions per site. Dividing by the "
    "alignment's site count does not recover a SNP count."
)


def units_sidecar_path(matrix_path: Path) -> Path:
    """`similarity.tsv` -> `similarity.units.json`, beside it."""
    matrix_path = Path(matrix_path)
    return matrix_path.with_name(matrix_path.name + UNITS_SIDECAR_SUFFIX)


def write_units_sidecar(matrix_path: Path, **facts: Any) -> Path:
    """Write the units record beside a matrix. Best effort, and never fatal.

    The matrix is the result; this is what makes it interpretable. Losing it must
    not lose the run, but losing it must be logged, because a matrix with no
    sidecar is exactly the ambiguous case this exists to remove.
    """
    payload: Dict[str, Any] = {
        "matrix": Path(matrix_path).name,
        "quantity": "patristic distance between two tips",
        "unit": DISTANCE_UNIT,
        "definition": (
            "the sum of the branch lengths along the path between two tips, "
            "i.e. depth(a) + depth(b) - 2*depth(LCA). The branch above the LCA "
            "is not part of the path and is not counted."
        ),
        "is_a_snp_count": False,
        "note": NOT_A_SNP_COUNT,
    }
    payload.update(facts)
    path = units_sidecar_path(matrix_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:  # pragma: no cover - filesystem dependent
        LOGGER.warning(
            "Stage 10: could not write the units sidecar at %s (%s). The matrix "
            "is intact, but a matrix with no sidecar is a matrix whose units a "
            "reader has to guess.",
            path,
            exc,
        )
    return path


def read_units_sidecar(matrix_path: Path) -> Optional[Dict[str, Any]]:
    """The sidecar's contents, or ``None`` when there is no sidecar."""
    path = units_sidecar_path(matrix_path)
    if not path.is_file():
        return None
    try:
        return dict(json.loads(path.read_text(encoding="utf-8")))
    except (OSError, ValueError) as exc:
        LOGGER.warning("Stage 10: unreadable units sidecar at %s: %s", path, exc)
        return None


class _Node:
    """One node of a parsed tree.

    **Identity is :attr:`node_id` and never :attr:`label`.** `node_id` is the
    node's position in the parse - a structural fact about the string, assigned
    in pre-order, unique by construction. `label` is carried as data and is used
    for nothing except (a) naming tip rows in the returned matrix and (b)
    reporting a duplicate *tip* label, which is a real error.

    That separation is the whole point of this class. The previous implementation
    used ``name`` as a dictionary key for the parent map, the depth map and the
    ancestor chains, which silently assumed node labels are unique. They are not:
    IQ-TREE writes its support values as internal labels, and on a real 10-tip
    tree ``100/100`` appears **three times**, at three different depths. One key,
    three writers, and every value downstream of it wrong. On that tree it
    raised; on a tree where the collision sits away from the root path it raised
    nothing at all and returned 15 of 45 pairs wrong, some of them negative.

    Holding identity in a field rather than in a dict key makes the collision
    unrepresentable: two nodes that share a label are still two distinct objects
    with two distinct ids, and nothing downstream compares them by name.
    """

    __slots__ = ("node_id", "label", "length", "children", "parent")

    def __init__(
        self,
        node_id: int,
        label: str,
        length: float,
        children: List["_Node"],
    ) -> None:
        self.node_id = node_id
        self.label = label
        self.length = length
        self.children = children
        self.parent: Optional["_Node"] = None

    @property
    def is_tip(self) -> bool:
        """A tip is a **leaf**, decided structurally."""
        return not self.children

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"_Node(id={self.node_id}, label={self.label!r})"


def _parse_newick(text: str) -> _Node:
    """Parse a Newick string into a tree of :class:`_Node`, linked to parents.

    Recursive descent over the actual grammar rather than a regular expression,
    so nested clades and quoted labels are handled correctly. The phylogeny stage
    has a parser for *tips*; this one keeps the branch lengths, which that one
    discards - a distance needs the path, not the label set.

    Two structural facts are established here and nowhere else:

    * every node gets a unique pre-order ``node_id``, which is its identity; and
    * every non-root node gets its ``parent`` pointer, so ancestry is a walk over
      objects rather than a lookup in a name-keyed table.

    The previous parser also *invented* a name for every unnamed internal node
    (``n1``, ``n2``, ...) because the parent chain was keyed on names and an
    unnamed node would otherwise have no entry. That workaround is gone: an
    unnamed internal node is now a perfectly ordinary node with an empty label
    and a real identity. IQ-TREE writes an unnamed root on every tree it
    produces, so this is the common case, not an edge case.
    """
    source = text.strip().rstrip(";").strip()
    position = 0
    counter = 0

    def parse_node() -> _Node:
        nonlocal position, counter
        counter += 1
        node_id = counter
        if position < len(source) and source[position] == "(":
            position += 1
            children: List[_Node] = []
            while True:
                children.append(parse_node())
                if position < len(source) and source[position] == ",":
                    position += 1
                    continue
                if position < len(source) and source[position] == ")":
                    position += 1
                    break
                raise DataContractError(
                    "Unbalanced parentheses in Newick string",
                    position=position,
                )
            label = _read_label()
            node = _Node(node_id, label, _read_length(), children)
        else:
            node = _Node(node_id, _read_label(), _read_length(), [])

        for child in node.children:
            child.parent = node
        return node

    def _read_label() -> str:
        nonlocal position
        start = position
        while position < len(source) and source[position] not in "(),:;[]":
            position += 1
        return source[start:position].strip()

    def _read_length() -> float:
        nonlocal position
        if position >= len(source) or source[position] != ":":
            return 0.0
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
                token=token,
            ) from exc

    return parse_node()


def _tips(root: _Node) -> List[_Node]:
    """Every leaf, in pre-order.

    Decided by structure - ``not node.children`` - and never by label. An
    earlier rule keyed on the label not matching ``^n\\d+$``, which held only for
    the committed fixture, whose internal nodes are all named ``n<N>``.
    """
    found: List[_Node] = []
    stack: List[_Node] = [root]
    while stack:
        node = stack.pop()
        if node.is_tip:
            found.append(node)
            continue
        for child in reversed(node.children):
            stack.append(child)
    found.sort(key=lambda n: n.node_id)
    return found


def _require_unique_tip_labels(tips: Sequence[_Node]) -> None:
    """Refuse two tips sharing a label.

    This is the *only* place a label is allowed to stand for a node, and it is
    the public matrix's boundary: rows are keyed by sample id. A duplicated tip
    label would collapse two rows into one, so it has to be an error rather than
    a silent drop - the same reasoning as the sample-id mismatch checks
    elsewhere in the pipeline.

    It cannot be papered over by a synthetic name, because inventing one would
    put a label in the matrix that is not a sample.
    """
    seen: Dict[str, _Node] = {}
    duplicates: Dict[str, List[int]] = {}
    for node in tips:
        if node.label in seen:
            duplicates.setdefault(node.label, [seen[node.label].node_id]).append(
                node.node_id
            )
        else:
            seen[node.label] = node
    if duplicates:
        label, ids = sorted(duplicates.items())[0]
        raise DataContractError(
            f"{len(duplicates)} tip label(s) occur more than once in the Newick "
            f"tree, starting with {label!r} at nodes {ids}. Two tips cannot share "
            "a name and still be two rows of a distance matrix, so this stops "
            "rather than collapsing them into one.",
            label=label,
            node_ids=",".join(str(i) for i in ids),
            n_duplicated_labels=len(duplicates),
        )


def _path_to_root(node: _Node) -> List[int]:
    """``node_id`` of every node from ``node`` up to the root, inclusive.

    **A cycle raises.** The previous version of this walk carried a `seen` guard
    that silently *stopped* when it revisited a node, which is what turned a
    label collision into a truncated ancestor chain and a wrong answer with no
    error. A parent chain built by :func:`_parse_newick` cannot cycle, because
    every ``parent`` pointer is assigned once, by a single pre-order walk from the
    root. So a cycle here means the structure was corrupted after parsing, and
    that is a fact to report rather than a length to clamp.
    """
    chain: List[int] = []
    visited: set = set()
    current: Optional[_Node] = node
    while current is not None:
        if current.node_id in visited:
            raise DataContractError(
                f"The tree's parent links contain a cycle at node "
                f"{current.node_id!r} (label {current.label!r}). Ancestry cannot "
                "be computed on a structure that loops, and truncating the walk "
                "would report a distance measured on part of the tree as though "
                "it were the whole tree.",
                node_id=current.node_id,
                label=current.label,
            )
        visited.add(current.node_id)
        chain.append(current.node_id)
        current = current.parent
    return chain


def _depths_by_node_id(root: _Node) -> Dict[int, float]:
    """Summed branch length from the root to every node, keyed by ``node_id``.

    Keyed by the structural id, so two nodes that share a label get two distinct
    entries. The root's own edge length is included as the base offset; it cancels
    out of every ``depth(tip) - depth(lca)`` and is kept so the values match what
    the module has always produced.
    """
    depths: Dict[int, float] = {}
    stack: List[Tuple[_Node, float]] = [(root, float(root.length))]
    while stack:
        node, depth = stack.pop()
        depths[node.node_id] = depth
        for child in node.children:
            stack.append((child, depth + float(child.length)))
    return depths


def patristic_distances(tree_text: str) -> Dict[str, Dict[str, float]]:
    """The full square matrix of patristic distances.

    Args:
        tree_text: A Newick string.

    Returns:
        ``{tip: {tip: distance}}``, including ``tip`` against itself at 0.0.

    Notes:
        Distances are ``depth(a) + depth(b) - 2 * depth(lca)``, which requires
        knowing the lowest common ancestor's depth. The alternative -
        walking up from each tip and adding until the paths meet - is what
        produces a **non-zero diagonal**, because a tip "meets itself" at the
        root rather than at itself. That error keeps the matrix symmetric and
        every off-diagonal value plausible, so it survives any check that reads
        a single entry.

        The lowest common ancestor is found by intersecting two *node-id* paths,
        so two internal nodes sharing a label cannot be confused for one another
        and a chain cannot be cut short.
    """
    root = _parse_newick(tree_text)
    tips = _tips(root)
    if not tips:
        raise DataContractError("Newick tree contains no tips")
    _require_unique_tip_labels(tips)

    depths = _depths_by_node_id(root)
    # Cached once per node rather than per pair: the LCA walk is the hot path,
    # and a path is a property of a node, not of a pair.
    paths = {node.node_id: _path_to_root(node) for node in tips}

    def distance(a: _Node, b: _Node) -> float:
        if a.node_id == b.node_id:
            return 0.0
        chain_a, chain_b = paths[a.node_id], paths[b.node_id]
        shared = set(chain_b)
        lca_id = next((i for i in chain_a if i in shared), None)
        if lca_id is None:
            raise DataContractError(
                f"Tips {a.label!r} and {b.label!r} are not in the same tree; the "
                "Newick string does not describe a single connected tree",
                tip_a=a.label,
                tip_b=b.label,
            )
        # Every node's depth is the summed length of the path from the root, so
        # subtracting gives the path from the tip up to that ancestor.
        return (depths[a.node_id] - depths[lca_id]) + (
            depths[b.node_id] - depths[lca_id]
        )

    return {
        a.label: {b.label: distance(a, b) for b in tips}
        for a in tips
    }


def matrix_rows(
    matrix: Mapping[str, Mapping[str, float]],
    order: Sequence[str],
) -> List[Dict[str, Any]]:
    """The matrix as rows in ``order``, each carrying its full distance vector.

    Row order follows the manifest order given, so the file is keyed on the
    cohort the run actually used rather than on dictionary iteration order -
    which would make the output differ between runs for no stated reason.
    """
    rows: List[Dict[str, Any]] = []
    for sample_id in order:
        distances = matrix.get(sample_id)
        if distances is None:
            raise DataContractError(
                f"{sample_id!r} is in the manifest but not in the tree",
                sample_id=sample_id,
                n_tips=len(matrix),
            )
        rows.append({ID_COLUMN: sample_id, DISTANCES: dict(distances)})
    return rows


def write_matrix(
    rows: Sequence[Mapping[str, Any]], path: Path
) -> Path:
    """Write the square matrix to ``path``.

    The header is ``sample_id`` followed by one column per sample, in the same
    order as the rows, so the file is a matrix a reader can index rather than a
    list of vectors. Written here rather than in the caller so the shape the
    contract describes and the shape on disk cannot drift.
    """
    order = [str(r[ID_COLUMN]) for r in rows]
    if not order:
        raise DataContractError("Cannot write an empty similarity matrix")

    path = Path(path)
    header = [ID_COLUMN] + order
    body: List[List[str]] = []
    for row in rows:
        distances = row[DISTANCES]
        body.append(
            [str(row[ID_COLUMN])]
            + [_format_distance(distances.get(other)) for other in order]
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join(r) for r in body]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _format_distance(value: Optional[float]) -> str:
    """Six significant figures, or an empty cell for a missing pair.

    A missing pair is left empty rather than filled with a zero: a zero would
    say "these two are identical", which is the one thing a distance matrix must
    never claim by accident.
    """
    if value is None:
        return ""
    return f"{value:.6g}"


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    tree_path: Optional[Path] = None,
    out_path: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    """Stage 10 entry point.

    Args:
        config: Loaded configuration. Unused today; taken so the eventual REAL
            caller's options have somewhere to come from.
        manifest: The cohort. Every sample must be a tip in the tree.
        mode: ``TEST`` reads the committed tree. ``REAL`` has no caller yet.
        tree_path: The stage-9 Newick tree.
        out_path: Where to write the matrix. Omit to skip writing.

    Returns:
        One row per manifest sample: ``{sample_id, distances}``.

    Raises:
        NotImplementedError: In ``REAL`` mode. There is no caller yet, and
            returning an empty matrix would claim every isolate is identical.
        DataContractError: A manifest sample is not a tip in the tree, or the
            tree has no tips.
    """
    if mode is RunMode.TEST:
        if tree_path is None:
            raise DataContractError(
                "TEST-mode similarity needs `tree_path`: it reads the committed "
                "stage-9 fixture tree. There is no default, because resolving one "
                "from configuration would let a TEST run measure a tree nobody "
                "chose."
            )
        tree = Path(tree_path)
    else:
        # REAL: the tree already exists. Stage 9 built it, so this reads it and
        # infers nothing - re-running IQ-TREE would produce a *different* tree,
        # and this matrix is a measurement of one specific tree, so a
        # re-inference would silently change the numbers being reported.
        #
        # The gate lives here rather than in a script because there is no
        # `scripts/similarity/run_similarity.py` - the Snakefile rule points at a
        # runner nobody wrote, which is the same reason this stage sits in
        # UNBUILT_STAGES. Every other REAL caller checks `allow_real_mode` in its
        # entry-point script; until that script exists, the stage is the only
        # place the check can be.
        if not bool(config.runtime.get("allow_real_mode", False)):
            raise NotImplementedError(
                "REAL-mode similarity is gated: runtime.allow_real_mode is false. "
                "Set it in the machine overlay, or open it for one session with "
                "PIPELINE_ALLOW_REAL_MODE=1. The committed overlays keep it shut "
                "so a REAL run cannot begin by accident. "
                f"(mode={getattr(mode, 'value', mode)})"
            )
        tree = config.phylogeny_dir(mode) / "tree.nwk"

    if not tree.is_file():
        raise DataContractError(
            f"Stage-9 tree not found at {tree}. The similarity matrix is read off "
            "that tree; without it there is nothing to measure. On a REAL run that "
            "means stage 9 produced no tree - and it cannot, because the core "
            "alignment it needs comes from panaroo, which has no installable "
            "osx-arm64 build. A missing tree is a missing input, not an absence of "
            "similarity.",
            path=str(tree),
            mode=getattr(mode, "value", mode),
        )

    newick = tree.read_text(encoding="utf-8")
    matrix = patristic_distances(newick)
    tips = set(matrix)
    declared = set(extract_tips(newick))

    missing = [s for s in manifest.sample_ids if s not in tips]
    if missing:
        raise DataContractError(
            f"{len(missing)} manifest sample(s) are not tips in the stage-9 "
            f"tree: {', '.join(missing[:10])}"
            + (" ..." if len(missing) > 10 else "")
            + ". The manifest decides cohort membership and the tree is evidence "
            "about it, so a sample in one and not the other is a mismatch to "
            "stop on - not a sample to drop. Dropping it would make every "
            "distance in the matrix describe a different cohort than the study "
            "claims.",
            missing=",".join(missing[:10]),
            n_missing=len(missing),
            n_tips=len(tips),
        )

    extra = sorted(tips - declared)
    if extra:
        # A tip the parser and the matrix disagree about would be a parser bug,
        # and it would put a label in the matrix that is not a sample.
        raise DataContractError(
            f"The tree has {len(extra)} tip(s) the matrix did not produce: "
            f"{', '.join(extra[:10])}",
            extra=",".join(extra[:10]),
        )

    order = [s for s in manifest.sample_ids]
    rows = matrix_rows(matrix, order)
    LOGGER.info(
        "Stage 10: %d x %d %s distance matrix from %s",
        len(rows), len(rows), METHOD, tree_path,
    )
    if out_path is not None:
        matrix_path = write_matrix(rows, Path(out_path))
        write_units_sidecar(
            matrix_path,
            method=METHOD,
            tree=str(tree),
            n_tips=len(tips),
            n_pairs=len(tips) * (len(tips) - 1) // 2,
        )
    return rows


__all__ = [
    "DISTANCE_UNIT",
    "DISTANCES",
    "ID_COLUMN",
    "METHOD",
    "NOT_A_SNP_COUNT",
    "matrix_rows",
    "patristic_distances",
    "read_units_sidecar",
    "run",
    "units_sidecar_path",
    "write_matrix",
    "write_units_sidecar",
]
