"""UI-D6 — node identity is the node's position in the tree, never its label.

Two traps, both measured, both data corruption rather than display:

1. **dendropy's `preserve_underscores` defaults to `False`.** With the default
   every `TEST_PA_001` comes back as `TEST PA 001`, which silently breaks the
   `sample_id` join on every page. A test oracle built on the default is
   therefore worthless, so the parser here is hand-written and preserves bytes
   exactly; dendropy is available to a *test* to compare against, with the flag
   set explicitly.

2. **Internal Newick labels are unusable as identities.** Measured on
   `test_data/similarity/real_tree.nwk`: `100/100` appears **three times**,
   alongside `99.9/92`, `100/92`, `83.4/88`, `99/100`, and one node with no
   label at all. So `node_id` is the **postorder index** from the root — one
   depth-first pass, children left to right, counter incremented after the
   subtree — and the label is carried for display only.

The `identity` string the endpoint returns says so in the response, so a
client cannot mistake "the node the server called 17" for "clade 17".

Distances come from branch lengths. `stages.similarity.patristic_distances` is
imported as the definition, and rule 11 of `docs/scientific_rules.md` holds: a
branch length is expected substitutions per site, not a SNP count.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

#: What `node_id` is. Returned verbatim in the `identity` field.
IDENTITY = "postorder index from the root, assigned server-side"

#: Support labels are `100/100`, `83.4/88`, `99.9/92` — never a bare number,
#: and never unique. A label that parses as a number is offered as `support`;
#: one that does not is display-only text.
_SUPPORT_SEPARATORS = ("/", "|")


class NewickError(Exception):
    """The Newick could not be parsed. The message names what was expected."""


@dataclass
class Node:
    """One node. `node_id` is assigned after the whole parse, not before."""

    label: Optional[str] = None
    length: Optional[float] = None
    children: List["Node"] = field(default_factory=list)
    node_id: Optional[int] = None
    parent_id: Optional[int] = None
    _line: int = 0

    @property
    def is_tip(self) -> bool:
        return not self.children

    @property
    def support(self) -> Optional[float]:
        """The support value, only when the label parses as one.

        `100/100` is a support pair in some dialects and gubbins' own
        bookkeeping in others; a slash-joined label is reported as text, not
        as a number, because treating `83.4/88` as `83.4` would be a different
        claim.
        """
        if not self.label or self.is_tip:
            return None
        for sep in _SUPPORT_SEPARATORS:
            if sep in self.label:
                return None
        try:
            return float(self.label)
        except (TypeError, ValueError):
            return None


class _Parser:
    """Recursive descent over the Newick grammar.

        tree    := subtree ';'
        subtree := leaf | '(' subtree (',' subtree)* ')' [label] [length]
        leaf    := label [length]
        length  := ':' number

    Labels may be single- or double-quoted, in which case whitespace is
    significant and `''` is an escaped quote. Nothing is substituted, trimmed
    beyond the Newick rule of stripping unquoted surrounding whitespace, or
    case-folded: a tip label is returned byte-identical to the file.
    """

    def __init__(self, text: str) -> None:
        self.text = text
        self.pos = 0

    def parse(self) -> Node:
        self._skip_ws()
        if self._peek() is None:
            raise NewickError("the Newick is empty")
        root = self._subtree()
        self._skip_ws()
        if self._peek() == ";":
            self.pos += 1
        self._skip_ws()
        if self._peek() is not None:
            raise NewickError(
                f"trailing content after the tree at character {self.pos}: "
                f"{self.text[self.pos:self.pos + 20]!r}"
            )
        if not isinstance(root, Node):
            raise NewickError("the Newick holds no subtree")
        return root

    # -- characters -----------------------------------------------------
    def _peek(self) -> Optional[str]:
        return self.text[self.pos] if self.pos < len(self.text) else None

    def _skip_ws(self) -> None:
        while self.pos < len(self.text) and self.text[self.pos] in " \t\r\n":
            self.pos += 1

    def _expect(self, char: str) -> None:
        self._skip_ws()
        if self._peek() != char:
            raise NewickError(
                f"expected {char!r} at character {self.pos}, found "
                f"{self._peek()!r}"
            )
        self.pos += 1

    # -- grammar --------------------------------------------------------
    def _subtree(self) -> Node:
        self._skip_ws()
        char = self._peek()
        if char is None:
            raise NewickError(
                f"unexpected end of the Newick at character {self.pos}"
            )
        if char == "(":
            return self._clade()
        node = Node(label=self._label())
        node.length = self._length()
        return node

    def _clade(self) -> Node:
        self._expect("(")
        children = [self._subtree()]
        while True:
            self._skip_ws()
            char = self._peek()
            if char == ",":
                self.pos += 1
                children.append(self._subtree())
            elif char == ")":
                self.pos += 1
                break
            else:
                raise NewickError(
                    f"expected ',' or ')' at character {self.pos}, found "
                    f"{char!r}"
                )
        self._skip_ws()
        label = self._label() if self._peek() not in (None, ":", ",", ")", ";") else None
        node = Node(label=label, children=children)
        node.length = self._length()
        return node

    def _label(self) -> Optional[str]:
        self._skip_ws()
        char = self._peek()
        if char is None:
            return None
        if char in ("'", '"'):
            return self._quoted(char)
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in "(),:;[] \t\r\n":
            self.pos += 1
        if self.pos == start:
            return None
        return self.text[start:self.pos]

    def _quoted(self, quote: str) -> str:
        self.pos += 1  # opening quote
        out: List[str] = []
        while self.pos < len(self.text):
            char = self.text[self.pos]
            if char == quote:
                if self.pos + 1 < len(self.text) and self.text[self.pos + 1] == quote:
                    out.append(quote)
                    self.pos += 2
                    continue
                self.pos += 1
                return "".join(out)
            out.append(char)
            self.pos += 1
        raise NewickError(
            f"unterminated {quote}-quoted label at character {self.pos}"
        )

    def _length(self) -> Optional[float]:
        if self._peek() != ":":
            return None
        self.pos += 1
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in "(),; \t\r\n":
            self.pos += 1
        raw = self.text[start:self.pos]
        if not raw:
            raise NewickError(f"branch length with no number at character {start}")
        try:
            return float(raw)
        except ValueError:
            raise NewickError(
                f"branch length {raw!r} at character {start} is not a number"
            ) from None


def parse(text: str) -> Node:
    """Parse a Newick string and assign postorder ids from the root."""
    root = _Parser(text).parse()
    counter = _assign_postorder(root, [0])
    return root


def _assign_postorder(node: Node, counter: List[int]) -> Optional[int]:
    """One depth-first pass, children left to right, counter after the subtree.

    Postorder means a node's id is larger than every id in its subtree, so a
    child always has a smaller id than its parent and the ordering is stable
    for the same bytes in the same order — which is what makes it reproducible
    across two parses of one file.
    """
    for child in node.children:
        _assign_postorder(child, counter)
    node_id = counter[0]
    counter[0] += 1
    node.node_id = node_id
    for child in node.children:
        child.parent_id = node_id
    return node_id


def assign_postorder(node: Node) -> Node:
    """Public entry, for a node built by other means."""
    _assign_postorder(node, [0])
    return node


def tip_labels(root: Node) -> List[str]:
    """Every tip label, in the tree's own left-to-right order."""
    labels: List[str] = []

    def walk(node: Node) -> None:
        if node.is_tip:
            labels.append(node.label or "")
            return
        for child in node.children:
            walk(child)

    walk(root)
    return labels


def to_json(root: Node, *, tip_metadata: Optional[Mapping[str, Mapping[str, Any]]] = None) -> Dict[str, Any]:
    """Structural JSON: `node_id`, `parent_id`, `children`, `length`, `is_tip`,
    `label`, `support`.

    `label` is display only. `tip_metadata` is keyed on the tip's own label,
    which is byte-identical to the file.
    """
    nodes: List[Dict[str, Any]] = []
    tips: List[str] = []
    lengths_seen = False

    def walk(node: Node) -> None:
        nonlocal lengths_seen
        if node.length is not None:
            lengths_seen = True
        payload: Dict[str, Any] = {
            "node_id": node.node_id,
            "parent_id": node.parent_id,
            "children": [c.node_id for c in node.children],
            "length": node.length,
            "is_tip": node.is_tip,
            "label": node.label,
            "support": node.support,
        }
        if node.is_tip and tip_metadata:
            extra = tip_metadata.get(node.label or "")
            if extra:
                payload["tip_metadata"] = dict(extra)
        nodes.append(payload)
        if node.is_tip:
            tips.append(node.label or "")
        for child in node.children:
            walk(child)

    walk(root)
    return {
        "identity": IDENTITY,
        "n_tips": len(tips),
        "n_internal": len(nodes) - len(tips),
        "has_branch_lengths": lengths_seen,
        "nodes": nodes,
    }


def tip_check(tree_tips: Sequence[str], cohort: Sequence[str]) -> Dict[str, Any]:
    """The same set difference `phylogeny.validate_tree_samples` reports.

    A tree whose tips cannot be joined is reported, not drawn silently. The
    underscore trap is exactly what this catches: with dendropy's default the
    tip set is `TEST PA 001` and every sample reads as missing.
    """
    tree_set = set(tree_tips)
    cohort_set = set(cohort)
    missing = sorted(cohort_set - tree_set)
    extra = sorted(tree_set - cohort_set)
    duplicates = len(tree_tips) - len(tree_set)
    reason_bits: List[str] = []
    if missing:
        reason_bits.append(f"{len(missing)} cohort sample(s) are not tips in this tree")
    if extra:
        reason_bits.append(f"{len(extra)} tip(s) in this tree are not in the cohort")
    if duplicates:
        reason_bits.append(f"{duplicates} tip label(s) appear more than once")
    return {
        "matched": not (missing or extra or duplicates),
        "n_missing_in_tree": len(missing),
        "n_extra_in_tree": len(extra),
        "missing": missing,
        "extra": extra,
        "reason": (
            "the tip set equals the cohort set"
            if not reason_bits
            else "; ".join(reason_bits)
            + " - this tree's tips cannot be joined to the isolates by "
            "sample_id, so a per-isolate join against it would be a guess"
        ),
    }


def patristic_distances(root: Node) -> Dict[str, Dict[str, float]]:
    """Every tip-to-tip patristic distance, summed along branch lengths.

    The MRCA's own branch is **not** in the path: a patristic distance is
    `depth(a) + depth(b) - 2*depth(LCA)`, and adding the LCA's stem overstates
    one sister pair by 128 % (rule 11). The pipeline's own
    `stages.similarity.patristic_distances` is the authority and is used by
    the similarity endpoint; this is the tree-level primitive that endpoint and
    the tests share.
    """
    depth: Dict[int, float] = {}
    labels: Dict[int, str] = {}

    def walk(node: Node, accumulated: float) -> None:
        total = accumulated + (node.length or 0.0)
        if node.node_id is not None:
            depth[node.node_id] = total
            labels[node.node_id] = node.label or ""
        for child in node.children:
            walk(child, total)

    walk(root, 0.0)

    # Ancestor path to the root for each node, used for the LCA by climbing.
    parent: Dict[int, Optional[int]] = {}

    def record(node: Node) -> None:
        if node.node_id is not None:
            parent[node.node_id] = node.parent_id
        for child in node.children:
            record(child)

    record(root)

    def ancestor_path(node_id: int) -> List[int]:
        chain = [node_id]
        cursor = parent.get(node_id)
        while cursor is not None:
            chain.append(cursor)
            cursor = parent.get(cursor)
        return chain

    tip_ids = [
        payload["node_id"] for payload in to_json(root)["nodes"] if payload["is_tip"]
    ]
    result: Dict[str, Dict[str, float]] = {}
    for left in tip_ids:
        left_chain = ancestor_path(left)
        left_set = set(left_chain)
        row: Dict[str, float] = {}
        for right in tip_ids:
            distance = depth[left] + depth[right]
            for candidate in left_chain:
                if candidate in set(ancestor_path(right)):
                    distance -= 2 * depth[candidate]
                    break
            row[labels[right]] = distance
        result[labels[left]] = row
    return result


def downsample_by_single_linkage(
    tips: Sequence[str],
    vectors: Mapping[str, Sequence[float]],
    k: int,
) -> List[str]:
    """Pick `k` tips so the picture keeps its shape.

    Greedy farthest-point sampling on the cohort's own distance vectors: start
    from the tip with the smallest total distance to the rest, then repeatedly
    add the tip furthest from everything already chosen. That keeps the spread
    of the cohort rather than its first `k` in file order, which is what a
    silent truncation does.

    The vectors are positional — index `i` is the distance to `tips[i]` — which
    is how `similarity.tsv` stores them: one row per sample holding a whole
    `distances` cell, keyed on the cohort rather than on columns.

    Near-linear in `n` for a fixed `k`. The naive greedy is O(k^2 * n) because
    every candidate re-measures its distance to every chosen tip on every round
    (300 x 900 x 300 = 81M lookups at the default). This keeps each
    not-yet-chosen tip's distance to its *nearest* chosen tip in an array, so
    adding a tip is one pass over the cohort and the next pick is an argmax
    over that array: O(k * n) total, ~270K lookups at the default.
    """
    ordered = list(tips)
    n = len(ordered)
    if n <= k or k <= 0:
        return ordered

    rows: List[Sequence[float]] = [vectors.get(name) or () for name in ordered]

    def cell(i: int, j: int) -> float:
        """``distance(ordered[i], ordered[j])``: row ``i``, column ``j``."""
        row = rows[i]
        return float(row[j]) if j < len(row) else 0.0

    # Seed: the tip with the smallest total distance to the rest of the cohort.
    first = min(range(n), key=lambda i: sum(rows[i][:n]))

    chosen: List[int] = [first]
    chosen_set = {first}

    # `nearest[i]` is ordered[i]'s distance to its nearest chosen tip.
    nearest = [cell(i, first) for i in range(n)]
    nearest[first] = 0.0

    while len(chosen) < k:
        nxt = -1
        best = -1.0
        for i in range(n):
            if i not in chosen_set and nearest[i] > best:
                best = nearest[i]
                nxt = i
        if nxt < 0:
            break
        chosen.append(nxt)
        chosen_set.add(nxt)
        for i in range(n):
            if i in chosen_set:
                continue
            d = cell(i, nxt)
            if d < nearest[i]:
                nearest[i] = d

    return [ordered[i] for i in chosen]


__all__ = [
    "IDENTITY",
    "NewickError",
    "Node",
    "assign_postorder",
    "downsample_by_single_linkage",
    "parse",
    "patristic_distances",
    "tip_check",
    "tip_labels",
    "to_json",
]