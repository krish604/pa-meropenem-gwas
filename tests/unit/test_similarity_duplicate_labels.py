"""Duplicate internal labels must not be able to corrupt a patristic distance.

**The bug this file exists for.** `patristic_distances` used to build its parent
map, its depth map and its ancestor chains out of dictionaries keyed on the node
**label**. That silently assumes internal labels are unique, and on any real
tree they are not. IQ-TREE writes branch-support values as internal labels, and
support values saturate: on the tree in `test_data/similarity/real_tree.nwk`,
`100/100` appears **three times**, at depths 0.2639616385, 0.5696136823 and
0.9612000964.

Three nodes, one key, last writer wins. The parent map then held
`parent['100/100'] = '99.9/92'` and `parent['99.9/92'] = '100/100'` - a 2-cycle -
and a `seen` guard in the ancestor walk turned that cycle into a silently
truncated chain. On this tree it raised `DataContractError`; on the variant below
it raised nothing at all and returned **15 of 45 pairs wrong, some negative**
(patristic distance cannot be negative, and nothing complained).

**Why the committed fixtures never caught it.** `test_data/phylogeny/tree.nwk`
names its internal nodes `n1`..`n18`, which are unique by construction and so
cannot collide however wrong the algorithm is. The oracle in
`test_similarity_real_caller.py` is sound - it compares all pairs against
dendropy - but its own 8-taxon fixture happens to draw three *distinct* support
values (`0/65`, `93.7/80`, `85.2/81`). The checker was right and the input was
lucky. So these tests point it at a real tree.

**Provenance of the fixture.** `test_data/similarity/real_tree.nwk` is
byte-identical to `iqtree.treefile` from a REAL run of the phylogeny stage on the
10-isolate smoke cohort (`iqtree -m GTR+G -T 1 --seed 20240617 -B 1000
--alrt 1000`), copied in unmodified. It is committed rather than regenerated
because a test fixture that depends on a search finding a particular topology is
a test that fails when the search changes, not when the code regresses.

**What "correct" is measured against.** dendropy, following the pattern already
established at `test_similarity_real_caller.py:63`. Every pairwise distance is
compared, not a sample of them: a 45-pair matrix has room for one bad LCA to look
plausible, and sampling pairs is how a 33% error rate gets reported as a pass.

One trap worth naming: dendropy's `newick` schema rewrites an unquoted `_` in a
label to a space, so `TEST_PA_001` arrives as `TEST PA 001`. The real tree uses
dotted accessions and is unaffected, but any oracle over the synthetic fixtures
has to undo that or it will report a missing tip rather than a distance error.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.stages import similarity as stage

REPO = Path(__file__).resolve().parents[2]
REAL_TREE = REPO / "test_data" / "similarity" / "real_tree.nwk"

#: Read straight off the fixture; asserted below so it cannot drift silently.
THE_COLLIDING_LABEL = "100/100"
THE_COLLISION_COUNT = 3

dendropy = pytest.importorskip("dendropy")


@pytest.fixture(scope="module")
def real_tree_text() -> str:
    if not REAL_TREE.is_file():
        pytest.fail(
            f"{REAL_TREE} is missing. It is the real IQ-TREE tree this file "
            "exists to test against; without it every other test here is vacuous."
        )
    return REAL_TREE.read_text(encoding="utf-8").strip()


def _ground_truth(newick: str) -> dict:
    """All pairwise patristic distances from dendropy.

    Keyed on `Taxon` objects inside dendropy and only converted to labels at the
    boundary, so the oracle does not commit the mistake under test - deciding
    node identity from a label string.
    """
    tree = dendropy.Tree.get_from_string(newick, schema="newick")
    matrix = tree.phylogenetic_distance_matrix()
    taxa = {t.label: t for t in matrix.taxon_namespace}
    return {
        a: {b: float(matrix.distance(taxa[a], taxa[b])) for b in taxa}
        for a in taxa
    }


def _all_pairs(tips):
    for i, a in enumerate(tips):
        for b in tips[i + 1:]:
            yield a, b


class TestTheFixtureStillHasTheCollision:
    """A test for a bug needs a tree that has the bug's precondition."""

    def test_the_internal_label_really_is_duplicated_three_times(self, real_tree_text):
        # Labels are what follows a closing paren up to the length separator.
        # Counted from the text rather than from the parser so the assertion is
        # independent of the code under test.
        labels = re.findall(r"\)(\S+?):", real_tree_text)
        assert labels.count(THE_COLLIDING_LABEL) == THE_COLLISION_COUNT, (
            f"expected {THE_COLLIDING_LABEL} to appear {THE_COLLISION_COUNT} "
            f"times, found {labels.count(THE_COLLIDING_LABEL)}. If the fixture was "
            "replaced with a collision-free tree, these tests would pass for the "
            "wrong reason and the bug would be untested again."
        )

    def test_the_colliding_nodes_sit_at_three_different_depths(self, real_tree_text):
        """Same label, different depth - which is what made it destructive."""
        root = stage._parse_newick(real_tree_text)

        colliding = []
        stack = [root]
        while stack:
            node = stack.pop()
            if node.label == THE_COLLIDING_LABEL:
                colliding.append(node)
            stack.extend(node.children)

        assert len(colliding) == THE_COLLISION_COUNT
        assert len({n.node_id for n in colliding}) == THE_COLLISION_COUNT, (
            "the colliding nodes share a node_id, so they are not distinct nodes "
            "and this fixture no longer reproduces the collision"
        )

        depths = stage._depths_by_node_id(root)
        distinct = {depths[n.node_id] for n in colliding}
        assert len(distinct) == THE_COLLISION_COUNT, (
            f"the {THE_COLLISION_COUNT} nodes labelled {THE_COLLIDING_LABEL} share "
            f"depth(s) {sorted(distinct)}; a label collision only destroys the "
            "distance when the colliding nodes are at different depths"
        )


class TestTheRealTreeGivesAWholeMatrix:
    """The old parser raised on this tree. A matrix is the deliverable."""

    def test_it_produces_a_full_square_matrix(self, real_tree_text):
        matrix = stage.patristic_distances(real_tree_text)
        assert len(matrix) == 10
        for tip, row in matrix.items():
            assert sorted(row) == sorted(matrix), f"row {tip} is not one entry per tip"

    def test_every_off_diagonal_distance_is_positive(
        self, real_tree_text
    ):
        """The silent failure produced NEGATIVE patristic distances.

        A distance cannot be negative. Asserting positivity does not prove
        correctness - a wrong positive passes - but it pins the specific
        nonsense the old parser emitted here, so a regression to it cannot be
        waved through by a matrix that "looks like distances".
        """
        matrix = stage.patristic_distances(real_tree_text)
        negatives = [
            (a, b, v)
            for a, b, v in (
                (a, b, matrix[a][b]) for a, b in _all_pairs(sorted(matrix))
            )
            if v <= 0.0
        ]
        assert not negatives, (
            f"{len(negatives)} pair(s) have a non-positive patristic distance, "
            f"first {negatives[0]}. A path length is a sum of branch lengths and "
            "cannot be <= 0 unless the tree was walked incorrectly."
        )


class TestDendropyOracle:
    """Every pair, not a sample of them."""

    def test_all_pairwise_distances_match_dendropy(self, real_tree_text):
        ours = stage.patristic_distances(real_tree_text)
        truth = _ground_truth(real_tree_text)

        assert sorted(ours) == sorted(truth), (
            f"tip sets differ: ours={sorted(ours)} dendropy={sorted(truth)}"
        )

        checked = 0
        wrong = []
        for a, b in _all_pairs(sorted(truth)):
            expected = truth[a][b]
            actual = ours[a][b]
            checked += 1
            if abs(actual - expected) > 1e-9:
                wrong.append(f"d({a},{b}): ours={actual!r} dendropy={expected!r}")
        assert not wrong, (
            f"{len(wrong)} of {checked} pairs disagree with dendropy:\n  "
            + "\n  ".join(wrong[:10])
        )
        assert checked == 45, f"expected 45 pairs for 10 tips, checked {checked}"

    def test_the_diagonal_is_exactly_zero(self, real_tree_text):
        matrix = stage.patristic_distances(real_tree_text)
        assert [t for t in matrix if matrix[t][t] != 0.0] == []

    def test_it_is_symmetric(self, real_tree_text):
        matrix = stage.patristic_distances(real_tree_text)
        for a, b in _all_pairs(sorted(matrix)):
            assert matrix[a][b] == pytest.approx(matrix[b][a])


class TestTheSilentCaseIsCaughtToo:
    """The loud case is the lucky one; the quiet one is the dangerous one.

    Renaming the TOPMOST internal node's `100/100` to `100/101` - itself an
    ordinary IQ-TREE support value - leaves a tree that is exactly as real as
    before and exactly as connected. The outermost clade is the unnamed root and
    its label is never read, so what actually changes is that the two remaining
    `100/100` nodes are now both well below the root. Under the old parser no
    chain was then cut, nothing raised, and 15 of 45 pairs came back wrong with
    some of them negative.

    Derived here rather than committed, so the two trees cannot drift apart.
    """

    @pytest.fixture(scope="class")
    def quiet_tree_text(self, real_tree_text) -> str:
        head, sep, tail = real_tree_text.rpartition(f"){THE_COLLIDING_LABEL}:")
        assert sep, "the topmost label is not where the fixture says it is"
        return f"{head})100/101:{tail}"

    def test_the_variant_still_collides(self, quiet_tree_text):
        labels = re.findall(r"\)(\S+?):", quiet_tree_text)
        assert labels.count(THE_COLLIDING_LABEL) == THE_COLLISION_COUNT - 1
        assert labels.count("100/101") == 1

    def test_the_variant_raises_nothing_but_is_still_exact(self, quiet_tree_text):
        ours = stage.patristic_distances(quiet_tree_text)
        truth = _ground_truth(quiet_tree_text)

        wrong = [
            f"d({a},{b}): ours={ours[a][b]!r} dendropy={truth[a][b]!r}"
            for a, b in _all_pairs(sorted(truth))
            if abs(ours[a][b] - truth[a][b]) > 1e-9
        ]
        assert not wrong, (
            f"{len(wrong)} pair(s) wrong on a tree the old parser handled without "
            "complaining:\n  " + "\n  ".join(wrong[:10])
        )

    def test_the_variant_differs_from_the_original_tree(self, quiet_tree_text):
        """Guard the fixture: the rename must be the ONLY difference.

        If a future edit disturbed anything else, the test above would be
        measuring a different tree from the one the bug was found on, and a
        regression would be reported against the wrong input. The two trees must
        also stay the same LENGTH, so `zip` below cannot hide a shifted
        character by truncating.
        """
        original = REAL_TREE.read_text(encoding="utf-8").strip()
        assert len(original) == len(quiet_tree_text), "the rename changed the length"
        differing = [
            (i, o, n) for i, (o, n) in enumerate(zip(original, quiet_tree_text))
            if o != n
        ]
        assert len(differing) == 1, (
            f"expected exactly one character to differ, got {len(differing)}: "
            f"{differing[:5]}"
        )
        index, was, now = differing[0]
        assert (was, now) == ("0", "1"), f"unexpected edit {was!r} -> {now!r}"
        # The edited character is the last of the TOPMOST INTERNAL NODE's support
        # label - the child of the unnamed root, i.e. the final `)LABEL:LENGTH;`
        # in the string. Every branch length is untouched, so the variant differs
        # from the original in which label collides and in nothing else.
        topmost = f"){THE_COLLIDING_LABEL}:0.2639616385);"
        renamed = ")100/101:0.2639616385);"
        assert original.endswith(topmost), (
            "the fixture's topmost label/length tail changed; re-check the rename"
        )
        assert quiet_tree_text.endswith(renamed), (
            "the variant is not the original with only that label renamed"
        )
        # And the single edit is the final character of that label: the last
        # `)LABEL:` in the string is the topmost internal node's.
        label_start = original.rfind(f"){THE_COLLIDING_LABEL}:")
        assert label_start > 0
        assert index == label_start + len(f"){THE_COLLIDING_LABEL}") - 1


class TestLabelsAreNeverIdentity:
    """The invariant, asserted on the structure rather than in prose."""

    def test_colliding_nodes_get_distinct_ids_and_distinct_depths(self, real_tree_text):
        root = stage._parse_newick(real_tree_text)
        depths = stage._depths_by_node_id(root)

        colliding = [
            n
            for n in _walk(root)
            if n.label == THE_COLLIDING_LABEL
        ]
        ids = {n.node_id for n in colliding}
        assert len(ids) == len(colliding) == THE_COLLISION_COUNT
        assert len({depths[i] for i in ids}) == THE_COLLISION_COUNT

    def test_the_depth_map_has_one_entry_per_node_not_per_label(self, real_tree_text):
        """The old code's fatal shape: fewer keys than nodes.

        A depth map keyed by label is smaller than the tree. Counting is the
        cheapest possible statement of the invariant, and it fails loudly if
        anyone reintroduces label keying.
        """
        root = stage._parse_newick(real_tree_text)
        nodes = list(_walk(root))
        depths = stage._depths_by_node_id(root)
        assert len(depths) == len(nodes) > len({n.label for n in nodes}), (
            f"{len(nodes)} nodes but only {len(depths)} depths; identity has gone "
            "back to being the label"
        )

    def test_an_unnamed_internal_node_still_has_an_identity(self):
        """IQ-TREE writes an unnamed root on every tree it produces.

        The old parser invented a name (`n1`, `n2`, ...) for unnamed internal
        nodes purely so the name-keyed parent map would have an entry. With
        identity in `node_id` there is nothing to invent, and an unnamed node is
        an ordinary node.
        """
        matrix = stage.patristic_distances("((a:1,b:1):1,(c:1,d:1):1);")
        assert matrix["a"]["b"] == pytest.approx(2.0)
        assert matrix["a"]["d"] == pytest.approx(4.0)
        root = stage._parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
        assert root.label == "", "the root of that Newick carries no label"
        assert root.node_id and root.parent is None


class TestTruncationRaises:
    """The old `seen` guard truncated; a cycle must now be an error."""

    def test_a_cyclic_parent_chain_raises_instead_of_truncating(self):
        """The exact old failure mode, constructed deliberately.

        Two nodes pointing at each other is what a label collision used to build
        by accident. The old `while current not in seen` walk stopped and returned
        a short chain, and the caller measured a distance on part of the tree as
        though it were the whole tree. It must refuse instead.
        """
        a = stage._Node(1, THE_COLLIDING_LABEL, 1.0, [])
        b = stage._Node(2, THE_COLLIDING_LABEL, 1.0, [])
        a.parent = b
        b.parent = a

        with pytest.raises(DataContractError) as excinfo:
            stage._path_to_root(a)
        message = str(excinfo.value)
        assert "cycle" in message
        assert THE_COLLIDING_LABEL in message, (
            "the error must name the label involved, or it is not diagnosable"
        )

    def test_a_walk_to_the_root_is_not_truncated(self):
        """The positive counterpart: a real chain reaches the root intact."""
        root = stage._Node(1, "", 0.0, [])
        mid = stage._Node(2, THE_COLLIDING_LABEL, 1.0, [])
        leaf = stage._Node(3, "PDT000294805.1", 1.0, [])
        mid.parent, leaf.parent = root, mid
        assert stage._path_to_root(leaf) == [3, 2, 1]

    def test_duplicate_tip_labels_are_refused(self):
        """The one place a label may name a node is checked, not assumed.

        Rows of the matrix are keyed by sample id. Two tips sharing a label would
        collapse into one row, so this has to stop rather than silently drop one.
        """
        with pytest.raises(DataContractError) as excinfo:
            stage.patristic_distances("(a:1,a:1,b:1);")
        message = str(excinfo.value)
        assert "more than once" in message
        assert "'a'" in message


def _walk(root):
    """Every node, pre-order."""
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))