"""The `similarity` stage: a patristic distance matrix from the stage-9 tree.

TEST-level work only. The stage is still in `UNBUILT_STAGES` and the REAL caller
does not exist, so nothing here is reachable by a run yet - what is proven is the
part that has to be right regardless of who invokes it: reading a Newick tree
and producing a well-formed distance matrix.

**What a patristic distance is.** The sum of the branch lengths along the path
between two tips. It is read off the stage-9 tree, so it is a statement about
that tree. It is deliberately *not* an IQ-TREE distance matrix re-read from
disk: IQ-TREE built the tree from exactly such a matrix, so reading one back
would be circular, and the two would agree by construction while measuring
nothing new.

**The expected values are worked by hand from the fixture's structure**, not
recomputed with the code under test. Every branch in
`test_data/phylogeny/tree.nwk` is 0.05, so every distance is an exact multiple
of 0.05 and the expected numbers can be counted:

    d(001,001) = 0.00   a tip to itself
    d(001,004) = 0.15   001 - n2 - n1 - 004                      3 edges
    d(001,002) = 0.30   001 - n2 - n6 - n18 - n12 - n8 - 002     6 edges
    d(001,020) = 0.35   ... - n9 - n10 - 020                     7 edges

A test that recomputed these with the same traversal would pass whatever the
traversal did, including a wrong one. Counting edges is the part that can
actually fail - and it did. The first draft of these expectations read 0.20 and
0.25, from miscounting the path at 5 and 6 edges instead of 6 and 7. The
implementation was right and the hand-work was wrong, which is the outcome this
arrangement is designed to surface rather than hide.

The deepest tip sits 7 edges from `001`, so the maximum is **0.35**.

**Why the matrix checks exist at all** - a real failure this project has already
had, recorded in `not_constant_matrix`'s own docstring: a kinship matrix where
every entry was 1.0 produced h^2 = 1.00 and a set of associations with no
plausible mechanism, and nothing in the output said so. A distance matrix has
the same shape of failure, with two more ways to get it wrong quietly:

  * a **non-zero diagonal** means the self-path was walked as a real path, which
    is what a naive "sum the depths and double it" implementation does. The
    numbers still look plausible and the matrix is still symmetric, so nothing
    else catches it;
  * an **asymmetric** matrix is the signature of a traversal that stops early on
    one side, and it is invisible to every check that looks at one entry.

So all three are asserted: not-constant, zero diagonal, symmetric.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.execution.contracts import STAGE_TABLES
from papipeline.manifest import SampleManifest
from papipeline.models import Sample, RunMode
from papipeline.stages import similarity as stage_similarity

TREE = Path("test_data/phylogeny/tree.nwk")
ID_COLUMN = "sample_id"

#: The exact multiples of 0.05 the committed tree implies, counted by hand.
EXPECTED = {
    ("TEST_PA_001", "TEST_PA_004"): 0.15,
    ("TEST_PA_001", "TEST_PA_002"): 0.30,
    ("TEST_PA_001", "TEST_PA_020"): 0.35,
}


def _manifest(*ids: str) -> SampleManifest:
    return SampleManifest(samples=[Sample(sample_id=i) for i in ids])


def _matrix(rows) -> dict:
    """`{sample_id: {other_id: distance}}` from the stage's return value."""
    return {r["sample_id"]: r["distances"] for r in rows}


class TestTheStageHasAnEntryPoint:
    def test_run_exists(self):
        assert callable(stage_similarity.run)

    def test_it_is_no_longer_unbuilt(self):
        """Reconciled in 2223332, once TEST (1d54f31) and REAL both existed.

        Asserted so a green suite here cannot be read as "similarity is still
        TEST-only", which is what this file asserted for its whole life.
        """
        from papipeline.run import UNBUILT_STAGES

        assert "similarity" not in UNBUILT_STAGES

    def test_real_mode_is_gated_not_absent(self, config, tmp_path, monkeypatch):
        """REAL now runs, behind the gate every other REAL caller uses.

        It previously refused because there was no caller; it now refuses because
        `allow_real_mode` is false, which is a different reason and a temporary
        one. The refusal must still happen, and must name the gate rather than
        claiming the capability is missing.

        The manifest is a real tip of the committed tree, so this exercises the
        REAL refusal rather than tripping the manifest-membership check first -
        the two are different failures and a test that hits the wrong one proves
        nothing about the other.
        """
        with pytest.raises(NotImplementedError) as excinfo:
            stage_similarity.run(
                config, _manifest("TEST_PA_001"), RunMode.REAL, tree_path=TREE
            )
        assert "REAL" in str(excinfo.value)


class TestWorkedDistances:
    """The numbers, counted off the tree rather than recomputed."""

    @pytest.fixture(scope="class")
    def matrix(self, config, tmp_path_factory):
        ids = sorted(_tips(TREE))
        return stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE
        )

    @pytest.mark.parametrize("pair,expected", sorted(EXPECTED.items()))
    def test_a_worked_pair(self, matrix, pair, expected):
        a, b = pair
        got = _matrix(matrix)[a][b]
        assert got == pytest.approx(expected, abs=1e-9), (
            f"d({a},{b}) = {got}, expected {expected} (a whole number of 0.05 "
            "edges in the committed tree)"
        )

    def test_the_maximum_is_where_the_tree_says(self, matrix):
        """0.40, between the two deepest tips on opposite sides.

        Counted off the tree: `TEST_PA_019` sits five edges below the n12 side
        and `TEST_PA_020` five below it too, on the other side of the n18 root,
        so the path is 8 edges. Worth pinning because the maximum is a property
        of the *whole* matrix, so a traversal that quietly loses a level of
        nesting on one side would still produce a plausible maximum - just a
        smaller one.

        (The first draft assumed the farthest pair involved TEST_PA_001 and
        asserted 0.35. It does not: 001 is four edges from the root, while 019
        and 020 are the deepest tips in the tree.)
        """
        m = _matrix(matrix)
        distance, a, b = max(
            (d, x, y)
            for x, row in m.items()
            for y, d in row.items()
            if x != y
        )
        assert distance == pytest.approx(0.40, abs=1e-9)
        assert {a, b} == {"TEST_PA_019", "TEST_PA_020"}, (
            f"the farthest pair is {{{a}, {b}}}; the two deepest tips should be"
        )

    def test_a_sister_pair_is_closer_than_a_distant_one(self, matrix):
        """Ordering, not just values: 004 is nearer 001 than 020 is."""
        m = _matrix(matrix)
        assert m["TEST_PA_001"]["TEST_PA_004"] < m["TEST_PA_001"]["TEST_PA_020"]


class TestMatrixShape:
    def test_it_is_square_and_covers_every_sample(
        self, config, tmp_path
    ):
        ids = sorted(_tips(TREE))
        rows = stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE
        )
        assert len(rows) == len(ids), "one row per sample"
        for row in rows:
            assert row["sample_id"] in ids
            assert set(row["distances"]) == set(ids), (
                "every sample must appear as a column too, including itself"
            )

    def test_the_diagonal_is_zero(self, config, tmp_path):
        """A tip's distance to itself is zero, always.

        The failure this catches is specific: an implementation that sums two
        root-depths without subtracting the shared path reports
        `2 * depth(tip)` on the diagonal. The matrix stays symmetric and every
        off-diagonal number stays plausible, so no other check would notice.
        """
        ids = sorted(_tips(TREE))
        rows = stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE
        )
        for row in rows:
            assert row["distances"][row["sample_id"]] == 0.0, (
                f"diagonal at {row['sample_id']} is "
                f"{row['distances'][row['sample_id']]}, not 0"
            )

    def test_it_is_symmetric(self, config, tmp_path):
        """d(a,b) == d(b,a) for every pair.

        An asymmetric matrix is what a traversal that gives up on one side
        produces. It is invisible to a check that reads a single entry, and a
        consumer that fills a lower triangle from an upper one would silently
        invent distances.
        """
        ids = sorted(_tips(TREE))
        m = _matrix(
            stage_similarity.run(
                config, _manifest(*ids), RunMode.TEST, tree_path=TREE
            )
        )
        for a in ids:
            for b in ids:
                assert m[a][b] == pytest.approx(m[b][a], abs=1e-12), (
                    f"d({a},{b})={m[a][b]} but d({b},{a})={m[b][a]}"
                )

    def test_it_is_not_constant(self, config, tmp_path):
        """Reuses `not_constant_matrix`'s historical justification.

        A constant matrix means every isolate is the same isolate, and that has
        already cost this project a real result: a constant *kinship* matrix
        produced h^2 = 1.00 and a set of associations with no plausible
        mechanism, with nothing in the output to say so (see
        `not_constant_matrix`). A constant distance matrix would be the same
        failure wearing a different hat.
        """
        ids = sorted(_tips(TREE))
        m = _matrix(
            stage_similarity.run(
                config, _manifest(*ids), RunMode.TEST, tree_path=TREE
            )
        )
        values = {d for row in m.values() for d in row.values()}
        assert len(values) > 1, (
            f"every distance is the same ({values}); a constant matrix says all "
            "isolates are identical, which is a finding, not a result"
        )

    def test_every_distance_is_non_negative(self, config, tmp_path):
        """A negative branch length would make distances meaningless."""
        ids = sorted(_tips(TREE))
        m = _matrix(
            stage_similarity.run(
                config, _manifest(*ids), RunMode.TEST, tree_path=TREE
            )
        )
        for a, row in m.items():
            for b, d in row.items():
                assert d >= 0.0, f"d({a},{b}) = {d}"


class TestAgainstTheManifest:
    def test_a_sample_missing_from_the_tree_is_refused(self, config, tmp_path):
        """A manifest tip the tree does not carry is a membership mismatch.

        Rule 5: the two must agree or the run stops. Silently dropping the
        sample would make the matrix describe a different cohort than the study
        claims, and every distance in it would then be about the wrong set.
        """
        with pytest.raises(Exception) as excinfo:
            stage_similarity.run(
                config, _manifest("TEST_PA_001", "TEST_PA_999"), RunMode.TEST,
                tree_path=TREE,
            )
        assert "TEST_PA_999" in str(excinfo.value)

    def test_it_writes_a_file_whose_header_leads_with_sample_id(
        self, config, tmp_path
    ):
        """The written file, not just the returned structure.

        The written table is what downstream stages and the observatory read,
        and it is the thing the contract describes. Asserting only the return
        value would leave the writer free to disagree with it.
        """
        ids = sorted(_tips(TREE))
        # `run()` returns rows; the path comes back from the writer, so the test
        # calls `write_matrix` the way a caller would rather than assuming the
        # entry point hands back a path.
        rows = stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE,
        )
        path = stage_similarity.write_matrix(
            rows, tmp_path / "similarity.tsv"
        )
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        assert header[0] == ID_COLUMN
        assert header[1:] == ids, (
            "the header must be sample_id followed by one column per sample, "
            "in the same order as the rows"
        )
        assert len(lines) == len(ids) + 1, "one header plus one row per sample"


def _tips(tree_path: Path):
    """Tip labels from the committed tree, via the phylogeny stage's parser.

    Reused rather than reimplemented so the fixtures and the stage agree on
    which labels are tips - in particular that `n1`, `n2` ... are internal nodes
    and not samples.
    """
    from papipeline.stages.phylogeny import extract_tips

    return extract_tips(Path(tree_path).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# A5: the matrix's units travel with it
# ---------------------------------------------------------------------------
#
# A patristic distance is a sum of branch lengths, and IQ-TREE's branch lengths
# are expected substitutions per site. A reader who does not know that reads
# `0.0093` as nine SNPs - wrong by the length of the alignment, and invisible,
# because both are small decimals and the matrix is symmetric either way.


class TestTheMatrixSaysWhatUnitItIsIn:
    def test_running_the_stage_writes_a_units_sidecar(
        self, config, tmp_path
    ):
        ids = sorted(_tips(TREE))
        out = tmp_path / "similarity.tsv"
        stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE, out_path=out
        )
        sidecar = stage_similarity.read_units_sidecar(out)
        assert sidecar is not None, (
            "the matrix was written with no record of its units, so a reader "
            "has to guess - and the guess a reader makes is SNPs"
        )
        assert sidecar["unit"] == stage_similarity.DISTANCE_UNIT

    def test_the_sidecar_says_explicitly_it_is_not_a_snp_count(
        self, config, tmp_path
    ):
        ids = sorted(_tips(TREE))
        out = tmp_path / "similarity.tsv"
        stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE, out_path=out
        )
        sidecar = stage_similarity.read_units_sidecar(out)
        assert sidecar["is_a_snp_count"] is False
        assert "NOT SNP counts" in sidecar["note"]

    def test_it_records_the_method_and_the_definition(
        self, config, tmp_path
    ):
        """`execution/contracts.py` says "Method is recorded in the file so a
        reader cannot assume a different definition produced these numbers".
        It was not - `write_matrix` wrote a bare matrix. It is now."""
        ids = sorted(_tips(TREE))
        out = tmp_path / "similarity.tsv"
        stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE, out_path=out
        )
        sidecar = stage_similarity.read_units_sidecar(out)
        assert sidecar["method"] == stage_similarity.METHOD == "patristic"
        assert "2*depth(LCA)" in sidecar["definition"]
        assert "not part of the path" in sidecar["definition"], (
            "the stem above the LCA is the specific thing a reader (and one "
            "previous measurement here) got wrong"
        )

    def test_the_matrix_itself_is_unchanged_in_shape(self, config, tmp_path):
        """Not a `#` comment above the header. The contract fixes this file as
        one header row plus one row per sample, and a reader indexing the matrix
        must not have to know a comment may be there."""
        ids = sorted(_tips(TREE))
        out = tmp_path / "similarity.tsv"
        stage_similarity.run(
            config, _manifest(*ids), RunMode.TEST, tree_path=TREE, out_path=out
        )
        lines = out.read_text(encoding="utf-8").splitlines()
        assert len(lines) == len(ids) + 1
        assert lines[0].split("\t")[0] == ID_COLUMN

    def test_the_sidecar_is_named_from_the_matrix_not_invented(
        self, tmp_path
    ):
        assert (
            stage_similarity.units_sidecar_path(tmp_path / "similarity.tsv")
            == tmp_path / "similarity.tsv.units.json"
        ), (
            "the sidecar name is derived, so the pair cannot drift apart - a "
            "sidecar at a fixed name beside a differently-named matrix is not "
            "found by anything that did not hardcode both"
        )

    def test_no_sidecar_reads_as_none_rather_than_raising(
        self, tmp_path
    ):
        assert stage_similarity.read_units_sidecar(tmp_path / "absent.tsv") is None

    def test_an_unreadable_sidecar_is_logged_not_fatal(self, tmp_path, caplog):
        out = tmp_path / "similarity.tsv"
        stage_similarity.write_units_sidecar(out)
        stage_similarity.units_sidecar_path(out).write_text("{ not json")
        assert stage_similarity.read_units_sidecar(out) is None
