"""similarity's REAL caller: stage-9's tree in, a real distance matrix out.

`stages.similarity` has had exactly one commit in its life - 1d54f31, "Build the
similarity stage's TEST path" - and has refused every non-TEST mode ever since.
It is TEST-only, and the direction is settled: read the tree stage 9 produced and
compute patristic distances from it. No second IQ-TREE invocation, because the
tree already exists and re-inferring it would be a different tree.

**The computation is not reimplemented.** Everything after the guard - read the
tree, `patristic_distances`, check the tip set against the manifest, emit rows -
is already mode-agnostic and correct. The only thing REAL needs is to resolve the
tree from where phylogeny writes it rather than requiring the caller to hand one
over. So the REAL branch is one path resolution.

**The tree under test is a real one**, produced by running the phylogeny stage:
IQ-TREE on an alignment, through the same caller the pipeline uses. Not a
hand-written Newick, because a hand-written tree would only prove the matrix code
reads *that* shape of string - and the stage-9 tree is not shaped like anything
one would write by hand (IQ-TREE's own labels, support values, its ordering).

**Worked values are counted independently.** The distance check re-derives a
path length from the Newick with a different algorithm - a recursive split on
balanced parentheses, summing the branch lengths on the path between two named
tips - rather than calling the module's depth/LCA traversal. Agreeing with an
independent method is evidence; agreeing with the same traversal twice is not.
"""

from __future__ import annotations

import copy
import re
import shutil
from pathlib import Path

import pytest

from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]
HAS_IQTREE = shutil.which("iqtree") is not None

#: Six sequences with structure, so the tree is not arbitrary and distances vary.
ALIGNMENT = """\
>S01
ACGTACGTAAGGCCTTACGGAT
>S02
ACGTACGTAAGGCCTTACGGAT
>S03
ACGTACGTAAGGCCTTACGGAT
>S04
ACGTTCGTAAGGCCTTACGGAT
>S05
ACGTTCGTAAGGCCTTACGGAT
>S06
ACGTTCGTAAGGCCTTACGGAT
>S07
AAGTTCGTAAGGCCTTACGGAT
>S08
AAGTTCGTAAGGCCTTACGGAT
"""


dendropy = pytest.importorskip("dendropy")


def _independent_path_length(newick: str, a: str, b: str) -> float:
    """Patristic distance from dendropy - an independent implementation.

    Two earlier versions of this check were hand-rolled Newick scanners and both
    were wrong, in ways that made the comparison worthless: one folded a closed
    branch into its siblings, so two zero-length tips at the root came out 0.0479
    apart; the other lost internal depths and read a genuinely distant pair as
    0.0. A checker that is itself buggy agrees with a buggy implementation about
    as often as with a correct one, so the independence has to come from software
    that is not this repository's.

    dendropy parses the tree its own way and implements patristic distance from
    the definition, so agreement between it and `patristic_distances` is evidence
    about both.
    """
    tree = dendropy.Tree.get_from_string(newick, schema="newick")
    matrix = tree.phylogenetic_distance_matrix()
    # Keyed on Taxon objects, not label strings.
    taxa = {x.label: x for x in matrix.taxon_namespace}
    return float(matrix.distance(taxa[a], taxa[b]))


@pytest.fixture
def real_mode_open(monkeypatch):
    """Open REAL for the test through the real mechanism.

    `PIPELINE_ALLOW_REAL_MODE`, not a patched config value: the gate is supposed
    to be the session override added for exactly this, and a test that bypassed
    it would not prove the gate reads what an operator sets.
    """
    monkeypatch.setenv("PIPELINE_ALLOW_REAL_MODE", "1")
    from papipeline.config.loader import load_config

    return load_config(REPO / "config" / "science.yaml", machine="laptop")


#: Stage 3's contracted table. `lineage_label` is the MLST sequence type
#: (ruling R4), so the phylogeny producer reads it and refuses - naming the
#: path - when it is absent. This supplies it inside the isolated run root.
MLST_HEADER = "sample_id\tST\talleles\tMLST_status\tmlst_scheme\tallele_database"


@pytest.fixture(scope="module")
def real_tree(tmp_path_factory):
    """A tree produced by the phylogeny stage, not written by this test."""
    if not HAS_IQTREE:
        pytest.skip("iqtree is not on PATH")
    from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample
    from papipeline.stages import phylogeny as stage_phylo

    out = tmp_path_factory.mktemp("phylo")
    (out / "core_snp_alignment.fasta").write_text(ALIGNMENT, encoding="utf-8")
    ids = re.findall(r">(\S+)", ALIGNMENT)
    manifest = SampleManifest(samples=[Sample(sample_id=i) for i in ids])

    # `PIPELINE_RESULTS_ROOT` moves `tool_output_root(REAL)` for the duration of
    # the call only: the producer must read stage 3's table at its contracted
    # path, and that path is under the run root rather than beside the tree.
    redirect = pytest.MonkeyPatch()
    try:
        run_root = tmp_path_factory.mktemp("run")
        redirect.setenv(RESULTS_ROOT_ENV, str(run_root))
        mlst = (
            run_root / RunMode.REAL.value.lower()
            / "intermediate" / "mlst" / "mlst_results.tsv"
        )
        mlst.parent.mkdir(parents=True, exist_ok=True)
        mlst.write_text(
            "\n".join([MLST_HEADER] + [
                f"{i}\t{100 + n}\t" + ";".join(f"locus{k}:1" for k in range(7))
                + f"\ttyped\tpaerMLST\tpubmlst.org"
                for n, i in enumerate(ids)
            ]) + "\n",
            encoding="utf-8",
        )

        config = load_config(REPO / "config" / "science.yaml", machine="laptop")
        prepared = copy.copy(config)
        raw = copy.deepcopy(dict(config.raw))
        raw["phylogeny"] = dict(raw.get("phylogeny") or {})
        raw["phylogeny"]["require_exact_sample_match"] = False
        # ALIGNMENT has ONE polymorphic column, so the configured `GTR+G+ASC`
        # and the constant-after-gaps filter cannot both apply to it: filtered,
        # it is a single site, and IQ-TREE exits 2 on that. This fixture is
        # about the similarity matrix reading a real tree, not about the model;
        # the configured pair is covered by
        # `tests/unit/test_iqtree_adapter.py::TestAscIsUsableOnceTheFilterHasRun`.
        raw["phylogeny"]["models"] = "GTR+G"
        raw["phylogeny"]["asc_drop_partially_constant"] = False
        object.__setattr__(prepared, "raw", raw)

        stage_phylo.run(prepared, manifest, RunMode.REAL, out)
    finally:
        redirect.undo()
    tree = out / "tree.nwk"
    assert tree.is_file(), "the phylogeny stage wrote no tree"
    return tree, ids, prepared


class TestTheRealCallerProducesARealMatrix:
    def test_it_reads_the_tree_phylogeny_wrote(self, real_tree, monkeypatch, real_mode_open):
        """The seam: phylogeny's output directory, not a caller-supplied path."""
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        monkeypatch.setattr(
            type(real_mode_open), "phylogeny_dir",
            lambda self, mode: tree.parent, raising=True,
        )
        rows = stage.run(real_mode_open, _manifest(ids), RunMode.REAL)
        assert len(rows) == len(ids)

    def test_the_matrix_is_square_with_a_zero_diagonal(self, real_tree, monkeypatch, real_mode_open):
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        monkeypatch.setattr(
            type(real_mode_open), "phylogeny_dir",
            lambda self, mode: tree.parent, raising=True,
        )
        rows = stage.run(real_mode_open, _manifest(ids), RunMode.REAL)
        for row in rows:
            distances = row["distances"]
            assert sorted(distances) == sorted(ids), "not one entry per cohort member"
            assert distances[row["sample_id"]] == 0.0

    def test_it_is_symmetric(self, real_tree, monkeypatch, real_mode_open):
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        monkeypatch.setattr(
            type(real_mode_open), "phylogeny_dir",
            lambda self, mode: tree.parent, raising=True,
        )
        matrix = {r["sample_id"]: r["distances"] for r in
                  stage.run(real_mode_open, _manifest(ids), RunMode.REAL)}
        for a in ids:
            for b in ids:
                assert matrix[a][b] == pytest.approx(matrix[b][a]), (
                    f"d({a},{b}) != d({b},{a})"
                )

    def test_it_is_not_constant(self, real_tree, monkeypatch, real_mode_open):
        """A matrix of identical values would satisfy every other check here.

        This is the one that would catch a caller returning zeros - which reads
        as "every isolate is identical to every other", the failure the original
        refusal existed to prevent.
        """
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        monkeypatch.setattr(
            type(real_mode_open), "phylogeny_dir",
            lambda self, mode: tree.parent, raising=True,
        )
        rows = stage.run(real_mode_open, _manifest(ids), RunMode.REAL)
        off_diagonal = [
            row["distances"][other]
            for row in rows for other in ids if other != row["sample_id"]
        ]
        assert len(set(off_diagonal)) > 1, "every distance is identical"


class TestWorkedValuesCountedIndependently:
    def test_a_distance_matches_an_independent_path_walk(self, real_tree, monkeypatch, real_mode_open):
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        monkeypatch.setattr(
            type(real_mode_open), "phylogeny_dir",
            lambda self, mode: tree.parent, raising=True,
        )
        rows = stage.run(real_mode_open, _manifest(ids), RunMode.REAL)
        matrix = {r["sample_id"]: r["distances"] for r in rows}

        checked = 0
        for i, a in enumerate(ids):
            for b in ids[i + 1:]:
                expected = _independent_path_length(
                    tree.read_text(encoding="utf-8"), a, b
                )
                assert matrix[a][b] == pytest.approx(expected, abs=1e-6), (
                    f"d({a},{b}) = {matrix[a][b]} but an independent walk of the "
                    f"same Newick gives {expected}"
                )
                checked += 1
        assert checked >= 20, "too few pairs checked to be meaningful"

    def test_sibling_tips_at_the_root_are_zero_distance(self, real_tree):
        """A value read straight off the tree text, with no computation.

        IQ-TREE hangs the first and last sequences directly off an unnamed root
        with 0.0000000000 branches, so the path between them is empty and the
        distance is 0.0 by inspection. This pins the *checker*: both earlier
        hand-rolled versions got this pair wrong, and a checker that is itself
        wrong would have made the comparison above meaningless.
        """
        tree, _, _ = real_tree
        text = tree.read_text(encoding="utf-8")
        first = re.match(r"\(([A-Za-z0-9_.-]+):[0-9.]+", text).group(1)
        last = re.findall(r",([A-Za-z0-9_.-]+):[0-9.]+\);$", text)[-1]
        # Both are the root's immediate children, and both branches are zero.
        assert re.match(rf"\({re.escape(first)}:0\.0*[,)]", text), (
            f"{first} is not a zero-length child of the root"
        )
        assert re.search(rf",{re.escape(last)}:0\.0*\);", text), (
            f"{last} is not a zero-length child of the root"
        )
        assert _independent_path_length(text, first, last) == pytest.approx(0.0)


class TestRealIsGated:
    def test_it_refuses_when_real_mode_is_shut(self, real_tree, monkeypatch):
        """The gate has to be here, not in a script.

        Every other REAL caller's `allow_real_mode` check lives in its
        `scripts/<stage>/run_<stage>.py`. There is no `scripts/similarity/` - the
        Snakefile rule points at a runner nobody wrote - so the stage is the only
        place a gate can live until that script exists.
        """
        from papipeline.stages import similarity as stage

        tree, ids, config = real_tree
        shut = copy.copy(config)
        runtime = dict(config.runtime)
        runtime["allow_real_mode"] = False
        object.__setattr__(shut, "_runtime", runtime)
        monkeypatch.setattr(
            type(config), "phylogeny_dir", lambda self, mode: tree.parent, raising=True
        )
        monkeypatch.delenv("PIPELINE_ALLOW_REAL_MODE", raising=False)
        with pytest.raises(Exception) as excinfo:
            stage.run(shut, _manifest(ids), RunMode.REAL)
        assert "allow_real_mode" in str(excinfo.value)

    def test_test_mode_still_reads_a_supplied_tree(self, real_tree):
        """The guard is a branch, not a replacement."""
        from papipeline.stages import similarity as stage

        tree, ids, _ = real_tree
        rows = stage.run(_laptop(), _manifest(ids), RunMode.TEST, tree_path=tree)
        assert len(rows) == len(ids)


def _laptop():
    from papipeline.config.loader import load_config

    return load_config(REPO / "config" / "science.yaml", machine="laptop")


def _manifest(ids):
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    return SampleManifest(samples=[Sample(sample_id=i) for i in ids])
