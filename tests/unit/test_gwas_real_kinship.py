"""Stage-10 kinship: the pyseer similarity matrix, built from the IQ-TREE tree.

What is under test is :mod:`papipeline.gwas_real.kinship` - the transformation
of the stage-9 tree's patristic distances into the similarity (kinship) matrix
pyseer's ``--lmm`` consumes (``--similarity``), plus a stage-shaped entry point
that writes it.

Three decisions are asserted here rather than described:

* **The distances are stage 10's, not a second traversal of the tree.** Spec
  user story 41 asks for a distance matrix computed from the same tree so that
  the kinship matrix and the tree cannot disagree. The way that is enforced is
  by calling ``stages.similarity.patristic_distances`` itself, and one test
  below spies on that call - a comment saying "we reuse it" would pass even if
  a second parser had crept in.

* **The transformation is Gower centring of the squared distances.** pyseer
  consumes a similarity matrix (``pyseer/lmm.py::initialise_lmm`` reads it with
  ``pd.read_table(index_col=0)`` and rescales so its trace equals the sample
  count), and it ships **no tree-to-similarity helper**: the installed package
  declares exactly five console scripts (``pyseer``, ``scree_plot_pyseer``,
  ``square_mash``, ``phandango_mapper``, ``annotate_hits_pyseer``) and
  ``python -m pyseer.similarity --help`` exists but takes
  ``--kmers/--vcf/--pres`` variants, not a tree. So the conversion is this
  pipeline's own, and it is the classical one:

      K = -1/2 * (D^2 - row_mean - col_mean + grand_mean)

  where ``D`` is the matrix of patristic distances. A tree metric is of
  negative type, so this ``K`` is positive semi-definite - the property
  ``pyseer/fastlmm/lmm_cov.py`` documents for its ``K`` ("positive semi-
  definite Kernel matrix") - and the PSD check in the module is the fail-closed
  assertion of that fact rather than an assumption. The worked example below is
  derived by hand in fractions, independently of the code.

* **A degenerate tree is refused, not passed through.** With every branch
  length zero every distance is zero, ``K`` is the zero matrix, and pyseer
  computes ``factor = float(len(p)) / np.diag(K.values).sum()``
  (``pyseer/lmm.py``) - a division by zero that surfaces as a tool traceback
  rather than as a statement about the tree.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.models import RunMode
from papipeline.stages import similarity as stage_similarity

from papipeline.gwas_real import kinship as kg

REPO = Path(__file__).resolve().parents[2]

#: The committed stage-9 fixture tree. All its branch lengths are 0.05, so
#: every patristic distance is a multiple of 0.05 (see test_similarity_stage).
FIXTURE_TREE = REPO / "test_data" / "phylogeny" / "tree.nwk"

#: Three tips with structure, small enough to work the transform out by hand.
THREE_TIP_TREE = "((A:0.1,B:0.1):0.2,C:0.3);"

#: Distances of that tree: d(A,B) = 0.1+0.1, d(A,C) = d(B,C) = 0.1+0.2+0.3.
THREE_TIP_IDS = ("A", "B", "C")


def _load_config(pipeline_root):
    """The project config, loaded fresh so a gate env var is re-read."""
    from papipeline.config.loader import load_config

    return load_config(pipeline_root / "config" / "science.yaml")


# --------------------------------------------------------------------------
# The transformation
# --------------------------------------------------------------------------


class TestTheTransformation:
    def test_the_distances_come_from_the_stage10_module(
        self, monkeypatch, pipeline_root
    ):
        """Reuse, asserted by spying on the call rather than by reading prose.

        Two traversals of the same Newick that agree today can drift apart the
        moment one is "optimised", and nothing downstream would notice: the
        kinship matrix would simply start describing a slightly different tree
        than stage 10 reports. The call itself is the contract.
        """
        real = stage_similarity.patristic_distances
        seen = {"calls": 0}

        def spy(tree_text):
            seen["calls"] += 1
            return real(tree_text)

        monkeypatch.setattr(stage_similarity, "patristic_distances", spy)
        kg.similarity_from_tree(THREE_TIP_TREE, THREE_TIP_IDS)
        assert seen["calls"] == 1, (
            "similarity_from_tree did not call stages.similarity."
            "patristic_distances exactly once; a second distance computation "
            "here is how the kinship matrix and the stage-10 matrix disagree"
        )

    def test_a_hand_computed_three_tip_tree(self):
        """The Gower transform, checked against arithmetic done by hand.

        With D^2 = [[0, 1/25, 9/25], [1/25, 0, 9/25], [9/25, 9/25, 0]]:

            row_mean(A) = (0 + 1/25 + 9/25) / 3 = 2/15
            row_mean(C) = (9/25 + 9/25 + 0)   / 3 = 6/25
            grand       = (2*(1/25) + 4*(9/25)) / 9 = 38/225

            K_AA = -1/2 * (0    - 2/15 - 2/15 + 38/225) = 11/225
            K_AB = -1/2 * (1/25 - 2/15 - 2/15 + 38/225) = 13/450
            K_AC = -1/2 * (9/25 - 2/15 - 6/25 + 38/225) = -7/90
            K_CC = -1/2 * (0    - 6/25 - 6/25 + 38/225) = 7/45

        The means run over ALL n entries including the diagonal zero, which is
        what makes -1/2 * D^2 double-centred rather than merely demeaned; the
        row sums come out to exactly zero, which is a property of the transform
        and not a coincidence of this example.

        ``K_AC`` is negative on purpose: a similarity is not required to be
        non-negative, only positive semi-definite, and two lineages far apart
        on the tree are less related than the cohort average.
        """
        K = kg.similarity_from_tree(THREE_TIP_TREE, THREE_TIP_IDS)
        assert K["A"]["A"] == pytest.approx(11 / 225, rel=1e-12)
        assert K["A"]["B"] == pytest.approx(13 / 450, rel=1e-12)
        assert K["A"]["C"] == pytest.approx(-7 / 90, rel=1e-12)
        assert K["C"]["C"] == pytest.approx(7 / 45, rel=1e-12)
        for row in K.values():
            assert sum(row.values()) == pytest.approx(0.0, abs=1e-12), (
                "double centring makes every row sum to zero; a row that does "
                "not means the transform is not the one documented"
            )

    def test_it_is_symmetric(self):
        K = kg.similarity_from_tree(THREE_TIP_TREE, THREE_TIP_IDS)
        for a in THREE_TIP_IDS:
            for b in THREE_TIP_IDS:
                assert K[a][b] == pytest.approx(K[b][a])

    def test_the_fixture_tree_gives_a_positive_semidefinite_matrix(self):
        """The property fastlmm documents for its K, checked on the real tree.

        A negative eigenvalue beyond float noise would make
        ``la.eigh`` in ``pyseer/fastlmm/lmm_cov.py`` see an indefinite kernel
        and the fit would be measuring the transform's bug rather than the
        cohort's relatedness.
        """
        import numpy as np

        newick = FIXTURE_TREE.read_text(encoding="utf-8")
        ids = [f"TEST_PA_{i:03d}" for i in range(1, 21)]
        K = kg.similarity_from_tree(newick, ids)
        matrix = np.array([[K[a][b] for b in ids] for a in ids])
        eigenvalues = np.linalg.eigvalsh(matrix)
        scale = max(1.0, float(np.abs(eigenvalues).max()))
        assert eigenvalues.min() >= -1e-8 * scale, (
            f"the kinship matrix has a negative eigenvalue {eigenvalues.min()}; "
            "Gower centring of a tree metric must be PSD"
        )

    def test_a_manifest_sample_missing_from_the_tree_is_refused(self):
        """A sample in the cohort and not in the tree is a mismatch, not a gap.

        Same reasoning as the stage's own tip-set check: silently dropping the
        sample would hand pyseer a kinship matrix describing a different cohort
        than the phenotype file does, and pyseer would quietly analyse the
        intersection.
        """
        with pytest.raises(DataContractError) as excinfo:
            kg.similarity_from_tree(THREE_TIP_TREE, ["A", "B", "C", "D"])
        message = str(excinfo.value)
        assert "D" in message
        assert "tree" in message.lower()

    def test_a_tree_with_only_zero_branch_lengths_is_refused(self):
        """Zero trace would divide by zero inside pyseer, not here.

        ``pyseer/lmm.py`` computes ``factor = float(len(p)) /
        np.diag(K.values).sum()``; with an all-zero ``K`` that is 0/0 for the
        caller, and the honest statement - this tree measures no relatedness -
        has to come from this pipeline rather than from a numpy traceback.
        """
        zero_tree = "(A:0.0,B:0.0,C:0.0);"
        with pytest.raises(DataContractError) as excinfo:
            kg.similarity_from_tree(zero_tree, THREE_TIP_IDS)
        message = str(excinfo.value)
        assert "trace" in message.lower()
        assert "pyseer" in message.lower()

    def test_the_matrix_follows_the_requested_order(self):
        """Row order is the caller's cohort order, not dictionary iteration."""
        forward = kg.similarity_from_tree(THREE_TIP_TREE, list(THREE_TIP_IDS))
        backward = kg.similarity_from_tree(
            THREE_TIP_TREE, list(reversed(THREE_TIP_IDS))
        )
        assert list(forward) == list(THREE_TIP_IDS)
        assert list(backward) == list(reversed(THREE_TIP_IDS))


# --------------------------------------------------------------------------
# The written file
# --------------------------------------------------------------------------


class TestTheWrittenMatrix:
    def test_pyseers_own_reader_reads_it_back(self, tmp_path):
        """The exact read pyseer performs, not a reader written for this test.

        ``pyseer/lmm.py::initialise_lmm`` does
        ``pd.read_table(K_in, index_col=0)`` then ``K.index.astype(str)``. If
        the file were not shaped for that, pyseer would fail deep inside its
        own fit with a pandas error instead of this pipeline saying so.
        """
        pd = pytest.importorskip("pandas")
        ids = list(THREE_TIP_IDS)
        K = kg.similarity_from_tree(THREE_TIP_TREE, ids)
        path = tmp_path / "kinship.tsv"
        kg.write_matrix(path, ids, K, tmp_path / "tree.nwk")

        frame = pd.read_table(path, index_col=0)
        frame.index = frame.index.astype(str)
        assert list(frame.index) == ids
        assert list(frame.columns) == ids
        assert float(frame.loc["A", "B"]) == pytest.approx(13 / 450, rel=1e-6)

    def test_the_sidecar_records_what_the_numbers_are(self, tmp_path):
        """A bare matrix of small decimals is exactly the ambiguous case.

        The convention comes from ``stages/similarity.py``: a matrix gets a
        sidecar saying its quantity, its units and what it is not, and the
        sidecar's name is derived from the matrix's so the pair cannot drift.
        """
        ids = list(THREE_TIP_IDS)
        K = kg.similarity_from_tree(THREE_TIP_TREE, ids)
        path = tmp_path / "kinship.tsv"
        kg.write_matrix(path, ids, K, tmp_path / "tree.nwk")

        sidecar = json.loads(
            tmp_path.joinpath("kinship.tsv.meta.json").read_text(encoding="utf-8")
        )
        assert sidecar["matrix"] == "kinship.tsv"
        assert sidecar["method"] == kg.METHOD
        assert sidecar["input_unit"] == stage_similarity.DISTANCE_UNIT
        assert sidecar["is_a_snp_count"] is False
        assert sidecar["source_tree"].endswith("tree.nwk")
        # The raw entries are squared substitutions per site; pyseer makes them
        # dimensionless itself with its trace rescaling.
        assert "trace" in sidecar["note"].lower()


# --------------------------------------------------------------------------
# The stage-shaped entry point
# --------------------------------------------------------------------------


class TestTheEntryPoint:
    def test_test_mode_writes_the_matrix_for_the_committed_tree(
        self, config, manifest, tmp_path
    ):
        """TEST runs on the committed fixtures, in manifest order."""
        out = tmp_path / "kinship.tsv"
        rows = kg.run(
            config,
            manifest,
            RunMode.TEST,
            tree_path=FIXTURE_TREE,
            out_path=out,
        )
        assert [r[kg.ID_COLUMN] for r in rows] == list(manifest.sample_ids)
        assert out.is_file()
        assert kg.units_sidecar_path(out).is_file()

        pd = pytest.importorskip("pandas")
        frame = pd.read_table(out, index_col=0)
        assert list(frame.index) == list(manifest.sample_ids)
        # Zero diagonal: a sample's distance to itself is zero even after the
        # transform, because the transform is a function of the distance matrix
        # and D_ii = 0 survives centring as a positive self-similarity - so the
        # *distance* file has a zero diagonal and this one does not. What must
        # hold is that the diagonal is the largest entry in its row for a tree
        # with distinct tips... which is not generally true either. What is
        # true and worth asserting: the diagonal is positive.
        assert (frame.values.diagonal() > 0).all()

    def test_test_mode_without_a_tree_path_is_refused(
        self, config, manifest, tmp_path
    ):
        """No default tree: resolving one silently would measure a tree
        nobody chose. Same rule as the stage's own TEST path."""
        with pytest.raises(DataContractError) as excinfo:
            kg.run(config, manifest, RunMode.TEST, out_path=tmp_path / "k.tsv")
        assert "tree_path" in str(excinfo.value)

    def test_real_mode_refuses_while_the_gate_is_closed(
        self, pipeline_root, manifest, tmp_path, monkeypatch
    ):
        from papipeline.config import loader as loader_module

        monkeypatch.delenv(loader_module.ALLOW_REAL_MODE_ENV, raising=False)
        config = _load_config(pipeline_root)
        assert config.runtime.get("allow_real_mode") is False, (
            "the committed overlays keep REAL shut; this test is what says so"
        )
        with pytest.raises(NotImplementedError) as excinfo:
            kg.run(config, manifest, RunMode.REAL, out_path=tmp_path / "k.tsv")
        message = str(excinfo.value)
        assert "allow_real_mode" in message
        assert "runtime.allow_real_mode" in message

    def test_real_mode_past_the_gate_reads_the_configured_tree(
        self, pipeline_root, manifest, tmp_path, monkeypatch
    ):
        """Open the real gate through the real mechanism, then hit the next
        fact: the tree is where configuration says stage 9 writes it, and a
        missing tree is a missing intermediate."""
        from papipeline.config import loader as loader_module

        monkeypatch.setenv(loader_module.ALLOW_REAL_MODE_ENV, "1")
        # Results redirected into the test's own directory so the assertion
        # does not depend on what a previous run left behind.
        monkeypatch.setenv(loader_module.RESULTS_ROOT_ENV, str(tmp_path / "results"))
        config = _load_config(pipeline_root)
        assert config.runtime.get("allow_real_mode") is True
        with pytest.raises(DataContractError) as excinfo:
            kg.run(config, manifest, RunMode.REAL, out_path=tmp_path / "k.tsv")
        message = str(excinfo.value)
        assert "tree.nwk" in message
        assert "similarity" in message.lower() or "kinship" in message.lower()
