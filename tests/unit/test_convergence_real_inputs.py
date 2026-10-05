"""REAL-mode convergence and cooccurrence must refuse on INCOMPLETE input.

`46f5005` refused both stages unconditionally, because neither had a REAL
caller. That was right about the defect and wrong about the shape of the fix: it
made the stages unreachable even once their inputs are real.

The property this file pins is narrower and more useful than "it refuses":

* it refuses **by input name**, for the specific input that is absent, with the
  path it was expected at and what is holding it back;
* an input that is present but **degenerate** counts as absent - a lineages map
  where every value is ``unknown`` is a claim that every sample has a lineage
  called ``unknown``, not an absence of data, and treating it as data
  classifies every determinant ``UNKNOWN`` and reads as "no convergence found"
  (the shared predicate itself is tested in ``test_real_inputs.py``);
* it **does not fire** when every declared input has real-shaped content, and
  the stage then produces a real result. That last one matters most: a refusal
  that has never been observed to stand down is not known to be a refusal rather
  than a permanent refusal.

Both stages read their REAL inputs from **disk**, not from the in-memory
dictionaries a full run would pass them, because a stage run standalone or
under Snakemake has no upstream result to inherit. The existing pattern for
that is ``stages.amr.load_amr_table`` and ``stages.phylogeny.load_tree_metadata``,
both reused here rather than re-parsed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.errors import StageError
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import convergence as conv
from papipeline.stages import cooccurrence as cooc


# --------------------------------------------------------------------------
# Complete, real-shaped inputs. These are what "everything present" looks like.
# --------------------------------------------------------------------------

SAMPLES = ("S1", "S2", "S3", "S4", "S5", "S6")


def _manifest() -> SampleManifest:
    return SampleManifest(samples=[Sample(sample_id=s) for s in SAMPLES])


def _complete_lineages() -> dict:
    """Two lineages, so convergence has something to count independent of."""
    return {s: ("L1" if i < 3 else "L2") for i, s in enumerate(SAMPLES)}


AMR_COLUMNS = (
    "sample_id\tantibiotic\tdeterminant\tgene\tvariant\tdeterminant_type\t"
    "mechanism\tevidence_source\tdatabase\tdatabase_version\tconfidence\t"
    "claim_status\tidentity_pct\tcoverage_pct"
)


def _write_amr(path: Path, antibiotic: str = "imipenem") -> Path:
    """A stage-4 table with acquired determinants spread over two lineages."""
    rows = [AMR_COLUMNS]
    # blaOXA-1 in every sample of L1 and one of L2 -> two independent lineages.
    for s in ("S1", "S2", "S3", "S4"):
        rows.append(
            f"{s}\t{antibiotic}\tblaOXA-1\tblaOXA-1\t.\tAMR\t."
            f"\tamrfinderplus\tAMRFinderPlus\t2026-08-07.1\t.\tDETECTED\t99.0\t100.0"
        )
    # fosA in S4-S6, all L2: three carriers, so it clears `rare_max_samples`
    # (2) and reaches the lineage test rather than stopping at `rare_isolated`.
    # That is the branch worth exercising.
    for s in ("S4", "S5", "S6"):
        rows.append(
            f"{s}\t{antibiotic}\tfosA\tfosA\t.\tAMR\t."
            f"\tamrfinderplus\tAMRFinderPlus\t2026-08-07.1\t.\tDETECTED\t99.5\t100.0"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


REGULATOR_COLUMNS = (
    "sample_id\tgene\tvariant\treference\talternate\teffect\tmechanism\t"
    "confidence\tclaim_status"
)


def _write_regulators(path: Path) -> Path:
    rows = [REGULATOR_COLUMNS]
    # A real variant, not the `.` sentinel: `read_tsv` maps `.` to None, and a
    # row whose `variant` is None is skipped by `load_variant_names`. Writing
    # `.` here would produce a table that *looks* populated and yields an empty
    # mutation namespace - which is the degenerate case, not a passing fixture.
    for i, s in enumerate(SAMPLES, start=1):
        rows.append(
            f"{s}\toprD\tA{i}00T\tC\tT\ttruncating\tloss_of_function\t.\tDETECTED"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


MECHANISM_COLUMNS = (
    "sample_id\tantibiotic\tdeterminant\tmechanism\tevidence_level\tgene\t"
    "gene_type\tnotes"
)


def _write_mechanisms(path: Path, antibiotic: str = "imipenem") -> Path:
    rows = [MECHANISM_COLUMNS]
    for s in ("S1", "S2", "S3", "S4"):
        rows.append(f"{s}\t{antibiotic}\tblaOXA-1\tenzymatic_inhibition\tDETECTED\tblaOXA-1\tAMR\t.")
    # fosA mirrors the AMR fixture: S4-S6, all L2. Three carriers clears
    # `rare_max_samples` (2), so it reaches the lineage test rather than
    # stopping at `rare_isolated` - that is the branch worth exercising.
    for s in ("S4", "S5", "S6"):
        rows.append(f"{s}\t{antibiotic}\tfosA\tenzymatic_inactivation\tDETECTED\tfosA\tAMR\t.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


TREE_METADATA_COLUMNS = "sample_id\ttree_tip_label\tlineage_label\tst\tsource"


def _write_tree_metadata(path: Path, lineages: dict) -> Path:
    rows = [TREE_METADATA_COLUMNS]
    for s in SAMPLES:
        rows.append(f"{s}\t{s}\t{lineages[s]}\tST1\tiqtree")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def phylo_dir(tmp_path: Path) -> Path:
    """A stage-9 directory inside ``tmp_path``, not the configured one.

    ``config.phylogeny_dir(REAL)`` resolves into the worktree's real ``results/``
    tree. The first version of this fixture wrote there, so running the suite
    left a `tree_metadata.tsv` behind in `results/` - and then a "does this
    refuse?" check read that leftover and concluded lineages were present. Every
    file these tests write goes under ``tmp_path``.
    """
    return tmp_path / "phylogeny"


@pytest.fixture
def complete_intermediate(tmp_path: Path, phylo_dir: Path) -> Path:
    """An intermediate root where **every** declared REAL input is present.

    This is the fixture the "does not fire" tests build on. If the refusal were
    ever unconditional, these tests fail - which is the point: a guard that has
    never been observed to stand down is not known to be a guard.
    """
    root = tmp_path / "intermediate"
    _write_amr(root / "stages" / "04_amr.tsv")
    _write_regulators(root / "regulators" / "regulator_variants.tsv")
    _write_mechanisms(root / "stages" / "05_mechanisms.tsv")
    _write_tree_metadata(phylo_dir / "tree_metadata.tsv", _complete_lineages())
    return root


# --------------------------------------------------------------------------
# The shared predicate
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# convergence
# --------------------------------------------------------------------------


class TestConvergenceRefusesByInputName:
    def test_missing_regulator_variants_is_named(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        (complete_intermediate / "regulators" / "regulator_variants.tsv").unlink()
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "regulator_variants" in str(excinfo.value)

    def test_missing_lineages_is_named(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        (phylo_dir / "tree_metadata.tsv").unlink()
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        message = str(excinfo.value)
        assert "lineages" in message
        assert "tree_metadata.tsv" in message

    def test_an_all_unknown_lineage_map_counts_as_absent(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        """The degenerate case, and the one a naive check would miss.

        A file full of ``unknown`` is a stage *claiming* lineage data and having
        none. Every determinant would classify ``UNKNOWN`` and the table would
        read as "no convergence found".
        """
        _write_tree_metadata(
            phylo_dir / "tree_metadata.tsv", {s: "unknown" for s in SAMPLES}
        )
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "unknown" in str(excinfo.value)

    def test_missing_amr_table_is_named(self, config, tmp_path: Path):
        root = tmp_path / "intermediate"
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=root, antibiotic="imipenem",
            )
        assert "amr_calls" in str(excinfo.value)

    def test_no_intermediate_root_is_refused_not_defaulted(
        self, config
    ):
        """No root means nowhere to read from; falling back to the passed
        dictionaries is exactly the in-memory coupling being removed."""
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                amr_calls={"S1": ["blaOXA-1"]},
                regulator_variants={"S1": ["oprD:."]},
                lineages={"S1": "L1"},
            )
        assert "intermediate_root" in str(excinfo.value)

    def test_it_reports_every_missing_input_not_just_the_first(
        self, config, tmp_path: Path
    ):
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=tmp_path / "intermediate", antibiotic="imipenem",
            )
        assert excinfo.value.context["missing"].split(",") == [
            "amr_calls", "lineages", "regulator_variants",
        ]


class TestConvergenceRunsWhenInputsAreComplete:
    def test_it_produces_calls_and_does_not_refuse(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        calls = conv.run(
            config, _manifest(), RunMode.REAL,
            intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
        )
        assert calls, "complete inputs produced no convergence calls at all"
        by_name = {c.determinant: c for c in calls}
        # blaOXA-1 is in L1 (S1-S3) and L2 (S4): two independent lineages.
        assert by_name["blaOXA-1"].convergence_category.value == "recurrent_convergent"
        assert by_name["blaOXA-1"].independent_lineages == 2

    def test_the_amr_table_is_read_from_disk_not_the_argument(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        """The decoupling this whole change exists for.

        Passing in-memory mappings that contradict the table on disk must change
        nothing, because a standalone or Snakemake invocation has no in-memory
        upstream result to pass in the first place.
        """
        calls = conv.run(
            config, _manifest(), RunMode.REAL,
            amr_calls={"S1": ["a_gene_not_on_disk"]},
            regulator_variants={"S1": ["also:not:on:disk"]},
            lineages={"S1": "L9"},
            intermediate_root=complete_intermediate,
            antibiotic="imipenem",
            phylogeny_dir=phylo_dir,
            
        )
        names = {c.determinant for c in calls}
        assert "a_gene_not_on_disk" not in names
        assert "also:not:on:disk" not in names
        assert "blaOXA-1" in names

    def test_it_refuses_when_only_one_lineage_carries_a_determinant(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        """Present and complete, but the classification must still be right.

        Guards against the completeness check being the only thing tested: this
        is the input the *logic* has to get right once the guard stands down.
        """
        calls = conv.run(
            config, _manifest(), RunMode.REAL,
            intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
        )
        fos = next(c for c in calls if c.determinant == "fosA")
        # fosA is S5+S6, both L2 -> one lineage.
        assert fos.independent_lineages == 1
        assert fos.convergence_category.value == "lineage_associated"


# --------------------------------------------------------------------------
# cooccurrence
# --------------------------------------------------------------------------


class TestCooccurrenceRefusesByNamespace:
    def test_missing_mechanisms_is_named(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        (complete_intermediate / "stages" / "05_mechanisms.tsv").unlink()
        with pytest.raises(StageError) as excinfo:
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "mechanism" in str(excinfo.value)

    def test_missing_regulator_variants_is_named(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        (complete_intermediate / "regulators" / "regulator_variants.tsv").unlink()
        with pytest.raises(StageError) as excinfo:
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "mutation" in str(excinfo.value)

    def test_missing_lineages_is_named(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        (phylo_dir / "tree_metadata.tsv").unlink()
        with pytest.raises(StageError) as excinfo:
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "lineages" in str(excinfo.value)

    def test_an_all_unknown_lineage_map_counts_as_absent(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        """`test_pair` accepts lineages=None, so this would not raise itself.

        The consequence is not a crash but a quietly degraded table: the
        lineage-stratified counts the module docstring promises would be absent
        and nothing would say so.
        """
        _write_tree_metadata(
            phylo_dir / "tree_metadata.tsv", {s: "unknown" for s in SAMPLES}
        )
        with pytest.raises(StageError) as excinfo:
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert "lineages" in str(excinfo.value)

    def test_it_names_every_configured_namespace_that_is_absent(
        self, config, tmp_path: Path
    ):
        with pytest.raises(StageError) as excinfo:
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=tmp_path / "intermediate", antibiotic="imipenem",
            )
        missing = excinfo.value.context["missing"].split(",")
        # Config declares five tests spanning gene, mutation and mechanism, so
        # all three namespaces are required.
        assert set(missing) >= {"gene", "mutation", "mechanism"}


class TestCooccurrenceRunsWhenInputsAreComplete:
    def test_it_produces_pairs_and_does_not_refuse(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        results = cooc.run(
            config, _manifest(), RunMode.REAL,
            intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
        )
        assert results, "complete inputs produced no co-occurrence pairs at all"
        types = {r.feature_type for r in results}
        # genes, gene_mutation, gene_mechanism are all reachable now.
        assert {"gene_gene", "gene_mutation", "gene_mechanism"} <= types

    def test_the_tables_are_read_from_disk_not_the_arguments(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        results = cooc.run(
            config, _manifest(), RunMode.REAL,
            amr_genes={"S1": ["a_gene_not_on_disk"]},
            variants={"S1": ["nope:1"]},
            mechanisms={"S1": ["no_such_mechanism"]},
            lineages={"S1": "L9"},
            intermediate_root=complete_intermediate,
            antibiotic="imipenem",
            phylogeny_dir=phylo_dir,
            
        )
        features = {r.feature_a for r in results} | {r.feature_b for r in results}
        assert "a_gene_not_on_disk" not in features
        assert "no_such_mechanism" not in features
        assert "blaOXA-1" in features


# --------------------------------------------------------------------------
# The refusal must not leave a partial artefact
# --------------------------------------------------------------------------


class TestRefusalLeavesNothingBehind:
    def test_a_refused_convergence_run_writes_no_table(
        self, config, tmp_path: Path
    ):
        from papipeline.execution.contracts import table_path
        from papipeline.io.tsv import write_tsv

        stage_dir = tmp_path / "intermediate" / "stages"
        table = table_path(stage_dir, "convergence")
        assert not table.exists()
        with pytest.raises(StageError):
            conv.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=tmp_path / "intermediate", antibiotic="imipenem",
            )
        # run.py writes only after run() returns, so a raise means no file. This
        # asserts that rather than trusting it: the whole point of refusing is
        # that the artefact is absent.
        assert not table.exists()

    def test_the_stage_tables_written_before_the_refusal_are_untouched(
        self, config, complete_intermediate: Path, phylo_dir: Path
    ):
        """A refusal must not damage inputs it already read successfully."""
        amr = complete_intermediate / "stages" / "04_amr.tsv"
        before = amr.read_bytes()
        (complete_intermediate / "stages" / "05_mechanisms.tsv").unlink()
        with pytest.raises(StageError):
            cooc.run(
                config, _manifest(), RunMode.REAL,
                intermediate_root=complete_intermediate, antibiotic="imipenem",
                phylogeny_dir=phylo_dir,
            )
        assert amr.read_bytes() == before


# --------------------------------------------------------------------------
# Modes that must keep refusing
# --------------------------------------------------------------------------


class TestStubStillRefuses:
    @pytest.mark.parametrize(
        "stage,kwargs",
        [
            (conv, dict(amr_calls={}, regulator_variants={}, lineages={})),
            (cooc, dict(amr_genes={}, variants={}, mechanisms={})),
        ],
    )
    def test_stub_is_not_a_caller(self, config, stage, kwargs):
        """STUB fabricates these tables in `papipeline.stub`.

        Letting STUB reach the real computation would mean STUB had stopped
        stubbing, which would be a much harder thing to notice than a refusal.
        """
        with pytest.raises(NotImplementedError):
            stage.run(config, _manifest(), RunMode.STUB, **kwargs)