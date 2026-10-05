"""similarity is reconciled: both modes proven, so it leaves the unbuilt taxonomy.

Two commits established that, in this order and neither before the other:

* `1d54f31` built the TEST path - a committed tree in, a square matrix out;
* `2223332` built the REAL caller and, in the course of it, fixed a Newick
  parser that could not read a tree from a real tree builder.

Until now `similarity` was listed in two places - `run.UNBUILT_STAGES`, which
refuses it outside STUB mode, and `standins.UNBUILT_WITH_STANDIN`, which fed the
DAG a fabricated table. Both existed to describe the same fact, so a stage could
be half-reconciled: built in reality, still declared absent. That is the state
this test exists to prevent.

**The taxonomy is asserted against the code, not restated.** A test that
hard-codes `("recombination", "similarity")` passes today and lies the moment a
stage is reconciled - it would keep asserting a stage is unbuilt after it has a
caller. `test_the_taxonomy_is_not_hardcoded_anywhere` pins that: the lists the
tests iterate are derived from `UNBUILT_STAGES`, so adding a caller and dropping
the name is the whole change.

The stand-in fixture goes with it. `test_data/standins/unbuilt_stages/similarity.tsv`
exists so the DAG has *something* to read while the stage is absent; leaving it
behind would leave a fabricated matrix on disk for a stage that now computes a
real one, which is precisely the confusion the stand-in mechanism exists to
prevent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.run import UNBUILT_STAGES
from papipeline.standins import UNBUILT_WITH_STANDIN

REPO = Path(__file__).resolve().parents[2]
STANDIN_DIR = REPO / "test_data" / "standins" / "unbuilt_stages"


class TestSimilarityIsNoLongerUnbuilt:
    def test_it_is_not_in_the_run_taxonomy(self):
        assert "similarity" not in UNBUILT_STAGES, (
            "similarity has a TEST path (1d54f31) and a REAL caller (2223332); "
            "leaving it in UNBUILT_STAGES refuses a stage that now works"
        )

    def test_it_is_not_stood_in_for(self):
        assert "similarity" not in UNBUILT_WITH_STANDIN, (
            "similarity is stood in for only while it has no implementation; the "
            "DAG must now get a real matrix"
        )

    def test_the_two_lists_agree_with_each_other(self):
        """They describe one fact, so they must not drift.

        A stage in one and not the other is refused by `UNBUILT_STAGES` yet still
        handed a fabricated table - the DAG would read a stand-in for a stage the
        pipeline says does not exist, and the refusal would never be reached.
        """
        assert set(UNBUILT_WITH_STANDIN) == set(UNBUILT_STAGES), (
            f"stood in for {sorted(set(UNBUILT_WITH_STANDIN))} but unbuilt "
            f"{sorted(set(UNBUILT_STAGES))}"
        )

    def test_the_standin_is_declared_but_never_read(self):
        """The file stays; the *use* of it does not.

        Worth being explicit, because it looks like an oversight and is not.
        `variants` was reconciled first, and its stand-in file is still on disk
        and still named as a rule input in the Snakefile - declared so the DAG
        resolves, never read, because the stage reads real output. Deleting the
        file instead breaks the TEST DAG, which is exactly what happened when I
        tried it. So this asserts the declared/used distinction, not absence.
        """
        assert (STANDIN_DIR / "similarity.tsv").is_file(), (
            "the Snakefile still names this path as a rule input; deleting it "
            "breaks the TEST DAG, as it did for `variants`"
        )
        source = (
            REPO / "papipeline" / "stages" / "similarity.py"
        ).read_text(encoding="utf-8")
        assert "unbuilt_stages" not in source, (
            "the stage must not read a stand-in; it resolves the matrix from "
            "stage 9's tree"
        )


class TestSimilarityIsWantedNotMerelyPresent:
    def test_the_stage_is_enabled_in_the_science_config(self, config):
        """Reconciliation alone leaves the stage switched off.

        `UNBUILT_STAGES` bypassed the `analysis.*` flag check, so a stage dropped
        from it falls back to that flag - and a stage disabled there executes
        nothing. The `similarity` Snakemake rule failed in STUB with "Stage
        similarity is disabled in configuration" until this key was added, which
        is the same trap `variants` recorded when it was reconciled.
        """
        assert config.raw["analysis"].get("similarity") is True

    def test_it_is_still_a_stage_in_its_own_right(self, config):
        """Not demoted to an internal step of another stage.

        spec.md:351 folds `mechanisms`, `regulators`, `structural_variants`
        and `integration` into their parents. Reconciliation must not quietly do
        the same to similarity - it has its own rule and its own output.
        """
        from papipeline.run import STAGE_ORDER

        assert "similarity" in STAGE_ORDER


class TestTheTaxonomyIsNotHardcoded:
    def test_it_is_derived_where_it_is_iterated(self):
        """The pattern that makes reconciliation a one-line change.

        `test_standins.py` and `test_stage_taxonomy.py` both iterate a list. If
        either restates it, the next stage to be reconciled leaves a test
        asserting it is still absent.
        """
        for name in ("test_stage_taxonomy.py", "test_standins.py"):
            source = (REPO / "tests" / "integration" / name).read_text(
                encoding="utf-8"
            ) if (REPO / "tests" / "integration" / name).is_file() else (
                REPO / "tests" / "unit" / name
            ).read_text(encoding="utf-8")
            assert (
                'UNBUILT = ("recombination", "similarity")' not in source
            ), (
                f"{name} hard-codes the unbuilt taxonomy; derive it from "
                "UNBUILT_STAGES so reconciliation cannot leave it stale"
            )

    def test_recombination_is_no_longer_unbuilt_either(self):
        """The control, and it moved when `recombination` was reconciled.

        It was the standing control for this file: while `similarity` had been
        reconciled, `recombination` was the proof that reconciliation removes one
        name and not the concept. `dag-resolve` promoted `recombination` too, so
        the control has been retired rather than left asserting a stage is
        unbuilt after it has a caller - which is the exact failure mode this
        file's own module docstring warns about.

        What replaces it is not a restatement but a check on the *pair*: with
        both names gone, `UNBUILT_STAGES` and `UNBUILT_WITH_STANDIN` are both
        empty, and the taxonomy has one state left for a stage that is declared
        but not built. That is the state the next such stage will arrive in, and
        it is asserted below rather than assumed.
        """
        assert "recombination" not in UNBUILT_STAGES
        assert "recombination" not in UNBUILT_WITH_STANDIN
        # The stand-in *file* stays on disk, unreferenced. See the module
        # docstring of `papipeline/standins.py`: deleting it before the tuple
        # stops naming it breaks the TEST DAG.
        assert (STANDIN_DIR / "recombination.tsv").is_file()


class TestTheUnbuiltTaxonomyIsNowEmpty:
    """The state every stage is in, stated once so a change is deliberate."""

    def test_no_stage_is_declared_unbuilt(self):
        assert UNBUILT_STAGES == {}, (
            "every stage in STAGE_ORDER is built and dispatched. A name here "
            "means a stage `run.py` refuses outside STUB, and `tests/"
            "integration/test_stage_taxonomy_is_self_verifying.py` asserts the "
            "dispatch chain and this mapping agree in both directions - so "
            "adding one here without a branch (or with one, without removing "
            "the name) fails the suite."
        )

    def test_no_stage_is_stood_in_for(self):
        assert UNBUILT_WITH_STANDIN == (), (
            "a stand-in is a fabricated table satisfying a declared edge for a "
            "stage with no implementation. None is unbuilt, so none may be "
            "stood in for: a stage in this tuple that is runnable is feeding "
            "the DAG a table that was never computed."
        )
