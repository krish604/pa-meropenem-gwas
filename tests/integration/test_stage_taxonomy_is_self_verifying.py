"""The stage taxonomy must describe the wiring that actually exists.

**Why this file exists.** Six bugs of one shape have appeared in this codebase:
a stage that was built but still declared absent, or declared present while its
reality was something else. `similarity` had a TEST path (1d54f31) and a REAL
caller (2223332) and was still in `UNBUILT_STAGES`, so two places described one
fact and the stage could be half-reconciled. `virulence` read a table nothing
wrote, so a REAL run had nothing to read (7e28840). The tests that noticed were
`sixteen_stages_complete` and `test_stage_taxonomy.py`, both of which hard-coded
counts and passed while the code was wrong.

A hard-coded count is a tripwire on the wrong thing. It fires when a stage is
promoted and says nothing about whether the promoted stage works. So the checks
here derive everything from `run.py`'s own dispatch chain and from the stage
modules' AST, and assert only invariants.

**The taxonomy has three states, not two.** Treating "runnable" as a single
state is what let the bugs hide:

1. *runnable with a REAL caller* - `similarity`. Its `run()` does raise
   `NotImplementedError`, but only when `runtime.allow_real_mode` is false. That
   is a gate on the run, not a refusal of the stage, and the stage proceeds once
   the gate is open.
2. *runnable in TEST, refusing REAL by design* - `gwas`, `convergence`,
   `cooccurrence`. Each refuses rather than silently producing a
   `ReferenceEngine` result that cannot detect the lineage confounding a real
   cohort (f90d6c9, 46f5005). The refusal is the correct behaviour and they are
   correctly runnable in TEST.
3. *unbuilt* - a stage with no dispatch branch at all. **Empty as of the
   `dag-resolve` reconciliation**, which promoted `recombination` (the last
   member) by adding `run.derive_recombination_tables` and dropping the name in
   the same commit. The set and the mechanism are kept: they are the right
   answer for the next declared-but-unbuilt stage, and the invariants below are
   written against the mapping rather than against its current contents, so they
   still say something when it is non-empty again.

States 2 and 3 were both reported as "runnable", which is why "runnable" could
not be tested. `REAL_REFUSING_STAGES` in `run.py` now names state 2 explicitly,
and the invariant below is that the declared set equals the set of stage
modules that refuse REAL - derived from the AST, not restated.

**What this does and does not catch.** It catches the taxonomy-shaped bugs:
a stage promoted without a dispatch branch, a stage dispatched without being
promoted, a stage declared runnable while its `run()` refuses REAL, a stand-in
taxonomy that drifts, and any count that stops matching the wiring. It does not
catch correctness bugs inside a stage - the blastn `qlen`/`slen` mix-up, the
`HSI-2` header split and the factor-versus-locus count all returned *plausible
numbers* and no amount of taxonomy introspection would have seen them. Those
needed running against the real database.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from papipeline.run import (
    EXECUTION_ORDER,
    REAL_REFUSING_STAGES,
    STAGE_ORDER,
    UNBUILT_STAGES,
)
from papipeline.standins import UNBUILT_WITH_STANDIN

RUN_PY = Path(__file__).resolve().parents[2] / "papipeline" / "run.py"
STAGES_DIR = Path(__file__).resolve().parents[2] / "papipeline" / "stages"

#: Stage modules that `run.py`'s dispatch chain names. `structural_variants`,
#: `mechanisms`, `regulators` and `integration` are internal steps folded into
#: the stage that owns them (spec.md:351) and have no place in STAGE_ORDER.
FOLDED_INTO_OWNER = {
    "mechanisms",
    "regulators",
    "structural_variants",
    "integration",
}


def _dispatched_stages() -> set:
    """Stages `run.py` actually dispatches, read from its AST.

    Parsed rather than grepped: a mention of a stage name in a comment, a
    prerequisite map or a log message is not a dispatch, and a string search
    cannot tell those apart. Only comparisons against the loop variable named
    `stage` count.
    """
    tree = ast.parse(RUN_PY.read_text(encoding="utf-8"))
    found = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        if not (
            isinstance(node.left, ast.Name) and node.left.id == "stage"
        ):
            continue
        for comparator in node.comparators:
            if isinstance(comparator, ast.Constant) and isinstance(comparator.value, str):
                found.add(comparator.value)
    return found


def _run_function(tree: ast.Module):
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "run":
            return node
    return None


def _enclosing_ifs(node, function) -> list:
    """Every `if` statement inside `function` that encloses `node`."""
    found: list = []

    def walk(current, stack):
        for child in ast.iter_child_nodes(current):
            if isinstance(child, ast.If):
                if any(grandchild is node for grandchild in ast.walk(child)):
                    found.append(child)
                if walk(child, stack + [child]):
                    return True
            elif walk(child, stack):
                return True
        return False

    walk(function, [])
    return found


def _mentions_mode(test: ast.AST) -> bool:
    """True when a guard's condition is about `RunMode`.

    This is the distinction the whole classification rests on. A raise inside
    `if mode is not RunMode.TEST:` refuses exactly REAL, which is state 2 and
    belongs in `REAL_REFUSING_STAGES`. A raise inside
    `if not config.runtime.get("allow_real_mode", False):` is a *gate on the
    run*: open the gate and the stage proceeds, so listing the stage as refusing
    REAL would be a false positive that invites "fixing" a stage that works.
    """
    for child in ast.walk(test):
        if isinstance(child, ast.Attribute) and child.attr in {"TEST", "REAL"}:
            return True
        if isinstance(child, ast.Name) and child.id in {"TEST", "REAL"}:
            return True
    return False


def _refuses_real(module: Path) -> bool:
    """True when `run()` refuses REAL.

    Two shapes count: an unconditional `raise NotImplementedError`, and one
    guarded by the mode. Anything guarded by something else is a gate.

    Parsed rather than searched for the string. A text search would classify
    `similarity` as refusing REAL - it raises `NotImplementedError`, inside the
    `allow_real_mode` gate - and the resulting false positive is worse than no
    check, because it points at a working stage.
    """
    tree = ast.parse(module.read_text(encoding="utf-8"))
    run_fn = _run_function(tree)
    if run_fn is None:
        return False
    for node in ast.walk(run_fn):
        if not isinstance(node, ast.Raise):
            continue
        exc = node.exc
        if not (
            isinstance(exc, ast.Call)
            and getattr(exc.func, "id", "") == "NotImplementedError"
        ):
            continue
        guards = _enclosing_ifs(node, run_fn)
        if not guards:
            return True  # unconditional
        # The INNERMOST guard decides, not any enclosing one. Nearly every
        # REAL body sits inside `if mode is not RunMode.TEST:`, so an
        # "any guard mentions mode" rule classifies every raise in a REAL
        # branch as a refusal - including the one inside similarity's
        # `allow_real_mode` gate, which is the opposite.
        innermost = max(guards, key=lambda guard: guard.lineno)
        if _mentions_mode(innermost.test):
            return True  # refuses exactly REAL
    return False


@pytest.fixture(scope="module")
def dispatched() -> set:
    return _dispatched_stages()


@pytest.fixture(scope="module")
def runnable() -> set:
    return set(STAGE_ORDER) - set(UNBUILT_STAGES)


class TestTheDispatchChainIsTheSourceOfTruth:
    def test_every_runnable_stage_is_dispatched(self, runnable, dispatched):
        """Requirement 1.

        A stage in `STAGE_ORDER` and outside `UNBUILT_STAGES`, with no branch in
        the dispatch chain, is declared runnable and never called. This is the
        `similarity` failure mode: promoted in the taxonomy, and 22 tests failed
        because nothing wrote `similarity.tsv`.
        """
        assert not runnable - dispatched, (
            "declared runnable but never dispatched: "
            f"{sorted(runnable - dispatched)}"
        )

    def test_no_stage_is_dispatched_without_being_declared_runnable(
        self, runnable, dispatched,
    ):
        """Requirement 2, the inverse.

        A dispatched stage absent from the taxonomy is executing outside it.
        `UNBUILT_STAGES` bypasses both the `analysis.*` flag check and the
        dispatch chain, so such a stage can run when the taxonomy says it is
        absent - the mirror image of the case above, and the reason both
        directions are asserted.
        """
        assert not dispatched - runnable, (
            "dispatched but not declared runnable: "
            f"{sorted(dispatched - runnable)}"
        )

    def test_no_folded_in_step_is_dispatched_as_a_stage(self, dispatched):
        """`structural_variants` and friends are steps inside their owner.

        They have no place in STAGE_ORDER, so a branch for one would mean the
        spec's folding had been undone in code while the docs still described it.
        """
        assert not (dispatched & FOLDED_INTO_OWNER), (
            f"folded-in steps dispatched as stages: {sorted(dispatched & FOLDED_INTO_OWNER)}"
        )

    def test_the_dispatch_covers_the_execution_order_exactly(self, dispatched):
        """`EXECUTION_ORDER` is documented as identical to `STAGE_ORDER` minus
        the unbuilt stages. Derived here rather than restated."""
        assert set(EXECUTION_ORDER) - set(UNBUILT_STAGES) == dispatched


class TestUnbuiltStagesHaveNoWiring:
    def test_an_unbuilt_stage_has_no_dispatch_branch(self, dispatched):
        assert not (set(UNBUILT_STAGES) & dispatched), (
            "declared unbuilt but dispatched: "
            f"{sorted(set(UNBUILT_STAGES) & dispatched)}"
        )

    def test_an_unbuilt_stage_has_no_working_implementation(self):
        """The check that actually catches the `similarity` bug.

        `UNBUILT_STAGES` and the dispatch chain agreed with each other while
        `similarity` was built: at 2223332 it had a TEST path (1d54f31) and a
        REAL caller, and `run.py` had no branch for it - so every check about
        the two agreeing passed. Both were consistently describing a stage that
        was in fact half-built.

        So this looks at the module itself. A stage declared unbuilt must either
        have no module or have a `run()` that refuses, because "no caller" and
        "a caller nobody invokes" are the same thing observed from two places.
        """
        for stage in sorted(UNBUILT_STAGES):
            module = STAGES_DIR / f"{stage}.py"
            if not module.exists():
                continue
            assert _refuses_real(module), (
                f"{stage} is declared unbuilt but papipeline/stages/{stage}.py "
                "has a run() that does not refuse - it has an implementation the "
                "taxonomy says is absent. Reconcile both, or drop the stage "
                "from UNBUILT_STAGES and wire the dispatch."
            )

    def test_an_unbuilt_stage_names_what_is_missing(self):
        """A bare tuple would leave the refusal unable to tell a reader what
        will build it."""
        for stage, reason in UNBUILT_STAGES.items():
            assert reason and reason.strip(), f"{stage} has no reason recorded"

    def test_the_standin_taxonomy_agrees(self):
        """`standins.UNBUILT_WITH_STANDIN` is a second description of the same
        fact, which is how `similarity` could be half-reconciled - present in one
        mapping and absent from the other.

        A stage may have no stand-in once it is built, so the stand-in set is a
        subset of the unbuilt set rather than equal to it; what must hold is that
        nothing claims a stand-in for a stage it considers runnable.
        """
        assert not (set(UNBUILT_WITH_STANDIN) - set(UNBUILT_STAGES)), (
            "runnable stages still claim a stand-in: "
            f"{sorted(set(UNBUILT_WITH_STANDIN) - set(UNBUILT_STAGES))}"
        )


class TestTheRealRefusalSetIsDerived:
    """State 2 of the taxonomy: runnable in TEST, refusing REAL by design."""

    def test_it_matches_the_stage_modules_that_refuse(self):
        """The declared set must equal what the code does.

        Each module in `UNBUILT_STAGES` is exempt: it has no caller and is not
        dispatched, so it has no `run()` to inspect.
        """
        modules = {}
        for module in sorted(STAGES_DIR.glob("*.py")):
            name = module.stem
            # `sv` is a folded-in step inside `amr` and has no place in the
            # taxonomy; it refuses REAL, and that refusal is amr's problem
            # to surface, not a 17th stage.
            if name in UNBUILT_STAGES or name.startswith("_"):
                continue
            if name in FOLDED_INTO_OWNER:
                continue
            if name not in set(STAGE_ORDER):
                continue
            if _refuses_real(module):
                modules[name] = True
        declared = set(REAL_REFUSING_STAGES)
        assert declared == set(modules), (
            "REAL_REFUSING_STAGES does not match reality: "
            f"declared {sorted(declared)}, actual {sorted(modules)}"
        )

    def test_a_refusing_stage_is_still_runnable_in_test(self):
        """Otherwise the refusal would be better expressed as unbuilt.

        `gwas`, `convergence` and `cooccurrence` each refuse REAL rather than
        silently returning a `ReferenceEngine` result that cannot detect the
        lineage confounding a real cohort. That refusal is the correct
        behaviour *and* they still have real TEST implementations, so they
        belong in state 2 - named - rather than being lumped with the unbuilt
        stages under a single "runnable" label.
        """
        assert set(REAL_REFUSING_STAGES) <= set(STAGE_ORDER)
        assert not (set(REAL_REFUSING_STAGES) & set(UNBUILT_STAGES)), (
            "a stage cannot both refuse REAL and be unbuilt; one states the "
            "capability is absent, the other that it was never built"
        )

    def test_recombination_is_not_treated_as_refusing(self):
        """The stage `dag-resolve` promoted, and the shape it now has.

        `recombination.run()` raises `NotImplementedError`, but only when
        `runtime.allow_real_mode` is false - the `similarity` shape. It raised
        unconditionally while the stage was unbuilt, which is why it *was* in
        `REAL_REFUSING_STAGES`; the dispatch branch in `run.py` and the change
        to an `allow_real_mode` gate moved it out. Asserted by name because this
        is the transition that file exists to make impossible to half-do, and a
        future reader who restores the unconditional raise will otherwise be
        pushed to add it back here.
        """
        assert "recombination" not in REAL_REFUSING_STAGES

    def test_similarity_is_not_treated_as_refusing(self):
        """The gate/refusal distinction, asserted on the stage most likely to
        be misread.

        `similarity.run()` raises `NotImplementedError`, but only when
        `runtime.allow_real_mode` is false. Listing it in `REAL_REFUSING_STAGES`
        would be a false positive that pushed the next reader to "fix" a stage
        that works.
        """
        assert "similarity" not in REAL_REFUSING_STAGES


class TestTheReportedCountsMatch:
    """Requirement 3. Derived, never hard-coded."""

    def test_the_taxonomy_has_sixteen_stages(self):
        """spec.md D1 (amended 2026-09-29) says sixteen. This one *is* a
        spec-derived number rather than a count of the implementation, so it is
        asserted literally; everything derived from it is not."""
        assert len(STAGE_ORDER) == 16

    def test_runnable_plus_unbuilt_is_the_whole_taxonomy(self):
        assert len(set(STAGE_ORDER) - set(UNBUILT_STAGES)) + len(UNBUILT_STAGES) == len(
            STAGE_ORDER
        )

    def test_the_runnable_count_equals_the_dispatch_count(self, dispatched):
        """The number that used to be written down.

        `test_pipeline.py` asserted `len(UNBUILT_STAGES) == 2` and
        `len(runnable) == 14`. Both held while `similarity` was built in reality
        and declared absent, and both failed the moment it was not.
        """
        assert len(set(STAGE_ORDER) - set(UNBUILT_STAGES)) == len(dispatched)

    def test_no_stage_is_declared_twice(self):
        assert len(STAGE_ORDER) == len(set(STAGE_ORDER))

    def test_execution_order_matches_stage_order(self):
        assert EXECUTION_ORDER == STAGE_ORDER

    def test_every_stage_module_that_is_dispatched_exists(self, dispatched):
        """A dispatch branch naming a module that does not exist would be an
        ImportError on the first run to reach it."""
        for stage in dispatched:
            assert (STAGES_DIR / f"{stage}.py").exists() or stage in {
                "reporting",
                "phenotype",
            }, f"{stage} is dispatched but papipeline/stages/{stage}.py is absent"