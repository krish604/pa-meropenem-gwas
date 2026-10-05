"""convergence and cooccurrence must refuse REAL rather than half-run it.

Both accept `mode` and **never branch on it** - zero occurrences of
`RunMode.REAL` in either file. So on a REAL run they compute from whatever
intermediate tables happen to exist and return a result indistinguishable from
TEST mode's. That is worse than refusing, because:

* convergence is a *scientific claim* - that resistance arises independently in
  multiple lineages - and computing it over partial data reports a conclusion
  the data does not support;
* an empty or thin result reads as "no convergence found", which is a finding,
  where the truth is "this was never computed".

`gwas` and `similarity` already raise `NotImplementedError` on REAL for exactly
this reason, and `similarity` documents the principle: "an empty one says every
isolate is identical to every other, which is a finding rather than an absence."
These two were simply never given the guard.

**Updated: the refusal is now conditional, and the error type changed.**
`46f5005` refused unconditionally with `NotImplementedError`, which was correct
while neither stage had a REAL caller but left both unreachable once their
inputs became real. They now read their inputs from disk and refuse with
`StageError` naming the specific input that is absent - see
`test_convergence_real_inputs.py`, which is the substantive test of that
behaviour. What this file still guards is the property the original ticket
cared about: **TEST still computes, and neither stage has lost its mode
branch.** The `NotImplementedError` assertions below became `StageError`
assertions for the same reason the production code changed, not because the
refusal stopped happening - it did not.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.errors import StageError
from papipeline.models import RunMode
from papipeline.stages import convergence as stage_convergence
from papipeline.stages import cooccurrence as stage_cooccurrence


def _empty_inputs():
    return {
        "amr_calls": {"S1": ["blaOXA-1"]},
        "regulator_variants": {"S1": []},
        "amr_genes": {"S1": ["blaOXA-1"]},
        "variants": {"S1": []},
        "mechanisms": {"S1": []},
        "lineages": {"S1": "LINEAGE_A"},
    }


class TestConvergenceRefusesReal:
    def test_it_raises_on_real(self, config):
        """`StageError`, not `NotImplementedError`.

        The type changed with the contract: the stage now has a REAL path, so it
        is no longer unimplemented - it is refusing on the state of its inputs,
        which is what `StageError` is for. `NotImplementedError` here would claim
        there is no caller, which is no longer true.
        """
        data = _empty_inputs()
        with pytest.raises(StageError) as excinfo:
            stage_convergence.run(
                config, _manifest(), RunMode.REAL,
                amr_calls=data["amr_calls"],
                regulator_variants=data["regulator_variants"],
                lineages=data["lineages"],
                intermediate_root=Path("does/not/exist"),
            )
        assert "convergence" in str(excinfo.value).lower()

    def test_it_names_the_mode_and_says_why(self, config):
        """A refusal a reader can act on: what was asked, and what is missing."""
        data = _empty_inputs()
        with pytest.raises(StageError) as excinfo:
            stage_convergence.run(
                config, _manifest(), RunMode.REAL,
                amr_calls=data["amr_calls"],
                regulator_variants=data["regulator_variants"],
                lineages=data["lineages"],
                intermediate_root=Path("does/not/exist"),
            )
        message = str(excinfo.value)
        assert "REAL" in message
        # "no caller" is no longer the reason; naming the input is.
        assert "amr_calls" in message or "regulator_variants" in message

    def test_test_mode_is_untouched(self, config):
        """The guard is a branch, not a replacement - TEST still computes."""
        data = _empty_inputs()
        result = stage_convergence.run(
            config, _manifest(), RunMode.TEST,
            amr_calls=data["amr_calls"],
            regulator_variants=data["regulator_variants"],
            lineages=data["lineages"],
        )
        assert isinstance(result, list)


class TestCooccurrenceRefusesReal:
    def test_it_raises_on_real(self, config):
        data = _empty_inputs()
        with pytest.raises(StageError) as excinfo:
            stage_cooccurrence.run(
                config, _manifest(), RunMode.REAL,
                amr_genes=data["amr_genes"], variants=data["variants"],
                mechanisms=data["mechanisms"], lineages=data["lineages"],
                intermediate_root=Path("does/not/exist"),
            )
        assert "cooccurrence" in str(excinfo.value).lower()

    def test_test_mode_is_untouched(self, config):
        data = _empty_inputs()
        result = stage_cooccurrence.run(
            config, _manifest(), RunMode.TEST,
            amr_genes=data["amr_genes"], variants=data["variants"],
            mechanisms=data["mechanisms"], lineages=data["lineages"],
        )
        assert isinstance(result, list)


class TestTheRefusalIsNotJustAnErrorMessage:
    def test_neither_stage_branches_silently_in_real_mode(self):
        """A tripwire against the guard being deleted.

        Checks the guard *form* each file now uses - `if mode is not
        RunMode.TEST`, matching `stages.similarity` - rather than a particular
        token. A stage can still compute identically in both modes and merely
        mention RunMode in a comment, so the behavioural tests above carry the
        weight; this one only notices the branch disappearing.
        """
        root = Path(__file__).resolve().parents[2] / "papipeline" / "stages"
        for name in ("convergence.py", "cooccurrence.py"):
            source = (root / name).read_text(encoding="utf-8")
            assert "if mode is not RunMode.TEST:" in source, (
                f"{name} no longer guards on mode, so a REAL run would compute "
                "silently over partial intermediates again"
            )


def _manifest():
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    return SampleManifest(samples=[Sample(sample_id="S1")])
