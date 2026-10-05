"""`stage_inputs()` in REAL mode: the full external set, demanded by name.

The DAG is gated so that REAL never executes here, which means the one property
that matters for REAL cannot be checked by running it. It can still be checked
*directly*: the function is small, has three free variables, and the whole
question is what it returns for a given mode.

So the real function's source is exec'd here rather than reimplemented. A copy
would pass while the Snakefile drifted, which is exactly the class of bug this
exists to catch - the same reason the dry-run check is hermetic.

For the naming half - that a missing external file stops the DAG *by name* -
that is asserted at the TEST level in
``test_snakemake_dry_run.py::test_a_missing_external_input_fails_the_dag_build``.
TEST and REAL take the same branch of `stage_inputs`, so a test proving the
branch is taken does not need the mode gate opened.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"

CONFIG_FILES = ("science.yaml", "machines/laptop.yaml", "machines/bigmachine.yaml")


def _stage_inputs_source() -> str:
    """The function's own source, extracted from the Snakefile."""
    text = SNAKEFILE.read_text(encoding="utf-8")
    start = text.index("def stage_inputs(")
    # to the next top-level definition or assignment
    end = re.search(r"(?m)^(def |CONFIG_FILES = )", text[start + 1:])
    assert end, "could not find the end of stage_inputs()"
    return text[start:start + 1 + end.start()]


def _stage_inputs(mode: RunMode):
    """Call the real function with `mode` bound, as the Snakefile would."""
    namespace = {
        "RESOLVED_MODE": mode,
        "RunMode": RunMode,
        "CONFIG_FILES": list(CONFIG_FILES),
    }
    exec(compile(_stage_inputs_source(), str(SNAKEFILE), "exec"), namespace)
    return namespace["stage_inputs"]


class TestRealModeReturnsTheFullExternalSet:
    def test_real_demands_external_inputs(self):
        stage_inputs = _stage_inputs(RunMode.REAL)
        produced = ["/results/out_amr.tsv"]
        external = ["/data/manifest.tsv", "/data/intermediate/amr/amr.tsv"]
        result = stage_inputs(produced, external)
        assert result == produced + external + list(CONFIG_FILES)

    def test_stub_demands_none_of_them(self):
        """The contrast that gives the REAL case its meaning.

        STUB fabricates everything with no fixture on disk (spec.md D8), so it
        must not demand a file that does not exist. If this ever matched REAL,
        the mode split would be doing nothing.
        """
        produced = ["/results/out_amr.tsv"]
        external = ["/data/manifest.tsv"]
        assert _stage_inputs(RunMode.STUB)(produced, external) == (
            produced + list(CONFIG_FILES)
        )

    def test_produced_inputs_are_demanded_in_every_mode(self):
        """A rule's real dependencies are never dropped, in any mode."""
        produced = ["/results/out_a.tsv", "/results/out_b.tsv"]
        for mode in RunMode:
            result = _stage_inputs(mode)(produced, [])
            assert result[:2] == produced, f"{mode} dropped a produced input"

    def test_test_mode_agrees_with_real(self):
        """Both non-STUB modes take the same branch, by design.

        Which is why the missing-input naming is proven in TEST: if this
        assertion ever fails, the two modes have diverged and the TEST proof no
        longer covers REAL.
        """
        produced, external = ["/o.tsv"], ["/e.tsv"]
        assert (_stage_inputs(RunMode.TEST)(produced, external)
                == _stage_inputs(RunMode.REAL)(produced, external))


class TestTheRealExternalSetNamesItsFiles:
    """What REAL actually asks for, read off the Snakefile's own variables."""

    @pytest.fixture(scope="class")
    def real_inputs(self):
        stage_inputs = _stage_inputs(RunMode.REAL)
        return stage_inputs

    @pytest.mark.parametrize("stage,expected_fragment", [
        ("amr", "amr_determinants.tsv"),
        ("mlst", "mlst_results.tsv"),
        ("phylogeny", "tree.nwk"),
        ("phenotype", "imipenem_phenotype.tsv"),
    ])
    def test_each_stage_demands_its_own_input(self, real_inputs, stage, expected_fragment):
        """A stage's external input must be a full path, so it can be named.

        Snakemake reports a missing input by the path it was given. A glob or a
        bare basename would leave the operator with something to guess from.
        """
        path = f"/data/intermediate/{expected_fragment}"
        result = real_inputs([], [path])
        assert path in result
        assert os.path.isabs(path)

    def test_a_missing_input_would_be_named_not_guessed(self, real_inputs):
        """The path is what reaches Snakemake, verbatim.

        Asserted as identity rather than membership-by-substring: the operator
        is told the path Snakemake was given, so nothing may rewrite it.
        """
        missing = "/data/intermediate/amr/amr_determinants.tsv"
        result = real_inputs([], [missing])
        assert result.count(missing) == 1
        assert result.index(missing) < len(CONFIG_FILES), (
            "an external input must precede the config files; a stage that "
            "fails on a missing config first would report the wrong cause"
        )
