"""`regulators` reads a stage INPUT and writes a stage OUTPUT. Both. Neither is drift.

**A prior agent called this a path drift**: `run.py` writes `06_regulators.tsv`
while `stages/regulators.run` reads `<root>/regulators/regulator_variants.tsv`,
"so stage 6 emits an EMPTY table in REAL". **That claim is FALSE, on both halves,
and the evidence is below.** This file is the standing refutation, because the
"fix" that claim implies - handing `regulators.run` the run's own output root -
is a real regression that a future reader would otherwise apply.

**What `docs/data_contract.md` says.**

* `:114` lists `regulators/regulator_variants.tsv` under the heading
  "### Stage input tables (under `intermediate/`)", and `:97-99` says of that
  whole section: "In TEST mode these are checked-in synthetic fixtures. **In REAL
  mode the external tools write them into the run's intermediate directory.**"
* `:177` lists `06_regulators.tsv` under "## Output files", grain "one row per
  variant call".

Two different files, two different sections, two different grains. The reader
is reading an input; the writer is writing an output. Calling that a drift calls
the contract wrong.

**Why "emits an EMPTY table" is also false.**
`stages/regulators.run` -> `load_regulator_variants` -> `io.tsv.read_tsv`, and
`read_tsv` raises `DataContractError` on a missing file (`io/tsv.py:70`, "Raises:
DataContractError: File missing, empty, or violates the contract"). So in REAL
with no provisioned input the stage **refuses**, loudly, naming the path. It does
not write an empty table. An empty table is the worse outcome and it is not the
one that happens.

**The real gap, named honestly.** No code in this repository writes
`<intermediate>/regulators/regulator_variants.tsv` in REAL. That is a
provisioning gap - the external screen the contract calls for has no caller -
and it is the *same shape* `stages/sv.py` closed with `_run_real`
(`stages/sv.py:200-235`: REAL calls nucmer and writes
`structural_variants/structural_variants.tsv`; TEST reads the committed
fixture). `regulators.run` has no such branch. That is a missing REAL caller, not
a wrong path, and it is a scientific feature rather than a wiring fix - out of
scope for `run.py`.

**The regression this file guards.** `regulators.run`'s fourth parameter is
*named* `intermediate_root` but is correctly handed `tool_output_root`. Reading
the parameter name as the contract is the mistake; in TEST the two roots are
different paths, so the "fix" breaks every TEST run with a missing-input error.
`tests/unit/test_stage_root_invariant.py:36-44` records the same warning in
prose. This asserts it in code, on **distinct temporary roots**, so the two
candidate paths cannot coincide and a swapped-root bug is invisible.
"""

from __future__ import annotations

import inspect
import shutil
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import PipelineError
from papipeline.run import run_pipeline

REPO = Path(__file__).resolve().parents[2]

#: The two candidate paths the prior claim conflated. Asserted distinct before
#: anything else, because if they were equal this file could not tell a correct
#: wiring from a swapped one.
READER_SUFFIX = Path("regulators") / "regulator_variants.tsv"
WRITER_NAME = "06_regulators.tsv"


@pytest.fixture()
def project(tmp_path, monkeypatch):
    """A pipeline whose every root is a distinct path under ``tmp_path``.

    ``config/`` is copied so the pipeline root moves and every accessor becomes a
    temporary path; ``test_data`` is symlinked because ``tool_output_root(TEST)``
    has no redirect hook (see ``tests/unit/test_stage_root_invariant.py``, which
    establishes this shape and why it is necessary). ``PIPELINE_RESULTS_ROOT``
    moves the run's own output tree, which is what makes the two roots distinct.
    """
    root = tmp_path / "project"
    root.mkdir()
    shutil.copytree(REPO / "config", root / "config")
    (root / "test_data").symlink_to(REPO / "test_data", target_is_directory=True)
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "written"))
    return load_config(root / "config" / "science.yaml")


@pytest.fixture()
def recorded(project, monkeypatch) -> Tuple[List[Tuple[str, Any]], str]:
    """Run the whole TEST pipeline with a spy on ``stage_regulators.run``.

    A stage failure is captured rather than propagated: a mis-wired root usually
    makes the run *refuse*, which is the correct behaviour and an unhelpful test
    failure, because the reader learns that "variants" broke and not which
    argument was wrong.
    """
    import papipeline.run as run_module

    seen: List[Tuple[str, Any]] = []
    original = run_module.stage_regulators.run
    signature = inspect.signature(original)

    def spy(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        for name, value in bound.arguments.items():
            if isinstance(value, Path):
                seen.append((name, value))
        return original(*args, **kwargs)

    monkeypatch.setattr(run_module.stage_regulators, "run", spy, raising=True)
    failure = ""
    try:
        run_pipeline(config=project, mode="TEST", write_html=False)
    except PipelineError as exc:
        failure = f"{type(exc).__name__}: {exc}"
    return seen, failure


class TestTheTwoPathsAreDifferentFiles:
    """The premise. If these coincided, nothing below could fail."""

    def test_the_reader_path_is_not_the_writers_path(self, project):
        from papipeline.models import RunMode

        reader = project.tool_output_root(RunMode.TEST) / READER_SUFFIX
        writer = project.intermediate_root(RunMode.TEST) / "stages" / WRITER_NAME
        assert reader != writer
        assert reader.is_file(), (
            f"the committed stage-input fixture is absent at {reader}; without it "
            "this file cannot distinguish the two paths and the whole file is void"
        )

    def test_the_two_roots_are_distinct_in_test(self, project):
        from papipeline.models import RunMode

        assert project.intermediate_root(RunMode.TEST) != project.tool_output_root(
            RunMode.TEST
        ), (
            "intermediate_root and tool_output_root resolved to the same path, so "
            "a swapped wiring would be invisible below"
        )


class TestTheReaderIsHandedTheInputRoot:
    def test_regulators_reads_the_stage_input_root_not_the_output_root(
        self, project, recorded,
    ):
        """The assertion the prior agent's "fix" would break.

        ``regulators.run``'s parameter is named ``intermediate_root``. Handing it
        ``intermediate`` instead looks like a fix and is a regression: in TEST
        the input lives under ``test_data/intermediate`` and nothing writes
        ``regulator_variants.tsv`` into the results tree, so the run would stop
        with a missing-input error on every TEST run.
        """
        from papipeline.models import RunMode

        seen, failure = recorded
        assert not failure, f"the run stopped first: {failure}"
        assert seen, "stage_regulators.run took no path argument and was not checked"
        name, handed = seen[0]
        assert name == "intermediate_root"
        expected = project.tool_output_root(RunMode.TEST)
        assert handed == expected, (
            f"regulators.run was handed {handed}, but the contract "
            "(docs/data_contract.md:114) says its input is read from "
            f"tool_output_root = {expected}"
        )
        assert handed != project.intermediate_root(RunMode.TEST), (
            "regulators.run was handed this run's own output directory. Nothing "
            "writes regulator_variants.tsv there, so this is the path drift the "
            "prior agent reported - reached from the other direction."
        )


class TestTheWriterIsHandedTheOutputRoot:
    def test_the_declared_output_lands_in_this_runs_stage_directory(
        self, project, recorded,
    ):
        """`06_regulators.tsv` is an OUTPUT (`docs/data_contract.md:177`).

        Asserted to exist in the results tree and to carry the declared grain -
        one row per variant call, never a header-only stub, which in REAL would
        be indistinguishable from "no regulator variants exist".
        """
        from papipeline.models import RunMode

        seen, failure = recorded
        assert not failure, f"the run stopped first: {failure}"
        output = (
            project.intermediate_root(RunMode.TEST) / "stages" / WRITER_NAME
        )
        assert output.is_file(), (
            f"{output.name} was not written; docs/data_contract.md:177 declares it"
        )
        lines = output.read_text(encoding="utf-8").splitlines()
        data = [line for line in lines if line and not line.startswith("#")]
        assert len(data) >= 2, (
            f"{output.name} has a header and no rows. In REAL that reads as "
            "'no regulator variants were called', which is a finding this stage "
            "has not made."
        )


class TestAMissingInputRefusesRatherThanEmittingAnEmptyTable:
    """The second half of the prior claim, which is also false.

    The claim was that REAL "emits an EMPTY table". It refuses. Asserted against
    a directory with no provisioned regulator table, which is what a REAL run
    without the external screen looks like.
    """

    def test_a_missing_input_raises_naming_the_path(self, tmp_path):
        from papipeline.models import RunMode
        from papipeline.stages import regulators

        empty = tmp_path / "intermediate"
        (empty / "regulators").mkdir(parents=True)
        with pytest.raises(Exception) as excinfo:
            regulators.run(
                _stub_config(), _stub_manifest(), RunMode.REAL, empty,
            )
        message = str(excinfo.value)
        assert "regulator_variants.tsv" in message, (
            "a missing stage input must name the file. A refusal that only said "
            "'input missing' leaves the reader guessing between a wrong root and "
            "an unrun tool."
        )
        assert str(READER_SUFFIX.parent) in message, (
            "the refusal must name the directory it looked in, or a wrong "
            "tool_output_root is indistinguishable from an unprovisioned cohort"
        )


def _stub_config():
    class _C:
        raw: Dict[str, Any] = {}

    return _C()


def _stub_manifest():
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    return SampleManifest(
        samples=(Sample(sample_id="S1", assembly_path=None, source="test"),)
    )