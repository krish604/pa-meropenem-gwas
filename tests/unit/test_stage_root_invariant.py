"""Every root-consuming stage is handed the root its body actually reads.

**Why this file exists.** The pipeline has two roots that are easy to confuse and
that are *the same path in REAL and different paths in TEST*:

    intermediate_root(mode) = results_root(mode) / "intermediate"   # this run WRITES here
    tool_output_root(mode)  = data_root(mode)  / "intermediate"    # stage INPUTS live here
                              ... except in REAL, where it is intermediate_root

`loader.py:537-540` says why the accessor exists: "Keeping this distinct from
`intermediate_root` is what stops a run from reading its own outputs as if they
were inputs, and stops synthetic fixtures from being written into `results/`."

Stage 7 was handed the wrong one once. `tests/unit/test_pangenome_caller_root.py`
records it: it read its own output directory, found no annotations, and reported
`0 genes | 0 core | 0 accessory` with a success exit over a fixture declaring
10/5/5. Nothing refused. That test covers **one** stage, and only the two ways
that one stage can be reached.

**What this adds.** A whole-pipeline invariant over every stage entry point that
takes a root, in one run, with the roots made distinct first. Two halves:

* every path handed to a stage is inside the root that stage is contracted to
  read or write - nothing else, in particular not the other root;
* the run leaves the fixture tree byte-identical, which is the "writes ONLY its
  contracted root" half, checked on disk rather than inferred.

**Why the fixture root is a symlink and not a copy.** `tool_output_root(TEST)` is
`data_root(TEST) / "intermediate"` and `data_root(TEST)` is
`<config root>/test_data`, with no override hook short of rewriting the loader -
and `papipeline/config/loader.py` is off-limits. So `config/` is copied into
`tmp_path` (making the pipeline root temporary and every accessor a temporary
path) and `test_data` is symlinked. Path identity, which is what a mis-wiring is
about, is still distinct; the content is shared read-only, and the
no-writes-into-fixtures assertion is what proves it stayed read-only.

A note recorded because it will bite the next person: **the regulators call is
correct.** `regulators.run`'s fourth parameter is *named* `intermediate_root`
and is given `tool_output_root`; in REAL that is the same path, and in TEST the
fixtures are under `test_data/intermediate`. Reading the parameter name as the
contract produces a "fix" that breaks TEST, because nothing writes
`regulator_variants.tsv` into the results tree. Eight stage signatures name that
parameter `intermediate_root` for exactly this reason; the table below is the
authority, not the parameter name.
"""

from __future__ import annotations

import inspect
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import PipelineError
from papipeline.models import RunMode
from papipeline.run import run_pipeline

REPO = Path(__file__).resolve().parents[2]

#: (dotted attribute on `papipeline.run`, stage label) for every entry point that
#: is handed a root. A stage missing from this list is a stage whose wiring is
#: unchecked, so adding one is the way to extend the invariant.
ROOTED_ENTRY_POINTS: Tuple[Tuple[str, str], ...] = (
    ("stage_annotation.run", "annotation"),
    ("stage_mlst.run", "mlst"),
    ("stage_amr.run", "amr"),
    ("stage_sv.run", "structural_variants"),
    ("stage_variants.run", "variants"),
    ("stage_regulators.run", "regulators"),
    ("stage_virulence.run", "virulence"),
    ("stage_pangenome.run", "pangenome"),
    ("stage_phylo.run", "phylogeny"),
    ("stage_similarity.run", "similarity"),
    ("stage_phenotype.load_phenotype", "phenotype"),
    ("stage_gwas.run", "gwas"),
    ("stage_convergence.run", "convergence"),
    ("stage_cooccurrence.run", "cooccurrence"),
    ("stage_mechanisms.run", "mechanisms"),
    ("stage_integration.write_master_table", "master_table"),
    ("stage_reporting.write_report", "reporting"),
)

#: The contracted root for each parameter of each entry point, as an attribute
#: name on `PipelineConfig` (or one of the derived keys below). Read from the
#: accessors, never hard-coded, so a change to `loader.py` moves the expectation
#: with it instead of leaving a stale literal behind.
CONTRACT: Mapping[str, Mapping[str, str]] = {
    # These read stage INPUT TABLES and, in REAL, write them.
    "annotation": {"intermediate_root": "tool_output"},
    "mlst": {"intermediate_root": "tool_output", "data_root": "assembly"},
    "amr": {"intermediate_root": "tool_output", "data_root": "assembly"},
    "structural_variants": {"intermediate_root": "tool_output"},
    "regulators": {"intermediate_root": "tool_output"},
    "virulence": {"intermediate_root": "tool_output"},
    "pangenome": {"intermediate_root": "tool_output"},
    "gwas": {"intermediate_root": "tool_output"},
    # `data_root` is the assembly root; `workdir` is scratch inside this run's
    # own output tree, which is the one place that is correct.
    "variants": {"data_root": "assembly", "workdir": "intermediate/variants"},
    # A dedicated directory, not a stage-input root.
    "phylogeny": {"phylogeny_dir": "phylogeny"},
    "similarity": {"tree_path": "phylogeny/tree.nwk", "out_path": "stage_dir/*"},
    "phenotype": {"phenotype_dir": "phenotype"},
    # These re-read THIS RUN'S OWN tables, so the output root is correct here.
    #
    # `phylogeny_dir` on 13 and 14 is new in round 12 Phase 3 and is the `phylogeny`
    # root, not a stage-input root: `convergence.real_input_paths` and
    # `cooccurrence.real_input_paths` both spell the crosswalk
    # `<phylogeny_dir>/tree_metadata.tsv`, because stage 9 writes outside the
    # intermediate root. Recording it here is what makes "the reader and the
    # writer agree" checkable rather than a comment - and it is the same contract
    # stage 9 itself is held to two rows above, so a change that moved one
    # without the other fails.
    "convergence": {
        "intermediate_root": "intermediate", "phylogeny_dir": "phylogeny",
    },
    "cooccurrence": {
        "intermediate_root": "intermediate", "phylogeny_dir": "phylogeny",
    },
    # Stage 5 (mechanisms), folded into stage 14. In REAL it re-reads its four
    # inputs from its own intermediate root and ignores the mappings run.py hands
    # it, so the root is load-bearing there in a way it is not for 13/14 above.
    "mechanisms": {"intermediate_root": "intermediate"},
    "master_table": {"path": "stage_dir/*"},
    "reporting": {"out_dir": "reports"},
}

#: Contract values ending in `/*` mean "somewhere strictly inside this
#: directory"; a bare value means "exactly this path". Both appear, because a
#: scratch directory is named outright while a table is named by its parent.
CONTAINMENT = "/*"


def _expected(key: str, roots: "_Roots") -> Path:
    directory, _, tail = key.partition("/")
    path = _realpath(roots.__dict__[directory])
    for part in tail.split("/"):
        if part and part != "*":
            path = path / part
    return path


def _realpath(path: Any) -> Any:
    """Normalise for comparison: `/var` and `/private/var` are one directory."""
    if isinstance(path, Path):
        return Path(os.path.realpath(path))
    return path


class _Roots:
    """Every root this pipeline knows about, for one resolved configuration."""

    def __init__(self, config, mode: RunMode) -> None:
        self.tool_output = config.tool_output_root(mode)
        self.assembly = config.assembly_root(mode)
        self.intermediate = config.intermediate_root(mode)
        self.stage_dir = self.intermediate / "stages"
        self.phylogeny = config.phylogeny_dir(mode)
        self.phenotype = config.phenotype_dir(mode)
        self.reports = config.reports_root(mode)

    def get(self, key: str) -> Path:
        return _expected(key, self)

    @property
    def named(self) -> Dict[str, Path]:
        return {
            name: _realpath(value)
            for name, value in (
                ("tool_output", self.tool_output),
                ("assembly", self.assembly),
                ("intermediate", self.intermediate),
                ("stage_dir", self.stage_dir),
                ("phylogeny", self.phylogeny),
                ("phenotype", self.phenotype),
                ("reports", self.reports),
            )
        }


@pytest.fixture()
def project(tmp_path, monkeypatch):
    """A pipeline whose every root is a path under ``tmp_path``.

    ``config/`` is copied so the pipeline root moves; ``test_data`` is symlinked
    because `tool_output_root(TEST)` has no redirect hook (see module docstring).
    `PIPELINE_RESULTS_ROOT` moves the results, intermediate and reports roots.
    """
    root = tmp_path / "project"
    root.mkdir()
    shutil.copytree(REPO / "config", root / "config")
    (root / "test_data").symlink_to(REPO / "test_data", target_is_directory=True)
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "written"))
    config = load_config(root / "config" / "science.yaml")
    return config, _Roots(config, RunMode.TEST)


@pytest.fixture()
def fixture_tree() -> Dict[str, Tuple[int, int]]:
    """Size and mtime of every committed fixture file, for the no-writes check."""
    return _snapshot(REPO / "test_data")


def _snapshot(tree: Path) -> Dict[str, Tuple[int, int]]:
    out: Dict[str, Tuple[int, int]] = {}
    for path in sorted(tree.rglob("*")):
        if path.is_file():
            stat = path.stat()
            out[str(path)] = (stat.st_size, stat.st_mtime_ns)
    return out


def _record(
    monkeypatch, recorded: List[Tuple[str, str, Any]]
) -> None:
    """Wrap every rooted entry point so it records its path arguments."""
    import papipeline.run as run_module

    for dotted, label in ROOTED_ENTRY_POINTS:
        target = run_module
        for part in dotted.split(".")[:-1]:
            target = getattr(target, part)
        name = dotted.split(".")[-1]
        original = getattr(target, name)
        signature = inspect.signature(original)

        def make(original=original, signature=signature, label=label):
            def spy(*args, **kwargs):
                bound = signature.bind(*args, **kwargs)
                for param, value in bound.arguments.items():
                    if isinstance(value, Path):
                        recorded.append((label, param, value))
                return original(*args, **kwargs)

            return spy

        monkeypatch.setattr(target, name, make(), raising=True)


def _run_and_record(project, monkeypatch) -> Tuple[List[Tuple[str, str, Any]], _Roots, str]:
    """Run the whole TEST pipeline with spies installed.

    A stage failure is captured rather than propagated. A mis-wired root usually
    makes the run *refuse* - which is the correct behaviour and a useless test
    failure, because the reader learns that "variants" broke and not which
    argument was wrong. Capturing it lets the assertions below name the argument,
    and the captured text is reported alongside so the original refusal is not
    lost either.
    """
    config, roots = project
    recorded: List[Tuple[str, str, Any]] = []
    _record(monkeypatch, recorded)
    try:
        run_pipeline(config=config, mode="TEST", write_html=False)
    except PipelineError as exc:
        return recorded, roots, f"{type(exc).__name__}: {exc}"
    return recorded, roots, ""


class TestTheRootsAreDistinct:
    """The premise. If two roots were the same path, nothing below could fail."""

    def test_the_run_output_root_is_not_the_stage_input_root(self, project):
        _, roots = project
        assert roots.intermediate != roots.tool_output, (
            "intermediate_root and tool_output_root resolved to the same path, so "
            "this test could not distinguish a correct wiring from a swapped one"
        )

    def test_every_root_is_a_distinct_path(self, project):
        _, roots = project
        named = roots.named
        for left, left_path in named.items():
            for right, right_path in named.items():
                if left >= right:
                    continue
                assert left_path != right_path, f"{left} and {right} share {left_path}"

    def test_the_written_roots_are_disposable_and_the_input_roots_are_not(
        self, project, tmp_path
    ):
        """The premise, stated as the invariant that actually matters.

        `tool_output_root(TEST)` has no redirect hook, so it necessarily resolves
        through the symlink back into the repository (see the module docstring).
        That is acceptable only because the run never writes there - which is the
        assertion in `TestTheRunWritesOnlyItsOwnTree`. What must hold here is the
        separation: the roots this run writes into are temporary, and none of
        them contains a root it reads from.
        """
        _, roots = project
        real_tmp = _realpath(tmp_path)
        for name in ("intermediate", "stage_dir", "reports"):
            assert _realpath(roots.named[name]).is_relative_to(real_tmp), (
                f"{name} is not under {tmp_path}, so the run would have written "
                "into the repository"
            )
        writable = [roots.named[n] for n in ("intermediate", "stage_dir", "reports")]
        for name in ("tool_output", "assembly", "phylogeny", "phenotype"):
            read_root = roots.named[name]
            for out in writable:
                assert not _within(read_root, out), (
                    f"{name} ({read_root}) sits inside the run's own output tree "
                    f"({out}), so a stage could read its own outputs as inputs"
                )


class TestEachStageIsHandedItsContractedRoot:
    """The invariant proper: one row per root-consuming entry point."""

    def test_every_rooted_entry_point_was_reached(self, project, monkeypatch):
        recorded, _, failure = _run_and_record(project, monkeypatch)
        seen = {label for label, _, _ in recorded}
        missing = sorted({label for _, label in ROOTED_ENTRY_POINTS} - seen)
        assert not missing and not failure, (
            f"these entry points took no path and so were not checked: {missing}. "
            "Either the stage no longer runs in TEST, or its root argument was "
            "removed - update ROOTED_ENTRY_POINTS and CONTRACT deliberately."
            + (f"\nThe run stopped first: {failure}" if failure else "")
        )

    def test_each_argument_equals_its_contracted_root(
        self, project, monkeypatch
    ):
        recorded, roots, failure = _run_and_record(project, monkeypatch)
        problems: List[str] = []
        for label, param, value in recorded:
            expected_key = CONTRACT.get(label, {}).get(param)
            if expected_key is None:
                problems.append(
                    f"{label}.{param}={value} is not in CONTRACT - the invariant "
                    "does not describe this argument"
                )
                continue
            expected = _expected(expected_key, roots)
            if expected_key.endswith(CONTAINMENT):
                ok = _realpath(value) != expected and _within(_realpath(value), expected)
            else:
                ok = _realpath(value) == expected
            if not ok:
                problems.append(
                    f"{label}.{param} was handed {value}, but its contract says "
                    f"{expected_key} = {expected}"
                )
        assert not problems and not failure, "\n".join(
            problems + ([failure] if failure else [])
        )

    def test_no_stage_reads_the_other_root(
        self, project, monkeypatch
    ):
        """Named separately, because it is the specific regression.

        `intermediate_root` and `tool_output_root` are one path in REAL, so a
        caller that swapped them would pass every root-equality assertion above
        in REAL. This asserts the mapping rather than the value's identity.
        """
        recorded, _, failure = _run_and_record(project, monkeypatch)
        swapped = [
            f"{label}.{param}={value}"
            for label, param, value in recorded
            if CONTRACT.get(label, {}).get(param) == "tool_output"
            and _realpath(value) == _realpath(
                Path(project[1].intermediate)
            )
        ]
        assert not swapped, (
            "these stages were handed this run's own output directory instead of "
            f"the tool output root: {swapped}"
            + (f"\nThe run stopped first: {failure}" if failure else "")
        )

    def test_the_output_root_stages_are_not_handed_the_fixture_root(
        self, project, monkeypatch
    ):
        """The converse trap: `convergence` re-reads this run's own tables."""
        recorded, roots, failure = _run_and_record(project, monkeypatch)
        wrong = [
            f"{label}.{param}={value}"
            for label, param, value in recorded
            if CONTRACT.get(label, {}).get(param) == "intermediate"
            and _realpath(value) == _realpath(roots.tool_output)
        ]
        assert not wrong, (
            "these stages re-read this run's tables and must be handed the "
            f"output root, not the fixture root: {wrong}"
            + (f"\nThe run stopped first: {failure}" if failure else "")
        )

    def test_every_path_handed_to_a_stage_is_inside_a_known_root(
        self, project, monkeypatch
    ):
        """Nothing reaches a stage from outside the roots the config declares."""
        recorded, roots, failure = _run_and_record(project, monkeypatch)
        containers = list(roots.named.values())
        strays = []
        for label, param, value in recorded:
            real = _realpath(value)
            if not any(_within(real, c) for c in containers):
                strays.append(f"{label}.{param}={value}")
        assert not strays and not failure, (
            f"these paths are outside every configured root: {strays}"
            + (f"\nThe run stopped first: {failure}" if failure else "")
        )


class TestTheRunWritesOnlyItsOwnTree:
    """The write half of the invariant, checked on disk."""

    def test_the_committed_fixtures_are_untouched(
        self, project, monkeypatch, fixture_tree
    ):
        _, _, failure = _run_and_record(project, monkeypatch)
        after = _snapshot(REPO / "test_data")
        changed = sorted(
            name
            for name in set(fixture_tree) | set(after)
            if fixture_tree.get(name) != after.get(name)
        )
        assert not changed and not failure, (
            "a TEST run modified committed fixtures, so a stage wrote to its "
            f"input root: {changed}"
            + (f"\nThe run stopped first: {failure}" if failure else "")
        )

    def test_the_four_written_roots_are_all_under_the_results_redirect(
        self, project, monkeypatch, tmp_path
    ):
        _, roots = project
        _, _, failure = _run_and_record(project, monkeypatch)
        real_tmp = _realpath(tmp_path)
        for name in ("intermediate", "stage_dir", "reports"):
            path = _realpath(roots.named[name])
            assert path.is_relative_to(real_tmp)
            assert path.exists() and not failure, (
                f"{name} was never created: {path}"
                + (f"\nThe run stopped first: {failure}" if failure else "")
            )


def _within(path: Path, container: Path) -> bool:
    return path == container or container in path.parents
