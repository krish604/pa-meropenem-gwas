"""Every stage call site in `run.py` supplies every parameter the stage requires.

**What this is.** `papipeline/run.py` is the only place the pipeline's stage
entry points are invoked, and it does so by attribute call:
``stage_amr.run(config, manifest, mode, intermediate_root, ...)``. A stage that
grows a new *required* parameter - a root, a path, an id set - and a `run.py`
that is not updated to pass it produce two failures that look like one:

* Python raises `TypeError: run() missing 1 required positional argument` at
  the moment the stage is reached, which is **after** the run has already spent
  the cost of every earlier stage; and
* worse, when the new parameter is added **with a default of `None`**, nothing
  raises at all. The stage receives `None`, and whether that is a refusal or a
  silent wrong answer is decided inside the stage body, per parameter, forever.

So the first failure is loud and late, and the second is quiet and permanent.
Both are the same defect - a required input that the call site does not supply -
and this test is the check that catches it in a second rather than in a run.

**Why `inspect.Signature.bind` and not running the pipeline.** `bind` is a pure
function of the signature and the supplied argument names. It answers exactly
the question - "is every parameter this stage declares as required actually
supplied *here*?" - with no cohort, no tools, no mode gate, and no 400 fixtures.
Running a TEST pipeline to find out would also be *unable* to find out: the paths
that matter are the REAL-only ones, and a TEST run takes the other branch.

The arguments are replaced with sentinels because `bind` checks arity and names,
never types or values. A sentinel per positional argument and per keyword is
therefore a faithful stand-in for the real call, and substituting them keeps the
test from accidentally depending on what `run.py` happens to pass.

**The two halves, and why "no default" is the whole of the declaration.** There
is no `required_outside_test=` decorator in this codebase and there is not going
to be one: the codebase expresses the requirement in the signature itself. A
parameter with no default is unsuppliable-optional, so `bind` fails if the call
site omits it - that is the whole mechanism, and it is checked here against the
live signature. The other form is a parameter *with* a default, which `bind`
cannot help with: the call site omitting it is legal Python. For those,
`PartThree` asserts that the stage body either refuses on `None` or falls back to
a named config accessor - and that the set which does which is written down, so
adding a fourth defaulted root is a deliberate edit rather than a silent
widening of what may be `None` in a real run.

`PartFour` widens that to every defaulted parameter, not only the roots: a
parameter with an `Optional[...] = None` default that is not a root is supplied
by nobody and refused by nobody, and the stage body decides alone what its
absence means. Those omissions are enumerated and must each carry a reason in
`tests/unit/seam_call_site_defaulted_allowlist.txt`. An omission recorded there
may be a deliberate fallback, a TEST-only path, or a genuine defect - the file
distinguishes them, and round 12 found three real ones.

**Anti-vacuity.** A test that finds no call sites passes. Three guards prevent
that: `PartOne::test_every_stage_module_reached_by_run_py_is_covered` fails if
`run.py` imports a stage this file does not account for; the recorded
`REQUIRED_PARAMETERS` table is compared against the live signatures, so a stage
that stops requiring something cannot pass unnoticed; and `PartOne` asserts the
table is non-empty for every stage that has any required parameter at all.
"""

from __future__ import annotations

import ast
import importlib
import inspect
import re
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

import pytest

REPO = Path(__file__).resolve().parents[2]
RUN_PY = REPO / "papipeline" / "run.py"


class _Sentinel:
    """Stands in for a real argument. `bind` never inspects the value."""

    def __init__(self, label: str) -> None:
        self.label = label

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<sentinel {self.label}>"


def _import_map(tree: ast.Module) -> Dict[str, str]:
    """Local alias -> dotted module, over module-level *and* function-local imports.

    `run.py` imports most stages at module scope but imports
    ``stages.recombination`` *inside* `derive_recombination_tables`, so a
    module-scope-only map would silently drop stage 8 from this file's
    coverage - and stage 8 is the one stage whose call site passes four
    keyword arguments and no positional ones after `mode`.
    """
    aliases: Dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            # `run.py` lives in the `papipeline` package, so its own relative
            # imports resolve against `papipeline`. Without the prefix below,
            # `from .stages import amr as stage_amr` would import `stages`,
            # which is not a package - and the failure would look like a
            # packaging break rather than a bug in this file.
            prefix = ".".join(["papipeline"] * node.level)
            base = f"{prefix}.{node.module}" if prefix else node.module
            for alias in node.names:
                name = alias.asname or alias.name
                aliases[name] = f"{base}.{alias.name}"
        elif isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.asname or alias.name.split(".")[0]
                aliases[name] = alias.name
    return aliases


def _stage_calls(tree: ast.Module) -> list:
    """Every ``<stage_alias>.run(...)`` call, in source order.

    An attribute call on a name that looks like a stage alias. The alias set is
    taken from the import map, not from a name pattern, so `stage_` is a
    consequence of the imports rather than the definition of the scope - a stage
    imported as `from .stages import amr` (no alias) is still found.
    """
    aliases = _import_map(tree)
    calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not isinstance(func, ast.Attribute) or func.attr != "run":
            continue
        base = func.value
        if not isinstance(base, ast.Name) or base.id not in aliases:
            continue
        calls.append(
            {
                "line": node.lineno,
                "alias": base.id,
                "module": aliases[base.id],
                "n_positional": len(node.args),
                "keywords": sorted(
                    kw.arg for kw in node.keywords if kw.arg is not None
                ),
                "starargs": any(kw.arg is None for kw in node.keywords),
                "starargs_positional": node.args
                and any(
                    isinstance(a, ast.Starred) for a in node.args
                ),
                "source": ast.unparse(node),
            }
        )
    return sorted(calls, key=lambda c: c["line"])


def _callables() -> Mapping[str, Any]:
    calls = _stage_calls(ast.parse(RUN_PY.read_text(encoding="utf-8")))
    by_alias: Dict[str, Any] = {}
    for call in calls:
        by_alias.setdefault(call["alias"], call["module"])
    return by_alias


def _signature(module: str) -> inspect.Signature:
    return inspect.signature(importlib.import_module(module).run)


#: Parameters each stage's signature declares as required, with no default.
#:
#: Read from the live signature and compared against this table, so the table is
#: a record of the requirement rather than a substitute for it. Adding a required
#: parameter to a stage fails `PartTwo` until it is written down here, which is
#: the point: the caller in `run.py` is what has to change with it, and this
#: table is where that obligation becomes visible.
REQUIRED_PARAMETERS: Mapping[str, Tuple[str, ...]] = {
    "papipeline.stages.amr": ("config", "manifest", "mode", "intermediate_root", "antibiotic"),
    "papipeline.stages.annotation": ("config", "manifest", "mode", "intermediate_root"),
    "papipeline.stages.cohort_variants": ("config", "manifest", "calls_by_isolate"),
    "papipeline.stages.convergence": ("config", "manifest", "mode"),
    "papipeline.stages.cooccurrence": ("config", "manifest", "mode"),
    "papipeline.stages.gwas": (
        "config", "manifest", "mode", "intermediate_root", "phenotype_calls",
    ),
    "papipeline.stages.integration": (
        "config", "manifest", "mode", "antibiotic", "phenotype", "amr",
        "mechanisms", "regulator_variants", "structural", "mlst", "lineages",
        "virulence", "gwas_results", "convergence",
    ),
    "papipeline.stages.mechanisms": (
        "config", "manifest", "mode", "antibiotic", "amr_calls",
    ),
    "papipeline.stages.mlst": ("config", "manifest", "mode", "intermediate_root"),
    "papipeline.stages.pangenome": ("config", "manifest", "mode", "intermediate_root"),
    "papipeline.stages.phylogeny": ("config", "manifest", "mode", "phylogeny_dir"),
    "papipeline.stages.recombination": ("config", "manifest", "mode"),
    "papipeline.stages.regulators": ("config", "manifest", "mode", "intermediate_root"),
    "papipeline.stages.similarity": ("config", "manifest", "mode"),
    "papipeline.stages.sv": ("config", "manifest", "mode", "intermediate_root"),
    "papipeline.stages.validation": ("config", "manifest", "mode"),
    "papipeline.stages.variants": ("config", "manifest", "mode", "data_root", "workdir"),
    "papipeline.stages.virulence": ("config", "manifest", "mode", "intermediate_root"),
}

#: Stage modules `run.py` imports that it does **not** call `.run` on, and the
#: entry point it reaches instead. Both of these expose their work as a named
#: function rather than a `run`, so `PartOne` would otherwise report them as
#: uncovered stages - which is the failure mode this table exists to make
#: explicit rather than to suppress.
ALTERNATIVE_ENTRY_POINTS: Mapping[str, Tuple[str, ...]] = {
    "papipeline.stages.phenotype": ("load_phenotype",),
    "papipeline.stages.reporting": ("write_report",),
}

#: Parameter names that denote one of the pipeline's roots or the pinned
#: reference. `bind` cannot help with these when they carry a default, so
#: `PartThree` handles them separately.
ROOT_PARAMETERS = frozenset(
    {
        "intermediate_root",
        "tool_output_root",
        "data_root",
        "assembly_root",
        "workdir",
        "phylogeny_dir",
        "phenotype_dir",
        "reports_root",
        "out_dir",
        "db_root",
        "reference_fasta",
        "reference_gff",
    }
)

#: Stage modules whose `.run` `run.py` reaches **only from a TEST branch**.
#:
#: `stages/recombination.run` is one. The REAL path for stage 8 is
#: `papipeline.run.derive_recombination_tables`, which drives
#: `adapters.gubbins.run_gubbins` and `adapters.gubbins.build_outputs` and calls
#: `stages.recombination.block_rows` for the arithmetic - it never calls
#: `stages.recombination.run`. So `run()`'s two defaulted root parameters are
#: not REAL-mode inputs at all: `intermediate_root` is never read by the body,
#: and `phylogeny_dir` is guarded by `if phylogeny_dir is not None:` and the
#: TEST path deliberately omits it so a run cannot replace a committed fixture.
#:
#: That is a fact about the *call site*, not about the signature, so it is
#: checked at the call site: `test_a_test_only_stage_stays_test_only` fails if
#: `run.py` ever calls one of these modules outside a TEST branch. At that point
#: the two parameters above become REAL inputs that are neither required nor
#: refused, which is the defect this file exists to catch.
TEST_ONLY_ENTRY_POINTS = frozenset({"papipeline.stages.recombination"})

#: Root parameters that carry a default and are **refused** when None, rather
#: than defaulted. Read off the live source by `PartThree`, which fails if a
#: parameter is neither in here nor in `CONFIG_DEFAULTED_ROOT_PARAMETERS`.
#:
#: `Optional[Path] = None` is the only reason this file cannot simply read
#: `required_parameters` and stop: `convergence` and `cooccurrence` both take
#: `intermediate_root: Optional[Path] = None` and both refuse in REAL when it is
#: absent (`stages/convergence.py:365`, `stages/cooccurrence.py:450`), through
#: `stages/real_inputs.refuse_incomplete` rather than a bare `if ... is None`.
CONFIG_DEFAULTED_ROOT_PARAMETERS: Mapping[str, Mapping[str, str]] = {
    "papipeline.stages.convergence": {"phylogeny_dir": "phylogeny_dir"},
    "papipeline.stages.cooccurrence": {"phylogeny_dir": "phylogeny_dir"},
    "papipeline.stages.regulators": {"reference_fasta": "reference_fasta"},
}

#: Every defaulted parameter omitted at a call site, and why that is acceptable,
#: read from a file rather than written here.
#:
#: A file, because the list is a record of *reviewed* omissions and it has to be
#: reviewable on its own: a reader should be able to read the reasons without
#: reading this test, and a reviewer should be able to change one reason without
#: touching code. The format and the reasoning behind it are documented in the
#: file's own header, which is the authoritative description - this comment says
#: only what the file is FOR.
#:
#: `PartFour` fails on an omission with no entry AND on an entry for an omission
#: that no longer exists, so the file cannot rot in either direction.
#: Beside this file rather than in a `data/` subdirectory: `.gitignore` line 12 is
#: a bare `data`, deliberately without a trailing slash so it matches a
#: directory of that name at ANY depth (the rule exists to keep the real dataset
#: and its symlink untracked). A `tests/unit/data/` directory would have been
#: silently ignored, and the allowlist is load-bearing - a test that reads a file
#: git never tracks fails only on a clean checkout.
ALLOWLIST_PATH = REPO / "tests" / "unit" / "seam_call_site_defaulted_allowlist.txt"

#: The line separating the allowlist file's prose header from its entries.
ENTRIES_MARKER = "# The entries."


def _allowlist() -> Mapping[str, str]:
    """``"<module>.<parameter>" -> reason``, from the file.

    Everything above the ``ENTRIES_MARKER`` line is prose addressed to a human
    and is skipped; the file is documented in prose *and* in data, and the prose
    must not have to be phrased as an entry to be allowed to explain itself.

    Below the marker, a malformed line is a hard error rather than a silently
    skipped one. Skipping it would be worse than failing: the entry is the
    justification for an omission, so a line that fails to parse is an omission
    that has quietly lost the only record that anyone looked at it.
    """
    text = ALLOWLIST_PATH.read_text(encoding="utf-8")
    assert ENTRIES_MARKER in text, (
        f"{ALLOWLIST_PATH.name} has no {ENTRIES_MARKER!r} line, so this file "
        "cannot tell its prose from its entries. See the header."
    )
    body = text.split(ENTRIES_MARKER, 1)[1]

    entries: Dict[str, str] = {}
    for number, raw in enumerate(body.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "::" not in line:
            raise AssertionError(
                f"{ALLOWLIST_PATH.name} entry {number} is not an entry: {raw!r}. "
                "An entry is '<module> <parameter> :: <reason>' - the reason is "
                "the point of the file, so a line without one is an omission "
                "that has lost its justification."
            )
        key, _, reason = line.partition("::")
        parts = key.split()
        if len(parts) != 2:
            raise AssertionError(
                f"{ALLOWLIST_PATH.name} entry {number} names {len(parts)} "
                f"field(s), expected '<module> <parameter>': {raw!r}"
            )
        module, parameter = parts
        entries[f"{module}.{parameter}"] = reason.strip()
    return entries


def _omitted_defaulted(calls) -> Dict[str, list]:
    """Every defaulted parameter a call site does not supply.

    Keyed ``"<module>.<parameter>"``, valued by the list of `run.py` lines that
    omit it. The supplied set is worked out the way Python would: positional
    arguments fill the positional parameters in declaration order, then the
    keywords, and an unpacked ``*args``/``**kwargs`` is treated as supplying
    everything - which is why `PartOne` separately forbids star-arguments at a
    stage call site, so this function is never asked to guess.
    """
    omitted: Dict[str, list] = {}
    for call in calls:
        signature = _signature(call["module"])
        parameters = list(signature.parameters.values())
        supplied = set()
        if call["starargs_positional"]:
            supplied |= {
                p.name
                for p in parameters
                if p.kind
                in (
                    inspect.Parameter.POSITIONAL_ONLY,
                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                )
            }
        else:
            positional = [
                p
                for p in parameters
                if p.kind
                in (
                    inspect.Parameter.POSITIONAL_ONLY,
                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                )
            ]
            supplied |= {p.name for p in positional[: call["n_positional"]]}
        supplied |= set(call["keywords"])
        if call["starargs"]:
            supplied |= {p.name for p in parameters}

        for parameter in parameters:
            if parameter.kind in (
                inspect.Parameter.VAR_POSITIONAL,
                inspect.Parameter.VAR_KEYWORD,
            ):
                continue
            if parameter.default is inspect.Parameter.empty:
                continue
            if parameter.name in supplied:
                continue
            omitted.setdefault(
                f"{call['module']}.{parameter.name}", []
            ).append(call["line"])
    return omitted


@pytest.fixture(scope="module")
def calls():
    found = _stage_calls(ast.parse(RUN_PY.read_text(encoding="utf-8")))
    assert found, "no stage_X.run(...) call site found in run.py at all"
    return found


@pytest.fixture(scope="module")
def aliases():
    return _callables()


class TestEveryCallSiteBinds:
    """Part one: the cheap half. Arity and names, at every call site."""

    def test_every_stage_module_reached_by_run_py_is_covered(self, calls, aliases):
        """The premise. An unaccounted import means unchecked wiring.

        Without this, deleting a call site - or importing a stage and calling it
        through a spelling this file does not recognise - shrinks the set under
        test and the remaining assertions get easier, which is the wrong
        direction for a test whose whole job is to notice a stage gaining a
        requirement.
        """
        reached = {c["module"] for c in calls}
        for module, alternatives in ALTERNATIVE_ENTRY_POINTS.items():
            reachable = importlib.import_module(module)
            missing = [
                name
                for name in alternatives
                if not hasattr(reachable, name)
            ]
            assert not missing, (
                f"{module} is listed as reached through {alternatives}, but "
                f"{missing} no longer exist on it. Update "
                "ALTERNATIVE_ENTRY_POINTS deliberately."
            )
        unchecked = sorted(reached - set(REQUIRED_PARAMETERS))
        assert not unchecked, (
            "run.py calls these stage modules and this file has no record of "
            f"their required parameters: {unchecked}. Add them to "
            "REQUIRED_PARAMETERS."
        )

    def test_the_call_site_count_matches_what_run_py_contains(self, calls):
        """`calls` is the whole denominator, stated so a reader can check it.

        `run.py` reaches 18 stage entry points today. If this number changes,
        the change is a wiring change and belongs in this file's table too.
        """
        assert len(calls) == len(REQUIRED_PARAMETERS), (
            f"{len(calls)} call sites for {len(REQUIRED_PARAMETERS)} recorded "
            "stages: a stage is called from more than one place, or a call site "
            "was added or removed. Both are deliberate edits to this file."
        )

    def test_each_call_site_supplies_every_required_parameter(self, calls):
        """The assertion. `bind` fails on a missing required argument.

        Sentinels, because `bind` checks names and counts and never values.
        """
        problems = []
        for call in calls:
            signature = _signature(call["module"])
            args = [
                _Sentinel(f"{call['alias']}#{i}")
                for i in range(call["n_positional"])
            ]
            kwargs = {name: _Sentinel(name) for name in call["keywords"]}
            try:
                signature.bind(*args, **kwargs)
            except TypeError as exc:
                problems.append(
                    f"{RUN_PY.name}:{call['line']} {call['alias']}.run(...) - {exc}"
                )
        assert not problems, (
            "these stage call sites do not supply every parameter their stage "
            "declares as required:\n  " + "\n  ".join(problems)
        )

    def test_no_call_site_passes_a_name_the_stage_does_not_declare(self, calls):
        """The other direction: a keyword the signature does not have.

        `bind` rejects this too, and it is checked separately so the message
        names the cause. A renamed parameter that leaves the call site behind
        raises `TypeError` at run time; `**kwargs` on the stage would swallow
        it instead, which is why this is asserted rather than assumed.
        """
        problems = []
        for call in calls:
            parameters = _signature(call["module"]).parameters
            for name in call["keywords"]:
                if name not in parameters and not any(
                    p.kind is inspect.Parameter.VAR_KEYWORD
                    for p in parameters.values()
                ):
                    problems.append(
                        f"{RUN_PY.name}:{call['line']} passes {name}=... but "
                        f"{call['alias']}.run declares {sorted(parameters)}"
                    )
        assert not problems, "\n".join(problems)

    def test_no_call_site_hides_arguments_behind_a_star(self, calls):
        """`*args` / `**kwargs` at a stage call site defeats this file.

        `bind` would be handed one sentinel for a whole unpacked sequence, so a
        required parameter could go missing and the bind would still succeed on
        arity alone. Nothing in `run.py` does this today; the assertion is what
        keeps it that way.
        """
        hiding = [
            f"{RUN_PY.name}:{c['line']} {c['alias']}.run{'(...)' if c['starargs_positional'] else ''}"
            f"{'**' if c['starargs'] else ''}"
            for c in calls
            if c["starargs"] or c["starargs_positional"]
        ]
        assert not hiding, (
            "these stage call sites unpack arguments, so this file cannot "
            f"prove what they supply: {hiding}"
        )


class TestTheRequirementIsDeclaredInTheSignature:
    """Part two: what "required" means here, pinned against the live code."""

    def test_the_recorded_requirements_match_the_live_signatures(self):
        problems = []
        for module, recorded in sorted(REQUIRED_PARAMETERS.items()):
            live = tuple(
                name
                for name, parameter in _signature(module).parameters.items()
                if parameter.default is inspect.Parameter.empty
                and parameter.kind
            )
            if live != recorded:
                problems.append(
                    f"{module}: signature requires {list(live)}, "
                    f"REQUIRED_PARAMETERS records {list(recorded)}"
                )
        assert not problems, (
            "a stage's required parameters changed:\n  " + "\n  ".join(problems)
            + "\n\nUpdate REQUIRED_PARAMETERS in the same commit, and check "
            "that run.py's call sites pass the new ones."
        )

    def test_the_table_is_not_empty(self):
        """Anti-vacuity. Every stage in it must require something.

        A table that had emptied to `()` everywhere would make `PartOne`
        trivially true - `bind` has nothing to be missing. Each stage in this
        pipeline does take a cohort and a mode, so an empty one is a defect in
        the extraction, not in the stage.
        """
        empty = sorted(m for m, r in REQUIRED_PARAMETERS.items() if not r)
        assert not empty, (
            f"these stages are recorded as requiring nothing: {empty}. Either "
            "the signature genuinely lost its parameters or the record is "
            "wrong; either way `bind` has nothing to check."
        )

    def test_every_recorded_stage_is_a_real_importable_module(self):
        missing = []
        for module in REQUIRED_PARAMETERS:
            try:
                loaded = importlib.import_module(module)
            except ImportError as exc:  # pragma: no cover - packaging break
                missing.append(f"{module}: {exc}")
                continue
            if not callable(getattr(loaded, "run", None)):
                missing.append(f"{module} has no callable run()")
        assert not missing, "\n".join(missing)


class _RootGuard:
    """Whether a stage body refuses a None root, or defaults it to config."""

    def __init__(self, module: str) -> None:
        self.module = module
        self.source = Path(inspect.getfile(importlib.import_module(module)))
        self.tree = ast.parse(self.source.read_text(encoding="utf-8"))
        self.refuses: set = set()
        self.defaults: Dict[str, str] = {}
        self._scan()

    def _scan(self) -> None:
        for node in ast.walk(self.tree):
            if isinstance(node, ast.If):
                names = _none_tests(node.test)
                if not names:
                    continue
                if _raises(node.body):
                    self.refuses.update(names)
                    continue
                for name in names:
                    accessor = _config_accessor(node)
                    if accessor:
                        self.defaults.setdefault(name, accessor)
            elif isinstance(node, ast.IfExp):
                # `phylo = Path(phylogeny_dir) if phylogeny_dir else
                # config.phylogeny_dir(mode)` - the same contract as an
                # `if/else`, and the idiom both `convergence` and
                # `cooccurrence` use for `phylogeny_dir`.
                if not isinstance(node.test, ast.Name):
                    continue
                accessor = _config_accessor(node.orelse)
                if accessor:
                    self.defaults.setdefault(node.test.id, accessor)

    def verdict(self, parameter: str, test_only: bool = False) -> str:
        if parameter in self.refuses:
            return "REFUSED"
        if parameter in self.defaults:
            return f"DEFAULTED:{self.defaults[parameter]}"
        if test_only:
            return "TEST_PATH_ONLY"
        return "UNDECLARED"


def _parents(tree: ast.Module) -> Dict[int, ast.AST]:
    out: Dict[int, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            out[id(child)] = node
    return out


def _inside_mode_test(
    lineno: int, tree: ast.Module, parents: Dict[int, ast.AST], mode: str = "TEST"
) -> bool:
    """Whether the call at `lineno` sits inside an ``if ... RunMode.<mode>``.

    Ascends from the call node through every enclosing ``if``, and requires one
    of them to compare against ``RunMode.<mode>``. Ascending through *all* of
    them, rather than stopping at the first, is deliberate: a call nested inside
    a non-mode ``if`` that is itself inside a mode branch is still gated by the
    mode branch, and stopping early would miss it.
    """
    node = next(
        (
            n
            for n in ast.walk(tree)
            if isinstance(n, ast.Call) and n.lineno == lineno
            and isinstance(n.func, ast.Attribute)
            and n.func.attr == "run"
        ),
        None,
    )
    if node is None:  # pragma: no cover - a call site the scanner cannot place
        return False
    while id(node) in parents:
        node = parents[id(node)]
        if not isinstance(node, ast.If):
            continue
        if any(
            isinstance(sub, ast.Attribute)
            and isinstance(sub.value, ast.Name)
            and sub.value.id == "RunMode"
            and sub.attr == mode
            for sub in ast.walk(node.test)
        ):
            return True
    return False


def _none_tests(node: ast.AST) -> set:
    """Parameter names compared against None by `if <name> is None`."""
    names = set()
    if not isinstance(node, ast.Compare):
        return names
    if len(node.ops) != 1 or not isinstance(node.ops[0], (ast.Is, ast.IsNot)):
        return names
    if not isinstance(node.comparators[0], ast.Constant):
        return names
    if node.comparators[0].value is not None:
        return names
    for side in (node.left, node.comparators[0]):
        if isinstance(side, ast.Name):
            names.add(side.id)
    return names


#: The `PipelineConfig` accessors that denote a root or the reference. A
#: defaulted root parameter is only accounted for if it falls back to one of
#: these; anything else is a value nobody chose.
CONFIG_ROOT_ACCESSORS = frozenset(
    {
        "assembly_root", "data_root", "db_root", "intermediate_root",
        "phylogeny_dir", "phenotype_dir", "reference_fasta", "reference_gff",
        "reports_root", "tool_output_root",
    }
)


def _config_accessor(node) -> str:
    """The root accessor called in `node`, or `""`. Takes a node or a list."""
    nodes = node if isinstance(node, list) else [node]
    for sub in (n for root in nodes for n in ast.walk(root)):
        if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Attribute):
            if sub.func.attr in CONFIG_ROOT_ACCESSORS:
                return sub.func.attr
    return ""


def _raises(body) -> bool:
    for node in body:
        for sub in ast.walk(node):
            if isinstance(sub, ast.Raise):
                return True
    return False


class TestDefaultedRootsAreEitherRefusedOrConfigDefaulted:
    """Part three: the quiet half of the same defect.

    `bind` proves a required parameter is supplied. It cannot see a parameter
    the stage declared **with a default**, because omitting it at the call site
    is legal Python and the stage then receives `None`. For the roots and the
    reference that is not a cosmetic difference: `None` means "nowhere to read
    from", and a stage that treats it as anything other than a refusal has a
    defect that no arity check will ever report.
    """

    @pytest.fixture(scope="class")
    def defaulted(self):
        found = {}
        for module in REQUIRED_PARAMETERS:
            guard = _RootGuard(module)
            test_only = module in TEST_ONLY_ENTRY_POINTS
            for name, parameter in _signature(module).parameters.items():
                if parameter.default is inspect.Parameter.empty:
                    continue
                if name not in ROOT_PARAMETERS:
                    continue
                found.setdefault(module, {})[name] = guard.verdict(name, test_only)
        return found

    def test_the_defaulted_root_set_is_not_empty(self, defaulted):
        """Anti-vacuity.

        If every root parameter were required outright this whole part would be
        vacuous, and would keep being vacuous while a future edit added one more
        defaulted root. Asserting the set is non-empty means the day that
        happens, this file says so instead of going quietly quiet.
        """
        undeclared = {
            m: {k: v for k, v in d.items() if v == "UNDECLARED"}
            for m, d in defaulted.items()
        }
        undeclared = {m: d for m, d in undeclared.items() if d}
        assert not undeclared, (
            "these defaulted root parameters are neither refused when None nor "
            f"defaulted to a config accessor: {undeclared}. Every defaulted "
            "root must do one or the other - see the module docstring."
        )

    def test_every_config_defaulted_root_is_recorded_in_this_file(self, defaulted):
        problems = []
        for module, verdicts in sorted(defaulted.items()):
            recorded = CONFIG_DEFAULTED_ROOT_PARAMETERS.get(module, {})
            for name, verdict in sorted(verdicts.items()):
                if not verdict.startswith("DEFAULTED:"):
                    if name in recorded:
                        problems.append(
                            f"{module}.{name} is recorded as defaulting to "
                            f"{recorded[name]!r}, but its body {verdict}. A "
                            "parameter that is refused is not defaulted, and "
                            "recording it as both hides which one it is."
                        )
                    continue
                actual = verdict.split(":", 1)[1]
                if recorded.get(name) != actual:
                    problems.append(
                        f"{module}.{name} defaults to {actual} in the body, but "
                        "CONFIG_DEFAULTED_ROOT_PARAMETERS records "
                        f"{recorded.get(name)!r}. Write it down in the same "
                        "commit that introduced the default."
                    )
        assert not problems, "\n".join(problems)

    def test_a_test_only_stage_stays_test_only(self):
        """The premise behind every ``TEST_PATH_ONLY`` verdict above.

        A defaulted root on a stage `run.py` only ever calls from a TEST branch
        is not a REAL-mode input, so there is nothing for it to refuse. The
        moment that stops being true the verdict is wrong, and silently so -
        which is why it is asserted directly rather than trusted.
        """
        tree = ast.parse(RUN_PY.read_text(encoding="utf-8"))
        parents = _parents(tree)
        calls = _stage_calls(tree)
        offenders = []
        for call in calls:
            if call["module"] not in TEST_ONLY_ENTRY_POINTS:
                continue
            if not _inside_mode_test(call["line"], tree, parents):
                offenders.append(
                    f"{RUN_PY.name}:{call['line']} calls "
                    f"{call['alias']}.run outside a TEST branch"
                )
        assert not offenders, (
            "these entry points were TEST-only and are now reached from another "
            f"mode: {offenders}. Their defaulted roots "
            "(stages/recombination.py:301-302) are neither required nor refused, "
            "so promoting them to REAL makes those parameters load-bearing and "
            "silently nullable."
        )

    def test_the_test_only_verdict_is_not_the_only_one_available(self, defaulted):
        """Anti-vacuity for `TEST_PATH_ONLY`.

        If every defaulted root resolved to `TEST_PATH_ONLY` then `PartThree`'s
        first assertion would pass for the wrong reason: nothing would have to
        be refused and nothing would have to be config-defaulted. Both of those
        verdicts have to occur, and they do - `amr`/`mlst` refuse `data_root`,
        `convergence`/`cooccurrence` default `phylogeny_dir` from config.
        """
        verdicts = {v for d in defaulted.values() for v in d.values()}
        for required in ("REFUSED", "TEST_PATH_ONLY"):
            assert any(v == required for v in verdicts), (
                f"no defaulted root parameter is {required} any more, so that "
                "verdict is untested and may no longer describe the code"
            )
        assert any(v.startswith("DEFAULTED:") for v in verdicts), (
            "no defaulted root parameter falls back to a config accessor any "
            "more; CONFIG_DEFAULTED_ROOT_PARAMETERS may be describing a shape "
            "the stages no longer have"
        )

    def test_a_refused_root_is_not_also_recorded_as_defaulted(self, defaulted):
        """Named separately, because the two are opposite contracts.

        `amr.run` and `mlst.run` both take `data_root: Optional[Path] = None`
        and both **refuse** a REAL run without it. `convergence.run` takes
        `phylogeny_dir: Optional[Path] = None` and **defaults** it from config.
        Those are different promises to the operator and the table records only
        the second, so a parameter quietly changing from one to the other is a
        behaviour change with nothing to fail.
        """
        refused = {
            module: sorted(n for n, v in d.items() if v == "REFUSED")
            for module, d in defaulted.items()
        }
        refused = {m: n for m, n in refused.items() if n}
        assert refused, (
            "no defaulted root parameter is refused when None any more. Either "
            "the stages changed or the `_RootGuard` scan stopped matching the "
            "shape they use - see stages/convergence.py:365 for the idiom."
        )

    def test_no_root_parameter_is_left_optional_without_a_verdict(self, defaulted):
        """Named separately: the enumeration itself must be complete.

        `PartThree::test_the_defaulted_root_set_is_not_empty` is about the
        roots that are declared; this is about making sure no root parameter
        escaped the declaration because its name is not in `ROOT_PARAMETERS`.
        A root that stops being recognised is a root nobody is checking.
        """
        unrecognised = {
            m: sorted(
                name
                for name, parameter in _signature(m).parameters.items()
                if name not in ROOT_PARAMETERS
                and "root" in name
                and parameter.default is not inspect.Parameter.empty
            )
            for m in REQUIRED_PARAMETERS
        }
        unrecognised = {m: names for m, names in unrecognised.items() if names}
        assert not unrecognised, (
            "these defaulted parameters contain 'root' but are not in "
            f"ROOT_PARAMETERS, so nothing checks them: {unrecognised}"
        )


#: The six omissions round 12's WIRE.md A3 found to change behaviour, and the
#: two more the same round's cooccurrence and annotation-reuse branches recorded
#: the same way. Every one is now passed at its call site.
#:
#: Kept as a list rather than deleted with the entries, because the assertion
#: below is the point: an omission that was a defect must not reappear as a
#: *defended* omission. Naming them here is what lets that be checked without
#: re-deriving the findings from prose.
#:
#: `convergence.n_samples` is the one whose stage body was fixed first
#: (`stages.convergence.analysed_sample_count`, replacing `len(manifest)` as the
#: `widespread_fraction` denominator - CONV measured the inversion that published
#: 8-of-9 as `recurrent_convergent`). The call site now passes
#: `_assemblies_analysed` as well, so the denominator is stated rather than
#: inherited.
FIXED_FINDINGS = (
    "papipeline.stages.annotation.database_version",
    "papipeline.stages.annotation.event_sink",
    "papipeline.stages.annotation.run_key",
    "papipeline.stages.annotation.store",
    "papipeline.stages.annotation.tool_version",
    "papipeline.stages.convergence.n_samples",
    "papipeline.stages.convergence.phylogeny_dir",
    "papipeline.stages.cooccurrence.phylogeny_dir",
    "papipeline.stages.gwas.engine",
    "papipeline.stages.mechanisms.intermediate_root",
    "papipeline.stages.regulators.locus_coverage",
    "papipeline.stages.regulators.min_locus_coverage",
)


class TestEveryOmissionIsAccountedFor:
    """Part four: the defaulted parameters `bind` cannot see.

    `PartOne` proves a stage's *required* parameters are supplied, and
    `PartThree` proves each defaulted **root** is either refused when `None` or
    defaulted to a config accessor. Neither covers the general case: a defaulted
    parameter that is not a root, whose default is an empty container or a small
    constant, and whose body treats the default as a real value.

    The shape of the hazard is the one this file exists for, at its quietest.
    Omitting such a parameter is legal Python, so nothing raises. The stage
    receives the default and carries on. Whether that is a refusal, a documented
    fallback or a silent wrong answer is settled inside the stage body, per
    parameter, with no caller in the loop - and the scientific tables are written
    either way. Three real instances were found by this part in round 12 and are
    recorded in the allowlist file as FINDING rather than defended.

    So the assertion is not "no parameter is omitted" - that would be wrong, and
    most of these omissions are deliberate. It is that every omission is written
    down with a reason. An omission nobody has looked at is the defect; a
    reviewed one is a decision.
    """

    @pytest.fixture(scope="class")
    def omitted(self, calls):
        found = _omitted_defaulted(calls)
        assert found, (
            "no defaulted parameter is omitted at any call site, so this part "
            "is vacuous. Either `run.py` now passes every defaulted parameter "
            "explicitly - in which case the allowlist file should be deleted "
            "deliberately - or the extraction stopped matching the call sites."
        )
        return found

    def test_every_omission_has_an_allowlist_entry(self, omitted):
        recorded = _allowlist()
        undocumented = sorted(set(omitted) - set(recorded))
        assert not undocumented, (
            "these defaulted parameters are not passed at a call site in "
            "run.py and have no entry in "
            f"{ALLOWLIST_PATH.relative_to(REPO)}: {undocumented}\n\n"
            "Either pass the argument, or add a line "
            "'<module> <parameter> :: <reason>' with a reason specific to that "
            "call site. An omission nobody has looked at is the defect this "
            "file exists to catch."
        )

    def test_no_entry_refers_to_an_omission_that_no_longer_exists(self, omitted):
        recorded = _allowlist()
        stale = sorted(set(recorded) - set(omitted))
        assert not stale, (
            f"{ALLOWLIST_PATH.relative_to(REPO)} justifies omissions that no "
            f"longer happen: {stale}. Either run.py stopped passing the "
            "parameter, or the stage stopped declaring it, or the stage is no "
            "longer called. Delete the entry in the same commit, so the file "
            "cannot keep vouching for a claim that is no longer being made."
        )

    def test_every_entry_carries_a_reason_worth_reading(self, omitted):
        """A reason that would fit any entry is not a reason.

        Checked mechanically rather than by taste: a reason must name the
        parameter or the stage it is justifying, or cite the line of code that
        settles it. That is the minimum for it to be checkable by the next
        reader rather than merely present.
        """
        recorded = _allowlist()
        thin = []
        for key, reason in sorted(recorded.items()):
            if key not in omitted:
                continue  # staleness is its own assertion above
            module, _, parameter = key.rpartition(".")
            stem = parameter.rstrip("_")
            words = {w for w in re.findall(r"[a-z]+", reason.lower()) if len(w) > 3}
            if stem not in reason and not words and len(reason) < 40:
                thin.append(f"{key}: {reason!r}")
            elif len(reason) < 25:
                thin.append(
                    f"{key}: {reason!r} - {len(reason)} chars, too short to say "
                    "why an omission is acceptable"
                )
        assert not thin, (
            "these allowlist reasons do not say why the omission is acceptable:\n  "
            + "\n  ".join(thin)
        )


    def test_a_finding_is_now_passed_at_its_call_site(self, omitted):
        """The six findings were fixed, not re-justified.

        **R10, round 12 Phase 3: this replaces
        `test_the_findings_are_still_recorded_as_findings` and
        `test_a_fixed_finding_says_which_side_fixed_it`.** Both asserted that a
        FINDING marker survived in the allowlist file. That was the right
        assertion while the omissions were unfixed and the *wrong* one once they
        were fixed: an allowlist entry is a claim that an omission is benign, so
        keeping one for a parameter the call site now passes would have the file
        vouching for a claim that is no longer being made - which
        `test_no_entry_refers_to_an_omission_that_no_longer_exists` also rejects.

        So the assertion is now the stronger, direct one: each of these keys
        must be **absent from the derived omission set**. A later edit that
        deletes one of these arguments fails here, and it fails with a message
        naming the finding rather than reporting a generic staleness.
        """
        reintroduced = sorted(key for key in FIXED_FINDINGS if key in omitted)
        assert not reintroduced, (
            "these omissions were found to change behaviour and were fixed at "
            f"their call sites; they are omitted again: {reintroduced}. Either "
            "the fix was reverted, or a stage's signature grew the parameter "
            "back with a default after the fix. Restore the argument rather "
            "than adding an allowlist entry - see "
            "pa-artifacts/round12/omissions.md for what each fix was."
        )

    def test_no_allowlist_entry_vouches_for_a_fixed_finding(self):
        """No entry may exist for a parameter the call site now passes.

        Without this, adding one back would satisfy the entry-existence checks
        and fail only the staleness check, with a message about bookkeeping
        rather than about a finding being laundered into a defence. The two are
        checked separately because they fail for different reasons and a reader
        needs to know which one they are looking at.
        """
        recorded = _allowlist()
        vouched = sorted(key for key in FIXED_FINDINGS if key in recorded)
        assert not vouched, (
            "the allowlist justifies omissions that no longer happen: "
            f"{vouched}. An entry is a claim that an omission is benign, and "
            "these parameters are now passed at their call sites, so the entry "
            "is vouching for a claim nobody is making. Delete it; do not "
            "reword it."
        )
