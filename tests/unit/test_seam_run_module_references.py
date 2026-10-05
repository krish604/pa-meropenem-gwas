"""Seam G1: every ``module.attr`` ``run.py`` reaches must resolve.

**The hazard, precisely.** ``papipeline/run.py`` imports its collaborators two
ways. ``from .stages import amr as stage_amr`` binds a MODULE, and the call site
then reads ``stage_amr.run(...)``. A rename of that ``run`` is caught by nothing:
the import still succeeds, the module still imports, every test that calls
``stages.amr.run`` directly still passes, and the failure appears only when the
orchestrator reaches that line - i.e. at the point a real run would be halfway
through stage 4. The other import form, ``from .x import name``, does NOT have
this hazard: a rename there raises ``ImportError`` at import time, which is why
this file is scoped to the module-attribute shape and not to "every name
``run.py`` imports". (That scoping is what a previous round of this work got
wrong by trying to cover both and then deleting the result.)

**Scope, and why it is exactly this.** The list of attributes under test is
DERIVED from ``run.py``'s own source by AST walk, not typed in here. A typed-in
list rots silently: it survives the rename that made it wrong and reports green.
Deriving it means the test's subject is "whatever ``run.py`` currently reaches",
which is the claim worth making.

**Non-vacuity is asserted, because a scan of nothing passes.** Three guards:
the total number of checks must exceed zero; the number of *first-party* checks
must exceed zero (stdlib modules like ``json`` would otherwise satisfy the first
guard alone); and the specific attributes WIRE-ALL wired in round 12 must be in
the derived set, which fails loudly if the extraction stops working rather than
silently shrinking. A test that scans nothing must not be able to report green.

**Not asserted:** that each attribute is *called*, that it is called on every
path, or that the call is correct. That is W5's subject - call-site tests that
drive ``run_pipeline`` with an injected runner - and duplicating it here would
produce two tests that both rot.
"""

from __future__ import annotations

import ast
import inspect
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pytest

RUN_MODULE = "papipeline.run"

#: Attributes this module must find in the derived set. Round 12's W1 wired
#: these; they are the anchors that turn "the scan found something" into "the
#: scan found the THING", because a broken extraction empties the list and an
#: empty list cannot contain a named member.
W1_WIRED_ATTRIBUTES: Tuple[Tuple[str, str], ...] = (
    ("tblastn_adapter", "run_tblastn"),
    ("stage_gwas_features", "produce_gwas_features"),
    ("stage_regulators", "produce_regulator_variants"),
    ("stage_gwas", "PyseerEngine"),
)


@dataclass(frozen=True)
class ModuleAttributeReference:
    """One ``module.attr`` read in ``run.py``, and where it is read."""

    binding: str
    attribute: str
    line: int
    module_name: str

    def describe(self) -> str:
        return (
            f"{self.binding}.{self.attribute} at "
            f"{Path(RUN_MODULE.replace('.', '/')).name}:{self.line} "
            f"(resolves on {self.module_name})"
        )


def _run_namespace() -> Dict[str, object]:
    import importlib

    return vars(importlib.import_module(RUN_MODULE))


def _module_bindings(namespace: Dict[str, object]) -> Dict[str, object]:
    """Top-level names in run.py's namespace that are MODULE objects.

    Decided by asking the objects, not by reading the import statements: a
    binding is in scope for this test if ``run.<name>`` is a module, whatever
    syntax brought it in. That is the same question the hazard is about.
    """
    return {
        name: value
        for name, value in namespace.items()
        if inspect.ismodule(value)
    }


def _references_in_source(tree: ast.AST, modules: Dict[str, object]):
    """Every ``<module binding>.<attr>`` in the tree, with its line.

    Only attribute access on a bare Name is collected. ``a.b.c`` contributes
    ``b`` on ``a`` (and ``c`` on the *result* of ``a.b``, which is not a module
    binding and is therefore out of scope), and ``getattr(m, "x")`` is not
    attribute syntax at all - both are honest exclusions, stated so a reader
    knows what is not covered rather than discovering it.
    """
    found: Dict[Tuple[str, str], int] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute):
            continue
        if not isinstance(node.value, ast.Name):
            continue
        binding = node.value.id
        if binding not in modules:
            continue
        key = (binding, node.attr)
        # Keep the FIRST line: the earliest read is the one that fails first.
        found.setdefault(key, node.lineno)
    return found


@pytest.fixture(scope="module")
def run_source_tree():
    import importlib

    module = importlib.import_module(RUN_MODULE)
    return ast.parse(Path(module.__file__).read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def checked_attributes(run_source_tree) -> List[ModuleAttributeReference]:
    namespace = _run_namespace()
    modules = _module_bindings(namespace)
    references = _references_in_source(run_source_tree, modules)
    return sorted(
        (
            ModuleAttributeReference(
                binding=binding,
                attribute=attribute,
                line=line,
                module_name=getattr(modules[binding], "__name__", "?"),
            )
            for (binding, attribute), line in references.items()
        ),
        key=lambda r: (r.binding, r.attribute),
    )


class TestTheScanIsNotVacuous:
    """A test that scans nothing passes. These are the guards against that."""

    def test_the_scan_found_something(self, checked_attributes):
        assert checked_attributes, (
            "no `module.attr` reference was derived from run.py. Either run.py "
            "no longer binds any module it reads attributes from - in which case "
            "this file's subject has gone - or the AST extraction is broken, "
            "which is the more likely answer and the worse one, because the "
            "test would then report green while checking nothing."
        )

    def test_the_scan_found_first_party_modules_not_just_the_standard_library(
        self, checked_attributes
    ):
        """``json.dumps`` is not the hazard. ``stage_amr.run`` is.

        Without this, a regression that emptied the first-party part of the list
        while leaving the four stdlib modules in it would still satisfy
        `test_the_scan_found_something`.
        """
        first_party = [
            ref
            for ref in checked_attributes
            if ref.module_name.startswith("papipeline")
        ]
        assert first_party, (
            "no first-party module attribute was checked; the only references "
            f"found are {[r.describe() for r in checked_attributes]}, which are "
            "standard library and cannot be renamed by a refactor in this repo"
        )

    @pytest.mark.parametrize("binding,attribute", W1_WIRED_ATTRIBUTES)
    def test_a_specifically_wired_attribute_is_in_the_scanned_set(
        self, checked_attributes, binding, attribute
    ):
        """The anchors. If these are absent the extraction failed, not the code.

        Each of these was wired into run.py by round 12's W1 and is reached
        through a module binding, so each is exactly what this file exists to
        check. Their absence means the scan is broken; their presence means the
        scan is looking at the code that was actually written.
        """
        scanned = {(ref.binding, ref.attribute) for ref in checked_attributes}
        assert (binding, attribute) in scanned, (
            f"{binding}.{attribute} is reached in run.py but was not derived "
            "by the scan. The extraction has regressed - check "
            "_references_in_source before concluding anything about run.py."
        )

    def test_the_scanned_set_covers_every_stage_module_run_calls(self):
        """The bulk of the risk is the stage modules; count them, do not list them.

        Not an exhaustive list of stages, because adding a stage would then fail
        here for a reason unrelated to this file's subject. A floor instead: the
        orchestrator's whole job is calling stage modules, so a set of scanned
        attributes that has lost most of them is broken, not merely smaller.
        """
        import papipeline.run as run

        namespace = vars(run)
        modules = _module_bindings(namespace)
        stage_modules = {
            name for name, mod in modules.items()
            if name.startswith("stage_")
            and getattr(mod, "__name__", "").startswith("papipeline.stages")
        }
        assert stage_modules, "no stage module bindings found in run.py"

        import papipeline.run as run_module_instance

        references = _references_in_source(
            ast.parse(Path(run_module_instance.__file__).read_text(encoding="utf-8")),
            modules,
        )
        scanned_stage_bindings = {
            binding for (binding, _) in references if binding in stage_modules
        }
        missing = stage_modules - scanned_stage_bindings
        assert not missing, (
            f"run.py imports these stage modules but reaches no attribute on "
            f"them: {sorted(missing)}. Either they are unused imports - which is "
            "its own finding - or the scan is broken."
        )


class TestEveryReachedModuleAttributeResolves:
    """The assertion itself: the attribute exists, at the point of use."""

    def test_every_reference_resolves(self, checked_attributes):
        namespace = _run_namespace()
        missing: List[str] = []
        for reference in checked_attributes:
            module = namespace[reference.binding]
            if not hasattr(module, reference.attribute):
                missing.append(reference.describe())
        assert not missing, (
            f"{len(missing)} of {len(checked_attributes)} module attribute(s) "
            "run.py reaches do not resolve. Each is a failure that raises only "
            "when the orchestrator reaches the line, never at import time and "
            "never in a test that calls the stage directly:\n  "
            + "\n  ".join(missing)
        )

    def test_the_check_covers_one_reference_per_attribute_not_per_line(
        self, checked_attributes
    ):
        """One entry per (module, attribute), however many times it is read.

        Stated because the count is the non-vacuity evidence: if a future change
        collapsed this to one attribute the other guards would still pass, and
        this is the assertion that makes the number mean something.
        """
        pairs = [(ref.binding, ref.attribute) for ref in checked_attributes]
        assert len(pairs) == len(set(pairs))
        assert len(pairs) > 30, (
            f"only {len(pairs)} module attributes are checked, which is below "
            "the floor this file sets for itself; run.py is an orchestrator "
            "whose body is almost entirely module calls, so a set this small "
            "means the extraction regressed"
        )

    def test_every_reference_names_a_real_line_in_run_py(
        self, checked_attributes
    ):
        """A line number a reader cannot navigate to is not evidence.

        Cheap, and it catches the failure where the scan is pointed at a
        different file than the one it names.
        """
        import papipeline.run as run

        source_lines = Path(run.__file__).read_text(encoding="utf-8").splitlines()
        for reference in checked_attributes:
            assert 1 <= reference.line <= len(source_lines), (
                f"{reference.describe()} names a line outside "
                f"{Path(run.__file__).name} (which has {len(source_lines)})"
            )

    def test_the_module_bindings_are_first_party_or_stdlib_not_something_else(
        self, checked_attributes
    ):
        """Guards against the scan picking up an injected global.

        run.py is imported into a test process that other fixtures have
        monkeypatched before. A binding that is neither a stdlib module nor a
        `papipeline` one is something this test should refuse to reason about
        rather than assert about.
        """
        for reference in checked_attributes:
            name = reference.module_name
            assert name.split(".")[0] in (
                "papipeline",
                "json",
                "platform",
                "shutil",
                "sys",
            ), (
                f"{reference.describe()} resolves on {name}, which is neither "
                "a papipeline module nor one of the standard library modules "
                "run.py imports"
            )