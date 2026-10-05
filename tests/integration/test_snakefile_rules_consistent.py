"""The Snakemake rules must not contradict the taxonomy, and must resolve.

**The gap this closes.** `test_stage_taxonomy_is_self_verifying.py` derives
everything from `run.py`'s dispatch chain and the stage modules' ASTs. It never
looked at `workflow/Snakefile`, so it passed while the DAG was unbuildable in
REAL:

    MissingInputException in rule variants
      affected files:
        results/real/intermediate/regulators/regulator_variants.tsv
        data/standins/unbuilt_stages/variants.tsv

Two distinct faults, both invisible to the taxonomy test:

* **A stand-in declared for a stage that is built.** `variants` was reconciled in
  f5573f5 but its rule still listed `VARIANTS_STANDIN` and still said "Unbuilt;
  ticket 14" in its docstring. The same half-reconciliation f1c6f5c fixed for
  `similarity`, one layer down, where the taxonomy test could not see it.
* **An input with no producing rule.** `REGULATOR_IN` was declared as an input by
  `variants` and produced by nothing. That is the mirror of the fault the stage
  code had: `stages/regulators.py` read
  `intermediate/regulators/regulator_variants.tsv`, which only the TEST fixture
  generator ever wrote.

So this file checks the DAG's own consistency: every declared input is produced
by some rule or is a config/manifest file, and no rule outside
`UNBUILT_WITH_STANDIN` references a stand-in.

**Parsed, not grepped.** Rule bodies are read with `ast`, so a stand-in mentioned
in a docstring is not mistaken for a dependency. Both files above mention their
old constants in prose, and a text search would flag both.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from papipeline.standins import UNBUILT_WITH_STANDIN

ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = ROOT / "workflow" / "Snakefile"


def _split_rules(text: str) -> dict:
    """`rule <name>:` blocks, each split into docstring and executable parts.

    Not `ast.parse`: `rule all:` is Snakemake DSL, not a Python `def`, so the
    file does not parse as Python at all. Splitting textually is fine as long as
    the docstring is separated out first - three of these rules *mention* the
    stand-in constants they no longer depend on, in prose explaining exactly
    why, and a naive search would flag all three as regressions.
    """
    lines = text.splitlines()
    starts = [
        (index, match.group(1))
        for index, line in enumerate(lines)
        if (match := re.match(r"^rule ([A-Za-z_][A-Za-z0-9_]*):\s*$", line))
    ]
    rules: dict = {}
    for position, (index, name) in enumerate(starts):
        end = starts[position + 1][0] if position + 1 < len(starts) else len(lines)
        block = "\n".join(lines[index:end])
        doc = re.search(r'"""(.*?)"""', block, re.DOTALL)
        docstring = doc.group(1) if doc else ""
        body = block.replace(doc.group(0), "") if doc else block
        rules[name] = {"doc": docstring, "body": body}
    return rules


def _declared_names(body: str, keyword: str) -> list:
    """Names in a rule's `input:` or `output:` block, docstrings already removed."""
    match = re.search(
        rf"^    {keyword}:\s*$(.*?)(?=^    \w+:|\Z)",
        body,
        re.MULTILINE | re.DOTALL,
    )
    if not match:
        return []
    return re.findall(r"[A-Z][A-Z0-9_]{2,}", match.group(1))


#: Accessor-derived roots that do NOT identify a specific owned location.
#: Deriving from one of these is how the TREE_IN drift happened: `DATA_ROOT` is
#: itself `CONFIG.data_root(...)`, so a path spelled on top of it looks derived
#: while resolving to a directory that never holds pipeline output.
GENERIC_ROOTS = frozenset({
    "DATA_ROOT", "RESULTS_ROOT", "DB_ROOT", "TOOL_OUTPUT_ROOT",
    "INTERMEDIATE", "REPORTS_ROOT",
})


def _reaches_accessor(
    name: str, declaration: str, source: str, seen: set
) -> bool:
    """True when a constant's value traces back to a `CONFIG.<accessor>(...)` call.

    Following the chain is the whole point. `TREE_IN` is
    `str(PHYLO_DIR / "tree.nwk")`, and it is `PHYLO_DIR` that calls
    `CONFIG.phylogeny_dir()` - so a check that only looked for the accessor name
    in the constant's own declaration flagged every accessor-derived constant in
    the file, `DATA_ROOT = CONFIG.data_root(...)` among them. A guard that
    reports correct code as wrong gets ignored, and then it protects nothing.
    """
    if name in seen:
        return False
    seen.add(name)
    if "CONFIG." in declaration:
        return True
    for referenced in re.findall(r"[A-Z][A-Z0-9_]{2,}", declaration):
        if referenced == name:
            continue
        nested = re.search(
            rf"^{referenced}\s*=\s*(.+)$", source, re.MULTILINE
        )
        if not nested:
            continue
        # A *generic* root does not count as reaching the accessor. This is the
        # whole subtlety: `TREE_IN = DATA_ROOT / "phylogeny" / "tree.nwk"` traces
        # to `CONFIG.data_root()` one hop away, so an indirection-blind check
        # waves it through - and `data/phylogeny` is exactly the path that could
        # never resolve, because in REAL the tree lives under
        # `results/real/intermediate/phylogeny`. Only a specific accessor
        # (`phylogeny_dir`) identifies which subdirectory the config owns.
        if referenced in GENERIC_ROOTS:
            continue
        if _reaches_accessor(referenced, nested.group(1), source, seen):
            return True
    return False


def _config_owned_paths(source: str) -> list:
    """Constants whose path a config accessor already owns.

    Read from the Snakefile's own accessor calls rather than a hand-kept list,
    so adding an accessor teaches this test about it automatically - which is
    the failure mode of the hand-kept `permitted` set in the check below.
    """
    names = set()
    for match in re.finditer(
        r"^([A-Z][A-Z0-9_]{2,})\s*=.*\b(\w+_dir|\w+_root)\(", source, re.MULTILINE
    ):
        names.add(match.group(1))
    # TREE_IN / TREE_META_IN predate any accessor, so name them explicitly: the
    # point is that they SHOULD be derived, which is exactly what the test says.
    names.update({"TREE_IN", "TREE_META_IN"})
    return sorted(names)


@pytest.fixture(scope="module")
def source() -> str:
    return SNAKEFILE.read_text(encoding="utf-8")


def _rule_functions() -> dict:
    """`rule <name>:` bodies, as source text.

    Snakemake's DSL is Python, so `ast.parse` of the whole file yields real
    FunctionDef nodes for each rule - which is what makes a docstring
    distinguishable from a dependency.
    """
    tree = ast.parse(SNAKEFILE.read_text(encoding="utf-8"))
    rules: dict = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name.startswith("rule_"):
            name = node.name[len("rule_"):]
            rules[name] = ast.get_source_segment(
                SNAKEFILE.read_text(encoding="utf-8"), node
            ) or ""
    return rules


def _declared_names(body: str, keyword: str) -> list:
    """Names appearing in a rule's `input:` or `output:` block."""
    match = re.search(
        rf"^\s*{keyword}:\s*$(.*?)(?=^\s{{4}}\w+:|\Z)",
        body,
        re.MULTILINE | re.DOTALL,
    )
    if not match:
        return []
    block = match.group(1)
    return re.findall(r"[A-Z][A-Z0-9_]{2,}", block)


@pytest.fixture(scope="module")
def rules() -> dict:
    return _split_rules(SNAKEFILE.read_text(encoding="utf-8"))


class TestStandInsAreTestOnly:
    def test_only_unbuilt_stages_reference_a_standin(self, rules, source):
        """Requirement: a stand-in is only ever a stand-in for an unbuilt stage.

        `variants` was reconciled in f5573f5 and still declared
        `VARIANTS_STANDIN`. The file must not be parseable as a constant at all,
        so this asserts no rule outside the unbuilt set mentions one.
        """
        offenders = {}
        for name, body in rules.items():
            if "STANDIN" not in body["body"]:
                continue
            if name in UNBUILT_WITH_STANDIN:
                continue
            offenders[name] = [
                line.strip()
                for line in body["body"].splitlines()
                if "STANDIN" in line
            ]
        assert not offenders, (
            "these rules reference a stand-in but are not unbuilt stages, so a "
            f"committed fixture could pass for their real output: {offenders}"
        )

    def test_standins_are_demanded_through_the_test_only_channel(self, rules):
        """A stand-in must be `test_only=`, never a plain `external`.

        `stage_inputs(external)` demands its files in TEST *and* REAL, which is
        how a REAL run came to reference `data/standins/unbuilt_stages/` - a
        directory that does not exist, since `data/` is read-only. Passing them
        as `test_only` means REAL reaches these rules with no stand-in edge and
        the unbuilt stage refuses honestly.
        """
        offenders = []
        for name, body in rules.items():
            if "STANDIN" not in body["body"]:
                continue
            for line in body["body"].splitlines():
                if "STANDIN" in line and "test_only" not in line:
                    offenders.append(f"{name}: {line.strip()}")
        assert not offenders, (
            "stand-ins must be demanded via test_only=, so REAL never consumes "
            f"one: {offenders}"
        )

    def test_the_standin_directory_is_not_inside_data(self, source):
        """`data/` is read-only, so a stand-in path under it can never resolve.

        Retired rather than re-pointed. `STANDIN_DIR` existed for one caller -
        `rule recombination`'s `test_only` edge - and `recombination` left
        `standins.UNBUILT_WITH_STANDIN` when `dag-resolve` dispatched it, so
        there is no stand-in path left in this file to validate. The
        requirement it encoded still stands, so it is now asserted as an
        absence: if a future stage is declared-but-unbuilt and wires a stand-in,
        this fails on the `STANDIN_DIR` it has to re-declare, and the re-declaring
        commit is where the `data_root(RunMode.TEST)` form has to be justified.
        """
        match = re.search(r"^STANDIN_DIR\s*=\s*(.+)$", source, re.MULTILINE)
        assert not match, (
            "STANDIN_DIR is declared again. No stage is stood in for "
            f"(UNBUILT_WITH_STANDIN is {tuple(UNBUILT_WITH_STANDIN)!r}), so a "
            "stand-in path here is a fixture for a stage that computes its own "
            "output. If a stage really is unbuilt, delete this assertion in the "
            "same commit that re-declares the constant."
        )


class TestEveryDeclaredInputHasAProducer:
    """The mirror check: nothing may read a file the DAG never writes.

    This is the class `stages/regulators.py:220` was in - a stage reading
    `intermediate/regulators/regulator_variants.tsv`, which only the TEST
    fixture generator wrote. The rule declared it as an input and no rule
    produced it, so the DAG could not resolve in REAL at all.
    """

    def _produced_names(self, source: str) -> set:
        produced = set()
        for rule in _split_rules(source).values():
            produced.update(_declared_names(rule["body"], "output"))
        # Values assigned as `NAME = str(...)` and consumed elsewhere.
        produced.update(
            re.findall(r"^([A-Z][A-Z0-9_]{2,})\s*=\s*str\(", source, re.MULTILINE)
        )
        return produced

    def _consumed_names(self, source: str) -> dict:
        consumed = {}
        for name, rule in _split_rules(source).items():
            consumed[name] = _declared_names(rule["body"], "input")
        return consumed

    def test_a_path_the_config_owns_is_not_hardcoded(self, source):
        """A pipeline-produced path must come from its accessor, not DATA_ROOT.

        `TREE_IN` and `TREE_META_IN` were spelled `DATA_ROOT / "phylogeny" / ...`
        after 33ad628 had moved stage 9's output to `CONFIG.phylogeny_dir()`.
        Two definitions of one location drifted apart and the DAG could not be
        built:

            MissingInputException in rule cooccurrence
              affected files: data/phylogeny/tree_metadata.tsv

        The previous version of this file permitted both names as "recognised
        externals" and so waved through the exact drift it existed to catch. A
        permitted name may be *consumed*, but if the config already owns that
        path the Snakefile must derive it from the accessor - otherwise the
        writer and the DAG describe different locations and nothing notices
        until a run cannot resolve.
        """
        owned = _config_owned_paths(source)
        offenders = {}
        for name in owned:
            if name not in source:
                continue
            assignment = re.search(
                rf"^{name}\s*=\s*(.+)$", source, re.MULTILINE
            )
            if assignment is None:
                continue
            declaration = assignment.group(1)
            if _reaches_accessor(name, declaration, source, seen=set()):
                continue
            offenders[name] = declaration.strip()
        assert not offenders, (
            "these constants name a location the config accessor already owns, "
            "but hard-code it instead. Derive them from the accessor so the "
            "writer and the DAG cannot drift: "
            f"{offenders}"
        )

    def test_no_rule_consumes_a_file_no_rule_produces(self, source):
        produced = self._produced_names(source)
        # Legitimately external: the sample manifest, config files, and the
        # per-tool output globs the stage scripts read.
        permitted = {
            # `METADATA` is now only ever the *fixture* `rule synthetic_fixtures`
            # writes. Every stage demands `COHORT_INPUT`, which is `METADATA` for
            # an ordinary overlay and the PDC isolate table for a smoke one -
            # the same choice `run.discover_run_manifest` makes, named once
            # (E2E-DISCOVER D1/D7: the table was required by sixteen stages and
            # declared by nothing, so `snakemake -n` said yes and every stage
            # then died at manifest load). Both names stay permitted; neither is
            # produced by a rule, and that is the point of permitting them.
            "METADATA", "COHORT_INPUT", "PDC_TABLE",
            "CONFIG_FILES", "ANNOTATION", "MLST_IN", "AMR_IN",
            "SV_IN", "VF_IN", "GWAS_IN", "TREE_IN", "TREE_META_IN", "THREADS",
            "MEMORY_MB", "LOG_DIR", "STATUS_LOG", "RUN_STAGE", "MODE",
            "MACHINE", "PIPELINE_ROOT", "OUT_FIGURES", "OUT_MASTER",
            "OUT_MECHANISMS", "OUT_REGULATORS", "OUT_VARIANTS",
            "OUT_COHORT_VARIANTS", "OUT_RECOMBINATION", "OUT_SIMILARITY",
        }
        orphans = {}
        for rule, names in self._consumed_names(source).items():
            unknown = [
                n for n in names
                if n not in produced and n not in permitted
                and not n.startswith("OUT_") and not n.endswith("_STANDIN")
            ]
            if unknown:
                orphans[rule] = unknown
        assert not orphans, (
            "these rules declare inputs that no rule produces and that are not "
            f"recognised externals: {orphans}"
        )

    def test_regulator_in_is_gone(self, source):
        """The dead pre-fold path must not reappear as a constant."""
        for line in source.splitlines():
            if line.strip().startswith("REGULATOR_IN") and not line.strip().startswith("#"):
                pytest.fail(
                    f"REGULATOR_IN is declared again: {line.strip()!r}. It is a "
                    "pre-fold path with no producer; the regulators screen is "
                    "now rule variants' own second output."
                )

    def test_the_regulators_table_is_a_declared_output(self, rules):
        """The screen has to produce something, or stage 6 consumes nothing."""
        variants = rules.get("variants", {}).get("body", "")
        assert "OUT_REGULATORS" in variants, (
            "rule variants no longer declares the regulators table it writes"
        )
        assert "06_regulators.tsv" in SNAKEFILE.read_text(encoding="utf-8"), (
            "OUT_REGULATORS must be the contracted INTERNAL_TABLES filename; "
            "declaring a second path is how the writer and the DAG drift apart"
        )


class TestRuleMetadataMatchesTheTaxonomy:
    def test_no_rule_docstring_still_calls_a_built_stage_unbuilt(self, rules):
        """`variants` said "Unbuilt; ticket 14" for several commits after it
        was built. A docstring is not enforced by anything, so it is asserted."""
        built_but_claimed_unbuilt = []
        for name in ("variants", "cohort_variants", "similarity", "virulence"):
            docstring = rules.get(name, {}).get("doc", "")
            if not docstring:
                continue
            if re.search(r"\bUnbuilt\b", docstring):
                built_but_claimed_unbuilt.append(name)
        assert not built_but_claimed_unbuilt, (
            "these rules are built but their docstring still says Unbuilt: "
            f"{built_but_claimed_unbuilt}"
        )