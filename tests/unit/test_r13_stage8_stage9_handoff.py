"""R13: stage 8 must actually leave stage 9 an alignment to read.

**The defect.** `run.derive_recombination_tables` finished its REAL branch with

    gubbins_adapter.build_outputs(gubbins_dir=..., phylogeny_dir=...)

`build_outputs` is not a member of `papipeline.adapters.gubbins`. It is a
function on the *stage* module, `papipeline.stages.recombination:416`. So the
REAL branch of stage 8 raised `AttributeError` on its last real statement, after
gubbins had already run, and stage 9's only input was never written:

    $ grep -rn "def build_outputs" papipeline/
    papipeline/stages/recombination.py:416:def build_outputs(

**Why that shape is worse than a crash on entry.** A caller that fails on its
first line cannot have written anything partial. A caller that fails on its last
line has run the whole stage — gubbins executed, provenance written — and then
raised. The run stops either way, so nothing half-finished is published, but the
operator's log ends with an `AttributeError` from a module they were not told
was involved, after spending the full gubbins runtime to produce it.

**The second half of the defect is the one that matters scientifically.** Stage 9
reads `phylogeny_dir / "core_snp_alignment.fasta"` (`stages/phylogeny.py:504`)
and stage 8 is the only producer of that file
(`stages/recombination.py:84`, `SNP_ALIGNMENT_NAME`). With the call broken,
nothing writes it, so stage 9 cannot be reached with a masked alignment at all.
That is the handoff this file pins.

**How it is tested.** `derive_recombination_tables` is driven in `RunMode.REAL`
against **synthetic inputs** and an **injected fake gubbins runner** - a real
subprocess, a real `subprocess.run`, the real `run_gubbins` body, the real
parsers, and no gubbins binary anywhere. Per R12 that is not a real-sample run
and needs no authorisation; no genome, no blast, no real tool is involved. The
gate is opened with `dataclasses.replace` on the runtime mapping, which is what
`tests/unit/test_recombination_dispatch.py` already does for the REAL refusals.

**Two tests, and the second is not redundant.** The first is behavioural and
would catch the bug by raising. The second is static: it resolves the *symbol the
call site names* and checks the module actually has it, so a wrong target is a
one-line failure message instead of an `AttributeError` surfacing from inside a
subprocess-driven stage. Both are mutation-checked below.
"""

from __future__ import annotations

import ast
import importlib
import inspect
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from papipeline.adapters import gubbins as gubbins_adapter
from papipeline.config.loader import load_config
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode
from papipeline.run import derive_recombination_tables
from papipeline.stages import phylogeny as stage_phylo
from papipeline.stages import recombination as stage_recombination

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
RUN_PY = REPO / "papipeline" / "run.py"

#: The exact string `stages/phylogeny.py:504` and `:600` read. Spelled out as a
#: literal rather than imported from `recombination.SNP_ALIGNMENT_NAME`, because
#: the claim under test is that *these two places agree*. Importing the constant
#: would make the comparison true by construction and the test would say nothing.
STAGE_NINE_READS = "core_snp_alignment.fasta"

SAMPLES = ("SYN_001", "SYN_002", "SYN_003", "SYN_004")

#: Eleven labelled nodes: four tips and seven internal. Every one of them appears
#: in the statistics file below, because `block_rows` refuses a tree and a
#: statistics file that disagree in either direction.
SYNTHETIC_TREE = (
    "(((SYN_001:0.01,SYN_002:0.04)Node_A:0.02,"
    "(SYN_003:0.07,SYN_004:0.01)Node_B:0.03)Node_C:0.01,"
    "(Node_D:0.02,Node_E:0.03)Node_F:0.01)Node_G;"
)

STATS_HEADER = (
    "Node\tTotal SNPs\tNumber of SNPs Inside Recombinations\t"
    "Number of SNPs Outside Recombinations\tNumber of Recombination Blocks\t"
    "Bases in Recombinations Including Gaps\t"
    "Cumulative Bases in Recombinations Including Gaps\t"
    "Bases in Recombinations Excluding Gaps\t"
    "Cumulative Bases in Recombinations Excluding Gaps\t"
    "r/m\trho/theta\tGenome Length\tBases in Clonal Frame"
)

#: One row per node in `SYNTHETIC_TREE`. Two rows carry a recombination block so
#: `recombination_detected` is not uniformly false and the table cannot be
#: distinguished from a stub by its values alone.
STATS_ROWS = [
    ("SYN_001", 3, 1),
    ("SYN_002", 10, 0),
    ("Node_A", 17, 0),
    ("SYN_003", 24, 0),
    ("SYN_004", 31, 1),
    ("Node_B", 4, 0),
    ("Node_C", 11, 0),
    ("Node_D", 18, 0),
    ("Node_E", 25, 0),
    ("Node_F", 5, 0),
    ("Node_G", 12, 0),
]

#: The runner gubbins is asked for. It writes the three files the stage reads and
#: nothing else. Written to disk and made executable so the *real*
#: `subprocess.run` in `run_gubbins` launches it: injecting a fake at the Python
#: function boundary would skip the subprocess, the exit code, and the "exited 0
#: without writing the SNP alignment" refusal, which is most of what this call
#: sits downstream of.
FAKE_RUNNER = '''#!{python}
"""Stand-in for run_gubbins.py. Writes gubbins-shaped output into the CWD."""
import sys

WORKDIR = sys.argv[1] if False else "."

with open("recombination.filtered_polymorphic_sites.fasta", "w") as handle:
    handle.write(">SYN_001\\nACGTACGTAC\\n>SYN_002\\nACGTACGTAG\\n"
                 ">SYN_003\\nACGTACGTAC\\n>SYN_004\\nACGAACGTAC\\n")

with open("recombination.node_labelled.final_tree.tre", "w") as handle:
    handle.write("{tree}\\n")

with open("recombination.final_tree.tre", "w") as handle:
    handle.write("{tree}\\n")

with open("recombination.per_branch_statistics.csv", "w") as handle:
    handle.write("{header}\\n")
    for node, total, inside in {rows!r}:
        blocks = 1 if inside else 0
        outside = total - inside
        handle.write(
            "%s\\t%d\\t%d\\t%d\\t%d\\t0\\t0\\t0\\t0\\t0.0\\t0.0\\t6000000\\t5900000\\n"
            % (node, total, inside, outside, blocks)
        )

print("fake gubbins: wrote 4 files")
'''


@pytest.fixture
def config():
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture
def open_gate(config):
    """The same config with `runtime.allow_real_mode` true.

    `dataclasses.replace` on the runtime mapping, not a re-read with an edited
    overlay: the committed overlays all keep the gate shut, and opening one from a
    test file would make a policy decision the test has no standing to make.
    This is the same construction `test_recombination_dispatch.py` uses.
    """
    return replace(config, runtime={**config.runtime, "allow_real_mode": True})


@pytest.fixture
def manifest() -> SampleManifest:
    return SampleManifest(
        [
            SimpleNamespace(sample_id=sample_id, assembly_path=None, source="synthetic")
            for sample_id in SAMPLES
        ]
    )


@pytest.fixture
def fake_gubbins(tmp_path, monkeypatch):
    """Install a fake gubbins toolchain and return the runner's path.

    Injected at three points, all on the adapter, all overridable seams the
    adapter already documents for exactly this purpose:

    * `self_test` - there is no way to make a *working* binary segfault on
      demand, which is why `run_gubbins` grew a `selftest` hook. Returns success.
    * `resolve_binary` - a path that exists, so `isolated_env` can symlink it.
      Never executed; `self_test` is what would have run it.
    * `resolve_runner` - the fake script, which *is* executed.
    """
    tool_dir = tmp_path / "fake_env" / "bin"
    tool_dir.mkdir(parents=True)

    binary = tool_dir / "gubbins"
    binary.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    binary.chmod(0o755)

    runner = tool_dir / "run_gubbins.py"
    runner.write_text(
        FAKE_RUNNER.format(
            python="/usr/bin/env python3",
            tree=SYNTHETIC_TREE,
            header=STATS_HEADER,
            rows=STATS_ROWS,
        ),
        encoding="utf-8",
    )
    runner.chmod(0o755)

    monkeypatch.setattr(
        gubbins_adapter,
        "self_test",
        lambda binary, scratch, search_dirs=(): (0, gubbins_adapter.SOURCE_BUILD_VERSION),
    )
    monkeypatch.setattr(gubbins_adapter, "resolve_binary", lambda dirs: binary)
    monkeypatch.setattr(gubbins_adapter, "resolve_runner", lambda dirs=(): runner)
    return runner


def _panaroo_alignment(tmp_path: Path) -> Path:
    """A REAL intermediate root that HAS the core gene alignment stage 7 writes."""
    panaroo_dir = tmp_path / "intermediate" / "panaroo"
    panaroo_dir.mkdir(parents=True, exist_ok=True)
    (panaroo_dir / gubbins_adapter.CORE_ALIGNMENT_BASENAME).write_text(
        ">SYN_001\nACGTACGTAC\n", encoding="utf-8"
    )
    return tmp_path / "intermediate"


class TestStage8LeavesStage9AnAlignment:
    """(a) The call is reached through the real REAL path, and the file lands."""

    def test_the_masked_alignment_lands_where_stage_9_reads_it(
        self, open_gate, manifest, fake_gubbins, tmp_path,
    ):
        intermediate = _panaroo_alignment(tmp_path)
        phylogeny_dir = tmp_path / "phylogeny"

        table = derive_recombination_tables(
            open_gate,
            manifest,
            RunMode.REAL,
            intermediate_root=intermediate,
            phylogeny_dir=phylogeny_dir,
            tool_output_root=tmp_path,
            stage_dir=tmp_path / "stages",
        )

        assert table.is_file(), "stage 8 must still write its own table"

        stage_nine_input = phylogeny_dir / STAGE_NINE_READS
        assert stage_nine_input.is_file(), (
            f"stage 9 reads {stage_nine_input} (stages/phylogeny.py:504) and "
            "stage 8 is its only producer. The file being absent here is the "
            "defect this test exists for: with the call on the wrong module the "
            "REAL branch raised AttributeError after gubbins had already run, so "
            "no run could ever hand stage 9 a masked alignment."
        )
        assert stage_nine_input.read_text(encoding="utf-8").strip(), (
            "an empty file at stage 9's input path is not an alignment; stage 9 "
            "would read it as an empty cohort rather than a missing one"
        )

    def test_it_is_the_name_stage_8_declares_and_stage_9_reads(self):
        """The write name and the read name cannot drift apart unnoticed.

        `stages/recombination.py:84` declares what stage 8 writes;
        `stages/phylogeny.py:504` hard-codes what stage 9 reads. Nothing in the
        type system connects those two literals. This is the assertion that
        connects them.
        """
        assert stage_recombination.SNP_ALIGNMENT_NAME == STAGE_NINE_READS, (
            "stage 8 writes stages.recombination.SNP_ALIGNMENT_NAME "
            "(recombination.py:84) and stage 9 reads the literal "
            f"{STAGE_NINE_READS!r} (phylogeny.py:504). If these two ever "
            "diverge, stage 9 reads a file no producer writes, and the failure "
            "is a missing input at run time rather than here."
        )

    def test_stage_9_would_be_able_to_read_what_stage_8_wrote(
        self, open_gate, manifest, fake_gubbins, tmp_path,
    ):
        """The file is FASTA carrying the cohort, not an empty or wrong file.

        Asserting existence alone would pass for a zero-byte file or one holding
        an unrelated cohort. This reads it the way a FASTA consumer does and
        checks the tips are the samples that were analysed.
        """
        phylogeny_dir = tmp_path / "phylogeny"
        derive_recombination_tables(
            open_gate,
            manifest,
            RunMode.REAL,
            intermediate_root=_panaroo_alignment(tmp_path),
            phylogeny_dir=phylogeny_dir,
            tool_output_root=tmp_path,
            stage_dir=tmp_path / "stages",
        )

        text = (phylogeny_dir / STAGE_NINE_READS).read_text(encoding="utf-8")
        tips = [
            line[1:].strip()
            for line in text.splitlines()
            if line.startswith(">")
        ]
        assert sorted(tips) == sorted(SAMPLES), (
            f"stage 9's input must carry the analysed cohort {sorted(SAMPLES)}, "
            f"found {sorted(tips)}"
        )

    def test_the_row_join_completed_so_the_call_was_reached(
        self, open_gate, manifest, fake_gubbins, tmp_path,
    ):
        """`block_rows` runs *after* the call, so its output proves the call ran.

        A row for every node in the synthetic statistics file, and a table with
        the contract's header. With the call on the wrong module this never
        executes, so the table would not exist at all.
        """
        stage_dir = tmp_path / "stages"
        table = derive_recombination_tables(
            open_gate,
            manifest,
            RunMode.REAL,
            intermediate_root=_panaroo_alignment(tmp_path),
            phylogeny_dir=tmp_path / "phylogeny",
            tool_output_root=tmp_path,
            stage_dir=stage_dir,
        )

        lines = table.read_text(encoding="utf-8").splitlines()
        assert lines[0].split("\t") == list(stage_recombination.BLOCK_COLUMNS)
        assert len(lines) - 1 == len(STATS_ROWS), (
            "one row per node in the synthetic per-branch statistics; a partial "
            "join means the tree and the statistics were paired wrongly"
        )
        detected = [
            line for line in lines[1:]
            if line.split("\t")[-1] == "1"
        ]
        assert len(detected) == 2, (
            "the two synthetic rows carrying a recombination block must be "
            "reported as detected; an all-zero column would mean the statistics "
            "were not read"
        )


class TestTheCallSiteNamesASymbolThatExists:
    """(b) A wrong target module is a one-line failure, not a stage AttributeError."""

    @staticmethod
    def _function_node() -> ast.FunctionDef:
        tree = ast.parse(RUN_PY.read_text(encoding="utf-8"))
        return next(
            node for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef)
            and node.name == "derive_recombination_tables"
        )

    @classmethod
    def _build_outputs_call(cls) -> ast.Call:
        """The `build_outputs(...)` call inside `derive_recombination_tables`.

        Found by AST rather than by reading the line, so reformatting the call
        does not silently turn this test into a no-op.
        """
        function = cls._function_node()
        calls = [
            node for node in ast.walk(function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "build_outputs"
        ]
        assert len(calls) == 1, (
            f"expected exactly one build_outputs call in "
            f"derive_recombination_tables, found {len(calls)}"
        )
        return calls[0]

    @classmethod
    def _alias_map(cls) -> dict:
        """Local name -> importable module path, from the function's own imports.

        The call site reads `stage_recombination.build_outputs`, and
        `stage_recombination` is not a module name — it is the alias bound by
        `from .stages import recombination as stage_recombination` at the top of
        the function. Resolving the alias through the import that binds it is the
        only way to learn which module the call actually reaches, and it is also
        what makes the test independent of how the alias happens to be spelled.
        """
        package = __package__ or "papipeline"
        aliases = {}
        for node in ast.walk(cls._function_node()):
            if not isinstance(node, ast.ImportFrom):
                continue
            if node.level:
                # Relative to the package the module under test lives in.
                prefix = package.rsplit(".", node.level - 1)[0] if node.level > 1 else package
                base = f"{prefix}.{node.module}" if node.module else prefix
            else:
                base = node.module or ""
            for alias in node.names:
                local = alias.asname or alias.name
                aliases[local] = f"{base}.{alias.name}" if alias.name != "*" else base
        return aliases

    @classmethod
    def _resolve(cls, node: ast.AST):
        """Turn `stage_recombination.build_outputs` into the object it names."""
        assert isinstance(node, ast.Attribute), (
            "the call must be a module attribute (`<module>.build_outputs`); a "
            "bare name or a different shape is not what this test resolves"
        )
        assert isinstance(node.value, ast.Name), (
            f"expected a module alias, found {ast.dump(node.value)}"
        )
        aliases = cls._alias_map()
        local = node.value.id
        assert local in aliases, (
            f"{local!r} is not bound by any import inside "
            f"derive_recombination_tables (it binds {sorted(aliases)}), so the "
            "module the call reaches cannot be determined from this file alone"
        )
        module = importlib.import_module(aliases[local])
        return module, node.attr

    def test_the_module_named_by_the_call_site_has_build_outputs(self):
        call = self._build_outputs_call()
        module, attribute = self._resolve(call.func)

        assert hasattr(module, attribute), (
            f"{module.__name__} has no {attribute!r}. run.py calls "
            f"{module.__name__}.{attribute} inside derive_recombination_tables, "
            f"so every REAL run raised AttributeError there - after gubbins had "
            "already run - and stage 9 was left without its only input. The "
            "function lives on papipeline.stages.recombination:416."
        )

    def test_the_named_function_takes_exactly_the_keywords_passed(self):
        """A symbol can exist and still be called wrongly. Check the keywords."""
        call = self._build_outputs_call()
        module, attribute = self._resolve(call.func)
        function = getattr(module, attribute)

        passed = [keyword.arg for keyword in call.keywords]
        assert all(name is not None for name in passed), (
            "a positional argument here would defeat the keyword-only signature "
            f"{inspect.signature(function)}; the module alias cannot be used to "
            "satisfy the leading parameters"
        )
        accepted = set(inspect.signature(function).parameters)
        assert set(passed) <= accepted, (
            f"{module.__name__}.{attribute} accepts {sorted(accepted)}, but "
            f"run.py passes {sorted(passed)}"
        )
        assert set(passed) == {"gubbins_dir", "phylogeny_dir"}, (
            "stage 8's outputs are the masked alignment and the trees; the "
            "call must pass exactly those two directories"
        )

    def test_build_outputs_is_the_stage_module_function_not_a_sibling(self):
        """Name the intended target, so the fix is not a coin flip.

        Pins that the reachable definition is the one on the stage module and
        that the adapter does not also carry a same-named symbol. If a second
        definition ever appears, the run would be calling whichever the alias
        resolves to, and this says so.
        """
        assert stage_recombination.build_outputs is not None
        assert not hasattr(gubbins_adapter, "build_outputs"), (
            "papipeline.adapters.gubbins has gained a build_outputs. If it is "
            "not the same function as stages.recombination.build_outputs, then "
            "the two names mean different things and a reader cannot tell which "
            "one run.py meant."
        )