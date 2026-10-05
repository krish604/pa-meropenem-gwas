"""`run_pipeline` dispatches stage 8, and refuses REAL for two named reasons.

**What this file is.** `papipeline/stages/recombination.py` was written and
committed with 448 lines, a REAL refusal and a self-test, and it had *no
caller*. `recombination` sat in `run.UNBUILT_STAGES`, and
`config/science.yaml` had no `analysis.recombination` key, so
`PipelineConfig.stage_enabled("recombination")` was False, `run.py`'s
`enabled()` returned False, the stage logged at INFO and was skipped - and
stage 9 then built its tree from whatever alignment it found, silently
unmasked. Two descriptions of one stage, both saying "absent", and the
pipeline reported success.

So the reconciliation is three changes that only make sense together:

* `run.derive_recombination_tables` dispatches the stage;
* `recombination` leaves `run.UNBUILT_STAGES` and
  `standins.UNBUILT_WITH_STANDIN`;
* `analysis.recombination: true` is in `config/science.yaml`.

**What this file asserts, and why each assertion is a behaviour.** The two
REAL refusals are separate faults and are tested separately, because they have
separate remedies and an operator told the wrong one rebuilds a binary that is
fine:

1. `runtime.allow_real_mode` is false - the gate on the run. Raised before
   anything else, so an operator whose gate is shut is not also told their
   alignment is missing.
2. The resolved gubbins binary exits 139 (SIGSEGV) on its own probe - the
   broken conda ``osx-arm64`` build. Raised by the adapter's
   `refuse_segfaulting_binary`, which names 139 and both versions.

**No REAL run happens here.** `gubbins` cannot be installed on this machine, so
the REAL tests inject the adapter's own `selftest` hook - the seam
`run_gubbins(selftest=...)` exists for, and the only reason it exists: there is
no way to make a good binary segfault on demand, so a refusal that could only be
tested by shipping a broken binary would be a refusal nobody tests.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from papipeline.adapters import gubbins as gubbins_adapter
from papipeline.config.loader import load_config
from papipeline.errors import (
    ModeNotAllowedError,
    StageError,
    ToolExecutionError,
    ToolNotAvailableError,
)
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode
from papipeline.run import (
    STAGE_ORDER,
    UNBUILT_STAGES,
    derive_recombination_tables,
)
from papipeline.stages import recombination as stage_recombination
from papipeline.stages.recombination import BLOCK_COLUMNS

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"

#: The committed gubbins-shaped output TEST reads. Not a stand-in: real output
#: of the real parser's input format, written once by a deterministic script so
#: TEST exercises stage 8 with no external tool.
TEST_GUBBINS_DIR = REPO / "test_data" / "intermediate" / "gubbins"


@pytest.fixture
def manifest() -> SampleManifest:
    return SampleManifest(
        [SimpleNamespace(sample_id=f"TEST_PA_{i:03d}", assembly_path=None, source="test")
         for i in range(1, 21)]
    )


@pytest.fixture
def config():
    return load_config(SCIENCE, machine="laptop")


def _open_gate(config):
    """The same config with `runtime.allow_real_mode` true.

    Built by `dataclasses.replace` on the runtime mapping rather than by
    re-reading an overlay: the committed overlays all keep the gate shut, and
    opening one in a test file is not a decision a test gets to make.
    """
    from dataclasses import replace

    return replace(config, runtime={**config.runtime, "allow_real_mode": True})


# ---------------------------------------------------------------------------
# 1. The stage is built, and the taxonomy says so.
# ---------------------------------------------------------------------------


class TestTheStageIsBuilt:
    def test_it_is_in_the_stage_order(self):
        assert "recombination" in STAGE_ORDER

    def test_it_is_not_unbuilt(self):
        assert "recombination" not in UNBUILT_STAGES, (
            "declared unbuilt means `enabled()` refuses it outside STUB, so "
            "analysis.recombination: true would raise NotImplementedError in TEST "
            "instead of running the stage"
        )

    def test_the_config_key_is_present_and_true(self, config):
        assert config.stage_enabled("recombination") is True, (
            "the dispatch branch exists, so the flag must enable it. Without the "
            "key the stage is skipped in silence - the exact defect this "
            "reconciliation closes."
        )


# ---------------------------------------------------------------------------
# 2. TEST: the real parser over the committed fixture, no external tool.
# ---------------------------------------------------------------------------


class TestTestModeReadsTheCommittedFixture:
    def test_it_writes_the_declared_contract_table(self, config, manifest, tmp_path):
        table = derive_recombination_tables(
            config,
            manifest,
            RunMode.TEST,
            intermediate_root=tmp_path / "intermediate",
            phylogeny_dir=tmp_path / "phylogeny",
            tool_output_root=REPO / "test_data" / "intermediate",
            stage_dir=tmp_path,
        )
        assert table.is_file()
        header = table.read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header == list(BLOCK_COLUMNS), (
            "the header must be the contract's, column for column; "
            "execution/contracts.py checks it and stages/recombination.py "
            "re-declares it so the two cannot drift"
        )

    def test_it_measures_the_fixture_rather_than_inventing_rows(self, config, manifest, tmp_path):
        """Rows come from the fixture, so a row count proves the parser ran."""
        table = derive_recombination_tables(
            config,
            manifest,
            RunMode.TEST,
            intermediate_root=tmp_path / "intermediate",
            phylogeny_dir=tmp_path / "phylogeny",
            tool_output_root=REPO / "test_data" / "intermediate",
            stage_dir=tmp_path,
        )
        statistics = gubbins_adapter.parse_per_branch_statistics(
            TEST_GUBBINS_DIR
            / f"{gubbins_adapter.PREFIX}"
              f"{gubbins_adapter.PER_BRANCH_STATISTICS_SUFFIX}"
        )
        written = table.read_text(encoding="utf-8").splitlines()[1:]
        assert len(written) == len(statistics), (
            "one row per node in per_branch_statistics.csv. A fixture read that "
            "dropped or invented nodes would change the count."
        )

    def test_it_does_not_write_into_the_committed_phylogeny_fixture(self, config, manifest, tmp_path):
        """TEST must not replace a fixture with this run's output.

        `phylogeny_dir` under TEST is the committed `test_data/phylogeny/`.
        Copying the masked alignment there in TEST would overwrite the fixture
        the phylogeny stage reads, which is a different kind of silent damage:
        it would pass locally and fail for the next reader.
        """
        from papipeline.config.loader import REPO_ROOT

        committed = REPO_ROOT / "test_data" / "phylogeny"
        before = {p.name: p.read_bytes() for p in sorted(committed.glob("*"))}
        derive_recombination_tables(
            config,
            manifest,
            RunMode.TEST,
            intermediate_root=tmp_path / "intermediate",
            # Deliberately the committed directory, not a tmp_path one: the
            # point is that this code path does not write here.
            phylogeny_dir=committed,
            tool_output_root=REPO / "test_data" / "intermediate",
            stage_dir=tmp_path,
        )
        after = {p.name: p.read_bytes() for p in sorted(committed.glob("*"))}
        assert after == before, (
            f"stage 8 in TEST modified the committed phylogeny fixture: "
            f"{sorted(set(after) ^ set(before)) or 'changed contents'}"
        )

    def test_it_refuses_when_the_fixture_directory_is_absent(self, config, manifest, tmp_path):
        """A missing input is a refusal, not an empty table."""
        from papipeline.errors import DataContractError

        with pytest.raises(DataContractError):
            derive_recombination_tables(
                config,
                manifest,
                RunMode.TEST,
                intermediate_root=tmp_path / "intermediate",
                phylogeny_dir=tmp_path / "phylogeny",
                tool_output_root=tmp_path / "no_such_tool_output",
                stage_dir=tmp_path,
            )


# ---------------------------------------------------------------------------
# 3. REAL refusal 1: the gate on the run.
# ---------------------------------------------------------------------------


class TestRealRefusesWhileTheGateIsShut:
    def test_it_names_the_flag_the_override_and_the_overlay(self, config, manifest, tmp_path):
        with pytest.raises(ModeNotAllowedError) as caught:
            derive_recombination_tables(
                config,
                manifest,
                RunMode.REAL,
                intermediate_root=tmp_path / "intermediate",
                phylogeny_dir=tmp_path / "phylogeny",
                tool_output_root=tmp_path,
                stage_dir=tmp_path,
            )
        message = str(caught.value)
        assert "allow_real_mode" in message
        assert "PIPELINE_ALLOW_REAL_MODE" in message, (
            "a one-session override has to be named, or the only remedy offered "
            "is editing a tracked overlay"
        )
        assert "laptop" in message, (
            "the overlay that currently holds the gate shut, by name: an "
            "operator with three machines needs to know which one"
        )

    def test_the_gate_fires_before_the_alignment_is_complained_about(
        self, config, manifest, tmp_path,
    ):
        """Ordering, which is the whole point of checking the gate first.

        `alignment` does not exist in this tmp_path. If the alignment check ran
        first the operator would be told to run panaroo - a long, expensive
        command - and the gate would still shut the stage when they came back.
        """
        with pytest.raises(ModeNotAllowedError):
            derive_recombination_tables(
                config,
                manifest,
                RunMode.REAL,
                intermediate_root=tmp_path / "absent",
                phylogeny_dir=tmp_path / "phylogeny",
                tool_output_root=tmp_path,
                stage_dir=tmp_path,
            )

    def test_the_stages_own_entry_point_gates_the_same_way(self, config):
        """`stages.recombination.run` is the second of the two checks.

        Not the dispatch path - `derive_recombination_tables` drives
        `run_gubbins` in REAL - but it is the entry point a caller reaches
        directly, so it must not be the weaker of the two. Asserted because the
        two were deliberately given the same shape when the stage was promoted.
        """
        with pytest.raises(NotImplementedError) as caught:
            stage_recombination.run(
                config,
                SampleManifest(()),
                RunMode.REAL,
                gubbins_dir=TEST_GUBBINS_DIR,
            )
        message = str(caught.value)
        assert "allow_real_mode" in message
        assert "PIPELINE_ALLOW_REAL_MODE" in message


# ---------------------------------------------------------------------------
# 4. REAL refusal 2: the resolved binary fails its self-test (exit 139).
# ---------------------------------------------------------------------------


class TestRealRefusesASegfaultingBinary:
    """The conda ``osx-arm64`` build, verified on this machine to exit 139."""

    @pytest.fixture
    def real_inputs(self, tmp_path):
        """A REAL intermediate root that HAS a panaroo alignment."""
        panaroo_dir = tmp_path / "intermediate" / "panaroo"
        panaroo_dir.mkdir(parents=True)
        (panaroo_dir / gubbins_adapter.CORE_ALIGNMENT_BASENAME).write_text(
            ">TEST_PA_001\nACGT\n", encoding="utf-8"
        )
        return tmp_path

    def test_exit_139_is_refused_by_name(self, config, manifest, real_inputs, tmp_path, monkeypatch):
        """The refusal names 139, the binary, and both versions.

        Asserted on the numbers rather than on the prose: the message is what
        tells an operator this is a packaging defect and not a bad cohort, and
        the numbers are what makes it checkable. `SOURCE_BUILD_VERSION` and
        `CONDA_BUILD_VERSION` are compared because the whole diagnosis is "the
        binary self-reports the wrong version" - if either moved without the
        other, the message would stop being true.
        """
        monkeypatch.setattr(gubbins_adapter, "resolve_binary", lambda dirs: Path("/opt/envs/gubbins/bin/gubbins"))
        monkeypatch.setattr(gubbins_adapter, "resolve_runner", lambda dirs=(): Path("/opt/envs/gubbins/bin/run_gubbins.py"))
        monkeypatch.setattr(
            gubbins_adapter, "self_test",
            lambda binary, scratch, search_dirs=(): (
                gubbins_adapter.SEGFAULT_EXIT, gubbins_adapter.CONDA_BUILD_VERSION,
            ),
        )

        with pytest.raises(ToolNotAvailableError) as caught:
            derive_recombination_tables(
                _open_gate(config),
                manifest,
                RunMode.REAL,
                intermediate_root=real_inputs / "intermediate",
                phylogeny_dir=real_inputs / "phylogeny",
                tool_output_root=tmp_path,
                stage_dir=tmp_path / "stages",
            )

        message = str(caught.value)
        assert str(gubbins_adapter.SEGFAULT_EXIT) in message, (
            "the refusal must state the exit code, so an operator can "
            "distinguish this from a usage error"
        )
        assert gubbins_adapter.CONDA_BUILD_VERSION in message
        assert gubbins_adapter.SOURCE_BUILD_VERSION in message
        assert "build_gubbins_from_source.sh" in message, (
            "the remedy is a specific script, and a refusal that only says "
            "'gubbins is broken' sends the reader to rebuild a binary that is "
            "fine"
        )

    def test_the_segfault_check_is_not_folded_into_the_generic_one(
        self, config, manifest, real_inputs, tmp_path, monkeypatch,
    ):
        """139 gets the segfault message; any other non-zero exit does not.

        `run_gubbins` folds a non-139 failure into a *different* error on
        purpose: a missing tree builder and a broken dylib are different problems,
        and telling an operator their tool is mis-packaged when it is merely
        mis-invoked sends them to rebuild a binary that is fine. This is the
        test that says the two were not merged.
        """
        for name, patch in (
            ("self_test", lambda binary, scratch, search_dirs=(): (2, "3.4.2")),
        ):
            monkeypatch.setattr(gubbins_adapter, name, patch)
        monkeypatch.setattr(gubbins_adapter, "resolve_binary", lambda dirs: Path("/opt/envs/gubbins/bin/gubbins"))
        monkeypatch.setattr(gubbins_adapter, "resolve_runner", lambda dirs=(): Path("/opt/envs/gubbins/bin/run_gubbins.py"))

        with pytest.raises(ToolExecutionError) as caught:
            derive_recombination_tables(
                _open_gate(config),
                manifest,
                RunMode.REAL,
                intermediate_root=real_inputs / "intermediate",
                phylogeny_dir=real_inputs / "phylogeny",
                tool_output_root=tmp_path,
                stage_dir=tmp_path / "stages",
            )
        message = str(caught.value)
        assert "not the exit-" in message, (
            "a non-139 failure must say so explicitly; the segfault remedy "
            "will not fix a missing tree builder"
        )


# ---------------------------------------------------------------------------
# 5. REAL refusal 3: the panaroo alignment is absent. Named, actionable.
# ---------------------------------------------------------------------------


class TestRealNamesThePanarooCommandItNeeds:
    def test_the_refusal_names_the_command_and_both_output_paths(
        self, config, manifest, tmp_path,
    ):
        """The pipeline never runs panaroo, so the operator has to be told how.

        Asserted on the pieces an operator needs: the tool, the flags that are
        actually required (not cosmetic), the output directory, and the
        presence-table basename stage 7 reads from it - because one panaroo run
        produces both and being told about only one sends them round again.
        """
        intermediate = tmp_path / "intermediate"
        with pytest.raises(StageError) as caught:
            derive_recombination_tables(
                _open_gate(config),
                manifest,
                RunMode.REAL,
                intermediate_root=intermediate,
                phylogeny_dir=tmp_path / "phylogeny",
                tool_output_root=tmp_path,
                stage_dir=tmp_path / "stages",
            )
        message = str(caught.value)
        assert "panaroo -i" in message, (
            "the refusal must carry the invocation, not just the tool name. The "
            f"command is the one verified in docs/environment-arm64.md section "
            f"3. Got: {message}"
        )
        assert "--remove-invalid-genes" in message, (
            "--remove-invalid-genes is REQUIRED: panaroo raises ValueError on "
            "any CDS whose length is not a multiple of 3, is under 34 bp, or "
            "contains an internal stop. On real Bakta data that is 10 genes out "
            "of 63,359. A command without it fails."
        )
        assert gubbins_adapter.CORE_ALIGNMENT_BASENAME in message
        from papipeline.adapters.panaroo import PRESENCE_CSV

        assert PRESENCE_CSV in message, (
            "stage 7 reads the presence table from the same directory; naming "
            "both means one command unblocks two stages"
        )
        assert str(intermediate / "panaroo") in message, (
            "the output directory must be the resolved path the parser reads, "
            "not a placeholder, so an operator can check it"
        )
        assert "not an absence of recombination" in message, (
            "a missing input and an absent finding are different, and only one "
            "of them is what a missing file means"
        )

    def test_the_documented_command_is_the_one_quoted(self):
        """A refusal that quotes a command nobody has run is a guess.

        AGENTS.md rule 1. `PROVISION_COMMAND` is asserted against the invocation
        recorded in `docs/environment-arm64.md` section 3, so a future edit to
        either has to be made deliberately.
        """
        from papipeline.adapters.panaroo import PROVISION_COMMAND

        doc = (REPO / "docs" / "environment-arm64.md").read_text(encoding="utf-8")
        assert (
            "panaroo -i *.gff3 -o out -t 4 --clean-mode moderate "
            "--remove-invalid-genes"
        ) in doc, (
            "docs/environment-arm64.md section 3 no longer records the "
            "verified invocation; PROVISION_COMMAND's provenance comment needs "
            "updating before this test can mean anything"
        )
        for flag in ("-i", "-o", "-t", "--clean-mode", "--remove-invalid-genes"):
            assert flag in PROVISION_COMMAND, (
                f"{flag} is in the verified invocation but not in the refusal. "
                "AGENTS.md rule 1: a wrong or missing flag is worse than no code."
            )
        assert "panaroo build" not in PROVISION_COMMAND, (
            "panaroo's verified invocation here has no subcommand word. "
            "`panaroo build` was written from memory and then corrected; it is "
            "recorded so nobody re-adds it."
        )

    def test_the_pangenome_adapter_never_invokes_panaroo(self):
        """The refusal's premise, asserted rather than assumed.

        `adapters/panaroo.py` imports no `subprocess`. If a future change adds
        an invocation, the refusal above stops being true and the operator is
        told to run a command the pipeline already ran. Checked from the AST so a
        mention of the word in a docstring or a comment cannot satisfy it.
        """
        import ast

        tree = ast.parse((REPO / "papipeline" / "adapters" / "panaroo.py").read_text(encoding="utf-8"))
        imported = {
            alias.name.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in getattr(node, "names", [])
        }
        assert "subprocess" not in imported, (
            "adapters/panaroo.py now imports subprocess. Either the stage-8 "
            "refusal is wrong or the pipeline runs panaroo, and the two cannot "
            "both be true."
        )