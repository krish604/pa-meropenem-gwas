"""The Linux runner scripts, the Linux machine overlay and the Linux docs.

These are the guarantees that hold on *any* machine, including this one, with
no Linux host, no Snakemake run and no Bakta invocation:

* every script under ``scripts/linux/`` parses under ``bash -n`` -- and the
  parser is the bash 3.2 macOS still ships, so a bash-4-only construct fails
  here rather than on the analysis machine;
* ``run_meropenem_real.sh`` **refuses** unless ``CONFIRM_REAL=yes``. It is
  executed in this test with an empty ``PATH`` on purpose: if the gate were
  ever moved below the first tool invocation, the command would fail with
  "command not found" instead of silently running the pipeline;
* the gate precedes ``PIPELINE_ALLOW_REAL_MODE=1`` and the snakemake call, and
  the overlay it names keeps ``runtime.allow_real_mode: false`` -- AGENTS.md
  rule 6, asserted on the file the runner actually loads;
* ``bootstrap_linux.sh`` creates an environment and records versions, and
  contains nothing that could start a pipeline;
* ``config/machines/linux.yaml`` loads as a machine overlay, declares the
  tools its stages invoke, and carries no machine-specific absolute path in
  ``paths`` (an absolute path belongs in an argument, not in a tracked file);
* the two documents exist, in order, and the flowchart is a mermaid block.

Nothing here reads ``data/``, ``db/`` or ``PDC_essential.tsv``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
import yaml

from papipeline.config.loader import load_machine_config

REPO = Path(__file__).resolve().parents[2]
LINUX_DIR = REPO / "scripts" / "linux"
OVERLAY = REPO / "config" / "machines" / "linux.yaml"
LINUX_ENV = REPO / "environment" / "environment-linux.yml"
LAPTOP_ENV = REPO / "environment" / "environment.yml"
LINUX_RUN_DOC = REPO / "docs" / "LINUX_RUN.md"
FLOWCHART_DOC = REPO / "docs" / "PIPELINE_FLOWCHART.md"

SCRIPTS = (
    LINUX_DIR / "bootstrap_linux.sh",
    LINUX_DIR / "run_bakta.sh",
    LINUX_DIR / "run_meropenem_real.sh",
)

#: The absolute bash, used so the PATH below cannot hide it.
BASH = "/bin/bash"

#: Benign commands the scripts may need *before* their gate: path handling,
#: dates, text filtering. Nothing here can start an analysis or a bioinformatics
#: tool, which is the point: a script that reached past its gate would find no
#: `bakta`, no `snakemake`, no `python3` and die with "command not found"
#: instead of doing the thing these tests exist to forbid.
SAFE_COMMANDS = (
    "dirname",
    "basename",
    "cat",
    "date",
    "grep",
    "sed",
    "tr",
    "awk",
)

_SANDBOX: Path | None = None


def _sandbox_bin() -> Path:
    """One directory of symlinks to SAFE_COMMANDS, built on first use."""
    global _SANDBOX
    if _SANDBOX is None:
        sandbox = Path(tempfile.mkdtemp(prefix="linux-scripts-test-bin-"))
        for name in SAFE_COMMANDS:
            found = shutil.which(name)
            if found:
                os.symlink(found, sandbox / name)
        _SANDBOX = sandbox
    return _SANDBOX


def _isolated_env(**overrides: str) -> dict:
    """The child's environment: no confirmation, and a PATH with no tools.

    An empty PATH would do it too, but then even `dirname` and `sed` would be
    missing and the scripts could not reach their own `--help`. The sandbox
    gives them exactly those and nothing else.
    """
    env = dict(os.environ)
    env.pop("CONFIRM_REAL", None)
    env["PATH"] = str(_sandbox_bin())
    env.update(overrides)
    return env


def _run(script: Path, *args: str, env: dict | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [BASH, str(script), *args],
        capture_output=True,
        text=True,
        cwd=REPO,
        env=env or _isolated_env(),
        timeout=60,
    )


class TestEveryLinuxScriptParses:
    @pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
    def test_bash_n_accepts_it(self, script: Path):
        assert script.is_file(), f"{script} is missing"
        result = subprocess.run(
            [BASH, "-n", str(script)],
            capture_output=True,
            text=True,
            cwd=REPO,
            env=_isolated_env(),
            timeout=60,
        )
        assert result.returncode == 0, (
            f"{script.name} does not parse under `bash -n` (the bash 3.2 on this "
            f"machine): {result.stderr.strip()}"
        )

    @pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
    def test_it_is_executable(self, script: Path):
        assert os.access(script, os.X_OK), (
            f"{script.name} is not executable; `docs/LINUX_RUN.md` tells the "
            "operator to run it as ./scripts/linux/..."
        )


class TestTheRealRunnerRefusesWithoutTheConfirmation:
    def test_no_confirmation_means_refusal(self):
        result = _run(LINUX_DIR / "run_meropenem_real.sh")
        assert result.returncode != 0, (
            "run_meropenem_real.sh ran without CONFIRM_REAL=yes. AGENTS.md rule 6 "
            "requires the operator to say `run real samples` first."
        )
        assert "CONFIRM_REAL" in (result.stderr + result.stdout), (
            f"the refusal must name the variable it wants: {result.stderr!r}"
        )

    @pytest.mark.parametrize("value", ["", "1", "true", "YES", "yes please"])
    def test_only_the_exact_word_opens_the_gate(self, value: str):
        result = _run(
            LINUX_DIR / "run_meropenem_real.sh", env=_isolated_env(CONFIRM_REAL=value)
        )
        assert result.returncode != 0, (
            f"CONFIRM_REAL={value!r} was accepted; only `yes` may open the gate"
        )
        assert "CONFIRM_REAL" in (result.stderr + result.stdout)

    def test_help_works_without_the_confirmation(self):
        result = _run(LINUX_DIR / "run_meropenem_real.sh", "--help")
        assert result.returncode == 0, result.stderr
        assert "CONFIRM_REAL" in (result.stdout + result.stderr)

    def test_the_gate_precedes_the_real_mode_variable_and_snakemake(self):
        """Static, because executing this path would start the pipeline.

        The order is the guarantee: the confirmation is checked before
        ``PIPELINE_ALLOW_REAL_MODE`` is set and before snakemake is named, so a
        run that skipped the gate could not have reached either.
        """
        text = (LINUX_DIR / "run_meropenem_real.sh").read_text(encoding="utf-8")
        gate = text.find("CONFIRM_REAL")
        allow = text.find("PIPELINE_ALLOW_REAL_MODE=1")
        snake = text.find("snakemake")
        overlay = text.find("config/machines/linux.yaml")
        assert gate != -1 and allow != -1 and snake != -1, (
            "the script must name CONFIRM_REAL, PIPELINE_ALLOW_REAL_MODE=1 and "
            "snakemake"
        )
        assert gate < allow < snake, (
            "the confirmation gate must come before PIPELINE_ALLOW_REAL_MODE=1 and "
            "before the snakemake invocation"
        )
        assert overlay != -1 and overlay < snake, (
            "the runner must hand snakemake machine=config/machines/linux.yaml"
        )

    def test_it_never_enables_real_mode_in_the_committed_overlay(self):
        raw = yaml.safe_load(OVERLAY.read_text(encoding="utf-8"))
        runtime = raw.get("runtime") or {}
        assert runtime.get("allow_real_mode") is False, (
            "config/machines/linux.yaml must keep runtime.allow_real_mode false; "
            "only PIPELINE_ALLOW_REAL_MODE=1 for one session opens it"
        )


class TestTheBaktaRunnerIsBatchOnly:
    def test_it_refuses_without_a_manifest(self):
        result = _run(LINUX_DIR / "run_bakta.sh")
        assert result.returncode != 0, "run_bakta.sh ran with no manifest"
        assert "--manifest" in (result.stderr + result.stdout), (
            f"the refusal must name the argument it needs: {result.stderr!r}"
        )

    def test_help_is_available(self):
        result = _run(LINUX_DIR / "run_bakta.sh", "--help")
        assert result.returncode == 0, result.stderr
        assert "--manifest" in result.stdout

    def test_it_names_the_output_layout_the_readme_expects(self):
        text = (LINUX_DIR / "run_bakta.sh").read_text(encoding="utf-8")
        assert "results" in text and "intermediate" in text and "bakta" in text, (
            "the runner must write where README.md says annotation lives: "
            "results/<mode>/intermediate/bakta/<sample_id>"
        )


class TestBootstrapRecordsVersionsAndRunsNothing:
    def test_help_is_available(self):
        result = _run(LINUX_DIR / "bootstrap_linux.sh", "--help")
        assert result.returncode == 0, result.stderr

    def test_it_calls_the_versions_script(self):
        text = (LINUX_DIR / "bootstrap_linux.sh").read_text(encoding="utf-8")
        assert "environment/versions.sh" in text, (
            "bootstrap must record tool and database versions after creating the "
            "environment; a run whose provenance is not recorded is an anecdote"
        )

    def test_it_cannot_start_a_pipeline(self):
        text = (LINUX_DIR / "bootstrap_linux.sh").read_text(encoding="utf-8")
        for forbidden in ("snakemake", "full_run", "run_meropenem_real", "--mode REAL"):
            assert forbidden not in text, (
                f"bootstrap_linux.sh mentions {forbidden!r}; it creates an "
                "environment and records versions, nothing else"
            )


class TestTheLinuxMachineOverlay:
    @pytest.fixture(scope="class")
    def machine(self):
        return load_machine_config(OVERLAY)

    def test_it_loads_and_names_itself(self, machine):
        assert machine.name == "linux"
        assert machine.root == REPO, (
            "a path-loaded overlay must root itself at the repository, not at the "
            "overlay's own directory"
        )

    def test_real_mode_is_closed(self, machine):
        assert machine.allow_real_mode is False

    def test_resources_are_declared(self, machine):
        assert machine.threads, "runtime.threads is required and has no default"
        assert machine.memory_mb, "runtime.memory_mb is required and has no default"

    def test_it_declares_the_tools_its_stages_invoke(self, machine):
        """stage 3 (mlst) and stage 4 (amrfinder + blast).

        Same invariant as tests/unit/test_overlay_tool_declarations.py, asserted
        here on the new overlay so a later edit cannot drop the block.
        """
        for tool in ("mlst", "amrfinder", "blast"):
            assert isinstance(machine.tool_available(tool), bool), (
                f"`{tool}` must be declared with a boolean `available`; an "
                "undeclared tool reads as unknown, not as missing"
            )

    def test_it_declares_no_machine_specific_absolute_path(self, machine):
        """Every `paths` entry is repository-relative.

        An absolute path belongs in an argument to the runner, not in a tracked
        file: the file would be correct only for the one machine its author was
        standing on.
        """
        offending = {
            key: value
            for key, value in (machine.paths or {}).items()
            if isinstance(value, str) and value.startswith("/")
        }
        assert not offending, (
            f"config/machines/linux.yaml hard-codes absolute paths: {offending}"
        )

    def test_it_declares_no_science(self):
        raw = yaml.safe_load(OVERLAY.read_text(encoding="utf-8"))
        for section in ("reference", "antibiotics", "cohort", "annotation"):
            assert section not in raw, (
                f"the linux overlay declares {section!r}; science lives in "
                "config/science.yaml so two machines cannot disagree about it"
            )


class TestTheLinuxEnvironmentFile:
    def test_it_parses_and_pins_the_same_python_as_the_laptop(self):
        linux = yaml.safe_load(LINUX_ENV.read_text(encoding="utf-8"))
        laptop = yaml.safe_load(LAPTOP_ENV.read_text(encoding="utf-8"))
        assert linux["name"] == laptop["name"], "both files must build one env name"
        assert linux["channels"] == ["conda-forge", "bioconda"]
        assert "python=3.11.16" in linux["dependencies"], (
            "the scientific stack is identical on both machines on purpose; a "
            "divergence between the laptop result and the analysis result is "
            "exactly the failure that invalidates a GWAS"
        )
        assert "python=3.11.16" in laptop["dependencies"]

    def test_every_pin_carries_an_explicit_version(self):
        linux = yaml.safe_load(LINUX_ENV.read_text(encoding="utf-8"))
        for dep in linux["dependencies"]:
            if not isinstance(dep, str):
                continue
            if dep.startswith("pip"):
                continue
            assert "=" in dep, (
                f"{dep!r} is unpinned; every version in this file must be a "
                "measured one, not a floating range"
            )

    def test_it_documents_the_dry_run_check_it_requires(self):
        text = LINUX_ENV.read_text(encoding="utf-8")
        assert "--dry-run" in text, (
            "the file's own contract is that every version is verified with a "
            "dry-run solve for the target platform; the command must be written "
            "where the next reader will find it"
        )


class TestTheLinuxDocuments:
    def test_linux_run_lists_the_steps_in_order(self):
        text = LINUX_RUN_DOC.read_text(encoding="utf-8")
        assert LINUX_RUN_DOC.is_file()
        positions = [
            text.find(name)
            for name in (
                "bootstrap_linux.sh",
                "run_bakta.sh",
                "run_meropenem_real.sh",
            )
        ]
        assert all(p != -1 for p in positions), (
            "docs/LINUX_RUN.md must name all three runner scripts"
        )
        assert positions == sorted(positions), (
            "the steps must be ordered: bootstrap, then Bakta, then the run"
        )

    def test_the_flowchart_is_mermaid_and_covers_the_stage_order(self):
        text = FLOWCHART_DOC.read_text(encoding="utf-8")
        assert "```mermaid" in text, "the flowchart must be a mermaid block"
        # 6a is the amendment (spec.md:363): a flowchart without it is the
        # old fifteen-stage taxonomy.
        for stage in ("validation", "cohort_variants", "gwas", "reporting"):
            assert stage in text, f"the flowchart is missing stage {stage!r}"
