"""`snakemake -n` must resolve the DAG, for both machine configurations.

The DAG is checked against a throwaway results root supplied via
``PIPELINE_RESULTS_ROOT``. That isolation is the point: run against the real
`results/`, a workflow with orphaned rules still passes, because Snakemake
finds the target file already on disk and never looks for a producer. That is
not hypothetical - it is how this file was green while `pangenome` and
`phenotype` were unreachable from `all`.

The first version of this test asserted only that the Snakefile loads. It now
also asserts that every stage rule is reachable from `all`, which is what the
spec's 15 stages actually mean.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV
from papipeline.run import STAGE_ORDER

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = PIPELINE_ROOT / "workflow" / "Snakefile"

pytestmark = pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake is not installed; this check cannot run here",
)


def _dry_run(
    machine: str,
    results_root: Path,
    mode: str = "TEST",
    quiet: bool = True,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--config", f"mode={mode}", f"machine={machine}",
            "--cores", "1",
            "--dry-run",
            *(["--quiet"] if quiet else []),
        ],
        cwd=PIPELINE_ROOT,
        # The redirect is an env var, not --config: the rules shell out to a
        # separate process that recomputes its paths from config, so a --config
        # value would reach this DAG and not the stage that writes the file.
        env={**os.environ, RESULTS_ROOT_ENV: str(results_root)},
        capture_output=True,
        text=True,
        timeout=600,
    )


#: The stages this codebase actually declares, taken from STAGE_ORDER.
#:
#: NOT the spec's canonical list, and deliberately not named "canonical": the
#: two disagree on membership, not just on count. The spec's D1 lists fifteen
#: stages including `variants` (minimap2/samtools vs PAO1), `recombination`
#: (gubbins) and `similarity`, none of which exist in the code. The code has
#: `structural_variants` in place of `variants`, and adds `mechanisms` and
#: `regulators`, which the spec's list does not mention. Reconciling that is a
#: design decision, not a rename, so this constant records what exists.
DECLARED_STAGES = frozenset(STAGE_ORDER)


def _scheduled(result: subprocess.CompletedProcess) -> set:
    """Rule names Snakemake actually scheduled for this dry run."""
    combined = result.stdout + result.stderr
    return {m.group(1) for m in re.finditer(r"(?m)^\s*rule (\w+):", combined)}


@pytest.mark.parametrize("machine", ["laptop", "bigmachine"])
def test_the_dag_resolves_for_a_machine(machine: str, tmp_path: Path):
    """The Snakefile must load and its DAG must build from nothing on disk."""
    result = _dry_run(machine, tmp_path / machine)
    assert result.returncode == 0, (
        f"snakemake --dry-run failed for machine={machine}\n"
        f"--- stdout ---\n{result.stdout[-3000:]}\n"
        f"--- stderr ---\n{result.stderr[-3000:]}"
    )


@pytest.mark.parametrize("machine", ["laptop", "bigmachine"])
def test_every_stage_is_reachable_from_all(machine: str, tmp_path: Path):
    """No stage may be an orphan.

    A rule that nothing depends on is never scheduled, so a green dry run says
    nothing about it. This is the check that catches it.
    """
    result = _dry_run(machine, tmp_path / machine, quiet=False)
    combined = result.stdout + result.stderr
    if result.returncode != 0:
        pytest.fail(f"DAG did not resolve; cannot assess reachability:\n{combined[-3000:]}")

    scheduled = _scheduled(result)
    orphans = sorted(DECLARED_STAGES - scheduled)
    assert not orphans, (
        f"these stages are never scheduled from `all` for machine={machine}: {orphans}. "
        "A stage nothing depends on does not run, and a dry run will not say so."
    )


def test_a_rule_may_not_declare_two_execution_keywords(tmp_path: Path):
    """The original Phase 2 defect, kept so it cannot come back.

    Snakemake allows one of run/shell/script/notebook/wrapper/
    template_engine/cwl per rule.
    """
    result = _dry_run("laptop", tmp_path / "exec")
    combined = result.stdout + result.stderr
    if "Multiple" in combined and "keywords" in combined:
        offenders = sorted({
            line.strip()
            for line in combined.splitlines()
            if "Multiple" in line or "keywords in rule" in line
        })
        pytest.fail(
            "Every rule must declare exactly one execution keyword; the "
            f"Snakefile declares both `script:` and `shell:` in: {offenders}"
        )
    assert result.returncode == 0, combined[-3000:]


# --------------------------------------------------------------------------
# A missing upstream file must fail the DAG build, not a run
# --------------------------------------------------------------------------

#: An external input that TEST declares and that no rule in the workflow
#: produces. Removing it is the cleanest way to provoke a MissingInput, because
#: the file exists in the committed fixtures and nothing can regenerate it.
PROBE = Path("intermediate") / "amr" / "amr_determinants.tsv"


def _temp_root_with_input_removed(tmp_path: Path) -> Path:
    """A throwaway overlay whose test_data is a copy missing one input.

    Redirects via a machine overlay *path*, which `load_config` already accepts,
    rather than a new seam: `--config machine=<path>` is an existing feature.
    The committed fixtures are copied, never edited, so `data/` and
    `test_data/` stay untouched.
    """
    import shutil

    machine_dir = tmp_path / "machine"
    machine_dir.mkdir()
    copied_data = tmp_path / "test_data"
    shutil.copytree(PIPELINE_ROOT / "test_data", copied_data)
    (copied_data / PROBE).unlink()

    source = (PIPELINE_ROOT / "config" / "machines" / "laptop.yaml").read_text(
        encoding="utf-8"
    )
    patched = source.replace(
        "test_data_root: test_data", f"test_data_root: {copied_data}"
    )
    assert patched != source, "laptop.yaml no longer declares test_data_root"
    overlay = machine_dir / "laptop.yaml"
    overlay.write_text(patched, encoding="utf-8")
    return overlay


def _dry_run_with_overlay(overlay: Path, results_root: Path, mode: str):
    return subprocess.run(
        [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--config", f"mode={mode}", f"machine={overlay}",
            "--cores", "1", "--dry-run",
        ],
        cwd=PIPELINE_ROOT,
        env={**os.environ, RESULTS_ROOT_ENV: str(results_root)},
        capture_output=True, text=True, timeout=600,
    )


def test_a_missing_external_input_fails_the_dag_build(tmp_path: Path):
    """An absent upstream file must stop the workflow before anything runs.

    This is the property `stage_inputs` exists to provide: TEST and REAL
    declare their external inputs, so an upstream tool that produced nothing
    stops the run immediately. Without it the stage would fail mid-run, hours
    in, after the expensive stages had already been paid for.
    """
    overlay = _temp_root_with_input_removed(tmp_path)
    result = _dry_run_with_overlay(overlay, tmp_path / "results", mode="TEST")
    combined = result.stdout + result.stderr

    assert result.returncode != 0, (
        "removing a required external input still resolved the DAG, so a "
        "missing upstream tool result would not be caught until a stage ran"
    )
    assert "MissingInput" in combined or "MissingOutput" in combined, (
        f"expected a DAG-build failure, got:\n{combined[-3000:]}"
    )
    assert "amr_determinants.tsv" in combined, (
        f"the failure must name the file that is missing:\n{combined[-3000:]}"
    )


def test_real_mode_is_refused_before_the_dag_is_even_built(tmp_path: Path):
    """REAL is gated, and the gate is the *first* thing that refuses.

    Deliberately asserted rather than worked around: exercising REAL's
    DAG-build path would mean setting `runtime.allow_real_mode`, which is the
    user's decision and not a test's. So the reachable behaviour is checked -
    the gate fires, and it names REAL and the flag - and the missing-input
    path above is covered in TEST, where the same `stage_inputs` call decides.
    """
    overlay = _temp_root_with_input_removed(tmp_path)
    result = _dry_run_with_overlay(overlay, tmp_path / "results", mode="REAL")
    combined = result.stdout + result.stderr

    assert result.returncode != 0
    assert "REAL" in combined, f"the refusal must name REAL:\n{combined[-3000:]}"
    assert "allow_real_mode" in combined, (
        f"the refusal must name the flag that would change it:\n{combined[-3000:]}"
    )
    assert "MissingInput" not in combined, (
        "the mode gate should fire before any input is checked"
    )
