"""Tool discovery must not depend on hard-coded machine paths.

Phase 2, item 1. Three call sites fell back to a literal
``/opt/homebrew/bin/<tool>`` when a tool was not on ``PATH``. On the analysis
machine that directory does not exist, so the fallback silently resolved to
nothing and a *missing* tool was reported as present - the exact failure mode
that makes a pre-flight check worthless.

The version-awareness these call sites genuinely need (an ancient samtools on
``PATH`` shadowing a working one) is preserved; only the hard-coded locations
go away. Extra locations are configuration.

Written before the fix.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

from papipeline.config.loader import load_config, load_machine_config
from papipeline.errors import ConfigError

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"
LAPTOP = PIPELINE_ROOT / "config" / "machines" / "laptop.yaml"

HARDCODED = ("/opt/homebrew", "/usr/local/bin", "/opt/local", "C:\\\\")


# ---------------------------------------------------------------------------
# the repository's own source must not name a machine path
# ---------------------------------------------------------------------------


def test_no_script_hard_codes_a_machine_bin_directory():
    """The literal that caused the bug must not reappear anywhere."""
    offenders = []
    for path in sorted(PIPELINE_ROOT.glob("scripts/**/*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for lineno, line in enumerate(text.splitlines(), start=1):
            for needle in HARDCODED:
                if needle in line and "homebrew" in line.lower():
                    offenders.append(f"{path.relative_to(PIPELINE_ROOT)}:{lineno}")
    assert not offenders, f"hard-coded machine paths remain: {offenders}"


# ---------------------------------------------------------------------------
# candidates come from PATH plus configuration, and nothing else
# ---------------------------------------------------------------------------


def test_candidates_are_path_first_and_unique():
    laptop = load_machine_config(LAPTOP)
    candidates = laptop.tool_candidates("minimap2")
    assert candidates, "at least the PATH result should be present when installed"
    assert len(candidates) == len(set(candidates)), "candidates must be de-duplicated"
    found = [c for c in candidates if c == shutil.which("minimap2")]
    assert found, "the PATH result must be among the candidates"


def test_configured_search_dirs_are_appended_after_path(monkeypatch, tmp_path):
    laptop = load_machine_config(LAPTOP)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    real_dir = tmp_path / "tools"
    real_dir.mkdir()
    # the candidate must actually exist: a declared directory whose tool is
    # absent is correctly skipped, which is a different behaviour entirely
    (real_dir / "minimap2").write_text("#!/bin/sh\n", encoding="utf-8")
    laptop_with_dir = load_machine_config(
        _overlay_with_search_dirs(LAPTOP, [str(real_dir)])
    )
    assert laptop_with_dir.tool_search_dirs == (str(real_dir),)
    assert laptop_with_dir.tool_candidates("minimap2") == (
        str(real_dir / "minimap2"),
    )
    # The shipped overlays DO declare a search directory, and there are two
    # independent reasons for one on the laptop. mummer is homebrew-only on
    # osx-arm64 (`mummer4` does not solve - libcxx >=18 is absent). gubbins is
    # declared because conda's osx-arm64 build is broken (it links no zlib and
    # SIGSEGVs) and scripts/build_gubbins_from_source.sh installs a working one
    # under tools/gubbins/, which is not on PATH. Hard rule 2 forbids burying
    # either path in code, so they are declared here. What must hold is that the
    # declarations are deliberate and their tools are real - not that the list
    # is empty, which was true only before stage 7 gained a caller.
    #
    # Compared as a SET, not a tuple: this is an exact-equality tripwire on the
    # declaration itself, so it should fail when a directory is added or
    # removed, but not merely because YAML lists them in a different order.
    #
    # gubbins is asserted RELATIVE. An absolute entry would encode the absolute
    # location of one checkout and be wrong in every other worktree or clone.
    assert set(laptop.tool_search_dirs) == {
        "/opt/homebrew/opt/mummer/bin",
        "tools/gubbins/bin",
    }, f"tool_search_dirs changed: {list(laptop.tool_search_dirs)}"
    assert Path("/opt/homebrew/opt/mummer/bin").is_dir(), (
        "the declared mummer directory does not exist; install mummer or correct "
        "tool_search_dirs in the overlay"
    )


def test_a_relative_search_dir_resolves_against_the_repo_root_not_the_cwd():
    """The bug this guards, stated as the test that would have caught it.

    `tool_candidates` used to join a relative entry against the process working
    directory. That makes discovery depend on where you are standing: it works
    under pytest (which runs from the repo root) and resolves to nothing from a
    Snakemake run directory or a git worktree checked out elsewhere. Nothing
    raises - the tool is simply reported absent, which reads as "this machine
    does not have gubbins" rather than "this overlay points somewhere wrong".

    So run from a deliberately foreign CWD and require the same answer.
    """
    laptop = load_machine_config(LAPTOP)
    elsewhere = Path(tempfile.mkdtemp())
    original = Path.cwd()
    os.chdir(elsewhere)
    try:
        from_repo_root = laptop.tool_candidates("gubbins")
    finally:
        os.chdir(original)

    expected = str(PIPELINE_ROOT / "tools" / "gubbins" / "bin" / "gubbins")
    # Only assert the gubbins entry when the build artifact actually exists;
    # tools/ is gitignored, so a fresh clone legitimately has no binary and
    # tool_candidates correctly skips a directory whose tool is absent.
    if Path(expected).exists():
        assert expected in from_repo_root, (
            f"a repo-relative tool_search_dirs entry was not resolved against the "
            f"repo root; from a foreign CWD it produced {from_repo_root}"
        )
    assert laptop.resolved_tool_search_dirs()[-1] == str(
        PIPELINE_ROOT / "tools" / "gubbins" / "bin"
    ), "resolved_tool_search_dirs must anchor a relative entry to the repo root"


def test_an_absolute_search_dir_is_not_reinterpreted():
    """A Homebrew prefix is a machine fact and must be passed through as-is.

    Resolving is for repository-relative entries only. Rewriting an absolute
    prefix against the repo root would silently produce a directory that does
    not exist, turning a correct declaration into a missing tool.
    """
    laptop = load_machine_config(LAPTOP)
    resolved = laptop.resolved_tool_search_dirs()
    assert "/opt/homebrew/opt/mummer/bin" in resolved, (
        f"the absolute mummer prefix was rewritten; resolved: {list(resolved)}"
    )


def test_a_declared_search_dir_that_does_not_exist_is_skipped(tmp_path):
    """A configured location is still subject to being real."""
    overlay = _overlay_with_search_dirs(LAPTOP, ["/definitely/not/here"])
    laptop = load_machine_config(overlay)
    import shutil as _shutil

    original = _shutil.which
    _shutil.which = lambda name: None
    try:
        assert laptop.tool_candidates("minimap2") == ()
    finally:
        _shutil.which = original


def test_an_undeclared_tool_yields_only_its_path_result(monkeypatch):
    """No search-dir guessing for a tool nobody configured a location for."""
    laptop = load_machine_config(LAPTOP)
    # The shipped overlays DO declare a search directory, for the two reasons given
    # in `test_configured_search_dirs_are_appended_after_path` above: mummer is
    # homebrew-only on osx-arm64, and gubbins is source-built into tools/gubbins/
    # because the conda package is broken there.
    assert set(laptop.tool_search_dirs) == {
        "/opt/homebrew/opt/mummer/bin",
        "tools/gubbins/bin",
    }, f"tool_search_dirs changed: {list(laptop.tool_search_dirs)}"
    assert Path("/opt/homebrew/opt/mummer/bin").is_dir(), (
        "the declared mummer directory does not exist; install mummer or correct "
        "tool_search_dirs in the overlay"
    )
    monkeypatch.setattr(shutil, "which", lambda name: None)
    assert laptop.tool_candidates("never_installed") == ()


# ---------------------------------------------------------------------------
# wiring
# ---------------------------------------------------------------------------


def test_a_composed_config_exposes_the_candidates():
    config = load_config(SCIENCE, machine="laptop")
    assert isinstance(config.machine, object)
    assert config.tool_candidates("minimap2") == config.machine.tool_candidates("minimap2")


def test_a_self_contained_config_has_no_search_dirs():
    config = load_config(SCIENCE, machine=None)
    assert config.machine is None
    assert config.tool_candidates("minimap2") == ()


# ---------------------------------------------------------------------------
# the scripts must actually use the helper rather than their own list
# ---------------------------------------------------------------------------


def _overlay_with_search_dirs(source: Path, dirs) -> Path:
    import shutil as _shutil
    import tempfile

    tmp = Path(tempfile.mkdtemp()) / "machines" / "with-dirs.yaml"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    _shutil.copy(source, tmp)
    text = tmp.read_text(encoding="utf-8")
    # Drop any existing declaration first. Inserting a second
    # `tool_search_dirs:` key made a duplicate YAML key, and the overlay's own
    # value won - so the directory under test never appeared at all. The helper
    # has to *replace* the declaration, not add to it.
    text = re.sub(
        r"(?:  tool_search_dirs:\n(?:    - .*\n)*)",
        "",
        text,
    )
    text = text.replace(
        "  status_dir: status\n",
        "  status_dir: status\n"
        + "".join(f"  tool_search_dirs:\n    - {d}\n" for d in dirs),
    )
    tmp.write_text(text, encoding="utf-8")
    return tmp


@pytest.mark.parametrize(
    "module", ["scripts.run_downstream_stages", "scripts.preflight"]
)
def test_the_scripts_no_longer_build_their_own_candidate_list(module: str):
    """Importing them must not reintroduce a hard-coded location."""
    if str(PIPELINE_ROOT) not in sys.path:
        sys.path.insert(0, str(PIPELINE_ROOT))
    try:
        source = (PIPELINE_ROOT / (module.replace(".", "/") + ".py")).read_text(
            encoding="utf-8"
        )
    except FileNotFoundError:
        pytest.skip(f"{module} is not present")
    assert "/opt/homebrew" not in source, f"{module} still hard-codes a location"
