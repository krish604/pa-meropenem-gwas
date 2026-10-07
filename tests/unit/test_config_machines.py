"""Tests for the science/machine config split.

Ticket 03. The property that matters here is not that a config loads; it is
that a machine overlay **cannot express a scientific constant at all**, so the
laptop and the analysis machine cannot drift apart on the science. That is the
whole reason for the split, and it is what these tests pin down.

Written before the implementation.
"""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pytest

from papipeline.config.loader import (
    load_config,
    load_machine_config,
)
from papipeline.errors import ConfigError

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"
LAPTOP = PIPELINE_ROOT / "config" / "machines" / "laptop.yaml"
BIGMACHINE = PIPELINE_ROOT / "config" / "machines" / "bigmachine.yaml"

# The rule is not a fixed allowlist of machine keys. It is: an overlay may not
# contain any top-level key that science.yaml also declares. That is derived
# from the science file rather than hard-coded here, so adding a new scientific
# section automatically closes the door behind it. A fixed allowlist would
# instead have to be updated by hand each time, which is exactly the drift this
# split exists to prevent.
MACHINE_OWNED = {"machine", "machine_class", "paths", "runtime", "tools", "reference"}


def _science_keys() -> set[str]:
    import yaml

    return set(yaml.safe_load(SCIENCE.read_text(encoding="utf-8")))


# ---------------------------------------------------------------------------
# the files exist and are shaped as the spec requires
# ---------------------------------------------------------------------------


def test_science_file_exists_and_has_no_paths_or_runtime():
    """Science is machine-independent, so it holds no paths and no resources."""
    assert "paths" not in _science_keys()
    assert "runtime" not in _science_keys()


@pytest.mark.parametrize("overlay", [LAPTOP, BIGMACHINE], ids=["laptop", "bigmachine"])
def test_no_overlay_declares_a_section_that_science_owns(overlay: Path):
    import yaml

    assert overlay.is_file(), f"{overlay} is missing"
    raw = yaml.safe_load(overlay.read_text(encoding="utf-8"))
    collisions = set(raw) & _science_keys()
    assert not collisions, (
        f"{overlay.name} redeclares scientific section(s) {sorted(collisions)}; "
        "the two machines could disagree about them"
    )


def test_every_key_an_overlay_does_use_is_machine_owned():
    import yaml

    for overlay in (LAPTOP, BIGMACHINE):
        raw = yaml.safe_load(overlay.read_text(encoding="utf-8"))
        unknown = set(raw) - MACHINE_OWNED
        assert not unknown, f"{overlay.name} has unrecognised keys: {sorted(unknown)}"


def test_laptop_overlay_declares_a_sample_cap_and_bigmachine_does_not():
    laptop = load_machine_config(LAPTOP)
    assert laptop.max_samples == 20
    assert load_machine_config(BIGMACHINE).max_samples is None


def test_both_machines_declare_the_same_antimicrobial_gate_closed():
    """A safety gate is not a per-machine preference."""
    laptop = load_machine_config(LAPTOP)
    big = load_machine_config(BIGMACHINE)
    assert laptop.allow_real_mode is False
    assert big.allow_real_mode is False


# ---------------------------------------------------------------------------
# the load-bearing property: overlays cannot carry science
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "section",
    ["organism", "gwas", "phenotype", "convergence", "antibiotics", "cooccurrence"],
)
def test_an_overlay_may_not_declare_a_scientific_section(section: str, tmp_path: Path):
    """If an overlay can name a scientific section, the two machines can disagree."""
    overlay = tmp_path / "machines" / "rogue.yaml"
    overlay.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(LAPTOP, overlay)
    with overlay.open("a", encoding="utf-8") as handle:
        handle.write(f"\n{section}:\n  invented_value: 1\n")
    with pytest.raises(ConfigError) as excinfo:
        load_machine_config(overlay)
    assert section in str(excinfo.value)
    assert "science" in str(excinfo.value).lower()


def test_a_scientific_value_changed_in_one_overlay_cannot_affect_the_other(tmp_path: Path):
    """The two overlays are loaded independently; neither sees the other's content."""
    rogue = tmp_path / "machines" / "rogue.yaml"
    rogue.parent.mkdir(parents=True, exist_ok=True)
    rogue.write_text(
        textwrap.dedent(
            """
            machine: rogue
            paths:
              data_root: /somewhere/else
            runtime:
              threads: 999
              memory_mb: 4096
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )
    rogue_cfg = load_machine_config(rogue)
    laptop_cfg = load_machine_config(LAPTOP)
    assert rogue_cfg.threads == 999
    assert laptop_cfg.threads != 999


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------


def test_load_config_composes_science_with_a_named_machine():
    config = load_config(SCIENCE, machine="laptop")
    assert config.threads == 4
    assert config.max_samples == 20


def test_load_config_accepts_a_machine_by_path():
    by_name = load_config(SCIENCE, machine="laptop")
    by_path = load_config(SCIENCE, machine=LAPTOP)
    assert by_name.threads == by_path.threads
    assert by_name.max_samples == by_path.max_samples


def test_the_same_science_yields_the_same_science_on_both_machines():
    """The point of the split: identical science, different resources."""
    laptop = load_config(SCIENCE, machine="laptop")
    big = load_config(SCIENCE, machine="bigmachine")
    for field in ("antibiotics", "allowed_phenotypes", "analysis", "organism",
                  "mechanisms", "regulators", "references", "gwas",
                  "phylogeny", "convergence", "qc"):
        assert getattr(laptop, field) == getattr(big, field), field
    assert laptop.threads != big.threads, "the machines should differ in resources"


def test_both_machines_agree_on_where_the_test_fixtures_live():
    """Fixtures are shared, so the two machines must point at them identically."""
    from papipeline.models import RunMode

    laptop = load_config(SCIENCE, machine="laptop")
    big = load_config(SCIENCE, machine="bigmachine")
    assert laptop.data_root_for(RunMode.TEST) == big.data_root_for(RunMode.TEST)


def test_derived_paths_land_at_the_repository_root_not_inside_config():
    """A derived path must be the real one, not merely the same on both machines.

    Comparing the two machines to each other is not enough: both were once
    wrong in the same way, because the overlay root was computed as
    `parent.parent` while the overlays live in `config/machines/`.
    """
    from papipeline.models import RunMode

    laptop = load_config(SCIENCE, machine="laptop")
    assert laptop.machine.root == PIPELINE_ROOT
    assert laptop.data_root_for(RunMode.TEST) == PIPELINE_ROOT / "test_data"
    assert laptop.data_root_for(RunMode.REAL) == PIPELINE_ROOT / "data"
    assert laptop.results_root(RunMode.TEST) == PIPELINE_ROOT / "results" / "test"


# ---------------------------------------------------------------------------
# tool availability (feeds the laptop gating in ticket 11)
# ---------------------------------------------------------------------------


def test_panaroo_is_available_on_the_laptop():
    """Verified in docs/environment-arm64.md section 3; this makes it executable.

    Available, but not through conda: the bioconda recipe cannot solve on
    osx-arm64 at any version, and panaroo is never actually invoked against
    prokka (its prokka.py is a GFF3 reader with no subprocess call), so it is
    pip-installed from an exact commit instead. Kept separate from the gubbins
    assertion because the two tools are unblocked for unrelated reasons and
    merging them into one "laptop tools" claim would hide that.
    """
    laptop = load_machine_config(LAPTOP)
    assert laptop.tool_available("panaroo") is True


def test_gubbins_is_available_on_the_laptop_by_source_build_only():
    """Available, but not by conda - and the difference is load-bearing.

    conda's gubbins does not run on this machine. The bioconda osx-arm64 build
    of 3.4.3 (py310hdfa5cb7_1, the newest published) links no zlib while
    referencing `_gzopen`, so the lazy-binding stub jumps to a null pointer and
    the binary SIGSEGVs (exit 139) on any existing alignment file - including an
    empty one. `scripts/build_gubbins_from_source.sh` builds the same v3.4.3
    sources into `tools/gubbins/` instead, which does work, and which produces
    byte-identical output once zlib is supplied to the conda binary.

    Distinct from the panaroo assertion above on purpose: both tools are
    available, but for unrelated reasons, and "available" must not be read as
    "the conda package works". If the source build is dropped and only conda
    remains, this must go back to false.

    See https://github.com/nickjcroucher/gubbins/issues/452
    """
    laptop = load_machine_config(LAPTOP)
    assert laptop.tool_available("gubbins") is True


def test_the_laptop_declares_where_its_source_built_gubbins_is():
    """The binary is not on PATH, so the overlay has to say where it is.

    `tool_search_dirs` is the only mechanism the loader consults - PATH first,
    then these directories - so a binary that lives only in `tools/gubbins/` is
    undiscoverable unless it is declared here.

    Declared RELATIVE, and that is the point. An absolute entry would encode the
    absolute location of one checkout: correct in the main worktree, wrong in
    every other clone and in every `git worktree`, and wrong silently -
    `tool_candidates` skips a directory whose tool is absent, so a stale path
    reports the tool as MISSING rather than the config as wrong.
    `MachineConfig.resolved_tool_search_dirs` anchors a relative entry to the
    repository root, which means the same thing in every checkout.
    """
    laptop = load_machine_config(LAPTOP)
    declared = [d for d in laptop.tool_search_dirs if "tools/gubbins" in d]
    assert declared, (
        "the laptop overlay no longer declares tools/gubbins; the source-built "
        f"gubbins will not be found. declared: {list(laptop.tool_search_dirs)}"
    )
    for directory in declared:
        assert not Path(directory).is_absolute(), (
            f"{directory!r} is absolute, so it only resolves in the one checkout "
            "it was written in; declare it relative to the repository root"
        )
    assert str(PIPELINE_ROOT / "tools" / "gubbins" / "bin") in (
        laptop.resolved_tool_search_dirs()
    ), "a relative entry must resolve to <repo root>/tools/gubbins/bin"


def test_the_big_machine_still_uses_the_conda_gubbins():
    """The conda defect is macOS-specific, so bigmachine must not copy the laptop.

    On Linux `libz` is an ordinary system library and the same link that
    produces a segfault here produces a working binary. Repointing bigmachine at
    `tools/gubbins/` would make a Linux machine depend on a checkout-local source
    build for no reason. Asserted separately from the laptop's declaration above
    because the two files are required to disagree.
    """
    big = load_machine_config(BIGMACHINE)
    assert not any("tools/gubbins" in d for d in big.tool_search_dirs), (
        "bigmachine must not declare the laptop's source-built gubbins"
    )


def test_panaroo_availability_on_the_laptop_is_backed_by_the_environment_pin():
    """The flag must agree with how panaroo is actually installed.

    `available: true` is only true because environment.yml installs panaroo from
    a pinned git source rather than from the unsolvable conda recipe. If that
    pin is reverted, this flag is a lie and stage 7 would be gated open onto a
    tool that cannot be installed.
    """
    import yaml

    env = yaml.safe_load(
        (Path(__file__).resolve().parents[2] / "environment" / "environment.yml")
        .read_text(encoding="utf-8")
    )
    pip_entries = [
        entry
        for dep in env["dependencies"]
        if isinstance(dep, dict) and "pip" in dep
        for entry in dep["pip"]
    ]
    assert any("panaroo" in entry for entry in pip_entries), (
        "laptop.yaml says panaroo is available, but environment.yml no longer "
        "installs it via pip. The conda recipe cannot solve on osx-arm64 - see "
        "docs/environment-arm64.md section 3."
    )


def test_biota_and_the_rest_are_marked_available_on_the_laptop():
    laptop = load_machine_config(LAPTOP)
    for tool in ("bakta", "mlst", "amrfinder", "blast", "minimap2", "samtools",
                 "bcftools", "iqtree", "seqkit", "snp-sites", "pyseer", "mafft"):
        assert laptop.tool_available(tool) is True, tool


def test_a_tool_absence_is_reported_rather_than_guessed():
    """An unlisted tool is unknown, not 'available' and not 'missing'."""
    laptop = load_machine_config(LAPTOP)
    assert laptop.tool_available("something_never_heard_of") is None


def test_the_big_machine_declares_panaroo_and_gubbins_available():
    big = load_machine_config(BIGMACHINE)
    assert big.tool_available("panaroo") is True
    assert big.tool_available("gubbins") is True


# ---------------------------------------------------------------------------
# failures
# ---------------------------------------------------------------------------


def test_an_unknown_machine_name_fails_with_the_known_names_listed():
    with pytest.raises(ConfigError) as excinfo:
        load_config(SCIENCE, machine="the-sun-laptop")
    assert "laptop" in str(excinfo.value)
    assert "bigmachine" in str(excinfo.value)


def test_a_missing_overlay_file_fails_clearly():
    with pytest.raises(ConfigError):
        load_config(SCIENCE, machine=Path("/nonexistent/machines/ghost.yaml"))


def test_an_overlay_missing_its_machine_name_is_rejected():
    bad = Path("/tmp") / "nameless.yaml"
    bad.write_text("runtime:\n  threads: 2\n", encoding="utf-8")
    with pytest.raises(ConfigError) as excinfo:
        load_machine_config(bad)
    assert "machine" in str(excinfo.value)


def test_an_overlay_with_no_resources_is_rejected():
    """An overlay that sets no thread count would silently fall back to a default."""
    bad = Path("/tmp") / "resourceless.yaml"
    bad.write_text("machine: ghost\n", encoding="utf-8")
    with pytest.raises(ConfigError) as excinfo:
        load_machine_config(bad)
    assert "threads" in str(excinfo.value)


# ---------------------------------------------------------------------------
# the repository's own configuration
# ---------------------------------------------------------------------------


def test_the_repository_configuration_loads_cleanly():
    config = load_config(SCIENCE, machine="laptop")
    # Both antibiotics the meropenem GWAS build enables; the old assertion
    # `== ("imipenem",)` encoded the pre-build state. imipenem stays first:
    # config.antibiotics[0] is the documented fallback.
    assert config.antibiotics == ("imipenem", "meropenem")
    assert config.antibiotics[0] == "imipenem"
    assert config.raw["organism"]["taxid"] == 287
