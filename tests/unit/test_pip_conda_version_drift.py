"""panaroo is pip-installed; the rest of the scientific stack is conda-pinned.

**Why this test exists.** `environment/environment.yml` installs panaroo from a
git source because the bioconda recipe cannot solve on osx-arm64 (see
`docs/environment-arm64.md` section 3). pip resolves independently of conda, and
panaroo's `setup.py` declares:

    install_requires=['networkx', 'gffutils', 'BioPython', 'joblib', 'tqdm',
                      'edlib', 'scipy', 'numpy', 'matplotlib', 'scikit-learn',
                      'plotly', 'dendropy', 'intbitset', 'biocode']

Five of those - `numpy`, `scipy`, `networkx`, `matplotlib`, `scikit-learn` - are
packages this project also pins as conda packages. If pip installs its own copy
or upgrades one, the pinned stack silently stops being what actually runs. That
stack is what keeps laptop and big-machine results comparable, so drift here is
a correctness problem, not a tidiness problem.

`environment/environment.yml` carries a comment saying the same thing. A comment
is not a safeguard: it cannot fail, and nobody reads it at the moment a `pip
install` moves a version. This test can.

**Scope.** Only packages that are BOTH conda-pinned here AND declared by panaroo
are checked. `scikit-learn` is declared by panaroo but is not conda-pinned, so
it is an unpinned transitive dependency - real, and recorded in
`_PANAROO_DECLARES_BUT_WE_DO_NOT_PIN` below rather than silently ignored.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version as installed_version
from pathlib import Path

import pytest
import yaml

ENVIRONMENT_YML = Path(__file__).resolve().parents[2] / "environment" / "environment.yml"

# Declared by panaroo's setup.py AND pinned as a conda package in
# environment/environment.yml. These are the ones where pip can contradict conda.
_GUARDED = {
    "numpy": "numpy",
    "scipy": "scipy",
    "networkx": "networkx",
    # conda package is `matplotlib-base`; the import/distribution name is
    # `matplotlib`. Same thing - conda-forge splits the base install out.
    "matplotlib": "matplotlib-base",
    # scikit-learn is panaroo's, not ours - asserted separately, not here.
}

# Declared by panaroo, installed by pip, and NOT pinned by us. These arrive with
# whatever version panaroo's dependencies resolve to. Listed so the gap is
# visible in the test suite rather than only in a comment.
_PANAROO_DECLARES_BUT_WE_DO_NOT_PIN = {
    "scikit-learn", "gffutils", "biopython", "joblib", "tqdm", "edlib",
    "plotly", "dendropy", "intbitset", "biocode",
}


def _environment() -> dict:
    return yaml.safe_load(ENVIRONMENT_YML.read_text(encoding="utf-8"))


def _conda_pins(env: dict) -> dict[str, str]:
    """Name -> version for every ``name=version`` conda dependency."""
    pins: dict[str, str] = {}
    for dep in env["dependencies"]:
        if not isinstance(dep, str) or "=" not in dep or dep.startswith(("git+", "pip")):
            continue
        name, _, ver = dep.partition("=")
        pins[name.strip().lower()] = ver.strip()
    return pins


def _pip_entries(env: dict) -> list[str]:
    entries: list[str] = []
    for dep in env["dependencies"]:
        if isinstance(dep, dict) and "pip" in dep:
            entries.extend(dep["pip"])
    return entries


class TestTheGuardedSetIsReal:
    """A guard over an empty set passes forever and protects nothing."""

    def test_the_guarded_set_is_not_empty(self):
        assert _GUARDED, "the guarded set has been emptied out"

    def test_every_guarded_distribution_is_actually_conda_pinned(self):
        pins = _conda_pins(_environment())
        missing = {
            dist: pkg for dist, pkg in _GUARDED.items() if pkg not in pins
        }
        assert not missing, (
            f"guarded distributions no longer have a conda pin in "
            f"environment.yml: {missing}. Either re-pin them or drop them from "
            f"the guarded set - do not leave them unchecked."
        )

    def test_panaroo_is_still_installed_via_pip(self):
        """If this fails, panaroo went back to the conda recipe.

        That recipe cannot solve on osx-arm64, so this is the signal that the
        pin has been reverted to something that does not install.
        """
        entries = _pip_entries(_environment())
        panaroo = [e for e in entries if "panaroo" in e]
        assert panaroo, (
            "panaroo is no longer declared in the pip: section of "
            "environment.yml. The bioconda recipe is unsolvable on osx-arm64 "
            "because of its hard prokka dependency - see "
            "docs/environment-arm64.md section 3."
        )

    def test_the_panaroo_pin_names_an_exact_commit(self):
        """A branch or tag pin lets upstream change what runs the science."""
        entry = next(e for e in _pip_entries(_environment()) if "panaroo" in e)
        # The URL itself contains no "@"; the last one separates the ref.
        ref = entry.rpartition("@")[2]
        assert len(ref) == 40 and all(c in "0123456789abcdef" for c in ref), (
            f"panaroo is not pinned to a full 40-character commit: {entry!r}. "
            f"Pin an exact commit so the science cannot move under us."
        )


class TestNoDrift:
    """The actual check: conda pin == what is importable right now."""

    @pytest.mark.parametrize("dist,conda_pkg", sorted(_GUARDED.items()))
    def test_installed_version_matches_the_conda_pin(self, dist, conda_pkg):
        pinned = _conda_pins(_environment())[conda_pkg]
        try:
            actual = installed_version(dist)
        except PackageNotFoundError:  # pragma: no cover - environment problem
            pytest.fail(f"{dist} is not installed at all")

        assert actual == pinned, (
            f"VERSION DRIFT: {dist} is pinned to {pinned} by conda in "
            f"environment.yml, but {actual} is what is installed.\n"
            f"panaroo is pip-installed and declares {dist} as a dependency, so "
            f"a `pip install` has almost certainly shadowed or upgraded the "
            f"conda package. Laptop and big-machine results stop being "
            f"comparable while this differs.\n"
            f"Fix by reinstalling the conda package "
            f"(`micromamba install -n pa-amr {conda_pkg}={pinned}`) and "
            f"reinstalling panaroo with --no-deps, or by constraining the pip "
            f"side. Do not 'fix' this by loosening the conda pin."
        )


class TestKnownUnpinnedTransitiveDependencies:
    """Recorded so the gap is visible, not because it is currently a failure."""

    def test_the_set_is_still_a_subset_of_what_panaroo_declares(self):
        # Guards against this list drifting into fiction. If panaroo's setup.py
        # changes, this is the test that notices.
        assert "scikit-learn" in _PANAROO_DECLARES_BUT_WE_DO_NOT_PIN
        assert not (_PANAROO_DECLARES_BUT_WE_DO_NOT_PIN & set(_GUARDED)), (
            "a package is listed as unpinned and guarded at once"
        )

    def test_scikit_learn_is_really_absent_from_our_pins(self):
        pins = _conda_pins(_environment())
        assert "scikit-learn" not in pins, (
            "scikit-learn is now conda-pinned - move it into the guarded set so "
            "drift in it starts failing the suite"
        )
