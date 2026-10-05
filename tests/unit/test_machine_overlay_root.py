"""A machine overlay loaded by path must resolve exactly like a named one.

Found while fixing an unrelated thing: `PipelineConfig.data_root_for` returns
`/test_data` when the overlay is passed as a path, where `data_root` returns the
real `test_data/`. `MachineConfig.root` comes from the overlay file's own
ancestry (`path.parent.parent.parent`), which is the repository root for
`config/machines/laptop.yaml` and `/` for an overlay anywhere else.

That split is dormant only because no production caller uses the delegating
accessors. It is not a cosmetic difference: the bounded REAL smoke run is meant
to use a *custom* overlay pointing at `db/smoke_genomes/`, so a path-loaded
overlay is about to become the normal case rather than the exception.

The rule asserted here: an overlay declares machine facts, not a location. Where
those facts are rooted is the pipeline's business, so a path-loaded overlay
resolves against the same root a named one does.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import load_config, load_machine_config
from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"

#: Delegating accessor -> the non-delegating one it must agree with. Any pair
#: where the two disagree is the bug, whichever way round they disagree.
DELEGATING_PAIRS = (
    ("data_root_for", "data_root"),
    ("results_root_for", "results_root"),
    ("intermediate_root_for", "intermediate_root"),
    ("reports_root_for", "reports_root"),
)


def _overlay_copy(tmp_path: Path) -> Path:
    """A byte-identical copy of the laptop overlay, somewhere else entirely."""
    target = tmp_path / "elsewhere" / "laptop.yaml"
    target.parent.mkdir(parents=True)
    shutil.copy(REPO / "config" / "machines" / "laptop.yaml", target)
    return target


class TestAPathLoadedOverlayResolvesLikeANamedOne:
    def test_the_root_is_the_pipeline_root_not_the_overlay_location(self, tmp_path):
        named = load_machine_config("laptop")
        by_path = load_machine_config(_overlay_copy(tmp_path))
        assert by_path.root == named.root, (
            f"a path-loaded overlay rooted at {by_path.root}, a named one at "
            f"{named.root}. The overlay declares machine facts, not a location."
        )

    def test_it_resolves_the_same_through_load_config(self, tmp_path):
        named = load_config(SCIENCE, machine="laptop")
        by_path = load_config(SCIENCE, machine=_overlay_copy(tmp_path))
        assert by_path.machine is not None
        assert by_path.machine.root == named.machine.root

    @pytest.mark.parametrize("mode", list(RunMode))
    def test_every_delegating_accessor_agrees_with_its_twin(self, tmp_path, mode):
        """The general assertion, so adding an accessor cannot reintroduce it.

        Checked in both directions: the pair must be equal, and it must be equal
        to what a *named* overlay produces. A pair that agreed with each other
        but not with the named overlay would be the same bug wearing a
        different hat.
        """
        overlay = _overlay_copy(tmp_path)
        # Per-pair reference, captured from the named overlay. A single variable
        # here would be compared against the wrong pair's value and report a
        # disagreement that is not one.
        reference = {}
        for label, machine in (("named", "laptop"), ("path", overlay)):
            config = load_config(SCIENCE, machine=machine)
            for delegating, direct in DELEGATING_PAIRS:
                got = getattr(config.machine, delegating)(mode)
                want = getattr(config, direct)(mode)
                if label == "named":
                    reference[(delegating, direct)] = want
                    continue
                assert got == want, (
                    f"{delegating} and {direct} disagree for a "
                    f"{label}-loaded overlay in {mode}: {got} vs {want}"
                )
                assert got == reference[(delegating, direct)], (
                    f"{delegating} resolves differently for a {label}-loaded "
                    f"overlay: {got} vs {reference[(delegating, direct)]}"
                )


class TestTheOverlayStillDeclaresItsOwnFacts:
    """Rooting at the pipeline must not flatten what the overlay says."""

    def test_paths_still_come_from_the_overlay(self, tmp_path):
        config = load_config(SCIENCE, machine=_overlay_copy(tmp_path))
        # laptop declares results_root: results
        assert config.results_root(RunMode.TEST).name == "test"
        assert config.results_root(RunMode.TEST).parent.name == "results"

    def test_a_custom_results_root_is_still_honoured(self, tmp_path):
        target = tmp_path / "custom" / "laptop.yaml"
        target.parent.mkdir(parents=True)
        text = (REPO / "config" / "machines" / "laptop.yaml").read_text(encoding="utf-8")
        text = text.replace("  results_root: results", "  results_root: custom_results")
        target.write_text(text, encoding="utf-8")

        config = load_config(SCIENCE, machine=target)
        assert config.results_root(RunMode.TEST).parent.name == "custom_results"
        assert config.results_root(RunMode.TEST).is_absolute()
        # and rooted at the pipeline, not at tmp_path
        assert config.results_root(RunMode.TEST) == (
            REPO / "custom_results" / "test"
        )
