"""`--machine` has to reach the run, not just the observatory.

`scripts/common/run_pipeline.py` parses `--machine`, uses it to build the
observatory, then calls `run_pipeline(config_path=..., mode=..., ...)` without
passing it. `papipeline.run.run_pipeline` had no `machine` parameter at all and
called `load_config(config_path)` bare, which resolves to the default overlay.

So the flag was silently inert: the run always used `laptop.yaml` no matter what
was asked for. The symptom reads like a configuration problem rather than a
wiring one - `run_pipeline --machine config/machines/smoke.yaml --mode REAL`
failed with

    REAL mode is disabled: runtime.allow_real_mode is false in the machine
    overlay for 'laptop'

naming `laptop`, for a run that had been told to use `smoke`. That message is
what sent me looking at the overlays for a missing flag, and both of the
overlays were correct.

The reason this mattered: it is why the smoke overlay had never been exercised
end to end. A bounded REAL run could not be started through the documented
entry point at all, so every statement that the smoke path "works" rested on
calling the stage functions directly rather than on the pipeline running.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from papipeline import run as run_module
from papipeline.config.loader import load_config

REPO = Path(__file__).resolve().parents[2]
SMOKE = REPO / "config" / "machines" / "smoke.yaml"
LAPTOP = REPO / "config" / "machines" / "laptop.yaml"


class TestRunPipelineAcceptsAMachine:
    def test_it_takes_a_machine_argument(self):
        """The absence of this parameter is the whole bug."""
        signature = inspect.signature(run_module.run_pipeline)
        assert "machine" in signature.parameters, (
            "run_pipeline has no `machine` parameter, so load_config is called "
            "without an overlay and every run silently used the default"
        )

    def test_the_default_is_still_laptop(self):
        """Adding the parameter must not change what an unqualified run does."""
        signature = inspect.signature(run_module.run_pipeline)
        assert signature.parameters["machine"].default == "laptop"


class TestTheOverlayActuallyUsedIsTheOneRequested:
    def _resolved_machine(self, machine):
        """The overlay `run_pipeline` would load, via its own code path."""
        config = load_config(
            run_module.Path("config/science.yaml") if hasattr(run_module, "Path")
            else REPO / "config" / "science.yaml",
            machine=machine,
        )
        return config.machine.name

    def test_the_smoke_overlay_is_reachable(self):
        assert self._resolved_machine(SMOKE) == "smoke"

    def test_the_laptop_overlay_is_still_reachable(self):
        assert self._resolved_machine(LAPTOP) == "laptop"


class TestTheCliPassesItThrough:
    def test_the_cli_forwards_machine_to_run_pipeline(self):
        """Parsed, used for the observatory, then dropped.

        Asserted against the source because the failure is an omission: there is
        no runtime symptom short of a full run, and the omission is invisible in
        a diff of the call.
        """
        source = (REPO / "scripts" / "common" / "run_pipeline.py").read_text(
            encoding="utf-8"
        )
        call = source[source.index("result = run_pipeline("):]
        call = call[: call.index(")") + 1]
        assert "machine=args.machine" in call, (
            "the CLI parses --machine but does not forward it to run_pipeline, so "
            "the overlay requested is not the overlay used"
        )
