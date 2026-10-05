"""`--dry-run` must plan the run you asked for, not the laptop's.

`scripts/common/run_pipeline.py` parses `--machine`, forwards it to
`run_pipeline(...)` (fixed in an earlier round, see
`test_run_pipeline_machine_flag.py`) and uses it to build the observatory - but
`print_plan` called `load_config(config_path)` bare. `load_config` defaults
`machine="laptop"`, so the branch was:

    if args.dry_run:
        return print_plan(args.config, args.mode, args.only, args.skip)

`machine` parsed, was used two lines later for something else, and never reached
the one call that decides what the plan says. Provenance therefore depended on
which flag you passed: `--dry-run` validated against `laptop.yaml` for every
machine on earth.

Why `--dry-run` is the worst place for that bug, and not a harmless one:

* `--dry-run` is the command a user runs **before** committing to a run, to ask
  "what would this do?". It is the only surface that reports a plan without
  executing one.
* The two overlays differ in exactly the values that decide whether a run is
  viable: `laptop.yaml` caps the cohort at 20 samples (`max_samples: 20`) and
  allows 4 threads / 8192 MB; `bigmachine.yaml` has no cap and allows 32 /
  65536. A bigmachine plan that reported the laptop's 4 threads and 8 GB would
  understate the resources by 8x and imply a 20-sample ceiling the bigmachine
  does not have.
* Nothing in the printed plan distinguished the two. Both overlays declare
  identical `paths`, so `data root`, `results` and `sources` came out the same
  either way. The wrong plan was not merely wrong, it was *indistinguishable*
  from the right one - which is why this went unnoticed and why the fix also
  prints `machine`, `threads` and `memory_mb`.

The tests below therefore use **two machines whose resolved values genuinely
differ** and assert on those values, not on the mere presence of a `machine`
key. A test that only asserted `--machine bigmachine` prints the word
`bigmachine` would be satisfied by an echo of the flag.
"""

from __future__ import annotations

import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "common" / "run_pipeline.py"
SCIENCE = REPO / "config" / "science.yaml"

#: Resolved from config/machines/*.yaml. Read here rather than hard-coded into
#: the assertions below so a change to an overlay shows up as a changed
#: expectation rather than a mysteriously passing test.
EXPECTED = {
    "laptop": {"machine": "laptop", "threads": "4", "memory_mb": "8192"},
    "bigmachine": {"machine": "bigmachine", "threads": "32", "memory_mb": "65536"},
}


def _load_script():
    """Import the script by path; `scripts/` is not a package."""
    spec = importlib.util.spec_from_file_location("_run_pipeline_script", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def script():
    return _load_script()


def _plan(script, machine: str) -> str:
    """The exact stdout of `--dry-run --machine <machine>`."""
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        code = script.main(
            ["--config", str(SCIENCE), "--machine", machine,
             "--mode", "TEST", "--dry-run"]
        )
    assert code == 0, f"--dry-run --machine {machine} exited {code}"
    return buffer.getvalue()


def _field(plan: str, key: str) -> str:
    for line in plan.splitlines():
        if line.startswith(f"{key}"):
            return line.split(":", 1)[1].strip()
    raise AssertionError(
        f"the plan has no {key!r} field, so --machine cannot be verified from "
        f"its own output:\n{plan}"
    )


class TestTheOverlaysActuallyDiffer:
    """Precondition. Without this the tests below would prove nothing."""

    def test_laptop_and_bigmachine_resolve_to_different_values(self):
        from papipeline.config.loader import load_config

        laptop = load_config(SCIENCE, machine="laptop")
        big = load_config(SCIENCE, machine="bigmachine")
        assert laptop.runtime.get("threads") != big.runtime.get("threads")
        assert laptop.runtime.get("memory_mb") != big.runtime.get("memory_mb")
        # And the cap differs structurally: laptop caps, bigmachine omits it.
        assert (laptop.runtime.get("max_samples") or 0) == 20
        assert big.runtime.get("max_samples") is None


class TestDryRunHonoursTheMachine:
    @pytest.mark.parametrize("machine", sorted(EXPECTED))
    def test_the_plan_reports_that_machines_resolved_values(self, script, machine):
        plan = _plan(script, machine)
        for key, value in EXPECTED[machine].items():
            assert _field(plan, key) == value, (
                f"--machine {machine} planned with {key}={_field(plan, key)!r}, "
                f"but that overlay resolves {key}={value!r}"
            )

    def test_the_two_machines_do_not_produce_the_same_plan(self, script):
        """The regression itself, stated as a single fact.

        Fails on the old code because `print_plan` called `load_config` bare:
        both invocations resolved `laptop.yaml` and produced byte-identical
        plans for two machines whose overlays differ by 8x in both threads and
        memory.
        """
        assert _plan(script, "laptop") != _plan(script, "bigmachine")


class TestPrintPlanSignature:
    def test_print_plan_takes_a_machine(self, script):
        import inspect

        assert "machine" in inspect.signature(script.print_plan).parameters, (
            "print_plan has no machine parameter, so the dry-run branch cannot "
            "pass one and --dry-run silently plans against the default overlay"
        )

    def test_the_default_is_still_laptop(self, script):
        """Threading the flag must not change what an unqualified dry-run does."""
        import inspect

        signature = inspect.signature(script.print_plan)
        assert signature.parameters["machine"].default == "laptop"

    def test_the_dry_run_branch_forwards_it(self, script):
        """Asserted against the source: the defect was an omission.

        There is no runtime symptom distinct from the plan values above; this
        only notices the forwarding being deleted again.
        """
        source = SCRIPT.read_text(encoding="utf-8")
        branch = source[source.index("if args.dry_run:"):]
        branch = branch[: branch.index("try:")]
        assert "machine=args.machine" in branch, (
            "the --dry-run branch parses --machine and does not forward it to "
            "print_plan, so the plan is resolved from the default overlay"
        )
