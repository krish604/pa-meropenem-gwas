"""Q6 — launcher tests with a fake runner.

Nothing here executes snakemake, bakta, panaroo, gubbins, iqtree or pyseer.
`launcher.Runner` is the seam; a `FakeRunner` records every call and the
default `Runner.execute` is asserted to refuse.
"""

from __future__ import annotations

import pytest

from conftest import REPO_ROOT, SMALL_LIVE

from dashboard.server import launcher as launcher_module


class FakeRunner(launcher_module.Runner):
    """Records calls; executes nothing."""

    def __init__(self) -> None:
        super().__init__()
        self.dry_calls = []
        self.execute_calls = []

    def dry_run(self, repo_root, extra=()):
        self.dry_calls.append((str(repo_root), tuple(extra)))
        return launcher_module.RunnerResult(
            argv=(*self.executable, "--dry-run", launcher_module.FULL_RULE, *extra),
            returncode=0,
            stdout="FAKE PLAN",
            stderr="",
            dry_run=True,
        )

    def execute(self, argv, env, cwd):
        self.execute_calls.append(tuple(argv))
        raise launcher_module.LauncherRefused("the fake runner never executes anything")


def _typed_phrase_check(payload):
    for check in payload.get("real_gates", []):
        if check["name"] == "typed_phrase":
            return check
    raise AssertionError("no typed_phrase check in the preflight payload")


def test_launcher_is_disabled_by_default(app_factory):
    client, _state = app_factory(results_root=SMALL_LIVE)
    capabilities = client.get("/api/launcher/capabilities").json()
    assert capabilities["can_launch"] is False
    assert capabilities["reason"].strip()

    payload = client.post("/api/launcher/preflight", json={}).json()
    assert payload["would_start"] is False
    assert payload["plan"] is None
    checks = {check["name"]: check for check in payload["checks"]}
    assert checks["launch_enabled"]["passed"] is False
    for tool in ("bakta", "panaroo", "gubbins", "iqtree", "pyseer"):
        assert tool in capabilities["expensive_tools"]


def test_default_action_is_a_dry_run_and_nothing_is_executed(app_factory):
    fake = FakeRunner()
    client, _state = app_factory(results_root=SMALL_LIVE, allow_launch=True, runner=fake)
    payload = client.post("/api/launcher/preflight", json={}).json()

    assert payload["would_start"] is False
    plan = payload["plan"]
    assert plan["kind"] == "dry_run"
    assert "full_run" in plan["argv"]
    assert "--dry-run" in plan["argv"] or "-n" in plan["argv"]
    assert fake.dry_calls, "the dry run must go through the runner seam"
    assert fake.execute_calls == [], "nothing may be executed"
    # The plan is a plan: the runner returns stdout but the dashboard says the
    # command was not run.
    assert "NOT executed" in plan["note"]


def test_default_runner_execute_refuses():
    with pytest.raises(launcher_module.LauncherRefused):
        launcher_module.Runner().execute(("snakemake", "full_run"), {}, REPO_ROOT)


@pytest.mark.parametrize(
    "phrase",
    [
        "run real",
        "run real sample",
        "yes",
        "Run Real Samples",
        "RUN REAL SAMPLES",
        "runrealsamples",
        "run  real samples",
        "run real samples.",
        "",
        " run real samples ",
        "\trun real samples\n",
    ],
)
def test_real_refuses_everything_but_the_exact_phrase(phrase):
    payload = launcher_module.preflight(
        allow_launch=False,
        repo_root=REPO_ROOT,
        results_root=None,
        log_path=None,
        mode="real",
        phrase=phrase,
        overlay="bigmachine",
    )
    check = _typed_phrase_check(payload)
    assert check["passed"] is False, (
        f"phrase {phrase!r} passed the REAL gate. DESIGN requires the phrase "
        f"typed exactly; dashboard/server/launcher.py:302 compares "
        f"`phrase.strip()`, which accepts surrounding whitespace."
    )


def test_real_accepts_only_the_exact_phrase():
    payload = launcher_module.preflight(
        allow_launch=False,
        repo_root=REPO_ROOT,
        results_root=None,
        log_path=None,
        mode="real",
        phrase="run real samples",
        overlay="bigmachine",
    )
    assert _typed_phrase_check(payload)["passed"] is True


def test_unknown_mode_is_refused_and_no_command_is_interpolated(app_factory):
    fake = FakeRunner()
    client, _state = app_factory(results_root=SMALL_LIVE, allow_launch=True, runner=fake)

    response = client.post("/api/launcher/preflight", json={"mode": "rm -rf /"})
    assert response.status_code == 400
    assert "unknown mode" in response.json()["error"]

    # A command in the body cannot become the plan.
    payload = client.post(
        "/api/launcher/preflight",
        json={"mode": "dry_run", "command": "rm -rf /", "phrase": "rm -rf /"},
    ).json()
    assert payload["would_start"] is False
    assert "rm" not in " ".join(payload["plan"]["argv"])
    assert fake.execute_calls == []

    # A REAL attempt with a destructive phrase is refused by the phrase gate.
    real = client.post(
        "/api/launcher/preflight",
        json={"mode": "real", "phrase": "rm -rf /", "overlay": "bigmachine"},
    ).json()
    assert real["would_start"] is False
    assert _typed_phrase_check(real)["passed"] is False
    assert real["real_gates_all_passed"] is False
