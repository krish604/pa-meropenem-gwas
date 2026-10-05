"""A tool's version is read by RUNNING it only when that tool is about to be run.

**The defect this closes.** ``adapters.external.detect_all()`` used to iterate
every tool in ``VERSION_FLAGS`` - including ``bakta`` and ``panaroo`` - and
``subprocess.run`` each one's version flag. It is reached from
``run_pipeline`` before stage 1, unconditionally, over a list of every tool the
project has ever heard of. A REAL 10-isolate run's tripwire shims logged exactly
two invocations for it::

    bakta   pid=93671  argv=--version
    panaroo pid=93679  argv=--version

Neither annotates anything and neither touches a genome, so the scientific
constraint held - but the standing rule is that Bakta is not executed at all,
and it was broken twice per run by a version probe. It also broke the rule on a
run configured with ``annotation.reuse_tool_output: require`` PRECISELY so Bakta
would not be executed, because reuse mode decides whether a tool is INVOKED and
a version probe never consults it.

**What is asserted, and how it is observed.** Executable shims named ``bakta``
and ``panaroo`` are placed first on ``PATH``. Each appends its own name and argv
to a log and exits 97. So "was this tool executed" is a file's contents, not an
inference from a log line, and the failure is loud as well as recorded. The
non-vacuity guard invokes a shim directly and requires the log to grow.

**Non-vacuity in both directions**, because "the log is empty" is satisfied by a
detection path that does nothing at all:

* the sweep reports PRESENCE for a shimmed tool and says its version is
  ``UNKNOWN`` - not dropped, not invented;
* a tool that genuinely is about to be invoked IS probed, and its version comes
  back from the shim's own output. ``detect_tools`` defaults to NOT executing, so
  this is asserted against both the default and the explicit ``execute=True``;
* under ``reuse_tool_output: require`` the annotation path returns ``None``
  without probing, and under ``off`` it probes - so the rule keys off reuse mode
  rather than being broken once for both.

**Nothing is executed here.** The only binaries run are the shell shims this
file writes. No real Bakta, Panaroo, or any other bioinformatics tool is
executed by this file, in either direction of the fix.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List, Tuple

import pytest
import yaml

from papipeline.adapters import external
from papipeline.config.loader import REPO_ROOT, RESULTS_ROOT_ENV, load_config

SCIENCE = REPO_ROOT / "config" / "science.yaml"
SMOKE = REPO_ROOT / "config" / "machines" / "smoke.yaml"

#: The two tools this project must never execute, and the two the REAL run's
#: tripwire caught being executed by a version probe.
FORBIDDEN: Tuple[str, ...] = ("bakta", "panaroo")

#: What a shim prints for its version flag. Distinct per tool so a test can tell
#: whose output came back.
SHIM_VERSION = "tripwire {name} 0.0.0"

#: Non-zero, so an accidental invocation is recorded AND fails loudly.
SHIM_EXIT = 97


@pytest.fixture()
def tripwire(tmp_path, monkeypatch):
    """Executable shims for the forbidden tools, first on ``PATH``.

    Returns ``(log_path, shim_dir)``. The log starts EMPTY; every invocation
    appends one ``name<TAB>argv`` line. Nothing here is a bioinformatics tool -
    each shim is a two-line shell script that writes a line and exits 97.
    """
    bindir = tmp_path / "tripwire-bin"
    bindir.mkdir()
    log = tmp_path / "calls.log"
    log.write_text("", encoding="utf-8")

    for name in FORBIDDEN:
        shim = bindir / name
        shim.write_text(
            "#!/bin/sh\n"
            f'printf "%s\\t%s\\n" "{name}" "$*" >> "{log}"\n'
            f'printf "{SHIM_VERSION.format(name=name)}\\n"\n'
            f"exit {SHIM_EXIT}\n",
            encoding="utf-8",
        )
        shim.chmod(0o755)

    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    # The shims must be the ones that resolve, or the test would be measuring the
    # real tools' absence instead of its own shims' silence.
    for name in FORBIDDEN:
        resolved = external.shutil.which(name)
        assert resolved == str(bindir / name), (
            f"`which {name}` resolved to {resolved}, not the tripwire shim; the "
            "test would be measuring something other than its own shims"
        )
    return log, bindir


def _invocations(log: Path) -> List[str]:
    return [line for line in log.read_text(encoding="utf-8").splitlines() if line]


class TestTheShimsThemselvesWork:
    """The guard on the guard: the tripwire records, and this proves it."""

    def test_a_shim_records_its_invocation_and_exits_non_zero(self, tripwire):
        log, bindir = tripwire
        assert _invocations(log) == [], "the log must start empty"
        completed = subprocess.run(
            [str(bindir / "bakta"), "--version"], capture_output=True, text=True
        )
        assert completed.returncode == SHIM_EXIT, (
            f"the shim exited {completed.returncode}; it is supposed to fail loudly "
            "as well as record"
        )
        assert _invocations(log) == ["bakta\t--version"], (
            f"the shim did not record its invocation: {_invocations(log)}"
        )

    def test_the_shim_answers_a_version_question(self, tripwire):
        """So a passing test cannot be explained by the shim refusing to speak."""
        _log, bindir = tripwire
        completed = subprocess.run(
            [str(bindir / "panaroo"), "--version"], capture_output=True, text=True
        )
        assert completed.stdout.strip() == SHIM_VERSION.format(name="panaroo")


class TestTheSweepExecutesNothing:
    """`detect_all()`: presence for every known tool, no subprocess at all."""

    def test_no_forbidden_tool_is_executed(self, tripwire):
        log, _bindir = tripwire
        statuses = external.detect_all()
        assert _invocations(log) == [], (
            "detect_all() executed a tool it must not execute:\n"
            + "\n".join(_invocations(log))
            + "\n\nThis is the defect: detect_all iterates every tool the project "
            "has heard of and runs each one's version flag. Presence comes from "
            "shutil.which, which resolves a name without running anything."
        )
        assert set(statuses) == set(external.VERSION_FLAGS), (
            f"detect_all reported {sorted(statuses)}, not every tool in "
            f"VERSION_FLAGS {sorted(external.VERSION_FLAGS)}"
        )

    def test_a_shimmed_tool_is_still_reported_present_and_unprobed(self, tripwire):
        """Not dropped from the report. Presence yes, version UNKNOWN."""
        _log, bindir = tripwire
        statuses = external.detect_all()
        for name in FORBIDDEN:
            status = statuses[name]
            assert status.available, (
                f"{name} resolves to {bindir / name} but detect_all reported it "
                "absent; the sweep is supposed to report presence without "
                "executing, and reporting it absent is the same information loss"
            )
            assert status.executable == str(bindir / name)
            assert status.version == external.UNPROBED_VERSION, (
                f"{name} reported version {status.version!r}; the sweep must not "
                "measure a version by executing the tool, and must say so rather "
                "than inventing one"
            )

    def test_an_absent_tool_is_still_reported_absent(self, tripwire, monkeypatch):
        """The other half: presence detection must not start claiming everything."""
        monkeypatch.setattr(external.shutil, "which", lambda name: None)
        statuses = external.detect_all()
        assert all(not status.available for status in statuses.values())
        assert all(status.executable is None for status in statuses.values())


class TestAToolAboutToBeUsedIsStillProbed:
    """The other half of the fix: detection is not simply broken."""

    def test_detect_tools_does_not_execute_by_default(self, tripwire):
        log, bindir = tripwire
        statuses = external.detect_tools(FORBIDDEN)
        assert _invocations(log) == [], (
            f"detect_tools executed by default: {_invocations(log)}. A named "
            "detection may probe, but only when the caller says it is about to "
            "run the tool; the default is presence only."
        )
        assert statuses["bakta"].available, (
            "presence detection is broken: a shimmed tool resolved on PATH was "
            "reported absent"
        )
        assert statuses["bakta"].version == external.UNPROBED_VERSION

    def test_execute_true_probes_and_returns_the_tools_own_version(self, tripwire):
        """A stage that is about to invoke the tool gets the real version."""
        log, _bindir = tripwire
        statuses = external.detect_tools(FORBIDDEN, execute=True)
        assert sorted(_invocations(log)) == sorted(
            [f"{name}\t--version" for name in FORBIDDEN]
        ), (
            f"execute=True did not run each named tool's version flag: "
            f"{_invocations(log)}"
        )
        for name in FORBIDDEN:
            assert statuses[name].version == SHIM_VERSION.format(name=name), (
                f"{name} probed but reported version {statuses[name].version!r}, "
                "not the first line of its own output - so the probe is not "
                "actually reading anything"
            )

    def test_probe_on_demand_updates_only_the_named_tool(self, tripwire):
        log, _bindir = tripwire
        statuses = external.detect_all()
        assert _invocations(log) == []
        measured = external.probe_on_demand("bakta", statuses)
        assert measured.version == SHIM_VERSION.format(name="bakta")
        assert _invocations(log) == ["bakta\t--version"], (
            f"probe_on_demand ran more or less than the one tool it was given: "
            f"{_invocations(log)}"
        )
        assert statuses["panaroo"].version == external.UNPROBED_VERSION, (
            "probe_on_demand measured a tool it was not asked about"
        )

    def test_probe_on_demand_refuses_to_invent_a_tool(self, tripwire):
        """An absent tool is not resolved by trying to run something."""
        log, _bindir = tripwire
        statuses = {"nonesuch": external.ToolStatus("nonesuch", None, None, False)}
        result = external.probe_on_demand("nonesuch", statuses)
        assert result.available is False
        assert result.version is None
        assert _invocations(log) == [], (
            f"probing an absent tool executed something: {_invocations(log)}"
        )


class TestReuseModeDecidesWhetherBaktaIsEvenAsked:
    """`require` means Bakta is never invoked, so it is never probed either."""

    @staticmethod
    def _config(tmp_path, mode: str):
        """The committed smoke overlay with `reuse_tool_output` set to `mode`.

        The overlay is read and reloaded rather than hand-written, so the loader
        and the overlay's own exemption rules are exercised. Only the one key is
        changed; `require` is the shipped value and is asserted, not assumed.
        """
        overlay = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
        assert overlay["annotation"]["reuse_tool_output"] == "require", (
            "the committed smoke overlay is supposed to ship `require`; it says "
            f"{overlay['annotation']['reuse_tool_output']!r}"
        )
        overlay["annotation"]["reuse_tool_output"] = mode
        overlay["runtime"]["allow_real_mode"] = True
        path = tmp_path / f"reuse-{mode}.yaml"
        path.write_text(yaml.safe_dump(overlay), encoding="utf-8")
        return load_config(SCIENCE, machine=path)

    def test_under_require_bakta_is_never_probed(self, tripwire, tmp_path):
        import papipeline.run as run

        log, _bindir = tripwire
        config = self._config(tmp_path, "require")
        assert config.reuse_tool_output() == "require"
        tools = external.detect_all()
        version = run._bakta_tool_version(tools, config)
        assert _invocations(log) == [], (
            f"Bakta was probed under reuse_tool_output=require: "
            f"{_invocations(log)}. `require` exists so Bakta is never executed; "
            "asking it for its version executes it."
        )
        assert version is None, (
            f"a Bakta version {version!r} was forwarded to the annotation stage "
            "on a run that never invokes the tool - a version for an execution "
            "that did not happen"
        )

    def test_under_off_bakta_is_probed_because_it_is_about_to_be_run(
        self, tripwire, tmp_path
    ):
        import papipeline.run as run

        log, _bindir = tripwire
        config = self._config(tmp_path, "off")
        assert config.reuse_tool_output() == "off"
        tools = external.detect_all()
        version = run._bakta_tool_version(tools, config)
        assert _invocations(log) == ["bakta\t--version"], (
            "`off` annotates every genome with Bakta, so the version is asked "
            f"for by a real execution; it was not: {_invocations(log)}"
        )
        assert version == SHIM_VERSION.format(name="bakta"), (
            f"the probed version came back as {version!r}"
        )
        # And the map the run manifest records is updated, so `tools_detected`
        # does not go on saying UNKNOWN about a tool this run has measured.
        assert tools["bakta"].version == SHIM_VERSION.format(name="bakta")

    def test_an_absent_bakta_reports_no_version_without_probing(
        self, tripwire, tmp_path, monkeypatch
    ):
        import papipeline.run as run

        log, _bindir = tripwire
        monkeypatch.setattr(external.shutil, "which", lambda name: None)
        config = self._config(tmp_path, "off")
        tools = external.detect_all()
        assert run._bakta_tool_version(tools, config) is None
        assert _invocations(log) == []


class TestTheSweepIsWhatRunPipelineActuallyCalls:
    """The call site, not just the function it calls."""

    def test_run_pipeline_stub_executes_no_forbidden_tool(
        self, tripwire, tmp_path, monkeypatch
    ):
        """``run_pipeline`` reaches detection before stage 1, on every mode.

        STUB rather than REAL: detection happens before any stage runs, so STUB
        reaches the same line without needing a cohort, a genome or a database.
        The tripwire is the whole assertion.
        """
        log, _bindir = tripwire
        monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "stub-results"))
        from papipeline.run import run_pipeline

        run_pipeline(config_path=SCIENCE, config=None, mode="STUB", machine=None)
        offenders = [
            line for line in _invocations(log) if line.split("\t")[0] in FORBIDDEN
        ]
        assert offenders == [], (
            "run_pipeline executed a forbidden tool during detection:\n"
            + "\n".join(offenders)
        )

    def test_the_run_manifest_still_lists_every_tool(self, tripwire, tmp_path, monkeypatch):
        """Nothing was dropped from the record to achieve that."""
        import json

        log, _bindir = tripwire
        monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "manifest-results"))
        from papipeline.run import run_pipeline

        result = run_pipeline(
            config_path=SCIENCE, config=None, mode="STUB", machine=None
        )
        payload = json.loads(Path(result.run_manifest).read_text(encoding="utf-8"))
        detected = payload["tools_detected"]
        assert set(detected) == set(external.VERSION_FLAGS), (
            f"the manifest lists {sorted(detected)}, not every tool the sweep "
            f"knows: {sorted(external.VERSION_FLAGS)}"
        )
        for name in FORBIDDEN:
            assert name in detected
            assert detected[name]["version"] == external.UNPROBED_VERSION, (
                f"{name} appears in the manifest with version "
                f"{detected[name]['version']!r}; a sweep that does not execute a "
                "tool cannot know its version, and must say UNKNOWN"
            )
        assert _invocations(log) == []


class TestUnprobedVersionsDoNotOverwriteThePin:
    """`build_provenance` reads this map, and must not mistake UNKNOWN for one."""

    def test_an_unpinned_reference_keeps_its_own_marker(self, config):
        import papipeline.run as run

        tools = external.detect_all()
        rows = run.build_provenance(config, tools)
        by_tool = {row["tool"]: row for row in rows}
        unpinned = by_tool.get("panaroo")
        assert unpinned is not None, "the reference table pins no panaroo row"
        assert unpinned["version_status"] != "pinned"
        assert unpinned["tool_version"] == unpinned_spec_version(config, "panaroo"), (
            f"panaroo's tool_version became {unpinned['tool_version']!r}; the "
            f"detected value is {external.UNPROBED_VERSION!r}, which is the ABSENCE "
            "of a measurement. Overwriting the pin's own marker with it would "
            "replace 'we have not pinned this' with a string that reads like a "
            "version and hide the gap version_status exists to show."
        )


def unpinned_spec_version(config, tool: str) -> str:
    """The pin's own `tool_version` for `tool`, read from the reference table."""
    for spec in config.references.values():
        if spec.tool == tool:
            return spec.tool_version
    raise AssertionError(f"no reference row for {tool}")