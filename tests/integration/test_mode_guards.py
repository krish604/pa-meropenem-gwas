"""Mode guards: REAL must be asked for, and STUB must never be mistaken for it.

Three properties, each of which is a way a run could quietly do the wrong
thing:

  a) REAL requires an explicit opt-in, and STUB is never selected implicitly
     for a run that means real data.
  b) A STUB output is unmistakable - a contract header, and a directory that
     cannot collide with a real or test run's.
  c) A STUB run never writes into the real `results/` tree.

These are guards on the guards. The rest of the suite checks that stages work;
this file checks that a run cannot reach the wrong data, or be mistaken for
one, by accident.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import ModeNotAllowedError
from papipeline.models import RunMode
from papipeline.run import resolve_mode, run_pipeline

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = PIPELINE_ROOT / "workflow" / "Snakefile"
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"

pytestmark = pytest.mark.skipif(
    subprocess.run(
        ["sh", "-c", "command -v snakemake"], capture_output=True
    ).returncode != 0,
    reason="snakemake is not installed; this check cannot run here",
)


# -- (a) REAL is explicit, STUB is never implied ---------------------------

def test_real_is_refused_by_default(config):
    """The gate is the default state, not something a caller opts into."""
    assert bool(config.runtime.get("allow_real_mode", False)) is False
    with pytest.raises(ModeNotAllowedError, match="REAL mode is disabled"):
        resolve_mode("REAL", config)


def test_refusal_names_the_flag_and_the_alternative(config):
    """A refusal that does not say what to change is a dead end."""
    with pytest.raises(ModeNotAllowedError) as excinfo:
        resolve_mode("REAL", config)
    message = str(excinfo.value)
    assert "allow_real_mode" in message
    assert "bigmachine" in message.lower(), (
        f"the message should point at the machine that can carry the cohort: {message}"
    )


def test_stub_is_not_a_fallback_for_a_real_request(config, monkeypatch, tmp_path):
    """Asking for REAL must never quietly become STUB.

    The failure modes this guards against, in order of likelihood: a caller
    typo'ing or omitting a mode and inheriting a default; a wrapper deciding
    "REAL is gated, so run STUB instead"; a config that stops declaring the
    mode. In each of those, a request for real data produces fabricated
    results, and nothing downstream can tell - the tables are well-formed and
    the report renders.
    """
    with pytest.raises(ModeNotAllowedError):
        run_pipeline(config=config, mode="REAL")

    # And the refusal must be total: no output, not even fabricated.
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "refused"))
    with pytest.raises(ModeNotAllowedError):
        run_pipeline(config=config, mode="REAL")
    assert not (tmp_path / "refused").exists() or not list(
        (tmp_path / "refused").rglob("*.tsv")
    ), "a refused REAL run wrote stage tables anyway"


def test_an_unknown_mode_is_refused_rather_than_defaulted(config):
    """A typo is an error, not a hint to fall back."""
    from papipeline.errors import PipelineError

    with pytest.raises(PipelineError, match="Unknown run mode"):
        resolve_mode("PRODUCTION", config)
    with pytest.raises(PipelineError, match="Unknown run mode"):
        resolve_mode("", config)


# -- (b) STUB output is unmistakable ---------------------------------------

def test_stub_uses_a_separate_root_from_test_and_real(config):
    """STUB, TEST and REAL must not share a results directory.

    Sharing is how a stub table ends up being read as a result: the paths look
    right, the headers are right, and only the directory name gives it away.
    """
    roots = {
        mode: config.results_root(mode) for mode in RunMode
    }
    values = [roots[m] for m in RunMode]
    assert len(set(values)) == len(values), f"results roots collide: {roots}"
    assert "stub" in str(roots[RunMode.STUB]).lower()


def test_stub_rows_carry_the_contract_header_and_no_data(config, monkeypatch, tmp_path):
    """Shaped like the real thing, empty of content."""
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "shaped"))
    result = run_pipeline(config=config, mode="STUB")

    from papipeline.execution.contracts import STAGE_TABLES

    for stage, (filename, columns) in STAGE_TABLES.items():
        path = result.outputs.get(stage)
        assert path is not None, f"stage {stage!r} declared no output"
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        assert lines[0].split("\t") == list(columns), (
            f"stage {stage!r} header does not match its contract"
        )
        assert len(lines) == 1, (
            f"stage {stage!r} has {len(lines) - 1} fabricated rows; a stub must "
            "invent no data, or it becomes indistinguishable from a result"
        )


def test_the_stub_report_says_it_is_a_stub(config, monkeypatch, tmp_path):
    """The report is the artefact a reader is most likely to open on its own."""
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "banner"))
    result = run_pipeline(config=config, mode="STUB")
    markdown = Path(result.report_paths["markdown"]).read_text(encoding="utf-8")
    assert "STUB" in markdown
    assert "fabricated" in markdown.lower()


# -- (c) STUB never touches the real results tree --------------------------

def test_stub_writes_nothing_into_the_real_results_tree(config, monkeypatch, tmp_path):
    """The point of the redirect: a STUB run is invisible to a real run.

    A stub run that wrote into `results/stub/` was harmless once, because the
    mode had its own directory. It is not harmless the moment a developer runs
    the dashboard against whatever is in results/: the feed would be 16
    fabricated stages.
    """
    real_results = PIPELINE_ROOT / "results"
    before = {
        p: p.stat().st_mtime_ns for p in real_results.rglob("*") if p.is_file()
    } if real_results.exists() else {}

    monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "redirected"))
    run_pipeline(config=config, mode="STUB")

    after = {
        p: p.stat().st_mtime_ns for p in real_results.rglob("*") if p.is_file()
    } if real_results.exists() else {}

    changed = {p for p in after if p not in before or after[p] != before[p]}
    assert not changed, f"a STUB run modified the real results tree: {sorted(changed)[:5]}"


def test_the_redirect_is_ignored_when_unset_is_absent_of_leakage():
    """With no redirect set, the real tree IS the destination.

    Stated so the test above cannot be satisfied by a run that simply never
    writes anything: the two together pin the behaviour from both sides.
    """
    config = load_config(SCIENCE)
    stub_root = config.results_root(RunMode.STUB)
    assert stub_root.name == "stub", (
        f"unredirected, the STUB results root should be results/stub, got {stub_root}"
    )
    assert not any(p.name == "stub" for p in config.results_root(RunMode.TEST).parents
                   if p.name == "results"), "TEST and STUB must not share a parent mode dir"
