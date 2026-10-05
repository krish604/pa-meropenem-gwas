"""A stage that does not run must say so, in the manifest and in the log.

**The defect.** `run.py`'s `enabled()` returned `False` for four different
reasons and distinguished none of them to anyone reading the run:

  * named in `--skip`;
  * outside `effective` (not requested via `--only`);
  * unbuilt and not requested;
  * **`analysis.<stage>` false in configuration** - and it logged that at
    **INFO**, the same level as `--- stage X ---`.

The fourth is the one that did damage. `config/science.yaml` had no
`analysis.recombination` key, so `stage_enabled("recombination")` was False, the
stage was skipped, and stage 9 then built its tree from whatever alignment it
found - silently unmasked. `run_manifest.json` recorded nothing at all: not the
skip, not the reason, not even the stage's name. A reader of the manifest could
not tell stage 8 had been asked for and not done. The whole run reported
success.

**What this file asserts.**

  1. Every skipped stage is in `run_manifest.json` under `stages_skipped`, with
     the reason that skipped it - in every mode, not only REAL, because the
     mechanism that hid the recombination skip is mode-independent.
  2. The skip is logged at WARNING, not INFO. Asserted through `caplog`,
     because a log level is a promise to whoever reads the log.
  3. Skipping a stage that a REQUESTED stage needs is a refusal naming both
     stages. `--skip gwas` on a full run used to produce a report that claimed a
     complete run without ever running GWAS; `run.py`'s own note on
     `PREREQUISITES["reporting"]` says the report must have run every stage.

**Why the tests fail on the pre-change code.** Assertion 1 raises `KeyError` on
`run_manifest.json` because the key does not exist. Assertion 3 fails because
`run_pipeline(skip=["gwas"])` returned a result rather than raising. Both are
demonstrated in the round-11 report for `dag-resolve`.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import PipelineError
from papipeline.run import PREREQUISITES, STAGE_ORDER, resolve_prerequisites, run_pipeline

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture
def results_root(tmp_path: Path, monkeypatch) -> Path:
    """A throwaway results tree, so nothing reads a previous run's manifest."""
    root = tmp_path / "results"
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(root))
    return root


@pytest.fixture
def config():
    return load_config(SCIENCE, machine="laptop")


def _manifest(results_root: Path) -> dict:
    return json.loads((results_root / "test" / "run_manifest.json").read_text(encoding="utf-8"))


class TestASkippedStageIsRecorded:
    def test_the_manifest_has_a_place_for_it(self, config, results_root):
        """The key exists, so a reader is not left inferring from an absence."""
        run_pipeline(config=config, mode="TEST", skip=["gwas"], only=["gwas"])
        payload = _manifest(results_root)
        assert "stages_skipped" in payload, (
            "run_manifest.json has no `stages_skipped`. A stage that was asked "
            "for and did not run is currently indistinguishable from a stage "
            "that was never asked for, and from a pipeline that has no such "
            "stage. See the module docstring."
        )

    def test_it_names_the_stage_and_the_reason(self, config, results_root):
        run_pipeline(config=config, mode="TEST", skip=["gwas"], only=["gwas"])
        skipped = _manifest(results_root)["stages_skipped"]
        assert "gwas" in skipped, (
            f"gwas was named in --skip and is absent from stages_skipped: {skipped}"
        )
        reason = skipped["gwas"]
        assert isinstance(reason, str) and reason.strip(), (
            "the reason must be a non-empty string. A null or an empty value is "
            "the same silence with more structure."
        )
        assert "skip" in reason.lower(), (
            f"the reason must say it was skipped on request: {reason!r}"
        )

    def test_a_stage_disabled_in_configuration_is_recorded_too(
        self, config, results_root, monkeypatch,
    ):
        """The reason that actually bit: a stage switched off by a missing key.

        This is the recombination case. `enabled()` returned False, logged at
        INFO, and wrote nothing - so the skip was invisible to both the log's
        default reader and the manifest.
        """
        # A config whose `analysis.recombination` is false, built by patching the
        # raw mapping rather than by editing science.yaml: a test must not edit a
        # tracked config to make itself pass.
        from dataclasses import replace

        patched = replace(
            config,
            raw={
                **config.raw,
                "analysis": {**config.raw["analysis"], "recombination": False},
            },
            analysis={**config.analysis, "recombination": False},
        )
        run_pipeline(config=patched, mode="TEST", only=["recombination"])
        skipped = _manifest(results_root)["stages_skipped"]
        assert "recombination" in skipped, (
            "a stage disabled by configuration is a skip, and it must be "
            f"recorded as one. stages_skipped={skipped}"
        )
        assert "recombination" in skipped["recombination"].lower()
        assert "analysis.recombination" in skipped["recombination"], (
            "the reason must name the key that was false. 'disabled' alone does "
            "not tell a reader which of the two gates closed - this one is the "
            "config flag, not runtime.allow_real_mode"
        )

    def test_an_actually_run_stage_is_not_in_the_skipped_map(
        self, config, results_root,
    ):
        run_pipeline(config=config, mode="TEST")
        skipped = _manifest(results_root)["stages_skipped"]
        executed = _manifest(results_root)["stages"]
        assert not (set(skipped) & set(executed)), (
            f"a stage appears as both run and skipped: "
            f"{sorted(set(skipped) & set(executed))}. A reader cannot tell which "
            "happened, which is the ambiguity this key exists to remove."
        )

    def test_a_clean_run_records_nothing_skipped(self, config, results_root):
        """The default TEST run enables every stage, so the map is empty.

        Asserted so the key cannot become a place where every stage is listed
        'skipped' to satisfy the first test.
        """
        run_pipeline(config=config, mode="TEST")
        skipped = _manifest(results_root)["stages_skipped"]
        assert skipped == {}, (
            f"a default TEST run runs all sixteen stages; nothing is skipped. "
            f"stages_skipped={skipped}"
        )


class TestASkipIsLoggedAtWarning:
    def test_it_is_not_info(self, config, results_root, caplog):
        """A log level is a promise about who reads the log.

        `--- stage X ---` is INFO, so a skip at INFO is invisible to anyone
        running at WARNING, which is the level an operator watching a long run
        actually sees.
        """
        with caplog.at_level(logging.DEBUG, logger="papipeline.run"):
            run_pipeline(config=config, mode="TEST", skip=["gwas"], only=["gwas"])
        warnings = [
            r for r in caplog.records
            if r.levelno >= logging.WARNING and "gwas" in r.getMessage()
        ]
        assert warnings, (
            "a stage named in --skip must be logged at WARNING or above. At "
            "INFO it is the same level as the '--- stage X ---' lines it was "
            "meant to contrast with, so it does not read as an omission at all."
        )

    def test_the_reason_is_in_the_log_line(self, config, results_root, caplog):
        with caplog.at_level(logging.DEBUG, logger="papipeline.run"):
            run_pipeline(config=config, mode="TEST", skip=["gwas"], only=["gwas"])
        messages = [
            r.getMessage() for r in caplog.records
            if r.levelno >= logging.WARNING and "gwas" in r.getMessage()
        ]
        assert any("skip" in m.lower() for m in messages), (
            f"the WARNING must carry the reason, not just the stage name: {messages}"
        )


class TestSkippingAStageSomethingRequestedIsARefusal:
    """The half that is not a recording problem but a correctness one."""

    def test_skipping_a_stage_the_report_needs_refuses(self, config, results_root):
        """`reporting` requires `gwas`. A report without GWAS is not a report.

        `run.py` says so itself, on `PREREQUISITES["reporting"]`: "a report that
        claims a complete run should have run every stage". Skipping `gwas` on a
        full run used to produce exactly that report.
        """
        assert "gwas" in PREREQUISITES["reporting"]
        with pytest.raises(PipelineError) as caught:
            run_pipeline(config=config, mode="TEST", skip=["gwas"])
        message = str(caught.value)
        assert "gwas" in message, (
            "the refusal must name the stage that was skipped: " + message
        )
        assert "reporting" in message, (
            "the refusal must name the stage that needed it. A refusal naming "
            "only the skipped stage tells the reader nothing they did not know."
        )

    def test_the_refusal_happens_before_any_stage_runs(self, config, results_root):
        """Failing after annotation has cost the operator hours for nothing.

        The check is on disk: `results/test/intermediate/stages/` must not exist,
        or a stage already wrote output before the refusal fired.
        """
        with pytest.raises(PipelineError):
            run_pipeline(config=config, mode="TEST", skip=["gwas"])
        stages = results_root / "test" / "intermediate" / "stages"
        assert not stages.exists() or not any(stages.iterdir()), (
            "a stage ran before the refusal fired. This is the same ordering "
            "mistake the sample cap has: enforced as soon as the cohort is "
            "known, not after the expensive stages have been paid for."
        )

    @pytest.mark.parametrize("downstream,upstream", sorted(
        (downstream, upstream)
        for downstream, ups in PREREQUISITES.items()
        for upstream in ups
        if downstream in set(STAGE_ORDER) and upstream in set(STAGE_ORDER)
    ))
    def test_every_declared_edge_refuses_when_its_upstream_is_skipped(
        self, config, results_root, downstream: str, upstream: str,
    ):
        """The rule, derived from `PREREQUISITES` rather than restated.

        A hard-coded pair is a tripwire that fires on the next prerequisite to be
        added, which is what `PREREQUISITES` itself was written to avoid.
        """
        with pytest.raises(PipelineError) as caught:
            run_pipeline(config=config, mode="TEST", skip=[upstream], only=[downstream])
        message = str(caught.value)
        assert upstream in message and downstream in message, (
            f"skipping {upstream!r} with {downstream!r} requested must name both. "
            f"Got: {message}"
        )

    def test_a_stage_outside_the_requested_set_is_not_recorded_as_skipped(
        self, config, results_root,
    ):
        """Not an omission, so not a skip.

        `--only convergence` does not schedule `phenotype`. Recording that as a
        skip would make `stages_skipped` a list of every stage the run did not
        do, which is every run's whole stage list minus sixteen - noise, and
        noise is what this key was added to remove. The distinction is
        *scheduled and did not run* versus *never asked for*.
        """
        assert "phenotype" not in PREREQUISITES["convergence"]
        result = run_pipeline(config=config, mode="TEST", only=["convergence"])
        assert result.stage_status["convergence"] == "completed"
        payload = _manifest(results_root)
        assert "phenotype" not in payload["stages"]
        assert "phenotype" not in payload["stages_skipped"], (
            "phenotype was never scheduled, so it is not a skip. If this is "
            "failing because `--only` stopped narrowing the schedule, that is a "
            "different and worse bug."
        )

    def test_the_transitive_case_is_covered_too(self, config, results_root):
        """`reporting` needs `phylogeny`, which needs `recombination`.

        The refusal has to see through `resolve_prerequisites`, or a two-hop gap
        is exactly the one that slips through: `--skip recombination` on a full
        run would otherwise produce a report built on an unmasked tree.
        """
        assert "phylogeny" in resolve_prerequisites({"reporting"})
        assert "recombination" in resolve_prerequisites({"reporting"})
        with pytest.raises(PipelineError) as caught:
            run_pipeline(config=config, mode="TEST", skip=["recombination"])
        message = str(caught.value)
        assert "recombination" in message
        assert "reporting" in message