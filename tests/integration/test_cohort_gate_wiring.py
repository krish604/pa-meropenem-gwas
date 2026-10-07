"""The cohort gate, wired: evaluated before any stage, recorded wherever it ran.

`papipeline/cohort_gate.py` was built, unit-verified and dormant — nothing in
the pipeline called `evaluate_cohort_gate` (`.build/meropenem-config.registry-
notes.md`). This file pins the wiring, and the three mode decisions it had to
make:

* **REAL enforces `cohort_gate.min_resistant`.** The floor exists because
  "an association scan below this size finds lineage, not resistance"; a scan
  over real isolates is the only scan whose numbers are reported as findings.
  Enforcing it is the whole point of the gate, so the default of
  `evaluate_cohort_gate` stays `enforce=True` and is asserted here on the
  committed fixture's own counts — the same cohort a REAL run would see, and
  one REAL would refuse.
* **TEST evaluates and records but does not enforce.** The committed fixtures
  carry 7 R of 20 by construction (`generator_seed=20240617`), against a floor
  of 100. Lowering the floor to fit them would be tuning configuration to make
  a test pass; refusing every TEST run would make the floor a statement about
  synthetic data. So the counts are recorded and the refusal is not raised —
  and the assertion above is what keeps that from becoming "the gate is off".
* **STUB reads no cohort, so it records "not evaluated".** spec.md D8: a stub
  run discovers no manifest and needs no fixture. It cannot count R isolates
  it never looked at, and inventing a count would be the fabrication STUB
  exists to avoid — so it says so instead, and the reason names the mode.

  (This departs from a reading of the wiring brief that says *every* mode
  records counts. STUB cannot; recording "not evaluated" is the honest form of
  recording counts in a mode with no cohort.)
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.cohort_gate import (
    CohortGateError,
    evaluate_cohort_gate,
)
from papipeline.models import RunMode
from papipeline.run import run_pipeline

PIPELINE_ROOT = Path(__file__).resolve().parents[2]

#: The key a reader of a refusal has to be able to edit.
MIN_RESISTANT_KEY = "cohort_gate.min_resistant"

#: The committed fixture's own counts: 20 rows, 7 of them R. See the module
#: docstring for why the gate must not enforce these.
FIXTURE_RESISTANT = 7
FIXTURE_TESTED = 20


def _read_manifest(results_root: Path, mode: str) -> dict:
    path = Path(results_root) / mode / "run_manifest.json"
    assert path.is_file(), f"no run manifest was written at {path}"
    return json.loads(path.read_text(encoding="utf-8"))


class TestTestModeEvaluatesAndRecords:
    """TEST: the join runs, the counts land in the result and the manifest."""

    def test_the_run_records_five_counts_from_the_committed_fixture(
        self, config, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        result = run_pipeline(config=config, mode="TEST")

        gate = result.cohort_gate
        assert gate is not None, (
            "the cohort gate did not run in TEST, so nothing was counted"
        )
        assert gate["status"] == "evaluated", gate
        assert gate["n_tested"] == FIXTURE_TESTED, gate
        assert gate["n_resistant"] == FIXTURE_RESISTANT, gate
        assert gate["n_with_assembly"] == FIXTURE_TESTED, (
            f"every fixture row names a manifest genome: {gate}"
        )
        assert gate["antibiotic"] == "imipenem", gate
        assert gate["min_resistant"] == 100, gate

    def test_the_run_completes_below_the_floor(self, config, tmp_path, monkeypatch):
        """The floor is a REAL guard. TEST records and continues.

        7 < 100, so if this raises, enforcement has leaked out of REAL.
        """
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        result = run_pipeline(config=config, mode="TEST")
        assert result.cohort_gate["enforced"] is False, result.cohort_gate

    def test_the_manifest_carries_the_same_counts(self, config, tmp_path, monkeypatch):
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        run_pipeline(config=config, mode="TEST")

        payload = _read_manifest(tmp_path / "results", "test")
        gate = payload.get("cohort_gate")
        assert gate, "run_manifest.json has no `cohort_gate` key"
        assert gate["status"] == "evaluated"
        assert gate["n_resistant"] == FIXTURE_RESISTANT
        # The manifest is the reproducibility record: it has to say whether the
        # floor was applied, or a reader cannot tell a gated run from a soft
        # one.
        assert gate["enforced"] is False

    def test_the_cohort_membership_is_recorded_not_inferred(
        self, config, tmp_path, monkeypatch
    ):
        """Who is in the cohort is the gate's answer, not a re-count elsewhere."""
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        result = run_pipeline(config=config, mode="TEST")
        ids = result.cohort_gate["cohort_sample_ids"]
        assert len(ids) == 17, (
            f"20 manifest samples, I excluded by policy (default `exclude`), "
            f"SDD and ND never enter: got {len(ids)}"
        )
        assert len(set(ids)) == len(ids), "a sample appears twice in the cohort"
        assert all(str(s).startswith("TEST_PA_") for s in ids), ids


class TestRealModeIsTheOneThatEnforces:
    """Enforcement is the default, asserted on the fixture's own counts.

    No REAL run happens: `resolve_mode` refuses one without
    `runtime.allow_real_mode`, and rule 6 stands. What is asserted is the gate
    itself, called the way a REAL run would call it — so "TEST does not raise"
    above cannot be true because the raise was deleted.
    """

    def test_the_default_is_to_enforce(self, config, tmp_path):
        phenotype_dir = PIPELINE_ROOT / "test_data" / "phenotype"
        manifest_ids = [
            f"TEST_PA_{n:03d}" for n in range(1, 21)
        ]
        with pytest.raises(CohortGateError) as excinfo:
            evaluate_cohort_gate(
                config,
                manifest_ids=manifest_ids,
                phenotype_dir=phenotype_dir,
                antibiotic="imipenem",
            )
        message = str(excinfo.value)
        assert MIN_RESISTANT_KEY in message, message
        assert str(FIXTURE_RESISTANT) in message, message
        assert excinfo.value.context["n_resistant"] == FIXTURE_RESISTANT
        assert excinfo.value.context["min_resistant"] == 100

    def test_enforce_false_returns_the_report_instead(
        self, config, tmp_path, monkeypatch
    ):
        """The TEST path, at the seam: same call, the refusal not raised."""
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        report = evaluate_cohort_gate(
            config,
            manifest_ids=[f"TEST_PA_{n:03d}" for n in range(1, 21)],
            phenotype_dir=PIPELINE_ROOT / "test_data" / "phenotype",
            antibiotic="imipenem",
            enforce=False,
        )
        assert report.n_resistant == FIXTURE_RESISTANT
        assert report.n_resistant < report.min_resistant, (
            "the fixture is no longer below the floor; the two tests above "
            "would stop proving the mode distinction they exist for"
        )


class TestStubRecordsThatItDidNotEvaluate:
    """STUB cannot count a cohort it is designed not to read."""

    def test_the_result_says_not_evaluated_and_names_the_mode(
        self, config, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        result = run_pipeline(config=config, mode="STUB")

        gate = result.cohort_gate
        assert gate is not None, "STUB should still record the gate's absence"
        assert gate["status"] == "not_evaluated", gate
        assert "STUB" in gate["reason"], gate
        assert "n_resistant" not in gate, (
            "a not-evaluated gate must not carry counts; a fabricated count is "
            "the failure mode this key exists to prevent"
        )

    def test_the_manifest_carries_the_same_verdict(
        self, config, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        run_pipeline(config=config, mode="STUB")
        gate = _read_manifest(tmp_path / "results", "stub")["cohort_gate"]
        assert gate["status"] == "not_evaluated"
        assert "STUB" in gate["reason"]


class TestAPhenotypelessScheduleDoesNotFakeOne:
    """`--only validation` has no phenotype to gate on.

    Evaluating anyway would read a table the run does not otherwise need and
    record counts for a cohort nobody is analysing; skipping silently would
    leave a reader of the manifest to infer the gate ran from its absence.
    So it records why it did not run.
    """

    def test_the_gate_is_not_evaluated_and_says_why(
        self, config, tmp_path, monkeypatch
    ):
        monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "results"))
        result = run_pipeline(config=config, mode="TEST", only=["validation"])

        gate = result.cohort_gate
        assert gate is not None, gate
        assert gate["status"] == "not_evaluated", gate
        assert "phenotype" in gate["reason"], gate
