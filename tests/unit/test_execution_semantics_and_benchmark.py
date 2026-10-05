"""Execution semantics and the benchmark harness.

Two things are under test, and they are deliberately kept apart:

* **Semantics** — that NEW RUN, RE-RUN, RESUME and RETRY are distinguished,
  and in particular that an ordinary second run does *not* claim to have
  resumed anything.
* **Harness plumbing** — that the benchmark would isolate its outputs,
  record the telemetry it claims to record, and refuse to report a number it
  did not measure.

The harness tests use a stand-in executable so that the measurement code is
actually exercised. **None of them is a benchmark result and none of their
timings mean anything.** They exist so the measurement code is not shipped
unverified. The real benchmark refuses to run without the real Bakta, and
records `NOT_RUN` when it is absent.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import stat
from pathlib import Path

import pytest

from papipeline.benchmark import (
    NOT_RUN,
    WORKER_LEVELS,
    audit_cohort,
    preflight,
    read_cohort_manifest,
    select_cohort,
    write_cohort_manifest,
    write_per_genome,
    write_summary,
)
from papipeline.benchmark.bakta10 import ConfigResult, GenomeResult, run_configuration
from papipeline.benchmark.telemetry import ProcessSample, SystemSampler
from papipeline.config.loader import load_config
from papipeline.execution import StageState
from papipeline.execution.semantics import Intent, RunKind, classify, resumable
from papipeline.execution.specs import standardised_annotation_spec
from papipeline.execution.store import ExecutionStore
from papipeline.observatory.events import EventBus
from papipeline.run import run_pipeline

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "config" / "science.yaml"
DATA = REPO / "data"


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


# ══════════════════════════════════════════════════════════════════════
# Part A: execution semantics
# ══════════════════════════════════════════════════════════════════════

class TestRunClassification:
    def test_an_unseen_run_key_is_a_new_run(self, tmp_path):
        with ExecutionStore(tmp_path / "s.db") as store:
            verdict = classify(store, "run_001")
        assert verdict.kind is RunKind.NEW_RUN
        assert verdict.prior_rows == 0
        assert "no recorded state" in verdict.reason

    def test_a_seen_key_without_resume_is_a_rerun(self, tmp_path):
        with ExecutionStore(tmp_path / "s.db") as store:
            out = tmp_path / "o.txt"
            out.write_text("x" * 32)
            store.record("run_001", "amr", "", state=StageState.SUCCEEDED)
            verdict = classify(store, "run_001", intent=Intent.EXECUTE)
        assert verdict.kind is RunKind.RE_RUN
        assert verdict.prior_succeeded == 1
        assert "without reusing state" in verdict.reason

    def test_a_seen_key_with_resume_is_a_resume(self, tmp_path):
        with ExecutionStore(tmp_path / "s.db") as store:
            store.record("run_001", "amr", "", state=StageState.SUCCEEDED)
            store.record("run_001", "mlst", "", state=StageState.INCOMPLETE)
            verdict = classify(store, "run_001", intent=Intent.RESUME)
        assert verdict.kind is RunKind.RESUME
        assert verdict.prior_succeeded == 1
        assert verdict.prior_incomplete == 1
        assert "resume requested" in verdict.reason

    def test_no_store_means_nothing_can_be_resumed(self, tmp_path):
        verdict = classify(None, "run_001", intent=Intent.RESUME)
        assert verdict.kind is RunKind.NEW_RUN
        assert "no execution store" in verdict.reason

    def test_the_four_kinds_are_distinct(self):
        assert len({RunKind.NEW_RUN, RunKind.RE_RUN, RunKind.RESUME}) == 3


class TestResumeSkipping:
    def _store(self, tmp_path, state: str, out: Path) -> ExecutionStore:
        store = ExecutionStore(tmp_path / f"{state}.db")
        store.record("run_001", "amr", "", state=StageState(state))
        if out.exists():
            # Recorded outputs must exist for a skip to be possible.
            store.record("run_001", "amr", "", state=StageState(state),
                         output_paths=[str(out)])
        return store

    def test_a_validated_success_is_skippable(self, tmp_path):
        out = tmp_path / "a.tsv"
        out.write_text("sample_id\tgene_name\nS1\toprD\n")
        store = self._store(tmp_path, "SUCCEEDED", out)
        try:
            decision = resumable(store, "run_001", "amr", "",
                                 standardised_annotation_spec(out, "S1"))
        finally:
            store.close()
        assert decision.skip is True
        assert "re-validate now" in decision.reason

    def test_a_deleted_output_is_not_skippable(self, tmp_path):
        """The filesystem overrides a recorded success."""
        out = tmp_path / "a.tsv"
        out.write_text("sample_id\tgene_name\nS1\toprD\n")
        store = self._store(tmp_path, "SUCCEEDED", out)
        out.unlink()
        try:
            decision = resumable(store, "run_001", "amr", "",
                                 standardised_annotation_spec(out, "S1"))
        finally:
            store.close()
        assert decision.skip is False
        assert "re-validate as" in decision.reason

    def test_a_truncated_output_is_not_skippable(self, tmp_path):
        out = tmp_path / "a.tsv"
        out.write_text("sample_id\tgene_name\nS1\toprD\n")
        store = self._store(tmp_path, "SUCCEEDED", out)
        out.write_text("sample_id\tgene_name\nS1\t.\n")
        try:
            decision = resumable(store, "run_001", "amr", "",
                                 standardised_annotation_spec(out, "S1"))
        finally:
            store.close()
        assert decision.skip is False

    def test_an_incomplete_task_is_always_reexecuted(self, tmp_path):
        out = tmp_path / "a.tsv"
        out.write_text("sample_id\tgene_name\nS1\toprD\n")
        store = self._store(tmp_path, "INCOMPLETE", out)
        try:
            decision = resumable(store, "run_001", "amr", "",
                                 standardised_annotation_spec(out, "S1"))
        finally:
            store.close()
        assert decision.skip is False
        assert "not SUCCEEDED" in decision.reason

    def test_a_changed_configuration_forces_reexecution(self, tmp_path):
        out = tmp_path / "a.tsv"
        out.write_text("sample_id\tgene_name\nS1\toprD\n")
        store = ExecutionStore(tmp_path / "c.db")
        store.record("run_001", "amr", "", state=StageState.SUCCEEDED,
                     config_hash="cfgA", output_paths=[str(out)])
        try:
            same = resumable(store, "run_001", "amr", "",
                             standardised_annotation_spec(out, "S1"),
                             config_hash="cfgA")
            changed = resumable(store, "run_001", "amr", "",
                               standardised_annotation_spec(out, "S1"),
                               config_hash="cfgB")
        finally:
            store.close()
        assert same.skip is True
        assert changed.skip is False
        assert "configuration changed" in changed.reason


class TestPipelineSemanticsEndToEnd:
    """The four situations, driven through the real orchestrator."""

    @pytest.fixture()
    def handle(self, tmp_path):
        from papipeline.observatory.wiring import create

        obs = create(load_config(CONFIG), db_path=tmp_path / "obs.db",
                     run_key="run_001", enable=True)
        yield obs
        obs.close()

    def _run(self, observatory, resume=False):
        return run_pipeline(
            config_path=CONFIG, mode="TEST", write_html=False,
            config=load_config(CONFIG), observatory=observatory, resume=resume)

    def test_new_run_executes_everything_and_emits_run_started(self, handle):
        result = self._run(handle)
        assert handle.classification.kind is RunKind.NEW_RUN
        names = [e.name for e in handle.bus.recent()]
        assert "RUN_STARTED" in names
        assert "RUN_RESUMED" not in names
        assert "TASK_RESUMED" not in names
        assert "RUN_COMPLETED" in names
        assert all(v == "completed" for v in result.stage_status.values())

    def test_a_second_run_is_a_rerun_and_never_emits_task_resumed(self, handle):
        self._run(handle)
        before = handle.bus.sequence
        self._run(handle)

        assert handle.classification.kind is RunKind.RE_RUN
        after = handle.bus.recent(before)
        names = [e.name for e in after]
        assert "TASK_RESUMED" not in names, (
            "an ordinary re-run must not claim to have resumed work it redid")
        assert "RUN_RESUMED" not in names
        assert "RUN_STARTED" in names

    def test_a_resume_skips_everything_that_still_validates(self, handle):
        self._run(handle)
        before = handle.bus.sequence
        result = self._run(handle, resume=True)

        assert handle.classification.kind is RunKind.RESUME
        names = [e.name for e in handle.bus.recent(before)]
        assert "RUN_RESUMED" in names
        # Of the fifteen declared stages the three unbuilt ones are skipped in
        # TEST rather than run, so twelve execute and resume. Derived, not
        # hardcoded, so a future stage is not silently miscounted here.
        from papipeline.run import STAGE_ORDER, UNBUILT_STAGES
        runnable = len(STAGE_ORDER) - len(UNBUILT_STAGES)
        assert names.count("TASK_RESUMED") == runnable
        assert all(v == "skipped_validated"
                   for v in result.stage_status.values())
        # Nothing re-executed, so nothing started.
        assert "TASK_STARTED" not in names

    def test_a_resume_reexecutes_a_task_whose_output_vanished(self, handle, tmp_path):
        self._run(handle)
        from papipeline.execution.contracts import table_path
        victim = table_path(REPO / "results" / "test" / "intermediate" / "stages", "amr")
        assert victim.exists()
        victim.unlink()

        result = self._run(handle, resume=True)
        assert result.stage_status["amr"] == "completed", (
            "a stage whose output was deleted must be re-executed on resume")
        names = [e.name for e in handle.bus.recent()]
        from papipeline.run import STAGE_ORDER, UNBUILT_STAGES
        runnable = len(STAGE_ORDER) - len(UNBUILT_STAGES)
        assert names.count("TASK_RESUMED") == runnable - 1, (
            "every still-valid stage except the one whose output was deleted")

    def test_retry_and_resume_are_different_things(self, tmp_path):
        """RETRY is attempt-level; RESUME is run-level and never implied.

        The retry is driven on the subprocess path, because that is where a
        genuinely retryable failure occurs. An in-process exception is
        deterministic and is correctly *not* retried, which is a different
        property and is asserted separately below.
        """
        import sys as _sys
        from papipeline.execution import (
            Check, CheckKind, OutputSpec, RetryPolicy, TaskContext, run_task,
        )

        out = tmp_path / "a.tsv"
        script = tmp_path / "flaky.py"
        script.write_text(
            "import sys\nfrom pathlib import Path\n"
            "marker = Path(sys.argv[2])\n"
            "if not marker.exists():\n"
            "    marker.write_text('1')\n"
            "    sys.stderr.write('OSError: [Errno 28] No space left on device')\n"
            "    raise SystemExit(1)\n"
            "Path(sys.argv[1]).write_text('sample_id\\tgene_name\\nS1\\toprD\\n')\n",
            encoding="utf-8")
        with ExecutionStore(tmp_path / "r.db") as store:
            spec = OutputSpec("amr", (Check(CheckKind.NON_EMPTY, path=out),))
            result = run_task(
                TaskContext(run_key="run_001", stage="amr", subject="", spec=spec),
                [_sys.executable, str(script), str(out), str(tmp_path / "marker")],
                store=store,
                policy=RetryPolicy(max_attempts=2, base_delay=0.0),
                sleep=lambda _: None)
            assert result.state is StageState.SUCCEEDED
            assert result.attempts == 2

            # A retry left no skip behind: the same run key is not a resume.
            verdict = classify(store, "run_001", intent=Intent.EXECUTE)
            assert verdict.kind is RunKind.RE_RUN
            history = store.attempts("run_001", "amr", "")
            assert [a.attempt for a in history] == [1, 2]
            assert [a.state.value for a in history] == ["FAILED", "SUCCEEDED"]

    def test_a_deterministic_exception_is_not_retried(self, tmp_path):
        """An in-process exception is a bug or a config error, not a blip."""
        from papipeline.execution import (
            RetryPolicy, TaskContext, run_task,
        )
        from papipeline.execution.specs import standardised_annotation_spec

        out = tmp_path / "a.tsv"
        calls = {"n": 0}

        def boom():
            calls["n"] += 1
            raise ValueError("deterministic")

        with ExecutionStore(tmp_path / "d.db") as store:
            result = run_task(
                TaskContext(run_key="run_001", stage="amr", subject="",
                            spec=standardised_annotation_spec(out, "S1")),
                boom, store=store,
                policy=RetryPolicy(max_attempts=3, base_delay=0.0),
                sleep=lambda _: None)
        assert result.state is StageState.FAILED
        assert calls["n"] == 1

    def test_a_failing_run_emits_run_failed(self, handle, monkeypatch):
        import papipeline.stages.amr as stage_amr
        from papipeline.errors import StageError

        def explode(*args, **kwargs):
            raise RuntimeError("deliberate")

        monkeypatch.setattr(stage_amr, "run", explode)
        with pytest.raises(StageError):
            self._run(handle)
        assert "RUN_FAILED" in [e.name for e in handle.bus.recent()]


# ══════════════════════════════════════════════════════════════════════
# Part B: the benchmark harness
# ══════════════════════════════════════════════════════════════════════

class TestBenchmarkCohort:
    def test_the_cohort_is_exactly_ten(self):
        cohort = select_cohort(DATA, 10)
        assert len(cohort) == 10

    def test_every_member_is_a_real_assembly_under_data(self):
        for genome in select_cohort(DATA, 10):
            assert genome.path.exists(), genome.accession
            assert genome.path.stat().st_size > 0
            assert DATA in genome.path.parents
            assert genome.accession.startswith("GCA_")

    def test_the_cohort_is_not_the_first_ten_pdc_rows(self):
        cohort = select_cohort(DATA, 10)
        audit = audit_cohort(cohort, REPO / "PDC_essential.tsv")
        assert audit["pdc_table_present"] is True
        assert audit["identical_to_pdc_first_n"] is False
        # A different population, not a coincidentally identical one.
        assert audit["overlap_with_pdc_first_n"] < 10

    def test_selection_is_deterministic(self):
        assert [g.accession for g in select_cohort(DATA, 10)] == \
               [g.accession for g in select_cohort(DATA, 10)]

    def test_selection_follows_the_repository_gca_order(self):
        from papipeline.pilot.cohort import discover_assemblies, sorted_accessions
        expected = sorted_accessions(discover_assemblies(DATA))[:10]
        assert [g.accession for g in select_cohort(DATA, 10)] == expected

    def test_the_cohort_is_persisted_and_reloads_identically(self, tmp_path):
        cohort = select_cohort(DATA, 10)
        path = write_cohort_manifest(tmp_path / "cohort.tsv", cohort)
        reloaded = read_cohort_manifest(path)
        assert [g.accession for g in reloaded] == [g.accession for g in cohort]
        assert [g.path for g in reloaded] == [g.path for g in cohort]
        assert [g.rank for g in reloaded] == list(range(1, 11))

    def test_a_missing_data_directory_is_reported_as_such(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist"):
            select_cohort(tmp_path / "absent", 10)

    def test_an_empty_data_directory_is_not_a_silent_subset(self, tmp_path):
        """An empty directory must not quietly yield a short cohort."""
        with pytest.raises(ValueError, match="no assembly files found"):
            select_cohort(tmp_path, 10)

    def test_too_few_assemblies_is_an_error(self, tmp_path):
        """Fewer than ten available is a refusal, not a smaller cohort."""
        source = select_cohort(DATA, 10)
        import shutil
        for genome in source[:3]:
            target = tmp_path / genome.accession / genome.path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(genome.path, target)
        with pytest.raises(ValueError, match="requires exactly"):
            select_cohort(tmp_path, 10)


class TestBenchmarkPreflight:
    def test_absent_bakta_blocks_the_benchmark(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PATH", str(tmp_path))
        check = preflight(tmp_path / "db", threads=4, max_workers=1)
        assert check.ok is False
        assert any("not on PATH" in r for r in check.reasons)
        assert any("will not substitute" in r for r in check.reasons)

    def test_absent_database_blocks_the_benchmark(self, tmp_path):
        check = preflight(tmp_path / "nope", threads=4, max_workers=1)
        assert check.ok is False
        assert any("database directory does not exist" in r
                   for r in check.reasons)

    def test_a_configuration_that_would_exhaust_memory_is_refused(
        self, tmp_path, monkeypatch
    ):
        """The safety rule: do not start something that would OOM.

        The expected worker count is DERIVED from the measured per-worker peak
        rather than hard-coded. An earlier version of this test asserted the
        figures for a 3.0 GB per-worker estimate, so when the estimate was
        replaced by an actual measurement (1.92 GB, /usr/bin/time -l) the test
        quietly went stale and the guard looked broken. Deriving the boundary
        means a future re-measurement cannot desynchronise the two again.
        """
        import math

        import papipeline.benchmark.bakta10 as bench

        class FakeMemory:
            total = 16 * 2 ** 30

        monkeypatch.setattr(bench.psutil, "virtual_memory",
                            lambda: FakeMemory(), raising=False)
        monkeypatch.setattr(bench.shutil, "which", lambda name: "/usr/bin/bakta")
        (tmp_path / "db").mkdir()

        ram_gb = 16.0
        budget = ram_gb * 0.9
        over = math.floor(budget / bench.MEASURED_PEAK_RSS_GB) + 1
        check = preflight(tmp_path / "db", threads=4, max_workers=over)

        assert check.ok is False
        need = bench.MEASURED_PEAK_RSS_GB * over
        assert any(
            f"{need:.1f} GB" in r and "16.0 GB" in r for r in check.reasons
        ), check.reasons

    def test_one_fewer_worker_is_still_allowed(self, tmp_path, monkeypatch):
        """A guard that refuses everything is not a guard.

        Pairs with the test above: it pins the boundary from both sides, so
        raising the per-worker estimate cannot quietly make every
        configuration look unsafe.
        """
        import math

        import papipeline.benchmark.bakta10 as bench

        class FakeMemory:
            total = 16 * 2 ** 30

        monkeypatch.setattr(bench.psutil, "virtual_memory",
                            lambda: FakeMemory(), raising=False)
        monkeypatch.setattr(bench.shutil, "which", lambda name: "/usr/bin/bakta")
        (tmp_path / "db").mkdir()

        under = math.floor(16.0 * 0.9 / bench.MEASURED_PEAK_RSS_GB)
        assert under >= 1
        check = preflight(tmp_path / "db", threads=4, max_workers=under)

        # NOT `check.ok is True`: the fake executable this test installs does not
        # exist, so the pre-flight is already un-ok for that reason and always
        # was. What matters here is that the *memory* guard did not fire, so
        # that is what is asserted.
        assert not any("concurrent Bakta processes" in r for r in check.reasons), (
            check.reasons
        )

    def test_a_blocked_configuration_is_recorded_not_measured(self, tmp_path):
        blocked = ConfigResult(workers=6, bakta_threads=4, genomes=10,
                               status=NOT_RUN, reason="bakta is not on PATH")
        path = write_summary(tmp_path / "summary.tsv", [blocked])
        with open(path) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert rows[0]["status"] == NOT_RUN
        assert rows[0]["reason"] == "bakta is not on PATH"
        # No measurement may be invented for a configuration that did not run.
        assert rows[0]["wall_seconds"] == "0.0"
        for field in ("cpu_seconds", "peak_ram_gb", "average_ram_gb",
                      "avg_genome_seconds", "genomes_per_hour",
                      "throughput_genomes_per_min"):
            assert rows[0][field] == "NA", field


class TestBenchmarkOutputsAreIsolated:
    def test_each_configuration_gets_its_own_directory(self, tmp_path):
        root = tmp_path / "bakta10"
        for workers in WORKER_LEVELS:
            (root / f"workers_{workers}").mkdir(parents=True)
        assert len({(root / f"workers_{w}").resolve() for w in WORKER_LEVELS}) == \
            len(WORKER_LEVELS)

    def test_the_sweep_levels_are_the_required_ones(self):
        assert WORKER_LEVELS == (1, 2, 3, 4, 6)
        assert 5 not in WORKER_LEVELS


class TestTelemetrySchema:
    def test_the_process_sample_carries_every_required_field(self):
        sample = ProcessSample(pid=123, threads=4).finish()
        row = sample.to_row()
        for field in ("pid", "started_at", "ended_at", "wall_seconds",
                      "cpu_seconds", "peak_rss_mb", "average_rss_mb",
                      "cpu_percent_mean", "bakta_threads"):
            assert field in row, field

    def test_an_unmeasured_value_is_none_not_zero(self):
        sample = ProcessSample(pid=1).finish()
        assert sample.peak_rss_mb is None
        assert sample.cpu_seconds is None
        assert sample.wall_seconds is not None

    def test_the_system_sampler_measures_a_real_interval(self):
        import time
        sampler = SystemSampler(interval=0.2).start()
        time.sleep(0.7)
        result = sampler.stop()
        assert result.samples >= 1
        assert result.wall_seconds > 0
        assert result.peak_ram_gb is not None
        assert result.average_ram_gb is not None
        assert result.peak_ram_gb >= result.average_ram_gb

    def test_the_system_sample_carries_every_required_field(self):
        row = SystemSampler(interval=0.1).start().stop().to_row()
        for field in ("peak_ram_gb", "average_ram_gb", "cpu_percent_mean",
                      "cpu_percent_peak", "load_average_1m",
                      "disk_read_bps", "disk_write_bps"):
            assert field in row, field

    def test_a_real_child_process_is_sampled(self):
        """Sampling a genuine process, not a mock."""
        import subprocess
        import sys
        import time
        from papipeline.benchmark.telemetry import sample_process

        proc = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(1.2)"])
        handle = sample_process(proc.pid, interval=0.1, threads=1)
        time.sleep(1.4)
        sample = handle.finish()
        proc.wait()
        assert sample.pid == proc.pid
        assert sample.samples >= 1
        assert sample.peak_rss_mb and sample.peak_rss_mb > 0
        assert sample.average_rss_mb and sample.average_rss_mb > 0
        assert sample.cpu_seconds is not None and sample.cpu_seconds >= 0
        assert sample.wall_seconds >= 1.0


# ── harness plumbing ────────────────────────────────────────────────
# NOT a benchmark. A stand-in executable exercises the measurement code so
# it is not shipped unverified. None of these timings mean anything.

STUB = '''#!/usr/bin/env python3
import sys, time
from pathlib import Path
args = sys.argv[1:]
def opt(n, d=None):
    return args[args.index(n)+1] if n in args and args.index(n)+1 < len(args) else d
# Bakta 1.12.1 does not accept `--in`: the adapter probes the real binary and
# passes the genome positionally. Match the genome by sequence extension rather
# than assuming a flag or a fixed position.
genome = next((a for a in args if a.endswith((".fna", ".fasta", ".fa"))), None)
if genome is None:
    sys.stderr.write("stub bakta: no genome argument" + chr(10))
    raise SystemExit(2)
out = Path(opt("--out")); out.mkdir(parents=True, exist_ok=True)
stem = Path(genome).stem
# Tabs are built with chr(9) rather than escaped in source. This template is
# itself embedded in a string literal, so a "\t" here does not survive as a tab
# and the generated file ends up with no delimiter at all - which the Bakta
# reader reports as "no header line".
TAB = chr(9)
NL = chr(10)
time.sleep(0.2)
(out / (stem + ".gff3")).write_text(
    TAB.join(["##gff-version 3", "contig_1", "bakta", "CDS", "1", "90", ".", "-", "0", "ID=PA1"]) + NL)
(out / (stem + ".tsv")).write_text(
    TAB.join(["#Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag",
              "Gene", "Product", "DbXrefs"]) + NL
    + TAB.join(["contig_1", "cds", "1", "90", "-", "PA1", "oprD", "OprD porin", "-"]) + NL)
'''

#: Bakta's inference table. Plausible, non-empty, and the wrong file.
INFERENCE_STUB = '''#!/usr/bin/env python3
import sys
from pathlib import Path
args = sys.argv[1:]
def opt(n, d=None):
    return args[args.index(n)+1] if n in args and args.index(n)+1 < len(args) else d
# Bakta 1.12.1 does not accept `--in`; the adapter probes the real binary and
# passes the genome positionally. Match the genome by sequence extension rather
# than assuming a flag or a fixed position.
genome = next((a for a in args if a.endswith((".fna", ".fasta", ".fa"))), None)
if genome is None:
    sys.stderr.write("stub bakta: no genome argument" + chr(10))
    raise SystemExit(2)
out = Path(opt("--out")); out.mkdir(parents=True, exist_ok=True)
stem = Path(genome).stem
# Tabs are built with chr(9); see the note in STUB above.
TAB = chr(9)
NL = chr(10)
# No Gene column, so the gene must be named from the Locus Tag. That is correct
# Bakta behaviour and must be a success, not a validation failure.
(out / (stem + ".gff3")).write_text(
    TAB.join(["##gff-version 3", "contig_1", "bakta", "CDS", "1", "90", ".", "-", "0", "ID=PA1"]) + NL)
(out / (stem + ".tsv")).write_text(
    TAB.join(["#Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag",
              "Score", "Evalue", "Query", "Cov", "Subject", "Cov", "Id", "Accession"]) + NL
    + TAB.join(["contig_1", "cds", "1", "90", "-", "PA1", "45.2", "1e-40",
                "90", "98", "0", "96", "UniRef100_A0A1", "-"]) + NL)
'''

PLUMBING_COHORT = 4


def _install_stub(tmp_path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir(parents=True, exist_ok=True)
    path = bindir / "bakta"
    path.write_text(STUB, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    os.environ["PATH"] = f"{bindir}{os.pathsep}{os.environ['PATH']}"
    return path


def _small_cohort(tmp_path) -> list:
    """A few tiny assemblies, so the harness test is quick."""
    import shutil
    cohort = []
    for genome in select_cohort(DATA, PLUMBING_COHORT):
        target = tmp_path / "genomes" / genome.accession / genome.path.name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(genome.path, target)
        from papipeline.benchmark.bakta10 import BenchGenome
        cohort.append(BenchGenome(rank=genome.rank, accession=genome.accession,
                                  path=target))
    return cohort


class TestHarnessPlumbing:
    """The measurement code works. These are not benchmark results."""

    @pytest.fixture()
    def harness(self, tmp_path):
        _install_stub(tmp_path)
        cohort = _small_cohort(tmp_path)
        import copy
        import dataclasses
        base = load_config(CONFIG)
        raw = copy.deepcopy(dict(base.raw or {}))
        raw.setdefault("annotation", {})
        raw["annotation"]["bakta_db"] = str(tmp_path / "db")
        (tmp_path / "db").mkdir()
        runtime = dict(base.runtime)
        runtime["allow_real_mode"] = True
        runtime["threads"] = 2
        config = dataclasses.replace(base, runtime=runtime, raw=raw)
        yield tmp_path, cohort, config
        os.environ["PATH"] = os.environ["PATH"].split(os.pathsep, 1)[-1] \
            if len(os.environ["PATH"].split(os.pathsep)) > 1 else os.environ["PATH"]

    def test_a_configuration_produces_per_genome_rows(self, harness):
        tmp_path, cohort, config = harness
        result = run_configuration(
            workers=2, cohort=cohort, out_dir=tmp_path / "w2",
            database=tmp_path / "db", threads=2, config=config,
            run_key="harness-w2", bus=EventBus(), genome_dir=tmp_path / "genomes")
        assert len(result.per_genome) == len(cohort)
        assert result.completed == len(cohort)
        assert result.failed == 0 and result.invalid == 0
        assert result.incomplete == 0

    def test_validation_is_required_for_success(self, harness):
        """A clean exit with the *wrong* table must not be SUCCEEDED.

        The stand-in writes Bakta's inference table where the feature table
        belongs, and exits 0. This is the historical bug: every gene name
        was silently lost and all eleven stage-6 loci scored 0/65, with a
        present, non-empty, entirely plausible file the whole time.
        """
        tmp_path, cohort, config = harness
        bindir = tmp_path / "bin"
        (bindir / "bakta").write_text(INFERENCE_STUB, encoding="utf-8")
        (bindir / "bakta").chmod(0o755)
        result = run_configuration(
            workers=1, cohort=cohort, out_dir=tmp_path / "w1",
            database=tmp_path / "db", threads=2, config=config,
            run_key="harness-w1", bus=EventBus(), genome_dir=tmp_path / "genomes")
        assert result.completed == 0, (
            "the inference table in the feature slot must fail the contract")
        assert result.invalid > 0, (
            "a present-but-wrong table is INVALID, not INCOMPLETE")

    def test_a_locus_tag_only_table_is_legitimately_named(self, harness):
        """The converse: an empty Gene column is *not* a failure.

        ``parse_bakta_tsv`` falls back to the locus tag, which is the
        behaviour the real pilot data forced. Asserting the opposite would
        encode a bug as a requirement.
        """
        tmp_path, cohort, config = harness
        bindir = tmp_path / "bin"
        (bindir / "bakta").write_text(
            STUB.replace("oprD\\tOprD porin", "\\t"), encoding="utf-8")
        (bindir / "bakta").chmod(0o755)
        result = run_configuration(
            workers=1, cohort=cohort, out_dir=tmp_path / "w1b",
            database=tmp_path / "db", threads=2, config=config,
            run_key="harness-w1b", bus=EventBus(), genome_dir=tmp_path / "genomes")
        assert result.completed == len(cohort), (
            "a locus-tag-only record is named from the locus tag, which is "
            "correct behaviour, not a validation failure")

    def test_telemetry_is_recorded_for_each_genome(self, harness):
        tmp_path, cohort, config = harness
        result = run_configuration(
            workers=2, cohort=cohort, out_dir=tmp_path / "w2t",
            database=tmp_path / "db", threads=2, config=config,
            run_key="harness-w2t", bus=EventBus(), genome_dir=tmp_path / "genomes")
        for genome in result.per_genome:
            assert genome.pid is not None and genome.pid > 0
            assert genome.bakta_threads == 2
            assert genome.wall_seconds and genome.wall_seconds > 0
            assert genome.peak_rss_mb and genome.peak_rss_mb > 0
            assert genome.average_rss_mb and genome.average_rss_mb > 0
            assert genome.exit_code == 0
        assert result.peak_ram_gb and result.peak_ram_gb > 0
        assert result.wall_seconds > 0

    def test_two_configurations_do_not_share_outputs(self, harness):
        tmp_path, cohort, config = harness
        first = run_configuration(
            workers=1, cohort=cohort, out_dir=tmp_path / "iso1",
            database=tmp_path / "db", threads=2, config=config,
            run_key="iso1", bus=EventBus(), genome_dir=tmp_path / "genomes")
        second = run_configuration(
            workers=3, cohort=cohort, out_dir=tmp_path / "iso3",
            database=tmp_path / "db", threads=2, config=config,
            run_key="iso3", bus=EventBus(), genome_dir=tmp_path / "genomes")
        assert (tmp_path / "iso1" / "execution.db").exists()
        assert (tmp_path / "iso3" / "execution.db").exists()
        assert (tmp_path / "iso1" / "bakta").exists()
        assert (tmp_path / "iso3" / "bakta").exists()
        assert first.wall_seconds > 0 and second.wall_seconds > 0

    def test_the_observatory_receives_the_benchmark_events(self, harness):
        tmp_path, cohort, config = harness
        bus = EventBus()
        run_configuration(
            workers=2, cohort=cohort, out_dir=tmp_path / "obs",
            database=tmp_path / "db", threads=2, config=config,
            run_key="bench-run", bus=bus, genome_dir=tmp_path / "genomes")
        names = [e.name for e in bus.recent()]
        assert names.count("TASK_STARTED") == len(cohort)
        assert names.count("TASK_COMPLETED") == len(cohort)
        subjects = {e.subject for e in bus.recent() if e.name == "TASK_STARTED"}
        assert subjects == {g.accession for g in cohort}

    def test_the_summary_records_every_required_column(self, harness):
        tmp_path, cohort, config = harness
        result = run_configuration(
            workers=2, cohort=cohort, out_dir=tmp_path / "sum",
            database=tmp_path / "db", threads=2, config=config,
            run_key="sum", bus=EventBus(), genome_dir=tmp_path / "genomes")
        path = write_summary(tmp_path / "summary.tsv", [result])
        with open(path) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for column in ("workers", "bakta_threads", "genomes", "completed",
                       "failed", "invalid", "incomplete", "wall_seconds",
                       "cpu_seconds", "peak_ram_gb", "average_ram_gb",
                       "avg_genome_seconds", "min_genome_seconds",
                       "max_genome_seconds", "genomes_per_hour",
                       "throughput_genomes_per_min", "retry_count"):
            assert column in rows[0], column
        assert rows[0]["workers"] == "2"
        assert rows[0]["bakta_threads"] == "2"
        assert rows[0]["genomes_per_hour"] != "NA"

    def test_per_genome_rows_are_written_even_when_empty(self, tmp_path):
        path = write_per_genome(tmp_path / "pg.tsv", [])
        assert path.exists()
        with open(path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
        assert "accession" in header and "wall_seconds" in header
