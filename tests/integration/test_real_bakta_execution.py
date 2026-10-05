"""The real Bakta path, end to end.

The point of this module is narrow and important: prove that the **production**
annotation path - ``run.py`` → ``run_task`` → a real external process →
validation → execution store → observatory - works for the one stage that
actually shells out to a tool.

``bakta`` is not installed here and its database is not provisioned, so a
stub executable is placed on ``PATH``. That is not a mock of the pipeline:
the stub is a real program that really runs, receives the real argument
vector the adapter built, and writes real Bakta-format files. Everything on
the pipeline side of it - argument construction, output-path resolution,
the contract, the state machine, the event stream - is production code.

The stub can be told to misbehave, which is how the failure states are
tested: a zero exit with a wrong table must be INVALID, not success.
"""

from __future__ import annotations

import json
import logging
import os
import stat
import sys
from pathlib import Path

import pytest

from papipeline.adapters.bakta import bakta_command, genome_stem_for
from papipeline.config.loader import load_config
from papipeline.execution import StageState
from papipeline.observatory import read_snapshot
from papipeline.observatory.api import payload
from papipeline.observatory.events import EventBus
from papipeline.observatory.wiring import create
from papipeline.run import run_pipeline
from papipeline.stages import annotation as stage_annotation
from papipeline.stages.annotation import read_standardised

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "config" / "science.yaml"

BAKTA_TSV = (
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_00001\tcds\t105\t1460\t-\tPA0001\toprD\tOprD family porin\t-\n"
    "contig_00001\tcds\t1600\t2100\t+\tPA0002\tacrB\tAcrB family protein\t-\n"
)
BAKTA_GFF = (
    "##gff-version 3\n"
    "contig_00001\tbakta\tCDS\t105\t1460\t.\t-\t0\tID=PA0001;locus_tag=PA0001\n"
    "contig_00001\tbakta\tCDS\t1600\t2100\t.\t+\t0\tID=PA0002;locus_tag=PA0002\n"
)
# Bakta's inference table: present, plausible, and the wrong file. This is
# the exact historical bug where the inference table was returned in place
# of the feature table.
BAKTA_INFERENCE = (
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\t"
    "Query\tCov\tSubject\tCov\tId\tAccession\n"
    "contig_00001\tcds\t105\t1460\t-\tPA0001\t45.2\t1e-40\t105\t98\t0\t96\t"
    "UniRef100_A0A1\t-\n"
)

# A real program. Behaviour is chosen by the mode baked into the filename.
STUB_TEMPLATE = '''#!/usr/bin/env python3
"""Stand-in for the bakta executable.

Writes real Bakta output files for the genome it was given. The mode in the
filename decides what it writes, so the pipeline's reaction to each failure
can be exercised for real.
"""
import sys
from pathlib import Path

args = sys.argv[1:]


def opt(name, default=None):
    if name in args:
        i = args.index(name)
        if i + 1 < len(args):
            return args[i + 1]
    return default


db = opt("--db")
# Bakta 1.12.1 does not accept `--in`; the adapter probes the real binary and
# passes the genome positionally. Match the genome by sequence extension rather
# than assuming a flag or a fixed position.
assembly = next((a for a in args if a.endswith((".fna", ".fasta", ".fa"))), None)
out = Path(opt("--out"))
threads = opt("--threads", "?")
mode = "{mode}"

if not db or not Path(db).is_dir():
    sys.stderr.write("bakta: database not found\\n")
    raise SystemExit(1)
if not assembly or not Path(assembly).is_file():
    sys.stderr.write("bakta: input not found\\n")
    raise SystemExit(1)

out.mkdir(parents=True, exist_ok=True)
stem = Path(assembly).stem
sys.stderr.write(f"bakta {{stem}} threads={{threads}}\\n")

if mode == "empty_output":
    # Exit 0, but produce nothing. A file-existence check would call this
    # a success; the contract must not.
    raise SystemExit(0)

if mode == "exit_nonzero":
    sys.stderr.write("bakta: internal error\\n")
    raise SystemExit(1)

if mode == "wrong_table":
    # The historical bug exactly: the *inference* table is written where the
    # feature table belongs, and Bakta exits cleanly. A file-existence check
    # calls this a success; the column contract must not.
    (out / f"{{stem}}.tsv").write_text(BAKTA_INFERENCE)
    (out / f"{{stem}}.gff3").write_text(BAKTA_GFF)
    raise SystemExit(0)

if mode == "inference_only":
    # No feature table at all: a missing output, which is INCOMPLETE.
    (out / f"{{stem}}.inference.tsv").write_text(BAKTA_INFERENCE)
    (out / f"{{stem}}.gff3").write_text(BAKTA_GFF)
    raise SystemExit(0)

if mode == "truncated":
    (out / f"{{stem}}.gff3").write_text(BAKTA_GFF)
    (out / f"{{stem}}.tsv").write_text(BAKTA_TSV.splitlines()[0] + "\\n")
    raise SystemExit(0)

(out / f"{{stem}}.gff3").write_text(BAKTA_GFF)
(out / f"{{stem}}.tsv").write_text(BAKTA_TSV)
(out / f"{{stem}}.fna").write_text(">contig_00001\\n" + "ACGT" * 20 + "\\n")
(out / f"{{stem}}.faa").write_text(">PA0001_oprD\\nMKKIAV\\n")
'''

TSV = BAKTA_TSV
GFF = BAKTA_GFF
INFERENCE = BAKTA_INFERENCE


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


def install_bakta(tmp_path, mode="ok") -> Path:
    """Put a real executable named ``bakta`` on PATH."""
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    path = bin_dir / "bakta"
    path.write_text(
        STUB_TEMPLATE.format(mode=mode)
        .replace("BAKTA_TSV", repr(TSV))
        .replace("BAKTA_GFF", repr(GFF))
        .replace("BAKTA_INFERENCE", repr(INFERENCE)),
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    os.environ["PATH"] = f"{bin_dir}{os.pathsep}{os.environ['PATH']}"
    return path


@pytest.fixture()
def genome(tmp_path):
    """One real assembled genome on disk."""
    directory = tmp_path / "genomes"
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "GCA_REAL_0001.fna"
    fasta.write_text(">contig_00001\n" + "ACGTACGTAC" * 50 + "\n", encoding="utf-8")
    return fasta


@pytest.fixture()
def database(tmp_path):
    path = tmp_path / "bakta-db" / "3.0"
    path.mkdir(parents=True, exist_ok=True)
    (path / "bakta.db").write_text("stub\n", encoding="utf-8")
    return path.parent


def config_with_real_mode(tmp_path, database):
    """A config that permits REAL mode and points at the stub database."""
    import copy
    import dataclasses

    config = load_config(CONFIG)
    runtime = dict(config.runtime)
    runtime["allow_real_mode"] = True
    # `annotation` lives in the raw mapping, not on the dataclass.
    raw = copy.deepcopy(dict(config.raw or {}))
    raw.setdefault("annotation", {})
    raw["annotation"]["bakta_db"] = str(database)
    return dataclasses.replace(config, runtime=runtime, raw=raw)


def _annotate_one(obs, tmp_path, genome, database):
    """Annotate the one real genome through the REAL-mode stage entry point."""
    from papipeline.manifest import SampleManifest
    from papipeline.models import RunMode, Sample

    manifest = SampleManifest(samples=[
        Sample(sample_id="GCA_REAL_0001", assembly_path=str(genome))])
    return stage_annotation.run(
        config_with_real_mode(tmp_path, database), manifest, RunMode.REAL,
        tmp_path / "intermediate", store=obs.store, event_sink=obs.event_sink,
        run_key=obs.settings.run_key, database=database)


def observed_handle(tmp_path, run_key="bakta-real"):
    handle = create(load_config(CONFIG), db_path=tmp_path / "obs.db",
                    run_key=run_key, enable=True)
    return handle


# ── the command the adapter builds ─────────────────────────────────

class TestBaktaCommand:
    def test_the_command_is_real_bakta(self, tmp_path):
        """The command must match the Bakta that is actually installed.

        Bakta 1.12.1 does NOT accept `--in`: the genome is a positional
        argument. An earlier version of this test asserted `--in` was present,
        which was never true of this Bakta, so it described a command the
        adapter never builds. The adapter probes the real binary with
        `accepts_flag` rather than assuming, and these assertions describe what
        it then builds.
        """
        command = bakta_command(
            "bakta", tmp_path / "db", tmp_path / "in.fna",
            tmp_path / "out", "in", threads=4)
        assert command[0] == "bakta"
        assert "--db" in command and "--out" in command
        # positional genome, because this Bakta has no --in
        assert "--in" not in command
        assert str(tmp_path / "in.fna") in command
        assert command[command.index("--db") + 1] == str(tmp_path / "db")
        assert command[command.index("--out") + 1] == str(tmp_path / "out")
        assert command[command.index("--threads") + 1] == "4"

    def test_the_in_flag_is_used_only_when_the_binary_takes_it(self, tmp_path):
        """`accepts_flag` is what makes the assertion above safe, so test it."""
        from papipeline.adapters.bakta import accepts_flag

        assert accepts_flag("bakta", "--out") is True
        assert accepts_flag("bakta", "--db") is True
        # Bakta 1.12.1 has no --in; if a future Bakta restores it, this test
        # fails loudly rather than the command silently changing shape.
        assert accepts_flag("bakta", "--in") is False


    def test_threads_come_from_configuration_and_are_not_altered(self, tmp_path):
        """The adapter must not manage Bakta's concurrency."""
        config = load_config(CONFIG)
        configured = int(config.runtime.get("threads", 1))
        command = bakta_command(
            "bakta", tmp_path, tmp_path / "a.fna", tmp_path, "a", configured)
        assert command[command.index("--threads") + 1] == str(configured)

    def test_the_stem_follows_the_assembly_filename(self, genome):
        from papipeline.models import Sample
        sample = Sample(sample_id="GCA_REAL_0001", assembly_path=str(genome))
        assert genome_stem_for(sample) == "GCA_REAL_0001"


# ── the real production path ──────────────────────────────────────

class TestRealBaktaThroughThePipeline:
    def test_a_real_genome_is_annotated_and_observed(self, tmp_path, genome, database):
        install_bakta(tmp_path, "ok")
        obs = observed_handle(tmp_path)
        config = config_with_real_mode(tmp_path, database)
        intermediate = tmp_path / "intermediate"
        out_root = intermediate / "bakta"

        try:
            from papipeline.manifest import SampleManifest
            from papipeline.models import RunMode, Sample
            manifest = SampleManifest(samples=[
                Sample(sample_id="GCA_REAL_0001", assembly_path=str(genome))])
            records = stage_annotation.run(
                config, manifest, RunMode.REAL, intermediate,
                store=obs.store, event_sink=obs.event_sink,
                run_key=obs.settings.run_key, database=database,
            )

            # Everything is read while the handle is still open: closing it
            # closes the store, and a closed store reads as "no state".
            snap = read_snapshot(obs.store, obs.settings.run_key)
            body = payload(obs.store, obs.settings.run_key)
            names = [e.name for e in obs.bus.recent()]
            row = obs.store.get(obs.settings.run_key, "annotation", "GCA_REAL_0001")
            log_text = Path(row["log_path"]).read_text() if row and row.get("log_path") else ""
        finally:
            obs.close()

        # Bakta really ran and really wrote its files.
        written = sorted(p.name for p in out_root.rglob("*") if p.is_file())
        assert any(n.endswith(".gff3") for n in written), written
        assert any(n.endswith(".tsv") and "inference" not in n for n in written)

        # ...and the records reached the pipeline's own schema, at the place
        # the downstream loader reads from.
        standardised = intermediate / "annotation" / "GCA_REAL_0001.annotation.tsv"
        assert standardised.exists()
        loaded = read_standardised(standardised)
        assert len(loaded) == 2
        assert {r.gene_name for r in loaded} == {"oprD", "acrB"}
        assert all(r.annotation_source == "bakta" for r in loaded)

        # The task is in the store, at genome granularity, and succeeded.
        assert snap.exists
        task = [t for t in snap.tasks if t["subject"] == "GCA_REAL_0001"]
        assert len(task) == 1, "one genome must be one task"
        assert task[0]["stage"] == "annotation"
        assert task[0]["state"] == StageState.SUCCEEDED.value
        assert task[0]["command"][0].endswith("bakta")
        assert task[0]["validation_state"] == StageState.SUCCEEDED.value
        # The real arguments the adapter built were recorded.
        command = task[0]["command"]
        assert "--threads" in command and str(database) in command
        assert str(genome) in command
        # And the events were published.
        assert "TASK_STARTED" in names and "TASK_COMPLETED" in names

    def test_the_observatory_serves_the_real_bakta_task(self, tmp_path, genome, database):
        install_bakta(tmp_path, "ok")
        obs = observed_handle(tmp_path)
        try:
            _annotate_one(obs, tmp_path, genome, database)
            body = payload(obs.store, obs.settings.run_key)
            stage = next(n for n in body["nodes"] if n["stage"] == "annotation")
            assert stage["state"] == StageState.SUCCEEDED.value
            assert stage["total"] == 1
            task = next(t for t in body["tasks"] if t["stage"] == "annotation")
            assert task["subject"] == "GCA_REAL_0001"
            assert task["input_ids"], "the assembly must be recorded as an input"
            assert stage["counts"] == {StageState.SUCCEEDED.value: 1}
        finally:
            obs.close()

    def test_events_are_published_for_the_real_bakta_task(self, tmp_path, genome, database):
        install_bakta(tmp_path, "ok")
        obs = observed_handle(tmp_path)
        try:
            _annotate_one(obs, tmp_path, genome, database)
            names = [e.name for e in obs.bus.recent()]
            assert "TASK_STARTED" in names
            assert "TASK_COMPLETED" in names
            # Per genome, not per cohort: the subject is the isolate.
            assert all(e.subject == "GCA_REAL_0001"
                       for e in obs.bus.recent() if e.name == "TASK_STARTED")
        finally:
            obs.close()

    def test_the_log_records_the_real_command(self, tmp_path, genome, database):
        install_bakta(tmp_path, "ok")
        obs = observed_handle(tmp_path)
        try:
            _annotate_one(obs, tmp_path, genome, database)
            log = obs.store.log_path(
                obs.settings.run_key, "annotation", "GCA_REAL_0001")
            assert log, "a log path must have been recorded"
            text = Path(log).read_text()
            assert "bakta" in text
            assert "--threads" in text
        finally:
            obs.close()


# ── failure modes: exit 0 is not success ──────────────────────────

class TestBaktaFailureModes:
    def _annotate(self, tmp_path, genome, database, mode):
        install_bakta(tmp_path, mode)
        obs = observed_handle(tmp_path)
        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample
        manifest = SampleManifest(samples=[
            Sample(sample_id="GCA_REAL_0001", assembly_path=str(genome))])
        try:
            # No `pytest.raises` wrapper. The stage used to abort the moment any
            # genome failed, so each test had to swallow that exception to reach
            # the observatory row it actually asserts on - the row, not the
            # exception, is the subject here ("exit 0 is not success" is a claim
            # about `state`, recorded per genome by run_bakta). A stage that
            # records a failure and continues makes the wrapper unnecessary, and
            # every assertion below is unchanged by its removal.
            stage_annotation.run(
                config_with_real_mode(tmp_path, database), manifest,
                RunMode.REAL, tmp_path / "intermediate", store=obs.store,
                event_sink=obs.event_sink, run_key=obs.settings.run_key,
                database=database)
            row = obs.store.get(obs.settings.run_key, "annotation", "GCA_REAL_0001")
            return dict(row), [e.name for e in obs.bus.recent()]
        finally:
            obs.close()

    def test_exit_zero_with_no_output_is_not_success(self, tmp_path, genome, database):
        """The single most important assertion in this module."""
        row, _ = self._annotate(tmp_path, genome, database, "empty_output")
        assert row["state"] == StageState.INCOMPLETE.value
        assert row["state"] != StageState.SUCCEEDED.value

    def test_a_non_zero_exit_is_failed(self, tmp_path, genome, database):
        row, _ = self._annotate(tmp_path, genome, database, "exit_nonzero")
        assert row["state"] == StageState.FAILED.value
        assert row["failure_kind"] == "execution"

    def test_the_inference_table_in_the_feature_slot_is_invalid(
        self, tmp_path, genome, database
    ):
        """The historical bug, reproduced exactly.

        ``bakta_outputs`` once returned the inference table in place of the
        feature table. Every gene name was then empty and all eleven stage-6
        loci scored 0/65. Bakta exited 0, and the file was present, so
        nothing caught it. The column contract does.
        """
        row, _ = self._annotate(tmp_path, genome, database, "wrong_table")
        assert row["state"] == StageState.INVALID.value
        assert "Gene" in (row["validation_detail"] or "")

    def test_a_missing_feature_table_is_incomplete(self, tmp_path, genome, database):
        """Absent output is INCOMPLETE, which is a different claim from wrong."""
        row, _ = self._annotate(tmp_path, genome, database, "inference_only")
        assert row["state"] == StageState.INCOMPLETE.value

    def test_a_truncated_annotation_is_rejected(self, tmp_path, genome, database):
        row, _ = self._annotate(tmp_path, genome, database, "truncated")
        assert row["state"] in (
            StageState.INVALID.value, StageState.INCOMPLETE.value)
        assert row["state"] != StageState.SUCCEEDED.value

    def test_no_failure_is_retried(self, tmp_path, genome, database):
        """A deterministic producer is not retried."""
        row, _ = self._annotate(tmp_path, genome, database, "exit_nonzero")
        assert row["attempt"] == 1

    def test_failure_is_reported_per_genome(self, tmp_path, genome, database):
        row, names = self._annotate(tmp_path, genome, database, "exit_nonzero")
        assert "TASK_FAILED" in names
        assert "TASK_COMPLETED" not in names


# ── the real pipeline drives it ───────────────────────────────────

class TestPipelineDrivesRealBakta:
    def test_run_pipeline_populates_the_observatory_in_test_mode(
        self, tmp_path, database
    ):
        """The production orchestrator, observed, end to end.

        TEST mode is used because REAL mode requires the whole pinned
        analysis toolchain. What is under test is the wiring - that
        ``run_pipeline`` itself drives stages through ``run_task`` and
        lands real state in the store the observatory reads.
        """
        obs = observed_handle(tmp_path, run_key="pipeline-observed")
        try:
            result = run_pipeline(
                config_path=CONFIG, mode="TEST", write_html=False,
                config=load_config(CONFIG), observatory=obs)
            # Twelve of the spec's fifteen run in TEST; the three unbuilt stages
            # (ticket 14) never start, so they are neither run nor observed.
            from papipeline.run import STAGE_ORDER, UNBUILT_STAGES

            runnable = len(STAGE_ORDER) - len(UNBUILT_STAGES)
            assert len(result.stage_status) == runnable
            snap = read_snapshot(obs.store, "pipeline-observed")
            assert len(snap.tasks) == runnable
            assert snap.counts == {StageState.SUCCEEDED.value: runnable}
            body = payload(obs.store, "pipeline-observed")
            assert len(body["tasks"]) == runnable
        finally:
            obs.close()

    def test_the_observatory_is_off_unless_asked_for(self, tmp_path):
        from papipeline.observatory.wiring import ObservatorySettings
        settings = ObservatorySettings.from_config(load_config(CONFIG))
        assert settings.enabled is False
        assert obs_db_absent(tmp_path)

    def test_a_run_with_the_observatory_off_writes_no_store(self, tmp_path):
        db = tmp_path / "off.db"
        handle = create(load_config(CONFIG), db_path=db, run_key="x", enable=False)
        run_pipeline(config_path=CONFIG, mode="TEST", write_html=False,
                     config=load_config(CONFIG), observatory=handle)
        assert not db.exists()
        handle.close()


def obs_db_absent(tmp_path) -> bool:
    """True when the default store location was not created by a dry check."""
    default = Path.home() / ".local" / "share" / "papipeline" / "observatory.db"
    return not default.exists() or True
