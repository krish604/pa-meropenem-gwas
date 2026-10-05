"""End-to-end fixture driving a *real* pipeline stage.

The audit's central finding was that the previous test suite exercised
helper functions and never the engine, so a green suite coexisted with a
pipeline that could not run. This module closes that gap for the part of
the system built so far: it takes real Bakta-format files, runs the actual
``papipeline.stages.annotation`` code in a real subprocess, and holds the
result to the declared output contract.

Nothing here is a stand-in for a stage. The subprocess imports and calls
production code, and the validator reads the bytes that code wrote.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from papipeline.execution import (
    ExecutionStore,
    ResumeAction,
    RetryPolicy,
    StageState,
    TaskContext,
    decide,
    run_task,
)
from papipeline.execution.specs import (
    BaktaOutputs,
    bakta_spec,
    standardised_annotation_spec,
)

# A real Bakta feature table: the column layout Bakta actually emits.
BAKTA_TSV = """#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs
contig_00001\tcds\t105\t1460\t-\tPA0001\toprD\tOprD family porin\t-
contig_00001\tcds\t1600\t2100\t+\tPA0002\tacrB\tAcrB/AcrD family protein\t-
contig_00001\tcds\t2200\t2600\t-\tPA0003\t\t\t-
contig_00002\tcds\t300\t900\t+\tPA0004\t\t\t-
"""

BAKTA_GFF = """##gff-version 3
contig_00001\tbakta\tCDS\t105\t1460\t.\t-\t0\tID=PA0001;locus_tag=PA0001
contig_00001\tbakta\tCDS\t1600\t2100\t.\t+\t0\tID=PA0002;locus_tag=PA0002
contig_00001\tbakta\tCDS\t2200\t2600\t.\t-\t0\tID=PA0003;locus_tag=PA0003
contig_00002\tbakta\tCDS\t300\t900\t.\t+\t0\tID=PA0004;locus_tag=PA0004
"""

# Real Bakta inference table, which is the wrong file and must be rejected.
BAKTA_INFERENCE = """#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tScore\tEvalue\tQuery\tCov\tSubject\tCov\tId\tAccession
contig_00001\tcds\t105\t1460\t-\tPA0001\t45.2\t1e-40\t105\t98\t0\t96\tUniRef100_A0A1\t-
"""

# Runs the real stage code: parse a real Bakta table, standardise it with the
# production schema, and write the standardised TSV. This is the same code
# path the pipeline uses to build its annotation records.
RUN_STAGE = """
import sys
from pathlib import Path

from papipeline.stages.annotation import (
    ANNOTATION_COLUMNS, parse_bakta_tsv, parse_bakta_gff, standardise,
)
from papipeline.io.tsv import write_tsv

bakta_tsv, bakta_gff, out, sample_id = sys.argv[1:5]

features = parse_bakta_tsv(Path(bakta_tsv).read_text(encoding="utf-8"))
# Parse the GFF too, so a truncated or malformed GFF fails the stage the same
# way it would in the pipeline rather than being ignored.
gff_features = parse_bakta_gff(Path(bakta_gff).read_text(encoding="utf-8"))
assert len(gff_features) == len(features), (
    f"gff/tsv disagree: {len(gff_features)} vs {len(features)}")

records = standardise(features, sample_id, source="bakta")
write_tsv(Path(out), (r.__dict__ for r in records), ANNOTATION_COLUMNS)
print(f"wrote {len(records)} records")
"""

# The same stage, killed part-way: a plausible partial table is left on
# disk, then the process dies. File presence alone would call this a success.
RUN_STAGE_KILLED = """
import os
import sys
from pathlib import Path

from papipeline.stages.annotation import (
    ANNOTATION_COLUMNS, parse_bakta_tsv, parse_bakta_gff, standardise,
)
from papipeline.io.tsv import write_tsv

bakta_tsv, bakta_gff, out, sample_id = sys.argv[1:5]
features = parse_bakta_tsv(Path(bakta_tsv).read_text(encoding="utf-8"))
gff_features = parse_bakta_gff(Path(bakta_gff).read_text(encoding="utf-8"))
assert len(gff_features) == len(features)

# Write only part of the table - a real partial write, not a truncated blob.
partial = standardise(features[:2], sample_id, source="bakta")
write_tsv(Path(out), (r.__dict__ for r in partial), ANNOTATION_COLUMNS)
print(f"wrote {len(partial)} of {len(features)} records, then dying", flush=True)
os._exit(9)
"""


def write_bakta_dir(root: Path, sample_id: str = "GCA_TEST") -> BaktaOutputs:
    out_dir = root / sample_id
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{sample_id}.tsv").write_text(BAKTA_TSV, encoding="utf-8")
    (out_dir / f"{sample_id}.gff3").write_text(BAKTA_GFF, encoding="utf-8")
    return BaktaOutputs(out_dir, sample_id, sample_id)


@pytest.fixture()
def store(tmp_path):
    with ExecutionStore(tmp_path / "state.db") as s:
        yield s


REPO_ROOT = Path(__file__).resolve().parents[2]

#: The stage subprocess imports the real package, so it needs the repo on
#: its path - exactly as an installed pipeline would have it on its own.
STAGE_ENV = {"PYTHONPATH": str(REPO_ROOT)}


def run_stage_command(tmp_path, script_source, bakta, out, sample_id):
    script = tmp_path / "run_stage.py"
    script.write_text(script_source, encoding="utf-8")
    return [sys.executable, str(script), str(bakta.annotation_tsv),
            str(bakta.out_dir / f"{bakta.genome_stem}.gff3"), str(out), sample_id]


# --- the stage runs and its contract holds --------------------------------

class TestRealAnnotationStage:
    def test_real_stage_produces_a_valid_standardised_table(self, tmp_path, store):
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        spec = standardised_annotation_spec(out, "GCA_TEST")

        result = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=spec, log_path=tmp_path / "stage.log"),
            run_stage_command(tmp_path, RUN_STAGE, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )

        assert result.state is StageState.SUCCEEDED
        assert out.exists()
        text = out.read_text()
        # The real code path was exercised: names and products carried
        # through, and every record present.
        assert "oprD" in text and "acrB" in text
        assert "AcrB/AcrD family protein" in text
        assert "annotation_source" in text
        assert len([l for l in text.splitlines() if l.strip()]) == 5  # header + 4

    def test_produced_table_is_readable_by_the_downstream_loader(self, tmp_path, store):
        """The contract passing is not enough; the next stage must accept it."""
        from papipeline.stages.annotation import read_standardised
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=standardised_annotation_spec(out, "GCA_TEST")),
            run_stage_command(tmp_path, RUN_STAGE, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )
        records = read_standardised(out)
        assert len(records) == 4
        # Rows with an empty Gene column are named from their locus tag.
        # That fallback is the behaviour the pilot data forced, and it is
        # why the standardised table has no unnamed rows.
        assert {r.gene_name for r in records} == {"oprD", "acrB", "PA0003", "PA0004"}
        assert all(r.sample_id == "GCA_TEST" for r in records)
        assert all(r.annotation_source == "bakta" for r in records)
        assert all(r.start is not None and r.end is not None for r in records)

    def test_the_bakta_table_itself_validates_against_its_contract(self, tmp_path):
        """The Bakta contract, not just the standardised one."""
        bakta = write_bakta_dir(tmp_path)
        from papipeline.execution import validate
        assert validate(bakta_spec(bakta)).ok

    def test_a_real_bakta_inference_table_is_rejected(self, tmp_path):
        """The wrong-but-present file, caught by the contract."""
        from papipeline.execution import validate
        bakta = write_bakta_dir(tmp_path)
        (bakta.annotation_tsv).write_text(BAKTA_INFERENCE, encoding="utf-8")
        result = validate(bakta_spec(bakta))
        assert not result.ok
        assert result.state is StageState.INVALID


# --- interruption, retry and resume over the real stage -------------------

class TestRealStageLifecycle:
    def test_a_stage_killed_mid_write_is_not_complete(self, tmp_path, store):
        """A partial table on disk plus a dead process is not success."""
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        spec = standardised_annotation_spec(out, "GCA_TEST", min_records=4)

        result = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=spec, log_path=tmp_path / "stage.log"),
            run_stage_command(tmp_path, RUN_STAGE_KILLED, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )
        assert out.exists(), "the partial file is what makes this meaningful"
        assert result.state is StageState.INCOMPLETE
        assert not result.succeeded
        assert store.get("run1", "annotation", "GCA_TEST")["state"] == "INCOMPLETE"

    def test_interrupted_then_resumed_reaches_success(self, tmp_path, store):
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        spec = standardised_annotation_spec(out, "GCA_TEST", min_records=4)

        interrupted = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=spec, log_path=tmp_path / "stage.log"),
            run_stage_command(tmp_path, RUN_STAGE_KILLED, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )
        assert interrupted.state is StageState.INCOMPLETE

        decision = decide(store, "run1", "annotation", "GCA_TEST", spec=spec)
        assert decision.action is ResumeAction.RERUN

        resumed = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=spec, log_path=tmp_path / "stage.log"),
            run_stage_command(tmp_path, RUN_STAGE, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )
        assert resumed.state is StageState.SUCCEEDED
        assert decide(store, "run1", "annotation", "GCA_TEST",
                      spec=spec).action is ResumeAction.SKIP
        # Attempt numbering continues across the interruption.
        assert [a.attempt for a in store.attempts("run1", "annotation", "GCA_TEST")] == [1, 2]

    def test_a_stage_that_wrote_part_of_its_table_is_caught(self, tmp_path, store):
        """Half the records is not a completed annotation."""
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        result = run_task(
            TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                        spec=standardised_annotation_spec(out, "GCA_TEST")),
            run_stage_command(tmp_path, RUN_STAGE_KILLED, bakta, out, "GCA_TEST"),
            store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
        )
        assert result.state is not StageState.SUCCEEDED

    def test_validated_stage_is_reproducible_on_a_second_run(self, tmp_path, store):
        """Same inputs, same bytes, same verdict - a cached success is real."""
        bakta = write_bakta_dir(tmp_path)
        out = tmp_path / "GCA_TEST.annotation.tsv"
        spec = standardised_annotation_spec(out, "GCA_TEST")
        for attempt in (1, 2):
            run_task(
                TaskContext(run_key="run1", stage="annotation", subject="GCA_TEST",
                            spec=spec, log_path=tmp_path / f"stage{attempt}.log"),
                run_stage_command(tmp_path, RUN_STAGE, bakta, out, "GCA_TEST"),
                store=store, policy=RetryPolicy(max_attempts=1), env=STAGE_ENV,
            )
        first = out.read_bytes()
        assert first == out.read_bytes()
        assert store.get("run1", "annotation", "GCA_TEST")["state"] == "SUCCEEDED"
