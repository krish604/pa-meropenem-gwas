"""`require` on the smoke overlay: the tool is unreachable, in BOTH directions.

**The claim under test.** ``config/machines/smoke.yaml`` sets
``annotation.reuse_tool_output: "require"``. That mode is the only one of the
three in which Bakta is never invoked, and the property that matters is a
property of the *control flow*: a genome whose curated output is absent must be
recorded as failed, by name, with the tool never reached. ``prefer`` would reuse
verified curated output and then EXECUTE BAKTA for the rest, which is a second
way to get the answer this overlay exists to avoid.

Round 12 set this mode, reverted it, and gave a reason that was measured rather
than assumed. This file is the measurement, and it is deliberately written in
two directions because a one-directional test of a guard is satisfied by a stage
that annotates nothing at all:

* **curated output present for all ten** - the run annotates all ten, reuses
  every one, and the injected runner that raises on call is never called.
* **curated output absent for one** - that isolate is refused WITH ITS NAME, the
  other nine are unaffected, and the runner is still never called. A cohort that
  is mostly unprepared must not be ended by the isolate that is.

**Nothing here is real and nothing here is a mock of the pipeline.** The cohort
is ten synthetic samples with synthetic sequences (R12). The curated output is
synthetic Bakta-format output written to ``<intermediate_root>/bakta/<id>/``,
which is where ``run_bakta`` would have written it - so the *path contract* is
real and exercised, while the bytes are ours. ``run_bakta`` is replaced by a
function that raises, so "no Bakta invocation occurred" is an observed fact and
not an inference from a log line.

**Ten, not 967.** The cohort size here mirrors the smoke overlay's
``max_samples: 10`` deliberately. A previous round read ``discover_pdc_manifest``'s
unfiltered 967-member membership as the smoke cohort and concluded that
``require`` refused "all 967 roster members"; the smoke run's manifest is
narrowed to ten by ``cohort.subset_file`` at ``papipeline/run.py:501``. The
number matters for the cost of ``require`` (ten curated trees must be staged,
not 967), not for whether the tool is reachable.
"""

from __future__ import annotations

import copy
import logging
from pathlib import Path

import pytest

from papipeline.config.loader import REPO_ROOT, load_config
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import annotation as stage

SCIENCE = REPO_ROOT / "config" / "science.yaml"
SMOKE = REPO_ROOT / "config" / "machines" / "smoke.yaml"

#: The overlay's own cap. Ten is the cohort this file builds, and the number the
#: staging cost of `require` is quoted against.
SMOKE_COHORT_SIZE = 10

#: The database string `expected_bakta_database_string` derives from
#: `version.json` below. The curated output's header must match it byte for
#: byte or `decide_reuse` refuses on `database_mismatch`, which would make
#: every test here pass for the wrong reason.
DATABASE_STRING = "v6.0, light"

FEATURE_TSV = (
    f"# Annotated with Bakta\n"
    f"# Software: v1.12.1\n"
    f"# Database: {DATABASE_STRING}\n"
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "contig_00001\tcds\t1\t300\t+\tTESTPA_00001\toprD\tOprD family porin\t-\n"
)

GFF3 = (
    "##gff-version 3\n"
    "contig_00001\tBakta\tCDS\t1\t300\t.\t+\t0\t"
    "ID=TESTPA_00001;Name=TESTPA_00001\n"
)

INFERENCE_TSV = (
    f"# Software: v1.12.1\n"
    f"# Database: {DATABASE_STRING}\n"
    "Sequence Id\tLocus Tag\tScore\tEvalue\n"
    "contig_00001\tTESTPA_00001\t45.2\t1e-40\n"
)


def isolate_ids(n: int = SMOKE_COHORT_SIZE) -> list[str]:
    """Cohort member ids, manifest order preserved."""
    return [f"SMOKE_{i:03d}" for i in range(n)]


def stage_curated_output(out_root: Path, sample_id: str) -> Path:
    """Write synthetic curated Bakta output where `decide_reuse` looks for it.

    ``out_root / sample_id``, with the three files ``decide_reuse`` requires:
    the feature table, the GFF3, and the inference table. That is the real
    contract - see ``annotation.decide_reuse``'s ``required`` tuple - so a
    staged set that this function writes is a set that production would accept.
    """
    from papipeline.adapters.bakta import BaktaOutputs

    directory = Path(out_root) / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    outputs = BaktaOutputs(
        out_dir=directory, sample_id=sample_id, genome_stem=sample_id
    )
    outputs.annotation_tsv.write_text(FEATURE_TSV, encoding="utf-8")
    outputs.gff.write_text(GFF3, encoding="utf-8")
    outputs.inference_tsv.write_text(INFERENCE_TSV, encoding="utf-8")
    return directory


@pytest.fixture()
def smoke_run(tmp_path):
    """Everything one smoke-shaped REAL annotation run needs, under tmp_path.

    The config is the COMMITTED smoke overlay - read, patched to redirect three
    keys, and reloaded - not a hand-built stand-in. So the reuse mode under test
    is the one the repository ships, taken through the same accessor stage 2
    reads. The three patches are the minimum a test needs and each is a
    redirection rather than a behaviour change: ``runtime.allow_real_mode`` (the
    overlay is closed by policy and this is a test), ``paths.smoke_genome_dir``
    (temp sequences), ``annotation.bakta_db`` (a temp database whose
    ``version.json`` yields the database string the curated output's header
    declares).
    """
    import dataclasses

    import yaml

    genomes = tmp_path / "smoke_genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    for sample_id in isolate_ids():
        (genomes / f"{sample_id}.fna").write_text(
            ">contig_00001\nACGTACGTAC\n", encoding="utf-8"
        )

    database = tmp_path / "bakta-db" / "db-light"
    database.mkdir(parents=True, exist_ok=True)
    (database / "bakta.db").write_text("stub\n", encoding="utf-8")
    (database / "version.json").write_text(
        '{"major": 6, "minor": 0, "type": "light"}', encoding="utf-8"
    )

    overlay = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
    # The mode is NOT asserted here. It is asserted in
    # `TestTheCommittedModeIsTheOneUnderTest`, and deliberately not in this
    # fixture: a guard here would abort every behavioural test in the file when
    # the mode regresses, which would hide the thing they exist to catch. Under
    # `prefer` they must FAIL on the runner that raises - that is the proof the
    # property is detected and not merely declared.
    overlay["paths"]["smoke_genome_dir"] = str(genomes)
    overlay["runtime"]["allow_real_mode"] = True
    overlay_path = tmp_path / "smoke-under-test.yaml"
    overlay_path.write_text(yaml.safe_dump(overlay), encoding="utf-8")

    base = load_config(SCIENCE, machine=overlay_path)
    raw = copy.deepcopy(dict(base.raw or {}))
    raw.setdefault("annotation", {})
    raw["annotation"]["bakta_db"] = str(database)
    config = dataclasses.replace(base, raw=raw)

    return {
        "config": config,
        "database": database,
        "genomes": genomes,
        "intermediate": tmp_path / "intermediate",
        "curated_root": tmp_path / "intermediate" / "bakta",
    }


@pytest.fixture()
def forbidden_bakta(monkeypatch):
    """A `run_bakta` that FAILS THE TEST if it is ever called.

    Not a mock that returns something plausible: the whole assertion is that
    this function is never reached, so the honest implementation of "reached" is
    to raise with a message that names the sample.
    """
    calls: list[str] = []

    def exploding_run_bakta(sample, **kwargs):
        calls.append(sample.sample_id)
        raise AssertionError(
            f"Bakta was invoked for {sample.sample_id} under "
            f"reuse_tool_output=require; `require` must never reach the tool"
        )

    monkeypatch.setattr(
        "papipeline.adapters.bakta.run_bakta", exploding_run_bakta
    )
    return calls


def _manifest() -> SampleManifest:
    return SampleManifest(
        samples=[
            Sample(sample_id=sample_id, assembly_path=f"{sample_id}.fna")
            for sample_id in isolate_ids()
        ]
    )


def _provenance_rows(intermediate: Path) -> dict[str, dict[str, str]]:
    """The per-sample reuse decision the stage wrote, keyed by sample_id."""
    from papipeline.io.tsv import read_tsv

    path = Path(intermediate) / "bakta" / stage.REUSE_PROVENANCE_FILENAME
    assert path.is_file(), f"no reuse provenance at {path}"
    return {row["sample_id"]: row for row in read_tsv(path)}


class TestTheCommittedModeIsTheOneUnderTest:
    def test_the_overlay_ships_require(self, smoke_run):
        """Guards the rest of this file against testing a stand-in config.

        If the committed mode stops being `require`, every assertion below
        would still pass - `prefer` with all ten curated also annotates ten -
        while the property they exist to pin (no edge to the tool) would be
        gone. So the mode is asserted, from the accessor, once.
        """
        assert smoke_run["config"].reuse_tool_output() == "require"


class TestCuratedOutputPresentForTheWholeCohort:
    """Direction one: the tool is not reached, and all ten are annotated."""

    def test_no_bakta_invocation_occurs(
        self, smoke_run, forbidden_bakta, caplog
    ):
        for sample_id in isolate_ids():
            stage_curated_output(smoke_run["curated_root"], sample_id)

        with caplog.at_level(logging.INFO, logger="stages.annotation"):
            result = stage.run(
                smoke_run["config"], _manifest(), RunMode.REAL,
                smoke_run["intermediate"], database=smoke_run["database"],
            )

        assert forbidden_bakta == [], (
            "the runner that raises on call was reached for "
            f"{forbidden_bakta}; `require` must not invoke Bakta when curated "
            "output is present"
        )
        # Non-vacuity: the run actually annotated the cohort rather than
        # producing an empty stage, which would also never call the tool.
        assert len(result) == SMOKE_COHORT_SIZE
        assert sum(1 for rows in result.values() if rows) == SMOKE_COHORT_SIZE

    def test_every_sample_is_recorded_as_reused_with_its_digest(
        self, smoke_run, forbidden_bakta
    ):
        """Reuse is auditable after the fact, not just true during the run."""
        for sample_id in isolate_ids():
            stage_curated_output(smoke_run["curated_root"], sample_id)

        stage.run(
            smoke_run["config"], _manifest(), RunMode.REAL,
            smoke_run["intermediate"], database=smoke_run["database"],
        )

        rows = _provenance_rows(smoke_run["intermediate"])
        assert set(rows) == set(isolate_ids())
        assert all(row["reused"] == "true" for row in rows.values()), rows
        assert all(row["reason"] == "verified" for row in rows.values()), rows
        # A reuse with no digest would be a reuse nobody could check later.
        assert all(len(row["sha256_feature_tsv"]) == 64 for row in rows.values())
        assert all(len(row["sha256_gff3"]) == 64 for row in rows.values())

    def test_the_records_reach_the_pipeline_schema(
        self, smoke_run, forbidden_bakta
    ):
        """Reuse is not a side channel: the genes land where the loader reads."""
        for sample_id in isolate_ids():
            stage_curated_output(smoke_run["curated_root"], sample_id)

        stage.run(
            smoke_run["config"], _manifest(), RunMode.REAL,
            smoke_run["intermediate"], database=smoke_run["database"],
        )

        table = (
            Path(smoke_run["intermediate"]) / "annotation"
            / f"{isolate_ids()[0]}.annotation.tsv"
        )
        assert table.is_file(), f"no standardised table at {table}"
        loaded = stage.read_standardised(table)
        assert len(loaded) == 1
        assert loaded[0].gene_name == "oprD"
        assert loaded[0].annotation_source == "bakta"


class TestCuratedOutputAbsentForOneSample:
    """Direction two: the refusal names the sample and the cohort survives."""

    @staticmethod
    def _stage_without_one(smoke_run, absent: str):
        for sample_id in isolate_ids():
            if sample_id == absent:
                continue
            stage_curated_output(smoke_run["curated_root"], sample_id)

    def test_the_run_refuses_and_names_the_sample(
        self, smoke_run, forbidden_bakta, caplog
    ):
        absent = isolate_ids()[3]
        self._stage_without_one(smoke_run, absent)

        with caplog.at_level(logging.WARNING, logger="stages.annotation"):
            result = stage.run(
                smoke_run["config"], _manifest(), RunMode.REAL,
                smoke_run["intermediate"], database=smoke_run["database"],
            )

        assert forbidden_bakta == [], (
            f"`require` reached Bakta for {forbidden_bakta} because curated "
            "output was absent; absent curated output must be a refusal, not a "
            "fallback to the tool"
        )
        # Named, in the one place a reader of a failed run will look.
        text = caplog.text
        assert absent in text, (
            f"the refusal does not name {absent}: {text!r} - a refusal a reader "
            "cannot act on is indistinguishable from a sample that was skipped"
        )
        # ...and in the durable record, not only the log line.
        row = _provenance_rows(smoke_run["intermediate"])[absent]
        assert row["reused"] == "false", row
        assert row["reason"] == stage.REUSE_MISSING, row
        assert row["detail"], "a refusal must say what was looked for"

    def test_the_other_nine_are_unaffected(
        self, smoke_run, forbidden_bakta
    ):
        """One unverifiable isolate must not end the cohort.

        This is the whole reason `require` is usable on a bounded run: the
        cohort is mostly unprepared by construction, so a per-isolate refusal
        is a recording, and a fatal one would make the mode unusable.
        """
        absent = isolate_ids()[3]
        self._stage_without_one(smoke_run, absent)

        result = stage.run(
            smoke_run["config"], _manifest(), RunMode.REAL,
            smoke_run["intermediate"], database=smoke_run["database"],
        )

        # Every member keeps its key: an absent key and an empty list mean
        # different things to every downstream denominator.
        assert set(result) == set(isolate_ids())
        annotated = {k for k, rows in result.items() if rows}
        assert absent not in annotated
        assert annotated == set(isolate_ids()) - {absent}

    def test_absent_curated_output_writes_no_bakta_files(
        self, smoke_run, forbidden_bakta
    ):
        """The refusal leaves no tool output behind, so a rerun can tell."""
        absent = isolate_ids()[3]
        self._stage_without_one(smoke_run, absent)

        stage.run(
            smoke_run["config"], _manifest(), RunMode.REAL,
            smoke_run["intermediate"], database=smoke_run["database"],
        )

        staged = Path(smoke_run["curated_root"]) / absent
        assert not staged.exists() or not any(staged.iterdir()), (
            f"{staged} holds files after a refusal; nothing should have written "
            "there"
        )