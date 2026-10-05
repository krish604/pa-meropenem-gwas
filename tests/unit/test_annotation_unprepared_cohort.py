"""Stage 2 must survive a cohort that is mostly unprepared.

A smoke run's manifest is deliberately the **whole** PDC roster - 967 isolates -
with only ten carrying a prepared assembly. `discover_run_manifest` says so
explicitly: every isolate stays a cohort member and "refuses per sample later,
rather than silently disappearing". Later turned out to mean stages 1, 3 and 6,
which each record an unusable isolate and carry on. Stage 2 was the exception, so
the refusal never happened where it was documented to:

* `run_bakta` **raises** `PipelineError` from `locate_assembly` before it
  returns, so the loop's `result.state != "SUCCEEDED"` check never saw it;
* the stage then raised `did not succeed for every genome` anyway.

Either way one absent genome halted the cohort. A real end-to-end smoke run
died on `No assembly found for PDT000152708.2` - an isolate that was never
supposed to be annotated - having successfully annotated several others first.

**The two failure modes are recorded differently, though both mean
not-annotated.** "No assembly" and "Bakta crashed on a real assembly" are
different things: the first is a cohort-composition fact, the second is a tool or
data fault that may deserve attention. Collapsing them would mean a run in which
Bakta failed on every genome it was actually given looked the same as a bounded
run that correctly skipped 957. The existing `failed` list already encodes
`sample=reason`, so both go there - the loop just stops treating the reason as
fatal.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import PipelineError
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import annotation as stage

REPO = Path(__file__).resolve().parents[2]


ANNOTATION_TSV = (
    "# Annotated with Bakta\n"
    "Sequence Id\tLocus Tag\tStart\tStop\tStrand\tProduct\tQuery\n"
    "contig_00001\tTESTPA_00001\t1\t300\t+\thypothetical protein\t1\n"
)

#: `ID=` carries the locus tag, as real Bakta writes it - verified against
#: `results/real/intermediate/bakta/*/*.gff3`, where ID and the TSV's `Locus Tag`
#: are the same string. This fixture used `ID=1`, a bare number that matches
#: nothing, and it passed only because the old cross-check compared row *counts*
#: rather than identities. Containment compares identities, so the fixture had to
#: become realistic rather than the check being loosened back.
GFF_BODY = (
    "##gff-version 3\n"
    "contig_00001\tBakta\tCDS\t1\t300\t.\t+\t0\t"
    "ID=TESTPA_00001;Name=TESTPA_00001\n"
)


@pytest.fixture
def provisioned_config():
    """The real config, with a minimal database that satisfies the preflight."""
    import copy

    config = load_config(REPO / "config" / "science.yaml")
    database = REPO / "db" / "bakta_db" / "db-light"
    if not (database / "version.json").is_file():
        pytest.skip("bakta database not provisioned here; db/ is gitignored")
    prepared = copy.copy(config)
    raw = copy.deepcopy(dict(config.raw))
    raw["annotation"] = dict(raw.get("annotation") or {})
    raw["annotation"]["bakta_db"] = str(database)
    object.__setattr__(prepared, "raw", raw)
    return prepared


def manifest_of(prepared_ids, cohort_size=967):
    """A manifest shaped like the real smoke one: mostly without assemblies.

    Distinct prefixes, because a manifest is keyed on `sample_id` and duplicate
    ids are refused - the unprepared run and the prepared run must not overlap.
    """
    unprepared = [f"PDT_U{i:06d}" for i in range(cohort_size - len(prepared_ids))]
    return SampleManifest(
        samples=[
            Sample(sample_id=sample_id, assembly_path=f"{sample_id}.fna")
            if sample_id in prepared_ids
            else Sample(sample_id=sample_id)
            for sample_id in unprepared + list(prepared_ids)
        ]
    )


def prepared_ids(n):
    return [f"PDT_P{i:06d}" for i in range(n)]


@pytest.fixture
def cohort_root(tmp_path, monkeypatch):
    """A stand-in data root holding files for the prepared isolates only.

    ``data/`` is read-only (AGENTS.md rule 5 territory, and the user's standing
    instruction), so the cohort's sequences are written under tmp_path and
    ``data_root`` is pointed at it. Declaring an ``assembly_path`` is not enough
    - ``locate_assembly`` checks the file exists, which is the whole point of it
    refusing rather than guessing.
    """
    root = tmp_path / "data"

    def _use(ids):
        for sample_id in ids:
            directory = root / sample_id
            directory.mkdir(parents=True, exist_ok=True)
            (directory / f"{sample_id}.fna").write_text(
                ">contig_00001\nACGTACGTAC\n", encoding="utf-8"
            )
        # Patched on the class, not the instance: PipelineConfig is a frozen
        # dataclass, so an instance attribute cannot be set. monkeypatch undoes
        # it, and only this test's own stage call reads data_root.
        monkeypatch.setattr(
            PipelineConfig, "data_root", lambda self, mode: root, raising=True
        )
        return root

    return _use


@pytest.fixture
def stub_bakta(monkeypatch):
    """Replace `run_bakta` so no real annotation runs in a unit test.

    Succeeds for a sample with an assembly, and records what it was asked for so
    the test can assert the loop reached every prepared isolate.
    """
    from papipeline.adapters.bakta import BaktaOutputs, BaktaResult, assembly_for

    asked: list = []

    def fake_run_bakta(sample, **kwargs):
        # Resolve the assembly first, exactly as the real adapter does. Without
        # this the stub would "succeed" for isolates that have no sequence, and
        # the loop's missing-assembly branch would never be reached - a stub
        # more permissive than the thing it stands in for tests nothing.
        assembly_for(sample, Path(kwargs["genomes_dir"]))
        asked.append(sample.sample_id)
        out_dir = Path(kwargs["out_root"]) / sample.sample_id
        out_dir.mkdir(parents=True, exist_ok=True)
        outputs = BaktaOutputs(
            out_dir=out_dir, sample_id=sample.sample_id, genome_stem="s"
        )
        for path in (outputs.annotation_tsv, outputs.gff):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                ANNOTATION_TSV if path.suffix == ".tsv" else GFF_BODY,
                encoding="utf-8",
            )
        return BaktaResult(
            sample_id=sample.sample_id, state="SUCCEEDED", attempts=1,
            outputs=outputs, command=["bakta"],
        )

    monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", fake_run_bakta)
    return asked


class TestTheLoopToleratesAnUnpreparedCohort:
    def test_a_mostly_unprepared_manifest_completes(
        self, provisioned_config, stub_bakta, tmp_path, cohort_root
    ):
        """The end-to-end failure this exists to prevent.

        967 cohort members, 10 with a prepared assembly. The stage must finish
        and return a record for every isolate, with the 957 unprepared ones
        present and empty - the representation `load_from_intermediate` already
        uses, so downstream sees "not annotated" rather than a missing key.
        """
        prepared = prepared_ids(10)
        cohort_root(prepared)
        result = stage.run(
            provisioned_config, manifest_of(prepared), RunMode.REAL, tmp_path
        )
        assert len(result) == 967
        assert sorted(stub_bakta) == sorted(prepared), (
            "the loop must attempt exactly the isolates that have an assembly"
        )
        assert sum(1 for v in result.values() if v) == 10

    def test_every_isolate_is_present_not_only_the_prepared_ones(
        self, provisioned_config, stub_bakta, tmp_path, cohort_root
    ):
        """An absent key and an empty list mean different things downstream.

        Dropping the 957 would shrink every denominator the report computes,
        which is the specific harm `stages.variants` documents for its own
        merge. The cohort stays whole; the annotation is what is missing.
        """
        prepared = prepared_ids(10)
        cohort_root(prepared)
        manifest = manifest_of(prepared)
        result = stage.run(provisioned_config, manifest, RunMode.REAL, tmp_path)
        assert set(result) == set(manifest.sample_ids)
        assert all(isinstance(v, list) for v in result.values())


class TestTheTwoFailureModesAreToldApart:
    def test_a_missing_assembly_is_not_attributed_to_bakta(
        self, provisioned_config, stub_bakta, tmp_path, caplog, cohort_root
    ):
        """"No assembly" and "Bakta crashed" must not read alike.

        A cohort in which Bakta failed on every genome it was actually given is
        a very different event from a bounded run that correctly skipped 957
        isolates, and the log is where a reader would tell them apart.
        """
        import logging

        prepared = prepared_ids(3)
        cohort_root(prepared)
        with caplog.at_level(logging.WARNING, logger="stages.annotation"):
            stage.run(provisioned_config, manifest_of(prepared), RunMode.REAL, tmp_path)
        text = caplog.text
        assert "no_assembly" in text
        assert "no_assembly" not in stub_bakta, (
            "an isolate with no assembly must never reach the tool"
        )

    def test_a_bakta_failure_is_recorded_as_its_own_state(
        self, provisioned_config, monkeypatch, tmp_path, caplog, cohort_root
    ):
        """A present assembly whose annotation fails is a different reason."""
        import logging

        from papipeline.adapters.bakta import BaktaOutputs, BaktaResult

        def failing(sample, **kwargs):
            from papipeline.adapters.bakta import assembly_for

            assembly_for(sample, Path(kwargs["genomes_dir"]))
            out_dir = Path(kwargs["out_root"]) / sample.sample_id
            out_dir.mkdir(parents=True, exist_ok=True)
            return BaktaResult(
                sample_id=sample.sample_id, state="FAILED", attempts=3,
                outputs=BaktaOutputs(
                    out_dir=out_dir, sample_id=sample.sample_id, genome_stem="s"
                ),
                command=["bakta"], validation_detail="exit 1",
            )

        monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", failing)
        prepared = prepared_ids(2)
        cohort_root(prepared)
        with caplog.at_level(logging.WARNING, logger="stages.annotation"):
            result = stage.run(
                provisioned_config, manifest_of(prepared), RunMode.REAL, tmp_path
            )
        # Asserted on the counts, not the sample list: the warning names only
        # the first ten, which for this cohort are all `no_assembly`, so a
        # tool failure would be invisible by name while still being counted.
        text = caplog.text
        assert "965 no assembly" in text, text
        assert "2 annotated but unusable" in text, (
            "a tool failure must be counted separately from a missing assembly, "
            "or a run where Bakta failed on every genome it was given looks the "
            "same as a bounded run that correctly skipped the rest"
        )
        assert all(v == [] for v in result.values())


class TestAPreflightFailureStillStopsTheStage:
    def test_a_missing_database_is_still_fatal(
        self, provisioned_config, stub_bakta, tmp_path
    ):
        """Tolerance for absent assemblies must not extend to a broken database.

        A misconfigured database is a run-level fault, not a per-sample one:
        every genome would fail the same way, and continuing would turn one
        clear refusal into 967 recorded failures.
        """
        import copy

        config = copy.copy(provisioned_config)
        raw = copy.deepcopy(dict(provisioned_config.raw))
        raw["annotation"] = dict(raw["annotation"])
        raw["annotation"]["bakta_db"] = str(tmp_path / "no-such-db")
        object.__setattr__(config, "raw", raw)

        with pytest.raises(PipelineError):
            stage.run(config, manifest_of(prepared_ids(1)), RunMode.REAL, tmp_path)


class TestThePreflightRunsOnceNotOncePerGenome:
    def test_the_database_is_verified_once_for_the_cohort(
        self, provisioned_config, stub_bakta, tmp_path, monkeypatch, cohort_root
    ):
        """It sat inside the loop, so a 967-isolate cohort verified the database
        967 times - re-reading references.tsv each time to reach the same answer.
        """
        from papipeline.adapters import bakta_db as preflight

        calls: list = []
        real = preflight.preflight_database

        def counting(**kwargs):
            calls.append(kwargs)
            return real(**kwargs)

        monkeypatch.setattr(preflight, "preflight_database", counting)
        cohort_root(prepared_ids(10))
        stage.run(
            provisioned_config,
            manifest_of(prepared_ids(10)),
            RunMode.REAL,
            tmp_path,
        )
        assert len(calls) == 1, f"preflight ran {len(calls)} times for one cohort"
