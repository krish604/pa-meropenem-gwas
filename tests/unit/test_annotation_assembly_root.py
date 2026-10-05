"""Stage 2 must read the bounded directory, not `data/`.

A smoke overlay authorises ten assemblies. Stage 2 resolved them from
``data_root`` instead of the overlay's genome directory, so it found every
isolate in the 967-member roster that happens to have a sequence in ``data/`` -
**602** of them - and began annotating isolates the overlay never authorised:

    PDT000050773.2, PDT000084290.3, PDT000084304.3 ...

At the measured ~1.5 min/genome that is a ~15-hour run on sixty times the
authorised sample count. Stage 1 got this right all along: `run.py` builds the
smoke cohort from `smoke_genome_dir` and reports `pass=10`, so the log said ten
while stage 2 did six hundred. The two disagreed about what the run was.

The cohort is *deliberately* the whole roster, so an isolate with no prepared
assembly refuses per sample rather than silently disappearing. That design only
holds if the stages needing a sequence read from the bounded directory -
widening the cohort's membership is the manifest's job, and widening the *work*
is nobody's.

**Also fixed by the same change:** the curated copies in `db/smoke_genomes/` used
to be decorative. Stage 2 reached `data/GCA_*/` by way of the manifest's
`assembly` attribute and annotated an incidentally-matching copy from the full
cohort instead - the same sequence, but not the fixture the smoke overlay is
supposed to be exercising. `test_the_curated_copies_are_what_gets_read` pins
that by asserting on the path actually opened.
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
SMOKE = REPO / "config" / "machines" / "smoke.yaml"

AUTHORISED = 10
#: Cohort members with a real, non-empty assembly in `data/` that the overlay
#: does not authorise. Asserted as "> 0" rather than an exact figure, because the
#: figure moves as `data/` is populated; the point is that it is large.
UNAUTHORISED_IN_DATA = 100


@pytest.fixture(scope="module")
def smoke_config():
    return load_config(REPO / "config" / "science.yaml", machine=SMOKE)


@pytest.fixture
def cohort(tmp_path, smoke_config, monkeypatch):
    """A smoke-shaped cohort: 967 members, ten with a curated assembly.

    Two directories, both populated:

    * ``db/smoke_genomes/`` - the ten the overlay authorises;
    * ``data/`` - a further `UNAUTHORISED_IN_DATA` isolates that have sequences
      and must never be read.

    A stage that reached for the wrong one finds real files rather than an empty
    directory, so the bound is tested by what it *reads*, not by what happens to
    be absent.
    """
    authorised_dir = tmp_path / "smoke_genomes"
    authorised_dir.mkdir()
    authorised = []
    for index in range(AUTHORISED):
        sample_id = f"PDT_AUTH{index:04d}"
        (authorised_dir / f"{sample_id}.fna").write_text(
            ">c\nACGTACGTAC\n", encoding="utf-8"
        )
        authorised.append(sample_id)

    data_dir = tmp_path / "data"
    unauthorised = []
    for index in range(UNAUTHORISED_IN_DATA):
        sample_id = f"PDT_UNAUTH{index:04d}"
        accession = f"GCA_{index:09d}.1"
        directory = data_dir / accession
        directory.mkdir(parents=True)
        # Named the way a real assembly is - accession *in the filename*, inside a
        # directory named by accession. `locate_assembly._candidates` matches
        # either the sample id exactly or an accession pattern in the name, so a
        # file called `PDT_UNAUTH0000_genomic.fna` is invisible to it. Getting
        # this wrong makes the unauthorised isolates undiscoverable, the stage
        # skips them for a different reason, and the bound test passes whether or
        # not the fix is present.
        (directory / f"{accession}_genomic.fna").write_text(
            ">c\nACGTACGTAC\n", encoding="utf-8"
        )
        unauthorised.append(Sample(sample_id=sample_id, assembly=accession))

    remaining = [
        Sample(sample_id=f"PDT_NONE{index:04d}")
        for index in range(967 - AUTHORISED - UNAUTHORISED_IN_DATA)
    ]

    monkeypatch.setattr(
        PipelineConfig, "data_root", lambda self, mode: data_dir, raising=True
    )
    monkeypatch.setattr(
        PipelineConfig,
        "assembly_root",
        lambda self, mode: (
            authorised_dir if self.is_smoke_overlay() else self.data_root(mode)
        ),
        raising=True,
    )
    manifest = SampleManifest(
        samples=[
            Sample(sample_id=sid, assembly_path=str(authorised_dir / f"{sid}.fna"))
            for sid in authorised
        ]
        + unauthorised
        + remaining
    )
    return manifest, authorised, unauthorised


@pytest.fixture
def recording_bakta(monkeypatch):
    """Records the assembly path each call is given, and writes valid output."""
    from papipeline.adapters.bakta import (
        BaktaOutputs,
        BaktaResult,
        assembly_for,
    )

    opened: list = []

    def fake(sample, **kwargs):
        assembly = assembly_for(sample, Path(kwargs["genomes_dir"]))
        opened.append(Path(assembly))
        out_dir = Path(kwargs["out_root"]) / sample.sample_id
        out_dir.mkdir(parents=True, exist_ok=True)
        outputs = BaktaOutputs(
            out_dir=out_dir, sample_id=sample.sample_id, genome_stem="s"
        )
        outputs.annotation_tsv.write_text(
            "Sequence Id\tLocus Tag\tStart\tStop\tStrand\tProduct\tQuery\n"
            "c\tL1\t1\t9\t+\tprotein\t1\n",
            encoding="utf-8",
        )
        outputs.gff.write_text(
            "##gff-version 3\nc\tBakta\tCDS\t1\t9\t.\t+\t0\tID=1;Name=L1\n",
            encoding="utf-8",
        )
        return BaktaResult(
            sample_id=sample.sample_id, state="SUCCEEDED", attempts=1,
            outputs=outputs, command=["bakta"],
        )

    monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", fake)
    return opened


def _provisioned_config(tmp_path):
    """The smoke config with a database the preflight accepts.

    **And with `annotation.reuse_tool_output` pinned to `off`.**

    This file's subject is the ASSEMBLY ROOT - which directory stage 2 resolves
    assemblies from - and its observation point is `recording_bakta`, i.e. the
    assembly path handed to the tool. That observation point is only reachable in
    a mode where stage 2 actually looks for an assembly.

    It is not reachable under the smoke overlay's shipped `require`: with curated
    output present, reuse returns before any assembly is read, and with curated
    output absent, `require` refuses before one is. So `require` makes this
    file's subject unobservable, vacuously or otherwise - `test_the_curated_copies_are_what_gets_read`
    failed exactly that way when the mode moved to `require`, asserting
    "no assembly was read at all".

    So the key is pinned here rather than inherited from an overlay whose value
    this test does not care about. `off` is also the STRICTEST setting for this
    bound: every cohort member with a resolvable assembly reaches the tool, so
    "only the ten authorised were attempted" is at its strongest. Reuse is not
    under test here and has its own module,
    `tests/unit/test_smoke_require_never_reaches_bakta.py`.
    """
    import copy

    config = load_config(REPO / "config" / "science.yaml", machine=SMOKE)
    database = REPO / "db" / "bakta_db" / "db-light"
    if not (database / "version.json").is_file():
        pytest.skip("bakta database not provisioned here; db/ is gitignored")
    prepared = copy.copy(config)
    raw = copy.deepcopy(dict(config.raw))
    raw["annotation"] = dict(raw.get("annotation") or {})
    raw["annotation"]["bakta_db"] = str(database)
    object.__setattr__(prepared, "raw", raw)
    return _pin_reuse_mode(prepared, "off")


def _pin_reuse_mode(config, mode: str):
    """Force `config.reuse_tool_output()` to return `mode`.

    `annotation.reuse_tool_output` lives in the raw mapping for a config loaded
    WITHOUT a machine overlay, but the smoke overlay's value is read off
    `MachineConfig`, and that field is what `reuse_tool_output()` consults
    (`loader.py:855`). So the overlay's own field is the thing to set - patching
    `raw` alone would be silently overridden by the overlay.
    """
    import dataclasses

    machine = config.machine
    if machine is None:
        raise AssertionError("this fixture expects a machine overlay to patch")
    return dataclasses.replace(
        config, machine=dataclasses.replace(machine, reuse_tool_output=mode)
    )


class TestTheOverlayBoundIsHonoured:
    def test_only_the_authorised_assemblies_are_attempted(
        self, cohort, recording_bakta, tmp_path
    ):
        manifest, authorised, _ = cohort
        stage.run(
            _provisioned_config(tmp_path), manifest, RunMode.REAL, tmp_path
        )
        assert len(recording_bakta) == AUTHORISED, (
            f"stage 2 attempted {len(recording_bakta)} genomes; the smoke "
            f"overlay authorises {AUTHORISED}"
        )

    def test_nothing_under_data_is_ever_read(
        self, cohort, recording_bakta, tmp_path
    ):
        """The check that matters, on the path actually opened.

        Asserting a count alone would pass if the stage read 10 sequences from
        `data/` instead of the 10 it was meant to. The overlay's bound is about
        *which* directory, not how many.
        """
        manifest, _, _ = cohort
        stage.run(
            _provisioned_config(tmp_path), manifest, RunMode.REAL, tmp_path
        )
        offenders = [p for p in recording_bakta if "smoke_genomes" not in p.parts]
        assert not offenders, (
            f"stage 2 read {len(offenders)} assembly(ies) outside the smoke "
            f"genome directory, e.g. {offenders[:3]}. Those isolates have "
            "sequences in data/ but are not authorised by the overlay."
        )

    def test_the_curated_copies_are_what_gets_read(
        self, cohort, recording_bakta, tmp_path
    ):
        """The fixture stops being decorative.

        Stage 2 used to reach `data/GCA_*/` by way of the manifest's `assembly`
        attribute and annotate an incidentally-matching copy from the full
        cohort. Same sequence, but not the file the smoke overlay exists to
        exercise - so a problem with the curated copy would never have surfaced.
        """
        manifest, authorised, _ = cohort
        stage.run(
            _provisioned_config(tmp_path), manifest, RunMode.REAL, tmp_path
        )
        assert recording_bakta, "no assembly was read at all"
        for path in recording_bakta:
            assert path.parent.name == "smoke_genomes"
            assert path.stem in authorised

    def test_the_unauthorised_isolates_are_recorded_not_annotated(
        self, cohort, recording_bakta, tmp_path
    ):
        manifest, _, unauthorised = cohort
        result = stage.run(
            _provisioned_config(tmp_path), manifest, RunMode.REAL, tmp_path
        )
        assert len(result) == 967
        for sample in unauthorised:
            assert result[sample.sample_id] == [], (
                "an isolate the overlay does not authorise must be recorded as "
                "not-annotated, never annotated because a copy happens to exist"
            )


class TestTheAccessorItself:
    def test_a_smoke_overlay_points_at_its_genome_directory(self, smoke_config):
        assert smoke_config.assembly_root(RunMode.REAL).name == "smoke_genomes"

    def test_an_ordinary_run_is_unaffected(self):
        config = load_config(REPO / "config" / "science.yaml", machine="laptop")
        assert config.assembly_root(RunMode.REAL) == config.data_root(RunMode.REAL)

    def test_it_agrees_with_the_manifest_builder(self, smoke_config):
        """The stage and `discover_run_manifest` must not disagree.

        Stage 1 reported `pass=10` while stage 2 read 602 - the two had separate
        answers to "which assemblies does this run use". One accessor, so the
        next such divergence cannot be introduced by only one of them changing.
        """
        from papipeline.run import discover_run_manifest

        manifest = discover_run_manifest(smoke_config, RunMode.REAL)
        prepared = [
            s for s in manifest.samples if getattr(s, "assembly_path", None)
        ]
        root = smoke_config.assembly_root(RunMode.REAL)
        for sample in prepared:
            assert str(sample.assembly_path).startswith(str(root)), (
                f"{sample.sample_id} has an assembly_path outside {root}, so the "
                "manifest builder and assembly_root disagree"
            )