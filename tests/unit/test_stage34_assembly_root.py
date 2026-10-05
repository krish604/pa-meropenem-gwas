"""Stages 3 and 4 must read the bounded directory, not `data/`.

Stage 2 had this bug and was fixed in `83d6aba`. The same defect was latent in
stage 3 and stage 4, and in two different shapes:

* both stages took `data_root` as a parameter, and `run.py` supplied
  `config.data_root(resolved)` - so a smoke overlay pointed them at `data/` and
  at 602 assemblies it does not authorise;
* neither call site passed `data_root` at all, so on a REAL run both raised
  "needs a data_root". Stage 3's refusal would have been the *second* thing to
  stop the run; the first was stage 2's unbound annotation.

Neither had been reached, because the run never got past annotation - which is
why a defect this reachable survived. Fixing a stage without checking its
neighbours is how the same bug appears three times.

**The decoy filenames here are real.** The stage 2 test named its unauthorised
assemblies `PDT_UNAUTH0000_genomic.fna` and passed against the *buggy* code,
because `locate_assembly._candidates` matches either the sample id exactly or an
accession pattern in the filename - and `PDT_...` matched neither, so the decoys
were invisible and the stage skipped those isolates for an unrelated reason. A
bound test whose decoys cannot be found proves no bound.

So the fixture writes what NCBI actually writes: `data/GCA_000937495.2/
GCA_000937495.2_ASM93749v2_genomic.fna`. `test_the_decoys_are_actually_findable`
asserts that directly, so this mistake cannot be repeated silently.

**Which tests actually guard the regression.** Only the two in
`TestTheOrchestratorSuppliesTheBoundedRoot`. They were verified to fail against
the unfixed `run.py` and pass against the fixed one.

The eight stage-level tests pass either way, and that is by construction rather
than by weakness: each one *injects* the authorised directory as `data_root`, so
they verify the property that matters locally - a stage reads only from the root
it is handed - but they cannot detect `run.py` handing it the wrong root. Both
halves are kept: the stage tests stop a stage from quietly reaching for
`config.data_root` itself, and the orchestrator tests stop `run.py` from
supplying `data/`. Reading the eight as regression guards for the run.py bug
would be a mistake, which is the same category of mistake as the invisible decoys.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import PipelineError
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import amr as stage_amr
from papipeline.stages import mlst as stage_mlst

REPO = Path(__file__).resolve().parents[2]
SMOKE = REPO / "config" / "machines" / "smoke.yaml"

AUTHORISED = 10
UNAUTHORISED = 60

#: Exactly the shape NCBI writes, from the real tree: an accession-named
#: directory containing `<accession>_<label>_genomic.fna`.
UNAUTHORISED_ACCESSIONS = tuple(
    f"GCA_0009{37495 + index:05d}.2" for index in range(UNAUTHORISED)
)


@pytest.fixture(scope="module")
def smoke_config():
    return load_config(REPO / "config" / "science.yaml", machine=SMOKE)


@pytest.fixture
def cohort(tmp_path, monkeypatch):
    """A smoke-shaped cohort with *findable* decoys.

    Returns the manifest, the authorised ids, and the unauthorised accessions.
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
    unauthorised_samples = []
    for index, accession in enumerate(UNAUTHORISED_ACCESSIONS):
        directory = data_dir / accession
        directory.mkdir(parents=True)
        (directory / f"{accession}_ASM_test_genomic.fna").write_text(
            ">c\nACGTACGTAC\n", encoding="utf-8"
        )
        unauthorised_samples.append(
            Sample(sample_id=f"PDT_UNAUTH{index:04d}", assembly=accession)
        )

    filler = [
        Sample(sample_id=f"PDT_NONE{i:04d}")
        for i in range(967 - AUTHORISED - UNAUTHORISED)
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
    return (
        SampleManifest(
            samples=[
                Sample(
                    sample_id=sid,
                    assembly_path=str(authorised_dir / f"{sid}.fna"),
                )
                for sid in authorised
            ]
            + unauthorised_samples
            + filler
        ),
        authorised,
        list(UNAUTHORISED_ACCESSIONS),
    )


class TestTheDecoysAreActuallyFindable:
    """Guards the guard.

    The stage 2 test passed against buggy code because its decoys could not be
    found. This asserts the decoys resolve, so a bound test cannot pass for that
    reason again.
    """

    def test_a_decoy_resolves_through_locate_assembly(self, cohort, tmp_path):
        from papipeline.assemblies import locate_assembly

        manifest, _, accessions = cohort
        sample = next(
            s for s in manifest.samples if getattr(s, "assembly", None) == accessions[0]
        )
        resolved = locate_assembly(sample, tmp_path / "data")
        assert resolved.is_file()
        assert accessions[0] in resolved.name, (
            "the decoy filename must carry its accession; a file named after the "
            "isolate is invisible to locate_assembly and proves nothing"
        )

    def test_an_authorised_assembly_also_resolves(self, cohort, tmp_path):
        from papipeline.assemblies import locate_assembly

        manifest, authorised, _ = cohort
        sample = next(s for s in manifest.samples if s.sample_id == authorised[0])
        assert locate_assembly(sample, tmp_path / "smoke_genomes").is_file()


class TestStage3IsBounded:
    def test_only_authorised_isolates_are_called(self, cohort, monkeypatch, tmp_path):
        from papipeline.adapters import mlst as mlst_tool

        manifest, authorised, _ = cohort
        opened: list = []

        def fake(sample_id, *, assembly, **kwargs):
            opened.append(Path(assembly))
            return mlst_tool.MlstResult(
                sample_id=sample_id, st="235", alleles={"acsA": "1"}, found=True
            )

        monkeypatch.setattr(mlst_tool, "call_isolate", fake)
        calls = stage_mlst.run(
            _config(tmp_path), manifest, RunMode.REAL, tmp_path,
            data_root=_authorised_dir(tmp_path),
        )
        assert len(opened) == AUTHORISED
        assert len(calls) == 967

    def test_nothing_under_data_is_read(self, cohort, monkeypatch, tmp_path):
        from papipeline.adapters import mlst as mlst_tool

        manifest, _, _ = cohort
        opened: list = []
        monkeypatch.setattr(
            mlst_tool, "call_isolate",
            lambda sample_id, *, assembly, **kw: (
                opened.append(Path(assembly)),
                mlst_tool.MlstResult(sample_id=sample_id, st="1", found=True),
            )[1],
        )
        stage_mlst.run(
            _config(tmp_path), manifest, RunMode.REAL, tmp_path,
            data_root=_authorised_dir(tmp_path),
        )
        assert not [p for p in opened if "smoke_genomes" not in p.parts]

    def test_it_refuses_a_real_run_with_no_data_root(self, tmp_path):
        """The second half of the bug, still worth pinning.

        The call site passed no `data_root`, so a REAL run raised instead of
        running. Better to refuse than to fall back to `data/` silently.
        """
        manifest = SampleManifest(samples=[Sample(sample_id="PDT_A")])
        with pytest.raises(PipelineError):
            stage_mlst.run(_config(tmp_path), manifest, RunMode.REAL, tmp_path)


class TestStage4IsBounded:
    def test_only_authorised_isolates_are_screened(self, cohort, monkeypatch, tmp_path):
        manifest, authorised, _ = cohort
        screened: list = []

        class Backend:
            name = "stub"

            def detect(self, sample_id, assembly):
                screened.append(Path(assembly))
                return []

            def provenance(self):
                return {}

        monkeypatch.setattr(
            stage_amr, "build_amr_adapter", lambda *a, **k: Backend()
        )
        result = stage_amr.run(
            _config(tmp_path), manifest, RunMode.REAL, tmp_path, "imipenem",
            data_root=_authorised_dir(tmp_path),
        )
        assert len(screened) == AUTHORISED
        assert len(result) == 967

    def test_nothing_under_data_is_screened(self, cohort, monkeypatch, tmp_path):
        manifest, _, _ = cohort
        screened: list = []

        class Backend:
            name = "stub"

            def detect(self, sample_id, assembly):
                screened.append(Path(assembly))
                return []

            def provenance(self):
                return {}

        monkeypatch.setattr(
            stage_amr, "build_amr_adapter", lambda *a, **k: Backend()
        )
        stage_amr.run(
            _config(tmp_path), manifest, RunMode.REAL, tmp_path, "imipenem",
            data_root=_authorised_dir(tmp_path),
        )
        assert not [p for p in screened if "smoke_genomes" not in p.parts]

    def test_it_refuses_a_real_run_with_no_data_root(self, tmp_path):
        manifest = SampleManifest(samples=[Sample(sample_id="PDT_A")])
        with pytest.raises(PipelineError):
            stage_amr.run(
                _config(tmp_path), manifest, RunMode.REAL, tmp_path, "imipenem"
            )


class TestTheOrchestratorSuppliesTheBoundedRoot:
    def test_run_pipeline_resolves_assemblies_through_assembly_root(self):
        """The fix has to be at the call site, not only in the stages.

        Both stages accept `data_root` as a parameter, so a stage-level default
        would have left `run.py` free to keep supplying `data/`.
        """
        source = (REPO / "papipeline" / "run.py").read_text(encoding="utf-8")
        assert "data_root = config.assembly_root(resolved)" in source, (
            "run.py still resolves assemblies from data_root, so every stage it "
            "feeds reads the unbounded directory"
        )
        assert "data_root = config.data_root(resolved)" not in source

    def test_both_stages_receive_it(self):
        source = (REPO / "papipeline" / "run.py").read_text(encoding="utf-8")
        for stage in ("stage_mlst.run", "stage_amr.run"):
            index = source.index(stage)
            window = source[index: index + 260]
            assert "data_root=data_root" in window, (
                f"{stage} is called without data_root, so a REAL run raises "
                "rather than running"
            )


def _authorised_dir(tmp_path: Path) -> Path:
    return tmp_path / "smoke_genomes"


def _config(tmp_path: Path) -> PipelineConfig:
    config = load_config(REPO / "config" / "science.yaml", machine=SMOKE)
    return config