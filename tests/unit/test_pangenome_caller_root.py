"""Stage 7's caller passes the wrong root, and the empty pan-genome hides it.

**The bug.** `papipeline/run.py` computed two roots side by side -
`intermediate = config.intermediate_root(resolved)` (line 617) and
`source_dir = config.tool_output_root(resolved)` (line 618) - and handed every
stage except one `source_dir`. The pangenome stage got `intermediate`.

The two are not interchangeable, and in TEST they are two different
directories:

    intermediate_root(TEST) = results/test/intermediate   <- this run's OUTPUTS
    tool_output_root(TEST)  = test_data/intermediate      <- the committed fixtures

`config/lint`... `loader.tool_output_root` exists for exactly this reason
(loader.py:537-539): "Keeping this distinct from `intermediate_root` is what
stops a run from reading its own outputs as if they were inputs."

So the stage was pointed at a directory this run was writing into, asked for
annotations there, found none, and returned an empty partition - `0 genes | 0
core | 0 accessory` - with a success exit and a well-formed
`pangenome_summary.tsv`. Nothing refused. The committed fixture says 10 / 5 / 5.

**Why the fallback path is the one that bites.** `run.py` also passes
`annotations or None`, and in a full run the annotation stage has already
populated that, which makes the root argument inert. The wrong root only becomes
observable when the argument is consulted: `only=["pangenome"]` (it has no
prerequisites in `PREREQUISITES`), or any run where annotation produced nothing.
That is why these tests drive the caller rather than the stage - a stage-level
test passes today and would not have caught this.

**The empty result was the second half of the failure.** Even reading the right
root, a cohort with no genes in it is not a finding - it is a run that looked in
the wrong place. `docs/scientific_rules.md` §9: "An empty input file is a
`DataContractError`, so a truncated or failed run is never mistaken for a cohort
with no findings." Zero genes now refuses in both modes.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import panaroo as adapter
from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.errors import PipelineError, StageError
from papipeline.models import RunMode
from papipeline.run import run_pipeline
from papipeline.stages import pangenome as stage

REPO = Path(__file__).resolve().parents[2]

#: The committed fixture's own numbers, read rather than restated. If the
#: fixture is regenerated these follow it; a hard-coded 10/5/5 would instead
#: fail and teach nothing about the code.
FIXTURE_SUMMARY = REPO / "test_data" / "intermediate" / "pangenome" / "pangenome_summary.tsv"


def _fixture_metrics() -> dict:
    metrics = {}
    for line in FIXTURE_SUMMARY.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        metric, _, value = line.partition("\t")
        metrics[metric.strip()] = value.strip()
    return metrics


def _written_metrics(path: Path) -> dict:
    metrics = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        metric, _, value = line.partition("\t")
        metrics[metric.strip()] = value.strip()
    return metrics


@pytest.fixture()
def redirected(monkeypatch, tmp_path):
    """Point the run's outputs at tmp_path so the repo tree is never written."""
    root = tmp_path / "results"
    monkeypatch.setenv(RESULTS_ROOT_ENV, str(root))
    return root


class TestTheCallerHandsOverTheToolOutputRoot:
    """The stage must be told where its INPUTS are, not where its outputs go."""

    def test_it_is_given_the_tool_output_root(self, config, monkeypatch, redirected):
        seen = {}

        def spy(cfg, manifest, mode, root, annotations=None):
            seen["root"] = root
            return stage.build_from_annotations({}, ())

        monkeypatch.setattr("papipeline.run.stage_pangenome.run", spy)
        run_pipeline(config=config, mode="TEST", only=["pangenome"], write_html=False)

        assert seen, "the pangenome stage was never called"
        assert Path(seen["root"]) == config.tool_output_root(RunMode.TEST), (
            f"stage 7 was handed {seen['root']}, but its inputs are the fixtures "
            f"under {config.tool_output_root(RunMode.TEST)}"
        )

    def test_it_is_not_given_the_run_intermediate_root(self, config, monkeypatch, redirected):
        """Named explicitly, because in REAL the two are the same path.

        In TEST they differ and the wrong one resolves to an empty directory, so
        the bug is loud there. In REAL `tool_output_root` *is*
        `intermediate_root`, so an assertion that only checked the TEST value
        would pass on a caller that was wrong in both modes. This asserts the
        run's output directory is not what gets handed over.
        """
        seen = {}

        def spy(cfg, manifest, mode, root, annotations=None):
            seen["root"] = root
            return stage.build_from_annotations({}, ())

        monkeypatch.setattr("papipeline.run.stage_pangenome.run", spy)
        run_pipeline(config=config, mode="TEST", only=["pangenome"], write_html=False)

        assert Path(seen["root"]) != config.intermediate_root(RunMode.TEST), (
            "stage 7 was handed this run's own output directory, so it reads its "
            "own outputs as inputs"
        )


class TestTestModeReportsTheCommittedFixture:
    """The end-to-end consequence: 10 / 5 / 5, not 0 / 0 / 0."""

    def test_the_counts_match_the_fixture(self, config, redirected):
        expected = _fixture_metrics()
        result = run_pipeline(
            config=config, mode="TEST", only=["pangenome"], write_html=False
        )
        written = _written_metrics(result.outputs["pangenome"])

        assert written["n_genes_total"] == expected["n_genes_total"], (
            f"TEST reported {written['n_genes_total']} genes against a fixture "
            f"declaring {expected['n_genes_total']}"
        )
        assert written["n_core_genes"] == expected["n_core_genes"]
        assert written["n_accessory_genes"] == expected["n_accessory_genes"]
        assert (written["n_genes_total"], written["n_core_genes"], written["n_accessory_genes"]) == (
            "10",
            "5",
            "5",
        ), "the committed fixture is 10 genes / 5 core / 5 accessory"

    def test_zero_genes_would_have_been_written_silently(self, config, redirected):
        """Why this is a bug and not a shrug: 0/0/0 is indistinguishable from a result.

        `write_outputs` emits a well-formed summary for an empty partition. Before
        the fix that file was produced and the run reported success, so nothing
        downstream could tell an empty cohort from a wrong directory.
        """
        empty = stage.build_from_annotations({}, ("S1", "S2"))
        paths = stage.write_outputs(empty, redirected / "empty")
        metrics = _written_metrics(paths["pangenome_summary"])
        assert metrics["n_genes_total"] == "0"
        assert paths["pangenome_summary"].exists(), (
            "an empty partition still produces a valid-looking summary - which is "
            "the whole reason it has to be refused upstream"
        )


class TestAnEmptyPangenomeIsRefusedInBothModes:
    """scientific_rules §9: never mistaken for a cohort with no findings."""

    def test_test_mode_refuses_rather_than_reporting_zero(self, config, tmp_path):
        empty_root = tmp_path / "no-annotations"
        empty_root.mkdir()

        from papipeline.manifest import discover_manifest

        manifest = discover_manifest(config.metadata_dir(RunMode.TEST))

        with pytest.raises(StageError) as excinfo:
            stage.run(config, manifest, RunMode.TEST, empty_root, None)

        message = str(excinfo.value)
        assert "0 genes" in message.lower(), message
        assert "missing input" in message.lower(), (
            f"the refusal must distinguish this from a cohort with no findings: {message}"
        )
        assert "gene_name" in message, (
            "the message must say what to check, not only that the count is zero: "
            f"{message}"
        )
        assert str(empty_root) in message, (
            "the message must name the directory it actually looked in, so the "
            f"operator can see which root was wrong: {message}"
        )

    def test_real_mode_refuses_rather_than_reporting_zero(
        self, monkeypatch, tmp_path
    ):
        """A header-only panaroo table parses cleanly and yields zero genes.

        `parse_presence_csv` refuses a wholly empty file, but a file with the
        right header and no rows is a well-formed table describing an empty
        cohort - which is exactly what a truncated or mis-rooted panaroo run
        leaves behind.
        """
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        monkeypatch.setenv("PIPELINE_ALLOW_REAL_MODE", "1")
        monkeypatch.setattr(adapter, "preflight", lambda: None)
        config = load_config(REPO / "config" / "science.yaml", machine="laptop")

        cohort = ("S1", "S2")
        manifest = SampleManifest(
            samples=[Sample(sample_id=s, assembly_path=f"/nowhere/{s}.fna") for s in cohort]
        )
        panaroo_dir = tmp_path / adapter.OUTPUT_DIRNAME
        panaroo_dir.mkdir(parents=True)
        header = ",".join(["Gene", "Non-unique Gene name", "Annotation", *cohort])
        (panaroo_dir / adapter.PRESENCE_CSV).write_text(header + "\n", encoding="utf-8")

        with pytest.raises(StageError) as excinfo:
            stage.run(config, manifest, RunMode.REAL, tmp_path)

        message = str(excinfo.value)
        assert "0 genes" in message.lower(), message
        assert "missing input" in message.lower(), message
        assert "panaroo" in message, (
            f"a REAL refusal must point at the table panaroo wrote: {message}"
        )

    def test_the_refusal_is_a_pipeline_error(self, config, tmp_path):
        """So an operator's `except PipelineError` still catches it."""
        from papipeline.manifest import discover_manifest

        manifest = discover_manifest(config.metadata_dir(RunMode.TEST))
        with pytest.raises(PipelineError):
            stage.run(config, manifest, RunMode.TEST, tmp_path / "absent", None)