"""The wiring test for the layer feature tables (written first).

`papipeline/layers` shipped complete and unit-tested (`tests/unit/test_layers_*`)
and was called by nothing: no rule, no stage branch, no entry point. It is one
of the five dormant packages this integration pass exists to wire. This file
pins the call site, so the wiring cannot be deleted without a red test.

What the wiring is, in the words of `.build/layer-encoder.registry-notes.md`:

* the layer tables run inside the **`pangenome`** branch — the same place stage
  12's other feature table (`gwas_features`) is produced, because that is where
  the stage-12 feature input is produced today, and every source it reads
  (stage 4's determinant table, stage 6's regulator table, stage 7's gene
  matrix) is written by a stage that runs before stage 7;
* **REAL and TEST run it, STUB never does.** STUB fabricates declared outputs
  and runs no parser, and `encode_layers` reads four source tables and a
  phenotype table, none of which a stub run has;
* `input_root` is `tool_output_root` (in TEST that is `test_data/intermediate`,
  where the four source tables are committed fixtures) and `output_dir` is the
  **run's own results tree**, never `test_data/` — the fixtures are read-only,
  so a run must not write an eight-file directory into the directory the next
  run reads from.

The drop rule is asserted rather than assumed: on the committed twenty-isolate
fixtures the encoder builds 32 features and keeps 2 (the registry notes record
the same two survivors), which is what proves the wiring fed it the real
fixtures instead of an empty directory.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.io import read_tsv
from papipeline.run import run_pipeline

LAYER_FILES = {
    "all_layers.tsv",
    "feature_dictionary.tsv",
    "imipenem_feature_matrix.tsv",
    "layer1_acquired.tsv",
    "layer2_pdc.tsv",
    "layer3_oprD.tsv",
    "layer4_efflux.tsv",
    "layer5_targets.tsv",
}

#: What the committed fixtures leave standing after `layers.min_carriers`.
SURVIVORS = {"L3_oprD_LoF_tier1", "L3_oprD_off_any"}


@pytest.fixture(scope="module")
def test_run(pipeline_root: Path):
    """One full TEST run, shared by every assertion below."""
    return run_pipeline(
        config_path=pipeline_root / "config" / "science.yaml",
        mode="TEST",
        config=load_config(pipeline_root / "config" / "science.yaml"),
    )


@pytest.fixture(scope="module")
def stub_run(pipeline_root: Path):
    """One full STUB run, for the mode that must not build them."""
    return run_pipeline(
        config_path=pipeline_root / "config" / "science.yaml",
        mode="STUB",
        config=load_config(pipeline_root / "config" / "science.yaml"),
    )


class TestALayersCallSiteExists:
    def test_a_test_run_records_the_layer_matrix(self, test_run, pipeline_root):
        """The call site ran: `layers` is a declared output of the run."""
        assert "layers" in test_run.outputs, (
            "encode_layers is not called by any stage; the package is dormant "
            "again"
        )
        path = Path(test_run.outputs["layers"])
        assert path.is_file(), path
        assert path.name == "imipenem_feature_matrix.tsv", path

    def test_all_eight_outputs_land_together(self, test_run):
        layer_dir = Path(test_run.outputs["layers"]).parent
        written = {p.name for p in layer_dir.iterdir()}
        assert written == LAYER_FILES, (
            f"missing: {sorted(LAYER_FILES - written)}; "
            f"unexpected: {sorted(written - LAYER_FILES)}"
        )

    def test_the_pangenome_status_is_completed(self, test_run):
        """The branch the wiring hangs off is the one that actually ran."""
        assert test_run.stage_status.get("pangenome") == "completed"


class TestWhatWasBuilt:
    def test_the_drop_rule_ran_on_the_committed_fixtures(self, test_run):
        """32 features built, 2 kept — real fixtures, not an empty directory.

        If the wiring ever handed the encoder an empty `input_root`, no
        feature would be built and this would be an empty dictionary instead.
        """
        layer_dir = Path(test_run.outputs["layers"]).parent
        rows = read_tsv(
            layer_dir / "feature_dictionary.tsv",
            required_columns=["feature", "layer", "source_tokens", "rule"],
            unique_columns=["feature"],
        )
        assert {r["feature"] for r in rows} == SURVIVORS, (
            "the committed fixtures' two survivors changed — either a fixture "
            "moved or the wiring stopped reading the committed source tables"
        )

    def test_the_banner_records_the_survival_count(self, test_run):
        """The count travels with the file, so a reader sees what was dropped."""
        text = Path(test_run.outputs["layers"]).read_text()
        assert "# features: 2 kept of 32 built" in text, (
            "the survival banner is missing or has changed:\n"
            + "\n".join(line for line in text.splitlines() if line.startswith("#"))
        )

    def test_the_matrix_has_the_manifest_cohort_and_no_one_else(
        self, test_run, manifest
    ):
        """AGENTS.md rule 5: exactly the manifest's samples, no fuzzy join."""
        rows = read_tsv(
            Path(test_run.outputs["layers"]),
            required_columns=["sample_id", "phenotype"],
            unique_columns=["sample_id"],
        )
        assert {r["sample_id"] for r in rows} == set(manifest.sample_ids)

    def test_the_phenotype_column_is_last(self, test_run):
        header = next(
            line
            for line in Path(test_run.outputs["layers"]).read_text().splitlines()
            if line and not line.startswith("#")
        ).split("\t")
        assert header[-1] == "phenotype", header[-3:]
        assert header[0] == "sample_id", header[:3]


class TestWhereItWasWritten:
    def test_the_outputs_live_in_the_run_tree(self, test_run, pipeline_root):
        """`output_dir` is the run's own intermediate tree, not `test_data/`."""
        path = Path(test_run.outputs["layers"])
        expected = (
            pipeline_root
            / "results"
            / "test"
            / "intermediate"
            / "stages"
            / "layers"
            / "imipenem_feature_matrix.tsv"
        )
        assert path == expected, path

    def test_nothing_is_written_under_the_committed_fixtures(
        self, test_run, pipeline_root
    ):
        """The fixtures are inputs; a run that writes into them poisons the
        next run (and dirties a byte-stable tree, see docs/reproducibility.md).
        """
        fixture_root = pipeline_root / "test_data"
        assert not (fixture_root / "layers").exists()
        stray = sorted(p.name for p in fixture_root.rglob("*layers*"))
        assert stray == [], f"layers output found under test_data/: {stray}"


class TestStubNeverBuildsThem:
    def test_a_stub_run_has_no_layer_output(self, stub_run):
        assert "layers" not in stub_run.outputs, (
            "STUB fabricates declared outputs and runs no parser; "
            "encode_layers reads five tables and must not be called here"
        )

    def test_a_stub_run_leaves_no_layer_directory(self, stub_run, pipeline_root):
        layer_dir = pipeline_root / "results" / "stub" / "intermediate" / "stages" / "layers"
        assert not layer_dir.exists(), layer_dir
