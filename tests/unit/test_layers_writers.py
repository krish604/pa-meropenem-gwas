"""The layer writers: eight outputs, one isolate set, a logged drop list.

Integration over the committed synthetic fixtures in ``test_data/`` plus small
hand-built bundles in ``tmp_path``. Nothing under ``data/``, ``db/`` or
``PDC_essential.tsv`` is read and no real-mode run happens.

Decisions pinned here:

* **One isolate set, every output.** The writers refuse a bundle whose sample
  set differs from the manifest, naming the output file, rather than emitting
  eight tables that quietly disagree about the cohort.
* **The drop rule applies to every output the same way**, the drop list is
  logged *and* recorded in the provenance banner of ``all_layers.tsv``, and the
  feature dictionary describes exactly the features that survived.
* **The matrix is drug-parameterised**: its name follows
  ``project.primary_antibiotic`` (``gwas_features.target_antibiotic``), so the
  ticket's ``meropenem_feature_matrix.tsv`` appears exactly when the project is
  configured for meropenem - a file named for one drug while carrying another
  drug's phenotype would be a lie.
* **Phenotype is S/I/R only.** The matrix carries the category and no MIC
  column; an isolate with no phenotype record stays missing rather than being
  defaulted.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, List, Sequence

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import DataContractError
from papipeline.io.tsv import read_tsv
from papipeline.layers import encode_layers
from papipeline.layers.features import Feature, LayerBundle
from papipeline.layers.settings import LayerSettings, layer_settings
from papipeline.layers.writers import (
    ALL_LAYERS_FILENAME,
    DICTIONARY_FILENAME,
    FLAG_COLUMN,
    apply_drop,
    feature_matrix_filename,
    write_layers,
)
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode

PIPELINE_ROOT = Path(__file__).resolve().parents[2]

LAYER_FILENAMES = {
    1: "layer1_acquired.tsv",
    2: "layer2_pdc.tsv",
    3: "layer3_oprD.tsv",
    4: "layer4_efflux.tsv",
    5: "layer5_targets.tsv",
}


def _header_line(path: Path) -> List[str]:
    """The header row of a written TSV, skipping the provenance banner."""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        return line.split("\t")
    raise AssertionError(f"{path} has no header row")


def _banner(path: Path) -> List[str]:
    return [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("#")
    ]


def _config_under(root: Path, *, replace: tuple[str, str] | None = None) -> PipelineConfig:
    root = Path(root)
    shutil.copytree(PIPELINE_ROOT / "config", root / "config")
    science = root / "config" / "science.yaml"
    text = science.read_text(encoding="utf-8")
    if replace is not None:
        old, new = replace
        assert old in text, f"config/science.yaml no longer carries {old!r}"
        text = text.replace(old, new, 1)
    science.write_text(text, encoding="utf-8")
    return load_config(science, machine=None)


def _run(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
    phenotype_dir: Path,
    output_dir: Path,
):
    return encode_layers(
        config,
        manifest,
        mode=RunMode.TEST,
        input_root=intermediate_root,
        output_dir=output_dir,
        phenotype_dir=phenotype_dir,
    )


class TestTheCommittedFixtureRun:
    """Encode the 20 synthetic isolates once and audit every output."""

    @pytest.fixture()
    def result(self, config, manifest, intermediate_root, phenotype_dir, tmp_path):
        return _run(config, manifest, intermediate_root, phenotype_dir, tmp_path / "layers")

    @pytest.fixture()
    def expected_ids(self, manifest) -> set:
        return set(manifest.sample_ids)

    def test_all_eight_outputs_exist(self, result):
        assert set(result.paths) == {
            "layer1",
            "layer2",
            "layer3",
            "layer4",
            "layer5",
            "all_layers",
            "matrix",
            "dictionary",
        }
        names = {path.name for path in result.paths.values()}
        assert names == {
            LAYER_FILENAMES[1],
            LAYER_FILENAMES[2],
            LAYER_FILENAMES[3],
            LAYER_FILENAMES[4],
            LAYER_FILENAMES[5],
            ALL_LAYERS_FILENAME,
            feature_matrix_filename("imipenem"),
            DICTIONARY_FILENAME,
        }
        assert all(path.is_file() for path in result.paths.values())

    def test_every_isolate_table_covers_the_same_isolate_set(
        self, result, expected_ids
    ):
        # The dictionary is a feature table, not an isolate table; it is
        # checked against the feature set below.
        for name, path in result.paths.items():
            if name == "dictionary":
                continue
            rows = read_tsv(path)
            ids = [row["sample_id"] for row in rows]
            assert set(ids) == expected_ids, name
            assert len(ids) == len(set(ids)) == len(expected_ids), name

    def test_the_fixture_leaves_exactly_the_two_oprd_features(
        self, result, expected_ids
    ):
        # 32 features are built from this fixture (4 L1 + 7 L2 + 4 L3 + 12 L4
        # + 5 L5); only the two oprD tier features reach 5 carriers in 20
        # isolates, so the other 30 are dropped and must be listed as such.
        assert {feature.name for feature in result.features} == {
            "L3_oprD_LoF_tier1",
            "L3_oprD_off_any",
        }
        assert len(result.dropped) == 30
        assert len(result.features) + len(result.dropped) == 32

    def test_the_dropped_list_names_the_features_and_their_carriers(self, result):
        by_name = {drop.name: drop for drop in result.dropped}
        assert by_name["L1_KPC"].carriers == 2
        assert by_name["L3_oprD_absent"].carriers == 1
        assert "layers.min_carriers" in by_name["L1_KPC"].reason
        assert by_name["L5_gyrA"].carriers == 0
        assert by_name["L4_pump_high_MexXY"].layer == 4

    def test_the_drop_list_is_logged(self, config, manifest, intermediate_root,
                                     phenotype_dir, tmp_path, caplog):
        import logging

        with caplog.at_level(logging.WARNING, logger="papipeline.layers"):
            _run(config, manifest, intermediate_root, phenotype_dir, tmp_path / "layers")
        assert "L1_KPC" in caplog.text
        assert "dropped" in caplog.text.lower()

    def test_the_drop_list_is_recorded_in_the_banner(self, result):
        banner = "\n".join(_banner(result.paths["all_layers"]))
        assert "drop_list" in banner
        assert "L1_KPC" in banner
        assert "layers.min_carriers" in banner

    def test_each_layer_file_holds_only_its_own_prefix(self, result):
        for layer in (1, 2, 4, 5):
            columns = _header_line(result.paths[f"layer{layer}"])
            assert columns[0] == "sample_id"
            assert all(
                column.startswith(f"L{layer}_") for column in columns[1:]
            ), (layer, columns)
        columns = _header_line(result.paths["layer3"])
        assert columns == [
            "sample_id",
            "L3_oprD_LoF_tier1",
            "L3_oprD_off_any",
        ]

    def test_the_matrix_carries_the_phenotype_and_no_mic(self, result, phenotype_dir):
        columns = _header_line(result.paths["matrix"])
        assert "phenotype" in columns
        assert "MIC" not in columns
        assert "MIC_unit" not in columns
        assert FLAG_COLUMN in columns

        expected = {
            row["sample_id"]: row["phenotype"]
            for row in read_tsv(phenotype_dir / "imipenem_phenotype.tsv")
        }
        for row in read_tsv(result.paths["matrix"]):
            want = expected.get(row["sample_id"])
            assert row["phenotype"] == want

    def test_partial_call_flag_is_missing_not_zero_when_nothing_tracks_it(
        self, result
    ):
        rows = read_tsv(result.paths["matrix"])
        assert {row[FLAG_COLUMN] for row in rows} == {None}
        banner = "\n".join(_banner(result.paths["all_layers"]))
        assert FLAG_COLUMN in banner
        assert "call_state" in banner

    def test_the_dictionary_describes_exactly_the_surviving_features(self, result):
        rows = read_tsv(result.paths["dictionary"])
        assert [row["feature"] for row in rows] == [
            "L3_oprD_LoF_tier1",
            "L3_oprD_off_any",
        ]
        assert _header_line(result.paths["dictionary"]) == [
            "feature",
            "layer",
            "source_tokens",
            "rule",
        ]
        for row in rows:
            assert row["layer"] == "3"
            assert row["source_tokens"]
            assert row["rule"]

    def test_the_dictionary_matches_the_all_layers_columns(self, result):
        columns = _header_line(result.paths["all_layers"])
        assert FLAG_COLUMN in columns
        feature_columns = [
            column for column in columns if column not in ("sample_id", FLAG_COLUMN)
        ]
        assert feature_columns == [row["feature"] for row in read_tsv(result.paths["dictionary"])]

    def test_the_run_is_byte_stable(self, config, manifest, intermediate_root,
                                    phenotype_dir, tmp_path):
        first = _run(config, manifest, intermediate_root, phenotype_dir, tmp_path / "a")
        second = _run(config, manifest, intermediate_root, phenotype_dir, tmp_path / "b")
        for key, path in first.paths.items():
            assert path.read_bytes() == second.paths[key].read_bytes(), key


class TestTheMatrixNameFollowsTheConfiguredDrug:
    def test_primary_antibiotic_decides_the_file_name(
        self, manifest, intermediate_root, tmp_path
    ):
        config = _config_under(
            tmp_path,
            replace=("primary_antibiotic: imipenem", "primary_antibiotic: meropenem"),
        )
        phenotype_dir = tmp_path / "phenotype"
        phenotype_dir.mkdir(parents=True)
        (phenotype_dir / "meropenem_phenotype.tsv").write_text(
            "sample_id\tantibiotic\tphenotype\n"
            + "".join(
                f"{sample_id}\tmeropenem\tR\n" for sample_id in manifest.sample_ids
            ),
            encoding="utf-8",
        )
        result = encode_layers(
            config,
            manifest,
            mode=RunMode.TEST,
            input_root=intermediate_root,
            output_dir=tmp_path / "layers",
            phenotype_dir=phenotype_dir,
        )
        assert result.paths["matrix"].name == "meropenem_feature_matrix.tsv"
        assert result.antibiotic == "meropenem"
        rows = read_tsv(result.paths["matrix"])
        assert {row["phenotype"] for row in rows} == {"R"}


class TestTheDropRule:
    @staticmethod
    def _bundle() -> LayerBundle:
        sample_ids = tuple(f"S{i}" for i in range(1, 7))
        features = [
            Feature(name="L1_too_rare", layer=1, rule="r", source_tokens=("x",)),
            Feature(name="L1_at_the_floor", layer=1, rule="r", source_tokens=("x",)),
            Feature(name="L1_everywhere", layer=1, rule="r", source_tokens=("x",)),
            Feature(
                name="L5_qrdr_count",
                layer=5,
                value_kind="count",
                rule="r",
                source_tokens=("x",),
            ),
        ]
        values = {
            "L1_too_rare": {s: int(s in {"S1", "S2", "S3", "S4"}) for s in sample_ids},
            "L1_at_the_floor": {
                s: int(s in {"S1", "S2", "S3", "S4", "S5"}) for s in sample_ids
            },
            "L1_everywhere": {s: 1 for s in sample_ids},
            "L5_qrdr_count": {s: int(s == "S1") for s in sample_ids},
        }
        return LayerBundle(sample_ids=sample_ids, features=features, values=values)

    def test_the_boundaries_are_exactly_the_configured_numbers(self):
        kept, dropped = apply_drop(self._bundle(), LayerSettings(5, 0.98))
        assert [feature.name for feature in kept] == ["L1_at_the_floor"]
        reasons = {drop.name: drop for drop in dropped}
        assert reasons["L1_too_rare"].carriers == 4
        assert "layers.min_carriers" in reasons["L1_too_rare"].reason
        assert reasons["L1_everywhere"].carriers == 6
        assert "layers.max_prevalence" in reasons["L1_everywhere"].reason
        assert reasons["L5_qrdr_count"].carriers == 1

    def test_the_settings_come_from_the_configuration(self, config):
        assert layer_settings(config) == LayerSettings(min_carriers=5, max_prevalence=0.98)


class TestTheIsolateSetGuard:
    @pytest.fixture()
    def partial_bundle(self, manifest) -> LayerBundle:
        sample_ids = tuple(manifest.sample_ids)[:-1]
        feature = Feature(name="L1_KPC", layer=1, rule="r", source_tokens=("x",))
        return LayerBundle(
            sample_ids=sample_ids,
            features=[feature],
            values={"L1_KPC": {sample_id: 1 for sample_id in sample_ids}},
        )

    def test_a_bundle_missing_an_isolate_is_refused_naming_the_file(
        self, manifest, partial_bundle, tmp_path
    ):
        with pytest.raises(DataContractError) as excinfo:
            write_layers(
                partial_bundle,
                manifest=manifest,
                phenotype_by_sample={},
                partial_flags={},
                antibiotic="imipenem",
                settings=LayerSettings(5, 0.98),
                mode=RunMode.TEST,
                output_dir=tmp_path / "layers",
            )
        message = str(excinfo.value)
        assert LAYER_FILENAMES[1] in message
        assert set(manifest.sample_ids) - set(partial_bundle.sample_ids)

    def test_a_bundle_with_an_extra_isolate_is_refused(
        self, manifest, tmp_path
    ):
        sample_ids = tuple(manifest.sample_ids) + ("TEST_PA_999",)
        feature = Feature(name="L1_KPC", layer=1, rule="r", source_tokens=("x",))
        bundle = LayerBundle(
            sample_ids=sample_ids,
            features=[feature],
            values={"L1_KPC": {sample_id: 1 for sample_id in sample_ids}},
        )
        with pytest.raises(DataContractError) as excinfo:
            write_layers(
                bundle,
                manifest=manifest,
                phenotype_by_sample={},
                partial_flags={},
                antibiotic="imipenem",
                settings=LayerSettings(5, 0.98),
                mode=RunMode.TEST,
                output_dir=tmp_path / "layers",
            )
        assert LAYER_FILENAMES[1] in str(excinfo.value)
