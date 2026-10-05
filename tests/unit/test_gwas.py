"""Tests for GWAS input construction and statistics (stage 12).

The critical behaviours: I is excluded rather than folded into R or S,
multiple-testing correction is applied, and a feature confined to one
lineage is flagged as lineage-linked rather than reported as a resistance
association.
"""

from __future__ import annotations

import copy

import pytest

from papipeline.errors import DataContractError
from papipeline.models import GwasResult, Phenotype, PhenotypeCall
from papipeline.stages import gwas as stage


def _call(sample_id, value, mic=None):
    return PhenotypeCall(
        sample_id=sample_id,
        antibiotic="imipenem",
        phenotype=Phenotype(value),
        mic=mic,
    )


def _features(rows: str):
    header = "sample_id\tgene__a\tsnp__b\tunitig__c\n"
    parsed = []
    for line in (header + rows).strip().splitlines()[1:]:
        fields = line.split("\t")
        parsed.append(
            {
                "sample_id": fields[0],
                "gene__a": fields[1],
                "snp__b": fields[2],
                "unitig__c": fields[3],
            }
        )
    return parsed


def _cohort(n_pos=6, n_neg=6, excluded=()):
    """A cohort with n_pos R and n_neg S samples, plus optional exclusions.

    Returns ``(calls, feature_rows, manifest)``. The manifest is built from
    the same IDs so ``build_input`` sees the cohort it is given rather than
    the synthetic-fixture manifest.
    """
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    calls, rows, ids = [], [], []
    for i in range(n_pos):
        sid = f"T{i:03d}"
        calls.append(_call(sid, "R"))
        rows.append(f"{sid}\t1\t0\t0")
        ids.append(sid)
    for i in range(n_neg):
        sid = f"S{i:03d}"
        calls.append(_call(sid, "S"))
        rows.append(f"{sid}\t0\t1\t0")
        ids.append(sid)
    for j, value in enumerate(excluded):
        sid = f"X{j:03d}"
        calls.append(_call(sid, value))
        rows.append(f"{sid}\t0\t0\t0")
        ids.append(sid)
    manifest = SampleManifest([Sample(sid, None, "test") for sid in ids])
    return calls, _features("\n".join(rows) + "\n"), manifest


class TestParseFeatureName:
    def test_splits_type_and_label(self):
        assert stage.parse_feature_name("snp__POS_123") == ("snp", "POS_123")

    def test_unprefixed_feature(self):
        assert stage.parse_feature_name("POS_123") == ("unspecified", "POS_123")

    def test_label_may_contain_underscores(self):
        assert stage.parse_feature_name("unitig__utg_000001_len_500") == (
            "unitig",
            "utg_000001_len_500",
        )


class TestBuildInput:
    def test_binarises_correctly(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert gwas_input.n_positive == 6
        assert gwas_input.n_negative == 6

    def test_intermediate_categories_are_excluded(self, config):
        """I, SDD and ND are excluded, never reassigned to R or S."""
        calls, features, cohort = _cohort(n_pos=6, n_neg=6, excluded=("I", "SDD", "ND"))
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert gwas_input.n_excluded == 3
        assert set(gwas_input.sample_ids).isdisjoint({"X000", "X001", "X002"})

    def test_excluded_samples_are_recorded(self, config):
        calls, features, cohort = _cohort(n_pos=6, n_neg=6, excluded=("I",))
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert gwas_input.excluded_samples == ("X000",)

    def test_sample_without_phenotype_is_excluded(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert all(
            sid in {c.sample_id for c in calls} for sid in gwas_input.sample_ids
        )

    def test_too_few_in_a_group_raises(self, config, sample_manifest):
        calls, features, cohort = _cohort(n_pos=1, n_neg=8)
        with pytest.raises(DataContractError, match="Too few samples"):
            stage.build_input(config, cohort, calls, features)

    def test_duplicate_feature_row_raises(self, config):
        calls, features, cohort = _cohort()
        features.append(dict(features[0]))
        with pytest.raises(DataContractError, match="Duplicate sample"):
            stage.build_input(config, cohort, calls, features)

    def test_feature_row_without_sample_id_raises(self, config):
        calls, features, cohort = _cohort()
        features[0]["sample_id"] = ""
        with pytest.raises(DataContractError, match="no sample_id"):
            stage.build_input(config, cohort, calls, features)

    def test_non_binary_value_treated_as_absent(self, config):
        calls, features, cohort = _cohort()
        features[0]["gene__a"] = "maybe"
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert gwas_input.features["gene__a"]["T000"] == 0

    def test_lineages_are_carried(self, config):
        calls, features, cohort = _cohort()
        lineages = {c.sample_id: "L1" for c in calls}
        gwas_input = stage.build_input(
            config, cohort, calls, features, lineages=lineages
        )
        assert gwas_input.lineages["T000"] == "L1"


class TestTheAntibioticIsSelectedByName:
    """The project's antibiotic, not the first one in the list.

    `build_input` read `config.antibiotics[0]`, which is a *position*, not a
    choice. `science.yaml` declares `project.primary_antibiotic` for exactly
    this, and `run.py` already reads it with a fallback - so two places in one
    codebase disagreed about which antibiotic a run was about, and no
    end-to-end test could have caught it while the list held one entry.

    Which is why this is latent rather than active: with a single antibiotic
    configured, `antibiotics[0]` and `primary_antibiotic` coincide. So every
    test here **reorders** the list to put a different antibiotic first, which
    is the only way the difference is observable. A test written against the
    real single-antibiotic config would pass against the bug.

    The symptom is narrow but not harmless: the name reaches a
    `DataContractError` raised when an outcome group is too small, so a failed
    run reports the wrong drug as the one under analysis.
    """

    @pytest.fixture(autouse=True)
    def _restore_config(self, config):
        """Undo the reordering after each test.

        The `config` fixture is session-scoped, so a mutation here would
        otherwise persist into every later test in the run - and did, until
        this fixture was added. Restoring both the raw mapping and the captured
        tuple is necessary because the code reads each.
        """
        original_raw = copy.deepcopy(config.raw)
        original_tuple = config.antibiotics
        yield
        config.raw.clear()
        config.raw.update(original_raw)
        object.__setattr__(config, "antibiotics", original_tuple)

    @staticmethod
    def _reordered(config, *, first, primary):
        """A config whose antibiotic *order* disagrees with the project's choice.

        Both halves have to move. `PipelineConfig.antibiotics` is a tuple
        captured at load time, so editing `config.raw["antibiotics"]` alone
        leaves `antibiotics[0]` unchanged and the test would pass for the wrong
        reason - which is how the first attempt at this test came to be
        vacuous.
        """
        order = [first, "imipenem", "meropenem"]
        config.raw["antibiotics"] = list(order)
        config.raw["project"]["primary_antibiotic"] = primary
        object.__setattr__(config, "antibiotics", tuple(order))
        return config

    def test_the_reported_antibiotic_is_the_projects_not_the_first(
        self, config, sample_manifest
    ):
        """Reordered so the first entry is *not* the primary antibiotic.

        Fails against the positional pick, which would report `meropenem`.
        """
        self._reordered(config, first="meropenem", primary="imipenem")
        calls, features, cohort = _cohort(n_pos=1, n_neg=8)
        with pytest.raises(DataContractError) as excinfo:
            stage.build_input(config, cohort, calls, features)
        assert excinfo.value.context["antibiotic"] == "imipenem", (
            f"the error names {excinfo.value.context['antibiotic']!r}, which is "
            "the first antibiotic in the list rather than the project's "
            "primary_antibiotic"
        )

    def test_it_follows_the_configured_primary_rather_than_a_hardcoded_name(
        self, config, sample_manifest
    ):
        """Any primary must be honoured, not just the one this project uses.

        A fix that hardcoded `imipenem` would pass the test above and be wrong
        for the next project.
        """
        self._reordered(config, first="imipenem", primary="meropenem")
        calls, features, cohort = _cohort(n_pos=1, n_neg=8)
        with pytest.raises(DataContractError) as excinfo:
            stage.build_input(config, cohort, calls, features)
        assert excinfo.value.context["antibiotic"] == "meropenem"

    def test_the_single_antibiotic_case_is_unaffected(self, config, sample_manifest):
        """The real config: one antibiotic, so order cannot matter.

        Asserted so a fix for the reordered case cannot regress the ordinary one.
        """
        calls, features, cohort = _cohort(n_pos=1, n_neg=8)
        with pytest.raises(DataContractError) as excinfo:
            stage.build_input(config, cohort, calls, features)
        assert excinfo.value.context["antibiotic"] == "imipenem"

    def test_a_successful_run_is_unaffected_by_the_change(
        self, config, sample_manifest
    ):
        """The fix touches only the error path, so the happy path still works."""
        calls, features, cohort = _cohort(n_pos=6, n_neg=6)
        result = stage.build_input(config, cohort, calls, features)
        assert result.n_samples == len(cohort.sample_ids)


class TestTestableFeatures:
    def test_feature_present_in_all_samples_is_skipped(self, config):
        """A feature with no contrast carries no information."""
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        calls = [_call(f"T{i:03d}", "R") for i in range(6)] + [
            _call(f"S{i:03d}", "S") for i in range(6)
        ]
        manifest = SampleManifest([Sample(c.sample_id, None, "t") for c in calls])
        rows = "\n".join(f"{c.sample_id}\t1\t0\t0" for c in calls)
        gwas_input = stage.build_input(
            config, manifest, calls, _features(rows + "\n")
        )
        assert "gene__a" not in gwas_input.testable_features()

    def test_feature_absent_everywhere_is_skipped(self, config):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        calls = [_call(f"T{i:03d}", "R") for i in range(6)] + [
            _call(f"S{i:03d}", "S") for i in range(6)
        ]
        manifest = SampleManifest([Sample(c.sample_id, None, "t") for c in calls])
        rows = "\n".join(f"{c.sample_id}\t0\t0\t0" for c in calls)
        gwas_input = stage.build_input(
            config, manifest, calls, _features(rows + "\n")
        )
        assert gwas_input.testable_features() == []


class TestBenjaminiHochberg:
    def test_monotonic_and_bounded(self):
        adjusted = stage.benjamini_hochberg(
            {"a": 0.001, "b": 0.01, "c": 0.5, "d": 0.9}
        )
        assert all(0 <= v <= 1 for v in adjusted.values())
        assert adjusted["a"] <= adjusted["b"] <= adjusted["c"] <= adjusted["d"]

    def test_adjusted_is_never_below_raw(self):
        raw = {"a": 0.01, "b": 0.2, "c": 0.6}
        adjusted = stage.benjamini_hochberg(raw)
        for key, value in raw.items():
            assert adjusted[key] >= value - 1e-12

    def test_single_test_is_unchanged(self):
        assert stage.benjamini_hochberg({"a": 0.03})["a"] == pytest.approx(0.03)

    def test_empty_input(self):
        assert stage.benjamini_hochberg({}) == {}

    def test_known_worked_example(self):
        """p = [0.01, 0.02, 0.03, 0.04, 0.05] under BH."""
        raw = {str(i): p for i, p in enumerate([0.01, 0.02, 0.03, 0.04, 0.05])}
        adjusted = stage.benjamini_hochberg(raw)
        assert adjusted["0"] == pytest.approx(0.05)
        assert adjusted["4"] == pytest.approx(0.05)

    def test_multiple_testing_penalty_is_real(self):
        """A small p-value among many large ones is penalised.

        Identical p-values correctly incur no penalty, because they all rank
        last; the penalty comes from the number of tests ranked above.
        """
        many = stage.benjamini_hochberg(
            {"small": 0.001, **{str(i): 0.5 for i in range(49)}}
        )
        alone = stage.benjamini_hochberg({"small": 0.001})
        assert many["small"] == pytest.approx(0.05)
        assert many["small"] > alone["small"]

    def test_identical_p_values_incur_no_penalty(self):
        """All-equal p-values all rank last, so BH leaves them unchanged."""
        adjusted = stage.benjamini_hochberg({str(i): 0.02 for i in range(20)})
        assert all(v == pytest.approx(0.02) for v in adjusted.values())


class TestReferenceEngine:
    def test_finds_a_perfectly_separating_feature(self, config):
        calls, features, cohort = _cohort()
        lineages = {c.sample_id: f"L{i % 3}" for i, c in enumerate(calls)}
        gwas_input = stage.build_input(
            config, cohort, calls, features, lineages=lineages
        )
        results = stage.ReferenceEngine().run(gwas_input, config)
        by_feature = {r.feature: r for r in results}
        assert by_feature["gene__a"].p_value < 0.05

    def test_multiple_testing_correction_is_applied(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        results = stage.ReferenceEngine().run(gwas_input, config)
        for result in results:
            assert result.adjusted_p_value is not None
            assert result.adjusted_p_value >= result.p_value - 1e-12

    def test_model_records_the_missing_kinship_correction(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        results = stage.ReferenceEngine().run(gwas_input, config)
        assert all("NO_KINSHIP_CORRECTION" in r.model for r in results)

    def test_lineage_distribution_is_reported(self, config):
        calls, features, cohort = _cohort()
        lineages = {c.sample_id: f"L{i % 3}" for i, c in enumerate(calls)}
        gwas_input = stage.build_input(
            config, cohort, calls, features, lineages=lineages
        )
        results = stage.ReferenceEngine().run(gwas_input, config)
        assert all(r.lineage_distribution for r in results)

    def test_results_sorted_by_adjusted_p(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        results = stage.ReferenceEngine().run(gwas_input, config)
        values = [r.adjusted_p_value for r in results]
        assert values == sorted(values)


class TestLineageConfounding:
    def _result(self, distribution):
        return GwasResult(
            feature="unitig__x",
            feature_type="unitig",
            effect=1.0,
            p_value=0.001,
            adjusted_p_value=0.01,
            effect_size=1.0,
            frequency=0.3,
            lineage_distribution=distribution,
        )

    def test_single_lineage_feature_is_confounded(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        result = self._result({"L1": 6})
        assert stage.is_lineage_confounded(result, gwas_input, 0.9)

    def test_evenly_spread_feature_is_not_confounded(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        result = self._result({"L1": 2, "L2": 2, "L3": 2})
        assert not stage.is_lineage_confounded(result, gwas_input, 0.9)

    def test_no_lineage_info_is_not_confounded(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        result = self._result({})
        assert not stage.is_lineage_confounded(result, gwas_input, 0.9)

    def test_flagging_selects_confounded_features(self, config):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        results = [
            self._result({"L1": 6}),
            self._result({"L1": 2, "L2": 2, "L3": 2}),
        ]
        flagged = stage.flag_lineage_linked(results, gwas_input, config)
        assert len(flagged) == 1


class TestSignificant:
    def test_respects_threshold(self, config):
        results = [
            GwasResult("a", "snp", 1.0, 0.001, 0.01, 1.0, 0.5),
            GwasResult("b", "snp", 1.0, 0.4, 0.6, 1.0, 0.5),
            GwasResult("c", "snp", 1.0, 0.2, None, 1.0, 0.5),
        ]
        hits = stage.significant(results, config)
        assert [r.feature for r in hits] == ["a"]

    def test_missing_adjusted_p_is_never_significant(self, config):
        results = [GwasResult("a", "snp", 1.0, 1e-10, None, 1.0, 0.5)]
        assert stage.significant(results, config) == []


class TestRun:
    def test_runs_on_synthetic_fixtures(self, config, manifest, intermediate_root, phenotype_dir, phylogeny_dir):
        from papipeline.models import RunMode
        from papipeline.stages import phenotype as pheno_stage
        from papipeline.stages import phylogeny as phylo_stage

        calls = pheno_stage.load_phenotype(
            config, phenotype_dir, "imipenem", sample_ids=manifest.sample_ids
        )
        lineages = phylo_stage.load_tree_metadata(phylogeny_dir / "tree_metadata.tsv")
        results, gwas_input = stage.run(
            config, manifest, RunMode.TEST, intermediate_root, calls, lineages=lineages
        )
        assert results
        assert gwas_input.n_positive >= config.gwas.min_samples_per_group

    def test_latent_signal_is_recovered(self, config, manifest, intermediate_root, phenotype_dir, phylogeny_dir):
        """The generator plants a recoverable signal; the machinery must find it.

        This tests the machinery, not biology: the signal is a generator
        artefact, documented in papipeline/testing/synthetic.py.
        """
        from papipeline.models import RunMode
        from papipeline.stages import phenotype as pheno_stage
        from papipeline.stages import phylogeny as phylo_stage

        calls = pheno_stage.load_phenotype(
            config, phenotype_dir, "imipenem", sample_ids=manifest.sample_ids
        )
        lineages = phylo_stage.load_tree_metadata(phylogeny_dir / "tree_metadata.tsv")
        results, _ = stage.run(
            config, manifest, RunMode.TEST, intermediate_root, calls, lineages=lineages
        )
        hits = stage.significant(results, config)
        assert hits, "no feature passed the adjusted-p threshold"
        assert hits[0].feature == "snp__TESTPOS_00100"

    def test_lineage_features_are_flagged(self, config, manifest, intermediate_root, phenotype_dir, phylogeny_dir):
        from papipeline.models import RunMode
        from papipeline.stages import phenotype as pheno_stage
        from papipeline.stages import phylogeny as phylo_stage

        calls = pheno_stage.load_phenotype(
            config, phenotype_dir, "imipenem", sample_ids=manifest.sample_ids
        )
        lineages = phylo_stage.load_tree_metadata(phylogeny_dir / "tree_metadata.tsv")
        results, gwas_input = stage.run(
            config, manifest, RunMode.TEST, intermediate_root, calls, lineages=lineages
        )
        flagged = stage.flag_lineage_linked(results, gwas_input, config)
        assert flagged, "the planted lineage-confounded features were not flagged"
        assert all(
            r.feature.startswith(("unitig__", "kmer__")) for r in flagged
        )


class TestPyseerEngine:
    def test_writes_expected_input_files(self, config, tmp_path):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        engine = stage.PyseerEngine("pyseer", tmp_path)
        paths = engine.write_inputs(gwas_input)
        assert set(paths) == {"phenotype", "features", "lineages"}
        for path in paths.values():
            assert path.exists()

    def test_phenotype_file_is_binary(self, config, tmp_path):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        paths = stage.PyseerEngine("pyseer", tmp_path).write_inputs(gwas_input)
        content = paths["phenotype"].read_text()
        assert "\t0" in content or "\t1" in content

    def test_missing_output_returns_empty(self, config, tmp_path):
        calls, features, cohort = _cohort()
        gwas_input = stage.build_input(config, cohort, calls, features)
        assert stage.parse_pyseer_output(tmp_path / "absent.tsv", gwas_input) == []


class TestOptionalAdapters:
    def test_card_rgi_adapter_refuses_when_not_enabled(self):
        from papipeline.stages.amr import CardRgiAdapter

        adapter = CardRgiAdapter("CARD", "2023-01", "imipenem")
        with pytest.raises(NotImplementedError, match="not enabled"):
            adapter.detect("S1", None)
