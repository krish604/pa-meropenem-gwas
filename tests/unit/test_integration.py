"""Tests for convergence (13), co-occurrence (14) and integration (15)."""

from __future__ import annotations

import pytest

from papipeline.models import (
    ClaimStatus,
    ConvergenceCategory,
    GwasResult,
    MASTER_TABLE_COLUMNS,
    MasterRecord,
    Phenotype,
    PhenotypeCall,
    StructuralCallStatus,
    StructuralVariant,
)
from papipeline.stages import cooccurrence as cooc
from papipeline.stages import convergence as conv
from papipeline.stages import integration as integ


# --------------------------------------------------------------------------
# convergence
# --------------------------------------------------------------------------


def _carriers(name, samples):
    return conv.DeterminantCarriers(determinant=name, carriers=tuple(sorted(samples)))


class TestLineageDistribution:
    def test_counts_per_lineage(self):
        distribution = conv.lineage_distribution(
            ["a", "b", "c"], {"a": "L1", "b": "L1", "c": "L2"}
        )
        assert distribution == {"L1": 2, "L2": 1}

    def test_unknown_lineage_recorded(self):
        assert conv.lineage_distribution(["a"], {}) == {"unknown": 1}


class TestClassify:
    def test_widespread_background(self, config):
        record = _carriers("geneX", [f"s{i}" for i in range(18)])
        distribution = {"L1": 6, "L2": 6, "L3": 6}
        assert (
            conv.classify(record, distribution, config, n_samples=20)
            is ConvergenceCategory.WIDESPREAD_BACKGROUND
        )

    def test_rare_isolated(self, config):
        record = _carriers("geneX", ["s1", "s2"])
        assert (
            conv.classify(record, {"L1": 2}, config, n_samples=20)
            is ConvergenceCategory.RARE_ISOLATED
        )

    def test_single_is_rare_not_lineage_associated(self, config):
        """One carrier is too few to distinguish chance from convergence."""
        record = _carriers("geneX", ["s1"])
        assert (
            conv.classify(record, {"L1": 1}, config, n_samples=20)
            is ConvergenceCategory.RARE_ISOLATED
        )

    def test_lineage_associated(self, config):
        record = _carriers("geneX", [f"s{i}" for i in range(6)])
        assert (
            conv.classify(record, {"L1": 6}, config, n_samples=20)
            is ConvergenceCategory.LINEAGE_ASSOCIATED
        )

    def test_recurrent_convergent(self, config):
        record = _carriers("geneX", [f"s{i}" for i in range(6)])
        assert (
            conv.classify(record, {"L1": 3, "L2": 3}, config, n_samples=20)
            is ConvergenceCategory.RECURRENT_CONVERGENT
        )

    def test_unknown_without_lineage_information(self, config):
        record = _carriers("geneX", [f"s{i}" for i in range(6)])
        assert (
            conv.classify(record, {"unknown": 6}, config, n_samples=20)
            is ConvergenceCategory.UNKNOWN
        )


class TestCollectCarriers:
    def test_inverts_sample_to_determinant(self):
        carriers = conv.collect_carriers({"s1": ["gA"], "s2": ["gA", "gB"]})
        by_name = {c.determinant: c.carriers for c in carriers}
        assert by_name["gA"] == ("s1", "s2")
        assert by_name["gB"] == ("s2",)

    def test_duplicates_within_a_sample_collapse(self):
        carriers = conv.collect_carriers({"s1": ["gA", "gA"]})
        assert carriers[0].n_carriers == 1


class TestSupportedAssociations:
    def test_requires_both_significance_and_convergence(self, config):
        calls = [
            conv.ConvergenceCall(
                "gA", 3, 6, {"L1": 3, "L2": 3}, ConvergenceCategory.RECURRENT_CONVERGENT
            ),
            conv.ConvergenceCall(
                "gB", 1, 6, {"L1": 6}, ConvergenceCategory.LINEAGE_ASSOCIATED
            ),
        ]
        gwas = [
            GwasResult("gA", "gene", 5.0, 0.001, 0.01, 5.0, 0.3),
            GwasResult("gB", "gene", 5.0, 0.001, 0.01, 5.0, 0.3),
        ]
        supported = conv.supported_associations(calls, gwas, config)
        assert [c.determinant for c in supported] == ["gA"]

    def test_significant_but_not_convergent_is_excluded(self, config):
        calls = [
            conv.ConvergenceCall(
                "gB", 1, 6, {"L1": 6}, ConvergenceCategory.LINEAGE_ASSOCIATED
            )
        ]
        gwas = [GwasResult("gB", "gene", 5.0, 0.001, 0.01, 5.0, 0.3)]
        assert conv.supported_associations(calls, gwas, config) == []

    def test_convergent_but_not_significant_is_excluded(self, config):
        calls = [
            conv.ConvergenceCall(
                "gA", 3, 6, {"L1": 3, "L2": 3}, ConvergenceCategory.RECURRENT_CONVERGENT
            )
        ]
        gwas = [GwasResult("gA", "gene", 5.0, 0.4, 0.8, 5.0, 0.3)]
        assert conv.supported_associations(calls, gwas, config) == []


# --------------------------------------------------------------------------
# co-occurrence
# --------------------------------------------------------------------------


class TestFeatureSets:
    def test_builds_three_namespaces(self, sample_manifest):
        namespaces = cooc.build_feature_sets(
            sample_manifest, {"TEST_A_01": ["gA"]}, {}, {}
        )
        assert set(namespaces) == {"gene", "mutation", "mechanism"}

    def test_sample_absent_from_source_gets_empty_set(self, sample_manifest):
        namespaces = cooc.build_feature_sets(sample_manifest, {}, {}, {})
        assert namespaces["gene"]["TEST_A_01"] == set()


class TestAllPairs:
    def test_pairs_are_unordered_and_distinct(self):
        pairs = cooc.all_pairs(["b", "a", "c"])
        assert pairs == [("a", "b"), ("a", "c"), ("b", "c")]

    def test_single_name_yields_no_pairs(self):
        assert cooc.all_pairs(["a"]) == []

    def test_duplicates_collapse(self):
        assert len(cooc.all_pairs(["a", "a", "b"])) == 1


class TestTestPair:
    def _sets(self):
        set_a = {"s1": {"gA"}, "s2": {"gA"}, "s3": set(), "s4": set()}
        set_b = {"s1": {"gB"}, "s2": set(), "s3": {"gB"}, "s4": set()}
        return set_a, set_b

    def test_counts_are_correct(self):
        set_a, set_b = self._sets()
        result = cooc.test_pair(
            "gA", "gB", set_a, set_b, "gene_gene", ["s1", "s2", "s3", "s4"]
        )
        assert result.n_a == 2
        assert result.n_b == 2
        assert result.n_both == 1

    def test_interpretation_limit_is_recorded(self):
        set_a, set_b = self._sets()
        result = cooc.test_pair(
            "gA", "gB", set_a, set_b, "gene_gene", ["s1", "s2", "s3", "s4"]
        )
        assert result.interpretation_limit == "association_only_not_causal"

    def test_degenerate_pair_returns_none(self):
        """A feature present in every sample cannot be tested for association."""
        set_a = {f"s{i}": {"gA"} for i in range(4)}
        set_b = {"s1": {"gB"}, "s2": set(), "s3": set(), "s4": set()}
        assert (
            cooc.test_pair("gA", "gB", set_a, set_b, "gene_gene", list(set_a)) is None
        )

    def test_p_value_is_produced(self):
        set_a, set_b = self._sets()
        result = cooc.test_pair(
            "gA", "gB", set_a, set_b, "gene_gene", ["s1", "s2", "s3", "s4"]
        )
        assert 0.0 <= result.adjusted_p_value <= 1.0


class TestBenjaminiHochberg:
    def test_available_locally(self):
        adjusted = cooc.benjamini_hochberg({"a": 0.01, "b": 0.02, "c": 0.9})
        assert adjusted["a"] <= adjusted["b"] <= adjusted["c"]

    def test_empty(self):
        assert cooc.benjamini_hochberg({}) == {}


# --------------------------------------------------------------------------
# integration
# --------------------------------------------------------------------------


def _master(sample_id, **kwargs):
    defaults = dict(sample_id=sample_id, antibiotic="imipenem")
    defaults.update(kwargs)
    return MasterRecord(**defaults)


class TestJoin:
    def test_none_means_empty_not_empty_string(self):
        assert integ._join([]) is None
        assert integ._join(["a", "b"]) == "a,b"

    def test_values_are_sorted_and_deduped(self):
        assert integ._join(["b", "a", "b"]) == "a,b"

    def test_none_entries_dropped(self):
        assert integ._join(["a", None, "b"]) == "a,b"


class TestSampleConfidence:
    def test_supported_requires_significance_and_convergence(self, config):
        from papipeline.models import MechanismCall

        mechanisms = [
            MechanismCall("S1", "imipenem", "gA", "efflux", ClaimStatus.DETECTED)
        ]
        gwas = [GwasResult("gA", "gene", 5.0, 0.001, 0.01, 5.0, 0.3)]
        conv_calls = [
            conv.ConvergenceCall(
                "gA", 2, 6, {"L1": 3, "L2": 3}, ConvergenceCategory.RECURRENT_CONVERGENT
            )
        ]
        assert (
            integ.sample_confidence(config, "S1", mechanisms, gwas, conv_calls, [])
            is ClaimStatus.SUPPORTED
        )

    def test_associated_when_only_significant(self, config):
        from papipeline.models import MechanismCall

        mechanisms = [
            MechanismCall("S1", "imipenem", "gA", "efflux", ClaimStatus.DETECTED)
        ]
        gwas = [GwasResult("gA", "gene", 5.0, 0.001, 0.01, 5.0, 0.3)]
        assert (
            integ.sample_confidence(config, "S1", mechanisms, gwas, [], [])
            is ClaimStatus.ASSOCIATED
        )

    def test_detected_without_statistical_support(self, config):
        from papipeline.models import MechanismCall

        mechanisms = [
            MechanismCall("S1", "imipenem", "gA", "efflux", ClaimStatus.DETECTED)
        ]
        assert (
            integ.sample_confidence(config, "S1", mechanisms, [], [], [])
            is ClaimStatus.DETECTED
        )

    def test_candidate_sv_gives_predicted_not_detected(self, config):
        sv = StructuralVariant(
            sample_id="S1",
            variant_id="SV1",
            variant_type="deletion",
            position=1,
            affected_gene="oprD",
            size=100,
            evidence="weak",
            confidence="low",
            call_status=StructuralCallStatus.CANDIDATE,
        )
        assert (
            integ.sample_confidence(config, "S1", [], [], [], [sv])
            is ClaimStatus.PREDICTED
        )

    def test_no_evidence_is_unknown(self, config):
        assert (
            integ.sample_confidence(config, "S1", [], [], [], []) is ClaimStatus.UNKNOWN
        )

    def test_confidence_is_never_causal(self, config):
        """There is no CAUSAL member in the vocabulary at all."""
        assert "CAUSAL" not in {s.name for s in ClaimStatus}


class TestBuildMasterTable:
    def _run(self, config, manifest):
        phenotype = {
            manifest.sample_ids[0]: PhenotypeCall(
                manifest.sample_ids[0], "imipenem", Phenotype.R
            )
        }
        return integ.build_master_table(
            config,
            manifest,
            "imipenem",
            phenotype,
            {},
            {},
            {},
            {},
            {},
            {},
            {},
            [],
            [],
        )

    def test_one_record_per_manifest_sample(self, config, manifest):
        records = self._run(config, manifest)
        assert len(records) == len(manifest)
        assert [r.sample_id for r in records] == manifest.sample_ids

    def test_missing_phenotype_stays_none(self, config, manifest):
        records = self._run(config, manifest)
        without = [r for r in records if r.sample_id != manifest.sample_ids[0]]
        assert all(r.phenotype is None for r in without)

    def test_notes_state_detection_is_not_resistance(self, config, manifest):
        records = self._run(config, manifest)
        assert all(
            "detection_is_not_resistance" in (r.evidence_notes or "") for r in records
        )

    def test_candidate_sv_is_labelled(self, config, manifest):
        sample_id = manifest.sample_ids[0]
        sv = StructuralVariant(
            sample_id=sample_id,
            variant_id="SV_CAND",
            variant_type="deletion",
            position=1,
            affected_gene="oprD",
            size=100,
            evidence="weak",
            confidence="low",
            call_status=StructuralCallStatus.CANDIDATE,
        )
        records = integ.build_master_table(
            config, manifest, "imipenem", {}, {}, {}, {}, {sample_id: [sv]}, {}, {}, {}, [], []
        )
        record = [r for r in records if r.sample_id == sample_id][0]
        assert record.structural_variant == "candidate:SV_CAND"
        assert record.confidence is ClaimStatus.PREDICTED

    def test_confirmed_sv_has_no_candidate_prefix(self, config, manifest):
        sample_id = manifest.sample_ids[0]
        sv = StructuralVariant(
            sample_id=sample_id,
            variant_id="SV_CONF",
            variant_type="insertion_sequence",
            position=1,
            affected_gene="oprD",
            size=100,
            evidence="strong",
            confidence="high",
            call_status=StructuralCallStatus.CONFIRMED,
        )
        records = integ.build_master_table(
            config, manifest, "imipenem", {}, {}, {}, {}, {sample_id: [sv]}, {}, {}, {}, [], []
        )
        record = [r for r in records if r.sample_id == sample_id][0]
        assert record.structural_variant == "SV_CONF"
        assert "candidate" not in (record.structural_variant or "")


class TestMasterTableSchema:
    def test_declared_columns_match_specification(self):
        expected = [
            "sample_id",
            "antibiotic",
            "phenotype",
            "amr_gene",
            "amr_variant",
            "chromosomal_mutation",
            "regulator",
            "mechanism",
            "structural_variant",
            "mlst",
            "lineage",
            "virulence_profile",
            "gwas_feature",
            "gwas_status",
            "convergence_status",
            "confidence",
        ]
        assert list(MASTER_TABLE_COLUMNS)[: len(expected)] == expected

    def test_to_row_projects_all_columns(self):
        row = _master("S1").to_row()
        assert set(row) == set(MASTER_TABLE_COLUMNS)

    def test_write_and_read_round_trip(self, tmp_path):
        from papipeline.io.tsv import read_tsv

        records = [_master("S1", phenotype="R"), _master("S2")]
        path = tmp_path / "master.tsv"
        integ.write_master_table(records, path)
        rows = read_tsv(path, required_columns=list(MASTER_TABLE_COLUMNS))
        assert len(rows) == 2
        assert rows[0]["phenotype"] == "R"
        assert rows[1]["phenotype"] is None


# --------------------------------------------------------------------------
# A2: the completeness check on `lineage_label` (FINISH LINE item 5)
# --------------------------------------------------------------------------
#
# Shown RED first. At cf8476f, `build_master_table` had no `mode` argument and
# no completeness check, so this input produced three records - S2 among them
# with `lineage=None`:
#
#   PRE-FIX: build_master_table returned 3 records with no mode argument
#      S1  lineage='235'
#      S2  lineage=None          <- no refusal, and none was possible
#      S3  lineage='235'
#
# And under the producer as it then stood, EVERY sample read "unknown", so the
# table could not distinguish "we do not know" from "we know and it is
# nothing".


class TestRealRefusesOnAnIncompleteLineageColumn:
    SAMPLES = ("S1", "S2", "S3", "S4")

    def _manifest(self):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        return SampleManifest(samples=[Sample(sample_id=s) for s in self.SAMPLES])

    def _build(self, config, manifest, lineages, mode=None):
        return integ.build_master_table(
            config, manifest, "imipenem",
            {},          # phenotype
            {},          # amr
            {},          # mechanisms
            {},          # regulator_variants
            {},          # structural
            {},          # mlst
            lineages,
            {},          # virulence
            [],          # gwas_results
            [],          # convergence
            mode=mode,
        )

    def test_a_missing_label_refuses_and_names_the_sample(self, config):
        from papipeline.errors import StageError
        from papipeline.models import RunMode

        manifest = self._manifest()
        with pytest.raises(StageError) as excinfo:
            self._build(config, manifest, {"S1": "235", "S3": "164"}, RunMode.REAL)
        message = str(excinfo.value)
        assert "S2" in message and "S4" in message, (
            "the samples are the actionable part; a count is not"
        )
        assert excinfo.value.context["missing"] == "S2,S4"
        assert excinfo.value.context["n_missing"] == 2

    def test_the_unknown_sentinel_refuses_too(self, config):
        """It is the value the producer actually wrote for an untyped sample, so
        refusing only on absence would have missed every one of them."""
        from papipeline.errors import StageError
        from papipeline.models import RunMode

        manifest = self._manifest()
        with pytest.raises(StageError) as excinfo:
            self._build(
                config, manifest,
                {"S1": "235", "S2": "unknown", "S3": "164", "S4": "235"},
                RunMode.REAL,
            )
        assert "S2" in str(excinfo.value)
        assert excinfo.value.context["sentinel"] == "unknown"

    def test_a_blank_label_refuses(self, config):
        """`load_tree_metadata` maps a missing cell to the sentinel, and a
        whitespace-only cell would otherwise pass the truthiness test."""
        from papipeline.errors import StageError
        from papipeline.models import RunMode

        manifest = self._manifest()
        with pytest.raises(StageError):
            self._build(
                config, manifest,
                {"S1": "235", "S2": "   ", "S3": "164", "S4": "235"},
                RunMode.REAL,
            )

    def test_the_message_says_why_the_column_is_load_bearing(self, config):
        """`lineage` is one column of seventeen. A reader who has to decide
        whether to care needs to be told what breaks without it."""
        from papipeline.errors import StageError
        from papipeline.models import RunMode

        with pytest.raises(StageError) as excinfo:
            self._build(config, self._manifest(), {}, RunMode.REAL)
        message = str(excinfo.value)
        assert "convergence" in message and "co-occurrence" in message
        assert "mlst_results.tsv" in message, (
            "the message must point at the table a reader can act on"
        )

    def test_a_complete_crosswalk_passes(self, config):
        from papipeline.models import RunMode

        records = self._build(
            config, self._manifest(),
            {"S1": "235", "S2": "235", "S3": "164", "S4": "308"},
            RunMode.REAL,
        )
        assert [r.lineage for r in records] == ["235", "235", "164", "308"]

    def test_one_distinct_lineage_is_enough_to_build_but_not_to_conclude(self, config):
        """The check is completeness, not statistical adequacy. A cohort with a
        single lineage is a real limitation the report states; a cohort with a
        missing label is a broken input, and refusing it here is correct.
        Pinning the boundary so the check cannot quietly grow into a
        significance test."""
        from papipeline.models import RunMode

        records = self._build(
            config, self._manifest(), {s: "235" for s in self.SAMPLES}, RunMode.REAL
        )
        assert len(records) == len(self.SAMPLES)

    @pytest.mark.parametrize("mode_name", ["TEST", "STUB"])
    def test_the_other_modes_are_not_gated(self, config, mode_name):
        """The fixture cohort is synthetic and its crosswalk is whatever the test
        passed; gating TEST would break every existing caller for no scientific
        reason. The gate is on REAL, where the master table is a result."""
        from papipeline.models import RunMode

        records = self._build(
            config, self._manifest(), {}, getattr(RunMode, mode_name)
        )
        assert len(records) == len(self.SAMPLES)

    def test_no_mode_means_no_gate(self, config):
        """A caller that has not said which mode it is in has not claimed to be a
        production run. This is why `mode` is optional, and it is asserted so
        the default cannot quietly become REAL."""
        assert len(self._build(config, self._manifest(), {})) == len(self.SAMPLES)

    def test_run_passes_the_mode_through(self):
        """`run` is what `run.py` calls; a `build_master_table` gate that `run`
        forgot to forward would be a gate that never fires."""
        import inspect

        source = inspect.getsource(integ.run)
        assert "mode=mode," in source, (
            "integration.run must forward its mode to build_master_table, or "
            "the REAL completeness check is dead code"
        )
