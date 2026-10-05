"""Tests for figure data preparation (papipeline.viz).

These assert that the figure layer is genuinely data-driven: no hard-coded
genes, categories or results, and no fabricated numbers when a field is
missing.
"""

from __future__ import annotations

import pytest

from papipeline.models import (
    AmrDeterminant,
    GwasResult,
    MechanismCall,
    Phenotype,
    PhenotypeCall,
    RegulatorVariant,
    StructuralCallStatus,
    StructuralVariant,
    VirulenceFactor,
    VariantType,
)
from papipeline import viz


def _phenotype(**mapping):
    return {
        sample_id: PhenotypeCall(sample_id, "imipenem", Phenotype(value))
        for sample_id, value in mapping.items()
    }


def _amr(**mapping):
    return {
        sample_id: [
            AmrDeterminant(
                sample_id=sample_id,
                antibiotic="imipenem",
                determinant=gene,
                gene=gene,
                variant=None,
                determinant_type="acquired",
                mechanism=None,
                evidence_source="test",
                database="DB",
                database_version="v1",
                confidence="99",
            )
        ]
        for sample_id, gene in mapping.items()
    }


def _mechanisms(**mapping):
    return {
        sample_id: [
            MechanismCall(sample_id, "imipenem", determinant, mechanism, "DETECTED")
            for determinant, mechanism in pairs
        ]
        for sample_id, pairs in mapping.items()
    }


def _variant(sample_id, gene, variant_type=VariantType.SNV):
    return RegulatorVariant(
        sample_id=sample_id,
        gene=gene,
        variant=f"{gene}_v1",
        variant_type=variant_type.value,
        position=100,
        reference="A",
        alternate="T",
        effect="test",
        mechanism="efflux",
        confidence="high",
    )


class TestPhenotypeDistribution:
    def test_counts_and_order(self):
        categories, counts = viz.phenotype_distribution_table(
            _phenotype(a="R", b="S", c="S", d="I")
        )
        assert counts["S"] == 2
        assert categories == ["R", "I", "S"]

    def test_clinical_order_is_stable(self):
        categories, _ = viz.phenotype_distribution_table(
            _phenotype(a="S", b="R", c="ND", d="SDD")
        )
        assert categories == ["R", "SDD", "S", "ND"]

    def test_empty_input(self):
        categories, counts = viz.phenotype_distribution_table({})
        assert categories == []
        assert dict(counts) == {}


class TestDeterminantDistribution:
    def test_counts_samples_not_calls(self):
        """A gene carried twice by one sample counts once."""
        sample = "S1"
        determinants = [
            AmrDeterminant(
                sample_id=sample,
                antibiotic="imipenem",
                determinant="gA",
                gene="gA",
                variant=None,
                determinant_type="acquired",
                mechanism=None,
                evidence_source="t",
                database="D",
                database_version="v",
                confidence="1",
            )
            for _ in range(3)
        ]
        _categories, counts = viz.determinant_distribution_table({sample: determinants})
        assert counts["gA"] == 1

    def test_multiple_genes(self):
        _categories, counts = viz.determinant_distribution_table(
            _amr(s1="gA", s2="gA", s3="gB")
        )
        assert counts["gA"] == 2
        assert counts["gB"] == 1


class TestMechanismDistribution:
    def test_counts_mechanisms(self):
        _names, counts = viz.mechanism_distribution_table(
            _mechanisms(s1=[("d1", "efflux")], s2=[("d2", "efflux"), ("d3", "acquired_determinant")])
        )
        assert counts["efflux"] == 2
        assert counts["acquired_determinant"] == 1


class TestCrosstab:
    def test_carriers_and_totals(self):
        feature_sets = {"s1": {"gA"}, "s2": {"gA"}, "s3": set()}
        phenotype = {"s1": "R", "s2": "S", "s3": "S"}
        order, features, counts = viz.crosstab(feature_sets, "gA", phenotype)
        assert counts[("R", "gA")] == 1
        assert counts[("S", "gA")] == 1
        assert features == ["gA"]

    def test_missing_phenotype_excluded_not_binned(self):
        """A sample with no phenotype is not folded into a neighbour."""
        feature_sets = {"s1": {"gA"}, "s2": {"gA"}}
        phenotype = {"s1": "R", "s2": None}
        _order, _features, counts = viz.crosstab(feature_sets, "gA", phenotype)
        assert ("S", "gA") not in counts
        assert counts[("R", "gA")] == 1

    def test_gene_by_phenotype_shape(self):
        table = viz.gene_by_phenotype_table(
            _amr(s1="gA", s2="gA", s3="gB"), _phenotype(s1="R", s2="S", s3="R")
        )
        assert set(table) == {"gA", "gB"}
        assert table["gA"]["R"]["carriers"] == 1
        assert table["gA"]["R"]["total"] == 2

    def test_category_axis_spans_all_observed_phenotypes(self):
        """A category where the feature is absent must still appear.

        Regression test: deriving the axis from the observed cells made a
        phenotype category disappear from the figure entirely.
        """
        table = viz.gene_by_phenotype_table(
            _amr(s1="gA", s2="gA"), _phenotype(s1="R", s2="S")
        )
        assert set(table["gA"]) == {"R", "S"}
        assert table["gA"]["R"]["carriers"] == 1
        assert table["gA"]["S"]["carriers"] == 1

    def test_totals_are_cohort_wide_not_feature_wide(self):
        table = viz.gene_by_phenotype_table(
            _amr(s1="gA", s2="gA", s3="gB"), _phenotype(s1="R", s2="S", s3="R")
        )
        for gene in ("gA", "gB"):
            assert table[gene]["R"]["total"] == 2
            assert table[gene]["S"]["total"] == 1

    def test_crosstab_indexes_by_sample_not_by_feature(self):
        """Regression test: the gene/phenotype axes were transposed.

        Passing a gene->samples map where sample->genes was expected made
        every cell zero.
        """
        table = viz.gene_by_phenotype_table(
            _amr(s1="gA", s2="gB"), _phenotype(s1="R", s2="R")
        )
        assert table["gA"]["R"]["carriers"] == 1
        assert table["gB"]["R"]["carriers"] == 1

    def test_invert_helper(self):
        inverted = viz._invert({"gA": ["s1", "s2"], "gB": ["s2"]})
        assert inverted == {"s1": {"gA"}, "s2": {"gA", "gB"}}

    def test_mechanism_by_phenotype_shape(self):
        table = viz.mechanism_by_phenotype_table(
            _mechanisms(s1=[("d1", "efflux")], s2=[("d2", "efflux")]),
            _phenotype(s1="R", s2="S"),
        )
        assert table["efflux"]["R"]["carriers"] == 1
        assert table["efflux"]["R"]["total"] == 1


class TestOprdStatus:
    def test_absence_from_variant_screen(self):
        status = viz.oprd_status_per_sample(
            {"s1": [_variant("s1", "oprD", VariantType.GENE_ABSENCE)]}
        )
        assert status["s1"] == "absent"

    def test_disruption(self):
        status = viz.oprd_status_per_sample(
            {"s1": [_variant("s1", "oprD", VariantType.GENE_DISRUPTION)]}
        )
        assert status["s1"] == "disrupted"

    def test_intact_requires_annotation_evidence(self):
        """Absence of a variant call is not evidence of an intact gene."""
        status = viz.oprd_status_per_sample({"s1": []}, None)
        assert status == {}

        status = viz.oprd_status_per_sample({"s1": []}, {"s1": ["oprD", "mexR"]})
        assert status["s1"] == "intact"

    def test_absent_only_from_an_absent_verdict(self):
        """A missing annotation symbol is NOT evidence the gene is gone.

        This previously asserted ``absent``. That was the fabricated negative:
        ``"oprD" in genes`` reads Bakta's ``gene=`` symbol, which Bakta also
        writes on the OprP/OprQ paralog family and which it omits on isolates
        that do carry the true locus as an unlabelled CDS. Eight of the ten
        smoke isolates were called ``absent`` this way while all ten carry the
        gene. The only legitimate sources of ``absent`` are the variant screen's
        GENE_ABSENCE and a structural ``absent`` verdict.
        """
        status = viz.oprd_status_per_sample({"s1": []}, {"s1": ["mexR"]})
        assert status["s1"] == "not_assessed"

    def test_crosstab_lists_all_three_states(self):
        variants = {
            "s1": [_variant("s1", "oprD", VariantType.GENE_ABSENCE)],
            "s2": [_variant("s2", "oprD", VariantType.GENE_DISRUPTION)],
        }
        annotated = {"s3": ["oprD"], "s1": [], "s2": []}
        status = viz.oprd_status_per_sample(variants, annotated)
        phenotypes, states, counts = viz.oprd_by_phenotype_table(
            status, _phenotype(s1="R", s2="S", s3="S")
        )
        assert set(states) == {"intact", "disrupted", "absent"}
        assert counts[("R", "absent")] == 1
        assert counts[("S", "disrupted")] == 1

    def test_unknown_state_appears_but_is_not_dropped(self):
        status = viz.oprd_status_per_sample({"s1": []}, None)
        _phenotypes, states, _counts = viz.oprd_by_phenotype_table(
            status, _phenotype(s1="R")
        )
        assert states == ["not_assessed"]


class TestOprdStatusFromLocusResolution:
    """The seam where the paralog problem is fixed.

    `"oprD" in genes` asks whether Bakta wrote the symbol, and on real data Bakta
    writes it on the OprD/OprP/OprQ family: 34 of the 36 `gene=oprD` CDS across
    the ten smoke isolates are paralogs at 33.7-36.8% identity to PAO1 PA0958,
    and 8 of the 10 isolates carry the true locus unlabelled. So the symbol test
    reports `absent` where the gene is present. These tests pin the replacement.
    """

    def _resolution(self, sample_id, verdict):
        return type(
            "R", (), {"is_resolved": verdict == "resolved", "verdict": verdict}
        )()

    def test_a_resolved_locus_is_intact(self):
        status = viz.oprd_status_per_sample(
            {"s1": []}, None, {"s1": self._resolution("s1", "resolved")}
        )
        assert status["s1"] == "intact"

    def test_a_refusal_is_not_assessed_and_never_absent(self):
        """The property the whole change exists to protect.

        Not locating a locus is not evidence of deleting it. Reporting `absent`
        here would manufacture the study's central negative from a failed
        alignment.
        """
        for verdict in (
            "refused:no_orthologous_hit",
            "refused:insufficient_coverage",
            "refused:ambiguous_locus",
            "refused:ambiguous_second_locus",
        ):
            status = viz.oprd_status_per_sample(
                {"s1": []}, None, {"s1": self._resolution("s1", verdict)}
            )
            assert status["s1"] == "not_assessed", verdict

    def test_a_refusal_overrides_the_symbol_test(self):
        """An isolate whose paralogs made Bakta write `oprD` but whose own locus
        could not be located must not be reported intact on the symbol alone."""
        status = viz.oprd_status_per_sample(
            {"s1": []},
            {"s1": ["oprD", "oprD"]},
            {"s1": self._resolution("s1", "refused:insufficient_coverage")},
        )
        assert status["s1"] == "not_assessed"

    def test_a_paralog_only_isolate_is_not_absent(self):
        """The measured case: paralogs present, symbol `oprD` written, no locus."""
        status = viz.oprd_status_per_sample(
            {"s1": []},
            {"s1": ["oprD"]},
            {"s1": self._resolution("s1", "refused:no_orthologous_hit")},
        )
        assert status["s1"] != "absent"

    def test_a_resolution_does_not_override_the_variant_screen(self):
        status = viz.oprd_status_per_sample(
            {"s1": [_variant("s1", "oprD", VariantType.GENE_ABSENCE)]},
            None,
            {"s1": self._resolution("s1", "resolved")},
        )
        assert status["s1"] == "absent"

    def test_without_a_resolution_the_symbol_test_still_runs(self):
        """Callers that pass nothing still get a status, and it is never ``absent``.

        The second half of this previously asserted ``absent`` for a sample whose
        annotation lacks the ``oprD`` symbol. That expectation encoded the
        fabricated negative and is deliberately reversed: without a locus search
        there is no evidence of absence, so the honest state is
        ``not_assessed``. This IS a visible change to figure 6 for any caller
        that supplies no resolution, and it is the intended correction.
        """
        status = viz.oprd_status_per_sample({"s1": []}, {"s1": ["oprD"]}, None)
        assert status["s1"] == "intact"
        status = viz.oprd_status_per_sample({"s1": []}, {"s1": ["mexR"]}, None)
        assert status["s1"] == "not_assessed"

    def test_not_assessed_is_a_state_the_figure_already_knows(self):
        """If this were not in OPRD_STATE_ORDER the state would be appended as an
        unknown rather than ordered, which is a silent change to figure 6."""
        from papipeline.viz import OPRD_STATE_ORDER

        assert "not_assessed" in OPRD_STATE_ORDER


class TestEffluxRegulator:
    def test_lists_every_screened_locus(self):
        """Loci with no variants must still appear, with zero counts."""
        variants = {"s1": [_variant("s1", "mexR")]}
        _phenotypes, genes, counts = viz.efflux_regulator_by_phenotype_table(
            variants, ["mexR", "nfxB", "mexZ"], _phenotype(s1="R")
        )
        assert genes == ["mexR", "nfxB", "mexZ"]
        assert counts[("R", "mexR")] == 1
        assert ("R", "nfxB") not in counts

    def test_non_regulator_gene_is_ignored(self):
        variants = {"s1": [_variant("s1", "oprD")]}
        _phenotypes, _genes, counts = viz.efflux_regulator_by_phenotype_table(
            variants, ["mexR"], _phenotype(s1="R")
        )
        assert counts == {}


class TestPresenceMatrix:
    def test_rows_are_features(self):
        feature_sets = {"s1": {"gA"}, "s2": set()}
        features, samples, matrix = viz.presence_matrix(feature_sets, ["s1", "s2"])
        assert features == ["gA"]
        assert samples == ["s1", "s2"]
        assert matrix == [[1, 0]]

    def test_amr_heatmap(self):
        features, samples, matrix = viz.amr_heatmap_matrix(
            _amr(s1="gA", s2="gB"), ["s1", "s2"]
        )
        assert sorted(features) == ["gA", "gB"]
        assert len(matrix) == 2

    def test_mechanism_heatmap(self):
        features, _samples, matrix = viz.mechanism_heatmap_matrix(
            _mechanisms(s1=[("d1", "efflux")]), ["s1"]
        )
        assert features == ["efflux"]
        assert matrix == [[1]]

    def test_empty_cohort(self):
        features, samples, matrix = viz.amr_heatmap_matrix({}, [])
        assert features == []
        assert samples == []
        assert matrix == []

    def test_matrix_is_rectangular(self):
        _f, _s, matrix = viz.amr_heatmap_matrix(
            _amr(s1="gA", s2="gB", s3="gC"), ["s1", "s2", "s3"]
        )
        assert all(len(row) == 3 for row in matrix)


class TestTreeAnnotations:
    def test_one_row_per_sample(self):
        rows = viz.tree_annotation_table(
            ["s1", "s2"],
            _amr(s1="gA"),
            _mechanisms(s1=[("d1", "efflux")]),
            _phenotype(s1="R", s2="S"),
            {"s1": "L1", "s2": "L2"},
            {},
        )
        assert len(rows) == 2
        assert rows[0]["amr_genes"] == ["gA"]
        assert rows[0]["mechanisms"] == ["efflux"]

    def test_missing_phenotype_is_none(self):
        rows = viz.tree_annotation_table(["s1"], {}, {}, {}, {}, {})
        assert rows[0]["phenotype"] is None

    def test_empty_sets_are_lists_not_none(self):
        rows = viz.tree_annotation_table(["s1"], {}, {}, _phenotype(s1="R"), {}, {})
        assert rows[0]["amr_genes"] == []
        assert rows[0]["mechanisms"] == []


class TestGwasAssociationTable:
    def _results(self):
        return [
            GwasResult("gene__a", "gene", 5.0, 0.001, 0.01, 5.0, 0.3, {"L1": 1, "L2": 1, "L3": 1}),
            GwasResult("unitig__b", "unitig", 1.0, 0.2, 0.5, 1.0, 0.3, {"L1": 6}),
        ]

    def test_orders_by_adjusted_p(self, config):
        rows = viz.gwas_association_table(self._results(), config)
        assert rows[0]["feature"] == "gene__a"

    def test_lineage_linked_flag(self, config):
        rows = viz.gwas_association_table(self._results(), config)
        by_name = {r["feature"]: r for r in rows}
        assert by_name["unitig__b"]["lineage_linked"] is True
        assert by_name["gene__a"]["lineage_linked"] is False

    def test_threshold_flag(self, config):
        rows = viz.gwas_association_table(self._results(), config)
        by_name = {r["feature"]: r for r in rows}
        assert by_name["gene__a"]["passes_threshold"] is True
        assert by_name["unitig__b"]["passes_threshold"] is False

    def test_top_n_limit(self, config):
        results = [
            GwasResult(f"f{i}", "snp", 1.0, 0.001 * (i + 1), 0.01, 1.0, 0.5)
            for i in range(50)
        ]
        assert len(viz.gwas_association_table(results, config, top=5)) == 5

    def test_empty_results(self, config):
        assert viz.gwas_association_table([], config) == []


class TestConvergenceTable:
    def test_counts_categories(self):
        from papipeline.models import ConvergenceCall, ConvergenceCategory

        calls = [
            ConvergenceCall("gA", 2, 4, {"L1": 2, "L2": 2}, ConvergenceCategory.RECURRENT_CONVERGENT),
            ConvergenceCall("gB", 1, 4, {"L1": 4}, ConvergenceCategory.LINEAGE_ASSOCIATED),
        ]
        categories, counts = viz.convergence_table(calls)
        assert counts["recurrent_convergent"] == 1
        assert counts["lineage_associated"] == 1
        assert set(categories) == {"recurrent_convergent", "lineage_associated"}


class TestCooccurrenceNetwork:
    def _pair(self, a, b, p, n_both=2):
        from papipeline.models import Cooccurrence

        return Cooccurrence(
            feature_a=a,
            feature_b=b,
            feature_type="gene_gene",
            n_a=5,
            n_b=5,
            n_both=n_both,
            statistic="odds_ratio",
            statistic_value=3.0,
            adjusted_p_value=p,
        )

    def test_only_significant_pairs_become_edges(self):
        nodes, edges = viz.cooccurrence_network_table(
            [self._pair("gA", "gB", 0.001), self._pair("gA", "gC", 0.9)]
        )
        assert len(edges) == 1
        assert edges[0]["source"] == "gA"

    def test_nodes_derived_from_edges(self):
        nodes, edges = viz.cooccurrence_network_table([self._pair("gA", "gB", 0.001)])
        assert {n["name"] for n in nodes} == {"gA", "gB"}
        assert all(n["degree"] == 1 for n in nodes)

    def test_interpretation_limit_is_carried(self):
        _nodes, edges = viz.cooccurrence_network_table([self._pair("gA", "gB", 0.001)])
        assert edges[0]["interpretation_limit"] == "association_only_not_causal"

    def test_empty_input(self):
        nodes, edges = viz.cooccurrence_network_table([])
        assert nodes == []
        assert edges == []
