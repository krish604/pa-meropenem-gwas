"""The two configuration tables this package owns.

``config/known_determinants.tsv`` is DERIVED from ``config/mechanisms.tsv``.
The first test is the one that keeps the derivation honest: regenerate the rows
in memory and compare them with the committed file, so the two cannot drift.
``config/interaction_pairs.tsv`` carries the tier-1 pre-specified pairs, and
every endpoint must resolve to something the knowledge tables (or, for the two
layer composites this build specifies, the layer-encoder spec) declares - an
endpoint that resolves to nothing is refused rather than accepted on trust.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.downstream import interactions as ix
from papipeline.downstream import known_determinants as kd
from papipeline.errors import DataContractError

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"


@pytest.fixture(scope="module")
def config() -> PipelineConfig:
    return load_config(SCIENCE, machine="laptop")


class TestKnownDeterminants:
    def test_the_committed_file_matches_the_derivation(self, config):
        """The file on disk IS the derivation from mechanisms.tsv.

        If config/mechanisms.tsv changes and the derived table is not
        regenerated, this fails - which is the point: a stale novelty filter
        silently calls a newly-known determinant novel.
        """
        committed = kd.read_known_determinants(kd.known_determinants_path(config))
        assert committed == kd.derive_known_determinants(config)

    def test_every_gene_row_traces_to_the_mechanism_table(self, config):
        rows = kd.derive_known_determinants(config)
        gene_rows = [r for r in rows if r["kind"] == kd.KIND_GENE]
        assert {r["determinant"] for r in gene_rows} == set(config.mechanisms)
        for row in gene_rows:
            spec = config.mechanisms[row["determinant"]]
            assert row["mechanism"] == spec.mechanism
            assert row["mechanism_class"] == spec.mechanism_class
            assert row["claim_ceiling"] == spec.claim_ceiling.value
            assert row["source"].endswith("mechanisms.tsv")

    def test_the_mbl_composite_is_derived_from_roles_not_invented(self, config):
        rows = {r["determinant"]: r for r in kd.derive_known_determinants(config)}
        assert kd.MBL_COMPOSITE in rows
        composite = rows[kd.MBL_COMPOSITE]
        assert composite["kind"] == kd.KIND_COMPOSITE
        members = [m for m in composite["members"].split(",") if m]
        assert len(members) >= 2, "an MBL composite over one gene is not a composite"
        for gene in members:
            role = config.mechanisms[gene].biological_role.lower()
            assert "metallo-beta-lactamase" in role, gene
        non_mbl = {
            g
            for g, spec in config.mechanisms.items()
            if "metallo-beta-lactamase" not in spec.biological_role.lower()
        }
        assert not set(members) & non_mbl

    def test_the_known_set_recognises_the_shapes_a_scan_emits(self, config):
        known = kd.load_known_determinants(config)
        # Composite, with and without a layer prefix.
        assert known.is_known("any_MBL")
        assert known.is_known("L1_any_MBL")
        # A gene feature and an oprD layer feature.
        assert known.is_known("gene__oprD_absent")
        assert known.is_known("L3_oprd_off_any")
        # bla-stripped family name, as the layer encoder spells KPC.
        assert known.is_known("L1_KPC")
        # A gene not in the knowledge table at all.
        assert not known.is_known("gene__femA_novel")

    def test_members_are_exposed_for_the_control_gate(self, config):
        known = kd.load_known_determinants(config)
        members = known.members_for(kd.MBL_COMPOSITE)
        assert members and all(m in config.mechanisms for m in members)
        # A control named after a gene resolves to that gene, so
        # `oprD_burden` gets `oprD`'s constituents.
        assert known.members_for("oprD_burden") is None or all(
            m in config.mechanisms for m in known.members_for("oprD_burden")
        )
        assert known.is_known("oprD_burden")


class TestInteractionPairs:
    def test_the_committed_tier_one_pairs_are_all_present(self, config):
        pairs = ix.load_interaction_pairs(
            ix.interaction_pairs_path(config), known=kd.load_known_determinants(config)
        )
        assert pairs, "the tier-1 table is empty"
        assert all(p.tier == 1 for p in pairs)
        endpoints = {(p.feature_a, p.feature_b) for p in pairs}
        for expected in ix.REQUIRED_CARBAPENEM_PAIRS:
            assert expected in endpoints, f"missing pre-specified pair {expected}"

    def test_the_pairs_are_relevant_to_carbapenems_only(self, config):
        pairs = ix.load_interaction_pairs(ix.interaction_pairs_path(config))
        for pair in pairs:
            assert pair.relevant_to("imipenem")
            assert not pair.relevant_to("ciprofloxacin")

    def test_every_declared_endpoint_resolves(self, config):
        known = kd.load_known_determinants(config)
        pairs = ix.load_interaction_pairs(ix.interaction_pairs_path(config), known=known)
        for pair in pairs:
            for endpoint in (pair.feature_a, pair.feature_b):
                assert ix.endpoint_is_declared(
                    endpoint, known
                ), f"{endpoint} resolves to no known determinant and no layer feature"

    def test_tier_three_is_refused_rather_than_quietly_run(self, config, tmp_path):
        path = tmp_path / "pairs.tsv"
        path.write_text(
            "tier\tantibiotic\tfeature_a\tfeature_b\tnotes\n"
            "3\timipenem\tnalC\tmexR\tdeliberately skipped\n"
        )
        with pytest.raises(DataContractError) as exc:
            ix.load_interaction_pairs(path)
        assert "tier 3" in str(exc.value)

    def test_an_undeclared_endpoint_is_refused_by_name(self, config, tmp_path):
        path = tmp_path / "pairs.tsv"
        path.write_text(
            "tier\tantibiotic\tfeature_a\tfeature_b\tnotes\n"
            "1\timipenem\tnotARealLocusX\tnalC\tinvented\n"
        )
        known = kd.load_known_determinants(config)
        with pytest.raises(DataContractError) as exc:
            ix.load_interaction_pairs(path, known=known)
        assert "notARealLocusX" in str(exc.value)
