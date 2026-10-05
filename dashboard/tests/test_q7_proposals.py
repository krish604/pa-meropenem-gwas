"""Q7 — the Phase-1 proposals INTEGRATE accepted, each with a test.

These are the backend-vs-DESIGN mismatches DATA, VIZ and SHELL measured and
proposed. Every one is exercised against the committed 10-isolate fixture or a
copy of it in ``tmp_path``; none touches the network or a real results tree.
"""

from __future__ import annotations

import json

import pytest

from conftest import REPO_ROOT, SMALL_LIVE, copy_layout


@pytest.fixture
def client(app_factory):
    return app_factory(results_root=SMALL_LIVE)[0]


# ---------------------------------------------------------------------------
# DATA — /api/isolates projection and the three promised filters
# ---------------------------------------------------------------------------


def test_isolates_columns_projection_keeps_the_named_columns(client):
    """`?columns=sample_id,st` must return those cells, not `items: [{}]`.

    `app.py` built `keep = set(columns)` from the raw csv string, which is a
    set of *characters*, so every projection kept nothing.
    """
    payload = client.get(
        "/api/isolates", params={"columns": "sample_id,st", "limit": 5}
    ).json()
    assert payload["items"]
    for row in payload["items"]:
        assert set(row) == {"sample_id", "st", "basis"}, row
        assert row["sample_id"].startswith("TEST_PA_")


def test_isolates_columns_projection_rejects_an_unknown_column(client):
    response = client.get("/api/isolates", params={"columns": "sample_id,nope"})
    assert response.status_code == 400
    assert "nope" in response.json()["error"]


def test_isolates_has_amr_gene_filter_matches(client):
    payload = client.get(
        "/api/isolates", params={"filter.has_amr_gene": "1", "limit": 1000}
    ).json()
    assert payload["meta"]["total"] > 0
    for row in payload["items"]:
        assert row["has_amr_gene"] is True
    none = client.get(
        "/api/isolates", params={"filter.has_amr_gene": "0", "limit": 1000}
    ).json()
    assert none["meta"]["total"] == 0


def test_isolates_has_point_mutation_filter_matches(client):
    payload = client.get(
        "/api/isolates", params={"filter.has_point_mutation": "1", "limit": 1000}
    ).json()
    for row in payload["items"]:
        assert row["has_point_mutation"] is True


def test_isolates_virulence_min_filter_is_a_minimum(client):
    """`filter.virulence_min=N` means "at least N", not "equals N"."""
    all_rows = client.get("/api/isolates", params={"limit": 1000}).json()["items"]
    high = max(row["n_virulence"] for row in all_rows)
    payload = client.get(
        "/api/isolates", params={"filter.virulence_min": str(high), "limit": 1000}
    ).json()
    assert payload["meta"]["total"] > 0
    for row in payload["items"]:
        assert row["n_virulence"] >= high


# ---------------------------------------------------------------------------
# DATA — stage-row filters and the unknown table key
# ---------------------------------------------------------------------------


def test_stage_row_filter_on_a_non_key_column_is_a_400(client):
    """A filter the index cannot answer is refused, not silently dropped.

    `_rows_via_index` set `predicate = None` when a clause named a column other
    than the key, so the filter was echoed in `meta.filters` while `total` was
    the unfiltered count — a filter that quietly did not apply looks like a
    filter.
    """
    response = client.get(
        "/api/stages/amr/rows", params={"filter.variant": "AmpC:G183A"}
    )
    assert response.status_code == 400
    assert "variant" in response.json()["error"]


def test_stage_row_filter_on_the_key_column_still_works(client):
    response = client.get(
        "/api/stages/amr/rows", params={"filter.sample_id": "TEST_PA_001"}
    )
    assert response.status_code == 200


def test_unknown_table_key_is_a_404(client):
    response = client.get("/api/tables/not_a_real_table")
    assert response.status_code == 404
    assert "not_a_real_table" in response.json()["error"]


# ---------------------------------------------------------------------------
# DATA / VIZ — basis.n is the cohort, not the row count
# ---------------------------------------------------------------------------


def _basis(payload):
    return (payload.get("meta") or {}).get("basis") or payload.get("basis")


def test_pangenome_basis_n_is_the_cohort_not_the_metric_count(client):
    payload = client.get("/api/pangenome").json()
    basis = payload["basis"]
    assert basis["n"] == 10, "the pangenome was computed on the cohort, not on 5 metrics"
    assert basis["rows_total"] == 5, "the metric-row count is still reported"
    assert payload["summary"]["n_samples"] == "10"


@pytest.mark.parametrize("endpoint", ["/api/gwas", "/api/convergence", "/api/cooccurrence"])
def test_statistic_basis_n_is_the_cohort(client, endpoint):
    payload = client.get(endpoint).json()
    basis = _basis(payload)
    assert basis["n"] == 10, f"{endpoint} basis.n must be the cohort (10)"
    assert basis["rows_total"] > 0


def test_gwas_top_basis_n_is_the_cohort_and_reports_n_returned(client):
    payload = client.get("/api/gwas/top", params={"limit": 5}).json()
    assert payload["basis"]["n"] == 10
    assert payload["n_returned"] == len(payload["items"])
    assert payload["n_tests"] >= payload["n_returned"]


def test_power_flag_is_not_underpowered_when_the_cohort_is_adequate(client):
    """The pipeline's `power_flag` says "underpowered" for every n; the flag
    text here is conditional on `basis.n >= min_samples` (3)."""
    payload = client.get("/api/gwas").json()
    assert payload["power"]["underpowered"] is False
    assert "underpowered" not in payload["power"]["flag"]
    assert payload["power"]["flag"] == "n=10, not a finding"
    # The R16 meaning is always carried verbatim.
    assert "Ruling R16" in payload["power"]["meaning"]


def test_power_flag_is_underpowered_below_the_minimum(client):
    """A statistic computed on fewer than `min_samples` isolates still flags."""
    from dashboard.server.app import _power

    flag = _power(2, min_samples=3)
    assert flag["underpowered"] is True
    assert flag["flag"] == "n=2, underpowered, not a finding"


# ---------------------------------------------------------------------------
# VIZ — the tree tip phenotype is populated from 11_phenotype.tsv
# ---------------------------------------------------------------------------


def test_tree_tip_phenotype_is_joined_from_the_phenotype_table(client):
    payload = client.get("/api/tree").json()
    tips = {node["label"]: node for node in payload["nodes"] if node["is_tip"]}
    assert tips["TEST_PA_001"]["tip_metadata"]["phenotype"] == "R"
    assert tips["TEST_PA_002"]["tip_metadata"]["phenotype"] == "S"
    # The tips endpoint exposes the same join.
    first = client.get("/api/tree/tips", params={"limit": 1}).json()["items"][0]
    assert first["phenotype"] in ("R", "I", "S", "SDD", "ND")


# ---------------------------------------------------------------------------
# DATA — the oprD evidence carries every additional column (§4.3)
# ---------------------------------------------------------------------------


def test_oprd_evidence_carries_repeat_compensating_indel_and_tblastn(app_factory, tmp_path):
    root = copy_layout(SMALL_LIVE, tmp_path / "oprd")
    stages = root / "intermediate" / "stages"
    (stages / "oprd.tsv").write_text(
        "# SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS\n"
        "sample_id\tverdict\tlesion_type\tposition\ttruncation_aa\tidentity_pct\t"
        "coverage_pct\trepeat_flag\tcompensating_indel\ttblastn_summary\n"
        "TEST_PA_001\tresolved\ttruncation\t124500\t41\t99.5\t100\t"
        "no_repeat\tyes\thit_at_124500\n",
        encoding="utf-8",
    )
    client = app_factory(results_root=root)[0]
    payload = client.get("/api/oprd").json()
    assert payload["available"] is True
    item = payload["items"][0]
    evidence = item["evidence"]
    assert evidence["repeat_flag"] == "no_repeat"
    assert evidence["compensating_indel"] == "yes"
    assert evidence["tblastn_summary"] == "hit_at_124500"
