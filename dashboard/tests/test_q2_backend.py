"""Q2 — backend contract tests (pytest + FastAPI TestClient).

Covers both layouts, the full openapi surface, UI-D1 (read-only), UI-D2
(honest states), UI-D4 (path safety and bind/token), UI-D5 (huge table paging),
UI-D6 (duplicate-label tree), oprD (§4.3) and provenance (`not reported`, never
0).
"""

from __future__ import annotations

import json
import os
import re
import tracemalloc
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

from conftest import (
    FIXTURES_DIR,
    REPO_ROOT,
    SMALL_BUNDLE,
    SMALL_LIVE,
    copy_layout,
    edge_layout,
    snapshot_tree,
)

from dashboard.server import newick as newick_module
from dashboard.server import security as security_module
from dashboard.server.security import BindRefused, ContainmentError, resolve_bind, resolve_contained

OPENAPI = REPO_ROOT / "dashboard" / "openapi.yaml"


# ---------------------------------------------------------------------------
# The openapi surface
# ---------------------------------------------------------------------------


def openapi_operations() -> List[Tuple[str, str]]:
    """`(METHOD, path)` for every operation in `dashboard/openapi.yaml`."""
    lines = OPENAPI.read_text(encoding="utf-8").splitlines()
    operations: List[Tuple[str, str]] = []
    current = None
    for line in lines:
        match = re.match(r"^  (/api/\S+):\s*$", line)
        if match:
            current = match.group(1)
            continue
        match = re.match(r"^    (get|post):\s*$", line)
        if match and current:
            operations.append((match.group(1).upper(), current))
    return operations


SUBSTITUTIONS = {
    "/api/stages/{stage}": "/api/stages/mlst",
    "/api/stages/{stage}/rows": "/api/stages/mlst/rows",
    "/api/tables/{key}": "/api/tables/master_table",
    "/api/isolates/{sample_id}": "/api/isolates/TEST_PA_001",
    "/api/reports/content": "/api/reports/content?path=reports/test_pipeline_report.md",
}

REQUIRED_KEYS = {
    "/api/health": {"ok", "source_kind", "root", "manifest_writer", "mode"},
    "/api/run": {"manifest_writer", "stages"},
    "/api/run/counts": {"n_stages_completed", "n_stages_not_run"},
    "/api/stages": {"stage_order", "items"},
    "/api/stages/{stage}": {"name", "state", "reason", "present", "path", "table"},
    "/api/stages/{stage}/rows": {"stage", "present", "reason", "meta"},
    "/api/tables/{key}": {"present", "reason", "meta"},
    "/api/isolates": {"items", "meta", "membership_source"},
    "/api/isolates/{sample_id}": {"sample_id", "tables"},
    "/api/oprd": {"available", "probes"},
    "/api/pangenome": {"present", "reason"},
    "/api/tree": {"present", "reason"},
    "/api/tree/tips": {"items", "meta"},
    "/api/similarity": {"present", "reason"},
    "/api/gwas": {"present", "reason", "n_tests"},
    "/api/gwas/top": {"present", "reason", "items"},
    "/api/convergence": {"present", "reason", "category_order", "items"},
    "/api/cooccurrence": {"present", "reason", "items"},
    "/api/provenance": {"manifest_writer"},
    "/api/provenance/bakta": {"available"},
    "/api/provenance/logs": {"items"},
    "/api/provenance/digests": {"items"},
    "/api/reports": {"items"},
    "/api/reports/content": set(),
    "/api/events/state": {"transport", "cursor", "stages"},
    "/api/launcher/capabilities": {"can_launch", "reason"},
    "/api/launcher/status": {"running", "per_isolate_timings"},
    "/api/launcher/preflight": {"checks", "would_start"},
}


def test_openapi_lists_the_expected_operations():
    operations = openapi_operations()
    assert len(operations) == 29, operations
    paths = {path for _m, path in operations}
    assert paths == set(REQUIRED_KEYS) | {"/api/events"}


def test_every_openapi_endpoint_responds(live_client, monkeypatch):
    client, _state = live_client

    # The SSE generator tails forever by design, so the route is exercised with
    # a finite stub of the same shape; the real generator is asserted separately
    # by `test_sse_frames_starts_with_retry_then_snapshot`.
    async def finite_frames(log_path, *, initial, last_event_id=None, poll_interval=0.1):
        yield b"retry: 2000\n\n"
        yield b"event: snapshot\ndata: {}\n\n"

    monkeypatch.setattr("dashboard.server.events.sse_frames", finite_frames)

    operations = openapi_operations()
    for method, path in operations:
        if path == "/api/events":
            response = client.get(path)
            assert response.status_code == 200, path
            assert response.headers["content-type"].startswith("text/event-stream")
            assert response.text.startswith("retry: 2000")
            continue
        url = SUBSTITUTIONS.get(path, path)
        if method == "GET":
            response = client.get(url)
        else:
            response = client.post(url, json={})
        assert response.status_code == 200, f"{method} {url} -> {response.status_code}"
        if path == "/api/reports/content":
            assert response.headers["content-type"].split(";")[0] in ("text/plain", "text/markdown")
            continue
        payload = response.json()
        for key in REQUIRED_KEYS.get(path, set()):
            assert key in payload, f"{url} missing {key!r}"


def test_sse_frames_starts_with_retry_then_snapshot(live_client):
    import asyncio

    from dashboard.server import events as event_module

    _client, state = live_client

    async def collect():
        generator = event_module.sse_frames(
            state.log_path(), initial={"transport": "sse"}, last_event_id=0
        )
        frames = []
        async for frame in generator:
            frames.append(frame)
            if len(frames) >= 2:
                break
        return frames

    frames = asyncio.run(collect())
    assert frames[0] == b"retry: 2000\n\n"
    assert frames[1].startswith(b"event: snapshot")


# ---------------------------------------------------------------------------
# Layout auto-detection
# ---------------------------------------------------------------------------


def test_live_layout_is_detected_as_live(live_client):
    client, state = live_client
    assert state.kind == "live"
    assert state.manifest_writer == "run"
    assert client.get("/api/health").json()["source_kind"] == "live"


def test_bundle_layout_is_detected_as_bundle(bundle_client):
    client, state = bundle_client
    assert state.kind == "bundle"
    assert state.manifest_writer == "run"
    payload = client.get("/api/health").json()
    assert payload["source_kind"] == "bundle"
    assert client.get("/api/isolates").json()["meta"]["total"] == 10


def test_unmatched_source_names_every_probe_in_order(app_factory, tmp_path, monkeypatch):
    monkeypatch.delenv("PA_DASH_RESULTS_ROOT", raising=False)
    monkeypatch.delenv("PA_DASH_BUNDLE", raising=False)
    # Probes 4-6 resolve a live results root through `config.results_root(mode)`,
    # which honours PIPELINE_RESULTS_ROOT. Point it at a path that does not exist
    # so all three fail with "no such directory" regardless of what this
    # developer's checkout happens to contain. Without this the test only passed
    # on a pristine clone: any `results/test/` left by an earlier pipeline run
    # made probe 5 match, and the dashboard was right to open it.
    monkeypatch.setenv("PIPELINE_RESULTS_ROOT", str(tmp_path / "absent"))
    # No flag, no env: probe 1 is "not given" and probe 3 is "not set", which
    # is the exact §2.3 failure message.
    client, state = app_factory()
    assert state.source is None
    response = client.get("/api/health")
    assert response.status_code == 503
    payload = response.json()
    probes = payload["probes"]
    assert len(probes) == 6
    assert [p["probe"] for p in probes] == [
        "--results-root",
        "PA_DASH_RESULTS_ROOT",
        "bundle (PA_DASH_BUNDLE)",
        "live REAL",
        "live TEST",
        "live STUB",
    ]
    assert [p["n"] for p in probes] == [1, 2, 3, 4, 5, 6]
    # The message names every path, in order, and distinguishes not-configured
    # from configured-and-empty.
    message = payload["message"]
    assert message.index("--results-root") < message.index("PA_DASH_RESULTS_ROOT")
    assert "not given" in message
    assert "not set" in message


# ---------------------------------------------------------------------------
# UI-D1 — the results root is read-only; the state dir is the only write
# ---------------------------------------------------------------------------


def test_results_tree_is_byte_identical_after_a_full_sweep(live_client, tmp_path):
    client, state = live_client
    before = snapshot_tree(SMALL_LIVE)
    for _method, path in openapi_operations():
        if path == "/api/events":
            continue
        url = SUBSTITUTIONS.get(path, path)
        client.get(url) if _method == "GET" else client.post(url, json={})
    after = snapshot_tree(SMALL_LIVE)
    assert before == after, "a dashboard request modified the results root"

    # The state directory is written (it is the only permitted write location).
    state_files = list(Path(state.state_dir.root).rglob("*"))
    assert state_files
    assert not str(state.state_dir.root).startswith(str(SMALL_LIVE))


def test_state_dir_writes_do_not_land_under_the_results_root(live_client):
    client, state = live_client
    client.get("/api/isolates")
    assert state.state_dir.root.exists()
    assert SMALL_LIVE not in state.state_dir.root.parents
    assert state.state_dir.root != SMALL_LIVE


# ---------------------------------------------------------------------------
# UI-D2 — honest states
# ---------------------------------------------------------------------------


def test_six_badges_are_distinct_with_named_reasons(app_factory, tmp_path):
    root = edge_layout(tmp_path)
    client, _state = app_factory(results_root=root)
    items = {item["name"]: item for item in client.get("/api/stages").json()["items"]}

    assert items["validation"]["state"] == "completed"

    failed = items["gwas"]
    assert failed["state"] == "failed"
    assert failed["reason"].strip()

    from papipeline.run import REAL_REFUSING_STAGES

    refused_real = items["convergence"]
    assert refused_real["state"] == "refused"
    assert refused_real["reason"] == REAL_REFUSING_STAGES["convergence"]

    refused_skipped = items["annotation"]
    assert refused_skipped["state"] == "refused"
    assert "analysis.annotation" in refused_skipped["reason"]

    not_assessed = items["mlst"]
    assert not_assessed["state"] == "not_assessed"
    assert not_assessed["reason"].strip()
    assert not_assessed["present"] is False

    not_run = items["pangenome"]
    assert not_run["state"] == "not_run"
    assert not_run["reason"].strip()

    # `not_assessed` and `not_run` are kept apart, and neither is `absent`.
    assert not_assessed["state"] != not_run["state"]
    assert not_assessed["reason"] != not_run["reason"]
    assert "absent" not in (not_assessed["state"], not_run["state"])

    # A header-only table is `not produced`, never an empty result.
    header_only = items["amr"]
    assert header_only["present"] is False
    assert "header and no rows" in header_only["reason"]


def test_missing_values_are_named_not_zero_or_blank(app_factory, tmp_path):
    # A copy with the virulence table removed: its absence must read
    # `not assessed`, not 0.
    root = copy_layout(SMALL_LIVE, tmp_path / "no_virulence")
    (root / "intermediate" / "stages" / "08_virulence.tsv").unlink()
    client, _state = app_factory(results_root=root)
    row = client.get("/api/isolates").json()["items"][0]
    assert row["n_virulence"] == "not assessed"
    assert row["st"] not in ("", 0, None)
    for column in ("amr_genes", "oprd_state", "n_snv"):
        assert row[column] not in ("", None)


def test_present_table_with_no_row_is_a_measurement(app_factory, tmp_path):
    # Manifest names a sample that appears in no table: the present virulence
    # table holds no row for it, so 0 is a real measurement; a missing table
    # would be `not assessed`.
    root = copy_layout(SMALL_LIVE, tmp_path / "extra_sample")
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["samples"] = list(manifest["samples"]) + ["TEST_PA_999"]
    manifest["n_samples"] = len(manifest["samples"])
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    client, _state = app_factory(results_root=root)
    items = {row["sample_id"]: row for row in client.get("/api/isolates").json()["items"]}
    extra = items["TEST_PA_999"]
    assert extra["n_virulence"] == 0
    assert extra["amr_genes"] == "not assessed"


# ---------------------------------------------------------------------------
# UI-D4 — path safety and bind/token
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "requested",
    ["../../etc/passwd", "/etc/passwd", "~/etc/passwd", "sub/../../etc/passwd"],
)
def test_traversal_and_absolute_paths_are_refused(small_live, requested):
    with pytest.raises(ContainmentError):
        resolve_contained(small_live, requested)


def test_encoded_traversal_is_refused_over_http(live_client):
    client, _state = live_client
    response = client.get("/api/reports/content", params={"path": "../../etc/passwd"})
    assert response.status_code == 400
    response = client.get("/api/reports/content", params={"path": "%2e%2e%2f%2e%2e%2fetc%2fpasswd"})
    assert response.status_code in (400, 404)


def test_symlink_escape_is_refused(app_factory, tmp_path):
    root = copy_layout(SMALL_LIVE, tmp_path / "symlinked")
    link = root / "reports" / "escape.md"
    link.symlink_to("/etc/passwd")
    client, _state = app_factory(results_root=root)
    response = client.get("/api/reports/content", params={"path": "reports/escape.md"})
    assert response.status_code == 400
    with pytest.raises(ContainmentError):
        resolve_contained(root, "reports/escape.md")


def test_default_bind_is_loopback_and_needs_no_token(monkeypatch):
    monkeypatch.delenv("PA_DASH_TOKEN", raising=False)
    policy = resolve_bind(None, None)
    assert policy.host == "127.0.0.1"
    assert policy.loopback is True
    assert policy.token_required is False


def test_wildcard_bind_is_always_refused(monkeypatch):
    monkeypatch.setenv("PA_DASH_TOKEN", "secret")
    with pytest.raises(BindRefused):
        resolve_bind("0.0.0.0", None)


def test_non_local_bind_requires_a_token(monkeypatch):
    monkeypatch.delenv("PA_DASH_TOKEN", raising=False)
    with pytest.raises(BindRefused):
        resolve_bind("192.168.1.20", None)
    monkeypatch.setenv("PA_DASH_TOKEN", "s3cret")
    policy = resolve_bind("192.168.1.20", None)
    assert policy.loopback is False
    assert policy.token_required is True


def test_token_required_on_a_non_local_bind_and_header_only(app_factory, tmp_path, monkeypatch):
    from fastapi.testclient import TestClient

    monkeypatch.setenv("PA_DASH_TOKEN", "s3cret")
    client, _state = app_factory(results_root=SMALL_LIVE, host="192.168.1.20", state_dir=tmp_path / "s")
    # `raise_server_exceptions=False` so a refusal is a response, not a raise.
    guarded = TestClient(client.app, raise_server_exceptions=False)
    assert guarded.get("/api/stages").status_code == 401
    assert guarded.get("/api/stages", headers={"X-PA-Dash-Token": "wrong"}).status_code == 401
    assert guarded.get("/api/stages", headers={"X-PA-Dash-Token": "s3cret"}).status_code == 200
    # Never accepted as a query parameter.
    assert guarded.get("/api/stages", params={"X-PA-Dash-Token": "s3cret"}).status_code == 401


# ---------------------------------------------------------------------------
# UI-D5 — a huge table is not fully loaded
# ---------------------------------------------------------------------------


@pytest.fixture(scope="function")
def wide_client(app_factory, tmp_path):
    from dashboard.fixtures import generate as gen

    root = gen.build_wide_root(tmp_path / "wide", rows=3000, cols=2000)
    client, _state = app_factory(results_root=root, state_dir=tmp_path / "wide_state")
    return client, root / "intermediate" / "stages" / "variants.tsv"


def test_d5_seam_reports_a_page_far_smaller_than_the_file(wide_client):
    client, table = wide_client
    # Warm the byte-offset index once, then page.
    client.get("/api/stages/variants/rows", params={"limit": 20})
    response = client.get("/api/stages/variants/rows", params={"limit": 20})
    assert response.status_code == 200
    read = response.json()["meta"]["read"]
    file_size = table.stat().st_size
    assert read["file_size"] == file_size
    assert read["fully_loaded"] is False
    assert read["bytes_read"] < file_size * 0.02
    assert "byte-range" in read["mechanism"]


def test_d5_paged_fetch_does_not_fully_load_the_table(wide_client):
    """The honest memory check: the endpoint must not materialise the table.

    `meta.read` reports only the page. This test measures the whole request's
    peak allocation against a full `read_tsv` baseline, so a page that quietly
    reads the file through `report_tables.read_table` first is caught rather
    than papered over.
    """
    from papipeline.execution.contracts import required_columns
    from papipeline.io.tsv import read_tsv

    client, table = wide_client
    client.get("/api/stages/variants/rows", params={"limit": 20})  # warm

    tracemalloc.start()
    client.get("/api/stages/variants/rows", params={"limit": 20})
    _current, endpoint_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # The baseline is what the endpoint's own `read_table` call costs: the same
    # `read_tsv` with the same required columns.
    tracemalloc.start()
    try:
        read_tsv(table, required_columns=list(required_columns("variants")))
    except Exception:
        pass
    _current, full_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert endpoint_peak < full_peak * 0.5, (
        f"the paged request allocated {endpoint_peak/1e6:.1f} MB against a full "
        f"read's {full_peak/1e6:.1f} MB: the endpoint materialised the table "
        f"before paging (dashboard/server/app.py:1275 calls `stage_table` -> "
        f"`read_table` -> `read_tsv`, which loads every row)"
    )


# ---------------------------------------------------------------------------
# UI-D6 — duplicate-label tree
# ---------------------------------------------------------------------------


def test_tree_tips_and_distances_match_a_dendropy_oracle(small_live):
    text = (small_live / "intermediate" / "phylogeny" / "stage9.nwk").read_text(encoding="utf-8")
    root = newick_module.parse(text)
    ours_tips = newick_module.tip_labels(root)
    assert ours_tips == [f"TEST_PA_{i:03d}" for i in range(1, 11)]

    dendropy = pytest.importorskip("dendropy")
    tree = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)
    taxa = list(tree.taxon_namespace)
    assert sorted(taxon.label for taxon in taxa) == sorted(ours_tips)

    pdm = tree.phylogenetic_distance_matrix()
    ours = newick_module.patristic_distances(root)
    for left in taxa:
        for right in taxa:
            expected = pdm(left, right)
            assert ours[left.label][right.label] == pytest.approx(expected, abs=1e-9)


def test_underscore_default_is_a_trap_we_avoid(small_live):
    text = (small_live / "intermediate" / "phylogeny" / "stage9.nwk").read_text(encoding="utf-8")
    root = newick_module.parse(text)
    assert all("_" in label for label in newick_module.tip_labels(root))

    dendropy = pytest.importorskip("dendropy")
    mangled = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=False)
    assert any(" " in taxon.label for taxon in mangled.taxon_namespace)


def test_tree_endpoint_reports_the_tip_set(live_client):
    client, _state = live_client
    payload = client.get("/api/tree").json()
    assert payload["present"] is True
    assert payload["n_tips"] == 10
    assert payload["tip_check"]["matched"] is True
    # Duplicate support labels are display-only, never identities.
    node_ids = [node["node_id"] for node in payload["nodes"]]
    assert len(node_ids) == len(set(node_ids))


# ---------------------------------------------------------------------------
# oprD (§4.3)
# ---------------------------------------------------------------------------


def test_oprd_returns_not_produced_with_the_pipeline_reason(live_client):
    from papipeline.stages.reporting import OPRD_NOT_ON_DISK

    client, _state = live_client
    payload = client.get("/api/oprd").json()
    assert payload["available"] is False
    assert payload["reason"] == OPRD_NOT_ON_DISK
    assert payload["probes"]
    assert [probe["n"] for probe in payload["probes"]] == list(range(1, len(payload["probes"]) + 1))
    for item in payload["items"]:
        assert item["verdict"] != "absent"
        assert item["display_state"] != "absent"


def test_isolate_oprd_state_is_never_absent(live_client):
    client, _state = live_client
    for row in client.get("/api/isolates").json()["items"]:
        assert row["oprd_state"] in ("intact", "disrupted", "not_assessed")
        assert row["oprd_state"] != "absent"


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def test_bakta_executions_absent_reads_not_reported(live_client):
    client, _state = live_client
    provenance = client.get("/api/provenance").json()
    assert provenance["bakta_executions"] is None
    assert provenance["bakta_executions_label"] == "not reported"
    bakta = client.get("/api/provenance/bakta").json()
    assert bakta["bakta_executions"] is None
    assert bakta["bakta_executions_label"] == "not reported"
    assert bakta["bakta_executions"] != 0
