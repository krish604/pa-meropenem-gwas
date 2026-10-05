"""Q3 — performance gates on the 900-isolate set.

The 900-isolate fixture is generated on demand and never committed. A normal
run skips these tests; set ``PA_FIXTURES_LARGE=1`` to run them. The measured
timings are recorded in ``pa-artifacts/ui/QA.md``.
"""

from __future__ import annotations

import os
import statistics
import time

import pytest

pytestmark = pytest.mark.skipif(
    os.environ.get("PA_FIXTURES_LARGE") != "1",
    reason="set PA_FIXTURES_LARGE=1 to generate the 900-isolate set and run the perf gates",
)

ISOLATES_GATE_S = 2.0
TREE_GATE_S = 3.0
HEATMAP_GATE_S = 3.0


@pytest.fixture(scope="module")
def large_client(tmp_path_factory):
    from fastapi.testclient import TestClient

    from dashboard.fixtures import generate as gen
    from dashboard.server.app import build_state, create_app

    base = tmp_path_factory.mktemp("large")
    root = gen.build_large(base / "large", n=900)
    state = build_state(
        results_root=root,
        repo_root=gen.REPO_ROOT,
        state_dir=base / "state",
    )
    return TestClient(create_app(state))


def _median_seconds(client, url, params=None, runs=3, warm=True):
    if warm:
        client.get(url, params=params)
    samples = []
    for _ in range(runs):
        start = time.perf_counter()
        response = client.get(url, params=params)
        samples.append(time.perf_counter() - start)
        assert response.status_code == 200
    return statistics.median(samples)


def test_isolate_table_first_page_under_two_seconds(large_client):
    seconds = _median_seconds(large_client, "/api/isolates", {"limit": 200})
    print(f"\nQ3 /api/isolates 900: {seconds * 1000:.0f} ms")
    assert seconds < ISOLATES_GATE_S, f"{seconds:.3f}s >= {ISOLATES_GATE_S}s"


def test_tree_json_under_three_seconds(large_client):
    seconds = _median_seconds(large_client, "/api/tree")
    print(f"\nQ3 /api/tree 900: {seconds * 1000:.0f} ms")
    assert seconds < TREE_GATE_S, f"{seconds:.3f}s >= {TREE_GATE_S}s"


def test_heatmap_data_under_three_seconds(large_client):
    seconds = _median_seconds(large_client, "/api/similarity")
    print(f"\nQ3 /api/similarity 900 (default max_tips=300): {seconds * 1000:.0f} ms")
    assert seconds < HEATMAP_GATE_S, f"{seconds:.3f}s >= {HEATMAP_GATE_S}s"
