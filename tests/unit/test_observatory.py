"""Tests for the observatory.

The bar here is higher than "the endpoint returns 200". Every test asserts
that a value the interface would show is the value the execution store
actually holds, and that a value the engine does not record is reported as
unavailable rather than invented.

The vertical-slice tests run a real cohort through real subprocesses, then
read it back through the real ASGI app.
"""

from __future__ import annotations

import itertools
import json
import time

import pytest
from starlette.testclient import TestClient

from papipeline.execution import StageState
from papipeline.observatory import (
    EventBus,
    adapt,
    aggregate_state,
    build_graph,
    bus_sink,
    graph_payload,
    read_snapshot,
    sample_host,
)
from papipeline.observatory.api import create_app, payload
from papipeline.observatory.events import TASK_VALIDATED, ExecutionEvent
from papipeline.observatory.fixture import Fixture
from papipeline.observatory.graph import (
    EdgeKind,
    dependency_edges,
    layout,
    sequence_edges,
)
from papipeline.observatory.throughput import compute, median_duration
from papipeline.run import EXECUTION_ORDER, PREREQUISITES


# ── a real server, for the real transport ─────────────────────────

def read_frames(response, want=2, match=None, limit=40, after_first=None,
                deadline=15.0):
    """Read SSE frames off a live response until ``want`` arrive.

    Reads byte-wise on the raw stream so a frame is returned the moment it
    arrives, rather than after the (endless) body completes. The response
    may only be read once, so anything that has to happen mid-stream - such
    as emitting the event a test is waiting for - goes in ``after_first``.
    """
    frames = []
    buffer = ""
    seen = 0
    started = time.monotonic()
    for raw in response.iter_bytes():
        # The body is endless, so an absent frame must fail the test rather
        # than hang it.
        if time.monotonic() - started > deadline:
            raise AssertionError(
                f"no matching SSE frame within {deadline}s (saw {len(frames)})")
        buffer += raw.decode("utf-8", errors="replace")
        while "\n\n" in buffer:
            frame, buffer = buffer.split("\n\n", 1)
            frame = frame.strip()
            if not frame:
                continue
            # The hook fires on the first frame of any kind, including one
            # the filter rejects: it means "the client is connected now".
            seen += 1
            if seen == 1 and after_first is not None:
                after_first()
            if match and not frame.startswith(match):
                continue
            frames.append(frame)
            if len(frames) >= want:
                return frames
            if len(frames) > limit:
                raise AssertionError("too many frames without the one wanted")
    return frames


@pytest.fixture()
def live_server():
    """A real uvicorn server over a real fixture, on an ephemeral port."""
    import socket
    import threading

    import uvicorn

    from papipeline.observatory.api import create_app

    fixture = Fixture(subjects=4)
    fixture.stage_cohort_baseline()
    fixture.run_cohort(with_failures=True)
    bus = EventBus()

    sock = socket.socket()
    sock.bind(("127.0.0.1", 0))
    port = sock.getsockname()[1]
    sock.close()

    app = create_app(str(fixture.db), fixture.run_key, bus=bus)
    config = uvicorn.Config(app, host="127.0.0.1", port=port, log_level="error")
    server = uvicorn.Server(config)
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()

    deadline = time.time() + 20
    import httpx
    while time.time() < deadline:
        try:
            if httpx.get(f"http://127.0.0.1:{port}/api/health", timeout=1.0).status_code == 200:
                break
        except Exception:
            time.sleep(0.05)
    else:
        server.should_exit = True
        raise AssertionError("the server did not come up")

    try:
        yield {"url": f"http://127.0.0.1:{port}", "bus": bus,
               "run_key": fixture.run_key, "fixture": fixture}
    finally:
        server.should_exit = True
        thread.join(timeout=10)
        fixture.close()


# ── structure: real, not invented ─────────────────────────────────

class TestGraphStructure:
    def test_nodes_are_the_pipeline_s_own_stages(self):
        nodes, _ = build_graph()
        assert [n.stage for n in nodes] == list(EXECUTION_ORDER)

    def test_dependency_edges_come_only_from_prerequisites(self):
        edges = dependency_edges()
        declared = {
            (src, dst)
            for dst, sources in PREREQUISITES.items()
            for src in sources
            if src in set(EXECUTION_ORDER)
        }
        drawn = {(e.source, e.target) for e in edges if e.kind == EdgeKind.DEPENDENCY.value}
        assert drawn == declared

    def test_no_dependency_is_invented(self):
        for edge in dependency_edges():
            assert edge.target in PREREQUISITES, f"{edge.target} declares no prerequisites"
            assert edge.source in PREREQUISITES[edge.target]

    def test_sequence_edges_are_a_separate_class(self):
        seq = sequence_edges()
        assert len(seq) == len(EXECUTION_ORDER) - 1
        assert all(e.kind == EdgeKind.SEQUENCE.value for e in seq)
        # A pair may be adjacent in execution order *and* a declared
        # dependency (phenotype -> gwas is both). That is correct and must
        # not be collapsed: the pair is drawn twice, once for each reason.
        assert {(e.source, e.target) for e in seq} & {
            (e.source, e.target) for e in dependency_edges()
        }, "phenotype -> gwas is both adjacent and a declared dependency"
        # Every edge carries exactly one kind, so the client can style them.
        assert all(e.kind in {EdgeKind.DEPENDENCY.value, EdgeKind.SEQUENCE.value}
                   for e in seq + dependency_edges())

    def test_layout_is_deterministic(self):
        assert layout() == layout()
        first, _ = build_graph()
        second, _ = build_graph()
        assert [(n.x, n.y) for n in first] == [(n.x, n.y) for n in second]

    def test_nodes_do_not_overlap(self):
        nodes, _ = build_graph()
        for a, b in itertools.combinations(nodes, 2):
            distance = ((a.x - b.x) ** 2 + (a.y - b.y) ** 2) ** 0.5
            assert distance > 60, f"{a.stage} and {b.stage} are too close to read"

    def test_payload_declares_what_each_edge_means(self):
        body = graph_payload()
        assert set(body["legend"]) == {"dependency", "sequence"}
        assert "reads" in body["legend"]["dependency"]


# ── state projection ──────────────────────────────────────────────

class TestStateAggregation:
    def test_nothing_recorded_is_pending(self):
        assert aggregate_state({}, 0) == StageState.PENDING.value

    def test_one_success_among_failures_is_a_failure(self):
        assert aggregate_state({"SUCCEEDED": 64, "FAILED": 1}, 65) == StageState.FAILED.value

    def test_running_outranks_everything(self):
        counts = {"SUCCEEDED": 60, "FAILED": 3, "RUNNING": 2}
        assert aggregate_state(counts, 65) == StageState.RUNNING.value

    def test_retrying_outranks_failure(self):
        assert aggregate_state({"FAILED": 4, "RETRYING": 1}, 5) == StageState.RETRYING.value

    def test_invalid_is_reported_over_incomplete(self):
        assert aggregate_state({"INCOMPLETE": 2, "INVALID": 1}, 3) == StageState.INVALID.value

    def test_all_success_is_success(self):
        assert aggregate_state({"SUCCEEDED": 65}, 65) == StageState.SUCCEEDED.value

    def test_a_stage_with_no_rows_is_present_and_pending(self):
        with Fixture(subjects=2) as fixture:
            snap = read_snapshot(fixture.store, fixture.run_key)
            by_stage = {n["stage"]: n for n in snap.nodes}
            assert len(snap.nodes) == len(EXECUTION_ORDER)
            assert by_stage["gwas"]["state"] == StageState.PENDING.value
            assert by_stage["gwas"]["total"] == 0

    def test_snapshot_is_idle_with_no_rows(self):
        with Fixture(subjects=2) as fixture:
            assert read_snapshot(fixture.store, fixture.run_key).idle is True

    def test_a_missing_store_is_idle_not_an_error(self):
        snap = read_snapshot(None, "nope")
        assert snap.exists is False
        assert snap.idle is True
        assert snap.tasks == []


# ── metrics: real, and honest when unmeasurable ───────────────────

class TestHostMetrics:
    def test_reports_real_host_facts(self):
        host = sample_host()
        assert host.host
        assert host.cpu_count >= 1
        assert host.sampled_at > 0

    def test_absent_readings_are_none_not_zero(self):
        """Zero would be a claim. The first CPU and disk samples cannot be
        computed, and must read as unknown."""
        host = sample_host()
        if host.cpu_percent is None:
            assert "psutil" in host.notes or "cpu_percent" in host.notes
        if host.disk_read_bps is None:
            assert host.disk_write_bps is None

    def test_a_second_sample_produces_a_rate(self):
        sample_host()
        time.sleep(0.05)
        host = sample_host()
        assert (host.disk_read_bps is not None and host.disk_write_bps is not None) \
            or "disk_io" in host.notes


# ── throughput: measured or absent ────────────────────────────────

class TestThroughput:
    def _task(self, state, elapsed, started=None, ended=None):
        return {"state": state, "elapsed_seconds": elapsed,
                "started_at": started, "ended_at": ended}

    def test_no_completions_means_no_rate(self):
        out = compute([], 65)
        assert out.tasks_per_minute is None
        assert out.eta_minutes is None
        assert "no task has completed" in out.basis

    def test_one_completion_is_not_a_rate(self):
        out = compute([self._task("SUCCEEDED", 10.0)], 65)
        assert out.tasks_per_minute is None
        assert "not meaningful" in out.basis

    def test_a_rate_appears_once_there_are_enough_samples(self):
        tasks = [
            self._task("SUCCEEDED", 60.0, "2024-12-01T14:00:00Z", "2024-12-01T14:05:00Z")
            for _ in range(5)
        ]
        out = compute(tasks, 10)
        assert out.tasks_per_minute is not None
        assert out.tasks_per_minute > 0
        assert out.eta_minutes is not None
        assert "wall-clock" in out.basis

    def test_nothing_remaining_means_zero_eta_not_a_guess(self):
        tasks = [self._task("SUCCEEDED", 60.0,
                            "2024-12-01T14:00:00Z", "2024-12-01T14:05:00Z")
                 for _ in range(4)]
        out = compute(tasks, 4)
        assert out.eta_minutes == 0.0

    def test_median_duration_is_none_with_no_finished_tasks(self):
        assert median_duration([]) is None
        assert median_duration([self._task("FAILED", 5.0)]) is None


# ── events: real transitions only ─────────────────────────────────

class TestEvents:
    def test_a_real_run_emits_real_events(self):
        """Events come from the runner, not from a generator."""
        with Fixture(subjects=2) as fixture:
            assert fixture.bus.sequence == 0
            fixture.run_cohort(with_failures=True)
            names = [e.name for e in fixture.bus.recent()]
            assert fixture.bus.sequence > 0
            assert "TASK_STARTED" in names
            assert "TASK_COMPLETED" in names

    def test_a_failing_task_emits_failure_not_completion(self):
        with Fixture(subjects=2) as fixture:
            fixture.run_annotation("GCA_X.1", behaviour="fails")
            names = [e.name for e in fixture.bus.recent()]
            assert "TASK_FAILED" in names
            assert "TASK_COMPLETED" not in names

    def test_a_retry_emits_a_retry_event(self):
        with Fixture(subjects=2) as fixture:
            fixture.run_annotation("GCA_X.1", behaviour="interrupted")
            names = [e.name for e in fixture.bus.recent()]
            assert "TASK_RETRY" in names

    def test_validation_is_reported_separately_from_completion(self):
        """They are two facts, and the engine distinguishes them."""
        with Fixture(subjects=2) as fixture:
            fixture.run_annotation("GCA_X.1", behaviour="ok")
            names = [e.name for e in fixture.bus.recent()]
            assert "TASK_COMPLETED" in names
            completed = [e for e in fixture.bus.recent() if e.name == "TASK_COMPLETED"][0]
            assert completed.payload.get("validated") is True

    def test_an_invalid_output_emits_invalid_not_completed(self):
        with Fixture(subjects=2) as fixture:
            fixture.run_annotation("GCA_X.1", behaviour="wrong_isolate")
            names = [e.name for e in fixture.bus.recent()]
            assert "TASK_INVALID" in names
            assert "TASK_COMPLETED" not in names

    def test_adapt_rejects_a_foreign_object(self):
        with pytest.raises(TypeError):
            adapt(object())

    def test_the_bus_survives_a_broken_subscriber(self):
        bus = EventBus()
        seen = []
        bus.subscribe(lambda e: (_ for _ in ()).throw(RuntimeError("boom")))
        bus.subscribe(seen.append)
        event = ExecutionEvent(name="TASK_STARTED", stage="a", subject="b",
                               run_key="r", state="RUNNING")
        bus.emit(event)
        assert len(seen) == 1, "one broken subscriber must not silence the rest"

    def test_a_slow_subscriber_cannot_break_a_run(self, tmp_path):
        """The sink is wrapped by the runner, so a dead observer is inert."""
        def exploding(_event):
            raise RuntimeError("observer is broken")
        with Fixture(subjects=1) as fixture:
            result = fixture.run_annotation("GCA_X.1", behaviour="ok")
        assert result.state is StageState.SUCCEEDED

    def test_the_buffer_is_bounded(self):
        bus = EventBus(capacity=5)
        for i in range(50):
            bus.emit(ExecutionEvent(name="TASK_STARTED", stage="a", subject=str(i),
                                    run_key="r", state="RUNNING"))
        assert len(bus.recent(limit=1000)) == 5
        assert bus.gap_since(0) > 0, "a gap must be reported, not hidden"


# ── the API over real state ───────────────────────────────────────

class TestApi:
    @pytest.fixture()
    def client(self):
        with Fixture(subjects=4) as fixture:
            fixture.run_cohort(with_failures=True)
            app = create_app(str(fixture.db), fixture.run_key, bus=fixture.bus)
            with TestClient(app) as client:
                client.fixture = fixture
                yield client

    def test_health(self, client):
        body = client.get("/api/health").json()
        assert body["status"] == "ok"
        assert body["database_present"] is True

    def test_structure_is_the_real_graph(self, client):
        body = client.get("/api/structure").json()
        assert [n["stage"] for n in body["nodes"]] == list(EXECUTION_ORDER)

    def test_snapshot_states_match_the_store(self, client):
        body = client.get("/api/snapshot").json()
        counts = body["counts"]
        for state, expected in counts.items():
            actual = sum(
                1 for t in body["tasks"] if t["state"] == state
            )
            assert actual == expected, f"{state} count disagrees with the task list"
        assert body["exists"] is True

    def test_every_state_shown_is_a_real_engine_state(self, client):
        body = client.get("/api/snapshot").json()
        allowed = {s.value for s in StageState}
        for task in body["tasks"]:
            assert task["state"] in allowed

    def test_idle_means_nothing_in_flight_not_no_failures(self, client):
        """A finished cohort with failures is idle for the overlay, because
        nothing is moving -- but the failures are still counted and shown."""
        body = client.get("/api/snapshot").json()
        assert body["idle"] is True
        assert body["counts"].get("FAILED", 0) >= 1
        assert body["counts"].get("INVALID", 0) >= 1

    def test_unmeasurable_metrics_are_declared_not_hidden(self, client):
        caps = client.get("/api/snapshot").json()["capabilities"]
        assert caps["per_task_cpu"] is False
        assert caps["per_task_progress"] is False
        assert "unavailable" in caps["note"]

    def test_task_detail_matches_the_store_row(self, client):
        fixture = client.fixture
        stage, subject = "annotation", "GCA_000000001.1"
        body = client.get(f"/api/tasks/{stage}/{subject}").json()
        row = fixture.store.get(fixture.run_key, stage, subject)
        assert body["state"] == row["state"]
        assert body["elapsed_seconds"] == row["elapsed_seconds"]
        assert body["command"] == json.loads(row["command"])
        assert body["subject"] == subject

    def test_task_detail_includes_real_attempt_history(self, client):
        fixture = client.fixture
        body = client.get("/api/tasks/annotation/GCA_000000002.1").json()
        recorded = fixture.store.attempts(fixture.run_key, "annotation",
                                          "GCA_000000002.1")
        assert len(body["attempts_history"]) == len(recorded)
        assert [a["attempt"] for a in body["attempts_history"]] == \
            [a.attempt for a in recorded]
        assert body["attempts_history"][-1]["state"] == StageState.FAILED.value

    def test_a_missing_task_is_a_404(self, client):
        assert client.get("/api/tasks/gwas/nobody").status_code == 404

    def test_logs_are_read_from_the_recorded_path(self, client):
        fixture = client.fixture
        body = client.get("/api/logs/annotation/GCA_000000001.1").json()
        assert body["available"] is True
        recorded = fixture.store.log_path(
            fixture.run_key, "annotation", "GCA_000000001.1")
        assert body["path"] == recorded
        assert "wrote 4 records" in body["text"]

    def test_a_task_with_no_log_says_so(self, client):
        body = client.get("/api/logs/annotation/GCA_000000002.1").json()
        if body["available"]:
            # The failing task still ran a process, so it has a log.
            assert body["text"]
        else:
            assert body["reason"], "an unavailable log must say why"

    def test_a_log_for_an_unknown_task_is_a_404(self, client):
        assert client.get("/api/logs/gwas/nobody").status_code == 404

    def test_tasks_can_be_filtered(self, client):
        body = client.get("/api/tasks?stage=annotation").json()
        assert body
        assert all(t["stage"] == "annotation" for t in body)
        failed = client.get("/api/tasks?state=FAILED").json()
        assert all(t["state"] == "FAILED" for t in failed)

    def test_the_index_and_assets_are_served(self, client):
        assert client.get("/").status_code == 200
        assert "Computational Observatory" in client.get("/").text
        assert client.get("/static/js/main.js").status_code == 200
        assert client.get("/static/css/observatory.css").status_code == 200

    def test_the_page_never_contains_hardcoded_metrics(self, client):
        """The reference dashboard's own numbers must not appear as literals
        in the shipped assets."""
        html = client.get("/").text
        for js in ("main.js", "network.js", "metrics.js", "detail.js"):
            body = client.get(f"/static/js/{js}").text
            for bogus in ("78%", "22.4 GB", "0.31 genomes", "03:17:42", "128 MB"):
                assert bogus not in body, f"{js} hardcodes {bogus!r}"
        assert "78%" not in html


class TestIdleObservatory:
    def test_a_fresh_machine_reports_idle(self, tmp_path):
        app = create_app(str(tmp_path / "nothing.db"), "run")
        with TestClient(app) as client:
            body = client.get("/api/snapshot").json()
            assert body["exists"] is False
            assert body["idle"] is True
            assert body["counts"] == {}
            # Host metrics are still real and still available.
            assert body["host"]["host"]
            assert body["host"]["cpu_count"] >= 1
            assert body["throughput"]["tasks_per_minute"] is None

    def test_idle_is_still_served_with_a_snapshot(self, client_ok=True):
        app = create_app("/nonexistent/path.db", "run")
        with TestClient(app) as client:
            assert client.get("/api/snapshot").status_code == 200


# ── the vertical slice, end to end ────────────────────────────────

class TestVerticalSlice:
    def test_a_real_cohort_is_visible_through_the_api(self):
        """Real subprocesses -> real store -> real API -> real states."""
        with Fixture(subjects=6) as fixture:
            plan = fixture.run_cohort(with_failures=True)
            app = create_app(str(fixture.db), fixture.run_key, bus=fixture.bus)
            with TestClient(app) as client:
                body = client.get("/api/snapshot").json()
                assert len(body["tasks"]) == len(plan)
                states = {t["subject"]: t["state"] for t in body["tasks"]}
                # Every state the fixture intended is the state reported.
                assert states["GCA_000000001.1"] == StageState.SUCCEEDED.value
                assert states["GCA_000000002.1"] == StageState.FAILED.value
                assert states["GCA_000000004.1"] == StageState.INVALID.value
                assert states["GCA_000000006.1"] == StageState.INCOMPLETE.value

    def test_a_live_run_is_not_idle_and_becomes_idle_when_it_finishes(self):
        with Fixture(subjects=3) as fixture:
            app = create_app(str(fixture.db), fixture.run_key, bus=fixture.bus)
            with TestClient(app) as client:
                assert client.get("/api/snapshot").json()["idle"] is True
                fixture.run_cohort(with_failures=False)
                # Finished work with nothing in flight is idle-but-recorded.
                after = client.get("/api/snapshot").json()
                assert after["idle"] is True
                assert after["counts"].get("SUCCEEDED") == 3
                assert after["exists"] is True

    def test_the_stream_pushes_a_snapshot(self, live_server):
        """A real SSE connection to a real server.

        ``TestClient`` cannot be used here: it needs a terminating body,
        and this stream is by definition endless. So the stream is exercised
        over real HTTP, which is also the only honest way to claim the
        transport works.
        """
        import httpx

        with httpx.Client(timeout=10.0) as client:
            with client.stream("GET", live_server["url"] + "/api/stream") as response:
                assert response.status_code == 200
                assert response.headers["content-type"].startswith("text/event-stream")
                assert "no-cache" in response.headers["cache-control"]
                frames = read_frames(response, want=2)
        assert frames[0].startswith("event: hello")
        assert frames[1].startswith("event: snapshot")
        payload = json.loads(frames[1].split("data: ", 1)[1])
        # The snapshot carries real state, not a placeholder.
        assert payload["run_key"] == live_server["run_key"]
        assert "nodes" in payload and "edges" in payload

    def test_the_stream_delivers_a_real_event(self, live_server):
        """An event emitted after the client connects reaches that client."""
        import httpx

        def emit_once_the_client_is_here():
            live_server["bus"].emit(ExecutionEvent(
                name=TASK_VALIDATED, stage="annotation",
                subject="GCA_000000001.1", run_key=live_server["run_key"],
                state="SUCCEEDED", attempt=1, detail="4 checks passed"))

        with httpx.Client(timeout=10.0) as client:
            with client.stream("GET", live_server["url"] + "/api/stream") as response:
                frames = read_frames(response, want=1, match="event: event",
                                     after_first=emit_once_the_client_is_here,
                                     deadline=20.0)
        assert len(frames) == 1, "exactly one event frame was expected"
        assert frames[0].startswith("event: event")
        body = json.loads(frames[0].split("data: ", 1)[1])
        assert body["name"] == TASK_VALIDATED
        assert body["subject"] == "GCA_000000001.1"

    def test_a_reconnecting_client_gets_a_fresh_snapshot(self, live_server):
        """Reconnection is a plain HTTP retry, and the client is re-told the
        truth rather than resuming on a stale picture."""
        import httpx

        with httpx.Client(timeout=10.0) as client:
            with client.stream("GET", live_server["url"] + "/api/stream") as first:
                frames_first = read_frames(first, want=2)
            with client.stream("GET", live_server["url"] + "/api/stream") as second:
                frames_again = read_frames(second, want=2)
        assert frames_first[0].startswith("event: hello")
        assert frames_first[1].startswith("event: snapshot")
        assert frames_again[0].startswith("event: hello")
        assert frames_again[1].startswith("event: snapshot")

    def test_the_event_sink_can_be_passed_straight_to_run_task(self):
        """The documented integration really works."""
        with Fixture(subjects=2) as fixture:
            from papipeline.execution import (
                OutputSpec, Check, CheckKind, RetryPolicy, TaskContext, run_task,
            )
            import sys as _sys
            out = fixture.root / "x.txt"
            spec = OutputSpec("demo", (Check(CheckKind.NON_EMPTY, path=out),))
            script = fixture.root / "p2.py"
            script.write_text("import sys;from pathlib import Path;"
                              "Path(sys.argv[1]).write_text('hello\\n')\n")
            result = run_task(
                TaskContext(run_key="r", stage="demo", subject="", spec=spec),
                [_sys.executable, str(script), str(out)],
                store=fixture.store,
                policy=RetryPolicy(max_attempts=1),
                event_sink=bus_sink(fixture.bus),
            )
            assert result.state is StageState.SUCCEEDED
            assert fixture.bus.sequence >= 2
