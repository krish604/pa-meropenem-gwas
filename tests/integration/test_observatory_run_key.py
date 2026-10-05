"""The observatory must not render an empty pipeline when it is looking at the
wrong run key.

**The bug this locks down.** `python -m papipeline.observatory` defaulted to the
module constant `DEFAULT_RUN_KEY = "observatory"`, while every machine overlay
configures `observatory.run_key: pipeline` and every real run writes under that
key. A default-started server therefore asked a database full of real rows for a
run that never existed, and rendered all sixteen stages `PENDING` - with
`/api/health` reporting `status: ok`. Against the smoke run this looked exactly
like a pipeline that had not started. Nothing errored.

That is the failure worth a test: not a crash but a *confident wrong answer*.
A reader who saw sixteen PENDING nodes would conclude the run had not begun,
rather than that they were watching the wrong key.

Fixed in two places, and both are needed:

* `__main__` reads the overlay via `ObservatorySettings.from_config` instead of
  the constants, so a default start resolves the key the writer actually uses.
* The API reports `available_run_keys` / `run_key_found` (health) and
  `run_key_mismatch` (snapshot), so a *remaining* mismatch - a `--run` typo, or
  a server pointed at the wrong database - is stated rather than rendered.

The second half is the durable one. Reading the overlay fixes today's default;
nothing stops the next person passing `--run` a key that does not exist, so the
API has to be able to say so.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from papipeline.execution.store import ExecutionStore
from papipeline.execution.state import StageState
from papipeline.observatory.api import DEFAULT_RUN_KEY, create_app
from papipeline.observatory.__main__ import main


def _store_with_run(tmp_path: Path, run_key: str = "pipeline") -> Path:
    """A database holding one finished-looking run under `run_key`."""
    db = tmp_path / "observatory.db"
    store = ExecutionStore(db)
    store.record(
        run_key=run_key,
        stage="amr",
        state=StageState.SUCCEEDED,
        attempt=1,
        max_attempts=1,
        elapsed_seconds=1986.0,
    )
    store.close()
    return db


class TestTheDefaultComesFromTheOverlay:
    def test_the_default_run_key_is_not_the_module_constant(self, tmp_path):
        """A default start must resolve the overlay's key, not `observatory`.

        Asserted against the constant rather than by hard-coding the string, so
        that changing the constant cannot quietly make this test pass for the
        wrong reason.
        """
        captured = {}
        with mock.patch("uvicorn.run") as uvicorn:
            main([])
        app = uvicorn.call_args.args[0]
        captured["run"] = app.state.run_key
        assert captured["run"] != DEFAULT_RUN_KEY, (
            f"a default-started server resolved the module constant "
            f"{DEFAULT_RUN_KEY!r} instead of reading the overlay, so it cannot "
            "see any run a machine overlay writes"
        )
        assert captured["run"] == "pipeline", (
            f"expected the overlay's run_key 'pipeline', got {captured['run']!r}"
        )

    def test_an_explicit_run_still_wins(self, tmp_path):
        """The overlay is a default, not an override of the operator."""
        with mock.patch("uvicorn.run") as uvicorn:
            main(["--run", "explicit-key"])
        assert uvicorn.call_args.args[0].state.run_key == "explicit-key"


class TestARealRunIsVisibleByDefault:
    def test_a_default_server_finds_a_run_under_the_overlay_key(self, tmp_path):
        """The original symptom: real rows present, all-PENDING reported."""
        db = _store_with_run(tmp_path, run_key="pipeline")
        with mock.patch("uvicorn.run") as uvicorn:
            main([])
        app = uvicorn.call_args.args[0]
        app.state.db_path = str(db)

        body = _snapshot(app)
        states = {n["stage"]: n["state"] for n in body["nodes"]}
        assert states["amr"] != "PENDING", (
            "the server reported a completed run as PENDING - the run key it "
            "asked for does not match the one the run was written under"
        )
        assert body["run_key_mismatch"]["mismatch"] is False

    def test_health_reports_the_run_key_was_found(self, tmp_path):
        db = _store_with_run(tmp_path, run_key="pipeline")
        app = create_app(str(db), "pipeline")
        health = _get(app, "/api/health")
        assert health["available_run_keys"] == ["pipeline"]
        assert health["run_key_found"] is True


class TestAMismatchIsLoud:
    def test_a_wrong_run_key_is_reported_not_rendered(self, tmp_path):
        """The durable half of the fix.

        Reading the overlay fixes today's default but cannot catch a `--run`
        typo. A wrong key must be stated, not drawn as sixteen PENDING nodes.
        """
        db = _store_with_run(tmp_path, run_key="pipeline")
        app = create_app(str(db), "typo-key")

        health = _get(app, "/api/health")
        assert health["run_key_found"] is False
        assert health["available_run_keys"] == ["pipeline"]

        body = _snapshot(app)
        assert body["run_key_mismatch"] == {
            "requested": "typo-key",
            "available": ["pipeline"],
            "mismatch": True,
        }

    def test_an_empty_database_is_not_a_mismatch(self, tmp_path):
        """No runs yet is a normal state, not a wrong key.

        If this reported a mismatch, every fresh machine would show a warning
        for a pipeline that simply has not run - the same noise-to-signal
        problem in the other direction.
        """
        db = tmp_path / "empty.db"
        ExecutionStore(db).close()
        app = create_app(str(db), "pipeline")
        assert _get(app, "/api/health")["run_key_found"] is False
        assert _snapshot(app)["run_key_mismatch"]["mismatch"] is False

    def test_a_missing_database_is_not_a_mismatch(self, tmp_path):
        app = create_app(str(tmp_path / "absent.db"), "pipeline")
        assert _get(app, "/api/health")["database_present"] is False
        assert _snapshot(app)["run_key_mismatch"]["mismatch"] is False

    def test_a_secondary_key_still_finds_its_own_run(self, tmp_path):
        """`run_keys()` reports every key, so multiple runs coexist."""
        db = _store_with_run(tmp_path, run_key="pipeline")
        store = ExecutionStore(db)
        store.record(
            run_key="other", stage="mlst", state=StageState.SUCCEEDED, attempt=1
        )
        store.close()

        app = create_app(str(db), "other")
        assert _get(app, "/api/health")["available_run_keys"] == ["other", "pipeline"]
        assert _get(app, "/api/health")["run_key_found"] is True


class TestStoreRunKeys:
    def test_it_lists_distinct_keys_sorted(self, tmp_path):
        store = ExecutionStore(tmp_path / "k.db")
        for key in ("zulu", "alpha"):
            store.record(run_key=key, stage="amr", state=StageState.SUCCEEDED, attempt=1)
        store.record(run_key="alpha", stage="mlst", state=StageState.SUCCEEDED, attempt=1)
        assert store.run_keys() == ["alpha", "zulu"]
        store.close()

    def test_it_is_empty_on_a_fresh_store(self, tmp_path):
        store = ExecutionStore(tmp_path / "fresh.db")
        assert store.run_keys() == []
        store.close()


# --- helpers -------------------------------------------------------------
# The app is a plain ASGI app here rather than a live server: this is about
# which run key gets asked for, and serving HTTP would add nothing.


def _client(app):
    from fastapi.testclient import TestClient

    return TestClient(app)


def _get(app, path):
    return _client(app).get(path).json()


def _snapshot(app):
    return _client(app).get("/api/snapshot").json()