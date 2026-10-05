"""Shared fixtures for the dashboard's own test suite (QA-owned).

The suite reads the committed 10-isolate fixtures under
``dashboard/fixtures/small/`` and never touches the network or the real results
tree. Every app is built with a ``state_dir`` under ``tmp_path``, so UI-D1's
"the state directory is the only write" is testable without writing into the
developer's home directory.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

FIXTURES_DIR = REPO_ROOT / "dashboard" / "fixtures"
SMALL_LIVE = FIXTURES_DIR / "small" / "live"
SMALL_BUNDLE = FIXTURES_DIR / "small" / "bundle"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def small_live() -> Path:
    return SMALL_LIVE


@pytest.fixture(scope="session")
def small_bundle() -> Path:
    return SMALL_BUNDLE


def _build(results_root=None, bundle=None, state_dir=None, **kwargs):
    from dashboard.server.app import build_state, create_app

    state = build_state(
        results_root=Path(results_root) if results_root else None,
        bundle=Path(bundle) if bundle else None,
        repo_root=REPO_ROOT,
        state_dir=Path(state_dir) if state_dir else None,
        **kwargs,
    )
    return create_app(state), state


@pytest.fixture
def app_factory(tmp_path):
    """`make(results_root=..., bundle=..., ...) -> (TestClient, AppState)`."""
    from fastapi.testclient import TestClient

    created: List[str] = []

    def make(*, results_root=None, bundle=None, state_dir=None, **kwargs):
        directory = state_dir or (tmp_path / f"state_{len(created)}")
        created.append(str(directory))
        app, state = _build(
            results_root=results_root,
            bundle=bundle,
            state_dir=directory,
            **kwargs,
        )
        return TestClient(app), state

    return make


@pytest.fixture
def live_client(app_factory):
    return app_factory(results_root=SMALL_LIVE)


@pytest.fixture
def bundle_client(app_factory):
    return app_factory(bundle=SMALL_BUNDLE)


def snapshot_tree(root: Path) -> Dict[str, Tuple[int, int, str]]:
    """`{relative_path: (size, mtime_ns, sha256)}` for every file under root."""
    root = Path(root)
    out: Dict[str, Tuple[int, int, str]] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        stat = path.stat()
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        out[str(path.relative_to(root))] = (stat.st_size, stat.st_mtime_ns, digest)
    return out


def copy_layout(source: Path, destination: Path) -> Path:
    """Copy a committed fixture into ``tmp_path`` so a test may mutate it."""
    destination = Path(destination)
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(source, destination)
    return destination


def edge_layout(tmp_path: Path, name: str = "edge") -> Path:
    """A live copy with the six badges made distinct (DESIGN §4.2).

    `run_manifest.json` is rewritten and two tables are changed:

    * `mlst` stays `completed` but its table is removed -> `not_assessed`;
    * `amr` stays `completed` but its table is header-only -> `not produced`;
    * `convergence` is `completed` under `run_mode: REAL` -> `refused`
      (it is in `REAL_REFUSING_STAGES`);
    * `annotation` carries `skipped_*` and a `stages_skipped` cause ->
      `refused`;
    * `gwas` is `failed` -> `failed`;
    * `pangenome` is absent from `stages` -> `not_run`.
    """
    root = copy_layout(SMALL_LIVE, tmp_path / name)
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["run_mode"] = "REAL"
    manifest["stages"]["mlst"] = "completed"
    manifest["stages"]["amr"] = "completed"
    manifest["stages"]["convergence"] = "completed"
    manifest["stages"]["annotation"] = "skipped_no_phenotype"
    manifest["stages"]["gwas"] = "failed"
    manifest["stages"].pop("pangenome", None)
    manifest["stages_skipped"] = {
        "annotation": "disabled in configuration (analysis.annotation is false)"
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stages = root / "intermediate" / "stages"
    (stages / "03_mlst.tsv").unlink()
    (stages / "04_amr.tsv").write_text(
        "# SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS\nsample_id\tantibiotic\n",
        encoding="utf-8",
    )
    return root
