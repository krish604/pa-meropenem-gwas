"""Q5 — the pure-JS unit tests, run with `node --test` (no npm install).

The JS tests live in ``dashboard/tests/js/`` and use a minimal DOM shim. This
wrapper runs them from the repository root so the relative imports resolve, and
skips with the exact reason when ``node`` is not on PATH.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from conftest import REPO_ROOT

JS_DIR = REPO_ROOT / "dashboard" / "tests" / "js"


@pytest.mark.skipif(shutil.which("node") is None, reason="node is not on PATH; JS tests skipped")
def test_js_unit_tests_pass():
    result = subprocess.run(
        ["node", "--test", str(JS_DIR)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, (
        "node --test failed\n--- stdout ---\n" + result.stdout[-8000:] + "\n--- stderr ---\n" + result.stderr[-4000:]
    )
    assert "pass " in result.stdout
