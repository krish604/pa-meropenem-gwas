"""Deliver ``scripts/pyseer_compat.py``'s shims to *spawned* pyseer workers.

This file is a ``sitecustomize`` module: any interpreter started with this
directory on ``PYTHONPATH`` imports it during startup, before any user code.

Why it exists, stated exactly:

* ``papipeline/gwas_real/adapter.py`` runs pyseer through
  ``scripts/pyseer_compat.py``, which applies its shims (``smf.Logit`` ->
  ``sm.Logit``, ``scipy.arange`` -> ``numpy.arange``, ...) inside ``main()``.
* With ``--cpu > 1``, pyseer fits every variant in a ``multiprocessing.Pool``
  worker (``pyseer/__main__.py``: ``if options.cpu > 1: pool = Pool(...)``).
* On macOS the default start method is ``spawn`` (``python3 -c "import
  multiprocessing as mp; print(mp.get_start_method())"`` -> ``spawn``), so a
  worker does not inherit the parent's memory. It re-runs the main script's
  *top level* under the name ``__mp_main__`` - which does **not** call
  ``main()``, because that is behind ``if __name__ == "__main__"`` - and then
  imports ``pyseer.model`` fresh.
* ``pyseer/model.py`` calls ``smf.Logit`` inside ``fit_lineage_effect``, which
  runs in a worker when ``--lineage`` is given. Without the alias the worker
  dies with ``AttributeError: module 'statsmodels.formula.api' has no
  attribute 'Logit'`` while the parent, which *did* run ``main()``, fits the
  null model happily - the failure looks like a broken statsmodels rather than
  like a shim that never reached the process that needed it.

The logic itself is not duplicated here: this module loads
``scripts/pyseer_compat.py`` and calls its own ``apply_shims()``, so the two
can never drift. The path is derived from this file's own location, never
from the environment.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

#: ``scripts/pyseer_compat.py``, three directories up from this file:
#: ``papipeline/gwas_real/_pyseer_shim/sitecustomize.py`` -> repository root.
SHIM_PATH = Path(__file__).resolve().parents[3] / "scripts" / "pyseer_compat.py"


def apply_repository_shim() -> None:
    """Run the repository's own ``apply_shims()``, if the file is there.

    A missing file is not an error worth raising: an interpreter started from
    a partial checkout should still start, and the parent process applies the
    shim again from ``main()`` regardless.
    """
    if not SHIM_PATH.is_file():
        return
    spec = importlib.util.spec_from_file_location(
        "_papipeline_pyseer_compat", SHIM_PATH
    )
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        return
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.apply_shims()


try:
    apply_repository_shim()
except Exception as exc:  # pragma: no cover - a hook must never stop startup
    print(
        "gwas_real: could not apply scripts/pyseer_compat.py in this "
        f"interpreter: {exc}",
        file=sys.stderr,
    )
