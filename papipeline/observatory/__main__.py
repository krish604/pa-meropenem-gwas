"""Entry point: ``python -m papipeline.observatory``."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional, Sequence

from .api import create_app


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m papipeline.observatory",
        description="Serve the PAPipeline Computational Observatory.")
    parser.add_argument("--db", default=None,
                        help="execution store to read (default: the overlay's "
                             "observatory.db, else ~/.local/share/papipeline/observatory.db)")
    parser.add_argument("--run", default=None,
                        help="run key to display (default: the overlay's "
                             "observatory.run_key, else 'observatory')")
    parser.add_argument("--machine", default="laptop",
                        help="machine overlay whose observatory settings to read: "
                             "laptop, bigmachine, or a path to a machines/*.yaml")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--demo", action="store_true",
                        help="build a real fixture cohort first, then serve it")
    parser.add_argument("--subjects", type=int, default=8)
    args = parser.parse_args(argv)

    bus = None
    db = args.db
    run = args.run

    # Read the overlay rather than the module constants. `DEFAULT_RUN_KEY` is
    # "observatory", but every machine overlay configures `run_key: pipeline`,
    # so a default-started server asked for a key no run has ever used and
    # rendered sixteen PENDING stages against a database full of real rows.
    # That is a wrong answer with no error in it, which is worse than a failure.
    # The overlay is the config pattern already used for the *writer* side
    # (`ObservatorySettings.from_config`), so the reader now uses it too, and
    # --db/--run stay as overrides rather than being the only source of truth.
    if db is None or run is None:
        from ..config import load_config
        from .wiring import ObservatorySettings

        config = load_config(
            Path(__file__).resolve().parents[2] / "config" / "science.yaml",
            machine=args.machine,
        )
        settings = ObservatorySettings.from_config(
            config, db_path=db, run_key=run
        )
        db = db if db is not None else str(settings.db_path)
        run = run if run is not None else settings.run_key

    if args.demo:
        from .fixture import Fixture

        fixture = Fixture(subjects=args.subjects, run_key=run)
        fixture.stage_cohort_baseline()
        db, run, bus = str(fixture.db), fixture.run_key, fixture.bus
        print(f"demo fixture: {db} (run {run})")
        print("building a real cohort; open the UI to watch it execute")
        fixture.run_background()

    import uvicorn

    print(f"\nobservatory -> http://{args.host}:{args.port}")
    print(f"  store: {db}")
    print(f"  run:   {run}")
    uvicorn.run(create_app(db, run, bus=bus), host=args.host, port=args.port,
                log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
