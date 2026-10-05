#!/usr/bin/env python3
"""Run a single stage (plus its prerequisites).

Called by the per-stage Snakemake rules. The orchestrator expands the
requested stage to include everything it depends on, so this is safe to call
for a stage whose inputs are not yet in memory.

Examples::

    python3 scripts/common/run_stage.py --mode TEST --stage amr
    python3 scripts/common/run_stage.py --mode TEST --stage gwas
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.config.loader import load_config                     # noqa: E402
from papipeline.errors import PipelineError                          # noqa: E402
from papipeline.logging_utils import configure_logging               # noqa: E402
from papipeline.run import PREREQUISITES, STAGE_ORDER, resolve_mode, run_pipeline  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run one pipeline stage and its prerequisites.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Prerequisites:\n"
        + "\n".join(f"  {k}: {', '.join(sorted(v)) or '(none)'}" for k, v in PREREQUISITES.items()),
    )
    parser.add_argument("--config", type=Path, default=PIPELINE_ROOT / "config" / "science.yaml")
    parser.add_argument("--machine", default="laptop",
                        help="Machine overlay: laptop or bigmachine.")
    parser.add_argument("--mode", default="TEST")
    parser.add_argument(
        "--stage",
        required=True,
        choices=sorted(STAGE_ORDER),
        help="The stage to run.",
    )
    parser.add_argument("--log-file", type=Path, default=None)
    parser.add_argument("--no-html", action="store_true")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging("INFO", logfile=args.log_file)

    try:
        config = load_config(args.config, machine=args.machine)
        resolve_mode(args.mode, config)  # fail fast on a forbidden mode
        result = run_pipeline(
            config_path=args.config,
            mode=args.mode,
            only=[args.stage],
            write_html=not args.no_html,
            # `config`, not `machine=`. This script loads the overlay in order to
            # validate the mode against it, then used to drop it on the floor:
            # `run_pipeline` re-loads with its own default of "laptop", so
            # `--machine <path/to/smoke.yaml>` was silently ignored. Stage 11 in
            # a bounded REAL smoke run failed for want of
            # `data/metadata/sample_metadata.tsv` - a file the smoke overlay
            # never uses - while the overlay that would have supplied the
            # cohort sat unread. Passing the object already validated keeps the
            # gate above and the run below reading the same overlay.
            config=config,
        )
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    state = result.stage_status.get(args.stage, "not_run")
    print(f"stage={args.stage} status={state}")
    if state not in ("completed", "skipped_no_phenotype"):
        print(f"WARNING: requested stage did not complete cleanly ({state})", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
