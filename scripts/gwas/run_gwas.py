#!/usr/bin/env python3
"""Stage 12: imipenem GWAS. Requires the full upstream state; this wrapper runs the pipeline with --only gwas, which pulls in every prerequisite.

Associations are not causal claims. Features confined to a single lineage are reported as lineage-linked.

Runs the native orchestrator with ``--only gwas``; the orchestrator expands
that to the stage plus its prerequisites, so this wrapper is correct from a
cold start. It does not reimplement any pipeline logic.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.cli import summarise                            # noqa: E402
from papipeline.errors import PipelineError                     # noqa: E402
from papipeline.logging_utils import configure_logging          # noqa: E402
from papipeline.run import run_pipeline                         # noqa: E402

STAGE = "gwas"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Stage 12: imipenem GWAS. Requires the full upstream state; this wrapper runs the pipeline with --only gwas, which pulls in every prerequisite.")
    parser.add_argument("--config", type=Path, default=PIPELINE_ROOT / "config" / "science.yaml")
    parser.add_argument("--machine", default="laptop",
                        help="Machine overlay: laptop or bigmachine.")
    parser.add_argument("--mode", default="TEST", choices=["TEST", "REAL"])
    parser.add_argument("--no-html", action="store_true")
    parser.add_argument("--log-file", type=Path, default=None)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging("INFO", logfile=args.log_file)

    try:
        result = run_pipeline(
            config_path=args.config,
            mode=args.mode,
            only=[STAGE],
            write_html=not args.no_html,
        )
    except PipelineError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1

    print()
    for line in summarise("gwas", result.stage_status.get(STAGE, "not_run")):
        print(line)
    for line in summarise("samples", result.n_samples):
        print(line)
    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print("  - {}".format(warning))
    for kind, path in result.report_paths.items():
        print("\n{}: {}".format(kind, path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
