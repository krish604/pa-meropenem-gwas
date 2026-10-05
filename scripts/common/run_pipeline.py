#!/usr/bin/env python3
"""Run the whole pipeline in one process.

This is the primary entry point. It needs only the pinned Python
dependencies and does not require Snakemake.

Examples::

    python3 scripts/common/run_pipeline.py --mode TEST
    python3 scripts/common/run_pipeline.py --mode TEST --no-html
    python3 scripts/common/run_pipeline.py --mode TEST --only gwas
    python3 scripts/common/run_pipeline.py --mode TEST --skip gwas cooccurrence
    python3 scripts/common/run_pipeline.py --mode TEST --dry-run
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
from papipeline.observatory import wiring                           # noqa: E402
from papipeline.run import (                                         # noqa: E402
    EXECUTION_ORDER,
    STAGE_ORDER,
    resolve_mode,
    run_pipeline,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the P. aeruginosa imipenem AMR pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Execution order (differs from stage numbering):\n  "
        + " -> ".join(EXECUTION_ORDER),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=PIPELINE_ROOT / "config" / "science.yaml",
        help="Path to the science config (default: %(default)s).",
    )
    parser.add_argument(
        "--machine",
        default="laptop",
        help="Machine overlay: a name (laptop, bigmachine) or a path to a "
             "config/machines/*.yaml file. Default: %(default)s",
    )
    parser.add_argument(
        "--mode",
        default="STUB",
        help="Run mode: STUB (fabricated outputs), TEST (synthetic fixtures) "
             "or REAL. Default: %(default)s",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="STAGE",
        help="Run only these stages plus their prerequisites.",
    )
    parser.add_argument(
        "--skip",
        nargs="+",
        metavar="STAGE",
        help="Skip these stages.",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Write only the Markdown report.",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Also write logs to this file.",
    )
    parser.add_argument(
        "--observatory",
        action="store_true",
        help="Record stage state, validate stage outputs and publish live "
             "events for the Computational Observatory. Off by default; the "
             "pipeline does not require it.",
    )
    parser.add_argument(
        "--observatory-db",
        type=Path,
        default=None,
        help="SQLite execution store for the observatory "
             "(default: ~/.local/share/papipeline/observatory.db, or "
             "runtime.observatory.db in the machine overlay).",
    )
    parser.add_argument(
        "--observatory-run",
        default=None,
        help="Run key the observatory should display (default: "
             "runtime.observatory.run_key, or 'pipeline').",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Continue this run key rather than re-executing it. Stages "
             "whose recorded outputs still re-validate are skipped; the "
             "rest are re-executed. Without this flag a second invocation "
             "re-runs everything, which is a deliberate re-run, not a resume.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate configuration and print the stage plan without running.",
    )
    parser.add_argument(
        "--list-stages",
        action="store_true",
        help="Print the stage list and exit.",
    )
    return parser


def print_plan(config_path: Path, mode: str, only, skip, machine: str = "laptop") -> int:
    # `machine` is threaded in rather than defaulted away: `load_config` already
    # defaults to the laptop overlay, so calling it bare made `--dry-run` report
    # the laptop's resolved configuration for *every* machine, including
    # `--machine bigmachine`. The plan is a statement about the run that would
    # follow it, so a plan that ignored the machine was worse than no plan - it
    # was a correct-looking plan for a different run.
    config = load_config(config_path, machine=machine)
    try:
        resolved = resolve_mode(mode, config)
    except PipelineError as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2

    selected = set(only) if only else set(STAGE_ORDER)
    selected -= set(skip or [])
    from papipeline.run import resolve_prerequisites

    selected |= resolve_prerequisites(selected) if only else set()

    print(f"config     : {config_path}")
    # Which overlay the plan below was resolved from, plus the two runtime
    # values it actually changes. Without these the plan is unverifiable from
    # its own output: the two overlays share every `paths` entry, so a plan
    # printed for the wrong machine was indistinguishable from a right one.
    print(f"machine    : {config.machine_name}")
    print(f"mode       : {resolved.value}")
    print(f"organism   : {config.organism.get('name')}")
    print(f"antibiotics: {', '.join(config.antibiotics)}")
    print(f"data root  : {config.data_root(resolved)}")
    print(f"results    : {config.results_root(resolved)}")
    print(f"sources    : {config.tool_output_root(resolved)}")
    print(f"threads    : {config.runtime.get('threads')}")
    print(f"memory_mb  : {config.runtime.get('memory_mb')}")
    print(f"real mode  : {'ENABLED' if config.runtime.get('allow_real_mode') else 'disabled'}")
    print("plan       :")
    for stage in EXECUTION_ORDER:
        mark = "run" if stage in selected else "skip"
        print(f"  [{mark}] {stage}")
    unpinned = config.unpinned_references()
    if unpinned:
        print(f"warning    : {len(unpinned)} unpinned reference(s) in config/references.tsv")
    return 0


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_stages:
        print("Stages (numbered as in the specification):")
        for index, stage in enumerate(STAGE_ORDER, start=1):
            print(f"  {index:2d}. {stage}")
        print("\nExecution order (dependency order):")
        for index, stage in enumerate(EXECUTION_ORDER, start=1):
            print(f"  {index:2d}. {stage}")
        return 0

    configure_logging("INFO", logfile=args.log_file)

    if args.dry_run:
        return print_plan(
            args.config, args.mode, args.only, args.skip, machine=args.machine
        )

    try:
        observatory = wiring.create(
            load_config(args.config, machine=args.machine),
            db_path=args.observatory_db,
            run_key=args.observatory_run,
            # An explicit --observatory-db or --observatory-run implies the
            # user wants it on, without needing both flags.
            enable=True if (args.observatory
                            or args.observatory_db
                            or args.observatory_run) else None,
        )
        result = run_pipeline(
            config_path=args.config,
            mode=args.mode,
            only=args.only,
            skip=args.skip,
            write_html=not args.no_html,
            observatory=observatory,
            resume=args.resume,
            machine=args.machine,
        )
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print()
    print(f"Run complete: mode={result.mode.value} samples={result.n_samples}")
    for stage, state in result.stage_status.items():
        print(f"  {stage:24s} {state}")
        
    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print(f"  - {warning}")
    for kind, path in result.report_paths.items():
        print(f"\n{kind}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
