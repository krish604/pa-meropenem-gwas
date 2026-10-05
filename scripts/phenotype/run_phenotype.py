#!/usr/bin/env python3
"""Stage 11: imipenem phenotype interface.

Thin CLI wrapper. The logic lives in
``papipeline.stages.phenotype``; this script only resolves paths, loads the
configuration and the manifest, and prints a summary.

Stage 11 of 16. Scientific constraints for this stage are documented in
``docs/scientific_rules.md`` and enforced in the stage module, not here.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.config.loader import load_config                      # noqa: E402
from papipeline.errors import PipelineError                           # noqa: E402
from papipeline.logging_utils import configure_logging                # noqa: E402
from papipeline.models import RunMode                                 # noqa: E402
from papipeline.cli import summarise                            # noqa: E402
from papipeline.run import discover_run_manifest                       # noqa: E402
from papipeline.stages import phenotype as stage                          # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Stage 11: imipenem phenotype interface.")
    parser.add_argument("--config", type=Path, default=PIPELINE_ROOT / "config" / "science.yaml")
    parser.add_argument("--machine", default="laptop",
                        help="Machine overlay: laptop or bigmachine.")
    parser.add_argument("--mode", default="TEST", choices=["TEST", "REAL"])
    parser.add_argument("--antibiotic", default=None, help="Defaults to the configured primary antibiotic.")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging("INFO")

    try:
        config = load_config(args.config, machine=args.machine)
        mode = RunMode(args.mode)
        if mode is RunMode.REAL and not config.runtime.get("allow_real_mode", False):
            print(
                "ERROR: REAL mode is disabled (runtime.allow_real_mode is false).",
                file=sys.stderr,
            )
            return 2

        antibiotic = args.antibiotic or str(
            (config.raw.get("project") or {}).get("primary_antibiotic")
            or config.antibiotics[0]
        )
        config.require_antibiotic(antibiotic)

        data_root = config.data_root(mode)
        # `discover_run_manifest`, not `discover_manifest`. The smoke overlay
        # builds its cohort from `PDC_essential.tsv` plus the prepared assemblies
        # in `paths.smoke_genome_dir`; going straight to `discover_manifest` asks
        # for `data/<mode>/metadata/sample_metadata.tsv`, which a smoke overlay
        # never creates. The pair that holds for stage 11 is the cohort and the
        # phenotype table, so the cohort has to be discovered the same way the
        # run that consumes this stage discovers it.
        manifest = discover_run_manifest(config, mode)

        result = stage.load_phenotype(
            config,
            config.phenotype_dir(mode),
            antibiotic,
            filename=None,
            sample_ids=manifest.sample_ids,
        )
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    _summarise(result)
    return 0


def _summarise(result) -> None:
    """Print a short, factual summary of what the stage produced."""
    for line in summarise("phenotype", result):
        print(line)


if __name__ == "__main__":
    raise SystemExit(main())
