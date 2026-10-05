#!/usr/bin/env python3
"""Generate the synthetic TEST fixtures.

The generator is seeded, so repeated runs produce byte-identical files. It
refuses to write into a directory named ``data``.

Examples::

    python3 scripts/common/make_synthetic_data.py
    python3 scripts/common/make_synthetic_data.py --n-samples 20 --seed 20240617
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.errors import PipelineError                          # noqa: E402
from papipeline.logging_utils import configure_logging               # noqa: E402
from papipeline.testing import (                                     # noqa: E402
    DEFAULT_N_SAMPLES,
    DEFAULT_SEED,
    SYNTHETIC_BANNER,
    generate,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate synthetic test fixtures. NOT biological data.",
        epilog=SYNTHETIC_BANNER,
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PIPELINE_ROOT / "test_data",
        help="Destination directory (default: %(default)s)",
    )
    parser.add_argument(
        "--n-samples", type=int, default=DEFAULT_N_SAMPLES,
        help="Number of synthetic samples (default: %(default)s)",
    )
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help="RNG seed. Same seed gives byte-identical output (default: %(default)s)",
    )
    parser.add_argument("--antibiotic", default="imipenem")
    parser.add_argument(
        "--no-genomes",
        action="store_true",
        help="Skip writing FASTA files (faster; stage 1 then reports missing input).",
    )
    parser.add_argument("--log-file", type=Path, default=None)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging("INFO", logfile=args.log_file)

    print(SYNTHETIC_BANNER)
    try:
        written = generate(
            args.out_dir,
            n_samples=args.n_samples,
            seed=args.seed,
            antibiotic=args.antibiotic,
            write_genome_files=not args.no_genomes,
            regulators_table=PIPELINE_ROOT / "config" / "regulators.tsv",
        )
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"\nWrote {len(written)} files under {args.out_dir}:")
    for key, path in sorted(written.items()):
        print(f"  {key:24s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
