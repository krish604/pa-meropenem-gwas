#!/usr/bin/env python3
"""Write the run manifest: tool, version, database and command provenance.

Also reports which external tools the pipeline expects versus which are
actually installed, so a missing dependency is visible before an analysis is
attempted rather than halfway through it.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline import __version__                                   # noqa: E402
from papipeline.adapters import detect_all                           # noqa: E402
from papipeline.config.loader import load_config                     # noqa: E402
from papipeline.errors import PipelineError                          # noqa: E402
from papipeline.logging_utils import configure_logging               # noqa: E402
from papipeline.run import build_provenance, resolve_mode, utc_now   # noqa: E402

#: Tools each stage needs. Reported as a matrix so a reviewer can see exactly
#: which stages are currently runnable.
STAGE_TOOL_REQUIREMENTS = {
    "annotation": ["bakta"],
    "mlst": ["mlst"],
    "amr": ["amrfinder"],
    "amr_optional": ["rgi"],
    "virulence": ["makeblastdb", "blastn"],
    "pangenome": ["panaroo"],
    "phylogeny": ["snp-sites", "iqtree2"],
    "gwas": ["pyseer"],
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write the provenance manifest for a run."
    )
    parser.add_argument("--config", type=Path, default=PIPELINE_ROOT / "config" / "science.yaml")
    parser.add_argument("--machine", default="laptop",
                        help="Machine overlay: laptop or bigmachine.")
    parser.add_argument("--mode", default="TEST")
    parser.add_argument("--output", type=Path, default=None, help="Defaults to <results>/run_manifest.json")
    parser.add_argument("--log-file", type=Path, default=None)
    parser.add_argument(
        "--print-tools", action="store_true", help="Only print the tool matrix."
    )
    return parser


def tool_matrix(tools) -> dict:
    """Availability of every tool the stages depend on."""
    matrix = {}
    for stage, names in STAGE_TOOL_REQUIREMENTS.items():
        matrix[stage] = {
            name: {
                "available": bool(tools.get(name) and tools[name].available),
                "version": tools[name].version if tools.get(name) else None,
            }
            for name in names
        }
    return matrix


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging("INFO", logfile=args.log_file)

    try:
        config = load_config(args.config, machine=args.machine)
        mode = resolve_mode(args.mode, config)
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    tools = detect_all()
    matrix = tool_matrix(tools)

    if args.print_tools:
        print(f"{'stage':16s} {'tool':16s} {'available':10s} version")
        for stage, entries in matrix.items():
            for name, info in entries.items():
                print(
                    f"{stage:16s} {name:16s} "
                    f"{'yes' if info['available'] else 'NO':10s} {info['version'] or '-'}"
                )
        return 0

    results_root = config.results_root(mode)
    results_root.mkdir(parents=True, exist_ok=True)
    output = args.output or (results_root / "run_manifest.json")

    provenance = build_provenance(config, tools)
    payload = {
        "pipeline_version": __version__,
        "run_mode": mode.value,
        "generated_at_utc": utc_now(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "references": provenance,
        "stage_tool_requirements": matrix,
        "unpinned_references": [spec.reference_id for spec in config.unpinned_references()],
        "note": (
            "Databases are never downloaded or updated during a run. Pin every "
            "version in config/references.tsv before a real analysis."
        ),
    }
    output.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"Provenance written to {output}")

    unpinned = payload["unpinned_references"]
    if unpinned:
        print(f"\nWARNING: {len(unpinned)} unpinned reference(s): {', '.join(unpinned)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
