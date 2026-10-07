#!/usr/bin/env python3
"""Plan the batch Bakta annotation that ``scripts/linux/run_bakta.sh`` runs.

Two modes, both read-only. Nothing here executes Bakta:

``layout``
    Print where annotation goes, which database it uses and how many threads
    each process gets - three lines, one value each. Every value comes from the
    pipeline's own configuration loader, so this script and stage 2 agree about
    ``results/<mode>/intermediate/bakta/<sample_id>`` by construction rather
    than by two copies of the same path arithmetic.

``plan``
    Print one tab-separated row per sample in the manifest::

        sample_id <TAB> action <TAB> out_dir <TAB> command-or-reason

    with ``action`` one of ``run``, ``reuse``, ``missing`` or ``mismatch``.

The cohort, the assembly lookup, the file stem, the reuse verdict and the Bakta
command line are all taken from :mod:`papipeline` - ``manifest``,
``assemblies``, ``adapters.bakta`` and ``stages.annotation`` - because a second
implementation of any of them is how a batch annotation ends up somewhere stage
2 never looks.

``mismatch`` exists for one trap: Bakta names its outputs after the assembly
FILE, while stage 2 looks for ``<genome_stem_for(sample)>.gff3``, which is the
manifest's ``assembly_path`` stem or, with none, the ``sample_id``. Those agree
only when the manifest says where the file is. Annotating anyway would write
files the stage then refuses to reuse, so the sample is refused here, where the
message can say what to change.

Environment:
    ``BAKTA_DB``     overrides ``annotation.bakta_db`` (this pipeline's own rule)
    ``BAKTA_BIN``    overrides the Bakta executable (``adapters.bakta``)
    ``BAKTA_OUT_ROOT``  overrides the output root (``run_bakta.sh --out-root``)
    ``BAKTA_THREADS``   overrides the thread count (``run_bakta.sh --threads``)
"""

from __future__ import annotations

import argparse
import os
import shlex
import sys
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # so `python3 bakta_plan.py` works anywhere
    sys.path.insert(0, str(REPO_ROOT))

from papipeline.adapters.bakta import (  # noqa: E402
    bakta_command,
    genome_stem_for,
    resolve_executable,
)
from papipeline.assemblies import locate_assembly  # noqa: E402
from papipeline.config.loader import PipelineConfig, load_config  # noqa: E402
from papipeline.errors import PipelineError  # noqa: E402
from papipeline.io.tsv import read_tsv  # noqa: E402
from papipeline.manifest import (  # noqa: E402
    METADATA_FILENAME,
    SampleManifest,
    discover_manifest,
    manifest_from_rows,
)
from papipeline.models import RunMode  # noqa: E402
from papipeline.stages.annotation import _bakta_database, decide_reuse  # noqa: E402

#: Where stage 2 looks for per-genome Bakta output. Derived from configuration,
#: never written as a literal: see ``PipelineConfig.intermediate_root``.
BAKTA_SUBDIR = "bakta"


def load_manifest(path: Path, cfg: PipelineConfig) -> SampleManifest:
    """The cohort, through the same code the pipeline itself uses.

    Accepts the directory holding ``sample_metadata.tsv``, that file itself, or
    any other TSV with a ``sample_id`` column - a named file goes through the
    same reader and the same validation rather than a bespoke parser.
    """
    pattern = cfg.sample_id_pattern
    unique = cfg.sample_id_unique_attributes
    if path.is_dir():
        return discover_manifest(path, pattern=pattern, unique_attributes=unique)
    if path.name == METADATA_FILENAME:
        return discover_manifest(path.parent, pattern=pattern, unique_attributes=unique)
    rows = read_tsv(path, required_columns=("sample_id",), unique_columns=("sample_id",))
    return manifest_from_rows(rows, source=str(path), pattern=pattern, unique_attributes=unique)


def annotation_root(cfg: PipelineConfig, out_root: str | None = None) -> Path:
    if out_root:
        return Path(out_root)
    return cfg.intermediate_root(RunMode.REAL) / BAKTA_SUBDIR


def thread_count(cfg: PipelineConfig, threads: str | None = None) -> int:
    value = threads or os.environ.get("BAKTA_THREADS") or cfg.threads
    if not value or int(value) < 1:
        raise SystemExit(
            "no thread count: set runtime.threads in the machine overlay or "
            "pass --threads to run_bakta.sh"
        )
    return int(value)


def emit(sample_id: str, action: str, out_dir: Path, payload: str) -> None:
    """One plan row. Tabs and newlines cannot appear in any of the parts."""
    parts = (sample_id, action, str(out_dir), payload)
    for part in parts:
        if "\t" in part or "\n" in part:
            raise SystemExit(f"internal: unprintable character in plan row {part!r}")
    sys.stdout.write("\t".join(parts) + "\n")


def plan_rows(
    cfg: PipelineConfig,
    manifest: SampleManifest,
    out_root: Path,
    database: Path,
    threads: int,
    executable: str,
) -> Iterable[tuple[str, str, Path, str]]:
    """One row per sample, in manifest order."""
    for sample in manifest:
        out_dir = out_root / sample.sample_id
        try:
            assembly = locate_assembly(sample, cfg.assembly_root(RunMode.REAL))
        except PipelineError as exc:
            yield (sample.sample_id, "missing", out_dir, str(exc))
            continue

        stem = genome_stem_for(sample)
        written_by_bakta = Path(assembly).stem
        if stem != written_by_bakta:
            yield (
                sample.sample_id,
                "mismatch",
                out_dir,
                f"stage 2 expects {stem}.gff3 but Bakta names its files after "
                f"the assembly ({written_by_bakta}); add assembly_path="
                f"{assembly.name} to the manifest so the two agree",
            )
            continue

        decision = decide_reuse(sample, out_root=out_root, database_dir=database)
        if decision.reused:
            yield (sample.sample_id, "reuse", out_dir, decision.detail or decision.reason)
            continue

        command = bakta_command(
            executable, database, assembly, out_dir, stem, threads
        )
        yield (
            sample.sample_id,
            "run",
            out_dir,
            " ".join(shlex.quote(part) for part in command),
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=("layout", "plan"))
    parser.add_argument("--science", required=True)
    parser.add_argument("--machine", required=True)
    parser.add_argument("--manifest", help="sample_metadata.tsv, or its directory")
    parser.add_argument("--out-root")
    parser.add_argument("--threads")
    args = parser.parse_args(argv)

    cfg = load_config(Path(args.science), machine=args.machine)
    out_root = annotation_root(cfg, args.out_root or os.environ.get("BAKTA_OUT_ROOT"))
    database = _bakta_database(cfg)
    threads = thread_count(cfg, args.threads)

    if args.mode == "layout":
        # One value per line: a path with a space needs no quoting scheme.
        print(out_root)
        print(database)
        print(threads)
        return 0

    if not args.manifest:
        parser.error("--manifest is required for plan mode")
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        raise SystemExit(f"manifest not found: {manifest_path}")

    manifest = load_manifest(manifest_path, cfg)
    executable = resolve_executable()

    reuse_mode = cfg.reuse_tool_output()
    if reuse_mode == "off":
        # Loud, because it decides whether this batch is work or waste: with
        # reuse off, stage 2 invokes Bakta for every genome regardless of what
        # is on disk, so a pre-annotation would simply be repeated.
        print(
            "NOTE: annotation.reuse_tool_output is \"off\", so stage 2 will run "
            "Bakta for every genome and this batch's output will be re-run "
            "rather than reused. Set it to \"prefer\" or \"require\" in "
            f"{args.science} to reuse what is annotated here.",
            file=sys.stderr,
        )

    for row in plan_rows(cfg, manifest, out_root, database, threads, executable):
        emit(*row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
