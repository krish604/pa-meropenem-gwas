"""Run the 10-genome Bakta concurrency benchmark.

    python3 -m papipeline.benchmark --out results/benchmarks/bakta10

The benchmark refuses to produce numbers it cannot stand behind. If Bakta
is unavailable, or the machine cannot safely run a configuration, that
configuration is recorded as NOT_RUN with the reason, and no figure is
invented to fill the gap.
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import os
import sys
import time
from pathlib import Path
from typing import List, Optional, Sequence

from ..config.loader import load_config
from ..observatory.events import EventBus
from .bakta10 import (
    NOT_RUN,
    N_GENOMES,
    WORKER_LEVELS,
    ConfigResult,
    audit_cohort,
    preflight,
    read_cohort_manifest,
    run_configuration,
    select_cohort,
    write_cohort_manifest,
    write_environment,
    write_per_genome,
    write_summary,
)
from .report import write_report

REPO = Path(__file__).resolve().parents[2]


def config_with_threads(config, threads: int, database: Path):
    """A config that permits REAL mode and names the pinned database."""
    runtime = dict(config.runtime)
    runtime["allow_real_mode"] = True
    runtime["threads"] = int(threads)
    raw = copy.deepcopy(dict(config.raw or {}))
    raw.setdefault("annotation", {})
    raw["annotation"]["bakta_db"] = str(database)
    return dataclasses.replace(config, runtime=runtime, raw=raw)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python3 -m papipeline.benchmark",
        description="Controlled Bakta concurrency benchmark over 10 genomes.")
    parser.add_argument("--out", type=Path,
                        default=REPO / "results" / "benchmarks" / "bakta10",
                        help="output directory (default: %(default)s)")
    parser.add_argument("--data", type=Path, default=REPO / "data",
                        help="directory holding the assemblies")
    parser.add_argument("--genomes", type=int, default=N_GENOMES)
    parser.add_argument("--workers", type=int, nargs="+", default=list(WORKER_LEVELS),
                        help="concurrency levels to measure")
    parser.add_argument("--threads", type=int, default=None,
                        help="Bakta --threads. Defaults to the configured "
                             "value and is NOT varied by this sweep.")
    parser.add_argument("--database", type=Path, default=None,
                        help="pinned Bakta database directory. Falls back to "
                             "$BAKTA_DB, then annotation.bakta_db in "
                             "the machine overlay. The database is referenced where "
                             "it lives; it is never copied.")
    parser.add_argument("--bakta", type=Path, default=None,
                        help="path to the Bakta executable. Falls back to "
                             "$BAKTA_BIN, then PATH. Bakta's helper tools "
                             "are resolved from the executable's own "
                             "directory.")
    parser.add_argument("--preflight-only", action="store_true",
                        help="run the checks and report, then stop")
    parser.add_argument("--verify-one", action="store_true",
                        help="annotate only the first genome of the frozen "
                             "cohort, through the normal PAPipeline path, as "
                             "an end-to-end plumbing check. This is NOT a "
                             "benchmark configuration and is never written to "
                             "summary.tsv.")
    args = parser.parse_args(argv)

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    # ── the cohort, fixed once and written down ──────────────────
    cohort = select_cohort(args.data, args.genomes)
    cohort_path = write_cohort_manifest(out_root / "cohort.tsv", cohort)
    audit = audit_cohort(cohort, REPO / "PDC_essential.tsv")

    print("BENCHMARK COHORT")
    for genome in cohort:
        print(f"{genome.rank:2d}. {genome.accession}  {genome.path.name}")
    print(f"\ncohort manifest: {cohort_path}")
    print(f"overlap with the first {args.genomes} rows of PDC_essential.tsv: "
          f"{audit['overlap_with_pdc_first_n']} of {args.genomes} "
          f"(identical: {audit['identical_to_pdc_first_n']})")
    print("cohort is selected from the assemblies that exist under "
          f"{args.data}, in the repository's own GCA accession order.\n")

    base = load_config(REPO / "config" / "science.yaml", machine=args.machine)
    threads = args.threads if args.threads is not None else int(
        base.runtime.get("threads", 1))
    import os

    configured_db = (dict(base.raw or {}).get("annotation") or {}).get("bakta_db")
    database = (
        args.database
        or (Path(os.environ["BAKTA_DB"]) if os.environ.get("BAKTA_DB") else None)
        or (Path(configured_db) if configured_db else None)
    )
    if database is None:
        # Only a last resort. The default location is almost never right:
        # a real database is usually provisioned with its environment.
        database = REPO / ".bakta" / "db"

    results: List[ConfigResult] = []
    checks = {
        w: preflight(database, threads, w, bakta=args.bakta)
        for w in args.workers
    }

    print("PREFLIGHT")
    for workers, check in checks.items():
        verdict = "OK" if check.ok else NOT_RUN
        print(f"  workers={workers}: {verdict}")
        for reason in check.reasons:
            print(f"      - {reason}")
    write_environment(
        out_root / "environment.tsv", checks[max(checks)],
        extra={"cohort_manifest": str(cohort_path),
               "cohort_overlap_with_pdc_first_n":
                   audit["overlap_with_pdc_first_n"],
               "cohort_identical_to_pdc_first_n":
                   audit["identical_to_pdc_first_n"]},
    )
    if args.preflight_only:
        print("\n--preflight-only: stopping before any execution.")
        return 0 if all(c.ok for c in checks.values()) else 2

    if args.verify_one:
        # A single genome through the real adapter, contract and store. Its
        # purpose is to prove the recovered Bakta setup works *through
        # PAPipeline*, not to measure anything.
        genome = cohort[0]
        directory = out_root / "verify_one"
        if directory.exists():
            import shutil
            shutil.rmtree(directory)
        directory.mkdir(parents=True, exist_ok=True)
        print(f"\nverify-one: annotating {genome.accession} through PAPipeline "
              f"(workers=1, --threads {threads})")
        print("  this is a plumbing check, not a benchmark result\n")
        result = run_configuration(
            workers=1, cohort=[genome], out_dir=directory,
            database=database, threads=threads,
            config=config_with_threads(base, threads, database),
            run_key="bakta10-verify", bus=EventBus(),
            genome_dir=args.data,
        )
        row = result.per_genome[0] if result.per_genome else None
        print(f"  exit/state : {row.state if row else 'no result'}")
        print(f"  attempts   : {row.attempts if row else 'NA'}")
        print(f"  wall       : {row.wall_seconds if row and row.wall_seconds else 'NA'} s")
        print(f"  peak RSS   : {row.peak_rss_mb if row and row.peak_rss_mb else 'NA'} MB")
        print(f"  PID        : {row.pid if row else 'NA'}")
        print(f"  validation : {row.detail[:150] if row else 'NA'}")
        write_per_genome(directory / "per_genome.tsv", [result])
        print(f"\n  outputs: {directory / 'bakta'}")
        print(f"  per-genome row: {directory / 'per_genome.tsv'}")
        return 0 if result.completed == 1 else 2

    bus = EventBus()
    for workers in args.workers:
        directory = out_root / f"workers_{workers}"
        check = checks[workers]
        if not check.ok:
            print(f"\nworkers={workers}: {NOT_RUN}")
            for reason in check.reasons:
                print(f"  reason: {reason}")
            results.append(ConfigResult(
                workers=workers, bakta_threads=threads, genomes=args.genomes,
                status=NOT_RUN, reason="; ".join(check.reasons),
                run_key=f"bakta10-w{workers}",
            ))
            continue

        # A fresh directory per configuration: a result can never be
        # accepted from an earlier sweep.
        if directory.exists():
            import shutil
            shutil.rmtree(directory)
        directory.mkdir(parents=True, exist_ok=True)
        print(f"\nworkers={workers}: running {args.genomes} genomes "
              f"with --threads {threads} ...")
        started = time.time()
        result = run_configuration(
            workers=workers, cohort=cohort, out_dir=directory,
            database=database, threads=threads,
            config=config_with_threads(base, threads, database),
            run_key=f"bakta10-w{workers}", bus=bus,
            genome_dir=args.data,
        )
        print(f"  wall={result.wall_seconds:.1f}s  completed={result.completed}"
              f"  failed={result.failed}  invalid={result.invalid}"
              f"  incomplete={result.incomplete}")
        results.append(result)
        write_per_genome(directory / "per_genome.tsv", [result])

    write_summary(out_root / "summary.tsv", results)
    write_per_genome(out_root / "per_genome.tsv", results)
    (out_root / "cohort_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    write_report(out_root / "benchmark_report.md", cohort=cohort, results=results,
                 checks=checks, audit=audit, threads=threads,
                 database=database, data_dir=args.data)

    ran = [r for r in results if r.status == "OK"]
    print(f"\nsummary: {out_root / 'summary.tsv'}")
    print(f"report : {out_root / 'benchmark_report.md'}")
    print(f"configurations measured: {len(ran)} of {len(results)}")
    if len(ran) < len(results):
        print("NOT_RUN configurations carry their reason in summary.tsv; "
              "no figures were invented for them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
