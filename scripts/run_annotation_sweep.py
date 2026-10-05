#!/usr/bin/env python3
"""Stage 2 background runner: Bakta annotation over the prepared cohort.

This does not reimplement any pipeline logic. It provisions nothing and
invents nothing; it invokes Bakta on the assemblies the existing prepare
stage already selected, writing each result where
``papipeline.stages.annotation`` expects to find it, so the existing parser
consumes the output unchanged.

Output layout, per sample, under ``<out_dir>/<Pilot_ID>/``:
  bakta.gff      parsed by stages.annotation.parse_bakta_gff
  bakta.tsv      parsed by stages.annotation.parse_bakta_tsv

Run modes:
  --prepare   create per-sample output directories and report the plan
  --run       annotate every sample (resumable: existing results are reused)

Resumability matters because a full sweep is hours long. A completed
sample is never re-annotated, so the sweep can be stopped and restarted
without losing work and without changing any result.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.errors import PipelineError                        # noqa: E402
from papipeline.io.tsv import read_tsv, write_tsv                  # noqa: E402
from papipeline.logging_utils import configure_logging, get_logger  # noqa: E402

LOGGER = get_logger("pilot.annotation")

OUT_ROOT = PIPELINE_ROOT / "results" / "pilot100" / "annotation"
MANIFEST = PIPELINE_ROOT / "results" / "pilot100" / "pilot100_manifest.tsv"
STATE = PIPELINE_ROOT / "results" / "pilot100" / "annotation_state.json"

#: Bakta's own completion marker. Present only on a finished run.
BAKTA_GFF = "bakta.gff"
BAKTA_TSV = "bakta.tsv"


def database_ready(database: Path) -> Tuple[bool, str]:
    """Whether the Bakta database is fully provisioned.

    A partially downloaded archive is not a usable database, so the check
    looks for the extracted database layout rather than "the directory is
    non-empty": an in-progress ``db-light.tar.xz`` must not pass.
    """
    if not database.exists():
        return False, "database directory does not exist"
    entries = {p.name for p in database.iterdir()}
    if entries and all(name.endswith((".tar.xz", ".tar.gz", ".part")) for name in entries):
        return False, "database archive still downloading (no extracted db/)"
    # The database must be the directory that actually holds version.json.
    if not (database / "version.json").exists():
        return False, (
            f"no version.json in {database.name} (found: {sorted(entries)[:4]}); "
            "--db must point at the extracted db-light/ directory"
        )
    return True, f"database present, version.json readable in {database.name}/"


def default_db() -> Path:
    """Locate the extracted Bakta database.

    ``bakta_db install`` writes ``<output>/db-light/`` (or ``db-full/``) and
    puts ``version.json`` inside it. Bakta's ``--db`` must point at *that*
    directory, not its parent, or it fails with "version file not readable".
    """
    root = Path.home() / "micromamba_pilot/envs/pilot100/share/bakta_db_pa"
    if not root.exists():
        return root
    # Prefer a directory that actually holds version.json.
    if (root / "version.json").exists():
        return root
    for candidate in sorted(root.iterdir()):
        if candidate.is_dir() and (candidate / "version.json").exists():
            return candidate
    # Otherwise fall back to the first extracted directory.
    dirs = [p for p in sorted(root.iterdir()) if p.is_dir()]
    return dirs[0] if dirs else root


def load_cohort(manifest: Path) -> List[Tuple[str, Path]]:
    """(Pilot_ID, genome path) for every analysable sample, in list order."""
    rows = read_tsv(manifest, required_columns=("Pilot_ID", "Assembly", "Genome_path"))
    cohort: List[Tuple[str, Path]] = []
    for row in rows:
        if str(row.get("analysis_included") or "TRUE").upper() != "TRUE":
            continue
        genome = Path(str(row["Genome_path"]))
        if not genome.exists():
            LOGGER.warning("Skipping %s: genome not found at %s", row["Pilot_ID"], genome)
            continue
        cohort.append((str(row["Pilot_ID"]), genome))
    return cohort


def is_complete(out_dir: Path) -> bool:
    """A sample counts as done only when Bakta's own outputs are present.

    Bakta names its outputs after the input file, not ``bakta.gff``, so
    completion is detected by looking for the annotation GFF plus the main
    annotation table regardless of the input's basename. Checking for
    literally ``bakta.gff`` would report every successful run as failed.

    The main annotation table is required rather than the inference table:
    the inference table can exist while the annotation table is absent, and
    only the annotation table carries gene names.
    """
    if not out_dir.exists():
        return False
    gff = bakta_gff(out_dir)
    table = bakta_annotation_tsv(out_dir)
    return gff is not None and table is not None


def bakta_gff(out_dir: Path) -> Optional[Path]:
    """The annotation GFF, whatever the input was named."""
    for pattern in ("*.gff3", "*.gff"):
        found = sorted(p for p in out_dir.glob(pattern) if p.stat().st_size > 0)
        if found:
            return found[0]
    return None


def bakta_annotation_tsv(out_dir: Path) -> Optional[Path]:
    """The main annotation table: the only one carrying gene names.

    Bakta writes three TSVs side by side and they are NOT interchangeable:

    * ``<input>.tsv``             - full annotation, has ``Gene``/``Product``
    * ``<input>.inference.tsv``   - expert-system hits, columns are
                                    Sequence Id/Type/Start/Stop/Strand/
                                    Locus Tag/Score/Evalue/Query Cov/
                                    Subject Cov/Id/Accession
    * ``<input>.hypotheticals.tsv`` - hypothetical proteins with hit scores

    Returning the inference table where the annotation table was expected
    loses every gene name and silently reduces gene-level analysis to locus
    tags, so the suffix is excluded explicitly.
    """
    candidates = [
        p
        for p in out_dir.glob("*.tsv")
        if p.stat().st_size > 0
        and not p.name.endswith(".inference.tsv")
        and not p.name.endswith(".hypotheticals.tsv")
    ]
    return sorted(candidates)[0] if candidates else None


def bakta_inference_tsv(out_dir: Path) -> Optional[Path]:
    """The expert-system inference table, if Bakta produced one."""
    return next(iter(sorted(out_dir.glob("*.inference.tsv"))), None)


def bakta_outputs(out_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """The (GFF, annotation TSV) pair the existing parser consumes."""
    return bakta_gff(out_dir), bakta_annotation_tsv(out_dir)


def bakta_executable() -> str:
    import shutil

    path = shutil.which("bakta")
    if path is None:
        raise PipelineError(
            "bakta is not on PATH",
            hint="micromamba run -n pilot100 ... (or activate the env)",
        )
    return path


def run_one(
    pilot_id: str,
    genome: Path,
    out_dir: Path,
    database: Path,
    threads: int,
    species: Optional[str],
    timeout: int,
) -> Dict[str, object]:
    """Annotate one genome. Returns a state record; never raises.

    The whole body is guarded. With 200 genomes in a queue a single
    unhandled exception would abandon every genome not yet started, so
    every failure mode - missing binary, unwritable output directory, a
    genome that is empty or absent, a timeout - is converted into a state
    record and the sweep moves on.
    """
    started = time.monotonic()
    if is_complete(out_dir):
        return {
            "pilot_id": pilot_id,
            "status": "reused",
            "seconds": 0.0,
            "out_dir": str(out_dir),
        }

    def failure(status: str, reason: str) -> Dict[str, object]:
        return {
            "pilot_id": pilot_id,
            "status": status,
            "seconds": round(time.monotonic() - started, 1),
            "error": reason,
        }

    # The genome is checked before the tool: a missing or empty input is the
    # actionable error, and reporting "bakta not found" instead would send
    # whoever is reading the log after the wrong problem.
    if not Path(genome).exists():
        return failure("error", f"genome not found: {genome}")
    try:
        if Path(genome).stat().st_size == 0:
            return failure("error", f"genome is empty: {genome}")
    except OSError as exc:
        return failure("error", f"genome unreadable: {exc}")

    try:
        executable = bakta_executable()
    except Exception as exc:  # noqa: BLE001 - must not abort the sweep
        return failure("error", f"bakta unavailable: {exc}")

    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return failure("error", f"cannot create output directory: {exc}")

    # Bakta takes the genome as a POSITIONAL argument; there is no --in
    # flag in 1.12.
    command = [
        executable,
        str(genome),
        "--out",
        str(out_dir),
        "--db",
        str(database),
        "--threads",
        str(threads),
        "--force",
    ]
    if species:
        command += ["--species", species]

    LOGGER.info("bakta %s <- %s", pilot_id, genome.name)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired:
        return failure("timeout", f"exceeded {timeout}s")
    except OSError as exc:
        return failure("error", str(exc))
    except Exception as exc:  # noqa: BLE001 - must not abort the sweep
        return failure("error", f"{type(exc).__name__}: {exc}")

    elapsed = round(time.monotonic() - started, 1)
    if completed.returncode != 0 or not is_complete(out_dir):
        return {
            "pilot_id": pilot_id,
            "status": "failed",
            "seconds": elapsed,
            "returncode": completed.returncode,
            "error": (completed.stderr or "")[-600:],
        }
    gff, table = bakta_outputs(out_dir)
    return {
        "pilot_id": pilot_id,
        "status": "completed",
        "seconds": elapsed,
        "out_dir": str(out_dir),
        "gff": str(gff) if gff else None,
        "inference_tsv": str(table) if table else None,
    }


def save_state(records: Dict[str, Dict[str, object]]) -> None:
    payload = {
        "updated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "samples": records,
    }
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(payload, indent=2), encoding="utf-8")


#: Bakta output the pipeline actually consumes. Everything else is dropped:
#: the annotation GFF is parsed by stage 2, the main TSV supplies gene names,
#: and the inference table is recorded in the state file for provenance.
KEEP_SUFFIXES = (".gff3", ".gff", ".tsv", ".log")


def prune_outputs(out_dir: Path, keep_suffixes=KEEP_SUFFIXES) -> Tuple[int, int]:
    """Delete Bakta output the pipeline never reads. Returns (removed, bytes).

    Bakta writes ~84 MB per genome: json 20 MB, embl 15 MB, gbff 14 MB, svg
    8 MB, fna 6 MB, ffn 6 MB, faa 2 MB, png 2 MB, log 2 MB, against a GFF of
    7 MB and a TSV of 0.7 MB. Roughly 90% of the volume is unused. At 200
    genomes that is the difference between ~17 GB and ~2 GB, which matters
    when free space is the binding constraint.

    The GFF, the main annotation TSV and the inference TSV are all kept, so
    pruning never invalidates :func:`is_complete` and never forces a
    re-annotation. The per-genome Bakta log is kept because it is the only
    record of what Bakta did.
    """
    if not out_dir.exists():
        return 0, 0
    removed = 0
    freed = 0
    for path in sorted(out_dir.iterdir()):
        if not path.is_file():
            continue
        if path.suffix in keep_suffixes:
            continue
        try:
            size = path.stat().st_size
            path.unlink()
        except OSError:
            continue
        removed += 1
        freed += size
    return removed, freed


def prune_all(out_dir: Path) -> Dict[str, int]:
    """Prune every sample directory under ``out_dir``."""
    total_removed = 0
    total_freed = 0
    for sample_dir in sorted(out_dir.iterdir()) if out_dir.exists() else []:
        if sample_dir.is_dir():
            removed, freed = prune_outputs(sample_dir)
            total_removed += removed
            total_freed += freed
    return {"files_removed": total_removed, "bytes_freed": total_freed}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Bakta sweep for the pilot cohort.")
    parser.add_argument("--prepare", action="store_true", help="Report the plan only.")
    parser.add_argument("--run", action="store_true", help="Annotate the cohort.")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=OUT_ROOT)
    parser.add_argument("--db", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--species", default="Pseudomonas aeruginosa")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--limit", type=int, default=None, help="Annotate only the first N pending.")
    parser.add_argument("--prune", action="store_true",
                        help="Delete unused Bakta output and exit. Reclaims ~90%% of the volume.")
    parser.add_argument("--no-prune", action="store_true",
                        help="Keep unused Bakta output after each genome (json, embl, gbff, svg, ...).")
    parser.add_argument("--log-file", type=Path, default=None)
    args = parser.parse_args(argv)

    configure_logging("INFO", logfile=args.log_file)

    if args.prune:
        result = prune_all(args.out_dir)
        print(f"pruned {result['files_removed']} files, "
              f"freed {result['bytes_freed'] / 1e9:.2f} GB from {args.out_dir}")
        return 0

    if not args.prepare and not args.run:
        parser.error("choose --prepare, --run or --prune")

    database = args.db or default_db()
    cohort = load_cohort(args.manifest)
    if not cohort:
        raise PipelineError("Cohort is empty; run --prepare-only first")

    pending = [c for c in cohort if not is_complete(args.out_dir / c[0])]
    done = len(cohort) - len(pending)
    ready, why = database_ready(database)

    print(f"cohort        : {len(cohort)} analysable assemblies")
    print(f"already done  : {done}")
    print(f"pending       : {len(pending)}")
    print(f"bakta database: {database}")
    print(f"database ready: {ready}  ({why})")
    print(f"threads/job   : {args.threads}  jobs: {args.jobs}")

    if not ready:
        print(
            f"\nERROR: Bakta database is not usable: {why}",
            file=sys.stderr,
        )
        print(
            "Provision it first: micromamba run -n pilot100 bakta_db download "
            "--type light --output <dir>",
            file=sys.stderr,
        )
        return 1

    if args.prepare:
        print("\nReady. Launch with: python3 scripts/run_annotation_sweep.py --run")
        return 0

    if not pending:
        print("\nNothing to do: every sample is already annotated.")
        return 0

    if args.limit:
        pending = pending[: args.limit]
        print(f"limiting to {len(pending)} sample(s) via --limit")

    previous: Dict[str, Dict[str, object]] = {}
    if STATE.exists():
        try:
            previous = json.loads(STATE.read_text(encoding="utf-8")).get("samples", {})
        except (json.JSONDecodeError, OSError):
            LOGGER.warning("Could not read prior state; starting fresh")

    records: Dict[str, Dict[str, object]] = dict(previous)
    started = time.monotonic()

    def work(item):
        pilot_id, genome = item
        return run_one(
            pilot_id, genome, args.out_dir / pilot_id, database,
            args.threads, args.species, args.timeout,
        )

    if args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for record in pool.map(work, pending):
                records[record["pilot_id"]] = record
                print(
                    f"  {record['pilot_id']:14s} {record['status']:10s} "
                    f"{record.get('seconds', 0)}s",
                    flush=True,
                )
                if record["status"] in ("completed", "reused") and not args.no_prune:
                    prune_outputs(args.out_dir / record["pilot_id"])
                save_state(records)
    else:
        for item in pending:
            record = work(item)
            records[record["pilot_id"]] = record
            print(
                f"  {record['pilot_id']:14s} {record['status']:10s} "
                f"{record.get('seconds', 0)}s",
                flush=True,
            )
            if record["status"] in ("completed", "reused") and not args.no_prune:
                prune_outputs(args.out_dir / record["pilot_id"])
            save_state(records)

    elapsed = time.monotonic() - started
    completed = sum(1 for r in records.values() if r.get("status") == "completed")
    failed = [r for r in records.values() if r.get("status") in ("failed", "timeout", "error")]

    rows = [
        {
            "Pilot_ID": pilot_id,
            "status": record.get("status"),
            "seconds": record.get("seconds"),
            "out_dir": record.get("out_dir", "."),
            "gff": record.get("gff", "."),
            "inference_tsv": record.get("inference_tsv", "."),
            "error": str(record.get("error", "."))[:300],
        }
        for pilot_id, record in sorted(records.items())
    ]
    write_tsv(
        args.out_dir.parent / "annotation_status.tsv",
        rows,
        ["Pilot_ID", "status", "seconds", "out_dir", "gff", "inference_tsv", "error"],
        header_comment=[
            "Stage 2 Bakta annotation status.",
            f"sweep wall time {elapsed:.0f}s",
        ],
    )
    save_state(records)

    print()
    print("Bakta sweep finished")
    print(f"  completed : {completed}")
    print(f"  failed    : {len(failed)}")
    print(f"  wall time : {elapsed / 60:.1f} min")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
