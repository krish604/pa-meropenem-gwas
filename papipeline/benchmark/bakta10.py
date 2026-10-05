"""A controlled Bakta concurrency benchmark over a fixed 10-genome cohort.

This is an *experiment*, not a scheduler. Its only job is to measure how
wall time, CPU, RAM and throughput respond to the number of concurrent
Bakta processes, holding everything else fixed, so that a scheduler design
can later be argued from data rather than intuition.

Design constraints that shape the code:

* **The cohort is fixed and auditable.** It is the first ten genome
  assemblies that physically exist under ``data/``, ordered by the
  repository's own natural GCA sort (``papipeline.pilot.cohort``). It is
  explicitly *not* the first rows of ``PDC_essential.tsv``; that ordering
  selects a materially different population, and :func:`audit_cohort`
  records the overlap so the distinction is evidence rather than a claim.
* **Real Bakta, or nothing.** Every task goes through
  :func:`papipeline.adapters.bakta.run_bakta` and
  :func:`papipeline.execution.run_task`. There is no simulated Bakta and no
  benchmark-only success path. If Bakta cannot run, the benchmark records
  ``NOT_RUN`` with the reason and reports no numbers.
* **Worker count and Bakta thread count are separate variables.** The sweep
  varies only the number of concurrent processes. ``--threads`` is recorded
  and left at whatever the configuration says.
* **Outputs are isolated per configuration.** Each concurrency level writes
  to its own directory, so a result can never be accepted from an earlier
  configuration.
* **Exit code 0 is not success.** Every output goes through the existing
  ``bakta_spec`` contract.
"""

from __future__ import annotations

import csv
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    import psutil
except ImportError:  # pragma: no cover
    psutil = None  # type: ignore[assignment]

from ..execution.state import StageState
from ..execution.store import ExecutionStore
from ..observatory.events import EventBus, bus_sink
from .telemetry import ProcessSample, SampleSystem, SystemSampler, sample_process

#: The concurrency levels to measure. 5 is omitted deliberately: it is not
#: needed to see the shape of the curve, and every extra configuration is
#: minutes of Bakta time.
WORKER_LEVELS: Tuple[int, ...] = (1, 2, 3, 4, 6)

N_GENOMES = 10
NOT_RUN = "NOT_RUN"


# --------------------------------------------------------------------------
# Cohort
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class BenchGenome:
    """One member of the benchmark cohort."""

    rank: int
    accession: str
    path: Path

    @property
    def stem(self) -> str:
        return Path(self.path).stem


def select_cohort(data_dir: Path, n: int = N_GENOMES) -> List[BenchGenome]:
    """The first ``n`` assemblies under ``data/`` in the repository's order.

    Delegates to :mod:`papipeline.pilot.cohort`, which is the pipeline's own
    discovery and ordering logic, so the benchmark cohort is the same kind
    of population the pilot used and not a second, subtly different one.
    """
    from ..pilot.cohort import discover_assemblies, sorted_accessions

    if not Path(data_dir).exists():
        raise FileNotFoundError(f"assembly directory does not exist: {data_dir}")
    found = discover_assemblies(Path(data_dir))
    if not found:
        # The directory exists but holds no assembly: that is a different
        # problem from a missing directory, and conflating them hides which.
        raise ValueError(
            f"no assembly files found under {data_dir}; the benchmark "
            f"cohort is selected from the filesystem, so an empty data "
            f"directory cannot produce a cohort")
    chosen = sorted_accessions(found)[:n]
    if len(chosen) < n:
        raise ValueError(
            f"only {len(chosen)} assemblies available under {data_dir}; "
            f"the benchmark cohort requires exactly {n}")
    out: List[BenchGenome] = []
    for rank, accession in enumerate(chosen, start=1):
        paths = found[accession]
        if len(paths) > 1:
            # Deterministic: the lexicographically first path, and the fact
            # is recorded rather than silently resolved.
            paths = sorted(paths, key=lambda p: str(p))
        out.append(BenchGenome(rank=rank, accession=accession, path=paths[0]))
    return out


def audit_cohort(cohort: Sequence[BenchGenome], pdc_table: Path) -> Dict[str, Any]:
    """Evidence that this cohort is not the first rows of the PDC table.

    The requirement is explicit about this, and "we did not use PDC_essential"
    is exactly the kind of claim that should be checkable rather than
    trusted, so the overlap is computed and written into the run directory.
    """
    result: Dict[str, Any] = {
        "cohort_n": len(cohort),
        "pdc_table": str(pdc_table),
        "pdc_table_present": Path(pdc_table).exists(),
        "cohort_accessions": [g.accession for g in cohort],
    }
    if not Path(pdc_table).exists():
        result["pdc_first_n"] = []
        result["overlap_with_pdc_first_n"] = None
        result["identical_to_pdc_first_n"] = None
        return result
    with open(pdc_table, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        first_n: List[str] = []
        for row in reader:
            value = (row.get("Assembly") or "").strip()
            if value:
                first_n.append(value)
            if len(first_n) >= len(cohort):
                break
    cohort_ids = [g.accession for g in cohort]
    result["pdc_first_n"] = first_n
    result["overlap_with_pdc_first_n"] = len(set(first_n) & set(cohort_ids))
    result["identical_to_pdc_first_n"] = first_n == cohort_ids
    return result


def write_cohort_manifest(path: Path, cohort: Sequence[BenchGenome]) -> Path:
    """Persist the cohort, so every configuration provably used this one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["rank", "accession", "assembly_path", "size_bytes"])
        for genome in cohort:
            writer.writerow([genome.rank, genome.accession,
                             str(genome.path), genome.path.stat().st_size])
    return path


def read_cohort_manifest(path: Path) -> List[BenchGenome]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            BenchGenome(rank=int(row["rank"]), accession=row["accession"],
                        path=Path(row["assembly_path"]))
            for row in reader
        ]


# --------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------


@dataclass
class Preflight:
    """Whether the benchmark can run, and exactly why not."""

    ok: bool
    reasons: List[str] = field(default_factory=list)
    bakta_path: Optional[str] = None
    bakta_version: Optional[str] = None
    database: Optional[str] = None
    database_bytes: Optional[int] = None
    database_version: Optional[str] = None
    threads: int = 1
    ram_gb: Optional[float] = None
    cpu_count: int = 0
    missing_dependencies: List[str] = field(default_factory=list)

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


#: Measured peak RSS of one real Bakta 1.12.1 process annotating
#: GCA_000710625.1 against the db-light database at 4 threads on this
#: machine: 1.92 GB (/usr/bin/time -l). The safety check below uses this
#: rather than a guess, but still refuses a configuration with little
#: headroom, because a memory figure from one genome is an estimate for a
#: concurrent cohort and must not be trusted to the byte.
MEASURED_PEAK_RSS_GB = 1.92

#: Tools Bakta resolves through PATH. tRNAscan-SE is the important one: when
#: it is missing, Bakta aborts within seconds and a careless harness can
#: mistake that for "Bakta is fast".
BAKTA_DEPENDENCIES = ("tRNAscan-SE", "diamond", "blastn", "makeblastdb",
                      "aragorn")


def resolve_bakta(explicit: Optional[Path] = None) -> Optional[str]:
    """Locate Bakta: explicit argument, then ``BAKTA_BIN``, then ``PATH``."""
    if explicit is not None:
        return str(Path(explicit).expanduser())
    from ..adapters.bakta import resolve_executable

    try:
        return resolve_executable()
    except Exception:  # noqa: BLE001 - absence is the answer, not a crash
        return None


def database_version_of(database: Path) -> Optional[str]:
    """Read a Bakta database's own ``version.json``, if it has one."""
    manifest = Path(database) / "version.json"
    if not manifest.exists():
        return None
    try:
        import json
        data = json.loads(manifest.read_text(encoding="utf-8"))
        return (f"{data.get('major')}.{data.get('minor')} "
                f"({data.get('type')}, {data.get('date')})")
    except (OSError, ValueError):
        return None


def bakta_version(executable: Optional[str] = None) -> Optional[str]:
    if not executable:
        return None
    try:
        out = subprocess.run([executable, "--version"], capture_output=True,
                             text=True, timeout=30, check=False)
        text = (out.stdout or out.stderr or "").strip()
        return text.splitlines()[0] if text else None
    except (OSError, subprocess.SubprocessError):
        return None


def missing_dependencies(executable: Optional[str]) -> List[str]:
    """Bakta's helper tools, as the child process would see them.

    Resolved against the executable's own directory, because that is the
    only way they are found when the environment is not on ``PATH``.
    """
    if not executable:
        return list(BAKTA_DEPENDENCIES)
    bindir = str(Path(executable).resolve().parent)
    missing: List[str] = []
    for tool in BAKTA_DEPENDENCIES:
        if not (Path(bindir) / tool).exists() and not shutil.which(tool):
            missing.append(tool)
    return missing


def directory_size(path: Path) -> Optional[int]:
    try:
        return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())
    except OSError:
        return None


def preflight(
    database: Path,
    threads: int,
    max_workers: int,
    *,
    bakta: Optional[Path] = None,
) -> Preflight:
    """Check the things that would make a run meaningless or destructive.

    Reported as data, not raised, so a blocked benchmark still leaves a
    record explaining itself.
    """
    result = Preflight(ok=True, threads=threads, cpu_count=os.cpu_count() or 0)
    if psutil is not None:
        result.ram_gb = round(psutil.virtual_memory().total / 2 ** 30, 2)

    result.bakta_path = resolve_bakta(bakta)
    if not result.bakta_path:
        result.ok = False
        result.reasons.append(
            "bakta is not on PATH and BAKTA_BIN is unset; the "
            "benchmark requires the real executable and will not "
            "substitute a simulation")
    else:
        if not Path(result.bakta_path).exists():
            result.ok = False
            result.reasons.append(
                f"bakta executable does not exist: {result.bakta_path}")
        else:
            result.bakta_version = bakta_version(result.bakta_path)
        absent = missing_dependencies(result.bakta_path)
        if absent:
            result.missing_dependencies = absent
            result.ok = False
            result.reasons.append(
                "Bakta's helper tools are not resolvable next to the "
                f"executable or on PATH: {', '.join(absent)}. Bakta aborts "
                "within seconds when tRNAscan-SE is missing, which can be "
                "mistaken for a fast run"
            )

    if not Path(database).exists():
        result.ok = False
        result.reasons.append(
            f"bakta database directory does not exist: {database}")
    else:
        result.database = str(database)
        result.database_bytes = directory_size(Path(database))
        result.database_version = database_version_of(Path(database))

    # Safety: refuse a configuration that would obviously exhaust memory.
    # Bakta's database is resident, so concurrency multiplies the cost.
    if result.ram_gb:
        per_worker = MEASURED_PEAK_RSS_GB
        need = per_worker * max_workers
        if need > result.ram_gb * 0.9:
            result.ok = False
            result.reasons.append(
                f"{max_workers} concurrent Bakta processes need roughly "
                f"{need:.1f} GB ({per_worker} GB each, measured) but only "
                f"{result.ram_gb} GB is available; this configuration would "
                f"risk swapping or OOM")
    return result


# --------------------------------------------------------------------------
# One configuration
# --------------------------------------------------------------------------


@dataclass
class GenomeResult:
    """The measured outcome for one genome at one concurrency level."""

    rank: int
    accession: str
    state: str
    attempts: int
    wall_seconds: Optional[float] = None
    cpu_seconds: Optional[float] = None
    peak_rss_mb: Optional[float] = None
    average_rss_mb: Optional[float] = None
    cpu_percent: Optional[float] = None
    exit_code: Optional[int] = None
    pid: Optional[int] = None
    bakta_threads: Optional[int] = None
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    detail: str = ""
    log_path: Optional[str] = None

    def to_row(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ConfigResult:
    """One row of the sweep."""

    workers: int
    bakta_threads: int
    genomes: int
    completed: int = 0
    failed: int = 0
    invalid: int = 0
    incomplete: int = 0
    retried: int = 0
    wall_seconds: float = 0.0
    cpu_seconds: Optional[float] = None
    peak_ram_gb: Optional[float] = None
    average_ram_gb: Optional[float] = None
    cpu_percent_mean: Optional[float] = None
    cpu_percent_peak: Optional[float] = None
    load_average_1m: Optional[float] = None
    disk_read_bps: Optional[float] = None
    disk_write_bps: Optional[float] = None
    per_genome: List[GenomeResult] = field(default_factory=list)
    status: str = "OK"
    reason: str = ""
    run_key: str = ""

    def to_summary_row(self) -> Dict[str, Any]:
        walls = [g.wall_seconds for g in self.per_genome
                 if g.wall_seconds is not None and g.state == StageState.SUCCEEDED.value]
        return {
            "workers": self.workers,
            "bakta_threads": self.bakta_threads,
            "genomes": self.genomes,
            "completed": self.completed,
            "failed": self.failed,
            "invalid": self.invalid,
            "incomplete": self.incomplete,
            "retry_count": self.retried,
            "wall_seconds": round(self.wall_seconds, 2),
            "cpu_seconds": NA if self.cpu_seconds is None else round(self.cpu_seconds, 2),
            "peak_ram_gb": NA if self.peak_ram_gb is None else round(self.peak_ram_gb, 2),
            "average_ram_gb": NA if self.average_ram_gb is None else round(self.average_ram_gb, 2),
            "avg_genome_seconds": NA if not walls else round(statistics.mean(walls), 2),
            "min_genome_seconds": NA if not walls else round(min(walls), 2),
            "max_genome_seconds": NA if not walls else round(max(walls), 2),
            "genomes_per_hour": NA if not self.wall_seconds
            else round(self.completed / self.wall_seconds * 3600, 2),
            "throughput_genomes_per_min": NA if not self.wall_seconds
            else round(self.completed / self.wall_seconds * 60, 4),
            "cpu_percent_mean": NA if self.cpu_percent_mean is None else round(self.cpu_percent_mean, 1),
            "cpu_percent_peak": NA if self.cpu_percent_peak is None else round(self.cpu_percent_peak, 1),
            "status": self.status,
            "reason": self.reason,
        }


NA = "NA"

SUMMARY_COLUMNS = [
    "workers", "bakta_threads", "genomes", "completed", "failed", "invalid",
    "incomplete", "retry_count", "wall_seconds", "cpu_seconds", "peak_ram_gb",
    "average_ram_gb", "avg_genome_seconds", "min_genome_seconds",
    "max_genome_seconds", "genomes_per_hour", "throughput_genomes_per_min",
    "cpu_percent_mean", "cpu_percent_peak", "status", "reason",
]


# --------------------------------------------------------------------------
# Running one configuration
# --------------------------------------------------------------------------


def run_configuration(
    *,
    workers: int,
    cohort: Sequence[BenchGenome],
    out_dir: Path,
    database: Path,
    threads: int,
    config: Any,
    run_key: str,
    bus: EventBus,
    genome_dir: Path,
    max_attempts: int = 1,
) -> ConfigResult:
    """Run the whole cohort at one concurrency level.

    Uses a thread pool only to overlap *waiting on subprocesses*; each
    genome's annotation is a separate real Bakta process. The threads are
    not workers in any persistent sense, which is why the result reports a
    concurrency and a process count rather than worker identities.
    """
    from ..adapters.bakta import run_bakta
    from ..manifest import SampleManifest
    from ..models import Sample

    out_dir.mkdir(parents=True, exist_ok=True)
    store = ExecutionStore(out_dir / "execution.db")
    sampler = SystemSampler(interval=1.0)
    sampler.start()
    sink = bus_sink(bus)
    results: List[Optional[GenomeResult]] = [None] * len(cohort)

    def one(index: int) -> None:
        genome = cohort[index]
        sample = Sample(sample_id=genome.accession,
                        assembly_path=str(genome.path), source="benchmark")
        holder: Dict[str, Any] = {}

        def record_pid(pid: int) -> None:
            holder["pid"] = pid
            # Sample the real child while it runs.
            holder["sampler"] = sample_process(
                pid, threads=threads, sink=sink, run_key=run_key,
                stage="annotation", subject=genome.accession)

        try:
            bakta_result = run_bakta(
                sample, config=config, genomes_dir=genome_dir,
                out_root=out_dir / "bakta", database=database,
                run_key=run_key, store=store, event_sink=sink,
                config_hash=None, pid_sink=record_pid,
            )
        except Exception as exc:  # noqa: BLE001 - recorded, not hidden
            results[index] = GenomeResult(
                rank=genome.rank, accession=genome.accession,
                state=StageState.FAILED.value, attempts=0,
                detail=f"{type(exc).__name__}: {exc}",
                bakta_threads=threads,
            )
            return

        observed = holder.get("sampler") or ProcessSample()
        if hasattr(observed, "finish"):
            observed.finish()
        row = store.get(run_key, "annotation", genome.accession) or {}
        # exit_code lives in attempt_log, not on the execution row, for the
        # same reason log_path does: the runner records per-attempt facts
        # per attempt.
        history = store.attempts(run_key, "annotation", genome.accession)
        exit_code = history[-1].exit_code if history else None
        results[index] = GenomeResult(
            rank=genome.rank, accession=genome.accession,
            state=bakta_result.state, attempts=bakta_result.attempts,
            wall_seconds=observed.wall_seconds or row.get("elapsed_seconds"),
            cpu_seconds=observed.cpu_seconds,
            peak_rss_mb=observed.peak_rss_mb,
            average_rss_mb=observed.average_rss_mb,
            cpu_percent=observed.cpu_percent_mean,
            exit_code=exit_code,
            pid=observed.pid or holder.get("pid"),
            bakta_threads=threads,
            started_at=observed.started_at,
            ended_at=observed.ended_at,
            detail=bakta_result.validation_detail[:400],
            log_path=store.log_path(run_key, "annotation", genome.accession),
        )

    # The store must exist before any task writes to it.
    store.close()
    store = ExecutionStore(out_dir / "execution.db")
    started = time.time()
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(one, range(len(cohort))))
    finally:
        wall = time.time() - started
        system = sampler.stop()
        store.close()

    per_genome = [r for r in results if r is not None]
    completed = sum(1 for r in per_genome if r.state == StageState.SUCCEEDED.value)
    out = ConfigResult(
        workers=workers, bakta_threads=threads, genomes=len(cohort),
        completed=completed,
        failed=sum(1 for r in per_genome if r.state == StageState.FAILED.value),
        invalid=sum(1 for r in per_genome if r.state == StageState.INVALID.value),
        incomplete=sum(1 for r in per_genome if r.state == StageState.INCOMPLETE.value),
        retried=sum(1 for r in per_genome if r.attempts > 1),
        wall_seconds=wall,
        per_genome=per_genome,
        run_key=run_key,
    )
    for genome in per_genome:
        if genome.cpu_seconds:
            out.cpu_seconds = (out.cpu_seconds or 0.0) + genome.cpu_seconds
    out.peak_ram_gb = system.peak_ram_gb
    out.average_ram_gb = system.average_ram_gb
    out.cpu_percent_mean = system.cpu_percent_mean
    out.cpu_percent_peak = system.cpu_percent_peak
    out.load_average_1m = system.load_average_1m
    out.disk_read_bps = system.disk_read_bps
    out.disk_write_bps = system.disk_write_bps
    return out


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def write_summary(path: Path, rows: Sequence[ConfigResult]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row.to_summary_row())
    return path


def write_environment(path: Path, pre: Preflight, extra: Optional[Dict[str, Any]] = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "bakta_path": pre.bakta_path or NA,
        "bakta_version": pre.bakta_version or NA,
        "bakta_database": pre.database or NA,
        "bakta_database_bytes": pre.database_bytes if pre.database_bytes is not None else NA,
        "bakta_threads": pre.threads,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "cpu_count": pre.cpu_count,
        "ram_gb": pre.ram_gb if pre.ram_gb is not None else NA,
        "psutil": "yes" if psutil is not None else "no",
        "preflight_ok": "yes" if pre.ok else "no",
        "preflight_reasons": " | ".join(pre.reasons) if pre.reasons else NA,
    }
    if extra:
        payload.update(extra)
    with open(path, "w", encoding="utf-8") as handle:
        for key, value in payload.items():
            handle.write(f"{key}\t{value}\n")
    return path


def write_per_genome(path: Path, rows: Sequence[ConfigResult]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Derive the header from the dataclass rather than from an instance:
    # GenomeResult has required fields, so a bare construction is invalid.
    import dataclasses as _dc
    fields = [f.name for f in _dc.fields(GenomeResult)]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            for genome in row.per_genome:
                writer.writerow(genome.to_row())
    return path
