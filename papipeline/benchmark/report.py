"""The benchmark report.

Written whether or not the benchmark ran. A blocked experiment still owes a
record of what was attempted, on what cohort, and why nothing was measured.

Every number in the output is either MEASURED or NA. Nothing is estimated,
interpolated, or carried over from a previous run.
"""

from __future__ import annotations

import platform
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from .bakta10 import NOT_RUN, ConfigResult, Preflight


def _fmt(value: Any, digits: int = 2) -> str:
    if value is None or value == "" or value == "NA":
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_report(
    path: Path,
    *,
    cohort: Sequence[Any],
    results: Sequence[ConfigResult],
    checks: Mapping[int, Preflight],
    audit: Dict[str, Any],
    threads: int,
    database: Path,
    data_dir: Path,
) -> Path:
    ran = [r for r in results if r.status == "OK"]
    blocked = [r for r in results if r.status != "OK"]
    any_check = checks[max(checks)] if checks else Preflight(ok=False)

    lines: List[str] = []
    add = lines.append

    add("# Bakta 10-genome concurrency benchmark")
    add("")
    add("Controlled measurement of how concurrent Bakta processes affect wall")
    add("time, CPU, RAM and throughput over a fixed 10-genome cohort, holding")
    add("everything else constant.")
    add("")
    add("Every figure below is **MEASURED** or **NA**. Nothing here is")
    add("estimated, interpolated, or carried over. `NA` means the measurement")
    add("was not obtained.")
    add("")

    add("## 1. Outcome summary")
    add("")
    if ran:
        add(f"**{len(ran)} of {len(results)} configurations were measured.**")
    else:
        add("**No configuration was measured. The benchmark did not run.**")
        add("")
        add("The reasons are recorded verbatim in section 3. No numbers are")
        add("reported, because reporting any would mean inventing them.")
    if blocked:
        add("")
        for row in blocked:
            add(f"- `workers={row.workers}`: **{row.status}** — {row.reason}")
    add("")

    add("## 2. Cohort (MEASURED)")
    add("")
    add(f"Selected from the assemblies that physically exist under `{data_dir}`,")
    add("ordered by the repository's own natural GCA accession sort")
    add("(`papipeline.pilot.cohort.discover_assemblies` + `sorted_accessions`).")
    add("")
    add("| # | Accession | Assembly file | Size |")
    add("| --- | --- | --- | --- |")
    for genome in cohort:
        add(f"| {genome.rank} | `{genome.accession}` | `{genome.path.name}` | "
            f"{genome.path.stat().st_size / 1e6:.1f} MB |")
    add("")
    add("This is **not** the first 10 rows of `PDC_essential.tsv`. That")
    add(f"ordering selects a materially different population: the overlap is")
    add(f"**{audit.get('overlap_with_pdc_first_n')} of {len(cohort)}**, and the")
    add(f"two orders are identical: **{audit.get('identical_to_pdc_first_n')}**.")
    add("The same cohort file is reused by every configuration, so no sweep can")
    add("silently compare different populations.")
    add("")

    add("## 3. Environment and preflight (MEASURED)")
    add("")
    add("| Item | Value |")
    add("| --- | --- |")
    add(f"| Bakta executable | {_fmt(any_check.bakta_path)} |")
    add(f"| Bakta version | {_fmt(any_check.bakta_version)} |")
    add(f"| Bakta database | {_fmt(any_check.database or str(database))} |")
    db_bytes = any_check.database_bytes
    add(f"| Bakta database size | "
        f"{'NA' if db_bytes is None else f'{db_bytes / 2**30:.2f} GB'} |")
    add(f"| Bakta `--threads` | {threads} |")
    add(f"| Python | {platform.python_version()} |")
    add(f"| Platform | {platform.platform()} |")
    add(f"| Logical CPUs | {any_check.cpu_count} |")
    add(f"| System RAM | {_fmt(any_check.ram_gb)} GB |")
    add(f"| psutil available | {'yes' if _psutil() else 'no'} |")
    add("")
    if any_check.reasons:
        add("Preflight findings:")
        add("")
        for reason in any_check.reasons:
            add(f"- {reason}")
        add("")

    add("## 4. Worker count vs Bakta thread count")
    add("")
    add("These are separate variables and were not varied together. The sweep")
    add(f"varies only the number of concurrent Bakta processes. `--threads` was")
    add(f"held at **{threads}** throughout.")
    add("")
    add("Concurrency is a count of concurrently running processes. There is no")
    add("persistent worker pool in this build, so no worker identity is reported")
    add("and none is implied.")
    add("")

    add("## 5. Results (MEASURED)")
    add("")
    if not ran:
        add("No configuration produced measurements. See section 3.")
        add("")
    else:
        add("| workers | threads | genomes | completed | failed | invalid | "
            "incomplete | retries | wall s | CPU s | peak RAM GB | avg RAM GB |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in ran:
            add(
                f"| {row.workers} | {row.bakta_threads} | {row.genomes} | "
                f"{row.completed} | {row.failed} | {row.invalid} | "
                f"{row.incomplete} | {row.retried} | "
                f"{_fmt(row.wall_seconds)} | {_fmt(row.cpu_seconds)} | "
                f"{_fmt(row.peak_ram_gb)} | {_fmt(row.average_ram_gb)} |"
            )
        add("")
        add("| workers | avg genome s | min s | max s | genomes/hour | "
            "genomes/min | CPU % mean | CPU % peak |")
        add("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in ran:
            summary = row.to_summary_row()
            add(
                f"| {row.workers} | {_fmt(summary['avg_genome_seconds'])} | "
                f"{_fmt(summary['min_genome_seconds'])} | "
                f"{_fmt(summary['max_genome_seconds'])} | "
                f"{_fmt(summary['genomes_per_hour'])} | "
                f"{_fmt(summary['throughput_genomes_per_min'], 4)} | "
                f"{_fmt(summary['cpu_percent_mean'], 1)} | "
                f"{_fmt(summary['cpu_percent_peak'], 1)} |"
            )
        add("")

    add("## 6. Per-genome runtimes (MEASURED)")
    add("")
    if not ran or not any(r.per_genome for r in ran):
        add("No per-genome measurements were taken.")
        add("")
    else:
        for row in ran:
            add(f"### workers = {row.workers}")
            add("")
            add("| # | Accession | state | attempts | wall s | CPU s | "
                "peak RSS MB | avg RSS MB | exit |")
            add("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
            for genome in row.per_genome:
                add(
                    f"| {genome.rank} | `{genome.accession}` | "
                    f"{genome.state} | {genome.attempts} | "
                    f"{_fmt(genome.wall_seconds)} | {_fmt(genome.cpu_seconds)} | "
                    f"{_fmt(genome.peak_rss_mb, 1)} | "
                    f"{_fmt(genome.average_rss_mb, 1)} | "
                    f"{_fmt(genome.exit_code)} |"
                )
            add("")

    add("## 7. Validation")
    add("")
    add("Every genome's Bakta output was checked against the existing")
    add("`bakta_spec` contract: the feature table must exist, be non-empty,")
    add("parse, carry the `Gene` and `Locus Tag` columns, and the GFF and TSV")
    add("must agree in record count.")
    add("")
    add("**Exit code 0 is never counted as success on its own.** A task is")
    add("`SUCCEEDED` only when the process exited cleanly *and* the contract")
    add("passed.")
    add("")
    if ran:
        add("| workers | SUCCEEDED | FAILED | INVALID | INCOMPLETE |")
        add("| --- | --- | --- | --- | --- |")
        for row in ran:
            add(f"| {row.workers} | {row.completed} | {row.failed} | "
                f"{row.invalid} | {row.incomplete} |")
        add("")
    else:
        add("No configuration ran, so no validation outcomes were produced.")
        add("")

    add("## 8. Output isolation")
    add("")
    add("Each concurrency level writes to its own `workers_N/` directory,")
    add("including its own SQLite execution store. Directories are removed")
    add("before use. A result therefore cannot be accepted from an earlier")
    add("configuration, and a partial earlier run cannot be mistaken for a")
    add("current one.")
    add("")

    add("## 9. Scaling behaviour")
    add("")
    if len(ran) >= 2:
        base_row = min(ran, key=lambda r: r.workers)
        add(f"Relative to `workers={base_row.workers}` (wall {_fmt(base_row.wall_seconds)} s):")
        add("")
        add("| workers | wall s | speed-up | efficiency | peak RAM GB |")
        add("| --- | --- | --- | --- | --- |")
        for row in sorted(ran, key=lambda r: r.workers):
            if base_row.wall_seconds:
                speedup = base_row.wall_seconds / row.wall_seconds if row.wall_seconds else None
                add(f"| {row.workers} | {_fmt(row.wall_seconds)} | "
                    f"{_fmt(speedup)}x | "
                    f"{_fmt(None if speedup is None else speedup / (row.workers / base_row.workers) * 100, 0)}% | "
                    f"{_fmt(row.peak_ram_gb)} |")
        add("")
        add("Efficiency is speed-up divided by the concurrency ratio: 100% is")
        add("perfect linear scaling, lower means the machine is already")
        add("saturated.")
    else:
        add("**Insufficient measurements to characterise scaling.** At least two")
        add("configurations are needed, and none completed.")
    add("")

    add("## 10. Optimal concurrency")
    add("")
    if len(ran) >= 2:
        best = max(ran, key=lambda r: (r.completed / r.wall_seconds) if r.wall_seconds else 0)
        add(f"By throughput, `workers={best.workers}` was fastest "
            f"({best.completed / best.wall_seconds * 60:.4f} genomes/min) on the")
        add("data measured here.")
        add("")
        add("This is a statement about **this machine, this database, this")
        add("cohort and this `--threads` value**, not a general recommendation.")
        add("It is not extrapolated to other hardware.")
    else:
        add("**Not determined.** No configuration produced measurements, so no")
        add("optimal worker count can be claimed. Declaring one from fewer than")
        add("two measured points would be a guess dressed as a finding.")
    add("")

    add("## 11. Limitations")
    add("")
    add("- One machine, one database, one cohort of 10 related *P. aeruginosa*")
    add("  isolates. Transfer to other organisms or hardware is not claimed.")
    add("- `--threads` was held constant. Optimising worker count and thread")
    add("  count together is a separate experiment and was not attempted.")
    add("- Each configuration is a single run, not a set of replicates, so the")
    add("  spread is unknown and small differences should not be read as real.")
    add("- Peak and average RAM are host-wide, not per process; per-process RSS")
    add("  is in the per-genome table.")
    add("- Bakta's database is resident, so concurrency multiplies memory. The")
    add("  preflight refuses configurations that would obviously exhaust RAM")
    add("  rather than risking a swap or an OOM kill.")
    add("- The pipeline itself remains strictly sequential. This benchmark")
    add("  measures the tool; it does not change how the pipeline schedules.")
    add("")
    add("## 12. What was NOT done")
    add("")
    add("No adaptive scheduler, dynamic worker count, RAM-based scaling, queue")
    add("optimisation, database-sharing change or GPU work was implemented. This")
    add("milestone is measurement only.")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _psutil() -> bool:
    try:
        import psutil  # noqa: F401
        return True
    except ImportError:
        return False
