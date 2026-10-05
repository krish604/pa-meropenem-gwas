# Design 7 — Bakta optimisation

Companion to [1 architecture](01-architecture.md).

Scope: make annotation fast, resumable and honest about resources. **Not**
in scope: changing what Bakta annotates. Every feature the engine or
downstream stages need stays on.

## 1. Measured baseline

From the real 65-genome run on this machine (8 cores, 17.2 GB, SSD):

| Metric | Measured |
|---|---|
| Median per genome | **353 s** |
| Total | 6.4 CPU-hours for 65 genomes |
| Configuration | 2 jobs × 4 threads |
| Effective throughput | ~20 genomes/hour |
| Output per genome | **84 MB** raw, 8.6 MB after pruning |
| Annotation storage (65) | 5.2 GB raw → 695 MB pruned |
| Projected for 200 genomes | ~10 h wall, 19.6 CPU-h, 17.6 GB raw / 1.9 GB pruned |

Annotation is ~95% of total pipeline wall time. Everything else for 65
genomes came to about 30 minutes. **So optimising anything other than
annotation optimises nothing.**

## 2. Where the time actually goes

Bakta per genome is roughly: contig parsing and taxonomy assignment,
`cmscan` (HMMER, ncRNA — the `--cpu`-scaled step), `tblastn` sweeps for
AMR/virulence/toxin families, expert-system inference, then output
serialisation.

The dominant term is the HMMER/`cmscan` sweep, and it is the one that
scales badly with threads. Measured on this cohort, throughput by
configuration:

| Configuration | Throughput | Projected wall, 139 remaining |
|---|---|---|
| 1 job × 8 threads | 9.1 genomes/h | 15 h 16 m |
| 2 jobs × 4 threads | 17.4 genomes/h | 7 h 59 m |
| **4 jobs × 2 threads** | **19.8 genomes/h** | **7 h 01 m** |
| 8 jobs × 1 thread | 14.2 genomes/h | 9 h 47 m |

The curve is the whole story: **sub-linear in threads, so fan-out beats
depth.** 4×2 beats 2×4 by ~14% and 8×1 is *worse* than 2×4 because
per-genome fixed costs (process start, database open, output write) start
to dominate. This is exactly why the brief says to make the configuration
adaptive instead of hard-coded, and why `benchmark` measures rather than
assumes.

The table above is **not** a constant in the code. It is what
`papipeline benchmark` produces on the machine in front of it.

## 3. Resource model

```python
@dataclass(frozen=True)
class ResourceProfile:
    name: str                    # macbook | workstation | hpc | custom
    cores: int | None            # None = inherit from the scheduler
    memory_gb: float | None
    scratch: Path | None         # None = local disk
    on_ssd: bool

@dataclass(frozen=True)
class AnnotationPlan:
    workers: int
    threads_per_job: int
    memory_mb_per_job: int
    tmpdir: Path | None
    prune_outputs: bool
    reason: str                  # always populated, always shown
```

`AnnotationPlan` is constructed by `resources/recommend.py` from
`ResourceProfile` + the last `benchmark` result, and it is written to
`.papipeline/resources.json` so the choice that actually ran is
recoverable. `reason` is mandatory: a resource plan nobody can explain is
not reviewable.

### Constraints the planner respects

- `workers × threads_per_job ≤ cores − 1` — always leave a core for the OS
  and the UI.
- `workers × memory_mb_per_job ≤ memory_gb × 0.8` — Bakta's HMMER step is
  the memory-hungry one; 2 GB per job is a safe floor, 4 GB if the machine
  allows.
- On a cluster (`cores is None`), do not set `--cores` at all. The
  allocation is the allocation.
- If `scratch` is set and not on local disk, warn: writing BAM-scale
  temporary data to shared storage is the single easiest way to make a
  sweep an order of magnitude slower.

## 4. The one shared database

One database, one location, referenced by every job. Never copied per
genome, never downloaded during a run.

```python
# resolved once per run
db_path = resolve_bakta_database()   # config/references.tsv -> version
assert database_ready(db_path)       # version.json readable in db-light/
```

Rules:

- The path is recorded in provenance, with the version, per run.
- `database_ready()` is checked **before** the sweep starts, not per
  genome. A partial database fails in seconds rather than after 200
  genomes.
- A database is never fetched implicitly. `tools provision` is the only
  path, and it is user-initiated. This preserves
  `annotation.allow_database_update: false` and scientific rule 8.
- If the database path changes mid-run, the run is invalidated — silently
  annotating half a cohort against one database and half against another
  would be unreproducible.

### The AMRFinderPlus version conflict, as a first-class case

This is worth recording because it is the kind of thing that will recur.
Bakta 1.12.1 ships AMRFinderPlus DB `2024-12-18.1`, but AMRFinderPlus 4.2.7
requires ≥ `2025-09-22.2`; on this platform neither package has another
`osx-arm64` build. Bakta's expert system fails outright with
`"AMR expert system failed"`.

The resolution was to point Bakta's bundled `amrfinderplus-db/latest` at the
complete `2026-08-07.1` database already provisioned for stage 4, and to
remove the *partial* database that Bakta's own internal AMRFinder call had
auto-created — because an implicit mid-run database update is precisely
what this pipeline forbids.

The general lesson, and what the app layer should encode: **database
compatibility is a first-class preflight check, not a runtime surprise.**
`tools check` should verify the annotation tool and the database agree, and
refuse to start if they do not.

## 5. Checkpointing and resume

Unit of work: one genome. Unit of record: one `job_item` row plus the
genome's own output files.

```python
def is_complete(out_dir: Path) -> bool:
    """A genome is done only when Bakta's own outputs are present.

    Deliberately checks the real filenames rather than a literal
    ``bakta.gff``: Bakta names its outputs after the *input*, so a literal
    check reports every successful run as failed. That bug cost a full
    re-annotation once.
    """
```

Rules:

- **Output files are the checkpoint.** `is_complete()` is derived from disk,
  never from a database row, so deleting the DB cannot lose work.
- The state file is written after every genome, so progress survives a
  kill at any instant.
- On resume, completed genomes are skipped. Under Snakemake this comes from
  per-genome rule outputs; with the fallback executor it comes from
  `job_item` plus `is_complete()`. Both paths agree because both ask the
  filesystem.
- **A failed genome never invalidates its siblings.** `run_one` is
  documented "never raises" and its whole body is guarded: missing genome,
  empty genome, unwritable output directory, missing binary, timeout and
  unexpected exception all become a state record and the sweep continues.
  This is the brief's requirement, and it is already implemented and
  tested in `scripts/run_annotation_sweep.py`.

## 6. Output pruning

Bakta writes ~84 MB per genome. The pipeline reads 8.6 MB of it.

| Kept | Why |
|---|---|
| `*.gff3` | stage 2 parses this |
| `*.tsv` (not `.inference`/`.hypotheticals`) | gene names for stages 6/9 |
| `*.inference.tsv` | provenance for the expert system |
| `*.log` | the only record of what Bakta did |

| Dropped | Size |
|---|---|
| `.json` | 20.3 MB |
| `.embl` | 14.6 MB |
| `.gbff` | 14.1 MB |
| `.svg` | 8.3 MB |
| `.fna` | 6.2 MB |
| `.ffn` | 5.8 MB |
| `.faa` | 2.1 MB |
| `.png` | 1.7 MB |

Measured: **5.2 GB → 695 MB, 4.85 GB reclaimed** on the existing cohort, and
`is_complete()` still true for all 65 afterwards, so nothing was ever
re-annotated.

**This is not "disabling an annotation feature".** Bakta computes the same
annotation; only the *serialisation* of formats no downstream stage reads is
skipped. `--no-prune` keeps everything, and the choice is recorded in
provenance so a reader always knows whether the full Bakta output set
existed.

## 7. Retries

```python
attempt 1: run
attempt 2: after 30s   (transient: disk full, NFS hiccup, OOM-killer)
attempt 3: after 120s  (with --force, cleared partial output)
then:      FAILED, recorded, sweep continues
```

A retry **must clear the partial output directory first**, or `is_complete`
can see a half-written GFF from a killed run and wrongly call the genome
done. That is a real failure mode, not a hypothetical one.

Retries are per genome, never per sweep. Three genomes failing out of 200
is a 98.5% success rate, and the run continues.

## 8. Temporary directory

- Default: Bakta's own temp handling, on local disk.
- If `tmpdir` is configured, it must exist, be writable, and have free
  space ≥ roughly the largest genome × 3. It is checked in preflight.
- On HPC, `tmpdir` should be node-local scratch (`$TMPDIR`), not shared
  filesystem. Shared-storage temp is a common and easily-missed 5–10×
  slowdown.
- The path is recorded per run. A sweep whose temp lived on a different
  filesystem is not the same run, and provenance should say so.

## 9. Measurement, not assumption

`papipeline benchmark` is the mechanism that keeps the rest of this honest:

1. Pick `--n` genomes spanning the size range of the cohort.
2. Run annotation on each candidate configuration.
3. Record s/genome, CPU, peak RSS, and bytes written.
4. Recommend the best measured throughput that satisfies the constraints.
5. Persist the measurement to `.papipeline/benchmarks.json`.

The recommendation carries its `reason`, so a user can see
`"cmscan does not scale linearly; 4x2 measured fastest"`. If the machine
changes, the next benchmark re-measures; nothing is hard-coded.

## 10. What is deliberately *not* optimised

| Not done | Why |
|---|---|
| Disabling Bakta's `--addons` / expert system | AMR and virulence annotation depend on it. Turning it off would silently remove stage 3 and stage 8 inputs to save time. |
| A shared on-disk cache of Bakta results across projects | Cross-project cache reuse needs a content-addressed key including the database version. Worth it at scale, but it is a later milestone with its own invalidation rules. |
| Running genomes in containers per job | Container start cost exceeds the per-genome saving at 353 s/genome. Revisit if the median drops below ~60 s. |
| Parallelising other stages across samples | They are single-process over the cohort already; intra-process parallelism is the right shape, and their combined cost is ~30 min for 200 genomes. |
