# Design 5 — CLI design

Companion to [1 architecture](01-architecture.md).

## Command surface

The brief's required verbs, plus what a real operator needs. The grouping
is `noun verb` so that `papipeline project --help` reads as a noun and
completion stays predictable.

```
papipeline project create      papipeline run         papipeline tools
papipeline project list        papipeline status      papipeline report
papipeline project show        papipeline resume      papipeline benchmark
papipeline project validate    papipeline logs        papipeline stages
papipeline project delete      papipeline doctor
```

Ten top-level nouns, 16 verbs. `doctor` is added because a local
bioinformatics toolchain has many ways to be half-installed, and the first
question is almost always "why did it not run".

## 1. `project create`

```bash
papipeline project create \
  --name pa_imipenem_2026 \
  --genomes  /data/pa_assemblies/ \
  --metadata /data/pa_sample_metadata.tsv \
  --phenotype /data/pa_imipenem_phenotype.tsv \
  [--amr /data/precomputed_amr.tsv] \
  [--antibiotics imipenem meropenem] \
  [--source copy|symlink] \
  [--profile macbook|workstation|hpc] \
  [--root ~/projects] \
  [--dry-run] \
  [--yes]
```

`--dry-run` runs the full validation set and prints the report without
creating anything. This is the cheap way to answer "how many of my genomes
are usable" before committing to a long run — the question that matters
most when 200 FASTA files turn out to contain 60 corrupt ones.

Validation output is per-sample and never a single pass/fail:

```
$ papipeline project create --name pa_test --genomes /data/pa/ \
    --metadata md.tsv --phenotype ph.tsv --dry-run

Scanning /data/pa/ ...
  200 files found, 200 unique

Validating 200 assemblies ...
  ok            138   (4.0-8.0 Mb, decodes cleanly)
  undersized      9   below qc.min_assembly_size=4000000
  corrupt        41   FASTA decode error at byte 1048576
  empty          12   zero bytes
  WARNING: 62 of 200 assemblies are not analysable.
           They will be recorded as excluded, never replaced.

Metadata join
  200 metadata rows, 200 matched, 0 unmatched
  0 genomes without a metadata row

Phenotype
  imipenem  200 samples   R=163  S=29  I=8

  138 samples would be analysed.

dry run: nothing written. Re-run without --dry-run to create the project.
```

The excluded samples are **reported, not hidden**, and the run proceeds on
the analysable subset. This mirrors what the real cohort already does:
835 assemblies on disk yielded 65 analysable, each exclusion recorded with
a reason in `assembly_integrity.tsv`, and nothing was ever substituted.

## 2. `run`

```bash
papipeline run <project> \
  [--stages annotation amr gwas] \
  [--skip convergence cooccurrence] \
  [--profile macbook|workstation|hpc] \
  [--workers 3] [--threads 2] [--memory 4096] [--tmpdir /scratch] \
  [--annotation-only] \
  [--dry-run] [--resume] [--foreground]
```

Behaviour:

- **Default is background.** Prints the run id and returns, so the shell is
  not held by a ten-hour job.
- `--dry-run` prints the DAG, the resolved resources and the preflight
  verdict, and writes nothing.
- `--resume` is the same as a fresh `run` in intent but continues the most
  recent incomplete run for the project. Resuming is the default behaviour
  of a re-run over an existing project directory; `--resume` just makes it
  explicit and prints what will be skipped.
- `--foreground` streams logs to the terminal, for CI and for debugging.

`run` refuses to start — with a `doctor`-style explanation, not a stack
trace — when:

- a required tool is missing or below its version floor;
- the Bakta database is absent or incomplete;
- free disk is below the measured requirement;
- a previous run for this project is still active.

## 3. `status`

```bash
papipeline status [<run_id>] [--json] [--watch] [--stage <stage>]
```

Default is the latest run of the current project. `--watch` re-reads
`.papipeline/state.json` (see design 3 §5) rather than polling SQLite, so
it stays cheap and keeps working if the database is locked by a writer.

```
$ papipeline status

project  pa_imipenem_2026
run      run-20260927T112233Z          RUNNING   elapsed 01:47:12
mode     REAL                          profile  macbook (8 workers x 2 threads)

  #  stage                status    elapsed   detail
  1  validation           COMPLETE      0:04   200/200 assemblies
  2  annotation           RUNNING      38:19   61/200 genomes, 2 workers
  3  mlst                 PENDING        -     needs annotation
  ...
 15  integration           PENDING        -     needs amr, mechanisms, regulators
 16  reporting             PENDING        -     needs validation, phenotype

  annotation: 61/200 complete, 0 failed, 139 remaining
  ETA annotation: 5h 12m (measured 348 s/genome, 8x2)
  log: results/../logs/stages/annotation.log
```

The ETA comes from measured throughput on this machine, never a constant.

## 4. `resume`

```bash
papipeline resume <run_id> [--dry-run] [--force-rerun <stage>]
```

Prints the plan before doing anything:

```
$ papipeline resume run-20260927T112233Z --dry-run

  COMPLETE  validation       200 outputs verified present
  COMPLETE  annotation       61/200 genomes  -> 61 kept, 139 will run
  FAILED    mlst             attempt 2/3, error: mlst: command not found
  PENDING   amr              will run
  SKIPPED   convergence      config-disabled

  140 stage jobs will run, 61 genome jobs will be skipped.
  note: annotation was interrupted; its 61 completed genomes are reused
        from the checkpoint, not re-annotated.

resume? [y/N]
```

`--force-rerun <stage>` invalidates one stage and everything downstream of
it. It never silently invalidates upstream work.

## 5. `tools`

```bash
papipeline tools                 # the registry
papipeline tools check [name]    # preflight, GO/NO-GO
papipeline tools databases       # database versions and pin status
papipeline tools provision bakta_db --type light   # explicit, opt-in
```

```
$ papipeline tools

  tool          version   floor   stage          source     status
  bakta         1.12.1    1.12    annotation     env        ok
  amrfinder     4.2.7     4.0     amr            env        ok
  mlst          2.33.1    2.19    mlst           env        ok
  pyseer        1.1.2     1.1     gwas           env        ok
  minimap2      2.31      2.17    phylogeny      homebrew   ok
  samtools      1.24      1.10    phylogeny      homebrew   ok
  bcftools      1.24      1.10    phylogeny      homebrew   ok
  FastTree      2.2.0     2.1     phylogeny      env        ok
  panaroo       -         -       pangenome      -          unavailable
  snip-sites    -         -       phylogeny      -          unavailable

  8 of 10 available. 2 stages will use a documented fallback.
```

Version resolution is by **version, not PATH order** — the reason is
concrete: this machine's env ships `samtools` and `bcftools` 0.1.19 from
2014, which shadow the working Homebrew 1.24 builds and fail *silently* on
modern SAM. The registry checks each candidate and reports which it
skipped, and `run` refuses to schedule a stage whose tool resolves to a
version below the floor.

## 6. `benchmark`

```bash
papipeline benchmark <project> [--n 5] [--profile ...]
```

Runs annotation on `--n` representative genomes, measures, and recommends.
This is the brief's "recommend based on measured throughput" requirement,
and it is a real measurement rather than a table.

```
$ papipeline benchmark pa_test --n 5

  5 genomes, 8 cores, 17.2 GB RAM, SSD

  configuration        throughput      wall for 139 remaining
  workers=1 threads=8    9.1 g/h        15h 16m
  workers=2 threads=4   17.4 g/h         7h 59m
  workers=4 threads=2   19.8 g/h         7h 01m   <- recommended
  workers=8 threads=1   14.2 g/h         9h 47m

  recommended: --workers 4 --threads 2
  reason: cmscan does not scale linearly with threads, so more jobs with
          fewer threads each beats fewer jobs with more threads. 4x2 was
          the best measured point and leaves 0 cores for the OS.
  projected: annotation 7h 01m, whole pipeline ~7h 40m
```

The recommendation is a recommendation. It is recorded in
`.papipeline/resources.json` and shown in the GUI's Resource Configuration
screen, and the user can always override.

## 7. `report`

```bash
papipeline report <run_id> --kind full_analysis|amr|gwas|antibiotic_comparison \
                  [--format md|html|both] [--antibiotics imipenem meropenem]
```

## 8. `stages` and `logs`

```bash
papipeline stages                     # the 16 stages, order, prerequisites, tools
papipeline logs [<run_id>] [--stage <s>] [--sample <id>] [--tail N] [--follow]
```

`--sample` tails one genome's Bakta log, which is the level of detail
needed when a single genome in a 200-genome sweep is the thing that failed.

## 9. Global flags

```
--home PATH        papipeline state root (default ~/.local/share/papipeline)
--project PATH     current project, so most verbs need no name
--json             machine-readable output on every verb
--log-level LEVEL  default INFO
--no-color         for CI logs
--version
```

Every verb supports `--json`. That is what lets the GUI reuse the same code
paths the CLI uses, and it makes each verb testable without parsing
human-formatted text.

## 10. Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | the verb ran and the *analysis* reported failure |
| 2 | usage error — bad flags |
| 3 | validation failed — bad inputs |
| 4 | preflight failed — missing tool or database |
| 5 | resource insufficient — disk, memory, cores |
| 130 | interrupted |

Distinct codes matter for CI: "the pipeline failed" and "the tool is not
installed" are different problems and a script should be able to tell them
apart without scraping text.
