# Design 6 — Snakemake integration

Companion to [1 architecture](01-architecture.md).

## 1. Status: what exists, and why it cannot just be turned on

`workflow/Snakefile` is 500 lines with 20 rules and has **never been
executed** — `README.md:212` says so, and `environment.yml` pins
`snakemake-minimal=8.20.5` while the same README records that `snakemake`
was not present in the build environment. Its own header (`Snakefile:14-29`)
explains the omission:

> […] an untested caching layer would be a liability rather than an asset.

That judgement was right, and the reason it was untestable is the same
reason it is now fixable: the stages exchange **in-memory structures**, so
Snakemake's file-based dependency model cannot see between-stage edges. A
per-stage rule has to re-run its prerequisites in-process, which means a
full run recomputes shared prerequisites repeatedly.

`workflow/rules/` is empty. There are no `profiles/`. There are no
`resources:` directives and no `configfile:` — config is read by a direct
`load_config()` call at `Snakefile:65`, so Snakemake cannot see a config
change as a dependency. `THREADS` comes from `--config threads=`, and
`runtime.threads`/`runtime.memory_mb` in `config.yaml` are never read by
the workflow at all.

**So Snakemake becomes authoritative only after stage state is persisted
(design 8, M3).** Doing it earlier would put a scheduler in charge of
dependencies it cannot see.

## 2. Target shape

```
workflow/
├── Snakefile                     thin: imports rules/, declares configfile + profiles
├── rules/
│   ├── annotation.smk            per-genome Bakta, checkpointed, one shared DB
│   ├── stages.smk                one rule per scientific stage
│   ├── report.smk                the four report kinds
│   └── provenance.smk
├── profiles/
│   ├── macbook.yaml
│   ├── workstation.yaml
│   └── hpc.yaml
├── envs/
│   ├── constraints-osx-arm64.yml
│   └── constraints-linux.yml
└── config/
    └── project.schema.json       the config keys Snakemake reads
```

### `Snakefile` shrinks

```python
configfile: "config/project.schema.json"
include: "rules/annotation.smk"
include: "rules/stages.smk"
include: "rules/report.smk"
include: "rules/provenance.smk"
```

All path derivation and config reading moves into rules. Today the
preamble is 82 lines of Python; that belongs in a rule or in
`papipeline_app`, not in the Snakefile.

## 3. Per-stage rules, made genuinely incremental

The key change is that each stage rule declares **file inputs and file
outputs**, and each in-memory stage first persists its predecessors'
artefacts to the paths the engine already writes.

```python
rule amr:
    input:
        metadata   = f"{OUT}/intermediate/metadata.tsv",
        amr_in     = f"{OUT}/intermediate/amr/amr_determinants.tsv",
        config     = expand("{root}/config/{f}", root=config["project_root"],
                            f=CONFIG_FILES),
    output:
        amr        = f"{OUT}/intermediate/stages/04_amr.tsv",
    threads:    config["threads"]
    resources:
        mem_mb    = config["memory_mb"],
        runtime   = config["runtime_minutes"],
    log:        f"{OUT}/logs/stages/amr.log"
    shell:
        "papipeline-runner --stage amr --project {OUT} --log-file {log}"
```

Two fixes make this honest:

1. **`CONFIG_FILES` must include `gene_families.tsv`.** The current list at
   `Snakefile:84-93` omits it, so a change to the ESBL/MBL nomenclature
   would not invalidate the AMR stage. That is a real dependency-tracking
   bug today, independent of the app layer.
2. **Every rule declares every file it reads**, including the config
   knowledge tables, so Snakemake's DAG is complete.

### The prerequisite problem, solved properly

Today `PREREQUISITES` (`papipeline/run.py:139-148`) is a Python dict and
`run_stage.py` expands it in-process. That stays — it is the engine's
contract and it is correct. What changes is that the *Snakemake* DAG also
encodes the same edges, so Snakemake schedules them without re-running
anything:

```python
def prerequisite_files(stage):
    """The artefacts a stage reads, per the engine's PREREQUISITES map."""
    return {
        "mechanisms":   ["04_amr.tsv", "06_regulators.tsv", "07_structural_variants.tsv"],
        "convergence":  ["04_amr.tsv", "06_regulators.tsv", "10_phylogeny.tsv"],
        "gwas":         ["11_phenotype.tsv", "12_gwas_features.tsv"],
        ...
    }[stage]
```

The map is generated **from `papipeline.run.PREREQUISITES`**, not
hand-written, so the DAG and the engine cannot disagree. A unit test
asserts every stage in `STAGE_ORDER` has a rule and that each rule's
`input` files cover its engine-declared prerequisites.

That is the change that makes a full run compute each stage once.

## 4. Profiles

Profiles are the mechanism for "do not assume all machines have the same
CPU/RAM". They carry **Snakemake's** resources; the profile name is also
handed to `papipeline_app.resources` so the app can reason about the same
machine.

### `profiles/macbook.yaml`

```yaml
# Apple Silicon laptop. Measured on 8 cores / 17 GB.
cores: 8
memory_mb: 8192
threads_default: 2
max_threads: 8          # never use all cores; leave one for the OS
annotation:
  workers: 4
  threads_per_job: 2    # cmsscan scales sub-linearly; more jobs beats more threads
  tmpdir: null          # local SSD; do not force a tmpdir
  prune_outputs: true
allow_gpus: false
```

### `profiles/workstation.yaml`

```yaml
# Many-core desktop or a local server.
cores: 64
memory_mb: 131072
threads_default: 8
max_threads: 32
annotation:
  workers: 16
  threads_per_job: 4
  tmpdir: null
  prune_outputs: true
allow_gpus: false
```

### `profiles/hpc.yaml`

```yaml
# Scheduler-allocated. The cluster decides CPUs; Snakemake is told only
# the ceiling, and resources are a request, not a reservation.
cores: null              # inherit from the allocation
memory_mb: 1024          # per job, so a job cannot hog a node
threads_default: 4
max_threads: 16
annotation:
  workers: 8
  threads_per_job: 4
  tmpdir: "$TMPDIR"      # node-local scratch is faster than shared storage
  prune_outputs: true
allow_gpus: false
cluster:
  submit: "sbatch"
  partition: null
  account: null
```

### Selection and override order

Later wins:

1. profile defaults
2. `config/runtime.{threads,memory_mb}` from `config.yaml`
3. `--profile <name>` on the command line
4. explicit `--cores/--threads/--resources` flags
5. `benchmark` recommendations, applied only when the user accepts them

`hpc.yaml` is deliberately different in kind: `cores: null` and a per-job
memory *request*. A profile that hard-codes cores is wrong on a cluster,
which is exactly the assumption the brief warns against.

## 5. Annotation: the only genuinely parallel stage

Annotation is `N` independent genomes. Every other stage is one process
over the whole cohort, so `threads` on those rules is coarse parallelism
inside a single process, not fan-out.

```python
rule annotation_one:
    input:
        genome   = lambda w: GENOMES[w.sample],
        db       = rules.bakta_db.output.path,
    output:
        gff      = f"{OUT}/results/annotation/{{sample}}/*.gff3",
        tsv      = f"{OUT}/results/annotation/{{sample}}/*.tsv",
    threads:  config["annotation"]["threads_per_job"]
    resources:
        mem_mb = config["memory_mb"],
        runtime = 60,
    log:     f"{OUT}/logs/bakta/{{sample}}.log"
    shell:   "bakta-runner --sample {wildcards.sample} --threads {threads} ..."
```

Snakemake's own incomplete-file handling plus per-genome output files give
checkpointing for free: a re-run sees 61 completed genomes' outputs present
and schedules only the remaining 139. That is the mechanism the brief asks
for, and it is the same property the current `ThreadPoolExecutor` sweep
implements by hand — only now the scheduler owns it, so `papipeline resume`
and `snakemake --rerun-incomplete` agree.

The shared database is a single rule output, so every genome depends on one
database and it is never copied per genome.

## 6. Who calls whom

```
CLI ─┐
     ├─▶ papipeline_app.backend ─▶ papipeline_app.engine ─▶ papipeline stages
GUI ─┘            │                                            (subprocess)
                  └─▶ store (SQLite)
                  └─▶ snakemake --profile <p> --cores N <target>
```

- The app never re-implements scheduling. It shells out to
  `snakemake --dry-run` for the DAG preview and `snakemake` for execution.
- A snakemake `run:` directive calls `papipeline_app.backend` for status
  bookkeeping, so SQLite stays authoritative for the UI even when the
  workflow was launched directly from the command line. A run started with
  plain `snakemake` therefore still appears in the GUI.
- If Snakemake is absent, `engine.py` falls back to calling
  `scripts/common/run_stage.py` directly, in `EXECUTION_ORDER`. The results
  are identical; only concurrency and resume granularity differ. The GUI
  reports which path was used.

## 7. Execution-order correctness

`EXECUTION_ORDER` (`papipeline/run.py:119-136`) differs from `STAGE_ORDER`:
`mechanisms` runs *after* `regulators` and `structural_variants` because it
consumes them. The Snakemake rules must follow `EXECUTION_ORDER`, and a
test asserts the two never diverge — otherwise the DAG will cheerfully
build a plan that cannot execute.

## 8. Verification before this is trusted

1. `snakemake --dry-run` on TEST fixtures — DAG builds, mode gate refuses REAL.
2. `snakemake -j1` full TEST run; outputs compared byte-for-byte against
   `run_pipeline.py --mode TEST`.
3. Kill mid-annotation, re-run, assert completed genomes are not re-run
   (`--dry-run` lists only pending ones).
4. `snakemake --profile hpc --dry-run` with `cores: null`; assert no rule
   over-requests.
5. Touch `config/gene_families.tsv`; assert the AMR and mechanism rules are
   invalidated. This is the test that catches the current
   `CONFIG_FILES` omission.
