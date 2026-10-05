# AGENTS.md

Guidance for coding agents working in this repository.

## What this is

A *Pseudomonas aeruginosa* imipenem antimicrobial-resistance GWAS pipeline, plus
a live monitoring dashboard. Input is **assembled genomes (FASTA)**. Read QC and
assembly are out of scope — the pipeline never assembles or cleans reads.

## Hard rules

These are not preferences. Violating any of them makes the work wrong.

1. **Never invent package versions, tool flags, or database names.** Verify with
   `conda search` / `micromamba search` and `<tool> --help`. If you are unsure,
   say so explicitly rather than guessing. A wrong flag is worse than no code.
2. **Never hard-code anything environmental.** Paths, thread counts, memory
   limits, database locations, and breakpoints are read from
   `config/laptop.yaml` or `config/bigmachine.yaml`. If you find yourself typing a
   path, a thread count, or an `os.cpu_count()` into code, stop.
3. **Never edit a test to make it pass.** A failing test is information. Use
   `/diagnosing-bugs` to find the cause.
4. **Databases and genome files are never committed.** See `.gitignore`.
5. **Sample IDs must match** across assemblies, annotations, variant tables and
   the phenotype file. Any mismatch is a hard failure with a clear message, never
   a silent drop or a fuzzy join.
6. **Do not start a real-genome run** until the user has explicitly said
   `run real samples`. `runtime.allow_real_mode` guards this in config; respect it.
7. **The laptop config must refuse to run more than 20 samples**, and must do so
   with a message that explains the limit and what to change.

## Modes

- **STUB** — every rule produces tiny fake outputs: each stage's
  contract-declared header, zero rows. One real `snakemake` run exercises the
  whole DAG and writes the event stream, with no bioinformatics tool installed
  and no fixture on disk. This is the default mode for development and CI.
  **The dashboard does not read this stream at all.** `scripts/emit.py` writes
  `status/events.jsonl`, and nothing in `papipeline/` reads that file — the
  observatory serves `/api/events` from an in-process `EventBus` instead. There
  are two observability mechanisms and they do not meet. A STUB run therefore
  proves the event *format*, and proves the DAG; it does not exercise the
  dashboard. Do not claim otherwise until one of the two is retired or bridged.
- **TEST** — the 20 committed synthetic fixtures in `test_data/`, byte-stable.
- **REAL** — actual genomes. Gated; see rule 6.

## Agent skills

### Issue tracker

Issues and specs are markdown files under `.scratch/imipenem-gwas-dashboard/`.
There is no git remote and no `gh`/`glab` CLI, so local markdown is the tracker of
record. `TASKS.md` is a roll-up of it, not the tracker. See
`docs/agents/issue-tracker.md`.

### Triage labels

Default vocabulary, recorded on a `Status:` line in each issue file:
`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`. See
`docs/agents/triage-labels.md`.

### Domain docs

Single-context layout. `CONTEXT.md` and `docs/adr/` do not exist yet and are
created lazily by `/domain-modeling`. The **pre-existing** design docs in `docs/`
and `docs/design/` are authoritative — read them before changing behaviour, and
surface any contradiction rather than overriding silently. See
`docs/agents/domain.md`.

## Pre-flight reads

Before touching a stage, read the doc that owns it:

| Changing | Read first |
| --- | --- |
| Any stage's behaviour | `docs/scientific_rules.md` |
| A file format between stages | `docs/data_contract.md` |
| **What a stage IS** | **`.scratch/imipenem-gwas-dashboard/spec.md`** |
| Orchestration or structure | the spec's Implementation Decisions, then `docs/architecture.md` as history |
| Anything touching the fixtures | `docs/reproducibility.md` |
| Test strategy | `docs/design/09-test-strategy.md` |
| Environment or tool availability | `docs/environment-arm64.md` |

### Three stage taxonomies exist. The spec wins.

`docs/architecture.md` describes 16 stages, `workflow/Snakefile` defines 20
rules, and the spec defines **sixteen**. An agent that reads the first two and
not the third will build the wrong pipeline. The spec supersedes both; they are
retained as history, not as instructions.

Sixteen, not fifteen: `spec.md:363` records the amendment on 2026-09-29, which
added `cohort_variants` at `6a` after stage 6 was run against a real assembly.
The stage list is fifteen *numbered* rows plus `6a`, which is sixteen stages and
one rule each. `papipeline/run.py:126` carries the same correction.

Four modules named in the old taxonomy — `mechanisms`, `regulators`,
`structural_variants`, `integration` — are **not stages and have no stage
number**. `spec.md:351` folds each into the stage that owns it, and
`papipeline/run.py:150-152` says which: `structural_variants` runs inside
`amr` (stage 4), `mechanisms` inside `amr`, `regulators` alongside `amr`,
`integration` inside `reporting`. So a grep for "stage 7" that lands on
structural-variant code is finding a stale number, not a stage.

Stage 10 (`similarity`) emits a **patristic distance matrix in substitutions per
site**, summed along the stage-9 tree path with the LCA's own stem excluded. It
is not a SNP count and not a count of differences. The same tree expressed
against a different alignment length gives different numbers, so the unit has to
travel with the value; the stage writes it beside the matrix in
`similarity.units.json`. See `docs/scientific_rules.md` §11.

### Running things on this machine

`MAMBA_ROOT_PREFIX` is unset, so environments live under the micromamba binary's
own prefix and `micromamba run -n pa-amr ...` works with no setup. Two
consequences:

- A `pilot100` env exists under `~/micromamba_pilot/envs/` and is **unreachable**
  — it belongs to a root prefix that is not active. Do not try to use it.
- `micromamba run -n pa-amr pytest` **prepends the env's `bin` to `PATH`**, which
  defeats any test that shadows a tool by prepending a stub directory. Run
  `pytest` from an activated environment instead:
  `eval "$(micromamba shell hook -s bash)" && micromamba activate pa-amr`

`gubbins` cannot be installed here from conda, so stage 8 cannot run on the
laptop from the package; `scripts/build_gubbins_from_source.sh` builds the same
tag from source into `tools/gubbins/`, which is why the binary is not on `PATH`.
`panaroo` **can** be installed here, but not from conda: the bioconda recipe
carries a `prokka` dependency that cannot be satisfied on arm64, and
`panaroo/prokka.py` is a GFF3 *reader* that never shells out to prokka, so
installing from source is correct rather than a workaround.
`pip install git+https://github.com/gtonkinhill/panaroo.git` into a venv has
been run end-to-end on real Bakta GFF3 for the 10-isolate smoke cohort, exit 0.
Stage 7 and stage 8 both have dispatch branches now, so neither tool's absence
is what stops a run — `runtime.allow_real_mode` is. See
`docs/environment-arm64.md`.

