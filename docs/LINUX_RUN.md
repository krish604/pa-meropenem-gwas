# Running on a Linux machine

Ordered steps for the analysis machine: solve the environment, build it, record
what is in it, annotate, then — only when you have said so — run.

Nothing here has been executed on Linux. The scripts were written and unit-tested
on macOS, the environment was solved with `--dry-run --platform linux-64`, and
`sbatch` was not present on the writing machine. What each step claims is stated
next to it.

## Step 0 — what you need

* a Linux x86-64 machine (glibc ≥ 2.17), `bash`, and `micromamba` on `PATH`
  (or `MICROMAMBA=/path/to/micromamba`);
* this repository checked out;
* no pipeline data here yet: `data/`, `db/` and `PDC_essential.tsv` are
  provisioned by you, are never committed (AGENTS.md rule 4), and are read only
  from step 4 onwards.

## Step 1 — solve the environment (dry run, installs nothing)

```bash
micromamba create --dry-run --platform linux-64 -f environment/environment-linux.yml
```

**Current verdict: this FAILS, on purpose, as the file is pinned.** The failure
is recorded, not worked around:

| Blocker | Measured |
|---|---|
| `gubbins=3.4.1` | no linux-64 build exists for Python 3.11.16 — every published build is Python 3.10, and under `python=3.11.16` even a bare `gubbins` does not resolve |
| `samtools=0.1.19` | requires `openssl >=1.1.0,<=1.1.1`; Python 3.11.16 brings `openssl >=3.5.7`, so the two cannot coexist |

`mlst=2.33.1` and `panaroo=1.7.0` are also reported as blockers of the full
file; which of their own dependencies conflicts was not isolated, so no claim
is made about them.

Evidence (all from this repository, nothing reconstructed):

* `.build/linux-ready.solve.linux-64.fresh-full.txt` — the three top-level blockers
* `.build/linux-ready.solve.linux-64.core-nopy3tools.txt` — the same file minus
  `gubbins`/`mlst`/`samtools`/`panaroo` solves cleanly: **357 packages, 698 MB**
* `.build/linux-ready.probe.gubbins-python311.txt`,
  `.build/linux-ready.probe.samtools019-py311.txt` — the two blockers isolated
  one at a time

Two consequences:

1. **No pin in `environment/environment-linux.yml` has been changed.** The pins
   are what makes the laptop result and the analysis result comparable, and
   relaxing one to make a solve go green would be inventing a version. Re-solving
   is a deliberate decision that has not been taken: pick a gubbins/Python
   strategy and a samtools strategy, change the file, re-run this step, and keep
   the new evidence beside the old.
2. **Solving from macOS needs two overrides**, because the solver resolves glibc
   from the host and otherwise reports every linux-only package as missing:

   ```bash
   CONDA_OVERRIDE_GLIBC=2.31 CONDA_OVERRIDE_LINUX=5.10 \
     micromamba create --dry-run --platform linux-64 -f environment/environment-linux.yml
   ```

   On the Linux host itself, neither override is needed.

## Step 2 — build the environment and record what is in it

```bash
scripts/linux/bootstrap_linux.sh                     # create + record
scripts/linux/bootstrap_linux.sh --dry-run           # solve only
scripts/linux/bootstrap_linux.sh --versions-out environment/versions.linux.txt
```

Creates the environment named by the file itself (`pa-amr`) and then runs
`environment/versions.sh` **inside** it, so every tool and database version is
written down before anything is analysed. Skips creation if the environment
already exists. It never starts a run and never downloads a database.

If step 1 fails, this fails too — that is the same blocker, one step later.

## Step 3 — provision the databases, out of band

The pipeline never downloads a database, and a stage preflight refuses a release
that is not the one pinned in `config/references.tsv`:

```bash
bakta_db download --type light --output db/bakta_db
amrfinder --update
bash environment/versions.sh --out environment/versions.linux.txt   # record again
```

## Step 4 — inputs

1. `data/metadata/sample_metadata.tsv` — one row per assembly; only `sample_id`
   is required, `assembly_path` is optional but see step 5. Contract:
   `docs/data_contract.md`.
2. `data/phenotype/imipenem_phenotype.tsv` — `phenotype` is one of `S`, `I`,
   `R` (plus `SDD`/`ND` as the contract allows). **No MICs are derived, and no
   phenotype is inferred from one**; an absent phenotype stays missing and is
   never defaulted to `S`.
3. Every sample id must match across assemblies, annotations, variant tables and
   the phenotype file. A mismatch is a hard failure, never a silent drop.

Read `docs/reproducibility.md` before the first real cohort: no reference genome
is pinned yet, and that blocks variant calling in stages 6 and 7.

## Step 5 — annotate with Bakta (stage 2's input)

```bash
scripts/linux/run_bakta.sh --manifest data/metadata/sample_metadata.tsv --dry-run
scripts/linux/run_bakta.sh --manifest data/metadata/sample_metadata.tsv
```

Writes where `README.md` and `docs/INSTALL.md` say annotation lives:

```
results/real/intermediate/bakta/<sample_id>/<assembly stem>.gff3
                                                          .tsv
                                                          .inference.tsv
```

That root is read from `config/machines/linux.yaml` through the pipeline's own
configuration loader, so the batch and stage 2 cannot disagree about it. The
cohort, the assembly lookup, the file stem, the reuse verdict and the Bakta
command line all come from `papipeline` (`manifest`, `assemblies`,
`adapters.bakta`, `stages.annotation`) — the script adds no parsing of its own.

* **Database**: `annotation.bakta_db` (`db/bakta_db/db-light`), verified against
  the pin before anything is annotated. Override with `--db` or `BAKTA_DB`.
* **Threads**: the overlay's `runtime.threads`, or `--threads`.
* **Resume**: a sample whose existing output passes stage 2's own acceptance
  check is skipped; a sample with no located assembly is reported by name and
  the batch exits non-zero. It never silently drops a genome.
* **A stem trap, refused by name**: Bakta names its files after the assembly
  *file*, while stage 2 looks for `<genome_stem_for(sample)>.gff3`. Those agree
  only when the manifest says where the file is. If they do not, the sample is
  refused with the `assembly_path=...` line to add — annotating anyway would
  write files stage 2 could never reuse.
* **Reuse is off by default.** With `annotation.reuse_tool_output: "off"` in
  `config/science.yaml`, stage 2 runs Bakta for every genome regardless of what
  is on disk, so this batch would be redone. The script says so loudly. Set
  `prefer` or `require` there to reuse what this step produced (that file is
  not part of this work package).
* **SLURM, optional and UNVERIFIED here**: `--slurm --slurm-time HH:MM:SS
  --slurm-mem-gb N` writes one job per pending sample under
  `<out-root>/slurm/` and submits it. `sbatch` was not on the writing machine,
  so this path has never been executed; it fails with a clear message if
  `sbatch` is absent rather than pretending to submit. The directives used are
  the ones `scripts/hpc/submit_slurm.sh` already uses in this repository.
  Without `--slurm`, annotation runs serially in the calling shell.

## Step 6 — authorise and run

Two gates, both required: your word, and the pipeline's own check.

```bash
CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh --dry-run     # plan first
CONFIRM_REAL=yes scripts/linux/run_meropenem_real.sh               # the run
```

Anything other than `CONFIRM_REAL=yes` exits 3 before reading, executing or
writing anything (AGENTS.md rule 6). Only then does the script export the
session override that opens REAL mode, and run the repository's workflow as:

```
snakemake --snakefile workflow/Snakefile \
    --config mode=REAL machine=config/machines/linux.yaml \
    --cores <overlay runtime.threads> --rerun-incomplete --keep-going full_run
```

`--cores N` overrides the overlay's thread count. The overlay keeps
`runtime.allow_real_mode: false` permanently — the override is per session and
never written to a tracked file.

**What a real run is allowed to claim.** `docs/STATUS.md` is the record: eight
of sixteen stages have run on real data, on a 10-isolate smoke cohort, and
nothing there is a biological finding. Stages 7 and 8 were blocked by tooling on
macOS; stages 12, 13 and 14 (`gwas`, `convergence`, `cooccurrence`) refuse REAL
outright — see `papipeline.run.REAL_REFUSING_STAGES`. Read it before you read
any output.

## Step 7 — verify before believing anything

```bash
bash environment/versions.sh --out environment/versions.after.txt
diff environment/versions.linux.txt environment/versions.after.txt
ls results/real/intermediate/stages/        # one table per stage
cat results/real/run_manifest.json          # written only by a completed run
```

A run that fails before the manifest step writes no `run_manifest.json`; the
per-stage log is then the only record. Compare the two version files: a run
whose provenance cannot be reproduced is not a result.

## Also see

* `docs/PIPELINE_FLOWCHART.md` — the stage order, drawn
* `docs/data_contract.md` — input file formats
* `docs/scientific_rules.md` — what each stage is allowed to conclude
* `docs/environment-arm64.md` — why the laptop differs (gubbins, panaroo)
