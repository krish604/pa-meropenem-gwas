# Pseudomonas aeruginosa Imipenem AMR Pipeline

An assembly-based *Pseudomonas aeruginosa* antimicrobial resistance analysis
pipeline focused on **imipenem**, built to be extended to other antibiotics
through configuration alone. A read-only monitoring **dashboard** lives in
[`dashboard/`](dashboard/).

> **STATUS: WORK IN PROGRESS. Not a complete 16-stage result.**
>
> All 16 stages are implemented. They are exercised against **synthetic fixtures**
> in `test_data/`, and a REAL run has been done on a **10-isolate smoke subset**.
> That run reached **stage 6 only**. Stages 7–16 have never run on real data.
>
> **n = 10 is grossly underpowered. Nothing this repository produces is a finding
> about imipenem susceptibility in *P. aeruginosa*.** See
> [`docs/STATUS.md`](docs/STATUS.md).
>
> **Meropenem is now enabled as configuration** (2026-10-08):
> `config/science.yaml` `antibiotics:` lists `imipenem, meropenem`, the
> `config/antibiotics.tsv` meropenem row is enabled, and every
> `config/mechanisms.tsv` row carries `imipenem,meropenem`. Configuration only —
> **no analysis code changed**, and `project.primary_antibiotic` is **still
> `imipenem`**, so a meropenem run cannot find its phenotype file yet. See
> [`TASKS.md`](TASKS.md) Phase 8.

## What has actually run on real data

REAL-mode run, 10 isolates, `snakemake full_run`, one stage each:

| Stage | Outcome on real data |
|---|---|
| 1 genome validation | **ran** — 10/10 assemblies pass |
| 2 genome annotation | **ran** — 10/10 imported from verified Bakta output, **Bakta never executed** |
| 3 MLST | **ran** — 10/10 typed (*paeruginosa*) |
| 4 AMR detection | **ran** — 199 determinant rows, AMRFinderPlus 4.2.7 |
| 4a structural variants | **ran** — 409 candidate calls, labelled `candidate:` |
| 5 virulence | **ran** — VFDB 4573 factors |
| 6 variants + regulator screen | **ran** — 509,612 alleles; 100 locus-coverage values measured |
| 6a cohort variants | **ran** — 73,752 sites |
| 7 pangenome | **refused** — `panaroo` not importable and `cd-hit` not on `PATH` |
| 8 recombination | **not run** — `gubbins` cannot run on this platform |
| 9–16 phylogeny → report | **not run** — downstream of stages 7 and 8 |

Stage 7's refusal above is the round-12 record, not the current state of the
tool: `panaroo` has since been installed from source (pinned) with `cd-hit` from
conda and **verified end-to-end on real Bakta GFF3 for this same 10-isolate
cohort** — 10,019 gene families, 4,828 core — but **no REAL run has been
repeated since**, so stages 7–16 remain unrun on real data. Stage 8 is still
blocked on this machine: `gubbins` has no usable macOS arm64 build (see
[`docs/environment-arm64.md`](docs/environment-arm64.md) §3 — every conda build
is Python-3.10-only and the shipped ones crash with SIGSEGV), and stages 9–16
consume stage 8. Running all 16 needs a Linux/conda machine, and the Linux
environment file does not solve as pinned yet (`TASKS.md` Phase 8).

## Known issues in the current results

1. **`variants.tsv` holds 9 of 10 isolates.** One isolate's `bcftools mpileup`
   was SIGKILLed (`returncode=-9`) and produced no calls. The pipeline keeps it
   as a cohort member and says so, which is the defensible choice — but any
   variant-based statistic computed from these tables has **n = 9, not 10**. The
   OOM cause is **unproven**; it was not reproduced.
2. **`max_depth` is 250 and was deliberately not lowered.** It is the only lever
   that bounds mpileup RSS, but it subsamples reads, so lowering it moves
   `DP`/`AD`. That is a science decision and it is still **pending**. Output size
   was bounded instead (`-O z`), which shrank one call from 69,488,493 bytes to
   1,861,798 with a byte-identical variant set.
3. **Four reporting fixes are specified but not implemented** — carrying
   `call_status` beside the cohort denominator, among others. See `BACKLOG.md`
   in the round-12 artefacts.
4. **Tool versions are mostly `UNKNOWN` in the manifest** by design: probing is
   on-demand so the pipeline never executes a tool it does not need.

## Inputs that are NOT in this repository

Nothing here is clinical data. Supply it yourself, out of band:

| Input | What it is | Where it goes |
|---|---|---|
| `PDC_essential.tsv` | the PDC essential-genes table | worktree root; real clinical data, never committed |
| `data/` | genome assemblies (`GCA_*` / `PDT_*` FASTA) | `data/` |
| `db/` | AMRFinderPlus, VFDB, MLST, Bakta (`db-light`), PAO1 reference | `db/` |
| Bakta annotation | per-isolate `bakta` output trees | `results/<mode>/intermediate/bakta/<sample_id>` |

**Bakta is an INPUT to this pipeline and is never re-run.** Stage 2 imports
existing `bakta` output under `annotation.reuse_tool_output: require`, which
refuses loudly rather than executing the tool when curated output is missing.
Three guards enforce this in a REAL run: the config mode, executable tripwire
shims on `PATH`, and a background process watcher.

## License

**License: not yet chosen.** No license file is committed. Until one is, treat
the code as all rights reserved and do not redistribute it.

## Quick start

```bash
# 1. environment (optional for TEST mode; the core deps are already common)
micromamba env create -f environment/environment.yml
micromamba activate pa-amr

# 2. run the test suite
python3 -m pytest tests -q

# 3. run the pipeline against synthetic fixtures
python3 scripts/common/run_pipeline.py --mode TEST

# 4. see what it would do, without running anything
python3 scripts/common/run_pipeline.py --mode TEST --dry-run

# 5. run one stage
python3 scripts/common/run_stage.py --mode TEST --stage amr

# 6. start the read-only results dashboard (writes only ~/.pa_dashboard/)
scripts/run_dashboard.sh --results <results-dir> --port 8765
```

The run writes:

- `results/test/intermediate/stages/*.tsv` — the 16 stage outputs
- `results/test/run_manifest.json` — tool, version and database provenance
- `reports/test_pipeline_report.md` / `.html` — the test report

## Try to analyse real data

`REAL` mode is gated by `runtime.allow_real_mode` and requires an explicit opt-in,
so it is refused by default:

```console
$ python3 scripts/common/run_pipeline.py --mode REAL
ERROR: REAL mode is disabled: runtime.allow_real_mode is false in config.yaml.
       REAL mode is enabled only for the separate analysis phase, after this
       pipeline skeleton has been reviewed.
```

With a machine overlay that sets it, and with the inputs supplied, the production
path is Snakemake:

```bash
PIPELINE_ALLOW_REAL_MODE=1 snakemake --snakefile workflow/Snakefile \
  --config mode=REAL machine=config/machines/smoke.yaml --cores 4 full_run
```

Each stage runs exactly once. The laptop config refuses more than 20 samples.

## What it does

```
QC / sample metadata
      ↓
 1  Genome validation          assembly metrics, QC thresholds, tool hooks
 2  Genome annotation          Bakta parser → standardised records
 3  MLST                       sequence type + allele profile
 4  AMR detection              AMRFinderPlus adapter, CARD optional
 6  Regulator screen           oprD, mexR, nalC/D, mexZ, nfxB, mexT/S, ampD/R, dacB
 7  SV / MGE                   confirmed / candidate / not_assessable
 5  Mechanism interpretation   permeability, efflux, AmpC, acquired
 8  Virulence                  VFDB, independent of AMR
 9  Pan-genome                  core / accessory partition
10  Phylogenomics               core alignment → SNPs → tree
11  Phenotype interface        R / I / S / SDD / ND
12  GWAS                       pyseer adapter, BH correction
13  Convergence                independent lineages, not sample counts
14  Co-occurrence              association statistics only
15  Integration                the master genotype–mechanism–phenotype table
16  Figures + report           13 data-driven figure tables
```

Stages are numbered as in the specification and **executed** in dependency
order: stage 5 consumes stages 6 and 7.

## The one thing to understand

The pipeline distinguishes five claim levels and never exceeds what the
evidence supports:

| Status | Meaning |
|---|---|
| `DETECTED` | directly observed in this sample |
| `PREDICTED` | inferred without direct observation |
| `ASSOCIATED` | statistically associated in this cohort |
| `SUPPORTED` | associated **and** robust across lineages |
| `UNKNOWN` | nothing to report |

There is no `CAUSAL` level — not because it is discouraged, but because it
cannot be expressed in the output schema.

Concretely, and enforced in code:

- A detected `oprD` gene means an **intact locus**, not susceptibility. Only
  disruption or absence carries `reduced_permeability`.
- A resistant sample can have no detected determinant, and a susceptible sample
  can have one. The synthetic fixture is built to keep this true.
- A feature carried by one lineage is **lineage-linked**, not a resistance
  association, however small its p-value.
- A candidate structural variant is `PREDICTED` and labelled `candidate:` in the
  master table. It is never promoted to confirmed.
- `R`/`I`/`S` never acquire an MIC or a zone diameter.
- Missing values stay missing.

The ten rules and the code and tests that enforce each are in
[`docs/scientific_rules.md`](docs/scientific_rules.md).

## Layout

```
config/          config.yaml + 4 editable knowledge tables
  config.yaml              organism, antibiotics, thresholds, switches
  antibiotics.tsv          antibiotic classes and phenotype standards
  mechanisms.tsv           gene → mechanism, with a claim ceiling per gene
  regulators.tsv           loci screened by stage 6
  references.tsv           pinned tool and database versions
data/            REAL inputs (835 GCA_* assemblies, untouched by this build)
test_data/       TEST fixtures: 20 synthetic TEST_PA_* samples
workflow/        Snakefile — thin wrapper over the same stage functions
scripts/         16 per-stage CLI wrappers + 4 shared entry points
  hpc/               900-isolate HPC runner (phases, SLURM, consent gate)
                     merged in from pa-aeruginosa-900; see scripts/hpc/README.md
papipeline/      the library: all scientific logic
results/         run outputs, per mode
reports/         generated reports
tests/           3764 tests collected (see Tests)
environment/     environment.yml + versions.sh
docs/            architecture, data contract, scientific rules, reproducibility
```

## Configuration-driven

Three files, no code:

| To change | Edit |
|---|---|
| The antibiotic analysed | `config.yaml` `antibiotics:`, `antibiotics.tsv` |
| A resistance locus | `mechanisms.tsv`, `regulators.tsv` |
| A QC threshold | `config.yaml` `qc:` |
| Whether a stage runs | `config.yaml` `analysis:` |
| A tool or database version | `references.tsv` |

Adding meropenem is three edits: name it in `config.yaml`, add a row to
`antibiotics.tsv`, and extend the `antibiotic` column of the relevant
`mechanisms.tsv` rows (it accepts a comma-separated list, or `all`).

## External tools

None are needed for TEST mode. The pipeline detects what is present and
refuses a stage whose tool is missing rather than degrading silently. A tool
that is present but not needed is **not executed** to discover its version.

| Stage | Tool | osx-arm64 laptop | Notes |
|---|---|---|---|
| 1 validate | `seqkit` | yes | native arm64 build |
| 2 annotation | `bakta` | **input only** | never re-run; curated output is imported |
| 3 MLST | `mlst` | yes | |
| 4 AMR | `ncbi-amrfinderplus` 4.2.7 | yes | native arm64 build |
| 4a SV | `nucmer` | yes | candidate calls only |
| 5 virulence | `blast` 2.17.0 | yes | native arm64 build |
| 6 variants | `minimap2` 2.31, `bcftools` 1.23.1 | yes | native arm64 builds |
| 7 pangenome | `panaroo` + `cd-hit` | **yes, by source install** | TASKS BLOCKER 3: the conda recipe cannot solve on osx-arm64 (a spurious `prokka` dependency), so panaroo is installed from source via pip and pinned to a commit in `environment/environment.yml`, with `cd-hit=4.8.1` declared as a conda dependency. Verified end-to-end on real Bakta GFF3 for the 10-isolate smoke cohort: exit 0, 10,019 gene families, 4,828 core. Neither is on `PATH` in the `pa-amr` env as it stands — see environment-arm64.md §3 |
| 8 recombination | `gubbins` | **no** | no usable arm64 build — see environment-arm64.md §3 |
| 9 phylogeny | `snp-sites`, `iqtree` | yes | native arm64 builds |
| 12 GWAS | `pyseer` 1.1.2 | limited | opt-in suite; 21 tests skip unless enabled (measured 2026-10-08) |
| — | `snakemake` | yes | the production path for REAL runs |

`bash environment/versions.sh` prints the current status. Per-tool evidence is in
[`docs/environment-arm64.md`](docs/environment-arm64.md).

## Tests

```bash
python3 -m pytest tests dashboard/tests -q   # 3879 tests
python3 -m pytest tests/unit -q             # parsers and stage logic
python3 -m pytest tests/integration -q      # config, full run, mode gating
python3 -m pytest dashboard/tests -q        # dashboard server and sources
node --test dashboard/tests/js/             # dashboard page modules (no npm)
```

The dashboard suite is **112 passed, 3 skipped** by default and **114 passed,
1 skipped** with `PA_FIXTURES_LARGE=1` (the 900-isolate perf gates); `node --test`
adds 32. See [`dashboard/README.md`](dashboard/README.md) for the UI's own
security model and layout notes.

Latest full run, 2026-10-08 (re-run for `TASKS.md` Phase 8):

| Invocation | Result |
|---|---|
| `pytest tests --collect-only -q` | 3764 collected, **0 errors** |
| `pytest tests -q` | **33 failed, 3519 passed, 151 skipped, 61 errors** |
| `pytest dashboard/tests -q` | 112 passed, 3 skipped, 0 failed |
| `node --test dashboard/tests/js/` | 32 passed, 0 failed |

**The suite is not green and is not reported as green.** All 94 non-passing
node IDs are listed in `.build/BASELINE_PREEXISTING_FAILURES.txt`, recorded
from the untouched repository *before* the meropenem build: **0 new failures,
and 12 previously-failing STUB tests now pass** (the STUB data-root fix). The
failures are environmental — they need `db/` (`db/reference/…`,
`db/smoke_genomes`, the Bakta and AMRFinderPlus databases) or `data/`, and
databases and genomes are never committed, so they fail on any clean checkout.

The 151 skips are not all environmental either: 117 need uncommitted inputs
(`db/`, `PDC_essential.tsv`, the source-built `tools/gubbins`), 21 are the
opt-in `pyseer` binary suite (`PAPIPELINE_TEST_PYSEER=1` enables it), 12 are
empty-parametrize tests that assert nothing in a default run, and 1 needs
`python-pptx`.

They cover the parsers (FASTA, metadata, annotation, MLST, AMR, mutations),
sample-ID validation, mechanism mapping, phenotype validation, GWAS input
construction, tree/sample matching, integration and figure data preparation —
and the failure modes: missing values, duplicate sample IDs, invalid
phenotypes, unknown genes, unknown mechanisms, missing/duplicate tree samples,
invalid antibiotics and malformed TSV.

Three bugs the tests caught during the build, for the record:

1. `stages/cooccurrence.py` indexed feature sets by dict key instead of
   checking whether a sample carried the feature, so **every pair looked
   perfectly co-occurring**.
2. `viz.gene_by_phenotype_table` passed a gene→samples map where
   sample→genes was expected, so **figures 4 and 5 were all zeros**.
3. `stages/phylogeny._parse_label` absorbed stray characters into tip names,
   turning a malformed tree into a misleading "sample missing" error.

## Adding an antibiotic

See [`docs/architecture.md`](docs/architecture.md#adding-an-antibiotic). No
code change. Tests `test_enabling_meropenem_requires_no_code_change` and
`test_all_keyword_covers_any_antibiotic` cover the path.

## Reading the output

Start with the **master table** (`15_master_table.tsv`), one row per sample.
Its `confidence` column is derived, not inherited, and its `evidence_notes`
column states what was observed and what was not concluded.

`docs/data_contract.md` documents every column and every file.

## Not done yet

Honest list. Round 12 closed items 1, 2, 3, 6 and 7 below; what remains is here
with the reason.

1. **Stages 7–16 have never run on real data.** Stage 7's tooling has since
   become available (`panaroo` from a pinned source install, `cd-hit` from
   conda, verified end-to-end on real Bakta GFF3), but **no REAL run has been
   repeated since round 12**, so nothing has changed about what ran. Stage 8 is
   still blocked on this machine — `gubbins` installs and then segfaults — and
   stages 9–16 consume stage 8. Needs a Linux/conda machine whose environment
   file solves (`environment/environment-linux.yml` does not solve as pinned;
   see `TASKS.md` Phase 8).
2. **The GWAS reference engine has no kinship correction, and the real pyseer
   path is built but not wired.** `papipeline/gwas_real/` now holds a stage-10
   kinship matrix (Gower-centred patristic distances) and a two-pass pyseer
   adapter whose flags are each confirmed against the installed package —
    both **built and unit-tested, wired to nothing**. Three input resolutions
    are undecided, and a wrong `variants_path` would silently narrow which
    variant families are tested, so it is deliberately left unwired. A fourth
    precondition surfaced 2026-10-08: the `--distances` file this adapter
    would read from stage 10 is written **packed** by dispatch while the
    adapter parses **square** — the shape must be decided before wiring
    (issue 27, `.scratch/imipenem-gwas-dashboard/issues/`). Until it
    is wired, stage 12's `ReferenceEngine` remains the mechanics-testing
    stand-in, and says so on every result row. `pyseer` 1.1.2 is available but
    only its opt-in suite is exercised.
3. **No plotting layer.** All 13 figures have prepared, tested data tables
   written as JSON; the drawing layer is not implemented.
4. **`max_depth` is unreviewed.** 250 is inherited, not chosen; the decision to
   lower it is pending (see Known issues).
5. **`run_manifest.json` is not written when a run fails before the manifest
   step.** Stage status is then only in the per-stage log.
6. **Promoter-region assessment is disabled** pending a pinned reference.
7. **Tool provenance is thin.** `tools_detected` is mostly `UNKNOWN` by design;
   `config/references.tsv` is the only place versions live and it is largely
   unpinned.
8. **`project.primary_antibiotic` is still `imipenem`.** Meropenem is enabled
   as configuration (banner above), but the key decides which
   `<antibiotic>_phenotype.tsv` is loaded, so **a meropenem run cannot find its
   phenotype file yet**. Flipping it now was measured to break 96 tests
   (113 failed / 89 errors in a scratch flip, restored byte-exact), because the
   fixtures ship only `imipenem_phenotype.tsv`. It needs a meropenem phenotype
   fixture first — open work, not done.
9. **`unitig-caller` is absent and its flags were never verified with
   `--help`.** `papipeline/layers/unitigs.py` records that and REAL raises
   `UnverifiedFlagsError` rather than guess a flag, so the unitig screen is
   stubbed in TEST/STUB and deliberately produces nothing in REAL. Install the
   tool, run `--help`, record the real flags. The `oprD_absent` / `oprD_LoF`
   features are **not** what this adapter produces — they come from
   `regulators.oprd_feature_rows` via the stage-12 feature builder, and layer 3
   writes `L3_oprD_absent` / `L3_oprD_LoF_tier1` / `tier2`.

## Documentation

| Document | Contents |
|---|---|
| [`docs/STATUS.md`](docs/STATUS.md) | what has actually run on real data, and what has not |
| [`docs/INSTALL.md`](docs/INSTALL.md) | install on macOS arm64 or Linux, and the platform limits |
| [`docs/architecture.md`](docs/architecture.md) | layering, stage order, extension points, limitations |
| [`docs/data_contract.md`](docs/data_contract.md) | every file, column and validation rule |
| [`docs/scientific_rules.md`](docs/scientific_rules.md) | the 10 rules, the code enforcing each, the tests |
| [`docs/reproducibility.md`](docs/reproducibility.md) | environment, provenance, determinism, gaps |
| [`docs/environment-arm64.md`](docs/environment-arm64.md) | per-tool availability, and why `gubbins` cannot run on arm64 |
| [`dashboard/README.md`](dashboard/README.md) | the monitoring dashboard |
