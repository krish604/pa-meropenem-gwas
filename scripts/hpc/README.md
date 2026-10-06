# *Pseudomonas aeruginosa* 900-isolate GWAS runner

Batch driver that runs an antimicrobial-resistance GWAS across 900 assembled
genomes on an HPC cluster. It splits the work into two phases and **stops to
ask your permission before the GWAS step**. All output is written to a single
`results/` folder. Nothing in the cohort is ever analysed twice.

---

## Requirements

Linux x86_64 cluster node · 64 GB RAM · 32 cores · [miniconda](https://docs.conda.io/)
· [Snakemake](https://snakemake.readthedocs.io/) ≥ 8. Tools (minimap2, samtools,
bcftools, AMRFinderPlus, Bakta, panaroo, gubbins, IQ-TREE, pyseer) come from the
pipeline's pinned environment file — do not install them by hand. Bakta is an
**input** here: existing Bakta GFF/TSV output is reused and never re-run.

---

## Setup (about 15 minutes, one-time)

This runner lives **inside** the pipeline repository, at `scripts/hpc/`. One
clone, no sibling checkout, and the runner can never drift from the workflow it
drives.

```bash
# 1. clone the pipeline -- the runner is already in it
git clone https://github.com/krish604/pa-amr-gwas-pipeline.git
cd pa-amr-gwas-pipeline

# 2. environment -- the pipeline pins every tool; do not install by hand
conda env create -f environment/environment-linux.yml
conda activate pa-amr

# 3. reference databases -- never downloaded by the pipeline; every stage
#    preflights them against config/references.tsv and refuses a wrong release
bakta_db download --type light --output db/bakta_db
amrfinder --update

# 4. one FASTA per isolate, named after its ID:
#    scripts/hpc/assemblies/PDT000034111.1.fna

# 5. cohort roster, one isolate ID per line
printf '%s\n' PDT000034111.1 PDT000034112.1 > scripts/hpc/config/samples.txt

# 6. tune if needed (PIPELINE_DIR defaults to this repository's root)
$EDITOR scripts/hpc/config/run.conf   # CORES, MEM_GB, EXPECTED_ISOLATES, MACHINE
```

The repository is **private** — clone it with `gh auth login` or a token on the
target machine.

`scripts/hpc/` is the study directory: its own `config/`, `assemblies/`,
`results/` and `docs/`. Nothing outside it is written by a run. To drive a
pipeline checked out somewhere else, set `PIPELINE_DIR` in
`scripts/hpc/config/run.conf` or in the environment.

REAL mode is `false` in every committed machine overlay, by design and by tests.
`run_all.sh` exports `PIPELINE_ALLOW_REAL_MODE=1` for the session instead, which
authorises the run **without editing a tracked file** — so the working tree stays
clean for the days a run takes. An explicit value you set yourself is respected,
never overwritten. `MACHINE` in `config/run.conf` picks the overlay: `bigmachine`
is uncapped, `laptop` refuses above 20 isolates.

---

## Run

```bash
scripts/hpc/run_all.sh
```

That single command:

1. **validates the cohort** — fails if any isolate ID appears twice, if the
   count differs from `EXPECTED_ISOLATES`, or if a listed assembly is missing
2. **runs Phase 1** — stages 1–11 (full list in `docs/STEPS.md`)
3. **stops and asks you**

```
 PHASE 1 COMPLETE.  GWAS HAS NOT BEEN APPROVED.
 ...
 Proceed to GWAS? [y/N]
```

4. **runs Phase 2 only if you approve** — stages 12–16: `pyseer` association,
   convergence, co-occurrence, report

Pressing Enter (or `n`) exits cleanly with Phase 1 results kept, so you can
inspect them and resume later with the same command. Approval is written to
`results/.gwas-consent`; `scripts/hpc/run_gwas.sh` refuses to run without it, so the
gate cannot be bypassed by calling that step directly. Restart from scratch
with `scripts/hpc/run_all.sh --restart`.

---

## On a cluster

Compute nodes have no terminal to prompt you, so approval happens on the
login node between jobs:

```bash
scripts/hpc/submit_slurm.sh phase1     # submit stages 1-11
scripts/hpc/submit_slurm.sh status     # watch it
scripts/hpc/submit_slurm.sh consent    # inspect results/, then approve
scripts/hpc/submit_slurm.sh gwas       # submit stages 12-16
```

Resources are read from `config/run.conf` (`CORES`, `MEM_GB`, `TIME_PHASE*`).

---

## Where the output goes

Everything is written under the study directory, never outside it:

```
scripts/hpc/results/
├── artifacts/stage_tables/    one table per stage
├── results/intermediate/      alignments, VCFs, pangenome graph
├── logs/                      per-stage logs + slurm-*.out
├── .phase1-complete           resume markers
└── .gwas-consent              your approval token
```

`scripts/hpc/results/` is git-ignored entirely (as is `assemblies/`, the
roster and the phenotype file — the repo ships the layout and the validation,
never the cohort). Re-running is safe: completed phases are skipped via the
marker files.

## Assumptions and limitations

- **No isolate reuse is enforced, not recommended.** `check_inputs.sh` aborts
  on a duplicate ID, because pooling the same organism twice silently biases
  allele frequencies and the case/control denominator.
- **The co-occurrence universe is capped** at 2,500 features (`max_features` in
  the pipeline's `config/science.yaml`); above that the stage refuses rather
  than runs — uncapped at n=900 it projects to 101 M pairs, ~54 h and 74 GB.
  So **mutation-inclusive** co-occurrence does not run at this scale;
  gene-vs-mechanism does.
- **A 900-isolate cohort is sized for discovery, not confirmation.** Treat
  p-values as candidates and correct for multiple testing before quoting a
  locus.
- **Stages 1–6 are validated on real data at n=10.** Stages 7–16 have never
  run end-to-end on real data — panaroo and gubbins have no usable macOS
  arm64 build. Expect stages 7–8 to be the wall-clock bottleneck; verify
  memory on a smaller subset first.
- **`max_depth` = 250** in variant calling is inherited, not chosen for this
  organism. Lowering it changes DP/AD everywhere, so it is a science decision.
- **Estimated wall time for 900 isolates:** ~1.5–3 days for Phase 1 on 32
  cores (stages 4 and 6 dominate; the pipeline loops isolates serially), plus
  Phase 2. Extrapolated from a measured 10-isolate run, not measured at 900.

## References

[Snakemake](https://snakemake.readthedocs.io/) ·
[Bakta](https://github.com/oschwengers/bakta) ·
[panaroo](https://github.com/gtonkinhill/panaroo) ·
[gubbins](https://github.com/andrewjpage/gubbins) ·
[pyseer](https://github.com/pyseer/pyseer) ·
[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinderPlus/) ·
[minimap2](https://github.com/lh3/minimap2) ·
[bcftools](https://samtools.github.io/bcftools/)
