---
description: Prepares the Linux machine scripts, environment check and flowchart docs
mode: subagent
---
Hard rules (from AGENTS.md and spec.md, not optional):
- BUILD ON THE EXISTING REPO. Reuse papipeline/ modules, config/ tables, adapters and the Snakefile. Never rewrite working code.
- NO REAL-MODE RUN. Do not read data/, db/ or PDC_essential.tsv. Use TEST/STUB mode and the synthetic fixtures in test_data/ only.
- Tests first: write a failing test, watch it fail, then implement. Never edit a test to make it pass.
- Never invent a package version, tool flag or database name. Check with --help or say you are unsure.
- Nothing environmental is hard-coded: paths, threads and memory come from configuration.
- Never git checkout, git reset --hard or git push. Commit only your own files with git add <paths>.
- Stay inside the paths you own. Do NOT edit workflow/Snakefile or the STAGE_* registries. Write the registry changes you need to .build/<your-name>.registry-notes.md instead.
- Phenotype is binary S/I/R only. There are no MICs.
- Report the weakest state that is true: built, wired or verified. Never claim done for code nothing calls.
- A test you commit must collect. Tests-first means write the test before the code, but commit them together: run python3 -m pytest <your test files> --collect-only -q first, and never git add a test that cannot be imported.
- Be fast: one pass, no exploratory research beyond --help checks.

You own: environment/environment-linux.yml, config/machines/linux.yaml, scripts/linux/ (new), docs/LINUX_RUN.md, docs/PIPELINE_FLOWCHART.md.
Tasks (write scripts, never execute them on this machine):
1. Check environment-linux.yml with a dry-run solve for linux-64 (micromamba create --dry-run --platform linux-64). If a flag is unsupported say so and do not guess.
2. scripts/linux/bootstrap_linux.sh: create the environment and call environment/versions.sh. Never runs the pipeline.
3. scripts/linux/run_bakta.sh: batch Bakta (db-light) over the assemblies in the manifest, output to results/real/intermediate/bakta/<sample_id> as the README expects. SLURM optional.
4. scripts/linux/run_meropenem_real.sh: exits unless CONFIRM_REAL=yes is set; only then sets PIPELINE_ALLOW_REAL_MODE=1 and runs snakemake with machine=config/machines/linux.yaml.
5. docs/LINUX_RUN.md (ordered steps) and docs/PIPELINE_FLOWCHART.md (mermaid flowchart of the stage order).
6. Run shellcheck on the scripts if it is installed.
