---
description: Builds the kinship matrix and the real pyseer GWAS path
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

You own: papipeline/gwas_real/ (new package), papipeline/stages/similarity.py (new), tests/unit/test_gwas_real_*.py.
Tasks:
1. Stage 10 similarity: build the pyseer similarity (kinship) matrix from the IQ-TREE tree. Confirm the helper script names and flags from the installed pyseer package --help.
2. pyseer adapter for the real path: --lmm with --similarity, --kmers or --pres or --vcf, --covariates with --use-covariates, --lineage. Burden tests via --vcf and --burden. Binary phenotype.
3. Reuse the existing unique-pattern reduction (scripts/gwas/count_patterns.py) for the significance threshold. Do not re-implement it.
4. Keep REAL refusing unless the existing mode gates are opened. TEST mode must run on synthetic fixtures. The opt-in suite (PAPIPELINE_TEST_PYSEER=1) must pass if pyseer is installed.
