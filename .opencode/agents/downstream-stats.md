---
description: Builds control gate, conditional scan, interaction tiers, lineage meta-analysis and evidence tiers
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

You own: papipeline/downstream/ (new package), config/interaction_pairs.tsv (new), config/known_determinants.tsv (new), tests/unit/test_downstream_*.py.
Tasks:
1. Positive-control gate: confirm the baseline scan recovers oprD burden, any_MBL and other configured controls. Fail loudly, naming the missing control.
2. Conditional scan: known layer features as covariates, plus a stratified scan in isolates lacking all known determinants.
3. Interaction tiers: tier 1 pre-specified pairs from interaction_pairs.tsv (for carbapenems, for example PDC-high with oprD_off_any, KPC with oprD_off_any, nalC with mexR, mexZ with ampD), tier 2 anchored scan with FDR. Logistic regression with an interaction term. Skip tier 3.
4. Lineage-stratified meta-analysis (fixed effect, Cochran Q, I squared).
5. Novelty filter against known_determinants.tsv (derive it from mechanisms.tsv) and evidence tiers A to D.
6. Use the repo claim levels DETECTED, PREDICTED, ASSOCIATED, SUPPORTED, UNKNOWN. Never add a CAUSAL level.
