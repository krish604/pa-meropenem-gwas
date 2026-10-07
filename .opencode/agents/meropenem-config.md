---
description: Adds meropenem to the config tables and builds the cohort size gate
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

You own: config/config.yaml, config/antibiotics.tsv, config/mechanisms.tsv, papipeline/cohort_gate.py (new), tests/unit/test_meropenem_*.py and tests/unit/test_cohort_gate.py.
Tasks:
1. Enable meropenem exactly as docs/architecture.md "Adding an antibiotic" describes: name it in config.yaml, add a row to antibiotics.tsv, extend the antibiotic column of the relevant mechanisms.tsv rows (comma list or all). The existing test test_enabling_meropenem_requires_no_code_change must still pass.
2. Build papipeline/cohort_gate.py: join assemblies with meropenem S/I/R using the existing papipeline.join, report n tested, n with assembly, n R, n I, n S. Make the intermediate-isolate policy a config key (default exclude). Refuse with a message naming the config key when resistant isolates fall below a configurable minimum (default 100).
3. Tests with synthetic fixtures only, including the refusal path.
