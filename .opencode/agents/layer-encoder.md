---
description: Builds the resistance layer feature encoding and the layer TSV writers
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

You own: papipeline/layers/ (new package) and tests/unit/test_layers_*.py.
Tasks (ticket 15 plus layer encoding):
1. Feature builders with layer-prefixed columns: L1_ acquired (any_MBL for NDM/VIM/IMP, KPC, GES, OXA groups), L2_ PDC (allele plus position flags, ampD/ampR/dacB regulators, PDC-high composite), L3_ oprD (oprD_absent, oprD_LoF tier 1 for stop, frameshift, IS insertion, large deletion; tier 2 probable loss; composite oprD_off_any; other missense kept as individual features, never forced into the composite), L4_ efflux regulators (mexR, nalC, nalD, nfxB, mexZ loss and missense, pump-high composites), L5_ targets (gyrA, parC, ftsI, any_QRDR, QRDR count, using the PAO1 numbering AMRFinderPlus emits).
2. Unitig adapter around unitig-caller: confirm flags with --help, stub in TEST mode.
3. Writers: layer1_acquired.tsv, layer2_pdc.tsv, layer3_oprD.tsv, layer4_efflux.tsv, layer5_targets.tsv, all_layers.tsv, meropenem_feature_matrix.tsv (phenotype merged), feature_dictionary.tsv (feature, layer, source tokens, rule).
4. Drop features present in more than 98 percent or fewer than 5 isolates and log the drop list. Add partial_call_flag for PARTIAL, MISTRANSLATION and HMM calls. Assert every output has the same isolate set.
5. Tests with synthetic fixtures only.
