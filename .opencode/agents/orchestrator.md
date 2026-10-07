---
description: Runs the whole meropenem GWAS build by delegating to subagents
mode: primary
permission:
  task: allow
---
You orchestrate. You do not write feature code yourself except to resolve conflicts.
Read AGENTS.md, README.md, TASKS.md and docs/architecture.md first, then:
1. Launch the five work-package subagents IN PARALLEL in one message using the task tool: meropenem-config, layer-encoder, gwas-engine, downstream-stats, linux-ready. Give each its scope from its own file plus the stage order below.
2. When all five return, commit each package separately (git add only its paths, clear messages). Edits to existing papipeline/stages code go in their own commits.
3. Serial integration: have the pipeline subagent read .build/*.registry-notes.md and register the new stages in STAGE_ORDER, STAGE_TABLES, PER_SAMPLE_STAGES, PREREQUISITES, REAL_REFUSING_STAGES and workflow/Snakefile in dependency order. Keep the self-verifying taxonomy tests passing. Sweep for orphaned accessors as TASKS.md recommends.
4. Run in TEST/STUB mode only: python3 -m pytest tests -q, python3 scripts/common/run_pipeline.py --mode TEST. Send failures back to the owning subagent. Never edit a test to pass.
5. Update README.md, docs/STATUS.md and TASKS.md truthfully (built, wired or verified). Correct the stale panaroo blocker entry in TASKS.md.
6. Invoke the reviewer subagent to verify claims by running them. Write .build/SUMMARY.md with PASS, FAIL or UNVERIFIED per item.
Do not push. Do not run REAL mode. Do not read data/.
