# 12: snakemake dry run with both machine configs

**What to build:** Proof that the DAG resolves under both the laptop and the big-machine configuration, which is what catches a path, a missing input or a config key that only works on one machine.

**Blocked by:** 11

**Status:** ready-for-agent

- [ ] snakemake --dry-run completes cleanly under the laptop config.
- [ ] snakemake --dry-run completes cleanly under the big-machine config.
- [ ] The 20-sample cap fires correctly under the laptop config during a dry run.
- [ ] The laptop config resolves to no path that exists only on the other machine, and vice versa.

## BLOCKER found in Phase 2: the Snakefile has never loaded

`snakemake --dry-run` on the current workflow/Snakefile fails immediately with:

```
SyntaxError in file "workflow/Snakefile", line 152:
Multiple run/shell/script/notebook/wrapper/template_engine/cwl keywords in
rule provenance.
```

**All 19 rules declare both `script:` and `shell:`.** Snakemake allows exactly
one execution keyword per rule, so the file is not loadable. This predates the
imipenem work: it is present at commit `8b35f4e` and at the Phase 1 checkpoint,
and the Phase 2 diff to the Snakefile is configuration-only.

Consequence: the Snakefile has never been executed even once, so
"the DAG resolves" has never been demonstrated for any configuration. Ticket 12
cannot pass until this is fixed, and the fix is a design question rather than a
typo: either the rules keep `shell:` (and drop `script:`) or they become true
Snakemake `script:` rules. Each has consequences for how `threads` and `mem_mb`
are declared, which is the same work as ticket 11. Fix them together.
