# 11: 15-stage Snakefile with threads, memory and tool gating

**What to build:** The canonical fifteen-stage pipeline, one rule per stage, each declaring its own threads and memory, with the existing sixteen-stage code mapped onto it and with an immediate, clear failure for any stage needing a tool this machine does not have.

**Blocked by:** 10, 04, 05

**Status:** ready-for-agent

- [ ] Fifteen rules exist, one per stage, each declaring threads and memory read from configuration.
- [ ] The existing mechanisms, regulators, structural_variants and integration modules are internal steps within the rules that own them, not stages of their own.
- [ ] A rule requiring a tool marked unavailable on this machine fails at DAG-build time with a message naming the tool, not partway through a run.
- [ ] STUB and TEST still exercise the entire DAG regardless of tool availability.
- [ ] REAL mode is refused at DAG-build unless the user has authorised it.

## Tool availability to encode (from Phase 1)

Two tools cannot be installed in the laptop environment, for different reasons,
and the DAG-build gating must distinguish them in its message:

- `panaroo` - unsolvable on osx-arm64 at any version (prokka's dead
  `tbl2asn-forever` / `perl 5.26.2` pins).
- `gubbins` - HAS three osx-arm64 builds (3.4.3) but they are all py39/py310,
  and cannot coexist with Python 3.11.16. `gubbins=3.4.3` + `python=3.10`
  solves; + `python=3.11` does not.

Also relevant to stage 6: the installed `samtools` 0.1.19 has **no `consensus`
subcommand**, so variant calling must go through `bcftools mpileup` + `call`.
`samtools=1.8` is not an option because `mlst`'s Perl stack pins `samtools <0.2`.
