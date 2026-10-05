# 14: Variant calling against a verified PAO1 reference

**What to build:** SNP and indel calls for every isolate against a pinned PAO1 reference, restricted to the loci that matter for carbapenem resistance. The reference is verified against the actual file, because a truncated or substituted reference shifts every variant coordinate and stays invisible until the biology is wrong.

**Blocked by:** 11, 08

**Status:** ready-for-agent

- [ ] The reference's actual length and checksum are verified against the pinned values at run time, and a mismatch fails the run.
- [ ] Nothing is downloaded during a run; provisioning is a separate deliberate step.
- [ ] Calls are restricted to oprD and the pinned regulators (mexR, nalC, nalD, mexZ, ampD, ampR) and the restriction is enforced, not merely intended.
- [ ] PAO1 is pinned as GCF_000006765.1 with its accession and checksum recorded in the references table.
- [ ] Call output is in the format the GWAS consumes.
