# 08: Breakpoint layer, failing while unverified

**What to build:** A breakpoint standard held in configuration with its provenance, which refuses to interpret any MIC at all until real values are supplied from a licensed copy. CLSI M100 is paywalled and no free primary source publishes the P. aeruginosa imipenem S/I/R triple, so the honest state is unverified.

**Blocked by:** 07

**Status:** ready-for-agent

- [ ] The breakpoint file records standard, organism, agent, edition and table number.
- [ ] The threshold numbers are an explicit UNVERIFIED placeholder, not a guess.
- [ ] Phenotype interpretation hard-fails with a message naming the missing field while the values are unverified.
- [ ] MIC-to-R/S conversion runs only on rows carrying a genuinely measured MIC.
- [ ] A measured MIC the active standard does not cover is a hard error naming the row.
- [ ] I, SDD and ND are never folded into R or S.
