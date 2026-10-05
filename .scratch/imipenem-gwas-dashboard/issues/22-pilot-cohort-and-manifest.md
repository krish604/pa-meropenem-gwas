# 22: Pilot cohort and REAL-mode manifest builder

**What to build:** The artefacts for a small, balanced pilot that actually exercises the GWAS, built without running anything. An arbitrary first-N-by-accession cohort would be all-susceptible and would not test the analysis at all.

**Blocked by:** 05, 06

**Status:** ready-for-agent

- [ ] A cohort file of 10-20 isolates, seeded and phenotype-stratified, records its own selection rationale inside the file.
- [ ] A manifest builder joins PDC metadata to the 835 downloaded assemblies, carrying isolate ID and BioSample as validated attributes.
- [ ] The builder produces artefacts only. Running against them stays gated behind the user's explicit instruction.
- [ ] The 20-sample laptop cap is exercised for real by this cohort size.
