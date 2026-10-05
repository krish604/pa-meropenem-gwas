# 06: Directional manifest to phenotype join with exclusions

**What to build:** A join rule that is deliberately asymmetric, because the manifest is the authority on what exists while the phenotype table is external evidence about it and may legitimately describe isolates that were never downloaded.

**Blocked by:** 05

**Status:** ready-for-agent | built. Directional and asymmetric. A manifest genome with zero or
two phenotype rows is a hard failure; a phenotype row with no assembly is
excluded; I with a measured MIC is KEPT (continuous trait), SDD/ND are
excluded regardless. Exclusion table + counts. 17 tests.

- [ ] Every manifest genome must have exactly one phenotype row. Zero rows or two rows for the same sample is a hard failure naming the sample.
- [ ] A phenotype row with no corresponding assembly is excluded and reported, never failed.
- [ ] An assembly whose phenotype is I, SDD, ND or absent is excluded, never failed, and is written to an exclusion table with a reason.
- [ ] Exclusion counts are logged and reported at the end of the run.
- [ ] Inside the pipeline, any sample-ID mismatch between assemblies, annotations, variant tables and the final phenotype table is a hard failure. There is no fuzzy or prefix matching anywhere.

**Build state: `built`.** Phase 2. Directional and asymmetric, with an exclusion table and counts. **`I` and `SDD` both keep a measured MIC** - SDD is a determination, and under EUCAST it is derived from an MIC inside the susceptible range, so excluding it (as an earlier revision did) discarded real measurements. Only `ND` is excluded regardless. **`papipeline.join` has no caller yet**: stage 11 does not use it, so this is built but not wired.
