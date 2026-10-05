# 07: Phenotype parser with AST provenance

**What to build:** A phenotype reader that treats the source laboratory's R/S/I call as authoritative, preserves genuinely measured MICs without inventing them, and carries the AST method, standard and edition wherever the source provides them.

**Blocked by:** 03

**Status:** ready-for-agent | wired (parser) / built (model), with a change of approach. The model phenotype is now
log2(MIC), a continuous trait, not the binary R/S label. AST method/standard/
edition carried per row with missing-counts reported. Thresholds deliberately
empty: no numbers are written down, and a test asserts they stay that way.
17 tests.

- [ ] All five categories are accepted and normalised; an unknown category is rejected.
- [ ] A categorical call never gains an MIC or a zone diameter.
- [ ] A genuinely measured MIC and its unit are preserved exactly.
- [ ] A non-numeric MIC and an MIC on an ND row are rejected.
- [ ] AST method, standard and edition are carried per row where present, and the count of rows lacking each is reported.
- [ ] A sample with no phenotype record stays absent and is never defaulted to susceptible.

**Build state: `wired (parser) / built (model)`.** Phase 2. The parser is wired - AST provenance is carried per row with missing-counts reported, and `log2(MIC)` is computed on every record. Thresholds stay deliberately empty; a test asserts they stay that way. **The association model is not switched**: `gwas.outcome` is still the binary R-vs-S mapping, and `phenotype_matrix` has no production caller. That is ticket 16, and the spec says so explicitly rather than claiming it is done.
