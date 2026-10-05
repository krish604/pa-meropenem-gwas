# 05: Sample identity and manifest

**What to build:** A canonical sample identifier that is already the genome directory name and FASTA filename stem, so the assembly join needs no lookup table, with the PDC isolate and BioSample carried alongside as validated attributes.

**Blocked by:** 03

**Status:** ready-for-agent | wired. `sample_id` is the versioned GCA accession; the format
rule is `sample_id.pattern` in science.yaml, so real accessions and the 20
committed TEST_PA_* fixtures both validate; isolate/BioSample/isolate-name are
carried and required to be unambiguous. 26 tests.

- [ ] sample_id is the versioned GCA assembly accession.
- [ ] The ID format rule is configuration-driven, so both real GCA accessions and the synthetic TEST_PA_ identifiers validate while a malformed ID is still rejected.
- [ ] Isolate ID, BioSample and isolate name are carried as attributes and must map 1:1 onto the manifest or the run fails.
- [ ] Duplicate sample IDs and missing assembly files are hard failures.
- [ ] Contig headers are explicitly not treated as sample IDs, and no check ties them together.

**Build state: `verified`.** Phase 2. `sample_id` is the versioned GCA accession; the format rule is `sample_id.pattern` in science.yaml; isolate/BioSample/isolate-name are carried and must be unambiguous. **Now applied on the path a run actually takes** - `discover_manifest` applies the pattern and the uniqueness checks, which it previously did not, so the rules were inert in production.

**Promoted to `verified` at the start of Phase 3.** An end-to-end test
now drives `run_pipeline` on a cloned project and shows that a malformed
identifier, a duplicate identifier, and both directions of ambiguous
provenance all stop the run *and write nothing first*. The pattern and
provenance unit tests cannot show either half of that, because a unit
test on a function nothing reaches is not evidence of a feature.
