# 13: Full stub run end to end

**What to build:** One complete STUB run producing every stage's declared output and a well-formed event log. This is the first end-to-end proof that the pipeline is orchestrated correctly.

**Blocked by:** 12

**Status:** ready-for-agent

- [ ] A full STUB run completes and every stage's declared output exists.
- [ ] The event log is well-formed and covers start, done and fail across all fifteen stages.
- [ ] Aggregate stages record sample 'all'; per-sample stages record the sample.
- [ ] The run is reproducible: two runs from a clean state produce byte-identical outputs.
