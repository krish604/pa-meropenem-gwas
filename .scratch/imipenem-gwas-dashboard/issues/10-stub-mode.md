# 10: STUB mode

**What to build:** A mode in which every rule emits a tiny, plausible-shaped fake output, so the entire DAG, the event stream and the dashboard can be exercised in seconds on a machine where no bioinformatics tool is installed and no fixture is on disk.

**Blocked by:** 03, 09

**Status:** ready-for-agent

- [ ] A configuration flag selects STUB mode, and it is the default for development and CI.
- [ ] All fifteen rules produce their declared outputs without invoking any real tool.
- [ ] A full STUB run completes with no fixture present and no tool installed.
- [ ] The event log from a STUB run covers all three event values and every stage.
- [ ] STUB and TEST remain distinct: STUB bypasses the real parsers, TEST exercises them.
