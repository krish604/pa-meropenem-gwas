# 09: Event emitter

**What to build:** One script through which every rule announces that it started, succeeded or failed, writing a single append-only JSON Lines file that the dashboard will replay. This is the dashboard's only input, so it is specified first and blocks nothing.

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent | built. scripts/emit.py writes exactly {t, event, stage, sample}
with event in {start, done, fail}; aggregate rules use the literal 'all';
append-only; a truncated final line from a killed run does not corrupt the
next write. 15 tests.

- [ ] Exactly one JSON object per line with exactly the fields t, event, stage and sample.
- [ ] event is exactly one of start, done or fail; nothing else is ever written.
- [ ] Aggregate rules record the literal sample 'all' so event volume stays proportionate to the work.
- [ ] Writes are append-only and a re-run does not truncate.
- [ ] A malformed or truncated line is skipped with a warning rather than crashing the stream.
- [ ] Emitted identically by the Snakemake rules and by the native orchestrator, so the dashboard's source of truth does not depend on the entry point.

**Build state: `built`.** Phase 2. `scripts/emit.py` writes exactly four fields with exactly three event values, aggregate rules use the literal `all`, and the log is append-only and survivable across a killed run. **No rule calls it yet** - wiring the emitter into the rules is ticket 11, and the reader side is ticket 18.
