# 25: Bridge status/events.jsonl into the observatory EventBus

**What to build:** A reader that takes the `status/events.jsonl` the workflow
already writes and feeds it into `papipeline.observatory.events.EventBus`, so
the dashboard has one event stream instead of two that do not meet.

**Blocked by:** 13 (the STUB run that produces the event log)

**Status:** ready-for-agent

## Why this exists

There are two observability mechanisms and neither reads the other:

| | written by | read by |
|---|---|---|
| `status/events.jsonl` — one object per rule transition: `{"t", "event", "stage", "sample"}` | `scripts/emit.py`, called by `scripts/run_rule.py`, which every stage rule calls | **nothing.** No module in `papipeline/` opens this file. |
| `papipeline.observatory/events.py` — in-process `EventBus` plus a SQLite `ExecutionStore`; `/api/events` serves the bus | `papipeline/execution/runner.py`, in the same process as the run | the observatory API, and the dashboard if it uses them |

So a STUB run produces a well-formed event log that no consumer ever sees. The
dashboard can only be exercised against the in-process bus, which means only
from inside a run. That is why the spec's claim that STUB "exercises the
dashboard" could not be made true, and why `spec.md:62` and `:562` now carry
dated status notes.

## What to build

- A replay-then-tail reader for `events.jsonl`, so a dashboard started after a
  run still shows that run. Replay from offset, then follow.
- A mapping from the file's vocabulary (`start`/`done`/`fail`, `sample: "all"`)
  onto `ExecutionEvent` (`name`/`stage`/`subject`/`state`/`attempt`/`seq`).
  Decide explicitly what `"all"` becomes for `subject`; the file uses it for
  cohort-level rules and the bus expects a subject.
- A decision about the two mechanisms: bridge them, or retire one. Bridging
  keeps the file as the durable artefact and makes the bus a view over it. The
  file is the more robust of the two — it survives the process that wrote it.

## Definition of done

- [ ] A test feeds a real STUB `events.jsonl` to the dashboard's API and gets
      the run back. This is the check that does not exist today.
- [ ] A malformed line is skipped with a warning, not fatal. A run that dies
      mid-write is still replayable to its last complete line.
- [ ] Replaying a file twice does not double-count; the bus is idempotent over
      a replay, or the reader dedupes.
- [ ] `spec.md:62` and `:562` lose their status notes, or the notes are updated
      to say what is still missing.
- [ ] `.opencode/agents/dashboard-ui.md` drops the "one source of truth"
      ambiguity, whichever way the decision goes.

## Not this ticket

Do not build the dashboard. Do not add a fifth field or a fourth `event` value
to the file — the format is four fields and three values on purpose, and a
consumer that needs more justifies a format change, not an ad-hoc extension.
