# 18: Dashboard backend: replay and stream

**What to build:** A backend that rebuilds the entire dashboard state from the event log when it starts, then streams new lines as they are written, reachable only over an SSH tunnel. It is a separate application from the existing observatory, which structurally cannot back it: that store holds one cohort-level row per stage, so it cannot produce a per-isolate grid.

**Blocked by:** 09

**Status:** ready-for-agent

- [ ] State is rebuilt by replaying the whole log on start, with no live stream required.
- [ ] New lines are then streamed over SSE as they are appended.
- [ ] The server binds to 127.0.0.1 only; 0.0.0.0 appears nowhere.
- [ ] A truncated final line left by a run killed mid-write does not break replay.
- [ ] A malformed line is skipped with a warning and the stream survives it.
- [ ] The backend does not import from, modify or share state with the existing observatory.
