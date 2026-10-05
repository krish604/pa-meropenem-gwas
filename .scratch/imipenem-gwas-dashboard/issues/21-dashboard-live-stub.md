# 21: Dashboard shown live against a stub run

**What to build:** The demo: the dashboard tracking a running STUB run, stages lighting up and isolates filling in, which is what proves the event contract and the UI actually agree.

**Blocked by:** 13, 19

**Status:** ready-for-agent

- [ ] The dashboard is observed updating during a live STUB run, not against a pre-generated file.
- [ ] Restarting the server rebuilds identical state from the log.
- [ ] The 127.0.0.1 binding is confirmed and the dashboard is reachable over an SSH tunnel.
