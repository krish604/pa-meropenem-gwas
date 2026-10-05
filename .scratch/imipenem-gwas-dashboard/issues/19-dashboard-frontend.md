# 19: Dashboard frontend

**What to build:** One page showing the four numbers a researcher actually checks, where the run is, which isolates failed, and what just happened.

**Blocked by:** 18

**Status:** ready-for-agent

- [ ] Metric cards show annotated, running, failed and elapsed.
- [ ] One card per stage showing waiting, running, done or failed.
- [ ] One square per isolate per stage, grey when waiting, blue when running, green when done, red when failed.
- [ ] A live event log shows what just happened rather than requiring it to be inferred from the grid.
- [ ] Colour is defined once in CSS custom properties and derived from, so dark mode works and one state is never two colours.
- [ ] Status is not conveyed by hue alone.
- [ ] Where the event stream carries no fact the UI shows an explicit 'not reported', never a plausible zero.
