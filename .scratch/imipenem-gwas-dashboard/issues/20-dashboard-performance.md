# 20: Dashboard performance at 900 isolates

**What to build:** Proof that the dashboard is usable at the scale the final run will actually reach: about 900 isolates across 15 stages, or roughly 13,500 cells. Measured, not asserted.

**Blocked by:** 19

**Status:** ready-for-agent

- [ ] Replay time against a realistically-sized synthetic log is measured and reported, not claimed.
- [ ] A single event updates a single cell; the grid is not rebuilt per event.
- [ ] DOM writes are batched rather than issued once per event.
- [ ] The measured numbers are reported plainly, including if they are poor.
