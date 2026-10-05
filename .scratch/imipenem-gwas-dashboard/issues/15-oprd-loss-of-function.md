# 15: OprD loss-of-function features

**What to build:** Two separate features describing the OprD porin, because 'the porin is gone' and 'the porin is present and broken' are different observations with different mechanisms, and collapsing them destroys the distinction.

**Blocked by:** 14

**Status:** ready-for-agent

- [ ] oprD_absent is produced when the locus is wholly absent.
- [ ] oprD_LoF is produced when the locus is present but carries a frameshift, a premature stop or an internal deletion.
- [ ] The two features are kept distinct end to end and neither is derived from the other.
- [ ] Near-boundary indels are not silently classified as disruptive.
- [ ] Presence of the locus is never treated as evidence of susceptibility.
