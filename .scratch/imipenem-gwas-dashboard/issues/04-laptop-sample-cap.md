# 04: 20-sample laptop cap with an explanatory failure

**What to build:** The laptop configuration refuses to run more than 20 samples, and the refusal tells the user what the limit is, which config key caused it, and what to change.

**Blocked by:** 03

**Status:** ready-for-agent | verified

- [ ] A run whose manifest exceeds the cap stops before doing any work.
- [ ] The message names the limit, the offending count, the config key and the file to change.
- [ ] The cap is read from configuration, not hard-coded.
- [ ] A manifest at exactly the cap runs normally; one over it does not.

**Build state: `verified`.** Phase 2, and now **wired**: `run_pipeline` calls `enforce_sample_cap` as soon as the manifest is known and before any stage runs. Proven end to end by a 21-sample cohort that stops with the limit named, while 20 samples still runs - the boundary the 20 committed fixtures could not test, since every one of them sits exactly at the cap.
