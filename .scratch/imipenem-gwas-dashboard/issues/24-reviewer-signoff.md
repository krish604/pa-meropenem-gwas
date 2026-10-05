# 24: Reviewer sign-off on the definition of done

**What to build:** An independent, read-only review that verifies rather than trusts, covering the seven hard rules, the seven items of the definition of done, and scientific soundness, with UNVERIFIED reported honestly where something could not be checked.

**Blocked by:** All other tickets

**Status:** ready-for-agent

- [ ] The reviewer reports PASS, FAIL or UNVERIFIED with evidence for each item of the definition of done.
- [ ] The seven hard rules are each checked individually.
- [ ] It is confirmed that the GWAS input carries all three feature families and that the threshold comes from unique variant patterns.
- [ ] It is confirmed that no genome, database or secret is tracked by git.
- [ ] It is confirmed that no REAL-mode run was started without the user's explicit instruction.
- [ ] Anything unverifiable is reported as unverified with the reason, not quietly passed.
