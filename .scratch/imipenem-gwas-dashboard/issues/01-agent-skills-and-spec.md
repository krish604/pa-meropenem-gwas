# 01: Agent skills, issue tracker and spec

**What to build:** The repository is set up so that agents can work on it without inventing conventions: agent guidance with the project's hard rules, a declared issue tracker and triage vocabulary, four specialised subagents, and a written specification of the imipenem GWAS pipeline and dashboard.

**Blocked by:** None (can start immediately)

**Status:** done (Phase 0, committed)

- [ ] AGENTS.md records the seven hard rules and the three execution modes.
- [ ] docs/agents/ records the local-markdown tracker, the triage label vocabulary and the single-context domain-doc layout.
- [ ] Four subagents exist: pipeline, tdd-writer, dashboard-ui, reviewer. The reviewer has no edit or write permission.
- [ ] spec.md exists and is marked ready-for-agent.
- [ ] .gitignore covers status/, db/, databases/, .env and .env.* while keeping .env.example committable.
