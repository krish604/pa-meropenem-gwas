# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root, or
- **`CONTEXT-MAP.md`** at the repo root if it exists: it points at one `CONTEXT.md` per context. Read each one relevant to the topic.
- **`docs/adr/`**: read ADRs that touch the area you're about to work in. In multi-context repos, also check `src/<context>/docs/adr/` for context-scoped decisions.

If any of these files don't exist, **proceed silently**. Don't flag their absence; don't suggest creating them upfront. The `/domain-modeling` skill (reached via `/grill-with-docs` and `/improve-codebase-architecture`) creates them lazily when terms or decisions actually get resolved.

## File structure

Single-context repo (most repos):

```
/
├── CONTEXT.md
├── docs/adr/
│   ├── 0001-event-sourced-orders.md
│   └── 0002-postgres-for-write-model.md
└── src/
```

Multi-context repo (presence of `CONTEXT-MAP.md` at the root):

```
/
├── CONTEXT-MAP.md
├── docs/adr/                          ← system-wide decisions
└── src/
    ├── ordering/
    │   ├── CONTEXT.md
    │   └── docs/adr/                  ← context-specific decisions
    └── billing/
        ├── CONTEXT.md
        └── docs/adr/
```

## This repo

**Single-context.** There is no `CONTEXT-MAP.md` and no monorepo layout. `CONTEXT.md`
and `docs/adr/` do not exist yet; the `/domain-modeling` skill will create them when
the imipenem GWAS terms and decisions are actually resolved.

## Pre-existing design documentation

This repo predates the skills setup and already carries a substantial, authoritative
design record that agents must read before changing behaviour:

- `docs/architecture.md` — how the pipeline is put together
- `docs/data_contract.md` — the file formats exchanged between stages
- `docs/scientific_rules.md` — the scientific invariants (e.g. R is never inferred from I)
- `docs/reproducibility.md` — determinism guarantees for the synthetic fixtures
- `docs/design/` — numbered design documents (`01`–`10`), including `08-migration-plan.md`
  and `09-test-strategy.md`

**Treat these as binding.** If a change contradicts one, surface the contradiction
explicitly rather than quietly overriding it.

## Use the glossary's vocabulary

When your output names a domain concept (in an issue title, a refactor proposal, a hypothesis, a test name), use the term as defined in `CONTEXT.md`. Don't drift to synonyms the glossary explicitly avoids.

If the concept you need isn't in the glossary yet, that's a signal: either you're inventing language the project doesn't use (reconsider) or there's a real gap (note it for `/domain-modeling`).

## Flag ADR conflicts

If your output contradicts an existing ADR, surface it explicitly rather than silently overriding:

> _Contradicts ADR-0007 (event-sourced orders), but worth reopening because…_
