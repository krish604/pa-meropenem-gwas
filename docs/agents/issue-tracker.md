# Issue tracker: Local Markdown

Issues and specs for this repo live as markdown files in `.scratch/`.

## Why local markdown

This repository has **no git remote** and neither the `gh` nor the `glab` CLI is
installed on the development machine. GitHub/GitLab-backed tracking is therefore
not available, and local markdown is the tracker of record.

## Conventions

- One feature per directory: `.scratch/<feature-slug>/`
- The spec is `.scratch/<feature-slug>/spec.md`
- Implementation issues are one file per ticket at `.scratch/<feature-slug>/issues/<NN>-<slug>.md`, numbered from `01`, never a single combined tickets file
- Triage state is recorded as a `Status:` line near the top of each issue file (see `triage-labels.md` for the role strings)
- Comments and conversation history append to the bottom of the file under a `## Comments` heading

## This repo's feature slug

The current effort is the *P. aeruginosa imipenem AMR-GWAS pipeline and live
monitoring dashboard*. Its slug is:

```
.scratch/imipenem-gwas-dashboard/
```

`TASKS.md` at the repo root is a **human-readable roll-up** of that directory, not
the tracker itself. If the two disagree, the files under `.scratch/` win.

## When a skill says "publish to the issue tracker"

Create a new file under `.scratch/<feature-slug>/` (creating the directory if needed).

## When a skill says "fetch the relevant ticket"

Read the file at the referenced path. The user will normally pass the path or the issue number directly.

## Wayfinding operations

Used by `/wayfinder`. The **map** is a file with one **child** file per ticket.

- **Map**: `.scratch/<effort>/map.md` (the Notes / Decisions-so-far / Fog body).
- **Child ticket**: `.scratch/<effort>/issues/NN-<slug>.md`, numbered from `01`, with the question in the body. A `Type:` line records the ticket type (`research`/`prototype`/`grilling`/`task`); a `Status:` line records `claimed`/`resolved`.
- **Blocking**: a `Blocked by: NN, NN` line near the top. A ticket is unblocked when every file it lists is `resolved`.
- **Frontier**: scan `.scratch/<effort>/issues/` for files that are open, unblocked, and unclaimed; first by number wins.
- **Claim**: set `Status: claimed` and save before any work.
- **Resolve**: append the answer under an `## Answer` heading, set `Status: resolved`, then append a context pointer (gist + link) to the map's Decisions-so-far in `map.md`.

## PRs as a request surface

**Off.** There is no remote, so external PRs cannot be triaged here. Flip this
flag in this file if a remote is ever added and PR triage is wanted.
