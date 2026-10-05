# Handoff: folding the stage taxonomy to spec.md D1's fifteen

Written because a session ran out of context mid-change, twice. The base branch
is clean and green; the work is not lost. Read this before touching
`papipeline/run.py` or `workflow/Snakefile`.

## Where everything is

| What | Where |
|---|---|
| Working branch, clean and green | `feature/imipenem-pa-pipeline` at `7ad7fdc` |
| First attempt, red, committed as-is | `feature/taxonomy-wip` at `e65d17d` |
| Second attempt, red, committed as-is | `wip/contracts-additive` |
| Plain-text backup of the dirty files | `/tmp/taxonomy_wip_backup/` |
| Backup of the additive contracts slice | `/tmp/contracts_additive.py` |

Verified at `7ad7fdc`: **917 unit, 202 integration, 0 failures, 0 skips.**

The two WIP branches are recoverable, not authoritative. Both commit messages
record what they contain and why they are red. Cherry-pick from them for
reference; re-verify everything.

## What the change is

`spec.md:349-353` already decided this, so it is migration work and not a design
question:

> The existing sixteen-stage code is **mapped onto** these, not discarded and not
> allowed to define its own taxonomy. The existing `mechanisms`, `regulators`,
> `structural_variants` and `integration` modules become internal steps within
> the rule that owns them; their tests and scientific content are preserved.
> Existing `STAGE_TABLES` and stage ordering constants are rewritten to the
> fifteen.

Agreed fold: `integration` → `report` (15); `mechanisms` → `combination` (14),
with the report only rendering it; `regulators` → `variants` (6), content and
tests preserved, oprD reaching the GWAS feature family in D6;
`structural_variants` → `amr` (4).

`structural_variants` **is** foldable. `sv.py` normalises a pre-computed SV table
from the tool-output root and contains no reference to `minimap2`, `samtools`,
`PAO1`, alignment or reference. It needs no ticket-14 input.

Three of the fifteen have no implementation: `variants`, `recombination`,
`similarity`. They are declared in the DAG, fabricate in STUB, and raise
`NotImplementedError` naming ticket 14 outside STUB.

## The finding that shapes the work: it is one atomic change

`STAGE_TABLES`, `STAGE_ORDER`, `PREREQUISITES`, the four folded bodies, the
`internal_table_path` call sites, STUB fabrication and the Snakefile must land
together. Demonstrated twice, independently:

1. Switching `STAGE_ORDER` to fifteen while the bodies are still stages stops
   five Snakefile rules getting `05_mechanisms.tsv`, `06_regulators.tsv` and
   `15_master_table.tsv`.
2. Adding the three contracts while the stages are not yet in `STAGE_ORDER`
   breaks every test that iterates `STAGE_TABLES`, because a STUB run writes no
   table for a stage it does not run.

Every boundary between the pieces is a red commit. Land it as one commit.

### Gates before committing

- Full unit **and** integration suite green.
- The reachability mutation check fails when a spec edge is deleted.
- TEST dry run reaches all fifteen stages.
- TEST cannot run `variants`, `recombination` or `similarity`.
- STUB dry run **and** a real STUB run both pass.

Verify each sub-step (full suite after each, re-read each edited region) but do
not commit them separately. The commit message must say the pieces landed
atomically because every boundary between them is red, list what it contains,
and state the suite was green at that commit.

## The cascade, and the answer to the 0-row question

Spec edges are declared literally in every mode, so `gwas` consumes `variants`
and `phylogeny` consumes `recombination`. Those stages refuse in TEST/REAL, so in
TEST their tables come from committed, contract-shaped stand-in fixtures supplied
by `stage_inputs()` as external inputs.

**A 0-row table is not enough, and 1 row is the minimum.** Two independent
reasons, and only one of them is obvious:

- *Nothing reads them.* `gwas.run()` (`stages/gwas.py:769`) reads
  `intermediate_root/gwas/gwas_features.tsv`, not `variants.tsv`.
  `phylogeny` reads `data_root/phylogeny/`. No code opens
  `variants.tsv`, `recombination.tsv` or `similarity.tsv`. So the spec edges are
  aspirational today, and 0 rows breaks no parser.
- *The contract rejects them.* `stage_spec` applies `non_empty(table)` **and**
  `min_rows(table, 1)` to every stage's principal table
  (`execution/contracts.py:155,163`). A 0-row stand-in fails validation whenever
  the observatory is on. One row is the minimum the contract demands.

So: stand-ins with **one** row each, labelled as stand-ins in the file header and
in the report, with a TASKS.md note to delete them when ticket 14 lands. A test
must fail if a stand-in is ever treated as real output.

The three contract column sets are **provisional** placeholders. They are marked
as such in `contracts.py`; ticket 14 sets them against the real output.

## Outstanding, in order

1. The atomic commit above.
2. Stand-in fixtures for TEST (1 row each), labelled, with the TASKS.md note.
3. Test: the TEST dry run reaches all fifteen; TEST cannot run the three
   unbuilt stages; the reachability mutation check still bites.
4. `stage_inputs()` unit test for mode=REAL: the full external set is returned
   and a missing file is named. Do **not** open `allow_real_mode`.
5. `spec.md:62` and `:562`: dated status note beneath each, "Intent; not
   implemented; see TASKS.md". Do not change the claim itself.
6. `.opencode/agents/dashboard-ui.md`: correct to describe the in-process
   `EventBus` and the absence of an `events.jsonl` consumer.
7. New ticket: bridge `events.jsonl` into the observatory `EventBus`
   (tail/replay). Do not start it.
8. TASKS.md: built/wired/verified, ticket 14's three absent stages as not built,
   the stand-in fixtures, and the two observability mechanisms that do not meet.

## Working method — the rules that were violated

- **Never `git checkout` a modified file.** It silently reverted an hour of
  uncommitted work. Back up to `/tmp` first, or commit to a WIP branch first.
  Branch checkouts are fine.
- **No `sed` or anchor-string edits on `run.py` or the Snakefile.** `sed` prepends
  indentation, so anchors copied from its output do not match the file, and two
  such edits silently did nothing. Read a region, edit with a precise tool, then
  re-read the region.
- **Full suite after each sub-step, and do not start the next while red.**
- Commit a WIP branch before a risky rewrite.

## Also worth knowing

- `git checkout -b` then commit-then-switch-back is the pattern that preserved
  this work. Use it before any large refactor.
- `stage_inputs()` lives in the Snakefile and switches on `RESOLVED_MODE`;
  `results_root` is redirected by `PIPELINE_RESULTS_ROOT` in the config loader
  (not `--config`, because the rules shell out to a process that recomputes its
  own paths from config).
- Nothing in `papipeline/` reads `status/events.jsonl`. The observatory serves
  `/api/events` from an in-process `EventBus`. Two observability mechanisms that
  do not meet; item 7 above is the bridge.
- `data/` and `test_data/` must stay untouched. Redirect with
  `PIPELINE_RESULTS_ROOT`, and copy fixtures rather than editing them.
