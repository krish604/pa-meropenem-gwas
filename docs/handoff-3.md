# Handoff 3 — round 11 closeout on `integration-r3`

Supersedes nothing: `docs/handoff-taxonomy.md` and `docs/handoff-ticket14.md` describe
work that is now **merged**. Read those for the reasoning behind ticket 14 and the stage
taxonomy; read this one for the state of the tree today.

---

## 1. State

| | |
| --- | --- |
| Integration branch | `integration-r3` |
| Head at time of writing | see §6 |
| Full suite | **2832 passed / 45 skipped / 0 failed** (464 s) |
| Round-11 baseline | 2448 passed / 33 skipped / 0 failed, at `cf8476f` |
| Merges | nine, all `--no-ff`, all two-parented, zero unresolved conflicts |
| Repairs | one authorised test-fixture fix (`9d7ad59`) |

`git status --porcelain` is empty. `feature/imipenem-pa-pipeline` is still `7cc9ece` and
must stay there. `papipeline/testing/synthetic.py` has **zero** changes across the whole
round: `git diff --stat a9c9548..HEAD -- papipeline/testing/synthetic.py` is empty.

### The +12 skips over baseline are intentional, and here is why

All twelve are `got empty parameter set`:

* `tests/unit/test_standins.py:51,58,67,73,89,107,111` — 7
* `tests/integration/test_stage_taxonomy.py:113,119,139,157,285` — 5

Both files parametrise over `run.UNBUILT_STAGES` / `UNBUILT_WITH_STANDIN`, which
`dag-resolve` emptied when recombination gained a dispatch branch. Nothing was silenced to
make a test pass; the same fact is asserted positively by
`tests/integration/test_stage_taxonomy_is_self_verifying.py`. **These twelve nodes can never
execute again and should be deleted rather than left skipping.**

---

## 2. Branch list — what each contributed

| branch | merge | what it brought |
| --- | --- | --- |
| `dag-resolve` | `cd0ef8e` | Recombination dispatched with its config key in one commit; REAL DAG resolves; `stages_skipped` with reasons at WARNING; `cohort.subset_file` applied at manifest load; actionable panaroo refusals. **Sole owner of `run.py` and `workflow/Snakefile`.** |
| `variants-indels` | `ae74755` | `min_ireads` proven to reach the mpileup argv; per-isolate SNV/indel provenance sidecar; zero-indel warning in REAL. No consumer code changed. |
| `amr-organism` | `da3da3d` | `--organism` passed and validated against the tool's own organism list; refusal on an unlisted one; provenance on the stage-4 table header. |
| `tree-lineage-asc` | `59c6544` | `lineage_label` produced from stage 3; integration refuses an incomplete column by name; `+ASC` live **behind** the column filter; similarity matrix carries its units. |
| `pheno-provenance` | `e81c12c` | R7 category accounting (exclude **and** count); `provenance_report` given a production caller. |
| `ticket-16-gwas` | `e34e761` | 14 files of GWAS input work: unique-pattern reduction, `count_patterns.py`, and two now-**required** config keys. |
| `oprd-production` | `507cc45` | The oprD production structural resolver, decided from the alignment rather than a codon index. **Not yet on a pipeline code path.** |
| `seam-tests` | `8a0612c` | Call-site arity over every `stage_X.run(...)`; producer/consumer closure over all 23 REAL inputs. Test-only. |
| `docs-apply` | `f655045` | 8 docs files brought into line with merged code; machine overlays corrected to **comment-only**. |

---

## 3. Rulings R1–R10 and where each is honoured

| ruling | honoured where |
| --- | --- |
| **R1** ticket-16 off-limits list | **Lifted** for this round. Exactly one edit under it by me: `tests/unit/test_round11_config_keys.py:210`, a construction line, under the authorised 2c repair. No non-test production file was edited by me. |
| **R2**–**R3** (see agent reports) | `DAG.md`, `TREE.md`, `VARIANTS.md`, `AMR.md`, `PHENO.md` each carry their own section. R3's `--organism` narrowing is documented as *unscreened for point mutations*, never *screened and negative*. |
| **R4** lineage | `lineage.method` is read and recorded; `tree_cut` is **not** implemented and a second method raises rather than labelling anything. |
| **R5**–**R6** | R6's `min_ireads: 1` is written as **necessary and not sufficient**, with the 853-of-2447 over-10-nt caveat. |
| **R7** phenotype | `I`/`SDD`/`ND` excluded **and counted**, never merged; absent ≠ zero; unknown codes refused differently from excluded ones. |
| **R8** `+ASC` | Live behind `drop_partially_constant_columns`; the counterfactual (no filter → binary refuses) is asserted, not described. |
| **R9** recombination | Dispatched, key and `UNBUILT_STAGES` removal in the **same commit**. |
| **R10** assertions | **Zero** assertion changes by me. Across the nine merges the edited pre-existing test files and their rulings are listed in `DAG.md:177-204` and `TREE.md:225-255`. |
| **NEW RULE** config key ships with its consumer | The Phase-0 failure mode, twice: `analysis.recombination` alone broke 76 nodes; `models` alone broke 12. All four round-11 keys now land with their consumers in one commit. Verified live on all three overlays. |

---

## 4. BACKLOG

Full text, one line per finding with `file:line`:
**`/Users/raghavkrishnankv/Desktop/pa-artifacts/round11/BACKLOG.md`**

The entries that most affect code written next, in the order I would take them:

1. **`run.py` wiring, twice, one line each.** `_build_report_context` does not pass
   `phenotype_calls`, so a report `run.py` writes shows absence wording instead of the
   provenance section. And nothing in a run calls `structural_call*`, so the oprD production
   resolver is proven but unreachable. Both fields exist, both are tested; each needs one
   argument at its call site. These are the whole of CC-REPORT items 7 and 4.
2. **`docs/environment-arm64.md:137`** still asserts the withdrawn zlib diagnosis, and now
   contradicts the corrected docstring in `test_gubbins_source_build.py` that the same merge
   landed. A reader comparing the two will find them inconsistent.
3. **No `data_root` env override** (`config/loader.py:41` offers only
   `PIPELINE_RESULTS_ROOT` and `PIPELINE_ALLOW_REAL_MODE`), so a REAL dry run's
   `sample_metadata.tsv` resolves through the `data/` symlink into the linked checkout.
4. **No indel plausibility bound.** 853 of 2447 `-m 1` indels exceed 10 nt on one-read
   support, and 1 of 10 isolates was OOM-killed. Choosing the bound is a scientific decision.
5. **D1–D11 from SEAM are NOT PROVEN.** The design document is not on disk; the two seam
   tests are a reconstruction from the brief's prose.
6. **No REAL cohort run has ever executed end-to-end.** Everything above rests on
   authorised per-stage reruns against the 10 smoke isolates.

---

## 5. How to resume

### Verify where you are

```
cd /Users/raghavkrishnankv/Desktop/pa-integration
git rev-parse --abbrev-ref HEAD && git rev-parse HEAD && git status --porcelain
eval "$(micromamba shell hook -s bash)" && micromamba activate pa-amr
uptime && python -m pytest -q          # expect 2832 passed / 45 skipped / 0 failed
```

**`micromamba run -n pa-amr pytest` is wrong here** — it prepends the env's `bin` to
`PATH`, which defeats any test that shadows a tool by prepending a stub directory. Activate
the environment instead. Check `uptime` before any full run; wait 5 min if 1-min load > 10.

### If you are picking up a PARTIAL item

`CC-REPORT.md` names exactly what is missing for items 4, 7 and 8. Items 4 and 7 are each one
line in `run.py`; item 8 is a false doc row plus one deferred scientific ruling.

### If you are running the pipeline for real

```
export PIPELINE_RESULTS_ROOT=/tmp/somewhere-real
snakemake --snakefile workflow/Snakefile --config mode=REAL machine=bigmachine \
    --cores 4 --dry-run
```

* REAL is refused at DAG-build time unless `PIPELINE_ALLOW_REAL_MODE=1` is also set; the
  refusal names the cap and points at `bigmachine`. Both facts verified.
* Eight externals must be operator-provisioned before the REAL DAG builds — the manifest,
  the phenotype table, the MLST/SV/AMR/virulence caller tables, the GWAS feature matrix, and
  the panaroo core-gene alignment. `CC-REPORT.md` §6 gives the table and why each is the
  operator's.
* `results_root` is redirected by the **environment variable**, not `--config`, because every
  rule shells out to a process that recomputes its own paths from config.
* `data/` and `test_data/genomes` are **symlinks into the canonical read-only worktree**.
  Writing a placeholder under `data/` writes into that checkout. Copy fixtures; do not edit
  them in place. `DAG.md` records an agent having to undo exactly this.
* `panaroo` and `gubbins` cannot be installed on `osx-arm64`. `gubbins` has a source build
  (`scripts/build_gubbins_from_source.sh`); the conda build segfaults with exit 139 and its
  refusal names the remedy.

### Two traps that nearly produced a wrong call this round

1. **A diff against the wrong base reads as catastrophic loss.** Every round-11 branch forks
   from `7cc9ece`, so `git diff cf8476f..<branch>` reports all of round 11 as deletions
   (ticket-16 looked like −20,897 lines across 83 files; the truth is 14 files). **Always
   diff against the merge base**, and remember that a branch which lacks a file entirely
   appears to delete it.
2. **An empty env var is not a redirect.** `PIPELINE_RESULTS_ROOT=""` silently falls back to
   the repo's own results tree; a `snakemake -n` batch ran clean and green against the wrong
   paths because a `mktemp -d` had failed first. Check the resolved paths in the output.

### Standing rules that are not negotiable here

- `AGENTS.md` rules 1–7. Rule 1 is the one agents break: **never invent a package version, a
  tool flag, or a database name.** Verify with `conda search` and `<tool> --help`; if unsure,
  say so. A wrong flag is worse than no code. `DAG.md` caught itself writing
  `panaroo build --remove-invalid-files` from memory and had to correct it.
- `papipeline/testing/synthetic.py` is FROZEN. `PDC_essential.tsv` is real clinical data:
  never commit, move, rename or delete it.
- Stage by explicit path. Never `git add -A` or `.`.
- `data/`, `db/`, `test_data/genomes` and `local/` are untracked and must never be staged.
- The canonical worktree `ncbi_genomes_pdc_only_global` is READ-ONLY.

---

## 6. Provenance

* `CC-REPORT.md` — the nine FINISH LINE items with evidence commands: **six PASS, three
  PARTIAL, none FAIL**, each PARTIAL naming precisely what is absent.
* `round11/{DAG,VARIANTS,AMR,TREE,PHENO,OPRD,SEAMT,DOCS}.md` — the per-agent reports, each
  with its own R1/R10 declarations.
* `round11/I8.md` — the Phase-2 integration log: nine merges, the union verifications, and
  the 2c gate failure and its repair.
* `round11/I8-final.md` — the Phase-3 closeout log.
* `round11/BACKLOG.md` — every finding deferred, one line each with `file:line`.
* `round11/STATUS.md` — the orchestrator's gate history. Note its Gate table is stale: it
  still reads "FAIL — halted here" for Gate A, which the resume decision overturned.