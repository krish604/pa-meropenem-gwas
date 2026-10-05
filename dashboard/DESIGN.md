# Dashboard design — Phase 0

Scope of this document: the read-only results dashboard. It fixes the data
contract between the pipeline's artefacts and the browser, the endpoint
surface, the frontend module contract, and the seven design decisions UI-D1
through UI-D7. It implements nothing. Every code reference below was read
directly; nothing here is quoted from memory.

This document uses `UI-D<n>` for the dashboard's own decisions because
`spec.md` already uses `D1`–`D13` for the pipeline's, and `D4` means something
else in each. Where a decision lines up with a spec decision it says so.

---

## 0. Verdict on the existing observatory: EXTEND, do not replace

**`papipeline/observatory/` is not reusable as the dashboard's state store, and
the dashboard does not read it.** Re-verified in this session:

- `papipeline/execution/store.py:322` — `ON CONFLICT(run_key, stage, subject)
  DO UPDATE`, with the unique key being those three columns
  (`store.py:40-44`, `store.py:173-179`).
- `papipeline/run.py:370` and `papipeline/run.py:405` — every stage is announced
  with `subject=""`.

So one row per `(run_key, stage)`. A 900-isolate run writes **16 rows**, not the
~14,400 cells the grid needs. The store cannot express the thing the dashboard
exists to show, and widening its key is a change to the pipeline's execution
layer, which is out of scope and forbidden here.

What the dashboard **does** reuse, by import:

| reused | from | verified |
|---|---|---|
| the 16-stage ordered list | `papipeline/run.py:130 STAGE_ORDER` | read |
| stage table filenames + headers | `papipeline/execution/contracts.py:40 STAGE_TABLES` | read |
| folded-step table paths | `contracts.py:134 INTERNAL_TABLES` | read |
| the stage→table key mapping | `papipeline/stages/report_tables.py:231 STAGE_TABLE_KEYS` | read |
| paths the Snakefile declares but `contracts.py` does not | `report_tables.py:219 SNAKEFILE_ONLY_TABLES` | read |
| banner-aware TSV reading | `papipeline/io/tsv.py` `read_tsv`, `iter_tsv_dicts`, `MISSING_SENTINELS` | import-tested |
| absent-table semantics with a named reason | `report_tables.TableRead` / `read_table` | import-tested |
| execution-state enum | `papipeline/execution/state.py:26 StageState` | import-tested |
| patristic-distance unit constants | `papipeline/stages/similarity.py` `DISTANCE_UNIT`, `units_sidecar_path` | import-tested |
| R16 power-flag wording | `papipeline/stages/reporting.py:123-141` | import-tested |

**Import-time side effects: none.** All sixteen candidate modules were imported
under stdout/stderr capture with the working tree checked before and after;
every one imported silently, and `git status --porcelain` was unchanged
afterwards. No reimplementation of banner handling is therefore required.

---

## 1. The canonical stage list, and the names the UI displays

`papipeline/run.py:130 STAGE_ORDER` is authoritative. It is a 16-tuple, and
`run.py:152 EXECUTION_ORDER` is `STAGE_ORDER`, so the display order and the
execution order are the same sequence and the UI needs only one list.

| # | code name (`STAGE_ORDER`) | spec.md D1 name | differs |
|---|---|---|---|
| 1 | `validation` | `validate` | yes |
| 2 | `annotation` | `annotate` | yes |
| 3 | `mlst` | `mlst` | — |
| 4 | `amr` | `amr` | — |
| 5 | `virulence` | `virulence` | — |
| 6 | `variants` | `variants` | — |
| 6a | `cohort_variants` | `cohort_variants` | — |
| 7 | `pangenome` | `pangenome` | — |
| 8 | `recombination` | `recombination` | — |
| 9 | `phylogeny` | `phylogeny` | — |
| 10 | `similarity` | `similarity` | — |
| 11 | `phenotype` | `phenotype` | — |
| 12 | `gwas` | `gwas` | — |
| 13 | `convergence` | `convergence` | — |
| 14 | `cooccurrence` | `combination` | yes |
| 15 | `reporting` | `report` | yes |

Four names differ, and spec.md's own numbering is 15 rows carrying 16 names
(`cohort_variants` is row `6a`), so "spec.md has fifteen stages" is a row count
rather than a name count.

**The UI displays the code names.** `STAGE_ORDER` is what `run_manifest.json`,
the event log, and every `TableRead` key use, so displaying anything else means
translating in both directions on every screen. `GET /api/stages` returns
`name` (code), `spec_name`, and `spec_number`, and the stage page shows the code
name with the spec name beside it in smaller type. `cooccurrence` is therefore
labelled `cooccurrence` / `combination` — a reader who knows the spec finds
their word.

### 1.1 A correction to the carry-forward notes

`workflow/Snakefile` defines **22** rules, not 20: the 16 stages plus `all`,
`provenance`, `synthetic_fixtures`, `full_run`, `clean_results` and
`list_stages`. Verified with `grep -c "^rule " workflow/Snakefile` → `22`. The
six extras are right; the total was not.

---

## 2. `ResultsSource`: the adapter, and how it is found

### 2.1 The interface

```python
class ResultsSource(Protocol):
    kind: Literal["live", "bundle"]          # which probe matched
    root: Path                                # realpath'd, once, at open time
    run_manifest: Mapping[str, Any] | None    # parsed, or None + reason
    manifest_reason: str                      # why there is no manifest

    def stage_table(self, stage: str) -> TableRead          # never raises
    def internal_table(self, name: str) -> TableRead        # never raises
    def optional(self, key: str) -> Path | None            # non-contracted artefacts
    def newick_paths(self) -> list[Path]                    # every tree it can see
    def report_files(self) -> list[Path]                    # md/html/figures only
```

Two rules, both load-bearing:

1. **Nothing raises for an absent artefact.** Every getter returns the
   `report_tables.TableRead` shape — `present`, `rows`, `reason` — because
   `read_table` (`report_tables.py:111`) already established that a reader which
   crashes on a stage that did not run is worse than one that says the stage did
   not run. The dashboard adopts that shape rather than inventing a second one.
2. **Paths come from `contracts`, never from literals.** A path the dashboard
   hard-codes is a second answer to a question `contracts.py` already answers,
   which is the drift class `contracts.py`'s own docstring says it exists to
   prevent.

### 2.2 AUTO-DETECT probe order

Tried in order; the first that satisfies its acceptance test wins; **all** are
recorded so the failure message can name every one.

| # | probe | value | accepted when |
|---|---|---|---|
| 1 | `--results-root` | CLI flag | dir exists |
| 2 | `PA_DASH_RESULTS_ROOT` | env var | dir exists |
| 3 | bundle | `$PA_DASH_BUNDLE` or the flag's sibling | `<root>/02_stage_outputs/run_manifest.json` exists **or** `<root>/02_stage_outputs` exists |
| 4 | live REAL | `PipelineConfig.results_root(RunMode.REAL)` | `run_manifest.json` or `intermediate/stages/` exists |
| 5 | live TEST | `…(RunMode.TEST)` | as above |
| 6 | live STUB | `…(RunMode.STUB)` | as above |

Probes 4–6 come from `config/machine overlay` paths — `paths.results_root:
results` in `config/machines/laptop.yaml:81` — through
`config/loader.py:917 results_root`, which also honours
`PIPELINE_RESULTS_ROOT` and appends the lowercased mode. Nothing here is a
literal path.

The `kind` is decided by **where `run_manifest.json` was found**, not by which
probe matched: `live` if it is at the root, `bundle` if it is under
`02_stage_outputs/`. This is what makes the bundle layout tolerant — the probe
order is a search, but the two adapters are a function of where the manifest
landed, so a bundle that is missing `02_stage_outputs` and has the manifest at
its root is still served correctly.

### 2.3 The failure message

When nothing matched, the API returns `503` with the whole search in it, in
order, each with the reason it failed:

```
No results source could be opened. Six probes were tried, in order:
  1. --results-root                : not given
  2. PA_DASH_RESULTS_ROOT          : not set
  3. bundle (PA_DASH_BUNDLE)       : not set
  4. live REAL  results/real       : no run_manifest.json and no intermediate/stages/
  5. live TEST  results/test       : no run_manifest.json and no intermediate/stages/
  6. live STUB  results/stub       : no run_manifest.json and no intermediate/stages/
No run has been executed in this worktree, or its outputs are elsewhere. Point
the dashboard at one with --results-root or PA_DASH_RESULTS_ROOT.
```

Probes 1–3 report *absence of configuration*, which is a different fact from
probes 4–6 reporting *the directory exists and holds nothing*. The dashboard
never collapses those two into "not found".

---

## 3. The authoritative table: stage → file → columns → rows → producer

`read_tsv` was opened against every committed fixture in this session with
`micromamba run -n pa-amr python`; row counts below are what that reader
returned, not what `wc -l` reported. The two differ on every file that carries a
`#` banner, which is most of them.

**Two different kinds of file are in play and they must not be conflated:**

- **input fixtures** under `test_data/` — what a TEST run *reads*. Counted below.
- **stage output tables** under `results/<mode>/intermediate/stages/` — what a
  run *writes*. **No results tree exists in this worktree or in any worktree
  visible here** (`ls results` → `No such file or directory`), so every output
  table below is `not produced`, with the reason given.

### 3.1 Stage output tables — the ResultsSource schema

`root` is the opened results root; `stage_dir` is `<root>/intermediate/stages`.
Producers are `module:function`.

| # | stage | path | format | columns | TEST rows | producer |
|---|---|---|---|---|---|---|
| 1 | `validation` | `intermediate/stages/01_validation.tsv` | TSV+banner | `sample_id, assembly_size, contig_count, n50, n90, largest_contig, gc_content, ambiguous_bases, basic_quality_status, species_confirmation, contamination_status, completeness_status, duplicate_status, flag_reasons` (14) | not produced | `run.py:_run_one_stage("validation")` → `stages/validation.py` |
| 2 | `annotation` | `intermediate/stages/02_annotation_summary.tsv` | TSV+banner | `sample_id, n_records, n_named_genes` (3) | not produced | `run.py:_run_one_stage("annotation")` → `stages/annotation.py` |
| 3 | `mlst` | `intermediate/stages/03_mlst.tsv` | TSV+banner | `sample_id, ST, alleles, MLST_status, mlst_scheme, allele_database` (6) | not produced | `run.py:_run_one_stage("mlst")` → `stages/mlst.py` |
| 4 | `amr` | `intermediate/stages/04_amr.tsv` | TSV+banner | `sample_id, antibiotic, determinant, gene, variant, determinant_type, mechanism, evidence_source, database, database_version, confidence, claim_status, identity_pct, coverage_pct` (14) | not produced | `run.py:_run_one_stage("amr")` → `stages/amr.py` |
| 5 | `virulence` | `intermediate/stages/08_virulence.tsv` | TSV+banner | `sample_id, virulence_factor, gene, category, database, database_version, confidence, identity_pct` (8) | not produced | `run.py:_run_one_stage("virulence")` → `stages/virulence.py` |
| 6 | `variants` | `intermediate/stages/variants.tsv` | TSV+banner | `sample_id, chrom, pos, ref, alt, qual, filter, GT, AC, AN, DP4, MQ, MQ0F` (13) | not produced | `run.py:_run_one_stage("variants")` → `stages/variants.py` (provisional; `contracts.py:83` marks these PROVISIONAL) |
| 6a | `cohort_variants` | `intermediate/stages/cohort_variants.tsv` | TSV+banner | `chrom, pos, ref, alt, ac, an, af` (7) | not produced | `run.py:_run_one_stage("cohort_variants")` → `stages/cohort_variants.py` |
| 7 | `pangenome` | `intermediate/stages/pangenome_summary.tsv` | TSV+banner | `metric, value` (2) | not produced | `run.py:1203` records `paths["pangenome_summary"]` |
| 8 | `recombination` | `intermediate/stages/recombination.tsv` | TSV+banner | `node, n_snps, mean_branch_length, recombination_detected` (4) | not produced | `run.py:derive_recombination_tables` (`run.py:1814`) |
| 9 | `phylogeny` | `intermediate/stages/10_phylogeny.tsv` | TSV+banner | `tree, n_tips, n_internal_nodes, has_branch_lengths, matches_manifest` (5) | not produced | `stages/phylogeny.py:TreeSummary.as_row` (`phylogeny.py:~95`) |
| 10 | `similarity` | `intermediate/stages/similarity.tsv` | TSV+banner | `sample_id, distances` (2) | not produced | `stages/similarity.py` |
| 11 | `phenotype` | `intermediate/stages/11_phenotype.tsv` | TSV+banner | `sample_id, antibiotic, phenotype, MIC, MIC_unit, zone_diameter, zone_unit, source` (8) | not produced | `run.py:_run_one_stage("phenotype")` → `stages/phenotype.py` |
| 12 | `gwas` | `intermediate/stages/12_gwas.tsv` | TSV+banner | `feature, feature_type, effect, p_value, adjusted_p_value, effect_size, frequency, lineage_distribution, model` (9) | not produced | `run.py:1367` → `stages/gwas.py` |
| 13 | `convergence` | `intermediate/stages/13_convergence.tsv` | TSV+banner | `determinant, independent_lineages, branch_count, distribution, convergence_category` (5) | not produced | `run.py:1396` → `stages/convergence.py` |
| 14 | `cooccurrence` | `intermediate/stages/14_cooccurrence.tsv` | TSV+banner | `feature_a, feature_b, feature_type, n_a, n_b, n_both, statistic, statistic_value, adjusted_p_value, interpretation_limit` (10) | not produced | `run.py:1452` → `stages/cooccurrence.py` |
| 15 | `reporting` | `intermediate/stages/16_figure_manifest.tsv` | TSV+banner | `figure, kind, n_items` (3) | not produced | `run.py:1499 write_tsv`, manifest from `run.py:_figure_manifest` |

### 3.2 Folded steps — written, consumed, not stages

`contracts.py:134 INTERNAL_TABLES`. The UI lists these beside their owner stage
and never as a 17th–20th stage.

| folded step | path | columns | owner |
|---|---|---|---|
| `structural_variants` | `intermediate/stages/07_structural_variants.tsv` | `sample_id, variant_id, variant_type, position, affected_gene, size, evidence, confidence, call_status, mge` (10) | `amr` |
| `regulators` | `intermediate/stages/06_regulators.tsv` | `sample_id, gene, variant, variant_type, position, reference, alternate, effect, mechanism, confidence, evidence_source, call_status` (12) | `amr` |
| `mechanisms` | `intermediate/stages/05_mechanisms.tsv` | `sample_id, antibiotic, determinant, mechanism, evidence_level, gene, gene_type, notes` (8) | `cooccurrence` |
| `master_table` | `intermediate/stages/15_master_table.tsv` | `sample_id, antibiotic, phenotype, amr_gene, amr_variant, chromosomal_mutation, regulator, mechanism, structural_variant, mlst, lineage, virulence_profile, gwas_feature, gwas_status, convergence_status, confidence, evidence_notes` (17) | `reporting` |

### 3.3 Snakefile-declared paths `contracts.py` does not carry

From `report_tables.py:219`, with the file that declares each:

| key | path | declared by |
|---|---|---|
| `gene_presence_absence` | `intermediate/stages/gene_presence_absence.tsv` | `workflow/Snakefile:265` |
| `core_genes` | `intermediate/stages/core_genes.tsv` | `workflow/Snakefile:266` |
| `accessory_genes` | `intermediate/stages/accessory_genes.tsv` | `workflow/Snakefile:267` |
| `alignment_summary` | `intermediate/stages/10_alignment_summary.tsv` | `papipeline/run.py:1250`, conditional |

### 3.4 Non-contracted artefacts the UI may show, each optional

| artefact | path | status |
|---|---|---|
| run manifest | `<root>/run_manifest.json` | `run.py:2464 _write_run_manifest` |
| pre-run manifest | `<root>/run_manifest.json` | `scripts/common/write_provenance.py:101` — **different schema**, see § 4.1 |
| figure data | `intermediate/stages/figure_data/*.json` (13 files) | `run.py:1509 _write_figure_data` |
| tree | `intermediate/phylogeny/*.nwk` | `config/loader.py:890 phylogeny_dir` |
| similarity units | `intermediate/stages/similarity.units.json` | `similarity.py:109 units_sidecar_path`, suffix `similarity.py:93` |
| variant accounting | `<variants workdir>/variants_provenance.json` | `variants.py:68 PROVENANCE_NAME`; best-effort, may be absent |
| reports | `results/reports/<name>.md`, `.html`, `figures/*` | `run.py:812 reports_root`; names from `reporting.py:2081` |
| event log | `<root>/status/events.jsonl` or repo `status/events.jsonl` | `scripts/emit.py:120` default |

### 3.5 TEST input fixtures — measured, not estimated

`read_tsv` on every committed fixture, `micromamba run -n pa-amr python`, 2026-10-05.
`wc -l` is given beside it because the difference is the `#` banner, and a
reader who trusts `wc -l` will be wrong by exactly the banner count.

| fixture | `wc -l` | banner | **`read_tsv` rows** | columns |
|---|---|---|---|---|
| `test_data/metadata/sample_metadata.tsv` | 24 | 3 | **20** | `sample_id, assembly_path, lineage, data_class` |
| `test_data/phenotype/imipenem_phenotype.tsv` | 26 | 5 | **20** | `sample_id, antibiotic, phenotype, MIC, MIC_unit, source` |
| `test_data/phylogeny/tree_metadata.tsv` | 23 | 2 | **20** | `sample_id, tree_tip_label, lineage_label, st, source` |
| `test_data/intermediate/mlst/mlst_results.tsv` | 23 | 2 | **20** | `sample_id, ST, alleles, MLST_status, mlst_scheme, allele_database` |
| `test_data/intermediate/gwas/gwas_features.tsv` | 24 | 3 | **20** | `sample_id` + 8 `type__label` feature columns |
| `test_data/intermediate/amr/amr_determinants.tsv` | 9 | 3 | **5** | 14 (matches stage 4 contract + `claim_status`) |
| `test_data/intermediate/regulators/regulator_variants.tsv` | 19 | 3 | **15** | 12 (matches `INTERNAL_TABLES["regulators"]`) |
| `test_data/intermediate/virulence/virulence_factors.tsv` | 71 | 3 | **67** | 8 (matches stage 5 contract) |
| `test_data/intermediate/structural_variants/structural_variants.tsv` | 17 | 3 | **13** | 10 (matches `INTERNAL_TABLES`) |
| `test_data/intermediate/pangenome/gene_presence_absence.tsv` | 203 | 2 | **200** | `gene, sample_id, present` |
| `test_data/intermediate/pangenome/core_genes.tsv` | 8 | 2 | **5** | `gene` |
| `test_data/intermediate/pangenome/accessory_genes.tsv` | 8 | 2 | **5** | `gene, n_samples` |
| `test_data/intermediate/pangenome/pangenome_summary.tsv` | 8 | 2 | **5** | `metric, value` |
| `test_data/intermediate/annotation/TEST_PA_001.annotation.tsv` | 10 | 2 | **7** | 10 |
| `test_data/standins/unbuilt_stages/variants.tsv` | 12 | 10 | **1** | 6 |
| `test_data/standins/unbuilt_stages/cohort_variants.tsv` | 15 | 13 | **1** | 7 |
| `test_data/standins/unbuilt_stages/recombination.tsv` | 8 | 6 | **raises** | duplicate header |
| `test_data/standins/unbuilt_stages/similarity.tsv` | 15 | 13 | **raises** | duplicate header |

Non-TSV fixtures: `test_data/phylogeny/tree.nwk` (511 B, 20 tips),
`test_data/similarity/real_tree.nwk` (435 B, 10 tips),
`.../gubbins/recombination.node_labelled.final_tree.tre` (602 B, 20 tips),
`.../gubbins/recombination.filtered_polymorphic_sites.fasta` (1,080 B, 40 lines),
`.../gubbins/recombination.per_branch_statistics.csv` (2,353 B, 40 lines, 12
tab-separated columns under a `.csv` name),
`test_data/phylogeny/core_alignment.fasta` (5,120 B),
`test_data/phylogeny/core_snp_alignment.fasta` (380 B).

### 3.6 Three findings that change the reader design

1. **`read_tsv` raises on two committed fixtures.** The `recombination` and
   `similarity` stand-ins have `standin_not_computed` repeated as three column
   names, and `read_tsv` refuses a header with duplicates
   (`io/tsv.py:92-98`). A dashboard that calls `read_tsv` and lets the exception
   reach the client would 500 on two of the sixteen stages. **The dashboard
   therefore reads through `report_tables.read_table`'s shape, not `read_tsv`
   alone**: catch `PipelineError`, and answer `present=False` with the
   exception's text as the reason. It still *imports* `read_tsv` — it does not
   reimplement banner handling.
2. **A header-only file is not a result.** `read_tsv` refuses it
   (`io/tsv.py:135`) and every stage contract requires `min_rows(table, 1)`
   (`contracts.py:247`). `report_tables.read_table` pre-checks this
   (`report_tables.py:133`). The dashboard adopts the same pre-check and the
   same sentence, so the UI says *not produced* and never renders an empty grid
   that reads as a measurement of zero.
3. **The committed stand-in headers no longer match `contracts.py`.**
   `standins/unbuilt_stages/variants.tsv` carries 6 columns; `STAGE_TABLES`
   declares 13. `UNBUILT_STAGES` is now **empty** (`run.py:188`), so these
   fixtures are stale for stages that exist. The dashboard therefore reports the
   header **the file actually carries** and cross-checks it against the contract
   as *informational*, never as a hard filter — a file whose columns have drifted
   is still the run's record, and hiding it would be worse than flagging it.

---

## 4. `run_manifest.json` — two writers, two schemas

This matters and is easy to get wrong: **two different functions write the same
filename.**

`papipeline/run.py:2464 _write_run_manifest` (after a run) writes:

```
pipeline_version, run_mode, generated_at_utc, antibiotic, python, platform,
stages {stage -> state}, stages_skipped {stage -> reason}, outputs {key -> path},
tools_detected {name -> {executable, version, available}},
references [ {reference_id, tool, tool_version, database, database_version,
              version_status} ]
```

`scripts/common/write_provenance.py:101` (before a run) writes a **subset plus
three keys the other lacks**:

```
pipeline_version, run_mode, generated_at_utc, python, platform,
references,                      # same row shape, from the same build_provenance
stage_tool_requirements {stage -> {tool -> {available, version}}},
unpinned_references [reference_id],
note                             # a prose string
```

and it has **no** `stages`, **no** `stages_skipped`, **no** `outputs`, **no**
`tools_detected`, **no** `antibiotic`.

The dashboard reads whichever is there and says which: `manifest_writer` is
`run` or `provenance`, inferred from the presence of `stages`. A manifest with no
`stages` key does **not** mean sixteen stages did not run — it means the
pre-run writer produced it, and every stage is `not_run` with the reason *"no
stage record: this manifest was written by `write_provenance.py` before the run,
which records no stage outcomes"*. That is a different sentence from *"the run
did not record this stage"*, and the UI shows the first.

`references` rows come from `run.py:1607 build_provenance` in both writers, so
the provenance page has one shape to render.

### 4.1 The manifest's own state vocabulary

`run.py:886 mark(name, state="completed")` is the only writer of `stages`. Values
that appear in the code: `completed`, `skipped_validated`
(`run.py:1549`), `skipped_no_phenotype` (`run.py:1369`). `stages_skipped`
(`run.py:736`) carries free-text reasons such as
`"disabled in configuration (analysis.gwas is false)"` and
`"unbuilt (ticket 14: ...) and not requested"`.

### 4.2 The stage-state taxonomy — UI badge vocabulary

Six badges, and **the badge is the single source of truth for the wording**
(UI-D2). No page writes a status sentence from its own logic; it asks
`state_badges` for the label, the reason and the tone.

| badge | meaning | emitted when |
|---|---|---|
| `completed` | ran and its declared output is present and readable | `stages[s] == "completed"` **and** the table is `present` |
| `running` | announced started, no terminal event | a `start` event with no later `done`/`fail` |
| `failed` | the process did not succeed | an event with `event == "fail"` for the stage, or `stages[s] == "failed"` |
| `refused` | it declined to run, by name | `s` in `REAL_REFUSING_STAGES` and `run_mode == "REAL"`; or `s` in `UNBUILT_STAGES`; or `stages[s]` starts with `skipped_` and the run recorded a cause |
| `not_assessed` | ran, and its output cannot answer the question | `stages[s] == "completed"` but the table is absent or header-only |
| `not_run` | never ran | absent from `stages`, or in `stages_skipped` |

Every non-`completed` badge carries a **non-empty `reason`**, sourced in this
order:

1. `stages_skipped[s]`, verbatim.
2. `REAL_REFUSING_STAGES[s]`, verbatim — e.g. `"no REAL engine: ReferenceEngine
   cannot handle lineage confounding"`.
3. `TableRead.reason` from `report_tables.read_table` — e.g. *"no file at
   `<path>`"*, or the header-only sentence.
4. The event that explains it: the `fail` event's timestamp and subject.
5. `"absent from run_manifest.stages and no event explains it"` — the last
   resort, never the first.

`not_assessed` and `not_run` are kept apart on purpose. A stage that ran and
produced nothing answerable is not a stage that never ran, and the pipeline's
own reporting makes the same distinction at `reporting.py:465`.

`papipeline/execution/state.py:26 StageState` (`PENDING, RUNNING, SUCCEEDED,
FAILED, INCOMPLETE, INVALID, RETRYING`) is **not** the UI vocabulary: it is the
execution layer's vocabulary and its rows come from a store the dashboard does not
read (§ 0). The UI's six badges map onto it in the provenance page's "execution
layer" column and nowhere else.

### 4.3 The oprD verdict — a hole, and the UI must show it as one

`papipeline/stages/reporting.py:810 OPRD_NOT_ON_DISK` states it in the pipeline's
own words: `contracts.py` declares no oprD table in `STAGE_TABLES` or
`INTERNAL_TABLES`, `workflow/Snakefile` declares no oprD output, and
`run.py`'s in-memory `oprd_locus_resolution` is the only copy of those verdicts.

So there is **no file** holding per-isolate oprD verdicts and their evidence. The
UI's oprD page therefore has exactly two honest states, and both are in the
`Verdict` vocabulary from `adapters/oprd_locus.py:93`:

- *verdicts available* — only if a future run writes them; the endpoint probes
  for a contracted path and a bundle-supplied table, in that order;
- *not produced* — with the reason above, quoted, plus the list of probe paths.

It **must not** fall back to the master table's `chromosomal_mutation` column or
to `viz.oprd_status_per_sample`, because that function is the one the
`oprd_locus` docstring was written to correct: on the ten-isolate cohort it
reported `absent` for 8 isolates that all carry the gene, 34 of the 36
`gene=oprD` CDS features being OprD/OprP/OprQ paralogs. Surfacing that as an
oprD verdict would manufacture the study's central negative.

The verdict vocabulary the UI must render verbatim (`oprd_locus.py:93`):

```
resolved                     <- the only one that asserts the locus is present
refused:no_orthologous_hit
refused:insufficient_coverage
refused:ambiguous_locus
refused:ambiguous_second_locus
```

and the display states (`viz.py:58`): `intact, disrupted, absent, not_assessed`.
The structural verdicts `absent` and `disrupted` come only from a confirmed
lesion; `not_assessed` is what every refusal maps to, never `absent`.

---

## 5. REST endpoints

Every list endpoint pages, sorts and filters **server-side** and returns
`{items, total, offset, limit, sort, filters}`. `total` is the count over the
whole declared artefact, computed before paging, and it is the number the UI
shows next to a page — never `items.length`. Every list endpoint also returns
`basis`: `{n, artefact, path, rows_total}` so a statistic is never displayed
without the thing it was computed on (UI-D7).

| method | path | returns |
|---|---|---|
| GET | `/api/health` | `{ok, source_kind, root, manifest_writer, mode}` |
| GET | `/api/run` | run summary: mode, antibiotic, n_samples, versions, platform, counts |
| GET | `/api/run/counts` | aggregate counts for the metric cards |
| GET | `/api/stages` | all 16, with `name`, `spec_name`, `spec_number`, `state`, `reason`, `n_rows`, `path`, `present` |
| GET | `/api/stages/{name}` | one stage: contract header, actual header, drift, row count, events |
| GET | `/api/stages/{name}/rows` | paged rows of that stage's principal table |
| GET | `/api/tables/{key}` | paged rows of a folded-step or Snakefile-declared table |
| GET | `/api/isolates` | the joined per-isolate table (§ 5.1) |
| GET | `/api/isolates/{sample_id}` | per-isolate detail: every contributing row |
| GET | `/api/oprd` | oprD verdicts, or `not produced` with probes (§ 4.3) |
| GET | `/api/pangenome` | `pangenome_summary.tsv` as metric/value + the three gene tables' counts |
| GET | `/api/tree` | structural JSON, server-assigned node ids (UI-D6) |
| GET | `/api/tree/tips` | paged tip metadata |
| GET | `/api/similarity` | downsampled matrix + units sidecar (§ 5.2) |
| GET | `/api/gwas` | `12_gwas.tsv` paged, plus `n_tests` and the QQ-ready p-value array |
| GET | `/api/gwas/top` | top hits by `adjusted_p_value` |
| GET | `/api/convergence` | `13_convergence.tsv` paged |
| GET | `/api/cooccurrence` | `14_cooccurrence.tsv` paged |
| GET | `/api/provenance` | tool versions, references, unpinned list, manifest shape |
| GET | `/api/provenance/bakta` | annotation-reuse provenance + `"Bakta executions: N"` |
| GET | `/api/provenance/logs` | tripwire / watcher log tails, paged |
| GET | `/api/provenance/digests` | sha256 of each stage artefact |
| GET | `/api/reports` | report files: path, size, mtime, kind |
| GET | `/api/reports/content` | rendered text of one report (UI-D4 containment) |
| GET | `/api/events` | SSE stream (§ 6) |
| GET | `/api/events/state` | polling fallback: full state snapshot + `cursor` |
| GET | `/api/launcher/capabilities` | what the launcher may legally do here |
| GET | `/api/launcher/status` | is a run in progress; last event; observed timings |
| POST | `/api/launcher/preflight` | dry-run checks; **starts nothing** |

### 5.1 Query-parameter schema, stated once

Every list endpoint accepts and every one is honoured server-side:

| param | type | default | notes |
|---|---|---|---|
| `offset` | int ≥ 0 | 0 | |
| `limit` | int 1–1000 | 200 | clamped, and the clamp is reported in `limit_applied` |
| `sort` | string | the table's declared column order | `<column>` or `-<column>`; unknown column → `400` naming the column, never a silent fallback |
| `q` | string | — | substring over the endpoint's declared searchable columns |
| `filter.<column>` | repeated | — | exact match; `filter.<column>__ne`, `__in` (comma-separated), `__null` (`1`/`0`) also accepted |
| `columns` | csv | all | projection; unknown column → `400` |
| `include_total` | bool | `1` | `0` lets a grid skip the count at the cost of not showing it |

`total` is always the unfiltered count unless a filter is given, in which case it
is the filtered count and `total_unfiltered` carries the other number. A page
that shows "12 of 900" after filtering to 12 must still be able to say it.

### 5.2 The isolate table

One row per `sample_id` in the manifest — or, where no manifest is readable, in
the union of the per-sample tables, with `membership_source` naming which. A
sample present in one table and absent from another **still gets a row**: the
cells that have no source read `not assessed`, not `0` and not blank. This is
rule 9 of `docs/scientific_rules.md` on the client.

| column | source | absent-value behaviour |
|---|---|---|
| `sample_id` | manifest / union | never absent |
| `st` | `03_mlst.tsv` | `not_assessed` if `MLST_status` is `no_call` or the row is missing |
| `amr_genes` | `04_amr.tsv`, comma-joined | empty **list** (`[]`) is a finding; a missing row is `not_assessed` — the two differ |
| `amr_point_mutations` | `04_amr.tsv` rows with non-empty `variant` | `reporting.py:383` — `variant` non-empty means "the tool says this element is mutated", **not** "disruptive". The column is labelled that way in the UI. |
| `n_virulence` | `08_virulence.tsv` count | `0` when the table is present and holds no row for the sample; `not_assessed` when the table is absent |
| `oprd_state` | § 4.3 | `not_assessed` unless a verdict source exists |
| `oprd_verdict` | § 4.3 | as above |
| `phenotype_sir` | `11_phenotype.tsv` for the run's `antibiotic` | `ND` is a recorded value, distinct from absent |
| `lineage` | `15_master_table.tsv` `lineage`, else `tree_metadata.tsv` `lineage_label` | the master table **refuses** a sentinel lineage (`data_contract.md:259`), so a sentinel reaching the UI is a bug worth surfacing, not rendering |
| `n_snv` | `variants_provenance.json` `per_isolate[].n_snv` | every key is present at zero when the file exists, so `0` is a real measurement there; absent file → `not_assessed` |
| `n_indel` (non-SNV) | same, `n_indel` | as above |

Filters: `phenotype_sir`, `lineage`, `oprd_state`, `has_amr_gene`,
`has_point_mutation`, `virulence_min`. Sort: any of the above plus `sample_id`.

### 5.3 The similarity matrix

`similarity.tsv` is one row per sample holding a whole `distances` cell
(`contracts.py:118`). 900 rows × 900 cells. Above **300 tips** the server
subsamples by single-linkage on the stage-9 tree so the picture keeps its shape,
and returns:

```
{ n_total: 900, n_shown: 300, downsampled: true,
  method: "single-linkage on the stage-9 tree, k=300",
  units: <similarity.units.json verbatim>,
  order: [sample_id...], matrix: [[...]] }
```

`units` is `similarity.units.json` passed through **unmodified** —
`quantity`, `definition`, `is_a_snp_count: false`
(`similarity.py:96 DISTANCE_UNIT`, `:131`). The UI must render that the matrix
is substitutions per site and never SNPs, which rule 11 of
`scientific_rules.md` makes a claim about what may be written down.

### 5.4 Provenance contents

Tool versions and availability (`tools_detected` or `stage_tool_requirements`,
whichever the manifest writer produced) · the `references` rows with
`version_status` and the `unpinned_references` list · the annotation-reuse
record and the literal string `"Bakta executions: N"` · tripwire and watcher log
tails with a size cap · SHA-256 per stage artefact · git HEAD hashes and
`git status --porcelain` of the opened root · the config snapshot
(`science.yaml` plus the machine overlay actually loaded) · free disk on the
results volume and system load average.

**"Bakta executions: N" is reported as a count or not at all.** The annotation
stage may reuse a previous run's Bakta output rather than re-invoking the tool
(`HEAD` commit 19e82d1). The UI's claim is exactly the number, and the reuse
provenance is what makes the number checkable. If neither the run manifest nor
the reuse record carries it, the field reads `not reported` — never `0`.

---

## 6. Events: the source, and how the client knows which it got

### 6.1 The choice

**The dashboard reads `status/events.jsonl`, and it is the only mechanism it
uses.** The observatory's in-process `EventBus` is not read.

Justification, in the order the arguments actually hold:

1. **The file is the only mechanism with no consumer.** `scripts/emit.py:4-16`
   says in its own docstring that nothing reads it. `AGENTS.md` records the same
   two-mechanism problem. Adding a third consumer of the `EventBus` means
   importing a process that is off by default
   (`runtime.observatory.enabled: false`) and whose store cannot represent the
   grid (§ 0).
2. **The file survives the process.** The `EventBus` lives in the runner's
   process; a run that dies takes its buffer with it. The file is appended with
   `flush()` + `fsync()` on every record (`emit.py:110-112`), so a run killed
   mid-write leaves a replayable log — which is what `emit.py:85-93` was written
   to guarantee.
3. **It is the same source regardless of entry point.** `spec.md` D9 says the
   native orchestrator emits the same events as the Snakemake rules. A dashboard
   reading one file works for both; a dashboard reading the in-process bus only
   works for one.
4. **The cost is a reader, and the reader is small.** Replay-and-tail is one
   offset-tracked file reader. The observatory would have cost a store migration.

This does not retire the other mechanism; it records which one the dashboard
depends on, so the bridging ticket can retire either without ambiguity. The
observatory is untouched.

### 6.2 The event schema, verbatim

```json
{"t": "<iso8601 UTC, 'Z'>", "event": "start|done|fail", "stage": "<name>", "sample": "<id|'all'>"}
```

Four fields. Three values for `event` (`emit.py:38`); `"all"` for an aggregate
rule (`emit.py:41`). A fifth field or a fourth value is a contract change, not a
tolerance — the reader **skips such a line with a warning naming its line
number** rather than coercing it.

`stage` in the log is a **Snakemake rule name**, which is not always the code
stage name (`rule validate` vs `validation`). The reader maps rule→stage through
an explicit table and reports unmapped names in the stream's own diagnostics
rather than dropping them.

### 6.3 SSE and the polling fallback

```
GET /api/events                     -> text/event-stream
GET /api/events/state?cursor=<int>  -> application/json, 200
```

The SSE response sets `Content-Type: text/event-stream`, `Cache-Control:
no-cache`, `X-Accel-Buffering: no`, and emits:

- `retry: 2000` once, at connect;
- one `event: snapshot` carrying the full replayed state **before** any
  `event: delta`, so a client that connects mid-run is never blank;
- `event: delta` per new line, batched on a 100 ms timer rather than per line;
- a `data:` comment heartbeat every 15 s, so an idle run does not look dead;
- `id: <byte offset>` on every event, so `Last-Event-ID` resumes exactly.

**How the client tells which it got.** The frontend opens the stream with
`EventSource` and arms a 3 s timer simultaneously. Three outcomes, all
explicit:

1. `open` + any frame received → SSE. The client announces
   `transport: "sse"` in the header and stops the timer.
2. timeout, then `GET /api/events/state` succeeds → **polling**. The client
   announces `transport: "poll"`, polls every 2 s carrying `cursor`, and shows
   the same banner either way.
3. timeout *and* `/api/events/state` fails → the run summary is still served and
   the event panel shows *"event stream unavailable: `<reason>`"*. It does not
   show an empty log, because an empty log reads as "nothing has happened".

The header carries the transport in both cases, so a screenshot of the UI says
which mechanism it is showing.

---

## 7. Frontend module contract

Every page module is an ES module:

```js
export function mount(container, ctx) { /* ... */ return teardown }
```

`teardown` is `() => void`, must remove every listener, abort every in-flight
`AbortController`, and close any `EventSource` the module opened. A page that
returns nothing is a leak, and the router refuses it in development.

```js
ctx = {
  api,      // { get(path, params, signal), post(path, body, signal) } — every call abortable
  router,   // { go(path), current(), on(event, fn), off(event, fn) }
  theme,    // { tokens } — CSS custom properties, read not written
  banner,   // { mode, power(n), underpowered } — see below
  badges,   // { state(name) -> {key,label,reason,tone}, tone(key) }
  on,       // { runState(fn), delta(fn) } — the live transport, one subscription
}
```

`ctx.badges` is the **only** place a stage or isolate state is turned into a
label, a colour token and a reason string. No page may assemble a status
sentence itself (UI-D2). This is why the taxonomy lives in one backend object:
two pages cannot disagree about what `not_assessed` means, and neither can
disagree with the report.

`ctx.banner` carries the mode banner and the D3 component's parameters.
`ctx.theme.tokens` is read-only; a component that hard-codes a colour is a
component that breaks in dark mode.

### 7.1 The D3 banner is one component

`PowerBanner` is a single component, mounted wherever a statistic is shown, and
it takes its text from the backend rather than composing it. Its wording is
`reporting.py:123-131`, verbatim:

```
n=<N>, underpowered, not a finding
Ruling R16. This is a description of the isolates in this run and nothing more.
It is not evidence about Pseudomonas aeruginosa imipenem susceptibility in
general, and no association, rate or difference below may be cited as a finding.
```

`N` is the cohort size the statistic was **actually computed on**, taken from
`basis.n` on the response — not the cohort size the page was told about. A page
that shows a filtered count with the full-cohort flag is exactly the error R16
exists to prevent. Status is never carried by hue alone: each badge pairs its
colour with a glyph and a text label.

---

## 8. File ownership map

Nobody edits another owner's files. A change one owner needs is written as a
proposal in that owner's report, not applied.

| owner | paths |
|---|---|
| **BACKEND** | `dashboard/server/**` · `dashboard/__main__.py` · `scripts/run_dashboard.sh` |
| **SHELL** | `dashboard/static/index.html` · `dashboard/static/app.js` · `dashboard/static/router.js` · `dashboard/static/css/**` · `dashboard/static/pages/pipeline.js` · `dashboard/static/pages/provenance.js` |
| **DATA** | `dashboard/static/pages/isolates.js` · `.../isolate.js` · `.../tables.js` · `.../summaries.js` · `.../oprd.js` |
| **VIZ** | `dashboard/static/pages/tree.js` · `.../matrix.js` · `.../gwas.js` · `.../evolution.js` |
| **QA** | `dashboard/tests/**` · `dashboard/fixtures/**` |
| **DESIGN** (this doc) | `dashboard/DESIGN.md` · `dashboard/openapi.yaml` · `dashboard/static/vendor/**` · `pa-artifacts/ui/**` |

`dashboard/openapi.yaml` is DESIGN-owned so that the contract and its
description cannot drift; BACKEND changes it in the same commit as the route.

---

## 9. UI-D1 … UI-D7

### UI-D1 — the results root is read-only, and the state directory is the only write

The dashboard **opens** the results root once, `realpath()`s it, and thereafter
opens every file read-only. It never creates, truncates, renames or deletes
anything under it, and it never writes a cache there — a dashboard that
rewrites the artefact it is reporting on has destroyed the evidence.

Every write goes to the state directory (§ 10), which is outside the results
root by construction. The enforcement is not a convention: `ResultsSource`
exposes no write method, and the containment check of UI-D4 would reject a write
path anyway.

### UI-D2 — absence is stated, never rendered as a value

Every absent fact is one of `not produced`, `not_run`, `not_assessed`,
`not reported`, or a named refusal string — never `0`, `""`, `[]` or a blank
cell. The wording is owned by one place (the six badges, § 4.2) and one rule
decides which: **a table that exists and holds no row for a sample is a
measurement of zero; a table that does not exist is not a measurement.**

This is not a UI preference; it is `docs/scientific_rules.md` rule 9 and
`report_tables.py:30-35` restated at the presentation layer. The 2026-09-30 note
in `data_contract.md:104-120` is the cautionary case: the `ast_*` columns are
dropped at the write, so a reader of `11_phenotype.tsv` cannot tell a provenance
gap from a source that recorded nothing. The UI must not create a second such
gap by rendering "no value" as "absent".

Degradation is explicit everywhere. The launcher, for one: the per-isolate
timings in `RUNBOOK_900.md` are not available to this session, so the launcher
field accepts them when a runbook supplies them and reads **"timings not
recorded"** when it does not. It never invents a per-isolate estimate.

### UI-D3 — the underpowered flag is one component, injected everywhere

§ 7.1. `PowerBanner` is mounted by the shell wherever a statistic renders; a
page cannot opt out, because the injection point is the statistic-rendering
helper and not the page. `N` comes from `basis.n`.

### UI-D4 — file serving is contained under the opened root, and a non-local bind needs a token

Three independent checks, all required:

1. **Resolve, then contain.** `candidate = (root / requested).resolve()` (which
   follows every symlink), then require
   `candidate == root or root in candidate.parents`. `os.path.commonpath` on the
   *resolved* pair, and a check that `candidate` is not `root` itself unless the
   endpoint serves the root. A traversal that escapes fails here even if it
   never touched `..` — e.g. a symlink inside the root pointing at `/etc`.
2. **Allowlist, then suffix.** Only `md`, `html`, `json`, `tsv`, `txt`, `csv`,
   `nwk`, `newick`, `tre`, `svg`, `png`, `pdf`, `log`, `fa`, `fasta`, `faa`,
   `gff`, `gff3`, `vcf`, `faa` are served, and `text/html` is served with
   `Content-Disposition: attachment` plus `X-Content-Type-Options: nosniff` so
   the report's own markup cannot execute in the dashboard's origin.
3. **Bind check.** The server refuses to start on a non-loopback host unless
   `PA_DASH_TOKEN` is set, and with `0.0.0.0` always. The token is compared with
   `secrets.compare_digest`, is required as `X-PA-Dash-Token` on every `/api/**`
   request, and is **never** accepted as a query parameter — a query string
   lands in shell history, in the browser's address bar and in every proxy log.

`127.0.0.1` and `::1` are loopback and need no token, which is the normal case:
the dashboard is reached through an SSH tunnel (`spec.md` D10).

### UI-D5 — a byte-offset index, cached in the state directory

Paging a 200,000-row variant table by re-reading it per request is quadratic
over a session. The index maps each key to the byte range of its row, so a page
is `seek` + `readline`.

File format, one file per table, `<state>/index/<stage>.<key16>.idx.json`:

```json
{
  "version": 1,
  "key": "<sha256 hex, see below>",
  "source": { "path": "<realpath>", "size": 12345, "mtime_ns": 1234567890 },
  "header": ["sample_id", "chrom", "..."],
  "key_column": "sample_id",
  "offsets": { "TEST_PA_001": [4096, 4241], "TEST_PA_002": [4241, 4389] },
  "row_count": 2
}
```

**Key** — `sha256(realpath + "\0" + str(size) + "\0" + str(mtime_ns))`. Size and
mtime together, so a truncated file and a rewritten file both miss.

**Invalidation** — the entry is discarded and rebuilt when the stored `key`
differs from the key computed now, when `version` differs, or when
`source.path` differs. There is no time-based expiry: mtime and size are a
statement about the file, and a poll of the file's stat is cheaper and more
honest than a guess about how long a cache lasts. A run still in progress
rewrites its tables, so its mtime moves and the index is rebuilt; that is the
intended behaviour, and it is why the index is never treated as a result.

Writing is atomic — temp file in the same directory, then `os.replace` — so a
kill mid-write leaves the previous index rather than a truncated one. The index
is **disposable by construction**: deleting the state directory costs a rebuild
and nothing else, which is the property that lets it be written at all under
UI-D1.

### UI-D6 — node identity is the node's position in the tree, never its label

`node_id` is assigned by the server as the node's **postorder index** from the
root: one depth-first pass, children left to right, counter incremented after
the subtree. Verified on `test_data/phylogeny/tree.nwk`: 38 nodes, 20 tips, ids
unique, and the mapping reproducible across two parses of the same bytes.

The internal Newick label is **never** an identity, and the fixtures show why:

- `test_data/similarity/real_tree.nwk` has eight internal nodes carrying
  `100/100` **three times**, `99.9/92`, `100/92`, `83.4/88`, `99/100` — and one
  with **no label at all** (`None`). Support values are neither unique nor
  present.
- `.../gubbins/recombination.node_labelled.final_tree.tre` has `Node_1` …
  `Node_19`, which are gubbins' own bookkeeping and mean nothing to a browser.

So `/api/tree` returns `{node_id, parent_id, children, length, is_tip, label,
support, tip_metadata}` where `label` is present **for display only** and every
client-side reference is by `node_id`. A postorder index is a position in *this*
file, not a biological identity, and the endpoint says so in its `identity`
field — so a client cannot mistake "the node the server called 17" for "clade 17".

One trap, and it is a data-corruption trap rather than a display one:
**dendropy's `preserve_underscores` defaults to `False`**, and with the default
every `TEST_PA_001` comes back as `TEST PA 001`. That silently breaks the
`sample_id` join on every page. The parse therefore passes
`preserve_underscores=True`, and the endpoint **asserts** that the tip set
equals the manifest set the way `phylogeny.validate_tree_samples`
(`phylogeny.py:266`) does, reporting the differences rather than drawing a tree
whose tips cannot be joined.

### UI-D7 — every number is computed server-side over the whole artefact, and arrives with its basis

The 900-isolate grid cannot be aggregated in the browser: a client that counts
what is on the page reports the page, not the run. So **no page computes a
count, a rate or a total.** Each response carries:

```
basis: { n: <isolates the statistic used>, artefact: "<stage>", path: "<realpath>",
         rows_total: <int>, filtered: <bool> }
```

and the UI renders `n` next to the statistic and passes it to the D3 banner. A
page that needs a count asks the server for it with the same filters and gets a
number with its own `basis`. This is what makes the D3 flag correct: `n` is the
denominator the statistic actually used.

The grid itself is a **custom virtualised canvas** component, not a table
library — 14,400 cells is beyond a DOM table's budget, and no general grid
handles the six-badge vocabulary, which is data the library would have to be
told about anyway.

---

## 10. State directory

Default `~/.pa_dashboard/`, overridable with `PA_DASH_STATE_DIR`.

```
~/.pa_dashboard/
  index/          <stage>.<key16>.idx.json    UI-D5, disposable
  logs/           dashboard.log                the dashboard's own log
  session.json    opened root, opened-at, transport   UI-D1: outside the root
```

**Contents are exactly these three kinds, and all three are derived.** Nothing in
the state directory is a result, a setting that changes what is reported, or
anything a run needs. Deleting the whole directory costs one rebuild and no
information — which is what makes writing here consistent with UI-D1.

**Writes.** Only `index/` and `logs/`. Never inside the opened results root, and
the containment check would reject it if a caller tried. `session.json` is
written once at open and rewritten only when the root changes.

**Invalidation.** `index/` entries invalidate on key mismatch (UI-D5), and the
whole directory may be deleted at any time. A stale index is a performance
regression, never a wrong answer, because the index stores offsets and the
`key` check is what makes a stale one unreachable.

---

## 11. Threat model (UI-D4)

| threat | defence | residual |
|---|---|---|
| `GET /api/reports/content?path=../../../../etc/passwd` | resolve-then-contain; `etc/passwd` is not under the root | none |
| symlink inside the root pointing at `/etc/passwd` | `resolve()` follows it, so the containment check sees the real target and fails | a symlink to another file *inside* the root is served; that is intended |
| encoded traversal (`%2e%2e%2f`, double-encoding) | the framework decodes before the handler sees it; the handler decodes nothing itself and re-checks after `resolve()` | none |
| a report containing `<script>` | `Content-Disposition: attachment` + `nosniff`; the UI renders report text with `textContent`, never `innerHTML` | none |
| a huge file exhausting memory | responses are streamed; every table read is capped at `limit ≤ 1000` rows | a 10 GB `.log` is streamed, not buffered |
| bind on `0.0.0.0` without a token | startup refuses; `PA_DASH_TOKEN` required for any non-loopback host | a token on a shared host is still a shared secret — SSH tunnel is the answer, and the token is the fallback |
| token in a URL | rejected: header only, `secrets.compare_digest` | none |
| CSRF from a page on another origin | the API is `GET`-only except the two launcher endpoints, which take no state-changing action in Phase 0; `POST /api/launcher/preflight` starts nothing | a later launcher POST must add an origin check — noted, not built |
| the dashboard writing to the results | UI-D1: no write method exists on the adapter | none |
| a symlinked state directory | the state dir is `resolve()`d at startup and logged | none |

**The read-only guarantee** is structural: the results root is opened once, every
file handle under it is opened `'r'`, and the adapter's interface has no write
method. Nothing the dashboard does can modify what it reports.

---

## 12. ASSUMPTIONS

Every unverified assumption, with what breaks if it is wrong.

| # | assumption | why | blast radius if wrong |
|---|---|---|---|
| A1 | **The delivery bundle layout** is `01_bakta_input/`, `02_stage_outputs/`, `03_report/`, `04_run_info/`, `05_validation/`, `06_for_900_isolates/RUNBOOK_900.md` | given in the brief; **not described in any doc or branch visible here** — grepped `docs/`, `.scratch/`, `TASKS.md` | probe 3 fails to find a bundle that exists. Mitigated: the probe order is tolerant and probes 4–6 do not depend on it; the `kind` is decided by where `run_manifest.json` lands, not by which probe matched. Fix: edit one probe, nothing else. |
| A2 | `02_stage_outputs/` is the bundle's copy of `intermediate/stages/` | inferred from the name | stage tables are not found in a bundle. The live adapter is unaffected. |
| A3 | `03_report/` holds the `.md`/`.html` | inferred from the name | report listing is empty; the endpoint returns `not produced` with probes. |
| A4 | `05_validation/` holds the validation-stage artefacts | inferred from the name | same as A3. |
| A5 | **`RUNBOOK_900.md` per-isolate timings are unavailable** | no delivery bundle exists (`PA_AMR_pipeline_real_run_10_isolates` is absent) | the launcher shows "timings not recorded" (UI-D2) — which is the designed behaviour, not a failure |
| A6 | **`run_manifest.json` carries no `bakta_executions` field** | not in either writer's payload; the annotation-reuse commit changed the stage, not the manifest | the provenance page shows `not reported` for "Bakta executions: N". The count must then be derived from the reuse provenance, or the field is added — a BACKEND proposal, not a silent `0`. |
| A7 | No tripwire/watcher log path is contracted | none found in `contracts.py` or the Snakefile | the log panel is empty with `not produced` and a probe list. |
| A8 | `variants_provenance.json` lands in stage 6's workdir, whose location is not contracted | `variants.py:562` passes `workdir / PROVENANCE_NAME`; the workdir path is not in `contracts.py` | per-isolate `n_snv`/`n_indel` read `not_assessed`. Probed at several candidate paths, all listed in the error. |
| A9 | Report filenames are `test_pipeline_report.{md,html}` (TEST/STUB) and `Pseudomonas_aeruginosa_Imipenem_AMR_Report.{md,html}` (REAL) | `reporting.py:2081` and `:2087`, confirmed in source | the report listing is discovered by extension, not by name, so it degrades to an empty list rather than a wrong link |
| A10 | The stage→stage-table mapping is stable at the commit this dashboard was designed against | read at 19e82d1 | a new stage appears in `STAGE_ORDER` and `STAGE_TABLES`; the UI renders it `not_run` with a reason. It cannot render a stage the pipeline does not have. |
| A11 | `dendropy` 5.1.0 parses all three committed Newick files | **verified** — all three parse with `string=`/`path=` and `force-unrooted` | none at the current commit; a future Newick dialect would need the fallback noted below |
| A12 | The observatory stays off by default and untouched | `runtime.observatory.enabled: false`; the dashboard does not import it | none. If it were enabled, the dashboard is still unaffected — that is the point of § 6.1 |
| A13 | `similarity.units.json` is written beside `similarity.tsv` | `similarity.py:109` | the matrix renders with `units: null` and the UI says the units are not recorded — **not** "substitutions per site" as an assumption |
| A14 | `stage_tool_requirements` and `tools_detected` are both rendered from whichever the manifest writer produced | both writers verified | the tool matrix shows `not reported` rather than an empty matrix |

**Fallback if A11 ever fails:** `papipeline/stages/phylogeny.py` already contains
a hand-rolled recursive-descent Newick parser (`_NewickParser`) that extracts tip
labels correctly and handles quoted labels. It collects tips only — no internal
structure — so a fallback would have to be extended to build nodes, in the
dashboard's own code, not in `papipeline/`. It is named here so that the
fallback is a decision rather than an emergency.

---

## 13. What this design does not decide

- **The launcher.** Phase 0 designs `POST /api/launcher/preflight`, which
  **starts nothing**. Which command a launch would run, with what arguments, is
  Phase 1's, and it must not be inferred from the pipeline's CLI by a dashboard.
- **Which mechanism gets retired.** § 6.1 picks one for the dashboard to read.
  Retiring the other is a pipeline change and belongs to the bridging ticket.
- **Colour values.** The tokens are read from CSS custom properties; the palette
  itself is SHELL's, and status must be distinguishable without hue alone.
- **The `analyses/` ordering of the grid.** STAGE_ORDER fixes the stage axis.
  Whether isolates are ordered by lineage, by phenotype or by `sample_id` is a
  DATA/VIZ decision that this document leaves open, and it must not be decided by
  insertion order, which is what makes a grid unreadable at 900.