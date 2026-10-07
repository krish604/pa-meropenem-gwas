# Data contract

Every file the pipeline reads or writes, its schema, and the rules that make it
valid. The rules are enforced in code by `papipeline/io/tsv.py` and
`papipeline/manifest.py`, not only described here.

## Directory layout

```
config/           configuration and knowledge tables (editable inputs)
data/             REAL mode inputs        (not read in TEST mode)
test_data/        TEST mode inputs        (synthetic fixtures)
results/test/     TEST mode outputs
results/real/     REAL mode outputs
reports/          reports
```

`data/` and `test_data/` are never mixed. `PipelineConfig.data_root(mode)`
resolves them separately, and the synthetic generator refuses to write into a
directory named `data`.

## Generic TSV rules

Applied by `papipeline.io.tsv.read_tsv`:

| Rule | Behaviour on violation |
|---|---|
| File must exist | `DataContractError: TSV file not found` |
| File must be non-empty | `DataContractError: TSV file is empty` |
| Header must have unique names | `DataContractError: duplicate column names` |
| Required columns must be present | `DataContractError: missing required columns` |
| Row must not exceed the header width | `DataContractError: more fields than the header` |
| Row must meet `min_columns` | `DataContractError: too few fields` |
| At least one data row | `DataContractError: no data rows` |
| Unique columns must be unique | `DuplicateSampleError` (reports both line numbers) |
| Composite keys must be unique and fully populated | `DuplicateSampleError` / `DataContractError` |

**Missing values.** These strings become `None`, never `""` and never `0`:
`.`, `NA`, `na`, `nan`, `None`, `none`, `null`, `-`, and the empty string. On
write, `None` is emitted as `.`.

**Comments.** Lines beginning with `#` are skipped anywhere in the file, so
outputs written with a provenance header round-trip through the reader.

**Short rows** are padded with `None` rather than rejected, because a trailing
empty cell is ambiguous rather than wrong. A partially populated *composite
key* is a hard error, because it cannot be deduplicated safely.

## Sample identifiers

`sample_id` is the join key for every stage. Validated by
`manifest.validate_sample_id`:

- 2–128 characters;
- alphanumeric plus `-`, `_`, `.` only;
- no whitespace, and no `/\'":;|,()[]{}#`;
- leading/trailing whitespace is an error, not silently stripped.

A sample ID must be unique within any single-source table. A sample may appear
once per antibiotic in the phenotype table, which is why that table uses a
composite key `(sample_id, antibiotic)`.

## Input files

### `metadata/sample_metadata.tsv`

| Column | Required | Notes |
|---|---|---|
| `sample_id` | yes | unique; see above |
| `assembly_path` | no | resolved relative to the metadata dir |
| `lineage` | no | free text; the pipeline does not trust it — lineages come from the tree |
| `data_class` | no | provenance label |

A sample with no `assembly_path` yields a QC record with status `missing_input`.
It is not dropped, and it is not given a fabricated value.

### `phenotype/<antibiotic>_phenotype.tsv`

| Column | Required | Notes |
|---|---|---|
| `sample_id` | yes | |
| `antibiotic` | yes | must be configured in `config.yaml` **and** `antibiotics.tsv` |
| `phenotype` | yes | one of `R`, `I`, `S`, `SDD`, `ND` |
| `MIC` | no | numeric, **only from a measurement** |
| `MIC_unit` | no | |
| `zone_diameter` | no | numeric, **only from a measurement** |
| `zone_unit` | no | |
| `source` | no | provenance label |

Composite key `(sample_id, antibiotic)`.

**Hard rules.** `R`/`I`/`S` are never converted into an MIC or a zone diameter.
An MIC reported against `ND` is an error, because `ND` means not determined.
An absent phenotype leaves the sample missing; it is never defaulted to `S`.

The three `ast_*` columns are deliberately **not** in the table above, and that
is a decision rather than an omission — see
`tests/unit/test_phenotype_real.py::TestPerIsolateMetadataIsNotAPhenotypeColumn`,
which asserts their absence so it cannot be reversed by accident. They are also
**not** in the written `11_phenotype.tsv`, which is a separate fact from the table
above and is easy to confuse with it.

They are nevertheless **read, held and counted**, at three different layers, and
only two of the three survive:

| layer | holds `ast_*`? |
|---|---|
| the phenotype table on disk | only if the source laboratory wrote them |
| the in-memory `PhenotypeCall` | **yes** — stored and exported (`models.py:473-475`, `:493-495`) |
| **`11_phenotype.tsv` on disk** | **no** — dropped at the write |

The drop point is this contract, not the writer: the write is handed
`STAGE_TABLES["phenotype"]`, which names no `ast_*` column, and `write_tsv`
passes `extrasaction="ignore"` to the underlying writer, so an unlisted key is
dropped without error. Had it raised, the stage would have failed loudly instead
of writing a narrower table — **silent is the hazard**. The consequence worth
knowing: `provenance_report` can audit the provenance gap only because stage 11
holds the values in memory; a consumer reading `11_phenotype.tsv` back cannot
distinguish a gap from a source that recorded nothing, because on disk they are
the same thing.

`phenotype.py:233-235` populates them and `phenotype.py:742-744` reports how many
records lack each (`lacking_method` / `lacking_standard` / `lacking_edition`), so
the gap is counted rather than invisible. Only `sample_id`, `antibiotic` and
`phenotype` are required (`phenotype.py:50`), which is why an absent AST method
is a provenance gap to be reported and not a parse failure.

It matters because an absent `ast_standard`/`ast_edition` means the `R`/`I`/`S`
call cannot be traced to a breakpoint standard, which is a weaker claim than one
that can — so those counts exist to surface exactly that, and must not be dropped
from any run summary. This is the same `docs/reproducibility.md` obligation as §8
(database version) applied to the phenotype side. Note also that
`config/science.yaml` naming a standard (CLSI M100, edition `UNSUPPLIED`) is an
**intention held by this repository**, not a provenance fact about any isolate,
and a report must not present it as the latter.

**`I` and `SDD` are excluded from the primary contrast and counted, never
merged.** The primary contrast is `R` vs `S`
(`gwas.outcome.positive: [R]`, `negative: [S]`). `I`, `SDD` and `ND` are real,
measured and outside that contrast, so they are **excluded and counted** rather
than folded into either arm: `EXCLUDED_CATEGORIES` is `("I", "SDD", "ND")`, it
restates `gwas.outcome.excluded`, and a test asserts the two agree against the
loaded config. Merging an intermediate into `R` would report a resistance call
the laboratory never made; merging an `SDD` into `S` would report a plain
susceptible. The exclusion line is logged **only when there is one**, so absent
is distinguishable from zero, and a report states the same counts in a table and
in prose — including saying explicitly that nothing was excluded when that is the
case.

An **unknown** code is refused differently from a known-but-excluded one, and
both refusals name the code: `RESISTANT` is not a spelling of `R`, and mapping
it would be a guess about a clinical record, whereas `I` is a valid AST code this
model declines to compare.

### Stage input tables (under `intermediate/`)

In TEST mode these are checked-in synthetic fixtures. In REAL mode the external
tools write them into the run's intermediate directory.

**`annotation/<sample_id>.annotation.tsv`**
`sample_id`, `contig_id`, `gene_id`, `gene_name`, `product`, `gene_type`,
`start`, `end`, `strand`, `annotation_source`

**`mlst/mlst_results.tsv`**
`sample_id` (unique), `ST`, `alleles` (`locus:allele;...`), `MLST_status`
(`typed`/`partial`/`allele_incomplete`/`no_call`), `mlst_scheme`,
`allele_database`

**`amr/amr_determinants.tsv`**
`sample_id`, `antibiotic`, `determinant`, `gene`, `variant`, `determinant_type`,
`mechanism`, `evidence_source`, `database`, `database_version`, `confidence`

- `variant` is populated for point-mutated elements and `.` otherwise. AMRFinderPlus
  spells a mutated element `<gene>_<variant>` in one column (`oprD_V359L`), so stage 4
  splits it: `gene` is the gene and `variant` is the substitution. `determinant` keeps
  the tool's own spelling unchanged.
- A REAL write carries a `#` provenance header naming `organism`, `tool`, `tool_version`,
  `database` and `database_version`. The organism is the input that decides whether the
  screen covers point mutations at all: **without a valid `--organism` AMRFinderPlus
  reports no POINT rows**, so a table with no variants and a table from a run that never
  screened for variants are otherwise indistinguishable. `read_tsv` skips `#` lines.

**Stage 4 passes `--organism`, and point-mutation screening was silently OFF
before it did.** The flag was not merely unset: the adapter's own docstring
claimed that omitting it "searches the whole database, which is the correct
behaviour for a cohort whose species the config does not assert" — and that claim
was false, because `config/science.yaml` has always recorded `organism.name`
(taxid 287) and `spec.md:343` names the invocation. The config carried the value,
the spec required the flag, and the adapter read neither.

Measured on the 10-isolate smoke cohort (AMRFinderPlus 4.2.7, database
`2026-08-07.1`, `--plus` both sides, only `--organism` differing): the run
**without** the flag reproduces the previous table exactly at 144 rows with
`variant` populated in **0 of 144**, and the run **with** it gives 199 rows —
**144 identical plus 55 POINT rows**, a strict superset with nothing lost. The
`variant` column of the real cohort went from entirely empty to populated.

The flag therefore **widens and narrows**. On this cohort it was purely additive,
but a cohort screened without it must be described as *unscreened for point
mutations*, never as *screened and negative*. An `oprD_V359L` call is the
practical difference.

**The organism value is read from the tool, not typed.** `derive_organism_flag`
converts the configured space-joined name to AMRFinderPlus' underscore-joined
spelling, and refuses an empty or absent `organism.name` by name. `list_organisms`
then runs `amrfinder --database <db> --list_organisms` and parses the
`Available --organism options:` line, because a hard-coded list would be a second
statement of what the database contains, free to drift from it. `resolve_organism`
refuses a value the database does not screen for. (`--list_organisms` **requires**
`--database`; without it the call exits non-zero.)

**`regulators/regulator_variants.tsv`**
`sample_id`, `gene`, `variant`, `variant_type`, `position`, `reference`,
`alternate`, `effect`, `mechanism`, `confidence`, `evidence_source`,
`call_status`

**`structural_variants/structural_variants.tsv`**
`sample_id`, `variant_id`, `variant_type`, `position`, `affected_gene`, `size`,
`evidence`, `confidence`, `call_status` (`confirmed`/`candidate`/
`not_assessable`), `mge`

**`virulence/virulence_factors.tsv`**
`sample_id`, `virulence_factor`, `gene`, `category`, `database`,
`database_version`, `confidence`, `identity_pct`

**`variants/variants.tsv`** (stage 6)
`sample_id`, `chrom`, `pos`, `ref`, `alt`, `qual`, `filter`, `GT`, `AC`, `AN`,
`DP4`, `MQ`, `MQ0F` — one row per variant call, joined across all isolates.

**`variants/variants_provenance.json`** (stage 6, alongside the table)
`per_isolate[]` carries `sample_id`, `call_status`, `n_snv`, `n_indel`,
`n_alleles`; `totals` carries `isolates`, `n_calls`, `n_snv`, `n_indel`,
`n_alleles`, and — only when status was tracked — `n_called`, `n_not_called`.

- `call_status` is one of `called`, `no_assembly`, `call_failed`,
  `call_failed_signal`. The last is distinct from `call_failed` because a
  negative return code is a signal, not an exit status: a tool that was killed
  (the round-12 `mpileup` SIGKILL) is a different failure from one that exited
  non-zero, and it points at a different response.
- **The key is omitted entirely when status was not tracked**, rather than
  emitted as `"unknown"`. This is `docs/scientific_rules.md` rule 9: missing
  data must remain missing, and a placeholder that reads like a verdict is the
  one thing a reader cannot distinguish from the truth.
- It is a *per-isolate* verdict and is deliberately **not** a column of
  `variants.tsv`. `tests/unit/test_variants_contract.py` rejects it there, and
  correctly: the table is one row per variant, and an isolate's failure to be
  read at all is not a property of any site.

**`cohort_variants/cohort_variants.tsv`** (stage 6a)
`chrom`, `pos`, `ref`, `alt`, `ac`, `an`, `an_calls`, `af` — one row per
cohort-polymorphic site.

- `an` is the cohort size: every isolate the study contains, read or not. It is
  unchanged from before `an_calls` existed and **stays the denominator of `af`**.
- `an_calls` is how many of those isolates produced a call set at all. It exists
  because `an` alone reported round 12's 10-isolate cohort as complete while one
  genome was never read, with nothing on the face of the table to say so.
- `af` remains exactly ``ac / an``. An unread isolate counts in the denominator
  as a non-carrier: it was never observed to carry the site, and excluding it
  would renumber every frequency already published. Reporting the shortfall and
  preserving the frequency are independent, and this table does both.

**`gwas/gwas_features.tsv`**
`sample_id` (unique) plus one column per feature, named `<type>__<label>`.
Values must be binary. The `<type>` prefix must appear in
`config.gwas.feature_types`.

**`phylogeny/tree.nwk`** — Newick. Tip labels must match the manifest sample
IDs **exactly**: no extras, none missing, no duplicates. Enforced by
`phylogeny.validate_tree_samples`.

**`phylogeny/tree_metadata.tsv`**
`sample_id` (unique), `tree_tip_label`, `lineage_label`, `st`, `source`

**All five columns are now written, in TEST and in REAL alike.** Stage 9 passes
`TREE_METADATA_COLUMNS` (all five) to `write_tsv`, not the two-element
`REQUIRED_TREE_METADATA`. That distinction was a real defect: the writer built a
three-key row including `source: "iqtree"` but handed `write_tsv` the required
two-column tuple, and `write_tsv` iterates only the columns it is given — so
`source`, `lineage_label` and `st` were all dropped, and every sample then read
back as the literal string `"unknown"`. `REQUIRED_TREE_METADATA` survives as the
**read-time** subset only: those two columns are what `load_tree_metadata` keys
on and the only ones it requires.

**`lineage_label` is the MLST sequence type from stage 3**, read from that
stage's own on-disk table at its contracted path
(`mlst/mlst_results.tsv` under the run's tool-output root) — not re-derived, and
not cut from the tree. The **method is recorded**, in the file's own `#`
provenance header alongside the stage-3 path, so a result always says how its
lineages were defined. Cutting the tree into clades is **deferred**: it needs a
support threshold, a minimum clade size and a rule for isolates the tree places
ambiguously, none of which has been chosen, so `lineage.method` accepts one value
and the loader refuses anything else at load time. Naming `tree_cut` as an
accepted value would be claiming a decision this project has not made.

**`"unknown"` is a sentinel, not a lineage.** A sample with no usable sequence
type carries the literal `UNKNOWN_LINEAGE`, and it does not reach an output: the
master-table builder **refuses** on an absent, whitespace-only or sentinel
label, naming every offending sample rather than a count, because convergence
counts independent lineages from this column, co-occurrence stratifies by it and
the GWAS lineage-confound check needs more than one distinct value. A table with
a sentinel entry would report those findings against a grouping that was never
established. This is rule 9 applied to lineage: "not determined" is not
"determined to be none".

`st` and `source` are advisory. A consumer must treat their absence as normal
rather than as a malformed file.

## Step 12a artefacts: unique variant patterns

These three are **not** stage input tables and do not live under
`intermediate/`. They are written into stage 12's *work directory* - the one
`PyseerEngine` is constructed with, alongside `phenotype.tsv`, `features.rtab`,
`kinship.tsv` and pyseer's own logs - so they are deliberately **absent from
`papipeline/execution/contracts.py`**, which declares the per-*stage* tables
under `intermediate/stages/`. Naming them there would imply a consumer that
does not exist.

They exist so the multiple-testing threshold can be read rather than trusted: a
number nobody can reconstruct is a number nobody should trust.

**`patterns.txt`** — pyseer's own `--output-patterns` output. Read, and never
written by this repository. One line per variant pyseer tested: the base64 of
the MD5 of that variant's binary presence/absence vector
(`pyseer/input.py::hash_pattern`).

Two details of that encoding are easy to get wrong, and a fixture built the
obvious way is a pattern pyseer could never emit. First, `hash_pattern` hashes
`k.view(np.uint8)`, and `k` is an array of Python ints, so numpy gives it the
platform's default integer width — **8 bytes per sample** on a 64-bit machine,
not one byte each. Second, under Python 3 `binascii.b2a_base64` appends the
newline, which is what makes the file one pattern per line — and therefore what
makes `sort -u | wc -l` a count of patterns at all.

`tests/unit/test_gwas_unique_patterns.py` produces its fixtures by **calling**
`hash_pattern` and `pyseer.utils.format_output`, rather than reimplementing
either, so they cannot drift from the installed package. It checks the *dtype*
against pyseer's own `.Rtab` reader too, because "one byte per sample" and "the
platform's default integer width" are different digests and only one is right.

The digest strings in `unique_patterns.tsv` are therefore platform-dependent:
the same cohort yields different `pattern` values on a machine where numpy's
default integer width differs. The **counts** — `variants_tested`,
`unique_patterns`, `redundant_variants` — are not, because equal vectors hash
equally whatever the width. Only the recorded strings move.

Its line count is pyseer's own *tested* count — **not** the number of variants
in the input matrix. Variants excluded by `--min-af`/`--max-af`, and variants
stopped by the prefilter gate, never reach the file.

That line count is the **input to the reduction**, not the divisor. The
threshold is divided by `unique_patterns` — the number of *distinct* patterns,
which is the entire reason the reduction exists, since two variants with the
same presence/absence vector are one test. `variants_tested` and
`unique_patterns` are reported as separate rows for exactly this reason.

Note also that pyseer writes a pattern for every non-prefilter variant *before* it evaluates its own p-value gate
(`pyseer/__main__.py`), so `--filter-pvalue` and `--lrt-pvalue` do **not**
narrow this file; the frequency bounds are what bound it.

`PyseerEngine.run` cross-checks the line count against the
`%d tested variants` count pyseer writes to stderr in `pyseer_pass1.log`, and
refuses the run if the two disagree. A log without that line is not a
disagreement and leaves the reduction standing.

**`unique_patterns.tsv`**
`pattern`, `n_variants`. One row per surviving distinct pattern, in the
byte order `LC_ALL=C sort -u` produces, so a reader who sorts the patterns file
themselves gets a match rather than a puzzle. The `n_variants` column
reconciles: the row count is the unique-pattern count, and the column sums back
to the number of tested variants.

**`unique_patterns_summary.tsv`**
`metric`, `value`. The audit record for the threshold: `method`, `alpha`,
`variants_tested`, `unique_patterns`, `redundant_variants`, `threshold`,
`helper_threshold`, `patterns_file`, `survivors_file`, `count_helper`. Every
value needed to recompute `threshold = alpha / unique_patterns` is here,
including which file the count came from and which helper produced it.

`helper_threshold` is pyseer's own arithmetic as printed by the vendored
`count_patterns.py`, which formats to three significant digits. It is a
cross-check on this repository's arithmetic, not the value handed to pyseer;
the run fails if the two disagree at that resolution.

## Knowledge tables

Editable reference data. Adding a locus is a table edit, not a code change.

**`antibiotics.tsv`** — `antibiotic` (unique), `antibiotic_class`,
`mechanism_scope`, `phenotype_source_standard`, `notes`

**`mechanisms.tsv`** — `mechanism`, `gene`, `gene_type`, `antibiotic`,
`biological_role`, `evidence_level`, `mechanism_class`, `gene_role`,
`claim_ceiling`, `notes`

- `antibiotic` accepts a comma-separated list or the keyword `all`.
- `evidence_level` and `claim_ceiling` must be one of the five
  `ClaimStatus` values. **`CAUSAL` is rejected at load time** with a
  `ConfigError`.
- An unknown gene is **not** silently passed through. `mechanism_for_gene`
  raises `UnknownGeneError` naming the gene, so the table gets extended
  deliberately.

**`regulators.tsv`** — `gene` (unique), `mechanism`, `gene_role`,
`screen_for` (comma-separated variant classes), `promoter_screen`,
`reference_length`, `evidence_level`, `notes`

`regulators` is **not a stage**. `spec.md:351` folds it into the rule that owns
it, alongside `structural_variants`, `mechanisms` and `integration`; the owner is
`amr`, stage 4. Its output table is `regulators/regulator_variants.tsv`
(declared above), and `config/regulators.tsv` is the **input** knowledge table,
not a stage artefact. It carries no stage number, and the header of the file says
so. Stage 6 is `variants` and is a different thing.

**`references.tsv`** — `reference_id` (unique), `tool`, `tool_version`,
`database`, `database_version`, `version_status`
(`pinned`/`unpinned`/`na`), `source`, `notes`

**`known_determinants.tsv`** — `determinant` (unique), `kind`, `mechanism`,
`mechanism_class`, `claim_ceiling`, `members`, `source`

- **Derived, not edited.** Regenerated from `mechanisms.tsv` by
  `papipeline.downstream.known_determinants.write_known_determinants(config)`;
  `tests/unit/test_downstream_config_tables.py` fails if the committed file
  drifts from the derivation. `source` names the row it came from per row.
- Read by `read_known_determinants` with `required_columns=COLUMNS` and
  `unique_columns=("determinant",)`, so a missing column or a duplicate
  determinant is a typed refusal rather than a silent partial table.
- This is what "already known" means for the positive-control gate, the
  conditional/stratified covariates and the novelty filter. `kind` decides
  whether detecting a bare determinant counts as recovery (`gene`, `composite`,
  `family`) or whether only a feature *of that state* does (a stem, such as
  `oprD`).

**`interaction_pairs.tsv`** — `tier`, `antibiotic`, `feature_a`, `feature_b`,
`notes` (all five required by the reader; the first four used for the joins)

- `tier` is an integer; **tier 3 is refused at load time**, named, not skipped
  (it is skipped by design, and a tier-3 row is an instruction nobody gave).
- `feature_a`/`feature_b` must resolve against `known_determinants.tsv` or the
  two declared layer composites (`PDC_high`, `oprD_off_any`); an endpoint that
  resolves to nothing is refused by name.
- `antibiotic` is a comma-separated list — a pair runs only when it is relevant
  to the run's own drug (the run resolves that as
  `stages.gwas_features.target_antibiotic`).
- Only `tier = 1` lives in this file. Tier 2 is constructed per run as
  known × novel, capped by `downstream.tier2_max_pairs` in `science.yaml`.

## Output files

Written under `results/<mode>/intermediate/stages/`.

| File | Grain |
|---|---|
| `01_validation.tsv` | one row per sample |
| `03_mlst.tsv` | one row per typed sample |
| `04_amr.tsv` | one row per determinant call |
| `05_mechanisms.tsv` | one row per mechanism call |
| `06_regulators.tsv` | one row per variant call |
| `07_structural_variants.tsv` | one row per SV call |
| `08_virulence.tsv` | one row per virulence detection |
| `10_similarity.tsv` | a **square** matrix, `sample_id` leading; `distances` is a patristic matrix in **substitutions per site**, not SNPs |
| `12_gwas.tsv` | one row per tested feature |
| `13_convergence.tsv` | one row per determinant |
| `14_cooccurrence.tsv` | one row per feature pair |
| `15_master_table.tsv` | one row per sample |

### Stage 6: `min_ireads` is 1, and why the tool's default of 2 was wrong here

`variants.bcftools.mpileup.min_ireads` is **1**, and it reaches the command as
`-m 1`. This is not a tuning preference. The tool's own documented default is 2,
and this repository verified that from the installed binary rather than from
memory — `bcftools mpileup` has no `--help`, so the authoritative text is the
usage block a bare invocation prints (`bcftools 1.23.1`, usage line 53):

```
  -m, --min-ireads INT   Minimum number gapped reads for indel candidates [2]
```

**It is not a read-depth filter.** It counts *gapped reads* supporting an indel
candidate, and 2 is a **reads** threshold — which no assembly input can ever
satisfy, because a draft assembly is a single consensus sequence. There is
exactly one gapped read per indel, not two. Left at the default, every assembly
indel is discarded before it is ever considered and stage 6 then reports the
locus as clean. Assembly alignments are depth 1, so the threshold is 1.

The same depth-1 fact governs the rest of the mpileup settings, and points the
other way: `-q`/`-Q` at 20 is **permissive**, because there is no depth to filter
on and a strict threshold would return nothing. QUAL is consequently constant
(30.4183 across every called site on the measured cohort) and carries no
per-site information, so no threshold may be derived from it.

`min_ireads: 1` is **necessary and not sufficient**. It makes every
single-read length difference an indel candidate, and on an assembly alignment
most of what it admits is contig breaks and repeat misresolution rather than
biology. No plausibility filter is applied, and choosing one is a scientific
decision for the study.

**Provenance is written beside the calls**, in `variants_provenance.json`: one
entry per isolate with `n_snv` / `n_indel` / `n_alleles` — **every key present at
zero**, so "called nothing" and "never measured" cannot look alike — plus a
cohort total recomputed from the per-isolate entries so the two cannot drift. An
unwritable sidecar path is a WARNING, not a stage failure: losing the receipt
must not lose the result. And a REAL run whose **cohort total** is zero indels
logs a WARNING naming `min-ireads` as the thing to check, because ten of ten
isolates with no indels on input that is by construction ten draft assemblies
each carrying its own indels is the signature of a filter that cannot see them.
The condition is the cohort total, not per isolate: one isolate with no indels is
ordinary.

### The caller was compared, not chosen on faith

Stage 6 calls with **`bcftools mpileup` + `bcftools call`**. The assembly-native
alternative, `paftools.js call`, was run over the same alignment as a check, and
the concordance is low enough to be worth stating rather than hiding:

| comparison (PDT000292995.1, allele pairs) | shared | mine only | theirs only |
|---|---|---|---|
| indels: paftools.js vs bcftools `-m 1` | **7** | 133 | 279 |
| indels: paftools.js vs bcftools default | 0 | 140 | 6 |
| SNVs: paftools.js vs bcftools `-m 1` | **8781** | 0 | 51816 |

**Concordance on indels is 7 of 140 (5%).** The SNVs are a perfect nested
subset — all 8,781 of paftools.js's are in bcftools's 60,597, none of its own.
Two independent callers agreeing on 100% of substitutions and 5% of indels is a
statement about indels, not about either tool.

Two limits on that number, both of which matter more than the headline:

- **`paftools.js called nothing at all in `oprD`** — 0 records in
  1043800–1045500, SNVs included. It reports 958,498 reference bases covered by
  exactly one contig out of 6,264,404, about 15%. Its zero there is a **coverage
  artefact, not a verdict**: it cannot adjudicate the oprD indels at all, and
  reporting "paftools.js also finds no oprD indel" would be false.
- The concordance is **genome-wide, not oprD**, for the same reason.

The 7 shared indels are the large ones (a 919 nt and a 794 nt insertion among
them). The two callers disagree almost completely on 1–2 nt indels, which is
exactly the size class the `oprD` findings are in.

### Cohort subsets

`cohort.subset_file` names a plain-text list of `sample_id`s, one per line, `#`
comments permitted. It is `null` in `config/science.yaml` — "the whole roster" —
and every production overlay inherits that, because **which** isolates are in the
cohort is a scientific claim rather than a machine preference. A member list is
sample-level information, which is why the value is nullable and why the file it
names lives outside configuration entirely. `config/machines/smoke.yaml` is the
one overlay that names its own, and it says why on the record.

A relative `subset_file` is resolved against the **repository root**, not the
working directory and not through the `data/`/`db/` symlinks — because a
Snakemake rule's shell command runs with whatever CWD the scheduler chose, and
because in a dev worktree those symlinks point at the canonical checkout, which
would make the list this worktree reads a property of another worktree's disk.
Absolute paths are returned unchanged. The smoke overlay's list lives under
`local/`, which is gitignored, so it cannot be committed:
`git check-ignore -v local/smoke_isolates.txt`.

Plus `run_manifest.json` (provenance) and `results/<mode>/intermediate/stages/figure_data/*.json`
(the prepared data for each of the thirteen figures).

### The master table

Columns, in order:

```
sample_id, antibiotic, phenotype, amr_gene, amr_variant,
chromosomal_mutation, regulator, mechanism, structural_variant,
mlst, lineage, virulence_profile, gwas_feature, gwas_status,
convergence_status, confidence, evidence_notes
```

The many-valued columns (`amr_gene`, `regulator`, `mechanism`,
`structural_variant`, `virulence_profile`) hold a comma-separated list, not a
single "winner". Choosing one determinant per sample would be an interpretive
act with no defensible basis.

`confidence` is **derived**, not inherited: `SUPPORTED` requires a
GWAS-significant feature *and* carriage across independent lineages;
`ASSOCIATED` requires either; otherwise it falls back to the strongest direct
observation. It is never `CAUSAL`.

`evidence_notes` records what was observed and what was not concluded,
including the fixed marker `detection_is_not_resistance`.

A candidate structural variant appears in `structural_variant` with the prefix
`candidate:` and drives `confidence` to `PREDICTED`. It is never presented as
confirmed.

### The downstream evidence table: `reports/downstream_evidence.tsv`

Written by `papipeline.downstream.runner.run_downstream_analyses`, which
`run.py` calls from its `convergence` branch — **REAL mode only**. This is the
one output in this file that does *not* live under
`results/<mode>/intermediate/stages/`: it is not a stage artefact, has no stage
number (`spec.md:351` folds no such step into any stage), and sits beside the
pipeline report. The path is resolved from `config.reports_root(mode)`.

Grain: **one row per feature graded by either scan** — the union of the baseline
(stage 12) and conditional scans, baseline order first.

Columns, in order (`papipeline.downstream.report.REPORT_COLUMNS`):

```
feature, novelty, evidence_tier, claim_status, basis,
conditional_adjusted_p, baseline_adjusted_p, independent_lineages, known_matches
```

- `claim_status` is a member of `ClaimStatus` — `DETECTED`, `PREDICTED`,
  `ASSOCIATED`, `SUPPORTED`, `UNKNOWN` — and **never `CAUSAL`**; `build_report`
  refuses anything else.
- `evidence_tier` is `A`–`D`. A feature no scan tested grades `D` / `UNKNOWN`
  rather than borrowing a p-value from anywhere; a missing adjusted p is not
  significance.
- `independent_lineages` is stage 13's count when it classified that
  determinant, otherwise the count derived from the feature's own lineage
  distribution, excluding `unknown`.
- `known_matches` is a `;`-separated list (`write_tsv` joins sequences).
- Written through `write_tsv`, so `read_tsv` with `REPORT_COLUMNS` round-trips
  it; `tests/integration/test_downstream_wiring.py` reads it back.

**Not written, by design:** the tier-1, tier-2 and meta-analysis sections
travel on the `DownstreamReport` object. They have their own tables and no
column order has been fixed for them; emitting them would be a format decision,
not a wiring one. That gap is recorded in
`.build/downstream-stats.registry-notes.md` §6.

**Gated:** the file is written only if the positive-control gate recovered every
control in `downstream.positive_controls`, and the gate runs before any of the
seven steps. On the committed TEST fixtures the gate cannot pass (no MBL
feature exists in that cohort and `gene__oprD_LoF` is at adjusted p = 0.901),
which is why the wiring is REAL-only — the measurement is pinned by
`tests/integration/test_downstream_wiring.py::TestTheGateIsFirst`.

## Adding a data contract change

1. Update the table above.
2. Update `read_tsv(...)` call sites to pass the new `required_columns`.
3. Add a test in `tests/unit/` that feeds a file violating the new rule and
   asserts the typed error.
