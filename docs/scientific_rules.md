# Scientific rules

These are the interpretive constraints this pipeline enforces. They are not
aspirational: each one is implemented in code and covered by a test, and the
file names the enforcing module and the test that would fail if it regressed.

The vocabulary of claims is `papipeline.models.ClaimStatus`: `DETECTED`,
`PREDICTED`, `ASSOCIATED`, `SUPPORTED`, `UNKNOWN`. **There is no `CAUSAL`
member**, so a causal claim cannot be expressed in the output schema at all.

---

## 1. Gene presence is not phenotypic resistance

A detected determinant is reported as `DETECTED` and nothing more. Detection
establishes that a sequence is present; it does not establish that the
sequence functions, is expressed, or confers resistance in this background.

**Enforced by** `stages/amr.py` (every call is constructed with
`claim_status=DETECTED`), `knowledge.ceiling_for()` (clamps any request to a
gene's configured `claim_ceiling`), and `stages/integration.py`
(`detection_is_not_resistance` is written into every record's
`evidence_notes`).

**Tested by** `test_input_claim_is_clamped`,
`test_all_ceilings_are_detected`, `test_amr_rows_are_only_detected`.

### The OprD case, specifically

OprD is the principal outer-membrane porin for carbapenem entry. Its *loss* is
the relevant observation, not its presence. The pipeline therefore treats the
two as different things:

| Observation | Reported as |
|---|---|
| `oprD` annotated (stage 2) | mechanism `locus_intact`, evidence `DETECTED` |
| `oprD` disrupted (stage 6) | mechanism `reduced_permeability` |
| `oprD` absent (stage 6) | mechanism `reduced_permeability` |
| insertion in `oprD` (stage 4, confirmed) | mechanism `reduced_permeability` |

Stage numbers here are `spec.md` D1's, which `papipeline/run.py:STAGE_ORDER`
mirrors. `structural_variants` is not a stage number of its own: `spec.md:351`
folds it into the rule that owns it, stage 4 (`amr`), and `run.py:150-152` says
so. Pangenome is stage 7, recombination stage 8, phylogeny stage 9 and
similarity stage 10.

**Enforced by** `stages/mechanisms.py`: `LOSS_OF_FUNCTION_REQUIRED` and
`annotate_intact_locus()`, whose note reads "presence indicates an intact
locus and is not evidence of reduced permeability or of susceptibility".

**Tested by** `test_oprd_presence_is_locus_intact_not_permeability`.

#### What `disrupted` means, and what it does not

`disrupted` is a claim about the **coding sequence**, and it has exactly two
admissible grounds:

- a **frameshift** — an indel whose length is not a multiple of three, so
  everything downstream is read out of frame; or
- a **premature stop** — a stop codon inside the open reading frame, before the
  native stop.

That is the whole definition. In particular:

- **An in-frame deletion that preserves the native stop is NOT loss of
  function.** Removing a multiple of three bases leaves the reading frame intact
  and leaves the terminator where it was, so the remaining sequence is still
  translated normally. It is a shorter protein, not a broken one. It must not be
  classified `disrupted`. This is the single easiest way to manufacture `oprD`
  loss of function from real data, because in-frame deletions are common in
  *P. aeruginosa* `oprD` and routinely leave the stop codon in place.
- **There is no ORF-fraction cut in this classification.** No threshold on "how
  much of the ORF is missing" participates in deciding `disrupted` vs intact. A
  locus retaining 95% of its coding sequence with its native stop is intact; a
  locus retaining 99% with a frameshift is disrupted. Fraction-of-ORF is not the
  quantity being classified, and introducing a cut on it would be a new
  scientific threshold, not a clarification of this one. The old banding could
  not have this property: under a 0.90 cut, PDT000292995.1 at 0.8378 was `intact`
  and at 0.85 `disrupted`.

**Aligned span is used only for presence versus absence.** How much of the locus
the alignment covers decides one thing: whether the locus is there. A span
covering the gene means present; no span means absent. Span length does not feed
the `disrupted` decision and must not be read as a severity. The two thresholds
that do exist — minimum identity to be considered the orthologue at all, and
minimum reference coverage before "is it whole?" is answerable — are both
**presence** thresholds, and `sensitivity_table` sweeps those two and nothing
else.

**Truncation length is reported as an attribute, never used as the test.** Where a
truncating variant is detected, its length is carried through to the output so a
reader can see what was lost. It is descriptive. It is not the classifier, and it
is not compared against any cutoff to decide the mechanism.

Note that this is a statement about what the classification does *not* do. It is
deliberately silent on intermediate cases — an in-frame deletion that *removes*
the native stop, or a start-codon loss — because resolving them requires a
threshold this rule does not set.

**Enforced by** `adapters/oprd_locus.py`: `DisruptionKind.DISRUPTING` is
`{FRAMESHIFT, PREMATURE_STOP, BOTH}` and nothing else;
`LocusEvidence.stop_is_reference_terminator` is what separates a stop that is the
gene's own ending from one that is a lesion; `orf_length_fraction` exists as a
reported attribute and is documented as unable to move a verdict.

**Tested by** `test_the_orf_length_cannot_move_a_verdict`, which sets
`orf_aa_length` to fourteen values from 1 aa to 1000 aa and asserts the verdict
and lesion type never change.

The synthetic fixture deliberately contains resistant samples with no
detected determinant and susceptible samples with one, so the distinction is
exercised rather than assumed.

---

## 2. Association is not causation

The pipeline reports statistical association. It does not report, imply or
infer a causal relationship, and the word "causal" does not appear in any
generated finding.

**Enforced by** the absence of a `CAUSAL` member in `ClaimStatus`; and by
`config/mechanisms.tsv`, which rejects `CAUSAL` in the `evidence_level` column
at load time with a `ConfigError`.

**Tested by** `test_causal_evidence_level_rejected`,
`test_confidence_never_causal`.

---

## 3. A GWAS association is not a causal claim

A significant GWAS p-value is evidence of an association in this cohort and
under this model. It is not evidence of a mechanism, and it is not evidence
that the feature would still be associated in another cohort.

Two further limits are reported rather than hidden:

- The reference engine performs **no kinship correction**, so its p-values are
  anti-conservative in a structured cohort. Every result it emits is stamped
  `reference_fisher:NO_KINSHIP_CORRECTION`.
- `I`, `SDD` and `ND` are **excluded**, not folded into `R` or `S`. The number
  excluded is reported as a run warning.

**Enforced by** `stages/gwas.py`: `GwasConfig.positive/negative/excluded`,
`benjamini_hochberg()` applied to every feature, and the `model` field on every
`GwasResult`.

**Tested by** `test_intermediate_categories_are_excluded_not_guessed`,
`test_model_records_the_missing_kinship_correction`,
`test_multiple_testing_correction_is_applied`.

### 3a. The multiple-testing threshold is derived in its own step, from pyseer's own helper

The number of independent tests is the number of **distinct variant
patterns**, not the number of variants: two variants with identical
presence/absence vectors are one test, not two. The threshold is

```
threshold = alpha / n_unique_patterns
```

This is Bonferroni, and it is Bonferroni because that is the one correction
pyseer documents. Its own `scripts/count_patterns.py` computes exactly this,
from the file `--output-patterns` writes, and its entire logic is
`LC_ALL=C sort -u <patterns> | wc -l` followed by `alpha / count`.

Two properties follow, and both are enforced rather than intended:

- **The count comes from pyseer's helper, not from this repository.** pyseer
  1.1.2 does not package that script — the wheel's `entry_points.txt` declares
  exactly five console scripts and it is not among them — so it is **vendored**
  at `scripts/gwas/count_patterns.py`, byte-identical to upstream, with its
  provenance, the SHA-256 of the upstream file and the SHA-256 of its own logic
  in the file's header. The logic digest is recorded *in the file* as well as
  in the test suite, so editing the correction and updating the recorded digest
  are two separately visible changes rather than one. Writing a local
  correction instead would have been the alternative, and it would have been a
  silent change of the science for a number the report presents as pyseer's.
  `config.gwas.unique_patterns.method` therefore accepts only `bonferroni`.
- **The step is separate from the association and its inputs are on disk.**
  Step 12a runs pyseer once with both p-value gates stated as open, reduces
  the patterns, and writes `patterns.txt`, `unique_patterns.tsv` and
  `unique_patterns_summary.tsv`; step 12b runs pyseer again with
  `--lrt-pvalue` set to the threshold 12a derived. Splitting them is what makes
  the number auditable: a user can read how many patterns were dropped and check
  the count against `sort -u | wc -l` themselves.

Three corollaries that are easy to get wrong, and are therefore checked:

- `alpha` is a scientific constant and lives in `config/science.yaml`. It is
  **required**, with no default — the helper has its own `0.05` to fall back
  on, and inheriting it would let a run report a threshold as though someone
  had chosen it. `method` is required for the same reason.
- **What the divisor is, and what bounds it.** The divisor in
  `alpha / n_unique_patterns` is `n_unique_patterns` - the number of
  *distinct* patterns. That is the whole point: two variants with the same
  presence/absence vector are one test. The count of variants *behind* it is
  pyseer's *tested* count, and that is the input to the reduction, not the
  divisor. Confusing the two is easy and changes the threshold, so both are
  reported separately in the audit record and neither is called the other.
  Variants outside the frequency bounds, and variants stopped by the prefilter
  gate, never reach the patterns file. Note also that pyseer writes a pattern
  for every non-prefilter variant *before* it evaluates its own p-value gate
  (`pyseer/__main__.py`), so `--filter-pvalue` and `--lrt-pvalue` do **not**
  narrow the file; those flags are stated explicitly so a future change to
  pyseer's defaults cannot either, and the frequency bounds in the shared
  command are what bound the pattern set.

  `PyseerEngine.run` also cross-checks the file's line count against the
  `%d tested variants` count pyseer writes to stderr. **This check can only
  refute, never establish**: if the line is absent - pyseer emitting it
  somewhere other than stderr, or a caller capturing stderr elsewhere - the
  check logs and stands down, and the reduction proceeds on the other two
  cross-checks. `PYSEER_TESTED_COUNT_PATTERN` is pinned against the installed
  package by `test_the_tested_count_string_is_pyseers_own`, because a check
  that fails open is only honest if you know it is one.
- The helper's `--memory`, `--cores` and `--temp` bound a **separate `sort`
  process**, so they are not this pipeline's `runtime.memory_mb`/`threads`.
  They come from the machine overlay, and reading them raises rather than
  falling back, because each has an upstream default (1024 Mb, 1 core, `/tmp`)
  that would otherwise be inherited silently. The scratch directory and the
  patterns path are both validated before being handed over, because the helper
  interpolates both into one unquoted `shell=True` string.

Review history for this section — including two claims it made that turned out
to be wrong — is on ticket 16 rather than here. A normative document should
state what is true now; where a claim was contested and how it was settled
belongs with the work.

**Enforced by** `stages/gwas.py`: `reduce_unique_patterns()`,
`count_unique_patterns()`, `tally_patterns()`, `parse_count_patterns_output()`,
`PyseerEngine.run()` (two passes, plus `_cross_check_pyseer_tested_count()`),
`_require_helper_safe_path()`, `_require_parallel_sort()`,
`GwasConfig.unique_patterns`,
`PipelineConfig.unique_patterns_memory_mb` / `_cores` / `_temp_dir`.

**Tested by** (all in `tests/unit/test_gwas_unique_patterns.py`):

*Formats checked against the installed pyseer rather than reimplemented*
`test_the_encoding_is_the_one_pyseer_builds_for_an_rtab`,
`test_the_pattern_encoding_is_pyseers_not_a_packed_byte_vector`,
`test_the_pattern_is_one_line_and_carries_the_newline_pyseer_depends_on`,
`test_the_tested_count_string_is_pyseers_own`,
`test_a_different_vector_gives_a_different_pattern`,
`test_the_synthetic_association_table_has_pyseers_own_column_count`,
`test_a_header_with_no_rows_is_an_empty_result`,
`test_the_synthetic_table_rows_come_from_pyseers_formatter`,
`test_every_column_the_stage_renames_is_one_pyseer_emits`,
`test_the_columns_the_stage_requires_are_the_ones_it_renames`

*The vendored helper is upstream's*
`test_the_logic_is_byte_identical_to_the_revision_it_was_taken_from`,
`test_the_header_records_the_same_digest_the_test_checks`,
`test_the_slicing_the_digest_is_taken_over_is_pinned`,
`test_it_is_vendored_at_the_path_the_stage_calls`,
`test_the_header_names_the_upstream_source_and_its_whole_file_digest`,
`test_it_actually_runs_as_a_script_and_reads_its_flags_from_its_own_parser`,
`test_the_helper_really_shells_out_to_sort_u_wc_l`,
`test_it_reproduces_sort_u_wc_l_on_a_synthetic_pyseer_patterns_file`

*The arithmetic and the audit record*
`test_the_threshold_is_alpha_over_unique_patterns`,
`test_alpha_is_read_from_config_not_hardcoded`,
`test_both_keys_are_required_rather_than_defaulted`,
`test_it_reports_the_survivor_count_and_what_was_dropped`,
`test_the_survivors_file_lists_one_row_per_unique_pattern`,
`test_the_survivors_file_reconciles_against_sort_u_wc_l`,
`test_the_survivors_file_records_how_the_threshold_was_derived`,
`test_the_summary_is_written_and_is_readable_tsv`,
`test_it_reruns_cleanly_over_its_own_outputs`

*Cross-checks, and the runs they stop*
`test_a_helper_that_disagrees_with_the_independent_tally_stops_the_run`,
`test_a_helper_whose_own_threshold_disagrees_with_ours_stops_the_run`,
`test_a_disagreement_stops_the_run`,
`test_a_disagreement_further_down_a_full_log_is_still_caught`,
`test_a_loaded_count_that_differs_does_not_trip_the_check`,
`test_a_log_without_the_line_leaves_the_reduction_standing`,
`test_a_missing_log_is_not_a_disagreement`,
`test_agreement_is_the_normal_case`,
`test_the_count_is_among_the_lines_pyseer_reports_not_derived_from_them`

*Refusals*
`test_a_missing_patterns_file_is_named_as_a_missing_intermediate`,
`test_an_empty_patterns_file_is_refused_rather_than_divided`,
`test_a_path_the_helper_cannot_safely_receive_is_refused_on_the_public_path`,
`test_a_scratch_directory_the_helper_cannot_safely_receive_is_refused`,
`test_a_shell_metacharacter_in_the_scratch_directory_is_refused`,
`test_a_memory_below_the_helpers_reserve_is_refused_naming_the_key`,
`test_more_cores_than_this_machines_sort_supports_is_refused`,
`test_more_cores_is_accepted_where_sort_supports_it`,
`test_a_missing_scratch_directory_is_refused_with_the_key_to_fix`,
`test_a_truncated_pyseer_output_is_refused_not_read_as_no_associations`,
`test_a_completely_empty_pyseer_output_is_refused_not_read_as_no_associations`,
`test_every_branch_of_normalise_rewrites_its_target`

*Configuration, read rather than inherited*
`test_every_committed_overlay_declares_all_three_helper_limits`,
`test_the_helper_limits_come_from_the_overlay`,
`test_an_overlay_that_omits_a_helper_limit_refuses_when_it_is_read`,
`test_a_minimal_overlay_still_loads`,
`test_a_config_with_no_machine_overlay_refuses_rather_than_defaulting`,
`test_a_config_without_the_unique_patterns_block_refuses_to_load`,
`test_only_the_correction_pyseer_documents_is_accepted`,
`test_an_alpha_outside_zero_to_one_is_refused`

*Helper-output parsing*
`test_it_reads_the_two_line_report`,
`test_a_zero_count_is_refused_not_divided`,
`test_a_changed_output_shape_is_refused_rather_than_guessed_at`,
`test_a_missing_line_names_which_one`,
`test_an_extra_keyed_line_is_refused_not_ignored`,
`test_a_repeated_key_is_refused_rather_than_one_of_the_two_taken`,
`test_a_warning_line_before_the_report_is_refused_not_swallowed`,
`test_a_non_numeric_count_is_refused`,
`test_the_helpers_own_arithmetic_agrees_with_ours`,
`test_duplicate_patterns_are_the_ones_that_collapse`,
`test_a_blank_lines_only_patterns_file_is_also_refused`

*The threshold reaches pyseer*
`test_pyseer_is_asked_for_patterns_and_no_p_value_filter_in_pass_one`,
`test_pass_two_carries_the_threshold_pass_one_derived`,
`test_the_frequency_bounds_not_the_p_value_gates_bound_the_pattern_set`,
`test_a_reduction_failure_stops_before_any_association_table_is_parsed`,
`test_no_variant_passing_the_threshold_is_an_empty_result_not_a_failure`,
`test_the_result_records_that_the_threshold_was_bonferroni_on_patterns`,
`test_pyseer_column_names_are_mapped_off_a_synthetic_pyseer_table`,
`test_the_parsed_feature_is_not_one_pyseer_never_reported`

---

## 4. Lineage association must be distinguished from resistance association

A feature carried only by samples from one lineage cannot be separated from
that lineage by any association test, no matter how small its p-value. Such a
feature is **lineage-linked**, not resistance-associated.

The pipeline computes a `dominant_lineage_share` for every GWAS result, flags
any feature above `config.gwas.lineage_confound_threshold` (default 0.9), and
reports those features in a dedicated report section, separate from the
association table. A flagged feature is not deleted; it is relabelled.

**Enforced by** `stages/gwas.py`: `is_lineage_confounded()` and
`flag_lineage_linked()`, which emit a warning naming the affected features;
`stages/reporting.py` renders them under "Lineage-linked features".

**Tested by** `test_single_lineage_feature_is_confounded`,
`test_evenly_spread_feature_is_not_confounded`,
`test_lineage_features_are_flagged`.

---

## 5. Convergence must be evaluated in phylogenetic context

Whether a determinant is a candidate for convergent selection depends on how
many *independent* lineages carry it, not on how many samples do. Stage 13
classifies each determinant into exactly one category, in this order:

| Order | Category | Condition |
|---|---|---|
| 1 | `widespread_background` | carried by > `widespread_fraction` (0.75) of samples |
| 2 | `rare_isolated` | carried by ≤ `rare_max_samples` (2) |
| 3 | `lineage_associated` | all carriers in one lineage |
| 4 | `recurrent_convergent` | ≥ `min_independent_lineages` (2) independent lineages |
| 5 | `unknown` | no lineage information available |

The order matters and is documented so the categories are mutually exclusive.
Every call carries its full lineage distribution, so a reader can check the
basis for the category.

Convergence is a **pattern in the data**, not a demonstration of a shared
selective cause.

**Enforced by** `stages/convergence.py`: `classify()`, with the ordering in the
module docstring.

**Tested by** `test_recurrent_convergent`, `test_lineage_associated`,
`test_widespread_background`, `test_rare_isolated`,
`test_single_is_rare_not_lineage_associated`.

---

## 6. R/I/S must not be converted into a quantitative MIC or zone diameter

A categorical susceptibility call carries no magnitude. Inventing an MIC from
`R` would fabricate a measurement that was never made.

- `R`, `I`, `S`, `SDD`, `ND` never acquire an MIC or a zone diameter.
- `MIC` is read only from an actual measurement column.
- An unparsable `MIC` is an error, not a silent `None`.
- An `MIC` reported against `ND` is an error: `ND` means not determined, so it
  cannot carry a value.
- A zone diameter is accepted only from a real measurement column, never
  derived from a category.

**Enforced by** `stages/phenotype.py`: `parse_phenotype_row()` performs no
derivation of any kind, and raises on an MIC with `ND`.

**Tested by** `test_ric_s_carry_no_mic`, `test_mic_with_nd_raises`,
`test_unparsable_mic_raises`, `test_zone_diameter_only_from_measurement`,
`test_no_mic_fabricated_from_categories`.

---

## 7. Candidate structural variants must not be treated as confirmed

Stage 4 keeps three states distinct:

| `call_status` | Meaning | Mechanism claim | Master table |
|---|---|---|---|
| `confirmed` | meets the configured evidence requirement | `DETECTED` | plain `variant_id` |
| `candidate` | suggestive only | `PREDICTED` | `candidate:<variant_id>` |
| `not_assessable` | region could not be evaluated | **no call emitted** | absent |

`not_assessable` is **not** "no variant present", and the two are never
conflated: a sample with an unassessable region is reported as unassessable,
and the count is surfaced as a run warning.

An unrecognised `call_status` **downgrades to `candidate`**, never upgrades to
`confirmed`.

`config.structural_variants.promote_candidate_calls` is deliberately ignored by
the code. Promotion requires orthogonal evidence and is a human decision.

**Enforced by** `stages/sv.py`: `parse_call_status()` (conservative default),
`confirmed_only()`; `stages/mechanisms.py`: `interpret_structural()` returns
`None` for `not_assessable` and `PREDICTED` for `candidate`;
`stages/integration.py` labels candidates in the output and derives
`confidence=PREDICTED`.

**Tested by** `test_candidate_is_predicted_not_detected`,
`test_not_assessable_yields_no_call`, `test_unknown_downgrades_to_candidate`,
`test_candidate_sv_is_labelled`, `test_candidate_svs_are_not_promoted`.

---

## 8. The database version must always be recorded

Every AMR and virulence record carries `database` and `database_version`. Both
are mandatory columns; a row missing either fails the data contract.

Databases are **never downloaded or updated during an analysis run**. All
`allow_database_update` flags are `false` in `config.yaml`, and
`adapters.ToolAdapter.ensure_database()` raises rather than honouring a
version change mid-run. Provisioning is an out-of-band step.

An unpinned reference is a **reported failure, not a silent pass**: every
unpinned row in `config/references.tsv` raises a run warning and appears in
`run_manifest.json` with `version_status: unpinned`. Currently all 10 are
unpinned, and this is stated in every test report.

**Enforced by** `stages/amr.py` (`REQUIRED` includes both columns),
`adapters/external.py`: `ensure_database()`, `probe_version()`;
`run.build_provenance()`; `environment/versions.sh`.

**Tested by** `test_amr_records_database_version`,
`test_no_implicit_database_updates`, `test_unpinned_references_are_surfaced`.

---

## 9. Missing data must remain missing

A missing value stays missing. It is never defaulted, imputed, zero-filled or
inferred from a related field.

- Sentinel strings (`.`, `NA`, `nan`, `None`, `null`, `-`, empty) become
  `None`, never `0` or `""`.
- An unmeasurable QC metric is `None`, not `0`.
- A sample with no phenotype record is absent from the phenotype matrix, not
  defaulted to `S`.
- A sample with no MLST record maps to `None`, which is distinct from a
  partial profile and from a `no_call`.
- "Nothing detected" is a list of length zero; "not screened" is an absent key.
  These are kept distinct throughout.
- A hook that has not been run reports `None`, not `pass`.
- An empty input file is a `DataContractError`, so a truncated or failed run is
  never mistaken for a cohort with no findings.

**Enforced by** `io/tsv.py` (`_clean_cell`), `io/fasta.py` (`_nxx` returns
`None`), and the explicit empty-list / `None` distinctions in each stage.

**Tested by** `test_sentinels_become_none`, `test_short_row_pads_with_none`,
`test_metrics_are_none_for_no_records`, `test_missing_samples_are_absent_not_defaulted`,
`test_empty_tsv_raises`.

---

## 10. Unsupported biological conclusions must not be generated automatically

The pipeline reports measurements, counts, associations and their stated
uncertainty. It does not generate narrative biological conclusions.

Concretely, the pipeline will not state that a gene *causes* resistance, that a
mechanism *explains* a phenotype, that a convergence pattern *demonstrates*
selection, or that a co-occurrence pair *is* an interaction. Those
interpretations belong to a human, in a document that cites this pipeline's
outputs and their limits.

The synthetic fixtures are labelled
`SYNTHETIC TEST DATA — NOT BIOLOGICAL RESULTS` in every file header and in
every report. The statistical relationships between the generated genotype,
mechanism and phenotype columns are **artefacts of the generator seed**, stated
as such in `papipeline/testing/synthetic.py`. They exist so stages 12–15 have a
signal to process and so assertions are stable. They are not evidence about
*Pseudomonas aeruginosa*.

**Enforced by** the absence of conclusion-generating code; the test-report
banner in `stages/reporting.py`; and the generator's refusal to write into
`data/`.

**Tested by** `test_report_is_marked_synthetic`,
`test_report_states_the_scientific_rules`, `test_generator_refuses_to_write_into_data`.

---

## 11. A branch length is a number of substitutions per site, not a count of SNPs

IQ-TREE's branch lengths, and therefore every patristic distance this pipeline
derives from them, are **substitutions per site**. They are not SNP counts, and
the two are not interchangeable: a branch length depends on the length of the
alignment it was estimated from, so the same tree expressed against a different
alignment yields different numbers.

Measured on the 10-isolate cohort (`pa-artifacts/phy2/patristic_report.txt`),
the same tree's patristic distances expressed against the SNP alignment (70,660
sites) and against the full-length core alignment (4,642,925 sites) have an
**identical** `max/min` of 2.3192, while the median ratio moves from 1.6316 to
107.2110:

```
  n_sites basis n70660     =     70660
    ratio  min = 0.9496   median = 1.6316   max = 2.2023
    spread (max-min) = 1.2527   max/min = 2.3192
  n_sites basis n4642925   =   4642925
    ratio  min = 62.3956   median = 107.2110   max = 144.7058
    spread (max-min) = 82.3103   max/min = 2.3192
```

Only the scale moves with the denominator; the spread does not. That is the
signature of a per-site quantity. A SNP count would not move at all.

Nor is the patristic distance a rescaled SNP distance. Against the masked SNP
matrix the ratio ranges 0.9496 to 2.2023 (median 1.6316, Spearman 0.988407). High
rank agreement is not scale agreement, and quoting a patristic distance as
"SNPs" would misstate its magnitude by a factor that varies per pair.

**The input alignment is variable-sites-only, so ascertainment bias is real and
`+ASC` is required.** The alignment given to IQ-TREE is the core-SNP alignment,
and IQ-TREE reports on it (`pa-artifacts/phy2/b_asc/iqtree.log:19-20`):

```
Alignment has 10 sequences with 70660 columns, 894 distinct patterns
44955 parsimony-informative, 25705 singleton sites, 0 constant sites
```

`0 constant sites` is not a property of the biology; it is what
variable-sites-only means. The model is therefore estimating site frequencies
from a sample in which the invariant class cannot appear, which biases them, and
biases the branch lengths derived from them. Running the same data as
`GTR+G+ASC` instead of `GTR+G` moves the median patristic-to-SNP ratio from
1.6316 to 0.8099 — `|median - 1|` from 0.6316 to 0.1901 — at a correlation with
the uncorrected run's ratios of 0.970535
(`pa-artifacts/phy2/scaling_check.txt`). `-fconst` is **not** the equivalent: it
collapses the fit (lnL −6800510.031, ~19× worse) because it strips the
informative positions inside gapped columns along with the uninformative ones.

**This is why `+ASC` sits behind a column filter rather than being passed
straight to IQ-TREE.** `+ASC` asserts that *every* alignment column is a
potential site; a column that is constant once gaps and N are ignored makes that
assertion false and the pinned binary refuses:

```
ERROR: Invalid use of +ASC because of 20 invariant sites in the alignment
```

So `phylogeny.asc_drop_partially_constant` removes those columns **before** the
search, and `TreeResult.column_filter` records how many went, so that `None` and
"dropped nothing" stay distinguishable. On the measured masked alignment the
filter is a correctly inert no-op (70,660 in, 70,660 kept, 0 dropped) — that
alignment has no constant sites — but it is what makes `+ASC` applicable **in
general**, and it was not optional: the committed fixture's 20 invariant sites
were fatal before the filter existed.

**The path sum excludes the LCA's own stem.** A patristic distance is the sum of
the branch lengths from each tip down to their MRCA. The MRCA's own branch is
**not** in it. Adding it is a plausible-looking error with a large effect: on one
sister pair it overstates by 0.2391, i.e. **128 %**, and turns a 0.8998 ratio
into 2.052 — manufacturing exactly the `GTR+G`-like dispersion this correction
set out to remove. Three implementations are checked against each other to
1e-12 on all 44 pairs: explicit edge traversal, the depth form
`d(a)+d(b)−2·d(MRCA)`, and the pipeline's own
`stages.similarity.patristic_distances`.

**Consequences for what may be written down.**

- Stage 10's `distances` column is in substitutions per site. The stage writes
  the unit, the definition and `is_a_snp_count: false` beside the matrix in
  `similarity.units.json`, derived from the matrix's own name so the pair cannot
  drift.
- Stage 9's alignment summary must not present a variable-sites-only alignment as
  though its length were the genome's.
- Any figure that puts a patristic distance next to a SNP count must label both
  units. They are different quantities that correlate; they do not measure the
  same thing.

**Enforced by** `adapters/iqtree.py`: `drop_partially_constant_columns()` and
the `+ASC` precondition check; `stages/similarity.py`:
`DISTANCE_UNIT`, `write_units_sidecar()`, `patristic_distances()`.

**Tested by** `test_without_the_filter_the_same_alignment_is_refused` (the
counterfactual is asserted, not described), and the similarity units tests.

---

## What a reader may and may not take from an output

**May** conclude: which determinants, variants, structural variants and
virulence factors were detected, in which samples; the mechanism each maps to
under the configured knowledge table; the observed phenotype distribution; which
features are statistically associated with which category, after correction, and
with what lineage distribution; how widely and in how many lineages a
determinant occurs; which feature pairs co-occur more than expected.

**May not** conclude: that any determinant confers resistance in any sample; that
any association is causal; that any lineage-linked feature is a resistance
determinant; that convergence indicates a shared selective cause; that
co-occurrence indicates interaction; that a candidate SV is present; or that
any category corresponds to a specific MIC.
