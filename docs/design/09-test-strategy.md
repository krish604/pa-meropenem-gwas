# Design 9 — Test strategy

Companion to [1 architecture](01-architecture.md).

## 1. Two rules, non-negotiable

1. **No existing test is deleted, skipped, or weakened.** Baseline is
   **641 passing** and the full suite runs after every migration step.
2. **Every claim in the design docs has a test.** A design document that
   says "a failed genome must not invalidate its siblings" is a claim, and
   claims are only worth what their tests are worth.

The second rule matters more than usual here, because this repo has already
been bitten by the failure mode it is guarding against: a stage that looked
complete while returning degenerate data. Two real examples from this
codebase, both found only by running the real thing:

- `bakta_outputs` returned Bakta's `.inference.tsv` instead of the main
  annotation table, so every gene name was lost and all 11 stage-6 loci
  scored 0/65 — a uniform false negative that read exactly like a
  biological finding.
- A leaked loop variable made all 65 per-sample VCFs byte-identical,
  producing a core-SNP alignment in which every sequence was the same and
  a tree with no resolution.

Neither was caught by a unit test, because both only manifest against real
data at scale. That is the argument for the integration layer below.

## 2. Layers

| Layer | Location | Scope | Speed |
|---|---|---|---|
| L0 existing | `tests/unit/`, `tests/integration/` | the engine, unchanged | ~14 s |
| L1 app unit | `tests/app/test_*.py` | store, resources, tools, registry, provenance | fast, no subprocess |
| L2 engine↔app | `tests/app/test_engine_*.py` | stage loading, resume, disk-only equality | seconds |
| L3 CLI | `tests/app/test_cli_*.py` | `main(argv=[...])`, exit codes, `--json` | seconds |
| L4 Snakemake | `tests/app/test_snakemake_integration.py` | DAG, profiles, invalidation, resume | minutes, marked |
| L5 web | `tests/app/test_web_api.py` | every route, parity with the CLI | seconds |
| L6 real-data | `tests/app/test_real_cohort.py` | the 65-genome cohort, opt-in | hours, opt-in |

L6 exists because of the two bugs above. It is skipped unless
`PAPIPELINE_REAL_DATA` points at a cohort, because it needs real data and
real tools.

## 3. L0 — the existing suite

Untouched. Conventions to match, taken from `tests/conftest.py`:

- class-based (`class Test*`), no top-level test functions;
- session-scoped fixtures for `pipeline_root`, `config`, `test_data_root`,
  `manifest`, `intermediate_root`, `phylogeny_dir`, `phenotype_dir`;
- factory fixtures `write_tsv_file` and `tsv_header` for on-disk inputs;
- `test_data_root` **auto-generates** fixtures via
  `papipeline.testing.generate()` if absent, so the committed `test_data/`
  tree is deletable and reproducible;
- synthetic fixtures carry `# SYNTHETIC TEST DATA - NOT BIOLOGICAL RESULTS`
  banners; new fixtures must too.

Reuse `papipeline.testing.synthetic` for every app test that needs data.
Never hand-write a fixture the generator already produces.

## 4. Coverage the brief requires

| Requirement | Test | Asserts |
|---|---|---|
| project creation | `test_projects.py::test_create_writes_the_documented_layout` | all seven directories exist; `project.yaml` valid |
| | `::test_create_is_idempotent_on_name` | second create with the same name is a `ConflictError` |
| | `::test_delete_requires_confirm` | without `confirm`, nothing is removed |
| input validation | `test_validation.py::test_empty_fasta_rejected` | 0 bytes → per-sample error, project not created |
| | `::test_truncated_fasta_reported_not_fatal` | corrupt sample is *excluded with a reason*; the rest proceed |
| | `::test_duplicate_sample_id_rejected` | `validate_sample_id` rules still enforced |
| | `::test_selection_precedes_metadata` | sample set derivable from files alone; phenotype cannot influence it |
| | `::test_phenotype_values_restricted` | a value outside `R,I,S,SDD,ND` is rejected |
| | `::test_no_mic_invented` | phenotype TSV with no MIC leaves `mic is None` |
| | `::test_unmatched_phenotype_sample_reported` | orphan phenotype row is a warning, not a silent drop |
| project persistence | `test_store.py::test_project_roundtrip` | create → read → update → delete |
| | `::test_cascade_delete_removes_samples_and_runs` | FK cascade |
| | `::test_database_is_rebuildable_from_disk` | delete the DB, rebuild, same state |
| tool discovery | `test_tools.py::test_version_floor_rejects_old_tool` | 0.1.19 samtools → unavailable |
| | `::test_resolution_prefers_version_over_path_order` | the Homebrew 1.24 wins over the env's 2014 stub |
| | `::test_missing_tool_is_reported_not_raised` | a missing tool is a registry row, not a crash |
| resource detection | `test_resources.py::test_detect_cores_and_memory` | reads the real machine, no hard-coded values |
| | `::test_plan_never_exceeds_cores` | `workers × threads ≤ cores − 1` |
| | `::test_plan_never_exceeds_memory` | `workers × mem ≤ 80% of RAM` |
| | `::test_hpc_profile_does_not_set_cores` | `cores is None` on a cluster profile |
| | `::test_every_plan_carries_a_reason` | `AnnotationPlan.reason` is never empty |
| Bakta scheduling | `test_bakta_scheduler.py::test_shared_database_used_once` | one DB path across all jobs; never copied per genome |
| | `::test_completion_detected_from_real_filenames` | input-named outputs count as complete — the `bakta.gff` bug |
| | `::test_prune_keeps_what_the_pipeline_reads` | gff/tsv/inference/log survive; json/embl/gbff/svg do not |
| | `::test_prune_does_not_invalidate_a_sample` | `is_complete()` true after prune |
| | `::test_retry_clears_partial_output` | a killed run's half-written GFF does not read as complete |
| | `::test_database_incompatible_with_tool_is_refused` | the AMRFinderPlus 2024-vs-2026 case |
| resume | `test_resume.py::test_completed_stages_are_skipped` | resume runs only PENDING/FAILED |
| | `::test_complete_but_missing_output_reruns` | the filesystem, not the DB, decides |
| | `::test_finished_run_resume_is_a_noop` | |
| checkpointing | `test_checkpointing.py::test_completed_genomes_not_rerun` | 61/200 done → 139 scheduled |
| | `::test_state_survives_kill` | state file valid after SIGKILL mid-sweep |
| failed jobs | `test_failed_jobs.py::test_one_bad_genome_does_not_stop_the_sweep` | the brief's explicit requirement |
| | `::test_missing_binary_is_a_record_not_an_exception` | `run_one` never raises |
| | `::test_failed_job_does_not_invalidate_siblings` | |
| | `::test_run_fails_if_any_job_failed` | and is `COMPLETE` when all are `COMPLETE` or `SKIPPED` |
| Snakemake | `test_snakemake_integration.py::test_dag_builds_in_test_mode` | |
| | `::test_dag_refuses_real_mode` | mode gate holds at DAG-build time |
| | `::test_full_run_matches_native_orchestrator` | byte-for-byte |
| | `::test_touching_gene_families_invalidates_amr` | the current `CONFIG_FILES` omission |
| | `::test_kill_mid_annotation_resumes` | only pending genomes rescheduled |
| | `::test_execution_order_matches_engine` | `EXECUTION_ORDER` and the DAG agree |
| GUI/API | `test_web_api.py::test_every_verb_has_a_route` | |
| | `::test_route_set_equals_cli_verb_set` | **the parity guarantee** |
| | `::test_error_codes_map_to_http` | `ValidationError`→422, `NotFoundError`→404, … |
| | `::test_concurrent_run_is_409` | |
| provenance | `test_provenance.py::test_run_records_everything_required` | inputs, tool versions, DB versions, config, resources, commit, timestamps, commands |
| | `::test_provenance_explains_the_report` | a report figure traces to a provenance row |
| | `::test_provision_is_the_only_network_verb` | no verb downloads during a run |
| result discovery | `test_results_registry.py::test_every_declared_type_discoverable` | all 14 `ResultType`s |
| | `::test_registry_fingerprints_headers` | schema change is detectable |
| | `::test_report_omits_nothing_the_run_produced` | a report cannot silently skip a stage |

## 5. Tests for the engine changes (migration step 2)

The highest-risk change, so the most explicit tests:

```python
@pytest.mark.parametrize("stage", [
    "mechanisms", "gwas", "convergence", "cooccurrence", "integration",
])
def test_stage_from_disk_equals_stage_in_process(stage, ...):
    """The disk path and the in-process path must agree exactly.

    Step 2 adds Optional disk loaders to six stages. If they disagree with
    the in-process path, the Snakemake DAG will silently produce different
    results from the native orchestrator - which is the exact failure the
    DAG's own header warns about.
    """
```

Plus: every stage is runnable **alone** from disk-only inputs, and
`run_pipeline(mode="TEST")` output is unchanged byte-for-byte.

## 6. L6 — the real-cohort regression suite

Opt-in via `PAPIPELINE_REAL_DATA=/path/to/cohort`. Slow, and the only layer
that would have caught the two degenerate-data bugs.

The assertions are deliberately about *data shape*, not about biology:

- a core-SNP alignment has **more than one distinct sequence** (catches the
  identical-VCF bug);
- per-sample VCFs have **more than one distinct fingerprint** (catches the
  leaked loop variable before it reaches the alignment);
- the identity-by-state kinship matrix is **not all 1.0** (catches the
  called/not-called mask bug);
- annotation output for N genomes is **not byte-identical** across
  samples;
- `oprD` is detected in a plausible fraction of genomes, and a locus
  reported 0/N is **flagged rather than silently accepted**.

None of these assert a biological conclusion. They assert that the
machinery is producing distinct, non-degenerate data — which is the class
of bug this repository has actually hit, three times.

## 7. CI

No CI exists today. Proposed, in priority order:

1. `pytest tests/ -q` on every push — 641 existing + L1–L3, L5.
2. `pytest tests/ --profile-snakemake` — L4, on a machine with Snakemake.
3. `python -m papipeline_app.tools.preflight` as a separate job; its output
   is an artifact, so drift in the toolchain is visible in the diff.
4. L6 nightly, opt-in, on a self-hosted runner with the real cohort.

Version pins: `environment.yml` for `osx-arm64` reality, plus the
constraints files from design 8 step 8.

## 8. What is deliberately not tested

| Not tested | Why |
|---|---|
| GUI visual appearance | No screenshot tests. The parity test covers behaviour; appearance is not a correctness property. |
| Snakemake internals | Snakemake's own test suite covers its scheduler. We test *our* DAG: that it builds, that it matches the engine, and that it resumes. |
| Real biology | Asserted nowhere. The pipeline's job is to report what the data says, and the L6 shape assertions are the boundary of what this repo claims. |
| Network behaviour | `tools provision` is tested against a local fixture server, never a live repository. |
