# Handoff: ticket 14 — variant calling, recombination masking, similarity matrix

Current state, not a history. Prior reasoning is in `git log`; the per-decision
back-and-forth behind anything here is recoverable from the commits named below.

## 1. State

- Branch `feature/imipenem-pa-pipeline`, tree clean. The commit that last
  rewrote this file is the newest on the branch — read `git log -1` rather than
  a hash recorded here, which is stale the moment this file is edited.
- **992 unit, 233 integration, 0 failures, 0 skips**
- `SkillOpt/` is untracked, unrelated, and unresolved — do not touch

**Green now:** D5 reference verifier + pin; `_overlay_root`; the single
`locate_assembly`; the fifteen-stage taxonomy; the STUB runtime; the join;
`consolidate_results.py` no longer states an unmeasured census.

**Proven working against real data:** `locate_assembly` resolves the first 25
real accession directories under `data/`.

**Counted but never constructed:** `PDC_essential.tsv` — 967 rows, 835 with
assemblies. No manifest has been built from it, because the function that would
does not exist.

## 2. THE ONE BLOCKER: two stages have parsers but no executor

**`variants` and `cohort_variants` can read and filter data but cannot do the
job.** They sit in `UNBUILT_STAGES` and refuse outside STUB, and removing that
refusal makes the suite fall to 30 failures because there is nothing behind it.

`run.py::_execute_stage` is the dispatcher. Its `variants` branch (`run.py:721`)
only writes the regulators internal table — its own comment reads *"Unreachable
until ticket 14 builds the stage"* — and it never calls the per-isolate caller.
There is **no `cohort_variants` branch at all**.

| piece | state |
|---|---|
| `stages/variants.py` — VCF → rows | done, contract derived from a real run |
| `stages/cohort_variants.py` — merge + symmetric filter | done, validated at N=4/12/34 |
| `run.py` dispatch calling either | **absent** |
| minimap2 + bcftools invocation | **absent everywhere** |

### Build this first

1. **Two `_execute_stage` branches**, red-first: align+call per sample writing
   `variants.tsv`, and the merge writing `cohort_variants.tsv`. Decide what each
   does in STUB, where stages fabricate rather than execute. This is real work,
   not bookkeeping.
2. **Then** the reconciliation, which is small and already proven to work:
   drop both from `UNBUILT_STAGES`; drop `variants` from
   `standins.UNBUILT_WITH_STANDIN`; add both to `run.py`'s stage imports. The
   fact currently lives in three places and they all disagree.
3. **Then** the narrow-scope `allow_real_mode` test, and only then the flag.
4. **Then** the real run on the 10 prepared genomes.

TEST fixtures already exist (`test_data/variants/`, 5 VCFs covering every filter
outcome) so step 3 is safe once 1-2 land. `db/smoke_genomes/` holds 10 assemblies
already proven usable end-to-end this session — gitignored, not committed.

## 2b. Then, in order

5. **gubbins** — not installed; solve into a separate py310 env (269 pkgs)
6. **Similarity matrix** — IQ-TREE *is* installed, so cheapest of the three
7. **Wire `gwas.run()` and `phylogeny`** to the real files
8. **The bounded REAL smoke run**

## 3. What each item needs

### 1 — pilot reports read the real census

- The census **is already persisted**: `run_pilot100.py:306` writes
  `pilot100_prepare_summary.json` with `exclusion_reasons` (the
  `exclusion_breakdown` output) plus verdict and analysable counts. The reports
  need to *read* it; no plumbing is required.
- The prose must be **restructured, not constant-swapped**. It currently walks
  through how 835 became 65 (233 zero-byte, 51 truncated, 426 corrupt, 93
  intact, 28 removed). A different run may have a different *set* of categories,
  so the paragraph has to describe whatever the map actually contains.
- Red-first with a fabricated `pilot100_prepare_summary.json` whose
  `exclusion_reasons` has a deliberately different shape — different counts, and
  ideally a category absent from the 835 run. No real pilot run is needed.
- **The keys are machine identifiers, not prose.** `exclusion_breakdown`
  (`pilot/cohort.py:556`) returns `Dict[str, int]` keyed by a reason-type prefix:
  `corrupt_binary_data`, `undersized`, `oversized`, `zero_byte_file`,
  `file_not_found`, `not_fasta_no_header`, and `other` for anything unmatched.
  A fixture using "truncated"/"corrupt" can never match real output. Presenting
  the map means deciding how identifiers become readable text, and handling a
  key absent from that list — which is exactly what `other` is for.

### 2/3 — the manifest

- Cohort membership is **all 967 isolates**, not the 835 with assemblies. An
  isolate with no assembly is an ordinary member; a stage needing a sequence
  refuses per sample via `locate_assembly`, and that refusal is correct, not a
  bug to work around.
- `PDC_essential.tsv` is read **directly** — do not populate
  `data/metadata/sample_metadata.tsv` as a copy.
- Seam: `papipeline/pdc.py` (`load_pdc_index(path, key)`, `clean_value(value)`,
  generic key) plus
  `papipeline/manifest.py::discover_pdc_manifest(pdc_path, *, assembly_paths=None)`.
  `assembly_paths` is how the smoke run supplies explicit per-sample paths —
  rule 1 of `locate_assembly`, used with no pattern matching.
- `pilot/cohort.py` imports the promoted functions; it keeps working unchanged.
- **`load_pdc_index` cannot be reused as-is**: it is keyed on `Assembly` and
  skips rows without one, so it holds 835 entries. Index by `Isolate` and keep
  `Assembly` as an attribute.
- The small-N test proves the same logic path serves 3–5 and 967. It is the
  first thing written and everything else is red against it.

### 4 — `config/machines/smoke.yaml` (answers settled)

- `smoke_genome_dir` as a config value in the overlay.
- Genome-source resolution **requires** it and refuses by key name if absent —
  no fallback to `data/` or any default.
- `max_samples: 10`, enforced by the existing `enforce_sample_cap`. No new cap
  logic.
- Report marker, **asserted** so it fails loudly if it cannot be attached:
  `SMOKE TEST — REAL data, N=<actual count> assemblies, not the full cohort
  analysis.`

### 5 — bcftools

- `bcftools` 1.23.1, `minimap2` 2.31, `samtools`, `mafft`, `iqtree` 3.1.3 are
  installed. `snippy` and `freebayes` are not.
- Design the contract from actual `mpileup`+`call` output. The current
  `variants.tsv` columns are marked `PROVISIONAL` — replace them, do not match
  them.
- Pin the exact invocation in `science.yaml`; not "latest".

### 6 — gubbins

- A solve into a **separate** env succeeds (`gubbins 3.4.3 py310hdfa5cb7_1`,
  269 packages). Nothing is installed.
- Record the environment requirement in **`science.yaml`**, not a machine
  overlay: whether masking must happen for the result to mean anything is the
  science's, while availability is the machine's. A run that silently skips
  masking produces a phylogeny that looks valid and is not.

## 4. Standing rules

- **Red-first**, via the `tdd` skill. Agree the seam before writing its test.
- Small commits; **full suite green after each**; never start the next item
  while red.
- **No `git checkout` on a modified file.** Branch checkouts are fine. After a
  verified fix, delete or rename the pre-fix backup immediately — stale backups
  have been nearly restored twice, and timestamp naming did not prevent it.
- `data/`, `PDC_essential.tsv` and the committed fixtures in `test_data/` are
  read-only.
- `allow_real_mode` stays closed. Item 4 *builds* the overlay; it is not
  exercised until the smoke run.
- Tickets 16, 08 and 25 are not started.
- No invented checksums, versions or tool flags. A value must be derived or
  downloaded.
- The pipeline must stay cohort-size-agnostic: nothing may assume an isolate
  count. Verified dynamic so far — `discover_manifest` iterates every row and
  logs `len(manifest)`; no slice or expected count anywhere in `manifest.py`.

## 5. History

`git log` for the reasoning behind each decision. The commits this document
supersedes: `3246562` (this rewrite's predecessor), `2dc53cd`, `1764d1f`,
`fa04663` (census removed from `consolidate_results.py`), `b112f01` (single
assembly lookup), `df4cb32` (overlay root), `a0e8bdd` (D5 verifier and pin).
