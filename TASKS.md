# TASKS

Roll-up of `.scratch/imipenem-gwas-dashboard/issues/`, which is the tracker of
record. If this file and the ticket files disagree, **the ticket files win**.

Spec: `.scratch/imipenem-gwas-dashboard/spec.md` (read it before starting any
ticket — the decisions it records are not optional).

Status vocabulary is on the `Status:` line of each ticket file:
`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`.

## Global constraints — these apply to every ticket

- Nothing environmental is hard-coded. Paths, threads, memory, database
  locations and breakpoints come from configuration. No `os.cpu_count()`.
- Never invent a package version, tool flag or database name. Verify by search
  and by `--help`, or say you are unsure.
- Never edit a test to make it pass.
- Never commit genomes, databases or secrets.
- Sample-ID mismatches fail loudly. No fuzzy joins, no silent drops.
- **No REAL-mode run** and no reading from `data/` until the user says
  `run real samples`.
- Tests are written before the implementation for every parser, validator,
  config loader and event function.
- Changes to existing `papipeline/stages/` code go in **separate commits**,
  never folded into another ticket's commit.
- **When a data model or a key changes, check for orphaned accessors.** This is
  the repeated way a model change silently breaks consumers: an attribute is
  carried correctly, a method exists and is correct, and nothing reads it — so
  the breakage is invisible in a diff, because nothing is missing where you are
  looking. Found live five times: `MachineConfig.reference_fasta`,
  `PipelineConfig.genomes_dir`, `intermediate_root`/`reports_root` not
  delegating to the overlay, `locate_assembly` never reading `Sample.assembly`,
  and `gwas.py` never reading a configured target (issue 26). The last one
  shipped a cohort where **0 of 835 real assemblies resolved**, with every
  refusal a legitimate "assembly not found" — a quiet failure, not a crash.
- **A fact declared in one place that the code never reads or classifies by** is
  now the **sixth** occurrence of the same defect class, and it deserves a sweep
  rather than a seventh mid-task discovery. The instances: `reference_fasta`,
  `genomes_dir`, `intermediate_root_for`/`reports_root_for` never delegating to
  the overlay, the assembly-lookup break, `gwas.py` never reading a configured
  target (issue 26), and `PER_SAMPLE_STAGES` omitting `variants` even though
  spec.md D1 declares stage 6 per-sample — so the taxonomy did not know the
  stage it had just been given a per-isolate contract.
  The shape is always the same: spec, contract or config says X, a
  classification set/dict/registry never learned X, and the code is then
  *correct with respect to its own stale list*. Nothing is missing where you
  are looking, which is why a diff never shows it.
  **Recommended sweep** (not done — deliberately deferred): grep for
  classification sets, dicts and registries, then cross-check every entry
  against what the spec and the contracts actually declare. `STAGE_FLAGS`,
  `PER_SAMPLE_STAGES`, `STAGE_TABLES`, `UNBUILT_STAGES`, `PREREQUISITES` and the
  Snakefile's own rule set are known instances of the same shape and drifted
  apart from each other while all looked fine.
- **Never `git checkout` a file, for any reason, on this ticket.** Both
  incidents this session happened during mutation testing, where the instinct is
  to undo a mutation. `cp` to `/tmp` before any edit that will be mutated or
  reverted, always. A third occurrence is a process failure to raise and discuss,
  not to self-correct and continue.

## Agents

`.opencode/agents/` holds four subagents: `pipeline` (backend/stages/config),
`tdd-writer` (tests only), `dashboard-ui` (dashboard only), `reviewer`
(read-only, run at the end of every phase). The reviewer verifies claims by
running them; it never trusts a summary.

---

## Phase 0 — setup (complete)

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 01 | Agent skills, issue tracker and spec | none | **done** |

Delivered: `AGENTS.md`; `docs/agents/{issue-tracker,domain,triage-labels}.md`;
`.opencode/agents/{pipeline,tdd-writer,dashboard-ui,reviewer}.md`;
`spec.md`; this file. `.gitignore` closed for `status/`, `db/`, `databases/`,
`.env`.

---

## Phase 1 — environment

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 02 | Build the environment on osx-arm64 and report the gaps | none | **done** — `374f352`; gap list corrected by `8b382a7` |

**The single most load-bearing unknown in the project.** Everything about
whether the pipeline can run on this machine is decided here.

Delivers: a corrected `environment/environment.yml` that resolves on
`osx-arm64`; a stated list of which required tools have no arm64 build, with the
evidence; correction of the suspect pins (AMRFinderPlus package name, BLAST
version line, iqtree package); gubbins and mafft added **only if** they resolve.

Acceptance: `micromamba env create` succeeds; a plain-language table of
available / unavailable tools exists; no pin in the file is unverified.

---

## Phase 2 — configuration, identity, phenotype

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 03 | Config resolver: science + machine overlays | none | **done** — `e554b97`; path-loaded overlay resolves like a named one, `df4cb32` |
| 04 | 20-sample laptop cap | 03 | **done** — cap enforced in `run.py`, `max_samples: 20` in `config/machines/laptop.yaml`, tests in `tests/unit/test_sample_cap.py`; mode gates `379fc49` |
| 05 | Sample identity and manifest | 03 | **verified** — `b98f300`, verified end to end `19bb7ab` |
| 06 | Directional manifest↔phenotype join | 05 | **verified** — implemented `b98f300`, wired and asserted against a STUB run `a48f976`. Superseded in practice by `b112f01` (one manifest-driven lookup) |
| 07 | Phenotype parser with AST provenance | 03 | **wired** — `b98f300`; smoke phenotype given its own directory `70a54c6` |
| 08 | Breakpoint layer - **DEFERRED; thresholds stay empty** | 07 | **deferred by decision** — see "Judgment calls settled" below. `max_maf` resolved `3c48442` |

**03** splits machine-independent science from per-machine settings so laptop and
big-machine results cannot diverge on science. Deliverable: a resolver plus a
test that proves the two overlays cannot disagree about a scientific constant.

**04** is the laptop safety rail. Deliverable: exceeding the cap stops the run
with a message naming the limit, the config key and file, and what to change.

**05** makes GCA accession the canonical `sample_id`, makes the ID format rule
configuration-driven so both real and synthetic IDs validate, and carries PDC
isolate/BioSample/isolate-name as validated attributes.

**06** implements the asymmetric join. A manifest genome with zero or two
phenotype rows is a hard failure. A phenotype row with no assembly is excluded.
An assembly whose phenotype is `I`/`SDD`/`ND`/absent is excluded with a reason.
Exclusion counts are reported. This is the ticket most likely to be got wrong by
accident, so its tests are the most important in the project.

**07** keeps source R/S/I authoritative, preserves genuinely measured MICs, and
carries AST method/standard/edition per row while reporting how many rows lack
each.

**08** writes the CLSI M100 structure with an explicit `UNVERIFIED` placeholder
and makes stage 11 **hard-fail naming the missing field**. It must be impossible
for an unsourced threshold to reach a result. Blocked pending real values from a
licensed M100 copy.

---

## Phase 3 — orchestration, STUB mode, events

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 09 | Event emitter | none | **wired** — `b98f300`, exercised by a STUB run `424b030`. **Not consumed by the dashboard**: `ff8f94a`, `f184cae`, `7ad7fdc` |
| 10 | STUB mode | 03, 09 | **done** — `424b030` |
| 11 | 15-stage Snakefile, threads, memory, tool gating | 10, 04, 05 | **done** — Snakefile fixed `efba315`, taxonomy folded to D1 `86ff513`, cascade asserted in the suite `0e615a2` |
| 12 | `snakemake -n` with both machine configs | 11 | **done** — DAG check made hermetic and proven non-vacuous, `697ff52` |
| 13 | Full stub run end to end | 12 | **done** — `424b030`, asserted `a48f976` |

**09** is the event contract, deliberately early and with no blockers: the
dashboard depends on it and nothing else does. One JSON object per line with
exactly `t`, `event`, `stage`, `sample`; `event` is one of `start`/`done`/`fail`;
aggregate rules use the literal sample `all`; append-only; a malformed line is
skipped with a warning rather than crashing the stream.

**10** makes the whole DAG runnable with no tools installed and no fixture on
disk. This is the mode that lets the orchestration and the dashboard be developed
on this machine at all.

**11** writes the fifteen rules, one per stage, each declaring threads and
memory, mapping the existing sixteen-stage code onto them, and adding the
DAG-build check that a rule needing a tool unavailable on this machine fails
immediately with a message naming the tool.

**13** is the first end-to-end proof: a complete STUB run producing every stage's
output and a well-formed event log covering all three event values.

---

## Phase 4 — real-tool logic that STUB mode skips

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 14 | Variant calling vs PAO1, verified reference, pinned loci | 11, 08 | **wired, not verified on real data** — reference verified and pinned `a0e8bdd`, index guard `55db12f`, invocation pinned from a measured run `d5063c3`, per-isolate caller `24ff924`, dispatch `f5573f5`. A REAL run has not reached stage 6 |
| 15 | OprD loss-of-function features | 14 | **not built** — no `oprD_absent` / `oprD_LoF` producer exists. `stages/mechanisms.py` interprets an already-detected determinant and documents that LoF must arrive from stage 6 or 7; nothing produces it |
| 16 | GWAS input: three families + unique-pattern reduction | 14, 15 | **unique-pattern reduction (12a) built**; families still partly built — `GwasInput` with families, lineages and exclusions exists in `stages/gwas.py`; step 12a reduces the variant matrix to unique patterns via pyseer's vendored helper, writes survivors and an audit record, and derives the threshold handed to `--lrt-pvalue`. Refuses REAL rather than use a mechanics-test engine, `f90d6c9`. Blocked on 15 for the OprD family |
| 17 | Convergence and combination analysis | 16 | **wired in TEST, refuses REAL by design** — `46f5005`. Now declared in `REAL_REFUSING_STAGES`, `239db00` |

**14** is the first stage that genuinely needs a tool. It verifies the
reference's actual length and checksum against the pinned values and fails on
mismatch — a truncated reference silently shifting every variant coordinate is
invisible until the biology is wrong. Calls are restricted to `oprD` and the
pinned regulators.

**15** produces **two** features, not one: `oprD_absent` (locus wholly absent)
and `oprD_LoF` (present but carrying a frameshift, premature stop or internal
deletion). They are mechanistically distinct; collapsing them destroys the
distinction between a porin that is gone and a porin that is present and broken.

**16** is the ticket the reviewer will scrutinise hardest. It must prove the GWAS
input carries all three families — SNP/indel, gene presence/absence, OprD
features — and that the multiple-testing threshold is derived from unique variant
patterns in a separate inspectable step. Unitigs are **not** among the families.
pyseer's helper script name and flags are confirmed from the installed package's
`--help`, not from memory.

**The 12a half of 16 is built.** pyseer's own `scripts/count_patterns.py` is
**vendored verbatim** at `scripts/gwas/count_patterns.py` — pyseer 1.1.2's wheel
declares five console scripts and that file is not among them, so it was never
installable. Its header records the upstream URL, retrieval date, the upstream
file's SHA-256 **and** the SHA-256 of the vendored logic; the test suite
recomputes the second and asserts it matches, so an edit to the correction fails
CI and an edit to the recorded digest without the correction fails CI too.
`alpha` and `method` are both **required** in `science.yaml` — the helper has
its own `0.05` to fall back on, and inheriting it would let a run report a
threshold as though someone had chosen it. The helper's
`--memory`/`--cores`/`--temp` bound a separate `sort` process and come from the
machine overlay, raising rather than defaulting when absent; both the patterns
path and the scratch directory are validated, because the helper interpolates
both into one unquoted `shell=True` string. The threshold is
`alpha / n_unique_patterns`, and pyseer now runs **twice**: once to record a
pattern for every variant it tested, then again with `--lrt-pvalue` set to the
derived threshold. Three outputs land in the GWAS work directory —
`patterns.txt`, `unique_patterns.tsv`, `unique_patterns_summary.tsv` — so the
count reconciles against `sort -u | wc -l` by hand. They are deliberately **not**
in `execution/contracts.py`: they are work-directory artefacts, not a stage's
principal table, and naming them there would imply a consumer that does not
exist.

**Four consequences worth recording.** First, the count is cross-checked three
ways — against an independent tally of the same file, against pyseer's own
printed threshold, and against the `%d tested variants` count pyseer writes to
stderr. A disagreement stops the run rather than picking a winner. Second, the
divisor is the **unique-pattern** count, and pyseer's **tested** count is the
reduction's *input* — two numbers it is easy to confuse, which changes the
threshold. The tested count is pyseer's own, not the input matrix's variant
count: pyseer writes a pattern for every non-prefilter variant *before* it
evaluates its own p-value gate (`pyseer/__main__.py`), so `--filter-pvalue` and
`--lrt-pvalue` do **not** narrow the patterns file — an earlier draft of this
file and of `docs/scientific_rules.md` §3a claimed they were what protected it,
and both have been corrected. What protects it is the frequency bounds, plus
the stderr cross-check. Third, a header-only association table is now an empty
result rather than a `DataContractError`, since pass 12b runs at a derived
cutoff and "no variant passes" is the ordinary null outcome — but a single
*non-header* line is refused, because reading it as "no associations" would turn
a broken run into a reported negative result. Fourth, **the families half of 16
is not done**: ticket 15 still has no `oprD_absent` / `oprD_LoF` producer, and
the SNP/indel family reaches pyseer through an `.Rtab` synthesised from the
gene presence/absence matrix, which is narrower than pyseer's `--vcf` route.
Both are upstream of the reduction and were left alone deliberately.

**17** is the stage that makes hits interpretable. It must require a determinant
to appear in at least two independent lineages before calling it convergent,
because a public cohort over-represents a few high-risk clones.

---

## Phase 5 — dashboard

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 18 | Dashboard backend: replay and stream | 09 | **built** — `39462da`, enabled for smoke `2829df8`, run-key read from the overlay and mismatches reported `791599e` |
| 19 | Dashboard frontend | 18 | **built** — `papipeline/observatory/static/{index.html,css,js}` |
| 20 | Performance at 900 isolates | 19 | **built** — `papipeline/observatory/throughput.py`; the 900-isolate measurement is **not recorded** |
| 21 | Dashboard shown live against a stub run | 13, 19 | **not done** — this is open item 5 of the definition of done |

**18** replays the whole event log on start to rebuild state, then tails and
streams over SSE, binds to `127.0.0.1` only, and survives a truncated final line.
Deliberately separate from the existing observatory — see the spec's evidence
section, where the decisive fact is that the observatory's store holds one
cohort-level row per stage, so it structurally cannot back a 900×15 grid.

**19** delivers metric cards, stage cards, one square per isolate per stage, and
a live event log, with colour defined once in CSS custom properties. Explicit
"not reported" wherever the stream carries no fact — never a plausible zero.

**20** is measured, not asserted: replay time and per-event update cost against
a realistically-sized synthetic log. One event must update one cell, not 13,500.

**21** is the demo: the dashboard tracking a live STUB run.

**OPEN — observatory `event_sequence` is stuck at 0 across process boundaries.**
Not a defect and not fixable inside the observatory alone, but it means
"observatory verified" is only half true. The state machine persists and reads
back correctly — `/api/tasks/{stage}` returns real state, timestamps and
elapsed time from the SQLite store, and a REAL smoke run recorded `amr
INCOMPLETE` at 1986s, verified by 791599e. The *event stream* does not cross the
process boundary: `EventBus` is in-process, so a separately launched server
reports `event_sequence: 0` and an empty `/api/events`, and `/api/snapshot`
reports `elapsed_seconds: 0.0` for a stage whose store row carries 1986.0.
`event_sequence` can only exceed zero when the server runs inside the pipeline's
own process.

This stays open until the dashboard is verified against a **live cross-process
run**, not just the state machine: a REAL run whose events a separately launched
server actually receives and streams. Until then, treat any observatory
verification as covering persistence only. See AGENTS.md on the two
observability mechanisms not meeting.

---

## Phase 6 — documentation, lock, sign-off

| # | Ticket | Blocked by | Status |
|---|---|---|---|
| 22 | Pilot cohort and REAL-mode manifest builder | 05, 06 | **done** — smoke overlay with three refusals `76312e5`, one manifest-driven lookup `b112f01`, cohort from `PDC_essential.tsv` `c791d1d`, Isolate keying fixed `8f6dc4a` |
| 23 | README verified from a clean folder; Linux lock | 13, 21, 02 | **partly built** — `README.md` and `environment/environment-linux.yml` exist; **not verified from a clean checkout**. Blocked on 21 |
| 24 | Reviewer sign-off on the definition of done | all | **not done** — open item 7 of the definition of done |

**22** produces the artefacts for a 10–20 genome, seeded, phenotype-stratified
cohort that records its own selection rationale, plus the manifest builder that
joins PDC metadata to the 835 downloaded assemblies. It builds the artefacts
**only**; running against them stays gated.

**23** documents the laptop and big-machine paths and is verified by following
it from a clean checkout, not by reading it. Also produces the locked Linux
environment file.

**24** is the final gate. The reviewer checks all seven items of the definition
of done and reports PASS / FAIL / **UNVERIFIED** with evidence for each.

---

## Blockers, as of 239db00

Three, and they are not the same kind of problem. Two are environment walls
that no amount of pipeline code gets past on this machine.

### BLOCKER 1 — structural variants: no REAL caller. **Unblocked, next up.**

Stage 4 folds `structural_variants` in (spec.md:351). `stages/sv.py` refuses REAL
with `NotImplementedError`, `a50bc0f`. A REAL smoke run therefore completes
stages 1-4 and stops: `Stage amr ended INCOMPLETE`, with `04_amr.tsv` written
(16,505 bytes) before the refusal.

This is the current stopping point and it is ordinary code work. It is also the
highest-value unblock available, because it is the only one of the three that
can be finished on the laptop. Nothing else stands between here and a REAL run
reaching stage 5.

### BLOCKER 2 — gubbins segfaults on osx-arm64. **Parked.**

`pa-gubbins` has gubbins 3.4.3 installed (an arm64 Mach-O binary), so the build
succeeds; the binary is what fails. It dies with **SIGSEGV (exit 139)** on a
3-sequence, 10-column alignment, with a zero-byte log, so neither input size nor
memory is involved.

Note for the next reader: `run_gubbins.py` reports this as `Gubbins crashed,
please ensure you have enough free memory`, because `subprocess.check_call`
raises `CalledProcessError` on a signal and that is caught as a
`SubprocessError`. That message is wrong for a segfault and will mislead anyone
who trusts it.

Parked pending a production-machine OS answer or an explicit go-ahead. Blocks
`recombination` (the last stage in `UNBUILT_STAGES`).

### BLOCKER 3 — panaroo on osx-arm64. **RESOLVED. Not an architecture wall.**

This entry previously read "panaroo has no osx-arm64 build. Parked, pending a
bigmachine decision." **That was wrong**, and it is corrected here rather than
quietly removed.

1. panaroo **does** have native `osx-arm64` builds (1.6.0, 1.7.0, 1.8.0), and
   `intbitset` has a working arm64 build too. It was never an architecture wall.
2. What blocked it was **packaging metadata**. The recipe carries a hard `prokka`
   dependency, and prokka cannot solve on arm64: `tbl2asn-forever` is in no
   channel, the `perl =5.26.2` pin is dead, and older panaroo wants `mkl`, which
   has no arm64 build. A `--platform osx-arm64` solve pulls **zero** osx-64
   packages, which is what proves it is metadata and not architecture.

The prokka dependency is **spurious**. `panaroo/prokka.py` is a GFF3 *reader* with
no subprocess call to prokka. panaroo's only external binary is cd-hit, which
installs cleanly on arm64 (4.8.1, verified).

**Resolution:** panaroo is installed **from source via pip**, pinned to an exact
commit in `environment/environment.yml`, with cd-hit declared as a conda
dependency. Verified end-to-end against unmodified real Bakta GFF3 for the
10-isolate smoke cohort — exit 0, 10,019 gene families, 4,828 core, 5,191
accessory, 0 genome-specific. `config/machines/laptop.yaml` now declares panaroo
available. Full evidence in `docs/environment-arm64.md` section 3.

**What this unblocked:** stage 7 (pangenome) runs on this machine. It was the
wrong blocker to have on record — it was recorded as an architecture ceiling that
would need a bigmachine decision, when the actual fix was a source install
available all along.

**What it did not unblock:** stage 8 (recombination) is still impossible here.
gubbins installs on arm64 and then segfaults — exit 139, `EXC_BAD_ACCESS
address=0x0` — so the downstream ceiling described in the old text still holds,
for a different and now accurately-recorded reason. See BLOCKER 2.

---

## Judgment calls settled

Recorded, not re-opened.

| Question | Resolution | Commit |
|---|---|---|
| `max_maf` ceiling | Measured to N=34; ceiling not resolvable by sampling, so made symmetric with the floor: two isolates on each side. Floor is evidence, ceiling is a design choice | `48a5f7d` → `3c48442` |
| Stage taxonomy | Three states, not two: runnable-with-a-REAL-caller, runnable-in-TEST-but-refusing-REAL (`gwas`, `convergence`, `cooccurrence`), and unbuilt. `REAL_REFUSING_STAGES` names the second; the taxonomy is now self-verifying against the dispatch chain and stage ASTs | `239db00` |
| Join / phenotype source rules | A bounded run gets its own phenotype directory; the full cohort's `data/phenotype/` is never written by a smoke run. One manifest-driven assembly lookup replaces three that disagreed | `70a54c6`, `b112f01` |

---

## What is left, stated plainly

Code written: 12 of 16 stages are REAL-capable. Verified end to end on real
data in a single pass: **4** (validation, annotation, mlst, amr). Stages 5-16
have never run together on real data - which is exactly where the last three
defects lived, so treat the remaining wiring as unproven rather than
mechanical.

Never started: a real smoke report containing the marker, and the full-cohort
run.

---

## Definition of done: met / not met

1. ~~The environment builds.~~ **MET** — `374f352`.
2. ~~All pytest tests pass.~~ **MET** — 1771 passed, 1 skipped at `239db00`.
   Must be run from an activated environment; invoking the env's `python`
   directly leaves external tools off `PATH` and fails 14 tests spuriously.
3. ~~`snakemake -n` passes with both machine configs.~~ **MET** — `697ff52`.
4. ~~A full STUB-mode run passes.~~ **MET** — `424b030`, asserted `a48f976`.
5. The dashboard shows that STUB run live. **NOT MET** — ticket 21. See the
   `event_sequence` entry above: persistence is verified, the cross-process
   event stream is not.
6. ~~The README is verified from a clean folder.~~ **NOT MET** — ticket 23. The
   README exists but has not been followed from a clean checkout; blocked on 21.
7. The reviewer agent signs off. **NOT MET** — ticket 24.

**4 of 7 met, not 5.** Items 6 and 7 are both open, not just 7.

---

---

## Phase 7 — every declared external has a real REAL-mode producer

**Logged before the work, deliberately.** Not a queue item: a project phase with
its own estimate, because the work is the test *plus* draining what it finds.

### Why this exists

Four separate defects in this project were the same shape: a path that a rule or
a stage depended on, which had a producer in TEST and **none in REAL**.

| Path | Consumers | Symptom |
|---|---|---|
| `intermediate/regulators/regulator_variants.tsv` | rule `variants` (`REGULATOR_IN`) | DAG unbuildable in REAL |
| `annotation/<sid>.annotation.tsv` | `load_from_intermediate` | `n_records = 0` for all 967 cohort members, reported as success |
| `intermediate/amr/amr_determinants.tsv` | rule `amr` (`AMR_IN`), `load_amr_table` | DAG unbuildable; ten real reports unread |
| `data/standins/unbuilt_stages/variants.tsv` | rule `variants` | pointed into read-only `data/` |

Each was found by accident, while chasing the previous one. That is the argument
for the check: not that four bugs existed, but that **nothing was looking for
this class**, so they surfaced one at a time through failed runs.

`tests/integration/test_snakefile_rules_consistent.py` closes the rule-declaration
half - stand-in references, unproduced rule inputs, stale docstrings. It cannot
close this half: the remaining defects live in **Python stage code**, which the
rule layer does not inspect. `amr_determinants.tsv` passed every rule-level check
while no stage wrote it.

### The work

**The test, ~1 day.** For each `intermediate/` path a rule declares as an
external, resolve the owning stage module and assert a REAL code path writes it.
Most of the effort is the resolver: mapping a declared external to the stage that
owns it, since `INTERNAL_TABLES` and the Snakefile constants are two independent
spellings of the same location.

**Draining what it finds, 2-3 days.** This is the larger half and the reason
this is its own phase. Known instances already fixed (`annotation`, `amr`) are
done; the remaining candidates are `structural_variants.tsv`, `mlst_results.tsv`,
`gwas_features.tsv` and `virulence_factors.tsv`, none of which has been
confirmed either way. Each finding is behavioural work, not a test fix - the
test lands red on purpose.

**Do not fold this into unblocking the DAG.** Two of the four defects above were
found while unblocking it, and fixing them inline is correct. The systematic
sweep is a different activity: it produces a work list rather than progress on
the run.

### Acceptance

Every `intermediate/` path any rule declares as an external is written by some
REAL stage, or the declaration is removed. No new instance of the four-row table
above can be introduced without the suite failing.

---

---

## Build state is not "done"

A ticket's `Status:` line now carries a build state next to the triage label:
**`built`** (implemented and unit-tested, but nothing calls it), **`wired`**
(reachable from the path a run actually takes), **`verified`** (wired, and an
end-to-end test exercises it on realistic input). See
`docs/agents/triage-labels.md`.

This was added because the review found three units that were finished,
well-tested and unwired while their tickets read "done":

| Unit | State at the end of Phase 2 | Wired by |
|---|---|---|
| sample cap | **verified** - now enforced before stage 1, proven by a 21-sample run | Phase 2, this pass |
| sample-id pattern + provenance checks | **verified** - `discover_manifest` applies them, and an end-to-end test shows a violation stops the run before anything is written | Phase 2 wired, Phase 3 verified |
| directional join (`papipeline.join`) | **verified** - wired into the phenotype branch | Phase 3 |

| continuous trait in the GWAS model | **built** - the model still takes binary R-vs-S | ticket 16, not started |
| SNP calling vs PAO1 (`variants`) | **not built** | ticket 14; spec D1 stage 6. No stage in `STAGE_ORDER`. |
| recombination masking (`recombination`) | **not built** | spec D1 stage 8 (gubbins). No stage in `STAGE_ORDER`; gubbins is also unrunnable on this machine (py39/py310 vs 3.11). |
| similarity matrix (`similarity`) | **not built** | spec D1 stage 10 (tree distance). No stage in `STAGE_ORDER`. |

### Phase 3, first pass: the workflow loads, and has never run

**Superseded by the table below.** Kept because it records the state this
phase started from: at that commit the workflow had never executed, and the STUB
runtime did not exist. Two rows in it are now false - do not read them as
current.


| Unit | State | Note |
|---|---|---|
| `workflow/Snakefile` | **built** | Loads, and the DAG shape is **verified** by a hermetic dry run against a throwaway results root. **The workflow has never actually executed** - every check is `snakemake -n`. |
| stage input roots | **verified** | The rules pointed at `intermediate_root`; the stages are passed `tool_output_root` (`loader.py:357`), a deliberately different directory. Corrected, and the dry run proves it. |
| `stage_inputs()` | **verified** | Splits rule-produced deps from on-disk files. TEST/REAL keep external inputs so a missing file fails at DAG-build. |
| stage reachability | **verified** | `pangenome` and `phenotype` were orphaned - nothing depended on them, so a green dry run said nothing. Both rewired; a mutation test proves the check bites. |
| `scripts/run_rule.py` | **built** | Every stage rule routes through it and its unit tests pass, but no Snakemake *run* has invoked it. |
| STUB runtime | **not implemented** | `papipeline/run.py` has no `RunMode.STUB` branch. `resolve_mode` raises `NotImplementedError` rather than accept a mode that would read real files. |

### Phase 3, later: STUB runs, and what it exposed

| Unit | State | Note |
|---|---|---|
| STUB runtime | **verified** | One real `snakemake` run produces 20 stage outputs, 2 reports, and a 32-record event log (16 start + 16 done). `papipeline/stub.py` fabricates; no parser, tool or fixture is touched. |
| `scripts/emit.py` | **wired** | Has a production caller (`run_rule.py`, called by every stage rule) and its output is asserted on a real run. But **it has no consumer**: nothing in `papipeline/` reads `events.jsonl`. The dashboard serves `/api/events` from an in-process `EventBus`. Verified as a writer, not as anything observed. |
| `papipeline.join` | **verified** | Wired into the phenotype branch. A manifest genome with no phenotype row is now a hard failure, not a warning. Exclusions are written to `11_phenotype_exclusions.tsv`. |
| `scripts/run_rule.py` | **verified** | Every stage rule routes through it, proven by a real run. |

### Phase 3, this pass: the taxonomy fold, and the stand-ins

| Unit | State | Note |
|---|---|---|
| stage taxonomy (fifteen) | **verified** | `spec.md:349-353` mapped the sixteen onto the fifteen; landed as one atomic commit because every boundary between the pieces is red. `STAGE_TABLES` 15, `INTERNAL_TABLES` 4, `STAGE_ORDER` 15, membership matching spec D1. |
| the four folded bodies | **verified** | `structural_variants` -> `amr`, `mechanisms` -> `cooccurrence`, `integration` -> `reporting`, `regulators` -> `variants`. Modules and tests unchanged; their outputs moved to `INTERNAL_TABLES` so the writer and the contract cannot disagree. |
| the three unbuilt stages | **verified** | `variants`, `recombination`, `similarity`: declared, fabricated in STUB, refused outside it with `NotImplementedError` naming ticket 14. No file is written when refused - asserted, because an empty table is a finding. |
| spec edges | **verified** | `gwas` <- `variants` (D6), `phylogeny` <- `recombination`. Asserted one by one, because reachability could not catch their loss: `variants` stays reachable via `reporting`, so a reachability-only check passes on a DAG that has quietly lost a required dependency. |
| stand-in fixtures | **wired** | `test_data/standins/unbuilt_stages/`, one row each, banner + `standin_*` values. `papipeline/standins.py` refuses to read one. Wired, not verified: nothing consumes them, which is the point. |
| `stage_inputs()` in REAL | **verified** | The real function's source is exec'd and called with mode=REAL. The mode gate is never opened; the missing-input naming is proven at TEST level, and a test pins that TEST and REAL take the same branch. |
| dashboard from a STUB log | **not built** | **Ticket 25.** Nothing in `papipeline/` reads `status/events.jsonl`; the observatory serves `/api/events` from an in-process `EventBus`. Two mechanisms that do not meet. |

#### DELETE WITH TICKET 14

`test_data/standins/unbuilt_stages/{variants,recombination,similarity}.tsv` stand in
for the three unbuilt stages. When ticket 14 builds them:

- **delete the three files and the `standins/` directory**, and the
  `VARIANTS_STANDIN` / `RECOMBINATION_STANDIN` / `SIMILARITY_STANDIN` inputs in
  `workflow/Snakefile`;
- **delete `papipeline/standins.py`**, and replace it with whatever the real
  stages need;
- **expect `tests/unit/test_standins.py` to fail**, in
  `TestNoOtherModuleReadsTheseTables`. That failure is the signal that real
  consumers now read those tables and the stand-in assertion has to be replaced,
  not deleted. The other tests in that file go with the fixtures.
- **replace the three contract column sets** in `STAGE_TABLES`, which are marked
  PROVISIONAL placeholders, with the real ones.

#### Two observability mechanisms that do not meet

`scripts/emit.py` writes `status/events.jsonl`, one object per rule transition,
and **nothing in `papipeline/` reads it**. The observatory serves `/api/events`
from an in-process `EventBus` fed by the runner. A STUB run therefore proves the
event *format* and the DAG; it does not exercise the dashboard. `spec.md:62` and
`:562` carry dated status notes saying so, and ticket 25 is the bridge. Not
started.

#### Two defects that only a real execution could find

Both were invisible to `snakemake -n`, which is the argument for the execution test:

1. `scripts/common/run_pipeline.py` read `args.machine` but never defined
   `--machine`, so it raised `AttributeError` from argparse before doing any
   work. The `reporting` rule has never worked.
2. No rule passed `--machine` to `run_stage.py` or `run_pipeline.py`, so the
   overlay chosen on the snakemake command line never reached the stage runner.

#### Doc claims the code does not support

- **The stage list is not the spec's.** The spec's D1 lists *fifteen* stages;
  the code declares *sixteen*, and they are not the same set. The spec's
  `variants` (minimap2/samtools vs PAO1), `recombination` (gubbins) and
  `similarity` have no implementation - `variants` is the subject of ticket 14.
  The code has `structural_variants` in its place, and adds `mechanisms` and
  `regulators`, which the spec's list does not mention. AGENTS.md already warns
  that three taxonomies exist and that the spec wins; this is that conflict,
  unrecorded. The reachability test now derives its stage list from
  `STAGE_ORDER` and calls it `DECLARED_STAGES` rather than "canonical".
- **"The whole DAG and the dashboard are exercised"** was claimed in AGENTS.md
  and is now corrected there. The deeper problem is that the dashboard reads no
  event log at all: `emit.py` writes `events.jsonl`, and the observatory serves
  `/api/events` from an in-process `EventBus`. Two mechanisms that do not meet.
  The same claim survives in three places that were **not** edited, because they
  are either authoritative or describe intended work:
  - `spec.md:62` and `spec.md:562` - the spec states the intent, and the code
    is behind it. Do not "fix" the spec to match the code here.
  - `.opencode/agents/dashboard-ui.md:2,38` - describes the dashboard agent's
    job as "replays and streams status/events.jsonl" and calls the file "one
    source of truth". False of the current code; it would send an agent looking
    for a replay path that does not exist. Worth correcting, flagged rather than
    done because that file is an agent role definition, not a status doc.
- **`run_pipeline.py --mode` help** listed only TEST and REAL; STUB is now
  listed and is the default, per the spec.

**Report the weakest state that is true.** A unit test on a function nothing
calls is evidence the function works, not that the feature exists.

## The 16 pre-existing Bakta failures: RESOLVED

All 16 now pass. They were not caused by this project and they were not caused
by the environment either; they were three separate defects in the tests, all of
which asserted things that were never true.

1. **The test stubs read the genome from `--in`** (5 unit + 4 integration +
   1 command assertion). Bakta 1.12.1 does not accept `--in`, so the adapter
   probes the real binary and passes the genome positionally. Every stub
   therefore got `assembly = None` and exited non-zero, and the stage reported
   `FAILED` where the test expected a classification. **The production code was
   right in every one of these cases.** The stubs now match the genome by
   sequence extension; the adapter was not changed.

2. **The stubs wrote no tab characters at all** (a second, independent defect in
   the same templates). The stub source is itself embedded in a string literal,
   so an escaped tab decayed twice on the way into the generated file and the
   output had no delimiter - which the Bakta reader reports as "no header line".
   The stubs now build rows with `chr(9)`, so no escape layer is involved.

3. **A message-substring mismatch.** The test looked for `"not on PATH"`; the
   pre-flight said `"bakta was not found on PATH"`, which does not contain it.
   The message now reads "bakta is not on PATH and BAKTA_BIN is unset", which
   says the same thing and contains what a reader would search for.

4. **The stale OOM guard - and the constant was not the stale part.** The test
   asserted the figures for a 3.0 GB per-worker estimate. `MEASURED_PEAK_RSS_GB`
   is 1.92, documented with its method and the genome it was measured on. So the
   constant was a *measurement* and the *test* was a stale guess; reversing them
   would have replaced a measurement with a guess. The test now derives its
   boundary from the constant, so a future re-measurement cannot desynchronise
   them again, and a second test pins the boundary from below so the guard
   cannot silently become a blanket refusal.

   That test also revealed something: its `assert check.ok is False` was passing
   for the wrong reason - the fake executable it installs does not exist, so the
   pre-flight was already un-ok. The new pair asserts the *memory guard* fired
   or did not, which is the thing under test.

## Blocker for Phase 3, found in Phase 2

`workflow/Snakefile` **has never loaded.** All 19 rules declare both `script:`
and `shell:`, which Snakemake rejects as a syntax error, so `snakemake -n`
fails before it builds anything. It predates this work (present at `8b35f4e`).
Ticket 12 cannot pass until it is fixed, and the fix is entangled with ticket
11 because the choice of execution keyword determines how `threads` and
`mem_mb` are declared. **No real DAG has therefore ever been resolved here.**

## Phase 2 change of scientific approach

**The association model uses a continuous trait.** The modelled phenotype is
`log2(MIC)` where a measured MIC exists, instead of the binary R/S label. Under a
binary outcome an *intermediate* isolate carrying an MIC would be discarded, and
the only quantitative thing in the row thrown away with it. A row with no
measured MIC is absent from the trait, never zero.

Consequence for the join (ticket 06): `I` with a measured MIC now **enters** the
analysis; `SDD` is also kept when it carries an MIC. `SDD` is *susceptible-dose
dependent*, a distinct **CLSI** category, and it is a determination rather than
a missing result; its EUCAST equivalent is `I` again (Susceptible, increased
exposure). A CLSI SDD call is derived from an MIC inside the susceptible range,
so MIC+SDD is the normal combination. Only `ND` is excluded regardless, because
"not determined" carries no measurement by definition. (An earlier revision of
this phase excluded SDD, on a false premise about the parser; corrected after
review.)

**Breakpoints became optional** and their thresholds remain deliberately empty
pending a licensed CLSI M100 copy. No numbers are written down anywhere. The
source R/S/I call is authoritative and no MIC is reinterpreted. Ticket 08 is
deferred to Phase 3 by instruction.

## Critical path

```
02 ──┐
03 ──┼─ 04 ──┐
│  05 ── 06 ──┼─ 11 ── 12 ── 13 ──┬─ 21 ── 23 ── 24
│  │           │                  │
│  │           └──── 14 ── 15 ── 16 ── 17
│  │
│  07 ── 08 ─────────────── 14
│
09 ──┬─ 10 ──┘
     └─ 18 ── 19 ── 20
```

Tickets with no blockers and no dependents on the critical path can proceed in
parallel: **02, 03, 09, 22** (given 05 and 06).

## Definition of done

1. The environment builds.
2. All pytest tests pass.
3. `snakemake -n` passes with both machine configs.
4. A full STUB-mode run passes.
5. The dashboard shows that STUB run live.
6. The README is verified from a clean folder.
7. The reviewer agent signs off.

**No real-genome run starts until the user explicitly says `run real samples`.**

---

## OPEN — stage 7 GPA contract: long vs wide (next session)

**Status: decided in principle, not implemented. Zero code committed toward
either shape.** The only thing landed is `papipeline/adapters/panaroo.py` (the
preflight and the presence/absence parser), which is independent of this
question and stays.

### Where it stands

The decision is to move `gene_presence_absence.tsv` from the cohort-level
`gene, n_samples, frequency` shape to per-sample **long** format
(`gene, sample_id, present`), because spec D6 needs gene presence/absence as a
pyseer feature family and `n_samples`/`frequency` cannot reconstruct a per-sample
matrix.

pyseer consumes the **wide** `.Rtab` matrix — verified from its own help, not
recalled:

    --pres PRES   Presence/absence .Rtab matrix as produced by roary

panaroo emits `gene_presence_absence.Rtab` natively, so the transpose is
available twice over: keep panaroo's file, or derive it from the long table.
Pick one and say so in the commit. Recommendation: preserve panaroo's `.Rtab`
and assert it agrees with a transpose of the long table, rather than making the
long table the only source of truth.

### Required order — do not reorder this

1. **Write the round-trip test FIRST, and watch it fail.** A test that writes
   via `write_outputs` and reads back via `read_gene_presence_absence`, asserting
   the presence mapping survives. It must fail on today's tree. If it passes
   before any format code changes, it is testing nothing — check that first.
2. **Then change reader and writer in ONE commit.** They must move together.
   Today `write_outputs` writes `matrix_rows()` (wide-ish, cohort counts) while
   `read_gene_presence_absence` parses "first column gene, every remaining column
   is a sample". Changing one without the other silently breaks the round trip.
3. **Regenerate `synthetic.py` and the `test_data/` fixture in that same
   commit.** Leaving either on the old shape is how the inconsistency recurs.

### The trap that cost this session its last hour

A narrow selector gave false confidence:

    pytest -k "pangenome or panaroo or synthetic or gwas"   →  78 passed

That is green **because nothing asserts the contract**, not because the contract
is fine. A half-finished reshape — writer long, reader wide, fixtures old — was
sitting in the tree and that selector reported 78 passes.

**Run the full suite before claiming the reshape works.** A filtered subset is
for iterating; only the full run is evidence. If the full suite hangs, find out
why rather than narrowing the selector to make it finish.

### Also still open, for whoever picks it up

- `TASKS.md` "BLOCKER 3 — panaroo has no osx-arm64 build" is **stale**. panaroo
  works on arm64 via a pinned pip source; `docs/environment-arm64.md` §3 and
  `environment/environment.yml` are correct and were updated on
  `step3-regulators`. This file was not.
- The Snakefile's pangenome docstring said "Stage 9"; corrected to **Stage 7**
  per spec (AGENTS.md: the spec wins). `run()`'s docstring and log message in
  `stages/pangenome.py` still say "Stage 9" — fix them in the same pass as the
  reshape, or as their own commit.
- A REAL-mode test asserting panaroo's identifier shape (`group_<n>`) does not
  exist yet. It is required, because `gene` means a literal gene name in
  TEST/STUB and an ortholog group id in REAL — an accepted TEST-fidelity gap, but
  only if it is written down and asserted.
- Not yet done: the REAL caller wiring into `run()` (preflight → parse →
  `partition` → `write_outputs`), gated on `allow_real_mode` **and** the preflight.
- Not yet done: side-env validation against the proven real numbers —
  **10,019 families / 4,828 core / 5,191 accessory / 0 genome-specific**.
  Do this in a disposable env; `pa-amr` is out of scope.
