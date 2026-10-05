# Spec: *Pseudomonas aeruginosa* imipenem AMR-GWAS pipeline and live monitoring dashboard

Status: ready-for-agent
Created: 2026-09-28
Supersedes: the 16-stage taxonomy in `docs/architecture.md` and `workflow/Snakefile` for stage naming. Where the two disagree, this spec wins; the older documents are not deleted, they are annotated as historical.
Target tracker: `.scratch/imipenem-gwas-dashboard/`

---

## Problem Statement

A researcher has 966 *Pseudomonas aeruginosa* isolates in the PDC Pathogen
Detection and Genomics dataset, of which 835 assemblies are downloaded and
present on disk. They want to know which genomic determinants are associated
with imipenem resistance, and they want to watch the analysis as it runs.

Three things make this hard today.

**The pipeline does not yet call variants.** The existing code ingests variant
calls rather than producing them; the architecture document states plainly that
no variant caller is wired in. Variant calling exists only as an ad-hoc pilot
script. A GWAS cannot be run on ingested tables alone: the headline mechanism
for carbapenem resistance in this organism is *loss of the OprD porin*, which
is a gene-level and structural observation, and the resistance phenotype is
strongly confounded by clonality. Without SNP/indel calls against a pinned
reference, without gene presence/absence, and without a phylogeny-aware model,
any result is uninterpretable and the multiple-testing correction is wrong.

**The pipeline cannot be run at all on this machine.** The environment has never
been built. Several pinned packages are believed to be misnamed or
mis-versioned, and some of the tools required by the target design have no
`osx-arm64` build at all. The configuration system hard-codes the development
machine's paths, thread counts and memory limits into a single file, so the same
code cannot be pointed at a large Linux machine without editing source.

**Long runs are invisible.** Annotation across hundreds of genomes takes hours.
There is no way to see, from another machine through an SSH tunnel, which
isolates are done, which stage is grinding, what has failed, or how long the run
has been going. A run that dies at isolate 700 of 900 leaves the user with log
files and no state.

There is also a scientific constraint that shapes everything below. Public
genome archives are not a random sample of the species. They over-represent
clinical and outbreak isolates, and they over-represent a small number of
high-risk clones and of the *bla*OXA-*alleles* that travel with them. A naive
association scan on such a cohort will confidently rediscover clonal structure
and report it as a resistance mechanism. This is precisely why the target
design includes a recombination filter, a phylogeny-based mixed linear model,
and an explicit convergence-across-independent-lineages stage: those three are
not optional extras, they are the controls that make a hit interpretable.

## Solution

A Snakemake pipeline of sixteen stages, one rule per stage, each declaring its
own threads and memory, driven entirely by configuration and runnable in three
modes. A live dashboard that replays an append-only event log on start and then
streams new lines, so the state of a long run is visible from anywhere.

The pipeline runs in three modes:

- **STUB** — every rule emits a tiny, plausible-shaped fake output. The whole
    DAG, the event stream and the dashboard are exercised in seconds with no
    bioinformatics tool installed and no fixture on disk. This is the default mode
    for development and CI.
    > **Status (2026-09-29): intent, partly not implemented.** The DAG and the
    > event stream are exercised and asserted, by a real `snakemake` run. **The
    > dashboard is not fed this stream**: `scripts/emit.py` writes
    > `status/events.jsonl` and nothing in `papipeline/` reads it — the
    > observatory serves `/api/events` from an in-process `EventBus`. See
    > TASKS.md.

- **TEST** — the real code paths run over 20 committed, byte-stable synthetic
  fixtures. Real parsers, synthetic input.
- **REAL** — actual genomes, gated on both a configuration flag and the user
  saying `run real samples`.

Three specific outcomes the solution must deliver:

1. **The GWAS input carries all three required feature families** — SNP/indel
   variants from the PAO1 alignment, gene presence/absence, and a gene-level
   OprD loss-of-function feature. Gene presence/absence alone is explicitly not
   sufficient, because OprD is frequently present-but-broken, which is a
   different observation with a different mechanism.
2. **The multiple-testing threshold is derived from unique variant patterns**,
   in a separate, inspectable step, so the number of independent tests is a
   function of the data rather than of the size of the variant table.
3. **The dashboard shows a live run**, survives a restart by rebuilding from
   the log file, and stays responsive at 900 isolates × 15 stages.

---

## User Stories

### Configuration and portability

1. As a researcher on a laptop, I want the pipeline to read all paths, thread
   counts, memory limits, database locations and tool availability from a single
   machine config, so that I never edit source to run it somewhere else.
2. As a researcher moving to a large Linux machine, I want to switch machines by
   changing which config file is passed, so that no code changes are needed.
3. As a researcher, I want the laptop config to refuse to run more than 20
   samples, so that I cannot accidentally start an overnight job on a machine
   with 16 GB of RAM.
4. As a researcher hitting that limit, I want an error message that tells me the
   limit, the file that caused it, and what to change, so that I can fix it
   without reading the source.
5. As a reviewer, I want every scientific constant — organism, genome size
   expectations, breakpoints, gene lists, tool and database versions — to live
   in one machine-independent place, so that the laptop result and the
   big-machine result cannot silently disagree about the science.
6. As a reviewer, I want to grep for a hard-coded path or an `os.cpu_count()`
   and find nothing, so that "nothing is hard-coded" is verifiable rather than
   asserted.
7. As a researcher, I want to know which tools are unavailable on this machine,
   so that I understand why a stage is gated before I hit it.
8. As a researcher on the laptop, I want a stage that needs an unavailable tool
   to fail at DAG-build time with a message naming the tool, rather than
   crashing mid-run, so that I find out before waiting hours.
9. As a researcher, I want STUB mode to run the entire DAG with no tools
   installed, so that I can develop the dashboard and the orchestration on a
   machine where nothing is provisioned.

### Data integrity and identity

10. As a researcher, I want `sample_id` to be the versioned GCA assembly
    accession, so that it is already the genome directory name and the FASTA
    filename stem and the assembly join needs no lookup table.
11. As a reviewer, I want the PDC isolate ID, BioSample and isolate name carried
    as validated attributes alongside the assembly accession, so that the
    provenance chain from assembly back to isolate is inspectable.
12. As a researcher, I want a mismatch between the manifest and the phenotype
    table to fail loudly, so that a silent drop never shrinks my cohort without
    my noticing.
13. As a researcher, I want every genome in the manifest to have exactly one
    phenotype row, so that a missing or duplicated row is an error I see rather
    than a sample that quietly disappears.
14. As a researcher, I want phenotype rows with no corresponding assembly to be
    excluded and reported, not failed, so that a phenotype table describing more
    isolates than I downloaded is a normal situation and not a crash.
15. As a researcher, I want isolates whose phenotype is `I`, `SDD` or `ND`, or
    which have no phenotype row at all, to be excluded from the binary analysis
    and written to an exclusion table with a reason, so that I know exactly who
    dropped out and why.
16. As a researcher, I want the count of excluded isolates reported, so that the
    analysed cohort size is never a mystery.
17. As a researcher, I want `I` to be excluded from the binary outcome rather
    than folded into `R` or `S`, so that intermediate susceptibility is not
    silently converted into a claim.
18. As a researcher, I want a sample with no phenotype record to stay absent and
    never default to susceptible, so that absence of evidence is not evidence of
    susceptibility.
19. As a reviewer, I want any sample-ID mismatch discovered *inside* the
    pipeline — between assemblies, annotations, variant tables and the final
    phenotype table — to be a hard failure, so that no fuzzy join can ever
    mis-attribute a variant to a sample.
20. As a researcher, I want the sample-ID format rule to be configuration-driven,
    so that both the synthetic fixture IDs and real GCA accessions validate
    while a malformed ID is still rejected.

### Reference and provenance

21. As a researcher, I want a pinned PAO1 reference used for all variant
    calling, so that variant coordinates are interpretable and stable.
22. As a researcher, I want the reference's accession and checksum recorded in
    the references table, so that I can prove which reference produced a result.
23. As a researcher, I want the reference's actual length and checksum verified
    against the real file on disk at run time, so that a truncated or
    substituted download is caught rather than trusted from a summary.
24. As a researcher, I want nothing downloaded during a run, so that a run is
    reproducible offline and cannot be perturbed by a remote server.
25. As a researcher, I want database provisioning to be a separate deliberate
    step, so that a run never silently upgrades a database mid-analysis.

### Phenotype interpretation

26. As a researcher, I want the source R/S/I call to be authoritative, so that
    the analysis reflects the reported laboratory interpretation.
27. As a researcher, I want MIC-to-R/S re-interpretation to run only on rows
    carrying a genuinely measured MIC, so that no categorical call is ever
    manufactured from a number that does not exist.
28. As a researcher, I want the breakpoint standard, its edition and its table
    number recorded in configuration, so that a susceptibility call can be
    traced to a specific published table.
29. As a researcher, I want a measured MIC that the active standard does not
    cover to be a hard error naming the row, so that an unstandardised value is
    never guessed into a category.
30. As a researcher, I want the pipeline to refuse to interpret MICs at all
    while the breakpoint values are unverified, so that no unsourced threshold
    can reach a result.
31. As a researcher, I want each phenotype row to carry the AST method, the
    standard and the edition where the source provides them, so that a
    susceptibility call can be assessed for comparability.
32. As a researcher, I want the number of rows lacking that provenance reported,
    so that I know how much of my cohort is methodologically comparable.
33. As a researcher, I want a planned sensitivity analysis comparing source
    calls against re-interpreted calls, so that I can quantify how much my
    conclusions depend on the source laboratory's interpretation.

### Variant calling and the reference loci

34. As a researcher, I want variants called against PAO1 for every isolate, so
    that I have a SNP/indel feature family for the GWAS.
35. As a researcher, I want the call restricted to a pinned set of reference
    loci including `oprD` and the regulators `mexR`, `nalC`, `nalD`, `mexZ`,
    `ampD` and `ampR`, so that the analysis targets the known mechanisms
    without pretending to have surveyed everything.
36. As a researcher, I want OprD loss-of-function to be assessed, so that the
    principal carbapenem entry porin is examined for the mechanism that actually
    drives most imipenem resistance.
37. As a researcher, I want "OprD absent" and "OprD present but disrupted" to be
    two separate features, so that the model can distinguish a porin that is
    gone from a porin that is present and broken, which are mechanistically
    different.
38. As a researcher, I want OprD loss-of-function to be a gene-level binary
    feature in the GWAS input, so that the strongest known mechanism is
    directly testable rather than only indirectly represented.

### Phylogeny and confounding

39. As a researcher, I want a recombination filter applied before tree building,
    so that horizontally acquired resistance determinants do not distort the
    tree and inflate support for false convergence.
40. As a researcher, I want a maximum-likelihood tree built from the
    recombination-filtered alignment, so that the phylogeny used for the
    downstream model is the corrected one.
41. As a researcher, I want a phylogenetic distance matrix computed from that
    same tree, so that the kinship matrix and the tree cannot disagree.
42. As a researcher, I want the tree's tip set to equal the analysis manifest
    exactly, so that a tree containing a dropped isolate is caught.
43. As a researcher, I want the convergence stage to require a determinant to
    appear in at least two independent lineages, so that a determinant riding a
    single successful clone is not reported as convergent.
44. As a researcher, I want widespread determinants distinguished from rare
    ones, so that background carriage is not reported as a finding.

### GWAS

45. As a researcher, I want pyseer run with a mixed linear model, so that
    relatedness is accounted for rather than assumed away.
46. As a researcher, I want the mixed model's phylogeny to be the corrected
    maximum-likelihood tree, so that kinship is modelled on recombination-filtered
    relatedness.
47. As a researcher, I want the multiple-testing threshold computed from unique
    variant patterns, so that perfectly correlated variants do not inflate the
    number of tests.
48. As a researcher, I want that threshold step to be a separate, inspectable
    step whose output I can read, so that I can see how many patterns were
    dropped and why.
49. As a researcher, I want the GWAS input to contain all three feature families
    together, so that I can ask which family carries the signal.
50. As a researcher, I want each group's minimum sample size enforced, so that a
    "significant" hit supported by two isolates cannot be reported.
51. As a researcher, I want a feature confounded with a single lineage reported
    as lineage-linked, so that clonal artefact is visible in the output rather
    than buried in a p-value.
52. As a researcher, I want unitigs available as an optional input, so that I
    can add a k-mer-based family later, but I do not want them counted among the
    families the pipeline currently claims to support.

### Combination and reporting

53. As a researcher, I want co-occurrence tested across the determinant classes,
    so that I can see which mechanisms travel together.
54. As a researcher, I want co-occurrence corrected for multiple testing and to
    respect a minimum cell count, so that a two-by-two cell of one isolate is
    not reported as an association.
55. As a researcher, I want a final report that compares the hits against the
    known imipenem resistance mechanisms, so that the output is interpretable
    biology rather than a table of p-values.
56. As a researcher, I want no stage to state or imply that detecting a
    determinant establishes phenotypic resistance, so that the report does not
    overstate what association analysis can show.

### The dashboard

57. As a researcher, I want a live dashboard, so that I can watch a long run from
    another machine instead of waiting.
58. As me, I want the dashboard to bind to the loopback interface only, so that a
    port on a shared analysis machine is not exposed to the network.
59. As me, I want to reach the dashboard over an SSH tunnel, so that I can view
    a run on a remote Linux machine without opening a firewall port.
60. As me, I want metric cards for annotated, running, failed and elapsed, so
    that the four numbers I actually check are above the fold.
61. As me, I want a card per stage showing waiting, running, done or failed, so
    that I can see which stage the run is in.
62. As me, I want one square per isolate per stage, coloured grey when waiting,
    blue when running, green when done and red when failed, so that I can spot
    the failures and the stragglers at a glance.
63. As me, I want a live event log, so that I can see what just happened rather
    than inferring it from the grid.
64. As me, I want the dashboard to rebuild its state from the log file when it
    restarts, so that a browser refresh or a server restart loses nothing.
65. As me, I want the dashboard to stay responsive at 900 isolates, so that it
    is usable at the scale the final run will actually reach.
66. As me, I want a single event written when each rule starts, succeeds and
    fails, so that the grid can be driven from events alone.
67. As me, I want aggregate rules to record one event for the whole cohort, so
    that the event volume stays proportionate to the work.
68. As me, I want the dashboard to work in dark mode, so that it is readable in
    the conditions I actually work in.
69. As me, I want the dashboard to show an explicit "not reported" where the
    event stream carries no fact, so that a missing measurement is never
    displayed as a plausible zero.
70. As a developer, I want a malformed or truncated line in the event log to be
    skipped with a warning rather than crashing the stream, so that a run killed
    mid-write still leaves a usable dashboard.
71. As a developer, I want one event to update one square rather than rebuilding
    the whole grid, so that the UI does not degrade as the run proceeds.
72. As a researcher, I want to be able to test the whole pipeline and the whole
    dashboard with no real data and no real tools, so that development is fast
    and carries no risk to the real dataset.

### Environment and reproducibility

73. As a researcher, I want the environment to build on this machine, so that I
    can actually run anything.
74. As a researcher, I want to be told which required tools have no build for
    this architecture, so that I know what cannot run here before I try.
75. As a researcher, I want a locked environment file for the Linux machine, so
    that the final run is reproducible.
76. As a reviewer, I want no database or genome file in version control, so that
    the repository stays source-only.
77. As a reviewer, I want to know exactly which tool versions and database
    versions a result was produced with, so that a result is reproducible or
    honestly not.
78. As a researcher, I want the README's instructions to be verified from a
    clean checkout, so that I can trust them.

---

## Implementation Decisions

### D1. Canonical stage list — sixteen stages

These are the stages. One Snakemake rule each, each declaring `threads` and
`mem_mb`.

| # | Stage | Tool | Kind |
|---|---|---|---|
| 1 | `validate` | seqkit | per-sample |
| 2 | `annotate` | bakta | per-sample |
| 3 | `mlst` | mlst (scheme `paeruginosa`) | per-sample |
| 4 | `amr` | AMRFinderPlus (`--organism Pseudomonas_aeruginosa`) | per-sample |
| 5 | `virulence` | BLAST vs VFDB | per-sample |
| 6 | `variants` | minimap2 + samtools vs PAO1 | per-sample |
| 6a | `cohort_variants` | internal | cohort |
| 7 | `pangenome` | panaroo | cohort |
| 8 | `recombination` | gubbins (needs mafft) | cohort |
| 9 | `phylogeny` | iqtree | cohort |
| 10 | `similarity` | phylogenetic distance from the stage-9 tree | cohort |
| 11 | `phenotype` | internal | cohort |
| 12 | `gwas` | pyseer, with 12a = unique-pattern filter, 12b = association | cohort |
| 13 | `convergence` | internal | cohort |
| 14 | `combination` | internal | cohort |
| 15 | `report` | internal | cohort |

The existing sixteen-stage code is **mapped onto** these, not discarded and not
allowed to define its own taxonomy. The existing `mechanisms`, `regulators`,
`structural_variants` and `integration` modules become internal steps within the
rule that owns them; their tests and scientific content are preserved. Existing
`STAGE_TABLES` and stage ordering constants are rewritten to the fifteen.

**Amended 2026-09-29: fifteen to sixteen.** `cohort_variants` was added at 6a
after stage 6 was actually run against a real assembly. Per-isolate calls against
PAO1 conflate species-wide fixed differences with cohort-informative
polymorphisms — a position every isolate differs at is not polymorphic — and only
a multi-isolate comparison separates them. It is a separate stage rather than a
second phase of `variants` because the per-sample/cohort distinction is
load-bearing in this DAG, and one rule with two input cardinalities would
undercut it. Precedent: snippy / snippy-core make the same split.

Its output contract and filtering threshold are deliberately **not** fixed here,
and that is a deliberate gap rather than an oversight — see
`docs/design/cohort-variant-merge.md`. The stage exists in the taxonomy and
refuses outside STUB until it is built, so the DAG is honest about the missing
step instead of silently running per-isolate calls as though they were the
answer.

### D2. The phylogeny chain

`gubbins` produces the recombination-filtered alignment **and** an intermediate
tree. `iqtree` builds the **authoritative** maximum-likelihood tree from the
filtered alignment. Stage 10's distance matrix and stage 12's mixed linear model
both consume the stage-9 tree. Gubbins' own tree is an intermediate and is not
used for inference. Rationale: modelling relatedness on a tree built without
recombination correction would reintroduce exactly the confounding stages 8 and
12 exist to control.

### D3. Sample identity and the directional join rule

`sample_id` is the versioned GCA assembly accession.

**Directional, and asymmetric by design:**

- Every genome in the manifest must have **exactly one** phenotype row. Zero
  rows, or two rows for the same sample, is a **hard failure** naming the sample.
- A phenotype row with **no** corresponding assembly is **excluded**, not failed.
  The phenotype table legitimately describes more isolates than were downloaded.
- An assembly whose phenotype is unusable (`I`, `SDD`, `ND`, or absent) is
  **excluded**, not failed, and is written to an exclusion table with a reason.
- Exclusion counts are logged and reported at the end of the run.

This replaces a symmetric 1:1 assertion. The asymmetry is the point: the
manifest is the authority on what exists, the phenotype table is external
evidence about it, and an external file may legitimately speak about things you
do not have.

**Inside** the pipeline, any sample-ID mismatch across assemblies, annotations,
variant tables and the final phenotype table is a **hard failure**. The
directional rule governs the manifest-to-phenotype edge only. There is no fuzzy
matching, no prefix matching, no silent drop.

The ID *format* rule is configuration-driven, so both `TEST_PA_001` and
`GCA_000000000.1` validate while a malformed ID is still rejected. Contig
headers are **not** sample IDs — the synthetic fixtures use
`TEST_PA_001_contig1` and real assemblies use NCBI-style identifiers — and no
check ties a contig header to `sample_id`.

### D4. Phenotype interpretation, and the continuous trait (D4a added Phase 2)

Source categorical R/S/I is authoritative. MIC-to-R/S re-interpretation runs
**only** on rows carrying a genuinely measured MIC, against a configuration-pinned
breakpoint standard. `I`, `SDD` and `ND` are excluded from the binary outcome and
never folded into `R` or `S`. An absent record stays absent.

**The association model is to use a continuous trait, not the binary label.**
The modelled phenotype is to be `log_base(MIC)` wherever a measured MIC exists,
rather than the R/S category. This is a change from the earlier binary-outcome
decision, made deliberately:

* it uses the measurement the row actually carries, rather than a 2-level
  summary of it;
* it keeps rows that a binary outcome would discard. An *intermediate*
  isolate with a measured MIC is a real determination with a quantity attached,
  and pooling it into R or S destroys information that was collected;
* pyseer's continuous mode expects a numeric trait.

Until ticket 16 lands, `gwas.outcome` in the science configuration remains the
binary R-vs-S mapping, and that is what pyseer receives. `gwas.outcome` is
retained rather than deleted so the change is a switch, not a rewrite.

A row with **no** measured MIC is **absent** from the trait, not zero. Zero
would be a fabricated measurement sitting inside the trait's range, pulling the
model towards it.

**How `I` and `SDD` are handled in the two traits, side by side.** This is
recorded explicitly because the two coexist until ticket 16, and a silent
disagreement between them would be invisible in the output.

| Row | Continuous trait (built, not yet modelled) | Binary R-vs-S model (runs today) |
| --- | --- | --- |
| `R` with an MIC | kept, `log2(MIC)` | positive (R) |
| `S` with an MIC | kept, `log2(MIC)` | negative (S) |
| `I` with an MIC | **kept**, `log2(MIC)` | **excluded** |
| `I` without an MIC | excluded, reason `unusable_category` | excluded |
| `SDD` with an MIC | **kept**, `log2(MIC)` | **excluded** |
| `SDD` without an MIC | excluded, reason `unusable_category` | excluded |
| `ND`, `I`, `SDD` with no MIC | excluded | excluded |
| `R` or `S` with no MIC | **kept** in the cohort join (`require_trait=False`), but contributes no value to the continuous matrix — absent, never zero | positive / negative |

So today the binary model analyses **strictly fewer** isolates than the
continuous trait would: every `I` and `SDD` isolate carrying a measured MIC is
in the cohort and is excluded from the association. That is deliberate for now
— a binary model has no way to represent "more susceptible than S but less than
R" — and it is recorded here so the two cannot drift apart unnoticed.

`gwas.outcome` in the science configuration lists `I`, `SDD` and `ND` as
`excluded`, and that is the authority for the binary model. When ticket 16
switches the model, this table is what changes, and the join's behaviour does
not.

> **Status: BUILT, NOT WIRED.** Phase 2 built the continuous trait — it is
> computed, carried on every record, and is what the directional join uses to
> decide inclusion. **The model does not consume it yet.** The GWAS stage still
> builds a binary R-vs-S outcome from `gwas.positive` / `gwas.negative`, and
> `phenotype_matrix` has no production caller. Switching the model is ticket
> 16. Do not read this decision as though the association run already uses
> log2(MIC); it does not.

The categorical view is not lost: it is retained separately for reporting, so
the two views are never conflated.

**Breakpoints are optional, and currently absent.** CLSI M100 is a paywalled
standard and no free primary source publishes the *P. aeruginosa* imipenem S/I/R
triple. The breakpoint configuration therefore ships with the standard, the
edition, the organism, the agent and the log base, and an **empty** `thresholds`
map. No numbers are written down. Consequences:

* the source R/S/I call is authoritative and is never overwritten;
* no MIC is reinterpreted into a category, and the run says so explicitly,
  naming the standard it *intends* to use;
* a test asserts that `thresholds` stays empty, so an unsourced number cannot
  drift in unnoticed.

Supplying the S/I/R values from a licensed copy switches re-interpretation on;
nothing else changes. The continuous trait does not depend on breakpoints and
works without them.

**Provenance per row.** Where the source provides them, each phenotype row
carries the AST method, the standard, and the edition. The run reports how many
rows lack each, and a missing *column* is counted as missing on every row, which
is the honest reading. This directly feeds the comparability question — an
MIC-derived category from one laboratory's method is not automatically
comparable with another's.

**A sensitivity run is planned, not built.** Comparing source calls against
re-interpreted calls, to quantify how much conclusions depend on the reporting
laboratory's interpretation, is deferred until thresholds exist. It is recorded
here so it is not lost.

### D5. Reference genome

PAO1, assembly `GCF_000006765.1` (ASM676v1), RefSeq `NC_002516.2`, 6,264,404 bp,
BioProject PRJNA57945, taxid 208964. Verified against NCBI assembly records.

The accession and checksum are pinned in the references table. The file path
comes from machine configuration. Provisioning is **out of band** — the
pipeline never downloads during a run. At run time the pipeline verifies the
file's actual length and checksum against the pinned values and **fails on
mismatch**, rather than trusting a summary or an earlier note. A truncated or
substituted reference silently shifting every variant coordinate is exactly the
failure that is invisible until the biology is wrong.

### D6. GWAS input construction

Three feature families, all required:

1. **SNP/indel** variants from the PAO1 alignment, in the format pyseer
   consumes.
2. **Gene presence/absence** from the pangenome.
3. **OprD features**, as two separate binary features:
   - `oprD_absent` — the locus is wholly absent.
   - `oprD_LoF` — the locus is present but carries a truncating variant:
     frameshift, premature stop, or an internal deletion.

   These are mechanistically distinct and collapsing them destroys the
   distinction between a porin that is gone and a porin that is present and
   broken. Presence of the locus is **not** evidence of susceptibility.

**Unitigs are deferred.** They are recorded as an optional future input, not as
one of the three families this pipeline supports. Claiming a k-mer family the
pipeline does not run would be a false statement in the report.

**Unique variant patterns.** A separate step (12a) reduces the variant matrix to
unique patterns, writing the survivors to disk; 12b runs the association on that
set. Splitting them makes the multiple-testing threshold auditable: the user can
read how many patterns were dropped and inspect the file. pyseer's helper script
name and flags will be confirmed against the installed package's `--help` before
use — not from memory.

`--lmm` is used, with the stage-9 tree. Minimum samples per group is enforced. A
feature confounded with a single lineage is reported as lineage-linked rather
than as an association.

### D7. Configuration layout

- A shared, machine-independent `science.yaml`: organism, genome size
  expectations, breakpoints, gene and regulator tables, tool versions, database
  versions.
- A `laptop.yaml` overlay: paths, threads, memory, `max_samples`, per-tool
  availability.
- A `bigmachine.yaml` overlay: the same keys, no cap.

**No duplication.** The scientific constants exist in exactly one place, so the
laptop result and the big-machine result cannot diverge on science. This is the
main reason for the split; two standalone files would drift.

The laptop config carries `max_samples: 20`. Exceeding it stops the run with a
message stating the limit, naming the config key and file, and saying what to
change.

**Tool availability.** Each overlay declares per-tool availability. A rule
requiring an unavailable tool fails at **DAG-build** time in REAL mode with a
message naming the tool and the machine class it needs, so the failure is
immediate. STUB and TEST still exercise the entire DAG regardless of tool
availability, because they invoke no tools.

### D8. Modes

| Mode | Fixtures | Real tools | Real parsers | Purpose |
|---|---|---|---|---|
| STUB | none | no | no | DAG, events, dashboard; default for dev/CI |
| TEST | 20 committed | no | yes | parser coverage, regression safety |
| REAL | manifest | yes | yes | the analysis |

STUB exists so the DAG, the event stream and the dashboard can be exercised with
no bioinformatics tool installed and no fixture on disk. It is the default mode
for development and CI.
> **Status (2026-09-29): intent, not implemented.** The DAG and the event stream
> are exercised and asserted, by a real `snakemake` run that produces every
> stage's declared output and a well-formed event log. **The dashboard is not
> exercised**: no test feeds a STUB event log to it, and none could, because
> nothing reads `status/events.jsonl` — the observatory uses an in-process
> `EventBus`. See TASKS.md for the two mechanisms and the bridging ticket.

REAL mode is gated twice: a configuration flag, and the user's explicit
`run real samples`. Both must hold.

### D9. Events

Every rule appends one JSON object per line to `status/events.jsonl` through a
single emitter script. Fields:

| Field | Type | Meaning |
|---|---|---|
| `t` | ISO-8601 string | timestamp |
| `event` | `start` \| `done` \| `fail` | transition |
| `stage` | string | stage name |
| `sample` | string | sample ID, or the literal `all` for aggregate rules |

`event` has exactly three values. `sample` is `all` for aggregate rules, which
keeps event volume proportionate to the work done rather than to the cohort
size. The native orchestrator emits the same events as the Snakemake rules, so
the dashboard's source of truth is the same file regardless of entry point.

### D10. Dashboard architecture

A new, separate application. It does not import from, modify, or share state
with the existing observatory. See "Further Notes" for the evidence behind that
separation.

- **Backend**: FastAPI. On start, replays the whole event log to rebuild state.
  Then tails the file and streams new lines over SSE.
- **Binding**: `127.0.0.1` only. Reached over an SSH tunnel. An externally bound
  port on a shared analysis machine is a data leak.
- **Resilience**: a line that is not valid JSON, or is missing a field, is
  skipped with a logged warning. A run killed mid-write leaves a replayable log.
- **Frontend**: one page. Metric cards (annotated, running, failed, elapsed);
  stage cards (waiting/running/done/failed); one square per isolate per stage
  (grey/blue/green/red); a live event log. Colours are CSS custom properties
  defined once, so dark mode works without duplicating a palette. Status is not
  conveyed by hue alone.
- **Honesty**: where the event stream carries no fact, the UI shows an explicit
  "not reported". Never a plausible zero.
- **Scale**: 900 isolates × 16 stages is up to ~14,400 cells. One event updates
  one cell; the grid is not rebuilt per event; DOM writes are batched.
  Performance is **measured** against a realistically-sized synthetic log, not
  asserted.

### D11. Environment

`environment/environment.yml` contains only packages that actually resolve on
`osx-arm64`, so the environment always builds on this machine. Tools with no
build for this architecture are recorded as unavailable in `laptop.yaml` and
their rules are gated there — not removed from the specification, and not
silently substituted with a different tool that does a different job. A separate
locked file is produced for the Linux machine.

Suspect pins to be corrected rather than trusted: the AMRFinderPlus package name,
the BLAST version line, and the iqtree package. Gubbins and mafft are added only
if they resolve. Every version is confirmed by search and by `--help` before it
is written down.

### D12. Selection bias is a design input

Public archives over-represent clinical and outbreak isolates and a small number
of high-risk clones. Three of the sixteen stages exist specifically to control
for this, and none of them is optional:

- **Stage 8** (recombination filter) — stops horizontally transferred
  determinants from distorting the tree and manufacturing convergence.
- **Stage 12's mixed linear model** — models relatedness instead of assuming
  independence, so a clonal artefact is absorbed rather than reported.
- **Stage 13** (convergence across independent lineages) — a determinant must
  appear in at least two independent lineages before it is called convergent.

Without all three, a scan of this cohort will confidently report clonal structure
as a resistance mechanism. This is stated in the report so a reader knows the
controls that were applied.

### D13. Ingestion is not inference

No stage states or implies that detecting a determinant establishes phenotypic
resistance. A variant in a resistance-associated locus is reported as a variant
plus its mechanism class, with its evidence ceiling. Association is quantified
in stage 12 and its robustness assessed in stage 13. Promoting a candidate
structural call to a confirmed one is disabled by configuration.

---

## Testing Decisions

### What makes a good test here

A good test asserts **externally observable behaviour** through a stable
interface: a returned record, a written file, an emitted event, an error
message. A bad test re-implements the logic and compares it to itself, or
asserts an intermediate representation that exists only as a refactoring
artefact.

The specific failure this guards against: a validator that "passes" because it
reproduces the same bug twice. Every validator therefore needs a case that the
*current* implementation gets wrong, or a case that would be wrong if the
implementation were subtly loosened.

**A failing test is never edited to pass.** Not by loosening an assertion, not
by widening a tolerance, not by adding `xfail` or a skip, not by changing a
fixture. The failure is diagnosed. Tests are written before the implementation.

### Test seams, highest first

The preference is few seams, as high as possible.

**Seam 1 — the whole-DAG stub run.** The highest seam. One STUB-mode run must
produce every stage's declared output, a well-formed event log covering all
three event values, and no hard-coded path. This single seam exercises the DAG
shape, the event protocol and the config plumbing at once. Most orchestration
regressions are caught here.

**Seam 2 — the mode-level pipeline run.** Run TEST mode over the 20 committed
fixtures and assert the final table and report. This catches contract and
scientific-logic regressions across all stages without a real tool.

**Seam 3 — per-component seams**, test-first, for the parts that need
diagnosing in isolation:

- **Config loader** — machine overlay resolution; `max_samples` enforcement
  including the message content; unknown-key rejection; per-tool availability.
- **Manifest and sample identity** — GCA-format and fixture-format IDs both
  valid; malformed IDs rejected; duplicate IDs rejected; missing assembly files
  rejected; PDC attribute mapping 1:1.
- **Directional join** — the asymmetric cases specifically: one phenotype row per
  manifest genome (zero and two both fail); extra phenotype rows excluded, not
  failed; `I`/`SDD`/`ND`/absent excluded with reasons; counts reported.
- **Phenotype parser** — all five categories; lower-case normalisation; unknown
  category rejected; non-numeric MIC rejected; MIC with `ND` rejected; measured
  MIC preserved; unverified breakpoints hard-fail naming the field; AST method /
  standard / edition carried and counted when missing.
- **Variant and OprD feature construction** — LoF classification for frameshift,
  premature stop and internal deletion; `oprD_absent` versus `oprD_LoF` kept
  distinct; near-boundary indels not silently called disruptive.
- **Unique-pattern reduction** — that a matrix of perfectly correlated variants
  collapses to one pattern, and the count of survivors is reported.
- **Event emitter** — exact field set; the three event values; `all` for
  aggregate rules; append-only behaviour; that a malformed line is skipped and
  the stream survives it.
- **Dashboard replay** — that state is rebuilt correctly from a written event
  log with no live stream; that a truncated final line does not break replay;
  that restart-and-replay matches continuous streaming.
- **Environment and gitignore** — that no genome, database or secret is tracked;
  that the environment file resolves on this architecture.

### Prior art in this repository

The existing suite is the prior art for style: pytest, plain functions, no
fixtures framework beyond what exists, real temporary directories for
filesystem-touching tests, and explicit assertions on error message content
rather than on exception type alone. The 20 committed synthetic fixtures are the
prior art for deterministic input, and their byte-stability is load-bearing for
every assertion that depends on them.

### Fixtures

The 20 committed synthetic fixtures are **kept as they are**, including their
`TEST_PA_*` identifiers. Renaming them would break the byte-stability guarantee
that existing assertions rest on, for no analytical gain, given that the ID
format rule is configuration-driven. Real data is never read by the test suite.

---

## Out of Scope

- **Read QC and assembly.** The pipeline consumes assembled genomes. It never
  assembles reads and never cleans them. This is stated in the README.
- **Unitigs / k-mer features.** Deferred, not part of the three current families.
- **The MIC-reinterpretation sensitivity analysis.** Planned and recorded in D4,
  not built.
- **Any real-genome run.** Gated on the user's explicit instruction. This spec
  is delivered and verified entirely on stub and synthetic data.
- **Downloading the remaining assemblies.** The 835 on disk are the available
  source data. The eventual ~899–900 isolate analysis is a later, larger-machine
  effort.
- **Modifying the existing observatory.** It is a separate application with a
  separate store. It is not merged, replaced, or imported from.
- **Multi-antibiotic analysis.** Imipenem only. `config/antibiotics.tsv` already
  carries meropenem as a configured-but-disabled example; enabling it is an
  append, not a code change.
- **Reading, modifying or committing anything under `data/`.**
- **Promoter-region assessment.** Requires a reliable upstream window and pinned
  coordinates; disabled in configuration until both exist.
- **Deleting the legacy 16-stage documentation.** The older documents are
  annotated as historical, not removed.

---

## Further Notes

### Dashboard overlap with the observatory — evidence

The existing observatory was read in full before this design was written. The
decision to build a separate application is not a preference for novelty; it is
forced by a specific structural fact.

**The decisive fact.** The observatory's store records **one cohort-level row
per stage**. The pipeline's stage runner sets the subject to the empty string,
so a 900-isolate run produces on the order of sixteen rows — one per stage — not
13,500. The observatory's own snapshot model is built on that assumption: its
expected-total calculation counts distinct non-empty subject values observed
anywhere, and its table contracts are per-stage tables.

**Therefore the store cannot produce the required grid.** The dashboard needs one
cell per isolate per stage. Reading the observatory's SQLite yields *better state
per task*, not *more tasks*. The nearest thing to a subject-by-state map in the
existing frontend is a three-line accessor that builds exactly the needed key
and is then never read — a two-dimensional grid with a dead accessor, not a
rendering.

**Second fact: the store is empty unless the observatory was enabled.** It is off
by default, in four independent places. When off, each stage takes a different
code path that calls the stage directly and writes nothing. A dashboard pointed
at the observatory database after a default run would show an empty, idle
observatory.

**Third fact: incompatible event vocabularies, with no lossless mapping.** The
observatory emits at least sixteen distinct event names carrying run key, state,
attempt, sequence, detail and payload. This specification has three event values
and four fields. The observatory deliberately distinguishes "completed" from
"validated" as two separate facts; the new protocol collapses both into `done`.
The new protocol's aggregate convention — the literal sample `all` — has no
observatory equivalent. Adapting either direction is a lossy translation layer
that must be maintained forever.

**Fourth fact: the in-memory event bus cannot satisfy replay-after-restart.** The
bus holds a bounded buffer inside the pipeline's own process. If the server
restarts, or the pipeline exits, that history is gone. Durability in the existing
design comes entirely from SQLite, and the stream is only ever a notification
that the store changed. There is no durable file-shaped event history. The new
dashboard's central requirement — replay the log on start, then stream — has no
existing implementation to inherit.

**Fifth fact: a different rendering problem.** The observatory draws sixteen
circles and roughly thirty-two edges on a canvas; its own code comments say the
bottleneck would be layout maths rather than fill rate. A 900 × 15 matrix with
per-cell incremental updates is a different problem requiring different
techniques.

**Sixth fact: history semantics differ.** The store upserts on
`(run_key, stage, subject)`, so re-running overwrites the previous row for the
same triple while accumulating attempt history. An append-only log cannot do
that. The two media have genuinely different history semantics, and the
append-only property is precisely what makes the dashboard crash-tolerant.

### What the new dashboard deliberately does not rebuild

The observatory already has, and is competent at, all of the following. They are
out of scope for the new dashboard rather than reimplemented: host CPU, memory
and disk-rate sampling; per-task log tailing with search and pause; provenance
capture including tool version, database version and config hash; retry and
attempt history; output-validation status; the dependency-DAG canvas with
layout, pan, zoom and fit; throughput and ETA estimation; and a capability
honesty block.

**What the new dashboard does inherit, as policy rather than code: the
never-fabricate discipline.** The existing documentation is unusually good on
this — the principle that a dashboard which invents its own telemetry is worse
than no dashboard, that a null measurement means "not measurable here" and never
zero, that an absent measurement renders as an explicit marker rather than a
number, and that a capability the engine does not record is reported as
unavailable rather than guessed. That standard is adopted wholesale.

**Also inherited as a known weakness to avoid.** The existing frontend defines
colours in two independent places: CSS custom properties in one block, *and*
several parallel hard-coded hex tables in JavaScript for states, stage hues,
sparklines, attempt states and event names — which defeats its own stated
contract that one state is never two different colours. The new dashboard
defines colour once, in CSS custom properties, and derives from them.

### Data reality, stated plainly

- PDC metadata population: **966** isolates.
- Assemblies present on disk: **835**.
- Eventual analysis target: approximately **899–900** isolates on a large Linux
  machine.
- Pilot for demonstrating correctness: **10–20** genomes, selected by a
  configuration-held, seeded, phenotype-stratified cohort file that records its
  own selection rationale. The pilot is a REAL-mode run and is therefore gated.

The 20-sample laptop cap and the 10–20 pilot are consistent: the pilot sits
exactly at the cap, so the cap is exercised for real rather than being a
theoretical limit.

### Known gaps carried forward, not hidden

- The breakpoint numbers are unavailable pending a licensed CLSI M100 copy. The
  pipeline hard-fails rather than guessing. This is the single largest open item
  in D4.
- The reference file is not on disk yet; provisioning is a manual, out-of-band
  step, and the run-time length and checksum verification is untested until it
  exists.
- Which packages lack an `osx-arm64` build is not yet known. Phase 1 establishes
  it; the result determines which stages are laptop-gated.
- The MIC-reinterpretation sensitivity analysis is designed in prose only.
