# Design note: cohort variant merge (not yet implemented)

Prep for the scoping decision. **Nothing here is built**, and the stage count
question in §1 is unresolved and needs you.

## 1. The spec conflict, stated rather than resolved

`spec.md:333` is **D1. Canonical stage list — fifteen stages**, a numbered table
of fifteen rows. It is a fixed list, not a set of concerns:

- `spec.md:54` — "A Snakemake pipeline of fifteen stages, one rule per stage"
- `spec.md:618` — "900 isolates × 15 stages is up to ~13,500 cells"
- `spec.md:640` — "Three of the fifteen stages exist specifically to control…"

So a 16th stage is a **spec amendment**, not an implementation detail. The three
places that state the count would each need updating, and they state it for
different reasons (DAG size, per-stage event cost, a claim about which stages
control what).

Stage 6 `variants` is currently listed as `per-sample` (`spec.md:345`). The merge
is cohort-wide, so it does not fit under that row as it stands.

**Two ways this can go, and the choice is yours:**

- **(a) New stage, spec amended.** D1 becomes sixteen; `variants` stays
  per-sample; a new cohort row sits between `variants` and `pangenome`, or
  between `similarity` and `phenotype`. Matches snippy / snippy-core, and keeps
  the per-sample/cohort distinction crisp — which is your stated reason.
- **(b) Folded into stage 6 as a second phase.** No count change, but one Snakemake
  rule then has two different input cardinalities, which the taxonomy currently
  treats as a structural signal. `spec.md:360` shows D1 has already folded
  modules into internal steps, so the precedent exists — but those were folds,
  not new units of work.

`UNBUILT_STAGES` and the reachability mutation check in
`tests/integration/test_stage_taxonomy.py` will both need updating either way.

## 2. What the merge must actually do

The problem it solves, concretely: per-isolate calls conflate two different
things. A position where **every** isolate differs from PAO1 is a species-wide
fixed difference, not a polymorphism. A position where isolates **disagree with
each other** is cohort-informative and is what a GWAS can test. Nothing
downstream can tell those apart from the per-isolate table alone.

So the merge is a **set operation across isolates**, not a re-call:

1. Ingest every isolate's `variants.tsv` (835 files, ~1.9 M rows at 2,300
   variants/Mb × 6.3 Mb — measure it, do not trust the arithmetic).
2. Group by `(chrom, pos, ref, alt)`.
3. Keep only positions where **≥2 distinct isolate states** exist. A site called
   by every isolate, or by exactly one, carries no within-cohort contrast.
4. Emit the surviving sites with, per site: the alternate-allele frequency, the
   count of isolates carrying it, and the cohort size.

**The filter threshold is the open scientific question.** "≥2 isolates disagree"
is almost certainly too permissive — with 835 isolates, sequencing and assembly
error will make almost every position look variable. A minimum-carrier-count
(e.g. ≥5 or ≥2% of isolates) is probably right, and that number is a scientific
choice, not an implementation one. It does not exist yet and must be decided,
recorded in `science.yaml`, and justified from the allele-frequency spectrum
rather than picked.

## 3. Where it sits in the DAG

Downstream of `variants` (6). Upstream of `phylogeny` (9) and `gwas` (12), both
of which currently consume stand-ins and are listed in the handoff as still
unwired to real files. A reasonable position is between `variants` and
`pangenome`, since recombination masking (8) and phylogeny (9) both want a
cohort-level SNP set.

## 4. RESOLVED 2026-09-29: output shape is long/tidy, not wide

**`spec.md:511` is wrong, and this corrects it.** It says variants must be "in
the format pyseer expects", implying a wide genotype matrix. pyseer is installed,
so this was checked against the tool rather than against the sentence.

pyseer takes exactly one of three mutually-exclusive variant inputs:

| flag | what it wants | shape |
|---|---|---|
| `--kmers` | k-mer file | long |
| `--vcf` | a VCF; **filters any non-PASS site** | long |
| `--pres` | "Presence/absence .Rtab matrix as produced by roary and piggy" | **wide** |

The `--pres` path is the only wide one, and it is wide by construction —
`pyseer/__main__.py` reads it as:

```python
var_type = "Rtab"
header = infile.readline().rstrip()
sample_order = header.split()[1:]      # samples live in the HEADER
```

Samples are in the **header row**, so only `--pres` is wide, and the header
comment even says "rather than sample names accessible by row".

**Decision: `cohort_variants.tsv` stays long/tidy**, one row per
`(chrom, pos, ref, alt)` with per-site counts, and the pivot to Rtab belongs to
whatever builds pyseer's input. Three reasons:

- A long table survives a change of cohort. A wide matrix has one column per
  isolate, so 835 isolates is an 835-column file, and every downstream consumer
  inherits that shape.
- The merge is a set operation. Doing it in long form is where the
  "how many isolates disagree" question is answerable at all; pivoting first
  throws away the counts that make the filter auditable.
- pyseer's `--vcf` takes long input too, so if the Rtab route is abandoned
  later, the long table is still directly usable.

The pivot is a *separate* concern and should get its own stage or internal step,
declared when a consumer exists. `gwas.run()` reads `gwas_features.tsv` today, so
that consumer is not yet real and the pivot has no home yet.

## 5. RESOLVED 2026-09-29: the threshold, from a real spectrum

Measured, not picked. Four real assemblies (`GCA_000710625.1`,
`GCA_000937465.2`, `GCA_000937495.2`, `GCA_001874795.1`) aligned to PAO1 and
called with the pinned invocation. 31k–60k calls each over the whole genome:

| isolates carrying the site | count | share |
|---|---|---|
| 1/4 | 72,114 | **68.3%** |
| 2/4 | 18,120 | 17.2% |
| 3/4 | 9,401 | 8.9% |
| 4/4 (fixed) | 5,955 | 5.6% |

105,590 distinct sites total.

**The naive rule is wrong, and the numbers show exactly how.** "At least two
isolates disagree" keeps 33,476 sites — 31.7% of the genome. But:

- **68.3% of all sites are singletons**, carried by exactly one isolate. At four
  isolates that is a 25% minor-allele frequency, i.e. a private variant or an
  alignment artefact.
- The "≥2" set is dominated by its own weakest member: 18,120 of its 33,476
  sites (54%) are carried by exactly two of four.

So a threshold of two is not a *conservative* filter, it is the *opposite* — it
is defined by the singleton-adjacent cases and admits most of the tail. Note the
"fixed" column is only 5.6% at four isolates, but that number is not meaningful
yet: four isolates cannot represent a species-wide fixed difference. It would
rise sharply at 835.

**Decision: the filter is on minor-allele frequency, not a raw carrier count.**

- A site is retained when its minor-allele frequency is within
  `1/N_cohort <= maf <= max_maf`. The `>= 2 isolates` part is implicit and comes
  free, since `maf < 1/N` means exactly one carrier.
- This scales. A fixed carrier count does not: "at least 2" is 50% of a
  4-isolate cohort and 0.24% of an 835-isolate one.
- `max_maf` stays below the fixed-difference case, and the 4/4 column shows what
  that looks like.

**The number is still open, and deliberately so.** What the 4-isolate spectrum
establishes is the *form* of the rule and that a raw count is wrong. It cannot
give a `max_maf`, because a minor-allele frequency needs a real minor allele
and four isolates cannot show a rare one. **Re-measure at ~50 isolates and read
the spectrum before fixing the constant.** Guessing it from 4 would be the same
error as the PROVISIONAL columns, one level up.

For reference, pyseer already floors at `--min-af 0.01` / `--max-af 0.99`. That
is not a substitute for a filter here: 0.01 of 835 isolates is 8.35 carriers, so
pyseer's floor is *lower* than anything defensible from this data, and it applies
after the merge has already decided what a site is.

## 5b. MEASURED 2026-09-29: the floor is confirmed, the ceiling is not

**Sample.** Stratified round-robin on `Location` over `PDC_essential.tsv`, 100
isolates drawn (`papipeline.cohort_sampling.sample_stratified`), giving 47
represented strata in the first 50. Reproducible and deterministic.

**Yield.** Only **12 of 67** attempted isolates produced any calls; 55 yielded
none. That is not a pipeline failure: those assemblies are corrupt, and the
pilot census already characterised the download — of 835 assemblies, 233 were
zero-byte, 51 truncated and 426 otherwise corrupt. Two concrete shapes seen:
gzip-less binary headers prepended to FASTA, and files with 795 `>` headers and
41 bases of sequence, which minimap2 aligns to garbage. **82% of the assemblies
are unusable, so a 50-isolate usable cohort needs ~280 attempted** — which is why
this ran short.

**Spectrum at N=12**, 150,234 distinct sites:

| carriers | AF | sites |
|---|---|---|
| 1/12 | 8.3% | **63,496 (42.3%)** |
| 2/12 | 16.7% | 14,228 |
| 3/12 | 25.0% | 14,852 |
| 4/12 | 33.3% | 12,020 |
| 5/12 | 41.7% | 14,920 |
| 6/12 | 50.0% | 10,661 |
| 7/12 | 58.3% | 6,986 |
| 8/12 | 66.7% | 4,164 |
| 9/12 | 75.0% | 3,287 |
| 10/12 | 83.3% | 3,408 |
| 11/12 | 91.7% | 2,004 |
| 12/12 | 100% | 208 (0.1%) |

**The floor is now evidence, not assumption.** The spectrum is bimodal: a
singleton mode at 63,496, then a valley at k=2 (14,228), then a broad hump
peaking at k=5. The k=1 bar is **4.5x** the next and is separated from it by a
local minimum — that is the natural break, and the `1/N` MAF floor removes
exactly it. This is the same 68.3% dominance seen at N=4, so the finding
replicated across a 3x larger and a properly stratified cohort.

**The ceiling still cannot be set, and the reason is now quantified.** The
fixed-difference mode is 208 sites (0.1%) at N=12 — far too small to be
species-wide, because twelve isolates sample strain variation, not the species.
The upper tail is still rising into the 11/12 bar. A `max_maf` chosen here would
be fitted to a mode that does not yet exist. **It needs N large enough for the
fixed-difference mode to become visible**, which given 82% assembly failure
means attempting several hundred isolates.

**Both earlier bugs confirmed absent at this scale**, on real data rather than
fixtures: no row has `ac > an`; `an` equals the cohort size on every row; the
minimum `ac` is 2, so singletons are dropped; and duplicating every call of one
isolate leaves the output byte-identical, so the double-count guard holds when
counts are large.

## 5c. MEASURED to N=34: the ceiling is NOT resolvable. Defer it.

300 isolates attempted, **34 usable (11%** — lower than the 18% estimated at
N=12, so ~900 attempts would be needed for N=100). 200,501 distinct sites:

| carriers | AF | sites |
|---|---|---|
| 1/34 | 2.9% | **62,605 (31.2%)** |
| 2/34 | 5.9% | 21,909 |
| 3/34 | 8.8% | 14,924 |
| 4/34 | 11.8% | 10,360 |
| 5/34 | 14.7% | 11,010 |
| … | | broadly declining with a small hump at 12–13/34 |
| 30/34 | 88.2% | 807 |
| 31/34 | 91.2% | 989 |
| **32/34** | **94.1%** | **1,752** ← local max |
| 33/34 | 97.1% | 1,236 |
| **34/34** | **100%** | **137** |

**The floor break got sharper, and replicates a third time.** The singleton
share falls monotonically with N — 68.3% at N=4, 42.3% at N=12, **31.2% at
N=34** — and the k=1 bar remains 2.9x the k=2 bar. The `1/N` MAF floor is
confirmed at three cohort sizes.

**The ceiling did not separate, and the evidence says why.** Two observations:

1. **The 100% bar SHRANK as N grew**: 183 sites on the same 12 isolates, 137 at
   N=34. A genuine species-fixed set can only hold steady or grow as isolates
   are added. A set that shrinks is not species-fixed.
2. **The upper-tail mass sits at 94%, not 100%.** The local maximum is
   `32/34`, with 1,752 sites — 8x the 100% bar. Mass piled just *below* full
   presence is the signature of shared assembly/annotation artefacts or
   lineage structure, where each isolate carries idiosyncratic errors on top of a
   shared core, not of a clean species-wide difference.

Total mass at 94–100% is 1.6% of sites.

**Decision: `max_maf` stays unset. Documented as "not yet resolvable", not as
an oversight.** A ceiling fitted to the 94% bump would be fitting an artefact of
assembly method, and would discard genuine polymorphisms on the strength of a
mode that has not been shown to be biological. The distinction this stage exists
to make — fixed difference versus cohort polymorphism — is exactly the one this
data cannot yet make, because at N=34 no species-wide mode is visible at all.

**What would resolve it.** Either N large enough that the 100% bar stops
shrinking and grows past the 94% bump — which at 11% yield is ~900 attempts, ~1
hour — or a different route: take fixed differences from a curated source (the
reference itself, or a published P. aeruginosa core-genome set) rather than
inferring them from 34 noisy assemblies. The second is probably better science
and much cheaper, and is worth considering before spending an hour of compute.

`merge_calls` therefore continues to accept `max_maf` as an explicit argument
and to refuse when `require_max_maf=True`. Nothing downstream depends on it yet,
so deferring costs nothing today and is honest about what is unknown.

## 5d. RESOLVED: the ceiling is the floor, mirrored

`max_maf` was the wrong question. It asked the merge to identify *genuine
species-wide fixation*, and §5c showed that claim is not recoverable here: the
100% bar shrank as N grew (183 → 137) and the upper-tail mass sat at 94% rather
than 100%. Fitting a proportion to that would have been fitting assembly noise.

**Reframed, the ceiling needs no calibration at all.** The floor asks for
support on one side — `ac >= 2`, because a lone carrier is indistinguishable
from an artefact. The mirror asks for support on the other side:
**`an - ac >= 2`**, because a lone dissenter is indistinguishable from that
isolate's own error by the same argument. Same failure mode, same remedy, either
tail. Neither bound makes a claim about biology; both just refuse sites resting
on a single isolate.

Algebraically `af > 1 - 2/an`, so it reuses the frequency comparison rather than
adding a second code path, and `min_dissenters=None` disables it.

**Measured on the real N=34 cohort:** removes **1.00%** of floor-passing sites —
1,373 of 137,896 — being exactly the 33/34 bar (1,236) and the 34/34 bar (137).

**This is the right cut, and demonstrably so.** The 32/34 bar holds 1,752 sites
and is *left alone*. A ceiling fitted to the 94% bump would have deleted it
along with genuine polymorphisms. The symmetric rule removes the two near-fixed
bars, which is where the artefact-driven mass was observed, and stops there.

### Why this beats both alternatives

| alternative | why not |
|---|---|
| Curated core-genome set / published differences | population mismatch. Those differences are from *some other* population; this cohort is a specific clinical collection, and importing another population's fixation calls is a stronger assumption than the artefact this guards against. |
| ~900 attempts to N=100 | an hour of compute to answer a question that needed no calibration. The symmetric rule is correct at N=4, and re-measuring it at N=100 would confirm it rather than derive it. |

The floor was already validated by mutation testing and confirmed empirically at
N=4, 12 and 34. The ceiling is the same rule applied to the other tail, so it
inherits that validation instead of needing its own.

## 6. Contract sketch — superseded, see §4

Not pinned. The GWAS consumes it, and `gwas.run()` currently reads
`gwas_features.tsv` instead, so its real needs are not yet known. Fixing a schema
now would be the same mistake as the PROVISIONAL columns: designing from an
expectation instead of from what consumes it.

Superseded by §4 and §5. The long/tidy shape is settled; the frequency columns
are still open pending the ~50-isolate re-measurement, so no column list is
pinned yet.

`STAGE_TABLES["cohort_variants"]` currently declares
`chrom, pos, ref, alt, ac, an, af` — the minimum every plausible contract needs.
Those are provisional and no consumer reads them; §4 explains why nothing should
until a real consumer exists.

## 5. Red-first notes for whoever builds it

- The **small-N trap applies twice here.** A merge test on 3 isolates cannot
  distinguish a correct filter from one that keeps everything, because with 3
  isolates almost every site is "variable". Include a case where a site is
  called by *all* isolates and must be **dropped** — that is the whole point, and
  it is invisible unless asserted.
- Assert the **sum of allele counts equals cohort size × ploidy**, or whatever
  the true invariant is. That catches a filter that silently drops a row.
- Assert cohort-size-agnosticism by running the same fixture at 3, 25 and 835
  isolates. A hard-coded cap here would be invisible at the sizes used so far.
