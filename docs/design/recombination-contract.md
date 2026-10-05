# Design — stage 8, recombination masking (Gubbins)

Status: written **before** the implementation, as the reviewable artifact. Every
empirical claim in it is annotated with the command that produced it and what
that command does and does not prove.

Owner of the code: `papipeline/stages/recombination.py` (the contract),
`papipeline/adapters/gubbins.py` (the process, the binary, the file formats).

---

## 1. What this stage is for

`spec.md` D1 puts stage 8 between the pangenome (7) and the phylogeny (9). It
takes the core gene alignment, asks Gubbins where recombination happened, and
produces two things the pipeline can use:

1. a **recombination-masked polymorphic-sites alignment**, which stage 9
   (`phylogeny`) and stage 12 (`gwas`) read as their SNP input; and
2. a **per-node table** saying which branches of the inferred tree carry
   recombination, which is the stage's own scientific output.

It does **not** produce the tree. Gubbins builds one internally; that tree is
kept as an extra output for inspection and is not the contract output. Stage 9
runs its own inference on the masked alignment. See §9.

---

## 2. Inputs

### 2.1 The alignment

**`core_gene_alignment_filtered.aln`**, under
`intermediate_root / "panaroo" / …` — the path is owned by
`papipeline.adapters.panaroo.output_dir(intermediate_root)`, not restated here
or in the workflow.

**Why the *filtered* alignment and not `core_gene_alignment.aln`.** Panaroo
writes two alignments. The unfiltered one is every core-gene column it aligned;
the filtered one is what survives panaroo's own model-based column filtering.
The evidence that the filtering did something is in the file sizes:

```
$ for f in core_gene_alignment.aln core_gene_alignment_filtered.aln; do ... ; done
core_gene_alignment.aln            seqs=10  sites=4758772  len_per_seq=4758772
core_gene_alignment_filtered.aln   seqs=10  sites=4645774  len_per_seq=4645774
```

10 sequences in both; 4,758,772 columns unfiltered, 4,645,774 filtered —
**112,998 columns (2.37 %) removed**, recorded by panaroo as
`core_alignment_header.embl` versus `core_alignment_filtered_header.embl`.

This matters because of what Gubbins counts. Gubbins derives its SNP set from
the alignment it is given and then computes recombination statistics over
windows of those SNPs. Feeding it columns panaroo has already judged
unreliable would inflate the SNP count those statistics are computed over, and
the recombination call is a statement about SNP density. The pipeline applies
two filters in series, each for a different reason: panaroo removes
**alignment-quality** columns (its entropy/HMM model), Gubbins removes
**recombination** signal. Passing Gubbins the unfiltered alignment would apply
the first filter twice and never apply the second properly.

**What this does not prove:** that panaroo's filtering is *correct*, only that
it is non-trivial and runs in the direction that reduces the SNP count. It also
does not prove the filtered alignment is the one a user would want for a
different purpose — for a *gene-content* question the unfiltered alignment is
the more faithful record.

### 2.2 Threads

From `config.runtime.threads` (the machine overlay). Never a literal. Hard
rule 2.

---

## 3. Which Gubbins binary — and why this is the load-bearing decision

### 3.1 Two builds exist and one of them is broken

| Build | Path | Self-reported version | On `tiny.aln` |
|---|---|---|---|
| source | `/Users/raghavkrishnankv/Desktop/pa-gubbins/tools/gubbins/bin/gubbins` | **3.4.2** | **exit 0** |
| conda (`osx-arm64`) | `…/micromamba/2.9.0/envs/pa-gubbins/bin/gubbins` | **3.4.3** | **exit 139** |

```
$ /Users/…/pa-gubbins/tools/gubbins/bin/gubbins tiny.aln ; echo EXIT=$?
EXIT=0
$ /opt/homebrew/Cellar/micromamba/2.9.0/envs/pa-gubbins/bin/gubbins tiny.aln ; echo EXIT=$?
EXIT=139
```

The conda build's `libgubbins.0.dylib` references `_gzopen`/`_gzread`/`_gzclose`
and cannot resolve them, so it dies with SIGSEGV. That is a bioconda
`osx-arm64` recipe defect, not a Gubbins source defect: the same 3.4.3 sources
build and link correctly through
`scripts/build_gubbins_from_source.sh`. Filed upstream as
nickjcroucher/gubbins#452. `tests/integration/test_gubbins_source_build.py`
already guards this.

The two builds report **different version strings (3.4.2 vs 3.4.3)**, so the
string is a usable discriminator. Note that `run_gubbins.py`'s own banner
prints `--- Gubbins 3.4.3 ---` on *both* builds — that is the **Python
package's** version, not the C binary's, so it must never be used to decide
which binary ran.

### 3.2 An explicit path is not enough — PATH isolation is required

`gubbins/common.py:79-89`:

```python
os.environ["PATH"] = os.environ["PATH"] + ":/usr/lib/gubbins/"
gubbins_exec = 'gubbins'
if utils.which(gubbins_exec) is None: ...
```

and `utils.which` is a plain left-to-right scan of `os.environ["PATH"]`. So
`run_gubbins.py` accepts no flag that names a binary, and its only fallback is
*appending* `/usr/lib/gubbins/` — which can never win against an entry that is
already earlier in PATH. **Prepending a directory is therefore the only
reliable override.**

The adapter builds a per-run directory containing exactly one entry — a
symlink to the source binary — and puts it first:

```
<workdir>/bin/gubbins -> <resolved source binary>
PATH = <workdir>/bin : <conda env bin> : <rest>
```

`utils.is_executable` uses `os.path.isfile`, which follows symlinks, so a
symlink is accepted. The conda env's `bin` stays on PATH *after* the isolated
directory, which is required: Gubbins needs `raxml-ng`, `pyjar` and friends
from there, and the point is to shadow *one* executable, not to empty PATH.

### 3.3 The PATH must also contain `/usr/sbin`

A minimal isolated PATH breaks Gubbins's RAxML-NG tree builder:

```
File ".../gubbins/treebuilders.py", line 753, in select_executable_based_on_threads
    single_threaded_exec = utils.choose_executable_based_on_processor(
File ".../gubbins/utils.py", line 118, in choose_executable_based_on_processor
    flags = output.lower().split()
UnboundLocalError: local variable 'output' referenced before assignment
```

`choose_executable_based_on_processor` probes CPU features by running
`sysctl -a | grep machdep.cpu.features` through `subprocess.Popen(shell=True)`.
If `sysctl` is not on PATH the `elif` never fires, `cpu_info` stays `False`,
and the `output` variable is **never assigned** — but the very next
unconditional statement reads it. Adding `/usr/sbin` fixes it; verified by
re-running the same command with `/usr/sbin` present, which then exits 0.

This is a Gubbins bug (a local read outside the branch that sets it) that this
adapter has to route around, and it is worth knowing that the workaround is
"keep a normal system PATH", not "make the PATH more minimal".

---

## 4. Tool configuration

| Setting | Value | Why |
|---|---|---|
| builder | `raxmlng` | `--tree-builder` takes the literal `raxmlng`, no hyphen (verified against `run_gubbins.py --help`). The **default is `raxml`**, and RAxML is absent from the `pa-amr` environment. |
| threads | `config.runtime.threads` | hard rule 2 |
| prefix | `recombination` | `--prefix` renames outputs; it does **not** relocate them |
| working dir | per-run directory under `intermediate_root` | see §6 |

**Correction to a briefed fact, recorded because it was load-bearing:** the
briefing said `raxml` is absent from the environment. It is present in
`pa-amr` (`/opt/…/envs/pa-amr/bin/raxml`) but **absent from `pa-gubbins`**,
which is the environment that actually holds `run_gubbins.py` and
`raxml-ng`. The conclusion (use `raxmlng`) is unchanged; the reason is
different, and the difference matters because it identifies *which* environment
the adapter must run in.

**`run_gubbins.py` has no output-directory flag.** Verified: the only output
flag is `--prefix`, and `common.py:1543-1551` shows `--prefix` maps input names
to output names by *string substitution* while everything is still written
relative to the process CWD. So the adapter must set `cwd=` to a scratch
directory.

---

## 5. Outputs

### 5.1 `core_snp_alignment.fasta` — the masked alignment

Copied from Gubbins' own `<prefix>.filtered_polymorphic_sites.fasta`
(`common.py:1545`). Written to **`phylogeny_dir / "core_snp_alignment.fasta"`**.

That location is not a choice available to this stage — it is where the
consumers already look:

* `stages/phylogeny.py:368` — `alignment = phylogeny_dir / "core_snp_alignment.fasta"`,
  passed straight to `adapters.iqtree.build_tree`;
* `stages/gwas.py:502` — records `f"core_snp_alignment:{self.snp_alignment.name}"`.

A file written anywhere else would be an alignment no consumer reads.

Format, verified on both a 5-taxon probe and the real 10-isolate cohort:

* plain FASTA, one `>taxon` record per sequence;
* **variable sites only** — the probe's 400-site input yields 57 columns, and
  the set of 57 matches the 57 records of `<prefix>.summary_of_snp_distribution.vcf`
  exactly (checked position-by-position, no discrepancy);
* ungapped, unambiguous — the character set is exactly `{A, C, G, T}`;
* all sequences the same length.

### 5.2 `recombination.tsv` — the per-node table

The stage's declared contract table, `execution/contracts.py`
`STAGE_TABLES["recombination"]`:

```
node	n_snps	mean_branch_length	recombination_detected
```

Mapped from Gubbins' `per_branch_statistics.csv` and the node-labelled tree:

| Column | Source | Meaning |
|---|---|---|
| `node` | `Node` column of `per_branch_statistics.csv` | tip label for a tip, `Node_<n>` for an internal node. The two files join exactly: 19 rows for 10 tips + 9 internal nodes. |
| `n_snps` | `Total SNPs` | SNPs assigned to that branch. |
| `recombination_detected` | `Number of Recombination Blocks > 0` | `1`/`0`. |
| `mean_branch_length` | **derived — see below** | |

**`mean_branch_length` is not in Gubbins' output.** `per_branch_statistics.csv`
has thirteen columns and none of them is a branch length; the columns are SNP
counts, block counts, base counts, `r/m`, `rho/theta`, genome length and clonal
frame. So this stage computes it, from
`<prefix>.node_labelled.final_tree.tre`:

> **`mean_branch_length` = the arithmetic mean of the lengths of all edges in
> the clade subtended by the node, including the node's own edge.**

For a tip that is just its own edge length; for an internal node it is the mean
over its whole subtree. It is defined for every row, and it is the quantity you
would plot against `n_snps` to see whether long branches carry the
recombination.

*Alternative considered and rejected:* the mean of a node's **immediate
children's** edges. It is the more common phrasing of "mean branch length",
but it is undefined for the ten tip rows — a tip has no children — so half the
table would be empty. Empty cells in a branch-length column are exactly the
kind of hole a later reader fills with a zero.

**`node` is not `sample_id`.** Internal nodes are `Node_<n>`, so this table is
*not* a per-sample table and is deliberately absent from
`DENSE_PER_SAMPLE_STAGES`. A caller that wants per-sample recombination must
map through the tree; this stage does not do that mapping and does not pretend
the two are the same thing.

### 5.3 Gubbins' tree — an extra output, not the contract output

Kept at `<intermediate_root>/gubbins/recombination.final_tree.tre` and
`…node_labelled.final_tree.tre`. **Not** renamed to `tree.nwk` and **not**
copied to `phylogeny_dir`. Stage 9 writes `tree.nwk` from its own IQ-TREE
inference. Copying Gubbins' tree there would make two files claim the same
path, and a downstream reader could not tell which inference produced the tree
its patristic distances are measured on.

### 5.4 Provenance

`<intermediate_root>/gubbins/provenance.json`:

* the resolved binary path, and the *isolated PATH directory* used;
* the binary's self-reported version string, read from the binary;
* the full `argv`;
* the exit code;
* wall-clock seconds;
* the input alignment path and its length in columns.

Written **whether the run succeeded or failed**. Provenance that only exists on
success cannot answer the question you most want answered, which is what ran
when something went wrong.

---

## 6. gubbins writes to CWD

`run_gubbins.py` has no output-directory flag and `--prefix` only renames
(§4). So the adapter creates a fresh per-run directory and sets `cwd=` to it.
The directory lives under the run's own `intermediate_root`, so it is removed
with the rest of the run and is never inside the repository working tree.

This is also why the adapter never runs Gubbins with the repository as CWD: a
run would scatter ~40 files of intermediate products into a tracked directory.

---

## 7. Refusals — three, each with its own named reason

All three are REAL-mode refusals. Each has a **distinct** message; no message
is reused for two causes, because a message that covers two failures cannot
tell an operator which one they have.

### (a) The REAL gate is shut

`runtime.allow_real_mode` is false in every committed overlay. Checked **first**,
before anything else, because an operator whose gate is shut does not need to
be told their alignment is missing — that is not their problem and fixing it
would not unblock them.

Wording follows the house style at `stages/similarity.py:377`: names the flag,
names the environment-variable override, names the overlay.

> REAL-mode stage 8 (recombination) is gated: `runtime.allow_real_mode` is
> false in the machine overlay for 'laptop'. … Set it in the overlay, or open it
> for one session with `PIPELINE_ALLOW_REAL_MODE=1`.

### (b) The alignment is missing

Names the expected path and what writes it.

> no core gene alignment at `<path>`. It is written by stage 7 (panaroo) …

A missing alignment is a **missing input, not an absence of recombination.** The
distinction is the whole point of refusing: a stage that returned an empty
block table here would be reporting "no recombination was found" when in fact
it looked at nothing.

### (c) The resolved binary fails its self-test

This is the one refusal that is specific to this machine, so it is the one that
must name the actual defect.

Before the real run, the adapter executes the resolved binary on a tiny
generated alignment in a scratch directory. If it exits **139** (SIGSEGV) the
run is refused, and the message says explicitly that this is the **conda
`osx-arm64` build and exit 139**:

> the gubbins binary resolved to `<path>` exits 139 (SIGSEGV) on a 5-taxon
> 400-site probe alignment. That is the bioconda `osx-arm64` build of gubbins
> 3.4.3, whose `libgubbins.0.dylib` cannot resolve `_gzopen`/`_gzread`/
> `_gzclose`. … Use the source build.

A non-139 non-zero exit is a *different* failure and says so, rather than being
folded into the segfault message.

---

## 8. Do not re-root the tree

Gubbins **already midpoint-roots** its own output tree. Verified: the real
10-isolate tree's root edge is `:0.0`, and `dendropy.reroot_at_midpoint()`
succeeds on it.

Re-rooting a tree that is already rooted is not harmless. It rebuilds the seed
edge and the bipartitions, and on a tree whose root edge has zero length it
walks a path with no positive edge to stop on.

### 8.1 The failure mode, precisely

`dendropy` `Tree.reroot_at_midpoint` contains, in the walk that finds the
break point:

```
assert break_on_node is not None or target_edge is not None
```

* dendropy **5.1.0**: `dendropy/datamodel/treemodel/_tree.py:2699`, inside
  `reroot_at_midpoint` (defined at line 2632).
* dendropy **4.6.1** (the version bundled with Gubbins):
  `_tree.py:2685`, same function, same assertion.

So this is **not a 5.x regression**. If a future tree triggers it, downgrading
or upgrading dendropy is not the fix. What rescues the real tree is that the
root edge has **zero length**, so the walk finds a valid break node on the
first step and never reaches the assert.

When it does fail, the symptom is a **bare `AssertionError`** with no message,
no node name, and no indication of which tree or which edge caused it —
because `assert` carries no message here and Python's default `AssertionError`
str is empty. In a Snakemake run it surfaces as an opaque traceback ending in
`AssertionError` in `_tree.py`, which reads like a library bug rather than a
property of the input.

**Integer rounding is not an acceptable workaround.** Rounding branch lengths
to integers "fixes" float-residue comparisons by making every length exactly
representable — and it trades a loud crash for silently wrong branch lengths,
which then propagate into stage 10's patristic distances, where a wrong length
is indistinguishable from a real distance. A crash is recoverable; a distance
matrix that is wrong by an unknown amount is not.

---

## 9. Relationship to stage 9, and what the tree is not

Gubbins builds a tree internally because it needs one to assign SNPs to
branches. That tree is a by-product of the recombination analysis, not the
phylogeny this pipeline reports. Stage 9 re-infers on the masked alignment
(IQ-TREE, `GTR+G`). Two trees, two methods, two files, no shared path.

---

## 10. TEST mode

TEST does not run Gubbins. It runs the **real parser and the real writer** over
a directory of Gubbins-shaped outputs supplied by the caller, which is the same
shape as `stages/similarity.run(tree_path=…)`: the real code path, a committed
input, no external tool.

`gubbins_dir` is **required** in TEST and has **no default**. The reasoning is
the one already written at `similarity.py:356`: resolving a fixture from
configuration would let a TEST run measure something nobody chose.

---

## 11. The stage module's own REAL refusal — read this

`run.py` lists `recombination` in `UNBUILT_STAGES` and has **no dispatch
branch** for it. `papipeline/run.py` is off-limits for this ticket (another
agent owns a hunk of it), so the stage module cannot be promoted into the
dispatch chain here.

`stages/recombination.py::run()` therefore **raises `NotImplementedError` in
REAL**, naming the adapter as the working path. This is not a workaround of
`tests/integration/test_stage_taxonomy_is_self_verifying.py`; it is what that
test's own docstring says the state means — *"A stage declared unbuilt must
either have no module or have a `run()` that refuses, because 'no caller' and
'a caller nobody invokes' are the same thing observed from two places."* There
is genuinely no caller on the `run_pipeline` path, so `run()` genuinely refuses.

`papipeline/adapters/gubbins.py` is the working REAL path and is fully
exercised by the tests, including against the real 10-isolate cohort. Promoting
the stage needs a `run.py` dispatch branch; proposed text is in the handoff.

`UNBUILT_WITH_STANDIN` and the TEST stand-in fixture are left in place for the
same reason: retiring them is part of the same `run.py` change.

---

## 12. Cost

### 12.1 Measured (n = 10, this machine)

| Quantity | Value | How obtained |
|---|---|---|
| input | 10 taxa × 4,645,774 columns | the filtered alignment |
| polymorphic sites | 71,165 | counted by column over called bases |
| Gubbins iterations | 5 (the maximum) | `Maximum number of iterations (5) reached.` |
| wall clock | **31.64 s** | `/usr/bin/time -l`, 4 threads |
| peak RSS | **245,972,992 B = 234.6 MiB** | `/usr/bin/time -l`, `maximum resident set size` |
| Gubbins' own total | 30.43 s | its log |

### 12.2 Extrapolated to n = 967 — **extrapolation, not measurement**

Nothing below was measured at 967. The model, and its weak points:

The cost has three parts, and they scale differently.

1. **RAxML-NG tree building and likelihood**, repeated once per iteration.
   RAxML's likelihood is `O(N · L)` per evaluation and its search is
   super-linear in `N`; in practice RAxML-family tree building is conventionally
   treated as roughly `O(N^2)`–`O(N^3)` in taxon count. Taking the midpoint
   `(N/L²)` scaling — 967/10 ≈ 97× in `N`, `L` fixed — gives roughly
   **10³–10⁴ ×** the tree-building cost.
2. **Ancestral reconstruction (pyjar + RAxML-NG `--ancestral`)**, `O(N · L)`.
   ~97 × the per-iteration cost.
3. **The iteration count.** Gubbins stops on convergence and defaults to a
   maximum of 5. It **hit 5** at n=10 — i.e. it did not converge, so 5 is a
   floor here, not a ceiling. With more taxa the recombination signal is
   stronger and convergence is usually *faster* in iteration count, but that is
   an expectation, not a measurement, and the run above is not evidence for it.

Combining: **hours, not minutes, on a single node** — order 10⁴ s at the low
end of the model. Memory is the gentler axis: the alignment matrix is
`N × L` bases, so 967 × 4.6 Mbp ≈ 4.5 GB as text and ~0.45 GB as 2-bit packed,
which puts a single-machine run in the **tens of GB** range for RAxML's model
and bootstrap arrays.

**Honest statement of the uncertainty:** the dominant term is the tree search,
which is the term I modelled rather than measured, and RAxML's real scaling
departs from any power law on real data. **This estimate should be treated as
"it will not finish this afternoon", not as a number to plan against.** The
defensible next step is to measure the scaling directly — run n = 10, 20, 40, 80
from the same cohort and fit — not to trust the power law.

### 12.3 Compared with a reference-mapped alignment from stage 6

Stage 6's cohort merge (`stages/cohort_variants.py`, `cohort_variants.tsv`
with `chrom, pos, ref, alt, ac, an, af`) already yields every variable position
in the cohort, called against one reference. A SNP-distance matrix from it is a
gather and a count over an `N × M` table — no tree, no ancestral reconstruction,
no iteration, no convergence criterion.

For n = 967 and M ≈ 10⁵–10⁶ SNPs that is ~10⁸ cells: **seconds, in memory, on
one core**, with no external tool beyond the variant caller stage 6 already
runs.

The asymmetry is the finding:

| | Gubbins on a core-gene alignment | Reference-mapped from stage 6 |
|---|---|---|
| Cost at n = 967 | hours, tens of GB, one iteration loop | seconds, ~1 GB, no loop |
| Gives recombination blocks | yes | no |
| Gives a distance matrix | indirectly, via the tree | directly |
| Gives a tree | yes (but not the one we report) | no |
| Masking | removes recombinant signal | not applicable |
| Detects *within-lineage* recombination | yes | yes, but only against the reference's allele set |

**The recommendation this supports:** the two are not substitutes. Gubbins'
unique contribution is the recombination *blocks*, and that is a per-branch
quantity that a reference-mapped SNP table cannot produce, because it has no
tree to hang a branch on. What the reference-mapped route buys is a
distance matrix and GWAS features at a cost that is not the bottleneck. If the
967-isolate study needs recombination *blocks*, Gubbins on a core-gene alignment
is the tool, and the cost is the price — but it should be an explicit,
measured decision, not the default because the alignment happens to be lying
around.

There is also a **scientific** caveat on the Gubbins route that the cost
argument does not cover, and it is larger: the cohort the cost was measured on
does not behave like one population (§13).

---

## 13. What the real run showed about the data itself

Recorded because it constrains how any of this may be reported. Full matrix in
the handoff; the shape:

* 45 pairs. **min 0, median 24,803, max 38,474.**
* **Exactly one pair at zero**: `PDT000294805.1` / `PDT000311294.1`.
* Two tight groups (≤ 169 SNPs within) sitting ~19,000–38,500 SNPs from
  everything else.

A median pairwise distance of ~24,800 SNPs across a 4.6 Mb *core* alignment is
not a population of near-clones; it is **three divergent lineages** in one
file. A phylogeny over it is a statement about lineages, and recombination
blocks computed on it mix within-lineage and between-lineage signal.

### 13.1 The zero-SNP pair, explained

`PDT000294805.1` and `PDT000311294.1` differ at **zero** columns where both
carry a called base. They are **not** identical sequences: they differ at **77
columns, every one of them gap placement** (one has a gap where the other has a
base), and their gap counts differ (665 vs 650).

That is why `per_branch_statistics.csv` gives both of them `Total SNPs = 0`:
Gubbins assigns a SNP to a branch only where the two taxa differ by a base at
that column. There is no such column. The ~20,000-SNP difference between this
pair *as a clade* and the rest of the cohort sits entirely on their parent
branch `Node_1` (`Total SNPs = 19,992`).

The likely mechanism is that the two assemblies agree on sequence but not on
gene boundaries — a start or stop one base different, or a flanking segment
included on one side only. **That is a hypothesis, not a finding.** Confirming
it means diffing the two assemblies directly, which this ticket does not do.

**Consequence for the table:** a `recombination_detected = 0` on a branch is
not evidence that recombination was looked for and not found. For these two
taxa it means Gubbins had no SNP to assign. Any report that counts
"branches with no recombination" must exclude or annotate zero-SNP branches, or
it will report the absence of a measurement as the absence of recombination.

---

## 14. Summary of decisions and their justifications

| # | Decision | Justification |
|---|---|---|
| 1 | Input = panaroo's **filtered** core alignment | panaroo filters alignment *quality*, Gubbins filters *recombination*; two filters, two purposes. §2.1 |
| 2 | Binary resolved by **isolated PATH prepend**, not a path argument | `run_gubbins.py` accepts no binary-path flag and its fallback only *appends* to PATH. §3.2 |
| 3 | Self-test on a probe alignment; exit 139 ⇒ refuse | a segfault is otherwise discovered 30 s in, after the real run has started. §7c |
| 4 | Conda env `bin` kept **after** the isolated dir | Gubbins needs `raxml-ng`/`pyjar` from it; the goal is to shadow one executable. §3.2 |
| 5 | `/usr/sbin` retained in PATH | otherwise `choose_executable_based_on_processor` raises `UnboundLocalError`. §3.3 |
| 6 | Builder `raxmlng` | the default `raxml` is absent from the env holding `run_gubbins.py`. §4 |
| 7 | Run in a per-run working dir; `--prefix` only renames | Gubbins has no output-directory flag. §6 |
| 8 | Masked alignment written to `phylogeny_dir` | that is where `phylogeny.py:368` and `gwas.py:502` read it. §5.1 |
| 9 | `mean_branch_length` = mean over the node's whole clade | the only definition that is defined for tip rows *and* internal rows. §5.2 |
| 10 | Gubbins' tree kept as an extra, not renamed to `tree.nwk` | stage 9 writes that path; two files must not claim it. §5.3 |
| 11 | Provenance written on failure too | provenance that exists only on success cannot answer "what ran?". §5.4 |
| 12 | Three refusals, three distinct messages | one message covering two failures cannot identify which occurred. §7 |
| 13 | **No re-rooting** | Gubbins already midpoint-roots; re-rooting a zero-length root edge is what trips the bare assert. §8 |
| 14 | No integer rounding as a dendropy workaround | trades a crash for silently wrong branch lengths. §8.1 |
| 15 | TEST requires an explicit `gubbins_dir`, no default | otherwise a TEST run measures something nobody chose. §10 |
| 16 | `run()` refuses REAL; adapter is the working path | `run.py` has no dispatch branch and is off-limits; the refusal states the true reason. §11 |