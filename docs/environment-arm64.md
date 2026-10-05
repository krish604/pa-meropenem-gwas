# Environment and tool availability on osx-arm64

Phase 1 output, ticket 02. Every claim below was checked against a command, and
§7 was verified by executing each tool. Where a conclusion was reached wrongly
and later corrected, the correction is recorded in place with the evidence that
overturned it — §3 and §5 are the two places that happened, and both are
labelled. Read those two as warnings as much as findings.

Machine: MacBook Air M2, arm64, 16 GB RAM, macOS.
Env manager: `micromamba 2.9.0`. `MAMBA_ROOT_PREFIX` is unset, so environments
live under the binary's own prefix and `micromamba run -n <env>` works with no
extra setup.

> **Gotcha worth knowing.** A `pilot100` environment exists under
> `~/micromamba_pilot/envs/`, which belongs to a *different* root prefix that is
> no longer active. `micromamba run -n pilot100` therefore fails with "the given
> prefix does not exist" even though `micromamba env list` shows the env. That
> is historical and does not affect `pa-amr`, which lives in the active prefix.

---

## 1. Corrections to the previous `environment.yml`

Four pins were wrong. Three of these were flagged as suspects before this phase
began, and all three were confirmed wrong.

| Was | Now | Why it was wrong |
|---|---|---|
| `amrfinder=3.12.0` | `ncbi-amrfinderplus=4.2.7` | `amrfinder` does not exist in any channel. The package is `ncbi-amrfinderplus`. |
| `blast=4.14.1` | `blast=2.17.0` | BLAST+ is versioned 2.x. There is no 4.x. |
| `iqtree=2.2.6` | `iqtree=3.1.3` | `iqtree2` does not exist. The package is `iqtree`, now shipping the IQ-TREE 3 codebase. |
| `pyseer=0.4.5` | `pyseer=1.1.2` | `0.4.5` does not exist. 1.1.2 is the **highest conda-installable** version; see §4. |

The previous pins also described a Python stack that never existed:
`numpy=1.26.4`, `pandas=2.2.2`, `pytest=8.2.2` and `snakemake-minimal=8.20.5`.
The environment that has actually run this pipeline (`pilot100`) uses
`numpy 2.4.6`, `pandas 3.0.6`, `pytest 9.1.1`, `snakemake 9.x`. The new pins
match that known-good baseline, so the stack is now described accurately.

`defaults` was removed as a channel. It resolves to `repo.anaconda.com`, which
failed DNS repeatedly during the build, and every package resolves from
`conda-forge` or `bioconda` without it.

## 2. Availability summary

| Stage | Tool | osx-arm64 | Evidence |
|---|---|---|---|
| 1 validate | seqkit 2.14.0 | available | native osx-arm64 build |
| 2 annotate | bakta 1.12.1 | available | noarch |
| 3 mlst | mlst 2.33.1 | available | noarch |
| 4 amr | ncbi-amrfinderplus 4.2.7 | available | native osx-arm64 build |
| 5 virulence | blast 2.17.0 | available | native osx-arm64 build |
| 6 variants | minimap2 2.31 | available | native osx-arm64 build |
| 6 variants | samtools 0.1.19 | available **but constrained** | see §3 |
| 6 variants | bcftools 1.23.1, htslib 1.23.1 | available | native osx-arm64 build |
| 7 pangenome | **panaroo** | available **via pip only** | conda recipe unsolvable, see §3 |
| 8 recombination | **gubbins** | **UNAVAILABLE here** | py310-only builds, see §3 |
| 8 (dependency) | mafft 7.526 | available | native osx-arm64 build |
| 9 phylogeny | iqtree 3.1.3 | available | native osx-arm64 build |
| 10 similarity | derived from the stage-9 tree | available | internal |
| 11 phenotype | internal | available | internal |
| 12 gwas | pyseer 1.1.2 | available **but limited** | see §4 |
| — | snp-sites 2.5.1 | available | native osx-arm64 build |

**One tool cannot run on this machine: `gubbins`.** Stage 8 (recombination
filter) is therefore impossible here. Because stages 9–15 consume stage 8, **the
laptop still cannot perform a full REAL run of the 15-stage pipeline.** That is
a property of the platform.

**`panaroo` (stage 7, pangenome) now works on this machine**, but **not through
conda**. It is installed from source via pip; see §3. Stage 7 is therefore no
longer a blocker. Stage 8 remains the only one.

This is the situation the planned machine-overlay tool-availability key exists to
express. **That config does not exist yet** — no `config/laptop.yaml` is present
in the repository, and no DAG-build refusal has been implemented. It is the
deliverable of ticket 04 (and 11 for the DAG-build enforcement). Until then, the
availability table in §2 is the only record of what can run here.

## 3. Why the conda gubbins package is unusable, and how panaroo is installed anyway

**`gubbins`** — an earlier draft of this document claimed gubbins had *zero*
`osx-arm64` builds. **That was false**, and it is corrected here rather than
quietly removed. gubbins has three `osx-arm64` builds of version 3.4.3:

```
$ micromamba search -c conda-forge -c bioconda --platform osx-arm64 --json gubbins
3.4.3  py310hdfa5cb7_1   osx-arm64
3.4.3  py39h67fc90f_0    osx-arm64
3.4.3  py310h483ef05_0   osx-arm64
```

and `gubbins=3.4.3` alone solves for `osx-arm64`. The real constraint is the
Python version:

| Solve | Result |
|---|---|
| `gubbins=3.4.3` (alone) | solves, drags `python 3.10.21`, `perl 5.32.1`, `raxml` |
| `gubbins=3.4.3` + `python=3.10` | **solves** — 272 packages, 170 MB |
| `gubbins=3.4.3` + `python=3.11` | **fails** — `requires scipy =*, which does not exist` |
| `gubbins=3.4.3` in this env (`python=3.11.16`) | **fails to solve** |

So every `osx-arm64` gubbins build is Python-3.10-only, and this codebase is
Python 3.11.16. Gubbins cannot be installed here without downgrading the entire
scientific stack, which would break the laptop/big-machine comparability the
config split exists to protect.

**The Python constraint above is real but no longer the operative one.** gubbins
*is* installable in a separate `osx-arm64` environment, and when installed it
**crashes immediately**. Verified on the newest build,
`gubbins-3.4.3-py310hdfa5cb7_1`, in a clean env created with
`--platform osx-arm64`:

```
$ run_gubbins.py --prefix smoke --starting-tree smoke.tree.nwk smoke.fasta
  Running Gubbins to detect SNPs...
  Gubbins crashed, please ensure you have enough free memory

$ gubbins smoke.fasta ; echo $?
139                      # SIGSEGV
```

The only fault evidence recorded here is the **exit status**: `139`, i.e. killed
by SIGSEGV. **No debugger backtrace of this crash is stored in this
repository**, so nothing here claims which instruction faulted or what address
it read. What follows is the set of candidate causes that were tested, and what
each test showed — not an identification of the fault.

Three candidate causes were tested. Two are ruled out on the evidence recorded
beside them. The third — the binary's linkage — is recorded in the last row
together with what was observed and what that observation does not establish.

| Hypothesis | Verdict |
|---|---|
| Rosetta / x86_64 translation | **No.** `gubbins`, `raxml-ng`, `FastTree`, `iqtree2` and `node` are all native `arm64`. Nothing is translated. |
| numba / llvmlite JIT codegen bug | **No.** `NUMBA_DISABLE_JIT=1` reproduces the identical 139. numba 0.67.0 / llvmlite 0.49.0 are present. |
| Unresolved `gz*` symbols at runtime (e.g. libz) | **This is the observed signature, and it is the opposite of what this row previously said.** `nm -u` reports `_gzopen`, `_gzread` and `_gzclose` undefined in the conda build's `libgubbins.0.dylib`, and `libSystem` exports no `gz*` at all (`nm -gU /usr/lib/libSystem.B.dylib \| grep -c gzopen` → `0`), so nothing in that process can satisfy them. Note that `otool -L` is **not** evidence either way: it lists what a file *links*, not what its undefined symbols can *resolve to*. |

The crash is input-independent (a 3-genome, 10-SNP alignment is enough to
trigger it) and reproduces in a pristine environment. On **one** identical
alignment with **one** identical argv, differing only in the binary, the conda
3.4.3 build exits 139 and writes nothing while the **source build of the same
tag** exits 0. So the conda **package** is unusable on `osx-arm64`; gubbins
itself is not, because a source build of the same version runs on this machine
(`config/machines/laptop.yaml`). **Relinking the conda binary is untested** —
it has never been tried, so no claim is made here about whether it would help.
The source build is a rebuild from source, which is a different thing from a
relink. Filed upstream; see the URL in `scripts/build_gubbins_from_source.sh`.

An earlier draft of this investigation attributed the crash to a missing
`LC_LOAD_DYLIB` entry and proposed an `install_name_tool` relink. That proposal
was **withdrawn**, and the relink has never been carried out, so there is no
evidence about it in either direction and none is claimed here. What is on
record instead is the comparison above: on one identical input and one
identical argv, the conda binary exits 139 while the source build of the same
tag exits 0.

**How the error happened, because it is the interesting part.** An early
cross-platform sweep in this phase *did* report three `osx-arm64` builds for
gubbins. A later single `micromamba search` returned nothing, and that later
result was believed in preference to the earlier contradicting evidence. The
lesson is recorded in §8 as well: when two sources disagree, the disagreement is
the finding. A silently-empty search result is not evidence of absence.

**`panaroo`** — **this section previously said the package was `noarch` and
"present everywhere". That was wrong.** panaroo *does* have native `osx-arm64`
builds on bioconda (1.6.0, 1.7.0, 1.8.0), and `intbitset` has a working
`osx-arm64` build. What fails is the **solve**, because the recipe carries a
hard `prokka` dependency that cannot be satisfied on arm64:

```
$ micromamba create -n panaroo-test --platform osx-arm64 -c bioconda -c conda-forge panaroo
Could not solve for environment specs
  panaroo [1.6.0|1.7.0|1.8.0] would require
    prokka =* * but there are no viable options
      prokka [1.14.6|1.15.6] -> tbl2asn-forever >=25.7, which does not exist
      prokka [1.13.7|1.14.0|1.14.5|1.14.6] -> perl >=5.26.2,<5.26.3.0a0, which does not exist
  panaroo [1.0.0|1.1.0|...|1.5.2] would require
    mkl =* *, which does not exist
```

Zero `osx-64` packages are pulled — this is a dependency-metadata problem, not
an architecture problem.

**The prokka dependency is spurious, and bypassing it is correct, not a hack.**
`panaroo/prokka.py` is a **GFF3 *reader***. It parses GFF3 text and contains no
subprocess call to prokka. panaroo shells out to *cd-hit*, which installs
cleanly on `osx-arm64` (verified, 4.8.1). So the conda recipe hard-requires a
tool that panaroo never invokes at runtime, and that is the whole reason the
solve fails.

Installing from source avoids the recipe entirely:

```
$ python3 -m venv panaroo-venv
$ panaroo-venv/bin/pip install git+https://github.com/gtonkinhill/panaroo.git
Successfully built panaroo biocode intbitset
Successfully installed ... BioPython-1.85 ... dendropy ... edlib ...
```

`intbitset` and `edlib` build from source on arm64 without complaint.

**Verified end-to-end against real Bakta GFF3**, unmodified, for the 10-isolate
smoke cohort (Bakta v1.12.1, DB v6.0 light) — exit 0:

```
panaroo -i *.gff3 -o out -t 4 --clean-mode moderate --remove-invalid-genes
```

| | gene families |
|---|---|
| total | 10,019 |
| **core** (present in all 10 genomes) | **4,828** |
| accessory (2–9 genomes) | 5,191 |
| genome-specific | 0 |

4,828 core genes is the expected figure for *P. aeruginosa*, so the clustering
is biologically credible, not merely non-crashing.

Two operational details that matter and are easy to get wrong:

- **`--remove-invalid-genes` is required.** panaroo raises
  `ValueError: Invalid gene sequence!` on any CDS whose length is not a
  multiple of 3, is under 34 bp, or contains an internal stop. Without the flag
  it *raises*; with it, those genes are skipped. On this cohort that is 10 genes
  out of 63,359 — real Bakta data hits the same error a hand-built fixture does.
- **Bakta GFF3 needs no adapter.** panaroo requires the GFF3 to carry an
  embedded `##FASTA` section, and Bakta already emits one, so nothing has to be
  spliced in between stages 2 and 7.

Do not "correct" this back to the conda package; it cannot solve on arm64.

## 4. `pyseer` can be installed, but the unique-pattern helper cannot

This is a real limitation of the plan as written, found by inspecting the
installed package rather than assuming.

pyseer 1.1.2 ships exactly these console scripts:

```
annotate_hits_pyseer, phandango_mapper, pyseer, scree_plot_pyseer, square_mash
```

There is **no `pyseer_unique_patterns`, no `pyseer_gwas`, and no
`pyseer_combine_features`**, and the package directory contains no
`unique_patterns` module. Versions 1.2.0 through 1.4.2 all require `glmnet_py`,
which is PyPI-only and in no conda channel, so 1.1.2 is the ceiling.

Per pyseer's own documentation, the unique-pattern threshold is produced by
`scripts/count_patterns.py` **from the pyseer GitHub repository** — a script in
the source tree, never packaged into the distribution. Its whole logic is:

```
LC_ALL=C sort -u <patterns file> | wc -l      # number of unique patterns
threshold = alpha / count                      # Bonferroni
```

and its input is produced by pyseer's own `--output-patterns` flag, which 1.1.2
*does* have.

So the intended workflow is achievable, but the "helper script" must be
**vendored into this repository**, not installed.

**Resolved on ticket 16: it is vendored, at `scripts/gwas/count_patterns.py`.**
The file's body is byte-identical to upstream's, and its header records the
source URL, the retrieval date, the SHA-256 of the upstream file as retrieved
(`1081bd1c…e817643`, from `master`, 2026-10-02) **and** the SHA-256 of the
vendored logic itself (`2d1898a0…9a34b1a7`, taken from the copyright line to
the end, trailing newline normalised).
`tests/unit/test_gwas_unique_patterns.py` recomputes the second digest and
asserts it equals the one recorded in the header, so editing the logic fails CI
and editing the recorded digest without the logic fails CI too. The test suite
cannot fetch upstream, so the first digest is for a human to check the
vendoring once against the file it claims to come from.

Three things the vendoring forced into the open, none of which the upstream
file handles:

- **Its `--memory`, `--cores` and `--temp` defaults are environmental.**
  `1024`, `1` and `/tmp` are facts about a machine, not about pyseer, so all
  three are passed from the machine overlay (`runtime.unique_patterns_*` in
  `config/machines/*.yaml`) rather than inherited. They are declared on
  `PipelineConfig` as accessors that **raise** when absent, because a fallback
  would make an unconfigured run look configured. Note that `sort` is a
  separate process with its own address space, so these are deliberately *not*
  `runtime.memory_mb` / `runtime.threads`.
- **It interpolates both its file argument and its scratch directory into a
  `shell=True` command without quoting.** A space in either would change the
  command; a shell metacharacter would be executed. The vendored file must stay
  byte-identical, so the *caller* validates both resolved paths and refuses
  anything outside `[A-Za-z0-9_-./]`. Do not relax that check to work around a
  troublesome directory - move it.
- **It hands `sort` a value it computed, and a flag `sort` may not have.**
  `--memory - mem_adjust` (10) goes negative if `--memory` is under 11, which
  `sort` reports as a buffer-size complaint naming nothing anyone set; and
  `--parallel=<n>`, added above one core, is not in POSIX — `/usr/bin/sort`
  2.3-Apple on macOS accepts it, but that is a fact about that binary and not
  something to infer from the platform. The caller checks the memory floor and
  probes `sort --help` for the flag.

Also noted: `--pres` expects a roary/piggy `.Rtab` matrix, not a TSV. The
pangenome stage will have to emit an Rtab.

## 5. `samtools 0.1.19` looks wrong and is not

bioconda's newest `osx-arm64` samtools is `1.8`, so `0.1.19` reads like a
mistake. It is not. `samtools 1.8` is **unsolvable alongside `mlst`**, because
`mlst` pulls `perl-bioperl` → `perl-bio-samtools` → `samtools >=0.1.19,<0.2`.

The `0.1.19` build's **samtools** has **no `consensus` subcommand**.

An earlier draft of this document claimed this single package also provides
`bcftools`, and that both binaries report `0.1.19-96b5f2294a`. **That was
false.** `bcftools` and `htslib` are separate packages, pinned separately, and
report their own versions:

```
$ bcftools --version
bcftools 1.23.1
Using htslib 1.23.1
```

The confusion came from running `bcftools --version` before the correct
`bcftools=1.23.1` pin was applied, when the `bcftools` on `PATH` really was the
old combined build. The conclusion — call variants through bcftools, never
`samtools consensus` — is correct; the supporting explanation was not, and is
now fixed.

**Consequence for stage 6:** variant calling must go through
`bcftools mpileup` + `call`. `bcftools 1.23.1` and `htslib 1.23.1` are pinned
explicitly so the calling tools are current even though samtools is old.
`consensus` must be taken from bcftools, never from samtools.

## 6. Bakta 1.12.1 flag facts, verified by probing the real binary

The adapter probes flags at runtime rather than assuming them. Confirmed for
the installed binary:

| Flag | Accepted | Note |
|---|---|---|
| `--in` | **no** | not supported; input is positional |
| `--out` | yes | undocumented alias for `--output` |
| `--output` | yes | the documented form |
| `--db`, `--threads` | yes | |
| `--skip-crispr` | yes | required by the spec |
| `--skip-ori` | yes | required by the spec |
| `--skip-gap` | yes | the spec's flag is correct |
| `--gap` | **no** | does not exist; do not use |
| `--keep-contig-headers` | yes | |

So the adapter's behaviour — omit `--in`, pass the genome positionally, use
`--out` — is correct for Bakta 1.12.1.

## 7. Verified working binaries

Every one of these was executed and printed its version:

```
python 3.11.16            pytest 9.1.1           seqkit v2.14.0
bakta 1.12.1              mlst 2.33.1            amrfinder 4.2.7
blastn 2.17.0+            minimap2 2.31-r1302   mafft v7.526 (2024/Apr/26)
snp-sites (usage verified) bcftools 1.23.1       htslib 1.23.1
iqtree 3.1.3 (MacOS ARM 64-bit)  samtools 0.1.19-96b5f2294a
```

AMRFinderPlus's organism flag is `--organism` / `-O`, and the correct way to
list valid organisms is `--list_organisms` (with an underscore). The database
is deliberately **absent**: provisioning is out of band, and the probe confirmed
AMRFinderPlus fails cleanly when the database is missing.

## 8. The Linux environment cannot be validated from this machine

`environment/environment-linux.yml` is written but **not verified**. A
`--platform linux-64` solve from macOS fails on the `__glibc` virtual package,
which is resolved from the host:

```
bakta =1.12.1 is not installable because it requires
  diamond >=2.2.0, which requires __glibc >=2.17,<3.0.a0, which is missing on the system
```

This is not a dependency conflict; it is an inability to cross-solve. **The
Linux file must be verified on the Linux machine.**

Per-platform maxima were still established, because channel searches do honour
`--platform`:

- `panaroo` → `1.7.0` (not 1.8.0)
- `gubbins` → `3.4.1` (not 3.4.3)

This is a trap worth recording: **anaconda.org's `latest_version` field is the
maximum across all platforms, not the maximum for a given platform.** An early
sweep in this phase was wrong for exactly this reason.

### Known risk on the Linux machine

`gubbins 3.4.1` is a `py310` build, while the environment pins Python 3.11.16.
`mlst` forced `samtools` down to 0.1.19 for a comparable reason, so this may
also be unsolvable. It cannot be tested from here. If it fails, the options are a
Python 3.10 environment for the phylogeny stages, or installing Gubbins outside
conda. **Decide that on the Linux machine, not by guessing now.**

## 9. Reproducing this

```bash
micromamba env create -f environment/environment.yml      # builds
micromamba env create -f environment/environment.yml -n x --dry-run   # solve only

# is a package installable here?
micromamba create -n _probe --dry-run -y -c bioconda -c conda-forge \
  --platform osx-arm64 <package>

# which versions exist for a platform?
micromamba search -c bioconda -c conda-forge --platform osx-arm64 <package>

# what does a tool actually accept?
bakta --help ; samtools ; bcftools --version ; iqtree --version
```

`environment/versions.sh` prints what is installed for comparison against
`config/references.tsv`.
