# Install

Two routes, and the choice is forced by the platform rather than by preference.

- **macOS arm64 (laptop)** — everything except stages 7–16. `gubbins` has no
  usable build here. Fine for development, TEST mode and stages 1–6.
- **Linux + conda/bioconda** — required for all 16 stages.

Full per-tool evidence, including the fault analysis behind the `gubbins`
conclusion, is in [`environment-arm64.md`](environment-arm64.md).

## macOS arm64

```bash
micromamba env create -f environment/environment.yml
micromamba activate pa-amr
```

`micromamba run -n pa-amr pytest` also works, but it **prepends the env's `bin`
to `PATH`**, which defeats any test that shadows a tool by prepending a stub
directory. Activate the environment instead when running tests that do that.

Verify what is available:

```bash
bash environment/versions.sh
```

### Limits of this platform

| Stage | Tool | State | Why |
|---|---|---|---|
| 7 pangenome | `panaroo`, `cd-hit` | **absent** | `panaroo` is pip-only on arm64 (no conda recipe); neither was in the run env |
| 8 recombination | `gubbins` | **cannot run** | see below |

`gubbins` is the hard blocker. The evidence, briefly:

- Every conda `osx-arm64` build is **Python-3.10-only**; this codebase is
  **3.11.16**, so the package will not solve into the project env.
- Installed anyway in a separate 3.10 env, it **crashes on the first run**:
  `run_gubbins.py` exits **139 (SIGSEGV)**. Input-independent — a 3-genome,
  10-SNP alignment triggers it — and it reproduces in a pristine environment.
- `nm -u` reports `_gzopen`, `_gzread`, `_gzclose` **undefined** in
  `libgubbins.0.dylib`, and `libSystem` exports no `gz*` at all
  (`nm -gU /usr/lib/libSystem.B.dylib | grep -c gopen` → `0`), so nothing in the
  process can satisfy them. Ruled out: Rosetta translation (every binary is
  native arm64) and the numba/llvmlite JIT (`NUMBA_DISABLE_JIT=1` reproduces the
  same 139).
- **No debugger backtrace is stored in this repository**, so no claim is made
  about which instruction faulted.

Because stages 9–16 consume stage 8, **a full 16-stage REAL run is not possible
on macOS arm64.** That is a property of the platform.

## Linux + conda (all 16 stages)

`gubbins` solves cleanly on Linux once Python 3.10 is available:

```bash
micromamba create -n pa-amr-linux -c conda-forge -c bioconda \
  python=3.10 gubbins=3.4.3 panaroo cd-hit   # 272 packages
```

Then create the project env as above and add the tools that are not on the base
image — `gubbins`, `panaroo`, `cd-hit`, `pyseer`, `snp-sites`, `iqtree`.

**Do not** downgrade the project env's Python to 3.10 to get `gubbins` on macOS.
That would break the laptop/big-machine comparability the config split exists to
protect. Keep `gubbins` on Linux.

## Databases

All provisioned out of band; never committed. Under `db/`:

| Database | Used by |
|---|---|
| AMRFinderPlus | stage 4 |
| VFDB | stage 5 |
| MLST scheme (*paeruginosa*) | stage 3 |
| Bakta `db-light` | stage 2, only to **verify** curated output |
| PAO1 reference | stages 6, 7 |

## Inputs

Not in this repository, and never committed:

- `PDC_essential.tsv` — the PDC essential-genes table, at the worktree root. Real
  clinical data.
- `data/` — genome assemblies.
- Bakta annotation trees, per isolate, at
  `results/<mode>/intermediate/bakta/<sample_id>`.

**Bakta is an input and is never re-run.** Stage 2 imports existing `bakta`
output under `annotation.reuse_tool_output: require` and refuses loudly when
curated output is missing, rather than executing the tool.

## Running

```bash
# TEST mode, synthetic fixtures, no external tools needed
python3 scripts/common/run_pipeline.py --mode TEST

# REAL mode: gated, and the production path is Snakemake
PIPELINE_ALLOW_REAL_MODE=1 snakemake --snakefile workflow/Snakefile \
  --config mode=REAL machine=config/machines/smoke.yaml --cores 4 full_run
```

The laptop config refuses more than 20 samples and says which limit it hit.

## Tests

```bash
python3 -m pytest tests dashboard/tests -q
```

`PDC_essential.tsv` is absent from a fresh clone, so 21 tests skip. That is by
design — supply it to exercise them.