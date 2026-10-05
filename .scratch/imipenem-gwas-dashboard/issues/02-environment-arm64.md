# 02: Build the environment on osx-arm64 and report the gaps

**What to build:** A conda environment that actually builds on this MacBook Air, and a plain statement of which required bioinformatics tools cannot be installed here, and the precise reason for each. Today the environment has never been built and several pinned packages are believed to be misnamed or mis-versioned, so nothing about the pipeline can be run here.

**Blocked by:** None (can start immediately)

**Status:** done (Phase 1). Findings in `docs/environment-arm64.md`.

**Result:** the environment builds on osx-arm64. `panaroo` is unsolvable on
arm64 at any version; `gubbins` HAS three arm64 builds but they are all
Python-3.10-only and cannot coexist with this environment's Python 3.11.16.
So stages 7 and 8 cannot run on the laptop and no full REAL run is possible
here. The Linux environment is written but **unverified** - a linux-64 solve
cannot be run from a macOS host, because glibc is resolved from the host.

**Correction (Phase 1 review):** an earlier draft claimed gubbins had "zero
osx-arm64 builds at any version". That was false. The build count is three
(3.4.3). The conclusion is unchanged; the stated reason was wrong and has been
corrected in `docs/environment-arm64.md` with the evidence.

- [ ] environment.yml resolves and the environment creates successfully via micromamba.
- [ ] The AMRFinderPlus package name, the BLAST version line and the iqtree package are each corrected against a real search, not from memory.
- [ ] gubbins and mafft are added only if they resolve; if not, they are recorded as unavailable with the evidence.
- [ ] A table lists every required tool as available or unavailable on osx-arm64, with the command that established it.
- [ ] No version pin in the file is unverified.
