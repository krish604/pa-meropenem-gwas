"""The source-built gubbins must not have the conda package's zlib defect.

**The defect this guards.** The bioconda `osx-arm64` build of gubbins 3.4.3
leaves `_gzopen`, `_gzread` and `_gzclose` **unresolved** in the conda process -
`nm -u` reports them undefined in `libgubbins.0.dylib`, and `libSystem` exports
no `gz*` at all (`nm -gU /usr/lib/libSystem.B.dylib | grep -c gzopen` → `0`), so
nothing in that process can satisfy them. What that costs is measured, not
inferred: on **one** identical input file with **one** identical argv and only
the binary differing, the conda build is killed by SIGSEGV and writes **zero**
files, where the source build of the same tag exits **0** and writes **nine**.

`run_gubbins.py` masks this behind "Gubbins crashed, please ensure you have
enough free memory" (`gubbins/common.py:246` swallows the real error), which is
why it reads as a memory problem and sends people looking in the wrong place.

gubbins' own `configure.ac:40` requires zlib
(`PKG_CHECK_MODULES([zlib], [zlib])`), which is consistent with a conda recipe
failing to honour a dependency the upstream build system demands. What is
**NOT** claimed here, because none of it was done: that the two builds produce
identical *output* (no diff of their results exists), that a debugger backtrace
shows the null jump (no backtrace is stored), or that **relinking** the conda
binary would help. Relinking was proposed once and withdrawn; it has never been
tried. `scripts/build_gubbins_from_source.sh` rebuilds from source, which is a
different thing from relinking.

**Why this skips.** `tools/` is gitignored - it is a build artifact, not a
source file. A fresh clone has no binary, and a test that demanded one would
fail on every clone that has not run the build script. So the absence of the
binary skips; the *presence* of the defect, when the binary is there, fails.
The config-level declarations that do not need a binary are asserted in
`tests/unit/test_config_machines.py` and `tests/unit/test_tool_discovery.py`.

Filed upstream; see the URL in `scripts/build_gubbins_from_source.sh`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
GUBBINS_BIN = PIPELINE_ROOT / "tools" / "gubbins" / "bin" / "gubbins"

pytestmark = pytest.mark.skipif(
    not GUBBINS_BIN.is_file(),
    reason=(
        "tools/gubbins/bin/gubbins is absent; it is a gitignored build artifact. "
        "Run scripts/build_gubbins_from_source.sh to produce it."
    ),
)


def _three_taxon_alignment() -> str:
    """A tiny already-aligned FASTA. gubbins wants an alignment, not sequences.

    Three taxa is the minimum its SNP-detection step accepts. It is below
    RAxML's four-species floor, but this test exercises the core binary only -
    the tree-building half of the pipeline is out of scope while stage 8 is
    unwired, and the crash being guarded happens before any tree is built.

    Lengths are equal because that is what an alignment is; the rows carry a few
    differences so the SNP step has something to find.
    """
    width = 120
    base = "ACGTACGTACGGATTACACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    rows = []
    for index in range(3):
        offset = index + 1
        rows.append(
            f">seq{index + 1}\n"
            + "".join(
                "ACGT"[("ACGT".index(char) + offset) % 4] if position % 37 == 0 else char
                for position, char in enumerate(base)
            )
        )
    assert len(rows) == 3 and width > 0
    return "\n".join(rows) + "\n"


def _run(binary: Path, alignment: Path, cwd: Path) -> subprocess.CompletedProcess:
    """Run the binary with a deliberately minimal environment.

    The point is to prove the installed binary needs nothing from the caller's
    environment: no `DYLD_LIBRARY_PATH`, no `DYLD_INSERT_LIBRARIES`, no conda
    env vars. The conda binary's whole failure mode was an unlinked import that
    only a library preload hid, so a test that ran it with the developer's shell
    already set up would not have caught anything.
    """
    return subprocess.run(
        [str(binary), str(alignment)],
        cwd=str(cwd),
        env={"PATH": "/usr/bin:/bin", "HOME": str(cwd)},
        capture_output=True,
        text=True,
        timeout=120,
    )


def test_the_built_binary_links_zlib():
    """The single assertion that would have caught the conda defect.

    Checked directly on the binary rather than inferred from a run: a link line
    that loses zlib still starts, and only fails when it reaches `gzopen`.
    """
    otool = shutil.which("otool")
    if otool is None:
        pytest.skip("otool is unavailable, so the linkage cannot be inspected")
    out = subprocess.run(
        [otool, "-L", str(GUBBINS_BIN)], capture_output=True, text=True, timeout=60
    )
    assert out.returncode == 0, out.stderr
    assert "libz" in out.stdout, (
        "the built gubbins does not link zlib. This is the exact defect filed "
        "upstream against the conda package; the binary would SIGSEGV on any "
        "existing alignment file.\n" + out.stdout
    )


def test_the_built_binary_does_not_segfault_on_an_existing_alignment(tmp_path):
    """The reproduction, run against the binary the pipeline would actually use."""
    alignment = tmp_path / "tiny.alignment"
    alignment.write_text(_three_taxon_alignment(), encoding="utf-8")

    result = _run(GUBBINS_BIN, alignment, tmp_path)

    assert result.returncode != 139, (
        "gubbins died with SIGSEGV (exit 139) - the conda package's zlib defect "
        "is present in this build"
    )
    assert result.returncode == 0, (
        f"gubbins exited {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )


def test_the_built_binary_does_not_segfault_on_an_empty_file(tmp_path):
    """The input-independent half of the signature.

    An existing but empty file crashed the conda binary too, because the crash
    was on the path taken once the file-existence check passes and before any
    parsing. If this ever crashes while the non-empty case passes, something
    other than the zlib linkage has regressed.
    """
    empty = tmp_path / "empty.alignment"
    empty.write_text("", encoding="utf-8")

    result = _run(GUBBINS_BIN, empty, tmp_path)

    assert result.returncode != 139, (
        "gubbins died with SIGSEGV on an empty existing file - input-independent "
        "crash, same signature as the conda defect"
    )


def test_a_missing_file_is_reported_rather_than_crashing(tmp_path):
    """The control case: a nonexistent file must exit non-zero but not crash.

    gubbins validates the path before opening anything, so this is the one input
    the broken conda binary handled correctly. Keeping it under test is what
    makes the two assertions above meaningful rather than incidental.
    """
    absent = tmp_path / "not_here.alignment"
    result = _run(GUBBINS_BIN, absent, tmp_path)
    assert result.returncode == 1
    assert result.returncode not in (139, -11)


def test_the_built_binary_writes_the_outputs_the_next_stage_reads(tmp_path):
    """Non-crashing is not the same as working.

    A run that segfaults and a run that silently produces nothing would both
    satisfy a bare exit-code check on some failure paths, so assert the actual
    artifacts gubbins exists to produce.
    """
    alignment = tmp_path / "tiny.alignment"
    alignment.write_text(_three_taxon_alignment(), encoding="utf-8")

    result = _run(GUBBINS_BIN, alignment, tmp_path)
    assert result.returncode == 0, result.stderr

    produced = {path.name for path in tmp_path.iterdir()}
    for expected in (
        "tiny.alignment.vcf",
        "tiny.alignment.phylip",
        "tiny.alignment.snp_sites.aln",
    ):
        assert expected in produced, (
            f"gubbins exited 0 but did not write {expected}; wrote: {sorted(produced)}"
        )
    vcf = tmp_path / "tiny.alignment.vcf"
    assert vcf.stat().st_size > 0, "the VCF gubbins wrote is empty"


def test_the_installed_library_ships_alongside_the_binary():
    """libgubbins records an absolute install path, so it must be reachable.

    This is the trade-off of the source build: the binary carries an absolute
    `install_name` for `libgubbins.0.dylib`, so it only loads from the prefix it
    was built with. If the library moves, the failure is a dyld error rather than
    a crash - still worth catching here rather than three stages into a run.

    Asserted by *existence*, not by string comparison against the repository
    root. Comparing paths would fail in any git worktree or clone that has its
    own `tools/` while the binary was built in another, which is a supported
    layout rather than a defect.
    """
    otool = shutil.which("otool")
    if otool is None:
        pytest.skip("otool is unavailable")
    out = subprocess.run(
        [otool, "-L", str(GUBBINS_BIN)], capture_output=True, text=True, timeout=60
    )
    assert out.returncode == 0, out.stderr

    recorded = [
        line.strip().split(" (")[0]
        for line in out.stdout.splitlines()[1:]
        if "libgubbins" in line
    ]
    assert recorded, (
        "otool -L reports no libgubbins dependency at all\n" + out.stdout
    )
    for path in recorded:
        assert Path(path).is_file(), (
            f"the binary's recorded libgubbins path does not exist: {path}\n"
            "The build is stale - re-run scripts/build_gubbins_from_source.sh"
        )


def test_the_build_script_is_not_a_committed_artifact_itself():
    """`tools/` is gitignored; the script that fills it is not.

    A regression here means someone committed a compiled binary, which hard rule
    4 forbids in spirit even though the file is not a genome.
    """
    ignore = (PIPELINE_ROOT / ".gitignore").read_text(encoding="utf-8")
    assert os.path.join("", "tools") + "/" in ignore or "\ntools/\n" in ignore, (
        "tools/ is not gitignored; the built binary would be committed"
    )
    assert (PIPELINE_ROOT / "scripts" / "build_gubbins_from_source.sh").is_file()