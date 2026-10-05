"""Stage 8 REAL caller: the gubbins binary, and the formats it emits.

Separated from ``papipeline.stages.recombination`` for the reason every other
adapter in this package is separated from its stage: the stage owns the
contract, this owns the process and the file formats a third party happens to
produce.

**The binary is the load-bearing decision here, and it is not a config value.**
The bioconda ``osx-arm64`` build of gubbins 3.4.3 ships a
``libgubbins.0.dylib`` that references ``_gzopen``/``_gzread``/``_gzclose`` and
cannot resolve them, so it dies with SIGSEGV — exit 139 — on any real input.
The same sources build and link correctly from source, which is what
``scripts/build_gubbins_from_source.sh`` does and what ``laptop.yaml`` declares
under ``tool_search_dirs``. Upstream: nickjcroucher/gubbins#452.

The two builds self-report **different versions — 3.4.2 (source) and 3.4.3
(conda)** — and that string is the cheap discriminator. Note that
``run_gubbins.py``'s own banner prints ``--- Gubbins 3.4.3 ---`` on *both*
builds, because that is the Python package's version and not the C binary's.
Reading the banner to decide which binary ran is a mistake this module makes
impossible by never consulting it.

**Why an isolated PATH rather than a path argument.** ``run_gubbins.py`` has no
flag that names a binary. ``gubbins/common.py`` sets ``gubbins_exec = 'gubbins'``
and resolves it with ``utils.which``, a plain left-to-right scan of ``PATH``,
and its only fallback *appends* ``/usr/lib/gubbins/`` — which can never beat an
entry that is already earlier. So the reliable override is to prepend a
directory holding exactly one entry: a symlink to the resolved source binary.
The environment's own ``bin`` stays on PATH *after* it, because gubbins needs
``raxml-ng`` and ``pyjar`` from there and the point is to shadow one executable,
not to empty PATH.

**Why the self-test runs before the real run.** A segfault discovered 30 seconds
into a real run has already written a working directory full of intermediates
and looks like a gubbins bug rather than a packaging defect. :func:`self_test`
runs the resolved binary against a tiny generated alignment and reports the
defect by name if it segfaults.

**Why the process runs in a scratch directory.** ``run_gubbins.py`` has no
output-directory flag; ``--prefix`` only renames (``common.py:1543-1551`` maps
input names to output names by string substitution) while everything is written
relative to the process CWD. Running it in the repository would scatter ~40
intermediate files into a tracked tree.

**One environment quirk is reproduced deliberately.** gubbins' RAxML-NG builder
probes CPU features by shelling out to ``sysctl``; when ``sysctl`` is absent the
probe leaves a local variable unassigned and the next statement reads it
(``UnboundLocalError``). So the PATH handed to the subprocess must retain
``/usr/sbin``. That is a gubbins bug, and the workaround is "keep a normal
system PATH", not "make PATH more minimal".
"""

from __future__ import annotations

import json
import os
import random
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from ..errors import (
    DataContractError,
    ModeNotAllowedError,
    ToolExecutionError,
    ToolNotAvailableError,
)
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.gubbins")

#: gubbins' own names for the outputs this stage reads. Taken from
#: ``translation_of_filenames_to_final_filenames`` (``common.py:1543-1551``), not
#: guessed: the prefix is interpolated, so a wrong constant here is a stage that
#: reports "no recombination found" because it looked for a file under a
#: different name.
SNP_ALIGNMENT_SUFFIX = ".filtered_polymorphic_sites.fasta"
FINAL_TREE_SUFFIX = ".final_tree.tre"
NODE_LABELLED_TREE_SUFFIX = ".node_labelled.final_tree.tre"
PER_BRANCH_STATISTICS_SUFFIX = ".per_branch_statistics.csv"
LOG_SUFFIX = ".log"

#: gubbins' two terminal iteration-loop messages, read out of `common.py`.
#:
#: `common.py:585-586` prints the convergence line and `common.py:590` prints the
#: ceiling line; one or the other is all the loop ever says about stopping, and
#: both go to stdout, which the adapter already captures.
#:
#: **Why this has to be parsed at all.** gubbins exits 0 either way: a run that
#: hits its iteration ceiling and a run that converges produce the same exit
#: code and the same final tree. Exit status therefore cannot distinguish a
#: result from a truncated one, and the tree, the masked alignment and
#: `recombination.tsv` are all written either way — so without this the only way
#: to tell a converged run from a non-converged one is to read gubbins' console
#: output by hand.
#:
#: Anchored at the start of a line, because gubbins echoes its own commands and
#: a substring match would fire on a line that merely quotes the message.
CONVERGED_MESSAGE_RE = re.compile(
    r"^Convergence after (\d+) iterations:", re.MULTILINE
)
ITERATION_CEILING_MESSAGE_RE = re.compile(
    r"^Maximum number of iterations \((\d+)\) reached\.", re.MULTILINE
)

#: The `--prefix` this stage passes. Also the stem of every output name above.
PREFIX = "recombination"

#: The input alignment, named by panaroo.
#:
#: Declared here rather than at each call site so the workflow, the adapter and
#: the stage cannot disagree about the file's name: the Snakefile needs it to
#: declare a DAG edge, and a literal in the Snakefile would be a second
#: definition of a location the code already owns.
#:
#: The **filtered** alignment, not ``core_gene_alignment.aln``. Panaroo filters
#: alignment *quality*; gubbins filters *recombination*. Two filters, two
#: purposes, applied in series. See ``docs/design/recombination-contract.md``
#: section 2.1.
CORE_ALIGNMENT_BASENAME = "core_gene_alignment_filtered.aln"

#: Where this stage's gubbins working directory lives, relative to the run's
#: intermediate root. Layout, not environment: the root itself comes from
#: configuration, exactly as `adapters.panaroo.OUTPUT_DIRNAME` works.
WORK_DIRNAME = "gubbins"

#: ``run_gubbins.py``'s literal spelling of the tree builder. No hyphen — this
#: was read out of ``run_gubbins.py --help`` rather than assumed, and the
#: hyphenated form is not accepted.
TREE_BUILDER = "raxmlng"

#: Exit status of a process killed by SIGSEGV. ``subprocess`` reports a signal
#: death as ``128 + signum``, so SIGSEGV(11) surfaces as 139. The conda
#: ``osx-arm64`` gubbins fails this way and no other way worth conflating.
SEGFAULT_EXIT = 139

#: Version string reported by the source build; the conda ``osx-arm64`` build
#: reports :data:`CONDA_BUILD_VERSION`. Used only for diagnostics and tests —
#: the refusal decision is made by the self-test's exit code, because a version
#: string can be spoofed by a wrapper script and a segfault cannot.
SOURCE_BUILD_VERSION = "3.4.2"
CONDA_BUILD_VERSION = "3.4.3"

#: Directories kept on PATH *after* the isolated binary directory. ``sysctl`` is
#: not optional: without it gubbins' RAxML-NG builder raises UnboundLocalError
#: before it ever runs. See the module docstring.
PATH_TAIL_DIRS: Tuple[str, ...] = ("/usr/sbin", "/usr/bin", "/bin", "/usr/local/bin", "/sbin")

#: Shape of the probe alignment used by :func:`self_test`. Small enough to be
#: instant, large enough that the binary must do real work before it can crash.
PROBE_TAXA = 5
PROBE_SITES = 420  # divisible by 60, so every line is the same width
PROBE_SNP_RATE = 0.12

#: Per-branch statistics columns this module reads. Located **by name**; the
#: file has 13 columns and gubbins is free to reorder them.
NODE_COLUMN = "Node"
TOTAL_SNPS_COLUMN = "Total SNPs"
BLOCKS_COLUMN = "Number of Recombination Blocks"


@dataclass(frozen=True)
class Convergence:
    """What gubbins' iteration loop actually did.

    `converged` is tri-state on purpose. ``False`` is a *finding* — the loop ran
    and stopped at its ceiling, so the reported tree and recombination blocks
    are from an unfinished iteration. ``None`` is an *absence of a finding* —
    the log said neither, which usually means the run died before the loop
    reported. Collapsing the two into one boolean is what would let a crashed
    run read as a converged one.

    `iterations` is the iteration the loop stopped at; `iterations_max` is the
    ceiling it was given, which is the same number gubbins prints in the
    non-converged message. Both are ``None`` when the log says neither.
    """

    converged: Optional[bool] = None
    iterations: Optional[int] = None
    iterations_max: Optional[int] = None
    note: str = ""

    @property
    def recordable(self) -> bool:
        """Whether this says anything a downstream consumer can act on."""
        return self.converged is not None


#: Neither converged nor not. Distinct from ``Convergence(converged=False)``.
UNKNOWN_CONVERGENCE = Convergence(
    note=(
        "gubbins' log carried neither a convergence line nor an iteration-ceiling "
        "line, so this run's convergence is unknown rather than false. Treat "
        "the outputs as unverified."
    ),
)


def parse_convergence(stdout: str) -> Convergence:
    """Read gubbins' convergence verdict out of its own console output.

    Both terminal messages are read, and the ceiling message wins if both are
    present: a run that converged early does not also print the ceiling line,
    but a resumed run (``--resume``) can carry lines from an earlier attempt in
    the same stdout, and in that case the ceiling is the later, weaker claim.
    Preferring ``False`` makes the conservative outcome the one that survives.

    Args:
        stdout: The captured stdout+stderr of ``run_gubbins.py``.

    Returns:
        A `Convergence`. Never raises: an unparseable log is
        `UNKNOWN_CONVERGENCE`, because a provenance record is written on failure
        too and must not itself be a second failure.
    """
    text = stdout or ""
    ceiling = ITERATION_CEILING_MESSAGE_RE.search(text)
    converged = CONVERGED_MESSAGE_RE.search(text)

    if ceiling is not None:
        return Convergence(
            converged=False,
            iterations=int(ceiling.group(1)),
            iterations_max=int(ceiling.group(1)),
            note=(
                f"gubbins stopped at its iteration ceiling ({ceiling.group(1)}) "
                "without converging. The final tree, the masked alignment and "
                "recombination.tsv are all written by this run and all come from "
                "an unfinished iteration loop; a downstream consumer must be able "
                "to see that."
            ),
        )
    if converged is not None:
        return Convergence(
            converged=True,
            iterations=int(converged.group(1)),
            note=f"gubbins reported convergence after {converged.group(1)} iterations.",
        )
    return UNKNOWN_CONVERGENCE


@dataclass(frozen=True)
class GubbinsResult:
    """What one gubbins run produced, and how it was invoked.

    `binary_version` is the string the **binary** reported, read from the binary
    itself, not the ``run_gubbins.py`` banner — see the module docstring for why
    that distinction is load-bearing.
    """

    work_dir: Path
    prefix: str
    snp_alignment: Path
    final_tree: Path
    node_labelled_tree: Path
    per_branch_statistics: Path
    log: Path
    binary: Path
    binary_version: str
    isolated_bin_dir: Path
    argv: List[str] = field(default_factory=list)
    exit_code: int = 0
    seconds: float = 0.0
    stdout: str = ""
    #: Whether gubbins' iteration loop converged. ``None`` when the log says
    #: neither way, which is not the same as ``False``.
    #:
    #: A field rather than a derived property because a `GubbinsResult` is also
    #: built by hand — by the tests, and by any caller holding gubbins' output
    #: from an earlier run — and a property over `stdout` would then report on a
    #: string that was never captured. `run_gubbins` is the only production
    #: writer, and it sets this from the stdout it captured.
    convergence: "Convergence" = field(default_factory=lambda: UNKNOWN_CONVERGENCE)

    def provenance(self) -> Dict[str, object]:
        """The record written to ``provenance.json``.

        Includes the exit code and the elapsed time even when they say the run
        failed, because provenance that only exists on success cannot answer the
        question worth asking, which is what ran when something went wrong.
        """
        return {
            "binary": str(self.binary),
            "binary_version": self.binary_version,
            "converged": self.convergence.converged,
            "iterations": self.convergence.iterations,
            "iterations_max": self.convergence.iterations_max,
            "convergence_note": self.convergence.note,
            "isolated_bin_dir": str(self.isolated_bin_dir),
            "isolated_path_prefix": str(self.isolated_bin_dir),
            "argv": list(self.argv),
            "exit_code": self.exit_code,
            "seconds": round(self.seconds, 3),
            "prefix": self.prefix,
            "work_dir": str(self.work_dir),
            "tree_builder": TREE_BUILDER,
            "snp_alignment": str(self.snp_alignment),
            "final_tree": str(self.final_tree),
            "node_labelled_tree": str(self.node_labelled_tree),
            "per_branch_statistics": str(self.per_branch_statistics),
        }


def resolve_binary(search_dirs: Sequence[str]) -> Path:
    """Find the gubbins binary to use, by scanning `search_dirs` in order.

    `search_dirs` is ``MachineConfig.resolved_tool_search_dirs()``, so on the
    laptop it contains ``<repo>/tools/gubbins/bin`` — where
    ``scripts/build_gubbins_from_source.sh`` puts a working build — and on
    Linux it is empty, leaving ``PATH`` to decide.

    **What this does not do:** it does not verify the binary works. It can only
    report what exists; :func:`self_test` is what decides. A search that returns
    a path whose contents nobody has executed is a fact about the filesystem,
    not about the tool.
    """
    for directory in search_dirs:
        candidate = Path(directory) / "gubbins"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    found = shutil.which("gubbins")
    if found:
        return Path(found)
    raise ToolNotAvailableError(
        "Stage 8 (recombination) cannot find a gubbins binary. It looked in "
        f"{list(search_dirs)!r} and then on PATH. On this platform the conda "
        "osx-arm64 package is broken (see the message from the self-test), so "
        "the binary has to be built from source: run "
        "scripts/build_gubbins_from_source.sh, which writes it to "
        "tools/gubbins/bin/gubbins and which the laptop overlay already declares "
        "in tool_search_dirs."
    )


def binary_version(binary: Path) -> str:
    """The version string the **binary** reports about itself.

    Invoked with no arguments, which makes gubbins print its usage block and
    exit 1 — that is expected here and is not an error. The version is on the
    `Version:` line of that block.

    A non-zero exit is tolerated on purpose. What this function is for is the
    string; refusing to return one because the tool exited oddly would turn a
    diagnostic into a second failure.
    """
    completed = subprocess.run(
        [str(binary)], capture_output=True, text=True, check=False, timeout=120
    )
    text = (completed.stdout or "") + (completed.stderr or "")
    for line in text.splitlines():
        if line.strip().lower().startswith("version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def write_probe_alignment(path: Path, seed: int = 20240617) -> Path:
    """A tiny already-aligned FASTA for the self-test to run against.

    Generated rather than committed, so the self-test needs no fixture on disk —
    which is what lets it run on a machine where the gubbins toolchain is
    installed but nothing else is.

    Gaps are deliberately absent: a probe containing gaps exercises the
    alignment-filtering path rather than the one being tested, and a crash there
    would be reported as the packaging defect it is not.
    """
    rng = random.Random(seed)
    base = [rng.choice("ACGT") for _ in range(PROBE_SITES)]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for index in range(PROBE_TAXA):
            seq = list(base)
            for _ in range(int(PROBE_SITES * PROBE_SNP_RATE)):
                position = rng.randrange(PROBE_SITES)
                seq[position] = rng.choice([c for c in "ACGT" if c != seq[position]])
            handle.write(f">probe_{index + 1}\n")
            handle.write("\n".join("".join(seq[i:i + 60]) for i in range(0, PROBE_SITES, 60)))
            handle.write("\n")
    return path


def isolated_env(
    binary: Path,
    search_dirs: Sequence[str] = (),
    bin_dir: Optional[Path] = None,
) -> Tuple[Dict[str, str], Path]:
    """A PATH in which `binary` shadows every other gubbins.

    Args:
        binary: The binary that must win.
        search_dirs: Appended to the PATH *after* the system directories, so a
            machine's declared tool directories stay reachable.
        bin_dir: Where to put the shadowing directory. **Defaults to a sibling
            of `binary`**, which is the source-build tree on the laptop — and
            that would mean writing a symlink into a repository checkout on
            every self-test. Callers running for real must therefore pass a
            directory under the run's own scratch tree.

    Returns the environment to hand to the subprocess and the directory that was
    prepended, so provenance can record *how* the shadowing was done rather
    than only which file won.

    The directory holds exactly one entry. Two would defeat the point: the
    purpose is to control one executable, and a directory of convenience
    symlinks is a second place for a tool to come from unnoticed.
    """
    bin_dir = Path(bin_dir) if bin_dir is not None else Path(binary).parent / "isolated_bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    link = bin_dir / "gubbins"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(Path(binary).resolve())

    tail = [str(Path(d)) for d in PATH_TAIL_DIRS if Path(d).is_dir()]
    for directory in search_dirs:
        candidate = str(Path(directory))
        if candidate not in tail:
            tail.append(candidate)
    path = os.pathsep.join([str(bin_dir), *tail])

    environment = dict(os.environ)
    environment["PATH"] = path
    return environment, bin_dir


def self_test(
    binary: Path, scratch: Path, search_dirs: Sequence[str] = ()
) -> Tuple[int, str]:
    """Run `binary` against a tiny probe alignment and return (exit code, version).

    A segfault is caught here, before a real run has written anything, and is
    reported by name — see :func:`refuse_segfaulting_binary`.

    A non-zero exit that is **not** 139 is returned as-is rather than folded into
    the segfault verdict. A usage error, a missing tree builder and a broken
    dylib are three different problems and telling an operator their tool is
    mis-packaged when it is merely mis-invoked sends them to rebuild a binary
    that is fine.
    """
    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    probe = write_probe_alignment(scratch / "self_test_alignment.aln")
    environment, _ = isolated_env(
        binary, search_dirs, bin_dir=scratch / "isolated_bin"
    )
    completed = subprocess.run(
        [str(binary), str(probe)],
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
        cwd=str(scratch),
        env=environment,
    )
    return completed.returncode, binary_version(binary)


def refuse_segfaulting_binary(
    binary: Path, version: str, exit_code: int, scratch: Path
) -> None:
    """Raise the exit-139 refusal, naming the conda osx-arm64 build and 139.

    Named separately from the probe run so the message has one home and cannot
    drift from the condition that triggers it.
    """
    raise ToolNotAvailableError(
        f"the gubbins binary at {binary} exits {SEGFAULT_EXIT} (SIGSEGV) on a "
        f"{PROBE_TAXA}-taxon {PROBE_SITES}-site probe alignment it was given "
        f"just now (self-test directory: {scratch}). That is the bioconda "
        f"osx-arm64 build of gubbins {CONDA_BUILD_VERSION}, whose "
        f"libgubbins.0.dylib references _gzopen/_gzread/_gzclose and cannot "
        f"resolve them at load time - it segfaults before it reads a single "
        f"base. This report was generated because the binary self-reports "
        f"version {version}; the source build reports {SOURCE_BUILD_VERSION}. "
        f"Do not proceed on this binary: a real run would fail the same way "
        f"after writing a working directory full of intermediates. Build gubbins "
        f"from source instead (scripts/build_gubbins_from_source.sh, which "
        f"writes to tools/gubbins/bin/gubbins and which the laptop overlay "
        f"already declares in tool_search_dirs), or point the adapter at an "
        f"explicit path. Filed upstream: "
        f"https://github.com/nickjcroucher/gubbins/issues/452",
        binary=str(binary),
        binary_version=version,
        exit_code=exit_code,
    )


def resolve_runner(search_dirs: Sequence[str] = ()) -> Path:
    """Find ``run_gubbins.py``: the Python driver that wraps the binary.

    Needed because the driver is a *separate artefact* from the binary. The
    binary comes from a source build under ``tools/gubbins/bin``; the driver
    comes from the gubbins conda environment, because it is the ``gubbins``
    Python package plus ``pyjar`` and RAxML-NG. So there is no single directory
    holding both, and a resolver that only looked beside the binary would report
    "absent" on a machine where the whole toolchain is installed.

    Overridable with ``PAPIPELINE_GUBBINS_RUNNER``, following the repository's
    existing opt-in convention for environment-provided paths
    (``PAPIPELINE_REAL_DATA``). Overridable because hard-coding a micromamba
    prefix would encode this machine's layout into the code.
    """
    override = os.environ.get("PAPIPELINE_GUBBINS_RUNNER")
    if override:
        candidate = Path(override)
        if candidate.is_file():
            return candidate

    for directory in search_dirs:
        candidate = Path(directory) / "run_gubbins.py"
        if candidate.is_file():
            return candidate

    found = shutil.which("run_gubbins.py")
    if found:
        return Path(found)
    raise ToolNotAvailableError(
        "Stage 8 (recombination) cannot find run_gubbins.py. It looked in "
        f"{list(search_dirs)!r} and then on PATH. The driver is part of the "
        "gubbins conda environment and is separate from the binary, which comes "
        "from a source build under tools/gubbins/bin — so the two are in "
        "different places by design. Point at it with PAPIPELINE_GUBBINS_RUNNER, "
        "or install the gubbins environment so the driver is on PATH."
    )


def build_command(
    runner: Path, alignment: Path, work_dir: Path, threads: int
) -> List[str]:
    """The ``run_gubbins.py`` invocation.

    Every scientific choice arrives as an argument; nothing here restates it.
    The alignment is passed **by absolute path** because gubbins writes relative
    to its CWD and is given a different CWD than the caller.
    """
    return [
        str(runner),
        "--prefix", PREFIX,
        "--tree-builder", TREE_BUILDER,
        "--threads", str(int(threads)),
        str(Path(alignment).resolve()),
    ]


def require_real_gate(config) -> None:
    """Refuse a REAL run while ``runtime.allow_real_mode`` is false.

    **Checked before anything else**, because an operator whose gate is shut does
    not need to be told their alignment is missing — that is not their problem,
    and fixing it would not unblock them. The wording follows the house style at
    ``stages/similarity.py:377``: the flag, the one-session environment
    override, and the overlay that currently holds it shut.

    A *gate on the run*, not a refusal of the stage: open it and the stage
    proceeds. That distinction is load-bearing elsewhere — the taxonomy derives
    "refuses REAL" from the AST and treats a gate as the opposite.
    """
    if bool(config.runtime.get("allow_real_mode", False)):
        return
    raise ModeNotAllowedError(
        f"REAL-mode stage 8 (recombination) is gated: runtime.allow_real_mode "
        f"is false in the machine overlay for {config.machine_name!r}. gubbins "
        "is provisioned out-of-band and reads its core gene alignment from this "
        "run's intermediate directory; neither happens unless REAL mode is "
        "authorised. Set it in the overlay, or open it for one session with "
        "PIPELINE_ALLOW_REAL_MODE=1. The committed overlays keep it shut so a "
        "REAL run cannot begin by accident.",
        machine=config.machine_name,
    )


def run_gubbins(
    alignment: Path,
    work_dir: Path,
    *,
    runner: Path,
    binary: Path,
    threads: int,
    search_dirs: Sequence[str] = (),
    config=None,
    timeout: int = 86400,
    selftest: Optional[Callable[[], Tuple[int, str]]] = None,
) -> GubbinsResult:
    """Run gubbins on `alignment` inside `work_dir` and return what it wrote.

    `work_dir` is created if absent and is used as the process CWD. It must not
    be the repository: gubbins writes roughly forty files relative to its CWD.

    Args:
        alignment: The core gene alignment, named by panaroo.
        work_dir: Per-run scratch directory, under the run's intermediate root.
        runner: ``run_gubbins.py``.
        binary: The gubbins binary to force, by PATH shadowing.
        threads: From configuration. Never a literal.
        search_dirs: Passed through to the PATH construction.
        config: When given, the REAL gate is checked before anything else.
        selftest: Overrides :func:`self_test`. Exists so the refusal paths are
            testable on a machine that has only a *working* binary — there is no
            way to make a good binary segfault on demand, and a refusal that can
            only be tested by shipping a broken binary is a refusal nobody
            tests. Defaults to the real self-test.

    Raises:
        ModeNotAllowedError: ``config`` was given and ``allow_real_mode`` is false.
        DataContractError: the alignment is absent.
        ToolNotAvailableError: the binary segfaulted its self-test.
        ToolExecutionError: the self-test failed some other way, gubbins exited
            non-zero, or it exited zero without writing the SNP alignment.
    """
    if config is not None:
        require_real_gate(config)

    alignment = Path(alignment)
    if not alignment.is_file():
        raise DataContractError(
            f"no core gene alignment at {alignment}. It is the input this stage "
            "reads, written by stage 7 (panaroo) as core_gene_alignment_filtered."
            "aln. A missing alignment is a missing input, not an absence of "
            "recombination: this stage has not looked for anything.",
            alignment=str(alignment),
        )

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    probe = selftest or (lambda: self_test(binary, work_dir / "self_test", search_dirs))
    exit_code, version = probe()
    if exit_code == SEGFAULT_EXIT:
        refuse_segfaulting_binary(binary, version, exit_code, work_dir / "self_test")
    if exit_code != 0:
        raise ToolExecutionError(
            f"the gubbins self-test exited {exit_code} on a {PROBE_TAXA}-taxon "
            f"{PROBE_SITES}-site probe alignment (binary {binary}, version "
            f"{version}). That is not the exit-{SEGFAULT_EXIT} segfault of the "
            "broken conda osx-arm64 build, so this is a different fault - a "
            "missing tree builder, an unreadable probe, or a bad install - and "
            "the segfault remedy will not fix it. Run the binary by hand on the "
            "probe alignment in the self-test directory to see its own message.",
            binary=str(binary),
            binary_version=version,
            exit_code=exit_code,
        )

    # The isolated directory shadows ONE executable; everything else gubbins
    # shells out to must stay reachable. Those tools - `raxml-ng`, `iqtree`,
    # `pyjar` - all live in the same `bin` as `run_gubbins.py`, so the runner's
    # own directory is put on the PATH tail. Omitting it does not fail loudly:
    # gubbins gets as far as "No usable version of IQTree could be found." only
    # after building a first tree, by which point the run has already spent its
    # setup. It was caught by running the real tool, not by any unit test.
    environment, bin_dir = isolated_env(
        binary,
        [*search_dirs, str(Path(runner).parent)],
        bin_dir=work_dir / "isolated_bin",
    )
    command = build_command(runner, alignment, work_dir, threads)

    LOGGER.info("Stage 8: %s on %s", TREE_BUILDER, alignment.name)
    started = time.monotonic()
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
            cwd=str(work_dir),
            env=environment,
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolExecutionError(
            "gubbins timed out",
            alignment=str(alignment), timeout_s=timeout,
        ) from exc
    except OSError as exc:
        raise ToolExecutionError(
            "could not launch run_gubbins.py",
            command=" ".join(command), error=str(exc),
        ) from exc
    elapsed = time.monotonic() - started

    stdout = (completed.stdout or "") + (completed.stderr or "")
    (work_dir / f"{PREFIX}{LOG_SUFFIX}.stdout").write_text(stdout, encoding="utf-8")

    snp_alignment = work_dir / f"{PREFIX}{SNP_ALIGNMENT_SUFFIX}"
    result = GubbinsResult(
        work_dir=work_dir,
        prefix=PREFIX,
        snp_alignment=snp_alignment,
        final_tree=work_dir / f"{PREFIX}{FINAL_TREE_SUFFIX}",
        node_labelled_tree=work_dir / f"{PREFIX}{NODE_LABELLED_TREE_SUFFIX}",
        per_branch_statistics=work_dir / f"{PREFIX}{PER_BRANCH_STATISTICS_SUFFIX}",
        log=work_dir / f"{PREFIX}{LOG_SUFFIX}",
        binary=Path(binary),
        binary_version=version,
        isolated_bin_dir=bin_dir,
        argv=command,
        exit_code=completed.returncode,
        seconds=elapsed,
        stdout=stdout,
        convergence=parse_convergence(stdout),
    )

    if completed.returncode != 0:
        raise ToolExecutionError(
            f"gubbins exited {completed.returncode} for {alignment.name}",
            alignment=str(alignment),
            command=" ".join(command),
            work_dir=str(work_dir),
            stderr=stdout[-2000:],
        )
    if not snp_alignment.is_file():
        # Exit 0 with no SNP alignment is not success. A stage that read that as
        # "no recombination" would report an absence of signal where there was a
        # tool failure - the exact substitution this stage exists to avoid.
        raise ToolExecutionError(
            f"gubbins exited 0 but wrote no {snp_alignment.name} in {work_dir}, "
            "so the recombination analysis did not complete",
            alignment=str(alignment), command=" ".join(command),
        )
    return result


def write_provenance(result: GubbinsResult, path: Path) -> Path:
    """Write ``result.provenance()`` to `path`.

    Written from whatever the caller has, so a caller that wants to record a
    failed run can. The stage calls this on both paths.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(result.provenance(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def parse_per_branch_statistics(path: Path) -> List[Dict[str, object]]:
    """Read ``per_branch_statistics.csv`` into dicts keyed by column name.

    Columns located by name, not position: the file has thirteen and gubbins
    orders them by its own logic. Reading positionally would silently swap
    ``Total SNPs`` for ``Number of SNPs Inside Recombinations`` — which is a
    number that looks right and means something else entirely.

    One row per node, in file order. Tips appear alongside internal nodes under
    their own names; the two are not separated here because the tree join is
    what distinguishes them, and that happens in the stage.
    """
    path = Path(path)
    if not path.is_file():
        raise DataContractError(
            f"gubbins did not write {path.name} in {path.parent}. This file is "
            "written on every successful run; its absence means the tool failed, "
            "not that no recombination was found."
        )

    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if not lines:
        raise DataContractError(
            f"{path.name} is empty; gubbins wrote no per-branch statistics",
            path=str(path),
        )

    header = lines[0].split("\t")
    missing = [
        column
        for column in (NODE_COLUMN, TOTAL_SNPS_COLUMN, BLOCKS_COLUMN)
        if column not in header
    ]
    if missing:
        raise DataContractError(
            f"{path.name} is missing the {missing!r} column(s); header was "
            f"{header!r}. Reading it positionally instead would be worse than "
            "refusing: these columns have similar names and different meanings.",
            path=str(path),
        )

    node_i, snps_i, blocks_i = (
        header.index(NODE_COLUMN),
        header.index(TOTAL_SNPS_COLUMN),
        header.index(BLOCKS_COLUMN),
    )

    rows: List[Dict[str, object]] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) <= max(node_i, snps_i, blocks_i):
            raise DataContractError(
                f"{path.name} row has {len(fields)} fields, fewer than the "
                f"{len(header)} its own header declares. A short row means the "
                "file is truncated, not that the node had no SNPs.",
                path=str(path),
            )
        node = fields[node_i].strip()
        if not node:
            continue
        try:
            n_snps = int(fields[snps_i].strip() or 0)
            n_blocks = int(fields[blocks_i].strip() or 0)
        except ValueError as exc:
            raise DataContractError(
                f"{path.name} node {node!r} has a non-integer SNP or block count "
                f"({fields[snps_i]!r}, {fields[blocks_i]!r})",
                path=str(path),
            ) from exc
        rows.append(
            {"node": node, "n_snps": n_snps, "recombination_detected": int(n_blocks > 0)}
        )
    return rows


__all__ = [
    "BLOCKS_COLUMN",
    "CONDA_BUILD_VERSION",
    "CORE_ALIGNMENT_BASENAME",
    "Convergence",
    "CONVERGED_MESSAGE_RE",
    "FINAL_TREE_SUFFIX",
    "GubbinsResult",
    "ITERATION_CEILING_MESSAGE_RE",
    "NODE_COLUMN",
    "NODE_LABELLED_TREE_SUFFIX",
    "PREFIX",
    "PER_BRANCH_STATISTICS_SUFFIX",
    "PROBE_SITES",
    "PROBE_TAXA",
    "SEGFAULT_EXIT",
    "SNP_ALIGNMENT_SUFFIX",
    "SOURCE_BUILD_VERSION",
    "TOTAL_SNPS_COLUMN",
    "TREE_BUILDER",
    "UNKNOWN_CONVERGENCE",
    "WORK_DIRNAME",
    "binary_version",
    "build_command",
    "isolated_env",
    "parse_convergence",
    "parse_per_branch_statistics",
    "refuse_segfaulting_binary",
    "require_real_gate",
    "resolve_binary",
    "resolve_runner",
    "run_gubbins",
    "self_test",
    "write_probe_alignment",
    "write_provenance",
]