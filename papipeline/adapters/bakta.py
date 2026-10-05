"""The Bakta adapter: one genome in, one validated annotation out.

This is the stage-2 adapter the pipeline declared but had not written. It
exists to make REAL-mode annotation a real external-tool invocation whose
outputs are held to a contract, rather than a ``NotImplementedError``.

It is deliberately thin. It assembles the argument vector, resolves output
paths, and hands execution to :func:`papipeline.execution.run_task`. It does
not interpret a single annotation record - that stays in
:mod:`papipeline.stages.annotation`, which owns the schema. And it does not
change Bakta's concurrency: the thread count is whatever the configuration
already said, passed straight through.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence

from ..errors import PipelineError
from ..execution.specs import BaktaOutputs, bakta_spec
from ..assemblies import locate_assembly
from ..models import Sample

if TYPE_CHECKING:  # pragma: no cover
    from ..config import PipelineConfig


@dataclass(frozen=True)
class BaktaResult:
    """What one genome's Bakta run produced, after validation."""

    sample_id: str
    state: str
    attempts: int
    outputs: BaktaOutputs
    command: Sequence[str]
    validation_detail: str = ""


def genome_stem_for(sample: Sample) -> str:
    """The file stem Bakta will use for this genome.

    Derived from the assembly filename, because that is what the file is
    actually called and a mismatch here is how outputs end up in the wrong
    place. Falls back to the sample id for an assembly-less sample.
    """
    if sample.assembly_path:
        return Path(sample.assembly_path).stem
    return sample.sample_id


def assembly_for(sample: Sample, data_root: Path) -> Path:
    """Resolve a sample's assembly. Thin alias over the shared lookup.

    The lookup itself lives in :mod:`papipeline.assemblies` so that the
    workflow, the variant caller and this adapter cannot drift into three
    different notions of where an assembly is. Kept as a named function because
    it is the adapter's own vocabulary.
    """
    return locate_assembly(sample, data_root)


def resolve_executable() -> str:
    """Locate the Bakta executable.

    Order: the ``BAKTA_BIN`` environment variable, then ``PATH``. The
    override exists because a Bakta installed inside an environment is not
    on the login ``PATH``; hardcoding a personal absolute path into the
    repository would make the project unrunnable on any other machine.
    """
    import os

    configured = os.environ.get("BAKTA_BIN")
    if configured:
        candidate = Path(configured).expanduser()
        if not candidate.is_file():
            raise PipelineError(
                "BAKTA_BIN is set but is not a file",
                bakta_bin=str(candidate),
            )
        return str(candidate)
    found = shutil.which("bakta")
    if not found:
        raise PipelineError(
            "bakta was not found. Put it on PATH, or set BAKTA_BIN to its "
            "full path. REAL-mode annotation requires the real tool; it is "
            "never substituted or simulated.",
            looked_in="PATH", env_var="BAKTA_BIN",
        )
    return found


def dependency_path_env(executable: str) -> Dict[str, str]:
    """Environment additions that let Bakta find its own dependencies.

    Bakta resolves ``tRNAscan-SE``, DIAMOND and BLAST through ``PATH``.
    When Bakta lives inside an environment that is not on ``PATH`` -
    exactly the usual case for a conda or micromamba install - those
    dependencies are invisible to it and the run aborts within a second or
    two with "dependency not found". The directory holding the executable
    is therefore prepended for the child process only. Nothing is
    installed, moved or modified; the environment is simply made visible.
    """
    import os

    bindir = str(Path(executable).resolve().parent)
    existing = os.environ.get("PATH", "")
    if bindir in existing.split(os.pathsep):
        return {}
    return {"PATH": f"{bindir}{os.pathsep}{existing}" if existing else bindir}


def accepts_flag(executable: str, flag: str) -> bool:
    """Does this Bakta actually accept ``flag``?

    The installed Bakta and the pinned version are not necessarily the same:
    ``environment.yml`` pins 1.9.3, while a real installation may be a
    later release whose command line differs. Bakta 1.10 moved the input
    genome from ``--in`` to a positional argument, so building the command
    from an assumed syntax produced "unrecognized arguments" and an
    instant, misleading failure. The real executable is therefore asked.
    """
    try:
        out = subprocess.run([executable, "--help"], capture_output=True,
                             text=True, timeout=60, check=False)
        text = (out.stdout or "") + (out.stderr or "")
    except (OSError, subprocess.SubprocessError):
        return False
    return flag in text


def bakta_command(
    executable: str,
    database: Path,
    assembly: Path,
    out_dir: Path,
    stem: str,
    threads: int,
    *,
    species: Optional[str] = None,
) -> List[str]:
    """The real Bakta invocation, in the syntax this executable accepts.

    ``--out`` is the output *directory* and Bakta names its files after the
    input stem, which is why the contract's paths are built from
    :class:`BaktaOutputs` rather than guessed here. ``--threads`` is passed
    through from configuration unchanged; this adapter does not manage
    Bakta's concurrency.
    """
    command = [executable]
    if accepts_flag(executable, "--in"):
        command += ["--in", str(assembly)]
    else:
        # Modern Bakta takes the genome positionally.
        command.append(str(assembly))
    command += ["--out", str(out_dir), "--db", str(database),
                "--threads", str(threads)]
    if accepts_flag(executable, "--force"):
        # Without it Bakta refuses to write into a non-empty directory, so a
        # re-run would fail on its own leftovers.
        command.append("--force")
    if species and accepts_flag(executable, "--species"):
        command += ["--species", species]
    return command


def run_bakta(
    sample: Sample,
    *,
    config: "PipelineConfig",
    genomes_dir: Path,
    out_root: Path,
    database: Path,
    run_key: str,
    store: Any = None,
    event_sink: Any = None,
    tool_version: Optional[str] = None,
    database_version: Optional[str] = None,
    config_hash: Optional[str] = None,
    pid_sink: Optional[Any] = None,
    species: Optional[str] = None,
) -> BaktaResult:
    """Annotate one genome with Bakta, through the real task machinery.

    Runs Bakta once. The result is ``SUCCEEDED`` only if the process exited
    cleanly *and* the outputs satisfy :func:`bakta_spec`; a zero exit with a
    wrong or truncated table is ``INVALID`` or ``INCOMPLETE``, never a
    success.

    ``pid_sink`` receives the child's process id as soon as it is spawned,
    which is the only moment a caller can attach to it. The benchmark uses
    this to sample real per-process CPU and memory.
    """
    from ..execution import RetryPolicy, TaskContext, run_task

    executable = resolve_executable()
    if not Path(database).exists():
        raise PipelineError(
            "The Bakta database directory does not exist. Databases are "
            "provisioned out of band and pinned in config/references.tsv; "
            "they are never downloaded during a run.",
            database=str(database), sample_id=sample.sample_id,
        )

    stem = genome_stem_for(sample)
    assembly = assembly_for(sample, genomes_dir)
    out_dir = Path(out_root) / sample.sample_id
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = BaktaOutputs(out_dir=out_dir, sample_id=sample.sample_id, genome_stem=stem)
    command = bakta_command(
        executable, Path(database), assembly, out_dir, stem,
        int(config.runtime.get("threads", 1)), species=species,
    )

    ctx = TaskContext(
        run_key=run_key,
        stage="annotation",
        subject=sample.sample_id,
        spec=bakta_spec(outputs),
        tool_version=tool_version,
        database_version=database_version,
        config_hash=config_hash,
        input_ids=[str(assembly)],
        log_path=Path(out_root) / "logs" / f"annotation.{sample.sample_id}.log",
    )
    result = run_task(
        ctx,
        command,
        store=store,
        policy=RetryPolicy(max_attempts=1),
        event_sink=event_sink,
        pid_sink=pid_sink,
        # Layered over the parent environment, so the executable's own
        # dependencies become visible without altering anything on disk.
        env=dependency_path_env(executable) or None,
    )
    return BaktaResult(
        sample_id=sample.sample_id,
        state=result.state.value,
        attempts=result.attempts,
        outputs=outputs,
        command=tuple(command),
        validation_detail=(result.validation.detail if result.validation else ""),
    )
