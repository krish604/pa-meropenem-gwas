"""Stage 4 preflight: the AMRFinderPlus tool and its database must agree.

`docs/design/07-bakta-optimization.md` treats this as a first-class check rather
than a runtime surprise, and records what happens without one: Bakta 1.12.1
ships database `2024-12-18.1`, AMRFinderPlus 4.2.7 requires `>= 2025-09-22.2`,
neither has another `osx-arm64` build, and Bakta's expert system fails outright
with "AMR expert system failed". A *partial* database auto-created by Bakta's own
internal AMRFinder call also had to be deleted by hand, because an implicit
mid-run database update is precisely what `allow_database_update: false` exists
to prevent.

**The pin is read from `config/references.tsv`**, which describes itself as the
reproducibility contract and carries `database_version` and `version_status` for
every reference. It is deliberately not hard-coded here. A floor written into
this module would be a second statement of the same fact, free to drift from the
file that owns it, and the drift would be invisible until a run compared a
cohort against an unrecorded database.

Three states are refused, each for a distinct reason:

* **absent** - the tool is installed and cannot run at all;
* **unpinned** - `references.tsv` calls this "a reported failure, not a pass";
* **mismatched** - the dangerous one, because it looks healthy. Half a cohort
  annotated against one database and half against another is unreproducible,
  and nothing downstream would report it.

Nothing here downloads or updates anything. Provisioning is out of band by
design, so the refusal names the tool that does it rather than doing it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from ..errors import PipelineError, ToolExecutionError
from ..io.tsv import read_tsv
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.amr")

#: The reference row that governs stage 4.
AMR_REFERENCE_ID = "ref_amrfinderplus"

#: The tool that provisions the database, named in refusals so the reader has
#: somewhere to go. It is never invoked from here.
PROVISIONING_TOOL = "amrfinder_update"

#: Filenames AMRFinderPlus has used to record a database version. Reused from
#: the pilot module rather than restated, so there is one list.
_VERSION_FILENAMES = ("version.txt", "amrfinderplus.version")

#: A provisioned database always carries the indexed protein FASTA.
_DATABASE_MARKER = "AMRProt.fa*"

#: Placeholder the references table uses for "not recorded".
_UNPINNED = "UNPINNED"


@dataclass(frozen=True)
class ReferencePin:
    """What `references.tsv` claims about a reference."""

    reference_id: str
    tool_version: str
    database: str
    database_version: str
    version_status: str

    @property
    def is_pinned(self) -> bool:
        return (
            self.version_status.strip().lower() == "pinned"
            and self.database_version.strip().upper() != _UNPINNED
        )


def read_reference_pin(
    references_path: Path, reference_id: str = AMR_REFERENCE_ID
) -> ReferencePin:
    """Read one row of `references.tsv`.

    Raises:
        PipelineError: The file is missing, or has no row for `reference_id`. A
            missing row is not the same as an unpinned one - the first means the
            contract does not describe stage 4 at all.
    """
    references_path = Path(references_path)
    if not references_path.is_file():
        raise PipelineError(
            f"{references_path} is missing, so stage 4 has no pinned reference "
            "to check its database against. Without it a run cannot record "
            "which database produced its results."
        )

    rows = read_tsv(
        references_path,
        required_columns=("reference_id", "database_version", "version_status"),
    )
    for row in rows:
        if str(row.get("reference_id", "")).strip() == reference_id:
            return ReferencePin(
                reference_id=reference_id,
                tool_version=str(row.get("tool_version") or "").strip(),
                database=str(row.get("database") or "").strip(),
                database_version=str(row.get("database_version") or "").strip(),
                version_status=str(row.get("version_status") or "").strip(),
            )

    raise PipelineError(
        f"{references_path} has no row for {reference_id!r}, so the "
        "reproducibility contract does not describe stage 4's database. Add the "
        "row rather than defaulting it: a default is a second, invisible "
        "statement of what the cohort was annotated against.",
        references_path=str(references_path),
        reference_id=reference_id,
        present=",".join(sorted(str(r.get("reference_id")) for r in rows)),
    )


def database_version_on_disk(database_dir: Path) -> Optional[str]:
    """The version recorded inside `database_dir`, or None if it records none.

    Reads the file the tool wrote rather than trusting the directory name, so a
    directory renamed to look current cannot pass for a current database.
    """
    database_dir = Path(database_dir)
    for name in _VERSION_FILENAMES:
        path = database_dir / name
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="replace").strip()
            if text:
                return text
    return None


def looks_provisioned(database_dir: Path) -> bool:
    """Whether `database_dir` holds something AMRFinderPlus could use.

    A directory that merely exists is the *partial* database docs/07 records,
    so existence is not enough: the indexed protein FASTA has to be there.
    """
    database_dir = Path(database_dir)
    if not database_dir.is_dir():
        return False
    return any(database_dir.glob(_DATABASE_MARKER))


def preflight_database(
    *,
    database_dir: Path,
    references_path: Path,
    reference_id: str = AMR_REFERENCE_ID,
) -> str:
    """Verify tool and database agree, and return the agreed version.

    Returns:
        The database version, identical in the reference file and on disk.

    Raises:
        PipelineError: The database is absent, the reference is unpinned, or the
            two disagree. Each message names what to change.
    """
    database_dir = Path(database_dir)

    if not database_dir.is_dir():
        raise PipelineError(
            f"the AMRFinderPlus database directory {database_dir} does not "
            "exist, so stage 4 cannot run. The database is provisioned out of "
            f"band - run {PROVISIONING_TOOL} yourself, or point "
            "`amr.database_path` at an existing one. This pipeline never "
            "downloads it during a run, because an implicit mid-run update "
            "would split the cohort across two databases with nothing to say so.",
            database_dir=str(database_dir),
            provisioning_tool=PROVISIONING_TOOL,
        )

    if not looks_provisioned(database_dir):
        raise PipelineError(
            f"{database_dir} exists but holds no indexed protein FASTA "
            f"({_DATABASE_MARKER}), so it is an incomplete database rather than "
            "a usable one. This is the partial-database state that docs/design/"
            "07-bakta-optimization.md records having to be deleted by hand after "
            "a tool created it implicitly. Re-provision it rather than letting "
            f"stage 4 discover the gap per genome; {PROVISIONING_TOOL} rebuilds "
            "it.",
            database_dir=str(database_dir),
            expected_marker=_DATABASE_MARKER,
            provisioning_tool=PROVISIONING_TOOL,
        )

    pin = read_reference_pin(references_path, reference_id)
    if not pin.is_pinned:
        raise PipelineError(
            f"{reference_id} in {references_path} is unpinned "
            f"(database_version={pin.database_version!r}, "
            f"version_status={pin.version_status!r}), and that file's own header "
            "calls an unpinned reference 'a reported failure, not a pass'. Set "
            f"the real database_version and version_status=pinned for "
            f"{reference_id} before a REAL run, so a result can be traced to the "
            "database that produced it.",
            references_path=str(references_path),
            reference_id=reference_id,
            database_version=pin.database_version,
            version_status=pin.version_status,
        )

    on_disk = database_version_on_disk(database_dir)
    if on_disk is None:
        raise PipelineError(
            f"{database_dir} has no version file "
            f"({' or '.join(_VERSION_FILENAMES)}), so the version cannot be "
            f"checked against the {pin.database_version} that {references_path} "
            "records. An unverifiable database is treated as a mismatched one: "
            "the point of the check is that every genome is annotated against a "
            "database the run can name.",
            database_dir=str(database_dir),
            expected=pin.database_version,
        )

    if on_disk != pin.database_version:
        raise PipelineError(
            f"the AMRFinderPlus database on disk is {on_disk} but "
            f"{references_path} records {pin.database_version} for "
            f"{reference_id}. Refusing rather than running: annotating part of "
            "the cohort against one database and part against another is "
            "unreproducible, and no later stage would report it. Either point "
            "the run at the pinned database or update the pin deliberately.",
            database_dir=str(database_dir),
            on_disk=on_disk,
            pinned=pin.database_version,
            reference_id=reference_id,
        )

    LOGGER.info(
        "Stage 4 preflight: AMRFinderPlus database %s matches %s",
        on_disk, reference_id,
    )
    return on_disk


def execute_command(
    command: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    timeout: int = 3600,
) -> "CommandResult":
    """Run `command` and return a :class:`CommandResult`.

    This is the *executor* the adapter's `runner` seam defaults to. It used to
    default to :class:`CommandResult` itself - a result dataclass - and then be
    called as ``self._runner(command=...)``, which raised

        TypeError: CommandResult.__init__() missing 4 required
        positional arguments

    on the first real run. Every test injected a fake, so the default was never
    executed and the seam was never executable; `amr` had therefore never
    screened anything, which is consistent with its long-standing "no REAL
    caller" status.

    A non-zero exit is *returned*, not raised. The adapter decides what a failure
    means, because "the tool failed" and "the tool wrote no report" are different
    conditions and only the caller can tell them apart.
    """
    import subprocess
    import time

    from .external import CommandResult

    started = time.monotonic()
    try:
        completed = subprocess.run(
            list(command),
            capture_output=True,
            text=True,
            check=False,
            cwd=str(cwd) if cwd else None,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolExecutionError(
            "AMRFinderPlus timed out",
            command=" ".join(command), timeout_s=timeout,
        ) from exc
    except OSError as exc:
        raise ToolExecutionError(
            "could not launch AMRFinderPlus",
            command=" ".join(command), error=str(exc),
        ) from exc
    return CommandResult(
        command=list(command),
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        duration_s=round(time.monotonic() - started, 3),
    )


#: The line `amrfinder --list_organisms` prints its choices on.
_ORGANISM_LIST_MARKER = "Available --organism options:"


def derive_organism_flag(organism_name: str) -> str:
    """Turn a configured species name into AMRFinderPlus' `--organism` spelling.

    The database spells every choice with underscores (``Pseudomonas_aeruginosa``)
    while `config/science.yaml` spells it the way a biologist writes it
    (``Pseudomonas aeruginosa``), so whitespace is the only substitution made.
    Nothing else is guessed: a subspecies, a pathovar or a strain would each
    need a rule, and inventing one would produce a plausible-looking flag the
    tool does not accept. The result is validated against the tool's own list by
    :func:`resolve_organism`, which is what makes "only whitespace is
    substituted" safe rather than merely convenient.
    """
    text = " ".join(str(organism_name or "").split())
    if not text:
        raise PipelineError(
            "config/science.yaml has no `organism.name`, so stage 4 has nothing "
            "to pass to AMRFinderPlus --organism. Refusing rather than running "
            "without it: without the flag AMRFinderPlus reports no POINT rows at "
            "all, so an imipenem cohort would be screened for acquired genes only "
            "and would look negative for point mutations rather than unscreened "
            "for them.",
            organism_name=None if organism_name is None else str(organism_name),
        )
    return text.replace(" ", "_")


def list_organisms(
    executable: str,
    database_dir: Path,
    *,
    runner=None,
) -> Tuple[str, ...]:
    """The `--organism` values this AMRFinderPlus database actually accepts.

    Asked of the tool rather than written down here. A hard-coded copy would be a
    second statement of what the database contains, free to drift from the
    database, and the drift would be invisible until a run passed the tool a
    species it has never heard of.
    """
    runner = runner or execute_command
    result = runner(
        command=[
            str(executable), "--database", str(database_dir), "--list_organisms",
        ]
    )
    if int(getattr(result, "returncode", 0) or 0) != 0:
        raise PipelineError(
            f"{executable} --list_organisms failed, so stage 4 cannot tell "
            "whether the configured organism is one the database screens for. "
            f"Refusing rather than passing --organism unchecked. "
            f"stderr: {(getattr(result, 'stderr', '') or '').strip()[-500:]}",
            executable=str(executable),
            database_dir=str(database_dir),
            returncode=int(getattr(result, "returncode", 0) or 0),
        )
    for line in (getattr(result, "stdout", "") or "").splitlines():
        if line.strip().startswith(_ORGANISM_LIST_MARKER):
            _, _, choices = line.partition(_ORGANISM_LIST_MARKER)
            return tuple(
                choice.strip()
                for choice in choices.split(",")
                if choice.strip()
            )
    raise PipelineError(
        f"{executable} --list_organisms printed no "
        f"{_ORGANISM_LIST_MARKER!r} line, so the accepted organism names cannot "
        "be read and stage 4 will not guess one. This is a tool or version "
        "change, not a configuration problem.",
        executable=str(executable),
        database_dir=str(database_dir),
        stdout=(getattr(result, "stdout", "") or "")[-500:],
    )


def resolve_organism(
    organism_name: str,
    *,
    executable: str,
    database_dir: Path,
    runner=None,
) -> str:
    """The validated `--organism` value for a configured species name.

    Returns:
        The value to pass, guaranteed to be one this database lists.

    Raises:
        PipelineError: The name is empty, or the derived value is not one the
            tool's own list contains. The message names the value it tried, the
            name it came from, and how to see the alternatives - because a
            silently ignored `--organism` is exactly the failure this replaces.
    """
    flag = derive_organism_flag(organism_name)
    available = list_organisms(executable, database_dir, runner=runner)
    if flag not in available:
        raise PipelineError(
            f"config organism.name {organism_name!r} derives --organism "
            f"{flag!r}, which the AMRFinderPlus database at {database_dir} does "
            f"not list, so stage 4 refuses rather than running an unvalidated "
            f"screen. Run `amrfinder --database {database_dir} "
            f"--list_organisms` to see the {len(available)} names this database "
            "accepts, and set config/science.yaml `organism.name` to one of "
            "them. Note that a mismatch here is not cosmetic: without a valid "
            "--organism AMRFinderPlus reports no POINT rows, so the "
            "point-mutation screen silently produces nothing.",
            organism_name=None if organism_name is None else str(organism_name),
            derived=flag,
            database_dir=str(database_dir),
            n_available=len(available),
        )
    LOGGER.info(
        "Stage 4: --organism %s (from organism.name %r), validated against the "
        "database's own list of %d", flag, organism_name, len(available),
    )
    return flag


def tool_version(
    executable: str,
    *,
    runner=None,
) -> Optional[str]:
    """The installed AMRFinderPlus version, as the tool reports it.

    Read from the tool rather than from `config/references.tsv`, because the pin
    records what the cohort is *reproducible against* and this records what
    *ran* - the two are worth keeping apart, and `stages.amr.preflight_database`
    already refuses when they disagree.
    """
    runner = runner or execute_command
    try:
        result = runner(command=[str(executable), "--version"])
    except ToolExecutionError as exc:
        LOGGER.warning("could not read the AMRFinderPlus version: %s", exc)
        return None
    if int(getattr(result, "returncode", 0) or 0) != 0:
        LOGGER.warning(
            "`%s --version` exited %s; recording the tool version as unknown",
            executable, getattr(result, "returncode", None),
        )
        return None
    text = (getattr(result, "stdout", "") or "").strip()
    return text.splitlines()[0].strip() if text else None


def commands_for_isolate(
    *,
    executable: str,
    assembly: Path,
    database_dir: Path,
    threads: int,
    organism: str,
    out_path: Optional[Path] = None,
) -> List[List[str]]:
    """The AMRFinderPlus invocations for one assembly.

    Exposed as a function rather than hidden inside the adapter so the no-update
    guarantee is testable: every command this module can produce is inspectable,
    and none carries `--update` or `--force_update`.

    **`--nucleotide`, never `--input`.** AMRFinderPlus 4.x renamed the flag, and
    the installed 4.2.7 rejects the old spelling outright with
    `ERROR: "--input" is not a valid option`. Verified against the binary rather
    than taken from release notes.

    **`--output` names a file, not `-`.** The v3 adapter piped the report to
    stdout; the v4 parser reads a path, because a report large enough to
    interleave with a tool's progress output on one stream is not worth parsing.

    **`--organism` is passed, and omitting it was a silent loss of a whole
    class of result.**  `config/science.yaml` has always recorded
    `organism.name: "Pseudomonas aeruginosa"` (taxid 287), and the project spec
    (`.scratch/imipenem-gwas-dashboard/spec.md:343`) names the invocation as
    ``AMRFinderPlus (--organism Pseudomonas_aeruginosa)``.  An earlier version
    of this docstring claimed the opposite - that "nothing in this pipeline
    records a species to put in it" - which was false; the config recorded it
    and the code simply never read it.

    The consequence was not cosmetic.  `--organism` is what makes AMRFinderPlus
    report **POINT** rows - single amino-acid substitutions in a gene the
    database associates with resistance in that species - and omitting it
    searched the whole database instead, so every row came back
    `Element subtype=POINT` free and the `variant` column was empty for every
    determinant in the cohort.  For an imipenem study, whose whole mechanistic
    story is loss of function in `oprD` and friends, that is the screen that
    matters most being switched off without a warning.

    It also **narrows** the gene set, not only widens it: AMRFinderPlus reports
    the genes it associates with the named organism, so switching the flag on
    removes rows as well as adding them. A reader comparing two cohorts screened
    with and without the flag must expect differences in both directions, and a
    cohort screened without it must be described as unscreened for point
    mutations rather than as screened and negative.
    """
    return [
        [
            str(executable),
            "--nucleotide", str(assembly),
            "--database", str(database_dir),
            "--organism", organism,
            "--plus",
            "--threads", str(int(threads)),
            # The caller names the output. Deriving it from the assembly put the
            # report *beside the sequence* - for a smoke run, inside
            # db/smoke_genomes/ - while the adapter then looked for it in the
            # report directory and raised "wrote no report" on a successful
            # screen. Tool output must not accumulate in a genome directory.
            "--output", str(out_path) if out_path else f"{assembly}.amrfinder.tsv",
        ]
    ]
