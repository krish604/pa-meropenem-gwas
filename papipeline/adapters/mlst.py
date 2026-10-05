"""Real ``mlst`` invocation, and parsing of the output it produces.

The stage reads a table; this module is what fills that table on a REAL run.
It runs the tool and reports what the tool said. It does not decide what a
profile *means* - the typed/partial/no_call contract, and the
``min_typing_locus_match`` floor, live in :func:`papipeline.stages.mlst.parse_row`
and are applied there. A second copy of those rules here would be free to drift
from the stage's, and the drift would show up as a REAL run disagreeing with a
TEST run on identical data.

Three things this module is careful about, each because the alternative is a
wrong answer rather than an error.

**``--legacy`` is required, not cosmetic.** The installed tool's own help
records that it "requires --scheme", and without it the header row naming the
loci is absent. The loci are what make an ST traceable: ``MlstCall.alleles``
exists so a sequence type can always be traced back to the alleles behind it,
and the header is the only place those names come from. A hard-coded list of
the seven *P. aeruginosa* loci would be the alternative, and it would silently
mis-parse any other scheme the config ever names.

**A missing locus is not an allele.** ``mlst`` writes ``-`` for a locus it did
not match. Passing that through as an allele would inflate the matched-loci
count, and the count is what the stage uses to decide typed versus partial - so
a two-of-seven profile could be promoted to ``typed``.

**No ST is not ST ``0``.** An unresolved ST is written ``-`` and stays ``None``,
because a sequence type is defined only for a complete profile and emitting a
number would be a fabrication.

Output on stdout, progress on stderr: ``--quiet`` keeps the parse target clean.
"""

from __future__ import annotations

import csv
import io
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from ..errors import PipelineError, ToolExecutionError
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.mlst")

#: Sentinel alleles `mlst` writes for a locus it did not resolve. `-` is no
#: match; `?` is a partial match. Neither is an allele number.
_UNMATCHED_ALLELE = "-"

#: Columns `mlst` puts ahead of the per-locus allele columns.
_LEADING_COLUMNS = ("FILE", "SCHEME", "ST")


@dataclass(frozen=True)
class MlstResult:
    """One isolate's MLST call, as the tool reported it.

    `found` distinguishes "the tool ran and typed nothing" from "the tool
    failed". Collapsing them would fail a cohort over one untypable genome, so
    an empty result is a successful call that happens to carry no ST.
    """

    sample_id: str
    st: Optional[str] = None
    alleles: Dict[str, str] = field(default_factory=dict)
    partial_loci: List[str] = field(default_factory=list)
    scheme: Optional[str] = None
    found: bool = False

    def as_stage_row(self, allele_database: Optional[str]) -> Dict[str, Optional[str]]:
        """Render as a `stages.mlst` table row.

        The `alleles` field is the stage's `locus:allele;locus:allele` form.
        A partial `?` allele is kept: dropping it would understate how much of
        the profile was actually determined, and the stage counts matched loci
        from this string.

        `MLST_status` starts as `typed` - the optimistic claim - because
        :func:`~papipeline.stages.mlst.parse_row` is built to demote it. The
        stage re-checks the profile against `min_typing_locus_match` and
        against whether a usable ST was actually returned, so a four-of-seven
        profile arrives as `partial` and an empty one as `no_call`. Claiming
        here and verifying there keeps the status vocabulary in one place
        instead of restating it.
        """
        return {
            "sample_id": self.sample_id,
            "ST": self.st or "",
            "alleles": ";".join(
                f"{locus}:{allele}" for locus, allele in self.alleles.items()
            ),
            "MLST_status": "typed" if self.found else "no_call",
            "mlst_scheme": self.scheme or "",
            "allele_database": allele_database,
        }


def query_command(
    executable: str,
    *,
    scheme: str,
    assembly: Path,
    threads: int = 1,
    datadir: Optional[str] = None,
    outfile: Optional[Path] = None,
) -> List[str]:
    """The `mlst` invocation for one assembly.

    Every value arrives from configuration; nothing here restates the science.
    """
    command = [
        executable,
        "--scheme", str(scheme),
        "--legacy",
        "--quiet",
        "--threads", str(int(threads)),
    ]
    if datadir:
        command += ["--datadir", str(datadir)]
    if outfile:
        command += ["--outfile", str(outfile)]
    command.append(str(assembly))
    return command


def parse_legacy_output(
    stdout: str, sample_id: str, *, path: Optional[Path] = None
) -> MlstResult:
    """Parse `mlst --legacy` stdout into an :class:`MlstResult`.

    The FILE column is authoritative for matching a row to its sample, rather
    than line order: a batched call returns many rows, and position-based
    matching would silently mis-assign them.
    """
    lines = [ln for ln in stdout.splitlines() if ln.strip()]
    if not lines:
        raise PipelineError(
            "mlst produced no output at all, not even a header row; the tool "
            "did not run as expected",
            sample_id=sample_id,
            path=str(path) if path else None,
        )

    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    fieldnames = list(reader.fieldnames or [])
    if not fieldnames:
        raise PipelineError(
            "mlst output has no header row, so the loci cannot be named",
            sample_id=sample_id,
            columns=",".join(fieldnames),
        )
    for column in _LEADING_COLUMNS:
        if column not in fieldnames:
            raise PipelineError(
                f"mlst output is missing the {column!r} column; --legacy output "
                f"is required to read an allele profile",
                sample_id=sample_id,
                columns=",".join(fieldnames),
            )

    loci = [f for f in fieldnames if f not in _LEADING_COLUMNS]
    if not loci:
        raise PipelineError(
            "mlst reported no loci, so there is no profile to record",
            sample_id=sample_id,
            columns=",".join(fieldnames),
        )

    rows = list(reader)
    if not rows:
        # Header only: the tool ran and typed nothing. `no_call`, not a failure.
        LOGGER.info("mlst returned no row for %s", sample_id)
        return MlstResult(sample_id=sample_id, found=False)

    row = _row_for(rows, sample_id)
    if row is None:
        raise PipelineError(
            f"mlst returned {len(rows)} row(s) but none for {sample_id}; the "
            "sample's result cannot be identified",
            sample_id=sample_id,
            samples=",".join(sorted(str(r.get("FILE", "")) for r in rows)),
        )

    raw_st = (row.get("ST") or "").strip()
    alleles: Dict[str, str] = {}
    partial: List[str] = []
    for locus in loci:
        value = (row.get(locus) or "").strip()
        if value == _UNMATCHED_ALLELE:
            continue
        if value == "?":
            # A partial match is a determination, just an incomplete one. It
            # is kept so the profile states what was actually resolved.
            partial.append(locus)
        if value:
            alleles[locus] = value

    return MlstResult(
        sample_id=sample_id,
        st=None if raw_st in ("", _UNMATCHED_ALLELE) else raw_st,
        alleles=alleles,
        partial_loci=partial,
        scheme=(row.get("SCHEME") or "").strip() or None,
        found=True,
    )


def _row_for(rows: Sequence[Dict[str, str]], sample_id: str) -> Optional[Dict[str, str]]:
    """The row belonging to `sample_id`, matched on the FILE column's stem.

    `mlst` echoes the path it was given, so the sample id is the file stem
    rather than the whole value.
    """
    for row in rows:
        name = str(row.get("FILE", "")).strip()
        if name == sample_id or Path(name).stem == sample_id:
            return row
    return None


def call_isolate(
    sample_id: str,
    *,
    assembly: Path,
    scheme: str,
    workdir: Path,
    tools: Dict[str, str],
    threads: int = 1,
    datadir: Optional[str] = None,
) -> MlstResult:
    """Run `mlst` on one assembly and return its call.

    `tools` maps tool name to resolved executable, so the caller is responsible
    for having verified the tool exists - the same division as the aligner,
    where `require_callers` runs once for the cohort rather than per isolate.
    """
    executable = tools.get("mlst")
    if not executable:
        raise PipelineError(
            "no resolved path for `mlst`; the tool must be probed before the "
            "cohort is called",
            sample_id=sample_id,
        )
    if not Path(assembly).is_file():
        raise PipelineError(
            f"mlst needs an assembly to read and {assembly} does not exist",
            sample_id=sample_id,
            assembly=str(assembly),
        )

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    outfile = workdir / f"{sample_id}.mlst.tsv"

    command = query_command(
        executable,
        scheme=scheme,
        assembly=Path(assembly),
        threads=threads,
        datadir=datadir,
        outfile=outfile,
    )
    LOGGER.info("Stage 3: %s -> mlst --scheme %s", sample_id, scheme)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise ToolExecutionError(
            f"could not execute mlst for {sample_id}: {exc}",
            sample_id=sample_id,
            command=" ".join(command),
        ) from exc

    if completed.returncode != 0:
        raise ToolExecutionError(
            f"mlst failed for {sample_id} (exit {completed.returncode})",
            sample_id=sample_id,
            command=" ".join(command),
            stderr=(completed.stderr or "").strip()[-2000:],
        )

    # `--outfile` is used, but the tool echoes its report to stdout as well;
    # reading the file is the deterministic choice, and stdout is the fallback
    # for a tool that honours only one of the two.
    text = outfile.read_text(encoding="utf-8") if outfile.is_file() else (
        completed.stdout or ""
    )
    return parse_legacy_output(text, sample_id, path=outfile)
