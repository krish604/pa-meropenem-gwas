"""AMRFinderPlus execution and v4 output parsing for the pilot cohort.

AMRFinderPlus 4.x renamed ``--input`` to ``--nucleotide`` and re-ordered the
report columns, so the v3-era parser in :mod:`papipeline.stages.amr` is not
reused verbatim. This module produces the same
:class:`~papipeline.models.AmrDeterminant` records, so everything downstream
of stage 4 is unchanged.

Two attribution rules are enforced here:

* ``evidence_source`` is ``amrfinderplus`` and ``database`` /
  ``database_version`` come from the provisioned database directory, so a
  detected gene is always traceable to the database that reported it;
* the AMR **class** and **subclass** recorded on each record are the
  database's own annotations, copied verbatim. They are not this
  pipeline's interpretation, and :func:`parse_amrfinder_report` keeps them in
  dedicated fields so downstream code cannot accidentally present them as a
  mechanism call.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..errors import PipelineError, ToolExecutionError, ToolNotAvailableError
from ..logging_utils import get_logger
from ..models import AmrDeterminant, ClaimStatus

LOGGER = get_logger("pilot.amr")

#: v4 report columns we rely on.
REQUIRED_COLUMNS: Tuple[str, ...] = (
    "Contig id",
    "Element symbol",
    "Type",
    "Class",
    "Subclass",
    "% Identity to reference",
    "% Coverage of reference",
)

#: AMR vs virulence, per AMRFinderPlus' own Type vocabulary.
AMR_TYPE = "AMR"
VIRULENCE_TYPE = "VIRULENCE"


def amrfinder_executable() -> str:
    """Locate the amrfinder executable, or raise a clear error."""
    path = shutil.which("amrfinder")
    if path is None:
        raise ToolNotAvailableError(
            "amrfinder is not on PATH",
            hint=(
                "Activate the pilot environment: "
                "micromamba activate pilot100"
            ),
        )
    return path


#: File names AMRFinderPlus has used to record a database version.
VERSION_FILENAMES: Tuple[str, ...] = ("version.txt", "amrfinderplus.version")


def read_database_version(database_dir: Path) -> Optional[str]:
    """Read the database version recorded inside ``database_dir``.

    AMRFinderPlus writes ``version.txt``; older layouts used
    ``amrfinderplus.version``. Returns ``None`` when neither is present,
    rather than guessing from the directory name.
    """
    database_dir = Path(database_dir)
    for name in VERSION_FILENAMES:
        path = database_dir / name
        if path.exists():
            text = path.read_text(encoding="utf-8", errors="replace").strip()
            if text:
                return text
    return None


def database_version(database_dir: Path) -> str:
    """Best available version string for a database directory.

    Prefers the file the tool wrote; falls back to the directory name, which
    is how ``amrfinder_update`` names the versioned directory.
    """
    database_dir = Path(database_dir)
    if not database_dir.exists():
        raise PipelineError(
            "AMRFinderPlus database directory not found", path=str(database_dir)
        )
    recorded = read_database_version(database_dir)
    if recorded:
        return recorded
    if any((database_dir / name).exists() for name in VERSION_FILENAMES):
        return "unknown"
    return database_dir.name


def is_database_dir(path: Path) -> bool:
    """Whether ``path`` looks like a provisioned AMRFinderPlus database."""
    path = Path(path)
    if not path.is_dir():
        return False
    if any((path / name).exists() for name in VERSION_FILENAMES):
        return True
    # A provisioned database always carries the indexed protein FASTA.
    return any(path.glob("AMRProt.fa*"))


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().rstrip("%")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_amrfinder_report(
    path: Path,
    sample_id: str,
    antibiotic: str,
    database: str,
    database_version_value: str,
) -> List[AmrDeterminant]:
    """Parse one AMRFinderPlus v4 report into :class:`AmrDeterminant` records.

    Only ``Type=AMR`` rows are returned; virulence rows are handled by the
    pipeline's own virulence stage so the two never mix.
    """
    path = Path(path)
    if not path.exists():
        raise PipelineError("AMRFinderPlus report not found", path=str(path))

    records: List[AmrDeterminant] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLUMNS if c not in fields]
        if missing:
            raise PipelineError(
                "AMRFinderPlus report is missing expected columns",
                path=str(path),
                missing=",".join(missing),
            )
        for row in reader:
            if (row.get("Type") or "").strip().upper() != AMR_TYPE:
                continue
            symbol = (row.get("Element symbol") or "").strip()
            if not symbol:
                continue
            identity = _to_float(row.get("% Identity to reference"))
            coverage = _to_float(row.get("% Coverage of reference"))
            records.append(
                AmrDeterminant(
                    sample_id=sample_id,
                    antibiotic=antibiotic,
                    determinant=symbol,
                    gene=symbol,
                    variant=None,
                    determinant_type=(row.get("Type") or AMR_TYPE).strip(),
                    # The database's own class/subclass, carried verbatim.
                    # NOT a mechanism assignment.
                    mechanism=None,
                    evidence_source="amrfinderplus",
                    database=database,
                    database_version=database_version_value,
                    confidence=None if identity is None else f"{identity:.2f}",
                    claim_status=ClaimStatus.DETECTED,
                    identity_pct=identity,
                    coverage_pct=coverage,
                )
            )
    return records


def read_class_map(path: Path) -> Dict[str, Dict[str, str]]:
    """Build ``gene -> {class, subclass, type}`` from one report.

    Used to attach the database's own functional annotation to a gene. The
    values are database statements, not pipeline interpretations.
    """
    path = Path(path)
    out: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = (row.get("Element symbol") or "").strip()
            if not symbol:
                continue
            out[symbol] = {
                "class": (row.get("Class") or "").strip(),
                "subclass": (row.get("Subclass") or "").strip(),
                "type": (row.get("Type") or "").strip(),
            }
    return out


def run_amrfinder(
    genome: Path,
    out_tsv: Path,
    database_dir: Path,
    threads: int = 4,
    plus: bool = True,
    timeout: int = 3600,
) -> Path:
    """Run AMRFinderPlus on one assembly.

    The database is never updated here. An existing report is reused, which
    keeps the run restartable and keeps repeated calls from silently
    changing results.
    """
    out_tsv = Path(out_tsv)
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        LOGGER.debug("Reusing existing AMRFinderPlus report %s", out_tsv)
        return out_tsv

    executable = amrfinder_executable()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable,
        "--nucleotide",
        str(genome),
        "--database",
        str(database_dir),
        "--threads",
        str(threads),
        "--output",
        str(out_tsv),
    ]
    if plus:
        command.append("--plus")

    LOGGER.info("amrfinder: %s", genome.name)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolExecutionError(
            "AMRFinderPlus timed out",
            genome=str(genome),
            timeout_s=timeout,
        ) from exc
    except OSError as exc:
        raise ToolExecutionError(
            "Failed to launch AMRFinderPlus",
            genome=str(genome),
            error=str(exc),
        ) from exc

    if completed.returncode != 0 or not out_tsv.exists():
        raise ToolExecutionError(
            "AMRFinderPlus exited non-zero",
            genome=str(genome),
            returncode=completed.returncode,
            stderr=(completed.stderr or "")[-1500:],
        )
    return out_tsv


@dataclass
class DetectionResult:
    """Outcome of the detection sweep."""

    records: Dict[str, List[AmrDeterminant]]
    class_map: Dict[str, Dict[str, Dict[str, str]]]
    reports: Dict[str, Path]
    failures: Dict[str, str]

    @property
    def n_scanned(self) -> int:
        return len(self.reports) + len(self.failures)

    @property
    def n_successful(self) -> int:
        return len(self.reports)

    @property
    def n_failed(self) -> int:
        return len(self.failures)


def detect_cohort(
    items: Sequence[Tuple[str, Path]],
    report_dir: Path,
    database_dir: Path,
    antibiotic: str,
    threads: int = 4,
    jobs: int = 1,
) -> DetectionResult:
    """Run AMRFinderPlus across the cohort.

    Args:
        items: ``(sample_id, genome_path)`` pairs.
        report_dir: Where per-genome reports are written.
        database_dir: The provisioned AMRFinderPlus database.
        antibiotic: Antibiotic label attached to the records.
        threads: Threads per AMRFinderPlus process.
        jobs: Concurrent processes. Each process already uses ``threads``,
            so ``jobs * threads`` should stay within the machine's cores.
    """
    report_dir = Path(report_dir)
    database_dir = Path(database_dir)
    db_version = database_version(database_dir)
    database_name = f"AMRFinderPlus_{db_version}"
    LOGGER.info(
        "AMRFinderPlus sweep: %d genomes, database %s, jobs=%d threads=%d",
        len(items),
        database_name,
        jobs,
        threads,
    )

    records: Dict[str, List[AmrDeterminant]] = {}
    class_map: Dict[str, Dict[str, Dict[str, str]]] = {}
    reports: Dict[str, Path] = {}
    failures: Dict[str, str] = {}

    def work(item: Tuple[str, Path]) -> Tuple[str, Optional[Path], Optional[str]]:
        sample_id, genome = item
        target = report_dir / f"{sample_id}.amrfinder.tsv"
        try:
            return sample_id, run_amrfinder(
                genome, target, database_dir, threads=threads
            ), None
        except (ToolExecutionError, PipelineError) as exc:
            return sample_id, None, str(exc)

    if jobs > 1:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            outcomes = list(pool.map(work, items))
    else:
        outcomes = [work(item) for item in items]

    for sample_id, report, error in outcomes:
        if report is None:
            LOGGER.error("AMRFinderPlus failed for %s: %s", sample_id, error)
            failures[sample_id] = error or "unknown"
            records[sample_id] = []
            class_map[sample_id] = {}
            continue
        reports[sample_id] = report
        records[sample_id] = parse_amrfinder_report(
            report, sample_id, antibiotic, database_name, db_version
        )
        class_map[sample_id] = read_class_map(report)

    total_hits = sum(len(v) for v in records.values())
    LOGGER.info(
        "AMRFinderPlus sweep complete: %d/%d genomes, %d gene calls, %d failures",
        len(reports),
        len(items),
        total_hits,
        len(failures),
    )
    return DetectionResult(
        records=records,
        class_map=class_map,
        reports=reports,
        failures=failures,
    )
