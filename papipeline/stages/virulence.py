"""Stage 8 - virulence analysis.

Kept deliberately independent of the AMR stage: this module imports nothing
from :mod:`papipeline.stages.amr` or :mod:`papipeline.knowledge`, and it does
not consult the mechanism tables. Virulence factor carriage is reported as
carriage only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..adapters.virulence import (
    VirulencePreflightError,
    database_version_on_disk,
    ensure_blast_db,
    parse_vfdb_headers,
    screen_isolate,
)
from ..config.loader import PipelineConfig
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode, VirulenceFactor

LOGGER = get_logger("stages.virulence")

VIRULENCE_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "virulence_factor",
    "gene",
    "category",
    "database",
    "database_version",
    "confidence",
    "identity_pct",
)

REQUIRED = (
    "sample_id",
    "virulence_factor",
    "gene",
    "database",
    "database_version",
)


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(str(value).strip().rstrip("%"))
    except ValueError:
        return None


def load_virulence_factors(path: Path) -> List[VirulenceFactor]:
    """Load and normalise the stage 8 table."""
    rows = read_tsv(path, required_columns=REQUIRED)
    results: List[VirulenceFactor] = []
    for row in rows:
        results.append(
            VirulenceFactor(
                sample_id=str(row["sample_id"]),
                virulence_factor=str(row["virulence_factor"]),
                gene=str(row["gene"]),
                category=row.get("category"),
                database=str(row["database"]),
                database_version=str(row["database_version"]),
                confidence=row.get("confidence"),
                identity_pct=_to_float(row.get("identity_pct")),
            )
        )
    LOGGER.info("Stage 8: %d virulence factor detections", len(results))
    return results


def group_by_sample(
    records: Sequence[VirulenceFactor], manifest: SampleManifest
) -> Dict[str, List[VirulenceFactor]]:
    """Group virulence detections by sample, restricted to the manifest."""
    grouped: Dict[str, List[VirulenceFactor]] = {sid: [] for sid in manifest.sample_ids}
    orphans: List[str] = []
    for record in records:
        if record.sample_id in grouped:
            grouped[record.sample_id].append(record)
        else:
            orphans.append(record.sample_id)
    if orphans:
        LOGGER.warning(
            "%d virulence rows reference samples not in the manifest",
            len(orphans),
        )
    return grouped


def profile_string(records: Sequence[VirulenceFactor]) -> Optional[str]:
    """Comma-free profile string for the master table.

    Returns ``None`` when nothing was detected, so that "screened, nothing
    found" is distinguishable from "not screened".
    """
    if not records:
        return None
    return ",".join(sorted({r.virulence_factor for r in records}))


def category_counts(
    records: Sequence[VirulenceFactor],
) -> Dict[str, int]:
    """Counts per virulence category across a cohort."""
    counts: Dict[str, int] = {}
    for record in records:
        key = record.category or "uncategorised"
        counts[key] = counts.get(key, 0) + 1
    return counts


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
) -> Dict[str, List[VirulenceFactor]]:
    """Stage 5 entry point.

    REAL screens each prepared assembly with blastn against VFDB and writes
    ``virulence/virulence_factors.tsv``. TEST keeps reading the committed
    fixture, which is what every other stage does in TEST.
    """
    if mode is not RunMode.TEST:
        return _run_real(config, manifest, intermediate_root)

    path = Path(intermediate_root) / "virulence" / "virulence_factors.tsv"
    return group_by_sample(load_virulence_factors(path), manifest)


def _run_real(
    config: PipelineConfig,
    manifest: SampleManifest,
    intermediate_root: Path,
) -> Dict[str, List[VirulenceFactor]]:
    """blastn against VFDB, per isolate with an assembly.

    The guard this replaces (7e28840) was correct while there was no caller.
    abricate has no installable osx-arm64 build, but its own dependency list is
    `blast >=2.7`, so blastn is the same mechanism without the perl wrapper.
    """
    # `virulence:` is a top-level science.yaml section with no typed accessor,
    # so it is read from `raw` rather than inventing a dataclass for it.
    section = dict(config.raw.get("virulence") or {})
    db_root = config.machine.db_root() if config.machine is not None else config.root / "db"
    sequences = db_root / str(section.get("sequences_db") or "vfdb/sequences")
    if not sequences.is_file():
        # Not a cue to download. `allow_database_update` is false, so a missing
        # database is a provisioning failure to report, not something to fetch.
        raise VirulencePreflightError(
            f"VFDB sequences not found at {sequences}. This stage does not "
            "download databases; provision them under the overlay's db_root "
            "first. (allow_database_update is "
            f"{section.get('allow_database_update')})",
            path=str(sequences),
        )

    program = str(section.get("blast_program") or "blastn")
    threads = config.threads or 1
    entries = parse_vfdb_headers(sequences)
    index = ensure_blast_db(sequences, makeblastdb="makeblastdb")
    version = database_version_on_disk(sequences)
    database_name = str(section.get("database") or "VFDB")
    genome_root = config.assembly_root(RunMode.REAL)

    rows: List[Dict[str, Any]] = []
    for sample in manifest:
        # An isolate with no assembly is a normal outcome, not an error: it is
        # recorded as screened-with-nothing rather than silently dropped, which
        # would shrink the denominator for carriage rate.
        if not sample.assembly_path:
            LOGGER.warning(
                "Stage 5: no assembly for %s, recorded as screened with no "
                "virulence factor found", sample.sample_id,
            )
            continue
        found = screen_isolate(
            genome=Path(sample.assembly_path),
            entries=entries,
            database=index,
            program=program,
            threads=threads,
            evalue=float(section.get("evalue") or 1e-5),
            min_identity_pct=float(section.get("min_identity_pct") or 90.0),
            min_coverage_pct=float(section.get("min_coverage_pct") or 90.0),
            database_name=database_name,
            database_version=version,
        )
        for row in found:
            row["sample_id"] = sample.sample_id
        rows.extend(found)

    target = Path(intermediate_root) / "virulence"
    target.mkdir(parents=True, exist_ok=True)
    _write_virulence_tsv(target / "virulence_factors.tsv", rows)
    LOGGER.info(
        "Stage 5: screened %d assemblies against %s (%d factors, version %s)",
        sum(1 for s in manifest if s.assembly_path),
        database_name,
        len(entries),
        version,
    )
    return group_by_sample([_row_to_factor(r) for r in rows], manifest)


def _write_virulence_tsv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    import csv

    columns = VIRULENCE_COLUMNS + ("coverage_pct", "vfdb_accession")
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ("" if row.get(key) is None else row.get(key))
                    for key in columns
                }
            )


def _row_to_factor(row: Dict[str, Any]) -> VirulenceFactor:
    return VirulenceFactor(
        sample_id=str(row["sample_id"]),
        virulence_factor=str(row["virulence_factor"]),
        gene=str(row.get("gene") or ""),
        category=row.get("category") or None,
        database=str(row.get("database") or "VFDB"),
        database_version=str(row.get("database_version") or "unknown"),
        confidence=row.get("confidence") or None,
        identity_pct=_to_float(row.get("identity_pct")),
    )
