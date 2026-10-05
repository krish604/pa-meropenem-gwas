"""Strict, dependency-free tabular I/O.

The data contract (see ``docs/data_contract.md``) is enforced here rather than
in each stage, so that a malformed file fails once with a clear message
instead of producing silently wrong tables downstream.

Rules enforced by :func:`read_tsv`:

* the header row must contain the required columns;
* a row must not have more fields than the header (malformed TSV);
* a missing trailing field is a hard error, not a silent ``None``;
* a sentinel value (``.``, ``NA``, ``nan``, ``None``, ``null``) becomes
  ``None`` so that missing data stays missing;
* duplicate values in a declared unique column are a hard error.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence

from ..errors import DataContractError, DuplicateSampleError
from ..logging_utils import get_logger

LOGGER = get_logger("io.tsv")

#: Strings treated as "missing" rather than as data.
MISSING_SENTINELS = frozenset({"", ".", "na", "nan", "none", "null", "-"})

TAB = "\t"


def _clean_cell(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    stripped = value.strip()
    if stripped.lower() in MISSING_SENTINELS:
        return None
    return stripped


def read_tsv(
    path: Path,
    required_columns: Sequence[str] = (),
    unique_columns: Sequence[str] = (),
    unique_together: Sequence[Sequence[str]] = (),
    min_columns: int = 1,
    comment_prefix: Optional[str] = "#",
) -> List[Dict[str, Optional[str]]]:
    """Read a TSV file into a list of dicts.

    Comment lines starting with ``comment_prefix`` are skipped anywhere in the
    file, so outputs written by :func:`write_tsv` with a provenance header
    round-trip through this reader.

    Args:
        path: File to read. Must exist.
        required_columns: Columns that must be present in the header.
        unique_columns: Columns whose values must be unique across rows.
        unique_together: Column groups whose *combined* values must be unique,
            e.g. ``[("sample_id", "antibiotic")]``.
        min_columns: Minimum number of columns a row may contain.
        comment_prefix: Line prefix treated as a comment. ``None`` disables.

    Returns:
        One dict per data row. Absent/sentinel values are ``None``.

    Raises:
        DataContractError: File missing, empty, or violates the contract.
        DuplicateSampleError: A unique column contains a repeated value.
    """
    path = Path(path)
    if not path.exists():
        raise DataContractError("TSV file not found", path=str(path))
    if path.stat().st_size == 0:
        raise DataContractError("TSV file is empty", path=str(path))

    with path.open("r", encoding="utf-8", newline="") as handle:
        raw_rows: List[List[str]] = []
        for raw_line in handle:
            if comment_prefix and raw_line.startswith(comment_prefix):
                continue
            raw_rows.append(next(csv.reader([raw_line], delimiter=TAB)))

    if not raw_rows:
        raise DataContractError(
            "TSV file contains only comments and no header row", path=str(path)
        )

    header = [c.strip() for c in raw_rows[0]]
    if len(set(header)) != len(header):
        duplicates = sorted({c for c in header if header.count(c) > 1})
        raise DataContractError(
            "TSV header contains duplicate column names",
            path=str(path),
            duplicates=",".join(duplicates),
        )

    missing = [c for c in required_columns if c not in header]
    if missing:
        raise DataContractError(
            "TSV header is missing required columns",
            path=str(path),
            missing=",".join(missing),
            found=",".join(header),
        )

    rows: List[Dict[str, Optional[str]]] = []
    for line_number, raw_row in enumerate(raw_rows[1:], start=2):
        if not raw_row or all(cell.strip() == "" for cell in raw_row):
            continue
        if len(raw_row) > len(header):
            raise DataContractError(
                "Malformed TSV: row has more fields than the header",
                path=str(path),
                line=line_number,
                fields=len(raw_row),
                header_fields=len(header),
            )
        if len(raw_row) < min_columns:
            raise DataContractError(
                "Malformed TSV: row has too few fields",
                path=str(path),
                line=line_number,
                fields=len(raw_row),
                minimum=min_columns,
            )
        row: Dict[str, Optional[str]] = {}
        for index, column in enumerate(header):
            value = raw_row[index] if index < len(raw_row) else None
            row[column] = _clean_cell(value)
        rows.append(row)

    if not rows:
        raise DataContractError(
            "TSV file contains a header but no data rows", path=str(path)
        )

    def _enforce_uniqueness(
        key_columns: Sequence[str], label: str
    ) -> None:
        absent = [c for c in key_columns if c not in header]
        if absent:
            raise DataContractError(
                "Cannot enforce uniqueness on absent column(s)",
                path=str(path),
                column=",".join(absent),
            )
        seen: Dict[tuple, int] = {}
        for index, row in enumerate(rows, start=2):
            values = tuple(row.get(c) for c in key_columns)
            if all(v is None for v in values):
                continue
            if any(v is None for v in values):
                # A partially-populated composite key is a data error: it
                # cannot be deduplicated and would silently join wrong.
                raise DataContractError(
                    f"Incompletely populated {label} key",
                    path=str(path),
                    line=index,
                    column=",".join(key_columns),
                )
            if values in seen:
                raise DuplicateSampleError(
                    f"Duplicate {label} key",
                    path=str(path),
                    column=",".join(key_columns),
                    value=",".join(str(v) for v in values),
                    first_line=seen[values],
                    second_line=index,
                )
            seen[values] = index

    for column in unique_columns:
        _enforce_uniqueness([column], "single-column")
    for group in unique_together:
        _enforce_uniqueness(list(group), "composite")

    LOGGER.debug("Read %d rows from %s", len(rows), path)
    return rows


def write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    columns: Sequence[str],
    header_comment: Optional[Sequence[str]] = None,
) -> Path:
    """Write rows to a TSV file with a stable column order.

    Args:
        path: Destination path. Parent directories are created.
        rows: Iterable of mappings.
        columns: Column order to emit.
        header_comment: Optional ``#``-prefixed provenance lines written
            above the header. These are ignored by :func:`read_tsv` callers
            that skip comments; use :func:`read_tsv` on the raw header only.

    Returns:
        The written path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        if header_comment:
            for line in header_comment:
                handle.write(f"# {line}\n")
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter=TAB, extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({c: _format(row.get(c)) for c in columns})
    LOGGER.debug("Wrote %s", path)
    return path


def _format(value: Any) -> str:
    if value is None:
        return "."
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (list, tuple, set)):
        return ";".join(str(v) for v in value)
    if isinstance(value, Mapping):
        return ";".join(f"{k}={v}" for k, v in sorted(value.items()))
    return str(value)


def iter_tsv_dicts(
    path: Path, skip_comments: bool = True
) -> Iterator[Dict[str, Optional[str]]]:
    """Stream a TSV, skipping ``#`` comment lines. Used for large tables."""
    path = Path(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        header: Optional[List[str]] = None
        for raw in handle:
            if skip_comments and raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split(TAB)
            if header is None:
                header = [f.strip() for f in fields]
                continue
            if not fields or all(f.strip() == "" for f in fields):
                continue
            yield {
                header[i]: _clean_cell(fields[i] if i < len(fields) else None)
                for i in range(len(header))
            }


def write_lines(path: Path, lines: Iterable[str]) -> Path:
    """Write plain text lines (newick, fasta, logs) to disk."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(f"{line}\n")
    return path


def write_fasta(path: Path, records: Iterable[tuple]) -> Path:
    """Write ``(header, sequence)`` pairs as FASTA with 80-column wrapping."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")
    return path
