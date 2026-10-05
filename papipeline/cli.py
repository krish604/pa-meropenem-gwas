"""Helpers for the command-line wrappers in ``scripts/``.

Keeps the per-stage scripts to argument parsing and nothing else: the
reporting of "what did this stage produce" is defined once, here, so every
stage prints a consistent and readable summary.
"""

from __future__ import annotations

import dataclasses
from typing import Any, List, Mapping, Optional, Sequence

from .models import _clean


def _count(value: Any) -> Optional[int]:
    """Length of a sized value, or ``None`` if it has no length."""
    try:
        return len(value)
    except TypeError:
        return None


def describe(value: Any, depth: int = 0) -> str:
    """One-line, factual description of a stage result.

    Deliberately terse: it reports shapes and counts, never an
    interpretation of the numbers.
    """
    if value is None:
        return "(no result)"

    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        fields = [
            f"{f.name}={_short(getattr(value, f.name))}"
            for f in dataclasses.fields(value)
            if not f.name.startswith("_")
        ]
        return f"{type(value).__name__}({', '.join(fields[:6])})"

    if isinstance(value, Mapping):
        n = _count(value)
        non_empty = sum(1 for v in value.values() if v)
        return f"mapping of {n} key(s), {non_empty} non-empty"

    if isinstance(value, (list, tuple, set, frozenset)):
        n = len(value)
        if isinstance(value, tuple) and value and all(
            isinstance(v, (str, int, float)) for v in value
        ):
            return f"({', '.join(_short(v) for v in value[:4])})"
        if value and dataclasses.is_dataclass(value[0]):
            columns = list(value[0].to_row().keys()) if hasattr(value[0], "to_row") else []
            suffix = f", columns: {', '.join(columns)}" if columns else ""
            return f"{n} {type(value[0]).__name__} record(s){suffix}"
        return f"{n} item(s)"

    return _short(value)


def _short(value: Any, limit: int = 40) -> str:
    """Compact rendering of a scalar, truncated for readability."""
    text = "" if value is None else str(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def summarise(label: str, result: Any) -> List[str]:
    """Lines describing a stage result.

    Handles the shapes stages actually return: a record list, a
    sample-keyed mapping, a dataclass, or a tuple of those.
    """
    lines: List[str] = [f"{label}: {describe(result)}"]

    if isinstance(result, Mapping):
        for key in sorted(result)[:3]:
            lines.append(f"  {key}: {describe(result[key])}")
    elif isinstance(result, (list, tuple)) and result:
        for index, part in enumerate(result[:3]):
            lines.append(f"  [{index}] {describe(part)}")

    return lines


def write_rows(path, records: Sequence[Any], columns: Sequence[str]) -> None:
    """Write dataclass records to a TSV using their ``to_row`` projection."""
    from .io.tsv import write_tsv

    rows = [r.to_row() if hasattr(r, "to_row") else dict(r) for r in records]
    write_tsv(path, rows, list(columns))
