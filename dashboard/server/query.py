"""§5.1 — the query-parameter schema, stated once and honoured server-side.

Every list endpoint pages, sorts and filters **server-side** (UI-D5), and
returns `{items, total, offset, limit, sort, filters, basis}`. `total` is the
count over the whole declared artefact, computed before paging, and it is the
number the UI shows next to a page — never `items.length`.

| param | type | default |
|---|---|---|
| `offset` | int ≥ 0 | 0 |
| `limit` | int 1–1000 | 200 |
| `sort` | `<column>` or `-<column>` | the declared column order |
| `q` | substring over the declared searchable columns | — |
| `filter.<column>` | exact match; `__ne`, `__in`, `__null` also accepted | — |
| `columns` | csv projection | all |
| `include_total` | bool | 1 |

Two rules that are refusals rather than tolerances:

- an **unknown sort column** is a `400` naming the column, never a silent
  fallback to the declared order — a sort that quietly did not apply looks
  like a sort;
- an **unknown projection or filter column** is a `400` for the same reason.

`limit` is clamped, and the clamp is reported in `limit_applied`, so a client
that asked for 5000 rows learns that it got 1000.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

#: The documented default page size (openapi `parameters.Limit`).
DEFAULT_LIMIT = 200
MAX_LIMIT = 1000
DEFAULT_OFFSET = 0

#: Separators a client may send filters with. Comma is the documented one.
_SEPARATORS = (",",)


class QueryError(Exception):
    """A bad query parameter. `message` is the sentence the client sees."""

    def __init__(self, message: str, detail: Optional[Dict[str, Any]] = None) -> None:
        super().__init__(message)
        self.message = message
        self.detail = detail or {}


@dataclass
class Filter:
    """One parsed `filter.<column>[__op]` clause."""

    column: str
    op: str
    value: Any

    def matches(self, row: Mapping[str, Any]) -> bool:
        actual = row.get(self.column)
        if self.op == "eq":
            return _loose_eq(actual, self.value)
        if self.op == "ne":
            return not _loose_eq(actual, self.value)
        if self.op == "in":
            return any(_loose_eq(actual, candidate) for candidate in self.value)
        if self.op == "null":
            absent = actual is None or str(actual).strip() == ""
            return absent if self.value else not absent
        if self.op == "min":
            # `filter.virulence_min=3` — "at least 3 virulence factors". An
            # absent cell (no virulence table) matches no minimum: a minimum
            # over something not measured is not a measurement of zero.
            if actual is None or isinstance(actual, str):
                return False
            try:
                return float(actual) >= float(self.value)
            except (TypeError, ValueError):
                return False
        raise QueryError(f"unsupported filter operator {self.op!r}")


def _loose_eq(actual: Any, expected: Any) -> bool:
    """Exact match, compared as text when either side is text.

    A `filter.n_both=3` from a query string has to match the integer 3 in a
    row, or every numeric filter would return nothing — which reads as "no
    rows match" rather than as a type mismatch.

    A **list** cell matches when any of its members matches. `amr_genes` is
    many-valued on purpose, so `filter.amr_genes=oprD` means "carries oprD",
    not "equals the list `["oprD"]`".
    """
    if isinstance(expected, (list, tuple)):
        return any(_loose_eq(actual, item) for item in expected)
    if isinstance(actual, (list, tuple)):
        return any(_loose_eq(item, expected) for item in actual)
    if actual is None:
        return expected is None or str(expected) in ("", "null", "none")
    if isinstance(expected, bool) or isinstance(actual, bool):
        # A boolean cell matches "1"/"true"/"yes" and "0"/"false"/"no", so
        # `filter.has_amr_gene=1` finds the isolates that carry a gene.
        if isinstance(actual, bool):
            truthy = actual
        else:
            truthy = str(actual).strip().lower() in ("1", "true", "yes")
        if isinstance(expected, bool):
            wanted = expected
        else:
            wanted = str(expected).strip().lower() in ("1", "true", "yes")
        return truthy == wanted
    if isinstance(actual, (int, float)) or isinstance(expected, (int, float)):
        try:
            return float(actual) == float(expected)
        except (TypeError, ValueError):
            return False
    return str(actual).strip() == str(expected).strip()


def parse_filter_params(raw: Mapping[str, Sequence[str]], columns: Sequence[str]) -> Dict[str, Any]:
    """Turn the repeated `filter.*` query params into a dict for the response.

    `raw` is the multi-dict FastAPI hands over, keys being `filter.col`,
    `filter.col__ne`, `filter.col__in`, `filter.col__null`.
    """
    known = set(columns)
    parsed: Dict[str, Any] = {}
    clauses: List[Filter] = []
    for key in sorted(raw):
        if not key.startswith("filter."):
            continue
        remainder = key[len("filter."):]
        column, _, suffix = remainder.partition("__")
        op = "eq"
        if suffix:
            op = {"ne": "ne", "in": "in", "null": "null"}.get(suffix, "")
            if not op:
                raise QueryError(
                    f"unknown filter modifier {('__' + suffix)!r} on column "
                    f"{column!r}. Accepted: __ne, __in, __null.",
                )
        if column == "virulence_min":
            # `virulence_min` is a derived filter, not a column: it means
            # "at least N virulence factors" and is applied as a `>=` over the
            # row's own `n_virulence` cell.
            column = "virulence_min"
            op = "min"
        if column not in known:
            raise QueryError(
                f"unknown filter column {column!r}. This table carries: "
                f"{', '.join(sorted(known)) or '(no columns)'}. A filter on a "
                f"column the table does not have would silently match nothing.",
                detail={"column": column, "columns": sorted(known)},
            )
        values = list(raw[key])
        value: Any
        if op == "in":
            value = []
            for item in values:
                for part in item.split(_SEPARATORS[0]):
                    part = part.strip()
                    if part:
                        value.append(part)
            clauses.append(Filter(column, op, value))
        elif op == "null":
            value = str(values[0]).strip() in ("1", "true", "yes")
            clauses.append(Filter(column, op, value))
        else:
            value = values[0]
            clauses.append(Filter(column, op, value))
        parsed[key] = values if len(values) > 1 else values[0]
    parsed["_clauses"] = clauses
    return parsed


@dataclass
class PageRequest:
    """One parsed list request."""

    offset: int = DEFAULT_OFFSET
    limit: int = DEFAULT_LIMIT
    limit_applied: Optional[int] = None
    sort: Optional[str] = None
    q: Optional[str] = None
    columns: Optional[Tuple[str, ...]] = None
    include_total: bool = True
    filters: Dict[str, Any] = field(default_factory=dict)
    clauses: Tuple[Filter, ...] = ()

    @property
    def has_filter(self) -> bool:
        return bool(self.clauses)


def parse_page(
    *,
    offset: Optional[int] = None,
    limit: Optional[int] = None,
    sort: Optional[str] = None,
    q: Optional[str] = None,
    columns: Optional[str] = None,
    include_total: Optional[int] = None,
    filter_params: Optional[Mapping[str, Sequence[str]]] = None,
    available: Sequence[str],
    searchable: Sequence[str] = (),
) -> PageRequest:
    """Validate and clamp one list request against the table's own columns."""
    known = list(available)

    # -- offset ---------------------------------------------------------
    chosen_offset = DEFAULT_OFFSET if offset is None else int(offset)
    if chosen_offset < 0:
        raise QueryError(
            f"offset must be >= 0; got {chosen_offset}.",
        )

    # -- limit, clamped and the clamp reported --------------------------
    requested = DEFAULT_LIMIT if limit is None else int(limit)
    if requested < 1:
        raise QueryError(f"limit must be >= 1; got {requested}.")
    chosen_limit = min(requested, MAX_LIMIT)
    limit_applied = chosen_limit if chosen_limit != requested else None

    # -- sort -----------------------------------------------------------
    chosen_sort: Optional[str] = None
    if sort is not None and str(sort).strip():
        chosen_sort = str(sort).strip()
        bare = chosen_sort[1:] if chosen_sort.startswith("-") else chosen_sort
        if bare not in known:
            raise QueryError(
                f"unknown sort column {bare!r}. This table carries: "
                f"{', '.join(known) or '(no columns)'}. An unknown sort is a "
                f"400 rather than a silent fallback to the declared order - a "
                f"sort that quietly did not apply looks like a sort.",
                detail={"column": bare, "columns": known},
            )

    # -- projection -----------------------------------------------------
    chosen_columns: Optional[Tuple[str, ...]] = None
    if columns is not None and str(columns).strip():
        parts = tuple(p.strip() for p in str(columns).split(",") if p.strip())
        unknown = [p for p in parts if p not in known]
        if unknown:
            raise QueryError(
                f"unknown projection column(s): {', '.join(unknown)}. This "
                f"table carries: {', '.join(known) or '(no columns)'}.",
                detail={"columns": unknown, "available": known},
            )
        chosen_columns = parts

    # -- filters --------------------------------------------------------
    parsed_filters: Dict[str, Any] = {}
    clauses: Tuple[Filter, ...] = ()
    if filter_params:
        parsed = parse_filter_params(filter_params, known)
        clauses = tuple(parsed.pop("_clauses"))
        parsed_filters = parsed

    return PageRequest(
        offset=chosen_offset,
        limit=chosen_limit,
        limit_applied=limit_applied,
        sort=chosen_sort,
        q=str(q) if q else None,
        columns=chosen_columns,
        include_total=bool(1 if include_total is None else include_total),
        filters=parsed_filters,
        clauses=clauses,
    )


def apply_page(
    rows: Sequence[Mapping[str, Any]],
    request: PageRequest,
    *,
    searchable: Sequence[str] = (),
    key_fn: Optional[Callable[[Mapping[str, Any]], Any]] = None,
) -> Tuple[List[Dict[str, Any]], int, int]:
    """Filter, sort, page — all server-side.

    Returns `(items, total, total_unfiltered)`. `total` is the filtered count
    when a filter is active, so "12 of 900" is sayable from
    `total_unfiltered`.

    An absent cell matches no exact filter except `__null`, and is not found
    by `q`. Rendering it as `""` would make a missing value equal an empty one.
    """
    total_unfiltered = len(rows)
    selected: List[Mapping[str, Any]] = list(rows)

    for clause in request.clauses:
        selected = [row for row in selected if clause.matches(row)]

    if request.q:
        needle = request.q.strip().lower()
        if needle:
            cols = list(searchable) or list(selected[0].keys()) if selected else list(searchable)
            selected = [
                row
                for row in selected
                if any(
                    row.get(c) is not None and needle in str(row.get(c)).lower()
                    for c in cols
                )
            ]

    if request.sort:
        descending = request.sort.startswith("-")
        bare = request.sort[1:] if descending else request.sort
        selected = sorted(
            selected, key=lambda row: _sort_key(row.get(bare)), reverse=descending
        )
    elif key_fn is not None:
        selected = sorted(selected, key=key_fn)

    total = len(selected)
    window = selected[request.offset: request.offset + request.limit]
    items = [dict(row) for row in window]
    if request.columns:
        keep = set(request.columns)
        items = [{k: v for k, v in row.items() if k in keep} for row in items]
    return items, total, total_unfiltered


def _sort_key(value: Any, _unused: Any = None) -> Tuple[int, Any]:
    """``None`` last, numbers before text, both sorted naturally.

    A missing cell is an unknown, not the smallest thing in the column, so it
    goes last rather than reading as -1.
    """
    if value is None:
        return (2, "")
    if isinstance(value, bool):
        return (0, float(value))
    if isinstance(value, (int, float)):
        return (0, float(value))
    text = str(value)
    try:
        return (0, float(text))
    except (TypeError, ValueError):
        return (1, text)


def page_meta(
    request: PageRequest,
    *,
    total: int,
    total_unfiltered: int,
    basis: Mapping[str, Any],
) -> Dict[str, Any]:
    """The `meta` block every list response carries.

    `total` is over the whole declared artefact before paging. `limit_applied`
    appears only when the limit was clamped. `filters` echoes what was actually
    honoured, so a client can tell a filter it sent from a filter it got.
    """
    meta: Dict[str, Any] = {
        "total": total,
        "offset": request.offset,
        "limit": request.limit,
        "sort": request.sort,
        "filters": {k: v for k, v in request.filters.items()},
        "basis": dict(basis),
    }
    if request.limit_applied is not None:
        meta["limit_applied"] = request.limit_applied
    if request.has_filter or request.q:
        meta["total_unfiltered"] = total_unfiltered
    meta["include_total"] = request.include_total
    return meta


def basis(
    *,
    n: int,
    artefact: str,
    path: Optional[str],
    rows_total: int,
    filtered: bool = False,
    min_samples: Optional[int] = None,
) -> Dict[str, Any]:
    """UI-D7 — every statistic arrives with the artefact it came from.

    `min_samples` travels with it so the D3 power flag can fire without the
    client having to know the configured minimum.
    """
    payload: Dict[str, Any] = {
        "n": int(n),
        "artefact": artefact,
        "path": str(path) if path is not None else None,
        "rows_total": int(rows_total),
        "filtered": bool(filtered),
    }
    if min_samples is not None:
        payload["min_samples"] = int(min_samples)
        payload["underpowered"] = int(n) < int(min_samples)
    return payload


__all__ = [
    "DEFAULT_LIMIT",
    "MAX_LIMIT",
    "Filter",
    "PageRequest",
    "QueryError",
    "apply_page",
    "basis",
    "page_meta",
    "parse_filter_params",
    "parse_page",
]