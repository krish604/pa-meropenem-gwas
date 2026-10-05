"""UI-D5 — the byte-offset index, cached in the state directory.

Paging a 200,000-row variant table by re-reading it per request is quadratic
over a session. The index maps each key to the byte range of its row, so a
page is `seek` + `readline`.

The format, the key and the invalidation rule are DESIGN §9 UI-D5 verbatim:

```json
{"version": 1, "key": "<sha256>",
 "source": {"path": "...", "size": 0, "mtime_ns": 0},
 "header": [...], "key_column": "sample_id",
 "offsets": {"TEST_PA_001": [4096, 4241]}, "row_count": 2}
```

The index is **disposable by construction** — deleting the state directory
costs a rebuild and nothing else — which is what makes writing it at all
consistent with UI-D1. Nothing under the results root is ever created,
truncated or renamed, and the counting reader below exists so a test can prove
that a huge table is never fully loaded (Q2).
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

from papipeline.io.tsv import MISSING_SENTINELS

from .state import StateDir

#: Bumped when the on-disk shape changes. A version mismatch discards the
#: entry rather than reading it (UI-D5 invalidation).
INDEX_VERSION = 1

#: How many rows a page may ask for. The openapi clamps at 1000 and this is
#: the ceiling the server enforces before that.
HARD_ROW_CAP = 1000


def index_key(path: Path) -> str:
    """The cache key: sha256 of realpath, size and mtime_ns.

    Size and mtime together, so a truncated file and a rewritten file both
    miss. There is no time-based expiry: mtime and size are a statement about
    the file, and polling `stat` is cheaper and more honest than a guess about
    how long a cache lasts. A run still in progress rewrites its tables, so its
    mtime moves and the index is rebuilt — which is why the index is never
    treated as a result.
    """
    real = Path(path).resolve()
    st = real.stat()
    payload = f"{real}\0{st.st_size}\0{st.st_mtime_ns}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _short(key: str) -> str:
    """The `key16` in `<stage>.<key16>.idx.json`."""
    return key[:16]


@dataclass
class CountingReader:
    """A read-only text handle that records how much it read.

    This is the seam QA asserts on (Q2): a huge TSV with tens of thousands of
    columns must never be fully loaded, and `bytes_read` is the evidence. It is
    also the reason `iter_rows` reads line by line rather than through
    `read_tsv`, which materialises the whole file as a list of dicts.
    """

    path: Path
    encoding: str = "utf-8"
    bytes_read: int = 0
    rows_read: int = 0

    def lines(self) -> Iterator[str]:
        with open(self.path, "r", encoding=self.encoding, newline="") as handle:
            while True:
                line = handle.readline()
                if not line:
                    return
                self.bytes_read += len(line.encode(self.encoding, errors="replace"))
                yield line

    def record_row(self, nbytes: int) -> None:
        self.rows_read += 1


def _clean(value: Optional[str]) -> Optional[str]:
    """`MISSING_SENTINELS` handling, matching `io.tsv._clean_cell`.

    Reused rather than reimplemented because a dashboard that renders `0` for a
    sentinel creates a second provenance gap (UI-D2).
    """
    if value is None:
        return None
    stripped = value.strip()
    if stripped.lower() in MISSING_SENTINELS:
        return None
    return stripped


def iter_rows(path: Path, *, reader: Optional[CountingReader] = None) -> Iterator[Tuple[List[str], Dict[str, Any]]]:
    """Stream a TSV, yielding `(header, row)` per data row.

    The `#`-comment rule and the sentinel rule are the pipeline's, taken from
    `io.tsv`. What this does differently is that it never raises: a header with
    duplicate names raises from `read_tsv` (trap 1) and here the caller has
    already been told the table is not produced, so a tolerant second pass is
    only ever used for information.
    """
    reader = reader or CountingReader(path)
    header: Optional[List[str]] = None
    for raw in reader.lines():
        if raw.startswith("#"):
            continue
        stripped = raw.rstrip("\n").rstrip("\r")
        if header is None:
            header = [c.strip() for c in stripped.split("\t")]
            continue
        if not stripped or not stripped.strip():
            continue
        fields = stripped.split("\t")
        yield header, {
            header[i]: _clean(fields[i] if i < len(fields) else None)
            for i in range(len(header))
        }


@dataclass
class TableIndex:
    """One table's byte offsets, as found on disk.

    `offsets` maps a key value to the byte range of *the first* row carrying it,
    in the format DESIGN §9 UI-D5 declares. A duplicate key keeps the first
    range and `duplicate_keys` counts the rest, because the pipeline's own
    contract treats a duplicate sample id as an error and this reader reports it
    rather than picking a winner silently.

    `spans` is the ordered `[key, start, end]` list for every row, which is what
    paging and sorting need — `offsets` alone cannot express a sort or a page
    past the first occurrence of a key. It is stored alongside `offsets` rather
    than instead of them, so the on-disk shape stays the one DESIGN fixes.
    """

    path: Path
    key: str
    header: Tuple[str, ...]
    key_column: Optional[str]
    offsets: Mapping[str, Tuple[int, int]]
    duplicate_keys: int = 0
    row_count: int = 0
    spans: Sequence[Sequence[Any]] = ()

    def rows_for(self, key_value: str) -> List[Dict[str, Any]]:
        """Every row carrying ``key_value``, read by seek + readline."""
        span = self.offsets.get(str(key_value))
        if span is None:
            return []
        start, end = span
        rows: List[Dict[str, Any]] = []
        with open(self.path, "r", encoding="utf-8", newline="") as handle:
            handle.seek(start)
            remaining = max(0, end - start)
            while remaining > 0:
                line = handle.readline()
                if not line:
                    break
                remaining -= len(line.encode("utf-8", errors="replace"))
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").rstrip("\r").split("\t")
                rows.append(
                    {
                        self.header[i]: _clean(fields[i] if i < len(fields) else None)
                        for i in range(len(self.header))
                    }
                )
        return rows

    def nth_range(self, index: int) -> Optional[Tuple[int, int]]:
        """The byte range of the ``index``-th data row, in file order."""
        order = self._ordered_keys
        if index < 0 or index >= len(order):
            return None
        return self.offsets[order[index]]

    def as_spans(self) -> List[Tuple[str, int, int]]:
        """`(key, start, end)` per row, in file order.

        Rebuilt from `spans` when the index carries them, and from `offsets`
        otherwise so an index written before `spans` existed still pages.
        """
        if self.spans:
            return [(str(s[0]), int(s[1]), int(s[2])) for s in self.spans]
        return [
            (key, span[0], span[1])
            for key, span in self.offsets.items()
        ]

    def _recompute_order(self) -> List[str]:
        return list(self.offsets.keys())

    @property
    def _ordered_keys(self) -> List[str]:
        cached = getattr(self, "_order_cache", None)
        if cached is None:
            cached = self._recompute_order()
            object.__setattr__(self, "_order_cache", cached)
        return cached


class IndexStore:
    """The state directory's index directory, and the cache in front of it."""

    def __init__(self, state: StateDir) -> None:
        self.state = state

    # -- paths -----------------------------------------------------------
    def path_for(self, stage: str, key: str) -> Path:
        return self.state.index_dir / f"{stage}.{_short(key)}.idx.json"

    # -- build -----------------------------------------------------------
    def build(self, stage: str, path: Path, *, key_column: Optional[str] = None) -> TableIndex:
        """Index one table by streaming it once.

        A huge table is read line by line and only the offsets are retained, so
        a 200,000-row file with tens of thousands of columns costs one pass and
        never a full materialisation. The `CountingReader` passed to QA is the
        evidence.
        """
        path = Path(path)
        key = index_key(path)
        reader = CountingReader(path)
        header: List[str] = []
        offsets: Dict[str, Tuple[int, int]] = {}
        spans: List[Tuple[str, int, int]] = []
        duplicates = 0
        row_count = 0
        position = 0
        first = True
        for raw in reader.lines():
            length = len(raw.encode("utf-8", errors="replace"))
            start = position
            position += length
            if raw.startswith("#"):
                continue
            if first:
                header = [c.strip() for c in raw.rstrip("\n").rstrip("\r").split("\t")]
                first = False
                continue
            if not raw.strip():
                continue
            row_count += 1
            reader.record_row(length)
            value = ""
            if key_column is not None and key_column in header:
                fields = raw.rstrip("\n").rstrip("\r").split("\t")
                key_at = header.index(key_column)
                cleaned = _clean(fields[key_at] if key_at < len(fields) else None)
                value = cleaned or ""
                if cleaned is not None:
                    if cleaned in offsets:
                        # A duplicate key keeps the first range and the second is
                        # counted: the pipeline's own contract treats a repeated
                        # sample id as an error, so this reader reports it rather
                        # than picking a winner silently.
                        duplicates += 1
                    else:
                        offsets[cleaned] = (start, start + length)
            spans.append((value, start, start + length))

        index = TableIndex(
            path=path,
            key=key,
            header=tuple(header),
            key_column=key_column,
            offsets=offsets,
            duplicate_keys=duplicates,
            row_count=row_count,
            spans=[[v, s, e] for v, s, e in spans],
        )
        self._persist(stage, key, index)
        return index

    # -- read ------------------------------------------------------------
    def get(self, stage: str, path: Path, *, key_column: Optional[str] = None) -> TableIndex:
        """The index for ``path``, rebuilt when the key moved.

        Invalidation is exactly the three UI-D5 conditions: the stored `key`
        differs from the key computed now, the `version` differs, or
        `source.path` differs. There is no time-based expiry.
        """
        path = Path(path)
        try:
            key = index_key(path)
        except OSError:
            # A file that cannot be stat'd has no index. Rebuilt next time it
            # can be.
            return TableIndex(
                path=path, key="", header=(), key_column=key_column, offsets={}
            )
        cached = self._load(stage, key, path)
        if cached is not None:
            return cached
        return self.build(stage, path, key_column=key_column)

    def _load(self, stage: str, key: str, path: Path) -> Optional[TableIndex]:
        target = self.path_for(stage, key)
        if not target.exists():
            return None
        try:
            payload = json.loads(target.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            # A truncated or unreadable index is discarded, not repaired: the
            # rebuild is cheap and a half-read index would be a wrong answer.
            return None
        if not isinstance(payload, Mapping):
            return None
        if payload.get("version") != INDEX_VERSION:
            return None
        if payload.get("key") != key:
            return None
        source = payload.get("source") or {}
        if source.get("path") != str(Path(path).resolve()):
            return None
        offsets_raw = payload.get("offsets") or {}
        if not isinstance(offsets_raw, Mapping):
            return None
        offsets: Dict[str, Tuple[int, int]] = {}
        for name, span in offsets_raw.items():
            try:
                offsets[str(name)] = (int(span[0]), int(span[1]))
            except (TypeError, ValueError, IndexError):
                return None
        return TableIndex(
            path=Path(path),
            key=key,
            header=tuple(str(h) for h in (payload.get("header") or ())),
            key_column=payload.get("key_column"),
            offsets=offsets,
            duplicate_keys=int(payload.get("duplicate_keys") or 0),
            row_count=int(payload.get("row_count") or 0),
            spans=payload.get("spans") or (),
        )

    def _persist(self, stage: str, key: str, index: TableIndex) -> None:
        """Write the index atomically, into the state directory only.

        Temp file in the same directory then `os.replace`, so a kill mid-write
        leaves the previous index rather than a truncated one (UI-D5).
        """
        payload = {
            "version": INDEX_VERSION,
            "key": index.key,
            "source": {
                "path": str(index.path.resolve()),
                "size": index.path.stat().st_size,
                "mtime_ns": index.path.stat().st_mtime_ns,
            },
            "header": list(index.header),
            "key_column": index.key_column,
            "offsets": {k: [v[0], v[1]] for k, v in index.offsets.items()},
            # Every row's key and byte range, in file order. What paging and
            # sorting need; `offsets` alone cannot express either.
            "spans": [list(s) for s in index.spans],
            "row_count": index.row_count,
            "duplicate_keys": index.duplicate_keys,
        }
        target = self.path_for(stage, index.key)
        target.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(target.parent),
            prefix=f".{target.name}.",
            suffix=".tmp",
            delete=False,
        )
        tmp = Path(handle.name)
        try:
            with handle:
                json.dump(payload, handle, sort_keys=True)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(str(tmp), str(target))
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

    # -- the QA seam -----------------------------------------------------
    def counting_reader(self, path: Path) -> CountingReader:
        """A reader QA can inspect after a request.

        Q2 asserts that a huge TSV was never fully loaded; this is what it
        measures. Nothing here caches the bytes.
        """
        return CountingReader(Path(path))


class PagedReader:
    """Stream a TSV into a page without materialising the whole table.

    UI-D5's claim is that a page is `seek` + `readline`, and the openapi's
    threat row is that "a huge file exhaust[s] memory" is answered by
    "responses are streamed; every table read is capped at limit<= 1000 rows".
    Both hold only if the reader never builds a list of every row's every cell.

    So there are two paths, and which one runs is stated rather than guessed:

    - **keyed** — the table was indexed on `key_column`, so a page is
      `seek(offset)` then `readline()` for the rows in range.
    - **streamed** — no index. Rows are read one at a time; only the *sort key*
      (or the filter column's value) is retained, plus the byte range, and the
      page is materialised afterwards by seeking to those ranges.

    Either way `reader` is a `CountingReader` and `bytes_read` is the evidence
    QA asserts on: a table with tens of thousands of columns must never be
    fully loaded, and this is how that is demonstrated rather than asserted.
    """

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.reader = CountingReader(self.path)

    # -- the streaming path ---------------------------------------------
    def scan(self, key_column: Optional[str]) -> Tuple[List[str], List[Tuple[int, int]]]:
        """One pass, keeping only the header, the key column and the ranges.

        Returns `(header, [(key_value_or_empty, start, end), ...])`. Cells are
        parsed only for the key column; everything else is measured by its byte
        range and read later if the page needs it.
        """
        header: List[str] = []
        spans: List[Tuple[str, int, int]] = []
        position = 0
        key_index = -1
        first = True
        for raw in self.reader.lines():
            length = len(raw.encode("utf-8", errors="replace"))
            start = position
            position += length
            if raw.startswith("#"):
                continue
            if first:
                header = [c.strip() for c in raw.rstrip("\n").rstrip("\r").split("\t")]
                first = False
                if key_column is not None:
                    key_index = header.index(key_column) if key_column in header else -1
                continue
            if not raw.strip():
                continue
            value = ""
            if key_index >= 0:
                fields = raw.rstrip("\n").rstrip("\r").split("\t")
                cleaned = _clean(fields[key_index] if key_index < len(fields) else None)
                value = cleaned or ""
            spans.append((value, start, start + length))
        return header, spans

    def read_range(
        self, start: int, end: int, header: Sequence[str], *, max_lines: int = 0
    ) -> List[Dict[str, Any]]:
        """Materialise the rows in one byte range. The only place cells are built.

        `max_lines` caps how many data lines are read (`0` is "the whole range"),
        so a caller that wants one row's value for a sort key does not pull the
        rest of the file in behind it. Bytes are counted either way, because the
        count is the evidence a test asserts on.
        """
        rows: List[Dict[str, Any]] = []
        with open(self.path, "r", encoding="utf-8", newline="") as handle:
            handle.seek(start)
            # `end == start` still reads one line when `max_lines` allows it: a
            # caller that wants a single row's cell passes `start == end`, and
            # treating `remaining == 0` as "read nothing" would return nothing
            # for a row that exists.
            remaining = max(0, end - start)
            taken = 0
            while max_lines <= 0 or taken < max_lines:
                if remaining <= 0 and taken > 0:
                    break
                line = handle.readline()
                if not line:
                    break
                size = len(line.encode("utf-8", errors="replace"))
                remaining -= size
                self.reader.bytes_read += size
                if line.startswith("#") or not line.strip():
                    continue
                taken += 1
                fields = line.rstrip("\n").rstrip("\r").split("\t")
                rows.append(
                    {
                        header[i]: _clean(fields[i] if i < len(fields) else None)
                        for i in range(len(header))
                    }
                )
        return rows

    def read_row(self, start: int, header: Sequence[str]) -> Dict[str, Any]:
        """One row, by seek + readline."""
        rows = self.read_range(start, start + 1 << 30, header)
        return rows[0] if rows else {}

    def page(
        self,
        *,
        header: Sequence[str],
        spans: Sequence[Tuple[str, int, int]],
        predicate: Optional[Any] = None,
        key_of: Optional[Any] = None,
        descending: bool = False,
        offset: int = 0,
        limit: int = 200,
        count_only: bool = False,
    ) -> Tuple[List[Dict[str, Any]], int, int, int]:
        """Filter, sort and page from byte ranges.

        Returns `(items, total, total_unfiltered, bytes_read)`. `bytes_read` is
        the total across the scan and the page materialisation, so a caller can
        compare it against `path.stat().st_size` and see that a 200,000-row file
        was not fully loaded.
        """
        kept: List[Tuple[Any, int, int]] = []
        total_unfiltered = len(spans)
        for value, start, end in spans:
            if predicate is not None and not predicate(value):
                continue
            kept.append((value if key_of is None else key_of(value, start), start, end))
        total = len(kept)
        if key_of is not None:
            kept.sort(key=lambda entry: _nulls_last(entry[0]), reverse=descending)
        window = kept[offset: offset + limit]
        items: List[Dict[str, Any]] = []
        if not count_only:
            for _key, start, end in window:
                items.extend(self.read_range(start, end, header))
        return items, total, total_unfiltered, self.reader.bytes_read


@dataclass(frozen=True)
class Peek:
    """Whether a table is present, and why not, without loading its rows.

    `report_tables.read_table` answers this question by reading the whole file,
    which is right for a ten-row summary and wrong for a 200,000-row variant
    table. This is the same answer for one pass that keeps nothing: existence,
    emptiness, the header-only refusal, and a row count.

    The refusal sentences are the ones `report_tables` uses, because a client
    that branches on `present` and reads `reason` must not have to care which
    reader produced it.
    """

    present: bool
    reason: str
    path: Path
    header: Tuple[str, ...] = ()
    n_rows: int = 0
    duplicate_header: bool = False


def peek_table(path: Path, *, name: str = "") -> Peek:
    """Presence, header and row count in one streaming pass. Never raises."""
    path = Path(path)
    if not path.exists():
        return Peek(False, f"no file at {path}", path)
    try:
        if path.stat().st_size == 0:
            return Peek(False, f"{path} is present but empty", path)
    except OSError as exc:
        return Peek(False, f"{path} cannot be stat'd: {exc}", path)

    header: List[str] = []
    rows = 0
    first = True
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                if not raw.strip():
                    continue
                if first:
                    header = [c.strip() for c in raw.rstrip("\n").rstrip("\r").split("\t")]
                    first = False
                    if len(set(header)) != len(header):
                        duplicates = sorted({c for c in header if header.count(c) > 1})
                        # Trap 1: `read_tsv` refuses this, and the refusal is the
                        # right answer — but it is a `not produced` with a named
                        # reason, not a 500.
                        return Peek(
                            False,
                            f"{path} is unreadable: TSV header contains duplicate "
                            f"column names ({', '.join(duplicates)})",
                            path,
                            tuple(header),
                            duplicate_header=True,
                        )
                    continue
                rows += 1
    except OSError as exc:
        return Peek(False, f"{path} is unreadable: {exc}", path)

    if not header:
        return Peek(False, f"{path} carries no header row", path)
    if rows == 0:
        return Peek(
            False,
            f"{path} carries a header and no rows. That is not a result - it is "
            f"what the file looks like when nothing was written into it, and "
            f"{name or 'the table'}'s output contract requires at least one row "
            f"(`contracts.stage_spec` adds `min_rows(table, 1)`)",
            path,
            tuple(header),
        )
    return Peek(True, "", path, tuple(header), rows)


def _nulls_last(value: Any) -> Tuple[int, Any]:
    """`None` last, numbers before text — same rule as `query._sort_key`."""
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


def sortable(values: Sequence[Any]) -> List[Tuple[int, Any]]:
    """A sort key that puts ``None`` last and numbers before text.

    `None` sorts after every real value rather than before it: an absent cell
    is not the smallest thing in a column, it is an unknown one, and rendering
    it as -1 would be a value the data does not carry.
    """
    out: List[Tuple[int, Any]] = []
    for value in values:
        if value is None:
            out.append((2, ""))
        elif isinstance(value, bool):
            out.append((0, int(value)))
        elif isinstance(value, (int, float)):
            out.append((0, float(value)))
        else:
            text = str(value)
            try:
                out.append((0, float(text)))
            except (TypeError, ValueError):
                out.append((1, text))
    return out


__all__ = [
    "HARD_ROW_CAP",
    "INDEX_VERSION",
    "CountingReader",
    "IndexStore",
    "PagedReader",
    "Peek",
    "TableIndex",
    "index_key",
    "iter_rows",
    "peek_table",
    "sortable",
]