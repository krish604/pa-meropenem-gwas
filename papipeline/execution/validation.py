"""Declarative output validation.

Every stage declares, as data, what its output must look like. The
validator turns that declaration into a verdict. The point is that
"the file is there" is never the whole test, and that a stage which
produced the wrong thing says so instead of reporting success.

Three failure shapes are distinguished, because they need different
responses:

* **missing / empty / truncated** -> ``INCOMPLETE``. Retryable: often a
  process killed part-way through.
* **present but malformed** -> ``INVALID``. Not retryable: re-running the
  same command reproduces the same wrong bytes.
* **present, parseable, but degenerate** -> ``INVALID``. This is the class
  that has actually bitten this repository: 65 byte-identical VCFs, and a
  kinship matrix that was entirely 1.0. Both looked like success.

Validation reuses the pipeline's own strict readers
(:func:`papipeline.io.tsv.read_tsv`, :mod:`papipeline.io.fasta`) so that a
file the rest of the pipeline would reject is also rejected here. It adds
no biological thresholds of its own.
"""

from __future__ import annotations

import enum
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

from ..errors import DataContractError, PipelineError
from ..io.fasta import read_fasta
from ..io.tsv import read_tsv
from .state import StageState

#: Sentinel values ``read_tsv`` maps to ``None``. A column that is entirely
#: sentinel means the producer wrote a shape but no data.
_SENTINEL_TEXT = frozenset({".", "-"})


class CheckKind(str, enum.Enum):
    """The kinds of check a declaration can make."""

    EXISTS = "exists"
    NON_EMPTY = "non_empty"
    MIN_SIZE = "min_size"
    PARSEABLE = "parseable"
    COLUMNS = "columns"
    MIN_ROWS = "min_rows"
    NO_SENTINEL_COLUMN = "no_sentinel_column"
    SAMPLE_COUNT = "sample_count"
    IDENTIFIERS = "identifiers"
    MATRIX_NOT_CONSTANT = "matrix_not_constant"
    #: The file is named after the expected input genome. A raw tool output
    #: often carries no isolate column, so the file name is the only
    #: identifier available - and a mismatched name means the wrong genome
    #: was processed.
    IDENTIFIES_FILE = "identifies_file"


#: Which state a failed check implies. Missing/empty/truncated is
#: INCOMPLETE (retryable); everything else is INVALID (not retryable).
_INCOMPLETE_KINDS = frozenset(
    {CheckKind.EXISTS, CheckKind.NON_EMPTY, CheckKind.MIN_SIZE, CheckKind.MIN_ROWS}
)


@dataclass(frozen=True)
class Check:
    """One declarative assertion about one artefact.

    ``kind`` selects the assertion; the remaining fields are its arguments.
    Keeping them in one flat record (rather than a class hierarchy) is what
    lets a spec be written as data, serialised, and compared in tests.
    """

    kind: CheckKind
    path: Optional[Path] = None
    parser: Optional[str] = None
    #: Table dialect for row-reading checks: ``pipeline`` or ``bakta``.
    dialect: str = "pipeline"
    columns: Tuple[str, ...] = ()
    minimum: Optional[int] = None
    key: str = ""
    expected: Optional[int] = None
    values: Tuple[str, ...] = ()
    description: str = ""

    def describe(self) -> str:
        if self.description:
            return self.description
        target = self.path.name if self.path else "artefact"
        if self.kind is CheckKind.EXISTS:
            return f"{target} exists"
        if self.kind is CheckKind.NON_EMPTY:
            return f"{target} is non-empty"
        if self.kind is CheckKind.PARSEABLE:
            return f"{target} parses as {self.parser}"
        if self.kind is CheckKind.COLUMNS:
            return f"{target} has columns {', '.join(self.columns)}"
        if self.kind is CheckKind.MIN_ROWS:
            return f"{target} has at least {self.minimum} data rows"
        if self.kind is CheckKind.MIN_SIZE:
            return f"{target} is at least {self.minimum} bytes"
        if self.kind is CheckKind.SAMPLE_COUNT:
            return f"{target} covers {self.expected} samples"
        if self.kind is CheckKind.IDENTIFIERS:
            return f"{target} identifiers are all expected"
        if self.kind is CheckKind.NO_SENTINEL_COLUMN:
            return f"{target} column {self.key!r} contains real values"
        if self.kind is CheckKind.MATRIX_NOT_CONSTANT:
            return f"{target} is not a constant matrix"
        if self.kind is CheckKind.IDENTIFIES_FILE:
            return f"{target} is named after {self.values[0] if self.values else '?'}"
        return f"{target} satisfies its declared check"

    @property
    def failure_state(self) -> StageState:
        return (
            StageState.INCOMPLETE
            if self.kind in _INCOMPLETE_KINDS
            else StageState.INVALID
        )


# --- ergonomic constructors, so specs read as declarations ---------------

def exists(path: Path, description: str = "") -> Check:
    return Check(CheckKind.EXISTS, path=path, description=description)


def non_empty(path: Path, minimum: int = 1, description: str = "") -> Check:
    return Check(CheckKind.NON_EMPTY, path=path, minimum=minimum, description=description)


def parses(path: Path, parser: str, description: str = "") -> Check:
    return Check(CheckKind.PARSEABLE, path=path, parser=parser, description=description)


def has_columns(path: Path, columns: Sequence[str], dialect: str = "pipeline",
                description: str = "") -> Check:
    return Check(CheckKind.COLUMNS, path=path, columns=tuple(columns),
                 dialect=dialect, description=description)


def min_rows(path: Path, count: int, description: str = "") -> Check:
    return Check(CheckKind.MIN_ROWS, path=path, minimum=count, description=description)


def covers_samples(path: Path, count: int, key: str = "sample_id", description: str = "") -> Check:
    return Check(CheckKind.SAMPLE_COUNT, path=path, expected=count, key=key, description=description)


def identifiers_known(path: Path, expected: Sequence[str], key: str = "sample_id",
                      description: str = "") -> Check:
    return Check(CheckKind.IDENTIFIERS, path=path, values=tuple(expected), key=key,
                 description=description)


def has_real_values(path: Path, key: str, dialect: str = "pipeline",
                    description: str = "") -> Check:
    return Check(CheckKind.NO_SENTINEL_COLUMN, path=path, key=key,
                 dialect=dialect, description=description)


def not_constant_matrix(path: Path, description: str = "") -> Check:
    return Check(CheckKind.MATRIX_NOT_CONSTANT, path=path, description=description)


def identifies_file(path: Path, stem: str, description: str = "") -> Check:
    """The file name must carry the expected genome stem."""
    return Check(CheckKind.IDENTIFIES_FILE, path=path, values=(stem,),
                 description=description)


@dataclass(frozen=True)
class SiblingSpec:
    """Cross-sample degeneracy check.

    A per-sample task is validated in isolation, which cannot detect that
    every sample produced the *same* thing. This is the check that would
    have caught the bug where a leaked loop variable made all 65 per-sample
    VCFs byte-identical and the resulting core-SNP alignment had one
    distinct sequence.
    """

    root: Path
    pattern: str = "*"
    #: At least this many distinct contents must exist among the siblings.
    min_distinct: int = 2
    label: str = "output"
    description: str = ""


@dataclass(frozen=True)
class OutputSpec:
    """A stage's declared output contract."""

    stage: str
    checks: Tuple[Check, ...] = ()
    sibling: Optional[SiblingSpec] = None
    #: Free-form provenance recorded alongside the verdict.
    expectations: Mapping[str, Any] = field(default_factory=dict)

    def paths(self) -> Tuple[Path, ...]:
        seen: List[Path] = []
        for check in self.checks:
            if check.path is not None and check.path not in seen:
                seen.append(check.path)
        return tuple(seen)

    def describe(self) -> str:
        lines = [f"{self.stage}: {len(self.checks)} check(s)"]
        lines += [f"  - {c.describe()}" for c in self.checks]
        if self.sibling is not None:
            lines.append(
                f"  - {self.sibling.label}s under {self.sibling.pattern} "
                f"have >= {self.sibling.min_distinct} distinct contents"
            )
        return "\n".join(lines)


@dataclass(frozen=True)
class CheckResult:
    check: Check
    passed: bool
    detail: str = ""

    def to_row(self) -> Dict[str, Any]:
        return {
            "check": self.check.kind.value,
            "target": self.check.path.name if self.check.path else "",
            "passed": self.passed,
            "detail": self.detail,
        }


@dataclass(frozen=True)
class ValidationResult:
    """The verdict for one attempt's outputs."""

    stage: str
    results: Tuple[CheckResult, ...] = ()
    siblings: Optional[CheckResult] = None

    @property
    def failures(self) -> Tuple[CheckResult, ...]:
        return tuple(r for r in self.results if not r.passed)

    @property
    def ok(self) -> bool:
        return not self.failures and (self.siblings is None or self.siblings.passed)

    @property
    def state(self) -> StageState:
        """``SUCCEEDED`` when clean, else the state of the first failure.

        Using the *first* failure keeps the verdict stable and explainable:
        a missing file is reported as missing even if a later check also
        failed because of it.
        """
        if self.ok:
            return StageState.SUCCEEDED
        if self.failures:
            return self.failures[0].check.failure_state
        if self.siblings is not None and not self.siblings.passed:
            return StageState.INVALID
        return StageState.INVALID

    @property
    def detail(self) -> str:
        if self.ok:
            return f"{self.stage}: {len(self.results)} check(s) passed"
        parts = [f"{r.check.describe()}: {r.detail}" for r in self.failures]
        if self.siblings is not None and not self.siblings.passed:
            parts.append(f"sibling: {self.siblings.detail}")
        return "; ".join(parts)

    def to_json(self) -> str:
        payload = {
            "stage": self.stage,
            "ok": self.ok,
            "state": self.state.value,
            "checks": [r.to_row() for r in self.results],
        }
        if self.siblings is not None:
            payload["siblings"] = self.siblings.to_row()
        return json.dumps(payload, indent=2, sort_keys=True)


#: How to read a table. The distinction is not cosmetic.
#:
#: ``pipeline`` uses :func:`papipeline.io.tsv.read_tsv`, which skips ``#``
#: comment lines *anywhere* in the file. A Bakta annotation table puts its
#: header on a ``#``-prefixed line, so ``read_tsv`` would treat the real
#: header as a comment and promote the first data row to be the header -
#: silently reporting duplicate column names, or worse, a plausible wrong
#: answer. ``bakta`` uses the pipeline's own
#: :func:`papipeline.stages.annotation.parse_bakta_tsv`, which is the
#: reader that actually understands these files.
DIALECTS: Dict[str, Callable[[Path], List[Dict[str, Any]]]] = {}


def _read_pipeline_tsv(path: Path) -> List[Dict[str, Any]]:
    return read_tsv(path)


def _read_bakta_tsv(path: Path) -> List[Dict[str, Any]]:
    """Normalised Bakta feature rows: seqid, gene_id, gene_name, product, ..."""
    from ..stages.annotation import parse_bakta_tsv
    return parse_bakta_tsv(path.read_text(encoding="utf-8"))


def _read_bakta_raw(path: Path) -> List[Dict[str, Any]]:
    """Bakta rows with the *file's own* column names preserved.

    :func:`_read_bakta_tsv` normalises to the pipeline's schema, so a column
    check against it cannot tell a feature table from an inference table -
    both produce the same key set. This dialect keeps ``Gene``, ``Score``,
    ``Id`` and friends, which is what distinguishes them.
    """
    from ..stages.annotation import read_tsv_text
    with path.open(encoding="utf-8") as handle:
        return read_tsv_text(handle)


def _read_vcf_rows(path: Path) -> List[Dict[str, Any]]:
    """One dict per VCF record, keyed by the file's own column names.

    A VCF is not a TSV: its header is ``#CHROM``-prefixed and rows carry
    genotype columns, so reading it with the pipeline TSV reader promotes a
    data row to be the header.
    """
    fields: List[str] = []
    out: List[Dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                fields = cols
                continue
            if line.strip() and fields:
                out.append({fields[i]: (cols[i] if i < len(cols) else None)
                            for i in range(len(fields))})
    return out


DIALECTS["pipeline"] = _read_pipeline_tsv
DIALECTS["bakta"] = _read_bakta_tsv
DIALECTS["bakta_raw"] = _read_bakta_raw
DIALECTS["vcf"] = _read_vcf_rows


def _read_gff(path: Path) -> int:
    """A GFF is line-oriented text; verify the directive and count features."""
    n = 0
    saw_directive = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##gff-version"):
                saw_directive = True
            elif line.strip() and not line.startswith("#"):
                n += 1
    if not saw_directive:
        raise ValueError("no ##gff-version directive; not a GFF")
    return n


# --- parsers -------------------------------------------------------------

def _read_tsv(path: Path) -> List[Dict[str, Any]]:
    return read_tsv(path)


def _read_fasta(path: Path) -> int:
    return sum(1 for _ in read_fasta(path))


def _read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _read_vcf(path: Path) -> Tuple[int, List[str]]:
    """Minimal VCF reader: header lines plus one record per data line.

    Deliberately not a full VCF parser. It verifies the file is shaped like
    a VCF and has data, which is what a completion check needs; deeper
    interpretation belongs to the variant-calling stage, not here.
    """
    n_records = 0
    samples: List[str] = []
    saw_header = False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    raise ValueError(
                        f"#CHROM line has {len(cols)} field(s); a VCF needs at least 8")
                samples = cols[9:] if len(cols) > 9 else []
                saw_header = True
                continue
            if line.strip():
                n_records += 1
    if not saw_header:
        raise ValueError("no #CHROM column header line; not a VCF")
    return n_records, samples


def _read_newick(path: Path) -> int:
    text = path.read_text(encoding="utf-8").strip()
    if not text.endswith((";", "\n")) and ";" not in text:
        raise ValueError("not a Newick string (no terminating ';')")
    return text.rstrip().rstrip(";").count(":")


PARSERS: Dict[str, Callable[[Path], Any]] = {
    "tsv": _read_tsv,
    "bakta": _read_bakta_tsv,
    "gff": _read_gff,
    "fasta": _read_fasta,
    "json": _read_json,
    "vcf": _read_vcf,
    "newick": _read_newick,
}


# --- check execution -----------------------------------------------------

def _run_check(check: Check) -> CheckResult:
    path = check.path
    assert path is not None

    if not path.exists():
        return CheckResult(check, False, "file does not exist")

    if check.kind is CheckKind.EXISTS:
        return CheckResult(check, True, "present")

    if check.kind is CheckKind.NON_EMPTY:
        size = path.stat().st_size
        need = check.minimum if check.minimum is not None else 1
        if size < need:
            return CheckResult(check, False, f"{size} bytes, expected >= {need}")
        return CheckResult(check, True, f"{size} bytes")

    if check.kind is CheckKind.MIN_SIZE:
        size = path.stat().st_size
        if size < (check.minimum or 0):
            return CheckResult(check, False, f"{size} bytes, expected >= {check.minimum}")
        return CheckResult(check, True, f"{size} bytes")

    def rows() -> List[Dict[str, Any]]:
        reader = DIALECTS.get(check.dialect)
        if reader is None:
            raise ValueError(f"unknown table dialect {check.dialect!r}")
        return reader(path)

    # Every remaining check needs to read the file.
    if check.kind is CheckKind.PARSEABLE:
        parser = PARSERS.get(check.parser or "")
        if parser is None:
            return CheckResult(check, False, f"unknown parser {check.parser!r}")
        try:
            parser(path)
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        return CheckResult(check, True, f"parses as {check.parser}")

    if check.kind is CheckKind.COLUMNS:
        try:
            data = rows()
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        if not data:
            return CheckResult(check, False, "no data rows")
        missing = [c for c in check.columns if c not in data[0]]
        if missing:
            return CheckResult(check, False, f"missing columns: {', '.join(missing)}")
        return CheckResult(check, True, f"{len(check.columns)} required column(s) present")

    if check.kind is CheckKind.MIN_ROWS:
        try:
            n = len(rows())
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        if n < (check.minimum or 0):
            return CheckResult(check, False, f"{n} data rows, expected >= {check.minimum}")
        return CheckResult(check, True, f"{n} data rows")

    if check.kind is CheckKind.NO_SENTINEL_COLUMN:
        try:
            data = rows()
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        if not data:
            return CheckResult(check, False, "no data rows")
        real = sum(1 for r in data if r.get(check.key) not in (None,) and
                   str(r.get(check.key)).strip() not in _SENTINEL_TEXT)
        if real == 0:
            return CheckResult(check, False, f"column {check.key!r} is empty or all sentinel")
        return CheckResult(check, True, f"{real}/{len(data)} rows have a real {check.key!r}")

    if check.kind is CheckKind.SAMPLE_COUNT:
        try:
            data = rows()
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        ids = {str(r.get(check.key)) for r in data if r.get(check.key) is not None}
        if check.expected is not None and len(ids) != check.expected:
            return CheckResult(
                check, False, f"{len(ids)} distinct {check.key} values, expected {check.expected}"
            )
        return CheckResult(check, True, f"{len(ids)} distinct {check.key} values")

    if check.kind is CheckKind.IDENTIFIERS:
        try:
            data = rows()
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        allowed = set(check.values)
        seen = {str(r.get(check.key)) for r in data if r.get(check.key) is not None}
        unknown = sorted(seen - allowed)
        if unknown:
            preview = ", ".join(unknown[:5])
            more = f" (+{len(unknown) - 5} more)" if len(unknown) > 5 else ""
            return CheckResult(check, False, f"unexpected {check.key}: {preview}{more}")
        return CheckResult(check, True, f"{len(seen)} identifier(s), all expected")

    if check.kind is CheckKind.IDENTIFIES_FILE:
        expected = check.values[0] if check.values else ""
        stem = path.name
        for suffix in (".gff3", ".gff", ".tsv", ".gz", ".faa", ".ffn", ".fna"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        if expected and expected not in stem:
            return CheckResult(
                check, False,
                f"file is named for {stem!r}, expected it to name {expected!r}")
        return CheckResult(check, True, f"file names {expected}")

    if check.kind is CheckKind.MATRIX_NOT_CONSTANT:
        try:
            data = rows()
        except (DataContractError, PipelineError) as exc:
            return CheckResult(check, False, f"unreadable: {exc}")
        except Exception as exc:  # noqa: BLE001
            return CheckResult(check, False, f"{type(exc).__name__}: {exc}")
        numeric: List[float] = []
        for row in data:
            for key, value in row.items():
                if key in (check.key, "sample", "sample_id"):
                    continue
                try:
                    numeric.append(float(value))
                except (TypeError, ValueError):
                    continue
        if not numeric:
            return CheckResult(check, False, "no numeric matrix values found")
        distinct = {round(v, 9) for v in numeric}
        if len(distinct) == 1:
            only = next(iter(distinct))
            return CheckResult(
                check, False,
                f"every matrix value is {only}; a constant matrix carries no information",
            )
        return CheckResult(check, True, f"{len(distinct)} distinct value(s)")

    return CheckResult(check, False, f"unhandled check kind {check.kind}")


def _run_sibling(spec: SiblingSpec) -> CheckResult:
    probe = Check(CheckKind.NON_EMPTY, path=spec.root / "x")
    files = sorted(p for p in spec.root.glob(spec.pattern) if p.is_file())
    if not files:
        return CheckResult(probe, False,
                           f"no {spec.label}s found under {spec.root}/{spec.pattern}")
    digests = {}
    for path in files:
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                h.update(chunk)
        digests.setdefault(h.hexdigest(), []).append(path.name)
    if len(digests) < spec.min_distinct:
        groups = sorted(digests.values(), key=len, reverse=True)
        worst = groups[0]
        return CheckResult(
            probe, False,
            f"all {len(files)} {spec.label}s are byte-identical or there is only "
            f"{len(digests)} distinct content(s); {len(worst)} share one "
            f"(e.g. {', '.join(worst[:3])})",
        )
    return CheckResult(probe, True,
                       f"{len(files)} {spec.label}s, {len(digests)} distinct contents")


def validate(spec: OutputSpec) -> ValidationResult:
    """Run a stage's declared checks and return the verdict.

    Read-only: nothing here writes, repairs, or re-runs. The verdict is
    recorded by the caller alongside the execution state.
    """
    results = tuple(_run_check(c) for c in spec.checks)
    siblings = _run_sibling(spec.sibling) if spec.sibling is not None else None
    return ValidationResult(stage=spec.stage, results=results, siblings=siblings)
