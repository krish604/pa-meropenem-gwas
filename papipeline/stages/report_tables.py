"""The run's own output tables, read back from disk at their contracted paths.

**Why this module exists.** Until now the report was assembled purely from
in-memory aggregates the orchestrator handed it (`run.py`'s
`_build_report_context`), so it could not tell "this stage found nothing" from
"this stage never ran". Worse, nothing in the repository read the stage tables
back: the round-12 seam analysis counted **24 output artefacts with zero
in-repo readers** (`pa-artifacts/round12/seam-matrix.md` §e), so every number in
a report was a second-hand copy of something the run had already written down.

This module is the reader that closes that gap. It opens each stage's principal
table **at the path the repository already declares** - `contracts.table_path`
and `contracts.internal_table_path` for the stage contract, and
`similarity.units_sidecar_path` for the units sidecar - and records, for each
one, whether it was there and why not when it was not.

**Paths are never invented here.** Every location below comes from one of:

* `papipeline/execution/contracts.py` `STAGE_TABLES` / `INTERNAL_TABLES`,
  which is the single source the writer and the contract validator both use;
* `workflow/Snakefile`'s `OUT_*` constants, for the three pan-genome gene
  tables and the alignment summary, which the Snakefile declares and
  `contracts.py` does not.

A table this module cannot find is reported as **not produced, with the reason**,
never as a blank section and never as a zero-row table. Those are different
claims: a zero-row table says "the stage ran and found nothing", which for a
cohort of ten assemblies is usually a bug rather than a result.

**Zero rows is treated as not produced, on purpose.** A header-only table
cannot distinguish "screened ten isolates, found nothing" from "never opened the
file", and for every table read here the run's own contract
(`contracts.stage_spec`) requires `min_rows(table, 1)`. A header-only file has
already failed its contract, so reporting it as a result would report a
contract violation as science.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..logging_utils import get_logger
from ..models import RunMode

LOGGER = get_logger("stages.report_tables")

#: The exact word a section uses when its table is not there.
#:
#: Named because a reader may rely on it, and because it is the alternative to
#: the two failure modes this module exists to remove: a blank section, and a
#: zero-row table presented as a finding.
NOT_PRODUCED = "not produced"

Row = Dict[str, Any]


@dataclass(frozen=True)
class TableRead:
    """One stage table, as found on disk - or as found *not* to be.

    ``present`` is the only thing a caller should branch on. ``rows`` is empty
    whenever ``present`` is ``False``, and ``reason`` is the sentence that says
    why, so the two cannot drift: a caller that renders ``rows`` without
    checking ``present`` renders nothing, which is still better than a
    fabricated finding but is not what this module is for.
    """

    name: str
    path: Path
    rows: Tuple[Row, ...] = ()
    present: bool = False
    reason: str = ""

    def __bool__(self) -> bool:
        return self.present

    @property
    def n_rows(self) -> int:
        return len(self.rows)

    @property
    def absent_sentence(self) -> str:
        """The line a section prints instead of a table."""
        return f"{NOT_PRODUCED}: {self.reason}"

    def count(self, column: str) -> Counter:
        """Value counts for ``column``, skipping absent cells.

        Absent cells are skipped rather than counted as their own value: a row
        with no `ST` says the stage had nothing to record for that isolate, and
        counting "" as a category would report an isolate with no sequence type
        as a sequence type called "".
        """
        counter: Counter = Counter()
        for row in self.rows:
            value = row.get(column)
            if value is None or str(value).strip() == "":
                continue
            counter[str(value).strip()] += 1
        return counter

    def distinct(self, column: str) -> List[str]:
        return sorted(self.count(column))


def read_table(
    name: str,
    path: Path,
    *,
    required_columns: Sequence[str] = (),
) -> TableRead:
    """Read one table, or say precisely why it could not be read.

    Never raises for an absent or empty file: a report that crashed because a
    stage did not run would be worse than a report that says the stage did not
    run. A file that exists but violates its contract *is* reported as not
    produced, with the violation named - that is a real problem in the run, and
    naming it is more useful than propagating an exception out of the report.
    """
    from ..errors import PipelineError
    from ..io.tsv import read_tsv

    path = Path(path)
    if not path.exists():
        return TableRead(name, path, reason=f"no file at {path}")
    if path.stat().st_size == 0:
        return TableRead(name, path, reason=f"{path} is present but empty")
    if _is_header_only(path):
        # Checked before `read_tsv` rather than caught from it.
        # `read_tsv` refuses a header-only file with "contains a header but no
        # data rows", which is true and useless: it does not say that every
        # contract in this module requires at least one row, so a reader cannot
        # tell a stage that found nothing from a file nothing was written into.
        return TableRead(
            name,
            path,
            reason=(
                f"{path} carries a header and no rows. That is not a result - it "
                f"is what the file looks like when nothing was written into it, "
                f"and {name}'s output contract requires at least one row "
                f"(`contracts.stage_spec` adds `min_rows(table, 1)`)"
            ),
        )

    try:
        rows = read_tsv(path, required_columns=list(required_columns))
    except PipelineError as exc:
        return TableRead(name, path, reason=f"{path} is unreadable: {exc}")

    if not rows:
        return TableRead(
            name,
            path,
            reason=(
                f"{path} carries a header and no rows. That is not a result - it "
                f"is what the file looks like when nothing was written into it, "
                f"and {name}'s output contract requires at least one row"
            ),
        )
    return TableRead(name, path, rows=tuple(rows), present=True)


def _is_header_only(path: Path) -> bool:
    """Whether ``path`` holds a header row and no data rows.

    Counted here, with the same comment rule `read_tsv` applies, so the two
    agree on what a data row is. Stops at the second content line, so a
    500,000-row variant table costs one line of reading to answer ``False``.
    """
    content = 0
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.strip():
                content += 1
                if content > 1:
                    return False
    return content == 1


def read_json_sidecar(name: str, path: Path) -> TableRead:
    """Read a JSON sidecar as a one-row table, or say why it is not there.

    The units sidecar is the reason this exists in the shape it does: a bare
    distance matrix does not say what its numbers measure, and the sidecar is the
    only place that answer is written down.
    """
    path = Path(path)
    if not path.exists():
        return TableRead(name, path, reason=f"no sidecar at {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return TableRead(name, path, reason=f"{path} is unreadable: {exc}")
    if not isinstance(payload, Mapping) or not payload:
        return TableRead(name, path, reason=f"{path} carries no keys")
    return TableRead(
        name,
        path,
        rows=({str(k): v for k, v in payload.items()},),
        present=True,
    )


#: Path-only artefacts the Snakefile declares that `contracts.py` does not.
#:
#: `OUT_PANGENOME_GPA`, `OUT_PANGENOME_CORE` and `OUT_PANGENOME_ACC` are
#: `workflow/Snakefile:265-267`; the alignment summary is written at
#: `papipeline/run.py:1250` as `stage_dir / "10_alignment_summary.tsv"` and is
#: conditional on `alignment_info` being non-empty. Each is recorded here with
#: the file that declares it, so a reader of this module can tell a declared
#: path from an invented one.
SNAKEFILE_ONLY_TABLES: Tuple[Tuple[str, str, str], ...] = (
    ("gene_presence_absence", "gene_presence_absence.tsv", "workflow/Snakefile:265"),
    ("core_genes", "core_genes.tsv", "workflow/Snakefile:266"),
    ("accessory_genes", "accessory_genes.tsv", "workflow/Snakefile:267"),
    ("alignment_summary", "10_alignment_summary.tsv", "papipeline/run.py:1250"),
)

#: Stage name -> the `contracts.STAGE_TABLES` key holding its principal table.
#:
#: Keys are the stage names `run.STAGE_ORDER` uses, so a stage that is dispatched
#: but has no table here is visible as a missing entry rather than as a silently
#: absent section.
STAGE_TABLE_KEYS: Tuple[Tuple[str, str], ...] = (
    ("validation", "validation"),
    ("annotation", "annotation"),
    ("mlst", "mlst"),
    ("amr", "amr"),
    ("virulence", "virulence"),
    ("variants", "variants"),
    ("cohort_variants", "cohort_variants"),
    ("pangenome", "pangenome"),
    ("recombination", "recombination"),
    ("phylogeny", "phylogeny"),
    ("similarity", "similarity"),
    ("phenotype", "phenotype"),
    ("gwas", "gwas"),
    ("convergence", "convergence"),
    ("cooccurrence", "cooccurrence"),
    ("reporting", "reporting"),
)

#: Folded steps whose table is written, consumed and reported even though it is
#: not a stage's principal output (`spec.md:351`).
INTERNAL_TABLE_KEYS: Tuple[Tuple[str, str], ...] = (
    ("structural_variants", "structural_variants"),
    ("regulators", "regulators"),
    ("mechanisms", "mechanisms"),
    ("master_table", "master_table"),
)


@dataclass(frozen=True)
class RunTables:
    """Every table this run wrote, and the reason for each one it did not.

    ``tables`` is keyed on the *stage* name where one exists, so a caller asks
    for ``run.amr`` and gets the AMR table whether it came from
    `contracts.STAGE_TABLES` or `INTERNAL_TABLES`.
    """

    stage_dir: Path
    tables: Mapping[str, TableRead] = field(default_factory=dict)
    similarity_units: Optional[TableRead] = None

    def __getitem__(self, key: str) -> TableRead:
        # An explicit `is None` test, not `or`. `TableRead.__bool__` is
        # `present`, so `self.tables.get(key) or _absent(key)` would report
        # "no output is declared for X" for every table that is declared and
        # simply absent from disk - replacing the real reason (the file is not
        # there) with a false one (X was never declared).
        found = self.tables.get(key)
        return found if found is not None else _absent(key)

    def get(self, key: str) -> Optional[TableRead]:
        return self.tables.get(key)

    @property
    def n_present(self) -> int:
        return sum(1 for t in self.tables.values() if t.present)

    @property
    def n_absent(self) -> int:
        return sum(1 for t in self.tables.values() if not t.present)


def _absent(key: str) -> TableRead:
    return TableRead(
        key,
        Path(f"<{key}: not a table this pipeline declares>"),
        reason=(
            f"no output is declared for {key!r} in "
            f"papipeline/execution/contracts.py, so there is no path to read it "
            f"from. Add it to STAGE_TABLES or INTERNAL_TABLES rather than "
            f"inventing one here."
        ),
    )


def load_run_tables(
    config: PipelineConfig,
    mode: RunMode,
    *,
    stage_dir: Optional[Path] = None,
) -> RunTables:
    """Read every declared output table of one run, from disk.

    Args:
        config: The run's config. Supplies `intermediate_root(mode)`.
        mode: The mode whose results tree is read.
        stage_dir: Override for `<intermediate>/stages`. Only tests and the
            real-data render pass one; a pipeline run does not, because
            `intermediate_root` already puts every run in its own tree.

    Returns:
        A :class:`RunTables`. Never raises: a missing table is recorded, not
        raised, so one absent stage cannot take the whole report down.
    """
    from ..execution.contracts import (
        INTERNAL_TABLES,
        internal_table_path,
        required_columns,
        table_path,
    )

    root = Path(stage_dir) if stage_dir is not None else (
        Path(config.intermediate_root(mode)) / "stages"
    )

    tables: Dict[str, TableRead] = {}
    for stage, key in STAGE_TABLE_KEYS:
        tables[stage] = read_table(
            stage,
            table_path(root, key),
            required_columns=required_columns(key),
        )
    for stage, key in INTERNAL_TABLE_KEYS:
        _filename, declared = INTERNAL_TABLES[key]
        tables[stage] = read_table(
            stage,
            internal_table_path(root, key),
            required_columns=declared,
        )
    for stage, filename, declared_by in SNAKEFILE_ONLY_TABLES:
        tables[stage] = read_table(stage, root / filename)
        if not tables[stage].present:
            LOGGER.debug(
                "reporting: %s not produced (%s); path declared by %s",
                stage, tables[stage].reason, declared_by,
            )

    from .similarity import units_sidecar_path

    units = read_json_sidecar(
        "similarity_units",
        units_sidecar_path(table_path(root, "similarity")),
    )
    return RunTables(stage_dir=root, tables=tables, similarity_units=units)


__all__ = [
    "INTERNAL_TABLE_KEYS",
    "NOT_PRODUCED",
    "Row",
    "RunTables",
    "SNAKEFILE_ONLY_TABLES",
    "STAGE_TABLE_KEYS",
    "TableRead",
    "load_run_tables",
    "read_json_sidecar",
    "read_table",
]