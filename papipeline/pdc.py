"""Row reading for ``PDC_essential.tsv``, shared by the pilot and the pipeline.

This module exists to keep one implementation of "read the PDC table" in a
neutral place. `pilot.cohort` reads the same file with its own parser, and the
natural next step - letting the pipeline build its cohort from PDC - would
otherwise mean either a second parser that can drift or an import from the
pipeline core into a pilot tool, which is backwards: the pilot is the part most
likely to be deleted.

One deliberate difference from the pilot's own index, and it is the whole point
of this module. `pilot.cohort.load_pdc_index` is keyed on ``Assembly`` and skips
any row without one, so it holds 835 entries for a 967-row file. That is correct
for its job - it attaches metadata to *selected assemblies* - and wrong for
building a cohort, where **every isolate is a member** and an isolate with no
assembly is an ordinary member with no sequence. So the key is a parameter here,
and no row is skipped for lacking the key column; a row whose *key* is empty is a
data error, not something to drop.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Mapping, Optional

from .errors import PipelineError

MISSING = "MISSING"

#: Column naming the cohort key. ``Isolate`` rather than ``Assembly`` because
#: cohort membership is a property of an isolate, not of an accession.
COHORT_KEY = "Isolate"

REQUIRED_COLUMNS = (COHORT_KEY, "Assembly", "BioSample")


def clean_value(value: Optional[str]) -> str:
    """Normalise a cell: strip it, and mark empty as ``MISSING``."""
    text = (value or "").strip()
    return text if text else MISSING


def load_pdc_index(pdc_path: Path, key: str = COHORT_KEY) -> Dict[str, Dict[str, str]]:
    """Index ``PDC_essential.tsv`` by ``key``, keeping every row.

    Args:
        pdc_path: The TSV to read.
        key: Column to index on. Defaults to ``Isolate``.

    Returns:
        A mapping of key value to the row's cleaned cells, in file order.

    Raises:
        PipelineError: The file is missing, empty, lacks the key column, or
            repeats a key. A repeated key is reported rather than resolved: the
            key is the join key for every downstream stage, so a duplicate means
            the data is wrong and guessing which row wins would be worse.
    """
    pdc_path = Path(pdc_path)
    if not pdc_path.exists():
        raise PipelineError("PDC_essential.tsv not found", path=str(pdc_path))
    if pdc_path.stat().st_size == 0:
        raise PipelineError(
            "PDC_essential.tsv is empty; no metadata can be read",
            path=str(pdc_path),
        )

    index: Dict[str, Dict[str, str]] = {}
    duplicates: List[str] = []
    with pdc_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if key not in fields:
            raise PipelineError(
                f"PDC_essential.tsv has no {key!r} column",
                path=str(pdc_path),
                found=",".join(fields),
            )
        for row in reader:
            value = (row.get(key) or "").strip()
            if not value:
                raise PipelineError(
                    f"PDC_essential.tsv has a row with an empty {key!r}",
                    path=str(pdc_path),
                )
            if value in index:
                duplicates.append(value)
                continue
            index[value] = {
                name: clean_value(cell) if isinstance(cell, str) else cell
                for name, cell in row.items()
                if name
            }

    if duplicates:
        raise PipelineError(
            f"PDC_essential.tsv contains duplicate {key} values",
            n_duplicates=len(duplicates),
            examples=",".join(sorted(set(duplicates))[:5]),
        )
    return index


def row_value(row: Mapping[str, str], column: str) -> Optional[str]:
    """A cleaned cell, or ``None`` when the column is absent or ``MISSING``.

    A missing assembly accession is ``None`` rather than the ``MISSING`` marker:
    for cohort construction "no assembly" is a fact to carry, and substituting a
    placeholder string would put it into a path or a comparison.
    """
    value = row.get(column)
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == MISSING:
        return None
    return text
