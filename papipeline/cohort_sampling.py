"""Sampling isolates for the allele-frequency measurement, stratified not arbitrary.

The measurement asks what allele frequencies look like across a cohort, so the
sample has to be one where a frequency means something. Taking the first 50 by
directory order measures whichever clone or hospital the accession numbering
happens to start in: the first four real isolates turned out to span 18-40
contigs and 6.4-6.9 Mb, which is already a spread, but nothing guarantees that
and a homogeneous cluster would collapse the spectrum to two spikes and make any
threshold derived from it an artefact of the sample.

Stratification here is by `Location`, the field with the most real variation (46
distinct values among the 835 isolates that carry an assembly), and then a
round-robin walk so small strata contribute rather than being crowded out by
`USA: Texas` (240 isolates) and `USA: Nashville, TN` (137). A plain random sample
would give `USA: Texas` roughly 29% of 50 isolates; round-robin gives every
stratum a share, which is what keeps a rare-allele class present at all.

Isolates with no location are not discarded - 35 of them - but they form their
own stratum rather than being dropped, so a run restricted to isolates that
happened to record a country is not silently assumed.

The choice matters for the threshold: a spectrum measured on one clone cannot
show a site that is rare cohort-wide, and `max_maf` is exactly the quantity that
needs one.
"""

from __future__ import annotations

import collections
import csv
from pathlib import Path
from typing import Dict, List, Optional, Sequence


def load_located(pdc_path: Path) -> List[Dict[str, str]]:
    """Isolates that carry an assembly accession, with their metadata."""
    rows: List[Dict[str, str]] = []
    with Path(pdc_path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("Assembly") or "").strip():
                rows.append({k: (v or "").strip() for k, v in row.items() if k})
    return rows


def stratify(rows: Sequence[Dict[str, str]], key: str = "Location") -> Dict[str, List[Dict[str, str]]]:
    """Group rows by one metadata field, empties forming their own stratum.

    An empty value is a stratum, not a filter. Dropping the 35 isolates with no
    recorded location would quietly bias the sample toward isolates whose labs
    filled the field in.
    """
    groups: Dict[str, List[Dict[str, str]]] = collections.defaultdict(list)
    for row in rows:
        groups[row.get(key, "") or "(unrecorded)"].append(row)
    return dict(groups)


def sample_stratified(
    rows: Sequence[Dict[str, str]],
    n: int,
    *,
    key: str = "Location",
    seed: Optional[int] = None,
) -> List[Dict[str, str]]:
    """Take ``n`` rows, round-robin across strata.

    Deterministic when ``seed`` is None, so the measured spectrum is
    reproducible from this code rather than from a run nobody can repeat. Within
    a stratum rows are taken in file order, which is stable.

    Raises:
        ValueError: ``n`` exceeds the number of available isolates.
    """
    if n > len(rows):
        raise ValueError(f"asked for {n} isolates but only {len(rows)} carry an assembly")
    groups = stratify(rows, key)
    # Largest strata first, so the walk starts with the ones that would otherwise
    # be truncated; ties broken by name so the order is stable across runs.
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    chosen: List[Dict[str, str]] = []
    index = 0
    while len(chosen) < n:
        progressed = False
        for _name, members in ordered:
            if index < len(members):
                chosen.append(members[index])
                progressed = True
                if len(chosen) == n:
                    break
        if not progressed:
            break
        index += 1
    return chosen
