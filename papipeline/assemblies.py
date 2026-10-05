"""Locating a sample's assembly file, and refusing when it cannot.

The pipeline's manifest is authoritative for *which* samples are in the cohort.
This module is only about finding the file for one of them, and it does that in
three ordered steps:

1. **An explicit ``assembly_path`` is used directly**, with no pattern matching.
   This is how the bounded REAL smoke run points at ``db/smoke_genomes/`` - the
   manifest says where the file is and nothing else is consulted.
2. **Otherwise the sample's accession directory is searched** for a file whose
   name carries the accession. Real assemblies are
   ``GCA_000710625.1_MRSN18971scaf_genomic.fna`` inside ``data/GCA_000710625.1/``,
   so requiring an exact ``{sample_id}.fna`` cannot see them - that is the whole
   reason this function exists. The TEST layout, where files are named exactly
   by sample id, is served by the same rule.
3. **Otherwise refuse.** No repository-wide recursive scan: a miss that widens
   into a search can attribute an unrelated file to this sample, and that
   failure mode belongs to ``pilot/cohort.py`` and to nothing else.

A pure function with no config coupling, so the stage, the Bakta adapter and the
workflow can all call it with what they already hold.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from .errors import PipelineError

if TYPE_CHECKING:  # pragma: no cover
    from .models import Sample

#: An assembly accession, as NCBI writes it. The version is part of the
#: identifier - `GCA_000710625.1` and `GCA_000710625.2` are different assemblies -
#: so the pattern requires it rather than matching a bare prefix.
ACCESSION = re.compile(r"GCA_\d+\.\d+")

#: Filename extensions an assembly may carry, longest first so `.fna.gz` is
#: tried before `.fa`.
SUFFIXES = (".fna.gz", ".fa.gz", ".fna", ".fasta", ".fa")


def _candidates(directory: Path, sample_id: str) -> list:
    """Files in ``directory`` that plausibly are this sample's assembly.

    Two accepted shapes, in preference order:

    * named exactly by the sample id - the TEST layout, and the historical
      convention;
    * named with the accession somewhere in the name - the real layout, where a
      file is ``<accession>_<lab strain>_scaf_genomic.fna``.

    Everything matching is returned, not just the first, so the caller can
    refuse on ambiguity instead of silently taking a match.
    """
    directory = Path(directory)
    if not directory.is_dir():
        return []

    found = []
    for suffix in SUFFIXES:
        exact = directory / f"{sample_id}{suffix}"
        if exact.is_file():
            found.append(exact)
    if found:
        return found

    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        if not ACCESSION.search(path.name):
            continue
        if path.suffix.lower() in (".gz",) and not any(
            path.name.lower().endswith(s) for s in SUFFIXES
        ):
            continue
        found.append(path)
    return found


def locate_assembly(sample: "Sample", data_root: Path) -> Path:
    """Return the assembly file for ``sample``, or refuse.

    Args:
        sample: The sample to locate. Its ``sample_id`` is expected to be a
            versioned assembly accession, which is also the genome directory
            name.
        data_root: Directory the cohort's genomes live under.

    Returns:
        The single assembly file for this sample.

    Raises:
        PipelineError: No assembly was found, or more than one matched. Both are
            refusals on purpose - a run that guessed here would analyse a
            sequence the manifest never named, and every coordinate downstream
            would be about the wrong genome.
    """
    data_root = Path(data_root)
    sample_id = sample.sample_id

    # Rule 1 - the manifest said where it is. No pattern matching.
    declared = getattr(sample, "assembly_path", None)
    if declared:
        candidate = Path(declared)
        if candidate.is_absolute():
            if candidate.is_file():
                return candidate
            raise PipelineError(
                f"Assembly for {sample_id} is declared at a path that does not exist",
                sample_id=sample_id,
                assembly_path=str(declared),
                resolved=str(candidate),
            )
        # A bare filename may be relative to the data root, or to the sample's
        # own accession directory - the manifest usually records a name, not a
        # full path, and both are common ways to write one.
        for parent in (data_root, data_root / sample_id):
            probe = parent / candidate.name
            if probe.is_file():
                return probe
        raise PipelineError(
            f"Assembly for {sample_id} is declared as {declared!r} but no such "
            f"file exists under {data_root} or {data_root / sample_id}",
            sample_id=sample_id,
            assembly_path=str(declared),
            searched=str(data_root),
        )

    # Rule 2 - named directories only, never a scan. The cohort is keyed on
    # `Isolate` (see `discover_pdc_manifest`), so `sample_id` is a PDC isolate
    # name and the real directory is named by the *assembly accession*, which is
    # carried on the sample as an attribute. Both keys are tried: searching only
    # `data/{sample_id}/` works for the TEST layout and silently finds nothing
    # for every real sample, which refuses as if the genome were missing.
    for key in (sample_id, getattr(sample, "assembly", None)):
        if not key:
            continue
        for directory in (data_root / key, data_root / "genomes"):
            found = _candidates(directory, key)
            if found:
                break
        if found:
            break
    else:
        found = []

    # Rule 3 - refuse. Never widen into a repository-wide search.
    if not found:
        raise PipelineError(
            f"No assembly found for {sample_id} under {data_root / sample_id}. "
            "An isolate with no downloaded assembly is a normal outcome, not an "
            "error to work around: a stage that needs the sequence must refuse "
            "this sample rather than substitute another one.",
            sample_id=sample_id,
            searched=str(data_root / sample_id),
            genomes_dir=str(data_root),
        )

    if len(found) > 1:
        raise PipelineError(
            f"Assembly for {sample_id} is ambiguous: {len(found)} files match "
            f"in {data_root / sample_id}. Refusing rather than choosing one, "
            "because the wrong choice analyses a different sequence silently. "
            "Set an explicit assembly_path on the sample to disambiguate.",
            sample_id=sample_id,
            candidates="; ".join(str(p) for p in found),
            searched=str(data_root / sample_id),
        )

    return found[0]
