"""The reference genome, verified before anything is called against it.

spec.md D5: the run checks the reference file it finds against the declared
length and checksum, because a truncated or substituted download shifts every
variant coordinate silently. Nothing downstream of a wrong reference is wrong in
a way you can see - the numbers come out, they are just about a different
genome - so the check has to happen before the first coordinate is produced, not
in the report.

This module was specified in D5 and unimplemented: nothing in `papipeline/`
read `science.yaml`'s `reference.md5`, so the field recorded an intent nobody
enforced. `verify_reference_file` is the enforcement.

The split here is deliberate. Verification is a property of *this file on this
machine*, not of the configuration, so it takes explicit expectations rather
than reaching into config, and it can be tested against a deliberately corrupted
file without a real genome anywhere in sight. `verify_reference` binds it to the
loaded configuration and is what a run calls.

Two states, because a missing checksum is not the same problem as a wrong one:

* **unpinned** (`md5: UNPINNED`, the deliberate placeholder) - verified on
  length alone, and reported as unverified rather than claimed as pinned. A
  length check still catches truncation, which is the common failure.
* **pinned** - length and checksum must both match, or the run stops.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping, Optional

from .errors import DataContractError
from .logging_utils import get_logger

LOGGER = get_logger("reference")

#: The placeholder that means "not pinned yet", rather than a real digest.
UNPINNED = "UNPINNED"

#: Chunk size for the digest. Large enough to keep syscall overhead irrelevant
#: on a 6 MB file, small enough not to matter.
_CHUNK = 1024 * 1024


def _sequence_length(path: Path) -> int:
    """Total bases in a FASTA file, counting every record.

    Not the file size: headers, newlines and line wrapping all differ between
    a freshly downloaded file and a reformatted one, and none of them change the
    genome. The declared length is a scientific fact about the assembly, so it
    is checked against the bases.
    """
    total = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def file_md5(path: Path) -> str:
    """MD5 of the file's bytes, streamed."""
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_reference_file(
    path: Path,
    *,
    expected_length: Optional[int] = None,
    expected_md5: Optional[str] = None,
    accession: str = "reference",
) -> Mapping[str, Any]:
    """Check one reference file against its declared identity.

    Args:
        path: The FASTA file to check.
        expected_length: Declared total bases, or None to skip.
        expected_md5: Declared digest, ``UNPINNED``, or None to skip.
        accession: Used in messages, so the operator knows which genome failed.

    Returns:
        A record of what was checked, for the run manifest.

    Raises:
        DataContractError: The file is missing, its length disagrees with the
            declaration, or a pinned checksum does not match. All three stop
            the run, because all three make every downstream coordinate
            untrustworthy and none of them is visible in the output.
    """
    path = Path(path)
    if not path.exists():
        raise DataContractError(
            f"Reference genome {accession} is not present at {path}. Nothing "
            "downstream can be trusted without it: a run that produced no "
            "coordinates and a run that produced wrong ones look identical from "
            "the output alone.",
            accession=accession,
            path=str(path),
        )

    observed_length = _sequence_length(path)
    if expected_length is not None and observed_length != expected_length:
        raise DataContractError(
            f"Reference {accession} is {observed_length} bases; "
            f"{expected_length} were declared. A truncated or substituted "
            "reference shifts every variant coordinate, so the run stops here "
            "rather than reporting positions against the wrong genome.",
            accession=accession,
            expected_length=expected_length,
            observed_length=observed_length,
            path=str(path),
        )

    pinned = bool(expected_md5) and expected_md5 != UNPINNED
    if not pinned:
        # Deliberate, not a failure: the declared checksum is a placeholder
        # until someone derives it from the real file. Length is still checked,
        # so the common failure - a partial download - is caught.
        LOGGER.warning(
            "Reference %s is UNPINNED: length was verified (%d bases) but the "
            "checksum was not, because none is declared. Coordinates from this "
            "run are unverified against the intended assembly.",
            accession, observed_length,
        )
        return {
            "accession": accession,
            "path": str(path),
            "observed_length": observed_length,
            "expected_length": expected_length,
            "pinned": False,
            "observed_md5": None,
        }

    observed_md5 = file_md5(path)
    if observed_md5 != expected_md5:
        raise DataContractError(
            f"Reference {accession} has md5 {observed_md5}; {expected_md5} was "
            "declared. The file is not the genome the science was written "
            "against, and every coordinate it produces would be about a "
            "different sequence.",
            accession=accession,
            expected_md5=expected_md5,
            observed_md5=observed_md5,
            path=str(path),
        )

    return {
        "accession": accession,
        "path": str(path),
        "observed_length": observed_length,
        "expected_length": expected_length,
        "pinned": True,
        "observed_md5": observed_md5,
    }


def verify_reference(config: Any, mode: Any = None) -> Optional[Mapping[str, Any]]:
    """Verify the reference this run is configured to use.

    Args:
        config: A loaded :class:`~papipeline.config.loader.PipelineConfig`.
        mode: Unused today, accepted so the call site can state the mode and
            stay correct if verification later becomes mode-dependent.

    Returns:
        The record from :func:`verify_reference_file`, or None when the
        configuration declares no reference at all - a pipeline with no
        reference has nothing to verify, which is not the same as a broken one.
    """
    # The identity lives in the *science* config, so it is read from `raw` -
    # the same place `require_antibiotic` reads the project block. It is NOT on
    # PipelineConfig, and it is NOT the machine overlay's: a per-machine
    # `reference` would let two machines disagree about which genome the
    # science was written against, which is precisely what D5 forbids.
    # `MachineConfig` has a `reference` field, but no overlay declares one, so
    # reading it here would have found nothing and verified nothing.
    identity = dict((getattr(config, "raw", None) or {}).get("reference") or {})
    if not identity:
        identity = dict(getattr(config, "reference", None) or {})

    if not identity:
        LOGGER.warning(
            "No reference is declared in the science configuration; there is "
            "nothing to verify, and any stage that needs a reference will say so."
        )
        return None

    accession = str(identity.get("accession") or identity.get("name") or "reference")
    length = identity.get("expected_length")
    return verify_reference_file(
        config.reference_fasta(),
        expected_length=int(length) if length is not None else None,
        expected_md5=identity.get("md5"),
        accession=accession,
    )


def verify_fai_index(fai_path: Path, *, expected_bases: int) -> int:
    """Check a `.fai` index against the reference's real base count.

    The FASTA is verified separately by :func:`verify_reference_file`, and
    re-checking it here would prove nothing: the source file is byte-identical
    whichever way it is indexed. **The index is what lies.**

    Observed on this machine: indexing the reference through `bgzip` +
    `samtools faidx` produced an index reporting 64,400 bases for a 6,264,404
    base file. Every tool then read a reference two orders of magnitude too
    small and produced confident nonsense — 64 alignments for a 6.4 Mb genome,
    and thousands of "variants" all sharing one constant QUAL. Indexing the
    uncompressed FASTA gives the correct count. The two are not interchangeable,
    so whichever path a caller takes, it has to come through here.

    Args:
        fai_path: The index to check.
        expected_bases: Base count from the D5-verified reference.

    Returns:
        The verified base count.

    Raises:
        DataContractError: The index is missing, empty, reports zero bases, or
            disagrees with ``expected_bases``. Every message names both numbers,
            because a refusal the operator cannot act on is just an obstacle.
    """
    fai_path = Path(fai_path)
    if not fai_path.is_file():
        raise DataContractError(
            "Reference index not found; nothing can be aligned or called "
            "against an unindexed reference",
            path=str(fai_path),
        )

    lines = [line for line in fai_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise DataContractError(
            "Reference index is empty", path=str(fai_path)
        )

    observed = 0
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        try:
            observed += int(parts[1])
        except ValueError as exc:
            raise DataContractError(
                "Reference index has an unreadable length field",
                path=str(fai_path),
                line=line,
            ) from exc

    if observed == 0:
        raise DataContractError(
            "Reference index reports zero bases. A `.fai` built by indexing a "
            "bgzipped FASTA has been seen to report a fraction of the real "
            "length, so index the uncompressed FASTA instead.",
            path=str(fai_path),
            expected_bases=expected_bases,
        )
    if observed != expected_bases:
        raise DataContractError(
            f"Reference index reports {observed} bases but the verified "
            f"reference has {expected_bases}. A `.fai` built by indexing a "
            "bgzipped FASTA has been seen to report a fraction of the real "
            "length. Every tool would then read a much smaller reference and "
            "produce confident wrong coordinates, so this stops here. Index the "
            "uncompressed FASTA, or delete the index and rebuild it.",
            path=str(fai_path),
            observed_bases=observed,
            expected_bases=expected_bases,
        )
    return observed
