"""Sample manifest: the single source of truth for the set of samples.

Every stage joins on ``sample_id``. The manifest is built once per run from
``data/metadata/`` and is the reference against which the phylogenetic tree
tip set is validated (stage 10).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set

from .errors import DataContractError, SampleIdError
from .io.tsv import read_tsv
from .logging_utils import get_logger
from .models import Sample
from .pdc import load_pdc_index, row_value

LOGGER = get_logger("manifest")

#: Sample identifiers must be filesystem- and TSV-safe.
_FORBIDDEN = set(" \t\n/\\'\";:|,()[]{}#")

METADATA_FILENAME = "sample_metadata.tsv"
METADATA_REQUIRED = ("sample_id",)


def validate_sample_id(sample_id: Optional[str], pattern: Optional[str] = None) -> str:
    """Validate a single sample identifier.

    The structural checks below are always applied. ``pattern``, when given, is
    the *configured* format rule from ``science.yaml`` and is applied on top, so
    the rule can be tightened without editing code.

    Raises:
        SampleIdError: The identifier is empty, contains whitespace or
            forbidden characters, has an unusable length, or does not match the
            configured format.
    """
    if sample_id is None:
        raise SampleIdError("Sample ID is missing")
    value = sample_id.strip()
    if not value:
        raise SampleIdError("Sample ID is empty")
    if value != sample_id:
        raise SampleIdError(
            "Sample ID has leading or trailing whitespace",
            sample_id=sample_id,
            stripped=value,
        )
    if len(value) < 2 or len(value) > 128:
        raise SampleIdError(
            "Sample ID length is outside the permitted range 2-128",
            sample_id=value,
            length=len(value),
        )
    bad = sorted(_FORBIDDEN & set(value))
    if bad:
        raise SampleIdError(
            "Sample ID contains forbidden characters",
            sample_id=value,
            characters="".join(bad),
        )
    if not value.replace("-", "").replace("_", "").replace(".", "").isalnum():
        raise SampleIdError(
            "Sample ID must be alphanumeric with optional - _ . separators",
            sample_id=value,
        )
    if pattern is not None and not re.fullmatch(pattern, value):
        raise SampleIdError(
            "Sample ID does not match the configured format",
            sample_id=value,
            pattern=pattern,
        )
    return value




def _check_unique_attributes(
    samples: Sequence["Sample"], attributes: Sequence[str]
) -> None:
    """Every provenance attribute must identify at most one assembly.

    Isolate and BioSample are carried for traceability, not for matching, so the
    only thing that matters is that they are not ambiguous: two assemblies
    claiming one isolate, or one assembly claimed by two isolates, both mean the
    provenance chain is broken and the row cannot be trusted.
    """
    for attribute in attributes:
        seen: Dict[str, str] = {}
        for sample in samples:
            value = getattr(sample, attribute, None)
            if not value:
                continue
            previous = seen.get(value)
            if previous is not None:
                raise DataContractError(
                    f"Sample attribute {attribute!r} is claimed by more than one "
                    f"assembly, so provenance is ambiguous",
                    attribute=attribute,
                    value=value,
                    first_sample=previous,
                    second_sample=sample.sample_id,
                )
            seen[value] = sample.sample_id

@dataclass(frozen=True)
class SampleManifest:
    """An ordered, duplicate-free set of samples."""

    samples: Sequence[Sample]

    def __post_init__(self) -> None:
        seen: Set[str] = set()
        for sample in self.samples:
            validate_sample_id(sample.sample_id)
            if sample.sample_id in seen:
                raise DataContractError(
                    "Duplicate sample_id in the manifest", sample_id=sample.sample_id
                )
            seen.add(sample.sample_id)

    def __len__(self) -> int:
        return len(self.samples)

    def __iter__(self):
        return iter(self.samples)

    @property
    def sample_ids(self) -> List[str]:
        return [s.sample_id for s in self.samples]

    def get(self, sample_id: str) -> Optional[Sample]:
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        return None

    def require(self, sample_id: str) -> Sample:
        sample = self.get(sample_id)
        if sample is None:
            raise DataContractError(
                "Sample is not present in the manifest", sample_id=sample_id
            )
        return sample

    def index(self) -> Dict[str, Sample]:
        return {s.sample_id: s for s in self.samples}

    def with_source(self, source: str) -> "SampleManifest":
        return SampleManifest(
            [Sample(s.sample_id, s.assembly_path, source) for s in self.samples]
        )


def manifest_from_rows(
    rows: Iterable[Mapping[str, object]],
    source: str = "metadata",
    *,
    pattern: Optional[str] = None,
    unique_attributes: Sequence[str] = ("assembly", "isolate_id", "biosample"),
) -> SampleManifest:
    """Build a manifest from metadata rows, resolving assembly paths.

    Args:
        rows: Metadata rows keyed by column name.
        source: Recorded on every sample for provenance.
        pattern: The configured sample-id format rule, from ``science.yaml``.
        unique_attributes: Provenance attributes that must be unambiguous.
    """
    samples: List[Sample] = []
    for row in rows:
        raw_id = row.get("sample_id")
        sample_id = validate_sample_id(
            str(raw_id) if raw_id is not None else None, pattern=pattern
        )
        assembly_path = row.get("assembly_path") or row.get("assembly") or None
        samples.append(
            Sample(
                sample_id=sample_id,
                assembly_path=str(assembly_path) if assembly_path else None,
                source=source,
                assembly=_optional_str(row.get("assembly")),
                isolate_id=_optional_str(row.get("isolate_id")),
                biosample=_optional_str(row.get("biosample")),
                isolate_name=_optional_str(row.get("isolate_name")),
            )
        )
    _check_unique_attributes(samples, unique_attributes)
    manifest = SampleManifest(samples)
    LOGGER.info("Loaded %d samples from %s", len(manifest), source)
    return manifest


def _optional_str(value: object) -> Optional[str]:
    """A missing or sentinel value is absent, not the string "None"."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in {".", "-", "NA", "na", "none", "None"}:
        return None
    return text


def discover_manifest(
    metadata_dir: Path,
    require_files: bool = False,
    *,
    pattern: Optional[str] = None,
    unique_attributes: Sequence[str] = ("assembly", "isolate_id", "biosample"),
) -> SampleManifest:
    """Load the manifest from ``<metadata_dir>/sample_metadata.tsv``.

    This is the path a run actually takes, so the configured format rule and the
    provenance checks are applied here rather than only in
    :func:`manifest_from_rows`. A rule enforced on one path and not the other is
    not a rule.

    Args:
        metadata_dir: Directory holding ``sample_metadata.tsv``.
        require_files: When True, every sample must have an existing
            ``assembly_path``. Used by REAL mode; TEST mode tolerates absent
            sequence files because synthetic genomes are optional.
        pattern: The configured sample-id format rule.
        unique_attributes: Provenance attributes that must be unambiguous.
    """
    metadata_dir = Path(metadata_dir)
    path = metadata_dir / METADATA_FILENAME
    if not path.exists():
        raise DataContractError(
            "Sample metadata file not found",
            path=str(path),
            hint="Create data/<mode>/metadata/sample_metadata.tsv",
        )
    rows = read_tsv(
        path,
        required_columns=METADATA_REQUIRED,
        unique_columns=METADATA_REQUIRED,
    )

    samples: List[Sample] = []
    for row in rows:
        sample_id = validate_sample_id(str(row["sample_id"]), pattern=pattern)
        assembly = row.get("assembly_path")
        if require_files and not assembly:
            raise DataContractError(
                "Sample has no assembly_path but files are required",
                sample_id=sample_id,
            )
        if assembly:
            resolved = Path(str(assembly))
            if not resolved.is_absolute():
                # normpath, not resolve(): collapse the ".." without following
                # symlinks, so a stored path compares equal to the same path
                # written the other way round.
                resolved = Path(os.path.normpath(str(metadata_dir / resolved)))
            if require_files and not resolved.exists():
                raise DataContractError(
                    "Assembly file does not exist", sample_id=sample_id, path=str(resolved)
                )
            assembly = str(resolved)
        samples.append(
            Sample(
                sample_id=sample_id,
                assembly_path=assembly,
                source="manifest",
                assembly=_optional_str(row.get("assembly")),
                isolate_id=_optional_str(row.get("isolate_id")),
                biosample=_optional_str(row.get("biosample")),
                isolate_name=_optional_str(row.get("isolate_name")),
            )
        )

    _check_unique_attributes(samples, unique_attributes)
    manifest = SampleManifest(samples)
    LOGGER.info("Manifest: %d samples from %s", len(manifest), path)
    return manifest


def discover_pdc_manifest(
    pdc_path: Path,
    *,
    assembly_paths: Optional[Mapping[str, str]] = None,
) -> SampleManifest:
    """Build the cohort from ``PDC_essential.tsv``, keyed on ``Isolate``.

    Every isolate in the file is a cohort member. This is the point of the
    function: an isolate with no assembly accession is an ordinary member whose
    ``assembly_path`` is unset, not a broken row. A stage that needs its
    sequence refuses per sample through `papipeline.assemblies.locate_assembly`,
    and that refusal is the correct outcome. Building the cohort from assemblies
    instead would drop those isolates silently, and the run would then analyse a
    different set than it claims.

    Nothing here assumes a cohort size, compares against an expected count, or
    slices. Two isolates and 967 are built by the same code path.

    Args:
        pdc_path: The PDC table.
        assembly_paths: Explicit sequence path per isolate, used as-is with no
            pattern matching. This is how a bounded run points a small cohort at
            its own genome directory.

    Returns:
        A manifest in file order, one sample per isolate.
    """
    rows = load_pdc_index(pdc_path, key="Isolate")
    supplied = dict(assembly_paths or {})
    samples = [
        Sample(
            sample_id=isolate,
            assembly_path=supplied.get(isolate),
            source="PDC_essential.tsv",
            assembly=row_value(row, "Assembly"),
            isolate_id=isolate,
            biosample=row_value(row, "BioSample"),
            isolate_name=row_value(row, "Isolate"),
        )
        for isolate, row in rows.items()
    ]
    if not samples:
        raise DataContractError(
            "PDC_essential.tsv yielded no cohort members",
            path=str(pdc_path),
        )
    return SampleManifest(samples=samples)
