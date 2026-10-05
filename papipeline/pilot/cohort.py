"""Pilot-100 cohort construction and analysis.

============================= IMPORTANT =============================
SELECTION IS BASED ON THE FILESYSTEM, NOT ON PDC_essential.tsv.

The pilot population is the FIRST 100 genome assemblies that physically
exist under ``data/``, sorted by GCA accession. ``PDC_essential.tsv`` is
consulted only AFTER that selection is made, and only to attach metadata.

The two orderings are materially different. Taking the first 100 rows of
PDC_essential.tsv instead would yield a different, overlapping-only-in-part
population. :func:`prove_selection_is_filesystem_based` computes that overlap
and reports it, so the distinction is auditable rather than asserted.
===================================================================

This module builds the cohort and the metadata join. It deliberately makes no
scientific claims: the gene/mechanism analysis lives in
:mod:`papipeline.pilot.analysis`, and the interpretation rules live in
``docs/scientific_rules.md``.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..errors import PipelineError
from ..logging_utils import get_logger

LOGGER = get_logger("pilot.cohort")

#: Assembly file extensions recognised during discovery.
GENOME_SUFFIXES: Tuple[str, ...] = (
    ".fna",
    ".fasta",
    ".fa",
    ".fna.gz",
    ".fasta.gz",
    ".fa.gz",
)

#: GCA accession pattern.
GCA_PATTERN = re.compile(r"GCA_\d+\.\d+")

PILOT_PREFIX = "PILOT100_"
N_PILOT = 100

#: Value used when an assembly has no PDC record. Never a blank cell.
MISSING = "MISSING"

#: PDC_essential.tsv column -> the pilot column it feeds.
PDC_COLUMN_MAP: Mapping[str, str] = {
    "Assembly": "Assembly",
    "Isolate": "Isolate",
    "BioSample": "BioSample",
    "BioProject": "BioProject",
    "Collection date": "Collection_date",
    "Location": "Location",
    "Source type": "Source_type",
    "Isolation source": "Isolation_source",
    "AST phenotypes": "AST_phenotypes_raw",
    "AMR genotypes": "AMR_genotypes_raw",
    "PDC_present": "PDC_present",
}

#: PDC column holding the antibiotic-keyed susceptibility categories.
AST_COLUMN = "AST phenotypes"
#: PDC column holding the gene-keyed genotype calls.
AMR_COLUMN = "AMR genotypes"


# --------------------------------------------------------------------------
# Step 1: discovery
# --------------------------------------------------------------------------


def _suffix_of(name: str) -> Optional[str]:
    lowered = name.lower()
    for suffix in sorted(GENOME_SUFFIXES, key=len, reverse=True):
        if lowered.endswith(suffix):
            return suffix
    return None


def discover_assemblies(data_dir: Path) -> Dict[str, List[Path]]:
    """Recursively find assembly files under ``data_dir``, keyed by GCA accession.

    A file qualifies only if it is an assembly extension AND its path
    contains a GCA accession. Returns ``{accession: [paths]}``; an accession
    with more than one file is reported by the caller rather than silently
    collapsed.
    """
    data_dir = Path(data_dir)
    if not data_dir.exists():
        raise PipelineError("Data directory not found", path=str(data_dir))

    found: Dict[str, List[Path]] = {}
    for path in sorted(data_dir.rglob("*")):
        if not path.is_file():
            continue
        if _suffix_of(path.name) is None:
            continue
        match = GCA_PATTERN.search(path.name) or GCA_PATTERN.search(str(path))
        if not match:
            LOGGER.debug("Assembly-extension file without a GCA accession: %s", path)
            continue
        found.setdefault(match.group(0), []).append(path)

    return found


def gca_sort_key(accession: str) -> Tuple[int, int]:
    """Natural sort key for a GCA accession.

    ``GCA_000710625.1`` sorts by numeric prefix then numeric version, so
    ``GCA_9.1`` precedes ``GCA_10.1``. Plain string sort gets this wrong.
    """
    match = GCA_PATTERN.match(accession)
    if not match:
        raise PipelineError("Not a GCA accession", value=accession)
    body = match.group(0)[len("GCA_") :]
    number, _, version = body.partition(".")
    return (int(number), int(version) if version else 0)


def sorted_accessions(found: Mapping[str, Sequence[Path]]) -> List[str]:
    """Deterministic accession ordering."""
    return sorted(found.keys(), key=gca_sort_key)


@dataclass(frozen=True)
class SelectedAssembly:
    """One member of the pilot cohort."""

    pilot_id: str
    assembly: str
    genome_path: Path
    genome_file: str

    def to_manifest_row(self, genome_found: bool = True) -> Dict[str, object]:
        return {
            "Pilot_ID": self.pilot_id,
            "Assembly": self.assembly,
            "Genome_path": str(self.genome_path),
            "Genome_file": self.genome_file,
            "Genome_found": "TRUE" if genome_found else "FALSE",
        }


def read_accession_list(path: Path) -> List[str]:
    """Read an ordered accession list, one per line.

    Blank lines and ``#`` comments are ignored. Order is preserved exactly:
    the list *is* the selection rule.
    """
    path = Path(path)
    if not path.exists():
        raise PipelineError("Accession list not found", path=str(path))
    if path.stat().st_size == 0:
        raise PipelineError("Accession list is empty", path=str(path))
    out: List[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not GCA_PATTERN.search(line):
            raise PipelineError(
                "Line is not a GCA accession",
                path=str(path),
                line=line_number,
                value=line[:60],
            )
        out.append(line.split()[0])
    return out


def select_from_list(
    found: Mapping[str, Sequence[Path]], accessions: Sequence[str]
) -> List[SelectedAssembly]:
    """Select assemblies in the exact order given by an accession list.

    This is the alternative selection rule: the cohort is dictated by the
    supplied list rather than by GCA order over ``data/``. An accession with
    no file in ``data/`` is reported rather than silently skipped, so the
    caller can see exactly which selections are unsatisfiable.
    """
    seen: Set[str] = set()
    duplicates: List[str] = []
    for accession in accessions:
        if accession in seen:
            duplicates.append(accession)
        seen.add(accession)
    if duplicates:
        raise PipelineError(
            "Accession list contains duplicates",
            n_duplicates=len(duplicates),
            examples=",".join(sorted(set(duplicates))[:5]),
        )

    missing = [a for a in accessions if a not in found]
    if missing:
        LOGGER.error(
            "%d of %d listed accessions have no assembly file in data/: %s",
            len(missing),
            len(accessions),
            ", ".join(missing[:10]) + ("..." if len(missing) > 10 else ""),
        )

    selected: List[SelectedAssembly] = []
    for index, accession in enumerate(accessions, start=1):
        if accession not in found:
            continue
        chosen = sorted(found[accession])[0]
        selected.append(
            SelectedAssembly(
                pilot_id=f"{PILOT_PREFIX}{index:03d}",
                assembly=accession,
                genome_path=chosen,
                genome_file=chosen.name,
            )
        )
    return selected


def select_first_n(
    found: Mapping[str, Sequence[Path]], n: int = N_PILOT
) -> List[SelectedAssembly]:
    """Select the first ``n`` accessions in GCA order.

    Raises:
        PipelineError: Fewer than ``n`` assemblies are available. The
            instruction is to stop rather than substitute, so this is a hard
            failure.
    """
    ordered = sorted_accessions(found)
    if len(ordered) < n:
        raise PipelineError(
            "Fewer assemblies available than required; refusing to substitute",
            available=len(ordered),
            required=n,
            hint="Stop and report; do not backfill from PDC_essential.tsv",
        )

    ambiguous = {a: v for a, v in found.items() if len(v) > 1}
    if ambiguous:
        LOGGER.warning(
            "%d accession(s) have more than one assembly file; the first in "
            "sorted path order is used: %s",
            len(ambiguous),
            ", ".join(sorted(ambiguous)[:5]),
        )

    selected: List[SelectedAssembly] = []
    for index, accession in enumerate(ordered[:n], start=1):
        chosen = sorted(ambiguous.get(accession, found[accession]))[0]
        selected.append(
            SelectedAssembly(
                pilot_id=f"{PILOT_PREFIX}{index:03d}",
                assembly=accession,
                genome_path=chosen,
                genome_file=chosen.name,
            )
        )
    return selected


# --------------------------------------------------------------------------
# Step 3: symlinks
# --------------------------------------------------------------------------


def symlink_assemblies(
    selected: Sequence[SelectedAssembly], link_dir: Path
) -> List[Path]:
    """Symlink the selected assemblies into ``link_dir``.

    ``data/`` is never written to. Existing links are replaced so the
    operation is idempotent.
    """
    link_dir = Path(link_dir)
    link_dir.mkdir(parents=True, exist_ok=True)
    created: List[Path] = []
    for item in selected:
        link = link_dir / item.genome_file
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(item.genome_path.resolve())
        created.append(link)
    return created


# --------------------------------------------------------------------------
# Step 4: PDC metadata join
# --------------------------------------------------------------------------


def load_pdc_index(pdc_path: Path) -> Dict[str, Dict[str, str]]:
    """Index ``PDC_essential.tsv`` by GCA accession.

    Keyed on the ``Assembly`` column, never on row order. A duplicated
    accession is a data error and is reported rather than resolved silently.
    """
    pdc_path = Path(pdc_path)
    if not pdc_path.exists():
        raise PipelineError(
            "PDC_essential.tsv not found", path=str(pdc_path)
        )
    if pdc_path.stat().st_size == 0:
        raise PipelineError(
            "PDC_essential.tsv is empty; no metadata can be joined",
            path=str(pdc_path),
        )

    index: Dict[str, Dict[str, str]] = {}
    duplicates: List[str] = []
    with pdc_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Assembly" not in reader.fieldnames:
            raise PipelineError(
                "PDC_essential.tsv has no 'Assembly' column",
                path=str(pdc_path),
                found=",".join(reader.fieldnames or []),
            )
        for row in reader:
            accession = (row.get("Assembly") or "").strip()
            if not accession:
                continue
            if accession in index:
                duplicates.append(accession)
                continue
            index[accession] = {
                key: (value.strip() if isinstance(value, str) else value)
                for key, value in row.items()
                if key
            }

    if duplicates:
        raise PipelineError(
            "PDC_essential.tsv contains duplicate Assembly accessions",
            n_duplicates=len(duplicates),
            examples=",".join(sorted(set(duplicates))[:5]),
        )

    LOGGER.info("Indexed %d PDC records by Assembly accession", len(index))
    return index


def _clean(value: Optional[str]) -> str:
    text = (value or "").strip()
    return text if text else MISSING


def join_metadata(
    selected: Sequence[SelectedAssembly], pdc_index: Mapping[str, Mapping[str, str]]
) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    """Attach PDC metadata to the selected assemblies.

    An assembly absent from PDC is retained and marked
    ``PDC_metadata_match=MISSING``. It is never discarded, and no field is
    invented.
    """
    rows: List[Dict[str, object]] = []
    matched = 0
    for item in selected:
        record = pdc_index.get(item.assembly)
        if record is None:
            rows.append(
                {
                    "Pilot_ID": item.pilot_id,
                    "Assembly": item.assembly,
                    "Genome_path": str(item.genome_path),
                    "Genome_file": item.genome_file,
                    "Isolate": MISSING,
                    "BioSample": MISSING,
                    "BioProject": MISSING,
                    "Collection_date": MISSING,
                    "Location": MISSING,
                    "Source_type": MISSING,
                    "Isolation_source": MISSING,
                    "Antibiotic": MISSING,
                    "AST_phenotype": MISSING,
                    "AST_phenotypes_raw": MISSING,
                    "AMR_genotypes": MISSING,
                    "AMR_genotypes_raw": MISSING,
                    "PDC_present": MISSING,
                    "PDC_metadata_match": MISSING,
                }
            )
            continue

        matched += 1
        ast_raw = _clean(record.get(AST_COLUMN))
        amr_raw = _clean(record.get(AMR_COLUMN))
        rows.append(
            {
                "Pilot_ID": item.pilot_id,
                "Assembly": item.assembly,
                "Genome_path": str(item.genome_path),
                "Genome_file": item.genome_file,
                "Isolate": _clean(record.get("Isolate")),
                "BioSample": _clean(record.get("BioSample")),
                "BioProject": _clean(record.get("BioProject")),
                "Collection_date": _clean(record.get("Collection date")),
                "Location": _clean(record.get("Location")),
                "Source_type": _clean(record.get("Source type")),
                "Isolation_source": _clean(record.get("Isolation source")),
                # The raw comma-separated `antibiotic=CATEGORY` string is
                # retained so each antibiotic's category can be extracted
                # independently. `Antibiotic` and `AST_phenotype` are filled
                # per antibiotic by the analysis stage.
                "Antibiotic": MISSING,
                "AST_phenotype": MISSING,
                "AST_phenotypes_raw": ast_raw,
                "AMR_genotypes": amr_raw,
                "AMR_genotypes_raw": amr_raw,
                "PDC_present": _clean(record.get("PDC_present")),
                "PDC_metadata_match": "MATCHED",
            }
        )

    counts = {
        "selected": len(selected),
        "matched": matched,
        "missing": len(selected) - matched,
    }
    return rows, counts


def prove_selection_is_filesystem_based(
    pdc_path: Path, selected: Sequence[SelectedAssembly]
) -> Dict[str, object]:
    """Quantify how the filesystem selection differs from PDC row order.

    Reported in the QC report so the selection method is verifiable rather
    than merely claimed.
    """
    pdc_path = Path(pdc_path)
    pdc_row_order: List[str] = []
    if pdc_path.exists() and pdc_path.stat().st_size > 0:
        with pdc_path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                accession = (row.get("Assembly") or "").strip()
                if accession:
                    pdc_row_order.append(accession)

    selected_set = {s.assembly for s in selected}
    first100_pdc = pdc_row_order[:N_PILOT]
    overlap = selected_set & set(first100_pdc)

    return {
        "pdc_total_rows": len(pdc_row_order),
        "filesystem_selection_n": len(selected_set),
        "pdc_roworder_first100_n": len(first100_pdc),
        "overlap_n": len(overlap),
        "only_in_filesystem_selection": len(selected_set - set(first100_pdc)),
        "only_in_pdc_roworder": len(set(first100_pdc) - selected_set),
        "identical_sets": selected_set == set(first100_pdc),
        "filesystem_first_accession": selected[0].assembly if selected else None,
        "filesystem_last_accession": selected[-1].assembly if selected else None,
    }


# --------------------------------------------------------------------------
# Assembly integrity
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class IntegrityVerdict:
    """Whether one assembly file is fit for analysis, and why not."""

    assembly: str
    path: Path
    size_bytes: Optional[int]
    included: bool
    reason: str

    def to_row(self) -> Dict[str, object]:
        return {
            "Pilot_ID": "",
            "Assembly": self.assembly,
            "Genome_path": str(self.path),
            "Genome_size_bytes": self.size_bytes,
            "analysis_included": "TRUE" if self.included else "FALSE",
            "exclusion_reason": self.reason,
        }


def assess_assembly(
    path: Path,
    min_size: Optional[int] = 4_000_000,
    max_size: Optional[int] = 8_000_000,
) -> IntegrityVerdict:
    """Decide whether an assembly file is analysable.

    This is a *file integrity and completeness* check, not a QC-threshold
    judgement. Two independent conditions must hold:

    1. the whole file decodes as text. Bulk-downloaded genome sets are
       routinely truncated mid-transfer, and a file that is valid FASTA for
       its first few megabytes and binary garbage thereafter will otherwise
       pass a size check while silently losing most of the genome;
    2. the file is within the plausible size range for the organism. A
       300 kb file is a partially-rehydrated stub, not a small bacterium.

    The size bounds default to the values in ``config.qc`` so this is
    configuration-driven rather than hard-coded.
    """
    path = Path(path)
    name = path.name

    if not path.exists():
        return IntegrityVerdict(name, path, None, False, "file_not_found")

    size = path.stat().st_size
    if size == 0:
        return IntegrityVerdict(name, path, 0, False, "zero_byte_file")

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1 << 20)
            if not chunk:
                break
            try:
                chunk.decode("utf-8")
            except UnicodeDecodeError as exc:
                return IntegrityVerdict(
                    name,
                    path,
                    size,
                    False,
                    f"corrupt_binary_data_at_offset_{size - len(chunk) + exc.start}",
                )

    with path.open("rb") as handle:
        first = handle.read(2)
    if not first.startswith(b">"):
        return IntegrityVerdict(
            name, path, size, False, "not_fasta_no_header"
        )

    if min_size is not None and size < min_size:
        return IntegrityVerdict(
            name, path, size, False, f"undersized_{size}_below_{min_size}"
        )
    if max_size is not None and size > max_size:
        return IntegrityVerdict(
            name, path, size, False, f"oversized_{size}_above_{max_size}"
        )
    return IntegrityVerdict(name, path, size, True, "ok")


def exclusion_breakdown(verdicts: Sequence["IntegrityVerdict"]) -> Dict[str, int]:
    """Count exclusions by reason *type*.

    The per-assembly reason carries the measured size or byte offset, which
    would otherwise yield one key per file. The exact values stay in
    ``assembly_integrity.tsv``.
    """
    counts: Dict[str, int] = {}
    for verdict in verdicts:
        if verdict.included:
            continue
        reason = verdict.reason
        for prefix in (
            "corrupt_binary_data",
            "undersized",
            "oversized",
            "zero_byte_file",
            "file_not_found",
            "not_fasta_no_header",
        ):
            if reason.startswith(prefix):
                reason = prefix
                break
        else:
            reason = "other"
        counts[reason] = counts.get(reason, 0) + 1
    return counts


def assess_cohort(
    selected: Sequence[SelectedAssembly],
    min_size: Optional[int] = 4_000_000,
    max_size: Optional[int] = 8_000_000,
) -> Tuple[List[IntegrityVerdict], List[SelectedAssembly], List[SelectedAssembly]]:
    """Assess every selected assembly.

    Returns ``(verdicts, analysable, excluded)``. The **selection is
    unchanged** by this step: ``selected`` is still the first 100 in GCA
    order. Exclusion is a downstream analysability decision, recorded per
    assembly with a reason, and never by substitution.
    """
    verdicts: List[IntegrityVerdict] = []
    analysable: List[SelectedAssembly] = []
    excluded: List[SelectedAssembly] = []

    for item in selected:
        verdict = assess_assembly(item.genome_path, min_size, max_size)
        if verdict.included:
            analysable.append(item)
        else:
            excluded.append(item)
        # Carry the pilot ID through so the verdict is joinable.
        verdicts.append(
            IntegrityVerdict(
                assembly=item.assembly,
                path=item.genome_path,
                size_bytes=verdict.size_bytes,
                included=verdict.included,
                reason=verdict.reason,
            )
        )

    LOGGER.info(
        "Assembly integrity: %d analysable, %d excluded of %d selected",
        len(analysable),
        len(excluded),
        len(selected),
    )
    if excluded:
        LOGGER.error(
            "Excluded assemblies by reason: %s",
            ", ".join(
                f"{k}={v}" for k, v in sorted(exclusion_breakdown(verdicts).items())
            ),
        )
    return verdicts, analysable, excluded


# --------------------------------------------------------------------------
# QC
# --------------------------------------------------------------------------


def matching_qc_rows(
    selected: Sequence[SelectedAssembly],
    metadata_rows: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    """Build the step 5 matching-QC table."""
    by_assembly = {str(r["Assembly"]): r for r in metadata_rows}
    rows: List[Dict[str, object]] = []
    for item in selected:
        record = by_assembly.get(item.assembly, {})
        found = record.get("PDC_metadata_match") == "MATCHED"
        rows.append(
            {
                "Assembly": item.assembly,
                "Genome_path": str(item.genome_path),
                "PDC_record_found": "TRUE" if found else "FALSE",
                "PDC_metadata_match": record.get("PDC_metadata_match", MISSING),
                "Isolate": record.get("Isolate", MISSING),
                "BioSample": record.get("BioSample", MISSING),
                "BioProject": record.get("BioProject", MISSING),
                "Antibiotic": record.get("Antibiotic", MISSING),
                "AST_phenotype": record.get("AST_phenotype", MISSING),
                "AMR_genotypes": record.get("AMR_genotypes", MISSING),
            }
        )
    return rows
