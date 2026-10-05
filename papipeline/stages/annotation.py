"""Stage 2 - genome annotation.

Bakta is the preferred backend. This module provides:

* :func:`parse_bakta_gff` - a tolerant parser for Bakta's ``bakta.gff``;
* :func:`parse_bakta_tsv` - a parser for ``bakta.tsv`` when preferred;
* :func:`standardise` - conversion of either into
  :class:`~papipeline.models.AnnotationRecord`;
* :func:`run` - the orchestration entry point;
* :func:`decide_reuse` - whether an existing Bakta output may be reused
  instead of re-invoking the tool (see :data:`REUSE_TOOL_OUTPUT_MODES`).

The parsers are pure functions over text, so they are fully testable against
checked-in fixtures without Bakta installed and without any genome data.

Deliberate behaviour: an annotation row whose ``gene_name`` is absent is
retained with ``gene_name=None``. Genes without a name are still real loci
and dropping them would silently bias the pan-genome.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config.loader import PipelineConfig
from ..errors import DataContractError, PipelineError, StageError
from ..io.tsv import read_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import AnnotationRecord, RunMode

LOGGER = get_logger("stages.annotation")

ANNOTATION_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "contig_id",
    "gene_id",
    "gene_name",
    "product",
    "gene_type",
    "start",
    "end",
    "strand",
    "annotation_source",
)

#: Bakta CDS qualifier keys that carry a locus tag.
_LOCUS_TAG_KEYS = ("locus_tag", "ID")

_ATTR_RE = re.compile(r'(?P<key>[A-Za-z_][A-Za-z0-9_]*)=(?P<value>"[^"]*"|[^;]*)')


def gff_features_are_covered(
    gff: Sequence[Dict[str, Optional[str]]],
    features: Sequence[Dict[str, Optional[str]]],
) -> bool:
    """True when every GFF feature has a counterpart in the TSV table.

    **Not an equality check.** It used to be `len(gff) == len(features)`, which
    rejected every real genome: Bakta's `inference.tsv` carries rows the GFF
    filter legitimately excludes - derived features, and rows whose gene id is
    `-`. On the smoke cohort that is 5,912 GFF features inside 5,944 TSV rows,
    and the length check discarded all ten annotated assemblies, leaving
    `n_records = 0` for every cohort member.

    Containment is the invariant that actually catches a broken pair. If the two
    tables came from different genomes, or one were truncated, some GFF gene id
    would have no TSV counterpart. Equality proves nothing extra and rejects
    valid data besides.

    Rows with no gene id (`-`) are excluded from both sides: they cannot be
    matched, and treating them as present on both would make containment
    vacuously true for a pair with nothing in common.
    """
    gff_ids = {
        row.get("gene_id")
        for row in gff
        if row.get("gene_id") and row.get("gene_id") != "-"
    }
    if not gff_ids:
        return False
    feature_ids = {
        row.get("gene_id")
        for row in features
        if row.get("gene_id") and row.get("gene_id") != "-"
    }
    return gff_ids <= feature_ids


def parse_bakta_gff(text: str) -> List[Dict[str, Optional[str]]]:
    """Parse a Bakta GFF3 body into raw feature dicts.

    Only ``CDS`` and ``tRNA``/``rRNA`` features with a name are kept; the
    GFF is parsed leniently because Bakta has changed column conventions
    across versions and the standardised schema is what downstream stages
    actually depend on.

    Returns:
        Raw feature dicts with keys ``seqid``, ``gene_id``, ``gene_name``,
        ``product``, ``gene_type``, ``start``, ``end``, ``strand``.
    """
    features: List[Dict[str, Optional[str]]] = []

    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            LOGGER.debug("Skipping short GFF line %d", line_number)
            continue
        seqid, _source, feature_type, start, end, _score, strand, _phase, attrs = fields[:9]
        if feature_type not in ("CDS", "tRNA", "rRNA", "ncRNA"):
            continue

        attributes = _parse_gff_attributes(attrs)
        gene_name = (
            attributes.get("gene")
            or attributes.get("locus_tag")
            or attributes.get("ID")
            or attributes.get("Name")
        )
        if not gene_name:
            LOGGER.debug("Skipping unnamed feature at GFF line %d", line_number)
            continue

        features.append(
            {
                "seqid": seqid,
                "gene_id": attributes.get("locus_tag") or attributes.get("ID") or gene_name,
                "gene_name": gene_name,
                "product": attributes.get("product") or attributes.get("Name"),
                "gene_type": feature_type,
                "start": start,
                "end": end,
                "strand": strand,
            }
        )

    return features


def _parse_gff_attributes(attributes: str) -> Dict[str, str]:
    """Parse GFF3 attribute column into a dict, stripping quotes and URL escapes."""
    out: Dict[str, str] = {}
    for match in _ATTR_RE.finditer(attributes):
        key = match.group("key")
        value = match.group("value").strip().strip('"')
        out[key] = value.replace("%2C", ",").replace("%3B", ";").replace("%3D", "=")
    return out


def parse_bakta_tsv(text: str) -> List[Dict[str, Optional[str]]]:
    """Parse a Bakta ``bakta.tsv`` body.

    Expected header: ``Sequence Id, Type, Strand, Locus Tag, Start, Stop,
    Length, GC, Gene, Product, Contig``.
    """
    from io import StringIO

    rows = read_tsv_text(
        StringIO(text),
        required=("Sequence Id", "Locus Tag", "Start", "Stop"),
    )
    out: List[Dict[str, Optional[str]]] = []
    for row in rows:
        gene_name = row.get("Gene") or row.get("Locus Tag")
        if not gene_name:
            continue
        out.append(
            {
                "seqid": row.get("Contig") or row.get("Sequence Id"),
                "gene_id": row.get("Locus Tag") or gene_name,
                "gene_name": gene_name,
                "product": row.get("Product"),
                "gene_type": row.get("Type") or "CDS",
                "start": row.get("Start"),
                "end": row.get("Stop"),
                "strand": row.get("Strand"),
            }
        )
    return out


#: Bakta writes its annotation table tab-delimited with the header on a
#: ``#``-prefixed line, e.g. ``#Sequence Id<TAB>Type<TAB>Start<TAB>...``.
#: The human-readable preamble above it starts ``# `` (hash-space).
BAKTA_TSV_COLUMNS: Tuple[str, ...] = (
    "Sequence Id",
    "Type",
    "Start",
    "Stop",
    "Strand",
    "Locus Tag",
    "Gene",
    "Product",
    "DbXrefs",
)


def read_tsv_text(handle, required: Sequence[str] = ()) -> List[Dict[str, Optional[str]]]:
    """Reader for Bakta's annotation table.

    Two real formats are supported, because both occur in the wild:

    * **tab-delimited with a ``#``-prefixed header** (Bakta 1.12, verified
      against real output)::

          # Annotated with Bakta
          # Software: v1.12.1
          #Sequence Id<TAB>Type<TAB>Start<TAB>Stop<TAB>...

    * **comma-delimited with a plain header** (legacy ``bakta.tsv`)::

          Sequence Id,Type,Strand,Locus Tag,Start,Stop

    The delimiter is inferred from the header line rather than assumed, and
    ``# `` human-readable preamble is never mistaken for the header.
    """
    import csv

    rows: List[Dict[str, Optional[str]]] = []
    fieldnames: Optional[List[str]] = None
    delimiter = "\t"

    for raw in handle:
        line = raw.rstrip("\n").rstrip("\r")
        if not line:
            continue

        if line.startswith("#"):
            candidate_text = line.lstrip("#")
            # The header is the comment line that names a real column.
            if "Sequence Id" in candidate_text and (
                "\t" in candidate_text or "," in candidate_text
            ):
                delimiter = "\t" if "\t" in candidate_text else ","
                fieldnames = [
                    c.strip() for c in candidate_text.split(delimiter)
                ]
            continue

        if fieldnames is None:
            # No header comment: this line is the header.
            delimiter = "\t" if "\t" in line else ","
            fieldnames = [c.strip() for c in line.split(delimiter)]
            continue

        fields = next(csv.reader([line], delimiter=delimiter))
        row: Dict[str, Optional[str]] = {}
        for index, name in enumerate(fieldnames):
            value = fields[index] if index < len(fields) else None
            text = (value or "").strip()
            row[name] = text or None
        rows.append(row)

    if fieldnames is None:
        raise DataContractError("Bakta TSV has no header line")
    missing = [c for c in required if c not in fieldnames]
    if missing:
        raise DataContractError(
            "Bakta TSV is missing required columns", missing=",".join(missing)
        )
    return rows


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        LOGGER.debug("Non-integer coordinate: %r", value)
        return None


def standardise(
    raw_features: Iterable[Dict[str, Optional[str]]],
    sample_id: str,
    source: str = "bakta",
) -> List[AnnotationRecord]:
    """Convert raw parsed features into the internal annotation schema."""
    records: List[AnnotationRecord] = []
    for feature in raw_features:
        records.append(
            AnnotationRecord(
                sample_id=sample_id,
                contig_id=str(feature.get("seqid") or ""),
                gene_id=str(feature.get("gene_id") or ""),
                gene_name=feature.get("gene_name"),
                product=feature.get("product"),
                gene_type=str(feature.get("gene_type") or "unknown"),
                start=_to_int(feature.get("start")),
                end=_to_int(feature.get("end")),
                strand=feature.get("strand"),
                annotation_source=source,
            )
        )
    return records


def read_standardised(path: Path) -> List[AnnotationRecord]:
    """Read a standardised annotation TSV back into records."""
    rows = read_tsv(
        path,
        required_columns=("sample_id", "gene_name"),
    )
    return [
        AnnotationRecord(
            sample_id=str(row["sample_id"]),
            contig_id=str(row.get("contig_id") or ""),
            gene_id=str(row.get("gene_id") or ""),
            gene_name=row.get("gene_name"),
            product=row.get("product"),
            gene_type=str(row.get("gene_type") or "unknown"),
            start=_to_int(row.get("start")),
            end=_to_int(row.get("end")),
            strand=row.get("strand"),
            annotation_source=str(row.get("annotation_source") or "unknown"),
        )
        for row in rows
    ]


def annotation_dir_for(intermediate_root: Path, sample_id: str) -> Path:
    """Conventional per-sample annotation directory."""
    return Path(intermediate_root) / "annotation" / sample_id


def load_from_intermediate(
    intermediate_root: Path, manifest: SampleManifest
) -> Dict[str, List[AnnotationRecord]]:
    """Load standardised annotations for every sample in the manifest.

    A sample with no annotation file maps to an empty list, so downstream
    stages see "not annotated" rather than a missing key.
    """
    root = Path(intermediate_root)
    out: Dict[str, List[AnnotationRecord]] = {}
    for sample in manifest:
        path = root / "annotation" / f"{sample.sample_id}.annotation.tsv"
        out[sample.sample_id] = read_standardised(path) if path.exists() else []
    missing = [sid for sid, recs in out.items() if not recs]
    if missing:
        LOGGER.warning(
            "%d samples have no annotation records: %s",
            len(missing),
            ",".join(missing[:10]) + ("..." if len(missing) > 10 else ""),
        )
    return out


def gene_names(records: Sequence[AnnotationRecord]) -> List[str]:
    """Distinct, non-null gene names in a record set."""
    return sorted({r.gene_name for r in records if r.gene_name})


# --------------------------------------------------------------------------
# Reusing Bakta output instead of re-invoking Bakta
# --------------------------------------------------------------------------
#
# What this must never become: "trust whatever is already on disk". The reuse
# path runs the SAME conversion chain a real run does - `parse_bakta_tsv`,
# `parse_bakta_gff`, `gff_features_are_covered`, `standardise`,
# `_write_standardised` - and differs in exactly one respect: the tool is not
# executed. The consequence is that provenance can no longer be established by
# watching a process, so it is established from the output itself.
#
# Three facts, each from the output rather than from a directory name:
#
# 1. the file set is present and non-empty - the FEATURE table, not the
#    inference table, because `bakta_outputs` returning `*.inference.tsv` once
#    lost every gene name and produced a uniform false negative;
# 2. the GFF's gene ids are contained in the feature table's (the existing
#    rule, unchanged - reuse skips the subprocess, not the validation);
# 3. the `# Software:` / `# Database:` header names the database the run is
#    configured against.

#: The three settings `annotation.reuse_tool_output` accepts. Stated here as
#: well as in the loader because this module is what acts on them, and a reader
#: of the stage should not have to open the loader to know what the modes do.
REUSE_TOOL_OUTPUT_MODES: Tuple[str, ...] = ("off", "prefer", "require")

#: Column order of the per-sample provenance record this stage writes beside
#: the Bakta output. Written, not merely returned: a run that reused a genome's
#: annotation has to be able to say afterwards which bytes it trusted, and a
#: value that only existed inside the process that made the decision is gone by
#: the time anyone asks.
REUSE_PROVENANCE_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "reused",
    "reason",
    "detail",
    "source_dir",
    "genome_stem",
    "sha256_feature_tsv",
    "sha256_gff3",
    "bakta_version",
    "database_string",
)

#: Where that record lands, relative to the Bakta output root. Nothing reads
#: it; it exists so the decision is auditable after the run.
REUSE_PROVENANCE_FILENAME = "reuse_provenance.tsv"

#: Named refusal reasons. A reason token rather than prose, because the stage's
#: existing `failed` list is read by counts and a reader needs to tell "no
#: output there" from "the output there is from another database".
REUSE_MISSING = "no_reusable_output"
REUSE_UNREADABLE = "header_unreadable"
REUSE_HEADER_INCOMPLETE = "header_incomplete"
REUSE_DATABASE_MISMATCH = "database_mismatch"
REUSE_GFF_TSV_DISAGREE = "gff_tsv_disagree"

#: The refusal tokens above, as a set, so the stage's summary counts can tell a
#: reuse refusal from a Bakta state without string-matching the message.
REUSE_REFUSAL_REASONS = frozenset(
    {
        REUSE_MISSING,
        REUSE_UNREADABLE,
        REUSE_HEADER_INCOMPLETE,
        REUSE_DATABASE_MISMATCH,
        REUSE_GFF_TSV_DISAGREE,
    }
)

#: How many lines of an output are searched for the `# Software:` /
#: `# Database:` block. Bakta writes them as the first lines of the file, so a
#: bounded read keeps this off the multi-megabyte tables - and a bounded read
#: that found no header is a refusal, not a licence to scan the whole file
#: hoping.
_HEADER_PROBE_LINES = 32

_SOFTWARE_LINE = re.compile(r"^#\s*Software:\s*v?(?P<version>\S+)\s*$")
_DATABASE_LINE = re.compile(r"^#\s*Database:\s*(?P<database>.+?)\s*$")


class ReusedOutputRefused(PipelineError):
    """An existing Bakta output cannot be trusted. Named per isolate.

    Carrying this as a typed error is a statement about the CHECK, not about
    the caller: :func:`decide_reuse` catches it and turns it into a refusal on
    the returned :class:`ReusedOutput`. It exists as an exception because the
    checks are chained - the first thing wrong with a candidate is the useful
    answer, and an exception is what stops a caller from quietly reading past
    it. The stage never lets one escape to the cohort: this stage records an
    unusable genome and carries on, and reuse has to behave the same way or one
    unverifiable isolate would end a bounded run.
    """


@dataclass(frozen=True)
class BaktaHeader:
    """What an output file says about the run that produced it."""

    software: str
    database: str


@dataclass(frozen=True)
class ReusedOutput:
    """One isolate's reuse decision, and the provenance behind it.

    Every field is present whether or not the output was reused. A refusal with
    a blank provenance row is what distinguishes "we did not use it" from "we
    did not look", and a reuse with a blank source path or digest is not a reuse
    anybody could check later.
    """

    sample_id: str
    reused: bool
    reason: str
    detail: str
    source_dir: Path
    genome_stem: str
    feature_tsv: Optional[Path] = None
    gff: Optional[Path] = None
    sha256_feature_tsv: str = ""
    sha256_gff3: str = ""
    bakta_version: str = ""
    database_string: str = ""

    def as_row(self) -> Dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "reused": "true" if self.reused else "false",
            "reason": self.reason,
            "detail": self.detail,
            "source_dir": str(self.source_dir),
            "genome_stem": self.genome_stem,
            "sha256_feature_tsv": self.sha256_feature_tsv,
            "sha256_gff3": self.sha256_gff3,
            "bakta_version": self.bakta_version,
            "database_string": self.database_string,
        }


def expected_bakta_database_string(database_dir: Path) -> str:
    """The ``# Database:`` string the configured database makes Bakta write.

    Derived from that database's own ``version.json`` and formatted the way
    **Bakta** formats it, which is::

        v{major}.{minor}, {type}

    Verified against Bakta's installed source rather than assumed -
    ``bakta/main.py:627`` and ``bakta/json_io.py:195`` both write exactly that
    f-string into both the TSV and the GFF header.

    **Not** :func:`papipeline.adapters.bakta_db.database_version`. That
    normalises to ``"{major}.{minor} ({type}, {date})"`` because it is compared
    against ``config/references.tsv``, whose schema has no column for the
    edition and so carries the release date instead. Reusing that string here
    would compare ``v6.0, light`` against ``6.0 (light, 2025-02-24)``, which
    can never be equal, so every genuinely reusable output would be refused for
    a reason that is about string shapes rather than about the database. The
    two forms answer two different questions and are built separately.

    Raises:
        StageError: the directory records no usable version. That is an
            incomplete database - Bakta itself refuses one - not an unknown
            version, so it is refused rather than guessed at.
    """
    from ..adapters import bakta_db

    path = Path(database_dir) / bakta_db.VERSION_FILE
    if not path.is_file():
        raise StageError(
            f"{path} does not exist, so the database string Bakta writes into "
            f"its output cannot be derived from {database_dir}. Existing output "
            "can therefore not be checked against the database this run is "
            "configured for. Provision a complete database out of band; this "
            "pipeline never downloads one during a run.",
            stage="annotation",
            database_dir=str(database_dir),
        )
    try:
        import json

        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise StageError(
            f"{path} could not be read as JSON, so the expected database string "
            f"is unknown: {exc}",
            stage="annotation",
            database_dir=str(database_dir),
        ) from exc

    major = payload.get("major")
    minor = payload.get("minor")
    edition = payload.get("type")
    if major is None or minor is None or not edition:
        raise StageError(
            f"{path} records no major/minor/type, so the string Bakta writes "
            f"into its output cannot be derived from it (got major={major!r}, "
            f"minor={minor!r}, type={edition!r}).",
            stage="annotation",
            database_dir=str(database_dir),
        )
    return f"v{major}.{minor}, {edition}"


def read_bakta_header(path: Path) -> BaktaHeader:
    """Read the ``# Software:`` / ``# Database:`` block from a Bakta output.

    Raises:
        ReusedOutputRefused: the file cannot be read, or records no software
            version, or records no database.

            A missing header is a refusal and never a default. Assuming the
            configured database because nothing contradicts it is precisely how
            output annotated against a *different* database gets accepted as
            current, and the failure it causes is invisible: the gene calls
            differ and every downstream stage reports them as biology.
    """
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            head = "".join(
                line for _, line in zip(range(_HEADER_PROBE_LINES), handle)
            )
    except OSError as exc:
        raise ReusedOutputRefused(
            f"existing Bakta output is unreadable, so its database cannot be "
            f"verified: {path} ({exc})",
            path=str(path), reason=REUSE_UNREADABLE,
        ) from exc

    software: Optional[str] = None
    database: Optional[str] = None
    for line in head.splitlines():
        match = _SOFTWARE_LINE.match(line)
        if match is not None:
            software = match.group("version")
            continue
        match = _DATABASE_LINE.match(line)
        if match is not None:
            database = match.group("database")

    if software is None or database is None:
        absent = [
            name
            for name, value in (("Software", software), ("Database", database))
            if value is None
        ]
        raise ReusedOutputRefused(
            f"existing Bakta output records no {' and no '.join(absent)} line in "
            f"its first {_HEADER_PROBE_LINES} lines ({path}), so the tool and "
            "database it was produced with cannot be established from it. Refusing "
            "rather than assuming the configured ones: an assumed pin is exactly "
            "how stale output from another database is accepted as current, and "
            "the resulting gene calls differ silently.",
            path=str(path), reason=REUSE_HEADER_INCOMPLETE,
        )
    return BaktaHeader(software=software, database=database)


def _sha256_file(path: Path) -> str:
    """sha256 of a file, read in chunks so a 20 Mb table is not slurped."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _refusal(
    sample_id: str, reason: str, detail: str, source_dir: Path, genome_stem: str
) -> ReusedOutput:
    return ReusedOutput(
        sample_id=sample_id,
        reused=False,
        reason=reason,
        detail=detail,
        source_dir=source_dir,
        genome_stem=genome_stem,
    )


def decide_reuse(
    sample: Any,
    *,
    out_root: Path,
    database_dir: Path,
) -> ReusedOutput:
    """May this one genome's existing Bakta output be reused as-is?

    The candidate set is looked for exactly where :func:`run_bakta` puts it:
    ``out_root / sample.sample_id``, with the files named by
    :func:`papipeline.adapters.bakta.genome_stem_for`, which is the stem of the
    assembly filename. Directory and file stem coincide only when the assembly
    happens to be called after the sample - true for the smoke cohort and an
    accident, not a rule - so the stem is asked of the adapter rather than
    assumed equal to the id.

    Nothing is copied, moved, renamed or deleted, and Bakta is not invoked.

    Checks, and what each one is for:

    1. the feature TSV, the GFF3 and the inference TSV are present and
       non-empty. The feature table, not the inference table: they share a
       header preamble, and reading the inference table instead once lost every
       gene name in the cohort. The protein FASTA is **not** required - this
       stage never reads it (``bakta_spec``'s ``require_faa`` defaults off), so
       requiring it would reject a set this stage could actually use.
    2. the header names the configured database, byte for byte. A mismatch
       refusal quotes both strings, because "not reused" without them is a
       refusal a reader cannot act on.
    3. the GFF's gene ids are contained in the feature table's, by the existing
       rule in :func:`gff_features_are_covered`.

    The ``Software:`` version is recorded but not compared. The database's
    ``version.json`` carries ``software-min``, which is a *minimum* and a
    different shape from the ``v1.12.1`` Bakta writes; turning one into the
    other to compare them would be a translation between shapes, and it would
    be wrong - ``software-min`` is 1.11 for a database whose real output says
    v1.12.1.

    Returns:
        A :class:`ReusedOutput` that never raises: a refusal is a value here,
        for the same reason ``run`` records an unusable genome rather than
        raising.
    """
    from ..adapters.bakta import BaktaOutputs, genome_stem_for

    stem = genome_stem_for(sample)
    directory = Path(out_root) / sample.sample_id
    outputs = BaktaOutputs(
        out_dir=directory, sample_id=sample.sample_id, genome_stem=stem
    )

    def _refuse(reason: str, detail: str) -> ReusedOutput:
        return _refusal(sample.sample_id, reason, detail, directory, stem)

    required = (
        (outputs.annotation_tsv, "feature table"),
        (outputs.gff, "GFF3"),
        (outputs.inference_tsv, "inference table"),
    )
    absent: List[str] = []
    empty: List[str] = []
    for path, role in required:
        if not path.is_file():
            absent.append(f"{role} {path.name}")
        elif path.stat().st_size == 0:
            empty.append(f"{role} {path.name}")
    if absent or empty:
        parts = []
        if absent:
            parts.append("missing: " + ", ".join(absent))
        if empty:
            parts.append("empty: " + ", ".join(empty))
        return _refuse(
            REUSE_MISSING,
            f"no reusable Bakta output for {sample.sample_id} under {directory} "
            f"({'; '.join(parts)})",
        )

    try:
        header = read_bakta_header(outputs.annotation_tsv)
    except ReusedOutputRefused as exc:
        return _refuse(
            str(exc.context.get("reason", REUSE_UNREADABLE)), str(exc)
        )

    expected = expected_bakta_database_string(database_dir)
    if header.database != expected:
        return _refuse(
            REUSE_DATABASE_MISMATCH,
            f"existing Bakta output for {sample.sample_id} records database "
            f"{header.database!r} but the configured database at {database_dir} "
            f"produces {expected!r}. Not reused: gene calls from a different "
            "database are a different annotation, and reusing them would feed "
            "that difference into mlst, amr and virulence without saying so.",
        )

    _, disagreement = _records_from_bakta_tables(
        outputs.annotation_tsv, outputs.gff, sample.sample_id
    )
    if disagreement:
        return _refuse(
            REUSE_GFF_TSV_DISAGREE,
            f"existing Bakta output for {sample.sample_id} is not usable: "
            f"{disagreement} ({outputs.gff.name} vs "
            f"{outputs.annotation_tsv.name})",
        )

    return ReusedOutput(
        sample_id=sample.sample_id,
        reused=True,
        reason="verified",
        detail=(
            f"reused existing Bakta output for {sample.sample_id} from "
            f"{directory} without invoking Bakta"
        ),
        source_dir=directory,
        genome_stem=stem,
        feature_tsv=outputs.annotation_tsv,
        gff=outputs.gff,
        sha256_feature_tsv=_sha256_file(outputs.annotation_tsv),
        sha256_gff3=_sha256_file(outputs.gff),
        bakta_version=header.software,
        database_string=header.database,
    )


def write_reuse_provenance(
    out_root: Path, provenance: Mapping[str, ReusedOutput]
) -> Optional[Path]:
    """Write the per-sample reuse provenance beside the Bakta output.

    One row per sample the stage decided about, reused or not. Returns the
    path, or ``None`` when there was nothing to record - the function is a no-op
    under ``reuse_tool_output: off``, so the pre-existing run gains no file.
    """
    from ..io.tsv import write_tsv

    if not provenance:
        return None
    directory = Path(out_root)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / REUSE_PROVENANCE_FILENAME
    write_tsv(
        path,
        (provenance[sample_id].as_row() for sample_id in sorted(provenance)),
        list(REUSE_PROVENANCE_COLUMNS),
    )
    return path


def _probe_reused_output(
    sample: Any,
    *,
    out_root: Path,
    database_dir: Path,
) -> ReusedOutput:
    """Provenance for a genome this run annotated itself.

    Same fields, ``reused=False``: a reader asking "which bytes were these gene
    calls made from?" gets an answer for every genome, not only for the ones
    that happened to skip the tool.
    """
    from ..adapters.bakta import BaktaOutputs, genome_stem_for

    stem = genome_stem_for(sample)
    directory = Path(out_root) / sample.sample_id
    outputs = BaktaOutputs(
        out_dir=directory, sample_id=sample.sample_id, genome_stem=stem
    )
    detail = f"{sample.sample_id} was annotated by Bakta in this run"
    try:
        header = read_bakta_header(outputs.annotation_tsv)
    except ReusedOutputRefused:
        return ReusedOutput(
            sample_id=sample.sample_id, reused=False, reason="not_reused",
            detail=detail, source_dir=directory, genome_stem=stem,
        )
    return ReusedOutput(
        sample_id=sample.sample_id, reused=False, reason="not_reused",
        detail=detail, source_dir=directory, genome_stem=stem,
        feature_tsv=outputs.annotation_tsv, gff=outputs.gff,
        sha256_feature_tsv=_sha256_file(outputs.annotation_tsv),
        sha256_gff3=_sha256_file(outputs.gff),
        bakta_version=header.software, database_string=header.database,
    )


def _records_from_bakta_tables(
    feature_tsv: Path, gff_path: Path, sample_id: str
) -> Tuple[Optional[List[AnnotationRecord]], str]:
    """Convert one genome's Bakta tables, or say why they are not usable.

    One implementation, called by both the tool path and the reuse path, so
    "the reused annotation was converted the same way" is a property of the
    code rather than a promise about two branches agreeing.

    Returns:
        ``(records, "")`` or ``(None, reason)``. The reason names the count of
        unmatched GFF ids rather than comparing counts, because Bakta's TSV
        legitimately holds rows the GFF filter excludes.
    """
    features = parse_bakta_tsv(Path(feature_tsv).read_text(encoding="utf-8"))
    gff = parse_bakta_gff(Path(gff_path).read_text(encoding="utf-8"))
    if not gff_features_are_covered(gff, features):
        missing = len({
            row.get("gene_id") for row in gff
            if row.get("gene_id") and row.get("gene_id") != "-"
        } - {
            row.get("gene_id") for row in features
            if row.get("gene_id") and row.get("gene_id") != "-"
        })
        return None, (
            f"{missing} GFF feature(s) absent from the "
            "TSV table (tables disagree or one is truncated)"
        )
    return standardise(features, sample_id, source="bakta"), ""


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    *,
    store: Any = None,
    event_sink: Any = None,
    run_key: str = "pipeline",
    database: Optional[Path] = None,
    tool_version: Optional[str] = None,
    database_version: Optional[str] = None,
) -> Dict[str, List[AnnotationRecord]]:
    """Stage 2 entry point.

    In TEST mode the standardised fixtures written by the synthetic
    generator are read, exactly as before.

    In REAL mode each genome is annotated by Bakta through
    :func:`papipeline.adapters.bakta.run_bakta`, which invokes the real
    tool and holds its outputs to a contract. A genome is a success only
    when the process exited cleanly *and* its annotation satisfied that
    contract, so a zero exit that produced an empty or wrong table is
    ``INVALID`` or ``INCOMPLETE`` rather than a pass.

    Annotation is a per-genome task even though this function is called
    once: that is the granularity at which Bakta runs, and it is what makes
    a single genome's failure attributable.

    ``annotation.reuse_tool_output`` may substitute an existing Bakta output
    for a tool invocation, per genome, through :func:`decide_reuse`. Under
    ``require`` a genome whose output does not verify is recorded as failed
    with the refusal by name and Bakta is never invoked - which is what makes
    a run against pre-computed annotation possible at all. Every mode stays
    per-genome and record-and-continue: one unverifiable isolate must not end
    a cohort that is mostly unprepared by construction.
    """
    if mode is not RunMode.REAL:
        return load_from_intermediate(intermediate_root, manifest)

    from ..adapters.bakta import run_bakta
    from ..observatory.wiring import config_hash

    reuse_mode = config.reuse_tool_output()
    genomes_dir = Path(config.assembly_root(mode))
    # Bakta writes where the pipeline expects per-genome tool output, and the
    # standardised table lands where load_from_intermediate looks for it.
    out_root = Path(intermediate_root) / "bakta"
    standardised_root = Path(intermediate_root)
    failed: List[str] = []
    standardised: Dict[str, List[AnnotationRecord]] = {}
    provenance: Dict[str, ReusedOutput] = {}

    # Resolved once for the cohort, not per genome. It re-reads
    # config/references.tsv and compares the database's own version against the
    # pin, so calling it inside the loop meant verifying the same answer 967
    # times for a smoke cohort.
    database_dir = Path(database) if database else _bakta_database(config)

    for sample in manifest:
        if reuse_mode != "off":
            decision = decide_reuse(
                sample, out_root=out_root, database_dir=database_dir
            )
            provenance[sample.sample_id] = decision
            if decision.reused:
                records, disagreement = _records_from_bakta_tables(
                    decision.feature_tsv, decision.gff, sample.sample_id
                )
                if records is not None:
                    standardised[sample.sample_id] = records
                    continue
                # decide_reuse ran the same conversion and reached the same
                # verdict, so this is unreachable in practice. Recorded rather
                # than asserted, because the cost of being wrong here is
                # silently reusing a table this stage would have refused.
                failed.append(
                    f"{sample.sample_id}={REUSE_GFF_TSV_DISAGREE}: {disagreement}"
                )
                continue
            if reuse_mode == "require":
                # Named, and only for this sample. Bakta is NOT invoked: that is
                # the whole point of `require`, and the alternative - falling
                # back to the tool - would spend an hour of CPU re-annotating a
                # cohort the run was configured not to annotate.
                failed.append(
                    f"{sample.sample_id}={decision.reason}: {decision.detail}"
                )
                continue
            LOGGER.info(
                "Stage 2: %s - %s; running Bakta for it",
                sample.sample_id, decision.detail,
            )
        try:
            result = run_bakta(
                sample,
                config=config,
                genomes_dir=genomes_dir,
                out_root=out_root,
                database=database_dir,
                run_key=run_key,
                store=store,
                event_sink=event_sink,
                tool_version=tool_version,
                database_version=database_version,
                config_hash=config_hash(config),
            )
        except PipelineError as exc:
            # **A missing assembly is a cohort fact, not a stage failure.**
            #
            # A smoke run's manifest is deliberately the whole PDC roster, so an
            # isolate with no prepared assembly is expected rather than
            # exceptional - `discover_run_manifest` says it "refuses per sample
            # later, rather than silently disappearing". Stages 1, 3 and 6 all
            # honour that; this loop did not, because `run_bakta` *raises* from
            # `locate_assembly` before returning, so the `result.state` check
            # below never saw it and one absent genome ended the cohort.
            #
            # Recorded as `no_assembly` rather than as a Bakta state, because
            # "there was nothing to annotate" and "Bakta failed on a real
            # assembly" are different events. A cohort in which the tool failed
            # on every genome it was actually given should not read the same as
            # a bounded run that correctly skipped the rest.
            failed.append(f"{sample.sample_id}=no_assembly")
            continue
        if result.state != "SUCCEEDED":
            failed.append(f"{sample.sample_id}={result.state}")
            continue
        if reuse_mode != "off":
            provenance[sample.sample_id] = _probe_reused_output(
                sample, out_root=out_root, database_dir=database_dir
            )
        records, disagreement = _records_from_bakta_tables(
            result.outputs.annotation_tsv, result.outputs.gff, sample.sample_id
        )
        if records is None:
            failed.append(f"{sample.sample_id}={disagreement}")
            continue
        standardised[sample.sample_id] = records
    _write_standardised(standardised_root, standardised)
    write_reuse_provenance(out_root, provenance)
    if failed:
        # Recorded, not fatal. Every reason here means "not annotated"
        # downstream - `load_from_intermediate` already represents that as an
        # empty list, so every cohort member keeps its key and only the
        # annotation is missing. Raising here instead would make a bounded run
        # impossible, since a smoke cohort is *mostly* unprepared by
        # construction.
        #
        # The three counts partition `failed` and sum to it. They have to: a
        # count that overlaps another is a reader who cannot add the row up.
        no_assembly = sum(1 for f in failed if f.endswith("=no_assembly"))
        refused = sum(1 for f in failed if _failure_reason(f) in REUSE_REFUSAL_REASONS)
        LOGGER.warning(
            "Stage 2: %d of %d cohort members were not annotated (%d no "
            "assembly, %d annotated but unusable, %d reusable output refused): %s",
            len(failed), len(manifest),
            no_assembly,
            len(failed) - no_assembly - refused,
            refused,
            ",".join(failed[:10]) + ("..." if len(failed) > 10 else ""),
        )
    return load_from_intermediate(standardised_root, manifest)


def _failure_reason(entry: str) -> str:
    """The `sample_id=<reason>` token of one ``failed`` entry."""
    return entry.split("=", 1)[1].split(":", 1)[0]


def _write_standardised(
    root: Path, records: Dict[str, List[AnnotationRecord]]
) -> None:
    """Write the standardised annotation table for each genome.

    Placed exactly where :func:`load_from_intermediate` reads it, so the
    loader is unchanged and there is one code path downstream.
    """
    from ..io.tsv import write_tsv

    directory = Path(root) / "annotation"
    directory.mkdir(parents=True, exist_ok=True)
    for sample_id, rows in records.items():
        write_tsv(
            directory / f"{sample_id}.annotation.tsv",
            (r.to_row() for r in rows),
            list(ANNOTATION_COLUMNS),
        )


def _bakta_database(config: PipelineConfig) -> Path:
    """Resolve the pinned Bakta database directory from configuration.

    Read from ``annotation.bakta_db`` and never downloaded. A missing
    setting is an error, not a default, because silently running Bakta
    against the wrong database is exactly the failure this pipeline is
    built to prevent.
    """
    import os

    raw = dict(getattr(config, "raw", {}) or {})
    section = dict(raw.get("annotation") or getattr(config, "annotation", {}) or {})
    # The environment wins so a machine can point at its own database
    # without the repository carrying anyone's absolute path.
    configured = (
        os.environ.get("BAKTA_DB")
        or section.get("bakta_db")
        or section.get("database_path")
    )
    if not configured:
        raise StageError(
            "REAL-mode annotation needs a pinned Bakta database. Set "
            "annotation.bakta_db in the science config to a directory provisioned "
            "out of band; databases are never downloaded during a run.",
            stage="annotation",
        )

    # Resolving a path is not the same as it being the right database. The
    # preflight is what notices that the directory on disk is a different
    # release - or a different edition - from the one the contract pins, which
    # is the failure docs/design/07-bakta-optimization.md records as the lesson
    # to encode rather than rediscover at stage 2.
    from ..adapters import bakta_db as bakta_preflight

    bakta_preflight.preflight_database(
        database_dir=Path(configured),
        references_path=config.references_path(),
    )
    return Path(configured)
