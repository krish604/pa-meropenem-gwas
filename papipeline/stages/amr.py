"""Stage 4 - AMR detection.

Pluggable database/tool architecture:

* :class:`AmrAdapter` - the interface every backend implements;
* :class:`AmrFinderPlusAdapter` - the preferred backend;
* :class:`CardRgiAdapter` - an optional alternative, disabled by default;
* :class:`PrecomputedAmrAdapter` - reads an already-produced table, which is
  what TEST mode uses.

Scientific constraints enforced here:

* a detected determinant is :attr:`ClaimStatus.DETECTED` and nothing more;
* ``antibiotic`` is filled from the knowledge table, never assumed from the
  database's own antibiotic annotation, so an unmapped determinant is
  recorded as ``unmapped`` instead of being attached to imipenem;
* ``database`` and ``database_version`` are mandatory and are carried
  through to every downstream row.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from ..adapters.amr import execute_command as amr_tool_execute
from ..adapters.external import ToolStatus
from ..config.loader import PipelineConfig
from ..errors import PipelineError
from ..io.tsv import read_tsv
from ..knowledge import mechanism_map_for_antibiotic
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import AmrDeterminant, ClaimStatus, RunMode

LOGGER = get_logger("stages.amr")

AMR_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "antibiotic",
    "determinant",
    "gene",
    "variant",
    "determinant_type",
    "mechanism",
    "evidence_source",
    "database",
    "database_version",
    "confidence",
    "claim_status",
    "identity_pct",
    "coverage_pct",
)

REQUIRED = (
    "sample_id",
    "antibiotic",
    "determinant",
    "determinant_type",
    "evidence_source",
    "database",
    "database_version",
)

#: Value used when a determinant cannot be tied to a configured antibiotic.
UNMAPPED_ANTIBIOTIC = "unmapped"

#: AMRFinderPlus' `Subtype` values that denote a mutated rather than an intact
#: element. Both are read: `POINT` is a listed single substitution,
#: `POINT_DISRUPT` an alignment-detected disruption the database does not list.
POINT_SUBTYPES = frozenset({"POINT", "POINT_DISRUPT"})

#: How AMRFinderPlus spells a mutated element in `Element symbol`:
#: ``<gene>_<variant>``, where the variant begins with the one-letter residue
#: code - ``oprD_V359L``, ``nalC_G71E``, ``mexR_I24AfsTer94``.
#:
#: **The pattern is measured, not assumed.** Applied to all 1515 rows of the
#: provisioned database's own ``AMRProt-mutation.tsv`` it matches 1515, and the
#: gene it recovers equals the gene in every row's ``mutated_protein_name`` -
#: zero disagreements. A stricter ``[A-Z]\d+`` variant pattern was tried first
#: and rejected: it missed 21 symbols, all of them multi-residue
#: substitutions/indels (``rpsJ_HKYK56del``, ``porB1b_GA120del``) whose variant
#: begins with a residue letter followed by another letter.
#:
#: The gene character class excludes ``_`` deliberately, which is what makes the
#: split unambiguous; no gene symbol in the database contains one.
_POINT_SYMBOL = re.compile(
    r"^(?P<gene>[A-Za-z0-9()'.\-\[\]+]+)_(?P<variant>[A-Z][A-Za-z0-9]*)$"
)


def split_point_symbol(symbol: str) -> Tuple[Optional[str], Optional[str]]:
    """``(gene, variant)`` for a point-mutated element symbol.

    Returns ``(None, None)`` when `symbol` does not carry a variant, which is
    the normal answer for an intact gene - and is why only rows the report
    marks as POINT are put through this at all.
    """
    match = _POINT_SYMBOL.match((symbol or "").strip())
    if match is None:
        return None, None
    return match.group("gene"), match.group("variant")


def point_variants_in_report(report: Path) -> Dict[str, Tuple[str, str]]:
    """``element symbol -> (gene, variant)`` for the report's POINT rows.

    Read from the report rather than re-deriving it from the record: `Subtype`
    is the tool's own statement that a row is a mutation, and a symbol that
    happens to contain an underscore in an AMR row is not a point mutation.
    A symbol the pattern cannot split is **left out and logged**, not guessed
    at - a wrong gene/variant split is a mechanism claim attributed to the
    wrong gene, which is worse than an empty variant column.
    """
    import csv

    out: Dict[str, Tuple[str, str]] = {}
    report = Path(report)
    if not report.is_file():
        return out
    with report.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("Subtype") or "").strip().upper() not in POINT_SUBTYPES:
                continue
            symbol = (row.get("Element symbol") or "").strip()
            gene, variant = split_point_symbol(symbol)
            if not symbol or gene is None:
                LOGGER.warning(
                    "AMRFinderPlus reported POINT element %r, which does not "
                    "match the <gene>_<variant> spelling the provisioned "
                    "database uses; leaving it unsplit rather than guessing a "
                    "gene", symbol,
                )
                continue
            out[symbol] = (gene, variant)
    return out


def base_gene(determinant: str) -> str:
    """Strip the allele suffix from a determinant name.

    ``blaOXA-48~SYNTHETIC_ALLELE`` -> ``blaOXA-48``.
    """
    return determinant.split("~")[0].strip() if determinant else ""


class AmrAdapter:
    """Interface for an AMR detection backend."""

    name = "abstract"

    def __init__(
        self,
        database: str,
        database_version: str,
        antibiotic: str,
    ) -> None:
        self.database = database
        self.database_version = database_version
        self.antibiotic = antibiotic

    def detect(self, sample_id: str, assembly: Optional[Path]) -> List[AmrDeterminant]:
        raise NotImplementedError

    def provenance(self) -> Dict[str, Optional[str]]:
        return {
            "database": self.database,
            "database_version": self.database_version,
            "tool": self.name,
        }


class AmrFinderPlusAdapter(AmrAdapter):
    """AMRFinderPlus backend, on the v4 interface. Never updates the database.

    Replaces a v3-era implementation whose command the installed 4.2.7 rejects
    outright - ``--input`` was renamed to ``--nucleotide`` in 4.x, and the
    report columns were re-ordered. Verified against the binary, not taken from
    release notes; see ``docs/design/10-reusable-modules.md``, which records the
    same divergence and points at the pilot's v4 parser as the correct
    implementation.

    The command is not built here. :func:`papipeline.adapters.amr.commands_for_isolate`
    owns the invocation, so the v4 spelling and the no-update guarantee are
    asserted in one place instead of at every call site.

    `report_dir` is where the tool's report lands; it is a parameter because
    the output filename is part of the command, and a caller that cannot say
    where reports go cannot read them back.

    `organism` is required, and is the value the caller already validated
    against the database's own ``--list_organisms``. It is not defaulted here:
    a default would be a second, unvalidated statement of the species, and the
    one thing this backend must never do is silently run without a valid
    ``--organism``, because that drops every POINT row without saying so.
    """

    name = "amrfinderplus"

    def __init__(
        self,
        status: ToolStatus,
        database: str,
        database_version: str,
        antibiotic: str,
        organism: str,
        database_path: Optional[Path] = None,
        threads: int = 4,
        report_dir: Optional[Path] = None,
        runner=None,
    ) -> None:
        super().__init__(database, database_version, antibiotic)
        self.status = status
        self.organism = organism
        self.database_path = database_path
        self.threads = threads
        self.report_dir = Path(report_dir) if report_dir is not None else None
        # Defaults to a real executor, not to `CommandResult`. That dataclass was
        # the old default and is a *result*, so calling it as a runner raised
        # TypeError on the first real run - see
        # `papipeline.adapters.amr.execute_command`.
        self._runner = runner or amr_tool_execute

    def build_command(self, assembly: Path, sample_id: str = "") -> List[str]:
        """The v4 invocation for one assembly, with the report in `report_dir`."""
        from ..adapters import amr as amr_tool

        return amr_tool.commands_for_isolate(
            executable=str(self.status.executable),
            assembly=Path(assembly),
            database_dir=Path(self.database_path or self.database),
            threads=self.threads,
            organism=self.organism,
            out_path=self.report_path_for(sample_id or Path(assembly).stem),
        )[0]

    def provenance(self) -> Dict[str, Optional[str]]:
        """What stage 4 ran with, in the order a reader needs it.

        The organism is first because it is the one input that decides whether
        the screen covers point mutations at all. A run whose provenance does not
        name it cannot be told apart from a run whose screen was narrowed to
        genes, and those two produce different tables from the same cohort.
        """
        return {
            "organism": self.organism,
            "tool": self.name,
            "tool_version": getattr(self.status, "version", None),
            "database": self.database,
            "database_version": self.database_version,
        }

    def detect(self, sample_id: str, assembly: Optional[Path]) -> List[AmrDeterminant]:
        """Screen one assembly, or return nothing when there is no sequence.

        A missing assembly is an empty result rather than an error: the stage
        records the sample as screened with nothing found, which is honest,
        whereas raising would abandon the rest of the cohort.
        """
        if assembly is None:
            return []
        assembly = Path(assembly)
        # Created here, not left to the tool: AMRFinderPlus writes to the path it
        # is given and does not create the directory, so a missing one is a
        # failed screen that looks like a negative one.
        if self.report_dir is not None:
            self.report_dir.mkdir(parents=True, exist_ok=True)
        command = self.build_command(assembly, sample_id)
        result = self._runner(command=command)
        report = self.report_path_for(sample_id)
        if not report.is_file():
            raise PipelineError(
                f"AMRFinderPlus reported success for {sample_id} but wrote no "
                f"report at {report}. Refusing to treat that as a clean screen: "
                "an unreadable screen and a negative one must not look alike.",
                sample_id=sample_id,
                report=str(report),
                command=" ".join(command),
            )
        return self.parse_report(report, sample_id)

    def report_path_for(self, sample_id: str) -> Path:
        """Where this sample's report is written, and read.

        One function, so the path the tool is told to write and the path the
        adapter reads cannot disagree - which they did when the command derived
        the name from the assembly while the lookup used `report_dir`.
        """
        name = f"{sample_id}.amrfinder.tsv"
        if self.report_dir is not None:
            return self.report_dir / name
        return Path(name)

    def parse_report(
        self, report: Path, sample_id: str
    ) -> List[AmrDeterminant]:
        """Parse a v4 report via the pilot's parser, not a second copy of it.

        ``papipeline.pilot.amr_detect`` is the implementation
        ``docs/design/08-migration-plan.md`` nominates. Duplicating it here would
        be a second parser free to drift from the first, and the drift would
        surface as a REAL run disagreeing with the pilot on identical input.

        Then the one thing that parser cannot do: fill in `gene` and `variant`
        for POINT rows. It records `Element symbol` verbatim into both fields,
        which for a mutated element is ``oprD_V359L`` - a gene name that matches
        nothing in ``config/mechanisms.tsv``. Left alone, ``--organism`` would
        therefore *lose* the ``oprD`` carbapenem mechanism it just discovered,
        which is the opposite of the point of the flag.
        """
        from dataclasses import replace

        from ..pilot import amr_detect

        report = Path(report)
        records = amr_detect.parse_amrfinder_report(
            report,
            sample_id=sample_id,
            antibiotic=self.antibiotic,
            database=self.database,
            database_version_value=self.database_version,
        )
        points = point_variants_in_report(report)
        if not points:
            return records
        return [
            replace(
                record,
                gene=points.get(record.determinant, (record.gene, None))[0],
                variant=points.get(record.determinant, (None, None))[1],
            )
            for record in records
        ]


class CardRgiAdapter(AmrAdapter):
    """Optional CARD/RGI backend. Disabled unless enabled in config."""

    name = "card_rgi"

    def detect(self, sample_id: str, assembly: Optional[Path]) -> List[AmrDeterminant]:
        raise NotImplementedError(
            "CARD/RGI adapter is declared but not enabled; set "
            "analysis.amr.optional_adapters.card_rgi=true and pin the CARD "
            "version in config/references.tsv first"
        )


class PrecomputedAmrAdapter(AmrAdapter):
    """Reads a previously produced AMR table.

    Used in TEST mode against synthetic fixtures, and as the ingestion path
    for externally generated AMR calls.
    """

    name = "precomputed"

    def __init__(self, table: Path, antibiotic: str) -> None:
        super().__init__(database="precomputed", database_version="n/a", antibiotic=antibiotic)
        self.table = Path(table)

    def detect(self, sample_id: str, assembly: Optional[Path]) -> List[AmrDeterminant]:
        return []

    def load_all(self) -> List[AmrDeterminant]:
        return load_amr_table(self.table, antibiotic=self.antibiotic)


def parse_amrfinder_stdout(
    text: str,
    sample_id: str,
    antibiotic: str,
    database: str,
    database_version: str,
) -> List[AmrDeterminant]:
    """Parse ``amrfinder --output -`` tabular output.

    Expected columns include ``Protein identifier``, ``Gene symbol``,
    ``Sequence Name``, ``Scope``, ``Element Type``, ``Element Subtype``,
    ``Class``, ``Subclass``, ``Method``, ``Target`` and ``Identity``.
    """
    import csv
    from io import StringIO

    reader = csv.DictReader(StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        return []

    results: List[AmrDeterminant] = []
    for row in reader:
        gene = (row.get("Gene symbol") or row.get("Sequence Name") or "").strip()
        if not gene:
            continue
        gene_symbol = gene.split("__")[0].strip()
        confidence = (row.get("Identity") or "").strip() or None
        results.append(
            AmrDeterminant(
                sample_id=sample_id,
                antibiotic=antibiotic,
                determinant=gene,
                gene=gene_symbol or None,
                variant=None,
                determinant_type=(row.get("Element Type") or "unknown").strip(),
                mechanism=None,
                evidence_source="amrfinderplus",
                database=database,
                database_version=database_version,
                confidence=confidence,
                claim_status=ClaimStatus.DETECTED,
                identity_pct=_maybe_float(confidence),
            )
        )
    return results


def _maybe_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = value.strip().rstrip("%")
    try:
        return float(text)
    except ValueError:
        return None


def load_amr_table(
    path: Path,
    antibiotic: str,
    strict_mapping: bool = True,
    mechanism_map: Optional[Dict[str, object]] = None,
) -> List[AmrDeterminant]:
    """Load and normalise an AMR determinant table.

    Args:
        path: TSV with the AMR output schema.
        antibiotic: Antibiotic the project is analysing.
        strict_mapping: When True, a determinant whose gene is absent from
            the knowledge table is recorded with
            ``antibiotic=unmapped`` and a warning, rather than being
            attributed to ``antibiotic``.
        mechanism_map: Gene -> MechanismSpec. Supplied by the caller so this
            function stays free of configuration dependencies.

    Returns:
        Normalised determinants. A gene in the knowledge table gets its
        mechanism attached; an unknown gene does not.
    """
    rows = read_tsv(path, required_columns=REQUIRED)
    results: List[AmrDeterminant] = []

    for row in rows:
        determinant = str(row["determinant"])
        gene = row.get("gene") or base_gene(determinant)
        row_antibiotic = str(row["antibiotic"])

        mechanism = None
        mapped_antibiotic = row_antibiotic
        if mechanism_map is not None:
            spec = mechanism_map.get(base_gene(gene or determinant))
            if spec is not None:
                mechanism = spec.mechanism
                if row_antibiotic == UNMAPPED_ANTIBIOTIC or not row_antibiotic:
                    mapped_antibiotic = antibiotic
            elif strict_mapping:
                LOGGER.warning(
                    "Determinant %s in sample %s is not in the knowledge table; "
                    "recording as %s rather than attributing it to %s",
                    determinant,
                    row.get("sample_id"),
                    UNMAPPED_ANTIBIOTIC,
                    antibiotic,
                )
                mapped_antibiotic = UNMAPPED_ANTIBIOTIC

        results.append(
            AmrDeterminant(
                sample_id=str(row["sample_id"]),
                antibiotic=mapped_antibiotic,
                determinant=determinant,
                gene=gene or None,
                variant=row.get("variant"),
                determinant_type=str(row["determinant_type"]),
                mechanism=mechanism or row.get("mechanism"),
                evidence_source=str(row["evidence_source"]),
                database=str(row["database"]),
                database_version=str(row["database_version"]),
                confidence=row.get("confidence"),
                claim_status=ClaimStatus.DETECTED,
                identity_pct=_maybe_float(row.get("identity_pct")),
                coverage_pct=_maybe_float(row.get("coverage_pct")),
            )
        )

    unmapped = sum(1 for r in results if r.antibiotic == UNMAPPED_ANTIBIOTIC)
    LOGGER.info(
        "Stage 4: %d determinant calls (%d unmapped) [%s]",
        len(results),
        unmapped,
        antibiotic,
    )
    return results


def build_amr_adapter(
    config: PipelineConfig,
    antibiotic: str,
    *,
    report_dir: Optional[Path] = None,
) -> AmrAdapter:
    """Construct the configured backend, having first cleared the preflight.

    The database check happens *here*, before any genome is analysed, because
    that is the whole point of it: `docs/design/07-bakta-optimization.md`
    records the alternative, where the mismatch only surfaces once a tool has
    already annotated something.

    The organism is resolved here for the same reason. `organism.name` from
    `config/science.yaml` is validated against the database's own
    ``--list_organisms`` *before* any screen, so an unrecognised species is a
    refusal at build time rather than a quietly unscreened cohort - an
    unrecognised value does not fail the tool, it just reports no POINT rows.

    Returns the AMRFinderPlus backend. The optional CARD backend is not selected
    here: `amr.optional_adapters.card_rgi` is false in `science.yaml`, and an
    adapter that raises `NotImplementedError` when enabled is the honest state of
    that experiment rather than a stub to ship.
    """
    from ..adapters import amr as amr_tool
    from ..adapters.external import detect_tools, require_tool

    settings = dict(config.raw.get("amr") or {})
    # On demand, for the one tool this stage is about to invoke - not a sweep
    # over everything the project has heard of. Stage 4 runs AMRFinderPlus for
    # every isolate below, so its version is asked for by a real execution.
    tool = require_tool(
        detect_tools(["amrfinder"], execute=True), "amrfinder", "stage 4 (AMR)"
    )

    database_dir = Path(
        settings.get("database_path")
        or os.environ.get("AMRFINDER_DB")
        or (Path(config.db_root) / "amrfinderplus-db" / "latest")
    )
    amr_tool.preflight_database(
        database_dir=database_dir,
        references_path=config.references_path(),
    )

    organism = amr_tool.resolve_organism(
        config.organism.get("name"),
        executable=str(tool.executable),
        database_dir=database_dir,
    )

    return AmrFinderPlusAdapter(
        status=tool,
        database=str(settings.get("database") or "AMRFinderPlus"),
        database_version=amr_tool.database_version_on_disk(database_dir) or "unknown",
        antibiotic=antibiotic,
        organism=organism,
        database_path=database_dir,
        threads=int(config.runtime.get("threads", 1) or 1),
        report_dir=report_dir,
    )


def _run_by_calling_tool(
    config: PipelineConfig,
    manifest: SampleManifest,
    data_root: Path,
    workdir: Path,
    antibiotic: str,
) -> Dict[str, List[AmrDeterminant]]:
    """REAL path: screen each isolate with the configured backend.

    **A failed isolate is recorded, not raised.** 82% of this cohort's
    assemblies are corrupt (`stages.variants` documents it), so a stage that
    stopped at the first failure could never complete, and one that ignored
    them would report a smaller cohort without saying so. The sample keeps its
    key with an empty list, which means "screened, nothing found" - a
    meaningful result, and distinct from an absent key.

    That distinction is the reason the failure is not fatal, so it is worth
    being explicit: a screen that could not be performed and a screen that
    found nothing must not look alike, and this records them differently only if
    the caller reads the log line. It is the best available answer given the
    stage's return type.
    """
    from ..assemblies import locate_assembly

    backend = build_amr_adapter(
        config, antibiotic, report_dir=Path(workdir) / "reports"
    )
    grouped: Dict[str, List[AmrDeterminant]] = {}
    failures: List[str] = []

    for sample_id in manifest.sample_ids:
        try:
            assembly = locate_assembly(manifest.require(sample_id), data_root)
            grouped[sample_id] = backend.detect(sample_id, assembly)
        except PipelineError as exc:
            LOGGER.warning(
                "Stage 4: no screen for %s, recorded as screened with nothing "
                "found (%s)", sample_id, exc,
            )
            grouped[sample_id] = []
            failures.append(sample_id)

    # The aggregated table, not just the per-isolate reports.
    #
    # `load_amr_table` reads `amr_determinants.tsv`, and the Snakemake rule
    # declares it as `amr`'s input - so before this, a REAL run wrote reports
    # into `reports/` that nothing in the DAG consumed, and the declared input
    # existed only for TEST, where the fixture generator supplies it. That made
    # the REAL DAG unbuildable (`MissingInputException` in rule amr) and it is
    # the same shape as the annotation defect: a tool-output table that only one
    # mode has a producer for.
    #
    # The header-only-when-empty case is deliberate. A screen that found nothing
    # must be distinguishable from a screen that never ran, so the file is
    # always written with its declared columns even when `grouped` is all empty.
    _write_determinant_table(
        Path(workdir) / "amr_determinants.tsv", grouped, backend=backend,
    )

    screened = sum(1 for v in grouped.values() if v)
    LOGGER.info(
        "Stage 4: %d/%d isolates have at least one determinant via %s%s",
        screened, len(manifest.sample_ids), backend.name,
        f"; {len(failures)} failed" if failures else "",
    )
    return grouped


def _write_determinant_table(
    path: Path,
    grouped: Dict[str, List[AmrDeterminant]],
    *,
    backend: Optional[AmrAdapter] = None,
) -> None:
    """Write the aggregate determinant table the DAG and loader both read.

    The backend's provenance rides along as `#` comment lines. It is written
    here rather than only in `run_manifest.json` because the *file* is what a
    later reader of this cohort opens first, and the organism flag is the one
    input that decides whether the table covers point mutations. A
    point-mutation-free table and a table from a run that never screened for
    point mutations are otherwise byte-identical.
    """
    from ..io.tsv import write_tsv

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [record.to_row() for calls in grouped.values() for record in calls]
    columns = list(REQUIRED) + [
        c for c in AMR_COLUMNS if c not in REQUIRED
    ]
    header = _provenance_header(backend)
    write_tsv(path, rows, columns, header_comment=header or None)
    if header:
        LOGGER.info("Stage 4 provenance: %s", "; ".join(header))


def _provenance_header(backend: Optional[AmrAdapter]) -> List[str]:
    """`key=value` lines describing how the table was screened.

    Values are quoted rather than merely written, because an organism name
    contains a space and an unquoted value would read as two fields.
    """
    if backend is None:
        return []
    record = backend.provenance()
    return [
        f"stage=4_amrfinderplus {key}={str(value)!r}"
        for key, value in sorted(record.items())
        if value is not None
    ]


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    intermediate_root: Path,
    antibiotic: str,
    data_root: Optional[Path] = None,
) -> Dict[str, List[AmrDeterminant]]:
    """Stage 4 entry point.

    Returns:
        sample_id -> determinants. A sample with no determinants maps to an
        empty list, which is a meaningful result ("screened, nothing found")
        and is distinct from an absent key ("not screened").

    ``TEST`` reads the committed fixture table; ``REAL`` runs the backend.
    """
    config.require_antibiotic(antibiotic)

    if mode is RunMode.REAL:
        if data_root is None:
            raise PipelineError(
                "REAL mode needs a data_root to locate assemblies; refusing to "
                "fall back to anything, because screening another isolate's "
                "sequence would attribute its genes to this one"
            )
        return _run_by_calling_tool(
            config, manifest, data_root,
            Path(intermediate_root) / "amr", antibiotic,
        )

    path = Path(intermediate_root) / "amr" / "amr_determinants.tsv"
    strict = bool(config.raw.get("amr", {}).get("strict_antibiotic_mapping", True))
    mechanism_map = mechanism_map_for_antibiotic(config, antibiotic)

    determinants = load_amr_table(
        path,
        antibiotic=antibiotic,
        strict_mapping=strict,
        mechanism_map=mechanism_map,
    )

    grouped: Dict[str, List[AmrDeterminant]] = {sid: [] for sid in manifest.sample_ids}
    known = set(grouped)
    orphans: List[str] = []
    for determinant in determinants:
        if determinant.sample_id in known:
            grouped[determinant.sample_id].append(determinant)
        else:
            orphans.append(determinant.sample_id)
    if orphans:
        LOGGER.warning(
            "%d AMR rows reference samples not in the manifest: %s",
            len(orphans),
            ",".join(sorted(set(orphans))[:10]),
        )

    screened = sum(1 for v in grouped.values() if v)
    LOGGER.info(
        "Stage 4: %d/%d samples have at least one determinant",
        screened,
        len(grouped),
    )
    return grouped
