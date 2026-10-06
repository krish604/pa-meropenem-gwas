"""K1 — `ResultsSource`: the read-only adapter over a run's artefacts.

Two implementations, `LiveSource` and `BundleSource`, chosen by **where
`run_manifest.json` landed** rather than by which probe matched (DESIGN §2.2).
That is what makes the bundle layout tolerant: a bundle missing
`02_stage_outputs/` but carrying the manifest at its root is still served
correctly, because the adapter is a function of the manifest's location.

Two rules are load-bearing (DESIGN §2.1):

1. **Nothing raises for an absent artefact.** Every getter returns the
   `report_tables.TableRead` shape — `present`, `rows`, `reason` — because
   `read_table` already established that a reader which crashes on a stage that
   did not run is worse than one that says the stage did not run.
2. **Paths come from `contracts`, never from literals.** `STAGE_TABLES` and
   `INTERNAL_TABLES` supply every stage filename. The four Snakefile-declared
   paths come from `report_tables.SNAKEFILE_ONLY_TABLES`, which names the
   declaring file for each. A path typed here would be a second answer to a
   question `contracts.py` already answers.

**Trap 1 — `read_tsv` raises on two committed fixtures.** The `recombination`
and `similarity` stand-ins repeat `standin_not_computed` as three column names
and `io/tsv.py` refuses a duplicate header. Everything below therefore reads
through `read_table`'s shape, which catches `PipelineError` and returns
`present=False` with the exception's text as the reason. It surfaces as a
named `not produced`, never as a 500.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from papipeline.execution.contracts import (
    INTERNAL_TABLES,
    STAGE_TABLES,
    internal_table_path,
    required_columns,
    table_path,
)
from papipeline.models import RunMode
from papipeline.stages.report_tables import (
    SNAKEFILE_ONLY_TABLES,
    TableRead,
    read_table,
)
from papipeline.stages.similarity import UNITS_SIDECAR_SUFFIX, units_sidecar_path

#: Environment variables the probe order reads.
ENV_RESULTS_ROOT = "PA_DASH_RESULTS_ROOT"
ENV_BUNDLE = "PA_DASH_BUNDLE"

#: The bundle directory names (DESIGN §12, assumptions A1–A4). Written down
#: here rather than scattered, and every one of them is optional: a bundle that
#: omits `03_report/` still serves its stage tables.
BUNDLE_STAGE_DIR = "02_stage_outputs"
BUNDLE_REPORT_DIR = "03_report"
BUNDLE_RUN_INFO_DIR = "04_run_info"
BUNDLE_VALIDATION_DIR = "05_validation"
BUNDLE_RUNBOOK_DIR = "06_for_900_isolates"
RUNBOOK_NAME = "RUNBOOK_900.md"

#: The **actual delivery layout** a REAL run was handed over in, which differs
#: from DESIGN §12's A1–A5 assumption: the stage tables are under
#: `artifacts/stage_tables/`, the logs under `artifacts/logs/`, the folded-step
#: tables under `artifacts/intermediate/`, and the guard evidence under
#: `guards/`. There is no `run_manifest.json` (written only on success) and no
#: `status/events.jsonl`, so per-stage state is derived from the logs by
#: `stage_logs.derive`. Both layouts are probed; neither is guessed.
DELIVERY_STAGE_DIR = Path("artifacts") / "stage_tables"
DELIVERY_INTERMEDIATE_DIR = Path("artifacts") / "intermediate"
DELIVERY_LOG_DIR = Path("artifacts") / "logs"
DELIVERY_GUARDS_DIR = "guards"
DELIVERY_PROVENANCE_DIR = "provenance"

#: Inside a bundle's `02_stage_outputs/`, the stage tables may be directly there
#: (A2 says it is the bundle's copy of `intermediate/stages/`), under a `stages/`
#: directory, or under an `intermediate/stages/` that mirrors the live tree. All
#: three are probed, in that order, so an unmatched layout fails loudly rather
#: than half-reading.
BUNDLE_STAGE_SUBDIRS: Tuple[Path, ...] = (
    Path("."),
    Path("stages"),
    Path("intermediate") / "stages",
)

#: Where the run's status event log lives. `emit.py`'s default is the
#: repository-relative `status/events.jsonl`; a results-root-relative one is
#: also probed because a bundle may have collected it.
EVENT_LOG_RELATIVE = Path("status") / "events.jsonl"

#: The exact phrase a caller must type to open a REAL run (K4).
REAL_RUN_PHRASE = "run real samples"


@dataclass(frozen=True)
class Probe:
    """One entry in the probe search, and why it did or did not match.

    Every probe is recorded, not only the ones that nearly matched, because a
    failure message that names only the near misses is a failure message that
    hides the answer (DESIGN §2.3).

    `note` carries the extra sentence when `reason` is one of the short enum
    values and the detail would otherwise be lost. `reason` stays inside the
    openapi enum; `note` is additive and documented.
    """

    n: int
    probe: str
    path: Optional[Path]
    ok: bool
    reason: str
    note: str = ""

    def as_dict(self) -> Dict[str, Any]:
        return {
            "n": self.n,
            "probe": self.probe,
            "path": str(self.path) if self.path is not None else None,
            "ok": self.ok,
            "reason": self.reason,
            "note": self.note,
        }


#: The reasons a probe may carry. Declared so a typo is visible here rather
#: than in a response, and so the set can be checked against the openapi enum.
#:
#: `no run_manifest.json and no intermediate/stages/` and `present but holds no
#: pipeline output` are kept apart: the first names two specific absences under a
#: directory that exists, the second says the directory is there and holds
#: nothing recognisable. Neither is collapsed into "not found", which would hide
#: the difference between *not configured* (probes 1-3) and *configured and
#: empty* (probes 4-6).
PROBE_REASONS = (
    "not given",
    "not set",
    "no such directory",
    "no run_manifest.json and no intermediate/stages/",
    "present but holds no pipeline output",
    "matched",
)


@dataclass(frozen=True)
class NoResultsSource(Exception):
    """Nothing matched. Carries every probe, in order."""

    probes: Tuple[Probe, ...]
    message: str = ""

    def __str__(self) -> str:  # pragma: no cover - exercised through the API
        return self.message or self.render()

    def render(self) -> str:
        """The §2.3 message, verbatim in shape.

        Probes 1–3 report *absence of configuration*; probes 4–6 report *the
        directory exists and holds nothing*. Those are different facts and are
        never collapsed into "not found".
        """
        lines = [
            f"No results source could be opened. "
            f"{len(self.probes)} probes were tried, in order:"
        ]
        for probe in self.probes:
            where = str(probe.path) if probe.path is not None else "-"
            lines.append(f"  {probe.n}. {probe.probe:<31}: {probe.reason}")
            if probe.path is not None:
                lines.append(f"       {where}")
            if probe.note:
                lines.append(f"       note: {probe.note}")
        lines.append(
            "No run has been executed in this worktree, or its outputs are "
            "elsewhere. Point the dashboard at one with --results-root or "
            f"{ENV_RESULTS_ROOT}."
        )
        return "\n".join(lines)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "message": self.render(),
            "probes": [p.as_dict() for p in self.probes],
        }


# ---------------------------------------------------------------------------
# The adapter
# ---------------------------------------------------------------------------


@dataclass
class _BaseSource:
    """Shared behaviour. Subclasses supply the stage directory and the extras."""

    root: Path
    stage_dir: Path

    # -- run manifest ----------------------------------------------------
    @property
    def _manifest_path(self) -> Path:
        return self.root / "run_manifest.json"

    def _load_manifest(self) -> Tuple[Optional[Mapping[str, Any]], str]:
        """The parsed manifest, or None with the sentence saying why.

        **Two functions write this filename with different schemas**
        (DESIGN §4). `manifest_writer` is `run` when the payload carries
        `stages` (`papipeline/run.py:2464`) and `provenance` when it does not
        (`scripts/common/write_provenance.py:101`). A manifest with no
        `stages` key does not mean sixteen stages did not run — it means the
        pre-run writer produced it.
        """
        path = self._manifest_path
        if not path.exists():
            return None, (
                f"no run_manifest.json at {path}. Nothing records the run's "
                f"mode, its stage outcomes or its tool versions, so the "
                f"dashboard is reporting artefacts with no run record behind "
                f"them."
            )
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            return None, f"{path} is unreadable as JSON: {exc}"
        if not isinstance(payload, Mapping):
            return None, (
                f"{path} is not a JSON object. It was written by something "
                f"other than this pipeline's two manifest writers."
            )
        return payload, ""

    # -- contracted tables ----------------------------------------------
    def stage_table_path(self, stage: str) -> Path:
        """One stage's principal table path, without reading it.

        The paged endpoints (`/api/stages/{stage}/rows`,
        `/api/tables/{key}`) must never materialise a huge table just to learn
        its byte offsets (UI-D5), so they resolve the path here and presence
        through `peek_table` instead of `read_table`.
        """
        return table_path(self.stage_dir, stage)

    def table_path_for(self, key: str) -> Path:
        """Any declared table's path by key, without reading it."""
        if key in STAGE_TABLES:
            return table_path(self.stage_dir, key)
        if key in INTERNAL_TABLES:
            _filename, _declared = INTERNAL_TABLES[key]
            return internal_table_path(self.stage_dir, key)
        for name, filename, _declared_by in SNAKEFILE_ONLY_TABLES:
            if name == key:
                return self.stage_dir / filename
        return Path(f"<{key}: not a declared table>")

    def stage_table(self, stage: str) -> TableRead:
        """One stage's principal table, or the reason there is not one.

        Never raises: a stage name the pipeline does not have is a different
        answer from a stage whose table is missing, and both are answers rather
        than exceptions.
        """
        if stage not in STAGE_TABLES:
            return TableRead(
                stage,
                Path(f"<{stage}: not a stage in STAGE_ORDER>"),
                reason=(
                    f"{stage!r} is not one of the {len(STAGE_TABLES)} stages "
                    f"papipeline/run.py STAGE_ORDER declares. A typo'd stage "
                    f"name must not read as a stage that did not run."
                ),
            )
        return read_table(
            stage,
            self.stage_table_path(stage),
            required_columns=required_columns(stage),
        )

    def internal_table(self, name: str) -> TableRead:
        """One folded step's table. `structural_variants`, `regulators`,
        `mechanisms`, `master_table`."""
        if name not in INTERNAL_TABLES:
            return TableRead(
                name,
                Path(f"<{name}: not a folded step INTERNAL_TABLES declares>"),
                reason=(
                    f"no output is declared for {name!r} in "
                    f"papipeline/execution/contracts.py INTERNAL_TABLES, so "
                    f"there is no path to read it from."
                ),
            )
        _filename, declared = INTERNAL_TABLES[name]
        return read_table(
            name,
            internal_table_path(self.stage_dir, name),
            required_columns=declared,
        )

    def snakefile_table(self, key: str) -> TableRead:
        """One path the Snakefile declares and `contracts.py` does not."""
        for name, filename, declared_by in SNAKEFILE_ONLY_TABLES:
            if name == key:
                read = read_table(name, self.stage_dir / filename)
                if not read.present:
                    # The reason names the declaring file, so a reader can tell
                    # a declared path from an invented one.
                    return TableRead(
                        name,
                        read.path,
                        reason=f"{read.reason} (path declared by {declared_by})",
                    )
                return read
        declared = ", ".join(name for name, _f, _d in SNAKEFILE_ONLY_TABLES)
        return TableRead(
            key,
            Path(f"<{key}: not a Snakefile-declared table>"),
            reason=(
                f"no output is declared for {key!r}. The Snakefile-declared "
                f"paths this pipeline knows are: {declared}."
            ),
        )

    def table_by_key(self, key: str) -> TableRead:
        """Any declared table by key: a stage, a folded step, or a Snakefile path."""
        if key in STAGE_TABLES:
            return self.stage_table(key)
        if key in INTERNAL_TABLES:
            return self.internal_table(key)
        return self.snakefile_table(key)

    def declared_keys(self) -> Tuple[str, ...]:
        """Every table key this source can serve, in a stable order."""
        return tuple(STAGE_TABLES) + tuple(INTERNAL_TABLES) + tuple(
            name for name, _f, _d in SNAKEFILE_ONLY_TABLES
        )

    # -- non-contracted artefacts ---------------------------------------
    def optional(self, key: str) -> Optional[Path]:
        """A non-contracted artefact's path, or None.

        None means "the dashboard found no such path". It is never a
        fabricated zero and never a guess at a location that is not declared.
        """
        for candidate in self._optional_paths(key):
            if candidate.exists():
                return candidate
        return None

    def _optional_paths(self, key: str) -> Sequence[Path]:
        return ()

    def newick_paths(self) -> List[Path]:
        """Every tree the source can see."""
        return []

    def report_files(self) -> List[Path]:
        """Report-shaped files: `.md`, `.html`, and the figure extensions."""
        return []


#: Extensions that make a file a "report" for the listing. Discovered by
#: extension, not by name, so a bundle that renames its report still lists it
#: (DESIGN §12 A9).
REPORT_SUFFIXES = frozenset({".md", ".html", ".htm", ".svg", ".png", ".pdf"})

#: Newick suffixes, so a `.tre` file is a tree and not a stray text file.
NEWICK_SUFFIXES = frozenset({".nwk", ".newick", ".tre"})


def _discover(root: Path, suffixes: frozenset) -> List[Path]:
    """Every file under ``root`` whose suffix is in ``suffixes``.

    Read-only: `os.walk` and a name test. Sorted so two calls agree, because a
    listing whose order changes between requests is a listing nobody can
    screenshot.
    """
    found: List[Path] = []
    if not root.is_dir():
        return found
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        dirnames.sort()
        for name in sorted(filenames):
            path = Path(dirpath) / name
            if path.suffix.lower() in suffixes:
                found.append(path)
    return found


# ---------------------------------------------------------------------------
# Live layout
# ---------------------------------------------------------------------------


@dataclass
class LiveSource(_BaseSource):
    """A run's own results tree: `<root>/intermediate/stages/`."""

    kind: str = field(default="live", init=False)

    def _optional_paths(self, key: str) -> Sequence[Path]:
        stage = self.stage_dir
        intermediate = self.root / "intermediate"
        if key == "run_manifest":
            return (self.root / "run_manifest.json",)
        if key == "figure_data":
            return (stage / "figure_data",)
        if key == "phylogeny_dir":
            return (intermediate / "phylogeny",)
        if key == "similarity_units":
            # A MEASURED DISCREPANCY between two authorities, so both are
            # probed and neither is silently preferred:
            #
            #  - `similarity.units_sidecar_path(matrix)` derives the sidecar from
            #    the matrix's own *name*, giving `similarity.tsv.units.json`
            #    (`UNITS_SIDECAR_SUFFIX = ".units.json"`, `similarity.py:109`).
            #  - DESIGN §3.4 and the openapi both name
            #    `intermediate/stages/similarity.units.json`.
            #
            # Which one a run wrote depends on the writer that ran, so both are
            # probed in that order and the one that exists is used verbatim. The
            # units JSON is passed through unmodified either way, so a wrong
            # pick could not change a number - only the path reported with it.
            matrix = table_path(stage, "similarity")
            return (
                units_sidecar_path(matrix),
                matrix.with_suffix(UNITS_SIDECAR_SUFFIX),
                intermediate / "similarity.units.json",
            )
        if key == "variants_provenance":
            return (
                stage / "variants_provenance.json",
                intermediate / "variants" / "variants_provenance.json",
                self.root / "variants_provenance.json",
            )
        if key == "reports_dir":
            return (self.root / "reports", self.root.parent / "reports")
        if key == "events_log":
            return (self.root / EVENT_LOG_RELATIVE,)
        if key == "tripwire_log":
            return (
                self.root / "logs" / "tripwire.log",
                self.root / "intermediate" / "tripwire.log",
            )
        if key == "watcher_log":
            return (
                self.root / "logs" / "watcher.log",
                self.root / "intermediate" / "watcher.log",
            )
        if key == "runbook":
            return ()
        if key == "annotation_reuse":
            return (
                intermediate / "annotation" / "reuse_provenance.tsv",
                self.root / "intermediate" / "reuse_provenance.tsv",
                stage / "reuse_provenance.tsv",
            )
        return ()

    def newick_paths(self) -> List[Path]:
        phylo = self.optional("phylogeny_dir")
        return _discover(Path(phylo), NEWICK_SUFFIXES) if phylo else []

    def report_files(self) -> List[Path]:
        reports = self.optional("reports_dir")
        found = _discover(Path(reports), REPORT_SUFFIXES) if reports else []
        # A bundle or a redirected run may keep the report inside the results
        # tree instead of beside it (`reports_root` follows the redirect only
        # when one is set). Probed, never assumed.
        if not found:
            found = _discover(self.root / "reports", REPORT_SUFFIXES)
        return found


# ---------------------------------------------------------------------------
# Bundle layout (DESIGN §12, A1–A5)
# ---------------------------------------------------------------------------


@dataclass
class BundleSource(_BaseSource):
    """A delivery bundle: `01_bakta_input/`, `02_stage_outputs/`, … .

    Every directory is optional and every one is probed. A bundle that omits
    `03_report/` still serves its stage tables; a bundle whose
    `02_stage_outputs/` holds neither stage tables directly nor an
    `intermediate/stages/` mirror is **not** half-read — the probe fails and
    the adapter for that layout is refused (K1).
    """

    kind: str = field(default="bundle", init=False)
    stage_probe: Path = Path(".")
    #: `"assumed"` is DESIGN §12's A1–A5 layout (`02_stage_outputs/`); the
    #: `"delivery"` layout is the one a REAL run actually arrived in
    #: (`artifacts/stage_tables/`, `artifacts/logs/`, `guards/`). The two share
    #: `stage_dir` semantics but not the non-contracted artefact locations.
    layout: str = "assumed"

    def _optional_paths(self, key: str) -> Sequence[Path]:
        if self.layout == "delivery":
            return self._delivery_paths(key)
        root = self.root
        stages = self.stage_dir
        if key == "run_manifest":
            return (
                root / BUNDLE_STAGE_DIR / "run_manifest.json",
                root / BUNDLE_RUN_INFO_DIR / "run_manifest.json",
                root / "run_manifest.json",
            )
        if key == "figure_data":
            return (stages / "figure_data",)
        if key == "phylogeny_dir":
            return (
                root / "intermediate" / "phylogeny",
                stages / "phylogeny",
                root / BUNDLE_STAGE_DIR / "intermediate" / "phylogeny",
            )
        if key == "similarity_units":
            matrix = table_path(stages, "similarity")
            return (
                units_sidecar_path(matrix),
                matrix.with_suffix(UNITS_SIDECAR_SUFFIX),
            )
        if key == "variants_provenance":
            return (
                stages / "variants_provenance.json",
                root / BUNDLE_STAGE_DIR / "intermediate" / "variants" / "variants_provenance.json",
                root / BUNDLE_RUN_INFO_DIR / "variants_provenance.json",
            )
        if key == "reports_dir":
            return (
                root / BUNDLE_REPORT_DIR,
                root / BUNDLE_STAGE_DIR / "reports",
                root / "reports",
            )
        if key == "events_log":
            return (
                root / BUNDLE_RUN_INFO_DIR / EVENT_LOG_RELATIVE,
                root / EVENT_LOG_RELATIVE,
            )
        if key in ("tripwire_log", "watcher_log"):
            name = f"{key[:-4]}.log"
            return (
                root / BUNDLE_RUN_INFO_DIR / "logs" / name,
                root / "logs" / name,
            )
        if key == "runbook":
            return (
                root / BUNDLE_RUNBOOK_DIR / RUNBOOK_NAME,
                root / RUNBOOK_NAME,
            )
        if key == "annotation_reuse":
            return (
                root / BUNDLE_STAGE_DIR / "intermediate" / "annotation" / "reuse_provenance.tsv",
                root / BUNDLE_RUN_INFO_DIR / "reuse_provenance.tsv",
                stages / "reuse_provenance.tsv",
            )
        return ()

    def _delivery_paths(self, key: str) -> Sequence[Path]:
        """The non-contracted artefacts of the actual delivery layout.

        Every path is probed, never assumed: the bundle's `README.md` names what
        it holds and each probe here corresponds to one line of it.
        """
        root = self.root
        stages = self.stage_dir
        intermediate = root / DELIVERY_INTERMEDIATE_DIR
        logs = root / DELIVERY_LOG_DIR
        guards = root / DELIVERY_GUARDS_DIR
        if key == "run_manifest":
            # Written only on success; a failed run leaves none. Probed anyway.
            return (root / "run_manifest.json", intermediate / "run_manifest.json")
        if key == "figure_data":
            return (stages / "figure_data",)
        if key == "phylogeny_dir":
            return (intermediate / "phylogeny", root / "intermediate" / "phylogeny")
        if key == "similarity_units":
            matrix = table_path(stages, "similarity")
            return (
                units_sidecar_path(matrix),
                matrix.with_suffix(UNITS_SIDECAR_SUFFIX),
            )
        if key == "variants_provenance":
            return (
                intermediate / "variants_provenance.json",
                intermediate / "variants" / "variants_provenance.json",
                stages / "variants_provenance.json",
            )
        if key == "reports_dir":
            # The bundle's own prose lives at its root: README.md, GUARDS.md,
            # PER_STAGE_OUTCOMES.md. Discovered by extension, never by name.
            return (root,)
        if key == "events_log":
            # No `status/events.jsonl` in this bundle; the per-stage record is
            # the logs below, read by `stage_logs.derive`.
            return (root / EVENT_LOG_RELATIVE, logs / "events.jsonl")
        if key == "tripwire_log":
            # GUARD 2's tripwire log is `guards/calls.log` (0 bytes = 0
            # invocations). The prior run's non-empty log is a different file
            # and is never read here.
            return (guards / "calls.log", logs / "tripwire.log")
        if key == "watcher_log":
            return (guards / "watcher.log", logs / "watcher.log")
        if key == "runbook":
            return (root / RUNBOOK_NAME,)
        if key == "annotation_reuse":
            return (
                intermediate / "reuse_provenance.tsv",
                intermediate / "annotation" / "reuse_provenance.tsv",
                stages / "reuse_provenance.tsv",
            )
        if key == "full_run_log":
            return (logs / "full_run.log",)
        if key == "stage_log_dir":
            return (logs,)
        return ()

    def newick_paths(self) -> List[Path]:
        phylo = self.optional("phylogeny_dir")
        return _discover(Path(phylo), NEWICK_SUFFIXES) if phylo else []

    def report_files(self) -> List[Path]:
        reports = self.optional("reports_dir")
        return _discover(Path(reports), REPORT_SUFFIXES) if reports else []


# ---------------------------------------------------------------------------
# The probe order (DESIGN §2.2)
# ---------------------------------------------------------------------------


def _probe_flag(given: Optional[Path]) -> Probe:
    """Probe 1 — `--results-root`."""
    if given is None:
        return Probe(1, "--results-root", None, False, "not given")
    if not given.exists():
        return Probe(1, "--results-root", given, False, "no such directory")
    if not given.is_dir():
        return Probe(1, "--results-root", given, False, "present but holds no pipeline output")
    return Probe(1, "--results-root", given, True, "matched")


def _probe_env(name: str, n: int, label: str) -> Probe:
    """Probes 2–3 — an environment variable. `not set` is a distinct fact
    from "set and holds nothing"."""
    raw = os.environ.get(name)
    if not raw:
        return Probe(n, label, None, False, "not set")
    path = Path(raw).expanduser()
    if not path.exists():
        return Probe(n, label, path, False, "no such directory")
    if not path.is_dir():
        return Probe(n, label, path, False, "present but holds no pipeline output")
    return Probe(n, label, path, True, "matched")


def _probe_live(mode: RunMode, n: int, label: str, config: Any) -> Probe:
    """Probes 4–6 — a live results root for one mode.

    Acceptance, per §2.2: `run_manifest.json` or `intermediate/stages/`
    exists. Nothing here is a literal path; the root comes from
    `PipelineConfig.results_root`, which also honours `PIPELINE_RESULTS_ROOT`
    and appends the lowercased mode.
    """
    try:
        path = Path(config.results_root(mode))
    except Exception as exc:  # a config that will not resolve is a refusal
        return Probe(
            n, label, None, False, f"no such directory (config would not resolve: {exc})"
        )
    if not path.exists():
        return Probe(n, label, path, False, "no such directory")
    if not path.is_dir():
        return Probe(n, label, path, False, "present but holds no pipeline output")
    if (path / "run_manifest.json").exists() or (path / "intermediate" / "stages").is_dir():
        return Probe(n, label, path, True, "matched")
    return Probe(
        n, label, path, False, "no run_manifest.json and no intermediate/stages/"
    )


def _resolve_live_stage_dir(path: Path) -> Optional[Path]:
    candidate = path / "intermediate" / "stages"
    return candidate if candidate.is_dir() else None


def _resolve_delivery_stage_dir(path: Path) -> Optional[Path]:
    """Where the actual delivery layout keeps its stage tables.

    `artifacts/stage_tables/` holds the stage and folded-step tables directly
    (A2 is wrong about the directory name, not about the shape). Accepted when
    that directory exists and holds at least one declared stage table **or**
    the run's `artifacts/logs/full_run.log`, so a run that failed before writing
    a table is still openable and can be reported honestly.
    """
    candidate = (path / DELIVERY_STAGE_DIR)
    if not candidate.is_dir():
        return None
    for _stage, (filename, _cols) in STAGE_TABLES.items():
        if (candidate / filename).exists():
            return candidate.resolve()
    if (path / DELIVERY_LOG_DIR / "full_run.log").is_file():
        return candidate.resolve()
    return None


def _resolve_bundle_stage_dir(path: Path) -> Optional[Tuple[Path, Path]]:
    """Where the bundle's stage tables are, and which candidate matched.

    Tries `02_stage_outputs/` itself first (A2: it is the bundle's copy of
    `intermediate/stages/`), then `02_stage_outputs/stages/`, then the
    `02_stage_outputs/intermediate/stages/` mirror. Returns None when none holds
    a declared table, so an unmatched layout fails loudly rather than
    half-reading (K1).
    """
    for sub in BUNDLE_STAGE_SUBDIRS:
        stage_dir = (path / BUNDLE_STAGE_DIR / sub).resolve()
        if not stage_dir.is_dir():
            continue
        for _stage, (_filename, _cols) in STAGE_TABLES.items():
            if (stage_dir / _filename).exists():
                return stage_dir, sub
    return None


def detect(
    *,
    flag_root: Optional[Path] = None,
    bundle_flag: Optional[Path] = None,
    config: Any = None,
) -> Tuple[_BaseSource, List[Probe], Dict[str, Any]]:
    """Run the probe order and open the first source that matches.

    Returns:
        `(source, probes, detail)` where `probes` is every probe tried, in
        order, and `detail` carries the manifest and the extras the endpoints
        need.

    Raises:
        NoResultsSource: nothing matched. The exception renders every probe,
            in order, each with why it failed (DESIGN §2.3).
    """
    probes: List[Probe] = []

    # 1. --results-root
    p1 = _probe_flag(flag_root)
    probes.append(p1)

    # 2. PA_DASH_RESULTS_ROOT
    p2 = _probe_env(ENV_RESULTS_ROOT, 2, ENV_RESULTS_ROOT)
    probes.append(p2)

    # 3. bundle - PA_DASH_BUNDLE, or the flag's sibling directory.
    bundle_candidate: Optional[Path] = None
    if os.environ.get(ENV_BUNDLE):
        bundle_candidate = Path(os.environ[ENV_BUNDLE]).expanduser()
        p3 = _probe_env(ENV_BUNDLE, 3, f"bundle ({ENV_BUNDLE})")
    elif bundle_flag is not None:
        bundle_candidate = Path(bundle_flag).expanduser()
        p3 = (
            Probe(3, "bundle (--bundle)", bundle_candidate, True, "matched")
            if bundle_candidate.is_dir()
            else Probe(
                3,
                "bundle (--bundle)",
                bundle_candidate,
                False,
                "present but holds no pipeline output",
            )
        )
    elif flag_root is not None:
        # "the flag's sibling": a bundle handed over next to the results root.
        sibling = Path(flag_root).expanduser().parent
        if (sibling / BUNDLE_STAGE_DIR).is_dir():
            bundle_candidate = sibling
            p3 = Probe(3, "bundle (sibling of --results-root)", sibling, True, "matched")
        else:
            p3 = Probe(3, "bundle (sibling of --results-root)", None, False, "not given")
    else:
        p3 = Probe(3, f"bundle ({ENV_BUNDLE})", None, False, "not set")
    probes.append(p3)

    # 4-6. live roots, one per mode.
    #
    # Recorded even when no configuration could be loaded. The probe order is
    # six steps and the failure message has to name all of them: a message that
    # listed three because the other three could not be *attempted* would read
    # as "there were only three places to look", which is false.
    live_probes: List[Tuple[int, str, RunMode]] = [
        (4, "live REAL", RunMode.REAL),
        (5, "live TEST", RunMode.TEST),
        (6, "live STUB", RunMode.STUB),
    ]
    for n, label, mode in live_probes:
        if config is None:
            probes.append(
                Probe(
                    n,
                    label,
                    None,
                    False,
                    "not set",
                    note=(
                        "the machine overlay could not be loaded, so "
                        "`paths.results_root` was unavailable to resolve "
                        "against; point the dashboard at a root with "
                        "--results-root"
                    ),
                )
            )
        else:
            probes.append(_probe_live(mode, n, label, config))

    # --- acceptance, in probe order ------------------------------------
    # Probes 1 and 2 are roots given by a person: open them as a LIVE tree,
    # because a live tree is `<root>/intermediate/stages`. If they hold a
    # bundle the manifest-location rule below still decides the kind.
    for probe in probes:
        if not probe.ok or probe.path is None:
            continue
        path = probe.path
        if probe.probe.startswith("bundle"):
            # Both bundle layouts are accepted at the bundle probe: DESIGN
            # §12's assumed `02_stage_outputs/` and the actual delivery
            # `artifacts/stage_tables/`. The assumed one is tried first so an
            # existing bundle is served exactly as before.
            resolved = _resolve_bundle_stage_dir(path)
            if resolved is not None:
                stage_dir, sub = resolved
                source = BundleSource(root=path.resolve(), stage_dir=stage_dir)
                return _finish(source, probes, {"bundle_probe": str(sub)})
            delivery = _resolve_delivery_stage_dir(path)
            if delivery is not None:
                source = BundleSource(
                    root=path.resolve(), stage_dir=delivery, layout="delivery"
                )
                return _finish(source, probes, {"bundle_probe": "delivery"})
            # Configured as a bundle but not shaped like one. Recorded as
            # a failed probe rather than silently opened as a live tree.
            probes[probe.n - 1] = Probe(
                probe.n,
                probe.probe,
                path,
                False,
                "present but holds no pipeline output",
            )
            continue

        # A person-given root (--results-root / PA_DASH_RESULTS_ROOT). The
        # actual delivery layout is recognised before the live tree, because a
        # bundle handed over as `--results-root` is still a bundle.
        delivery = _resolve_delivery_stage_dir(path)
        if delivery is not None:
            source = BundleSource(
                root=path.resolve(), stage_dir=delivery, layout="delivery"
            )
            return _finish(source, probes, {"bundle_probe": "delivery"})

        stage_dir = _resolve_live_stage_dir(path)
        manifest_at_root = (path / "run_manifest.json").exists()
        if stage_dir is None and not manifest_at_root:
            probes[probe.n - 1] = Probe(
                probe.n,
                probe.probe,
                path,
                False,
                "present but holds no pipeline output",
            )
            continue
        if stage_dir is None:
            stage_dir = (path / "intermediate" / "stages").resolve()
        source = LiveSource(root=path.resolve(), stage_dir=stage_dir)
        return _finish(source, probes, {})

    raise NoResultsSource(tuple(probes))


def _finish(
    source: _BaseSource, probes: List[Probe], detail: Dict[str, Any]
) -> Tuple[_BaseSource, List[Probe], Dict[str, Any]]:
    """Attach the manifest and decide `kind` by where the manifest landed.

    §2.2: `live` if `run_manifest.json` is at the root, `bundle` if it is under
    `02_stage_outputs/`. A bundle that is missing `02_stage_outputs` and has
    the manifest at its root is therefore served as `live` — correctly, because
    that is where its tables are too.
    """
    manifest, reason = source._load_manifest()
    manifest_path = source._manifest_path
    if manifest is None and (source.root / BUNDLE_STAGE_DIR / "run_manifest.json").exists():
        # A bundle whose manifest sits under `02_stage_outputs/`. Read it there
        # and mark the source a bundle.
        alt = source.root / BUNDLE_STAGE_DIR / "run_manifest.json"
        try:
            payload = json.loads(alt.read_text(encoding="utf-8"))
            manifest = payload if isinstance(payload, Mapping) else None
            if manifest is not None:
                reason = ""
                manifest_path = alt
                source = BundleSource(root=source.root, stage_dir=source.stage_dir)
        except (OSError, ValueError) as exc:
            reason = f"{alt} is unreadable as JSON: {exc}"
    elif isinstance(source, BundleSource) and (source.root / "run_manifest.json").exists():
        # The other half of the §2.2 rule: a bundle-shaped tree whose manifest
        # is at the **root** is served as `live`. The probe order is a search,
        # but the kind is a function of where the manifest landed, so a bundle
        # that is missing `02_stage_outputs/` and has the manifest at its root
        # is still served correctly — and reported as `live`, which is where its
        # tables are too.
        source = LiveSource(root=source.root, stage_dir=source.stage_dir)

    writer = manifest_writer(manifest)
    # No manifest is not "no record". A failed run leaves its per-stage state
    # in `artifacts/logs/full_run.log` and one `stage_<name>.log` per stage it
    # reached; derive it there so a stage that ran is never rendered `not_run`
    # with the pre-run-manifest sentence. None when there is nothing to derive.
    from . import stage_logs

    derived = stage_logs.derive(
        source.root, source.stage_dir, layout=getattr(source, "layout", "")
    )
    detail.update(
        {
            "run_manifest": manifest,
            "manifest_reason": reason,
            "manifest_path": manifest_path,
            "manifest_writer": writer,
            "derived_stage_state": derived.as_dict() if derived is not None else None,
        }
    )
    return source, probes, detail


def manifest_writer(manifest: Optional[Mapping[str, Any]]) -> str:
    """Which of the two writers produced this manifest.

    `run` when it carries `stages` (`papipeline/run.py:2464`); `provenance`
    when it does not (`scripts/common/write_provenance.py:101`); `none` when
    there is no manifest at all.
    """
    if not isinstance(manifest, Mapping):
        return "none"
    if "stages" in manifest:
        return "run"
    if "references" in manifest or "stage_tool_requirements" in manifest:
        return "provenance"
    return "none"


def manifest_gaps(manifest: Optional[Mapping[str, Any]], writer: str) -> List[str]:
    """Keys this writer did not produce.

    Shown so a manifest without `stages` reads as "the pre-run writer produced
    this" rather than as "no stages ran".
    """
    if writer == "none" or not isinstance(manifest, Mapping):
        return [
            "no manifest is readable, so the run's mode, its stage outcomes and "
            "its tool versions are not recorded",
        ]
    expected = {
        "run": (
            "pipeline_version", "run_mode", "generated_at_utc", "antibiotic",
            "python", "platform", "stages", "stages_skipped", "outputs",
            "tools_detected", "references",
        ),
        "provenance": (
            "pipeline_version", "run_mode", "generated_at_utc", "python",
            "platform", "references", "stage_tool_requirements",
            "unpinned_references", "note",
        ),
    }[writer]
    return [f"no {key} key" for key in expected if key not in manifest]


__all__ = [
    "BUNDLE_REPORT_DIR",
    "BUNDLE_RUN_INFO_DIR",
    "BUNDLE_RUNBOOK_DIR",
    "BUNDLE_STAGE_DIR",
    "BUNDLE_VALIDATION_DIR",
    "DELIVERY_GUARDS_DIR",
    "DELIVERY_INTERMEDIATE_DIR",
    "DELIVERY_LOG_DIR",
    "DELIVERY_PROVENANCE_DIR",
    "DELIVERY_STAGE_DIR",
    "ENV_BUNDLE",
    "ENV_RESULTS_ROOT",
    "EVENT_LOG_RELATIVE",
    "BundleSource",
    "LiveSource",
    "NoResultsSource",
    "Probe",
    "REAL_RUN_PHRASE",
    "RUNBOOK_NAME",
    "detect",
    "manifest_gaps",
    "manifest_writer",
]