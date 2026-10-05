"""K2 — the FastAPI app: all 29 paths from `dashboard/openapi.yaml`.

Shape of the module:

- `AppState` holds the opened results source, the state directory, the bind
  policy and the config-derived minimum sample size. It is built once at
  startup and read by every handler.
- A `require_source` dependency turns "no results source" into the §2.3 503
  rather than letting each route invent its own absence.
- Paging, sorting and filtering go through `query.parse_page` /
  `query.apply_page`, so every list endpoint honours §5.1 identically.
- Every statistic that leaves the server carries a `basis` (UI-D7) and the
  configured minimum sample size, so the D3 flag fires from `basis.n` rather
  than from whatever the page was told.

Nothing here writes to the results root. `ResultsSource` has no write method,
and every path a client names goes through `security.resolve_contained`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from fastapi import FastAPI, Query, Request
from fastapi.responses import JSONResponse, PlainTextResponse, StreamingResponse

from papipeline.models import RunMode
from papipeline.stages.report_tables import NOT_PRODUCED, TableRead

from . import badges as badge_module
from . import events as event_module
from . import isolates as isolate_module
from . import launcher as launcher_module
from . import newick as newick_module
from . import oprd as oprd_module
from . import provenance as provenance_module
from . import source as source_module
from .index import IndexStore, peek_table
from .papipeline_refs import (
    POWER_FLAG_MEANING,
    POWER_FLAG_TEMPLATE,
    SPEC_NAMES,
    SPEC_NUMBERS,
    STAGE_ORDER,
    power_flag,
)
from .query import QueryError, apply_page, basis, page_meta, parse_page
from .security import (
    BindPolicy,
    ContainmentError,
    TOKEN_HEADER,
    check_suffix,
    check_token,
    content_headers,
    resolve_bind,
    resolve_contained,
)
from .source import NoResultsSource
from .state import StateDir

#: Above this many tips the similarity matrix is downsampled server-side and the
#: response says so (§5.3).
DEFAULT_MAX_TIPS = 300

#: The p-value array a QQ plot needs is capped here; the full set is 900 rows
#: of nine columns and the browser only draws the quantiles.
QQ_MAX = 5000


@dataclass
class AppState:
    """Everything a handler needs, resolved once."""

    state_dir: StateDir
    bind: BindPolicy
    repo_root: Optional[Path] = None
    allow_launch: bool = False
    machine: Optional[str] = None
    min_samples: Optional[int] = None
    antibiotic: Optional[str] = None

    source: Optional[source_module._BaseSource] = None
    probes: Tuple[source_module.Probe, ...] = ()
    detail: Dict[str, Any] = field(default_factory=dict)
    failure: Optional[NoResultsSource] = None
    index: Optional[IndexStore] = None
    runner: Optional[launcher_module.Runner] = None

    # -- convenience ----------------------------------------------------
    @property
    def root(self) -> Optional[Path]:
        return self.source.root if self.source is not None else None

    @property
    def stage_dir(self) -> Optional[Path]:
        return self.source.stage_dir if self.source is not None else None

    @property
    def manifest(self) -> Optional[Mapping[str, Any]]:
        return self.detail.get("run_manifest")

    @property
    def manifest_writer(self) -> str:
        return str(self.detail.get("manifest_writer") or "none")

    @property
    def kind(self) -> str:
        return self.source.kind if self.source is not None else "live"

    def log_path(self) -> Optional[Path]:
        if self.source is None:
            return None
        return self.source.optional("events_log") or event_module.find_log(
            self.source.root, self.repo_root
        )[0]


# ---------------------------------------------------------------------------
# App construction
# ---------------------------------------------------------------------------


def build_state(
    *,
    results_root: Optional[Path] = None,
    bundle: Optional[Path] = None,
    repo_root: Optional[Path] = None,
    state_dir: Optional[Path] = None,
    host: Optional[str] = None,
    port: Optional[int] = None,
    allow_launch: bool = False,
    machine: Optional[str] = None,
    runner: Optional[launcher_module.Runner] = None,
) -> AppState:
    """Open the results source, or record the whole failed search.

    A failure is not raised here: the API answers `503` with every probe named,
    and a server that refuses to start cannot report which paths it tried.
    """
    bind = resolve_bind(host, port)
    directory = StateDir.open(state_dir)
    config = None
    try:
        from papipeline.config.loader import load_config

        config = load_config(machine=machine) if (repo_root and (repo_root / "config" / "science.yaml").exists()) else None
    except Exception:
        # A config that will not load removes probes 4-6, not the server. The
        # failure message then names five probes rather than six, which is
        # honest: those three could not be attempted.
        config = None

    state = AppState(
        state_dir=directory,
        bind=bind,
        repo_root=Path(repo_root).resolve() if repo_root else None,
        allow_launch=bool(allow_launch),
        machine=machine,
        index=IndexStore(directory),
        runner=runner,
    )
    if config is not None:
        # Read from the loaded configuration rather than typed: the minimum
        # sample size is a scientific setting in `config/science.yaml`, and a
        # dashboard that hard-coded it would report a threshold the pipeline
        # did not use.
        try:
            state.min_samples = int(config.gwas.min_samples_per_group)
        except Exception:
            state.min_samples = None
        try:
            antibiotics = list(config.antibiotics)
            state.antibiotic = antibiotics[0] if antibiotics else None
        except Exception:
            state.antibiotic = None

    try:
        source, probes, detail = source_module.detect(
            flag_root=results_root, bundle_flag=bundle, config=config
        )
        state.source = source
        state.probes = tuple(probes)
        state.detail = detail
        directory.write_session(
            root=source.root,
            kind=source.kind,
            manifest_writer=state.manifest_writer,
        )
    except NoResultsSource as failure:
        state.failure = failure
        state.probes = failure.probes
        directory.write_session(
            root=None, kind="none", manifest_writer="none"
        )
    return state


def create_app(state: AppState) -> FastAPI:
    """The FastAPI application. No side effects beyond wiring."""
    app = FastAPI(
        title="PA AMR results dashboard API",
        version="0.1.0-phase1",
        description=(
            "Read-only HTTP surface for a Pseudomonas aeruginosa imipenem AMR "
            "pipeline run. Absence is stated, never rendered as a value; every "
            "statistic arrives with the artefact it was computed on."
        ),
        docs_url="/api/docs",
        openapi_url="/api/openapi.json",
    )
    app.state.pa = state

    _read_only_middleware(app)
    _mount_static(app, state)
    _routes(app)
    return app


def _mount_static(app: FastAPI, state: AppState) -> None:
    """Serve the frontend, if it has been written.

    SHELL owns `dashboard/static/`. The backend only mounts it; a missing
    directory is not an error, because the API is useful without a shell.
    """
    static = Path(__file__).resolve().parent.parent / "static"
    if static.is_dir():
        from starlette.staticfiles import StaticFiles

        app.mount("/static", StaticFiles(directory=str(static)), name="static")

        index = static / "index.html"

        @app.get("/", include_in_schema=False)
        def _index() -> Any:
            if index.exists():
                from starlette.responses import FileResponse

                return FileResponse(str(index))
            return JSONResponse(
                {
                    "error": (
                        "The frontend has not been built. This is the API; see "
                        "/api/docs for its own contract."
                    )
                },
                status_code=200,
            )


# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------


def _state(request: Request) -> AppState:
    return request.app.state.pa


def require_source(request: Request) -> AppState:
    """The opened source, or the §2.3 503 naming every probe.

    `check_token` runs here, before any handler, so a non-local bind refuses
    every `/api/**` request without the token — including the ones that would
    otherwise answer.
    """
    state = _state(request)
    try:
        check_token(request.headers.get(TOKEN_HEADER), required=state.bind.token_required)
    except ContainmentError as exc:
        return _error(exc.message, exc.status_code)
    if state.source is None:
        failure = state.failure or NoResultsSource(state.probes)
        return JSONResponse(status_code=503, content=failure.as_dict())  # type: ignore[return-value]
    return state


def _error(message: str, status_code: int = 400, **extra: Any) -> JSONResponse:
    payload: Dict[str, Any] = {"error": message, "detail": extra.get("detail")}
    if "columns" in extra:
        payload["columns"] = extra["columns"]
    if "probes" in extra:
        payload["probes"] = extra["probes"]
    return JSONResponse(status_code=status_code, content=payload)


def _ok_query(fn, *args, **kwargs):
    """Run a handler, turning a `QueryError` into its 400."""
    try:
        return fn(*args, **kwargs)
    except QueryError as exc:
        return _error(exc.message, 400, **exc.detail)
    except ContainmentError as exc:
        return _error(exc.message, exc.status_code)


def _read_only_middleware(app: FastAPI) -> None:
    """UI-D1, enforced rather than assumed.

    `ResultsSource` exposes no write method and every file under the opened
    root is opened `'r'`, so nothing the dashboard does can modify what it
    reports. This middleware is the belt to that braces: a write verb against
    any path other than the two launcher POSTs is refused with a 405 that says
    why, so a future route cannot acquire one by accident.
    """

    @app.middleware("http")
    async def _guard(request: Request, call_next):
        if request.method in ("GET", "HEAD", "OPTIONS"):
            return await call_next(request)
        if request.url.path in ("/api/launcher/preflight",):
            return await call_next(request)
        return JSONResponse(
            status_code=405,
            content={
                "error": (
                    "the dashboard is read-only over the results root (UI-D1). "
                    "The only POST endpoint is /api/launcher/preflight, which "
                    "starts nothing."
                ),
                "detail": None,
            },
        )


# ---------------------------------------------------------------------------
# Stage and table helpers
# ---------------------------------------------------------------------------


def _manifest_state(stage: str, manifest: Optional[Mapping[str, Any]]) -> Optional[str]:
    if not isinstance(manifest, Mapping):
        return None
    stages = manifest.get("stages")
    if not isinstance(stages, Mapping):
        return None
    value = stages.get(stage)
    return str(value) if value is not None else None


def _manifest_skipped(manifest: Optional[Mapping[str, Any]]) -> Dict[str, str]:
    if not isinstance(manifest, Mapping):
        return {}
    skipped = manifest.get("stages_skipped")
    if not isinstance(skipped, Mapping):
        return {}
    return {str(k): str(v) for k, v in skipped.items()}


def _manifest_records_stages(manifest: Optional[Mapping[str, Any]]) -> bool:
    return isinstance(manifest, Mapping) and "stages" in manifest


def _run_mode(manifest: Optional[Mapping[str, Any]]) -> Optional[str]:
    if isinstance(manifest, Mapping):
        value = manifest.get("run_mode")
        if value:
            return str(value)
    return None


def _event_signals(state: AppState) -> Dict[str, Any]:
    """Per-stage `start`/`fail` evidence from the log, keyed on stage name.

    Read once per request and passed to `badges.classify`, which is step 4 of
    the reason order.
    """
    path = state.log_path()
    if path is None or not Path(path).exists():
        return {}
    replay = event_module.read_from(path, 0)
    signals: Dict[str, Any] = {}
    for event in replay.events:
        stage = event.mapped_stage
        if stage is None:
            continue
        entry = signals.setdefault(
            stage,
            {"start": False, "terminal": False, "fail": False, "sentence": ""},
        )
        if event.event == "start":
            entry["start"] = True
        else:
            entry["terminal"] = True
        if event.event == "fail":
            entry["fail"] = True
            entry["sentence"] = (
                f"a fail event at {event.t} names this stage for sample "
                f"{event.sample}; the process did not succeed"
            )
    return signals


def _stage_payload(
    state: AppState, stage: str, *, signals: Optional[Mapping[str, Any]] = None
) -> Dict[str, Any]:
    """One stage's row: badge, reason, table presence, row count, folded steps."""
    source = state.source
    assert source is not None
    manifest = state.manifest
    read = source.stage_table(stage)
    signal = (signals or {}).get(stage, {})
    badge = badge_module.classify(
        stage,
        manifest_state=_manifest_state(stage, manifest),
        table_present=read.present,
        table_reason=read.reason,
        stages_skipped=_manifest_skipped(manifest),
        run_mode=_run_mode(manifest),
        has_start_event=bool(signal.get("start")),
        has_terminal_event=bool(signal.get("terminal")),
        has_fail_event=bool(signal.get("fail")),
        fail_event_sentence=str(signal.get("sentence") or ""),
        manifest_records_stages=_manifest_records_stages(manifest),
    )
    internal: List[Dict[str, Any]] = []
    for key, owner in (("structural_variants", "amr"), ("regulators", "amr"),
                       ("mechanisms", "cooccurrence"), ("master_table", "reporting")):
        if owner != stage:
            continue
        folded = source.internal_table(key)
        internal.append(
            {
                "key": key,
                "path": str(folded.path),
                "owner": owner,
                "present": folded.present,
                "reason": folded.reason,
                "n_rows": folded.n_rows if folded.present else None,
            }
        )
    return {
        "name": stage,
        "spec_name": SPEC_NAMES.get(stage, stage),
        "spec_number": SPEC_NUMBERS.get(stage),
        "order": list(STAGE_ORDER).index(stage) + 1,
        "state": badge.key,
        "reason": badge.reason,
        "badge": badge.as_dict(),
        "present": read.present,
        "path": str(read.path),
        "n_rows": read.n_rows if read.present else None,
        # Advisory. The header the FILE carries is reported separately at
        # `/api/stages/{stage}`, and the two are cross-checked as information
        # rather than as a filter (§3.6.3).
        "declared_columns": _declared_columns(stage),
        "internal_tables": internal,
    }


def _declared_columns(stage: str) -> List[str]:
    from papipeline.execution.contracts import STAGE_TABLES

    entry = STAGE_TABLES.get(stage)
    return list(entry[1]) if entry else []


def _actual_header(path: Path) -> Tuple[List[str], str]:
    """The header the file carries, or the reason it could not be read.

    Cross-checked against the contract as **information**, never as a filter: a
    file whose columns have drifted is still the run's record, and hiding it
    would be worse than flagging it. The committed stand-ins already disagree
    with `STAGE_TABLES`.
    """
    path = Path(path)
    if not path.exists():
        return [], f"no file at {path}"
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            for raw in handle:
                if raw.startswith("#"):
                    continue
                if not raw.strip():
                    continue
                return [c.strip() for c in raw.rstrip("\n").rstrip("\r").split("\t")], ""
    except OSError as exc:
        return [], f"{path} is unreadable: {exc}"
    return [], f"{path} carries no header row"


def _header_drift(declared: Sequence[str], actual: Sequence[str]) -> List[Dict[str, str]]:
    drift: List[Dict[str, str]] = []
    for column in declared:
        if column not in actual:
            drift.append({"column": column, "kind": "missing_from_file"})
    for column in actual:
        if column not in declared:
            drift.append({"column": column, "kind": "missing_from_contract"})
    return drift


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


def _routes(app: FastAPI) -> None:

    # -- 1. /api/health --------------------------------------------------
    @app.get("/api/health")
    def health(request: Request) -> Any:
        state = _state(request)
        base = {
            "ok": state.source is not None,
            "source_kind": state.kind if state.source is not None else None,
            "root": str(state.root) if state.root is not None else None,
            "manifest_writer": state.manifest_writer,
            "mode": _run_mode(state.manifest),
            "pipeline_version": (
                state.manifest.get("pipeline_version")
                if isinstance(state.manifest, Mapping)
                else None
            ),
            "state_dir": str(state.state_dir.root),
            "loopback": state.bind.loopback,
            "token_required": state.bind.token_required,
            "badges": badge_module.badge_map(),
        }
        if state.source is None:
            failure = state.failure or NoResultsSource(state.probes)
            return JSONResponse(status_code=503, content={**failure.as_dict(), **base})
        return base

    # -- 2. /api/run -----------------------------------------------------
    @app.get("/api/run")
    def run_summary(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        manifest = state.manifest
        writer = state.manifest_writer
        signals = _event_signals(state)
        stages = [_stage_payload(state, stage, signals=signals) for stage in STAGE_ORDER]
        n_samples = _n_samples(state, manifest)
        payload: Dict[str, Any] = {
            "manifest_writer": writer,
            "manifest_gaps": source_module.manifest_gaps(manifest, writer),
            "manifest_path": str(state.detail.get("manifest_path") or ""),
            "run_mode": _run_mode(manifest),
            "antibiotic": (
                manifest.get("antibiotic") if isinstance(manifest, Mapping) else state.antibiotic
            ),
            "n_samples": n_samples,
            "generated_at_utc": manifest.get("generated_at_utc") if isinstance(manifest, Mapping) else None,
            "pipeline_version": manifest.get("pipeline_version") if isinstance(manifest, Mapping) else None,
            "python": manifest.get("python") if isinstance(manifest, Mapping) else None,
            "platform": manifest.get("platform") if isinstance(manifest, Mapping) else None,
            "banner": _banner(state),
            "stages_skipped": _manifest_skipped(manifest),
            "outputs": (manifest.get("outputs") if isinstance(manifest, Mapping) else None) or {},
            "stages": stages,
            "probes": [p.as_dict() for p in state.probes],
        }
        return payload

    # -- 3. /api/run/counts ---------------------------------------------
    @app.get("/api/run/counts")
    def run_counts(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        signals = _event_signals(state)
        stages = [_stage_payload(state, stage, signals=signals) for stage in STAGE_ORDER]
        tally: Dict[str, int] = {key: 0 for key in badge_module.BADGES}
        for row in stages:
            tally[row["state"]] = tally.get(row["state"], 0) + 1
        return {
            "n_stages_completed": tally.get("completed", 0),
            "n_stages_failed": tally.get("failed", 0),
            "n_stages_refused": tally.get("refused", 0),
            "n_stages_not_run": tally.get("not_run", 0),
            "n_stages_not_assessed": tally.get("not_assessed", 0),
            "n_stages_running": tally.get("running", 0),
            "n_samples": _n_samples(state, state.manifest),
            "bases": {
                "stages": basis(
                    n=len(stages),
                    artefact="papipeline.run.STAGE_ORDER",
                    path=str(state.root),
                    rows_total=len(stages),
                    min_samples=state.min_samples,
                ),
                "manifest": basis(
                    n=_n_samples(state, state.manifest) or 0,
                    artefact="run_manifest.json",
                    path=str(state.detail.get("manifest_path")),
                    rows_total=len(stages),
                    min_samples=state.min_samples,
                ),
            },
        }

    # -- 4. /api/stages --------------------------------------------------
    @app.get("/api/stages")
    def list_stages(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        signals = _event_signals(state)
        items = [_stage_payload(state, stage, signals=signals) for stage in STAGE_ORDER]
        return {"stage_order": list(STAGE_ORDER), "items": items}

    # -- 5. /api/stages/{stage} -----------------------------------------
    @app.get("/api/stages/{stage}")
    def get_stage(stage: str, request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        if stage not in STAGE_ORDER:
            return _error(
                f"{stage!r} is not a stage in STAGE_ORDER. The {len(STAGE_ORDER)} "
                f"stages are: {', '.join(STAGE_ORDER)}. A typo'd stage name must "
                f"not read as a stage that did not run.",
                404,
                columns=list(STAGE_ORDER),
            )
        signals = _event_signals(state)
        payload = _stage_payload(state, stage, signals=signals)
        source = state.source
        assert source is not None
        read = source.stage_table(stage)
        declared = _declared_columns(stage)
        actual, header_reason = _actual_header(read.path)
        payload.update(
            {
                "declared_header": declared,
                "actual_header": actual,
                "actual_header_reason": header_reason,
                "header_drift": _header_drift(declared, actual),
                "table": {
                    "present": read.present,
                    "reason": read.reason,
                    "path": str(read.path),
                    "n_rows": read.n_rows if read.present else None,
                },
                "events": _events_for(state, stage),
            }
        )
        return payload

    # -- 6. /api/stages/{stage}/rows ------------------------------------
    @app.get("/api/stages/{stage}/rows")
    def stage_rows(
        stage: str,
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        if stage not in STAGE_ORDER:
            return _error(
                f"{stage!r} is not a stage in STAGE_ORDER; the {len(STAGE_ORDER)} "
                f"stages are: {', '.join(STAGE_ORDER)}.",
                404,
                columns=list(STAGE_ORDER),
            )
        return _ok_query(
            _stage_rows_impl, state, stage, request, offset, limit, sort, q,
            columns, include_total, _filter_dict(request),
        )

    # -- 7. /api/tables/{key} -------------------------------------------
    @app.get("/api/tables/{key}")
    def list_table(
        key: str,
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        source = state.source
        assert source is not None
        if key not in source.declared_keys():
            return _error(
                f"{key!r} is not a table this pipeline declares. Known keys: "
                f"{', '.join(source.declared_keys())}.",
                404,
                columns=list(source.declared_keys()),
            )
        return _ok_query(
            _table_rows_impl, state, key, request, offset, limit, sort, q,
            columns, include_total, _filter_dict(request),
        )

    # -- 8. /api/isolates ------------------------------------------------
    @app.get("/api/isolates")
    def list_isolates(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
        ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(
            _isolates_impl,
            state,
            request,
            offset,
            limit,
            sort,
            q,
            columns,
            include_total,
            _filter_dict(request),
        )

    # -- 9. /api/isolates/{sample_id} -----------------------------------
    @app.get("/api/isolates/{sample_id}")
    def get_isolate(sample_id: str, request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_isolate_detail_impl, state, sample_id)

    # -- 10. /api/oprd ---------------------------------------------------
    @app.get("/api/oprd")
    def get_oprd(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_oprd_impl, state, request, offset, limit, _filter_dict(request))

    # -- 11. /api/pangenome ----------------------------------------------
    @app.get("/api/pangenome")
    def get_pangenome(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_pangenome_impl, state)

    # -- 12. /api/tree ----------------------------------------------------
    @app.get("/api/tree")
    def get_tree(request: Request, source: str = Query("stage9")) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_tree_impl, state, source)

    # -- 13. /api/tree/tips ----------------------------------------------
    @app.get("/api/tree/tips")
    def list_tree_tips(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(
            _tree_tips_impl, state, request, offset, limit, sort, q, columns, include_total, _filter_dict(request)
        )

    # -- 14. /api/similarity ---------------------------------------------
    @app.get("/api/similarity")
    def get_similarity(request: Request, max_tips: int = Query(DEFAULT_MAX_TIPS)) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_similarity_impl, state, max_tips)

    # -- 15. /api/gwas ---------------------------------------------------
    @app.get("/api/gwas")
    def list_gwas(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(
            _gwas_impl, state, request, offset, limit, sort, q, columns, include_total, _filter_dict(request)
        )

    # -- 16. /api/gwas/top ----------------------------------------------
    @app.get("/api/gwas/top")
    def gwas_top(request: Request, limit: int = Query(25)) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_gwas_top_impl, state, limit)

    # -- 17. /api/convergence --------------------------------------------
    @app.get("/api/convergence")
    def list_convergence(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(
            _convergence_impl, state, request, offset, limit, sort, q, columns, include_total, _filter_dict(request)
        )

    # -- 18. /api/cooccurrence -------------------------------------------
    @app.get("/api/cooccurrence")
    def list_cooccurrence(
        request: Request,
        offset: int = Query(0),
        limit: int = Query(200),
        sort: Optional[str] = Query(None),
        q: Optional[str] = Query(None),
        columns: Optional[str] = Query(None),
        include_total: int = Query(1),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(
            _cooccurrence_impl, state, request, offset, limit, sort, q, columns, include_total, _filter_dict(request)
        )

    # -- 19. /api/provenance ---------------------------------------------
    @app.get("/api/provenance")
    def get_provenance(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_provenance_impl, state)

    # -- 20. /api/provenance/bakta ---------------------------------------
    @app.get("/api/provenance/bakta")
    def get_bakta(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_bakta_impl, state)

    # -- 21. /api/provenance/logs ---------------------------------------
    @app.get("/api/provenance/logs")
    def provenance_logs(
        request: Request,
        name: str = Query("any"),
        offset: int = Query(0),
        limit: int = Query(200),
        tail_lines: int = Query(500),
    ) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_logs_impl, state, name, offset, limit, tail_lines)

    # -- 22. /api/provenance/digests ------------------------------------
    @app.get("/api/provenance/digests")
    def digests(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_digests_impl, state)

    # -- 23. /api/reports ------------------------------------------------
    @app.get("/api/reports")
    def list_reports(request: Request) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        return _ok_query(_reports_impl, state)

    # -- 24. /api/reports/content ----------------------------------------
    @app.get("/api/reports/content")
    def report_content(request: Request, path: str = Query(...)) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        try:
            resolved = resolve_contained(state.root, path)
            suffix = check_suffix(resolved)
        except ContainmentError as exc:
            return _error(exc.message, exc.status_code)
        try:
            text = resolved.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            return _error(f"{resolved} could not be read: {exc}", 404)
        from .security import CONTENT_TYPES

        return PlainTextResponse(
            text,
            media_type=CONTENT_TYPES.get(suffix, "text/plain; charset=utf-8"),
            headers=dict(content_headers(suffix, resolved.name)),
        )

    # -- 25. /api/events -------------------------------------------------
    @app.get("/api/events")
    async def stream_events(request: Request) -> Any:
        state = _state(request)
        check_token(request.headers.get(TOKEN_HEADER), required=state.bind.token_required)
        if state.source is None:
            failure = state.failure or NoResultsSource(state.probes)
            return JSONResponse(status_code=503, content=failure.as_dict())
        last = request.headers.get("Last-Event-ID")
        try:
            cursor = int(last) if last and str(last).strip().isdigit() else 0
        except ValueError:
            cursor = 0
        snapshot = _event_snapshot(state)
        return StreamingResponse(
            event_module.sse_frames(state.log_path(), initial=snapshot, last_event_id=cursor),
            media_type="text/event-stream",
            headers={
                "Cache-Control": "no-cache",
                "X-Accel-Buffering": "no",
                # The header carries the transport in both cases, so a
                # screenshot of the UI says which mechanism it is showing.
                "X-PA-Dash-Transport": "sse",
            },
        )

    # -- 26. /api/events/state -------------------------------------------
    @app.get("/api/events/state")
    def events_state(request: Request, cursor: int = Query(0)) -> Any:
        state = require_source(request)
        if isinstance(state, JSONResponse):
            return state
        payload = _event_snapshot(state)
        log_path = state.log_path()
        replay = event_module.read_from(log_path, cursor)
        payload.update(
            {
                "transport": "poll",
                "cursor": replay.cursor,
                "n_events": len(replay.events),
                "delta": [e.as_dict() for e in replay.events],
                "n_lines_skipped": replay.n_skipped,
                "skipped_warnings": replay.warnings[:200],
                "log_probes": event_module.find_log(state.root, state.repo_root)[1],
                "truncated": replay.truncated,
            }
        )
        return payload

    # -- 27. /api/launcher/capabilities ---------------------------------
    @app.get("/api/launcher/capabilities")
    def launcher_capabilities(request: Request) -> Any:
        state = _state(request)
        check_token(request.headers.get(TOKEN_HEADER), required=state.bind.token_required)
        return launcher_module.capabilities(
            allow_launch=state.allow_launch, loopback=state.bind.loopback
        )

    # -- 28. /api/launcher/status ----------------------------------------
    @app.get("/api/launcher/status")
    def launcher_status(request: Request) -> Any:
        state = _state(request)
        check_token(request.headers.get(TOKEN_HEADER), required=state.bind.token_required)
        timings, source_name, probes = launcher_module.per_isolate_timings(state.root)
        payload = launcher_module.status(log_path=state.log_path(), timings_source=source_name)
        payload["per_isolate_timings"] = timings
        payload["runbook_probes"] = probes
        return payload

    # -- 29. /api/launcher/preflight -------------------------------------
    @app.post("/api/launcher/preflight")
    async def launcher_preflight(request: Request) -> Any:
        state = _state(request)
        check_token(request.headers.get(TOKEN_HEADER), required=state.bind.token_required)
        try:
            body = await request.json()
        except Exception:
            body = {}
        if not isinstance(body, dict):
            body = {}
        mode = str(body.get("mode") or "dry_run")
        # Only `dry_run` and `real` are recognised. An unknown mode is refused
        # rather than mapped to one, so a typo cannot become a different action.
        if mode not in ("dry_run", "real"):
            return _error(
                f"unknown mode {mode!r}. The launcher recognises 'dry_run' and "
                f"'real' only; nothing is executable besides those.",
                400,
            )
        phrase = body.get("phrase")
        overlay = body.get("overlay") or body.get("machine")
        if phrase is not None and not isinstance(phrase, str):
            return _error("phrase must be a string typed by the caller.", 400)
        payload = launcher_module.preflight(
            allow_launch=state.allow_launch,
            repo_root=state.repo_root,
            results_root=state.root,
            log_path=state.log_path(),
            runner=state.runner,
            mode=mode,
            phrase=phrase,
            overlay=overlay if isinstance(overlay, str) else None,
        )
        launcher_module.log_launch(
            state.state_dir,
            {
                "event": "preflight",
                "mode": mode,
                "would_start": False,
                "phrase_supplied": bool(phrase),
                "overlay": overlay if isinstance(overlay, str) else None,
                "root": str(state.root) if state.root else None,
            },
        )
        return payload


# ---------------------------------------------------------------------------
# Endpoint implementations
# ---------------------------------------------------------------------------


def _filter_dict(request: Request) -> Dict[str, List[str]]:
    """The repeated `filter.*` params, as the multi-dict the parser expects."""
    out: Dict[str, List[str]] = {}
    for key, value in request.query_params.multi_items():
        out.setdefault(key, []).append(value)
    return out


def _banner(state: AppState) -> Optional[str]:
    """The mode banner, verbatim from `stages/reporting.py`.

    A STUB or TEST run must not read as an analysis, so the banner travels with
    the run summary rather than being composed by the client.
    """
    mode = _run_mode(state.manifest)
    try:
        from papipeline.models import RunMode as _Mode

        from papipeline.stages.reporting import banner_for

        return banner_for(_Mode(mode)) if mode else None
    except Exception:
        return None


def _n_samples(state: AppState, manifest: Optional[Mapping[str, Any]]) -> Optional[int]:
    """The cohort size, or None when nothing records it.

    Not 0: `0` standing for unknown is the plausible zero this project refuses
    everywhere.
    """
    if isinstance(manifest, Mapping):
        for key in ("n_samples", "n_isolates", "samples"):
            value = manifest.get(key)
            if isinstance(value, int) and not isinstance(value, bool):
                return value
    # Derived from a per-sample table when the manifest says nothing.
    source = state.source
    if source is None:
        return None
    for stage in ("validation", "mlst", "phenotype"):
        read = source.stage_table(stage)
        if read.present:
            return read.n_rows
    return None


def _power(n: int, min_samples: Optional[int] = None) -> Dict[str, Any]:
    """The R16 flag for a statistic computed on `n` isolates.

    `N` is the cohort the statistic was **actually computed on**, so this is
    called with `basis.n` and never with a cohort size the page was told about.

    The pipeline's `power_flag` reads "n={n}, underpowered, not a finding" for
    every `n` — only the boolean distinguishes. Repeating "underpowered" on an
    adequately powered cohort is misleading, so the flag text is conditional
    here: the pipeline's verbatim template when `n < min_samples`, and the same
    template minus "underpowered" when it is not. The R16 meaning is always
    carried verbatim.
    """
    underpowered = bool(min_samples is not None and int(n) < int(min_samples))
    if underpowered:
        flag = power_flag(n)
    else:
        flag = f"n={n}, not a finding"
    return {
        "flag": flag,
        "template": POWER_FLAG_TEMPLATE,
        "meaning": POWER_FLAG_MEANING,
        "n": n,
        "underpowered": underpowered,
    }


def _page(
    state: AppState,
    rows: Sequence[Mapping[str, Any]],
    *,
    available: Sequence[str],
    searchable: Sequence[str],
    artefact: str,
    path: Optional[Path],
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
    key_fn=None,
    extra: Optional[Mapping[str, Any]] = None,
    basis_n: Optional[int] = None,
    projection_extra: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Parse, apply and wrap — the identical path every list endpoint takes.

    `basis_n` overrides the `n` the basis and power flag are computed on. The
    default is the row count, which is right for a table whose rows *are* the
    statistic. For gwas / convergence / cooccurrence the rows are tests or
    pairs and the statistic was computed on the cohort, so the caller passes
    the cohort size and the R16 flag reads `n=900`, not `n=240` (DESIGN §7.1).
    """
    request = parse_page(
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        available=available,
        searchable=searchable,
    )
    if projection_extra and request.columns:
        # `basis` is the provenance of every cell, not a data column, so a
        # projection must not drop it. It is not in `available` (it is added to
        # the row after parsing), so it is unioned in here rather than passed
        # through `parse_page`'s validation.
        extra_cols = tuple(
            c for c in projection_extra if c not in request.columns
        )
        request.columns = tuple(request.columns) + extra_cols
    items, total, total_unfiltered = apply_page(
        rows, request, searchable=searchable, key_fn=key_fn
    )
    computed_on = total if request.has_filter or request.q else total
    n_for_flag = basis_n if basis_n is not None else computed_on
    payload: Dict[str, Any] = {
        "items": items,
        "meta": page_meta(
            request,
            total=total,
            total_unfiltered=total_unfiltered,
            basis=basis(
                n=n_for_flag,
                artefact=artefact,
                path=path,
                rows_total=total_unfiltered,
                filtered=bool(request.has_filter or request.q),
                min_samples=state.min_samples,
            ),
        ),
        "power": _power(n_for_flag, min_samples=state.min_samples),
    }
    if extra:
        payload.update(dict(extra))
    return payload


def _rows_via_index(
    state: AppState,
    *,
    path: Path,
    key_column: Optional[str],
    page_request,
    searchable: Sequence[str],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Serve a page from byte ranges rather than from a materialised table.

    This is UI-D5's actual mechanism: one streaming pass to learn each row's
    key value and byte range, then `seek` + `readline` for the rows on the page
    only. A 200,000-row file with tens of thousands of columns never becomes a
    list of dicts, and `bytes_read` is returned so a test can prove it.

    Filters and sorting that need a column other than the key fall back to the
    indexed column's value, and say so: a silent wrong-column sort would look
    like a sort that applied.
    """
    from .index import PagedReader

    reader = PagedReader(path)
    # The scan happens once per file, not once per request: the index it
    # produces is cached in the state directory and re-read when its key still
    # matches. A second request for the same table therefore does zero scanning
    # and pays only `seek` + `readline` for the page.
    cached = state.index.get(str(path.parent.name or path.name), path, key_column=key_column)
    header = list(cached.header)
    spans = cached.as_spans()
    scan_bytes = 0
    if not spans:
        header, raw_spans = reader.scan(key_column)
        spans = raw_spans
        scan_bytes = reader.bytes_read
    clauses = page_request.clauses

    unsupported = [c.column for c in clauses if key_column is not None and c.column != key_column]
    if unsupported:
        # A filter on a column the index does not carry is a 400, not a silent
        # no-op. The index stores only the key column's value, so a filter on
        # any other column cannot be answered from it — and a filter that
        # quietly did not apply looks like a filter (the same rule as an
        # unknown sort column).
        raise QueryError(
            f"unsupported filter column(s): {', '.join(sorted(set(unsupported)))}. "
            f"This table is indexed on {key_column!r}, so only filters on that "
            f"column can be answered without a full scan. A filter that silently "
            f"matched nothing would read as 'no rows match'.",
            detail={"columns": sorted(set(unsupported)), "key_column": key_column},
        )
    predicate = None
    if clauses:
        predicate = lambda value: all(c.matches({key_column: value}) for c in clauses)

    sort_column: Optional[str] = None
    descending = False
    if page_request.sort:
        descending = page_request.sort.startswith("-")
        sort_column = page_request.sort[1:] if descending else page_request.sort

    # Filter and sort from the cached key values, then materialise only the
    # page. `key_of` reads a column out of a row by seeking to that row's range,
    # so a sort on a column other than the key costs one read per candidate
    # rather than a full materialisation.
    def key_of(value: str, start: int) -> Any:
        if sort_column == key_column:
            return value
        # A sort on another column needs that cell. One `seek` + `readline` for
        # one line per candidate row, not a materialisation of the table.
        rows = reader.read_range(start, start, header, max_lines=1)
        return rows[0].get(sort_column) if rows else None

    items, total, total_unfiltered, page_bytes = reader.page(
        header=header,
        spans=spans,
        predicate=predicate,
        # `key_of=None` keeps the file's own row order, which is what a request
        # with no `sort` asked for.
        key_of=key_of if sort_column is not None else None,
        descending=descending,
        offset=page_request.offset,
        limit=page_request.limit,
    )
    size = path.stat().st_size if path.exists() else 0
    read_bytes = scan_bytes + page_bytes
    meta = page_meta(
        page_request,
        total=total,
        total_unfiltered=total_unfiltered,
        basis=basis(
            n=total,
            artefact=str(path.name),
            path=path,
            rows_total=total_unfiltered,
            filtered=bool(page_request.has_filter or page_request.q),
            min_samples=state.min_samples,
        ),
    )
    meta["read"] = {
        "bytes_read": read_bytes,
        "file_size": size,
        # `fully_loaded` compares the bytes this request read against the file
        # size. A page that reads a small fraction is the mechanism working; a
        # page that reads the whole file is a scan, which happens once per file
        # and is then cached.
        "fully_loaded": bool(size and read_bytes >= size),
        "mechanism": (
            "byte-range page over a cached byte-offset index (UI-D5)"
            if not scan_bytes
            else "streaming scan (first read of this table), then a cached "
            "byte-offset index (UI-D5)"
        ),
        "key_column": key_column,
        "sort_column": sort_column,
        "unsupported_filter_columns": unsupported,
        "index_cached": not scan_bytes,
    }
    return items, meta


def _stage_rows_impl(
    state: AppState,
    stage: str,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    path = state.source.stage_table_path(stage)
    peek = peek_table(path, name=stage)
    payload: Dict[str, Any] = {
        "stage": stage,
        # `peek_table`'s verdict: presence, header and row count in one streaming
        # pass that keeps no rows. A huge table is never materialised just to
        # learn its byte offsets (UI-D5).
        "present": peek.present,
        "reason": peek.reason,
        "path": str(path),
        "header": list(peek.header) or list(_actual_header(path)[0]),
    }
    if not peek.present:
        payload["rows"] = []
        payload["items"] = []
        payload["meta"] = page_meta(
            parse_page(offset=offset, limit=limit, available=[]),
            total=0,
            total_unfiltered=0,
            basis=basis(
                n=0,
                artefact=stage,
                path=path,
                rows_total=0,
                min_samples=state.min_samples,
            ),
        )
        payload["power"] = _power(0)
        return payload

    parsed = parse_page(
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        available=list(peek.header),
        searchable=list(peek.header),
    )
    items, meta = _rows_via_index(
        state,
        path=path,
        key_column="sample_id" if "sample_id" in peek.header else peek.header[0],
        page_request=parsed,
        searchable=list(peek.header),
    )
    if parsed.columns:
        keep = set(parsed.columns)
        items = [{k: v for k, v in row.items() if k in keep} for row in items]
    payload.update(
        {
            "items": items,
            "rows": items,
            "meta": meta,
            "power": _power(meta["basis"]["n"]),
            "n_rows_total": peek.n_rows,
        }
    )
    return payload


def _table_rows_impl(
    state: AppState,
    key: str,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    path = state.source.table_path_for(key)
    peek = peek_table(path, name=key)
    payload: Dict[str, Any] = {
        "key": key,
        "present": peek.present,
        "reason": peek.reason,
        "path": str(path),
        "header": list(peek.header) or list(_actual_header(path)[0]),
    }
    if not peek.present:
        payload["rows"] = []
        payload["items"] = []
        payload["meta"] = page_meta(
            parse_page(offset=offset, limit=limit, available=[]),
            total=0,
            total_unfiltered=0,
            basis=basis(n=0, artefact=key, path=path, rows_total=0, min_samples=state.min_samples),
        )
        payload["power"] = _power(0)
        return payload
    parsed = parse_page(
        offset=offset, limit=limit, sort=sort, q=q, columns=columns,
        include_total=include_total, filter_params=filter_params,
        available=list(peek.header), searchable=list(peek.header),
    )
    key_column = "sample_id" if "sample_id" in peek.header else peek.header[0]
    items, meta = _rows_via_index(
        state, path=path, key_column=key_column, page_request=parsed,
        searchable=list(peek.header),
    )
    if parsed.columns:
        keep = set(parsed.columns)
        items = [{k: v for k, v in row.items() if k in keep} for row in items]
    payload.update(
        {"items": items, "rows": items, "meta": meta, "power": _power(meta["basis"]["n"]),
         "n_rows_total": peek.n_rows}
    )
    return payload


def _events_for(state: AppState, stage: str) -> List[Dict[str, Any]]:
    path = state.log_path()
    if path is None or not Path(path).exists():
        return []
    replay = event_module.read_from(path, 0)
    return [
        e.as_dict()
        for e in replay.events
        if e.mapped_stage == stage or e.stage == stage
    ]


def _event_snapshot(state: AppState) -> Dict[str, Any]:
    """The full replayed state, which both transports carry.

    Emitted as the SSE `snapshot` frame before any `delta`, and returned by
    `/api/events/state`, so a client that connects mid-run is never blank.
    """
    signals = _event_signals(state)
    stages = [_stage_payload(state, stage, signals=signals) for stage in STAGE_ORDER]
    path = state.log_path()
    replay = event_module.read_from(path, 0)
    stage_state, totals, unmapped = event_module.aggregate(
        replay.events, stage_order=STAGE_ORDER
    )
    return {
        "transport": "sse",
        "cursor": replay.cursor,
        "n_events": len(replay.events),
        "n_lines_skipped": replay.n_skipped,
        "skipped_warnings": replay.warnings[:200],
        "unmapped_stages": unmapped[:200],
        "totals": totals,
        "stage_event_state": {
            name: {
                "started": entry.started,
                "done": entry.done,
                "failed": entry.failed,
                "running_samples": len(entry.running_samples),
                "badge_hint": entry.badge_hint,
                "last_event": entry.last_event.as_dict() if entry.last_event else None,
            }
            for name, entry in stage_state.items()
        },
        "stages": stages,
    }


def _isolate_reads(state: AppState) -> Dict[str, Any]:
    """Every per-isolate table, read once."""
    source = state.source
    return {
        "validation": source.stage_table("validation"),
        "mlst": source.stage_table("mlst"),
        "amr": source.stage_table("amr"),
        "virulence": source.stage_table("virulence"),
        "phenotype": source.stage_table("phenotype"),
        "mechanisms": source.internal_table("mechanisms"),
        "regulators": source.internal_table("regulators"),
        "structural_variants": source.internal_table("structural_variants"),
        "master_table": source.internal_table("master_table"),
    }


def _variants_provenance(state: AppState) -> Tuple[Optional[Mapping[str, Any]], str]:
    payload, reason = provenance_module.read_json(state.source.optional("variants_provenance"))
    if isinstance(payload, Mapping):
        return payload, ""
    return None, reason or f"{NOT_PRODUCED}: variants_provenance.json is not readable"


def _tree_metadata(state: AppState) -> Any:
    """`tree_metadata.tsv` beside the tree, where a run wrote it.

    Absent is not fatal: it blocks only the per-tip lineage and ST labels, and
    the reason says so rather than leaving the tree unexplained.
    """
    from papipeline.stages.report_tables import read_table

    source = state.source
    candidates: List[Path] = []
    phylo = source.optional("phylogeny_dir")
    if phylo is not None:
        root = Path(phylo)
        candidates.extend([root / "tree_metadata.tsv", root.parent / "tree_metadata.tsv"])
    candidates.append(state.stage_dir / "tree_metadata.tsv")
    for candidate in candidates:
        read = read_table("tree_metadata", candidate)
        if read.present:
            return read
    return TableRead(
        name="tree_metadata",
        path=candidates[0] if candidates else Path("<no phylogeny dir>"),
        reason=(
            "no tree_metadata.tsv was found beside the stage-9 tree, in its "
            "parent, or in the stage directory. Its absence does not block the "
            "tree itself; it blocks only the per-tip lineage and ST labels."
        ),
    )


def _tree_impl(state: AppState, source: str) -> Dict[str, Any]:
    """The stage-9 tree as structural JSON, with node ids by postorder index."""
    paths = state.source.newick_paths()
    if not paths:
        return {
            "present": False,
            "reason": (
                f"{NOT_PRODUCED}: no Newick file was found. Probed "
                f"the stage-9 phylogeny directory declared by "
                f"`config.loader.phylogeny_dir`. A tree's absence is a gap in "
                f"the run's durable record, not an empty tree."
            ),
            "probes": [
                {"n": 1, "path": str(p), "ok": False, "reason": "no file at this path"}
                for p in ([state.source.optional("phylogeny_dir")] if state.source.optional("phylogeny_dir") else [])
            ],
        }
    chosen = _choose_tree(paths, source)
    if chosen is None:
        return {
            "present": False,
            "reason": (
                f"{NOT_PRODUCED}: no {source!r} tree was found among "
                f"{len(paths)} Newick file(s). Candidates: "
                + ", ".join(p.name for p in paths)
            ),
            "candidates": [{"n": i + 1, "path": str(p), "ok": False,
                            "reason": f"not a {source} tree"} for i, p in enumerate(paths)],
        }
    try:
        text = Path(chosen).read_text(encoding="utf-8")
        root = newick_module.parse(text)
    except (OSError, newick_module.NewickError) as exc:
        return {
            "present": False,
            "reason": f"{chosen} could not be parsed as Newick: {exc}",
            "path": str(chosen),
        }
    metadata = _tree_metadata(state)
    phenotype_read = state.source.stage_table("phenotype")
    phenotype_by_sample: Dict[str, str] = {}
    if phenotype_read.present:
        for row in phenotype_read.rows:
            sample = row.get("sample_id")
            value = row.get("phenotype")
            if sample and value is not None:
                phenotype_by_sample[str(sample)] = str(value)
    tip_metadata: Dict[str, Dict[str, Any]] = {}
    for row in (metadata.rows if metadata.present else ()):
        sample = row.get("sample_id")
        if sample:
            tip_metadata[str(sample)] = {
                "lineage_label": row.get("lineage_label"),
                "st": row.get("st"),
                "source": row.get("source"),
                # The SIR call from `11_phenotype.tsv`, joined by `sample_id`.
                # Absent stays absent: a tip with no phenotype row reads `None`,
                # never a fabricated category (UI-D2).
                "phenotype": phenotype_by_sample.get(str(sample)),
            }
    payload = newick_module.to_json(root, tip_metadata=tip_metadata)
    payload.update(
        {
            "present": True,
            "reason": "",
            "source": "stage9" if "gubbins" not in Path(chosen).name else "gubbins",
            "path": str(Path(chosen).resolve()),
            "bytes": len(text.encode("utf-8")),
            "tip_check": newick_module.tip_check(
                newick_module.tip_labels(root),
                list(tip_metadata) or newick_module.tip_labels(root),
            ),
            "metadata_reason": metadata.reason if not metadata.present else "",
        }
    )
    return payload


def _choose_tree(paths: Sequence[Path], source: str) -> Optional[Path]:
    """Pick the requested tree from the candidates."""
    if source == "all":
        return paths[0] if paths else None
    for path in paths:
        name = path.name.lower()
        if source == "gubbins" and ("gubbins" in name or "final_tree" in name):
            return path
        if source == "stage9" and "gubbins" not in name and "final_tree" not in name:
            return path
    return None


def _tree_tips_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    tree = _tree_impl(state, "stage9")
    if not tree.get("present"):
        return {
            "items": [],
            "present": False,
            "reason": tree.get("reason", ""),
            "meta": page_meta(
                parse_page(offset=offset, limit=limit, available=[]),
                total=0,
                total_unfiltered=0,
                basis=basis(n=0, artefact="tree", path=tree.get("path"), rows_total=0,
                            min_samples=state.min_samples),
            ),
        }
    tips = [
        {
            "node_id": node["node_id"],
            "sample_id": node["label"] or "",
            "label": node["label"] or "",
            "lineage_label": (node.get("tip_metadata") or {}).get("lineage_label"),
            "st": (node.get("tip_metadata") or {}).get("st"),
            "phenotype": (node.get("tip_metadata") or {}).get("phenotype"),
            "source": (node.get("tip_metadata") or {}).get("source"),
        }
        for node in tree["nodes"]
        if node["is_tip"]
    ]
    available = ["node_id", "sample_id", "label", "lineage_label", "st", "phenotype", "source"]
    paged = _page(
        state,
        tips,
        available=available,
        searchable=("sample_id", "label", "lineage_label", "st"),
        artefact="tree",
        path=Path(tree["path"]) if tree.get("path") else None,
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        key_fn=lambda row: row["node_id"],
    )
    paged["present"] = True
    paged["reason"] = ""
    paged["identity"] = newick_module.IDENTITY
    return paged


def _similarity_impl(state: AppState, max_tips: int) -> Dict[str, Any]:
    """The patristic matrix, downsampled server-side above `max_tips`.

    The downsample is **reported**: `n_total`, `n_shown`, `downsampled` and a
    `method` string. A silent truncation is not offered (UI-D2), because a
    matrix that quietly holds 300 of 900 rows reads as the whole cohort.
    """
    read = state.source.stage_table("similarity")
    units_path = state.source.optional("similarity_units")
    units, units_reason = provenance_module.read_json(units_path)
    if not read.present:
        return {
            "present": False,
            "reason": read.reason,
            "units": units if isinstance(units, Mapping) else None,
            "units_reason": "" if isinstance(units, Mapping) else units_reason,
            "n_total": 0,
            "n_shown": 0,
            "downsampled": False,
            "method": None,
            "order": [],
            "matrix": [],
            "basis": basis(n=0, artefact="similarity", path=read.path, rows_total=0,
                           min_samples=state.min_samples),
        }
    # `similarity.tsv` is one row per sample holding a whole `distances` cell
    # (`contracts.py:118`), keyed on the cohort rather than on columns, so the
    # cell's position i is the i-th sample **in the file's own row order**.
    # That order is what the vectors are indexed against.
    order: List[str] = []
    vectors: Dict[str, List[float]] = {}
    for row in read.rows:
        sample = row.get("sample_id")
        distances = row.get("distances")
        if not sample or distances is None:
            continue
        parsed = _parse_vector(str(distances))
        if parsed is None:
            continue
        order.append(str(sample))
        vectors[str(sample)] = parsed

    # A SQUARE matrix: every kept row's vector must cover every kept column.
    # A short vector is a truncated read, and padding it with zeros would
    # render distances nobody computed, so those samples are dropped and the
    # drop is reported.
    complete = [s for s in order if len(vectors[s]) >= len(order)]
    incomplete = [s for s in order if s not in set(complete)]
    order = complete
    n_total = len(order)

    shown = order
    downsampled = False
    method: Optional[str] = None
    if n_total > max_tips:
        shown = newick_module.downsample_by_single_linkage(order, vectors, max_tips)
        downsampled = True
        method = (
            f"single-linkage on the stage-9 tree, k={max_tips}; farthest-point "
            f"sampling on the stage-10 patristic matrix"
        )
    positions = {name: i for i, name in enumerate(order)}
    matrix = [
        [vectors[sample][positions[other]] for other in shown] for sample in shown
    ]
    return {
        "present": True,
        "reason": "",
        "n_total": n_total,
        "n_shown": len(shown),
        "downsampled": downsampled,
        "method": method,
        "order": shown,
        "matrix": matrix,
        # Passed through unmodified. `is_a_snp_count` is false by construction,
        # and the UI must say substitutions per site and never SNPs (rule 11).
        "units": units if isinstance(units, Mapping) else None,
        "units_reason": "" if isinstance(units, Mapping) else units_reason,
        "n_rows_in_file": read.n_rows,
        "n_incomplete_vectors": len(incomplete),
        "incomplete_samples": incomplete[:200],
        "basis": basis(n=n_total, artefact="similarity", path=read.path,
                       rows_total=read.n_rows, min_samples=state.min_samples),
    }


def _parse_vector(raw: str) -> Optional[List[float]]:
    """The `distances` cell: comma- or semicolon-separated numbers."""
    text = raw.strip().strip("[]()")
    if not text:
        return None
    for separator in (";", ",", " "):
        if separator in text:
            parts = [p for p in text.split(separator) if p.strip()]
            break
    else:
        parts = text.split()
    try:
        return [float(p) for p in parts]
    except ValueError:
        return None


def _gwas_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    read = state.source.stage_table("gwas")
    if not read.present:
        return {
            "present": False,
            "reason": read.reason,
            "n_tests": 0,
            "items": [],
            "meta": page_meta(
                parse_page(offset=offset, limit=limit, available=[]),
                total=0,
                total_unfiltered=0,
                basis=basis(n=0, artefact="gwas", path=read.path, rows_total=0,
                            min_samples=state.min_samples),
            ),
        }
    rows: List[Dict[str, Any]] = []
    p_values: List[float] = []
    for row in read.rows:
        model = str(row.get("model") or "")
        entry = dict(row)
        # Shown verbatim. Every reference-engine result is stamped
        # `reference_fisher:NO_KINSHIP_CORRECTION`, and a p-value that omits
        # the stamp is anti-conservative in a structured cohort, so a row whose
        # model does not name its correction is flagged rather than hidden.
        entry["model_states_kinship_correction"] = "NO_KINSHIP_CORRECTION" in model
        # The numeric columns are coerced where they parse. `read_tsv` returns
        # text, and a QQ plot needs numbers; a cell that does not parse stays
        # as its string rather than becoming a zero (UI-D2).
        for column in ("p_value", "adjusted_p_value", "effect_size", "frequency"):
            number = _number_or_none(entry.get(column))
            if number is not None:
                entry[column] = number
        rows.append(entry)
        value = entry.get("p_value")
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            p_values.append(float(value))
    n_tests = len(rows)
    available = sorted({k for row in rows for k in row})
    paged = _page(
        state,
        rows,
        available=available,
        searchable=("feature", "feature_type", "lineage_distribution", "model"),
        artefact="gwas",
        path=read.path,
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        key_fn=lambda row: (row.get("p_value") is None, row.get("p_value") or 0.0),
        basis_n=_cohort_size(state),
    )
    paged.update(
        {
            "present": True,
            "reason": "",
            # The row count, from which the QQ plot derives its expected
            # quantiles server-side rather than in the browser (UI-D7).
            "n_tests": n_tests,
            "p_values": sorted(p_values)[:QQ_MAX],
            "p_values_truncated": len(p_values) > QQ_MAX,
        }
    )
    return paged


def _gwas_top_impl(state: AppState, limit: int) -> Dict[str, Any]:
    read = state.source.stage_table("gwas")
    if not read.present:
        return {
            "present": False,
            "reason": read.reason,
            "items": [],
            "basis": basis(n=0, artefact="gwas", path=read.path, rows_total=0,
                           min_samples=state.min_samples),
        }
    threshold = None
    try:
        from papipeline.config.loader import load_config

        config = load_config()
        threshold = float(config.gwas.lineage_confound_threshold)
    except Exception:
        threshold = 0.9
    items: List[Dict[str, Any]] = []
    for row in read.rows:
        share = _dominant_lineage_share(row.get("lineage_distribution"))
        items.append(
            {
                "feature": row.get("feature"),
                "adjusted_p_value": _number_or_none(row.get("adjusted_p_value")),
                "p_value": _number_or_none(row.get("p_value")),
                "dominant_lineage_share": share,
                # A feature carried only by samples of one lineage cannot be
                # separated from that lineage by any association test. It is
                # RELABELLED, never deleted, and never presented as a
                # resistance determinant.
                "lineage_linked": bool(share is not None and share >= threshold),
            }
        )
    items.sort(key=lambda entry: (entry["adjusted_p_value"] is None, entry["adjusted_p_value"] or 0.0))
    shown = items[: max(1, min(limit, 200))]
    return {
        "present": True,
        "reason": "",
        "items": shown,
        "lineage_confound_threshold": threshold,
        # `n` is the cohort the GWAS was computed on, not the number of hits on
        # this page. `n_returned` and `n_tests` say how many items there are,
        # so a 25-hit top list on a 900-isolate run does not read as n=25.
        "basis": basis(n=_cohort_size(state), artefact="gwas", path=read.path,
                       rows_total=len(items), min_samples=state.min_samples),
        "n_returned": len(shown),
        "n_tests": len(items),
    }


def _dominant_lineage_share(raw: Any) -> Optional[float]:
    """The largest lineage's share, from the `lineage_distribution` cell."""
    if raw is None:
        return None
    text = str(raw)
    total = 0.0
    largest = 0.0
    saw_number = False
    for part in text.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            _label, _, count = part.rpartition(":")
            try:
                value = float(count)
            except ValueError:
                value = 1.0
        else:
            try:
                value = float(part)
            except ValueError:
                value = 1.0
        saw_number = True
        total += value
        largest = max(largest, value)
    if not saw_number or total <= 0:
        return None
    return largest / total


def _number_or_none(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _convergence_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    read = state.source.stage_table("convergence")
    categories = (
        "widespread_background",
        "rare_isolated",
        "lineage_associated",
        "recurrent_convergent",
        "unknown",
    )
    if not read.present:
        return {
            "present": False,
            "reason": read.reason,
            "category_order": list(categories),
            "items": [],
            "meta": page_meta(
                parse_page(offset=offset, limit=limit, available=[]),
                total=0,
                total_unfiltered=0,
                basis=basis(n=0, artefact="convergence", path=read.path, rows_total=0,
                            min_samples=state.min_samples),
            ),
        }
    rows = [dict(r) for r in read.rows]
    for row in rows:
        for column in ("independent_lineages", "branch_count", "statistic_value",
                       "n_a", "n_b", "n_both", "adjusted_p_value"):
            number = _number_or_none(row.get(column))
            if number is not None and row.get(column) is not None:
                row[column] = number
    available = sorted({k for row in rows for k in row})
    paged = _page(
        state,
        rows,
        available=available,
        searchable=("determinant", "convergence_category", "distribution"),
        artefact="convergence",
        path=read.path,
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        basis_n=_cohort_size(state),
    )
    paged.update(
        {
            "present": True,
            "reason": "",
            # The five categories in the order they are evaluated, which is what
            # makes them mutually exclusive. Part of the science; not re-sorted
            # for display.
            "category_order": list(categories),
            "limit_note": (
                "Convergence is a pattern in the data, not a demonstration of a "
                "shared selective cause."
            ),
        }
    )
    return paged


def _cooccurrence_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    read = state.source.stage_table("cooccurrence")
    if not read.present:
        return {
            "present": False,
            "reason": read.reason,
            "items": [],
            "meta": page_meta(
                parse_page(offset=offset, limit=limit, available=[]),
                total=0,
                total_unfiltered=0,
                basis=basis(n=0, artefact="cooccurrence", path=read.path, rows_total=0,
                            min_samples=state.min_samples),
            ),
        }
    rows = [dict(r) for r in read.rows]
    for row in rows:
        for column in ("n_a", "n_b", "n_both", "statistic_value", "adjusted_p_value"):
            number = _number_or_none(row.get(column))
            if number is not None and row.get(column) is not None:
                row[column] = number
    available = sorted({k for row in rows for k in row})
    paged = _page(
        state,
        rows,
        available=available,
        searchable=("feature_a", "feature_b", "feature_type", "statistic"),
        artefact="cooccurrence",
        path=read.path,
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        basis_n=_cohort_size(state),
    )
    paged.update(
        {
            "present": True,
            "reason": "",
            "limit_note": (
                "Co-occurrence is association only and is never presented as an "
                "interaction."
            ),
        }
    )
    return paged


def _pangenome_impl(state: AppState) -> Dict[str, Any]:
    """`pangenome_summary.tsv` as metric/value pairs, plus the gene-table counts.

    The three gene tables are **counted, not paged**: a 900-isolate
    presence/absence table is 200,000 rows and the counts are what the page
    shows.
    """
    read = state.source.stage_table("pangenome")
    summary: Dict[str, Any] = {}
    if read.present:
        for row in read.rows:
            metric = row.get("metric")
            if metric is not None:
                summary[str(metric)] = row.get("value")
    gene_tables: List[Dict[str, Any]] = []
    for key, _filename, declared_by in _snakefile_only():
        table = state.source.snakefile_table(key)
        gene_tables.append(
            {
                "key": key,
                "path": str(table.path),
                "declared_by": declared_by,
                "present": table.present,
                "reason": table.reason,
                "n_rows": table.n_rows if table.present else None,
            }
        )
    return {
        "present": read.present,
        "reason": read.reason,
        "summary": summary,
        "gene_tables": gene_tables,
        "basis": basis(
            # The cohort the pangenome was computed on, not the metric-row count.
            # A 900-isolate pangenome whose summary table has 5 rows is n=900.
            n=_cohort_size(state),
            artefact="pangenome",
            path=read.path,
            rows_total=read.n_rows if read.present else 0,
            min_samples=state.min_samples,
        ),
    }


def _snakefile_only() -> Tuple[Tuple[str, str, str], ...]:
    from papipeline.stages.report_tables import SNAKEFILE_ONLY_TABLES

    return SNAKEFILE_ONLY_TABLES


def _isolates_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    sort: Optional[str],
    q: Optional[str],
    columns: Optional[str],
    include_total: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    """The joined per-isolate table. Instant at 900 rows.

    Every contributing table is read once and folded into a per-`sample_id`
    index, so a page costs one pass per table rather than one pass per row
    (UI-D5).
    """
    reads = _isolate_reads(state)
    variants_payload, variants_reason = _variants_provenance(state)
    metadata = _tree_metadata(state)
    index = isolate_module.build_index(
        mlst=reads["mlst"],
        amr=reads["amr"],
        virulence=reads["virulence"],
        phenotype=reads["phenotype"],
        master=reads["master_table"],
        tree_metadata=metadata,
        variants_provenance=variants_payload,
    )
    oprd_items, oprd_probes, oprd_reason = oprd_module.collect(state.root)
    oprd_by_sample = oprd_module.by_sample(oprd_items)

    samples, membership = isolate_module.cohort(_manifest_samples(state), reads)
    rows = [isolate_module.row_for(sample, index, oprd_by_sample) for sample in samples]
    # Every row carries its own `basis`, so a client reading a single isolate's
    # cells knows which cohort each cell came from without a second request.
    per_row_basis = basis(
        n=len(samples),
        artefact="isolates",
        path=reads["amr"].path,
        rows_total=len(samples),
        min_samples=state.min_samples,
    )
    for row in rows:
        row["basis"] = dict(per_row_basis)

    available = list(isolate_module.COLUMNS) + [
        "has_amr_gene",
        "has_point_mutation",
        "virulence_min",
    ]
    paged = _page(
        state,
        rows,
        available=available,
        searchable=list(isolate_module.SEARCHABLE),
        artefact="isolates",
        path=reads["amr"].path,
        offset=offset,
        limit=limit,
        sort=sort,
        q=q,
        columns=columns,
        include_total=include_total,
        filter_params=filter_params,
        key_fn=lambda row: row["sample_id"],
        # `basis` is the provenance of every cell, not a data column, so a
        # projection must not drop it. `_page` parses the csv `columns` and
        # keeps this key alongside the requested ones.
        projection_extra=("basis",),
        extra={
            "membership_source": membership,
            "variants_reason": variants_reason,
            "oprd_available": bool(oprd_items),
            "oprd_reason": oprd_reason,
            "contributing_sources": isolate_module.contributing_sources(index),
            "absent_behaviour": isolate_module.absent_behaviour(),
        },
    )
    return paged


def _manifest_samples(state: AppState) -> Optional[List[str]]:
    """The sample list the manifest records, when it records one."""
    manifest = state.manifest
    if not isinstance(manifest, Mapping):
        return None
    for key in ("samples", "sample_ids", "n_samples_list", "cohort"):
        value = manifest.get(key)
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            out = [str(v) for v in value]
            if out:
                return out
    outputs = manifest.get("outputs")
    if isinstance(outputs, Mapping):
        value = outputs.get("samples") or outputs.get("sample_ids")
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            return [str(v) for v in value]
    return None


def _cohort_size(state: AppState) -> int:
    """The cohort a statistic was computed on, for `basis.n` and the R16 flag.

    The manifest's sample list when there is one, else the union of the
    per-sample tables the way `isolates.cohort` builds it. This is the
    denominator the D3 power flag must use (DESIGN §7.1): a pangenome computed
    over 900 isolates is not `n=5` because the summary table has 5 metric rows.
    """
    samples = _manifest_samples(state)
    if samples:
        return len(set(str(s) for s in samples))
    reads = _isolate_reads(state)
    seen: set = set()
    for read in reads.values():
        if read is None or not read.present:
            continue
        for row in read.rows:
            sample = row.get("sample_id")
            if sample:
                seen.add(str(sample))
    return len(seen)


def _isolate_detail_impl(state: AppState, sample_id: str) -> Dict[str, Any]:
    reads = _isolate_reads(state)
    metadata = _tree_metadata(state)
    everything = dict(reads)
    everything["tree_metadata"] = metadata
    variants_payload, variants_reason = _variants_provenance(state)
    samples, membership = isolate_module.cohort(_manifest_samples(state), reads)
    if sample_id not in samples:
        return _error(
            f"{sample_id!r} is not in the manifest and appears in no per-sample "
            f"table. The cohort is {len(samples)} sample(s) drawn from the "
            f"{membership!r}.",
            404,
        )
    oprd_items, oprd_probes, oprd_reason = oprd_module.collect(state.root)
    index = isolate_module.build_index(
        mlst=reads["mlst"],
        amr=reads["amr"],
        virulence=reads["virulence"],
        phenotype=reads["phenotype"],
        master=reads["master_table"],
        tree_metadata=metadata,
        variants_provenance=variants_payload,
    )
    return {
        "sample_id": sample_id,
        "membership_source": membership,
        "row": isolate_module.row_for(sample_id, index, oprd_module.by_sample(oprd_items)),
        "tables": isolate_module.detail_tables(sample_id, reads=everything),
        "variants_provenance": _variants_for(variants_payload, sample_id),
        "variants_reason": variants_reason,
        "oprd": {
            "available": bool(oprd_items),
            "reason": oprd_reason,
            "probes": oprd_probes,
            "items": [i for i in oprd_items if i["sample_id"] == sample_id],
        },
        "basis": basis(
            n=1,
            artefact="isolates",
            path=reads["amr"].path,
            rows_total=len(samples),
            min_samples=state.min_samples,
        ),
        "power": _power(1),
    }


def _variants_for(payload: Optional[Mapping[str, Any]], sample_id: str) -> Optional[Mapping[str, Any]]:
    if not isinstance(payload, Mapping):
        return None
    per_isolate = payload.get("per_isolate")
    if not isinstance(per_isolate, Sequence):
        return None
    for entry in per_isolate:
        if isinstance(entry, Mapping) and str(entry.get("sample_id")) == sample_id:
            return dict(entry)
    return None


def _oprd_impl(
    state: AppState,
    request: Request,
    offset: int,
    limit: int,
    filter_params: Mapping[str, List[str]],
) -> Dict[str, Any]:
    items, probes, reason = oprd_module.collect(state.root)
    available = [
        "sample_id", "verdict", "lesion_type", "position",
        "truncation_aa", "identity_pct", "coverage_pct", "display_state",
    ]
    parsed = parse_page(
        offset=offset, limit=limit, filter_params=filter_params, available=available
    )
    selected, total, total_unfiltered = apply_page(items, parsed)
    payload: Dict[str, Any] = {
        "available": bool(items),
        "reason": reason,
        "probes": [{"n": index, **record} for index, record in enumerate(probes, start=1)],
        "items": selected,
        "verdict_vocabulary": list(oprd_module.VERDICTS),
        "display_states": list(oprd_module.DISPLAY_STATES),
        "meta": page_meta(
            parsed,
            total=total,
            total_unfiltered=total_unfiltered,
            basis=basis(
                n=total,
                artefact="oprd",
                path=str(state.root),
                rows_total=total_unfiltered,
                min_samples=state.min_samples,
            ),
        ),
    }
    if not items:
        # The two honest states, and the refusal to derive one from the master
        # table's `chromosomal_mutation`. That derivation reported `absent`
        # for 8 gene-carrying isolates on the ten-isolate cohort, and
        # surfacing it would manufacture the study's central negative (§4.3).
        payload["not_produced_note"] = (
            "There is no contracted file for these verdicts. The master "
            "table's `chromosomal_mutation` column and "
            "`viz.oprd_status_per_sample` are NOT used as a fallback: the latter "
            "reported `absent` for 8 isolates that all carry the gene, because "
            "34 of the 36 gene=oprD CDS features are OprD/OprP/OprQ paralogs."
        )
    return payload


def _provenance_impl(state: AppState) -> Dict[str, Any]:
    reuse = _reuse_facts(state)
    return provenance_module.provenance(
        state.manifest,
        reuse,
        repo_root=state.repo_root,
        results_root=state.root,
        machine=state.machine,
        writer=state.manifest_writer,
    )


def _bakta_impl(state: AppState) -> Dict[str, Any]:
    reuse = _reuse_facts(state)
    count = provenance_module.bakta_executions(state.manifest, reuse)
    return {
        "available": bool(reuse.get("available")),
        "reason": reuse.get("reason", ""),
        # Null, not 0, when the count is not recorded. `0` would be a false
        # claim of proof.
        "bakta_executions": count,
        "bakta_executions_label": provenance_module.bakta_executions_label(count),
        "n_samples": state.antibiotic and _n_samples(state, state.manifest),
        "reused_from_run": reuse.get("reused_from"),
        "verified": (
            bool(reuse.get("n_reused")) if reuse.get("available") else None
        ),
        "n_reused": reuse.get("n_reused"),
        "n_not_reused": reuse.get("n_not_reused"),
        "rows": reuse.get("rows", []),
    }


def _reuse_facts(state: AppState) -> Dict[str, Any]:
    path = state.source.optional("annotation_reuse")
    if path is None:
        return {
            "available": False,
            "reason": (
                f"{NOT_PRODUCED}: no reuse_provenance.tsv was "
                f"found. `annotation.write_reuse_provenance` writes it beside "
                f"the Bakta output, so its absence means the record is not on "
                f"disk - not that nothing was reused."
            ),
            "rows": [],
            "n_rows": None,
            "n_reused": None,
            "n_not_reused": None,
            "bakta_versions": [],
        }
    from papipeline.stages.report_tables import read_table

    read = read_table("annotation_reuse", Path(path))
    return provenance_module.bakta_from_reuse(read)


def _logs_impl(state: AppState, name: str, offset: int, limit: int, tail_lines: int) -> Dict[str, Any]:
    """Tripwire and watcher log tails.

    No log path is contracted anywhere (A7), so an absent log is `not produced`
    with its probes named — never an empty log that reads as "nothing was
    reported".
    """
    wanted = [n for n in provenance_module.LOG_PROBE_NAMES if name in ("any", n)]
    items: List[Dict[str, Any]] = []
    for key in wanted:
        path = state.source.optional(f"{key}_log")
        tail = provenance_module.tail_text(Path(path) if path else None, lines=tail_lines)
        probe_list = _log_probes(state, key)
        items.append(
            {
                "name": key,
                "present": bool(tail["present"]),
                "reason": tail["reason"],
                "path": tail["path"],
                "n_lines": tail["n_lines"],
                "lines": tail["lines"],
                "probes": probe_list,
            }
        )
    parsed = parse_page(offset=offset, limit=limit, available=["name", "present", "n_lines"])
    selected, total, _unfiltered = apply_page(items, parsed)
    return {
        "items": selected,
        "meta": page_meta(
            parsed,
            total=total,
            total_unfiltered=len(items),
            basis=basis(n=total, artefact="provenance_logs", path=str(state.root),
                        rows_total=len(items), min_samples=state.min_samples),
        ),
    }


def _log_probes(state: AppState, key: str) -> List[Dict[str, Any]]:
    """Every path probed for one log, from the adapter's own declarations."""
    names = [f"{key}.log", f"{key}.txt", f"{key}.jsonl"]
    probes: List[Dict[str, Any]] = []
    for name in names:
        path = state.stage_dir / name
        probes.append(
            {
                "n": len(probes) + 1,
                "path": str(path),
                "ok": path.exists(),
                "reason": "matched" if path.exists() else "no file at this path",
            }
        )
    return probes


def _digests_impl(state: AppState) -> Dict[str, Any]:
    """SHA-256 of each stage artefact, streamed."""
    paths: Dict[str, Path] = {}
    source = state.source
    for stage in STAGE_ORDER:
        read = source.stage_table(stage)
        if Path(read.path).exists():
            paths[stage] = Path(read.path)
    for key in ("structural_variants", "regulators", "mechanisms", "master_table"):
        read = source.internal_table(key)
        if Path(read.path).exists():
            paths[key] = Path(read.path)
    return {"items": provenance_module.digests(paths)}


def _reports_impl(state: AppState) -> Dict[str, Any]:
    """Report files, discovered by extension rather than by name.

    Only `.md`, `.html` and the figure extensions are returned, and the path is
    relative to the opened root so it can be handed straight to
    `/api/reports/content`.
    """
    root = Path(state.root)
    items: List[Dict[str, Any]] = []
    for path in state.source.report_files():
        try:
            stat = path.stat()
        except OSError:
            continue
        suffix = path.suffix.lower()
        if suffix in (".md",):
            kind = "markdown"
        elif suffix in (".html", ".htm"):
            kind = "html"
        elif suffix in (".svg", ".png", ".pdf"):
            kind = "figure"
        else:
            kind = "other"
        try:
            relative = path.resolve().relative_to(root)
        except ValueError:
            # Outside the opened root after following a symlink. Not listed,
            # because the content endpoint would refuse it and a link that
            # always 400s is worse than an absent one.
            continue
        items.append(
            {
                "path": str(relative),
                "name": path.name,
                "kind": kind,
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    return {"items": items}


__all__ = ["AppState", "build_state", "create_app"]