# Design 4 — API design

Companion to [1 architecture](01-architecture.md).

## 1. One verb layer, two front doors

```
                       ┌──────────────────────────┐
   papipeline run ────▶│                          │
   papipeline status ─▶│  papipeline_app.backend  │──▶ store / engine / resources / tools
   papipeline tools ──▶│      (the verbs)         │
   GUI  "Run" button ─▶│                          │
                       └──────────────────────────┘
```

`backend.py` is the whole application surface. It is an ordinary Python
module of plain functions — no framework types, no request objects, no
`Response`. It returns plain dicts and dataclasses.

That constraint is what makes the CLI/GUI parity requirement mechanically
testable: `test_backend_parity.py` walks the verb registry, calls each verb
against a fixture project, and asserts both front doors expose the same set
with the same argument names. If the GUI ever grows a capability the CLI
lacks, that test fails.

```python
# papipeline_app/backend.py — the shape, not the implementation
def project_create(
    name: str,
    genome_dir: Path,
    metadata_tsv: Path,
    phenotype_tsv: Path,
    *,
    amr_tsv: Path | None = None,
    antibiotics: Sequence[str] = ("imipenem",),
    genome_source: str = "copy",
    profile: str | None = None,
) -> dict: ...

def run_start(project: str, *, stages=None, skip=None, profile=None,
              dry_run=False, resume=False) -> dict: ...
def run_status(run_id: str) -> dict: ...
def run_resume(run_id: str) -> dict: ...
```

## 2. Verb registry

One declaration drives the CLI parser, the HTTP routes, the OpenAPI schema
and the parity test. Adding a verb in one place adds it everywhere.

```python
# papipeline_app/backend.py
@dataclass(frozen=True)
class Verb:
    name: str
    fn: Callable[..., dict]
    args: tuple[ArgSpec, ...]
    summary: str
    mutates: bool          # False for read-only verbs → safe for GET

VERBS = (
    Verb("project.create", project_create, (...), "Create a project", True),
    Verb("project.list",   project_list,   (),       "List projects", False),
    ...
)
```

`mutates=False` is what lets the HTTP layer map read verbs to `GET` and
write verbs to `POST` without a hand-maintained table that can drift.

## 3. Verb catalogue

| Verb | Args | Returns |
|---|---|---|
| `project.create` | `name, genome_dir, metadata_tsv, phenotype_tsv, amr_tsv?, antibiotics?, genome_source?, profile?` | project record + validation report |
| `project.list` | — | list of project summaries |
| `project.get` | `project` | project + samples + counts |
| `project.delete` | `project, confirm` | deleted paths |
| `project.validate` | `project` | per-sample validation, no writes |
| `input.scan` | `genome_dir` | discovered files, format check, duplicates |
| `input.preview` | `project, n=20` | sample table for the UI's selection grid |
| `config.effective` | `project` | resolved config + which fields differ from repo default |
| `workflow.stages` | — | 16 stages, order, prerequisites, tool requirements |
| `run.start` | `project, stages?, skip?, profile?, dry_run?, resume?` | run record + DAG preview |
| `run.status` | `run_id` | run + per-stage status + elapsed + current job |
| `run.list` | `project?` | recent runs |
| `run.resume` | `run_id` | what will be skipped, what will re-run |
| `run.cancel` | `run_id` | cancelled |
| `run.logs` | `run_id, stage?, tail?` | log text + path |
| `benchmark` | `project, n=5, profile?` | measured s/genome + recommended workers×threads |
| `tools.list` | — | registry: name, version, path, available, source |
| `tools.check` | `name?` | preflight GO/NO-GO, per-tool detail |
| `tools.databases` | — | database registry with versions and pin status |
| `tools.provision` | `database, type` | explicit, user-initiated download |
| `results.list` | `run_id, type?` | registry rows |
| `results.get` | `run_id, type` | typed rows + column schema + summary |
| `results.export` | `run_id, type?, format` | path to written file |
| `report.generate` | `run_id, kind, format?` | written report paths |
| `report.list` | `run_id` | available report kinds |

`kind` for reports: `full_analysis | amr | gwas | antibiotic_comparison`.

## 4. HTTP mapping

Framework choice is open decision **D4**; the mapping is framework-neutral.

```
GET    /api/verbs                          the registry itself
GET    /api/projects                       project.list
POST   /api/projects                       project.create
GET    /api/projects/{name}                project.get
DELETE /api/projects/{name}                project.delete
POST   /api/projects/{name}/validate       project.validate
GET    /api/projects/{name}/samples        input.preview
GET    /api/projects/{name}/config         config.effective
GET    /api/stages                         workflow.stages
POST   /api/runs                           run.start
GET    /api/runs                           run.list
GET    /api/runs/{run_id}                  run.status
POST   /api/runs/{run_id}/resume           run.resume
POST   /api/runs/{run_id}/cancel           run.cancel
GET    /api/runs/{run_id}/logs             run.logs
POST   /api/projects/{name}/benchmark      benchmark
GET    /api/tools                          tools.list
GET    /api/tools/check                    tools.check
GET    /api/databases                      tools.databases
POST   /api/databases                      tools.provision
GET    /api/runs/{run_id}/results          results.list
GET    /api/runs/{run_id}/results/{type}   results.get
GET    /api/runs/{run_id}/reports          report.list
POST   /api/runs/{run_id}/reports          report.generate
GET    /healthz
```

A run is long-lived, so the UI does not hold a request open waiting for it
to finish. Polling is the default and the interval is the client's choice.

> **Superseded in part — 2026-09-27.** The paragraph below originally
> concluded that "the server exposes no SSE or websocket in the first
> version, because log text is served from disk and there is nothing
> push-shaped to push". That reasoning was correct about *log text* and
> wrong about the run. Once `papipeline.execution` existed, a running
> pipeline *does* produce push-shaped state transitions, and the
> observatory needs them: a task moving to `RUNNING`, `RETRYING`,
> `SUCCEEDED`, `FAILED`, `INVALID` or `INCOMPLETE` is exactly the event
> that a node in the network graph animates. Polling at 1 Hz makes the
> network visibly lag behind the work it is depicting.
>
> The observatory therefore adds `GET /api/stream` as server-sent events,
> alongside the polling routes. SSE rather than a websocket because the
> traffic is strictly server-to-client, which needs no socket library and
> makes a reconnect a plain HTTP retry — the right property for a
> long pipeline on a laptop whose lid may close.
>
> Polling is retained: `GET /api/snapshot` and `GET /api/events` are what a
> reconnecting client uses to re-establish the truth, and the stream
> re-sends a full snapshot on every heartbeat tick so a dropped event costs
> a frame of freshness and nothing else. The execution store, not the
> stream, remains authoritative.
>
> The original position is kept below because the reasoning about log text
> still holds: log content is still served from disk on demand and is not
> streamed.

A run is long-lived, so the UI polls rather than holding a socket open.
Polling interval is the client's choice; the server exposes no SSE or
websocket in the first version, because log text is served from disk and
there is nothing push-shaped to push.

## 5. Errors

```python
class AppError(Exception):
    code = "app_error"
    http_status = 500
    def __init__(self, message, **context): ...

class ValidationError(AppError):   code="validation_error";  http_status=422
class NotFoundError(AppError):     code="not_found";        http_status=404
class ConflictError(AppError):     code="conflict";         http_status=409
class ToolUnavailable(AppError):   code="tool_unavailable"; http_status=503
class ResourceError(AppError):     code="insufficient_resource"; http_status=409
class PipelineError(AppError):     code="pipeline_error";   http_status=500
```

`AppError` deliberately mirrors `papipeline.errors.PipelineError(message,
**context)`: context is keyword-only and rendered into the message. One
shape for errors on both sides of the seam means the CLI prints
`str(exc)` and gets something a human can act on, with no formatting
layer.

Scientific errors keep their own vocabulary. A `DataContractError` from
stage 11 is **not** translated into a `ValidationError` — it means the
engine rejected data the app layer thought was fine, which is a different
and more serious thing. The API returns it as `pipeline_error` with the
original type in the payload.

## 6. Input validation — before anything runs

`project.create` refuses to create a project directory unless all of these
pass. Each failure is reported per-sample, never as a single boolean.

| Check | Rule | Source of truth |
|---|---|---|
| Genome files exist and are non-empty | FASTA, ≥1 record, decodes | `papipeline.io.fasta.file_integrity` |
| Genome filenames unique | case-sensitive collision check | new |
| `sample_id` valid | existing rules — no whitespace, no forbidden chars, 2–128, alnum + `-_.` | `papipeline.manifest.validate_sample_id` |
| `sample_id` unique | within the project | `SampleManifest.__post_init__` |
| Metadata joins to genomes | every `assembly_path` resolves | new |
| Every genome has a metadata row | no silent drops | new |
| Phenotype has required columns | `sample_id`, `antibiotic`, `phenotype` | `papipeline.stages.phenotype.REQUIRED_COLUMNS` |
| Phenotype values allowed | `R, I, S, SDD, ND` from config | `config.phenotype.allowed_values` |
| Phenotype samples ⊆ project samples | no orphans | new |
| Antibiotic configured | present in `config.yaml` **and** `antibiotics.tsv` | `PipelineConfig.require_antibiotic` |
| **Selection precedes metadata** | see below | new |

The last one matters scientifically. The engine forbids letting phenotype
influence which samples enter a cohort. So `project.create` takes the
genome directory and metadata, derives the sample set **from the files
alone**, and only then attaches phenotype. The ordering is asserted in
code and covered by a test, mirroring the existing
`tests/unit/test_pilot.py::selection_precedes_metadata`.

Quantitative phenotype is honoured if present and never invented: if the
TSV carries `MIC`/`MIC_unit` they are read; if it does not, the field
stays `None`. `tests/unit/test_phenotype.py::TestNoFabricatedQuantitation`
already guards the engine side and the app layer must not reintroduce the
violation on the way in.

## 7. Concurrency and safety

- One active run per project. A second `run.start` is a `ConflictError`
  (409), not a queue — silently queueing a second run against a mutating
  intermediate directory would be a data race.
- `run.start` is idempotent on `(project, stage_selection, config
  fingerprint)`: re-issuing an identical request returns the existing run.
- All writes that create or delete a project directory require an explicit
  `confirm=True`. The CLI prompts; the API requires the field.
- `tools.provision` is the **only** verb permitted to touch the network, and
  only on explicit user request. This preserves scientific rule 8 and
  `annotation.allow_database_update: false`: nothing downloads implicitly
  during a run.
