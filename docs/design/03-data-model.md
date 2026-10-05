# Design 3 — Data model (SQLite)

Companion to [1 architecture](01-architecture.md). Eight tables, stdlib
`sqlite3`, no ORM. DDL lives in `papipeline_app/store/schema.sql` and is
versioned by `PRAGMA user_version`.

## Design rules

1. **The database is an index, not the source of truth.** Every scientific
   artefact stays a TSV on disk, in the layout the engine already defines.
   SQLite records *what exists, what state it is in, and what produced
   it*. Deleting the database loses the UI, never the science. This is
   what makes design 8's migration cheap and the whole thing recoverable.
2. **Status is an explicit enum, never inferred** from the presence of a
   file. A `COMPLETE` row with a missing file is a bug worth seeing.
3. **Foreign keys are declared and enforced** (`PRAGMA foreign_keys=ON`).
4. **Timestamps are ISO-8601 UTC strings**, matching what `run.py` already
   writes into `run_manifest.json`, so the two never disagree.
5. **Every mutating method takes an explicit transaction.** No implicit
   commits.

---

## 1. Enumerations

```
JobStatus      PENDING | RUNNING | COMPLETE | FAILED | SKIPPED
RunStatus      PENDING | RUNNING | COMPLETE | FAILED | CANCELLED
ProjectStatus  DRAFT | READY | RUNNING | COMPLETE | FAILED | ARCHIVED
ResultType     AMR_GENES | AMR_MECHANISMS | REGULATORS | STRUCTURAL_VARIANTS
               | MLST | VIRULENCE | PAN_GENOME | PHYLOGENY | GWAS
               | CONVERGENCE | CO_OCCURRENCE | INTEGRATED | QC | PHENOTYPE
```

`JobStatus` is exactly the five values the brief specifies. `SKIPPED` is
distinct from `COMPLETE` on purpose: a stage disabled by config is
`SKIPPED`, and conflating the two would let a disabled analysis read as a
successful one.

---

## 2. DDL

```sql
PRAGMA foreign_keys = ON;

CREATE TABLE project (
    id                INTEGER PRIMARY KEY,
    name              TEXT    NOT NULL UNIQUE,
    root              TEXT    NOT NULL,           -- absolute project directory
    organism          TEXT    NOT NULL DEFAULT 'Pseudomonas aeruginosa',
    taxid             INTEGER,
    antibiotics       TEXT    NOT NULL,           -- JSON array
    status            TEXT    NOT NULL DEFAULT 'DRAFT',
    genome_source     TEXT    NOT NULL,           -- 'copy' | 'symlink'
    schema_version    INTEGER NOT NULL DEFAULT 1,
    created_at        TEXT    NOT NULL,
    updated_at        TEXT    NOT NULL
);

CREATE TABLE sample (
    id                INTEGER PRIMARY KEY,
    project_id        INTEGER NOT NULL REFERENCES project(id) ON DELETE CASCADE,
    sample_id         TEXT    NOT NULL,           -- papipeline's identifier rules
    assembly_path     TEXT    NOT NULL,           -- absolute, inside the project
    assembly_sha256   TEXT    NOT NULL,
    assembly_size     INTEGER,
    contig_count      INTEGER,
    n50               INTEGER,
    included          INTEGER NOT NULL DEFAULT 1, -- 0 = QC-excluded, never substituted
    exclusion_reason  TEXT,
    UNIQUE (project_id, sample_id)
);
CREATE INDEX idx_sample_project ON sample(project_id, included);

CREATE TABLE run (
    id                 INTEGER PRIMARY KEY,
    project_id         INTEGER NOT NULL REFERENCES project(id) ON DELETE CASCADE,
    run_id             TEXT    NOT NULL UNIQUE,  -- e.g. run-20260927T112233Z
    status             TEXT    NOT NULL,
    mode               TEXT    NOT NULL,          -- TEST | REAL
    stage_selection    TEXT    NOT NULL,          -- JSON: only/skip
    resource_profile   TEXT    NOT NULL,          -- macbook|workstation|hpc|custom
    resource_config    TEXT    NOT NULL,          -- JSON, the resolved numbers
    pipeline_version   TEXT,
    pipeline_commit    TEXT,
    started_at         TEXT,
    ended_at           TEXT,
    elapsed_seconds    REAL,
    error              TEXT
);
CREATE INDEX idx_run_project ON run(project_id, started_at DESC);

CREATE TABLE job (
    id                INTEGER PRIMARY KEY,
    run_id            INTEGER NOT NULL REFERENCES run(id) ON DELETE CASCADE,
    stage             TEXT    NOT NULL,          -- one of STAGE_ORDER
    seq               INTEGER NOT NULL,          -- position in EXECUTION_ORDER
    status            TEXT    NOT NULL,
    attempt           INTEGER NOT NULL DEFAULT 0,
    workers           INTEGER,                   -- for annotation: genome fan-out
    threads_per_job   INTEGER,
    memory_mb         INTEGER,
    tmpdir            TEXT,
    tool_version      TEXT,
    database_version  TEXT,
    command           TEXT,                      -- exact argv
    started_at        TEXT,
    ended_at          TEXT,
    elapsed_seconds   REAL,
    error             TEXT,
    UNIQUE (run_id, stage)
);
CREATE INDEX idx_job_run_status ON job(run_id, status);

-- Per-genome rows. Only stage 2 (annotation) is per-sample, but the table
-- is generic so MLST/AMR can adopt the same shape later without a migration.
CREATE TABLE job_item (
    id             INTEGER PRIMARY KEY,
    job_id         INTEGER NOT NULL REFERENCES job(id) ON DELETE CASCADE,
    sample_id      TEXT    NOT NULL,
    status         TEXT    NOT NULL,
    attempt        INTEGER NOT NULL DEFAULT 0,
    output_path    TEXT,
    elapsed_seconds REAL,
    tool_version   TEXT,
    error          TEXT,
    UNIQUE (job_id, sample_id)
);
CREATE INDEX idx_job_item_status ON job_item(job_id, status);

-- What the results registry discovered. One row per (run, type).
CREATE TABLE result (
    id            INTEGER PRIMARY KEY,
    run_id        INTEGER NOT NULL REFERENCES run(id) ON DELETE CASCADE,
    type          TEXT    NOT NULL,               -- ResultType
    stage         TEXT    NOT NULL,
    path          TEXT    NOT NULL,               -- TSV/fasta/newick on disk
    format        TEXT    NOT NULL,               -- tsv | fasta | newick | json
    n_rows        INTEGER,
    n_columns     INTEGER,
    schema_fingerprint TEXT,                      -- sha256 of the header line
    summary       TEXT,                           -- JSON, small counts only
    created_at    TEXT    NOT NULL,
    UNIQUE (run_id, type, path)
);
CREATE INDEX idx_result_run_type ON result(run_id, type);

CREATE TABLE provenance (
    id           INTEGER PRIMARY KEY,
    run_id       INTEGER NOT NULL REFERENCES run(id) ON DELETE CASCADE,
    kind         TEXT    NOT NULL,                -- input|software|database|config
                                                   -- |resource|command|environment
    subject      TEXT    NOT NULL,                -- path, tool name, config key
    value        TEXT,
    recorded_at  TEXT    NOT NULL,
    UNIQUE (run_id, kind, subject)
);
CREATE INDEX idx_provenance_run_kind ON provenance(run_id, kind);

CREATE TABLE tool_registry (
    id             INTEGER PRIMARY KEY,
    name           TEXT    NOT NULL,
    stage          TEXT,                          -- NULL = general
    path           TEXT,
    version        TEXT,
    available      INTEGER NOT NULL DEFAULT 0,
    min_version    TEXT,                          -- the floor from tools/preflight.py
    source         TEXT,                          -- env | homebrew | system | missing
    last_checked   TEXT    NOT NULL,
    UNIQUE (name, stage)
);
CREATE TABLE database_registry (
    id             INTEGER PRIMARY KEY,
    name           TEXT    NOT NULL,              -- bakta_db | amrfinderplus | vfdb ...
    reference_id   TEXT,                          -- links to config/references.tsv
    path           TEXT,
    version        TEXT,
    pinned         INTEGER NOT NULL DEFAULT 0,    -- from references.tsv version_status
    schema_version TEXT,
    size_bytes     INTEGER,
    last_checked   TEXT    NOT NULL,
    UNIQUE (name)
);
```

Nine tables. The brief's required fields are all present:
`project`, `sample`, `stage` (as `job.stage`), `job`, `status`,
`start_time` (`started_at`), `end_time` (`ended_at`), `elapsed_time`
(`elapsed_seconds`), `tool_version`, `database_version`, `command`,
`error`.

---

## 3. Status transitions

```
JobStatus:
  PENDING ──▶ RUNNING ──┬──▶ COMPLETE
                        ├──▶ FAILED        (attempt < max_attempts ──▶ PENDING)
                        └──▶ SKIPPED       (config-disabled; terminal, never retried)

RunStatus:
  PENDING ──▶ RUNNING ──┬──▶ COMPLETE     (all jobs COMPLETE|SKIPPED)
                        ├──▶ FAILED        (any job FAILED after retries)
                        └──▶ CANCELLED
```

Two invariants, enforced in code and covered by tests:

- **`SKIPPED` and `COMPLETE` are both terminal and both successful.** A run
  is `COMPLETE` when every job is `COMPLETE` or `SKIPPED`. A run with one
  `FAILED` job is `FAILED`.
- **A `FAILED` job never invalidates its siblings.** This is the brief's
  resume requirement and mirrors the fix already made in
  `scripts/run_annotation_sweep.py::run_one`, which is documented "never
  raises" so one bad genome cannot abandon a queue.

---

## 4. Resumability

Resume reads `job` + `job_item`, then trusts the filesystem only to
*confirm*:

```
for each job in EXECUTION_ORDER:
    if job.status == COMPLETE and every declared output exists:
        skip                      # cheap, and re-verified
    elif job.status in (PENDING, FAILED) and prerequisites COMPLETE:
        run
    else:
        mark SKIPPED(reason='prerequisite not complete')
```

The "every declared output exists" clause matters. A `COMPLETE` row whose
output was deleted must re-run; without the check the UI would show a
green run with nothing behind it.

For annotation, resume granularity is `job_item`, so 60 of 100 genomes
already done means 60 are skipped and only 40 run.

---

## 5. Concurrency

SQLite in WAL mode, one writer at a time. Snakemake may run several stage
jobs concurrently, so:

- Job status updates are short transactions; a stage's scientific work
  happens entirely outside a transaction.
- `busy_timeout` is set to 30 s.
- **SQLite is not the channel for progress.** Stage logs are appended to
  `logs/` by the worker; the UI tails the file. Polling a database for
  log lines would serialise every worker behind the writer lock.
- A `status` snapshot is written to `.papipeline/state.json` on every job
  transition, so a UI read never has to touch SQLite and the project
  directory stays self-describing if the database is lost.

---

## 6. Migrations

`PRAGMA user_version` gates `migrate()`:

```
0 ──▶ 1  initial schema (the DDL above)
1 ──▶ 2  add result.schema_fingerprint
      3  add job.tmpdir
      …
```

Migrations are forward-only, each in its own function, each in one
transaction. There is no down-migration: a project directory plus its
database is reconstructible from provenance, so the database is treated as
disposable cache. That is a deliberate trade — it is what makes "delete
the DB, rebuild from disk" a supported recovery path rather than an
emergency.
