"""Q8 — the actual delivery bundle: layout, no-manifest stage state, oprD.

The REAL 10-isolate run was handed over in a layout that differs from DESIGN
§12's A1–A5 assumption and carries **no `run_manifest.json`** (it is written
only on success; the run failed at stage 7). Before this ticket the dashboard
could not open that bundle at all: `source.detect` raised `NoResultsSource`
with every probe failing.

Two layers of test, deliberately:

* a **hermetic synthetic delivery bundle** built under `tmp_path`, so the
  layout detection and the no-manifest derivation are proven without depending
  on a path outside the repository and without committing a real accession; and
* **real-bundle assertions**, gated on `PA_DASH_REAL_BUNDLE` (default: the
  path this session was given), which check the numbers the bundle's own
  `README.md` / `PER_STAGE_OUTCOMES.md` / `GUARDS.md` state.

No test writes into the bundle: the read-only sweep is asserted on a byte
snapshot. The synthetic sample ids are `SYN_*`; the bundle's accessions never
appear in this file.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict

import pytest

from conftest import REPO_ROOT, snapshot_tree

from dashboard.server import source as source_module

#: Where the REAL delivery bundle is. Overridable so the suite is not tied to
#: one machine; the default is the path this session was given.
REAL_BUNDLE = Path(
    os.environ.get(
        "PA_DASH_REAL_BUNDLE",
        "/Users/raghavkrishnankv/Desktop/PA_AMR_pipeline_real_run_10_isolates",
    )
)


# ---------------------------------------------------------------------------
# A hermetic synthetic delivery bundle
# ---------------------------------------------------------------------------


def _write_table(stage_dir: Path, stage: str, sample: str = "SYN_001") -> None:
    from papipeline.execution.contracts import STAGE_TABLES

    filename, columns = STAGE_TABLES[stage]
    row = ["SYN_001" if column == "sample_id" else "x" for column in columns]
    (stage_dir / filename).write_text(
        "\t".join(columns) + "\n" + "\t".join(row) + "\n", encoding="utf-8"
    )


def build_synthetic_delivery(root: Path) -> Path:
    """A minimal bundle in the ACTUAL delivery layout, no manifest."""
    root = Path(root)
    stages = root / "artifacts" / "stage_tables"
    logs = root / "artifacts" / "logs"
    intermediate = root / "artifacts" / "intermediate"
    guards = root / "guards"
    for directory in (stages, logs, intermediate, guards):
        directory.mkdir(parents=True, exist_ok=True)

    for stage in (
        "validation", "annotation", "mlst", "amr",
        "virulence", "variants", "cohort_variants",
    ):
        _write_table(stages, stage)

    # Folded-step tables, owned by stage 4 / stage 6.
    for key in ("regulators", "structural_variants"):
        from papipeline.execution.contracts import INTERNAL_TABLES

        filename, columns = INTERNAL_TABLES[key]
        row = ["SYN_001" if column == "sample_id" else "x" for column in columns]
        (stages / filename).write_text(
            "\t".join(columns) + "\n" + "\t".join(row) + "\n", encoding="utf-8"
        )

    (logs / "full_run.log").write_text(
        "\n".join(
            [
                "2026-01-01 00:00:00 INFO papipeline.run: === REAL mode | 2 samples | antibiotic=imipenem | data=/tmp/x ===",
                "2026-01-01 00:00:01 INFO papipeline.run: --- stage validation ---",
                "2026-01-01 00:00:02 INFO papipeline.run: --- stage annotation ---",
                "2026-01-01 00:00:03 INFO papipeline.run: --- stage mlst ---",
                "2026-01-01 00:00:04 INFO papipeline.run: --- stage amr ---",
                "2026-01-01 00:00:05 INFO papipeline.run: --- stage virulence ---",
                "2026-01-01 00:00:06 INFO papipeline.run: --- stage variants ---",
                "2026-01-01 00:00:07 INFO papipeline.run: --- stage cohort_variants ---",
                "2026-01-01 00:00:08 INFO papipeline.run: --- stage pangenome ---",
                "2026-01-01 00:00:09 ERROR papipeline.run: Stage pangenome ended FAILED: "
                "ToolNotAvailableError: Stage 7 (pangenome) cannot verify the toolchain: "
                "panaroo (python package not importable); cd-hit (not on PATH). Both are required.",
                "2026-01-01 00:00:10 INFO papipeline.run: oprD locus SYN_001: resolved (coverage 99.0%)",
                "2026-01-01 00:00:11 INFO papipeline.run: oprD locus SYN_002: refused:insufficient_coverage (coverage 60.0%)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    for stage in (
        "validation", "annotation", "mlst", "amr",
        "virulence", "variants", "cohort_variants", "pangenome",
    ):
        (logs / f"stage_{stage}.log").write_text("log line\n", encoding="utf-8")

    (intermediate / "oprd_structural_calls.tsv").write_text(
        "sample_id\tstructural_verdict\tlesion_type\tidentity_pct\tcoverage_pct\ttruncation_aa\n"
        "SYN_001\tintact\tno_lesion\t97.0\t100\t0\n"
        "SYN_002\tdisrupted\tframeshift_with_premature_stop\t95.0\t100\t42\n",
        encoding="utf-8",
    )
    (intermediate / "reuse_provenance.tsv").write_text(
        "sample_id\treused\treason\tdetail\tsource_dir\tbakta_version\n"
        "SYN_001\ttrue\tverified\treused without invoking Bakta\t/tmp/bakta/SYN_001\t1.12.1\n"
        "SYN_002\ttrue\tverified\treused without invoking Bakta\t/tmp/bakta/SYN_002\t1.12.1\n",
        encoding="utf-8",
    )
    # GUARD 2: 0 bytes = 0 invocations.
    (guards / "calls.log").write_text("", encoding="utf-8")
    (guards / "watcher.log").write_text("# heartbeat t=1s samples=20 detections=0\n", encoding="utf-8")
    return root


@pytest.fixture
def delivery_root(tmp_path: Path) -> Path:
    return build_synthetic_delivery(tmp_path / "delivery")


@pytest.fixture
def delivery_client(app_factory, delivery_root):
    return app_factory(results_root=delivery_root)


# ---------------------------------------------------------------------------
# Layout detection (fails before: `NoResultsSource`)
# ---------------------------------------------------------------------------


def test_delivery_layout_is_detected(delivery_root):
    opened, _probes, detail = source_module.detect(flag_root=delivery_root)
    assert opened.kind == "bundle"
    assert opened.layout == "delivery"
    assert opened.stage_dir == (delivery_root / "artifacts" / "stage_tables").resolve()
    assert detail["run_manifest"] is None
    assert detail["derived_stage_state"]["records"] is True


def test_delivery_layout_is_detected_via_bundle_flag(delivery_root):
    opened, _probes, _detail = source_module.detect(bundle_flag=delivery_root)
    assert opened.layout == "delivery"


def test_assumed_bundle_layout_still_opens(bundle_client):
    """The assumed `02_stage_outputs/` layout is not regressed."""
    _client, state = bundle_client
    assert state.kind == "bundle"
    assert getattr(state.source, "layout", "assumed") == "assumed"


# ---------------------------------------------------------------------------
# No manifest: per-stage state derived from the logs
# ---------------------------------------------------------------------------


def test_no_manifest_stage_state_is_derived(delivery_client):
    client, state = delivery_client
    assert state.manifest is None
    assert state.manifest_writer == "none"
    items = {item["name"]: item for item in client.get("/api/stages").json()["items"]}
    for stage in (
        "validation", "annotation", "mlst", "amr",
        "virulence", "variants", "cohort_variants",
    ):
        assert items[stage]["state"] == "completed", stage
        assert items[stage]["present"] is True, stage

    failed = items["pangenome"]
    assert failed["state"] == "failed"
    assert "ToolNotAvailableError" in failed["reason"]
    assert "cd-hit" in failed["reason"]

    for stage in (
        "recombination", "phylogeny", "similarity", "phenotype",
        "gwas", "convergence", "cooccurrence", "reporting",
    ):
        assert items[stage]["state"] == "not_run", stage
        assert items[stage]["reason"].strip(), stage
        assert "blocked by the failure of stage pangenome" in items[stage]["reason"]


def test_run_mode_and_cohort_come_from_the_log(delivery_client):
    client, _state = delivery_client
    run = client.get("/api/run").json()
    assert run["run_mode"] == "REAL"
    assert run["antibiotic"] == "imipenem"
    assert run["n_samples"] == 2
    assert client.get("/api/health").json()["mode"] == "REAL"
    counts = client.get("/api/run/counts").json()
    assert counts["n_stages_completed"] == 7
    assert counts["n_stages_failed"] == 1
    assert counts["n_stages_not_run"] == 8


def test_stage_after_failure_is_never_rendered_not_assessed(delivery_client):
    client, _state = delivery_client
    items = {item["name"]: item for item in client.get("/api/stages").json()["items"]}
    # gwas is a REAL-refusing stage, but it never ran, so it must not be
    # relabelled `refused` with a reason that never applied to it.
    assert items["gwas"]["state"] == "not_run"
    assert "REAL engine" not in items["gwas"]["reason"]


# ---------------------------------------------------------------------------
# oprD: the real structural-call shape and the separate coverage gate
# ---------------------------------------------------------------------------


def test_oprd_reads_the_structural_call_table(delivery_client):
    client, _state = delivery_client
    payload = client.get("/api/oprd").json()
    assert payload["available"] is True
    by_state: Dict[str, int] = {}
    for item in payload["items"]:
        by_state[item["display_state"]] = by_state.get(item["display_state"], 0) + 1
        assert item["verdict"] == "resolved"
        assert item["display_state"] != "absent"
    assert by_state == {"intact": 1, "disrupted": 1}
    # The structural-call table is the tenth probe, not the first: every
    # contracted path was tried and named first.
    assert payload["probes"][-1]["ok"] is True
    assert payload["probes"][-1]["path"].endswith("oprd_structural_calls.tsv")


def test_locus_coverage_gate_is_reported_separately(delivery_client):
    client, _state = delivery_client
    gate = client.get("/api/oprd").json()["locus_coverage"]
    assert gate["available"] is True
    assert {item["verdict"] for item in gate["items"]} == {
        "resolved",
        "refused:insufficient_coverage",
    }
    assert "different instrument" in gate["note"]


def test_display_state_for_no_lesion_is_intact():
    from dashboard.server.oprd import display_state_for
    from papipeline.adapters.oprd_locus import Verdict

    assert display_state_for(Verdict.RESOLVED, "no_lesion") == "intact"
    assert display_state_for(Verdict.RESOLVED, "frameshift_with_premature_stop") == "disrupted"
    assert display_state_for(Verdict.INSUFFICIENT) == "not_assessed"


# ---------------------------------------------------------------------------
# Bakta executions: 0, derived from the reuse record (never a bare zero)
# ---------------------------------------------------------------------------


def test_bakta_executions_zero_is_derived_from_reuse(delivery_client):
    client, _state = delivery_client
    bakta = client.get("/api/provenance/bakta").json()
    assert bakta["bakta_executions"] == 0
    assert bakta["bakta_executions_label"] == "Bakta executions: 0"
    assert bakta["n_reused"] == 2
    assert bakta["n_not_reused"] == 0
    # The tripwire log is the independent check, and it is empty.
    logs = {item["name"]: item for item in client.get("/api/provenance/logs").json()["items"]}
    assert logs["tripwire"]["present"] is True
    assert logs["tripwire"]["n_lines"] == 0


# ---------------------------------------------------------------------------
# Absent stages are named, never zeros or blanks
# ---------------------------------------------------------------------------


def test_phenotype_absent_reads_not_assessed(delivery_client):
    client, _state = delivery_client
    for row in client.get("/api/isolates").json()["items"]:
        assert row["phenotype_sir"] == "not assessed"
        assert row["lineage"] == "not assessed"


def test_tree_absence_names_stage_nine(delivery_client):
    client, _state = delivery_client
    tree = client.get("/api/tree").json()
    assert tree["present"] is False
    assert tree["reason"].startswith("not produced")
    assert "Stage 9 (phylogeny)" in tree["reason"]


def test_delivery_bundle_is_not_modified(delivery_client, delivery_root):
    client, state = delivery_client
    assert state.source is not None, "the bundle did not open, so nothing was proven"
    before = snapshot_tree(delivery_root)
    for path in (
        "/api/health", "/api/run", "/api/stages", "/api/isolates", "/api/oprd",
        "/api/provenance", "/api/provenance/bakta", "/api/tree", "/api/reports",
    ):
        assert client.get(path).status_code == 200, path
    assert snapshot_tree(delivery_root) == before


# ---------------------------------------------------------------------------
# The REAL bundle, when it is present
# ---------------------------------------------------------------------------


@pytest.fixture
def real_bundle() -> Path:
    if not (REAL_BUNDLE / "artifacts" / "logs" / "full_run.log").is_file():
        pytest.skip(
            "the REAL 10-isolate delivery bundle is not present; set "
            "PA_DASH_REAL_BUNDLE to its path"
        )
    return REAL_BUNDLE


@pytest.fixture
def real_client(app_factory, real_bundle):
    return app_factory(results_root=real_bundle, state_dir=None)


def test_real_bundle_opens_and_is_a_bundle(real_bundle):
    opened, _probes, detail = source_module.detect(flag_root=real_bundle)
    assert opened.kind == "bundle"
    assert opened.layout == "delivery"
    assert detail["run_manifest"] is None
    assert detail["derived_stage_state"]["run_mode"] == "REAL"


def test_real_bundle_stage_table(real_client):
    client, _state = real_client
    items = {item["name"]: item for item in client.get("/api/stages").json()["items"]}
    completed = [name for name, item in items.items() if item["state"] == "completed"]
    # `cohort_variants` left this set when `an_calls` joined the contract.
    # This bundle was written before the column existed, so the reader
    # refuses it rather than present a table that cannot state its own
    # denominator -- the refusal test_report_reads_run_tables.py asserts.
    # See docs/data_contract.md, "Adding a data contract change". A re-run
    # of the delivery bundle puts it back in this set.
    assert set(completed) == {
        "validation", "annotation", "mlst", "amr", "virulence",
        "variants",
    }
    assert items["cohort_variants"]["state"] == "not_assessed"
    assert "missing required columns" in items["cohort_variants"]["reason"]
    assert "an_calls" in items["cohort_variants"]["reason"]
    assert items["pangenome"]["state"] == "failed"
    assert "ToolNotAvailableError" in items["pangenome"]["reason"]
    for name in (
        "recombination", "phylogeny", "similarity", "phenotype",
        "gwas", "convergence", "cooccurrence", "reporting",
    ):
        assert items[name]["state"] == "not_run"
        assert "Stage 7 (pangenome)" in items[name]["reason"]


def test_real_bundle_oprd_two_intact_eight_disrupted(real_client):
    client, _state = real_client
    payload = client.get("/api/oprd").json()
    assert payload["available"] is True
    states: Dict[str, int] = {}
    lesions: Dict[str, int] = {}
    for item in payload["items"]:
        states[item["display_state"]] = states.get(item["display_state"], 0) + 1
        lesions[item["lesion_type"]] = lesions.get(item["lesion_type"], 0) + 1
        assert item["display_state"] != "absent"
    assert states == {"disrupted": 8, "intact": 2}
    assert lesions == {"frameshift_with_premature_stop": 8, "no_lesion": 2}


def test_real_bundle_locus_coverage_refused_four(real_client):
    client, _state = real_client
    gate = client.get("/api/oprd").json()["locus_coverage"]
    assert gate["available"] is True
    refused = [item for item in gate["items"] if item["verdict"].startswith("refused:")]
    assert len(refused) == 4
    assert all(item["verdict"] == "refused:insufficient_coverage" for item in refused)


def test_real_bundle_bakta_executions_zero(real_client):
    client, _state = real_client
    bakta = client.get("/api/provenance/bakta").json()
    assert bakta["bakta_executions"] == 0
    assert bakta["bakta_executions_label"] == "Bakta executions: 0"
    assert bakta["n_reused"] == 10
    assert bakta["n_not_reused"] == 0


def test_real_bundle_cohort_is_ten_and_killed_isolate_is_still_counted(real_client):
    client, _state = real_client
    iso = client.get("/api/isolates", params={"limit": 100}).json()
    cohort = {row["sample_id"] for row in iso["items"]}
    assert len(cohort) == 10
    assert iso["membership_source"] == "union_of_tables"
    # The isolate whose variant calling was killed has no regulator rows and no
    # per-isolate variant rows, yet remains a cohort member.
    reg_rows = []
    for offset in (0, 1000):
        reg_rows.extend(
            client.get(
                "/api/tables/regulators", params={"limit": 1000, "offset": offset}
            ).json()["items"]
        )
    in_regulators = {row["sample_id"] for row in reg_rows}
    assert len(in_regulators) == 9
    assert len(cohort - in_regulators) == 1
    for row in iso["items"]:
        assert row["phenotype_sir"] == "not assessed"
        assert row["oprd_state"] in ("intact", "disrupted")


def test_real_bundle_tree_is_not_produced(real_client):
    client, _state = real_client
    tree = client.get("/api/tree").json()
    assert tree["present"] is False
    assert tree["reason"].startswith("not produced")
    assert "Stage 9 (phylogeny)" in tree["reason"]


def test_real_bundle_is_not_modified(real_client, real_bundle):
    client, state = real_client
    assert state.source is not None, "the bundle did not open, so nothing was proven"
    before = snapshot_tree(real_bundle)
    for path in (
        "/api/health", "/api/run", "/api/stages", "/api/isolates",
        "/api/oprd", "/api/provenance", "/api/tree",
    ):
        assert client.get(path).status_code == 200, path
    assert snapshot_tree(real_bundle) == before
