"""Sample identity must hold on the path a run actually takes.

Ticket 05, promoted to `verified`. The pattern and the provenance checks existed
and had unit tests, but `discover_manifest` - the function `run_pipeline` really
calls - ignored both, so the rules were inert in production. Unit tests on a
function nothing reaches are not evidence of a feature.

These go through `run_pipeline` on a cloned project, and assert the two things
the unit tests could not:

* a manifest that violates the rules stops the run;
* nothing is written before it stops. A refusal after annotation is not a
  refusal.

Written before the wiring.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import DataContractError, SampleIdError
from papipeline.run import run_pipeline

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = PIPELINE_ROOT / "test_data"
METADATA = "sample_metadata.tsv"


def _project(tmp_path: Path) -> Path:
    """A self-contained project, so nothing reaches back into the repository."""
    shutil.copytree(PIPELINE_ROOT / "config", tmp_path / "config")
    shutil.copytree(FIXTURES, tmp_path / "test_data")
    return tmp_path / "test_data"


def _rewrite_metadata(test_data: Path, transform) -> None:
    path = test_data / "metadata" / METADATA
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(transform(lines)) + "\n", encoding="utf-8")


def _run(test_data: Path):
    return run_pipeline(
        mode="TEST",
        config=load_config(
            test_data.parent / "config" / "science.yaml", machine="laptop"
        ),
    )


def _assert_nothing_written(project: Path) -> None:
    results = project / "results"
    stage_files = list(results.rglob("stages/*")) if results.exists() else []
    assert not stage_files, f"output was produced before the refusal: {stage_files[:5]}"


def _with_provenance(lines, columns):
    """Return provenance-enriched rows as FIELD LISTS, one value per sample.

    The committed metadata carries `#` comment lines above the real header, so
    the header is located rather than assumed to be line 0. Field lists rather
    than joined strings, so a test can collide one cell without re-parsing.
    """
    head = next(i for i, line in enumerate(lines) if line and not line.startswith("#"))
    header = lines[head].split("\t")
    for column in columns:
        if column not in header:
            header.append(column)
    rows = [header]
    for line in lines[head + 1:]:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        fields += [""] * (len(header) - len(fields))
        sid = fields[0]
        for column in columns:
            fields[header.index(column)] = {
                "assembly": sid,
                "isolate_id": f"PDT_{sid}",
                "biosample": f"SAM_{sid}",
                "isolate_name": f"name_{sid}",
            }[column]
        rows.append(fields)
    return rows


def _render(rows):
    return ["\t".join(r) for r in rows]


def _data_row_indices(rows):
    return [i for i, r in enumerate(rows) if i and r[0]]


# ---------------------------------------------------------------------------
# the happy path still works
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_the_committed_cohort_runs(tmp_path: Path):
    """With only negative cases, a rule that refuses everything would look
    like a rule that works."""
    test_data = _project(tmp_path)
    assert _run(test_data).n_samples == 20


# ---------------------------------------------------------------------------
# the format rule, on the real path
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_a_sample_id_that_breaks_the_configured_pattern_stops_the_run(tmp_path: Path):
    test_data = _project(tmp_path)
    # A well-formed, filesystem-safe identifier the configured pattern rejects:
    # neither a GCA accession nor a fixture id.
    _rewrite_metadata(
        test_data,
        lambda lines: [lines[0]]
        + [line.replace("TEST_PA_001", "isolate-17", 1) for line in lines[1:]],
    )
    with pytest.raises(SampleIdError) as excinfo:
        _run(test_data)
    assert "isolate-17" in str(excinfo.value)
    _assert_nothing_written(tmp_path)


@pytest.mark.slow
def test_a_duplicate_sample_id_stops_the_run(tmp_path: Path):
    test_data = _project(tmp_path)
    _rewrite_metadata(
        test_data,
        lambda lines: lines
        + [next(l for l in lines[1:] if l.startswith("TEST_PA_001"))],
    )
    with pytest.raises(DataContractError):
        _run(test_data)
    _assert_nothing_written(tmp_path)


# ---------------------------------------------------------------------------
# provenance must map 1:1
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_two_isolates_claiming_one_assembly_stop_the_run(tmp_path: Path):
    """Provenance that cannot be attributed to one assembly is not provenance."""
    test_data = _project(tmp_path)
    columns = ("assembly", "isolate_id", "biosample", "isolate_name")

    def collide(lines):
        rows = _with_provenance(lines, columns)
        data = _data_row_indices(rows)
        column = rows[0].index("assembly")
        rows[data[1]][column] = rows[data[0]][column]
        return _render(rows)

    _rewrite_metadata(test_data, collide)
    with pytest.raises(DataContractError) as excinfo:
        _run(test_data)
    assert "assembly" in str(excinfo.value)
    _assert_nothing_written(tmp_path)


@pytest.mark.slow
def test_one_isolate_claiming_two_assemblies_stops_the_run(tmp_path: Path):
    test_data = _project(tmp_path)

    def collide(lines):
        rows = _with_provenance(lines, ("isolate_id",))
        data = _data_row_indices(rows)
        column = rows[0].index("isolate_id")
        rows[data[1]][column] = rows[data[0]][column]
        return _render(rows)

    _rewrite_metadata(test_data, collide)
    with pytest.raises(DataContractError) as excinfo:
        _run(test_data)
    assert "PDT" in str(excinfo.value)
    _assert_nothing_written(tmp_path)


# ---------------------------------------------------------------------------
# the contig-header rule, stated rather than assumed
# ---------------------------------------------------------------------------


def test_a_contig_header_is_not_a_sample_id():
    """The fixtures use <sample_id>_contig1; nothing may derive identity from one."""
    from papipeline.io.fasta import read_fasta

    records = list(read_fasta(FIXTURES / "genomes" / "TEST_PA_001.fna"))
    assert records, "the fixture should parse"
    record_ids = [r.record_id for r in records]
    assert "TEST_PA_001" not in record_ids
    assert all(rid.startswith("TEST_PA_001_contig") for rid in record_ids)
