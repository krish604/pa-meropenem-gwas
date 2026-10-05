"""Shared pytest fixtures.

The suite runs against the synthetic fixtures in ``test_data/`` and against
purpose-built malformed files written into ``tmp_path``. It contacts no network.

One documented exception to "synthetic only": a handful of tests read the REAL
``PDC_essential.tsv`` cohort table, which is real clinical data (966 isolates,
real BioSample accessions, real AST results) and is deliberately untracked --
``.gitignore`` excludes ``PDC_*.tsv`` under "real dataset: never commit". Those
tests are listed in ``PDC_TABLE_TESTS`` below and are skipped, with a reason,
where that file is absent. Everything else stays synthetic-only.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

PIPELINE_ROOT = Path(__file__).resolve().parents[1]
if str(PIPELINE_ROOT) not in sys.path:
    sys.path.insert(0, str(PIPELINE_ROOT))

from papipeline.config.loader import PipelineConfig, load_config  # noqa: E402
from papipeline.manifest import SampleManifest, discover_manifest  # noqa: E402
from papipeline.models import RunMode  # noqa: E402
from papipeline.testing import SYNTHETIC_BANNER, generate  # noqa: E402

TEST_MODE = RunMode.TEST

#: Tests that read the real, untracked ``PDC_essential.tsv`` cohort table.
#:
#: These assert properties of the REAL cohort -- the isolate count, the
#: location distribution, that the benchmark cohort is not the first N PDC
#: rows. None can be satisfied by synthetic data, and substituting a fixture
#: would make them assert nothing, so they are skipped rather than rewritten.
#:
#: Deliberately an exact node-id allow-list rather than a path or name pattern:
#: a pattern would also swallow the neighbours that do NOT need the table (for
#: example ``test_stratify_keeps_empty_values_as_their_own_group``, which is
#: pure logic and passes), and a skip that hides a passing test is worse than a
#: stale entry. If a test here starts failing for a real reason it will simply
#: be absent from this list and fail loudly, which is the intended failure mode.
PDC_TABLE_TESTS = frozenset(
    {
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsStratified::test_it_spans_many_locations",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsStratified::test_small_strata_are_not_crowded_out",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsStratified::test_it_does_not_collapse_onto_one_hospital",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsStratified::test_unrecorded_location_is_a_stratum_not_a_filter",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsStratified::test_first_fifty_by_order_would_be_worse",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsUsable::test_every_chosen_isolate_has_an_assembly",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsUsable::test_it_is_reproducible",
        "tests/unit/test_cohort_sampling.py::TestTheSampleIsUsable::test_asking_for_more_than_exists_is_refused",
        "tests/unit/test_cohort_sampling.py::TestAgainstRealMetadata::test_location_is_the_strongest_stratifier_available",
        "tests/unit/test_annotation_assembly_root.py::TestTheAccessorItself::test_it_agrees_with_the_manifest_builder",
        "tests/unit/test_execution_semantics_and_benchmark.py::TestBenchmarkCohort::test_the_cohort_is_not_the_first_ten_pdc_rows",
        "tests/integration/test_smoke_overlay.py::TestManifestDiscoveryHonoursTheSmokeGenomeDir::test_smoke_discovery_reads_the_smoke_genome_dir",
        "tests/integration/test_smoke_overlay.py::TestManifestDiscoveryHonoursTheSmokeGenomeDir::test_no_sample_metadata_file_is_written",
        "tests/integration/test_smoke_overlay.py::TestTheCapCountsPreparedAssembliesOnASmokeRun::test_membership_is_untouched_still_967",
    }
)

PDC_TABLE_NAME = "PDC_essential.tsv"

#: Why the tests above cannot run without the table.
PDC_ABSENT_REASON = (
    "{name} needs the real PDC isolate table, which is not in this checkout. "
    "It is real clinical data (966 isolates with real BioSample accessions and "
    "AST results) and .gitignore excludes PDC_*.tsv on purpose, so it is "
    "provisioned out of band like the databases. Provide it at the repository "
    "root to exercise these tests; do not commit it."
)


def _pdc_table_available() -> bool:
    """Is the real cohort table readable?

    Checked both relative to the current working directory and at the
    repository root, because the tests reference it both ways: some pass a bare
    ``Path("PDC_essential.tsv")`` (resolved against the CWD) and one passes
    ``REPO / "PDC_essential.tsv"``.
    """
    return any(
        candidate.is_file()
        for candidate in (Path(PDC_TABLE_NAME), PIPELINE_ROOT / PDC_TABLE_NAME)
    )


def pytest_collection_modifyitems(config, items):
    """Skip the real-cohort tests when the real cohort table is absent."""
    if _pdc_table_available():
        return
    reason = PDC_ABSENT_REASON.format(name=PDC_TABLE_NAME)
    skip = pytest.mark.skip(reason=reason)
    for item in items:
        if item.nodeid in PDC_TABLE_TESTS:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def pipeline_root() -> Path:
    return PIPELINE_ROOT


@pytest.fixture(scope="session")
def config() -> PipelineConfig:
    """The project's real configuration, loaded once per session."""
    return load_config(PIPELINE_ROOT / "config" / "science.yaml")


@pytest.fixture(scope="session")
def test_data_root() -> Path:
    root = PIPELINE_ROOT / "test_data"
    if not (root / "metadata" / "sample_metadata.tsv").exists():
        generate(root)
    return root


@pytest.fixture(scope="session")
def manifest(test_data_root: Path) -> SampleManifest:
    return discover_manifest(test_data_root / "metadata")


@pytest.fixture(scope="session")
def intermediate_root(test_data_root: Path) -> Path:
    return test_data_root / "intermediate"


@pytest.fixture(scope="session")
def phylogeny_dir(test_data_root: Path) -> Path:
    return test_data_root / "phylogeny"


@pytest.fixture(scope="session")
def phenotype_dir(test_data_root: Path) -> Path:
    return test_data_root / "phenotype"


@pytest.fixture
def write_tsv_file(tmp_path: Path):
    """Factory writing raw TSV text to a temp file and returning its path.

    Lets a test construct a deliberately malformed table without checking a
    broken fixture into the repository.
    """

    def _write(name: str, content: str) -> Path:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    return _write


@pytest.fixture
def tsv_header():
    """Factory producing a well-formed header line for the fake tables."""

    def _header(*columns: str) -> str:
        return "\t".join(columns) + "\n"

    return _header


@pytest.fixture
def sample_manifest() -> SampleManifest:
    """A small hand-built manifest, independent of the generated fixtures."""
    from papipeline.models import Sample

    return SampleManifest(
        [
            Sample("TEST_A_01", None, "test"),
            Sample("TEST_A_02", None, "test"),
            Sample("TEST_A_03", None, "test"),
        ]
    )
