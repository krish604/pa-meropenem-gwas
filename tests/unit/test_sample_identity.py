"""Sample identity is the join key for the whole pipeline.

Ticket 05. Every assembly, annotation, variant table and phenotype row is tied
together by ``sample_id``, so a mistake here is not a bug in one stage - it is a
cohort that silently means something else. The properties that matter:

* the identifier is the versioned assembly accession, which is already the genome
  directory name and the FASTA stem, so the assembly join needs no lookup table;
* the *format* rule is configuration-driven, so real accessions and the
  synthetic fixture IDs both validate while a malformed ID is still rejected;
* the PDC isolate / BioSample / isolate name travel with the sample and must
  agree with the manifest one-to-one;
* a contig header is not a sample identifier, and nothing pretends it is.

Written before the implementation.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import DataContractError, SampleIdError
from papipeline.manifest import discover_manifest, manifest_from_rows, validate_sample_id

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"
FIXTURE_METADATA = PIPELINE_ROOT / "test_data" / "metadata" / "sample_metadata.tsv"

GCA = "GCA_000710625.1"
FIXTURE = "TEST_PA_001"


def _science_sample_id() -> dict:
    return load_config(SCIENCE, machine="laptop").raw.get("sample_id") or {}


# ---------------------------------------------------------------------------
# the format rule is configuration
# ---------------------------------------------------------------------------


def test_the_science_config_declares_a_sample_id_format():
    section = _science_sample_id()
    assert section, "science.yaml must declare a sample_id section"
    assert section.get("pattern"), "the format must be a declared pattern"
    re.compile(section["pattern"])


def test_a_versioned_assembly_accession_is_the_canonical_form():
    assert validate_sample_id(GCA, pattern=_science_sample_id()["pattern"]) == GCA


def test_the_synthetic_fixture_ids_remain_valid_under_the_same_rule():
    """The 20 committed fixtures must not have to be renamed."""
    section = _science_sample_id()
    assert validate_sample_id(FIXTURE, pattern=section["pattern"]) == FIXTURE
    for n in range(1, 21):
        assert validate_sample_id(f"TEST_PA_{n:03d}", pattern=section["pattern"])


@pytest.mark.parametrize(
    "bad",
    [
        "",
        " ",
        "GCA_000710625.1 ",           # trailing whitespace
        " GCA_000710625.1",           # leading whitespace
        "GCA 000710625 1",            # spaces inside
        "GCA/000710625.1",            # a path separator
        "GCA_000710625.1;rm -rf /",   # shell metacharacters
        "a",                          # too short
        "x" * 129,                    # too long
    ],
)
def test_a_malformed_identifier_is_rejected_under_the_configured_rule(bad: str):
    with pytest.raises(SampleIdError):
        validate_sample_id(bad, pattern=_science_sample_id()["pattern"])


def test_the_rule_rejects_something_the_built_in_checks_would_allow():
    """A config-driven rule must be able to be stricter, not merely decorative."""
    strict = r"^GCA_\d{6}\.\d+$"
    # passes the built-in alnum/separator check
    validate_sample_id(GCA)
    with pytest.raises(SampleIdError):
        validate_sample_id(GCA, pattern=strict)


def test_a_pattern_that_disagrees_with_the_sample_is_an_error_not_a_pass():
    with pytest.raises(SampleIdError):
        validate_sample_id("GCA_000710625.1", pattern=r"^TEST_\d+$")


# ---------------------------------------------------------------------------
# uniqueness and existence
# ---------------------------------------------------------------------------


def test_a_duplicate_sample_id_is_a_hard_failure():
    rows = [
        {"sample_id": GCA, "assembly_path": "a.fna"},
        {"sample_id": GCA, "assembly_path": "b.fna"},
    ]
    with pytest.raises(DataContractError) as excinfo:
        manifest_from_rows(rows)
    assert GCA in str(excinfo.value)


def test_a_missing_assembly_file_is_a_hard_failure(tmp_path: Path):
    metadata = tmp_path / "metadata"
    metadata.mkdir()
    (metadata / "sample_metadata.tsv").write_text(
        f"sample_id\tassembly_path\n{GCA}\t../genomes/absent.fna\n",
        encoding="utf-8",
    )
    with pytest.raises(DataContractError) as excinfo:
        discover_manifest(metadata, require_files=True)
    assert "absent.fna" in str(excinfo.value)


def test_the_assembly_path_resolves_against_the_metadata_directory(tmp_path: Path):
    (tmp_path / "genomes").mkdir()
    fasta = tmp_path / "genomes" / f"{GCA}.fna"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")
    metadata = tmp_path / "metadata"
    metadata.mkdir()
    (metadata / "sample_metadata.tsv").write_text(
        f"sample_id\tassembly_path\n{GCA}\t../genomes/{GCA}.fna\n", encoding="utf-8"
    )
    manifest = discover_manifest(metadata, require_files=True)
    assert Path(manifest.require(GCA).assembly_path) == fasta


# ---------------------------------------------------------------------------
# PDC provenance travels with the sample
# ---------------------------------------------------------------------------


def test_isolate_and_biosample_are_carried_as_attributes():
    manifest = manifest_from_rows([
        {
            "sample_id": GCA,
            "assembly_path": f"{GCA}.fna",
            "isolate_id": "PDT000034122.1",
            "biosample": "SAMN02673308",
            "isolate_name": "MRSN18971",
        }
    ])
    sample = manifest.require(GCA)
    assert sample.isolate_id == "PDT000034122.1"
    assert sample.biosample == "SAMN02673308"
    assert sample.isolate_name == "MRSN18971"


def test_the_assembly_accession_is_recorded_alongside_the_sample_id():
    manifest = manifest_from_rows([
        {"sample_id": GCA, "assembly_path": f"{GCA}.fna", "assembly": GCA}
    ])
    assert manifest.require(GCA).assembly == GCA


def test_two_isolates_claiming_one_assembly_is_a_hard_failure():
    """1:1 is the point: otherwise provenance is ambiguous."""
    rows = [
        {"sample_id": GCA, "assembly_path": f"{GCA}.fna",
         "assembly": GCA, "isolate_id": "PDT000000001.1", "biosample": "SAMN1"},
        {"sample_id": GCA + "X", "assembly_path": f"{GCA}.fna",
         "assembly": GCA, "isolate_id": "PDT000000002.1", "biosample": "SAMN2"},
    ]
    with pytest.raises(DataContractError) as excinfo:
        manifest_from_rows(rows)
    assert GCA in str(excinfo.value)


def test_one_isolate_claiming_two_assemblies_is_a_hard_failure():
    rows = [
        {"sample_id": GCA, "assembly_path": f"{GCA}.fna", "assembly": GCA,
         "isolate_id": "PDT000000001.1"},
        {"sample_id": "GCA_000000002.1", "assembly_path": "b.fna",
         "assembly": "GCA_000000002.1", "isolate_id": "PDT000000001.1"},
    ]
    with pytest.raises(DataContractError) as excinfo:
        manifest_from_rows(rows)
    assert "PDT000000001.1" in str(excinfo.value)


def test_one_biosample_claiming_two_assemblies_is_a_hard_failure():
    rows = [
        {"sample_id": GCA, "assembly_path": "a.fna", "assembly": GCA,
         "biosample": "SAMN1"},
        {"sample_id": "GCA_000000002.1", "assembly_path": "b.fna",
         "assembly": "GCA_000000002.1", "biosample": "SAMN1"},
    ]
    with pytest.raises(DataContractError) as excinfo:
        manifest_from_rows(rows)
    assert "SAMN1" in str(excinfo.value)


# ---------------------------------------------------------------------------
# contig headers are not sample identifiers
# ---------------------------------------------------------------------------


def test_the_fixture_contig_header_is_not_a_sample_id():
    """The fixtures use <sample_id>_contig1; real assemblies use NCBI names.

    Nothing may build a manifest from a FASTA header, and a header must never be
    accepted where a sample identifier belongs.
    """
    header = f"{FIXTURE}_contig1 SYNTHETIC_TEST_DATA length=2000 topology=linear"
    record_id = header.split()[0]
    assert record_id != FIXTURE
    # the sample id is recoverable, but by documented stripping, not by equality
    assert record_id.startswith(FIXTURE + "_contig")


def test_discovering_a_manifest_never_reads_the_assemblies():
    """The manifest is metadata only; it must not depend on opening a FASTA."""
    import inspect

    source = inspect.getsource(discover_manifest)
    for forbidden in (".fna", ".fasta", "read_fasta", "seqkit"):
        assert forbidden not in source, (
            f"discover_manifest must not touch sequence files (found {forbidden!r})"
        )


# ---------------------------------------------------------------------------
# the committed fixtures
# ---------------------------------------------------------------------------


def test_the_committed_fixture_manifest_still_loads_and_is_unchanged_in_size():
    manifest = discover_manifest(FIXTURE_METADATA.parent)
    assert len(manifest) == 20
    assert FIXTURE in manifest.sample_ids


def test_no_sample_identifier_is_baked_into_the_configuration():
    """No *concrete* identifier may live in configuration.

    The configured pattern legitimately contains the strings "GCA_" and
    "TEST_", so the check is not for the prefix but for any value that is itself
    a valid identifier: configuration that named a real isolate would be a stray
    record of who is in the study.
    """
    config = load_config(SCIENCE, machine="laptop")
    pattern = re.compile(_science_sample_id()["pattern"])

    def walk(node, path="config"):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "pattern":
                    continue  # the rule, not an instance of it
                walk(value, f"{path}.{key}")
        elif isinstance(node, (list, tuple)):
            for index, value in enumerate(node):
                walk(value, f"{path}[{index}]")
        elif isinstance(node, str) and pattern.fullmatch(node):
            offenders.append(f"{path} = {node!r}")

    offenders: list[str] = []
    walk(config.raw)
    assert not offenders, f"configuration names concrete sample ids: {offenders}"
