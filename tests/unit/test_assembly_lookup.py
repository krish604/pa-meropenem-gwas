"""`locate_assembly` finds a sample's file, and refuses rather than guesses.

The seam is agreed: `papipeline.assemblies.locate_assembly(sample, data_root)`.
A pure function, no config coupling, so both `adapters/bakta.py` and the
workflow can call it with what they already have.

It replaces three things that disagreed with each other and with the real data:

  - `adapters/bakta.py:assembly_for`, which probes for an exact
    `{sample_id}{suffix}` and therefore cannot see
    `GCA_000710625.1_MRSN18971scaf_genomic.fna`;
  - `Snakefile:154`'s `DATA_ROOT / "genomes"` glob, which resolves to nothing
    against the real layout because `data/genomes/` is empty;
  - `PipelineConfig.genomes_dir`, which had no caller at all.

The three rules, in order: an explicit `assembly_path` is used directly; otherwise
the sample's accession directory is searched for a file matching the accession;
otherwise refuse. No repository-wide recursive scan - that belongs to
`pilot/cohort.py` and to nothing else.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.assemblies import locate_assembly
from papipeline.errors import PipelineError
from papipeline.models import Sample

#: The real filename for this accession, verbatim from `data/`. The stem is NOT
#: the accession, which is the whole reason an exact-filename probe fails.
REAL_DIR = "GCA_000710625.1"
REAL_FILE = "GCA_000710625.1_MRSN18971scaf_genomic.fna"


def _real_layout(tmp_path: Path) -> Path:
    """A faithful copy of the real per-accession layout."""
    data = tmp_path / "data"
    (data / REAL_DIR).mkdir(parents=True)
    (data / REAL_DIR / REAL_FILE).write_text(">c\nACGT\n", encoding="utf-8")
    return data


def _test_layout(tmp_path: Path) -> Path:
    """The TEST layout: a flat genomes/ dir named exactly by sample id."""
    data = tmp_path / "data"
    (data / "genomes").mkdir(parents=True)
    (data / "genomes" / "TEST_PA_001.fna").write_text(">c\nACGT\n", encoding="utf-8")
    return data


class TestTheRealLayout:
    def test_it_resolves_the_real_filename_for_a_real_sample_id(self, tmp_path):
        """The case no exact-filename probe can satisfy.

        `data/GCA_000710625.1/GCA_000710625.1_MRSN18971scaf_genomic.fna` for
        `sample_id=GCA_000710625.1`.
        """
        data = _real_layout(tmp_path)
        found = locate_assembly(Sample(sample_id=REAL_DIR), data)
        assert found == data / REAL_DIR / REAL_FILE
        assert found.exists()

    def test_the_old_glob_would_have_found_nothing(self, tmp_path):
        """Documents the bug ticket 14 would have inherited.

        The workflow globbed `DATA_ROOT / "genomes"`, which is empty in the real
        layout - so a REAL run resolved an empty assembly list and the failure
        would have looked like "no assemblies", not "no assemblies configured".
        """
        data = _real_layout(tmp_path)
        assert not (data / "genomes").exists()
        assert sorted(data.glob("genomes/*.fna")) == []

    def test_it_still_serves_the_test_layout(self, tmp_path):
        """The TEST fixtures are named exactly by sample id, and must keep working."""
        data = _test_layout(tmp_path)
        found = locate_assembly(Sample(sample_id="TEST_PA_001"), data)
        assert found == data / "genomes" / "TEST_PA_001.fna"


class TestAnExplicitPathIsUsedDirectly:
    def test_an_absolute_existing_path_wins(self, tmp_path):
        """Rule 1: no pattern matching when the manifest says where it is.

        This is the path the bounded smoke run will use, pointing at
        db/smoke_genomes/ - so it must not depend on any directory layout.
        """
        target = tmp_path / "smoke" / "TEST_PA_009.fna"
        target.parent.mkdir(parents=True)
        target.write_text(">c\nACGT\n", encoding="utf-8")
        found = locate_assembly(
            Sample(sample_id="TEST_PA_009", assembly_path=str(target)),
            tmp_path / "data",
        )
        assert found == target

    def test_a_bare_filename_is_resolved_against_the_data_root(self, tmp_path):
        data = _real_layout(tmp_path)
        found = locate_assembly(
            Sample(sample_id=REAL_DIR, assembly_path=REAL_FILE), data
        )
        assert found == data / REAL_DIR / REAL_FILE


class TestItRefusesRatherThanGuesses:
    def test_nothing_found_is_refused(self, tmp_path):
        data = tmp_path / "data"
        data.mkdir()
        with pytest.raises(PipelineError) as excinfo:
            locate_assembly(Sample(sample_id="GCA_999999999.9"), data)
        message = str(excinfo.value)
        assert "GCA_999999999.9" in message, (
            "the refusal must name the sample, or the operator cannot act on it"
        )

    def test_ambiguity_is_refused(self, tmp_path):
        """Two candidate files for one accession is a question, not a choice."""
        data = tmp_path / "data"
        directory = data / REAL_DIR
        directory.mkdir(parents=True)
        (directory / f"{REAL_DIR}_scaffold_one.fna").write_text(">c\nA\n", encoding="utf-8")
        (directory / f"{REAL_DIR}_scaffold_two.fna").write_text(">c\nC\n", encoding="utf-8")

        with pytest.raises(PipelineError) as excinfo:
            locate_assembly(Sample(sample_id=REAL_DIR), data)
        message = str(excinfo.value).lower()
        assert "ambiguous" in message or "more than one" in message, (
            f"an ambiguous lookup must say so, not pick one: {excinfo.value}"
        )

    def test_refusal_does_not_scan_the_repository(self, tmp_path):
        """A miss must not widen into a search.

        If a miss fell back to a recursive scan it could find an unrelated file
        somewhere else and attribute it to this sample - which is the failure
        mode the pilot's rglob has and the pipeline must not inherit.
        """
        data = tmp_path / "data"
        (data / "unrelated").mkdir(parents=True)
        (data / "unrelated" / "GCA_123456789.1_other.fna").write_text(">c\nA\n", encoding="utf-8")
        with pytest.raises(PipelineError):
            locate_assembly(Sample(sample_id="GCA_999999999.9"), data)
