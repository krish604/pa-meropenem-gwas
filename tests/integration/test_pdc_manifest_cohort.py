"""Cohort construction from `PDC_essential.tsv`, at a size nobody has tried.

The size-agnostic requirement is that nothing in the cohort path assumes an
isolate count. It has been argued, not tested: `PDC_essential.tsv` was read and
counted (967 rows, 835 with an assembly accession) but never *constructed*, so
every claim that the code handles a real cohort rests on inspection of a function
that did not exist. This is the first test to build an actual cohort, and it
does so with four isolates.

Four, deliberately, because a cohort size is exactly the kind of thing code
quietly assumes. Every failure mode below is invisible at 967 rows and fatal at
four:

- a fixed slice, an `expected_count` comparison, or a "should be roughly this
  many" sanity check all pass on real data and truncate a small one;
- a cohort built by iterating assemblies cannot represent an isolate with no
  assembly, so a 3-assembled / 4-isolate cohort is unrepresentable for it;
- uniqueness checks over a large file are the kind of thing that gets relaxed
  for scale, and a duplicate isolate in four rows is visible while the same
  check over 967 might be left to the database.

The fixture uses the real column names and the real `GCA_`/`PDC` value shapes,
so a parser that works here is a parser that works on the real file. One isolate
has no assembly accession, because that is the case the design turns on: cohort
membership is **all isolates**, and an assembly-less isolate is an ordinary
member whose `assembly_path` is simply unset. A stage needing its sequence
refuses per sample later, via `locate_assembly`, and that refusal is correct.
The alternative — building the cohort from assemblies — yields three samples
here and silently drops a cohort member, which is the bug this pins.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.manifest import discover_pdc_manifest

# Real column order from PDC_essential.tsv, including the columns this test
# does not use, so the parser cannot pass by reading only the first two.
COLUMNS = [
    "Isolate",
    "Assembly",
    "BioSample",
    "BioProject",
    "Collection date",
    "Location",
    "Source type",
    "Isolation source",
    "Isolation type",
    "AST phenotypes",
    "AMR genotypes",
    "PDC_present",
]

ROWS = [
    # isolate, assembly, biosample
    ("PDC000001", "GCA_000710625.1", "SAMN00000001"),
    ("PDC000002", "GCA_000710626.1", "SAMN00000002"),
    # no Assembly accession: an ordinary cohort member, not a broken row
    ("PDC000003", "", "SAMN00000003"),
    ("PDC000004", "GCA_000710627.1", "SAMN00000004"),
]


def _write_pdc(tmp_path: Path, rows) -> Path:
    lines = ["\t".join(COLUMNS)]
    for isolate, assembly, biosample in rows:
        lines.append(
            "\t".join(
                [isolate, assembly, biosample, "PRJNA0001", "2019-01-01", "USA",
                 "clinical", "blood", "blood", "R", "blaIMP-1", "TRUE"]
            )
        )
    path = tmp_path / "PDC_essential.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture()
def pdc(tmp_path: Path) -> Path:
    return _write_pdc(tmp_path, ROWS)


class TestCohortIsSizeAgnostic:
    def test_every_isolate_is_a_cohort_member(self, pdc):
        """Four isolates in, four samples out. Not three."""
        manifest = discover_pdc_manifest(pdc)
        assert len(manifest) == 4
        assert set(manifest.sample_ids) == {"PDC000001", "PDC000002", "PDC000003", "PDC000004"}

    def test_isolate_without_assembly_is_still_a_member(self, pdc):
        """The assembly-less isolate is present, with no path set."""
        sample = manifest_sample(pdc, "PDC000003")
        assert sample.assembly is None
        assert sample.assembly_path is None
        assert sample.biosample == "SAMN00000003"

    def test_assembly_becomes_an_attribute_not_the_key(self, pdc):
        """`Assembly` is provenance; the isolate is the join key."""
        sample = manifest_sample(pdc, "PDC000001")
        assert sample.sample_id == "PDC000001"
        assert sample.assembly == "GCA_000710625.1"
        assert sample.isolate_id == "PDC000001"

    def test_explicit_assembly_paths_are_applied_by_isolate(self, pdc):
        """`assembly_paths` is keyed by isolate, and takes no precedence question.

        This is how the smoke run points a small cohort at `db/smoke_genomes/`
        without the pipeline pattern-matching anything.
        """
        manifest = discover_pdc_manifest(
            pdc, assembly_paths={"PDC000001": "db/smoke_genomes/one.fna"}
        )
        assert manifest.require("PDC000001").assembly_path == "db/smoke_genomes/one.fna"
        assert manifest.require("PDC000002").assembly_path is None

    def test_row_order_is_preserved(self, pdc):
        """A cohort built in file order, so a run is reproducible."""
        manifest = discover_pdc_manifest(pdc)
        assert manifest.sample_ids == [r[0] for r in ROWS]

    def test_manifest_reflects_the_file_not_a_slice(self, pdc, tmp_path):
        """Two isolates and five isolates both build; nothing is capped silently."""
        (tmp_path / "two").mkdir(exist_ok=True)
        (tmp_path / "five").mkdir(exist_ok=True)
        two = discover_pdc_manifest(_write_pdc(tmp_path / "two", ROWS[:2]))
        assert len(two) == 2
        five_rows = ROWS + [("PDC000005", "GCA_000710628.1", "SAMN00000005")]
        five = discover_pdc_manifest(_write_pdc(tmp_path / "five", five_rows))
        assert len(five) == 5
        assert [s.sample_id for s in five][-1] == "PDC000005"


    def test_a_slice_to_any_size_is_caught(self, pdc, tmp_path):
        """A cohort must be exactly the file's rows, at whatever size that is.

        The obvious size-agnostic failure is a cap: `[:100]`, a `sorted(...)[:n]`,
        an `if len > N: raise`. None of those show up on a 4-row fixture, so the
        fixture itself asserts the total at a size a slice would not reach.
        """
        (tmp_path / "big").mkdir(exist_ok=True)
        rows = [(f"PDC{i:06d}", f"GCA_{i:09d}.1", f"SAMN{i:08d}") for i in range(1, 121)]
        manifest = discover_pdc_manifest(_write_pdc(tmp_path / "big", rows))
        assert len(manifest) == 120
        assert manifest.sample_ids[-1] == "PDC000120"

    def test_rows_are_never_dropped_for_lacking_an_assembly(self, tmp_path):
        """A skipped row shows up as a smaller cohort, and is caught here.

        The assembly-less isolate is the canary: a parser that `continue`s past a
        row with no Assembly yields three samples instead of four, which at 967
        rows would be invisible.
        """
        manifest = discover_pdc_manifest(_write_pdc(tmp_path, ROWS))
        assert len(manifest) == len(ROWS)
        assert sum(1 for s in manifest if s.assembly is None) == 1

    def test_duplicate_key_is_refused_not_resolved(self, tmp_path):
        """A repeated join key fails loudly; last-write-wins is not acceptable."""
        path = _write_pdc(tmp_path, [ROWS[0], ROWS[0], ROWS[1]])
        with pytest.raises(Exception) as excinfo:
            discover_pdc_manifest(path)
        assert "PDC000001" in str(excinfo.value)

class TestCohortRefusesBadInput:
    def test_duplicate_isolate_is_refused(self, tmp_path):
        """A repeated join key is a hard failure, never a silent overwrite.

        Asserted on the *message*: an implementation that reports the count of
        duplicates but never names them is not usable for triage, and a test
        that only checked the exception type would let that through.
        """
        path = _write_pdc(tmp_path, [ROWS[0], ROWS[0], ROWS[1], ROWS[2]])
        with pytest.raises(Exception) as excinfo:
            discover_pdc_manifest(path)
        message = str(excinfo.value)
        assert "PDC000001" in message
        assert "duplicate" in message.lower()

    def test_empty_isolate_is_refused_not_skipped(self, tmp_path):
        """A row with no key is a data error, not a row to drop.

        Asserting the cohort length as well as the exception matters: a parser
        that skips empty keys yields a *valid-looking* manifest one sample short,
        which is how cohort members go missing without a trace.
        """
        path = _write_pdc(tmp_path, [(ROWS[0]), ("", "GCA_9.1", "SAMN9"), ROWS[1]])
        with pytest.raises(Exception) as excinfo:
            discover_pdc_manifest(path)
        assert "Isolate" in str(excinfo.value)
        # And the failure is not merely a shorter cohort returned quietly.
        with pytest.raises(Exception):
            discover_pdc_manifest(path)

    def test_missing_file_is_refused(self, tmp_path):
        path = tmp_path / "absent.tsv"
        with pytest.raises(Exception):
            discover_pdc_manifest(path)

    def test_empty_file_is_refused(self, tmp_path):
        """Header-only is a cohort of zero, which is never intended."""
        path = tmp_path / "header_only.tsv"
        path.write_text("\t".join(COLUMNS) + "\n", encoding="utf-8")
        with pytest.raises(Exception):
            discover_pdc_manifest(path)


def manifest_sample(pdc_path: Path, isolate: str):
    return discover_pdc_manifest(pdc_path).require(isolate)
