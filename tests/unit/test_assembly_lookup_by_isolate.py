"""`locate_assembly` must follow the cohort's actual key, not an old convention.

`discover_pdc_manifest` keys the cohort on `Isolate`, so `sample_id` is a PDC
isolate name (`PDT000034122.1`). `locate_assembly` searched only
`data/{sample_id}/`, which is correct while sample ids are assembly accessions
and wrong now: the real directories are named by accession
(`data/GCA_000710625.1/`), so the search never finds anything and every real
sample refuses as if the genome were missing.

The `assembly` attribute is already carried on every `Sample` for exactly this
reason, but was never consulted - the accessor existed and nothing read it,
which is the "method with no caller" pattern this repo has now hit repeatedly.

Worth stating plainly: this fails *quietly*. Every refusal is a legitimate
"assembly not found", so nothing crashes and no test that only counts
refusals would notice. The cohort resolves, the cap passes, and the run reports
zero genomes. That is worse than an error.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.assemblies import locate_assembly
from papipeline.errors import PipelineError
from papipeline.manifest import discover_pdc_manifest


def _pdc(tmp_path: Path) -> Path:
    lines = ["\t".join(["Isolate", "Assembly", "BioSample", "PDC_present"])]
    lines.append("PDT000034122.1\tGCA_000710625.1\tSAMN00000001\tTRUE")
    lines.append("PDT000050773.2\tGCA_000710626.1\tSAMN00000002\tTRUE")
    path = tmp_path / "pdc.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _real_layout(root: Path) -> None:
    """The layout actually on disk: accession-named dirs, sample-named files."""
    d = root / "GCA_000710625.1"
    d.mkdir(parents=True)
    (d / "GCA_000710625.1_MRSN18971scaf_genomic.fna").write_text(">c\nACGT\n", encoding="utf-8")
    d2 = root / "GCA_000710626.1"
    d2.mkdir(parents=True)
    (d2 / "GCA_000710626.1.fna").write_text(">c\nACGT\n", encoding="utf-8")


class TestAssemblyIsFoundByItsAccession:
    def test_isolate_keyed_sample_resolves_through_its_assembly(self, tmp_path):
        """The real layout, with the real PDC isolate key.

        Neither the isolate id nor an exact `{sample_id}.fna` can see this file:
        the directory is the accession and the file carries a lab-strain suffix.
        """
        _real_layout(tmp_path)
        manifest = discover_pdc_manifest(_pdc(tmp_path))
        resolved = locate_assembly(manifest.require("PDT000034122.1"), tmp_path)
        assert resolved.name == "GCA_000710625.1_MRSN18971scaf_genomic.fna"
        assert resolved.parent.name == "GCA_000710625.1"

    def test_both_samples_resolve(self, tmp_path):
        _real_layout(tmp_path)
        manifest = discover_pdc_manifest(_pdc(tmp_path))
        for sample in manifest:
            assert locate_assembly(sample, tmp_path).is_file()

    def test_explicit_path_still_wins(self, tmp_path):
        """Rule 1 is unchanged: a declared path is used with no pattern matching."""
        _real_layout(tmp_path)
        declared = tmp_path / "smoke" / "one.fna"
        declared.parent.mkdir()
        declared.write_text(">c\nACGT\n", encoding="utf-8")
        manifest = discover_pdc_manifest(
            _pdc(tmp_path), assembly_paths={"PDT000034122.1": str(declared)}
        )
        assert locate_assembly(manifest.require("PDT000034122.1"), tmp_path) == declared

    def test_test_layout_by_sample_id_still_works(self, tmp_path):
        """The TEST convention is not sacrificed to fix the real one."""
        d = tmp_path / "PDT000050773.2"
        d.mkdir(parents=True)
        (d / "PDT000050773.2.fna").write_text(">c\nACGT\n", encoding="utf-8")
        manifest = discover_pdc_manifest(_pdc(tmp_path))
        resolved = locate_assembly(manifest.require("PDT000050773.2"), tmp_path)
        assert resolved.name == "PDT000050773.2.fna"


class TestRefusalsStillHold:
    def test_missing_assembly_refuses(self, tmp_path):
        """An isolate with no accession refuses, and says which sample."""
        lines = ["\t".join(["Isolate", "Assembly", "BioSample", "PDC_present"])]
        lines.append("PDT000099999.9\t\tSAMN00000009\tTRUE")
        path = tmp_path / "pdc_none.tsv"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        manifest = discover_pdc_manifest(path)
        with pytest.raises(PipelineError) as excinfo:
            locate_assembly(manifest.require("PDT000099999.9"), tmp_path)
        assert "PDT000099999.9" in str(excinfo.value)

    def test_ambiguity_still_refuses(self, tmp_path):
        """Two candidates in one accession directory is a refusal, not a pick."""
        d = tmp_path / "GCA_000710625.1"
        d.mkdir(parents=True)
        (d / "GCA_000710625.1_scaf1.fna").write_text(">c\nACGT\n", encoding="utf-8")
        (d / "GCA_000710625.1_scaf2.fna").write_text(">c\nACGT\n", encoding="utf-8")
        manifest = discover_pdc_manifest(_pdc(tmp_path))
        with pytest.raises(PipelineError):
            locate_assembly(manifest.require("PDT000034122.1"), tmp_path)

    def test_no_repository_wide_scan(self, tmp_path):
        """A genome elsewhere in the tree must not be attributed to a sample."""
        elsewhere = tmp_path / "somewhere" / "deep"
        elsewhere.mkdir(parents=True)
        (elsewhere / "GCA_000710625.1_x.fna").write_text(">c\nACGT\n", encoding="utf-8")
        manifest = discover_pdc_manifest(_pdc(tmp_path))
        with pytest.raises(PipelineError):
            locate_assembly(manifest.require("PDT000034122.1"), tmp_path)
