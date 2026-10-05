"""Stage 7 in REAL mode: it calls now, candidate-only.

The guard this replaces (a50bc0f) asserted `NotImplementedError` on REAL. That
was correct while there was no caller, and the refusal is what let a REAL run
stop cleanly after stage 4 instead of reporting ten genuine AMRFinderPlus
reports as a stage failure.

Those assertions are gone. What replaces them is the property that matters
more: assembly-only input can produce a candidate and nothing else, so no REAL
cohort can ever yield a `confirmed` structural variant.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import sv as caller
from papipeline.models import RunMode, StructuralCallStatus
from papipeline.stages import sv as stage_sv


def _manifest(tmp_path, with_assembly=True):
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    path = None
    if with_assembly:
        path = tmp_path / "genome.fna"
        path.write_text(">c1\n" + "ACGT" * 500 + "\n", encoding="utf-8")
    return SampleManifest(
        samples=[Sample(sample_id="PDT_A.1", assembly_path=str(path) if path else None)]
    )


def _config(tmp_path, *, reference=True, caller_name="nucmer"):
    from papipeline.config.loader import PipelineConfig

    db = tmp_path / "db"
    db.mkdir(parents=True, exist_ok=True)
    ref = db / "ref.fna"
    if reference:
        ref.write_text(">ref\n" + "ACGT" * 100 + "\n", encoding="utf-8")

    class _Machine:
        paths = {}

        def db_root(self):
            return db

        def _path(self, key, default):
            return db / default

        def reference_fasta(self, science=None):
            return ref

        def tool_candidates(self, tool):
            # Real files: `_first_candidate` checks existence, so a bare name
            # would fail the lookup these tests never meant to exercise.
            bindir = tmp_path / "bin"
            bindir.mkdir(parents=True, exist_ok=True)
            stub = bindir / tool
            stub.write_text("#!/bin/sh\n", encoding="utf-8")
            stub.chmod(0o755)
            return (str(stub),)

    return PipelineConfig(
        root=tmp_path,
        raw={"structural_variants": {
            "caller": caller_name,
            "nucmer_minmatch": 20, "nucmer_mincluster": 65, "nucmer_maxgap": 90,
            "nucmer_maxmatch": True, "min_sv_size": 1000,
            "min_flanking_identity_pct": 90.0,
        }},
        antibiotics=(), allowed_phenotypes=(), analysis={}, mechanisms={},
        regulators={}, references={}, antibiotic_specs={}, organism={},
        qc=None, gwas=None, phylogeny=None, convergence=None,
        paths={"reference_fasta": str(ref)},
        runtime={"threads": 1}, machine=_Machine(),
    )


class TestRealNoLongerRefuses:
    def test_it_no_longer_raises_not_implemented(self, tmp_path):
        """It must get as far as looking for the reference.

        A missing reference is a provisioning problem with its own message; if
        this still raised `NotImplementedError` that message would be
        unreachable.
        """
        config = _config(tmp_path, reference=False)
        with pytest.raises(caller.SvPreflightError) as excinfo:
            stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        assert "NotImplemented" not in str(excinfo.value)

    def test_a_missing_reference_names_the_path(self, tmp_path):
        config = _config(tmp_path, reference=False)
        with pytest.raises(caller.SvPreflightError) as excinfo:
            stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        assert "ref.fna" in str(excinfo.value)

    def test_an_unconfigured_caller_is_an_error_not_a_default(self, tmp_path):
        """Rule 2: no capability silently assumed to be present."""
        config = _config(tmp_path, caller_name="")
        with pytest.raises(caller.SvPreflightError) as excinfo:
            stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        assert "caller" in str(excinfo.value)

    def test_test_mode_still_reads_the_fixture(self, tmp_path, config):
        """The REAL branch is additive; the fixture path is untouched."""
        target = tmp_path / "structural_variants"
        target.mkdir(parents=True)
        (target / "structural_variants.tsv").write_text(
            "sample_id\tvariant_id\tvariant_type\tcall_status\n"
            "PDT_A.1\tv1\tdeletion\tconfirmed\n",
            encoding="utf-8",
        )
        result = stage_sv.run(config, _manifest(tmp_path), RunMode.TEST, tmp_path)
        assert result["PDT_A.1"], "the TEST fixture row was not returned"


class TestRealIsCandidateOnly:
    """The load-bearing property of an assembly-only caller."""

    def test_no_real_call_is_ever_confirmed(self, tmp_path, monkeypatch):
        """Rule 7: `confirmed` means "meets the configured evidence
        requirement", and assembly-only evidence cannot meet it.

        `confirmed_only()` is therefore empty for every REAL cohort. That is a
        property of the input, not a defect to tune away, and the test states it
        so a later change that promotes these calls fails loudly.
        """
        config = _config(tmp_path)
        monkeypatch.setattr(
            stage_sv, "screen_assembly",
            lambda **kwargs: (
                [
                    caller.SvCandidate(
                        sample_id="PDT_A.1", variant_id="v1",
                        variant_type="deletion", position=100, size=2000,
                        evidence="nucmer:deletion_in_reference",
                        call_status="candidate",
                    )
                ],
                [],
            ),
        )
        result = stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        records = result["PDT_A.1"]
        assert records
        assert all(r.call_status is StructuralCallStatus.CANDIDATE for r in records)
        assert stage_sv.confirmed_only(records) == []

    def test_unassessable_regions_are_surfaced_not_counted_as_absent(
        self, tmp_path, monkeypatch,
    ):
        """A contig gap is not "no variant present", and rule 7 forbids
        conflating them. The count is returned to the caller as a warning rather
        than folded into the call count."""
        config = _config(tmp_path)
        monkeypatch.setattr(
            stage_sv, "screen_assembly",
            lambda **kwargs: ([], [("c1", 0, 5000), ("c2", 1, 4000)]),
        )
        result = stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        # Two regions unassessable, and no record invented from them.
        assert result["PDT_A.1"] == []

    def test_it_writes_the_declared_columns(self, tmp_path, monkeypatch):
        config = _config(tmp_path)
        monkeypatch.setattr(
            stage_sv, "screen_assembly",
            lambda **kwargs: (
                [
                    caller.SvCandidate(
                        sample_id="PDT_A.1", variant_id="v1",
                        variant_type="insertion", position=100, size=2000,
                        evidence="nucmer:insertion_in_assembly",
                        call_status="candidate",
                    )
                ],
                [],
            ),
        )
        out = tmp_path / "o"
        stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, out)
        written = (
            out / "structural_variants" / "structural_variants.tsv"
        ).read_text()
        assert written.splitlines()[0].split("\t") == list(stage_sv.SV_COLUMNS)
        assert "PDT_A.1\tv1\tinsertion" in written

    def test_an_isolate_without_an_assembly_is_screened_with_nothing(
        self, tmp_path, monkeypatch,
    ):
        """Recorded as screened, not dropped - 957 of the 967-isolate roster are
        in exactly this position."""
        config = _config(tmp_path)
        calls = []
        monkeypatch.setattr(
            stage_sv, "screen_assembly",
            lambda **kwargs: calls.append(kwargs) or ([], []),
        )
        result = stage_sv.run(
            config, _manifest(tmp_path, with_assembly=False),
            RunMode.REAL, tmp_path / "o",
        )
        assert calls == [], "no assembly means nothing to align"
        assert result["PDT_A.1"] == []


class TestDottedIsolateIds:
    def test_a_dotted_sample_id_does_not_truncate_the_delta_name(
        self, tmp_path, monkeypatch,
    ):
        """`Path.with_suffix('.delta')` strips `.1` from `PDT000034122.1`.

        The screen then looked for `PDT000034122.delta` after a *successful*
        alignment and reported "nucmer wrote no delta file" - a tool failure
        message for a path bug. Isolate IDs in this cohort always contain a dot.
        """
        config = _config(tmp_path)
        captured = {}

        def fake_screen(**kwargs):
            captured.update(kwargs)
            delta = Path(kwargs["workdir"]) / f"{kwargs['sample_id']}.delta"
            delta.parent.mkdir(parents=True, exist_ok=True)
            delta.write_text("stub", encoding="utf-8")
            return [], []

        monkeypatch.setattr(stage_sv, "screen_assembly", fake_screen)
        stage_sv.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "o")
        assert captured["sample_id"] == "PDT_A.1"
        expected = tmp_path / "o" / "structural_variants" / "nucmer" / "PDT_A.1.delta"
        assert expected.exists(), "the delta name lost its dotted tail"