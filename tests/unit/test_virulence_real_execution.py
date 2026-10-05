"""Stage 5 in REAL mode: it screens now, and it refuses for the right reason.

The guard this replaces (7e28840) asserted `NotImplementedError` on REAL. Those
assertions are inverted here deliberately, and the refusal that remains is a
*different* one: a missing database, not a missing caller.
"""

from __future__ import annotations

import pytest

from papipeline.adapters import virulence as vf
from papipeline.models import RunMode
from papipeline.stages import virulence as stage_virulence


def _manifest(tmp_path, with_assembly=True):
    from papipeline.manifest import SampleManifest
    from papipeline.models import Sample

    path = None
    if with_assembly:
        path = tmp_path / "genome.fna"
        path.write_text(">contig\nATGTGTGCGCTGGATCGTAGAGAAAGGCCACTT\n", encoding="utf-8")
    return SampleManifest(samples=[Sample(sample_id="PDT_A", assembly_path=str(path) if path else None)])


def _config_with_vfdb(tmp_path, sequences=True):
    from papipeline.config.loader import PipelineConfig

    db_root = tmp_path / "db"
    (db_root / "vfdb").mkdir(parents=True, exist_ok=True)
    if sequences:
        (db_root / "vfdb" / "sequences").write_text(
            ">vfdb~~~fimV~~~X~~~ (fimV) fimbrial [Type IV pili (VF0082) - "
            "Adherence (VFC0001)] [Pseudomonas aeruginosa]\nACGTACGTACGT\n",
            encoding="utf-8",
        )

    class _Machine:
        paths = {}

        def db_root(self):
            return db_root

        def _path(self, key, default):
            return db_root / default

    return PipelineConfig(
        root=tmp_path,
        raw={"virulence": {
            "sequences_db": "vfdb/sequences",
            "database": "VFDB",
            "blast_program": "blastn",
            "min_identity_pct": 90.0,
            "min_coverage_pct": 90.0,
            "evalue": 1e-5,
            "allow_database_update": False,
        }},
        antibiotics=(), allowed_phenotypes=(), analysis={}, mechanisms={},
        regulators={}, references={}, antibiotic_specs={}, organism={},
        qc=None, gwas=None, phylogeny=None, convergence=None, paths={},
        runtime={"threads": 1}, machine=_Machine(),
    )


class TestRealNoLongerRefuses:
    def test_it_no_longer_raises_not_implemented(self, tmp_path):
        """The whole point of this change.

        It must get as far as *looking* for the database. A missing database is
        a provisioning problem with its own message, and is asserted separately
        below - if this still raised NotImplementedError, that message would be
        unreachable.
        """
        config = _config_with_vfdb(tmp_path, sequences=False)
        with pytest.raises(vf.VirulencePreflightError) as excinfo:
            stage_virulence.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "out")
        assert "NotImplemented" not in str(excinfo.value)

    def test_a_missing_database_names_the_path_and_does_not_download(
        self, tmp_path,
    ):
        """`allow_database_update` is false, so a missing database is a
        provisioning failure to report, never a cue to fetch."""
        config = _config_with_vfdb(tmp_path, sequences=False)
        with pytest.raises(vf.VirulencePreflightError) as excinfo:
            stage_virulence.run(config, _manifest(tmp_path), RunMode.REAL, tmp_path / "out")
        message = str(excinfo.value)
        assert "vfdb/sequences" in message
        assert "does not download" in message or "not download" in message
        assert "allow_database_update" in message

    def test_test_mode_still_reads_the_fixture(self, tmp_path, config):
        """The REAL branch is additive; TEST is unchanged."""
        target = tmp_path / "virulence"
        target.mkdir(parents=True)
        (target / "virulence_factors.tsv").write_text(
            "sample_id\tvirulence_factor\tgene\tdatabase\tdatabase_version\n"
            "PDT_A\tType IV pili\tfimV\tSYNTHETIC_VFDB\tv0-synthetic\n",
            encoding="utf-8",
        )
        result = stage_virulence.run(config, _manifest(tmp_path), RunMode.TEST, tmp_path)
        assert result["PDT_A"], "the TEST fixture row was not returned"


class TestRealWritesTheTable:
    def test_it_screens_and_writes_one_row_per_factor(self, tmp_path, monkeypatch):
        """End to end with blastn stubbed, so the stage's own writing is covered.

        `screen_isolate` is the seam; blastn itself is exercised against the real
        VFDB in `test_virulence_blast_adapter` and by the smoke run.
        """
        config = _config_with_vfdb(tmp_path)
        monkeypatch.setattr(
            stage_virulence, "ensure_blast_db", lambda *a, **k: tmp_path / "db"
        )
        monkeypatch.setattr(
            stage_virulence, "database_version_on_disk", lambda *a, **k: "sha256:test"
        )
        monkeypatch.setattr(
            stage_virulence,
            "screen_isolate",
            lambda **kwargs: [
                {"sample_id": "", "virulence_factor": "Type IV pili",
                 "gene": "fimV", "category": "Adherence",
                 "database": "VFDB", "database_version": "sha256:test",
                 "confidence": None, "identity_pct": 99.5,
                 "coverage_pct": 100.0, "vfdb_accession": "vfdb~~~fimV~~~X~~~"}
            ],
        )

        out = tmp_path / "out"
        result = stage_virulence.run(
            config, _manifest(tmp_path), RunMode.REAL, out
        )

        assert result["PDT_A"][0].virulence_factor == "Type IV pili"
        written = (out / "virulence" / "virulence_factors.tsv").read_text()
        assert "PDT_A\tType IV pili\tfimV" in written
        # Every declared column is present, or a reader parsing the contract's
        # header would trip over the row.
        header = written.splitlines()[0].split("\t")
        assert header[:8] == list(stage_virulence.VIRULENCE_COLUMNS)

    def test_an_isolate_without_an_assembly_is_screened_with_nothing(
        self, tmp_path, monkeypatch,
    ):
        """Recorded as screened, not dropped.

        Silently omitting it would shrink the denominator for a carriage rate.
        957 of the 967-isolate roster are in exactly this position.
        """
        config = _config_with_vfdb(tmp_path)
        monkeypatch.setattr(
            stage_virulence, "ensure_blast_db", lambda *a, **k: tmp_path / "db"
        )
        monkeypatch.setattr(
            stage_virulence, "database_version_on_disk", lambda *a, **k: "sha256:test"
        )
        calls = []

        def _screen(**kwargs):
            calls.append(kwargs)
            return []

        monkeypatch.setattr(stage_virulence, "screen_isolate", _screen)
        result = stage_virulence.run(
            config, _manifest(tmp_path, with_assembly=False), RunMode.REAL, tmp_path / "out"
        )
        assert calls == [], "no assembly means nothing to screen"
        assert result["PDT_A"] == []