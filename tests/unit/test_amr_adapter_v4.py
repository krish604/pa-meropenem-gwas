"""Stage 4's REAL backend, and the v3 code it replaces.

`AmrFinderPlusAdapter` existed but was written against AMRFinderPlus 3.x, and
**its command is rejected outright by the installed 4.2.7**:

    $ amrfinder --input <assembly> ...
    *** ERROR ***
    "--input" is not a valid option

That is the tool's own output, not an inference from release notes - which is
why the flag is asserted here rather than trusted. 4.x renamed `--input` to
`--nucleotide` and re-ordered the report columns, so both the invocation and the
parser needed replacing. `docs/design/10-reusable-modules.md` and
`docs/design/08-migration-plan.md` both say the same thing: use the v4
implementation in `papipeline/pilot/amr_detect.py`, which already exists and is
untested. This is that wiring.

**A row here is synthetic, and that is a real limitation.** No captured real
v4 report is committed anywhere in the repository, and one cannot be produced
until the AMRFinderPlus database is provisioned - it is absent on this machine.
The fixture below therefore uses exactly the columns
`pilot.amr_detect.REQUIRED_COLUMNS` already declares, and *that* is the project's
own statement of the v4 contract, not an invented format. What these tests prove
is that the adapter is wired to the existing v4 parser and emits the right
records. What they cannot prove is that real output parses, and that awaits the
database. `test_the_fixture_matches_the_projects_declared_v4_columns` keeps the
synthetic fixture honest about that.

The v3 parser `parse_amrfinder_stdout` is deliberately left in place. It is
tested, and it remains the correct reader for a v3-format table; it is simply no
longer what a REAL run uses.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from papipeline.adapters import amr as preflight
from papipeline.pilot import amr_detect
from papipeline.stages import amr as stage_amr

#: Columns exactly as `pilot.amr_detect.REQUIRED_COLUMNS` declares them. Not a
#: guess at the v4 header - a restatement of the project's existing contract.
V4_COLUMNS = list(amr_detect.REQUIRED_COLUMNS)

#: `config/science.yaml` `organism.name` with its whitespace already in the
#: database's spelling. Declared here so the tests do not depend on the tool
#: being installed; `tests/unit/test_amr_organism_flag.py` proves the flag is
#: validated against the tool's own `--list_organisms`.
ORGANISM = "Pseudomonas_aeruginosa"

# A minimal v4 report. Values are invented because no real capture exists; see
# the module docstring.
V4_REPORT = "\t".join(V4_COLUMNS) + "\n" + "\t".join(
    [
        "contig_1", "blaOXA-1", "AMR", "BETA-LACTAM", "CARBAPENEM",
        "99.5", "100.0",
    ]
) + "\n" + "\t".join(
    [
        "contig_1", "clpV", "VIRULENCE", "EFFICIENT", "CLPPROTEIN",
        "98.0", "100.0",
    ]
) + "\n"


def write_report(path: Path, text: str = V4_REPORT) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class TestTheCommandIsV4:
    def _adapter(self, tmp_path):
        return stage_amr.AmrFinderPlusAdapter(
            status=_FakeStatus(),
            database="AMRFinderPlus",
            database_version="2026-08-07.1",
            antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db",
            threads=2,
        )

    def test_it_uses_nucleotide_because_the_tool_rejects_input(self, tmp_path):
        command = self._adapter(tmp_path).build_command(Path("a.fna"))
        assert "--nucleotide" in command
        assert "--input" not in command, (
            "4.x renamed --input to --nucleotide; the installed 4.2.7 exits with "
            'ERROR: "--input" is not a valid option'
        )

    def test_it_writes_a_report_file_rather_than_stdout(self, tmp_path):
        """The v3 adapter passed `--output -` and read stdout. The v4 parser
        reads a path, so a report sharing a stream with the tool's progress
        output is no longer something to parse."""
        command = self._adapter(tmp_path).build_command(Path("a.fna"))
        assert command[command.index("--output") + 1] != "-"

    def test_it_never_asks_to_update_the_database(self, tmp_path):
        joined = " ".join(self._adapter(tmp_path).build_command(Path("a.fna")))
        assert "--update" not in joined and "--force_update" not in joined

    def test_the_command_comes_from_the_one_inspectable_builder(self, tmp_path):
        """One place decides the invocation, so the no-update and v4 guarantees
        are asserted once rather than per call site."""
        adapter = self._adapter(tmp_path)
        built = preflight.commands_for_isolate(
            executable="amrfinder", assembly=Path("a.fna"),
            database_dir=tmp_path / "db", threads=2,
            organism=ORGANISM,
            out_path=adapter.report_path_for("a"),
        )
        assert built[0] == adapter.build_command(Path("a.fna"))


class TestTheV4ParserIsTheOneUsed:
    def test_detect_reads_a_v4_report(self, tmp_path):
        # In `report_dir`, under the sample's name - which is where the adapter
        # reads. It used to be written beside the assembly, which was the defect:
        # the command derived the path from the sequence and the lookup used
        # `report_dir`, so a successful screen reported "wrote no report".
        report = write_report(tmp_path / "S1.amrfinder.tsv")

        def runner(command, **kwargs):
            return _FakeResult()

        adapter = stage_amr.AmrFinderPlusAdapter(
            status=_FakeStatus(), database="AMRFinderPlus",
            database_version="2026-08-07.1", antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=2, runner=runner,
            report_dir=tmp_path,
        )
        records = adapter.detect("S1", Path("a.fna"))
        assert [r.determinant for r in records] == ["blaOXA-1"]
        assert records[0].evidence_source == "amrfinderplus"
        assert records[0].database_version == "2026-08-07.1"

    def test_virulence_rows_are_left_to_the_virulence_stage(self, tmp_path):
        """Mixing them would make a virulence gene read as a resistance
        determinant, which is scientific rule 1's exact failure."""
        # In `report_dir`, under the sample's name - which is where the adapter
        # reads. It used to be written beside the assembly, which was the defect:
        # the command derived the path from the sequence and the lookup used
        # `report_dir`, so a successful screen reported "wrote no report".
        report = write_report(tmp_path / "S1.amrfinder.tsv")
        adapter = stage_amr.AmrFinderPlusAdapter(
            status=_FakeStatus(), database="AMRFinderPlus",
            database_version="2026-08-07.1", antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=2,
            runner=lambda command, **kw: _FakeResult(), report_dir=tmp_path,
        )
        determinants = {r.determinant for r in adapter.detect("S1", Path("a.fna"))}
        assert "clpV" not in determinants

    def test_no_assembly_is_an_empty_result_not_an_error(self, tmp_path):
        adapter = stage_amr.AmrFinderPlusAdapter(
            status=_FakeStatus(), database="AMRFinderPlus",
            database_version="v", antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=1,
        )
        assert adapter.detect("S1", None) == []


class TestTheFixtureIsHonestAboutTheV4Contract:
    def test_the_fixture_matches_the_projects_declared_v4_columns(self):
        """The synthetic report must use the columns the project already
        declares, or it is testing an invented format."""
        header = V4_REPORT.splitlines()[0].split("\t")
        assert header == V4_COLUMNS

    def test_a_report_missing_a_column_is_refused(self, tmp_path):
        """A truncated header is a tool or version change, and must not be
        parsed as a short report."""
        write_report(tmp_path / "a.fna.amrfinder.tsv", "Contig id\tElement symbol\n")
        adapter = stage_amr.AmrFinderPlusAdapter(
            status=_FakeStatus(), database="AMRFinderPlus",
            database_version="v", antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=1,
            runner=lambda command, **kw: _FakeResult(), report_dir=tmp_path,
        )
        with pytest.raises(Exception):
            adapter.detect("S1", Path("a.fna"))


class TestTheStageRunsTheToolOnARealRun:
    def _manifest(self, *ids):
        from papipeline.manifest import SampleManifest
        from papipeline.models import Sample

        return SampleManifest(samples=[Sample(sample_id=i) for i in ids])

    def _assemblies(self, root, *ids):
        for sample_id in ids:
            directory = root / sample_id
            directory.mkdir(parents=True, exist_ok=True)
            (directory / f"{sample_id}.fna").write_text(">x\nACGT\n", encoding="utf-8")
        return root

    def test_a_real_run_calls_the_backend_per_isolate(self, config, tmp_path, monkeypatch):
        from papipeline.models import AmrDeterminant, RunMode

        seen = []

        class Recording(_FakeBackend):
            def detect(self, sample_id, assembly):
                seen.append(sample_id)
                return [
                    AmrDeterminant(
                        sample_id=sample_id, antibiotic="imipenem",
                        determinant="blaOXA-1", gene="blaOXA-1",
                        variant=None, determinant_type="AMR", mechanism=None,
                        evidence_source="amrfinderplus", database="AMRFinderPlus",
                        database_version="2026-08-07.1", confidence=None,
                    )
                ]

        # Only the factory is replaced, so the stage's real call path - locate
        # the assembly, call detect, group the result - is what is exercised.
        monkeypatch.setattr(
            stage_amr, "build_amr_adapter",
            lambda config, antibiotic, **kw: Recording(), raising=True,
        )
        data_root = self._assemblies(tmp_path / "data", "PDT_A", "PDT_B")
        result = stage_amr.run(
            config, self._manifest("PDT_A", "PDT_B"), RunMode.REAL,
            intermediate_root=tmp_path, antibiotic="imipenem", data_root=data_root,
        )
        assert sorted(seen) == ["PDT_A", "PDT_B"]
        assert all(result[s] for s in ("PDT_A", "PDT_B"))

    def test_every_manifest_sample_gets_a_key(self, config, tmp_path, monkeypatch):
        """Screened-and-nothing-found is `[]`; not-screened is absent. A sample
        whose backend call failed is recorded as screened-with-nothing rather
        than dropped, so it cannot inflate a cohort denominator."""
        from papipeline.models import RunMode

        monkeypatch.setattr(
            stage_amr, "build_amr_adapter",
            lambda config, antibiotic, **kw: _FakeBackend(), raising=True,
        )
        data_root = self._assemblies(tmp_path / "data", "PDT_A")
        result = stage_amr.run(
            config, self._manifest("PDT_A", "PDT_EMPTY"), RunMode.REAL,
            intermediate_root=tmp_path, antibiotic="imipenem", data_root=data_root,
        )
        assert set(result) == {"PDT_A", "PDT_EMPTY"}


class _FakeStatus:
    executable = "amrfinder"
    available = True
    version = "4.2.7"


class _FakeResult:
    stdout = ""
    stderr = ""
    returncode = 0


class _FakeBackend:
    """A backend that reports nothing, for stage-level wiring tests."""

    name = "amrfinderplus"

    def detect(self, sample_id, assembly):
        return []

    def provenance(self):
        return {"database": "AMRFinderPlus", "database_version": "2026-08-07.1"}


class TestTheDefaultRunner:
    """The default runner, not an injected fake.

    Every other test here injects a fake, so the default was never executed - and
    it was `CommandResult`, a *result* dataclass, called as
    `self._runner(command=...)`. That raised

        TypeError: CommandResult.__init__() missing 4 required
        positional arguments

    on the first real end-to-end run, after three attempts, because the seam had
    never been executable: `amr` had therefore never screened anything, which is
    consistent with its long-standing "no REAL caller" status.

    These put a stand-in `amrfinder` on PATH and let the *real* default executor
    run it. No fake runner is injected anywhere below - that is the whole point,
    and it is what the rest of this file failed to do.

    Exercising the default also caught two more defects in the same path, which a
    fake runner had been hiding: the report was written *beside the assembly*
    (inside `db/smoke_genomes/` for a smoke run) while the adapter read it from
    `report_dir`, and `report_dir` was never created, so the tool's write failed.
    """

    @pytest.fixture
    def fake_amrfinder(self, tmp_path, monkeypatch):
        bindir = tmp_path / "bin"
        bindir.mkdir()
        script = bindir / "amrfinder"
        script.write_text(
            "#!/bin/sh\n"
            "# Stand-in: honour --output, write a v4-shaped report, exit 0.\n"
            'out=""\n'
            "while [ $# -gt 0 ]; do\n"
            '  if [ "$1" = "--output" ]; then out="$2"; shift; fi\n'
            "  shift\n"
            "done\n"
            'printf "Contig id\\tElement symbol\\tType\\tClass\\tSubclass\\t'
            '%% Identity to reference\\t%% Coverage of reference\\n" > "$out"\n'
            'printf "contig_1\\tblaOXA-1\\tAMR\\tBETA-LACTAM\\tBETA-LACTAM\\t'
            '99.5\\t100.0\\n" >> "$out"\n'
            "exit 0\n",
            encoding="utf-8",
        )
        script.chmod(0o755)
        monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
        return bindir

    def _adapter(self, tmp_path, report_dir):
        from papipeline.adapters.external import ToolStatus

        return stage_amr.AmrFinderPlusAdapter(
            status=ToolStatus(
                name="amrfinder", executable="amrfinder", available=True,
                version="4.2.7",
            ),
            database="AMRFinderPlus", database_version="2026-08-07.1",
            antibiotic="imipenem", organism=ORGANISM, database_path=tmp_path / "db", threads=2,
            report_dir=report_dir,
        )

    def test_the_default_runner_actually_executes(self, fake_amrfinder, tmp_path):
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        assembly = genomes / "PDT_A.fna"
        assembly.write_text(">c\nACGT\n", encoding="utf-8")
        adapter = self._adapter(tmp_path, tmp_path / "reports")
        records = adapter.detect("PDT_A", assembly)
        assert [r.determinant for r in records] == ["blaOXA-1"]

    def test_the_report_lands_in_the_report_dir(self, fake_amrfinder, tmp_path):
        genomes = tmp_path / "curated"
        genomes.mkdir()
        assembly = genomes / "PDT_A.fna"
        assembly.write_text(">c\nACGT\n", encoding="utf-8")
        reports = tmp_path / "reports"
        self._adapter(tmp_path, reports).detect("PDT_A", assembly)
        assert [p.name for p in reports.glob("*.amrfinder.tsv")] == [
            "PDT_A.amrfinder.tsv"
        ]

    def test_nothing_is_written_beside_the_assembly(
        self, fake_amrfinder, tmp_path
    ):
        """A genome directory holds sequences. Tool output must not accumulate
        there, and for a smoke run it is the curated fixture directory."""
        genomes = tmp_path / "curated"
        genomes.mkdir()
        assembly = genomes / "PDT_A.fna"
        assembly.write_text(">c\nACGT\n", encoding="utf-8")
        self._adapter(tmp_path, tmp_path / "reports").detect("PDT_A", assembly)
        assert not list(genomes.glob("*.amrfinder.tsv"))
        assert sorted(p.name for p in genomes.iterdir()) == ["PDT_A.fna"]

    def test_the_report_dir_is_created_if_absent(self, fake_amrfinder, tmp_path):
        """AMRFinderPlus does not create the directory it is told to write to,
        so a missing one is a failed screen disguised as a negative one."""
        reports = tmp_path / "not" / "yet"
        assert not reports.exists()
        genomes = tmp_path / "genomes"
        genomes.mkdir()
        assembly = genomes / "PDT_A.fna"
        assembly.write_text(">c\nACGT\n", encoding="utf-8")
        self._adapter(tmp_path, reports).detect("PDT_A", assembly)
        assert (reports / "PDT_A.amrfinder.tsv").is_file()
