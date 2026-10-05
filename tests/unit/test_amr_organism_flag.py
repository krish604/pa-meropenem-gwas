"""Stage 4's `--organism`, which was absent and cost the whole POINT screen.

**What was wrong.** `papipeline.adapters.amr.commands_for_isolate` built the
AMRFinderPlus command with no `--organism`. The docstring that explained the
omission claimed "nothing in this pipeline records a species to put in it" -
which was false. `config/science.yaml` has carried

    organism:
      name: "Pseudomonas aeruginosa"
      taxid: 287

for the life of the project, and `.scratch/imipenem-gwas-dashboard/spec.md:343`
names the invocation as ``AMRFinderPlus (--organism Pseudomonas_aeruginosa)``.
The code simply never read the key.

**What it cost, measured on the 10-isolate smoke subset** (AMRFinderPlus 4.2.7,
database 2026-08-07.1, `--plus` on both sides; reports under
`pa-artifacts/round11/amr/{without,with}/`):

* without the flag, **every** row came back `Subtype=AMR` and the `variant`
  column was empty for the whole cohort - 144 rows, 144 `.`, in the stage-4
  table that already existed;
* with the flag, rows with `Subtype=POINT` appear, e.g. `oprD_V359L`,
  `nalC_G71E`, `nalC_S209R`, `mexR_I24AfsTer94`.

For an imipenem study whose entire mechanistic story is loss of function in
`oprD`, that is the screen that matters most being switched off, with no error
and no warning: AMRFinderPlus treats an unknown organism and no organism the
same way, reporting nothing rather than complaining.

**The spelling is measured, not assumed.** `ORGANISM` below is validated
against the installed tool's own `--list_organisms` by
`test_the_configured_organism_is_one_the_installed_tool_lists`, and the
`<gene>_<variant>` split is measured against the provisioned database's own
`AMRProt-mutation.tsv` by `test_the_symbol_split_matches_every_database_row`.
Both skip when the tool or database is absent, because this repository does not
ship either; nothing here is checked against a captured file.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pytest

from papipeline.adapters import amr as adapter
from papipeline.adapters.external import CommandResult, ToolStatus
from papipeline.stages import amr as stage_amr

#: What `config/science.yaml` says, and what `--list_organisms` must contain.
CONFIG_ORGANISM_NAME = "Pseudomonas aeruginosa"
ORGANISM = "Pseudomonas_aeruginosa"

#: The configured database, when this machine has one. Untracked by design
#: (AGENTS.md rule 4), so every test that needs it skips rather than fails.
#: Resolved against the repository root rather than the working directory, so a
#: test run from anywhere finds the same database - or skips for the same reason.
CONFIGURED_DB = (
    Path(__file__).resolve().parents[2] / "db" / "amrfinderplus-db" / "2026-08-07.1"
)


def _list_organisms_stdout(choices):
    body = ", ".join(choices)
    return (
        "Running: amrfinder --database <db> --list_organisms\n"
        "Software version: 4.2.7\n"
        f"Available --organism options: {body}\n"
        "amrfinder took 0 seconds to complete\n"
    )


def _listing(*choices, returncode=0, stdout=None, stderr=""):
    return CommandResult(
        command=["amrfinder", "--list_organisms"],
        returncode=returncode,
        stdout=_list_organisms_stdout(choices) if stdout is None else stdout,
        stderr=stderr,
        duration_s=0.0,
    )


def _real_tool():
    import shutil

    return shutil.which("amrfinder")


class TestTheOrganismIsDerivedNotGuessed:
    def test_whitespace_becomes_the_databases_underscore_spelling(self):
        assert adapter.derive_organism_flag(CONFIG_ORGANISM_NAME) == ORGANISM

    def test_an_already_spelled_name_passes_through_unchanged(self):
        assert adapter.derive_organism_flag(ORGANISM) == ORGANISM

    def test_surrounding_and_repeated_whitespace_collapses(self):
        assert adapter.derive_organism_flag("  Pseudomonas   aeruginosa ") == ORGANISM

    def test_an_absent_name_refuses_and_names_the_key(self):
        """A missing key must not become a default organism. Defaulting would
        be a second, unvalidated statement of the species - and the species is
        what decides whether point mutations are screened at all."""
        with pytest.raises(adapter.PipelineError) as excinfo:
            adapter.derive_organism_flag(None)
        assert "organism.name" in str(excinfo.value)
        assert "POINT" in str(excinfo.value), (
            "the refusal has to say what is lost, or it reads as pedantry"
        )

    def test_a_blank_name_refuses_too(self):
        with pytest.raises(adapter.PipelineError):
            adapter.derive_organism_flag("   ")


class TestTheOrganismListComesFromTheTool:
    def test_the_list_is_parsed_out_of_the_tools_own_line(self, tmp_path):
        seen = {}

        def runner(command, **kwargs):
            seen["command"] = list(command)
            return _listing(ORGANISM, "Escherichia")

        assert adapter.list_organisms(
            "amrfinder", tmp_path / "db", runner=runner
        ) == (ORGANISM, "Escherichia")
        assert seen["command"] == [
            "amrfinder", "--database", str(tmp_path / "db"), "--list_organisms",
        ], "the database matters: the list is per-database, not per-tool"

    def test_a_failing_list_is_a_refusal_not_an_empty_list(self, tmp_path):
        """An empty list would make every organism look unknown, which is a
        confusing way to report "the database is not there"."""
        def runner(command, **kwargs):
            return _listing(returncode=1, stdout="", stderr="No valid database")

        with pytest.raises(adapter.PipelineError) as excinfo:
            adapter.list_organisms("amrfinder", tmp_path / "db", runner=runner)
        assert "No valid database" in str(excinfo.value)

    def test_output_without_the_marker_is_a_refusal(self, tmp_path):
        """A tool that stops printing the line is a version change, and guessing
        the names from memory is exactly the drift this avoids."""
        def runner(command, **kwargs):
            return _listing(stdout="Software version: 4.2.7\n")

        with pytest.raises(adapter.PipelineError) as excinfo:
            adapter.list_organisms("amrfinder", tmp_path / "db", runner=runner)
        assert "Available --organism options" in str(excinfo.value)


class TestAnUnknownOrganismRefusesAndNamesIt:
    def test_it_names_the_value_it_tried_and_the_name_it_came_from(self, tmp_path):
        def runner(command, **kwargs):
            return _listing(ORGANISM, "Escherichia")

        with pytest.raises(adapter.PipelineError) as excinfo:
            adapter.resolve_organism(
                "Pseudomonas aureolata", executable="amrfinder",
                database_dir=tmp_path / "db", runner=runner,
            )
        message = str(excinfo.value)
        assert "Pseudomonas_aureolata" in message
        assert "Pseudomonas aureolata" in message

    def test_it_names_the_organism_in_real_amrfinder_would_reject(self, tmp_path):
        """This is the whole failure this replaces: AMRFinderPlus does not
        complain about an organism it does not know. It just reports no POINT
        rows, so the refusal has to say so."""

        def runner(command, **kwargs):
            return _listing(ORGANISM)

        with pytest.raises(adapter.PipelineError) as excinfo:
            adapter.resolve_organism(
                "Escherichia coli K-12", executable="amrfinder",
                database_dir=tmp_path / "db", runner=runner,
            )
        assert "--list_organisms" in str(excinfo.value)
        assert "POINT" in str(excinfo.value)

    def test_a_known_organism_resolves_to_the_same_string(self, tmp_path):
        def runner(command, **kwargs):
            return _listing(ORGANISM, "Escherichia")

        assert adapter.resolve_organism(
            CONFIG_ORGANISM_NAME, executable="amrfinder",
            database_dir=tmp_path / "db", runner=runner,
        ) == ORGANISM


class TestTheCommandCarriesTheOrganism:
    def _command(self, tmp_path, organism=ORGANISM):
        return adapter.commands_for_isolate(
            executable="amrfinder", assembly=Path("a.fna"),
            database_dir=tmp_path / "db", threads=2, organism=organism,
            out_path=tmp_path / "reports" / "a.amrfinder.tsv",
        )[0]

    def test_it_passes_the_validated_flag(self, tmp_path):
        command = self._command(tmp_path)
        assert command[command.index("--organism") + 1] == ORGANISM

    def test_the_adapter_builds_the_same_command(self, tmp_path):
        stage = stage_amr.AmrFinderPlusAdapter(
            status=ToolStatus(
                name="amrfinder", executable="amrfinder", available=True,
                version="4.2.7",
            ),
            database="AMRFinderPlus", database_version="2026-08-07.1",
            antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=2,
            report_dir=tmp_path / "reports",
        )
        assert stage.build_command(Path("a.fna")) == self._command(tmp_path)

    def test_there_is_no_default_that_could_silently_drop_the_flag(self, tmp_path):
        """`organism` is keyword-only and required. Making it optional would
        reintroduce the exact failure, one layer down."""
        with pytest.raises(TypeError):
            adapter.commands_for_isolate(
                executable="amrfinder", assembly=Path("a.fna"),
                database_dir=tmp_path / "db", threads=2,
            )


class TestProvenanceNamesTheOrganism:
    def _adapter(self, tmp_path):
        return stage_amr.AmrFinderPlusAdapter(
            status=ToolStatus(
                name="amrfinder", executable="amrfinder", available=True,
                version="4.2.7",
            ),
            database="AMRFinderPlus", database_version="2026-08-07.1",
            antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=2,
            report_dir=tmp_path / "reports",
        )

    def test_it_records_organism_tool_version_and_database_version(self, tmp_path):
        record = self._adapter(tmp_path).provenance()
        assert record["organism"] == ORGANISM
        assert record["tool_version"] == "4.2.7"
        assert record["database_version"] == "2026-08-07.1"

    def test_the_table_written_by_a_real_run_carries_it(
        self, config, tmp_path, monkeypatch
    ):
        """The header is on the *file*, because the file is what a later reader
        of this cohort opens first. Without it, a point-mutation-free table and
        a table from a run that never screened for point mutations are
        byte-identical."""
        class _Empty:
            name = "amrfinderplus"

            def detect(self, sample_id, assembly):
                return []

            def provenance(self):
                return {
                    "organism": ORGANISM, "tool": "amrfinderplus",
                    "tool_version": "4.2.7", "database": "AMRFinderPlus",
                    "database_version": "2026-08-07.1",
                }

        from papipeline.manifest import SampleManifest
        from papipeline.models import RunMode, Sample

        monkeypatch.setattr(
            stage_amr, "build_amr_adapter",
            lambda config, antibiotic, **kw: _Empty(), raising=True,
        )
        genomes = tmp_path / "genomes"
        (genomes / "PDT_A").mkdir(parents=True)
        (genomes / "PDT_A" / "PDT_A.fna").write_text(">x\nACGT\n")
        stage_amr.run(
            config, SampleManifest(samples=[Sample(sample_id="PDT_A")]),
            RunMode.REAL, tmp_path, "imipenem", data_root=genomes,
        )
        text = (tmp_path / "amr" / "amr_determinants.tsv").read_text()
        header = [line for line in text.splitlines() if line.startswith("#")]
        joined = "\n".join(header)
        assert f"organism={ORGANISM!r}" in joined
        assert "tool_version='4.2.7'" in joined
        assert "database_version='2026-08-07.1'" in joined

    def test_a_header_only_table_still_carries_the_provenance(self, tmp_path):
        """A screen that found nothing and a screen that never ran must not look
        alike, and the provenance is part of that difference."""
        path = tmp_path / "amr_determinants.tsv"

        class _Backend:
            def provenance(self):
                return {"organism": ORGANISM, "database_version": "2026-08-07.1"}

        stage_amr._write_determinant_table(
            path, {"PDT_A": []}, backend=_Backend(),
        )
        text = path.read_text()
        assert f"organism={ORGANISM!r}" in text
        assert "sample_id" in text.splitlines()[-1]

    def test_the_tool_version_is_read_from_the_tool_not_the_pin(self, tmp_path):
        """`config/references.tsv` records what the cohort is reproducible
        against; this records what ran. `preflight_database` already refuses
        when they disagree, so reading the pin here would be a second, weaker
        check."""

        def runner(command, **kwargs):
            assert command == ["amrfinder", "--version"]
            return CommandResult(
                command=command, returncode=0, stdout="4.2.7\n",
                stderr="", duration_s=0.0,
            )

        assert adapter.tool_version("amrfinder", runner=runner) == "4.2.7"

    def test_an_unreadable_version_is_none_rather_than_a_guess(self, tmp_path):
        def runner(command, **kwargs):
            return CommandResult(
                command=command, returncode=1, stdout="", stderr="boom",
                duration_s=0.0,
            )

        assert adapter.tool_version("amrfinder", runner=runner) is None


class TestThePointMutationIsSplitIntoGeneAndVariant:
    """`oprD_V359L` must not stay a gene name.

    `pilot.amr_detect.parse_amrfinder_report` copies `Element symbol` into both
    `determinant` and `gene` and leaves `variant` empty. Left alone that makes
    the `oprD` carbapenem mechanism fail its lookup in
    `config/mechanisms.tsv`, so switching the flag *on* would lose the very
    signal it recovered.
    """

    HEADER = (
        "Contig id\tElement symbol\tType\tSubtype\tClass\tSubclass\t"
        "% Identity to reference\t% Coverage of reference"
    )

    def _report(self, path, symbol, subtype="POINT"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            self.HEADER + "\n"
            + f"contig_1\t{symbol}\tAMR\t{subtype}\tBETA-LACTAM\tCARBAPENEM\t"
              "99.0\t100.0\n",
            encoding="utf-8",
        )
        return path

    def _records(self, tmp_path, symbol, subtype="POINT"):
        report = self._report(
            tmp_path / "reports" / "S1.amrfinder.tsv", symbol, subtype
        )
        return stage_amr.AmrFinderPlusAdapter(
            status=ToolStatus(
                name="amrfinder", executable="amrfinder", available=True,
                version="4.2.7",
            ),
            database="AMRFinderPlus", database_version="2026-08-07.1",
            antibiotic="imipenem", organism=ORGANISM,
            database_path=tmp_path / "db", threads=1,
            report_dir=tmp_path / "reports",
            runner=lambda command, **kw: CommandResult(
                command=command, returncode=0, stdout="", stderr="", duration_s=0.0
            ),
        ).parse_report(report, "S1")

    def test_oprd_is_the_gene_and_v359l_the_variant(self, tmp_path):
        record, = self._records(tmp_path, "oprD_V359L")
        assert record.gene == "oprD"
        assert record.variant == "V359L"

    def test_the_determinant_keeps_the_tools_own_spelling(self, tmp_path):
        """The tool's string is the evidence; nothing is normalised away from
        it, so a reader can match it against the report."""
        record, = self._records(tmp_path, "oprD_V359L")
        assert record.determinant == "oprD_V359L"

    def test_a_point_disrupt_row_is_split_too(self, tmp_path):
        """`POINT_DISRUPT` is the alignment-detected case - `mexR_I24AfsTer94`
        from the smoke run - and it is not in the database's listed point
        mutations, so a lookup-only implementation would miss it."""
        record, = self._records(tmp_path, "mexR_I24AfsTer94", "POINT_DISRUPT")
        assert (record.gene, record.variant) == ("mexR", "I24AfsTer94")

    def test_an_intact_gene_is_left_alone(self, tmp_path):
        record, = self._records(tmp_path, "mexA", subtype="AMR")
        assert record.gene == "mexA"
        assert record.variant is None

    def test_a_gene_containing_a_hyphen_is_split_correctly(self, tmp_path):
        record, = self._records(tmp_path, "aph(3')-IIb_V103I")
        assert record.gene == "aph(3')-IIb"
        assert record.variant == "V103I"

    def test_an_unsplittable_point_symbol_is_left_unsplit_and_logged(
        self, tmp_path, caplog
    ):
        """A wrong gene/variant split attributes a mechanism to the wrong gene,
        which is worse than an empty variant column. So it refuses to guess."""
        record, = self._records(tmp_path, "noUnderscoreHere")
        assert record.gene == "noUnderscoreHere"
        assert record.variant is None
        assert any("noUnderscoreHere" in r.message for r in caplog.records)


class TestTheSplitIsMeasuredAgainstTheRealDatabase:
    """Not a captured file: the provisioned database, when this machine has one."""

    MUTATION_TABLE = CONFIGURED_DB / "AMRProt-mutation.tsv"

    def _symbols(self):
        rows = []
        with self.MUTATION_TABLE.open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#") or line.startswith("taxgroup"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 8:
                    rows.append(fields)
        return rows

    @pytest.mark.skipif(
        not MUTATION_TABLE.is_file(),
        reason="no provisioned AMRFinderPlus database on this machine",
    )
    def test_the_split_matches_every_mutation_row_in_the_database(self):
        """The `<gene>_<variant>` pattern is measured here rather than assumed.

        Every row's recovered gene is compared against the gene in that row's
        own `mutated_protein_name`, so the check is against a second statement
        of the gene made by the database rather than against itself.
        """
        rows = self._symbols()
        assert rows, f"{self.MUTATION_TABLE} parsed to no rows"
        mismatched = []
        for fields in rows:
            gene, variant = stage_amr.split_point_symbol(fields[3])
            if gene is None or variant is None:
                mismatched.append((fields[3], "did not split"))
                continue
            if gene.lower() != fields[7].rsplit("_", 1)[-1].lower():
                mismatched.append((fields[3], gene))
        assert not mismatched, (
            f"{len(mismatched)} of {len(rows)} mutation symbols do not split into "
            f"the gene the database names: {mismatched[:5]}"
        )

    @pytest.mark.skipif(
        not MUTATION_TABLE.is_file(),
        reason="no provisioned AMRFinderPlus database on this machine",
    )
    def test_no_gene_symbol_in_the_database_contains_an_underscore(self):
        """This is what makes the split unambiguous rather than a guess about
        where the gene ends."""
        genes = {
            stage_amr.split_point_symbol(fields[3])[0]
            for fields in self._symbols()
        }
        assert genes, "no genes parsed"
        assert not [g for g in genes if g and "_" in g]


class TestAgainstTheInstalledTool:
    """These run the real `amrfinder`, so they prove the spelling and not the
    project's belief about it."""

    @pytest.mark.skipif(
        _real_tool() is None or not CONFIGURED_DB.is_dir(),
        reason="amrfinder or the provisioned database is absent on this machine",
    )
    def test_the_configured_organism_is_one_the_installed_tool_lists(self):
        assert adapter.list_organisms(_real_tool(), CONFIGURED_DB) != ()
        assert ORGANISM in adapter.list_organisms(_real_tool(), CONFIGURED_DB)

    @pytest.mark.skipif(
        _real_tool() is None or not CONFIGURED_DB.is_dir(),
        reason="amrfinder or the provisioned database is absent on this machine",
    )
    def test_resolve_returns_the_flag_unchanged_for_the_real_config_value(self):
        assert adapter.resolve_organism(
            CONFIG_ORGANISM_NAME, executable=_real_tool(),
            database_dir=CONFIGURED_DB,
        ) == ORGANISM

    @pytest.mark.skipif(
        _real_tool() is None or not CONFIGURED_DB.is_dir(),
        reason="amrfinder or the provisioned database is absent on this machine",
    )
    def test_the_real_tool_reports_its_own_version(self):
        completed = subprocess.run(
            [_real_tool(), "--version"], capture_output=True, text=True, check=False
        )
        assert completed.returncode == 0
        recorded = adapter.tool_version(_real_tool())
        assert recorded == completed.stdout.strip()


class TestBuildAmrAdapterRefusesBeforeAnyGenomeIsScreened:
    def test_an_unlisted_organism_refuses_at_build_time(self, config, monkeypatch):
        """The refusal has to happen before any isolate is analysed, because an
        invalid `--organism` does not fail the tool - it just reports nothing."""
        from papipeline.adapters import external

        monkeypatch.setattr(
            # **kwargs, not `names` alone: stage 4 now says `execute=True`
            # explicitly, because it is about to run AMRFinderPlus for the
            # cohort. The double accepts the keyword so it still short-circuits
            # detection; `require_tool` is patched below to return the status.
            external, "detect_tools", lambda names, **kwargs: [], raising=True
        )
        monkeypatch.setattr(
            external, "require_tool",
            lambda tools, name, why: ToolStatus(
                name=name, executable="amrfinder", available=True, version="4.2.7",
            ),
            raising=True,
        )
        monkeypatch.setattr(
            adapter, "preflight_database",
            lambda **kw: "2026-08-07.1", raising=True,
        )
        monkeypatch.setattr(
            adapter, "list_organisms",
            lambda executable, database_dir, runner=None: ("Escherichia",),
            raising=True,
        )
        with pytest.raises(adapter.PipelineError) as excinfo:
            stage_amr.build_amr_adapter(config, "imipenem")
        assert "Pseudomonas_aeruginosa" in str(excinfo.value)

    def test_a_listed_organism_reaches_the_adapter(self, config, monkeypatch):
        from papipeline.adapters import external

        monkeypatch.setattr(
            # **kwargs, not `names` alone: stage 4 now says `execute=True`
            # explicitly, because it is about to run AMRFinderPlus for the
            # cohort. The double accepts the keyword so it still short-circuits
            # detection; `require_tool` is patched below to return the status.
            external, "detect_tools", lambda names, **kwargs: [], raising=True
        )
        monkeypatch.setattr(
            external, "require_tool",
            lambda tools, name, why: ToolStatus(
                name=name, executable="amrfinder", available=True, version="4.2.7",
            ),
            raising=True,
        )
        monkeypatch.setattr(
            adapter, "preflight_database",
            lambda **kw: "2026-08-07.1", raising=True,
        )
        monkeypatch.setattr(
            adapter, "list_organisms",
            lambda executable, database_dir, runner=None: (ORGANISM,),
            raising=True,
        )
        built = stage_amr.build_amr_adapter(config, "imipenem")
        assert built.organism == ORGANISM
        assert built.build_command(Path("a.fna"))[5] == "--organism"


def _config():
    from papipeline.config.loader import PipelineConfig

    return PipelineConfig(
        raw={
            "organism": {"name": CONFIG_ORGANISM_NAME, "taxid": 287},
            "amr": {"database": "AMRFinderPlus"},
            "runtime": {"threads": 2},
        }
    )