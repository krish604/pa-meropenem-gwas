"""The tool and its database must agree before stage 4 starts.

`docs/design/07-bakta-optimization.md` records why this is a first-class
preflight and not a runtime surprise: Bakta 1.12.1 ships AMRFinderPlus database
`2024-12-18.1`, AMRFinderPlus 4.2.7 requires `>= 2025-09-22.2`, neither has
another `osx-arm64` build, and Bakta's expert system then fails outright with
"AMR expert system failed". A *partial* database was also auto-created by
Bakta's own internal AMRFinder call and had to be deleted, because an implicit
mid-run database update is exactly what this pipeline forbids.

So the failure has to arrive before the first genome is analysed, naming the
database and the reason. Three states are refused, and each for its own reason:

* **no database** - the tool is installed and cannot run;
* **pin not recorded** - `references.tsv` says `unpinned`, which that file's own
  header calls "a reported failure, not a pass";
* **pin and disk disagree** - the most dangerous of the three, because it looks
  healthy. Half a cohort annotated against one database and half against
  another is unreproducible, and nothing downstream would ever say so.

The floor is read from `references.tsv`, the file that calls itself the
reproducibility contract. It is not hard-coded here, because a hard-coded floor
would be a second statement of the same fact, free to drift from the first.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import amr as adapter
from papipeline.config.loader import PipelineConfig
from papipeline.errors import PipelineError

REFERENCES = Path("config/references.tsv")
AMR_REFERENCE_ID = "ref_amrfinderplus"


@pytest.fixture
def provisioned(tmp_path):
    """A database directory that looks genuinely provisioned."""
    directory = tmp_path / "amrfinderplus-db" / "2026-08-07.1"
    directory.mkdir(parents=True)
    (directory / "version.txt").write_text("2026-08-07.1\n", encoding="utf-8")
    (directory / "AMRProt.fa").write_text(">x\nMKK\n", encoding="utf-8")
    return directory


@pytest.fixture
def references(tmp_path):
    """A references.tsv with a real pin, unless overridden."""
    def _write(database_version="2026-08-07.1", version_status="pinned"):
        path = tmp_path / "references.tsv"
        path.write_text(
            "reference_id\ttool\ttool_version\tdatabase\tdatabase_version\t"
            "version_status\tsource\tnotes\n"
            f"{AMR_REFERENCE_ID}\tamrfinder\t4.2.7\tAMRFinderPlus\t"
            f"{database_version}\t{version_status}\tamrfinder_update\ttest\n",
            encoding="utf-8",
        )
        return path

    return _write


class TestTheDatabaseMustBeThere:
    def test_an_absent_database_is_refused(self, tmp_path, references):
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=tmp_path / "not-provisioned",
                references_path=references(),
            )
        assert "not-provisioned" in str(excinfo.value)

    def test_the_refusal_says_how_to_provision_it(self, tmp_path, references):
        """"Add the key" is easy to act on; a bare "database not found" sends
        the reader looking for a download flag that does not exist.

        `references.tsv` states the database is provisioned out of band and
        never updated during a run, so the message has to say the same thing
        rather than implying the pipeline will fetch it.
        """
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=tmp_path / "absent", references_path=references()
            )
        message = str(excinfo.value)
        assert "amrfinder_update" in message or "out of band" in message

    def test_an_empty_directory_is_refused(self, tmp_path, references):
        """A directory that exists but holds no database is the *partial*
        database docs/07 records, and it is the state that fails confusingly."""
        empty = tmp_path / "amrfinderplus-db" / "2026-08-07.1"
        empty.mkdir(parents=True)
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=empty, references_path=references()
            )
        assert "2026-08-07.1" in str(excinfo.value)


class TestThePinMustBeRecorded:
    def test_an_unpinned_reference_is_refused(self, tmp_path, provisioned, references):
        """`references.tsv` calls unpinned "a reported failure, not a pass"."""
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=provisioned,
                references_path=references("UNPINNED", "unpinned"),
            )
        message = str(excinfo.value)
        assert "unpinned" in message.lower()
        assert AMR_REFERENCE_ID in message

    def test_it_names_the_file_to_edit(self, provisioned, references):
        """The pin lives in one file, and the reader has to be told which."""
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=provisioned,
                references_path=references("UNPINNED", "unpinned"),
            )
        assert "references.tsv" in str(excinfo.value)


class TestThePinAndTheDiskMustAgree:
    def test_a_mismatch_is_refused(self, tmp_path, provisioned, references):
        """The dangerous one: everything looks healthy and the cohort is
        silently split across two databases."""
        other = tmp_path / "amrfinderplus-db" / "2025-01-01.1"
        other.mkdir(parents=True)
        (other / "version.txt").write_text("2025-01-01.1\n", encoding="utf-8")
        (other / "AMRProt.fa").write_text(">x\nMKK\n", encoding="utf-8")
        with pytest.raises(PipelineError) as excinfo:
            adapter.preflight_database(
                database_dir=other, references_path=references()
            )
        message = str(excinfo.value)
        assert "2025-01-01.1" in message and "2026-08-07.1" in message

    def test_an_agreement_passes_and_returns_the_version(
        self, provisioned, references
    ):
        assert adapter.preflight_database(
            database_dir=provisioned, references_path=references()
        ) == "2026-08-07.1"

    def test_the_version_comes_from_the_files_not_a_literal(
        self, provisioned, references
    ):
        """A hard-coded expectation would pass here and be wrong the moment the
        pin is legitimately bumped."""
        (provisioned / "version.txt").write_text("2026-09-30.9\n", encoding="utf-8")
        assert adapter.preflight_database(
            database_dir=provisioned, references_path=references("2026-09-30.9")
        ) == "2026-09-30.9"


class TestTheRealContractIsPinned:
    """`references.tsv` now records a pinned AMRFinderPlus database.

    This class previously asserted the opposite - that the pin was unpinned and
    the preflight therefore refused. That was true when it was written, and it
    said so: "a future commit that pins the reference has to mean it
    deliberately - the test will change." The database has now been provisioned
    out of band and pinned, so the invariant is the forward-looking one: the
    contract records a real version, and the preflight accepts the database
    rather than refusing it.

    Skipped when the database is absent, because `db/` is gitignored
    (AGENTS.md rule 4). A fresh clone has the pin and not the files, which is
    the intended split: the contract is committed, the data is provisioned.
    """

    DATABASE = Path("db/amrfinderplus-db/2026-08-07.1")

    def test_the_pin_is_a_real_version(self):
        pin = adapter.read_reference_pin(REFERENCES, AMR_REFERENCE_ID)
        assert pin.is_pinned, (
            f"references.tsv records database_version={pin.database_version!r} "
            f"status={pin.version_status!r} for {AMR_REFERENCE_ID}. An unpinned "
            "reference is a reported failure, not a pass, per that file's header."
        )

    def test_the_tool_version_is_recorded_too(self):
        pin = adapter.read_reference_pin(REFERENCES, AMR_REFERENCE_ID)
        assert pin.tool_version not in ("", "UNPINNED"), (
            "the tool version is half the compatibility question - a database "
            "that satisfies one tool may not satisfy another"
        )

    def test_the_configured_path_agrees_with_the_pin(self, config: PipelineConfig):
        """The path and the pin are two statements of the same thing.

        A path pointing at a different release than the pin records would pass
        a file-exists check and fail the preflight at run time instead, which is
        the ordering that hides the mistake.
        """
        configured = str(config.raw["amr"]["database_path"])
        pin = adapter.read_reference_pin(REFERENCES, AMR_REFERENCE_ID)
        assert pin.database_version in Path(configured).name, (
            f"amr.database_path is {configured!r} but the pin is "
            f"{pin.database_version!r}; the two must name the same release"
        )

    def test_the_preflight_accepts_the_provisioned_database(self, config: PipelineConfig):
        if not adapter.looks_provisioned(self.DATABASE):
            pytest.skip("database not provisioned here; db/ is gitignored")
        assert adapter.preflight_database(
            database_dir=self.DATABASE, references_path=REFERENCES
        ) == "2026-08-07.1"

    def test_the_database_really_is_the_pinned_release(self):
        """Read the release's own record, not the pin.

        This is the check that would catch a database directory being replaced
        in place, or an `amrfinder -u` having been run against it - both of
        which leave a path that resolves and a version that disagrees.
        """
        if not self.DATABASE.is_dir():
            pytest.skip("database not provisioned here; db/ is gitignored")
        recorded = adapter.database_version_on_disk(self.DATABASE)
        pin = adapter.read_reference_pin(REFERENCES, AMR_REFERENCE_ID)
        assert recorded == pin.database_version


class TestNothingEverUpdatesTheDatabase:
    def test_no_command_the_adapter_builds_can_update_the_database(
        self, provisioned
    ):
        """`allow_database_update: false` everywhere, and docs/07 records an
        implicit update as the thing that had to be cleaned up by hand."""
        for command in adapter.commands_for_isolate(
            executable="/usr/bin/amrfinder",
            assembly=Path("a.fna"),
            database_dir=provisioned,
            threads=4,
            organism="Pseudomonas_aeruginosa",
        ):
            joined = " ".join(command)
            assert "--update" not in joined
            assert "--force_update" not in joined
