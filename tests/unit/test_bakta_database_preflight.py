"""The Bakta database must be present, pinned, and the pinned edition.

`papipeline.stages.annotation` already refuses a REAL run with no Bakta database,
but it refuses on *absence only*. It cannot notice that the database on disk is a
different release from the one `references.tsv` records, and it has no concept of
the light/full edition at all. So the failure modes docs/design/07 records - a
database that is present, usable, and wrong - pass straight through.

Three states, refused, matching `papipeline.adapters.amr.preflight_database`
exactly so there is one preflight to understand rather than two:

* **absent** - the one case the stage already caught;
* **unpinned** - `references.tsv`'s own header calls this "a reported failure,
  not a pass";
* **mismatched, including the edition** - the dangerous one. A full database
  where a light one was pinned, or a light one where a full one was pinned, is
  a different annotation from the one the science was reasoned about, and
  nothing downstream would say so.

**The edition is part of the version string, not a separate field.** The
reference schema is `reference_id, tool, tool_version, database,
database_version, version_status` and has nowhere to put light-vs-full. Rather
than widen a contract that `docs/data_contract.md` owns, the version is
normalised to the shape the repository's own `benchmark/bakta10.py` produces -
`"{major}.{minor} ({type}, {date})"` - which carries the edition as part of the
version. The alternative, a column that exists nowhere else, would be the second
place the edition is stated.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from papipeline.adapters import bakta_db as preflight

REFERENCES = Path("config/references.tsv")
BAKTA_REFERENCE_ID = "ref_bakta"


def write_version_json(directory: Path, **overrides) -> Path:
    """A database directory carrying a realistic `version.json`."""
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "date": "2025-02-24",
        "major": 6,
        "minor": 0,
        "type": "light",
        "doi": "10.5281/zenodo.14916843",
    }
    payload.update(overrides)
    (directory / "version.json").write_text(json.dumps(payload), encoding="utf-8")
    return directory


@pytest.fixture
def references(tmp_path):
    def _write(database_version="6.0 (light, 2025-02-24)", version_status="pinned"):
        path = tmp_path / "references.tsv"
        path.write_text(
            "reference_id\ttool\ttool_version\tdatabase\tdatabase_version\t"
            "version_status\tsource\tnotes\n"
            f"{BAKTA_REFERENCE_ID}\tbakta\t1.12.1\tBaktaDB\t{database_version}\t"
            f"{version_status}\tzenodo\ttest\n",
            encoding="utf-8",
        )
        return path

    return _write


class TestTheEditionIsMachineReadable:
    """`version.json` records the edition, and it is what makes the pin checkable."""

    def test_the_light_edition_is_recovered(self, tmp_path):
        directory = write_version_json(tmp_path / "db-light")
        assert preflight.database_version(directory) == "6.0 (light, 2025-02-24)"

    def test_the_full_edition_is_distinguished_from_light(self, tmp_path):
        full = write_version_json(tmp_path / "db-full", type="full")
        light = write_version_json(tmp_path / "db-light", type="light")
        assert preflight.database_version(full) != preflight.database_version(light)

    def test_a_missing_version_json_is_not_a_version(self, tmp_path):
        """Bakta refuses a database without one ("version file not readable"),
        so a directory lacking it is incomplete rather than unversioned."""
        empty = tmp_path / "db-light"
        empty.mkdir()
        assert preflight.database_version(empty) is None

    def test_it_matches_the_format_the_benchmark_already_produces(self):
        """One format, not two.

        `benchmark/bakta10.py` builds this string already, so the preflight
        reuses its shape rather than inventing a second spelling that would need
        translating between the two call sites.
        """
        directory = write_version_json(Path("/tmp") / "unused")
        assert preflight.database_version(directory).startswith("6.0 (")


class TestTheThreeRefusedStates:
    def test_absent_is_refused(self, tmp_path, references):
        with pytest.raises(Exception) as excinfo:
            preflight.preflight_database(
                database_dir=tmp_path / "absent", references_path=references()
            )
        assert "absent" in str(excinfo.value)

    def test_the_refusal_names_how_to_provision_it(self, tmp_path, references):
        with pytest.raises(Exception) as excinfo:
            preflight.preflight_database(
                database_dir=tmp_path / "absent", references_path=references()
            )
        assert "bakta_db download" in str(excinfo.value)

    def test_unpinned_is_refused(self, tmp_path, references):
        directory = write_version_json(tmp_path / "db-light")
        with pytest.raises(Exception) as excinfo:
            preflight.preflight_database(
                database_dir=directory,
                references_path=references("UNPINNED", "unpinned"),
            )
        message = str(excinfo.value)
        assert "unpinned" in message.lower()
        assert "references.tsv" in message

    def test_a_different_release_is_refused(self, tmp_path, references):
        """Present, usable, and wrong - the case the stage's own check misses."""
        directory = write_version_json(tmp_path / "db-light", date="2025-09-01")
        with pytest.raises(Exception) as excinfo:
            preflight.preflight_database(
                database_dir=directory, references_path=references()
            )
        message = str(excinfo.value)
        assert "2025-09-01" in message and "2025-02-24" in message

    def test_a_different_edition_is_refused(self, tmp_path, references):
        """A full database where a light one was pinned.

        This is the one that matters most for the science: a species-specific
        light database and a full one annotate differently, so a run that
        silently got the other edition would report genes the pinned science
        never reasoned about.
        """
        directory = write_version_json(tmp_path / "db-full", type="full")
        with pytest.raises(Exception) as excinfo:
            preflight.preflight_database(
                database_dir=directory, references_path=references()
            )
        message = str(excinfo.value)
        assert "full" in message and "light" in message

    def test_an_agreement_passes(self, tmp_path, references):
        directory = write_version_json(tmp_path / "db-light")
        assert preflight.preflight_database(
            database_dir=directory, references_path=references()
        ) == "6.0 (light, 2025-02-24)"


class TestTheCheckedInContract:
    """`ref_bakta` is pinned, and the pin carries the edition."""

    def test_the_pin_is_real(self):
        pin = preflight.read_reference_pin(REFERENCES, BAKTA_REFERENCE_ID)
        assert pin.is_pinned, (
            f"ref_bakta records database_version={pin.database_version!r} "
            f"status={pin.version_status!r}"
        )

    def test_the_pin_names_the_edition(self):
        """Because the schema has no edition column, the version string is the
        only place it can live. If it is missing, the pin cannot be checked."""
        pin = preflight.read_reference_pin(REFERENCES, BAKTA_REFERENCE_ID)
        assert "light" in pin.database_version

    def test_the_configured_path_agrees_with_the_pin(self, config):
        configured = str(config.raw["annotation"]["bakta_db"])
        pin = preflight.read_reference_pin(REFERENCES, BAKTA_REFERENCE_ID)
        assert "light" in Path(configured).name, (
            f"annotation.bakta_db is {configured!r} but the pin says "
            f"{pin.database_version!r}; a light pin pointed at a full database "
            "would fail the preflight at run time rather than here"
        )

    @pytest.mark.skipif(
        not Path("db/bakta_db/db-light/version.json").is_file(),
        reason="database not provisioned here; db/ is gitignored",
    )
    def test_the_provisioned_database_is_the_pinned_one(self):
        assert preflight.preflight_database(
            database_dir=Path("db/bakta_db/db-light"),
            references_path=REFERENCES,
        ) == "6.0 (light, 2025-02-24)"
