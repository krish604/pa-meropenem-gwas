"""Stage 2 preflight: the Bakta database must be present, pinned, and light.

`papipeline.stages.annotation` already refuses a REAL run with no Bakta database,
but it refuses on *absence only*. It cannot notice that the database on disk is a
different release from the one `config/references.tsv` records, and it has no
concept of the light/full edition. So the failure mode
`docs/design/07-bakta-optimization.md` records as the lesson to encode - tool and
database disagreeing - passes straight through, and the run proceeds against a
database the contract does not describe.

The three refused states match
:func:`papipeline.adapters.amr.preflight_database` exactly, so there is one
preflight to understand rather than two:

* **absent** - the single case the stage already caught;
* **unpinned** - `references.tsv`'s own header calls this "a reported failure,
  not a pass";
* **mismatched, edition included** - present, usable, and wrong.

**The edition is part of the version string.** The reference schema has no
column for light-vs-full, so rather than widen a contract `docs/data_contract.md`
owns, the version is normalised to the shape
`papipeline/benchmark/bakta10.py` already produces - `"{major}.{minor}
({type}, {date})"` - which carries the edition inside the version. A new column
would be a second place for the edition to be stated and a second to be wrong.

This reads the database's own `version.json` and never provisions, downloads or
updates anything, so `allow_database_update: false` keeps holding.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from ..adapters.amr import ReferencePin, read_reference_pin
from ..errors import PipelineError
from ..logging_utils import get_logger

LOGGER = get_logger("adapters.bakta_db")

#: The reference row that governs stage 2.
BAKTA_REFERENCE_ID = "ref_bakta"

#: The file Bakta itself reads to decide the database is usable, and refuses
#: without ("version file not readable").
VERSION_FILE = "version.json"

#: Named in refusals so the reader has somewhere to go. Never invoked from here.
PROVISIONING_COMMAND = "bakta_db download --type light --output <dir>"


def database_version(database_dir: Path) -> Optional[str]:
    """The version a Bakta database reports about itself.

    Normalised to ``"{major}.{minor} ({type}, {date})"`` - the same string
    `papipeline.benchmark.bakta10` builds, so the two call sites agree without
    translating between them. The `type` is the edition, and it is inside the
    version because the reference schema has no column for it.

    Returns ``None`` when there is no `version.json`. That is not an unknown
    version but an incomplete database: Bakta itself refuses such a directory,
    so it is treated the same way.
    """
    path = Path(database_dir) / VERSION_FILE
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PipelineError(
            f"{path} could not be read as JSON, so the Bakta database's version "
            f"cannot be established: {exc}",
            database_dir=str(database_dir),
            version_file=str(path),
        ) from exc

    major = payload.get("major")
    minor = payload.get("minor")
    edition = payload.get("type")
    date = payload.get("date")
    if major is None or minor is None or not edition or not date:
        raise PipelineError(
            f"{path} does not record the fields a version is built from "
            f"(major, minor, type, date); got major={major!r}, minor={minor!r}, "
            f"type={edition!r}, date={date!r}. Refusing rather than assembling a "
            "partial version, which would not match the pin for a reason nobody "
            "could later reconstruct.",
            database_dir=str(database_dir),
            version_file=str(path),
        )
    return f"{major}.{minor} ({edition}, {date})"


def preflight_database(
    *,
    database_dir: Path,
    references_path: Path,
    reference_id: str = BAKTA_REFERENCE_ID,
) -> str:
    """Verify the database exists and matches the pin. Returns the version.

    Raises:
        PipelineError: absent, unpinned, or mismatched - the last including the
            light/full edition, since a full database where a light one was
            pinned is a different annotation from the one the science was
            reasoned about.
    """
    database_dir = Path(database_dir)

    if not database_dir.is_dir():
        raise PipelineError(
            f"the Bakta database directory {database_dir} does not exist, so "
            "stage 2 cannot run. Provision it out of band with "
            f"`{PROVISIONING_COMMAND}`; this pipeline never downloads a database "
            "during a run, because an implicit mid-run update would split the "
            "cohort across two databases with nothing to say so.",
            database_dir=str(database_dir),
            provisioning_command=PROVISIONING_COMMAND,
        )

    pin: ReferencePin = read_reference_pin(references_path, reference_id)
    if not pin.is_pinned:
        raise PipelineError(
            f"{reference_id} in {references_path} is unpinned "
            f"(database_version={pin.database_version!r}, "
            f"version_status={pin.version_status!r}), and that file's own header "
            "calls an unpinned reference 'a reported failure, not a pass'. Set "
            f"the real database_version and version_status=pinned for "
            f"{reference_id} before a REAL run.",
            references_path=str(references_path),
            reference_id=reference_id,
            database_version=pin.database_version,
            version_status=pin.version_status,
        )

    on_disk = database_version(database_dir)
    if on_disk is None:
        raise PipelineError(
            f"{database_dir} has no {VERSION_FILE}, so its version - and "
            "whether it is the light or full edition - cannot be checked "
            f"against the {pin.database_version} that {references_path} records. "
            "Bakta itself refuses such a directory. An unverifiable database is "
            "treated as a mismatched one: the point of the check is that every "
            "genome is annotated against a database the run can name.",
            database_dir=str(database_dir),
            expected=pin.database_version,
        )

    if on_disk != pin.database_version:
        raise PipelineError(
            f"the Bakta database on disk is {on_disk!r} but {references_path} "
            f"records {pin.database_version!r} for {reference_id}. Refusing "
            "rather than running: annotating part of the cohort against one "
            "database and part against another is unreproducible, and if the "
            "difference is the edition then the annotation itself changes - a "
            "full database reports genes a species-specific light one does not. "
            "Either point the run at the pinned database or update the pin "
            "deliberately.",
            database_dir=str(database_dir),
            on_disk=on_disk,
            pinned=pin.database_version,
            reference_id=reference_id,
        )

    LOGGER.info(
        "Stage 2 preflight: Bakta database %s matches %s", on_disk, reference_id
    )
    return on_disk
