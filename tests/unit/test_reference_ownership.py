"""The reference genome and the tool versions belong to the science.

Ticket 03, corrected after review. The first version of the machine overlay
carried a `reference:` block and a per-tool `database:` path, which quietly
re-opened the exact hole the science/machine split exists to close: two machines
could load *different reference genomes* and *different AMR database releases*
without error, because the accessor of that information was per-machine.

The rule now:

* **Which reference** - accession, assembly, organism, taxid, expected length -
  is a scientific fact. It lives in science.yaml, and an overlay may not
  restate it.
* **Which database release** is governed by config/references.tsv, the
  existing reproducibility contract. It is not duplicated in an overlay, and
  every version in it is currently `UNPINNED`, which is a *reported failure*,
  not a pass.
* **Where things live** and **whether a tool exists on this machine** are the
  only per-machine facts. A tool entry may say `available` and `path`; it may
  not carry a version or a database release.
* No `latest`. A moving symlink is the opposite of a pin, and an analysis that
  reads it twice can read two different databases.

Written before the fix.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest
import yaml

from papipeline.config.loader import (
    MACHINE_OWNED_KEYS,
    load_config,
    load_machine_config,
)
from papipeline.errors import ConfigError

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SCIENCE = PIPELINE_ROOT / "config" / "science.yaml"
LAPTOP = PIPELINE_ROOT / "config" / "machines" / "laptop.yaml"
BIGMACHINE = PIPELINE_ROOT / "config" / "machines" / "bigmachine.yaml"
REFERENCES_TSV = PIPELINE_ROOT / "config" / "references.tsv"


def _overlay(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "machines" / "rogue.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(LAPTOP, path)
    path.write_text(path.read_text(encoding="utf-8") + text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# the reference identity is science
# ---------------------------------------------------------------------------


def test_science_owns_the_reference_identity():
    section = (load_config(SCIENCE, machine="laptop").raw.get("reference") or {})
    assert section, "science.yaml must declare the reference genome"
    assert section.get("accession") == "GCF_000006765.1"
    assert section.get("organism") == "Pseudomonas aeruginosa"
    assert section.get("taxid") == 208964


def test_the_reference_length_is_recorded_so_a_truncated_file_is_caught():
    section = load_config(SCIENCE, machine="laptop").raw["reference"]
    assert section.get("expected_length") == 6264404


def test_no_checksum_is_invented():
    """A checksum may only be present if it can be re-derived from the real file.

    This was originally unconditional - md5 had to stay UNPINNED, because a
    checksum cannot be derived without the file and an invented one verifies
    nothing. That guard did its job: it is why the pin below was derived rather
    than guessed.

    The file now exists, so the checksum is legitimate, and the protection is
    preserved rather than removed: a present checksum must match the digest of
    the provisioned genome. `test_reference_verification.py` proves the
    derivation independently; this keeps the ownership property - the value must
    be UNPINNED *or* real - so a future edit that types a digest from memory
    fails here even with no reference on the machine.
    """
    section = load_config(SCIENCE, machine="laptop").raw["reference"]
    value = section.get("md5")
    if value is None or str(value).upper() in {"UNPINNED", ""}:
        return  # legitimately unpinned, and the strongest state

    from papipeline.config.loader import load_config as _load
    from papipeline.reference import file_md5, verify_reference_file

    config = _load(SCIENCE, machine="laptop")
    path = config.reference_fasta()
    if not path.exists():
        # No file to check against, so the value cannot be shown to be derived.
        # Refuse rather than wave it through: an uncheckable checksum is the
        # thing this test exists to prevent.
        pytest.fail(
            f"md5 is {value!r} but there is no reference at {path} to derive it "
            "from. Either provision the genome, or set md5 back to UNPINNED."
        )

    observed = file_md5(path)
    assert observed == str(value).lower(), (
        f"md5 {value!r} does not match the provisioned genome ({observed}). A "
        "checksum that does not describe the file verifies nothing."
    )
    # And it must be internally consistent with the declared length, or the two
    # describe different genomes.
    verify_reference_file(
        path,
        expected_length=int(section["expected_length"]),
        expected_md5=str(value),
        accession=str(section.get("accession", "reference")),
    )


@pytest.mark.parametrize("overlay", [LAPTOP, BIGMACHINE], ids=["laptop", "bigmachine"])
def test_no_overlay_declares_the_reference_identity(overlay: Path):
    raw = yaml.safe_load(overlay.read_text(encoding="utf-8"))
    assert "reference" not in raw, (
        "which reference genome the science uses is a scientific fact, not a "
        "per-machine one"
    )
    assert "reference" not in MACHINE_OWNED_KEYS


def test_an_overlay_cannot_point_at_a_different_reference(tmp_path: Path):
    """The exploit the review demonstrated: two machines, two genomes, no error."""
    rogue = _overlay(tmp_path, '\nreference:\n  accession: "GCF_999999999.9"\n')
    with pytest.raises(ConfigError):
        load_machine_config(rogue)


# ---------------------------------------------------------------------------
# tool versions and database releases are governed by references.tsv
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("overlay", [LAPTOP, BIGMACHINE], ids=["laptop", "bigmachine"])
def test_an_overlay_tool_entry_carries_no_database_release(overlay: Path):
    raw = yaml.safe_load(overlay.read_text(encoding="utf-8"))
    for tool, entry in (raw.get("tools") or {}).items():
        assert set(entry) <= {"available", "path"}, (
            f"{overlay.name}: tool {tool!r} may declare only availability and a "
            f"path; found {sorted(set(entry) - {'available', 'path'})}. A "
            f"database release here would be per-machine, which is the hole "
            f"this split closes."
        )


def test_no_overlay_points_at_a_moving_latest_symlink():
    for overlay in (LAPTOP, BIGMACHINE):
        text = overlay.read_text(encoding="utf-8")
        assert "/latest" not in text, (
            f"{overlay.name} points at a moving 'latest'; a run must name the "
            f"pinned release from config/references.tsv"
        )


def test_the_version_contract_still_lives_in_references_tsv():
    """Single source of truth: the overlay must not become a second one."""
    rows = [
        line
        for line in REFERENCES_TSV.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    header = rows[0].split("\t")
    assert "database_version" in header
    assert "version_status" in header
    assert any("ref_reference_genome" in r for r in rows)


def test_an_unpinned_reference_is_reported_not_treated_as_pinned():
    """UNPINNED is a reported failure, so a run cannot quietly use whatever.

    `ref_amrfinderplus` used to be asserted here, as a second reference known to
    be unpinned. It is now pinned, deliberately, and a test that kept demanding
    it be unpinned would fail the moment the project did the right thing. The
    intent is that *unpinned references are reported*, so it is asserted against
    the set as a whole: several remain, and the list must not be empty.

    `ref_reference_genome` was asserted here for the same reason and was removed
    the day the reference was actually pinned: GCF_000006765.1 ASM676v1, with
    the GFF verified against the FASTA's own sequence-region and .fai length.
    That is the second time this test has had to drop a name because the project
    fixed something, which is the intended direction of failure - the assertion
    is a statement about the table, not a fixture.
    """
    from papipeline.config.loader import load_references

    refs = load_references(REFERENCES_TSV)
    unpinned = [r.reference_id for r in refs.values() if not r.is_pinned]
    assert len(unpinned) >= 2, (
        f"only {unpinned} remain unpinned. If the project has genuinely pinned "
        "everything, delete this assertion rather than keeping a reference in "
        "the table as a placeholder - an UNPINNED row is a reported failure, "
        "and keeping one just to satisfy a test makes it decorative."
    )


def test_a_pinned_reference_reports_its_real_version():
    """The counterpart: a pin that resolves is a pin, not a placeholder.

    `ref_amrfinderplus` is pinned to a real release, so it must now read as
    pinned. Without this, a pin that silently failed to register would satisfy
    the unpinned test above just as well as a real one.
    """
    from papipeline.config.loader import load_references

    refs = load_references(REFERENCES_TSV)
    amr = refs["ref_amrfinderplus"]
    assert amr.is_pinned
    assert amr.database_version == "2026-08-07.1"
    assert amr.tool_version not in ("", "UNPINNED")


def test_an_overlay_declaring_a_tool_version_is_rejected(tmp_path: Path):
    rogue = _overlay(tmp_path, "\n")
    text = rogue.read_text(encoding="utf-8")
    rogue.write_text(
        text.replace(
            "  amrfinder:   { available: true",
            "  amrfinder:   { available: true, version: '4.2.7'",
        ),
        encoding="utf-8",
    )
    with pytest.raises(ConfigError) as excinfo:
        load_machine_config(rogue)
    assert "version" in str(excinfo.value).lower()


# ---------------------------------------------------------------------------
# availability stays a machine fact
# ---------------------------------------------------------------------------


def test_availability_is_read_from_the_overlay_not_a_global_list(tmp_path: Path):
    """Availability is a per-machine fact, and this proves it mechanically.

    **This test used to lean on a coincidence in the production configs.** It
    asserted that `gubbins` was `false` on laptop and `true` on bigmachine, so
    the two overlays demonstrably differed. That example is gone: panaroo became
    available on both (pip, from source) and gubbins then became available on
    both (conda recipe links no zlib on osx-arm64, so the laptop builds it from
    source instead). The two `tools` blocks are now identical, and no tool
    differs.

    Asserting that they differ anyway would be worse than losing the assertion.
    It would couple a test about *mechanism* to the accident of which tools this
    project currently needs, so the next time a tool becomes genuinely available
    on both machines the suite would go red for a reason that has nothing to do
    with the property being protected. That is precisely how "availability is one
    global list" would slip back in unnoticed.

    So build two purpose-made overlays that differ, and require the loader to
    report the difference. The property under test is that `tool_available`
    consults the overlay it was given; the fixtures supply the difference rather
    than waiting for the configs to provide one.
    """
    false_side = _overlay(tmp_path, "\n")
    false_side.write_text(
        false_side.read_text(encoding="utf-8").replace(
            "  bakta:       { available: true }",
            "  bakta:       { available: false }",
        ),
        encoding="utf-8",
    )
    as_false = load_machine_config(false_side)
    as_shipped = load_machine_config(LAPTOP)

    assert as_shipped.tool_available("bakta") is True
    assert as_false.tool_available("bakta") is False, (
        "tool_available returned the same answer for two overlays that declare "
        "it differently - availability is being read from somewhere other than "
        "the overlay"
    )
    # And an undeclared tool stays unknown in both, rather than inheriting a
    # neighbour's value.
    assert as_shipped.tool_available("a_tool_nobody_declared") is None
    assert as_false.tool_available("a_tool_nobody_declared") is None


def test_both_overlays_now_declare_the_same_tools_available(tmp_path: Path):
    """Records why the test above no longer leans on a production difference.

    Asserted so that the next person to add a laptop-only tool knows this file
    deliberately stopped depending on one, and so that "the overlays agree" is a
    recorded fact rather than an accident waiting to be read as a bug.
    """
    laptop = yaml.safe_load(LAPTOP.read_text(encoding="utf-8"))
    big = yaml.safe_load(BIGMACHINE.read_text(encoding="utf-8"))
    left = {k: (v or {}).get("available") for k, v in (laptop.get("tools") or {}).items()}
    right = {k: (v or {}).get("available") for k, v in (big.get("tools") or {}).items()}
    differing = {k for k in set(left) | set(right) if left.get(k) != right.get(k)}
    assert not differing, (
        "the overlays now differ on "
        f"{sorted(differing)}; if that is intended, this file's per-overlay "
        "mechanism test still holds on purpose-made overlays, but the "
        "'availability must differ' framing no longer describes reality"
    )
    # panaroo and gubbins are available on both, for unrelated documented reasons.
    for tool in ("panaroo", "gubbins"):
        assert left.get(tool) is True and right.get(tool) is True, tool


def test_panaroo_is_available_on_both_overlays_and_that_is_recorded():
    """Both overlays now declare panaroo available, deliberately.

    Asserted rather than left implicit so that re-diverging the two overlays
    looks like a change someone made, not something that quietly happened.
    """
    laptop = load_machine_config(LAPTOP)
    big = load_machine_config(BIGMACHINE)
    assert laptop.tool_available("panaroo") is True
    assert big.tool_available("panaroo") is True


def test_tool_paths_may_not_be_absolute():
    """An absolute path is a fact about one machine wearing a machine's clothes."""
    for overlay in (LAPTOP, BIGMACHINE):
        raw = yaml.safe_load(overlay.read_text(encoding="utf-8"))
        for value in (raw.get("paths") or {}).values():
            if isinstance(value, str):
                assert not value.startswith("/"), f"{overlay.name}: {value!r} is absolute"


# ---------------------------------------------------------------------------
# State-agnostic invariants
# ---------------------------------------------------------------------------
#
# The assertion removed above named `ref_reference_genome` and had to be edited
# out the day the reference was actually pinned - the second time this file has
# had to drop a name because the project fixed something. A test that must be
# edited every time reality improves is a maintenance liability wearing the
# costume of a safety check.
#
# These assert what must hold *whatever* the current state is. They fail when the
# table becomes incoherent, not when it changes.


def test_pin_state_and_declared_status_never_disagree():
    """`is_pinned` is derived from `version_status`, so they cannot disagree.

    If that ever becomes two independent fields, this is the test that catches
    them drifting - and it holds whether the project has pinned everything or
    nothing.
    """
    from papipeline.config.loader import load_references

    for ref in load_references(REFERENCES_TSV).values():
        assert ref.is_pinned == (ref.version_status == "pinned"), (
            f"{ref.reference_id}: is_pinned={ref.is_pinned} but "
            f"version_status={ref.version_status!r}"
        )


def test_a_pinned_reference_states_a_version_and_an_unpinned_one_does_not():
    """The reporting distinction is only meaningful if the two differ.

    A row that claims `pinned` while carrying `UNPINNED` in its version would be
    reported as pinned and used as one - the exact outcome `is_pinned` exists to
    prevent, reachable by typo rather than by decision.
    """
    from papipeline.config.loader import load_references

    for ref in load_references(REFERENCES_TSV).values():
        if ref.is_pinned:
            assert ref.database_version and "UNPINNED" not in ref.database_version, (
                f"{ref.reference_id} claims pinned but its version is "
                f"{ref.database_version!r}"
            )
        else:
            # Not the converse: several unpinned rows carry no version at all
            # rather than the literal UNPINNED, and demanding one would be
            # asserting a convention the table does not follow. The direction
            # that matters is checked above - a row must never *claim* pinned
            # while carrying an absent or UNPINNED version.
            if ref.database_version:
                assert "UNPINNED" in ref.database_version, (
                    f"{ref.reference_id} is unpinned but reports version "
                    f"{ref.database_version!r}; an unpinned row carrying a "
                    "concrete version reads as pinned"
                )


def test_every_reference_declares_a_source():
    """Independent of pinned state.

    A reference with no source cannot be re-fetched, so its pin is a claim with
    nothing behind it - which is how a version silently stops matching the bytes
    on disk.
    """
    from papipeline.config.loader import load_references

    missing = [
        r.reference_id
        for r in load_references(REFERENCES_TSV).values()
        if not (r.source or "").strip()
    ]
    assert not missing, f"references with no declared source: {missing}"


def test_the_table_is_not_decorative():
    """At least one reference is pinned.

    Asserted about the *table* rather than any single row: an all-`UNPINNED`
    table is a placeholder, and this catches that without demanding any
    particular reference be in either state.
    """
    from papipeline.config.loader import load_references

    refs = load_references(REFERENCES_TSV)
    pinned = [r.reference_id for r in refs.values() if r.is_pinned]
    assert pinned, (
        "no reference is pinned. Either the project has pinned nothing yet, or "
        "the table has become a placeholder - and an UNPINNED row kept only to "
        "satisfy a test is decorative."
    )
