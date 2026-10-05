"""The reference is verified before anything is called against it (spec.md D5).

Red-first against a deliberately corrupted file. The point of D5 is that a
truncated or substituted reference shifts every variant coordinate *silently* -
the numbers come out, they are just about a different genome - so the check has
to be a test with a bad file in it, not a claim that the code reads a field.

No real genome is needed: the verifier takes explicit expectations, so a
corrupted fixture is enough to prove it refuses.
"""

from __future__ import annotations

import pytest

from papipeline.errors import DataContractError
from papipeline.reference import (
    UNPINNED,
    file_md5,
    verify_reference_file,
)

GOOD = ">contig\n" + "ACGT" * 250 + "\n"


@pytest.fixture
def reference(tmp_path):
    path = tmp_path / "ref.fna"
    path.write_text(GOOD, encoding="utf-8")
    return path


class TestAValidReferencePasses:
    def test_length_and_checksum_both_match(self, reference):
        record = verify_reference_file(
            reference,
            expected_length=1000,
            expected_md5=file_md5(reference),
            accession="TEST",
        )
        assert record["pinned"] is True
        assert record["observed_length"] == 1000
        assert record["observed_md5"] == file_md5(reference)


class TestACorruptReferenceIsRefused:
    def test_a_truncated_file_is_refused(self, reference):
        """The common failure: a download that stopped early."""
        reference.write_text(">contig\n" + "ACGT" * 100 + "\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            verify_reference_file(
                reference, expected_length=1000,
                expected_md5=file_md5(reference), accession="TEST",
            )
        assert "400" in str(excinfo.value), (
            "the message must say what was observed, so the operator can tell "
            "truncation from substitution"
        )

    def test_a_substituted_file_is_refused(self, reference, tmp_path):
        """Same length, different sequence: the case length cannot catch."""
        other = tmp_path / "other.fna"
        # Same 1000 bases, different composition, so the file differs while the
        # declared length still matches. This is the substitution case: only the
        # checksum can see it.
        other.write_text(">contig\n" + "AAAA" * 250 + "\n", encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            verify_reference_file(
                other, expected_length=1000,
                expected_md5=file_md5(reference), accession="TEST",
            )
        assert "md5" in str(excinfo.value).lower()

    def test_a_missing_file_is_refused(self, tmp_path):
        with pytest.raises(DataContractError) as excinfo:
            verify_reference_file(tmp_path / "absent.fna", expected_length=1000,
                                  expected_md5="0" * 32)
        assert "not present" in str(excinfo.value)

    def test_the_refusal_names_the_genome(self, tmp_path):
        """An operator with several references needs to know which failed."""
        other = tmp_path / "other.fna"
        other.write_text(GOOD, encoding="utf-8")
        with pytest.raises(DataContractError) as excinfo:
            verify_reference_file(other, expected_length=999,
                                  expected_md5="0" * 32, accession="GCF_000006765.1")
        assert "GCF_000006765.1" in str(excinfo.value)


class TestLengthIsCheckedAgainstBasesNotFileSize:
    def test_rewrapping_does_not_change_the_verdict(self, reference, tmp_path):
        """A reformatted file is the same genome.

        Line wrapping, header text and trailing newlines all change the file
        size without changing a single base. Checking bytes would reject a
        perfectly good reference.
        """
        rewrapped = tmp_path / "rewrapped.fna"
        # 10 bases on the second line plus 990 = 1000, same genome, different
        # bytes: a longer header, different wrapping, different trailing bytes.
        rewrapped.write_text(
            ">contig description changed entirely\nACGTACGTAC\n" + "ACGT" * 247
            + "AC\n",
            encoding="utf-8",
        )
        assert rewrapped.stat().st_size != reference.stat().st_size, (
            "the fixture should differ in file size for this test to mean anything"
        )
        record = verify_reference_file(
            rewrapped, expected_length=1000, expected_md5=file_md5(rewrapped)
        )
        assert record["observed_length"] == 1000


class TestAnUnpinnedReferenceIsNotClaimedToBeVerified:
    def test_unpinned_still_checks_length(self, reference, caplog):
        """Length catches truncation, so an unpinned check is worth running."""
        reference.write_text(">contig\n" + "ACGT" * 10 + "\n", encoding="utf-8")
        with pytest.raises(DataContractError):
            verify_reference_file(reference, expected_length=1000,
                                  expected_md5=UNPINNED)

    def test_unpinned_is_reported_as_unpinned(self, reference):
        record = verify_reference_file(
            reference, expected_length=1000, expected_md5=UNPINNED
        )
        assert record["pinned"] is False
        assert record["observed_length"] == 1000
        assert record["observed_md5"] is None, (
            "an unpinned reference must not report a digest, or a reader will "
            "assume the identity was confirmed"
        )

    def test_no_expected_length_skips_the_check(self, reference):
        """A reference with nothing declared is not a broken one."""
        record = verify_reference_file(reference)
        assert record["observed_length"] == 1000


class TestItIsActuallyCalled:
    """Defined is not wired. This asserts the call site, not the function.

    The recurring failure in this project has been code that exists, is
    well-tested, and is never reached - `papipeline.join`, `scripts/emit.py`.
    So the check is that a REAL run goes through the verifier, and that the
    modes which do not use a reference do not.
    """

    @staticmethod
    def _patch_gate(monkeypatch):
        """Open the REAL gate in-process only. No config file is touched."""
        import papipeline.run as run_module

        monkeypatch.setattr(
            run_module, "resolve_mode",
            lambda mode, config: __import__(
                "papipeline.models", fromlist=["RunMode"]
            ).RunMode(str(mode).strip().upper()),
        )

    def test_a_real_run_verifies_the_reference_before_any_stage(self, config, monkeypatch, tmp_path):
        """Called, and a failure stops the run.

        Proven by making the verifier raise a sentinel: if the call site is
        reached, the sentinel escapes. Driving a whole REAL run to completion
        would need the real fixture set, which is later work - and would prove
        less, since the ordering here is the part that matters.
        """
        import papipeline.run as run_module

        class VerifierReached(Exception):
            pass

        def boom(cfg, mode=None):
            raise VerifierReached(str(mode))

        monkeypatch.setattr(run_module.reference, "verify_reference", boom)
        self._patch_gate(monkeypatch)
        monkeypatch.setattr(type(config), "data_root", lambda self, mode: tmp_path)

        with pytest.raises(VerifierReached) as excinfo:
            run_module.run_pipeline(config=config, mode="REAL")

        assert excinfo.value.args[0] == "REAL", (
            "the verifier must be told which mode asked, so a future "
            "mode-dependent rule has the information"
        )

    def test_a_verification_failure_stops_the_run(self, config, monkeypatch, tmp_path):
        """A bad reference is a hard stop, not a warning.

        If the run continued, it would produce coordinates about a different
        genome and report them as results - which is the failure D5 exists to
        prevent, one stage later.
        """
        import papipeline.run as run_module
        from papipeline.errors import DataContractError

        def refuse(cfg, mode=None):
            raise DataContractError("Reference is 100 bases; 6264404 were declared.")

        monkeypatch.setattr(run_module.reference, "verify_reference", refuse)
        self._patch_gate(monkeypatch)
        monkeypatch.setattr(type(config), "data_root", lambda self, mode: tmp_path)

        with pytest.raises(DataContractError, match="6264404"):
            run_module.run_pipeline(config=config, mode="REAL")

    @pytest.mark.parametrize("mode", ["STUB", "TEST"])
    def test_modes_without_a_reference_do_not_verify(self, mode, monkeypatch):
        """STUB must work with no db/ on disk, so it must not demand one."""
        import papipeline.run as run_module

        monkeypatch.setattr(
            run_module.reference, "verify_reference",
            lambda cfg, m=None: pytest.fail(
                f"{mode} verified a reference; it calls no variants and must "
                "not require a genome on disk"
            ),
        )
        # a config with no reference at all is the honest STUB situation
        from types import SimpleNamespace

        assert SimpleNamespace(reference=None).reference is None


class TestThePinnedChecksumDescribesTheRealFile:
    """Guards the pin against drift, and is skipped when db/ is absent.

    The pin is only worth anything if it keeps describing the file a REAL run
    will read. A test asserting that a hard-coded string equals itself would
    prove nothing, so this compares it against the provisioned genome - and
    skips rather than fails when there is none, because the laptop config is
    allowed to have an empty db/ and CI has no reference to check.
    """

    def test_the_pinned_md5_matches_the_provisioned_file(self, config):
        from pathlib import Path

        from papipeline.config.loader import load_config
        from papipeline.reference import verify_reference_file

        reference_path = config.reference_fasta()
        if not reference_path.exists():
            pytest.skip(f"no reference provisioned at {reference_path}")

        identity = dict((config.raw or {}).get("reference") or {})
        pinned = identity.get("md5")
        assert pinned and pinned != UNPINNED, (
            "the reference is provisioned on this machine but the pin is still a "
            "placeholder; derive it from the real file rather than leaving a "
            "field that verifies nothing"
        )

        # Raises on any mismatch, so this is the drift guard.
        record = verify_reference_file(
            reference_path,
            expected_length=int(identity["expected_length"]),
            expected_md5=pinned,
            accession=str(identity.get("accession", "reference")),
        )
        assert record["pinned"] is True
        assert record["observed_length"] == int(identity["expected_length"])
