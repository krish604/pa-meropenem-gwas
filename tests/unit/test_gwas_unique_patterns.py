"""Tests for step 12a: unique variant patterns and the testing threshold.

What is being pinned here, in order of how badly it would hurt to get wrong:

1. The vendored helper is byte-identical to pyseer's own
   ``scripts/count_patterns.py``. A helper that has been quietly edited is a
   multiple-testing threshold from somewhere unknown, reported as if it were
   pyseer's.
2. The count comes from that helper, and it agrees with an independent tally
   of the same file. If the two disagree, the run stops rather than picking
   one.
3. The threshold is ``alpha / n_unique_patterns`` with ``alpha`` read from
   ``config/science.yaml`` - not a literal in this code.
4. The reduction is *inspectable*: the survivors and the audit record are on
   disk, and the survivors file reconciles against ``sort -u | wc -l``.

The inputs are synthetic pyseer output, in pyseer 1.1.2's own formats:

* ``patterns.txt`` - what ``--output-patterns`` writes. pyseer hashes each
  tested variant's binary vector and appends the base64 of the digest, one per
  line. The encoding is **not** the obvious one - it spans the platform's
  default integer width, not one byte per sample - so it is produced by
  **calling** ``pyseer.input.hash_pattern``. See
  ``TestSyntheticOutputMatchesTheInstalledPyseer``.
* ``pyseer_raw.tsv`` - the association table, with pyseer's own column names
  from ``pyseer/__main__.py``: ``variant af filter-pvalue lrt-pvalue beta
  beta-std-err variant_h2 notes`` for an ``--lmm`` run. The column **count** is
  verified against pyseer's ``format_output``; the **names** cannot be, because
  ``format_output`` emits values and not names. A wrong name fails loud - the
  parser requires ``variant`` and ``lrt-pvalue`` and refuses by name otherwise
  - but it is not verified here.

These formats were read off the installed package, not from memory; see
``docs/environment-arm64.md`` section 4. Where pyseer's own code can produce
them, these fixtures **call** it rather than reimplementing it - a first draft
reimplemented `hash_pattern` and got the encoding wrong. The tests do not
*invoke* pyseer, because spec.md D8 defines TEST mode as running no real tools,
and they skip where the package is absent rather than requiring it. The
consequence - that no test in this repository observes the real pyseer binary's
stdout or stderr - is recorded as still-open on ticket 16 rather than papered
over.
"""

from __future__ import annotations

import base64
import hashlib
import inspect
import re
import shutil
import subprocess
import sys

from pathlib import Path
from types import SimpleNamespace
from typing import Dict

import numpy as np
import pytest

from papipeline.config import load_config
from papipeline.config.loader import load_machine_config
from papipeline.errors import ConfigError, DataContractError
from papipeline.stages import gwas as stage


# --------------------------------------------------------------------------
# Synthetic pyseer output
# --------------------------------------------------------------------------

#: The digest of pyseer's `scripts/count_patterns.py` as retrieved from
#: https://github.com/mgalardini/pyseer/blob/master/scripts/count_patterns.py on
#: 2026-10-02, whole file including its shebang. Recorded in the vendored file's
#: header for a human to check the vendoring once; CI cannot verify it offline,
#: and does not pretend to.
UPSTREAM_WHOLE_FILE_SHA256 = (
    "1081bd1c1a2409cc115facd639af0d876ba0b01ab10ccafd1cabdf926e817643"
)

#: sha256 of the upstream *logic*: from the copyright line to the end, trailing
#: newline normalised. This is the digest CI enforces, and the vendored file's
#: header records the same value, so the file and the test cannot disagree.
VENDORED_BODY_SHA256 = "2d1898a029095c82a1ce6aa095b4fa190a1e4d646291fac1aecfae4e9a34b1a7"

#: Where the upstream logic starts. The digest is taken from here so that
#: documenting the vendoring - which this repository's header does, at length -
#: does not read as modifying it.
UPSTREAM_BODY_MARKER = "# Copyright 2017 Marco Galardini and John Lees"


def _require_pyseer():
    """pyseer's library, or skip with the reason stated.

    Two of this file's inputs cannot be constructed without it: pyseer's
    pattern encoding and its output format. Every fixture here is built from
    one of those two, so **most** of this file needs the package and skips
    without it - the reduction arithmetic and the two-pass wiring included.

    That is a deliberate trade. spec.md D7 has STUB and TEST exercise the whole
    DAG regardless of tool availability, and a test suite that hard-fails on a
    machine without pyseer is the same mistake one layer up. The survivors are
    whatever never reaches a fixture: the helper-output parser, the machine
    configuration checks, the vendored-helper provenance and flag tests, and the
    refusals that fire before a fixture is built. Measured, that is 30 of 75.

    The price is real and should not be glossed: on such a machine those 45 are
    **unverified rather than wrong**. The fixture builders delegate to pyseer,
    so they cannot fabricate a format - only fail - and a skip says which claim
    went unchecked.

    pyseer 1.1.2 is pinned in both `environment/*.yml`, so any environment built
    from this repository's own specification has it. The skip is for machines
    that did not.
    """
    pytest.importorskip(
        "pyseer",
        reason=(
            "pyseer's own pattern encoding and output format cannot be "
            "checked without the package; the formats are then unverified "
            "rather than wrong"
        ),
    )


def pyseer_pattern(vector) -> str:
    """One line of a ``--output-patterns`` file, as pyseer writes it.

    Calls **pyseer's own** ``hash_pattern`` rather than reimplementing it,
    because the encoding is not what it looks like. `pyseer.input.read_variant`
    accumulates a ``sample -> 0/1`` dict and hands ``hash_pattern`` an
    ``np.array`` of Python ints, so numpy gives it the platform's default
    *integer* width and the ``k.view(np.uint8)`` inside spans **8 bytes per
    sample** on a 64-bit machine, not one byte each. Hashing a packed
    ``b"\\x01\\x01\\x00"`` gives a different digest, and every fixture built that
    way would have described a pattern pyseer could never emit.

    A first draft of this file did exactly that, and it was wrong - which is
    why `test_the_encoding_is_the_one_pyseer_builds_for_an_rtab` exercises
    pyseer's reader rather than trusting this function. Delegating also means
    that if pyseer's encoding changes, the fixtures follow it instead of
    quietly diverging.

    This is an *import*, not an invocation: TEST mode is defined as running no
    real tools (spec.md D8), so nothing here executes ``pyseer``.
    """
    _require_pyseer()
    from pyseer.input import hash_pattern

    column = np.asarray([int(value) for value in vector])
    return hash_pattern(column).decode().strip()


def pyseer_lmm_table(rows) -> str:
    """A synthetic pyseer association table, tab-separated, header first.

    Built through **pyseer's own** ``format_output`` so the column count and
    order come from the installed package rather than from a hand-written
    string that could drift from it.

    Args:
        rows: One tuple per variant, in the column order below.

    Returns:
        A tab-separated table whose header is ``PYSEER_LMM_HEADER``.
    """
    _require_pyseer()
    from pyseer.classes import LMM
    from pyseer.utils import format_output

    lines = [PYSEER_LMM_HEADER]
    for variant, af, filter_p, lrt_p, beta, stderr, h2 in rows:
        item = LMM(
            kmer=variant,
            pattern=b"synthetic",
            af=af,
            prep=filter_p,
            pvalue=lrt_p,
            kbeta=beta,
            bse=stderr,
            frac_h2=h2,
            max_lineage=None,
            kstrains=np.array([]),
            nkstrains=np.array([]),
            notes=[],
            prefilter=False,
            filter=False,
        )
        lines.append(format_output(item, lmm=True))
    return "\n".join(lines) + "\n"


#: Eight sample columns, so a few vectors distinguish cleanly.
SAMPLE_COLUMNS = 8

VECTORS = (
    # A perfect positive: every one of the four R samples carries it.
    (1, 1, 1, 1, 0, 0, 0, 0),
    # Byte-identical to the above - carries no independent information, and
    # pyseer fits it separately. Counting both as distinct patterns would halve
    # the divisor - `alpha / 8` instead of `alpha / 4` - for no reason.
    (1, 1, 1, 1, 0, 0, 0, 0),
    # Perfect negative.
    (0, 0, 0, 0, 1, 1, 1, 1),
    # Present in half the cohort, unrelated to the outcome.
    (1, 0, 1, 0, 1, 0, 1, 0),
    (1, 0, 1, 0, 1, 0, 1, 0),
    # Present in a single sample: below min_samples_per_group, so it is not
    # testable, but pyseer still writes its pattern if it survives min-af.
    (0, 0, 0, 0, 0, 0, 0, 1),
)

EXPECTED_UNIQUE = 4
EXPECTED_TESTED = len(VECTORS)
EXPECTED_REDUNDANT = EXPECTED_TESTED - EXPECTED_UNIQUE


#: pyseer 1.1.2's own association-table header for an `--lmm` run with no
#: lineage, no covariates and no sample lists, read off `pyseer/__main__.py`
#: (the `header` list, plus the `variant_h2` and `notes` additions).
#:
#: `test_the_synthetic_association_table_has_pyseers_own_column_count` checks
#: this against the installed `format_output`, so a mismatch cannot pass.
PYSEER_LMM_HEADER = (
    "variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\tvariant_h2\tnotes"
)


def write_patterns(path: Path, vectors) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(pyseer_pattern(v) + "\n" for v in vectors))
    return path


def parse_helper(stdout: str):
    """The vendored helper's own report, read by the production parser.

    Deliberately the production parser, not a reimplementation: these tests
    assert the helper's output is readable by the code that reads it, and a
    second parser here would pass while the real one broke.
    """
    fields: Dict[str, str] = {}
    for line in stdout.splitlines():
        if line.strip():
            key, sep, value = line.partition(":")
            assert sep, f"unexpected helper output line: {line!r}"
            fields[key.strip()] = value.strip()
    return SimpleNamespace(
        patterns=int(fields["Patterns"]),
        threshold=float(fields["Threshold"]),
    )


def wc_lines(command: str) -> int:
    """``wc -l``'s answer, parsed.

    BSD ``wc`` pads its output to a fixed width, so ``len(stdout.split())``
    counts the padding as tokens and silently returns 1.
    """
    return int(
        subprocess.run(
            command, shell=True, capture_output=True, text=True, check=True
        ).stdout.strip()
    )


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


@pytest.fixture
def overlay(tmp_path: Path, pipeline_root: Path) -> Path:
    """A machine overlay that is the laptop's, with a writable scratch dir.

    Copied rather than referenced so ``unique_patterns_temp_dir`` points at a
    directory that exists and belongs to this test. Editing the committed
    overlay from a test would leave the working tree dirty, and two tests in
    this repository already assert that it stays clean.
    """
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    text = (pipeline_root / "config" / "machines" / "laptop.yaml").read_text()
    text = text.replace(
        "unique_patterns_temp_dir: /tmp",
        f"unique_patterns_temp_dir: {scratch}",
    )
    path = tmp_path / "overlay.yaml"
    path.write_text(text)
    return path


@pytest.fixture
def machine_config(overlay: Path, pipeline_root: Path):
    return load_config(pipeline_root / "config" / "science.yaml", machine=str(overlay))


@pytest.fixture
def patterns_file(tmp_path: Path) -> Path:
    return write_patterns(tmp_path / "patterns.txt", VECTORS)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    return tmp_path / "work"


# --------------------------------------------------------------------------
# The vendored helper
# --------------------------------------------------------------------------


class TestSyntheticOutputMatchesTheInstalledPyseer:
    """The fixtures are built from pyseer, not from a memory of pyseer.

    Everything else in this file is only as trustworthy as the notion of a
    "pyseer patterns file" and "pyseer association table" the fixtures encode.
    If that notion were wrong, the independent tally, the helper cross-check,
    the stderr cross-check and the parser would all agree with each other - and
    with this file's own model - while the threshold was wrong. So the two
    formats are checked against the installed package.

    These tests *import* pyseer; they do not invoke it. spec.md D8 defines
    TEST mode as running no real tools, and an end-to-end `PyseerEngine.run`
    against the binary would cross that line. What it would buy - observing
    pyseer's stdout and stderr rather than reading its formatter - is recorded
    as still-open on ticket 16 rather than smuggled in here.

    Nor is the import fully inside that line, and it is worth being precise
    rather than convenient: importing a real tool makes it a dependency of the
    test suite. That is why these tests use `importorskip` and skip, rather
    than erroring, on a machine without pyseer - so that "no real tools" means
    "the suite still runs", not "the suite needs one installed".

    What is *not* checked here, and cannot be from source: the column **names**
    in `PYSEER_LMM_HEADER`. `format_output` emits values, not names, so the
    names are read from `pyseer/__main__.py`'s `header` list. Only the column
    *count* is verified against pyseer. A wrong name fails loud rather than
    silent - `_normalise_columns` requires `variant` and `lrt-pvalue` and
    refuses by name otherwise - but it is not verified.
    """

    def setup_method(self):
        _require_pyseer()

    def test_the_pattern_encoding_is_pyseers_not_a_packed_byte_vector(self):
        """The fact a reimplementation gets wrong.

        `hash_pattern` MD5s `k.view(np.uint8)`, and `k` is an array of Python
        ints, so numpy gives it the platform's default integer width. On this
        machine that is 8 bytes per sample, and the digest is over all 64.
        Hashing one byte per value - the obvious reading - gives a different
        digest every time, so every fixture built that way describes a pattern
        pyseer could never emit.
        """
        from pyseer.input import hash_pattern

        vector = (1, 1, 1, 1, 0, 0, 0, 0)
        column = np.asarray([int(v) for v in vector])
        encoded = hash_pattern(column).decode().strip()

        assert encoded == pyseer_pattern(vector)
        # The naive reading, asserted to be different so the delegation above
        # is doing real work.
        packed = bytes(int(v) for v in vector)
        naive = base64.b64encode(hashlib.md5(packed).digest()).decode()
        assert naive != encoded, (
            "hashing one byte per sample agrees with pyseer on this platform, "
            "so the delegation to hash_pattern is not load-bearing here and the "
            "docstring's claim about 8 bytes per sample needs rechecking"
        )

    def test_a_different_vector_gives_a_different_pattern(self):
        """So that a pattern file is a usable proxy for distinct variants."""
        assert pyseer_pattern((1, 1, 1, 1, 0, 0, 0, 0)) != pyseer_pattern(
            (1, 1, 1, 1, 0, 0, 0, 1)
        )
        assert pyseer_pattern((1, 1, 1, 1, 0, 0, 0, 0)) == pyseer_pattern(
            (1, 1, 1, 1, 0, 0, 0, 0)
        )

    def test_the_pattern_is_one_line_and_carries_the_newline_pyseer_depends_on(self):
        """`sort -u | wc -l` counts *lines*.

        `hash_pattern` returns `binascii.b2a_base64(digest)`, which under
        Python 3 appends a newline. Without it the whole patterns file would be
        a single line and the count would be 1.
        """
        from pyseer.input import hash_pattern

        raw = hash_pattern(np.asarray([1, 0, 1, 0]))
        assert raw.endswith(b"\n"), "pyseer's pattern has no trailing newline"
        assert raw.count(b"\n") == 1

    def test_the_encoding_is_the_one_pyseer_builds_for_an_rtab(self, tmp_path: Path):
        """The load-bearing link, checked against pyseer's own reader.

        `pyseer_pattern` calls `hash_pattern` on an array this test file
        builds. That only means anything if pyseer hands `hash_pattern` an
        array of the same shape and dtype - and it does not build the column
        the obvious way: `read_variant` accumulates a `sample -> 0/1` dict and
        then does `np.array([d[x] for x in p.index if x in d])`, so numpy gives
        it the platform's default *integer* width and `hash_pattern`'s
        `k.view(np.uint8)` spans 8 bytes per sample here, not one.

        If pyseer ever packed the column to a single byte per value, every
        fixture in this file would describe a pattern pyseer could not emit -
        and asserting that `pyseer_pattern` agrees with `hash_pattern` would
        still pass, because both would be wrong the same way. So pyseer's
        reader is exercised here, on a real `.Rtab`.
        """
        import pandas as pd
        from pyseer.input import hash_pattern, read_variant

        vector = (1, 1, 1, 1, 0, 0, 0, 0)
        samples = [f"S{i:02d}" for i in range(len(vector))]

        rtab = tmp_path / "features.rtab"
        rtab.write_text(
            "Gene\t" + "\t".join(samples) + "\n"
            + "gene__ampC\t" + "\t".join(str(v) for v in vector) + "\n"
        )

        order = pd.DataFrame(index=samples)
        with rtab.open() as handle:
            # pyseer's own Rtab loader consumes the header to build `order`.
            handle.readline()
            _eof, column, name, _k, _nk, _af = read_variant(
                handle, order, "Rtab", None, None, False, set(samples), samples
            )

        assert name == "gene__ampC"
        assert list(column) == list(vector)
        # The dtype, stated rather than assumed: this is the whole point.
        assert column.dtype.kind == "i", f"column is {column.dtype}, not an int"
        assert column.dtype.itemsize > 1, (
            "pyseer's column is single-byte here, so hashing packed bytes WOULD "
            "have been correct and the docstring's claim about 8 bytes per "
            "sample needs rechecking"
        )
        assert column.view(np.uint8).size == column.size * column.dtype.itemsize

        # And the digest pyseer would write for that Rtab row is the one this
        # file's fixtures use.
        assert hash_pattern(column).decode().strip() == pyseer_pattern(vector)

    def test_the_tested_count_string_is_pyseers_own(self):
        """The stderr line 12a's cross-check parses, quoted from pyseer.

        That cross-check fails *open*: if the string ever stops matching, the
        line is simply not found and the check quietly stops firing — which
        would leave the docs claiming a guarantee that no longer exists. So the
        string is pinned against the installed package's own source, rather
        than left as a regex written from memory.
        """
        from pyseer import __main__ as pyseer_main

        source = inspect.getsource(pyseer_main)
        # `getsource` returns the *source text*, so pyseer's trailing newline
        # appears as the two characters backslash-n rather than as a newline.
        emitted = {
            text.replace("\\n", "")
            for text in re.findall(r"stderr\.write\('([^']*)' *%", source)
        }
        assert "%d tested variants" in emitted, (
            f"pyseer does not write that stderr line; it writes {sorted(emitted)}"
        )
        stage_pattern = re.compile(stage.PYSEER_TESTED_COUNT_PATTERN)

        # Render pyseer's real format string and check the stage's pattern
        # matches what the run's log will contain.
        assert stage_pattern.match("%d tested variants" % 1234).group(1) == "1234"

        # And it must not match the three neighbours, which pyseer writes to the
        # same stream: `loaded` is tested + prefiltered, so accepting it would
        # report a different population and, for most cohorts, never agree.
        for other in (
            "%d loaded variants",
            "%d filtered variants",
            "%d printed variants",
        ):
            assert other in emitted, f"pyseer no longer writes {other!r}"
            assert stage_pattern.match(other % 1234) is None, (
                f"the stage's pattern would match {other!r}"
            )

    def test_the_synthetic_association_table_has_pyseers_own_column_count(self):
        """The header constant must match what `format_output` emits.

        pyseer prints its header itself and builds each row with
        `format_output`; a row and its header disagreeing is how a parser ends
        up reading a beta as a p-value.
        """
        from pyseer.classes import LMM
        from pyseer.utils import format_output

        item = LMM(
            kmer="oprD_absent", pattern=b"synthetic", af=0.25, prep=1.0,
            pvalue=5e-5, kbeta=1.9, bse=0.5, frac_h2=0.7, max_lineage=None,
            kstrains=np.array([]), nkstrains=np.array([]), notes=[],
            prefilter=False, filter=False,
        )
        emitted = format_output(item, lmm=True)
        assert len(emitted.split("\t")) == len(PYSEER_LMM_HEADER.split("\t"))

    def test_the_synthetic_table_rows_come_from_pyseers_formatter(self):
        table = pyseer_lmm_table([
            ("oprD_absent", 0.25, 1.0, 5.0E-05, 1.9, 0.5, 0.7),
        ])
        body = table.splitlines()[1].split("\t")
        assert body[0] == "oprD_absent"
        # pyseer formats every numeric field in scientific notation, two
        # decimals; the parser reads them as floats, and a fixture written as
        # "0.25" would not be a row pyseer can produce.
        assert body[1] == "2.50E-01"
        assert float(body[3]) == pytest.approx(5.0e-05)

    def test_every_column_the_stage_renames_is_one_pyseer_emits(self):
        """The production rename map, checked against the header constant.

        Reads :data:`stage.PYSEER_COLUMN_RENAME` rather than restating it: a
        copy would pass even if the map changed. A rename for a column pyseer
        does not emit is dead code at best, and at worst a rename of the wrong
        column.
        """
        renamed = set(stage.PYSEER_COLUMN_RENAME)
        header = set(PYSEER_LMM_HEADER.split("\t"))
        assert renamed, "the rename map is empty - has it been restructured?"
        assert renamed <= header, (
            f"the stage renames columns pyseer's LMM header lacks: "
            f"{sorted(renamed - header)}"
        )
        # And the targets are the schema `parse_pyseer_output` requires.
        targets = set(stage.PYSEER_COLUMN_RENAME.values())
        assert {"feature", "pvalue"} <= targets

    def test_the_columns_the_stage_requires_are_the_ones_it_renames(self):
        """A header check and a rename map that disagree would be a dead gate.

        `_normalise_columns` refuses any file lacking `_PYSEER_HEADER_REQUIRED`.
        If one of those columns were not renamed, the file would be accepted and
        then read as if it already used the target names.
        """
        assert set(stage._PYSEER_HEADER_REQUIRED) <= set(
            stage.PYSEER_COLUMN_RENAME
        ), (
            "a column _normalise_columns requires is not one it renames, so the "
            "file is accepted and then read under the wrong schema"
        )


class TestVendoredHelperIsUpstreams:
    """The point of vendoring is that it *is* upstream's script."""

    def test_it_is_vendored_at_the_path_the_stage_calls(self, pipeline_root: Path):
        assert stage.COUNT_PATTERNS_HELPER == (
            pipeline_root / "scripts" / "gwas" / "count_patterns.py"
        )
        assert stage.COUNT_PATTERNS_HELPER.is_file(), (
            "pyseer 1.1.2 does not package scripts/count_patterns.py, so the "
            "repository vendors it. If it is missing, the multiple-testing "
            "threshold has no defined source."
        )

    def test_the_slicing_the_digest_is_taken_over_is_pinned(self):
        """The digest is only meaningful if the slice is defined.

        Taken from the copyright line so that this repository's provenance
        header - which is long, and exists precisely to document the
        vendoring - is excluded from the digest rather than pinning itself.
        If the marker moved, the digest would cover a different region and the
        comparison below would be meaningless rather than failing.
        """
        lines = stage.COUNT_PATTERNS_HELPER.read_text().splitlines()
        assert lines[0].startswith("#!"), "the vendored file must keep a shebang"
        assert sum(1 for line in lines if line.startswith("#!")) == 1, (
            "more than one shebang means the slice boundary is ambiguous"
        )
        assert lines.count(UPSTREAM_BODY_MARKER) == 1, (
            "the upstream copyright line must appear exactly once, at the "
            "start of the upstream body"
        )

    def test_the_logic_is_byte_identical_to_the_revision_it_was_taken_from(self):
        """An edited helper is a threshold from nowhere, wearing pyseer's name.

        The digest covers the body from the upstream copyright line onwards.
        The provenance header above it is this repository's and is excluded.
        """
        text = stage.COUNT_PATTERNS_HELPER.read_text()
        body = text[text.index(UPSTREAM_BODY_MARKER) :].rstrip("\n")
        assert (
            hashlib.sha256(body.encode()).hexdigest() == VENDORED_BODY_SHA256
        ), (
            "scripts/gwas/count_patterns.py no longer matches the pyseer "
            "revision its header names. Either restore it or re-vendor "
            "deliberately and update both digests - do not let them drift."
        )

    def test_the_header_records_the_same_digest_the_test_checks(self):
        """Header and test must not be able to disagree about the body.

        A digest held only in the test would be satisfied by any body, since
        whoever edits the logic would update the constant to match. Recording
        it in the file too means an edit to the logic and an edit to the
        recorded digest are two separate, individually visible changes.
        """
        text = stage.COUNT_PATTERNS_HELPER.read_text()
        assert VENDORED_BODY_SHA256 in text, (
            "the provenance header must record the sha256 of the logic below "
            "it, which is the digest the test suite enforces"
        )

    def test_the_header_names_the_upstream_source_and_its_whole_file_digest(self):
        """A human must be able to check the vendoring once, against upstream."""
        text = stage.COUNT_PATTERNS_HELPER.read_text()
        assert "github.com/mgalardini/pyseer" in text
        assert "Apache License 2.0" in text
        assert UPSTREAM_WHOLE_FILE_SHA256 in text, (
            "the header must record the sha256 of the upstream file as "
            "retrieved, so the vendoring can be checked against the original"
        )

    def test_it_actually_runs_as_a_script_and_reads_its_flags_from_its_own_parser(
        self, tmp_path: Path
    ):
        """The vendored file is executed, not imported.

        So its argparse block is the only thing defining its flags: if this
        fails, the names the pipeline passes on the command line are not the
        names the helper accepts. Run with `--help` and read the real output
        rather than grepping the source, because grepping the source cannot
        tell a declared flag from one that argparse rejected.
        """
        out = subprocess.run(
            [sys.executable, str(stage.COUNT_PATTERNS_HELPER), "--help"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for flag in ("--threshold", "--alpha", "--cores", "--memory", "--temp"):
            assert flag in out, f"{flag} is not in the helper's own --help"

    def test_the_helper_really_shells_out_to_sort_u_wc_l(self, tmp_path: Path):
        """Its whole logic is one shell pipeline; prove the pipeline, not the
        comment claiming it."""
        patterns = write_patterns(tmp_path / "patterns.txt", VECTORS)
        # If it were counting something else - lines, or uniq -c, or hashing
        # differently - the number would move. 24 distinct lines, 4 distinct
        # values, so a conflation of "lines" and "patterns" is visible.
        proc = subprocess.run(
            [
                sys.executable, str(stage.COUNT_PATTERNS_HELPER),
                str(patterns), "--alpha", "0.05", "--memory", "1024",
                "--cores", "1", "--temp", str(tmp_path),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        assert parse_helper(proc.stdout).patterns == EXPECTED_UNIQUE
        assert EXPECTED_TESTED != EXPECTED_UNIQUE, (
            "the fixture must contain redundant patterns, or this test cannot "
            "tell counting patterns from counting variants"
        )


class TestHelperAgreesWithSort:
    """The helper's output must be the number ``sort -u | wc -l`` gives.

    Not a restatement of the helper's own logic - a second implementation of
    the count, so that a wrong count is detectable.
    """

    def test_it_reproduces_sort_u_wc_l_on_a_synthetic_pyseer_patterns_file(
        self, tmp_path: Path
    ):
        patterns = write_patterns(tmp_path / "patterns.txt", VECTORS)
        expected = wc_lines("LC_ALL=C sort -u " + str(patterns) + " | wc -l")
        assert expected == EXPECTED_UNIQUE
        count, _ = stage.parse_count_patterns_output(
            subprocess.run(
                [
                    "python3",
                    str(stage.COUNT_PATTERNS_HELPER),
                    str(patterns),
                    "--alpha",
                    "0.05",
                    "--memory",
                    "1024",
                    "--cores",
                    "1",
                    "--temp",
                    str(tmp_path),
                ],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        )
        assert count == expected

    def test_duplicate_patterns_are_the_ones_that_collapse(self, tmp_path: Path):
        patterns = write_patterns(tmp_path / "patterns.txt", VECTORS)
        tally, total = stage.tally_patterns(patterns)
        assert total == EXPECTED_TESTED
        assert len(tally) == EXPECTED_UNIQUE
        # Two patterns are shared by two variants each; the rest are unique.
        assert sorted(tally.values()) == [1, 1, 2, 2]


# --------------------------------------------------------------------------
# The reduction
# --------------------------------------------------------------------------


class TestReduceUniquePatterns:
    def test_it_reports_the_survivor_count_and_what_was_dropped(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert reduction.n_variants_tested == EXPECTED_TESTED
        assert reduction.n_unique_patterns == EXPECTED_UNIQUE
        assert reduction.n_redundant_variants == EXPECTED_REDUNDANT
        assert reduction.n_variants_tested - reduction.n_redundant_variants == (
            reduction.n_unique_patterns
        )

    def test_the_threshold_is_alpha_over_unique_patterns(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert reduction.alpha == machine_config.gwas.unique_patterns.alpha
        assert reduction.threshold == pytest.approx(0.05 / EXPECTED_UNIQUE)

    def test_alpha_is_read_from_config_not_hardcoded(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        workdir: Path, tmp_path: Path,
    ):
        """Changing alpha in science.yaml must move the threshold.

        If the threshold were computed from a literal, this test would pass
        identically and the science would be unfalsifiable from the config.
        """
        # Copied, never symlinked. A symlinked `science.yaml` written through by
        # `write_text` rewrites the committed file, and a test that mutates the
        # repository's science while the rest of the suite reads it produces
        # failures in tests that have nothing to do with the one that caused
        # them.
        committed = (pipeline_root / "config" / "science.yaml").read_text()
        # `_resolve_root` takes the science file's grandparent as the pipeline
        # root and reads knowledge tables from `<root>/config/`, so the copy has
        # to reproduce that layout rather than just holding a science.yaml.
        alt_root = tmp_path / "altroot"
        shutil.copytree(pipeline_root / "config", alt_root / "config")
        altered_science = alt_root / "config" / "science.yaml"
        altered_science.write_text(
            committed.replace("    alpha: 0.05", "    alpha: 0.01")
        )
        assert altered_science.read_text() != committed

        original = load_config(
            pipeline_root / "config" / "science.yaml", machine=str(overlay)
        )
        changed = load_config(altered_science, machine=str(overlay))
        assert original.gwas.unique_patterns.alpha == 0.05
        assert changed.gwas.unique_patterns.alpha == 0.01

        base = stage.reduce_unique_patterns(original, patterns_file, workdir)
        moved = stage.reduce_unique_patterns(
            changed, patterns_file, tmp_path / "work2"
        )
        assert moved.threshold == pytest.approx(0.01 / EXPECTED_UNIQUE)
        assert moved.threshold < base.threshold
        # And the committed file is still what it was.
        assert (pipeline_root / "config" / "science.yaml").read_text() == committed

    def test_the_helpers_own_arithmetic_agrees_with_ours(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert f"{reduction.threshold:.2E}" == f"{reduction.helper_threshold:.2E}"

    def test_the_survivors_file_lists_one_row_per_unique_pattern(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        from papipeline.io.tsv import read_tsv

        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        rows = read_tsv(reduction.survivors_file, required_columns=("pattern",))
        assert len(rows) == EXPECTED_UNIQUE
        assert {r["pattern"] for r in rows} == {
            pyseer_pattern(v) for v in VECTORS
        }
        # Multiplicities reconcile against the input file: one row per distinct
        # value, and summing them gives back every tested variant.
        assert sum(int(r["n_variants"]) for r in rows) == EXPECTED_TESTED
        assert len({r["pattern"] for r in rows}) == len(rows), "duplicate rows"

    def test_the_survivors_file_reconciles_against_sort_u_wc_l(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        """The number a reader can check by hand must equal the number used."""
        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        by_hand = wc_lines("LC_ALL=C sort -u " + str(patterns_file) + " | wc -l")
        assert by_hand == reduction.n_unique_patterns

    def test_the_survivors_file_records_how_the_threshold_was_derived(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        text = reduction.survivors_file.read_text()
        assert "threshold = alpha / unique_patterns" in text
        assert f"{machine_config.gwas.unique_patterns.alpha} / {EXPECTED_UNIQUE}" in text

    def test_the_summary_is_written_and_is_readable_tsv(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        from papipeline.io.tsv import read_tsv

        reduction = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert reduction.summary_file == (
            workdir / stage.UNIQUE_PATTERNS_SUMMARY_FILENAME
        )
        rows = {
            r["metric"]: r["value"] for r in read_tsv(reduction.summary_file)
        }
        assert int(rows["variants_tested"]) == EXPECTED_TESTED
        assert int(rows["unique_patterns"]) == EXPECTED_UNIQUE
        assert int(rows["redundant_variants"]) == EXPECTED_REDUNDANT
        assert rows["method"] == "bonferroni"
        assert float(rows["alpha"]) == pytest.approx(0.05)
        assert float(rows["threshold"]) == pytest.approx(reduction.threshold)
        # The audit record has to name where the count came from, or it is a
        # number with no provenance.
        assert rows["patterns_file"] == str(patterns_file)
        assert rows["count_helper"] == str(stage.COUNT_PATTERNS_HELPER)

    def test_it_reruns_cleanly_over_its_own_outputs(
        self, machine_config, patterns_file: Path, workdir: Path
    ):
        """A rerun must produce the same files, not append to them.

        Snakemake re-runs rules; a step whose output grew on each run would make
        its own row count the thing under audit.
        """
        first = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        survivors_once = first.survivors_file.read_text()
        second = stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert second.survivors_file.read_text() == survivors_once
        assert second.n_unique_patterns == first.n_unique_patterns


class TestReduceRefuses:
    def test_a_missing_patterns_file_is_named_as_a_missing_intermediate(
        self, machine_config, tmp_path: Path, workdir: Path
    ):
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(
                machine_config, tmp_path / "absent.txt", workdir
            )
        assert "No patterns file" in str(exc.value)

    def test_an_empty_patterns_file_is_refused_rather_than_divided(
        self, machine_config, tmp_path: Path, workdir: Path
    ):
        """alpha / 0 is the failure mode; it gets a message, not a traceback."""
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(machine_config, empty, workdir)
        assert "tested no variant" in str(exc.value)

    def test_a_blank_lines_only_patterns_file_is_also_refused(
        self, machine_config, tmp_path: Path, workdir: Path
    ):
        blank = tmp_path / "blank.txt"
        blank.write_text("\n\n\n")
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(machine_config, blank, workdir)
        assert "tested no variant" in str(exc.value)

    def test_a_path_the_helper_cannot_safely_receive_is_refused_on_the_public_path(
        self, machine_config, tmp_path: Path, workdir: Path
    ):
        """The helper interpolates its file argument into a shell string.

        Upstream builds that command unquoted, and the vendored file must stay
        byte-identical, so the caller refuses rather than quotes. A space would
        change the command; a metacharacter would be executed.

        Driven through `reduce_unique_patterns`, not the private validator, so
        that the production path being unguarded would fail this.
        """
        hostile = tmp_path / "patterns with space.txt"
        hostile.write_text(pyseer_pattern((1, 1, 1, 1, 0, 0, 0, 0)) + "\n")
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(machine_config, hostile, workdir)
        assert "quoting" in str(exc.value)

    def test_a_scratch_directory_the_helper_cannot_safely_receive_is_refused(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path,
    ):
        """The scratch directory gets the same check, for the same reason.

        Upstream interpolates `options.temp` into the same unquoted shell
        string. It is lower risk than the patterns path - it comes from the
        overlay rather than the work directory - but it is the same defect, and
        guarding one of the two without a reason for the asymmetry would be
        worse than guarding neither.
        """
        hostile = tmp_path / "scratch dir"
        hostile.mkdir()
        text = overlay.read_text().replace(
            f"unique_patterns_temp_dir: {tmp_path / 'scratch'}",
            f"unique_patterns_temp_dir: {hostile}",
        )
        unsafe_overlay = tmp_path / "unsafe.yaml"
        unsafe_overlay.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(unsafe_overlay))
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(config, patterns_file, workdir)
        assert "scratch directory" in str(exc.value)
        assert "quoting" in str(exc.value)
        # The scratch directory has exactly one configuration key, and it is
        # the one a reader would go and edit.
        assert exc.value.context.get("key") == "runtime.unique_patterns_temp_dir"

    def test_a_shell_metacharacter_in_the_scratch_directory_is_refused(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path,
    ):
        """"; rm -rf ~" must never reach a `shell=True` string."""
        hostile = tmp_path / "scratch;touch"
        hostile.mkdir()
        text = overlay.read_text().replace(
            f"unique_patterns_temp_dir: {tmp_path / 'scratch'}",
            f"unique_patterns_temp_dir: {hostile}",
        )
        unsafe_overlay = tmp_path / "unsafe.yaml"
        unsafe_overlay.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(unsafe_overlay))
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(config, patterns_file, workdir)
        assert ";" in str(exc.value)
        assert not (tmp_path / "scratch;touch").joinpath("pwned").exists()

    def test_a_memory_below_the_helpers_reserve_is_refused_naming_the_key(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path,
    ):
        """The helper passes `--memory - 10` to `sort`.

        A value of 8 would reach `sort` as `-S -2M`, and sort's complaint would
        be about a buffer size, naming nothing that anyone set.
        """
        text = overlay.read_text().replace(
            "unique_patterns_memory_mb: 1024", "unique_patterns_memory_mb: 8"
        )
        tiny = tmp_path / "tiny.yaml"
        tiny.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(tiny))
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(config, patterns_file, workdir)
        assert "unique_patterns_memory_mb" in str(exc.value)

    def test_more_cores_than_this_machines_sort_supports_is_refused(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path, monkeypatch,
    ):
        """`--parallel` is not in POSIX, so it is probed rather than assumed.

        Where it is absent the helper would otherwise fail with a `sort` usage
        error naming nothing that anyone set. Stubbed to the absent case here
        because whether a given machine's sort has the flag is exactly the fact
        that varies - `/usr/bin/sort` 2.3-Apple on macOS does have it, and
        pretending otherwise in a test would be asserting a platform claim this
        repository has not verified.
        """
        import subprocess as _sp

        real_run = _sp.run

        def fake_run(cmd, *args, **kwargs):
            if cmd[:2] == ["sort", "--help"]:
                return SimpleNamespace(stdout="  -S, --buffer-size=SIZE\n", returncode=0)
            return real_run(cmd, *args, **kwargs)

        text = overlay.read_text().replace(
            "unique_patterns_cores: 1", "unique_patterns_cores: 4"
        )
        many = tmp_path / "many.yaml"
        many.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(many))
        monkeypatch.setattr(_sp, "run", fake_run)

        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(config, patterns_file, workdir)
        assert "unique_patterns_cores" in str(exc.value)
        assert "--parallel" in str(exc.value)

    def test_more_cores_is_accepted_where_sort_supports_it(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path, monkeypatch,
    ):
        """The cores>1 branch is exercised when the capability is present.

        Otherwise the only tested value is 1 and the branch `bigmachine.yaml`
        actually uses - 4 - is never run. Stubs `sort --help` to report the
        flag, so this test does not depend on which `sort` this machine has.
        """
        import subprocess as _sp

        real_run = _sp.run

        def fake_run(cmd, *args, **kwargs):
            if cmd[:2] == ["sort", "--help"]:
                return SimpleNamespace(
                    stdout="  --parallel=N  sort in parallel\n", returncode=0
                )
            return real_run(cmd, *args, **kwargs)

        text = overlay.read_text().replace(
            "unique_patterns_cores: 1", "unique_patterns_cores: 4"
        )
        many = tmp_path / "many.yaml"
        many.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(many))
        monkeypatch.setattr(_sp, "run", fake_run)

        reduction = stage.reduce_unique_patterns(config, patterns_file, workdir)
        # Same count either way - parallelism is `sort`'s business, not the
        # helper's - but the run completed rather than refusing.
        assert reduction.n_unique_patterns == EXPECTED_UNIQUE

    def test_a_missing_scratch_directory_is_refused_with_the_key_to_fix(
        self, overlay: Path, pipeline_root: Path, patterns_file: Path,
        tmp_path: Path, workdir: Path,
    ):
        text = overlay.read_text().replace(
            f"unique_patterns_temp_dir: {tmp_path / 'scratch'}",
            f"unique_patterns_temp_dir: {tmp_path / 'not-created'}",
        )
        absent = tmp_path / "no_scratch.yaml"
        absent.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(absent))
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(config, patterns_file, workdir)
        assert "unique_patterns_temp_dir" in str(exc.value)

    def test_a_helper_that_disagrees_with_the_independent_tally_stops_the_run(
        self, machine_config, patterns_file: Path, workdir: Path, monkeypatch
    ):
        """One count, or none.

        If the helper and the local tally see different numbers the run must
        stop. Choosing the more convenient one would mean the reported
        threshold and the reported reduction came from different files.
        """
        monkeypatch.setattr(
            stage,
            "count_unique_patterns",
            lambda config, path: (EXPECTED_UNIQUE + 1, 0.01),
        )
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert "disagreed" in str(exc.value)

    def test_a_helper_whose_own_threshold_disagrees_with_ours_stops_the_run(
        self, machine_config, patterns_file: Path, workdir: Path, monkeypatch
    ):
        """A count right and a threshold wrong means two different alphas.

        The helper's arithmetic is upstream's, so where it disagrees with ours,
        upstream's is the one the pipeline would be reporting as pyseer's - and
        one of the two is not the threshold this stage claims to use.
        """
        real = stage.count_unique_patterns
        monkeypatch.setattr(
            stage,
            "count_unique_patterns",
            lambda config, path: (EXPECTED_UNIQUE, 0.5),
        )
        assert real(machine_config, patterns_file)[0] == EXPECTED_UNIQUE
        with pytest.raises(DataContractError) as exc:
            stage.reduce_unique_patterns(machine_config, patterns_file, workdir)
        assert "does not agree" in str(exc.value)


class TestParseCountPatternsOutput:
    def test_it_reads_the_two_line_report(self):
        count, threshold = stage.parse_count_patterns_output(
            "Patterns:\t1234\nThreshold:\t4.05E-05\n"
        )
        assert count == 1234
        assert threshold == pytest.approx(4.05e-05)

    def test_a_zero_count_is_refused_not_divided(self):
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output("Patterns:\t0\nThreshold:\t0\n")
        assert "alpha / 0 is undefined" in str(exc.value)

    def test_a_changed_output_shape_is_refused_rather_than_guessed_at(self):
        """Guessing which line is which would put a number in the threshold."""
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output("1234\n4.05E-05\n")
        assert "no ':' separator" in str(exc.value)

    def test_a_missing_line_names_which_one(self):
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output("Patterns:\t1234\n")
        assert "Threshold" in str(exc.value)

    def test_an_extra_keyed_line_is_refused_not_ignored(self):
        """Two keys are expected. Silently dropping a third is how a renamed
        field becomes a silently-defaulted number."""
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output(
                "Patterns:\t1234\nThreshold:\t4.05E-05\nVerbose:\t1\n"
            )
        assert "unexpected Verbose" in str(exc.value)

    def test_a_repeated_key_is_refused_rather_than_one_of_the_two_taken(self):
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output(
                "Patterns:\t1234\nPatterns:\t99\nThreshold:\t4.05E-05\n"
            )
        assert "twice" in str(exc.value)

    def test_a_warning_line_before_the_report_is_refused_not_swallowed(self):
        """The helper writes nothing but its two lines. Anything else on stdout
        means it changed, or something else is writing to it."""
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output(
                "warning: something\nPatterns:\t1234\nThreshold:\t4.05E-05\n"
            )
        assert "unexpected warning" in str(exc.value)

    def test_a_non_numeric_count_is_refused(self):
        with pytest.raises(DataContractError) as exc:
            stage.parse_count_patterns_output("Patterns:\tmany\nThreshold:\t1E-05\n")
        assert "non-numeric" in str(exc.value)


class TestMachineConfigurationIsRequired:
    def test_the_helper_limits_come_from_the_overlay(self, machine_config):
        assert machine_config.unique_patterns_memory_mb == 1024
        assert machine_config.unique_patterns_cores == 1
        assert machine_config.unique_patterns_temp_dir.is_dir()

    def test_every_committed_overlay_declares_all_three_helper_limits(
        self, pipeline_root: Path
    ):
        """The accessors below are how the three keys get enforced.

        They raise at use rather than at load, because these three bound one
        `sort` in one stage and an overlay for a run that never reaches stage
        12 should still be valid. The price of that choice is that nothing
        stops an overlay quietly losing a key, so this pins all three in every
        committed overlay.
        """
        for name in ("laptop", "smoke", "bigmachine"):
            overlay = load_machine_config(
                pipeline_root / "config" / "machines" / f"{name}.yaml"
            )
            assert overlay.unique_patterns_memory_mb, f"{name}: memory_mb"
            assert overlay.unique_patterns_cores, f"{name}: cores"
            assert overlay.unique_patterns_temp_dir, f"{name}: temp_dir"

    def test_an_overlay_that_omits_a_helper_limit_refuses_when_it_is_read(
        self, tmp_path: Path, pipeline_root: Path
    ):
        """Absent keys are a refusal naming the key, never an inherited default.

        Inheriting the helper's own defaults instead would leave the run
        looking configured while `sort` used a 1-core, 1024 Mb, /tmp scratch
        nobody chose.
        """
        text = (pipeline_root / "config" / "machines" / "laptop.yaml").read_text()
        text = text.replace("  unique_patterns_memory_mb: 1024\n", "")
        path = tmp_path / "thin.yaml"
        path.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(path))
        with pytest.raises(ConfigError) as exc:
            config.unique_patterns_memory_mb
        assert "unique_patterns_memory_mb" in str(exc.value)

    def test_a_minimal_overlay_still_loads(
        self, tmp_path: Path, pipeline_root: Path
    ):
        """The keys are required at use, not at load, and that is deliberate.

        A bounded overlay for a run that stops before stage 12 has no business
        declaring that stage's `sort` limits, and refusing to load it would be a
        rule about a stage the run never reaches.
        """
        text = (pipeline_root / "config" / "machines" / "laptop.yaml").read_text()
        for key in (
            "  unique_patterns_memory_mb: 1024\n",
            "  unique_patterns_cores: 1\n",
            "  unique_patterns_temp_dir: /tmp\n",
        ):
            text = text.replace(key, "")
        path = tmp_path / "minimal.yaml"
        path.write_text(text)
        config = load_config(pipeline_root / "config" / "science.yaml", machine=str(path))
        assert config.threads == 4
        for attribute in (
            "unique_patterns_memory_mb",
            "unique_patterns_cores",
            "unique_patterns_temp_dir",
        ):
            with pytest.raises(ConfigError):
                getattr(config, attribute)

    def test_a_config_with_no_machine_overlay_refuses_rather_than_defaulting(
        self, pipeline_root: Path
    ):
        config = load_config(pipeline_root / "config" / "science.yaml", machine=None)
        for attribute in (
            "unique_patterns_memory_mb",
            "unique_patterns_cores",
            "unique_patterns_temp_dir",
        ):
            with pytest.raises(ConfigError):
                getattr(config, attribute)

    def test_only_the_correction_pyseer_documents_is_accepted(self):
        from papipeline.config.loader import UniquePatternConfig

        assert UniquePatternConfig.from_dict(
            {"method": "bonferroni", "alpha": 0.05}
        ).method == "bonferroni"
        with pytest.raises(ConfigError) as exc:
            UniquePatternConfig.from_dict(
                {"method": "benjamini-hochberg", "alpha": 0.05}
            )
        assert "does not implement" in str(exc.value)

    def test_both_keys_are_required_rather_than_defaulted(self):
        """Neither `alpha` nor `method` may be inherited from a default.

        `alpha` has the helper's own 0.05 to fall back on and `method` has the
        single supported correction, so a missing key would leave the run
        reporting a threshold as though someone had chosen it - the same
        failure the machine overlay's helper limits are built to refuse.
        """
        from papipeline.config.loader import UniquePatternConfig

        with pytest.raises(ConfigError) as exc:
            UniquePatternConfig.from_dict({"method": "bonferroni"})
        assert "alpha is required" in str(exc.value)

        with pytest.raises(ConfigError) as exc:
            UniquePatternConfig.from_dict({"alpha": 0.05})
        assert "method is required" in str(exc.value)

    def test_an_alpha_outside_zero_to_one_is_refused(self):
        from papipeline.config.loader import UniquePatternConfig

        for bad in (0.0, -0.1, 1.5):
            with pytest.raises(ConfigError):
                UniquePatternConfig.from_dict(
                    {"method": "bonferroni", "alpha": bad}
                )

    def test_a_config_without_the_unique_patterns_block_refuses_to_load(
        self, tmp_path: Path, pipeline_root: Path
    ):
        """`science.yaml` losing the block is a load-time error, not a default.

        The stage that reads it would not run in most modes, so a lazy check
        here would surface as a run that got as far as stage 12 and then
        refused - after everything upstream had already been done.
        """
        committed = (pipeline_root / "config" / "science.yaml").read_text()
        stripped = committed.split("  unique_patterns:")[0]
        alt_root = tmp_path / "altroot"
        shutil.copytree(pipeline_root / "config", alt_root / "config")
        (alt_root / "config" / "science.yaml").write_text(stripped)
        with pytest.raises(ConfigError) as exc:
            load_config(alt_root / "config" / "science.yaml", machine=None)
        assert "unique_patterns" in str(exc.value)
        assert (pipeline_root / "config" / "science.yaml").read_text() == committed


class TestThresholdReachesPyseer:
    """The derived threshold has to be the value pyseer is actually given."""

    def test_pyseer_is_asked_for_patterns_and_no_p_value_filter_in_pass_one(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """Pass 1's gates are stated open, and the reduction sees every pattern.

        The assertion is that the flags are passed with those values. What
        they *do* is a separate matter, and it is not what they do: pyseer
        writes a pattern before evaluating either gate, so they cannot change
        the patterns file. See `test_the_patterns_gate_flags_are_not_what_
        bounds_the_pattern_set` for the part that is load-bearing.
        """
        seen = {}

        def fake_invoke(self, cmd, table, log, phase):
            seen["phase"] = phase
            seen.setdefault("pass1", list(cmd))
            seen.setdefault("phases", []).append(phase)
            if "--output-patterns" in cmd:
                write_patterns(
                    Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS
                )
            Path(table).write_text(pyseer_lmm_table([]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        gwas_input = _gwas_input()
        results = engine.run(gwas_input, machine_config)

        cmd = seen["pass1"]
        assert any("12a" in phase for phase in seen["phases"])
        assert any("12b" in phase for phase in seen["phases"])
        # Stated open, so pyseer's defaults cannot change them later.
        assert cmd[cmd.index("--filter-pvalue") + 1] == "1"
        assert cmd[cmd.index("--lrt-pvalue") + 1] == "1"

        assert engine.reduction.n_unique_patterns == EXPECTED_UNIQUE
        assert engine.reduction.threshold == pytest.approx(0.05 / EXPECTED_UNIQUE)
        # Pass 1's association table is discarded, not reported.
        assert results == []

    def test_pass_two_carries_the_threshold_pass_one_derived(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        calls = []

        def fake_invoke(self, cmd, table, log, phase):
            calls.append((phase, list(cmd)))
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
                Path(table).write_text(pyseer_lmm_table([]))
            else:
                Path(table).write_text(pyseer_lmm_table([
                    ("gene__ampC", 0.5, 1.0, 1.0E-04, 2.5, 0.4, 0.6),
                ]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        results = engine.run(_gwas_input(), machine_config)

        assert len(calls) == 2
        pass_two = calls[1][1]
        threshold = float(pass_two[pass_two.index("--lrt-pvalue") + 1])
        assert threshold == pytest.approx(0.05 / EXPECTED_UNIQUE)
        assert threshold == pytest.approx(engine.reduction.threshold)

        # And the threshold reached pyseer agrees with what 12a wrote down.
        from papipeline.io.tsv import read_tsv

        summary = {
            r["metric"]: r["value"]
            for r in read_tsv(engine.reduction.summary_file)
        }
        assert float(summary["threshold"]) == pytest.approx(threshold)

        assert [r.feature for r in results] == ["gene__ampC"]

    def test_the_result_records_that_the_threshold_was_bonferroni_on_patterns(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        def fake_invoke(self, cmd, table, log, phase):
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
            Path(table).write_text(pyseer_lmm_table([
                ("oprD_absent", 0.25, 1.0, 5.0E-05, 1.9, 0.5, 0.7),
            ]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        results = engine.run(_gwas_input(), machine_config)
        assert "bonferroni_unique_patterns" in results[0].model

    def test_a_reduction_failure_stops_before_any_association_table_is_parsed(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """An undefined threshold must not degrade into an unthresholded run."""
        calls = []

        def fake_invoke(self, cmd, table, log, phase):
            calls.append(phase)
            if "--output-patterns" in cmd:
                Path(cmd[cmd.index("--output-patterns") + 1]).write_text("")

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        with pytest.raises(DataContractError):
            engine.run(_gwas_input(), machine_config)
        assert len(calls) == 1, "12b ran despite 12a failing"

    def test_no_variant_passing_the_threshold_is_an_empty_result_not_a_failure(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """A cohort with no association reports none; it does not crash.

        pyseer prints its header and no rows when nothing passes
        --lrt-pvalue. Before 12a that could not happen, because the stage ran
        with pyseer's default `--lrt-pvalue 1`, which passes everything. Now
        that pass 12b runs at a derived cutoff, "nothing passes" is the ordinary
        outcome of a null result and has to be representable.
        """
        def fake_invoke(self, cmd, table, log, phase):
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
            Path(table).write_text(pyseer_lmm_table([]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        results = engine.run(_gwas_input(), machine_config)
        assert results == []
        # The threshold that produced the empty table is still on disk.
        assert engine.reduction.summary_file.is_file()

    def test_a_truncated_pyseer_output_is_refused_not_read_as_no_associations(
        self, tmp_path: Path
    ):
        """One line that is not a header is neither of the two empty shapes.

        Three situations present as "a line or fewer" and only two are
        results: no lines at all, and a header with no rows. A single
        non-header line is a truncated write or a stray diagnostic, and reading
        it as "pyseer found nothing" would convert a broken run into a reported
        negative result - which is the opposite of what splitting the threshold
        into its own step is for.
        """
        raw = tmp_path / "pyseer_raw.tsv"
        raw.write_text("Traceback (most recent call last):\n")
        normalised = stage.PyseerEngine._normalise_columns(raw)
        with pytest.raises(DataContractError) as exc:
            stage.parse_pyseer_output(normalised, _gwas_input())
        assert "no 'feature' column" in str(exc.value)
        assert "Traceback" in str(exc.value)

    def test_a_header_with_no_rows_is_an_empty_result(
        self, tmp_path: Path
    ):
        """The contrast with the case above, asserted together so it holds."""
        raw = tmp_path / "pyseer_raw.tsv"
        raw.write_text(PYSEER_LMM_HEADER + "\n")
        normalised = stage.PyseerEngine._normalise_columns(raw)
        assert stage.parse_pyseer_output(normalised, _gwas_input()) == []

    def test_a_completely_empty_pyseer_output_is_refused_not_read_as_no_associations(
        self, tmp_path: Path
    ):
        """Distinct from header-only, which is a real empty result.

        With the real binary an exit-0 run always writes a header before the
        association loop, so a zero-byte table is a truncated or unwritten file.
        Reading it as "pyseer found nothing" would turn a broken run into a
        reported negative result - the same defect as the non-header single
        line, one branch over.
        """
        raw = tmp_path / "pyseer_raw.tsv"
        raw.write_text("")
        normalised = stage.PyseerEngine._normalise_columns(raw)
        with pytest.raises(DataContractError) as exc:
            stage.parse_pyseer_output(normalised, _gwas_input())
        assert "empty output file" in str(exc.value)

    def test_every_branch_of_normalise_rewrites_its_target(self, tmp_path: Path):
        """A stale `pyseer_results.tsv` must not survive into a later read.

        The parser is handed whatever `_normalise_columns` returns, so if any
        branch returned early without rewriting it, the previous run's results
        would be read as this run's.
        """
        cases = {
            "empty": "",
            "garbage": "not a table\n",
            "header_only": PYSEER_LMM_HEADER + "\n",
            "with_rows": pyseer_lmm_table([
                ("oprD_absent", 0.25, 1.0, 5.0E-05, 1.9, 0.5, 0.7),
            ]),
        }
        for name, body in cases.items():
            target = tmp_path / f"case_{name}"
            target.mkdir()
            raw = target / "pyseer_raw.tsv"
            raw.write_text(body)
            # A recognisable stale artefact from a previous run.
            stale = target / "pyseer_results.tsv"
            stale.write_text("feature\tpvalue\nSTALE_FROM_A_PREVIOUS_RUN\t1\n")

            normalised = stage.PyseerEngine._normalise_columns(raw)
            assert normalised == stale, f"{name}: returned a different path"
            assert "STALE_FROM_A_PREVIOUS_RUN" not in normalised.read_text(), (
                f"{name}: left the previous run's results in place"
            )

    def test_the_frequency_bounds_not_the_p_value_gates_bound_the_pattern_set(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """What actually determines the pattern set, asserted rather than assumed.

        The pass-1 command states `--filter-pvalue 1 --lrt-pvalue 1`, and it is
        tempting to read those as what guarantees every eligible variant is
        counted. They do not: `pyseer/__main__.py` writes a pattern for every
        non-prefilter variant *before* evaluating `x.filter`, so the p-value
        gates cannot remove one. `--min-af`/`--max-af` are what bound the set.

        So the invariant that matters is about the frequency bounds, and it is
        asserted here on the command that carries them.
        """
        seen = {}

        def fake_invoke(self, cmd, table, log, phase):
            seen.setdefault("pass1", list(cmd))
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
            Path(table).write_text(pyseer_lmm_table([]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        engine.run(_gwas_input(), machine_config)

        cmd = seen["pass1"]
        n_samples = len(_gwas_input().sample_ids)
        assert float(cmd[cmd.index("--min-af") + 1]) == pytest.approx(1.0 / n_samples)
        assert float(cmd[cmd.index("--max-af") + 1]) == pytest.approx(
            1.0 - 1.0 / n_samples
        )
        # And the two passes agree on them, since they come from one place.
        assert engine.reduction.n_variants_tested == EXPECTED_TESTED

    def test_pyseer_column_names_are_mapped_off_a_synthetic_pyseer_table(
        self, tmp_path: Path
    ):
        """`variant` and `lrt-pvalue` mapped onto the schema the parser reads.

        `af` is renamed to `freq` too, but deliberately not asserted against a
        value: `parse_pyseer_output` computes each feature's frequency from
        `gwas_input.features` rather than trusting the engine's, so a synthetic
        `af` here would be checked against something the parser ignores.
        """
        raw = tmp_path / "pyseer_raw.tsv"
        raw.write_text(pyseer_lmm_table([
            ("oprD_absent", 0.25, 1.0, 5.0E-05, 1.9, 0.5, 0.7),
            ("gene__ampC", 0.50, 1.0, 4.0E-01, 0.1, 0.9, 0.6),
        ]))
        normalised = stage.PyseerEngine._normalise_columns(raw)
        header = normalised.read_text().splitlines()[0].split("\t")
        assert header[0] == "feature"
        assert "pvalue" in header and "freq" in header
        assert "lrt-pvalue" not in header and "variant" not in header

        results = stage.parse_pyseer_output(normalised, _gwas_input())
        assert [r.feature for r in results] == ["oprD_absent", "gene__ampC"]
        assert results[0].p_value == pytest.approx(5.0e-05)
        # `oprD_absent` is carried by exactly half of the eight-sample cohort in
        # _gwas_input, so this is the locally computed frequency.
        assert results[0].frequency == pytest.approx(0.5)

    def test_the_parsed_feature_is_not_one_pyseer_never_reported(
        self, tmp_path: Path
    ):
        """A feature pyseer did not report gets no frequency, not a wrong one.

        `parse_pyseer_output` looks the feature up in the input matrix. An
        unknown name must not fall back to a neighbour's numbers.
        """
        raw = tmp_path / "pyseer_raw.tsv"
        raw.write_text(pyseer_lmm_table([
            ("a_feature_pyseer_never_saw", 0.25, 1.0, 5.0E-05, 1.9, 0.5, 0.7),
        ]))
        normalised = stage.PyseerEngine._normalise_columns(raw)
        results = stage.parse_pyseer_output(normalised, _gwas_input())
        assert len(results) == 1
        assert results[0].frequency is None
        assert results[0].lineage_distribution == {}


class TestPatternCountMatchesPyseer:
    """The patterns file is checked against pyseer's own tested-variant tally.

    pyseer writes `%d tested variants` to stderr. That is its own statement of
    how many variants reached the model, and it is the population the reduction
    runs over - so the two disagreeing means every count derived from the file
    rests on a population pyseer does not endorse.

    To be exact about the arithmetic, because this is the easiest thing in the
    stage to get backwards: the *divisor* in `alpha / n_unique_patterns` is the
    **unique**-pattern count. This check validates the *input* population, not
    the divisor.
    """

    @staticmethod
    def _run(machine_config, monkeypatch, tmp_path: Path, log_body: str):
        def fake_invoke(self, cmd, table, log, phase):
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
                Path(log).write_text(log_body)
            Path(table).write_text(pyseer_lmm_table([]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        return engine, engine.run(_gwas_input(), machine_config)

    def test_agreement_is_the_normal_case(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        _, results = self._run(
            machine_config, monkeypatch, tmp_path, f"{EXPECTED_TESTED} tested variants\n"
        )
        assert results == []

    def test_a_disagreement_stops_the_run(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """pyseer says one number, the file holds another: no threshold.

        Not a warning. The threshold would be divided by one of them, and
        reporting either without saying so would be the failure.
        """
        with pytest.raises(DataContractError) as exc:
            self._run(
                machine_config, monkeypatch, tmp_path, f"{EXPECTED_TESTED + 3} tested variants\n"
            )
        assert "tested" in str(exc.value)
        assert str(EXPECTED_TESTED + 3) in str(exc.value)

    def test_a_log_without_the_line_leaves_the_reduction_standing(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """No tally to compare against is not a disagreement.

        A caller who captured pyseer's stderr elsewhere leaves nothing to check,
        and inventing a second opinion would be worse than recording that there
        was none. The other two checks - the helper and the local tally - are
        unaffected.
        """
        engine, _ = self._run(
            machine_config, monkeypatch, tmp_path, "no counts here\n"
        )
        assert engine.reduction.n_unique_patterns == EXPECTED_UNIQUE

    def test_a_missing_log_is_not_a_disagreement(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        def fake_invoke(self, cmd, table, log, phase):
            if "--output-patterns" in cmd:
                write_patterns(Path(cmd[cmd.index("--output-patterns") + 1]), VECTORS)
                # No log written at all.
            Path(table).write_text(pyseer_lmm_table([]))

        monkeypatch.setattr(stage.PyseerEngine, "_invoke", fake_invoke)
        engine = stage.PyseerEngine("pyseer", tmp_path / "wd")
        engine.run(_gwas_input(), machine_config)
        assert engine.reduction.n_variants_tested == EXPECTED_TESTED

    def test_the_count_is_among_the_lines_pyseer_reports_not_derived_from_them(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """A full pyseer log carries loaded/filtered/tested/printed counts.

        Parsed out of a realistic log rather than a one-line stub, so that
        picking the wrong line - say `loaded`, which is tested plus prefiltered
        - would fail rather than pass.
        """
        _, results = self._run(
            machine_config,
            monkeypatch,
            tmp_path,
            f"{EXPECTED_TESTED + 2} loaded variants\n"
            "2 filtered variants\n"
            f"{EXPECTED_TESTED} tested variants\n"
            "0 printed variants\n",
        )
        assert results == []

    def test_a_loaded_count_that_differs_does_not_trip_the_check(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """`loaded` is `tested + prefiltered`, so it is always >= `tested`.

        Reading the wrong line would raise on every run of a cohort with any
        prefiltered variant, because `loaded` would never equal the patterns
        file's line count. Matching here proves the parser reads `tested`.
        """
        _, results = self._run(
            machine_config,
            monkeypatch,
            tmp_path,
            f"{EXPECTED_TESTED + 40} loaded variants\n"
            "40 filtered variants\n"
            f"{EXPECTED_TESTED} tested variants\n"
            "0 printed variants\n",
        )
        assert results == []

    def test_a_disagreement_further_down_a_full_log_is_still_caught(
        self, machine_config, monkeypatch, tmp_path: Path
    ):
        """The line is found wherever it is, not only in a one-line stub."""
        with pytest.raises(DataContractError) as exc:
            self._run(
                machine_config,
                monkeypatch,
                tmp_path,
                f"{EXPECTED_TESTED} loaded variants\n"
                "0 filtered variants\n"
                f"{EXPECTED_TESTED + 1} tested variants\n"
                "3 printed variants\n",
            )
        assert str(EXPECTED_TESTED + 1) in str(exc.value)


def _gwas_input() -> "stage.GwasInput":
    """A small cohort with two feature families and two lineages."""
    samples = tuple(f"S{i:02d}" for i in range(8))
    binary = {"S00": 1, "S01": 1, "S02": 1, "S03": 1,
              "S04": 0, "S05": 0, "S06": 0, "S07": 0}
    return stage.GwasInput(
        sample_ids=samples,
        binary_outcome=binary,
        features={
            "gene__ampC": {s: binary[s] for s in samples},
            "oprD_absent": {s: 1 - binary[s] for s in samples},
        },
        feature_types={"gene__ampC": "gene", "oprD_absent": "gene"},
        lineages={s: ("L1" if binary[s] else "L2") for s in samples},
        n_positive=4,
        n_negative=4,
        n_excluded=0,
    )