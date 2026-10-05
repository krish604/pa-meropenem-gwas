"""The oprD locus resolver is a declared input of the run, not an afterthought.

**The invariant this file exists to protect.** ``oprD_absent`` may come ONLY
from an ``absent`` verdict and ``oprD_LoF`` only from ``disrupted``. Everything
else - including *nothing at all* - must land on ``not_assessed``. A missing
input that produced ``absent`` would manufacture the study's central negative
from a file that was never read, which is the failure
:mod:`papipeline.adapters.oprd_locus` was written to prevent and the failure
that the annotation-symbol fallback in :func:`papipeline.viz.oprd_status_per_sample`
does produce on real data.

**Why the mechanism is a coverage assertion, not a code review.** A mapping
that omits a sample lets ``oprd_status_per_sample`` fall through to the symbol
test *for that sample*, silently reintroducing the fabricated negative. So the
assertion is: every manifest sample appears in the returned mapping, whatever
the resolver concluded about it.

**What is NOT tested here.** No blastp is run. ``blastp_command`` and
``run_blastp`` are covered in ``test_oprd_locus_resolver.py``, including a
non-zero exit. These tests cover the wiring in ``run.py``: which paths are
declared, what happens when one is absent, and what the verdict mapping is
allowed to contain.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.adapters.oprd_locus import Verdict
from papipeline.errors import StageError
from papipeline.manifest import SampleManifest
from papipeline.models import Sample, RunMode
from papipeline.run import OPRD_LOCUS_WORK_DIRNAME, resolve_oprd_loci
from papipeline.viz import oprd_status_per_sample

SAMPLES = ("isolate_A", "isolate_B", "isolate_C")


def _manifest(sample_ids=SAMPLES) -> SampleManifest:
    return SampleManifest(
        samples=tuple(
            Sample(sample_id=sid, assembly_path=None, source="test")
            for sid in sample_ids
        )
    )


class _FakeBlastp:
    """Stands in for the ``blastp`` on PATH.

    ``resolve_oprd_loci`` resolves the executable before it touches any input,
    so a test for the missing-input refusals has to get past that check. This
    replaces ``shutil.which`` for the name only, and points at this very file -
    a real, executable file, so the check under test is "did we resolve a
    program" rather than "is the string non-empty".
    """

    def __call__(self, name, *args, **kwargs):
        if name == "blastp":
            return str(Path(__file__).resolve())
        return shutil.which(name, *args, **kwargs)


def _patch_blastp(monkeypatch):
    monkeypatch.setattr("papipeline.run.shutil.which", _FakeBlastp())


def _write_inputs(tmp_path: Path, sample_ids=SAMPLES, *, assembly=True,
                   faa=True, gff3=True):
    """Materialise the declared per-sample inputs under the intermediate root.

    Assemblies go in ``genomes/<sample_id>/<sample_id>.fna`` because that is the
    layout `locate_assembly` rule 2 looks for. A flat ``genomes/<sample_id>.fna``
    is matched by neither rule, and every test below that wants a *present*
    assembly would then be asserting the missing-assembly refusal instead.
    """
    genomes = tmp_path / "genomes"
    intermediate = tmp_path / "intermediate"
    genomes.mkdir(parents=True, exist_ok=True)
    for sid in sample_ids:
        if assembly:
            (genomes / sid).mkdir(parents=True, exist_ok=True)
            (genomes / sid / f"{sid}.fna").write_text(
                ">x\nACGT\n", encoding="utf-8"
            )
        bakta = intermediate / "bakta" / sid
        bakta.mkdir(parents=True, exist_ok=True)
        if gff3:
            (bakta / f"{sid}.gff3").write_text("##gff-version 3\n", encoding="utf-8")
        if faa:
            (bakta / f"{sid}.faa").write_text(">gene_1\nMKV\n", encoding="utf-8")
    return genomes, intermediate


class _Config:
    """The minimum the resolver reads: threads, and the pinned reference.

    ``reference_protein`` is patched out in every test below, so the two
    reference accessors only have to exist and never have to resolve - the
    point of each test is what happens BEFORE or AROUND the search, not the
    search. A real ``PipelineConfig`` is not used because it would drag in
    ``db/reference/``, which is untracked, and a test that skips when a database
    is absent cannot protect an invariant.
    """

    def __init__(self, tmp_path=None):
        self.runtime = {"threads": 2}
        root = Path(tmp_path) if tmp_path is not None else Path("/nonexistent")
        self._reference = root / "reference.fna"

    def reference_fasta(self):
        return self._reference

    def reference_gff(self):
        return self._reference.with_suffix(".gff")


@pytest.fixture(autouse=True)
def no_reference_protein(monkeypatch):
    """Replace the PAO1 query derivation with a fixed 443-residue string.

    Autouse, because no test here is about the query: they are about what
    happens before and around the search.

    ``reference_protein`` is a pure function over two pinned files and is
    covered directly in ``test_oprd_locus_resolver.py``. Leaving it in would
    make every test here depend on ``db/reference/`` being populated, which is
    an untracked directory - a test that skips when a database is absent cannot
    protect an invariant.
    """
    monkeypatch.setattr(
        "papipeline.run.reference_protein",
        lambda gff, fasta, tag: "M" + "A" * 442,
    )


class TestAMissingInputIsRefusedByName:
    def test_a_missing_assembly_names_the_sample(self, tmp_path, monkeypatch,):
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path, faa=False, gff3=False)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        message = str(excinfo.value)
        assert "assembly" in message
        for sid in SAMPLES:
            assert sid in message, f"{sid} was not named in the refusal"

    def test_a_missing_assembly_names_the_directory_it_looked_in(
        self, tmp_path, monkeypatch,
    ):
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path, assembly=False)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        assert str(genomes) in str(excinfo.value), (
            "the refusal must say where the assembly was expected, or a reader "
            "cannot tell an absent cohort from a wrong genomes dir"
        )

    def test_a_missing_proteome_names_the_file_and_the_stage_that_writes_it(
        self, tmp_path, monkeypatch,
    ):
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path, faa=False)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        message = str(excinfo.value)
        assert ".faa" in message
        assert "Bakta" in message, (
            "a refusal that names the file but not what produces it sends the "
            "reader to look for a tool they do not have"
        )

    def test_a_missing_gff_names_the_file(self, tmp_path, monkeypatch,):
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path, gff3=False)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        assert ".gff3" in str(excinfo.value)

    def test_every_missing_input_is_reported_at_once(
        self, tmp_path, monkeypatch,
    ):
        """One refusal, not one per isolate.

        A loop that raised on the first absent assembly would report one sample
        of ten and send the reader round the loop nine more times.
        """
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path, assembly=False)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        message = str(excinfo.value)
        assert message.count("assembly for") == len(SAMPLES)

    def test_a_missing_blastp_is_refused_before_any_input_is_read(
        self, tmp_path, monkeypatch,
    ):
        """The tool check comes first, and the message says it is not `absent`.

        An operator whose blastp is missing should not be told their isolates
        lack oprD, and should not be walked through a cohort whose inputs are
        fine.
        """
        monkeypatch.setattr("papipeline.run.shutil.which", lambda name: None)
        genomes, intermediate = _write_inputs(tmp_path)
        with pytest.raises(StageError) as excinfo:
            resolve_oprd_loci(
                _Config(tmp_path), _manifest(), RunMode.REAL,
                genomes_dir=genomes, intermediate_root=intermediate,
            )
        message = str(excinfo.value)
        assert "blastp" in message
        assert "not the same as an isolate that lacks oprD" in message

    def test_the_work_directory_lives_under_the_intermediate_root(
        self, tmp_path, monkeypatch,
    ):
        """Named because it is the invariant AGENTS.md rule 2 rests on.

        A literal ``oprd_locus`` relative to the repository would put run output
        into a tracked tree, and the probe alignment written into it is the
        auditable record of what blastp was asked.
        """
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path)
        seen = {}

        def fake_resolve_isolate(sample_id, proteome, gff3, reference, *, workdir,
                                 **kwargs):
            seen[sample_id] = workdir
            from papipeline.adapters.oprd_locus import LocusResolution
            return LocusResolution(
                sample_id=sample_id, verdict=Verdict.RESOLVED, reason="test",
            )

        monkeypatch.setattr(
            "papipeline.run.resolve_isolate", fake_resolve_isolate
        )
        resolve_oprd_loci(
            _Config(tmp_path), _manifest(), RunMode.REAL,
            genomes_dir=genomes, intermediate_root=intermediate,
        )
        assert set(seen) == set(SAMPLES)
        for sample_id, workdir in seen.items():
            assert Path(workdir).is_relative_to(intermediate), (
                f"{sample_id}'s scratch at {workdir} is outside the run's own "
                "intermediate root"
            )
        assert seen[SAMPLES[0]].parent.name == OPRD_LOCUS_WORK_DIRNAME
        assert seen[SAMPLES[0]].parent.parent == intermediate


class TestTheVerdictMappingIsComplete:
    """The coverage invariant, which is what makes the refusal mapping safe."""

    def _run(self, tmp_path, monkeypatch, verdict_for):
        _patch_blastp(monkeypatch)
        genomes, intermediate = _write_inputs(tmp_path)
        from papipeline.adapters.oprd_locus import LocusResolution

        def fake_resolve_isolate(sample_id, *args, **kwargs):
            return LocusResolution(
                sample_id=sample_id, verdict=verdict_for(sample_id),
                reason="test",
            )

        monkeypatch.setattr(
            "papipeline.run.resolve_isolate", fake_resolve_isolate
        )
        return resolve_oprd_loci(
            _Config(tmp_path), _manifest(), RunMode.REAL,
            genomes_dir=genomes, intermediate_root=intermediate,
        )

    def test_every_manifest_sample_appears_whatever_the_verdict(
        self, tmp_path, monkeypatch,
    ):
        """Mixed outcomes, including two different refusals.

        The mixed case is the one that matters: a mapping that carried only the
        resolved isolates would look correct and would leave the other two to
        the symbol test.
        """
        verdicts = {
            "isolate_A": Verdict.RESOLVED,
            "isolate_B": Verdict.NO_HIT,
            "isolate_C": Verdict.INSUFFICIENT,
        }
        resolutions = self._run(
            tmp_path, monkeypatch, lambda sid: verdicts[sid]
        )
        assert set(resolutions) == set(SAMPLES)

    def test_a_mapping_missing_a_sample_reintroduces_the_symbol_test(self):
        """The failure this prevents, demonstrated.

        Same isolate, two ways of describing the evidence: one mapping that
        covers the cohort and one that omits it. Omitting a sample hands that
        sample back to the annotation-symbol test, and completeness of the
        mapping is therefore load-bearing and not tidiness.

        What changed is what the symbol test is now allowed to conclude. It used
        to report ``absent`` for an isolate with no ``oprD`` symbol, which is
        wrong on real data 34 times out of 36. Since ``fd38486`` ("stop the
        symbol test faking absence") it reports ``not_assessed`` instead:
        ``papipeline/viz.py`` reads
        ``"intact" if "oprD" in set(genes) else "not_assessed"``, where
        pre-``fd38486`` it read ``else "absent"``. So an incomplete mapping is
        now a *silent* degradation - the result stays safe, but the locus
        resolution stops being authoritative and nothing raises.
        """
        annotated = {"isolate_A": ["mexS"], "isolate_B": ["nalC"]}
        covered = {
            "isolate_A": _resolution("isolate_A", Verdict.NO_HIT),
        }
        status = oprd_status_per_sample({}, annotated, covered)
        assert status == {
            "isolate_A": "not_assessed",
            # Omitted from `covered`, so the symbol test decided - and it found
            # no `oprD` symbol. That is not evidence the gene is gone, so it
            # must be reported as unassessed rather than as `absent`.
            "isolate_B": "not_assessed",
        }
        # Stated separately so the safe-value contract is checked in its own
        # right, and not merely as a side effect of the dict above.
        assert "absent" not in status.values(), (
            "an incomplete locus mapping let the annotation-symbol test report "
            "`absent`. Not seeing an `oprD` symbol is not evidence of deleting "
            "the gene."
        )

    def test_a_complete_mapping_yields_no_absent_and_no_intact_from_the_symbol(
        self, tmp_path, monkeypatch,
    ):
        """The single most important assertion in this task.

        With the mapping covering the cohort, no status may come from the symbol
        test: the only states reachable are ``intact`` (from a resolved locus)
        and ``not_assessed`` (from a refusal). ``absent`` is unreachable, and
        that is the point.
        """
        resolutions = self._run(
            tmp_path, monkeypatch,
            lambda sid: (
                Verdict.RESOLVED if sid == "isolate_A" else Verdict.NO_HIT
            ),
        )
        # Annotated with every sample deliberately lacking an `oprD` symbol, so
        # the symbol test would say `absent` for all three if it ran.
        annotated = {sid: ["mexS", "nalC"] for sid in SAMPLES}
        status = oprd_status_per_sample({}, annotated, resolutions)
        assert status == {
            "isolate_A": "intact",
            "isolate_B": "not_assessed",
            "isolate_C": "not_assessed",
        }
        assert "absent" not in status.values(), (
            "a refused locus resolution produced `absent`. Not locating a locus "
            "is not evidence of deleting it."
        )

    def test_a_refused_resolution_is_never_absent(self, tmp_path, monkeypatch):
        """Every refusal verdict, asserted individually.

        A loop over ``Verdict.REFUSALS`` rather than one representative case,
        because a new refusal added to that frozenset would otherwise inherit
        this bug silently.
        """
        annotated = {sid: ["mexS"] for sid in SAMPLES}
        for verdict in sorted(Verdict.REFUSALS):
            resolutions = {
                sid: _resolution(sid, verdict) for sid in SAMPLES
            }
            status = oprd_status_per_sample({}, annotated, resolutions)
            assert set(status.values()) == {"not_assessed"}, (
                f"{verdict} produced {sorted(set(status.values()))}"
            )

    def test_outside_real_the_mapping_is_none_not_empty(self, tmp_path,
                                                        monkeypatch):
        """`None` and `{}` are different, and only one of them is safe.

        ``oprd_status_per_sample`` falls back to the symbol test when the
        argument is ``None``. Passing an empty dict instead would look like "the
        resolver ran and resolved nobody" and would still leave the symbol test
        in charge - so TEST would silently keep the behaviour this wiring exists
        to remove, without anything failing.
        """
        genomes, intermediate = _write_inputs(tmp_path)
        assert resolve_oprd_loci(
            _Config(tmp_path), _manifest(), RunMode.TEST,
            genomes_dir=genomes, intermediate_root=intermediate,
        ) is None
        assert resolve_oprd_loci(
            _Config(tmp_path), _manifest(), RunMode.STUB,
            genomes_dir=genomes, intermediate_root=intermediate,
        ) is None

    def test_outside_real_a_missing_assembly_is_not_an_error(self, tmp_path,
                                                             monkeypatch,):
        """TEST has no assemblies and no blastp, and must not start refusing.

        The inputs are absent *and* blastp is unresolvable, and neither stops a
        TEST run - otherwise this wiring would break every TEST run in the
        suite, which is the failure mode that makes a change like this get
        reverted rather than fixed.
        """
        monkeypatch.setattr("papipeline.run.shutil.which", lambda name: None)
        genomes, intermediate = _write_inputs(tmp_path, assembly=False)
        assert resolve_oprd_loci(
            _Config(tmp_path), _manifest(), RunMode.TEST,
            genomes_dir=genomes, intermediate_root=intermediate,
        ) is None


def _resolution(sample_id, verdict):
    from papipeline.adapters.oprd_locus import LocusResolution

    return LocusResolution(
        sample_id=sample_id, verdict=verdict, reason="test",
    )