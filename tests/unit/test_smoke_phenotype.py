"""A smoke run reads phenotype from a smoke-scoped path, and refuses to guess it.

`data/phenotype/` is the full cohort's directory and is empty. Stage 11 refuses
without phenotype data and stage 12 (GWAS) refuses below
`gwas.min_samples_per_group` per outcome group, so the first smoke cohort - ten
isolates, all imipenem-R - produced data that was entirely *correct* and still
could not be tested. That is stage 12 behaving correctly; the cohort composition
was what had to change. A test pins the floor so the next person does not
rebuild an untestable cohort and discover it at stage 12.

The calls are **measured**, read from the `AST phenotypes` column of
`PDC_essential.tsv`. Not assigned, not defaulted, not inferred - and the
refusal paths below exist because every one of those is a way to invent a call
the laboratory never made.

Two groups of tests, because the inputs are not all committed:

* `TestThePhenotypeDirectoryIsSmokeScoped` reads only committed config, so it
  always runs.
* `TestTheRealCohort` needs `PDC_essential.tsv` and the smoke assemblies, both
  of which are deliberately uncommitted (AGENTS.md rule 4 - databases and genome
  files are never committed). It skips when they are absent rather than
  pretending to pass, so a fresh clone is honest about what it has not checked.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
import yaml

from papipeline.config.loader import load_config
from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
MACHINES = REPO / "config" / "machines"
SMOKE_OVERLAY = MACHINES / "smoke.yaml"
PDC_TABLE = REPO / "PDC_essential.tsv"
SMOKE_GENOMES = REPO / "db" / "smoke_genomes"

BUILD = REPO / "scripts" / "smoke" / "build_smoke_phenotype.py"


def _load_builder():
    """Import the extractor by path, so the test does not need it on sys.path."""
    spec = importlib.util.spec_from_file_location("build_smoke_phenotype", BUILD)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


builder = _load_builder()


@pytest.fixture(scope="module")
def smoke_config():
    return load_config(SCIENCE, machine=SMOKE_OVERLAY)


@pytest.fixture
def pdc(tmp_path):
    """A synthetic PDC isolate table, for the refusal paths.

    The real one is uncommitted and 967 rows wide; testing refusals against it
    would mean editing it, which is exactly what must not happen.
    """

    def _write(rows, name="PDC_essential.tsv"):
        path = tmp_path / name
        lines = ["Isolate\tAssembly\t" + "AST phenotypes"]
        lines.extend(rows)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    return _write


@pytest.fixture
def cohort(tmp_path):
    def _make(*isolates):
        directory = tmp_path / "genomes"
        directory.mkdir(exist_ok=True)
        for isolate in isolates:
            (directory / f"{isolate}.fna").write_text(">x\nACGT\n", encoding="utf-8")
        return directory

    return _make


# --------------------------------------------------------------------------
class TestThePhenotypeDirectoryIsSmokeScoped:
    """Config only, so these run everywhere."""

    def test_a_smoke_overlay_has_its_own_phenotype_dir(self, smoke_config):
        directory = smoke_config.machine.smoke_phenotype_dir()
        assert directory != REPO / "data" / "phenotype", (
            "the smoke phenotype directory must not be data/phenotype/; that is "
            "the real cohort's directory, and a bounded run must not write into "
            "it - the table would then be a 10-isolate file sitting where the "
            "full cohort's belongs"
        )

    def test_phenotype_dir_returns_the_smoke_one(self, smoke_config):
        """Stage 11 must not be handed the full cohort's directory.

        `phenotype_dir` is what the stage actually calls, so this is the
        assertion that protects the run; the accessor test above covers only
        the accessor.
        """
        assert smoke_config.phenotype_dir(RunMode.REAL) == (
            smoke_config.machine.smoke_phenotype_dir()
        )

    def test_a_non_smoke_overlay_still_reads_the_real_directory(self):
        """Asserted positively, not as an absence.

        A laptop run must go on reading `data/phenotype/`, and the smoke branch
        must not leak into it.
        """
        config = load_config(SCIENCE, machine=MACHINES / "laptop.yaml")
        assert not config.is_smoke_overlay()
        assert config.phenotype_dir(RunMode.REAL) == REPO / "data" / "phenotype"

    def test_asking_a_non_smoke_overlay_is_refused(self):
        config = load_config(SCIENCE, machine=MACHINES / "laptop.yaml")
        with pytest.raises(Exception) as excinfo:
            config.machine.smoke_phenotype_dir()
        assert "smoke_phenotype_dir" in str(excinfo.value)

    def test_a_missing_key_is_refused_at_load_by_name(self, tmp_path):
        """No default, no fall-through - refused at *load*, like the genome key.

        A lazy check would let the overlay load and the run begin, then fail
        once something asked for the path - at which point it reads as a missing
        directory rather than a misconfigured overlay.
        """
        data = yaml.safe_load(SMOKE_OVERLAY.read_text(encoding="utf-8"))
        data["paths"].pop("smoke_phenotype_dir", None)
        path = tmp_path / "smoke_no_pheno.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        with pytest.raises(Exception) as excinfo:
            load_config(SCIENCE, machine=path)
        assert "smoke_phenotype_dir" in str(excinfo.value)

    def test_the_refusal_says_why_there_is_no_fallback(self, tmp_path):
        """"Add the key" is easy to act on; "we did not default it" is not,
        unless the refusal explains the reason."""
        data = yaml.safe_load(SMOKE_OVERLAY.read_text(encoding="utf-8"))
        data["paths"].pop("smoke_phenotype_dir", None)
        path = tmp_path / "smoke_no_pheno.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        with pytest.raises(Exception) as excinfo:
            load_config(SCIENCE, machine=path)
        assert "data/phenotype/" in str(excinfo.value)


# --------------------------------------------------------------------------
class TestTheExtractorRefusesToGuess:
    """Every way this script could invent a call, closed.

    Each of these is a fabrication dressed as convenience: an omitted isolate
    silently shrinks an outcome group until stage 12 refuses for a reason that
    looks like a pipeline fault, and an `I` folded into R or S is a call the
    laboratory never made.
    """

    def test_it_names_the_missing_isolate_and_writes_nothing(
        self, pdc, cohort, tmp_path, capsys
    ):
        """The message has to name the sample, or rule 5's diagnosis is lost."""
        source = pdc(["PDT_A\tGCA_1\t\"imipenem=R\""])
        genomes = cohort("PDT_A", "PDT_GHOST")
        out = tmp_path / "out"
        assert builder.main(
            ["--pdc", str(source), "--genomes", str(genomes), "--out", str(out)]
        ) == 2
        assert "PDT_GHOST" in capsys.readouterr().err
        assert not (out / "imipenem_phenotype.tsv").exists()

    def test_it_refuses_to_default_an_unmeasured_isolate(self, pdc, cohort, tmp_path):
        """No imipenem entry at all. Defaulting it to R or S would be a guess."""
        source = pdc(["PDT_A\tGCA_1\t\"amikacin=R\""])
        genomes = cohort("PDT_A")
        assert builder.main(
            ["--pdc", str(source), "--genomes", str(genomes),
             "--out", str(tmp_path / "out")]
        ) == 2

    def test_it_refuses_to_refold_an_intermediate(self, pdc, cohort, tmp_path):
        """`I` is a real AST category the model excludes; R or S would be a lie."""
        source = pdc(["PDT_A\tGCA_1\t\"imipenem=I\""])
        genomes = cohort("PDT_A")
        assert builder.main(
            ["--pdc", str(source), "--genomes", str(genomes),
             "--out", str(tmp_path / "out")]
        ) == 2

    def test_it_names_the_forbidden_directory(self, tmp_path):
        """Tested as a function, never by aiming a write at the real directory.

        Were the guard dropped, the alternative test would put a ten-isolate
        table into `data/phenotype/` - the corruption the guard prevents,
        caused by the test that checks it.
        """
        forbidden = tmp_path / "phenotype"
        assert builder.refuses_to_write_into(forbidden, forbidden) is not None

    def test_it_permits_the_smoke_directory(self, tmp_path):
        assert builder.refuses_to_write_into(
            tmp_path / "smoke_phenotypes", tmp_path / "phenotype"
        ) is None

    def test_the_refusal_survives_a_sloppy_path(self, tmp_path):
        """`resolve()` first, or `./phenotype/` and a trailing slash slip past."""
        forbidden = tmp_path / "phenotype"
        for sloppy in (tmp_path / "." / "phenotype", tmp_path / "phenotype" / "."):
            assert builder.refuses_to_write_into(sloppy, forbidden) is not None

    def test_it_refuses_an_empty_cohort(self, pdc, tmp_path):
        source = pdc(["PDT_A\tGCA_1\t\"imipenem=R\""])
        assert builder.main(
            ["--pdc", str(source), "--genomes", str(tmp_path / "nothing"),
             "--out", str(tmp_path / "out")]
        ) == 2


# --------------------------------------------------------------------------
class TestParsing:
    """The exported AST field is quoted and comma-separated; a split is not enough."""

    @pytest.mark.parametrize(
        "field, expected",
        [
            ('"imipenem=R"', "R"),
            ('"amikacin=R,imipenem=S"', "S"),
            ('"imipenem=S,amikacin=R"', "S"),
            ('"imipenem=I"', "I"),
            ('"amikacin=R"', None),
            ('"meropenem=R"', None),
            ("", None),
            ('"imipenem=R,ciprofloxacin=I"', "R"),
        ],
    )
    def test_the_call_is_read_correctly(self, field, expected):
        assert builder.imipenem_call(field) == expected

    def test_it_does_not_match_a_drug_that_merely_contains_the_name(self):
        assert builder.imipenem_call('"ertapenem=R"') is None


# --------------------------------------------------------------------------
needs_real_data = pytest.mark.skipif(
    not (PDC_TABLE.is_file() and SMOKE_GENOMES.is_dir()),
    reason="needs PDC_essential.tsv and db/smoke_genomes/, both uncommitted by "
           "AGENTS.md rule 4 - absent in a fresh clone",
)


@needs_real_data
class TestTheRealCohort:
    """The acceptance criterion: measured calls, both groups, floor cleared.

    Regenerated into a tmp directory rather than read from `db/`, so this
    exercises the script end to end and cannot pass on a stale hand-edit.
    """

    @pytest.fixture(scope="class")
    def built(self, tmp_path_factory):
        out = tmp_path_factory.mktemp("smoke_phenotype")
        code = builder.main(
            ["--genomes", str(SMOKE_GENOMES), "--out", str(out)]
        )
        assert code == 0
        return out / "imipenem_phenotype.tsv"

    def test_the_cohort_keeps_both_groups_and_clears_the_floor(
        self, built, smoke_config
    ):
        """The assertion that would have caught the original 10/10-R cohort."""
        from papipeline.io.tsv import read_tsv

        rows = read_tsv(
            built, required_columns=("sample_id", "antibiotic", "phenotype")
        )
        calls = [r["phenotype"] for r in rows]
        floor = smoke_config.gwas.min_samples_per_group
        for group in ("R", "S"):
            n = calls.count(group)
            assert n >= floor, (
                f"{n} isolates are {group}; GWAS needs >= {floor} per group. A "
                "single-group cohort cannot be tested however correct its data "
                "is - stage 12 refuses, and it is right to."
            )

    def test_stage_11_loads_it_and_stage_12_would_accept(
        self, built, smoke_config
    ):
        """End to end through the real loader, not a reimplementation of it."""
        from papipeline.stages.phenotype import load_phenotype

        records = load_phenotype(smoke_config, built.parent, "imipenem")
        assert len(records) == 10
        calls = [str(r.phenotype) for r in records]
        assert set(calls) == {"R", "S"}
        floor = smoke_config.gwas.min_samples_per_group
        assert min(calls.count("R"), calls.count("S")) >= floor

    def test_it_covers_exactly_the_cohort(self, built):
        """Both directions matter: a row for an absent sample inflates the join,
        a missing one drops the sample silently."""
        from papipeline.io.tsv import read_tsv

        rows = read_tsv(built, required_columns=("sample_id",))
        assert {r["sample_id"] for r in rows} == {
            p.stem for p in SMOKE_GENOMES.glob("*.fna")
        }

    def test_every_call_is_binary(self, built):
        from papipeline.io.tsv import read_tsv

        rows = read_tsv(built, required_columns=("phenotype",))
        assert {r["phenotype"] for r in rows} <= {"R", "S"}

    def test_every_row_names_imipenem(self, built):
        """The filename is the only thing tying the table to a drug, so a row
        for another antibiotic in this file would be read as an imipenem call."""
        from papipeline.io.tsv import read_tsv

        rows = read_tsv(built, required_columns=("antibiotic",))
        assert {r["antibiotic"] for r in rows} == {"imipenem"}

    def test_no_mic_is_invented(self, built):
        """The AST gives a categorical call only. Converting R/S into a number
        would fabricate a measurement, and the config forbids exactly that."""
        from papipeline.io.tsv import read_tsv

        rows = read_tsv(built, required_columns=("phenotype", "MIC"))
        assert all(r["MIC"] is None for r in rows)

    def test_the_provenance_is_recorded(self, built):
        """A phenotype table with no stated source is indistinguishable from one
        a person typed."""
        text = built.read_text(encoding="utf-8")
        assert "PDC_essential.tsv" in text
        assert "AST" in text

    def test_it_does_not_look_like_the_full_cohort_table(self, built):
        """It must not be mistakable for `data/phenotype/`, which is empty."""
        assert "data/phenotype/" in built.read_text(encoding="utf-8")
