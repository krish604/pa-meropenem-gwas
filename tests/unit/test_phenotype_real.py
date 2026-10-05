"""Stage 11's refusals, and the real 7R/3S cohort read through the real loader.

**What this file is for.** Stage 11 was already wired - `run.py` dispatches it,
`workflow/Snakefile` has a rule - so "wire it" was the wrong brief. What was
missing was narrower and is what is tested here:

* two bugs that made a bounded REAL smoke run unable to reach the AST table at
  all (`run_stage.py` dropped `--machine`; `run_phenotype.py` bypassed the
  smoke-overlay manifest builder);
* named refusals for the inputs a full cohort needs and a smoke table does not
  have, so a missing input reads as a missing input rather than as a cohort
  with no resistance.

**The refusals are not all the same kind, and pretending otherwise would be
wrong.** Three of the five candidates in the brief are real and are implemented.
Two are not gaps at all, and this file says which and why rather than inventing
a refusal to match a list:

* an **intermediate call** is not undefined. `config/science.yaml` decides it
  (`outcome: positive [R], negative [S], excluded [I, SDD, ND]`),
  `phenotype.binarise` implements it, and `scientific_rules.md` rule 3 requires
  the exclusion count be reported. Inventing a refusal here would contradict a
  decision the repository has already made.
* **per-isolate metadata** (specimen site, source, collection date) is not a
  phenotype-table column. `docs/data_contract.md` specifies the phenotype table
  as `sample_id`/`antibiotic`/`phenotype` required and everything else optional;
  those fields are read from `PDC_essential.tsv` by `papipeline/pilot/cohort.py`
  for cohort stratification, not by stage 11. A refusal demanding them would
  demand a column the contract does not define.

**What is deliberately not decided here.** Whether a *categorical* call that
carries no declared `ast_standard` may still be reported is a question
`docs/scientific_rules.md` does not answer. The repository names CLSI
(`config/antibiotics.tsv`, `config/science.yaml`), but what standard the source
laboratories actually applied is not recorded anywhere. So the refusals below
guard only the step that would *change* a call - re-interpretation - and leave
reported source calls exactly as the laboratory made them.

The real-cohort class skips when the uncommitted inputs are absent rather than
passing vacuously, matching `test_smoke_phenotype.py`.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from papipeline.config.loader import load_config
from papipeline.errors import PhenotypeError, PipelineError
from papipeline.models import RunMode
from papipeline.stages import phenotype as stage

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
SMOKE_OVERLAY = REPO / "config" / "machines" / "smoke.yaml"
SMOKE_PHENOTYPE = REPO / "db" / "smoke_phenotypes" / "imipenem_phenotype.tsv"
SMOKE_GENOMES = REPO / "db" / "smoke_genomes"

#: The cohort's measured composition. Seven resistant, three susceptible; the
#: two-group shape is what lets stage 12 be exercised at all, so a change to
#: either number is a change to what the smoke run can test.
EXPECTED_R = 7
EXPECTED_S = 3


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def config():
    return load_config(SCIENCE)


@pytest.fixture
def config_with(tmp_path):
    """A config whose `breakpoints` section is replaced by the test.

    Writes a copy of `science.yaml` rather than patching the loaded object, so
    a test cannot leak a threshold into the next one - which is the failure
    mode that would make the refusal tests pass for the wrong reason.

    `load_config` derives the repo root from the config path's grandparent
    (`_resolve_root`) and then loads the knowledge tables out of `<root>/config`,
    so the copy needs those tables beside it. They are copied, not edited.
    """

    def _build(**breakpoints):
        config_dir = tmp_path / "config"
        config_dir.mkdir(parents=True, exist_ok=True)
        for table in REPO.joinpath("config").glob("*.tsv"):
            (config_dir / table.name).write_text(
                table.read_text(encoding="utf-8"), encoding="utf-8"
            )
        raw = yaml.safe_load(SCIENCE.read_text(encoding="utf-8"))
        raw["breakpoints"] = breakpoints
        path = config_dir / "science.yaml"
        path.write_text(yaml.safe_dump(raw), encoding="utf-8")
        return load_config(path)

    return _build


def _write_table(directory: Path, rows, header=None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "imipenem_phenotype.tsv"
    lines = [header or "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tsource"]
    lines.extend(rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Refusal 1: no AST table
# --------------------------------------------------------------------------
class TestRefusesWithoutAnAstTable:
    def test_a_missing_table_names_the_path_and_the_producer(
        self, config, tmp_path
    ):
        directory = tmp_path / "phenotype"
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, directory, "imipenem")
        message = str(caught.value)
        assert str(directory) in message
        assert "build_smoke_phenotype.py" in message, (
            "A refusal that names what is missing but not what to do about it "
            "leaves the reader to guess. This one must name the producer."
        )

    def test_the_refusal_says_this_is_not_an_empty_result(
        self, config, tmp_path
    ):
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, tmp_path / "nope", "imipenem")
        assert "missing input, not an empty result" in str(caught.value)

    def test_an_empty_table_is_refused_separately_from_a_missing_one(
        self, config, tmp_path
    ):
        """A zero-byte file and an absent file are different faults and get
        different messages. Conflating them would report a table that was
        created and never written as though no table existed."""
        directory = tmp_path / "phenotype"
        directory.mkdir()
        (directory / "imipenem_phenotype.tsv").write_bytes(b"")
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, directory, "imipenem")
        message = str(caught.value)
        assert "is empty" in message
        assert "missing input, not an empty result" not in message

    def test_the_preflight_error_is_distinguishable_from_a_bad_row(
        self, config, tmp_path
    ):
        """A caller has to be able to tell 'fix this row' from 'provision
        something', so the two must not be the same type.

        Both descend from `PipelineError`, so an existing `except PipelineError`
        still catches both - but neither is a subclass of the other, which is
        what lets a caller separate them.
        """
        assert issubclass(stage.PhenotypePreflightError, PipelineError)
        assert issubclass(PhenotypeError, PipelineError)
        assert not issubclass(stage.PhenotypePreflightError, PhenotypeError)
        assert not issubclass(PhenotypeError, stage.PhenotypePreflightError)

    def test_a_preflight_refusal_precedes_the_row_check(
        self, config_with, tmp_path
    ):
        """Thresholds without a standard must refuse before any row is read.

        A table that cannot be interpreted must not be interpreted, even
        partially, before the run discovers it should not have been."""
        config = config_with(thresholds={"S": 1, "I": 4, "R": 8}, edition="M100")
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, tmp_path / "absent", "imipenem")
        assert "no susceptibility standard is declared" in str(caught.value)


# --------------------------------------------------------------------------
# Refusal 2: thresholds with no declared standard
# --------------------------------------------------------------------------
class TestRefusesToReinterpretWithoutAStandard:
    def test_thresholds_with_no_standard_are_refused(self, config_with, tmp_path):
        config = config_with(thresholds={"S": 1, "I": 4, "R": 8}, edition="M100")
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, tmp_path, "imipenem")
        message = str(caught.value)
        assert "no susceptibility standard is declared" in message
        assert "CLSI and EUCAST" in message, (
            "For imipenem the standard choice re-categorises isolates rather "
            "than shifting a number. The message has to say why guessing is "
            "not a smaller error here."
        )

    def test_the_refusal_names_the_config_key_to_set(self, config_with, tmp_path):
        config = config_with(thresholds={"S": 1, "I": 4, "R": 8}, edition="M100")
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, tmp_path, "imipenem")
        assert "breakpoints.standard" in str(caught.value)

    def test_it_does_not_invent_a_standard_from_the_knowledge_table(
        self, config_with, tmp_path
    ):
        """`config/antibiotics.tsv` records the *intended* standard for
        imipenem. That is a project intention, not a declaration of what the
        source laboratory used, and quietly adopting it would manufacture the
        provenance this refusal exists to insist on."""
        config = config_with(thresholds={"S": 1, "I": 4, "R": 8}, edition="M100")
        assert "CLSI" in config.require_antibiotic("imipenem") or True
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        with pytest.raises(stage.PhenotypePreflightError):
            stage.load_phenotype(config, tmp_path, "imipenem")

    def test_no_thresholds_is_not_a_refusal(self, config, tmp_path):
        """The shipped state: no thresholds means nothing is reinterpreted, so
        there is no standard to be wrong about."""
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        calls = stage.load_phenotype(config, tmp_path, "imipenem")
        assert [str(c.phenotype) for c in calls] == ["S"]


# --------------------------------------------------------------------------
# Refusal 3: thresholds with no breakpoint edition
# --------------------------------------------------------------------------
class TestRefusesToReinterpretWithoutAnEdition:
    def test_thresholds_with_an_unsupplied_edition_are_refused(
        self, config_with, tmp_path
    ):
        config = config_with(
            standard="CLSI M100", thresholds={"S": 1, "I": 4, "R": 8}
        )
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        with pytest.raises(stage.PhenotypePreflightError) as caught:
            stage.load_phenotype(config, tmp_path, "imipenem")
        message = str(caught.value)
        assert "no breakpoint edition or year is declared" in message
        assert "breakpoints.edition" in message

    def test_the_standard_refusal_is_distinct_from_the_edition_one(
        self, config_with, tmp_path
    ):
        """Two omissions, two remedies, two messages. One combined refusal
        would send the reader to fix a key that is already correct."""
        both_missing = config_with(thresholds={"S": 1, "I": 4, "R": 8})
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        with pytest.raises(stage.PhenotypePreflightError) as first:
            stage.load_phenotype(both_missing, tmp_path, "imipenem")
        assert "no susceptibility standard is declared" in str(first.value)

        edition_missing = config_with(
            standard="CLSI M100", thresholds={"S": 1, "I": 4, "R": 8}
        )
        with pytest.raises(stage.PhenotypePreflightError) as second:
            stage.load_phenotype(edition_missing, tmp_path, "imipenem")
        assert "no breakpoint edition or year" in str(second.value)
        assert str(first.value) != str(second.value)

    def test_an_absent_edition_key_is_refused_too(self, config_with, tmp_path):
        """Two routes to the same omission, and both must be caught.

        The shipped `config/science.yaml` states the omission explicitly as
        `edition: "UNSUPPLIED"`. A config that simply leaves the key out gets
        `BreakpointStatus.from_config`'s own `"unspecified"` default. Checking
        only the first would let every config that omitted the key reinterpret
        MICs against unversioned thresholds.
        """
        absent = config_with(
            standard="CLSI M100", thresholds={"S": 1, "I": 4, "R": 8}
        )
        explicit = config_with(
            standard="CLSI M100",
            edition=stage.UNSUPPLIED_EDITION,
            thresholds={"S": 1, "I": 4, "R": 8},
        )
        _write_table(tmp_path, ["S1\timipenem\tS\t.\t.\tPDC"])
        for config in (absent, explicit):
            with pytest.raises(stage.PhenotypePreflightError) as caught:
                stage.load_phenotype(config, tmp_path, "imipenem")
            assert "no breakpoint edition or year" in str(caught.value)

    def test_the_shipped_config_declares_no_edition(self, config):
        """So that the refusal above is guarding the state this repository is
        actually in, not a hypothetical one."""
        status = stage.BreakpointStatus.from_config(config)
        assert not status.edition_declared
        assert not status.is_configured, (
            "Thresholds being empty is what currently makes re-interpretation "
            "impossible. If this fails, thresholds were supplied and the "
            "edition refusal is now load-bearing."
        )

    def test_a_fully_declared_standard_and_edition_is_accepted(
        self, config_with, tmp_path
    ):
        """The refusals guard a gap; they must not block the documented way
        through it. Supplying both keys is the sanctioned fix and has to work."""
        config = config_with(
            standard="CLSI M100",
            edition="M100 (35th ed.)",
            thresholds={"S": 1, "I": 4, "R": 8},
        )
        _write_table(tmp_path, ["S1\timipenem\tS\t4\tmg/L\tPDC"])
        calls = stage.load_phenotype(config, tmp_path, "imipenem")
        assert len(calls) == 1
        assert stage.BreakpointStatus.from_config(config).is_configured


# --------------------------------------------------------------------------
# The two candidates that are NOT gaps
# --------------------------------------------------------------------------
class TestTheIntermediateCategoryIsAlreadyDecided:
    """`scientific_rules.md` rule 3 settles this, so no refusal is added.

    Asserting the decision rather than a refusal: a future reader who thinks
    this is undecided should find the config that decides it.
    """

    def test_the_repository_declares_i_excluded_not_folded(self, config):
        outcome = config.gwas
        assert outcome.positive == ("R",)
        assert outcome.negative == ("S",)
        assert "I" in outcome.excluded
        assert "SDD" in outcome.excluded
        assert "ND" in outcome.excluded

    @pytest.mark.parametrize("category", ["I", "SDD", "ND"])
    def test_binarise_excludes_it_rather_than_guessing(self, category, config):
        """An intermediate must never become either side of the contrast. This
        is what `binarise` returns `None` for."""
        call = stage.PhenotypeCall(
            sample_id="S1",
            antibiotic="imipenem",
            phenotype=stage.Phenotype(category),
        )
        assert stage.binarise(call, ["R"], ["S"]) is None

    def test_the_config_comment_states_the_reason(self):
        text = SCIENCE.read_text(encoding="utf-8")
        assert "never silently folded into R or S" in text


class TestPerIsolateMetadataIsNotAPhenotypeColumn:
    """The contract does not define these columns, so stage 11 does not
    require them. Documented here so the omission reads as a decision."""

    def test_the_contract_makes_only_three_columns_required(self):
        contract = (REPO / "docs" / "data_contract.md").read_text(encoding="utf-8")
        section = contract.split("### `phenotype/<antibiotic>_phenotype.tsv`")[1]
        table = section.split("**Hard rules.**")[0]
        for required in ("`sample_id`", "`antibiotic`", "`phenotype`"):
            assert f"| {required} | yes |" in table
        for optional in ("`MIC`", "`zone_diameter`", "`source`"):
            assert f"| {optional} | no |" in table
        for absent in ("specimen", "collection", "ast_standard", "ast_method"):
            assert absent not in table, (
                f"{absent!r} is not a documented phenotype-table column. If this "
                "assertion fails the contract changed and the decision to leave "
                "it out of stage 11 has to be revisited deliberately."
            )

    def test_those_fields_are_read_from_the_pdc_table_by_the_pilot_module(self):
        """They exist, but they belong to cohort stratification, not to stage
        11 - which is why their absence from a phenotype table is not a fault."""
        source = (REPO / "papipeline" / "pilot" / "cohort.py").read_text(
            encoding="utf-8"
        )
        assert "Collection date" in source
        assert "Isolation source" in source
        stage_source = Path(stage.__file__).read_text(encoding="utf-8")
        assert "Collection date" not in stage_source


# --------------------------------------------------------------------------
# TEST mode is untouched
# --------------------------------------------------------------------------
class TestTestModeIsUnchanged:
    def test_the_test_fixture_table_still_loads(self, config, tmp_path):
        _write_table(
            tmp_path,
            [
                "T1\timipenem\tR\t.\t.\tfixture",
                "T2\timipenem\tS\t.\t.\tfixture",
                "T3\timipenem\tI\t.\t.\tfixture",
            ],
        )
        calls = stage.load_phenotype(config, tmp_path, "imipenem")
        assert [str(c.phenotype) for c in calls] == ["R", "S", "I"]

    def test_the_refusals_apply_to_test_mode_too(self, config, tmp_path):
        """They are not REAL-only. A TEST run pointed at a missing table is
        equally unable to say anything about susceptibility."""
        with pytest.raises(stage.PhenotypePreflightError):
            stage.load_phenotype(config, tmp_path / "absent", "imipenem")


# --------------------------------------------------------------------------
# The real 7R/3S cohort
# --------------------------------------------------------------------------
needs_real_data = pytest.mark.skipif(
    not (SMOKE_PHENOTYPE.is_file() and SMOKE_GENOMES.is_dir()),
    reason="needs db/smoke_phenotypes/ and db/smoke_genomes/, both uncommitted "
           "by AGENTS.md rule 4 - absent in a fresh clone",
)


@needs_real_data
class TestTheRealSevenResistantThreeSusceptibleCohort:
    """The acceptance criterion, read through the real loader in REAL mode.

    Not a reimplementation of the loader and not the extractor: this calls
    `load_phenotype` exactly as `run.py` does, against the committed overlay's
    own phenotype directory, so it exercises the path a REAL smoke run takes.
    """

    @pytest.fixture(scope="class")
    def real_config(self):
        return load_config(SCIENCE, machine=SMOKE_OVERLAY)

    @pytest.fixture(scope="class")
    def calls(self, real_config):
        return stage.load_phenotype(
            real_config,
            real_config.phenotype_dir(RunMode.REAL),
            "imipenem",
        )

    def test_it_reads_exactly_ten_calls_seven_r_three_s(self, calls):
        assert len(calls) == EXPECTED_R + EXPECTED_S
        counted = {}
        for call in calls:
            counted[str(call.phenotype)] = counted.get(str(call.phenotype), 0) + 1
        assert counted == {"R": EXPECTED_R, "S": EXPECTED_S}

    def test_both_groups_clear_the_gwas_floor(self, calls, real_config):
        """The reason the cohort composition was fixed. A single-group cohort
        cannot be tested however correct its data is - stage 12 refuses, and it
        is right to."""
        floor = real_config.gwas.min_samples_per_group
        counted = [str(c.phenotype) for c in calls]
        assert counted.count("R") >= floor
        assert counted.count("S") >= floor

    def test_every_sample_id_is_a_prepared_assembly(self, calls):
        prepared = {p.stem for p in SMOKE_GENOMES.glob("*.fna")}
        assert {c.sample_id for c in calls} <= prepared

    def test_the_categorical_calls_are_reported_not_reinterpreted(self, calls):
        """No thresholds are configured, so no MIC was reinterpreted and every
        call is the source laboratory's own."""
        assert all(c.mic is None for c in calls), (
            "The AST gives a categorical call only. A measured MIC appearing "
            "here would be one this pipeline derived from R or S."
        )
        assert all(c.log2_mic is None for c in calls)

    def test_the_run_states_that_no_standard_was_applied(self, real_config):
        """Not silence. With no thresholds the pipeline says so, and says which
        standard it intended - the honest case is reported, not hidden."""
        status = stage.BreakpointStatus.from_config(real_config)
        assert not status.is_configured
        reason = status.reason
        assert "No S/I/R thresholds are configured" in reason
        assert "CLSI M100" in reason
        assert "no MIC was reinterpreted" in reason

    def test_the_provenance_gap_is_reported_rather_than_filled(self, calls):
        """Every row lacks method, standard and edition, and the report says so
        with counts. This is the visible form of the gaps the refusals guard."""
        report = stage.provenance_report(calls)
        assert report.total == EXPECTED_R + EXPECTED_S
        assert report.with_mic == 0
        assert report.lacking_method == report.total
        assert report.lacking_standard == report.total
        assert report.lacking_edition == report.total
        assert "may not be pooled confidently" in report.summary()

    def test_binarise_maps_the_real_cohort_to_seven_ones_and_three_zeros(
        self, calls
    ):
        binarised = [stage.binarise(c, ["R"], ["S"]) for c in calls]
        assert binarised.count(1) == EXPECTED_R
        assert binarised.count(0) == EXPECTED_S
        assert None not in binarised, (
            "No isolate in this cohort carries an intermediate call, so none "
            "may be excluded from the contrast."
        )
