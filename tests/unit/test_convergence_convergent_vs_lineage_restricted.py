"""A mutation planted in three lineages IS convergent; one confined to a single
lineage is NOT. Both directions, on real code, with synthetic inputs.

**Why this file exists.** Every other test of `stages.convergence` pins a
*guard*: that the stage refuses on an absent input, that it refuses on a
degenerate crosswalk, that a single-lineage determinant is `lineage_associated`.
None of them pins the answer the stage exists to give, on a cohort whose answer
was planted rather than hoped for. A stage that refused correctly and computed
nonsense would pass all of them.

So this file runs the real `run()` - the REAL path, reading three TSVs off
disk - over `tests/fixtures/convergence_real/`, a nine-sample cohort in three
lineages of three, and asserts the calls against a table derived from the
fixture's construction rather than from a previous run's output.

**R12: real code path, synthetic inputs.** Nothing here is mocked. `run()` is
the production entry point, `_run_real` is the production REAL branch,
`load_amr_table` and `load_tree_metadata` are the production readers, and the
only thing substituted is the data. That is the point - the earlier rounds'
complaint about this stage was that its guards were well tested and its
computation was not.

**Why the denominator is the whole fixture.** `convergence.widespread_fraction`
is 0.75. Over these nine samples `widespread_background` needs seven carriers,
which `blaOXA-1` (eight) reaches. Over the 967-isolate PDC roster - which is
what `len(manifest)` returns, because a manifest is the roster whatever subset
this run analysed - it would need 726, and the category would be unreachable.
This fixture therefore fails if the `widespread_fraction` denominator is ever
reverted to the manifest size, which no assertion on a refusal path could have
caught.

**Not a claim about the real cohort.** Nine synthetic samples in three lineages
is a fixture. Nothing here is evidence about *Pseudomonas aeruginosa*, and the
counts asserted below are properties of a constructed table.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.errors import StageError
from papipeline.stages import convergence as conv

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "convergence_real"

#: The planted cohort, restated here so the expected answers are visible in the
#: assertion rather than only in `fixtures/convergence_real/README.md`. If the
#: fixture changes, this is the list that has to change with it.
EXPECTED = {
    # determinant: (n_carriers, n_independent_lineages, category)
    "blaOXA-1": (8, 3, "widespread_background"),
    "fosA": (3, 1, "lineage_associated"),
    "oprD:A100T": (3, 3, "recurrent_convergent"),
    "mexR:I24A": (3, 1, "lineage_associated"),
}


@pytest.fixture(scope="module")
def config():
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture
def cohort(tmp_path: Path):
    """The committed fixture, copied into ``tmp_path`` and laid out on disk in
    the contracted positions.

    Copied rather than read in place so that nothing this test does - and
    nothing a future edit to this file does - can write into the repository or
    into ``results/``. The first version of the sibling test file wrote a
    ``tree_metadata.tsv`` into ``results/``, and a later "does this refuse?"
    check then read that leftover and concluded the lineages input was present.
    A refusal test that passes because of a stale artefact is worse than no
    refusal test.
    """
    intermediate = tmp_path / "intermediate"
    phylogeny = tmp_path / "phylogeny"
    (intermediate / "stages").mkdir(parents=True)
    (intermediate / "regulators").mkdir(parents=True)
    phylogeny.mkdir(parents=True)
    shutil.copyfile(FIXTURES / "04_amr.tsv", intermediate / "stages" / "04_amr.tsv")
    shutil.copyfile(
        FIXTURES / "regulator_variants.tsv",
        intermediate / "regulators" / "regulator_variants.tsv",
    )
    shutil.copyfile(
        FIXTURES / "tree_metadata.tsv", phylogeny / "tree_metadata.tsv"
    )
    return intermediate, phylogeny


#: The manifest is deliberately **967 isolates** - the whole PDC roster - while
#: only nine are analysed. That is not decoration: it is how a smoke run looks,
#: and it is the only manifest shape under which the `len(manifest)` fallback is
#: visible. A manifest sized to the fixture would make every assertion in this
#: file pass with the defect in place.
ROSTER = tuple(f"PDT{i:09d}.1" for i in range(967))

#: The nine samples the fixture cohort actually contains, crosswalk included.
ANALYSED = ("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")


@pytest.fixture(scope="module")
def manifest() -> SampleManifest:
    """The whole roster, not the nine. See `ROSTER`."""
    return SampleManifest(samples=[Sample(sample_id=s) for s in ROSTER])


@pytest.fixture
def calls(config, manifest, cohort):
    intermediate, phylogeny = cohort
    return conv.run(
        config,
        manifest,
        RunMode.REAL,
        intermediate_root=intermediate,
        antibiotic="imipenem",
        phylogeny_dir=phylogeny,
    )


class TestThePlantedAnswersAreTheCalls:
    def test_every_planted_determinant_appears_exactly_once(self, calls):
        names = [c.determinant for c in calls]
        assert sorted(names) == sorted(EXPECTED), (
            f"the fixture's four determinants did not survive the read: {names}"
        )

    @pytest.mark.parametrize("determinant", sorted(EXPECTED))
    def test_carrier_count_and_lineage_count(self, calls, determinant):
        want_n, want_lineages, _ = EXPECTED[determinant]
        call = next(c for c in calls if c.determinant == determinant)
        assert call.branch_count == want_n, (
            f"{determinant}: read {call.branch_count} carriers, planted {want_n}"
        )
        assert call.independent_lineages == want_lineages, (
            f"{determinant}: read {call.independent_lineages} lineages, "
            f"planted {want_lineages}"
        )

    @pytest.mark.parametrize("determinant", sorted(EXPECTED))
    def test_category(self, calls, determinant):
        _, _, want = EXPECTED[determinant]
        call = next(c for c in calls if c.determinant == determinant)
        assert call.convergence_category.value == want

    def test_the_distribution_column_shows_the_planting(self, calls):
        """The evidence for the call, not just the call.

        A category with no distribution beside it is a category a reader has to
        take on trust, and the module docstring promises the distribution is
        published for exactly that reason.
        """
        by_name = {c.determinant: c for c in calls}
        assert by_name["oprD:A100T"].distribution == {"L1": 1, "L2": 1, "L3": 1}
        assert by_name["mexR:I24A"].distribution == {"L1": 3}
        assert by_name["fosA"].distribution == {"L2": 3}


class TestTheConvergenceCallIsTheMutationAndNotTheGene:
    """The headline requirement, stated as its own class.

    `oprD:A100T` and `mexR:I24A` have the *same* three carriers' worth of
    rarity and differ only in which lineages those carriers sit in. They cannot
    both be right unless the stage read the crosswalk; they cannot both be
    wrong in the same direction unless the stage ignored it.
    """

    def test_a_mutation_in_three_lineages_is_convergent(self, calls):
        call = next(c for c in calls if c.determinant == "oprD:A100T")
        assert call.independent_lineages == 3, (
            "oprD:A100T was planted in A1(L1), B1(L2) and C1(L3) - three "
            "independent lineages"
        )
        assert call.convergence_category.value == "recurrent_convergent"

    def test_a_mutation_confined_to_one_lineage_is_not(self, calls):
        call = next(c for c in calls if c.determinant == "mexR:I24A")
        assert call.independent_lineages == 1, (
            "mexR:I24A was planted in A1, A2 and A3 - all L1"
        )
        assert call.convergence_category.value != "recurrent_convergent"
        assert call.convergence_category.value == "lineage_associated"

    def test_same_carrier_count_opposite_answers(self, calls):
        """The pair that makes a coincidence impossible."""
        by_name = {c.determinant: c for c in calls}
        convergent = by_name["oprD:A100T"]
        restricted = by_name["mexR:I24A"]
        assert convergent.branch_count == restricted.branch_count == 3
        assert convergent.independent_lineages == 3
        assert restricted.independent_lineages == 1
        assert convergent.convergence_category is not restricted.convergence_category

    def test_convergent_calls_returns_only_the_convergent_one(self, calls):
        found = [c.determinant for c in conv.convergent_calls(calls)]
        assert found == ["oprD:A100T"]


class TestTheDenominatorIsTheAnalysedCohort:
    """`widespread_fraction` is 0.75 of the samples *analysed*.

    This is the assertion that the `n_samples or len(manifest)` fallback could
    not survive. `tests/unit/seam_call_site_defaulted_allowlist.txt` records
    that fallback as a FINDING against the run.py call site, on the grounds that
    a manifest is the whole PDC roster - 967 isolates - whatever subset was
    analysed. Eight carriers out of nine is `widespread_background`; eight out of
    967 is not. The manifest this file builds IS the 967 roster, so every
    assertion above runs at the shape where the defect is visible - a manifest
    sized to the fixture would let the defect back in silently.
    """

    def test_eight_of_nine_is_widespread_background(self, calls):
        call = next(c for c in calls if c.determinant == "blaOXA-1")
        assert call.branch_count == 8
        assert call.convergence_category.value == "widespread_background", (
            "8/9 = 0.889 exceeds widespread_fraction 0.75, so this determinant "
            "must read as background; it does not, so the denominator is not "
            "the analysed cohort"
        )

    def test_the_old_fallback_inverted_this_call(self, config, calls):
        """What the `n_samples or len(manifest)` fallback actually did.

        Not "the category was rarely reached", which is how the allowlist put
        it. It **inverted** the call: eight carriers in nine samples is a
        background determinant, and at `widespread_fraction` 0.75 with 967 in
        the denominator it needs 726 carriers to read as background, so it fell
        through to the lineage test and was published as
        `recurrent_convergent` - a convergent-selection candidate for a
        determinant three quarters of the analysed cohort carries. That is a
        claim about convergent selection, made on arithmetic, and it is exactly
        the over-claim `widespread_background` exists to prevent.
        """
        distribution = next(
            c.distribution for c in calls if c.determinant == "blaOXA-1"
        )
        record = conv.DeterminantCarriers(
            "blaOXA-1", tuple(f"s{i}" for i in range(8))
        )
        assert conv.classify(record, distribution, config, n_samples=967).value == (
            "recurrent_convergent"
        )
        assert conv.classify(record, distribution, config, n_samples=9).value == (
            "widespread_background"
        )

    def test_the_roster_size_cannot_move_a_call(self, config, cohort, calls):
        """Manifest-independence, which is the property the fix buys.

        A manifest is the whole PDC roster - 967 isolates - whatever subset the
        run analyses, precisely so an unprepared isolate stays visible and can
        refuse. So a stage whose denominator came from the manifest inherited
        the roster's size in its arithmetic. Reading the denominator off the
        crosswalk instead makes the answer independent of it: the same three
        tables produce the same calls under a 967-sample and a 9-sample
        manifest, field for field.

        `calls` is already computed under the roster (`ROSTER`), so the first
        half of that comparison is what the rest of this file asserts. This test
        runs the nine-sample case and requires them to agree.
        """
        intermediate, phylogeny = cohort
        nine = conv.run(
            config,
            SampleManifest(samples=[Sample(sample_id=s) for s in ANALYSED]),
            RunMode.REAL,
            intermediate_root=intermediate,
            antibiotic="imipenem",
            phylogeny_dir=phylogeny,
        )
        assert len(ROSTER) == 967
        assert [
            (c.determinant, c.branch_count, c.independent_lineages,
             c.convergence_category, c.distribution)
            for c in nine
        ] == [
            (c.determinant, c.branch_count, c.independent_lineages,
             c.convergence_category, c.distribution)
            for c in calls
        ]

    def test_an_explicit_n_samples_from_the_caller_is_honoured(self, config, cohort):
        """The parameter still works; it is the *fallback* that was wrong.

        A caller that knows the analysed count better than the crosswalk does -
        a Snakemake rule handed a cohort subset, say - must be able to say so.
        """
        intermediate, phylogeny = cohort
        calls = conv.run(
            config, SampleManifest(samples=[]), RunMode.REAL,
            n_samples=4,
            intermediate_root=intermediate,
            antibiotic="imipenem",
            phylogeny_dir=phylogeny,
        )
        by_name = {c.determinant: c for c in calls}
        # 8 carriers over a stated 4: over 1, so background by a wide margin.
        assert by_name["blaOXA-1"].convergence_category.value == "widespread_background"

    def test_the_analysed_count_is_the_crosswalk_not_the_manifest(self, config):
        assert conv.analysed_sample_count(
            {s: f"L{i}" for i, s in enumerate("ABCDEFGHI")}
        ) == 9

    def test_an_unclassifiable_sample_is_excluded_from_the_denominator(self):
        """A sample whose lineage is the sentinel cannot be classified, so it
        is not analysed. Counting it would deflate every carrier fraction and
        push calls AWAY from `widespread_background` - the same over-claim in
        the opposite direction.
        """
        assert conv.analysed_sample_count({"A1": "L1", "A2": "unknown"}) == 1

    def test_a_detergent_crosswalk_yields_zero_rather_than_the_manifest(self):
        assert conv.analysed_sample_count({}) == 0
        assert conv.analysed_sample_count({"A1": "unknown"}) == 0
        # In REAL that state is already a refusal, so zero never becomes a
        # division; assert the guard that makes that true is still in place.
        assert conv.UNKNOWN_LINEAGE == "unknown"


class TestTheRefusalStaysFailClosedAndActionable:
    """C3: each missing input is named, with the path it was expected at AND the
    stage that would produce it.

    `tests/unit/test_convergence_real_inputs.py` establishes that the stage
    refuses *by input name*. This class is the stronger claim, and the one a
    reader acts on: a refusal that names the input but not the path leaves them
    to guess where to look, and one that names neither the path nor the producer
    leaves them with a wall. Both halves are asserted per input rather than as a
    substring sweep, so losing one input's path fails one test that says which
    input lost it.
    """

    #: input name -> (its contracted path, the producer stage named for it)
    EXPECTED_REFUSAL = {
        "amr_calls": ("04_amr.tsv", "stage 4"),
        "regulator_variants": (
            "regulator_variants.tsv", "produce_regulator_variants",
        ),
        "lineages": ("tree_metadata.tsv", "build_tree_outputs"),
    }

    #: input name -> (fixture file, contracted path relative to the cohort root)
    SOURCES = {
        "amr_calls": (
            "04_amr.tsv", lambda i, p: i / "stages" / "04_amr.tsv"
        ),
        "regulator_variants": (
            "regulator_variants.tsv",
            lambda i, p: i / "regulators" / "regulator_variants.tsv",
        ),
        "lineages": (
            "tree_metadata.tsv", lambda i, p: p / "tree_metadata.tsv"
        ),
    }

    def _refuse(self, config, cohort, drop: str) -> str:
        """Every input present except ``drop``, which must be refused by name."""
        intermediate, phylogeny = cohort
        dropped = None
        for name, (filename, place) in self.SOURCES.items():
            target = place(intermediate, phylogeny)
            target.parent.mkdir(parents=True, exist_ok=True)
            if name == drop:
                # `cohort` writes all three; this is what makes it absent.
                target.unlink(missing_ok=True)
                dropped = target
                continue
            shutil.copyfile(FIXTURES / filename, target)
        assert dropped is not None and not dropped.exists()
        with pytest.raises(StageError) as excinfo:
            conv.run(
                config, self._manifest(), RunMode.REAL,
                intermediate_root=intermediate,
                antibiotic="imipenem",
                phylogeny_dir=phylogeny,
            )
        assert excinfo.value.context["missing"] == drop
        return str(excinfo.value)

    @staticmethod
    def _manifest() -> SampleManifest:
        return SampleManifest(samples=[Sample(sample_id=s) for s in ROSTER])

    @pytest.mark.parametrize("drop", sorted(EXPECTED_REFUSAL))
    def test_it_names_the_input_its_path_and_its_producer(
        self, config, cohort, drop
    ):
        message = self._refuse(config, cohort, drop)
        want_path, want_producer = self.EXPECTED_REFUSAL[drop]
        assert drop in message
        assert want_path in message, (
            f"the refusal for {drop} does not say which path was expected"
        )
        assert want_producer in message, (
            f"the refusal for {drop} does not say which stage would produce it"
        )

    def test_naming_a_dead_producer_is_not_left_in_the_message(self, config, cohort):
        """The corrected text must not still argue a closed gap.

        The previous message claimed stage 6 had no REAL producer and that
        `lineage_label` "has no REAL producer at all", naming
        `papipeline.testing.synthetic` as the only writer of each. Both were
        false: `stages.regulators.produce_regulator_variants` writes the first,
        called from `run.derive_regulator_table`, and
        `phylogeny.build_tree_outputs` writes the second. A reader who resolved
        every gap the message offered would have arrived at the same wall with
        no explanation left, which is the failure mode a refusal exists to
        prevent. Asserting the absence of the dead claims, so the correction
        cannot be silently reverted.
        """
        message = self._refuse(config, cohort, "regulator_variants")
        assert "only papipeline.testing.synthetic writes this table" not in message
        assert "has no REAL producer yet" not in message

        message = self._refuse(config, cohort, "lineages")
        assert "writes no lineage_label" not in message
        assert "no REAL producer at all" not in message
        # The gap that IS live: the core alignment, panaroo and gubbins.
        assert "panaroo" in message and "gubbins" in message