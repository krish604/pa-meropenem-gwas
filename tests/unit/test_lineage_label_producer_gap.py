"""The `lineage_label` producer gap, and the two docstrings that hid it.

Three separate lies, all of which survive in the tree at `c730676` and are
fixed in `stages/convergence.py`, `stages/cooccurrence.py` and
`scripts/common/run_pipeline.py`:

**1. The refusal blamed the wrong thing.** `_CONVERGENCE_BLOCKERS["lineages"]`
and `_COOCCURRENCE_BLOCKERS["lineages"]` named *only* tool availability -
"blocked by panaroo ... and gubbins ...". That is true and it is not the whole
truth, in the worst possible direction: installing or fixing both tools does
**not** unblock these stages.

`stages.phylogeny.build_tree_outputs` (`phylogeny.py:383-393`) writes rows of
`sample_id` / `tree_tip_label` / `source` and nothing else, under
`REQUIRED_TREE_METADATA = ("sample_id", "tree_tip_label")`. So a REAL
`tree_metadata.tsv` has **no `lineage_label` column at all**, and
`load_tree_metadata` (`phylogeny.py:325`) then does
`row.get("lineage_label") or "unknown"` - every sample becomes the literal
string `"unknown"`. `is_effectively_empty(..., sentinel="unknown")` is `True`,
so the stage refuses with *"the sample -> lineage crosswalk carries no real
lineage label"* and no stated cause.

A reader who installed panaroo and fixed gubbins would arrive at that wall with
every explanation the message offered already exhausted.

**2. `cooccurrence`'s module docstring promised counts the module never
computed.** `cooccurrence.py:19-20` claimed "Lineage-stratified counts are
reported alongside the pooled counts". No such count exists:
`COOCCURRENCE_COLUMNS` (`cooccurrence.py:36-47`) has no lineage field, and
`test_pair`'s body never reads the `lineages` argument it accepts.

The fix is not to add stratification - that is a feature. It is to stop
promising it, and to make the crosswalk's real role (a completeness gate)
explicit.

**3. `TREE_METADATA_COLUMNS` is dead code that hides the gap.**
`phylogeny.py:38-44` declares the five-column schema the contract requires
(`docs/data_contract.md:138`), nothing imports it, and `build_tree_outputs`
passes the two-column `REQUIRED_TREE_METADATA` to `write_tsv` instead. So the
declared schema and the written schema disagree, and `write_tsv`'s
`extrasaction` behaviour means the extra keys are dropped silently - `source`
too, not only `lineage_label` and `st`. Reported, not fixed here: making the
producer write the declared columns is a producer change, and deciding *what a
lineage is* needs a rule `config/science.yaml` does not contain.

## Why the guards stay fail-closed

The point of these tests is not that the message is nicer. It is that a REAL run
whose crosswalk was produced by the pipeline's own REAL producer is **refused**,
and that the refusal now names the reason. Both properties are asserted here
against a `tree_metadata.tsv` written byte-for-byte the way
`build_tree_outputs` writes one - not a hand-made fixture that happens to carry
`lineage_label`, which is what made the gap invisible.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import StageError
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import cooccurrence as cc
from papipeline.stages import convergence as cv
from papipeline.stages import phylogeny as phy

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"

SAMPLES = ("S1", "S2", "S3", "S4", "S5", "S6")

#: Schemas copied from `test_convergence_real_inputs.py`, which is the
#: authoritative description of a real-shaped intermediate. Reusing the strings
#: rather than inventing them is deliberate: a hand-written header that is one
#: column short fails in `read_tsv` with a DataContractError about the *fixture*,
#: which tells you nothing about the producer gap being tested.
AMR = (
    "sample_id\tantibiotic\tdeterminant\tgene\tvariant\tdeterminant_type\t"
    "mechanism\tevidence_source\tdatabase\tdatabase_version\tconfidence\t"
    "claim_status\tidentity_pct\tcoverage_pct"
)
REGULATORS = (
    "sample_id\tgene\tvariant\treference\talternate\teffect\tmechanism\t"
    "confidence\tclaim_status"
)
MECHANISMS = (
    "sample_id\tantibiotic\tdeterminant\tmechanism\tevidence_level\tgene\t"
    "gene_type\tnotes"
)


def _amr_rows(antibiotic: str = "imipenem"):
    """blaOXA-1 across two lineages, fosA in three carriers of the second."""
    rows = [
        f"{s}\t{antibiotic}\tblaOXA-1\tblaOXA-1\t.\tAMR\t."
        f"\tamrfinderplus\tAMRFinderPlus\t2026-08-07.1\t.\tDETECTED\t99.0\t100.0"
        for s in ("S1", "S2", "S3", "S4")
    ]
    rows += [
        f"{s}\t{antibiotic}\tfosA\tfosA\t.\tAMR\t."
        f"\tamrfinderplus\tAMRFinderPlus\t2026-08-07.1\t.\tDETECTED\t99.5\t100.0"
        for s in ("S4", "S5", "S6")
    ]
    return rows


def _regulator_rows():
    return [
        f"{s}\toprD\tA{i}00T\tC\tT\ttruncating\tloss_of_function\t.\tDETECTED"
        for i, s in enumerate(SAMPLES, start=1)
    ]


def _mechanism_rows(antibiotic: str = "imipenem"):
    return [
        f"{s}\t{antibiotic}\tblaOXA-1\tenzymatic_inhibition\tDETECTED\t"
        f"blaOXA-1\tAMR\t."
        for s in SAMPLES
    ]


@pytest.fixture(scope="module")
def config():
    return load_config(SCIENCE, machine="laptop")


@pytest.fixture
def manifest():
    return SampleManifest(samples=[Sample(sample_id=s) for s in SAMPLES])


def _write(path: Path, header: str, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join([header.rstrip("\n")] + list(rows)) + "\n", encoding="utf-8"
    )
    return path


@pytest.fixture
def real_producer_inputs(tmp_path: Path):
    """An intermediate root that is complete except for the tree crosswalk.

    The crosswalk is written by `_write_as_build_tree_outputs_does`, not by a
    fixture that carries `lineage_label`.
    """
    root = tmp_path / "intermediate"
    _write(root / "stages" / "04_amr.tsv", AMR, _amr_rows())
    _write(root / "regulators" / "regulator_variants.tsv", REGULATORS,
           _regulator_rows())
    _write(root / "stages" / "05_mechanisms.tsv", MECHANISMS, _mechanism_rows())
    phylo = tmp_path / "phylogeny"
    phylo.mkdir(parents=True, exist_ok=True)
    (phylo / "tree_metadata.tsv").write_text(
        "sample_id\ttree_tip_label\n" + "".join(f"{s}\t{s}\n" for s in SAMPLES),
        encoding="utf-8",
    )
    return root, phylo


# --------------------------------------------------------------------------
# The producer really does omit lineage_label
# --------------------------------------------------------------------------


class TestTheProducerGapIsClosed:
    """Ruling R4 closed the gap this file documented, and these tests were the
    ones asserting it was open. They are inverted rather than deleted, because
    the shape of the mistake is the point: a five-column schema declared in the
    contract and in a module constant, a two-column schema actually written,
    and `write_tsv`'s `extrasaction="ignore"` closing the difference silently.
    """

    def test_build_tree_outputs_writes_lineage_label(self):
        source = (REPO / "papipeline" / "stages" / "phylogeny.py").read_text(
            encoding="utf-8"
        )
        body = source[source.index("def build_tree_outputs("):]
        body = body[: body.index("\ndef ")]
        assert "lineage_label" in body, (
            "build_tree_outputs stopped writing lineage_label, so the producer "
            "gap this file documented is open again - update the blockers in "
            "convergence.py and cooccurrence.py and this test"
        )

    def test_the_declared_schema_is_what_is_written(self):
        """`TREE_METADATA_COLUMNS` vs the columns passed to `write_tsv`.

        These two used to disagree - five declared, two written, extras dropped
        without a word - which is what let a producer with no `lineage_label`
        look like a producer that simply had nothing to say.
        """
        assert "lineage_label" in phy.TREE_METADATA_COLUMNS
        assert "lineage_label" not in phy.REQUIRED_TREE_METADATA, (
            "REQUIRED_TREE_METADATA is the subset read_tsv demands of an "
            "existing file; if it grows to require lineage_label, a "
            "two-column crosswalk from any older producer becomes unreadable"
        )
        source = (REPO / "papipeline" / "stages" / "phylogeny.py").read_text(
            encoding="utf-8"
        )
        body = source[source.index("def build_tree_outputs("):]
        body = body[: body.index("\ndef ")]
        assert "TREE_METADATA_COLUMNS,\n        header_comment=" in body, (
            "build_tree_outputs must write the contract's columns and record its "
            "provenance in the same call; writing REQUIRED_TREE_METADATA again "
            "drops lineage_label, st AND source without a word"
        )

    def test_load_tree_metadata_degrades_every_sample_to_unknown(self, tmp_path: Path):
        """The sentinel is still the sentinel - and that is now a DEFECT to be
        caught downstream, not the pipeline's own answer.

        R4 gives `lineage_label` a producer, so a crosswalk that arrives without
        the column means something upstream failed. `integration.build_master_table`
        refuses on this literal; `convergence`/`cooccurrence` refuse on it too.
        What must not happen is a result that reads as a finding while grouping
        isolates under the string "unknown".
        """
        path = tmp_path / "tree_metadata.tsv"
        path.write_text(
            "sample_id\ttree_tip_label\n" + "".join(f"{s}\t{s}\n" for s in SAMPLES),
            encoding="utf-8",
        )
        assert phy.load_tree_metadata(path) == {s: "unknown" for s in SAMPLES}
        assert phy.UNKNOWN_LINEAGE == "unknown", (
            "the sentinel the refusing stages match on is named in "
            "phylogeny.UNKNOWN_LINEAGE; if the two drift apart, the refusals "
            "stop firing on the value the reader actually sees"
        )


class TestTheProducerNowReadsStageThree:
    """R4: `lineage_label` IS the MLST sequence type, so the producer reads the
    table stage 3 wrote rather than deriving anything of its own."""

    def test_the_contracted_path_is_the_one_stage_three_writes(self):
        from papipeline.stages import mlst as mlst_stage

        source = (REPO / "papipeline" / "stages" / "mlst.py").read_text(
            encoding="utf-8"
        )
        assert '"mlst" / "mlst_results.tsv"' in source or (
            "mlst_results.tsv" in source
        )
        assert phy.MLST_RESULTS_RELATIVE == ("mlst", "mlst_results.tsv")
        assert mlst_stage.MLST_COLUMNS, "stage 3 declares the columns it writes"

    def test_a_missing_table_is_refused_naming_the_path(
        self, config, manifest, tmp_path, monkeypatch
    ):
        from papipeline.config.loader import RESULTS_ROOT_ENV
        from papipeline.errors import PipelineError

        monkeypatch.setenv(RESULTS_ROOT_ENV, str(tmp_path / "results"))
        with pytest.raises(PipelineError) as excinfo:
            phy.load_lineage_labels(config, RunMode.REAL, manifest)
        message = str(excinfo.value)
        # The path a reader can go and look at, not just "no MLST results".
        assert str(phy.mlst_results_path(config, RunMode.REAL)) in message
        assert "stage 3" in message
        # And the reason it is not a tooling problem, which is what the
        # convergence/cooccurrence messages used to imply by naming only
        # panaroo and gubbins.
        assert "panaroo" in message and "gubbins" in message


# --------------------------------------------------------------------------
# The refusals stay fail-closed, and name the real reason
# --------------------------------------------------------------------------


class TestRealIsStillRefusedWithTheProducersOwnCrosswalk:
    @pytest.mark.parametrize(
        "stage", [cv.run, cc.run], ids=["convergence", "cooccurrence"]
    )
    def test_it_refuses_rather_than_reporting_a_finding(
        self, stage, config, manifest, real_producer_inputs
    ):
        root, phylo = real_producer_inputs
        with pytest.raises(StageError) as excinfo:
            stage(
                config, manifest, RunMode.REAL,
                intermediate_root=root, antibiotic="imipenem", phylogeny_dir=phylo,
            )
        assert excinfo.value.context["missing"] == "lineages"

    @pytest.mark.parametrize(
        "stage", [cv.run], ids=["convergence"]
    )
    def test_the_refusal_names_both_gaps_not_just_the_tools(
        self, stage, config, manifest, real_producer_inputs
    ):
        """The actual fix. Naming only panaroo/gupbins is what was wrong.

        Asserting the *second* gap specifically: a message that mentioned the
        tools and stopped there is exactly the pre-fix behaviour, and it is the
        more dangerous half because it looks actionable.

**cooccurrence was removed from this parametrisation in round 12**, and
        the removal is the point rather than a narrowing. This test pins
        "build_tree_outputs writes no lineage_label", which was true when
        written and is now false: `phylogeny.build_tree_outputs` emits
        `lineage_label` from `load_lineage_labels` (`phylogeny.py:533`), and
        ruling R4 makes that label the stage-3 MLST sequence type. A test that
        asserts a false statement about the producer is worse than no test,
        because it is a tripwire against the fix. cooccurrence's corrected
        refusal is asserted in
        `TestCooccurrenceRefusalNamesTheProducerThatExists` below.

        **R10 assertion change, convergence only.** The ``"clade" in message``
        assertion was dropped for `cv.run` in round 12, because the claim it
        guarded became false rather than because the refusal weakened. The text
        it policed said a lineage needed "a clade-assignment rule (what counts
        as a lineage, and at what support)" and that `science.yaml` carried "no
        lineage rule at all". Ruling R4 chose the rule - a lineage IS the MLST
        sequence type - and `science.yaml` carries `lineage.method: st`. A
        message still arguing for a clade rule would send a reader to author a
        scientific decision this project has already made.

        What is asserted instead is the *corrected* text, which is the stronger
        property: it must still name the surviving tool gap, must still name
        the producer of `lineage_label` so a reader knows which stage to run, and
        must say the second gap is closed rather than silently dropping it. A
        message that merely stopped mentioning `clade` would pass a weaker test
        than the one it replaces.
        """
        root, phylo = real_producer_inputs
        with pytest.raises(StageError) as excinfo:
            stage(
                config, manifest, RunMode.REAL,
                intermediate_root=root, antibiotic="imipenem", phylogeny_dir=phylo,
            )
        message = str(excinfo.value)
        # Gap 1, the tools. Still the live blocker.
        assert "panaroo" in message and "gubbins" in message
        # The producer of lineage_label, by name.
        assert "build_tree_outputs" in message
        assert "lineage_label" in message
        # The corrected half: the closed gap is stated as closed, and the
        # obsolete argument is gone rather than merely softened.
        assert "has been removed" in message
        assert "upstream stage failed" in message

    def test_stub_still_raises_on_both(self, config, manifest):
        """The fabrication guard is untouched by this change."""
        for stage in (cv.run, cc.run):
            with pytest.raises(NotImplementedError):
                stage(config, manifest, RunMode.STUB)


class TestCooccurrenceRefusalNamesTheProducerThatExists:
    """The round-12 correction to cooccurrence's own blocker text.

    convergence's message is still right for convergence - see
    `stages/convergence.py` - but cooccurrence's asserted that
    `build_tree_outputs` "writes no lineage_label" and that "lineage_label has
    no REAL producer at all; only papipeline.testing.synthetic writes it".
    Both were true when written and both are now false:
    `phylogeny.build_tree_outputs` writes `lineage_label` (`phylogeny.py:533`)
    and `phylogeny.load_lineage_labels` builds it from stage 3's MLST sequence
    type under ruling R4. TRIAGE ranked this stale text as a class-(c) defect.

    A reader who followed the old message would have gone looking for a
    producer that does not need to be built, and would not have learned that
    the thing they need is stage 3's `mlst_results.tsv`.
    """

    def test_it_refuses_and_names_stage_three_not_a_missing_producer(
        self, config, manifest, real_producer_inputs
    ):
        root, phylo = real_producer_inputs
        with pytest.raises(StageError) as excinfo:
            cc.run(
                config, manifest, RunMode.REAL,
                intermediate_root=root, antibiotic="imipenem", phylogeny_dir=phylo,
            )
        message = str(excinfo.value)
        assert excinfo.value.context["missing"] == "lineages"
        # The producer that exists, named by the stage that writes it.
        assert "stage 3" in message and "mlst_results.tsv" in message
        assert "stage 9" in message and "tree_metadata.tsv" in message

    def test_it_no_longer_claims_the_label_has_no_real_producer(self):
        """The false sentence, asserted absent so it cannot come back."""
        blockers = cc._COOCCURRENCE_BLOCKERS["lineages"]
        assert "no REAL producer at all" not in blockers
        assert "writes no lineage_label" not in blockers
        # And the trade this round declined to make is stated where a reader
        # deciding whether to change it will find it.
        assert "backlog" in blockers




# --------------------------------------------------------------------------
# cooccurrence: no stratified counts are promised
# --------------------------------------------------------------------------


class TestCooccurrencePromisesNothingItDoesNotCompute:
    def test_the_output_schema_has_no_lineage_field(self):
        assert not any("lineage" in c for c in cc.COOCCURRENCE_COLUMNS)

    def test_the_module_docstring_makes_no_stratification_claim(self):
        source = (REPO / "papipeline" / "stages" / "cooccurrence.py").read_text(
            encoding="utf-8"
        )
        docstring = source[: source.index('"""', source.index('"""') + 3) + 3]
        for phrase in (
            "Lineage-stratified counts are reported",
            "stratified counts are reported alongside",
        ):
            assert phrase not in docstring, (
                f"the module docstring claims {phrase!r}, which no code in this "
                f"module computes"
            )
        # And it must say the opposite, so the fix cannot be a silent deletion.
        assert "pooled" in docstring.lower()

    def test_test_pair_ignores_the_lineages_it_is_given(self):
        """The dead parameter, pinned so it cannot drift into a half-feature.

        Documented as deliberately unread rather than removed (removal is
        proposed, not done). This asserts the *behaviour* - identical results
        with and without a crosswalk - so if someone later wires stratification
        in, this fails and the docstring has to be rewritten with it.
        """
        set_a = {"s1": {"gA"}, "s2": {"gA"}, "s3": set(), "s4": set()}
        set_b = {"s1": {"gB"}, "s2": set(), "s3": {"gB"}, "s4": set()}
        samples = ["s1", "s2", "s3", "s4"]
        without = cc.test_pair("gA", "gB", set_a, set_b, "gene_gene", samples)
        with_crosswalk = cc.test_pair(
            "gA", "gB", set_a, set_b, "gene_gene", samples,
            lineages={"s1": "L1", "s2": "L2", "s3": "L1", "s4": "L2"},
        )
        assert without == with_crosswalk

    def test_a_real_run_with_a_degenerate_crosswalk_is_refused_not_pooled(self):
        """Why an unread parameter is still a declared input.

        Without the gate, an all-`unknown` crosswalk would publish pooled counts
        carrying no stratification at all, with nothing in the schema to say so.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            root = tmp / "intermediate"
            phylo = tmp / "phylogeny"
            _write(root / "stages" / "04_amr.tsv", AMR, _amr_rows())
            _write(root / "regulators" / "regulator_variants.tsv", REGULATORS,
                   _regulator_rows())
            _write(root / "stages" / "05_mechanisms.tsv", MECHANISMS,
                   _mechanism_rows())
            _write(phylo / "tree_metadata.tsv",
                   "sample_id\ttree_tip_label\tlineage_label",
                   [f"{s}\t{s}\tunknown" for s in SAMPLES])

            config = load_config(SCIENCE, machine="laptop")
            manifest = SampleManifest(samples=[Sample(sample_id=s) for s in SAMPLES])
            with pytest.raises(StageError) as excinfo:
                cc.run(
                    config, manifest, RunMode.REAL,
                    intermediate_root=root, antibiotic="imipenem", phylogeny_dir=phylo,
                )
            assert "lineages" in str(excinfo.value)
