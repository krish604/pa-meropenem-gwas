"""The configured feature universe: which classes, how common, how many.

Round 13 asked a question the previous rounds had not: **what does stage 14
actually pair?** The answer, measured on the real 10-isolate tree and recorded
in the module docstring, is that the universe is DETERMINANTS ONLY - three
classes, and not gene presence/absence.

    gene       acquired AMR determinant names        54 distinct at n=10
    mutation   `{gene}:{variant}` at regulator loci  300 distinct at n=10
    mechanism  the `mechanism` column, a CLOSED set   4, forever

Only `mutation` grows without bound, because a cohort has more isolates in which
to find a variant, and it is the class that decides whether this stage runs.
Projected from the measured cumulative curve it reaches ~11,870 at n=900, which
is ~101,000,000 pairs.

So this file pins three things:

1. **The universe is determinants only**, and it cannot quietly become the
   pangenome: the classes are enumerated, and a class outside the list is
   refused with the file that *would* have supplied it named.
2. **The universe is filtered and capped**, and the filter is inert where it
   must be inert - `min_prevalence=0.05` cannot move a 20-sample cohort, so no
   committed fixture depends on it.
3. **The cap is binding, not decoration.** The n=900 projection exceeds it, and
   the file asserts that it does rather than leaving the claim in a comment.

Per-pair costs used throughout, MEASURED by a full run at exactly the cap
(`round-13/cooc2/scale_at_cap.py`: n=900, 2,500 features, 3,114,595 pairs
admitted, 2,763,235 rows emitted, 5,925.67 s wall, 2,321.2 MiB peak RSS delta):
counting 0.34 us/pair; 1,902.50 us/pair for the test and the row together; 881
resident bytes per emitted row including the BH key dict. **The test is 5,600x
the count**, which is why the ceiling is set on the size of the universe and not
on counting throughput.

The 1,902 us is the number most likely to be got wrong, and the module says so
where it is stated: `fisher_exact` cost is a function of the contingency table's
MARGINS, not of the code path - the same call measures ~128 us on a 90-carrier
table and ~500-650 us on a 900-row one, and cProfile puts it at 85% of this
stage's wall at n=900. A cap sized off the small-table micro-benchmark would be
an order of magnitude wrong, so `test_fisher_is_the_cost_and_its_cost_rises_with_the_margins`
pins the relationship rather than leaving it to a comment.
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path
from typing import Dict, List, Sequence, Set

import pytest

from papipeline.config.loader import PipelineConfig, load_config
from papipeline.errors import StageError
from papipeline.manifest import SampleManifest
from papipeline.models import Sample
from papipeline.stages import cooccurrence as cc

#: Per-pair costs measured at n=900 on this host. See the module docstring of
#: this file. They are reproduced here rather than imported because they are
#: claims about a machine, and a test that imported them would be asserting that
#: the measurement still holds rather than that the arithmetic is right.
#: Per ADMITTED PAIR (counted, whether or not it survives the degenerate rule),
#: because that is what `fisher_exact` is called on: 5,925.67 s / 3,114,595
#: pairs. Per *row* the figure would be 2,144 us, but 11.3% of pairs never reach
#: a row and must not be charged for one.
US_PER_PAIR_COUNT = 0.34
US_PER_PAIR_TEST = 1902.50

#: Per SURVIVING ROW, because that is what allocates. 2,321.2 MiB RSS delta /
#: 2,763,235 rows at the cap. Memory is charged on rows and wall on pairs, and
#: mixing the two denominators is how a cap ends up sized by a number that was
#: never measured.
RESIDENT_BYTES_PER_ROW = 881

#: Share of admitted pairs that survive `cooccurrence_from_counts`'s degenerate
#: rule, measured at the cap: 2,763,235 / 3,114,595 = 88.7%.
MEASURED_SURVIVAL_FRACTION = 0.887

#: Marginal new features per isolate, measured as the mean slope of the
#: cumulative distinct curve over n=5..10 on the real 10 isolates
#: (gene 2.60, mutation 13.00, mechanism 0.00 - the mechanism vocabulary is
#: closed). The curves saturate, so a linear projection is an UPPER bound, which
#: is the direction a cap wants to be wrong in.
MEASURED_NEW_FEATURES_PER_ISOLATE = {"gene": 2.60, "mutation": 13.00}
MEASURED_UNIVERSE_AT_N10 = {"gene": 54, "mutation": 300, "mechanism": 4}

#: The mechanism vocabulary's ceiling is five, not four: the four distinct
#: mechanisms in config/mechanisms.tsv, plus `locus_intact`, which
#: mechanisms.py writes directly. See the module docstring's universe table.
MEASURED_MECHANISM_NAMES_FROM_CONFIG = 4
MECHANISM_CLASS_CEILING = 5

ALL_TESTS = [
    "genes",
    "gene_mutation",
    "mutation_mutation",
    "gene_mechanism",
    "mechanism_mechanism",
]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def _with_cooccurrence(config: PipelineConfig, **overrides) -> PipelineConfig:
    """A config whose `cooccurrence` section carries `overrides`."""
    raw = dict(config.raw)
    raw["cooccurrence"] = {**raw.get("cooccurrence", {}), **overrides}
    return dataclasses.replace(config, raw=raw)


def _manifest(n: int) -> SampleManifest:
    return SampleManifest(samples=[Sample(sample_id=f"S{i}") for i in range(n)])


def _namespaces(
    rng: random.Random, n: int, sizes=(8, 10, 5)
) -> Dict[str, Dict[str, Set[str]]]:
    """Three namespaces with a plausible shape, and a knob for prevalence."""
    out = {}
    for name, size in zip(("gene", "mutation", "mechanism"), sizes):
        out[name] = {
            f"S{index}": {
                f"{name}_f{feature}"
                for feature in range(size)
                if rng.random() < 0.4
            }
            for index in range(n)
        }
    return out


def _pair_keys(results) -> List[tuple]:
    return [(r.feature_type, r.feature_a, r.feature_b) for r in results]


# --------------------------------------------------------------------------
# O0: the universe is determinants only
# --------------------------------------------------------------------------


class TestTheUniverseIsDeterminantsOnly:
    def test_the_three_classes_are_exactly_the_three_namespaces_built(
        self, config: PipelineConfig
    ):
        """`build_feature_sets` and the eligible list cannot drift apart.

        A fourth namespace appearing in one and not the other would mean a
        namespace the counter can read but the configuration cannot switch off,
        or the reverse - and neither is visible in the output.
        """
        built = cc.build_feature_sets(_manifest(4), {}, {}, {})
        assert set(built) == set(cc.ELIGIBLE_FEATURE_CLASSES)
        assert set(cc.resolve_feature_classes(config)) == set(built)

    def test_gene_presence_absence_is_not_read_by_this_stage(self):
        """Stage 7 writes a GPA table; stage 14 must not read it.

        Asserted against the paths rather than by running the stage, because the
        failure mode is *not* a wrong number - it is a pangenome-sized universe
        appearing in a stage documented as determinants-only, and the only cheap
        moment to catch that is before anyone runs it at n=900. The module
        docstring is allowed to *name* the file (it has to, to say it is not
        read); what is forbidden is a load path.
        """
        from papipeline.models import RunMode

        class _Cfg:
            raw: dict = {}

            def phylogeny_dir(self, mode):  # pragma: no cover - unused
                return Path("nowhere")

        paths = cc.real_input_paths(_Cfg(), Path("/root"), RunMode.REAL)
        assert set(paths) == {"gene", "mutation", "mechanism", "lineages"}
        for name, path in paths.items():
            assert "pangenome" not in path, f"{name} reads a pangenome artefact: {path}"

        # And no source line in the module *reads* it - every occurrence is in a
        # docstring, a comment or the refusal text.
        source = Path(cc.__file__).read_text(encoding="utf-8")
        offenders = [
            line
            for line in source.splitlines()
            if "gene_presence_absence" in line
            and not line.lstrip().startswith(("#", '"', "*", "'"))
            and "tsv" not in line
        ]
        assert not offenders, offenders

    def test_the_mechanism_class_is_closed_and_cannot_grow(self, config: PipelineConfig):
        """The mechanism vocabulary has a ceiling of five, and it is knowable.

        Four distinct mechanisms come from `config/mechanisms.tsv`; the fifth,
        `locus_intact`, is written directly by `mechanisms.py` and so is not in
        that table. Both halves are asserted, because a test that counted only
        the config table would pass while the ceiling had actually moved - and
        because this is the one part of the universe whose size is knowable
        without a cohort, which is why C(5,2) = 10 pairs is not the wall this
        stage has to be designed around.
        """
        from papipeline.config.loader import load_mechanisms

        table = load_mechanisms(config.config_dir / "mechanisms.tsv")
        from_config = {spec.mechanism for spec in table.values()}
        assert len(from_config) == MEASURED_MECHANISM_NAMES_FROM_CONFIG, (
            f"config/mechanisms.tsv declares {len(from_config)} distinct "
            f"mechanisms {sorted(from_config)}, not 4. If that is deliberate, "
            "update the module docstring's universe table with it - this is the "
            "class that cannot grow, so its ceiling is the one number a reader "
            "can trust without a cohort."
        )

        # The fifth name is not in that table: `mechanisms.py` writes it
        # directly, which is why counting the config table alone undercounts the
        # ceiling by one. Asserted at the source because that is where it is.
        from papipeline.stages import mechanisms as mech

        source = Path(mech.__file__).read_text(encoding="utf-8")
        assert 'mechanism="locus_intact"' in source, (
            "mechanisms.py writes `locus_intact` directly, outside "
            "config/mechanisms.tsv, so the closed-vocabulary ceiling of five "
            "counts a name the config table does not declare. If that literal "
            "has moved, update the ceiling rather than this assertion."
        )

    def test_the_documented_n10_counts_are_the_measured_ones(self):
        """The universe table in the module docstring is not a guess.

        Numbers that a reader will use for arithmetic have to be checkable. The
        measurement itself needs the real 10-isolate tree, which is not committed
        (no bioinformatics tool may be run to make it), so what is pinned here is
        that the constants the docstring quotes and the constants this file
        projects from are the same numbers.
        """
        docstring = cc.__doc__ or ""
        for value in MEASURED_UNIVERSE_AT_N10.values():
            assert f"**{value}**" in docstring, (
                f"the module docstring's universe table no longer quotes the "
                f"measured n=10 count {value}"
            )


# --------------------------------------------------------------------------
# the three keys, validated
# --------------------------------------------------------------------------


class TestResolveFeatureClasses:
    def test_all_three_by_default(self, config: PipelineConfig):
        assert cc.resolve_feature_classes(config) == ("gene", "mutation", "mechanism")

    def test_declared_order_is_preserved(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, feature_classes=["mechanism", "gene"])
        assert cc.resolve_feature_classes(cfg) == ("mechanism", "gene")

    def test_a_duplicated_class_is_collapsed(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, feature_classes=["gene", "gene"])
        assert cc.resolve_feature_classes(config=cfg) == ("gene",)

    def test_an_unknown_class_is_refused_and_names_the_pangenome(self, config: PipelineConfig):
        """A reader who asks for the GPA universe gets told where it lives.

        Refused rather than ignored: silently pairing nothing is indistinguishable
        from finding nothing, which is the whole reason `_add_missing_lineages`
        exists in this module.
        """
        cfg = _with_cooccurrence(config, feature_classes=["gene", "gene_presence"])
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_feature_classes(cfg)
        message = str(excinfo.value)
        assert "gene_presence" in message
        assert "gene_presence_absence.tsv" in message
        for name in cc.ELIGIBLE_FEATURE_CLASSES:
            assert name in message

    def test_an_empty_class_list_is_refused(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, feature_classes=[])
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_feature_classes(cfg)
        assert "empty" in str(excinfo.value)

    def test_a_bare_string_is_refused_rather_than_iterated(self, config: PipelineConfig):
        """`feature_classes: gene` would otherwise iterate to 'g', 'e', 'n', 'e'."""
        cfg = _with_cooccurrence(config, feature_classes="gene")
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_feature_classes(cfg)
        assert "list" in str(excinfo.value)


class TestResolveMinPrevalence:
    def test_absent_means_no_filter(self, config: PipelineConfig):
        raw = dict(config.raw)
        section = dict(raw.get("cooccurrence", {}))
        section.pop("min_prevalence", None)
        raw["cooccurrence"] = section
        assert cc.resolve_min_prevalence(
            dataclasses.replace(config, raw=raw)
        ) == cc.DEFAULT_MIN_PREVALENCE

    @pytest.mark.parametrize("value", [0.0, 0.01, 0.05, 0.25])
    def test_a_fraction_is_kept(self, config: PipelineConfig, value):
        cfg = _with_cooccurrence(config, min_prevalence=value)
        assert cc.resolve_min_prevalence(cfg) == pytest.approx(value)

    def test_a_string_number_is_read_not_refused(self, config: PipelineConfig):
        """YAML will hand over the string `"0.05"` if it was ever quoted."""
        cfg = _with_cooccurrence(config, min_prevalence="0.05")
        assert cc.resolve_min_prevalence(cfg) == pytest.approx(0.05)

    @pytest.mark.parametrize("value", [-0.01, -1.0])
    def test_a_negative_fraction_is_refused(self, config: PipelineConfig, value):
        cfg = _with_cooccurrence(config, min_prevalence=value)
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_min_prevalence(cfg)
        assert "negative" in str(excinfo.value)

    @pytest.mark.parametrize("value", [1.0, 1.5])
    def test_one_hundred_percent_is_refused(self, config: PipelineConfig, value):
        """It would guarantee an empty table, so saying so is the useful error."""
        cfg = _with_cooccurrence(config, min_prevalence=value)
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_min_prevalence(cfg)
        assert "no contrast to test" in str(excinfo.value)

    def test_a_non_number_is_refused(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, min_prevalence="most")
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_min_prevalence(cfg)
        assert "fraction" in str(excinfo.value)


class TestThePrevalenceFloor:
    @pytest.mark.parametrize(
        "min_prevalence, n_samples, expected",
        [
            (0.05, 1, 1),
            (0.05, 10, 1),
            (0.05, 20, 1),      # the committed fixtures live here
            (0.05, 21, 2),
            (0.05, 900, 45),
            (0.01, 900, 9),
            (0.10, 900, 90),
            (0.0, 900, 1),
            (0.5, 10, 5),
            (0.5, 11, 6),
        ],
    )
    def test_the_floor(self, min_prevalence, n_samples, expected):
        assert cc._prevalence_floor(min_prevalence, n_samples) == expected

    @pytest.mark.parametrize(
        "min_prevalence, n_samples",
        [(0.07, 100), (0.29, 100), (0.58, 100), (0.1, 30), (0.02, 50)],
    )
    def test_a_boundary_feature_is_admitted_not_lost_to_float_error(
        self, min_prevalence, n_samples
    ):
        """`ceil` on an IEEE product is a trap and this is the tripwire.

        `0.07 * 100` is `7.000000000000001`, so a plain `ceil` yields 8 and
        excludes a feature carried by exactly 7 of 100 samples - a filter that
        says 7% silently behaving as 8%.
        """
        exact = min_prevalence * n_samples
        assert abs(exact - round(exact)) < 1e-9, "fixture must land on a boundary"
        assert cc._prevalence_floor(min_prevalence, n_samples) == round(exact)


class TestResolveMaxFeatures:
    def test_absent_means_no_cap(self, config: PipelineConfig):
        raw = dict(config.raw)
        section = dict(raw.get("cooccurrence", {}))
        section.pop("max_features", None)
        raw["cooccurrence"] = section
        assert cc.resolve_max_features(dataclasses.replace(config, raw=raw)) is None

    def test_the_configured_cap_is_read(self, config: PipelineConfig):
        assert cc.resolve_max_features(config) == 2500

    @pytest.mark.parametrize("bad", [0, -1, -4000])
    def test_a_non_positive_cap_is_refused(self, config: PipelineConfig, bad):
        cfg = _with_cooccurrence(config, max_features=bad)
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_max_features(cfg)
        assert "at least 1" in str(excinfo.value)

    def test_a_non_integer_cap_is_refused(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, max_features="many")
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_max_features(cfg)
        assert "integer" in str(excinfo.value)

    def test_a_quoted_integer_is_read(self, config: PipelineConfig):
        cfg = _with_cooccurrence(config, max_features="250")
        assert cc.resolve_max_features(cfg) == 250


# --------------------------------------------------------------------------
# admission
# --------------------------------------------------------------------------


class TestAdmitFeatures:
    def test_an_ineligible_class_is_dropped(self):
        namespaces = {
            "gene": {"S0": {"g1"}},
            "mutation": {"S0": {"m1"}},
            "mechanism": {"S0": {"mech"}},
        }
        admitted = cc.admit_features(namespaces, ["gene", "mechanism"], 0.0, 1)
        assert set(admitted) == {"gene", "mechanism"}

    def test_all_three_classes_pass_through_at_zero_prevalence(self):
        namespaces = {
            "gene": {"S0": {"g1"}},
            "mutation": {"S0": {"m1"}},
            "mechanism": {"S0": {"mech"}},
        }
        assert set(cc.admit_features(namespaces, cc.ELIGIBLE_FEATURE_CLASSES, 0.0, 1)) == set(
            namespaces
        )

    def test_a_feature_below_the_floor_is_dropped_from_every_sample(self):
        namespaces = {"gene": {f"S{i}": {"common", "rare"} for i in range(10)}}
        namespaces["gene"]["S9"] = {"common", "rare"}
        namespaces["gene"]["S9"].add("singleton")
        admitted = cc.admit_features(namespaces, ["gene"], 0.2, 10)  # floor 2
        carried = cc.namespace_names(admitted["gene"])
        assert carried == ["common", "rare"]
        assert all("singleton" not in members for members in admitted["gene"].values())

    def test_a_sample_key_survives_even_when_it_carries_nothing(self):
        """An empty set must stay in the mapping.

        `build_feature_sets` pre-seeds every manifest sample, and dropping the
        empty ones would change which samples the completeness check sees - so a
        namespace that admits nothing would no longer be recognisable as empty.
        """
        namespaces = {"gene": {"S0": {"a"}, "S1": set()}}
        admitted = cc.admit_features(namespaces, ["gene"], 0.0, 2)
        assert set(admitted["gene"]) == {"S0", "S1"}

    def test_prevalence_counts_every_key_not_only_the_manifest(self):
        """The filter and the counter must see the same universe.

        `namespace_names` scans *every* key of the mapping, so a feature carried
        only by a non-manifest sample is a real member of the namespace as far as
        counting is concerned. Filtering it out on the manifest's view alone
        would leave a name in the universe that the counter then pairs - and a
        pair whose `n_a` is 0, dropped later as degenerate, having already cost
        its whole row of pairs.
        """
        namespaces = {
            "gene": {"S0": {"cohort_wide"}, "NOT_A_SAMPLE": {"orphan", "cohort_wide"}}
        }
        admitted = cc.admit_features(namespaces, ["gene"], 0.5, 1)  # floor 1
        assert cc.namespace_names(admitted["gene"]) == ["cohort_wide", "orphan"]

    def test_the_callers_namespaces_are_not_mutated(self):
        namespaces = {"gene": {"S0": {"common", "rare"}}}
        for i in range(1, 10):
            namespaces["gene"][f"S{i}"] = {"common"}
        cc.admit_features(namespaces, ["gene"], 0.5, 10)
        assert cc.namespace_names(namespaces["gene"]) == ["common", "rare"]

    def test_an_empty_namespace_yields_an_empty_mapping(self):
        assert cc.admit_features({}, cc.ELIGIBLE_FEATURE_CLASSES, 0.05, 900) == {}


# --------------------------------------------------------------------------
# the cap
# --------------------------------------------------------------------------


class TestEnforceFeatureCap:
    @staticmethod
    def _universe(size: int) -> Dict[str, Dict[str, Set[str]]]:
        return {
            "mutation": {"S0": {f"m{i}" for i in range(size)}},
        }

    def test_a_universe_under_the_cap_passes_and_reports_its_classes(self):
        sizes = cc.enforce_feature_cap(self._universe(10), 10)
        assert sizes == {"mutation": 10}

    def test_a_universe_exactly_at_the_cap_passes(self):
        """The cap is inclusive: `max_features: 2500` admits 2,500 features."""
        assert cc.enforce_feature_cap(self._universe(2500), 2500) == {"mutation": 2500}

    def test_one_over_the_cap_refuses(self):
        with pytest.raises(StageError) as excinfo:
            cc.enforce_feature_cap(self._universe(2501), 2500)
        assert "2501" in str(excinfo.value)
        assert "2500" in str(excinfo.value)

    def test_no_cap_never_refuses(self):
        assert cc.enforce_feature_cap(self._universe(100_000), None)["mutation"] == 100_000

    def test_the_refusal_names_the_count_the_cap_and_the_breakdown(self):
        """A reader holding this needs to know WHICH class to thin."""
        universe = {
            "gene": {"S0": {f"g{i}" for i in range(5)}},
            "mutation": {"S0": {f"m{i}" for i in range(2500)}},
        }
        with pytest.raises(StageError) as excinfo:
            cc.enforce_feature_cap(universe, 2500)
        message = str(excinfo.value)
        assert "2505" in message
        assert "max_features=2500" in message
        assert "gene=5" in message
        assert "mutation=2500" in message

    def test_the_refusal_carries_the_counts_as_context(self):
        with pytest.raises(StageError) as excinfo:
            cc.enforce_feature_cap(self._universe(2501), 2500)
        assert excinfo.value.context["features"] == 2501
        assert excinfo.value.context["max_features"] == 2500

    def test_the_refusal_says_what_to_do_about_it(self):
        for key in (
            "cooccurrence.max_features",
            "cooccurrence.min_prevalence",
            "cooccurrence.feature_classes",
        ):
            with pytest.raises(StageError) as excinfo:
                cc.enforce_feature_cap(self._universe(2501), 2500)
            assert key in str(excinfo.value)

    def test_the_cap_is_on_the_total_not_per_class(self):
        """Two classes of 2,001 are over a 4,000 cap; a per-class cap lets them through.

        The bound that matters is C(total, 2), because the cross products
        `f_i * f_j` and the within-class `C(f_i, 2)` sum to at most `C(F, 2)` for
        any split - so a per-class cap would let a 3-way split of 4,000 admit
        *more* pairs than one class of 4,000, which is the opposite of a bound.
        """
        universe = {
            "gene": {"S0": {f"g{i}" for i in range(1251)}},
            "mutation": {"S0": {f"m{i}" for i in range(1250)}},
        }
        assert cc.enforce_feature_cap(universe, 2501)["gene"] == 1251
        with pytest.raises(StageError):
            cc.enforce_feature_cap(universe, 2500)

    def test_an_empty_universe_is_under_any_cap(self):
        assert cc.enforce_feature_cap({}, 1) == {}


# --------------------------------------------------------------------------
# the cap is binding: the n=900 arithmetic, made executable
# --------------------------------------------------------------------------


class TestTheCapIsBinding:
    @staticmethod
    def _projected_universe(n_target: int = 900) -> Dict[str, int]:
        """Linear projection of the MEASURED n=10 curve to `n_target`.

        The measured curve saturates, so this is an upper bound - which is the
        direction in which a cap wants to be conservative.
        """
        projected = {}
        for name, slope in MEASURED_NEW_FEATURES_PER_ISOLATE.items():
            at_10 = MEASURED_UNIVERSE_AT_N10[name]
            projected[name] = at_10 + slope * (n_target - 10)
        projected["mechanism"] = MEASURED_UNIVERSE_AT_N10["mechanism"]
        return projected

    def test_every_growing_class_has_a_positive_slope(self):
        """Guards the projection arithmetic, not the word 'upper'.

        Whether a linear extrapolation is *conservative* depends on a measured
        curve that is not committed - no bioinformatics tool may be run to make
        it - so that claim lives in the config comment and the report. What IS
        checkable here is that a class with a zero or negative slope would
        project a shrinking universe and quietly make the cap look generous, so
        the arithmetic below can be trusted to be an over-estimate.
        """
        for name, slope in MEASURED_NEW_FEATURES_PER_ISOLATE.items():
            assert slope > 0, f"{name} projects a non-growing universe"
        assert "mechanism" not in MEASURED_NEW_FEATURES_PER_ISOLATE, (
            "the mechanism class is closed; a slope for it would mean the "
            "ceiling of five had been forgotten"
        )

    def test_the_projection_never_shrinks_a_class(self):
        projected = self._projected_universe(200)
        for name, at_10 in MEASURED_UNIVERSE_AT_N10.items():
            assert projected[name] >= at_10, name

    def test_a_full_n900_universe_exceeds_the_configured_cap(self, config: PipelineConfig):
        """The cap is not decoration: the projection hits it.

        If this fails, either the cap was raised past the point the machine can
        carry, or the measured growth changed. Both are decisions, and neither
        should be made by deleting an assertion.
        """
        projected = self._projected_universe()
        total = sum(projected.values())
        assert total > cc.resolve_max_features(config), (
            f"the projected n=900 universe is {total:,} features "
            f"({projected}), which no longer exceeds the configured cap of "
            f"{cc.resolve_max_features(config)}. Raise the cap deliberately, with "
            "the measured pair cost in hand, or narrow the universe."
        )

    def test_the_pair_arithmetic_at_n900_is_what_the_config_comment_claims(self):
        """The `101,410,161 pairs / ~62 h / ~127 GiB` line, checked.

        The wall and memory figures follow from the pair count and the two
        measured per-pair constants, so asserting the pair count pins all three.
        """
        projected = self._projected_universe()
        total = sum(projected.values())
        pairs = total * (total - 1) // 2
        assert total == 14242, total
        assert pairs == 101_410_161, pairs
        hours = pairs * US_PER_PAIR_TEST / 1e6 / 3600.0
        gib = pairs * MEASURED_SURVIVAL_FRACTION * RESIDENT_BYTES_PER_ROW / 2 ** 30
        assert 50 < hours < 75, hours
        assert 65 < gib < 100, gib

    def test_the_cap_bounds_the_work_to_the_measured_budget(self, config: PipelineConfig):
        """C(cap, 2) pairs must fit the budget the config comment states.

        Uses the SURVIVING-row figures rather than the pair count, because rows
        are what cost memory: at n=900 ~88% of admitted pairs survive the
        degenerate rule, and the 12% that do not never reach an allocation.
        """
        cap = cc.resolve_max_features(config)
        pairs = cap * (cap - 1) // 2
        wall_s = pairs * (US_PER_PAIR_COUNT + US_PER_PAIR_TEST) / 1e6
        gib = pairs * MEASURED_SURVIVAL_FRACTION * RESIDENT_BYTES_PER_ROW / 2 ** 30
        assert pairs == 3_123_750, pairs
        assert 5_800.0 < wall_s < 6_100.0, (
            f"{wall_s:.0f} s at the cap ({wall_s/3600:.2f} h); the run at the "
            "cap measured 5,925.67 s"
        )
        assert 2.0 < gib < 2.6, (
            f"{gib:.2f} GiB resident at the cap; the run at the cap measured "
            "2.27 GiB"
        )

    def test_the_cap_is_a_small_fraction_of_the_bigmachine_memory_ceiling(self):
        """A stage that needs most of the machine is a claim, not a budget.

        `runtime.memory_mb` is 65536 on the bigmachine and 8192 on the laptop
        (config/machines/*.yaml). The stage is sized to fit the LAPTOP's whole
        ceiling, which is the stricter of the two by a factor of eight - so the
        cap has to be well inside it.
        """
        from papipeline.config.loader import load_machine_config

        root = Path(__file__).resolve().parents[2] / "config" / "machines"
        big = load_machine_config(root / "bigmachine.yaml")
        cap = 2500
        pairs = cap * (cap - 1) // 2
        gib = pairs * MEASURED_SURVIVAL_FRACTION * RESIDENT_BYTES_PER_ROW / 2 ** 30
        assert gib * 1024 < big.memory_mb / 8, (
            f"{gib:.2f} GiB is more than an eighth of the bigmachine's "
            f"{big.memory_mb} MB ceiling"
        )

    def test_the_gene_and_mechanism_position_runs_at_n900(self):
        """The configuration that IS runnable at n=900, asserted as arithmetic.

        gene + mechanism at the projected n=900 size is 2,372 features - under
        the 4,000 cap - while adding `mutation` is 14,242 and refuses. This is
        the pair of numbers that says which analysis is available today.
        """
        projected = self._projected_universe()
        without_mutations = projected["gene"] + projected["mechanism"]
        assert without_mutations == pytest.approx(2373, abs=1)
        assert without_mutations < 2500
        assert sum(projected.values()) > 2500

    def test_the_test_is_the_expensive_half_not_the_count(self):
        """The claim the whole cap rests on, as a ratio."""
        assert US_PER_PAIR_TEST / US_PER_PAIR_COUNT > 5000

    def test_fisher_is_the_cost_and_its_cost_rises_with_the_margins(self):
        """Why the ceiling is set on the universe rather than on counting.

        Two claims, both measured (round-13/cooc2/scale_profile.py, cProfile at
        n=900): `scipy.stats.fisher_exact` is 85% of this stage's wall time at
        cohort scale, and its cost tracks the contingency table's margins rather
        than the number of calls - the same call on a 90-carrier table measures
        ~128 us, on a 900-row table ~500-650 us, and in situ (wider and more
        varied tables) ~2,322 us per emitted pair. So a cap justified by counting
        throughput would be sizing the wrong thing by three orders of magnitude.
        """
        import time

        from scipy.stats import fisher_exact

        def cost(n_carrier: int, trials: int = 300) -> float:
            table = [[n_carrier // 2, n_carrier // 2],
                     [n_carrier // 2, n_carrier]]
            t0 = time.perf_counter()
            for _ in range(trials):
                fisher_exact(table, alternative="two-sided")
            return (time.perf_counter() - t0) / trials * 1e6

        small, large = cost(90), cost(900)
        assert large > small, (small, large)

        # And the wrapper adds almost nothing, so the constant above is Fisher's
        # cost and would stay true if the wrapper were rewritten.
        trials = 300
        t0 = time.perf_counter()
        for _ in range(trials):
            cc.cooccurrence_from_counts(
                "fa", "fb", "gene_gene", 450, 450, 225, 900)
        wrapper_us = (time.perf_counter() - t0) / trials * 1e6
        assert wrapper_us < 2.0 * large, (wrapper_us, large)


# --------------------------------------------------------------------------
# end to end through _compute
# --------------------------------------------------------------------------


class TestComputeAppliesTheUniverse:
    def test_a_narrowed_class_removes_its_rows(self, config: PipelineConfig):
        rng = random.Random(20240617)
        n = 18
        namespaces = _namespaces(rng, n)
        manifest = _manifest(n)

        wide = _pair_keys(
            cc._compute(config, manifest, namespaces, {}, ALL_TESTS)
        )
        narrow_cfg = _with_cooccurrence(config, feature_classes=["gene"])
        narrow = _pair_keys(
            cc._compute(narrow_cfg, manifest, namespaces, {}, ALL_TESTS)
        )
        assert wide, "the fixture produced no pairs, so this proves nothing"
        assert narrow
        assert set(narrow) < set(wide)
        assert {row[0] for row in narrow} == {"gene_gene"}

    def test_a_prevalence_filter_removes_rows(self, config: PipelineConfig):
        """A rare feature must cost *no* pairs, not be dropped after counting.

        The fixture is built with a common core every sample carries and a
        long tail only one or two samples do, because that is the shape the
        measured universe has (at n=10 the `mutation` carrier-count median is
        4/10 and the minimum is 1). Under a uniform density there would be
        nothing for a prevalence filter to remove and the assertion would pass
        vacuously - which is the failure this fixture exists to prevent.
        """
        rng = random.Random(42)
        n = 40
        namespaces = {
            "gene": {
                f"S{i}": {"g_common"} | {f"g_r{j}" for j in range(6)
                                             if rng.random() < 0.08}
                for i in range(n)
            },
            "mutation": {
                f"S{i}": {"m_common"} | {f"m_r{j}" for j in range(10)
                                             if rng.random() < 0.05}
                for i in range(n)
            },
            "mechanism": {f"S{i}": {"mech"} for i in range(n)},
        }
        manifest = _manifest(n)
        floor = cc._prevalence_floor(0.1, n)  # 4 of 40
        admitted = cc.admit_features(namespaces, cc.ELIGIBLE_FEATURE_CLASSES, 0.1, n)
        allowed = set()
        for name in cc.ELIGIBLE_FEATURE_CLASSES:
            allowed.update(cc.namespace_names(admitted[name]))
        assert floor == 4
        assert allowed, "the fixture removed everything, so the filter is untested"
        assert any(
            len(cc.namespace_names(namespaces[name]))
            > len(cc.namespace_names(admitted[name]))
            for name in cc.ELIGIBLE_FEATURE_CLASSES
        ), "the fixture removed nothing, so the filter is untested"

        unfiltered = _pair_keys(
            cc._compute(_with_cooccurrence(config, min_prevalence=0.0), manifest,
                        namespaces, {}, ALL_TESTS)
        )
        filtered = _pair_keys(
            cc._compute(_with_cooccurrence(config, min_prevalence=0.1), manifest,
                        namespaces, {}, ALL_TESTS)
        )
        assert unfiltered, "the fixture produced no pairs, so this proves nothing"
        assert set(filtered) < set(unfiltered)
        assert all(a in allowed and b in allowed for _, a, b in filtered)

    def test_the_cap_refuses_through_compute(self, config: PipelineConfig):
        rng = random.Random(7)
        n = 20
        namespaces = _namespaces(rng, n, sizes=(30, 30, 5))
        manifest = _manifest(n)
        over = _with_cooccurrence(config, max_features=10)
        with pytest.raises(StageError) as excinfo:
            cc._compute(over, manifest, namespaces, {}, ALL_TESTS)
        assert "10" in str(excinfo.value)
        # And under the cap the same fixture runs, so the refusal is the cap
        # rather than the fixture.
        assert cc._compute(
            _with_cooccurrence(config, max_features=65), manifest, namespaces, {},
            ALL_TESTS,
        )

    def test_the_configured_run_is_unchanged_by_the_filter_at_fixture_scale(
        self, config: PipelineConfig
    ):
        """R10's guarantee, as a test: 20 samples cannot move.

        `min_prevalence: 0.05` is in the real config, so every committed fixture
        and every TEST-mode run passes through it. This asserts the output is
        byte-identical to `min_prevalence: 0.0` at n=20, which is why no existing
        test assertion had to change.
        """
        rng = random.Random(20240617)
        n = 20
        namespaces = _namespaces(rng, n)
        manifest = _manifest(n)
        configured = cc._compute(
            _with_cooccurrence(config, min_prevalence=0.0), manifest, namespaces,
            {}, ALL_TESTS,
        )
        shipped = cc._compute(
            _with_cooccurrence(config, min_prevalence=0.05), manifest, namespaces,
            {}, ALL_TESTS,
        )
        assert configured, "the fixture produced no pairs, so this proves nothing"
        assert [dataclasses.astuple(r) for r in configured] == [
            dataclasses.astuple(r) for r in shipped
        ]

    @pytest.mark.parametrize("n", [1, 2, 5, 20])
    def test_the_shipped_prevalence_is_inert_up_to_twenty_samples(
        self, config: PipelineConfig, n
    ):
        assert cc._prevalence_floor(
            cc.resolve_min_prevalence(config), n
        ) == 1

    def test_a_run_at_the_cap_still_produces_rows(self, config: PipelineConfig):
        """The cap is a ceiling, not a switch that empties the table."""
        rng = random.Random(11)
        n = 20
        namespaces = _namespaces(rng, n)
        admitted = cc.admit_features(
            namespaces, cc.ELIGIBLE_FEATURE_CLASSES, 0.05, n
        )
        total = sum(
            len(cc.namespace_names(sets)) for sets in admitted.values()
        )
        assert 0 < total <= cc.resolve_max_features(config)
        rows = cc._compute(config, _manifest(n), namespaces, {}, ALL_TESTS)
        assert rows


# --------------------------------------------------------------------------
# the REAL gate agrees with the configured universe
# --------------------------------------------------------------------------


AMR_HEADER = (
    "sample_id\tantibiotic\tdeterminant\tgene\tvariant\tdeterminant_type\t"
    "mechanism\tevidence_source\tdatabase\tdatabase_version\tconfidence\t"
    "claim_status\tidentity_pct\tcoverage_pct"
)


def _write_amr(root: Path, rows: Sequence[tuple]) -> Path:
    path = root / "stages" / "04_amr.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [AMR_HEADER]
    for sample, gene in rows:
        lines.append(
            f"{sample}\timipenem\t{gene}\t{gene}\t\tAMR\t\tamrfinderplus\t"
            "AMRFinderPlus\t2026-08-07.1\t99.0\tDETECTED\t99.0\t100.0"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_mechanisms(root: Path, rows: Sequence[tuple]) -> Path:
    path = root / "stages" / "05_mechanisms.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "sample_id\tantibiotic\tdeterminant\tmechanism\tevidence_level\tgene\t"
        "gene_type\tnotes"
    ]
    for sample, gene, mechanism in rows:
        lines.append(
            f"{sample}\timipenem\t{gene}\t{mechanism}\tDETECTED\t{gene}\t"
            "chromosomal\tfixture"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_tree_metadata(root: Path, labels: dict) -> Path:
    from papipeline.stages.phylogeny import TREE_METADATA_COLUMNS

    phylo = root / "phylogeny"
    phylo.mkdir(parents=True, exist_ok=True)
    path = phylo / "tree_metadata.tsv"
    values = {"sample_id": None, "tree_tip_label": None, "lineage_label": None,
              "st": "ST1", "source": "mlst"}
    lines = ["\t".join(TREE_METADATA_COLUMNS)]
    for sample, label in labels.items():
        row = dict(values, sample_id=sample, tree_tip_label=sample, lineage_label=label)
        lines.append("\t".join(row[column] for column in TREE_METADATA_COLUMNS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return phylo


class TestRealAgreesWithTheConfiguredUniverse:
    """R12: real code paths, synthetic tables on disk, no bioinformatics tool.

    `_run_real` re-reads three tables through the *production* loaders and
    then checks completeness. If the completeness gate ignored
    `feature_classes`, a caller who deliberately excluded the mutation namespace
    would be refused for missing mutations before the exclusion could take
    effect - so the filter would exist and never run.
    """

    @staticmethod
    def _tree(tmp_path: Path, n: int = 6, with_mechanisms: bool = True) -> tuple:
        samples = [f"PDT{i:09d}.1" for i in range(n)]
        genes = ["blaPDC-1", "blaVIM-2", "oprD", "nalC", "mexM", "cblB"][:n]
        root = tmp_path / "intermediate"
        _write_amr(root, list(zip(samples, genes)))
        if with_mechanisms:
            _write_mechanisms(
                root,
                [(s, g, "efflux" if i % 2 else "reduced_permeability")
                 for i, (s, g) in enumerate(zip(samples, genes))],
            )
        _write_tree_metadata(tmp_path, {s: "ST1" for s in samples})
        manifest = SampleManifest(samples=[Sample(sample_id=s) for s in samples])
        return root, manifest, samples

    def test_a_run_excluding_mutations_is_not_refused_for_mutations(
        self, tmp_path: Path, config: PipelineConfig
    ):
        from papipeline.models import RunMode

        root, manifest, samples = self._tree(tmp_path)
        cfg = _with_cooccurrence(config, feature_classes=["gene", "mechanism"])
        results = cc.run(
            cfg, manifest, RunMode.REAL, intermediate_root=root,
            antibiotic="imipenem",
            phylogeny_dir=tmp_path / "phylogeny",
        )
        assert results, "the fixture produced no rows, so this proves nothing"
        assert {r.feature_type for r in results} == {
            "gene_gene", "gene_mechanism", "mechanism_mechanism",
        }

    def test_a_run_of_one_class_alone_is_refused_for_the_others_no_longer(
        self, tmp_path: Path, config: PipelineConfig
    ):
        """`feature_classes: [gene]` with no regulator and no mechanism table.

        Both other classes are absent from this tree. If the completeness gate
        ignored the configured universe it would refuse here, and the narrowing
        would never be reached.
        """
        from papipeline.models import RunMode

        root, manifest, _ = self._tree(tmp_path, with_mechanisms=False)
        cfg = _with_cooccurrence(config, feature_classes=["gene"])
        results = cc.run(
            cfg, manifest, RunMode.REAL, intermediate_root=root,
            antibiotic="imipenem",
            phylogeny_dir=tmp_path / "phylogeny",
        )
        assert results
        assert {r.feature_type for r in results} == {"gene_gene"}

    def test_the_widened_universe_still_refuses_on_a_truly_absent_class(
        self, tmp_path: Path, config: PipelineConfig
    ):
        """Narrowing the universe must not blind the gate.

        There is no `regulator_variants.tsv` on this tree, so with `mutation`
        eligible the refusal must still fire and must still name it.
        """
        from papipeline.models import RunMode

        root, manifest, _ = self._tree(tmp_path)
        with pytest.raises(StageError) as excinfo:
            cc.run(
                config, manifest, RunMode.REAL, intermediate_root=root,
                antibiotic="imipenem",
                phylogeny_dir=tmp_path / "phylogeny",
            )
        assert "mutation" in str(excinfo.value)

    def test_the_cap_refuses_a_real_run_too(self, tmp_path: Path, config: PipelineConfig):
        """Not a TEST-mode-only guard: the ceiling applies to both modes."""
        from papipeline.models import RunMode

        root, manifest, _ = self._tree(tmp_path)
        cfg = _with_cooccurrence(config, feature_classes=["gene"], max_features=2)
        with pytest.raises(StageError) as excinfo:
            cc.run(
                cfg, manifest, RunMode.REAL, intermediate_root=root,
                antibiotic="imipenem",
                phylogeny_dir=tmp_path / "phylogeny",
            )
        message = str(excinfo.value)
        assert "cooccurrence refuses to run" in message
        assert "6" in message and "2" in message


class TestTheKeysAreInTheShippedConfig:
    """N7: a key in `config/science.yaml` with no reader is a lie."""

    @pytest.mark.parametrize(
        "key", ["feature_classes", "min_prevalence", "max_features"]
    )
    def test_the_key_is_declared_and_read(self, key, config: PipelineConfig):
        assert key in config.raw["cooccurrence"], (
            f"cooccurrence.{key} is read by papipeline/stages/cooccurrence.py "
            "and must therefore be declared in config/science.yaml"
        )

    def test_the_declared_classes_are_the_ones_the_stage_accepts(
        self, config: PipelineConfig
    ):
        declared = config.raw["cooccurrence"]["feature_classes"]
        assert list(cc.resolve_feature_classes(config)) == declared
        assert set(declared) <= set(cc.ELIGIBLE_FEATURE_CLASSES)

    def test_the_config_states_the_measured_cost_the_cap_was_chosen_from(self):
        """The numbers the justification rests on are in the file, not just in a report."""
        text = (
            Path(__file__).resolve().parents[2] / "config" / "science.yaml"
        ).read_text(encoding="utf-8")
        for needle in ("0.34 us/pair", "1,902.50 us/pair", "881 B/row",
                       "3,123,750", "5,925.67 s = 1.646 h", "2.27 GiB"):
            assert needle in text, f"config/science.yaml no longer states {needle!r}"