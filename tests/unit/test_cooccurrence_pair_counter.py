"""The blocked matrix product must agree with the naive pair loop, exactly.

`cooccurrence` counted C(F,2) pairs by rescanning the cohort for each one:
`test_pair` built `carriers_a` and `carriers_b` with a set comprehension over
every sample, then intersected them. That is O(F^2 * n) Python-level work, and
at F = 50,000 it is 1.25e9 pairs - the reason SCALE ranked stage 14 as a
harder wall than stage 2's 80 hours.

The replacement is `matrix_product_pair_counts`: `n_both(a, b)` is the dot
product of the two carrier vectors, so the table is `M.T @ M`, walked in
blocks. This file is the evidence that it is the *same* arithmetic.

Three things are pinned, and the third is the one that would have rotted
silently:

1. **Triple-for-triple agreement, property-based over several seeds.** Not one
   hand-built example. The naive loop is kept in the module precisely so it can
   be run beside the fast path on random inputs; a single fixture would agree
   just as well if the block walk skipped a whole diagonal, which is the bug
   that costs the most and shows up least.
2. **Chunking is invisible.** The same input must give the same triples for
   every chunk edge, including 1, which forces one feature per block and
   exercises every diagonal-block branch. If the walk reordered or dropped the
   lower triangle, a single chunk size would hide it.
3. **The order is unchanged.** Not because any number depends on it - BH sorts
   by p-value - but because the order is the order of the written table, so a
   change shows up as a diff nobody can explain.

Both counters are also asserted to produce the *same* `_compute` output, since
`_compute` is what run.py actually calls and a counter that agrees in
isolation could still be wired to the wrong arguments.
"""

from __future__ import annotations

import dataclasses
import random
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.manifest import SampleManifest
from papipeline.models import Sample
from papipeline.stages import cooccurrence as cc

#: Seeds, not one. 0 is in the list because a fixed seed that happens to work
#: is indistinguishable from a fixed seed that happens to be unrepresentative.
SEEDS = (0, 1, 7, 42, 1337, 20240617)


def _random_namespace(
    rng: random.Random, n_samples: int, n_features: int, density: float
) -> Dict[str, Set[str]]:
    return {
        f"S{index}": {
            f"f{feature}"
            for feature in range(n_features)
            if rng.random() < density
        }
        for index in range(n_samples)
    }


def _sample_ids(n: int) -> List[str]:
    return [f"S{index}" for index in range(n)]


def _pairs(
    counter, sample_ids, set_a, set_b, same_namespace, **kwargs
) -> List[Tuple[str, str, int, int, int]]:
    """The counter's triples, in ITS order.

    Returned unsorted on purpose. The blocked walk is tile-major, so its order
    differs from `all_pairs`' and cannot be made to match without buffering an
    a-block's worth of pairs - see `matrix_product_pair_counts`. Order is
    therefore asserted where it is actually guaranteed, in `_compute`, not here.
    """
    return list(
        counter(
            sample_ids, set_a, set_b, same_namespace=same_namespace, **kwargs
        )
    )


def _same_triples(left, right) -> None:
    """Set equality on the triples, with the multiset length checked too.

    `Counter` rather than `set` so a duplicated triple fails: a walk that
    emitted one pair twice and dropped another would compare equal as a set.
    """
    assert Counter(left) == Counter(right), (
        f"{len(left)} vs {len(right)} triples; "
        f"only in the first: {sorted(Counter(left) - Counter(right))[:5]}; "
        f"only in the second: {sorted(Counter(right) - Counter(left))[:5]}"
    )


class TestTheCountersAgree:
    @pytest.mark.parametrize("seed", SEEDS)
    def test_within_namespace(self, seed):
        """The C(F,2) case, including a sparse one where most features are absent."""
        rng = random.Random(seed)
        n_samples, n_features = 24, 11
        for density in (0.15, 0.5, 0.85):
            namespace = _random_namespace(rng, n_samples, n_features, density)
            sample_ids = _sample_ids(n_samples)
            _same_triples(
                _pairs(
                    cc.naive_pair_counts, sample_ids, namespace, namespace, True
                ),
                _pairs(
                    cc.matrix_product_pair_counts,
                    sample_ids,
                    namespace,
                    namespace,
                    True,
                    chunk_features=4,
                ),
            )

    @pytest.mark.parametrize("seed", SEEDS)
    def test_cross_namespace(self, seed):
        """The full cross product, where the two feature sets differ in size."""
        rng = random.Random(seed)
        n_samples = 19
        genes = _random_namespace(rng, n_samples, 7, 0.4)
        mutations = _random_namespace(rng, n_samples, 13, 0.25)
        sample_ids = _sample_ids(n_samples)
        _same_triples(
            _pairs(
                cc.naive_pair_counts, sample_ids, genes, mutations, False
            ),
            _pairs(
                cc.matrix_product_pair_counts,
                sample_ids,
                genes,
                mutations,
                False,
                chunk_features=3,
            ),
        )

    @pytest.mark.parametrize("seed", SEEDS[:3])
    def test_a_feature_absent_from_the_cohort_still_yields_its_pairs(self, seed):
        """A name carried only by a non-manifest sample reaches both counters.

        `namespace_names` ranges over every key of the mapping, not the
        manifest, so such a feature produces a pair with `n_a == 0`. That pair
        is dropped later as degenerate - but it must be *counted* identically
        first, or the two counters disagree on a case the production path hits
        whenever a stage hands over a sample the manifest does not list.
        """
        rng = random.Random(seed)
        namespace = _random_namespace(rng, 12, 6, 0.4)
        # A name NO cohort sample carries, so n_a == 0 for every pair it makes.
        namespace["NOT_IN_MANIFEST"] = {"f_outsider"}
        sample_ids = _sample_ids(12)
        naive = _pairs(cc.naive_pair_counts, sample_ids, namespace, namespace, True)
        fast = _pairs(
            cc.matrix_product_pair_counts,
            sample_ids,
            namespace,
            namespace,
            True,
            chunk_features=2,
        )
        _same_triples(naive, fast)
        outsider_pairs = [t for t in naive if "f_outsider" in (t[0], t[1])]
        assert outsider_pairs, "the outsider produced no pair at all"
        # "f_outsider" sorts last, so it is always feature_b and its count is 0.
        assert all(triple[3] == 0 for triple in outsider_pairs)

    @pytest.mark.parametrize("chunk", [1, 2, 3, 7, 64])
    def test_the_chunk_edge_is_invisible(self, chunk):
        """One feature per block is the case that finds a dropped diagonal."""
        rng = random.Random(99)
        namespace = _random_namespace(rng, 15, 9, 0.45)
        sample_ids = _sample_ids(15)
        reference = _pairs(
            cc.naive_pair_counts, sample_ids, namespace, namespace, True
        )
        assert reference
        _same_triples(
            _pairs(
                cc.matrix_product_pair_counts,
                sample_ids,
                namespace,
                namespace,
                True,
                chunk_features=chunk,
            ),
            reference,
        )

    def test_the_blocked_walk_skips_the_lower_triangle_rather_than_the_upper(self):
        """A guard on the guard: no pair is emitted twice, and none reversed."""
        rng = random.Random(5)
        namespace = _random_namespace(rng, 10, 6, 0.5)
        sample_ids = _sample_ids(10)
        pairs = _pairs(
            cc.matrix_product_pair_counts,
            sample_ids,
            namespace,
            namespace,
            True,
            chunk_features=2,
        )
        keys = [(a, b) for a, b, *_ in pairs]
        assert len(keys) == len(set(keys))
        assert all(a < b for a, b in keys)

    def test_an_empty_namespace_yields_nothing_rather_than_raising(self):
        assert _pairs(cc.matrix_product_pair_counts, ["S1"], {}, {}, True) == []
        assert (
            _pairs(
                cc.matrix_product_pair_counts, ["S1"], {"S1": {"f0"}}, {}, False
            )
            == []
        )

    def test_a_chunk_edge_below_one_is_refused_with_the_key_named(self):
        with pytest.raises(ValueError) as excinfo:
            list(
                cc.matrix_product_pair_counts(
                    ["S1"], {"S1": {"a"}}, {"S1": {"a"}},
                    same_namespace=True, chunk_features=0,
                )
            )
        assert "count_chunk_features" in str(excinfo.value)


class TestComputeIsIdenticalWhicheverCounterRuns:
    """`_compute` is what run.py calls, so the counter has to be wired right."""

    @staticmethod
    def _manifest(n: int) -> SampleManifest:
        return SampleManifest(
            samples=[Sample(sample_id=f"S{i}") for i in range(n)]
        )

    @staticmethod
    def _namespaces(rng: random.Random, n: int):
        return {
            "gene": _random_namespace(rng, n, 8, 0.45),
            "mutation": _random_namespace(rng, n, 10, 0.3),
            "mechanism": _random_namespace(rng, n, 5, 0.6),
        }

    @pytest.mark.parametrize("seed", SEEDS[:3])
    def test_all_five_configured_tests_agree(self, config: PipelineConfig, seed):
        rng = random.Random(seed)
        n = 18
        namespaces = self._namespaces(rng, n)
        manifest = self._manifest(n)
        configured = [
            "genes",
            "gene_mutation",
            "mutation_mutation",
            "gene_mechanism",
            "mechanism_mechanism",
        ]

        def run_with(counter: str):
            raw = dict(config.raw)
            raw["cooccurrence"] = {
                **raw.get("cooccurrence", {}),
                "pair_counter": counter,
                "count_chunk_features": 3,
            }
            return cc._compute(
                dataclasses.replace(config, raw=raw), manifest, namespaces, {},
                configured,
            )

        naive = run_with("naive")
        fast = run_with("matrix_product")
        assert naive, "the fixture produced no pairs, so this proves nothing"
        assert [_row(r) for r in naive] == [_row(r) for r in fast]

    def test_the_table_order_does_not_depend_on_the_counter(self, config):
        """`_compute` normalises order, so the written table does not move.

        Not cosmetic: `14_cooccurrence.tsv` is a contract artefact, and a row
        order that changed with `cooccurrence.pair_counter` or with
        `count_chunk_features` would make every diff of that file unreadable and
        would make two runs of the same cohort incomparable.
        """
        rng = random.Random(20240617)
        n = 16
        namespaces = {
            "gene": _random_namespace(rng, n, 7, 0.4),
            "mutation": _random_namespace(rng, n, 9, 0.3),
            "mechanism": _random_namespace(rng, n, 4, 0.55),
        }
        manifest = self._manifest(n)
        configured = ["genes", "gene_mutation", "gene_mechanism"]

        def order_of(counter: str, chunk: int):
            raw = dict(config.raw)
            raw["cooccurrence"] = {
                **raw.get("cooccurrence", {}),
                "pair_counter": counter,
                "count_chunk_features": chunk,
            }
            results = cc._compute(
                dataclasses.replace(config, raw=raw), manifest, namespaces, {},
                configured,
            )
            return [(r.feature_type, r.feature_a, r.feature_b) for r in results]

        naive = order_of("naive", 512)
        assert naive == order_of("matrix_product", 1)
        assert naive == order_of("matrix_product", 512)
        assert naive == order_of("matrix_product", 4096)
        # And it is sorted within each family, which is what the old
        # all_pairs/cross-product order produced.
        for family in {row[0] for row in naive}:
            keys = [row[1:] for row in naive if row[0] == family]
            assert keys == sorted(keys)

    def test_the_default_counter_is_the_matrix_product(self, config: PipelineConfig):
        assert cc.resolve_pair_counter(config) == "matrix_product"

    def test_an_unrecognised_counter_is_refused_rather_than_defaulted(
        self, config: PipelineConfig
    ):
        """A reader who asked for the oracle must not silently get production."""
        raw = dict(config.raw)
        raw["cooccurrence"] = {
            **raw.get("cooccurrence", {}),
            "pair_counter": "bitset",
        }
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_pair_counter(dataclasses.replace(config, raw=raw))
        assert "matrix_product" in str(excinfo.value)

    @pytest.mark.parametrize("bad", [0, -1])
    def test_a_non_positive_chunk_is_refused(self, config: PipelineConfig, bad):
        raw = dict(config.raw)
        raw["cooccurrence"] = {
            **raw.get("cooccurrence", {}),
            "count_chunk_features": bad,
        }
        with pytest.raises(ValueError) as excinfo:
            cc.resolve_count_chunk_features(dataclasses.replace(config, raw=raw))
        assert "count_chunk_features" in str(excinfo.value)

    def test_a_non_integer_chunk_is_refused(self, config: PipelineConfig):
        raw = dict(config.raw)
        raw["cooccurrence"] = {
            **raw.get("cooccurrence", {}),
            "count_chunk_features": "wide",
        }
        with pytest.raises(ValueError):
            cc.resolve_count_chunk_features(dataclasses.replace(config, raw=raw))

    def test_an_absent_chunk_key_falls_back_to_the_pinned_default(
        self, config: PipelineConfig
    ):
        raw = dict(config.raw)
        raw["cooccurrence"] = {
            k: v
            for k, v in (raw.get("cooccurrence", {}) or {}).items()
            if k != "count_chunk_features"
        }
        assert (
            cc.resolve_count_chunk_features(dataclasses.replace(config, raw=raw))
            == cc.DEFAULT_COUNT_CHUNK_FEATURES
            == 512
        )


def _row(result) -> Tuple:
    """Everything a row carries, so a difference anywhere shows up."""
    return (
        result.feature_a,
        result.feature_b,
        result.feature_type,
        result.n_a,
        result.n_b,
        result.n_both,
        result.statistic,
        result.statistic_value,
        result.adjusted_p_value,
    )


# ---------------------------------------------------------------------------
# The two new config keys, asserted the way round 11's are
# ---------------------------------------------------------------------------


class TestTheConfigKeysAreReal:
    """A key with only a YAML value is a typo; a key with only a default is a
    comment. Each is asserted to *reach* something.

    Rule 2 is why `count_chunk_features` exists at all: a memory bound written
    into a stage is a memory bound nobody can change without editing a stage.
    """

    def test_the_science_file_states_both(self, config: PipelineConfig):
        section = config.raw["cooccurrence"]
        assert section["pair_counter"] == "matrix_product"
        assert int(section["count_chunk_features"]) == 512

    def test_the_values_reach_the_counter(self, config: PipelineConfig):
        """Not merely present: `pair_counter` selects the implementation."""
        sample_ids = ["S1", "S2"]
        namespace = {"S1": {"a"}, "S2": {"b"}}
        def which(counter_name: str):
            raw = dict(config.raw)
            raw["cooccurrence"] = {
                **raw.get("cooccurrence", {}),
                "pair_counter": counter_name,
            }
            stream = cc.pair_counter(
                dataclasses.replace(config, raw=raw), sample_ids, namespace,
                namespace, same_namespace=True,
            )
            # The generator's own code object names the implementation, so this
            # asserts the key selected a *different function*, not just a
            # different answer.
            return stream.gi_code.co_name, list(stream)

        assert which("naive")[0] == "naive_pair_counts"
        assert which("matrix_product")[0] == "matrix_product_pair_counts"
        # Same answers either way - the config switches the route, not the answer.
        assert which("naive")[1] == sorted(which("matrix_product")[1])

    def test_the_chunk_key_reaches_the_block_walk(self, config: PipelineConfig):
        """A chunk edge of 1 forces one feature per block.

        Asserted through the answer rather than through a mock, because a mock
        would pass even if the key were threaded to the wrong argument.
        """
        sample_ids = [f"S{i}" for i in range(8)]
        # Deterministic, not `hash()`: string hashing is salted per process, so
        # a hash-derived fixture changes shape between runs and a benchmark of
        # "the pairs this produces" would not be reproducible.
        namespace = {
            f"S{index}": {f"f{feature}" for feature in range(8)
                          if (index >> feature) & 1 and feature != index}
            for index in range(8)
        }
        raw = dict(config.raw)
        raw["cooccurrence"] = {
            **raw.get("cooccurrence", {}),
            "count_chunk_features": 1,
        }
        tight = list(
            cc.pair_counter(
                dataclasses.replace(config, raw=raw), sample_ids, namespace,
                namespace, same_namespace=True,
            )
        )
        raw["cooccurrence"]["count_chunk_features"] = 64
        loose = list(
            cc.pair_counter(
                dataclasses.replace(config, raw=raw), sample_ids, namespace,
                namespace, same_namespace=True,
            )
        )
        assert tight, "the fixture produced no pairs"
        assert Counter(tight) == Counter(loose)

    def test_min_cell_count_is_still_dead_and_says_so(self, config):
        """Recorded rather than fixed, and pinned so it cannot be fixed quietly.

        `cooccurrence.min_cell_count: 2` is the contingency-table filter that
        would decide whether this stage is runnable at F=50,000, and nothing
        reads it. Wiring it changes which rows the stage emits, which is a
        behaviour change for Phase 4 rather than a bug fix for this round - but
        a key that is dead *and* undocumented is how the next reader misses it.
        """
        science = (
            Path(__file__).resolve().parents[2] / "config" / "science.yaml"
        ).read_text(encoding="utf-8")
        assert "min_cell_count: 2" in science
        assert "DEAD" in science

        stages = Path(cc.__file__).parent
        # A *read*, not a mention: the module docstring names the key on
        # purpose, so the search is for a subscript or a `.get`.
        readers = [
            path.name
            for path in sorted(stages.glob("*.py"))
            if any(
                token in path.read_text(encoding="utf-8")
                for token in ('["min_cell_count"]', '.get("min_cell_count")',
                              "'min_cell_count'")
            )
        ]
        assert readers == [], (
            f"{readers} now read cooccurrence.min_cell_count. That is the "
            "Phase 4 behaviour change, not a fix - land it as its own commit "
            "with the row-count delta it causes recorded."
        )


# ---------------------------------------------------------------------------
# O4: pooled counts only, and the module says so
# ---------------------------------------------------------------------------


class TestThePoolingClaimIsHonest:
    """O4: the claim at the module docstring, re-checked rather than assumed.

    Commit `0f92078` removed the sentence *"Lineage-stratified counts are
    reported alongside the pooled counts"* and replaced it with an explicit
    statement that nothing is stratified. That removal is load-bearing, so it is
    pinned here: a docstring can drift back to a promise the code does not
    keep, and this is the promise most likely to be restored by someone adding
    stratification halfway.

    The alternative - actually computing stratified counts - was weighed and
    declined. See `docs` in the round-12 COOC.md; the three reasons are that
    `COOCCURRENCE_COLUMNS` and `contracts.STAGE_TABLES` declare ten columns with
    no lineage field, so a stratified count has nowhere to go without a change
    to `run.py` that this round does not own; that BH would then have to cover
    the pooled family *and* one family per lineage, multiplying m by
    (1 + n_lineages) at exactly the point where the cohort is too small to
    afford it; and that at n=10 an MLST-ST stratum holds one to three samples,
    where a 2x2 Fisher test is arithmetically defined and scientifically
    meaningless.
    """

    def test_the_old_promise_is_not_in_the_docstring(self):
        source = Path(cc.__file__).read_text(encoding="utf-8")
        assert "reported alongside the pooled counts" not in source

    def test_the_docstring_states_the_pooling_and_the_reason(self):
        source = Path(cc.__file__).read_text(encoding="utf-8")
        assert "pooled across the cohort only" in source
        assert "COOCCURRENCE_COLUMNS`` carries no lineage field" in source

    def test_no_lineage_field_exists_to_receive_a_stratified_count(self):
        """The structural reason stratification is not a message-level fix."""
        from papipeline.execution.contracts import STAGE_TABLES

        assert not any("lineage" in c for c in cc.COOCCURRENCE_COLUMNS)
        _name, columns = STAGE_TABLES["cooccurrence"]
        assert not any("lineage" in c for c in columns)
        assert tuple(columns) == cc.COOCCURRENCE_COLUMNS

    def test_the_crosswalk_is_still_refused_on_because_it_is_unused(self):
        """A pin on the trade this round did *not* make.

        `_add_missing_lineages` refuses REAL on a degenerate crosswalk even
        though no test reads it. That is deliberate (commit `0f92078`) and it
        costs a refusal without changing an output value - the reason it is
        still open rather than settled is recorded in the round-12 backlog.
        """
        assert cc._add_missing_lineages({}, {"S1": "unknown"}) == {
            "lineages": "the sample -> lineage crosswalk carries no real lineage label"
        }
        assert cc._add_missing_lineages({}, {"S1": "ST1"}) == {}

    def test_a_real_crosswalk_lets_the_stage_through(self):
        """So the refusal is not unconditional - R4 makes this reachable."""
        assert cc._add_missing_lineages({}, {"S1": "ST1", "S2": "ST2"}) == {}
