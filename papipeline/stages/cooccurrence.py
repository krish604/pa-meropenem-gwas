"""Stage 14 - co-occurrence analysis.

Tests whether pairs of features (gene x gene, gene x mutation,
mutation x mutation, gene x mechanism, mechanism x mechanism) co-occur more
often than expected by chance.

The statistical framing matters more than the arithmetic. Each test reports
the odds ratio, Fisher's exact p-value and a BH-adjusted p-value, plus the
counts needed to recompute it. Every row carries
``interpretation_limit=association_only_not_causal``, because:

* co-occurrence in a cross-sectional cohort cannot establish that one
  feature influenced another;
* the tests are not adjusted for lineage or for clonal structure, so pairs
  that co-occur because both belong to the same clone will look associated;
* the comparison is against the cohort's own marginal frequencies, which
  absorbs background linkage.

Counts are **pooled across the cohort only**. Nothing here is stratified by
lineage: ``COOCCURRENCE_COLUMNS`` carries no lineage field and no test reads a
lineage label, so a row states how two features co-occur *in this cohort* and
says nothing about whether an association survives within a lineage. That is
the same limitation the ``interpretation_limit`` column already records for
causality, stated here for lineage.

The sample -> lineage crosswalk is still a declared input, but only as a
completeness gate: REAL refuses rather than publish pooled counts beside a
crosswalk a reader might assume stratified them (see
``_add_missing_lineages``).

The feature universe: DETERMINANTS ONLY
--------------------------------------
This stage pairs three feature classes and nothing else. Measured on the real
10-isolate tree, with each name counted the way its loader spells it:

=====================  =========================  ==============================
class                  read from                  distinct names at n=10
=====================  =========================  ==============================
``gene``               ``04_amr.tsv``             **54** acquired AMR genes
``mutation``           ``regulator_variants.tsv`` **300** ``{gene}:{variant}``
                                                   over 10 regulator loci
``mechanism``          ``05_mechanisms.tsv``      **4**, a closed vocabulary
=====================  =========================  ==============================

**There is no gene presence/absence here, and adding one would be a different
stage.** Stage 7 writes ``pangenome/gene_presence_absence.tsv`` and nothing in
this module - no loader, no config key, no path - reads it. A pangenome
universe would be the ten-thousand-feature case rather than the
determinant-only one, and reaching it is a decision about what question stage 14
answers, not a tuning change.

Two consequences of the three classes above, both measured rather than assumed:

* the ``gene`` class is **acquired** determinants only. ``load_amr_names``
  keeps ``determinant_type == "AMR"`` and drops chromosomal and intrinsic calls,
  so its size is bounded by the AMR gene content of a *P. aeruginosa* isolate
  (tens of genes), not by the ~6,000 coding sequences of the reference.
* the ``mechanism`` class is **closed, and the ceiling is five names**: the four
  distinct mechanisms in ``config/mechanisms.tsv`` plus ``locus_intact``, which
  ``mechanisms.py:301`` writes directly rather than from that table. The real
  10-isolate tree shows four of them, because no ``locus_intact`` row survives
  once oprD carries a variant (see the round-12 report). So
  C(5,2) = 10 pairs, ever. It cannot be the wall and is not treated as one.

Only ``mutation`` grows without bound, because a cohort has more isolates in
which to find a variant. Projected from the measured cumulative curve (13.0 new
``{gene}:{variant}`` names per isolate over n=5..10, extrapolated linearly
because the true curve saturates and a linear projection is therefore an upper
bound): **11,870 at n=900, against 300 at n=10.** ``gene`` projects to 2,368.

Scale
-----
Counting the pairs is the part that was measured and fixed. Every feature is
present or absent per sample, so the whole counting problem is one matrix
product, and :func:`matrix_product_pair_counts` does it as ``M.T @ M`` in
blocks of ``cooccurrence.count_chunk_features`` rather than as a Python loop
over C(F, 2) pairs. :func:`naive_pair_counts` is the same arithmetic written
the slow way and is kept as the **oracle**: ``pair_counter: naive`` selects it,
and a property-based test asserts the two agree triple-for-triple. Measured at
n=900: **0.34 us per pair counted**, against a projected ~260 us for the Python
pair loop, so counting is no longer where the time goes.

What the matrix product does **not** fix is where the stage still dies, and it
is stated here with its measurement rather than glossed. On the same n=900
cohort, one **admitted** pair costs **1,902 us** downstream of counting - one
``scipy.stats.fisher_exact`` call on a 900-row contingency table, plus one
``Cooccurrence`` - against 0.34 us to count it. **The test is 5,600x the count.**
Each *surviving* row also costs **881 resident bytes**, measured as RSS delta
over emitted rows and including the ``pvalues`` key dict BH needs. 11.3% of
admitted pairs are dropped as degenerate and reach neither.

Those two figures are measured by a run at exactly this stage's cap - n=900,
2,500 features, 3,114,595 pairs admitted, 2,763,235 rows emitted, **5,925.67 s =
1.646 h wall and 2.27 GiB peak RSS**, end to end through this module's own
``_compute``.

The 1,902 us is worth its own sentence because it is the number most likely to be
got wrong: the *same call* on a 90-carrier table measures 128 us, and on a
900-row one ~500-650 us. Fisher's cost is a function of the table margins, not of
the code path, so a per-pair cost micro-benchmarked on a small table is not a
smaller version of the n=900 cost - it is a different one, and ``cProfile`` puts
``fisher_exact`` at **85% of this stage's wall** at n=900.

So the cost of a run is set by the *number of pairs*, and the number of pairs is
C(F, 2) in the feature universe above. That is what the remaining keys are for:

* ``cooccurrence.min_prevalence`` - a feature carried by fewer than this
  fraction of the cohort cannot be tested for co-occurrence with anything, and
  each one it removes removes F-1 pairs. At 5% it is **inert below n=21**, so
  the committed 20-sample fixtures are untouched and a 900-isolate run gets a
  real filter.
* ``cooccurrence.max_features`` - the hard ceiling on the admitted universe,
  summing the configured classes. Exceeding it **refuses**, naming the count and
  the cap, because the alternative is a run that spends 54 hours and 74 GiB
  producing a table nobody can read. It is set to 2,500 because that is the
  smallest round number above the *gene + mechanism* projection at n=900 (2,373):
  the one cohort-scale configuration that runs today is genes against mechanisms,
  and the mutation-inclusive analysis is refused with instructions instead of
  attempted. The full three-class n=900 universe projects to 14,242 features and
  101,410,161 pairs - about 54 h and 74 GiB at the measured rates - which is what
  the refusal exists to prevent.
* ``cooccurrence.feature_classes`` - which of the three classes are eligible, so
  a caller can turn a namespace off without editing this module.

``cooccurrence.min_cell_count`` is parsed by nothing and remains the obvious
*next* filter: it would drop pairs on the contingency table rather than features
on prevalence, and it changes which rows the stage emits. Recorded in the
round-12 backlog rather than wired here.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

from ..config.loader import PipelineConfig
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import Cooccurrence, RunMode

LOGGER = get_logger("stages.cooccurrence")

COOCCURRENCE_COLUMNS: Tuple[str, ...] = (
    "feature_a",
    "feature_b",
    "feature_type",
    "n_a",
    "n_b",
    "n_both",
    "statistic",
    "statistic_value",
    "adjusted_p_value",
    "interpretation_limit",
)

FeatureSets = Mapping[str, Set[str]]

#: One counted pair: ``(feature_a, feature_b, n_a, n_b, n_both)``.
PairCount = Tuple[str, str, int, int, int]

#: The two ways the pairs are counted. ``matrix_product`` is the production
#: path; ``naive`` is the oracle the property-based test compares against, and
#: is selectable from config so the oracle is reachable in a real run and not
#: only from a test.
PAIR_COUNTERS: Tuple[str, ...] = ("matrix_product", "naive")

#: Features per block edge in :func:`matrix_product_pair_counts`. Bounds the
#: transient to ``chunk x chunk`` float32, so peak memory is a function of this
#: key and not of the feature count. Read from
#: ``cooccurrence.count_chunk_features``; this is the fallback for a config that
#: omits it, and a test pins it.
DEFAULT_COUNT_CHUNK_FEATURES = 512

#: A float32 accumulator is exact while the running sum stays below 2**24, so a
#: cohort at or above that many samples must not be counted this way. 16.7
#: million isolates is not a reachable cohort, but silently losing precision
#: above a threshold nobody wrote down is the failure mode worth naming.
_FLOAT32_EXACT_LIMIT = 2 ** 24

#: The feature classes this stage can pair, and nothing else. The module
#: docstring carries the measurement that fixes this list: the universe is
#: determinants only, and stage 7's `gene_presence_absence.tsv` is not read
#: here. A name outside this tuple is a configuration error rather than a class
#: that quietly pairs nothing.
ELIGIBLE_FEATURE_CLASSES: Tuple[str, ...] = ("gene", "mutation", "mechanism")

#: Fallbacks when a config omits a key. ``None`` for the cap means *no cap*, and
#: that default is deliberate and is the reason the key has to be in
#: `config/science.yaml`: an uncapped stage is a stage whose cost nobody chose.
#: `min_prevalence` defaults to 0.0 (no filter) for the same reason the committed
#: fixtures must not move.
DEFAULT_MIN_PREVALENCE = 0.0


def build_feature_sets(
    manifest: SampleManifest,
    amr_genes: Mapping[str, Sequence[str]],
    variants: Mapping[str, Sequence[str]],
    mechanisms: Mapping[str, Sequence[str]],
) -> Dict[str, FeatureSets]:
    """Build the four feature namespaces used by the configured tests."""
    sample_ids = manifest.sample_ids
    genes: Dict[str, Set[str]] = {s: set() for s in sample_ids}
    mutations: Dict[str, Set[str]] = {s: set() for s in sample_ids}
    mech: Dict[str, Set[str]] = {s: set() for s in sample_ids}

    for sample_id, names in amr_genes.items():
        genes.setdefault(sample_id, set()).update(names)
    for sample_id, names in variants.items():
        mutations.setdefault(sample_id, set()).update(names)
    for sample_id, names in mechanisms.items():
        mech.setdefault(sample_id, set()).update(names)

    return {
        "gene": genes,
        "mutation": mutations,
        "mechanism": mech,
    }


def _fisher(odds_table: List[List[int]]) -> Tuple[Optional[float], Optional[float]]:
    from scipy.stats import fisher_exact

    odds_ratio, p_value = fisher_exact(odds_table, alternative="two-sided")
    if odds_ratio is not None and not (odds_ratio == float("inf") or odds_ratio == float("-inf")):
        odds: Optional[float] = float(odds_ratio)
    else:
        odds = None
    return odds, float(p_value)


def namespace_names(sets: FeatureSets) -> List[str]:
    """Every feature name in a namespace, sorted.

    Taken over *all* keys of ``sets``, not just the manifest's samples, so this
    matches what :func:`all_pairs` was fed before. A name that only a
    non-cohort sample carries still produces a pair here, and that pair has
    ``n_a == 0``, which :func:`cooccurrence_from_counts` drops as degenerate -
    the same outcome the naive path reached by scanning the manifest.
    """
    return sorted({name for members in sets.values() for name in members})


def cooccurrence_from_counts(
    feature_a: str,
    feature_b: str,
    feature_type: str,
    n_a: int,
    n_b: int,
    n_both: int,
    n_samples: int,
) -> Optional[Cooccurrence]:
    """The one place a counted triple becomes a tested pair.

    Every counting strategy funnels through here, so the statistics, the
    degenerate-pair rule and the ``Cooccurrence`` shape cannot differ between
    the matrix product and the naive oracle. Split out of :func:`test_pair` for
    exactly that reason; it is not a public entry point because it takes counts
    on trust.

    Returns ``None`` when a feature is fixed across the cohort, which is the
    only case with no contrast to test.
    """
    # A feature carried by every sample (or by none) cannot be tested for
    # association with anything.
    if min(n_a, n_b, n_samples - n_a, n_samples - n_b) == 0:
        LOGGER.debug(
            "Skipping degenerate pair %s x %s (a feature is fixed across the cohort)",
            feature_a,
            feature_b,
        )
        return None

    n_neither = n_samples - n_a - n_b + n_both
    table = [[n_both, n_a - n_both], [n_b - n_both, n_neither]]
    odds, p_value = _fisher(table)
    return Cooccurrence(
        feature_a=feature_a,
        feature_b=feature_b,
        feature_type=feature_type,
        n_a=n_a,
        n_b=n_b,
        n_both=n_both,
        statistic="odds_ratio",
        statistic_value=odds,
        adjusted_p_value=p_value,
    )


def _presence_matrix(
    sample_ids: Sequence[str],
    sets: FeatureSets,
    names: Sequence[str],
    index: Mapping[str, int],
):
    """``(n_samples, n_features)`` float32 of carrier indicators.

    float32 rather than uint8 because the product below accumulates in it, and
    an integer accumulator in numpy has no BLAS path. Exact while the cohort is
    under :data:`_FLOAT32_EXACT_LIMIT` samples, which
    :func:`matrix_product_pair_counts` checks.
    """
    import numpy as np

    matrix = np.zeros((len(sample_ids), len(names)), dtype=np.float32)
    for column, name in enumerate(names):
        for sample_id, members in sets.items():
            if name in members:
                row = index.get(sample_id)
                if row is not None:
                    matrix[row, column] = 1.0
    return matrix


def matrix_product_pair_counts(
    sample_ids: Sequence[str],
    set_a: FeatureSets,
    set_b: FeatureSets,
    *,
    same_namespace: bool,
    chunk_features: int = DEFAULT_COUNT_CHUNK_FEATURES,
) -> Iterator[PairCount]:
    """Count every pair with a blocked matrix product, not a Python pair loop.

    ``n_both(a, b)`` is the dot product of the two features' carrier vectors, so
    the whole C(F, 2) table is ``M.T @ M``. That table is never materialised: it
    is walked in ``chunk_features``-sided blocks in both indices, so the
    transient is ``chunk x chunk`` float32 and peak memory is set by
    ``chunk_features`` rather than by F.

    The blocks are walked a-major, then in blocks of ``feature_b`` inside that,
    so the pairs come out **tile-major**: for one a-block, every pair with a
    b-index below the diagonal is emitted before any pair above it. That is not
    lexicographic, and it cannot be made lexicographic without buffering a
    whole a-block's worth of pairs - 512 features x 50,000 b-features is 25.6e6
    triples, which is the memory wall this function exists to avoid.

    So the counter does not promise an order and :func:`_compute` normalises
    one: it sorts each test family's rows before concatenating, which is what
    ``all_pairs`` produced and therefore leaves the written table byte-identical
    whichever counter ran. That is the property worth having - the table must not
    depend on a memory-bound tuning key.
    """
    import numpy as np

    if chunk_features < 1:
        raise ValueError(
            f"cooccurrence.count_chunk_features is {chunk_features}; it is a "
            "block edge, so it must be at least 1."
        )
    if len(sample_ids) >= _FLOAT32_EXACT_LIMIT:
        LOGGER.warning(
            "%d samples is at or above the %d at which a float32 accumulator "
            "stops being exact; falling back to the naive counter. Unreachable "
            "at any cohort this pipeline has been run on, and named so nobody "
            "has to rediscover it.",
            len(sample_ids),
            _FLOAT32_EXACT_LIMIT,
        )
        yield from naive_pair_counts(
            sample_ids, set_a, set_b, same_namespace=same_namespace
        )
        return

    index = {sample_id: row for row, sample_id in enumerate(sample_ids)}
    names_a = namespace_names(set_a)
    names_b = names_a if same_namespace else namespace_names(set_b)
    if not names_a or not names_b:
        return

    matrix_a = _presence_matrix(sample_ids, set_a, names_a, index)
    matrix_b = matrix_a if same_namespace else _presence_matrix(
        sample_ids, set_b, names_b, index
    )
    carriers_a = matrix_a.sum(axis=0)
    carriers_b = carriers_a if same_namespace else matrix_b.sum(axis=0)
    # The j-blocks that lie entirely at or below the diagonal contribute
    # nothing for a within-namespace test; skipping them halves the work.
    first_j = 0 if not same_namespace else None

    blocks_a = range(0, len(names_a), chunk_features)
    for start_a in blocks_a:
        stop_a = min(start_a + chunk_features, len(names_a))
        block_a = matrix_a[:, start_a:stop_a]
        lower = start_a if first_j is None else first_j
        for start_b in range(lower, len(names_b), chunk_features):
            stop_b = min(start_b + chunk_features, len(names_b))
            product = block_a.T @ matrix_b[:, start_b:stop_b]
            if same_namespace:
                rows = np.arange(start_a, stop_a)[:, None]
                columns = np.arange(start_b, stop_b)[None, :]
                keep = columns > rows
            else:
                keep = np.ones(product.shape, dtype=bool)
            if not keep.any():
                continue
            local_rows, local_columns = np.nonzero(keep)
            counts = product[local_rows, local_columns]
            for local_row, local_column, count in zip(
                local_rows.tolist(), local_columns.tolist(), counts.tolist()
            ):
                yield (
                    names_a[start_a + local_row],
                    names_b[start_b + local_column],
                    int(carriers_a[start_a + local_row]),
                    int(carriers_b[start_b + local_column]),
                    int(round(count)),
                )


def naive_pair_counts(
    sample_ids: Sequence[str],
    set_a: FeatureSets,
    set_b: FeatureSets,
    *,
    same_namespace: bool,
    chunk_features: int = DEFAULT_COUNT_CHUNK_FEATURES,
) -> Iterator[PairCount]:
    """The same counting, written as the Python pair loop it replaced.

    Kept deliberately naive and deliberately reachable: it is the oracle the
    matrix product is tested against, and ``cooccurrence.pair_counter: naive``
    selects it so the comparison is possible on real data and not only in a
    test. Its cost is the number this change was made to remove - O(F^2) pairs,
    each rescanning the cohort for carriers - so it is not a fallback anybody
    should run on a cohort-sized feature set.
    """
    del chunk_features  # Accepted so both counters share one signature.
    names_a = namespace_names(set_a)
    names_b = names_a if same_namespace else namespace_names(set_b)
    if same_namespace:
        pairs = all_pairs(names_a)
    else:
        pairs = list(itertools.product(names_a, names_b))
    for name_a, name_b in pairs:
        n_a = sum(1 for s in sample_ids if name_a in set_a.get(s, ()))
        n_b = sum(1 for s in sample_ids if name_b in set_b.get(s, ()))
        n_both = sum(
            1
            for s in sample_ids
            if name_a in set_a.get(s, ()) and name_b in set_b.get(s, ())
        )
        yield name_a, name_b, n_a, n_b, n_both


def resolve_pair_counter(config: PipelineConfig) -> str:
    """Which counting strategy ``cooccurrence.pair_counter`` asks for.

    An unrecognised value is a configuration error rather than a silent
    fallback: defaulting a reader who asked for the oracle to production
    counting, or the reverse, is exactly the kind of quiet substitution that
    makes a number untraceable.
    """
    raw = (config.raw.get("cooccurrence", {}) or {}).get("pair_counter")
    if raw is None:
        return "matrix_product"
    name = str(raw)
    if name not in PAIR_COUNTERS:
        raise ValueError(
            f"cooccurrence.pair_counter is {name!r}; it must be one of "
            f"{', '.join(PAIR_COUNTERS)}."
        )
    return name


def resolve_count_chunk_features(config: PipelineConfig) -> int:
    """``cooccurrence.count_chunk_features``, validated as a positive int."""
    raw = (config.raw.get("cooccurrence", {}) or {}).get("count_chunk_features")
    if raw is None:
        return DEFAULT_COUNT_CHUNK_FEATURES
    try:
        value = int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"cooccurrence.count_chunk_features is {raw!r}; it is a block edge "
            "and must be an integer."
        ) from exc
    if value < 1:
        raise ValueError(
            f"cooccurrence.count_chunk_features is {value}; it is a block edge, "
            "so it must be at least 1."
        )
    return value


def pair_counter(
    config: PipelineConfig,
    sample_ids: Sequence[str],
    set_a: FeatureSets,
    set_b: FeatureSets,
    *,
    same_namespace: bool,
) -> Iterator[PairCount]:
    """The configured counter, with one signature whichever it is."""
    chunk = resolve_count_chunk_features(config)
    if resolve_pair_counter(config) == "naive":
        return naive_pair_counts(
            sample_ids, set_a, set_b, same_namespace=same_namespace,
            chunk_features=chunk,
        )
    return matrix_product_pair_counts(
        sample_ids, set_a, set_b, same_namespace=same_namespace,
        chunk_features=chunk,
    )


def resolve_feature_classes(config: PipelineConfig) -> Tuple[str, ...]:
    """``cooccurrence.feature_classes``: which classes are eligible to be paired.

    Order is the declared order, and it is preserved so a refusal that names the
    admitted universe reads in the order the configuration asked for rather than
    in whatever order a set happened to iterate.

    An empty list is refused rather than run: with no eligible class there is
    nothing to pair, and a stage that returns an empty table on purpose is
    indistinguishable from a stage that found nothing.
    """
    raw = (config.raw.get("cooccurrence", {}) or {}).get("feature_classes")
    if raw is None:
        return ELIGIBLE_FEATURE_CLASSES
    if isinstance(raw, str) or not isinstance(raw, (list, tuple)):
        raise ValueError(
            f"cooccurrence.feature_classes is {raw!r}; it is a list of class "
            f"names, one of: {', '.join(ELIGIBLE_FEATURE_CLASSES)}."
        )
    names: List[str] = []
    for entry in raw:
        name = str(entry)
        if name not in ELIGIBLE_FEATURE_CLASSES:
            raise ValueError(
                f"cooccurrence.feature_classes contains {name!r}; this stage "
                f"can only pair {', '.join(ELIGIBLE_FEATURE_CLASSES)}. Gene "
                "presence/absence is not one of them - stage 7's "
                "`gene_presence_absence.tsv` is not read here."
            )
        if name not in names:
            names.append(name)
    if not names:
        raise ValueError(
            "cooccurrence.feature_classes is empty; with no eligible class "
            "there is nothing to pair, and an empty table would be "
            "indistinguishable from a stage that found nothing."
        )
    return tuple(names)


def resolve_min_prevalence(config: PipelineConfig) -> float:
    """``cooccurrence.min_prevalence``: the smallest share of the cohort a
    feature must be carried by to be eligible.

    A fraction, not a count, because the thing it has to be is scale-relative: a
    filter expressed in samples would mean something different at n=20 and n=900,
    which is exactly the confusion this key exists to remove.

    ``1.0`` is refused rather than honoured: it admits only features present in
    every sample, and :func:`cooccurrence_from_counts` already drops those as
    having no contrast to test, so it would guarantee an empty table.
    """
    raw = (config.raw.get("cooccurrence", {}) or {}).get("min_prevalence")
    if raw is None:
        return DEFAULT_MIN_PREVALENCE
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"cooccurrence.min_prevalence is {raw!r}; it is a fraction of the "
            "cohort and must be a number."
        ) from exc
    if value != value or value < 0.0:
        raise ValueError(
            f"cooccurrence.min_prevalence is {value}; it is a fraction of the "
            "cohort, so it must not be negative."
        )
    if value >= 1.0:
        raise ValueError(
            f"cooccurrence.min_prevalence is {value}; a feature carried by the "
            "whole cohort has no contrast to test and "
            "`cooccurrence_from_counts` already drops it, so this value would "
            "guarantee an empty table. Use a value below 1.0."
        )
    return value


def resolve_max_features(config: PipelineConfig) -> Optional[int]:
    """``cooccurrence.max_features``: the ceiling on the admitted universe.

    ``None`` (an absent key) means no ceiling, which is the honest default for a
    key nobody has chosen a value for yet - and the reason the key is required
    in `config/science.yaml`. Any present value must be a positive integer; a
    cap of 0 or a negative cap would refuse every run including an empty one,
    which is a config error wearing a refusal's clothes.
    """
    raw = (config.raw.get("cooccurrence", {}) or {}).get("max_features")
    if raw is None:
        return None
    try:
        value = int(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"cooccurrence.max_features is {raw!r}; it is a count of features "
            "and must be an integer."
        ) from exc
    if value < 1:
        raise ValueError(
            f"cooccurrence.max_features is {value}; it is a count of features, "
            "so it must be at least 1."
        )
    return value


def _prevalence_floor(min_prevalence: float, n_samples: int) -> int:
    """How many cohort samples a feature must be carried by to be eligible.

    ``ceil``, and never below 1: a filter that admitted nothing would empty the
    table, and one that admitted a feature carried by no cohort sample would let
    through a name that only a non-manifest sample holds - which the degenerate
    rule in :func:`cooccurrence_from_counts` then drops anyway, one wasted pair
    per such feature rather than a clean universe.

    At ``min_prevalence = 0.05`` this is 1 for every cohort up to n=20, which is
    why the committed fixtures do not move when the key is set.

    Rounded to 9 decimals before the ceiling, because ``ceil`` on a float
    product is a trap: ``0.07 * 100`` is ``7.000000000000001`` in IEEE 754, so a
    plain ``ceil`` gives 8 and a feature carried by exactly 7 of 100 samples is
    excluded by a filter that says 7%. The rounding is well inside the precision
    a config-authored fraction can carry, and it puts the boundary back where a
    reader computes it by hand.
    """
    return max(1, math.ceil(round(min_prevalence * n_samples, 9)))


def admit_features(
    namespaces: Mapping[str, FeatureSets],
    classes: Sequence[str],
    min_prevalence: float,
    n_samples: int,
) -> Dict[str, FeatureSets]:
    """The configured classes, minus the features too rare to be testable.

    Applied by ``_compute`` to both TEST and REAL, so the universe a pair is
    drawn from is a property of the computation and not of where the features
    were loaded from. It runs before :func:`pair_counter`, which is the point:
    a feature removed here never reaches the matrix product, and so never costs
    F-1 pairs.

    Prevalence is counted over **every** key of a namespace mapping, not only
    the manifest's samples, so it matches :func:`namespace_names` - the counter's
    view of the universe. Filtering on one view and counting on another would
    leave a rare feature in the universe that the counter then pairs anyway.

    Returns a fresh mapping; the caller's namespaces are not mutated, because
    ``_run_real`` builds them once and a test may want to inspect what was
    refused.
    """
    wanted = set(classes)
    floor = _prevalence_floor(min_prevalence, n_samples)
    admitted: Dict[str, FeatureSets] = {}
    for name, sets in namespaces.items():
        if name not in wanted:
            continue
        carriers: Dict[str, int] = {}
        for members in sets.values():
            for feature in members:
                carriers[feature] = carriers.get(feature, 0) + 1
        keep = {feature for feature, count in carriers.items() if count >= floor}
        admitted[name] = {
            sample_id: {f for f in members if f in keep}
            for sample_id, members in sets.items()
        }
    return admitted


def enforce_feature_cap(
    admitted: Mapping[str, FeatureSets],
    cap: Optional[int],
) -> Dict[str, int]:
    """Refuse a universe larger than ``cooccurrence.max_features`` allows.

    The count is the **total across the configured classes**, not a per-class
    one, because the work is not per-class: a total of F features admits at most
    C(F, 2) pairs however they are distributed (the cross products
    ``f_i * f_j`` and the within-class ``C(f_i, 2)`` together sum to at most
    ``C(F, 2)``), so a total is what bounds the run. A per-class cap would let
    three classes of 4,000 through to a worse position than one class of 12,000.

    Names the count, the cap and the per-class breakdown, because a reader
    holding the refusal needs to know which class to thin and not only that
    something is too big.
    """
    from ..errors import StageError

    sizes = {
        name: len(namespace_names(sets))
        for name, sets in sorted(admitted.items())
        if sets
    }
    total = sum(sizes.values())
    if cap is None or total <= cap:
        return sizes
    breakdown = ", ".join(f"{name}={size}" for name, size in sizes.items())
    raise StageError(
        f"cooccurrence refuses to run: the eligible feature universe is "
        f"{total} features after cooccurrence.min_prevalence and "
        f"cooccurrence.feature_classes ({breakdown}), which is over "
        f"cooccurrence.max_features={cap}. At {total} features this stage "
        f"tests C({total},2) = {total * (total - 1) // 2:,} pairs, and every "
        "one of them costs a scipy.stats.fisher_exact call and a Cooccurrence "
        "row - measured at 128.0 us and 244 resident bytes per pair at n=900, "
        "against 0.34 us to count it. To run: raise cooccurrence.max_features "
        "if the machine can carry it, raise cooccurrence.min_prevalence, "
        "narrow cooccurrence.feature_classes, or split the analysis. Refusing "
        "is the point; the alternative is an out-of-memory failure that reads "
        "as a bioinformatics bug rather than a configuration decision.",
        stage="cooccurrence",
        features=total,
        max_features=cap,
    )


def test_pair(
    feature_a: str,
    feature_b: str,
    set_a: FeatureSets,
    set_b: FeatureSets,
    feature_type: str,
    sample_ids: Sequence[str],
    lineages: Optional[Mapping[str, str]] = None,
) -> Optional[Cooccurrence]:
    """Test one feature pair.

    ``set_a``/``set_b`` map sample ID -> set of features carried. A sample
    counts toward a feature only if the feature is in that sample's set.

    This counts the pair itself rather than going through
    :func:`pair_counter`, so it stays usable on one pair without building a
    matrix. The arithmetic is shared with the blocked path through
    :func:`cooccurrence_from_counts`, which is what keeps the two in step.

    ``lineages`` is accepted and **deliberately not read**. The test is pooled
    across the cohort, so a per-sample lineage label has nothing to act on
    here, and giving it one would be lineage stratification - a feature, not a
    detail, and not this module's to add silently. The crosswalk is enforced
    earlier instead: ``_add_missing_lineages`` makes REAL reject a missing or
    degenerate one, so a pooled row is never published next to a crosswalk the
    reader could assume stratified it. See the module docstring.

    Returns ``None`` when the pair has no contrast to test, which keeps
    sparse tables from being dominated by noise.
    """
    a = set_a
    b = set_b
    n = len(sample_ids)
    n_a = sum(1 for s in sample_ids if feature_a in a.get(s, ()))
    n_b = sum(1 for s in sample_ids if feature_b in b.get(s, ()))
    n_both = sum(
        1 for s in sample_ids
        if feature_a in a.get(s, ()) and feature_b in b.get(s, ())
    )
    return cooccurrence_from_counts(
        feature_a, feature_b, feature_type, n_a, n_b, n_both, n
    )


def all_pairs(names: Iterable[str]) -> List[Tuple[str, str]]:
    """Unordered pairs of distinct names, in sorted order."""
    ordered = sorted(set(names))
    return [
        (ordered[i], ordered[j])
        for i in range(len(ordered))
        for j in range(i + 1, len(ordered))
    ]


def benjamini_hochberg(pvalues: Mapping[str, float]) -> Dict[str, float]:
    """BH FDR adjustment, keyed by an arbitrary string."""
    items = sorted(((k, v) for k, v in pvalues.items() if v is not None), key=lambda kv: kv[1])
    n = len(items)
    out: Dict[str, float] = {}
    previous = 1.0
    for rank in range(n, 0, -1):
        key, p_value = items[rank - 1]
        value = min(previous, p_value * n / rank)
        out[key] = min(1.0, value)
        previous = out[key]
    return out


#: What holds each namespace back, for the refusal message.
_COOCCURRENCE_BLOCKERS = {
    "gene": "stage 4 (amr); run with --only amr first",
    "mutation": (
        "stage 6 (variants) has no REAL producer yet - ticket 14/15; only "
        "papipeline.testing.synthetic writes regulator_variants.tsv"
    ),
    "mechanism": (
        "stage 5 (mechanisms), which derives from the same regulator_variants "
        "table and so is empty for the same reason"
    ),
    "lineages": (
        "stage 9 (phylogeny), which builds `lineage_label` from stage 3's MLST "
        "sequence type (`lineage.method: st`). A REAL run needs BOTH of those "
        "outputs, and this stage reads stage 9's `tree_metadata.tsv`: stage 3 "
        "must have written `mlst_results.tsv` under this run's tool-output "
        "root, and stage 9 must then have written a `tree_metadata.tsv` whose "
        "`lineage_label` column is populated. An all-`unknown` crosswalk - "
        "samples present in the manifest but absent from stage 3's table - is "
        "refused rather than reported, because `COOCCURRENCE_COLUMNS` carries no "
        "lineage field and a pooled row published beside a sentinel crosswalk "
        "reads as if it had been stratified. NOTE: no test in this module reads "
        "the crosswalk (see `test_pair`), so this gate costs a refusal and "
        "changes no output value; that trade is recorded in the round-12 "
        "backlog rather than changed here."
    ),
}


def real_input_paths(
    config: PipelineConfig,
    intermediate_root: Path,
    mode: RunMode,
    phylogeny_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """Where each REAL input is read from, for the refusal message.

    The stage-table locations come from ``contracts`` so this cannot name a
    stale path. ``phylogeny_dir`` is overridable for the reason given in
    ``convergence.real_input_paths``.
    """
    from ..execution.contracts import internal_table_path, table_path

    stage_dir = Path(intermediate_root) / "stages"
    phylo = Path(phylogeny_dir) if phylogeny_dir else config.phylogeny_dir(mode)
    return {
        "gene": str(table_path(stage_dir, "amr")),
        "mutation": str(
            Path(intermediate_root) / "regulators" / "regulator_variants.tsv"
        ),
        "mechanism": str(internal_table_path(stage_dir, "mechanisms")),
        "lineages": str(phylo / "tree_metadata.tsv"),
    }


def run(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    amr_genes: Optional[Mapping[str, Sequence[str]]] = None,
    variants: Optional[Mapping[str, Sequence[str]]] = None,
    mechanisms: Optional[Mapping[str, Sequence[str]]] = None,
    lineages: Optional[Mapping[str, str]] = None,
    *,
    intermediate_root: Optional[Path] = None,
    antibiotic: Optional[str] = None,
    phylogeny_dir: Optional[Path] = None,
) -> List[Cooccurrence]:
    """Stage 14 entry point.

    Args:
        amr_genes: sample -> acquired gene names. In REAL these are **re-read
            from disk**; a standalone or Snakemake invocation has no in-memory
            stage-4 result to inherit.
        variants: sample -> variant names.
        mechanisms: sample -> mechanism names.
        intermediate_root: REAL only. Where the inputs above are read from.
        antibiotic: REAL only. Selects the AMR rows for the drug under analysis.

    Returns:
        All testable pairs, with BH-adjusted p-values. Significant pairs are
        *not* filtered out of the return value; the report selects them, so
        that a reader can apply their own threshold.

    Raises:
        StageError: In REAL, when any declared namespace carries no data.
            Nothing is written, so no partial table is left behind.
        NotImplementedError: In STUB, which fabricates this table instead.
    """

    if mode is not RunMode.TEST:
        if mode is RunMode.REAL:
            return _run_real(
                config,
                manifest,
                mode,
                amr_genes,
                variants,
                mechanisms,
                lineages,
                intermediate_root=intermediate_root,
                antibiotic=antibiotic,
                phylogeny_dir=phylogeny_dir,
            )
        # STUB fabricates this table in `papipeline.stub`, so reaching the real
        # computation here would mean STUB stopped stubbing.
        raise NotImplementedError(
            "STUB-mode cooccurrence is fabricated by papipeline.stub and has "
            f"no caller here. (mode={getattr(mode, 'value', mode)})"
        )

    return _compute(
        config,
        manifest,
        build_feature_sets(manifest, amr_genes or {}, variants or {}, mechanisms or {}),
        lineages or {},
        (config.raw.get("cooccurrence", {}) or {}).get("tests") or [],
    )


def _compute(
    config: PipelineConfig,
    manifest: SampleManifest,
    namespaces: Mapping[str, FeatureSets],
    lineages: Mapping[str, str],
    configured: Sequence[str],
) -> List[Cooccurrence]:
    """Test every configured pair set over already-built namespaces.

    Split out of `run` so the TEST and REAL paths run *identical* logic and
    cannot drift: the only difference between them is where the namespaces
    came from and whether they were checked for completeness first.

    ``lineages`` is threaded through to ``test_pair`, which does not read it -
    see that function's docstring. It is accepted here so the two paths cannot
    diverge in *which* crosswalk they were handed.

    The **configured universe** is applied here, once, for both modes: which
    classes are eligible (``cooccurrence.feature_classes``), how common a
    feature has to be (``cooccurrence.min_prevalence``) and how many may be
    admitted (``cooccurrence.max_features``). It belongs at the top of the one
    function both modes run rather than in either caller, because a universe that
    differed by mode would make a TEST result and a REAL result different numbers
    for the same isolates.
    """
    sample_ids = manifest.sample_ids
    tests = set(configured)

    classes = resolve_feature_classes(config)
    min_prevalence = resolve_min_prevalence(config)
    namespaces = admit_features(
        namespaces, classes, min_prevalence, len(sample_ids)
    )
    cap = resolve_max_features(config)
    sizes = enforce_feature_cap(namespaces, cap)
    LOGGER.info(
        "Stage 14 universe: %d eligible features (%s), min_prevalence=%s, "
        "cap=%s",
        sum(sizes.values()),
        ", ".join(f"{name}={size}" for name, size in sizes.items()) or "none",
        min_prevalence,
        cap if cap is not None else "unset",
    )

    results: List[Cooccurrence] = []

    def add(namespace_a: str, namespace_b: str, label: str) -> None:
        """Test every pair drawn from the two namespaces.

        For a within-namespace test the pairs are distinct unordered
        combinations; for a cross-namespace test they are the full cross
        product, with the pair written in a stable order.

        The pairs are enumerated by the configured counter rather than
        materialised, so a namespace with 50,000 features does not first build a
        1.25e9-element list. Every triple goes through
        :func:`cooccurrence_from_counts`, which is also what
        :func:`test_pair` calls, so there is one implementation of the
        statistics and one of the degenerate-pair rule.

        This family's rows are sorted before they join the result, so the order
        of ``14_cooccurrence.tsv`` is a property of the *features* and not of
        the counter or the chunk edge. Without that, switching
        ``cooccurrence.pair_counter`` would reorder the table for no reason a
        reader could see.
        """
        set_a = namespaces.get(namespace_a) or {}
        set_b = namespaces.get(namespace_b) or {}
        if not set_a or not set_b:
            return

        same = namespace_a == namespace_b
        if same and len(namespace_names(namespaces.get(namespace_a) or {})) < 2:
            return

        counted = 0
        degenerate = 0
        family: List[Cooccurrence] = []
        for name_a, name_b, n_a, n_b, n_both in pair_counter(
            config, sample_ids, set_a, set_b, same_namespace=same
        ):
            counted += 1
            result = cooccurrence_from_counts(
                name_a, name_b, label, n_a, n_b, n_both, len(sample_ids)
            )
            if result is None:
                degenerate += 1
                continue
            family.append(result)
        family.sort(key=lambda r: (r.feature_a, r.feature_b))
        results.extend(family)
        LOGGER.debug(
            "Stage 14 %s: %d pairs counted, %d degenerate, counter=%s",
            label,
            counted,
            degenerate,
            resolve_pair_counter(config),
        )

    if "genes" in tests:
        add("gene", "gene", "gene_gene")
    if "gene_mutation" in tests:
        add("gene", "mutation", "gene_mutation")
    if "mutation_mutation" in tests:
        add("mutation", "mutation", "mutation_mutation")
    if "gene_mechanism" in tests:
        add("gene", "mechanism", "gene_mechanism")
    if "mechanism_mechanism" in tests:
        add("mechanism", "mechanism", "mechanism_mechanism")

    pvalues = {
        f"{r.feature_a}|{r.feature_b}": r.adjusted_p_value
        for r in results
        if r.adjusted_p_value is not None
    }
    adjusted = benjamini_hochberg(pvalues)
    # Rebuilt IN PLACE, and this is not a style preference. At n=900 the emitted
    # rows are ~1.4 kB each (measured), so holding the pre-BH and post-BH lists
    # simultaneously doubles the largest allocation this stage makes - and the
    # alternative to halving it is a cap set by a number that is half the cost of
    # the operation. `Cooccurrence` is frozen, so each row is a new object either
    # way; overwriting the slot frees the old one immediately instead of leaving
    # it for the garbage collector to find at the end of the comprehension.
    for index, r in enumerate(results):
        results[index] = Cooccurrence(
            feature_a=r.feature_a,
            feature_b=r.feature_b,
            feature_type=r.feature_type,
            n_a=r.n_a,
            n_b=r.n_b,
            n_both=r.n_both,
            statistic=r.statistic,
            statistic_value=r.statistic_value,
            adjusted_p_value=adjusted.get(f"{r.feature_a}|{r.feature_b}"),
        )
    finalised = results

    threshold = float(
        config.raw.get("cooccurrence", {}).get("significance_threshold", 0.05)
    )
    hits = [
        r for r in finalised if r.adjusted_p_value is not None and r.adjusted_p_value <= threshold
    ]
    LOGGER.info(
        "Stage 14: %d pairs tested across %d namespaces | %d pairs at adj.p<=%s "
        "(association only, not causal)",
        len(finalised),
        len(namespaces),
        len(hits),
        threshold,
    )
    return finalised

def _load_mechanism_names(path: Path) -> Dict[str, List[str]]:
    """Sample -> mechanism names, read from stage 5's internal table.

    Written to ``intermediate/stages/05_mechanisms.tsv`` by the step spec.md:351
    folds into this stage, so a standalone invocation reads it rather than
    re-deriving it.
    """
    from ..io.tsv import read_tsv

    table = Path(path)
    if not table.is_file():
        # See `convergence.load_amr_names`: a missing input is a refusal that
        # names it, not a `read_tsv` error from three frames down.
        return {}
    names: Dict[str, List[str]] = {}
    for row in read_tsv(table, required_columns=("sample_id",)):
        mechanism = row.get("mechanism")
        if not mechanism:
            continue
        names.setdefault(str(row["sample_id"]), []).append(str(mechanism))
    return {sample: sorted(set(values)) for sample, values in names.items()}


def _run_real(
    config: PipelineConfig,
    manifest: SampleManifest,
    mode: RunMode,
    amr_genes: Optional[Mapping[str, Sequence[str]]],
    variants: Optional[Mapping[str, Sequence[str]]],
    mechanisms: Optional[Mapping[str, Sequence[str]]],
    lineages: Optional[Mapping[str, str]],
    *,
    intermediate_root: Optional[Path],
    antibiotic: Optional[str],
    phylogeny_dir: Optional[Path] = None,
) -> List[Cooccurrence]:
    """The REAL path: read from disk, then require every namespace to be real.

    The completeness check runs against the **built namespaces**, not the raw
    mappings. Those differ in a way that matters: a sample present in the
    manifest with no determinant maps to an empty list, which
    `build_feature_sets` turns into an empty set. Checking the raw mapping would
    see a non-empty dict of empty lists and conclude the namespace has data.
    """
    from .real_inputs import is_effectively_empty, refuse_incomplete

    if intermediate_root is None:
        raise refuse_incomplete(
            "cooccurrence",
            {
                "intermediate_root": (
                    "REAL mode reads its inputs from disk, and no intermediate "
                    "root was supplied, so there is nowhere to read them from"
                )
            },
            blocked_by={
                "intermediate_root": (
                    "a caller that passes intermediate_root - run.py, or the "
                    "Snakemake rule for this stage"
                )
            },
        )

    root = Path(intermediate_root)
    paths = real_input_paths(config, root, mode, phylogeny_dir)
    drug = antibiotic or str(
        (getattr(config, "raw", None) or {}).get("project", {}).get(
            "primary_antibiotic"
        )
        or config.antibiotics[0]
    )

    from .convergence import load_amr_names

    amr_genes = load_amr_names(paths["gene"], drug)
    variants = _load_variant_names(Path(paths["mutation"]))
    mechanisms = _load_mechanism_names(Path(paths["mechanism"]))
    lineages = _load_lineages(Path(paths["lineages"]))

    namespaces = build_feature_sets(manifest, amr_genes, variants, mechanisms)

    configured = [
        t for t in (config.raw.get("cooccurrence", {}) or {}).get("tests") or []
    ]
    needed = _namespaces_for_tests(configured, resolve_feature_classes(config))

    missing = {}
    for name in needed:
        if is_effectively_empty(
            {s: sorted(v) for s, v in (namespaces.get(name) or {}).items() if v}
        ):
            missing[name] = "no sample carries a feature in this namespace"
    missing = _add_missing_lineages(missing, lineages)

    if missing:
        raise refuse_incomplete(
            "cooccurrence",
            missing,
            paths=paths,
            blocked_by=_COOCCURRENCE_BLOCKERS,
        )

    return _compute(config, manifest, namespaces, lineages, configured)


def _namespaces_for_tests(
    configured: Sequence[str],
    classes: Optional[Sequence[str]] = None,
) -> List[str]:
    """Which namespaces a set of configured tests actually reads.

    A test that is not configured needs no data, so requiring its namespace
    would refuse a run that has everything the configuration asks for.

    ``classes`` narrows the answer to the configured feature universe. Without
    it, ``feature_classes: [gene]`` would still demand a non-empty mutation
    namespace here - so a run that deliberately excludes mutations would be
    refused by this gate for data it has declared it will not use, and the
    narrowing that happens in `_compute` would never be reached.
    """
    eligible = set(classes) if classes is not None else None
    needed: List[str] = []
    for test_name in configured:
        if test_name == "genes":
            names = ["gene"]
        elif test_name == "gene_mutation":
            names = ["gene", "mutation"]
        elif test_name == "mutation_mutation":
            names = ["mutation"]
        elif test_name == "gene_mechanism":
            names = ["gene", "mechanism"]
        elif test_name == "mechanism_mechanism":
            names = ["mechanism"]
        else:
            continue
        for name in names:
            if eligible is not None and name not in eligible:
                continue
            if name not in needed:
                needed.append(name)
    return needed


def _add_missing_lineages(
    missing: Dict[str, str], lineages: Mapping[str, str]
) -> Dict[str, str]:
    """Lineages are required, though no test reads them.

    ``test_pair`` pools across the cohort and takes no lineage input, so a
    missing or all-``unknown`` crosswalk would not raise on its own - the stage
    would publish pooled counts carrying no lineage stratification whatsoever.
    Nothing in the output says so, since ``COOCCURRENCE_COLUMNS`` has no lineage
    field, so the reader has no way to notice. Refusing is the alternative to
    letting them assume a stratification that was never performed.
    """
    from .real_inputs import is_effectively_empty

    if is_effectively_empty(lineages, sentinel="unknown"):
        missing["lineages"] = (
            "the sample -> lineage crosswalk carries no real lineage label"
        )
    return missing


def _load_variant_names(path: Path) -> Dict[str, List[str]]:
    """Reused from `convergence` so both stages spell a variant the same way."""
    from .convergence import load_variant_names as _load

    return _load(Path(path))


def _load_lineages(path: Path) -> Dict[str, str]:
    """Stage 9's crosswalk. A missing file yields ``{}``, which the completeness
    check then refuses - see `convergence._load_lineages` for why this does not
    warn-and-continue."""
    if not Path(path).exists():
        return {}
    from .phylogeny import load_tree_metadata

    return load_tree_metadata(Path(path))
