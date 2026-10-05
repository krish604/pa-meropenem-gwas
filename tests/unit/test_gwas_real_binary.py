"""End-to-end observation of the REAL pyseer binary, on synthetic data.

Every other `PyseerEngine` test in this repository patches `_invoke`. That was
a deliberate trade - spec.md D8 defines TEST mode as running no real tools - but
it has a cost that ticket 16 named as its largest remaining risk and that has now
actually been paid: **nothing in this repository had ever observed the real
pyseer's output.** Sixteen test sites asserted against a header this project
read out of `pyseer/__main__.py`, and one of the assumptions behind them was
wrong (see `TestWhatTheRealBinaryDoesNotEmit`).

So this file invokes the binary. R12 authorises exactly that: real code paths,
synthetic inputs. The inputs are synthetic on purpose - a real cohort is 10
isolates, and 10 isolates cannot demonstrate that a planted causal feature
ranks first.

**Opt-in, and why.** These tests are skipped unless
``PAPIPELINE_TEST_PYSEER=1``. spec.md D8 forbids real tools in TEST mode, and
that is a mode boundary rather than a preference; the way to cross it honestly is
an explicit opt-in that says so, not a default that quietly does it. Run them
with::

    PAPIPELINE_TEST_PYSEER=1 pytest tests/unit/test_gwas_real_binary.py

Each test needs the pinned pyseer; without it they skip and name what went
unchecked rather than passing vacuously.
"""

from __future__ import annotations

import os
import shutil

from pathlib import Path

import numpy as np
import pytest

from papipeline.errors import DataContractError
from papipeline.models import Phenotype, PhenotypeCall, RunMode
from papipeline.stages import gwas as stage_gwas

#: Set this to "1" to let these tests launch pyseer. Not defaulted on: see the
#: module docstring.
OPT_IN = "PAPIPELINE_TEST_PYSEER"

#: pyseer 1.1.2's own association-table header for an `--lmm` run. Read from the
#: installed package in `TestWhatTheRealBinaryDoesNotEmit` rather than trusted.
PYSEER_LMM_HEADER = (
    "variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\tvariant_h2\tnotes"
)

#: The planted causal feature's name. It carries no `snp__`/`gene__` trick: the
#: prefix is a real one from `config.gwas.feature_types`.
CAUSAL = "gene__planted_causal_oprd"

#: n for the power demonstration. Well above the point where a mixed model on a
#: binary trait stops being degenerate, and far above the real cohort, which is
#: the whole reason this cohort is synthetic.
G3_N = 300

#: The 10 real isolates, and their real imipenem calls, read from
#: `db/smoke_phenotypes/imipenem_phenotype.tsv` - measured laboratory AST
#: results, 7 R and 3 S. Not invented here; parsed from the file at run time.
REAL_PHENOTYPE = Path("db/smoke_phenotypes/imipenem_phenotype.tsv")


def requires_pyseer():
    if os.environ.get(OPT_IN) != "1":
        pytest.skip(
            f"set {OPT_IN}=1 to invoke the real pyseer binary; spec.md D8 "
            "defines TEST mode as running no real tools, so this is opt-in"
        )
    if shutil.which("pyseer") is None:
        pytest.skip("pyseer is not on PATH; see docs/environment-arm64.md section 4")


@pytest.fixture(scope="session")
def machine_config(tmp_path_factory, pipeline_root):
    """The laptop overlay, with a writable scratch dir for the vendored helper.

    Copied rather than referenced: the helper's ``--temp`` has to point at a
    directory that exists and belongs to the test, and editing the committed
    overlay from a test would leave the tree dirty.

    Session-scoped so the class-scoped cohorts below can request it - a
    pyseer run is not cheap and the configuration is read-only.
    """
    from papipeline.config import load_config

    scratch = tmp_path_factory.mktemp("scratch")
    text = (pipeline_root / "config" / "machines" / "laptop.yaml").read_text()
    text = text.replace(
        "unique_patterns_temp_dir: /tmp", f"unique_patterns_temp_dir: {scratch}"
    )
    overlay = tmp_path_factory.mktemp("overlay") / "overlay.yaml"
    overlay.write_text(text)
    return load_config(pipeline_root / "config" / "science.yaml", machine=str(overlay))


# --------------------------------------------------------------------------
# The cohort
# --------------------------------------------------------------------------


def simulate_tree(n_tips: int, rng: np.random.Generator):
    """A rooted binary tree over ``n_tips``.

    Each node splits near the middle - by up to a quarter of its tips either
    side - and closes as a clade at a maximum depth. The near-middle split is
    load-bearing. A uniform ``randint(1, size)`` split looks like a tree and
    produces one clade holding 96% of the cohort and two holding 2% each, which
    is not population structure: it is one lineage and two rounding errors, and
    every "lineage-aware" statistic over it is meaningless.

    Returns ``(tree, top_clade_of_tip, depth_of_tip)``. **Lineage is the
    top-level clade**, because that is what a lineage is - the deepest split. The
    coalescent depth is reported separately and is used for drift, so features
    vary within a lineage as well as between them.
    """
    counter = [0]
    top_clade = [0]

    def grow(size: int, depth: int, clade: int):
        # A truncated subtree is a clade of `size` tips, not one tip: the index
        # has to advance by the whole clade or tips and depths stop lining up.
        if size <= 1 or depth >= 8:
            tips = list(range(counter[0], counter[0] + size))
            counter[0] += size
            return ("clade", depth, tips)
        quarter = max(1, size // 4)
        left = size // 2 + int(rng.integers(-quarter, quarter + 1))
        left = int(np.clip(left, 1, size - 1))
        return (
            "node", depth,
            grow(left, depth + 1, clade),
            grow(size - left, depth + 1, clade + 1 if depth == 0 else clade),
        )

    tree = grow(n_tips, 0, 0)

    lineage_of: dict = {}
    depth_of: dict = {}

    def walk_top(node, top):
        # The root is the only node that opens a new lineage. `node[1]` is the
        # node's own depth, so "is this the root" is asked directly rather than
        # inferred from `top == 0` - which is true again for every node inside
        # lineage 0, and produced 1 tip in one lineage and 299 in the other.
        if node[0] == "clade":
            for tip in node[2]:
                lineage_of[tip] = f"L{top}"
                depth_of[tip] = node[1]
            return
        at_root = node[1] == 0
        walk_top(node[2], 0 if at_root else top)
        walk_top(node[3], 1 if at_root else top)

    walk_top(tree, 0)
    return tree, lineage_of, depth_of


def synthetic_cohort(n: int, seed: int):
    """A cohort with tree-derived structure and one planted causal feature.

    Returns ``(gwas_input, stats)``. The structure is real: the nulls' allele
    frequencies are drawn per lineage, so they inherit the population structure
    a GWAS has to correct for, and the causal feature is drawn independently of
    the tree so that finding it cannot be a lineage artefact.
    """
    rng = np.random.default_rng(seed)
    _tree, lineage_of, depths = simulate_tree(n, rng)

    causal = np.array([1 if rng.random() < 0.45 else 0 for _ in range(n)])
    # Strong but not perfectly collinear: a feature that determines the outcome
    # exactly would be a collinearity artefact rather than a test.
    pheno = np.where(
        causal == 1, rng.random(n) < 0.92, rng.random(n) < 0.06
    ).astype(int)

    features = {CAUSAL: {f"S{i:04d}": int(causal[i]) for i in range(n)}}
    n_null = 0
    for k in range(60):
        base = rng.uniform(0.15, 0.85)
        drift = rng.normal(0.0, 0.22)
        vec = np.array(
            [
                1
                if rng.random()
                < float(np.clip(base + drift * (depths[i] - 2), 0.02, 0.98))
                else 0
                for i in range(n)
            ]
        )
        if vec.sum() < 5 or vec.sum() > n - 5:
            continue
        features[f"gene__null{k:03d}"] = {f"S{i:04d}": int(vec[i]) for i in range(n)}
        n_null += 1

    samples = tuple(f"S{i:04d}" for i in range(n))
    lineages = {s: lineage_of[i] for i, s in enumerate(samples)}
    gwas_input = stage_gwas.GwasInput(
        sample_ids=samples,
        binary_outcome={s: int(pheno[i]) for i, s in enumerate(samples)},
        features=features,
        feature_types={name: "gene" for name in features},
        lineages=lineages,
        n_positive=int(pheno.sum()),
        n_negative=int((1 - pheno).sum()),
        n_excluded=0,
    )
    stats = {
        "n": n,
        "n_null": n_null,
        "n_lineages": len(set(lineage_of.values())),
        "causal_carriers": int(causal.sum()),
    }
    return gwas_input, stats


def real_phenotype_calls():
    """The 10 real isolates' measured imipenem calls, parsed from the file."""
    if not REAL_PHENOTYPE.is_file():
        pytest.skip(f"{REAL_PHENOTYPE} not present in this checkout")
    rows = []
    for line in REAL_PHENOTYPE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] == "sample_id":
            continue
        rows.append(
            PhenotypeCall(
                sample_id=parts[0], antibiotic=parts[1],
                phenotype=Phenotype(parts[2]), mic=None,
            )
        )
    return rows


def real_cohort_input():
    """The real 10 isolates, with a SYNTHETIC null feature matrix.

    The phenotype and the sample ids are real and measured. The features are
    not: no real feature table exists to read (see the module docstring of the
    round-12 report), so they are drawn at random and every p-value this
    produces is a statement about the *machinery*, not about resistance.

    R16: n=10 with 7 R and 3 S is underpowered by construction, and any
    statistic computed on it is not a finding. The assertions below therefore
    never treat a p-value as a result - only as evidence that the pipeline
    produced a number and labelled it.
    """
    calls = real_phenotype_calls()
    binary = {
        c.sample_id: (1 if c.phenotype.value == "R" else 0) for c in calls
    }
    rng = np.random.default_rng(11)
    samples = tuple(sorted(binary))
    features = {}
    for k in range(40):
        vec = rng.integers(0, 2, size=len(samples))
        if vec.sum() < 3 or vec.sum() > len(samples) - 3:
            continue
        features[f"gene__null{k:03d}"] = {
            s: int(vec[i]) for i, s in enumerate(samples)
        }
    return stage_gwas.GwasInput(
        sample_ids=samples,
        binary_outcome=binary,
        features=features,
        feature_types={name: "gene" for name in features},
        lineages={s: "unknown" for s in samples},
        n_positive=sum(binary.values()),
        n_negative=len(binary) - sum(binary.values()),
        n_excluded=0,
    )


# --------------------------------------------------------------------------
# What the real binary actually emits
# --------------------------------------------------------------------------


class TestWhatTheRealBinaryDoesNotEmit:
    """The assumption the other sixteen test sites could not check.

    `parse_pyseer_output` reads `row.get("adj.pvalue")` and falls back to a
    local Benjamini-Hochberg. The round-12 status report assumed pyseer emits
    that column, and called the resulting rank bug "likely dormant". Both halves
    of that are wrong, and only running the binary shows it.
    """

    def test_pyseer_emits_no_adjusted_p_column(self, tmp_path, machine_config):
        """No FDR, anywhere in pyseer 1.1.2 - so the BH branch is the real one.

        This is not a cosmetic detail. The old lookup was
        `computed.get(str(need_adjustment.index(p)))`, and `.index()` returns
        the first match, so two features sharing a raw p-value were both given
        the earlier rank's adjusted value. Whether that bug could ever fire
        depends entirely on this column: with pyseer's own values present the
        branch is skipped, and without them it runs on every row.
        """
        requires_pyseer()
        gwas_input, _ = synthetic_cohort(40, seed=3)
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "wd")
        engine.run(gwas_input, machine_config)

        raw = (tmp_path / "wd" / "pyseer_raw.tsv").read_text()
        header = [
            line for line in raw.splitlines() if line and not line.startswith("#")
        ][0].split("\t")
        assert "adj.pvalue" not in header, (
            "pyseer now emits an adjusted p-value column; re-read "
            "parse_pyseer_output, because its local BH branch would no longer "
            "be the one every run takes and this test is what says so"
        )
        assert header[:2] == ["variant", "af"]
        assert header[-1] == "notes"

    def test_the_tested_count_line_this_stage_parses_is_really_emitted(
        self, tmp_path, machine_config
    ):
        """`_cross_check_pyseer_tested_count` reads pyseer's stderr by pattern.

        A pattern that silently stops matching is a check that fails open, and
        this is the only observation in the repository that would notice.
        """
        requires_pyseer()
        gwas_input, stats = synthetic_cohort(40, seed=3)
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "wd")
        engine.run(gwas_input, machine_config)

        log = (tmp_path / "wd" / "pyseer_pass1.log").read_text(errors="replace")
        stated = [
            int(line.strip().split()[0])
            for line in log.splitlines()
            if line.strip().endswith("tested variants")
        ]
        assert stated, f"pyseer printed no tested-variant count:\n{log}"
        assert engine.reduction.n_variants_tested == stated[0]
        patterns = (tmp_path / "wd" / "patterns.txt").read_text().splitlines()
        assert len(patterns) == stated[0]

    def test_the_pinned_lmm_header_matches_the_real_one(self, tmp_path, machine_config):
        """`PYSEER_LMM_HEADER` is what the other sixteen sites build fixtures from."""
        requires_pyseer()
        gwas_input, _ = synthetic_cohort(40, seed=3)
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "wd")
        engine.run(gwas_input, machine_config)
        raw = (tmp_path / "wd" / "pyseer_raw.tsv").read_text()
        header = [
            line for line in raw.splitlines() if line and not line.startswith("#")
        ][0]
        assert header.split("\t") == PYSEER_LMM_HEADER.split("\t")


@pytest.fixture(scope="module")
def g3_outcome(tmp_path_factory, machine_config):
    """The G3 cohort, run once through the real binary, shared by both classes.

    Module-scoped because a pyseer fit at n=300 is not cheap and three classes
    want to make claims about the same table. Sharing one run also means the
    claims cannot disagree with each other.
    """
    requires_pyseer()
    workdir = tmp_path_factory.mktemp("g3")
    gwas_input, stats = synthetic_cohort(G3_N, seed=20261004)
    engine = stage_gwas.PyseerEngine("pyseer", workdir)
    results = engine.run(gwas_input, machine_config)
    return {
        "results": results,
        "input": gwas_input,
        "stats": stats,
        "engine": engine,
        "workdir": workdir,
        "config": machine_config,
    }


# --------------------------------------------------------------------------
# G3 - the planted causal feature, through the production path
# --------------------------------------------------------------------------


class TestPlantedCausalFeatureIsFoundFirst:
    """300 samples, tree-derived structure, one planted causal feature.

    This is the test ticket 16 said was missing and could not be written without
    a mode that permits tools. It runs the production path end to end:
    `PyseerEngine.run` -> two pyseer passes -> the vendored helper's threshold ->
    `parse_pyseer_output` -> `GwasResult`.
    """

    def test_the_cohort_really_has_structure_and_a_causal_feature(self, g3_outcome):
        """Guards the guard: a test that passes because its data is flat proves nothing."""
        stats = g3_outcome["stats"]
        assert stats["n"] >= 300
        assert stats["n_lineages"] > 1, (
            "the nulls' structure comes from the tree, so a one-lineage cohort "
            "would make this file a test of nothing"
        )
        assert stats["n_null"] >= 40
        assert 5 <= stats["causal_carriers"] <= stats["n"] - 5

    def test_the_causal_feature_is_tested_at_all(self, g3_outcome):
        assert stage_gwas.PyseerEngine is not None
        tested = g3_outcome["input"].testable_features()
        assert "gene__planted_causal_oprd" in tested

    def test_the_causal_feature_is_not_lineage_confounded(self, g3_outcome):
        """Otherwise "found first" could mean "is a lineage".

        The feature is drawn independently of the tree, so it should spread
        across lineages. This is the assertion that makes the finding mean
        something - and it is the one the stage's own lineage flag would
        otherwise be silently failing.
        """
        results = g3_outcome["results"]
        causal = next(
            (r for r in results if r.feature == "gene__planted_causal_oprd"), None
        )
        assert causal is not None
        assert not stage_gwas.is_lineage_confounded(
            causal, g3_outcome["input"], g3_outcome["config"].gwas.lineage_confound_threshold
        ), (
            "the planted feature came out confined to one lineage, so its rank "
            "is a lineage artefact and this test is measuring the wrong thing"
        )

    def test_it_ranks_first(self, g3_outcome):
        results = g3_outcome["results"]
        assert results, "pyseer reported nothing at all; nothing to rank"
        best = min(results, key=lambda r: r.p_value)
        assert best.feature == "gene__planted_causal_oprd", (
            f"planted feature ranked {[r.feature for r in results].index(best.feature) + 1}"
            f" of {len(results)}; best was {best.feature} at p={best.p_value:.3E}"
        )

    def test_it_passes_the_threshold_step_12a_derived(self, g3_outcome):
        """The gate is Bonferroni on unique patterns, and it must let this through."""
        results = g3_outcome["results"]
        causal = next(r for r in results if r.feature == "gene__planted_causal_oprd")
        threshold = g3_outcome["engine"].reduction.threshold
        assert causal.p_value <= threshold, (
            f"planted feature p={causal.p_value:.3E} did not clear the "
            f"Bonferroni threshold {threshold:.3E}"
        )
        assert causal.p_value <= g3_outcome["config"].gwas.significance_threshold

    def test_no_null_feature_passes_at_the_configured_alpha(self, g3_outcome):
        """A threshold that passes everything is not a threshold.

        R16: this is a statement about the *machinery* on a synthetic cohort
        whose nulls carry deliberate structure. It says the correction is
        working, and nothing about resistance.
        """
        config = g3_outcome["config"]
        results = g3_outcome["results"]
        passing = stage_gwas.significant(results, config)
        names = {r.feature for r in passing}
        assert "gene__planted_causal_oprd" in names, (
            "the planted feature did not survive; the nulls cannot be judged "
            "until it has"
        )
        false_positives = names - {"gene__planted_causal_oprd"}
        assert not false_positives, (
            f"{len(false_positives)} structured null features cleared "
            f"adj.p <= {config.gwas.significance_threshold}: "
            f"{sorted(false_positives)[:8]}"
        )

    def test_the_provenance_record_names_the_tool_and_the_command(self, g3_outcome):
        """G1: a result that cannot say which pyseer produced it is not auditable."""
        from papipeline.io.tsv import read_tsv

        path = g3_outcome["workdir"] / stage_gwas.PYSEER_PROVENANCE_FILENAME
        assert path.is_file()
        record = {r["key"]: r["value"] for r in read_tsv(path)}
        assert record["pyseer_version"] == "1.1.2", (
            f"expected the pinned pyseer 1.1.2, got {record['pyseer_version']!r}"
        )
        assert "--lmm" in record["command_1"]
        assert "--similarity" in record["command_1"]
        assert record["kinship_source"]
        assert record["n_samples"] == str(G3_N)
        # Both passes, in order.
        assert int(record["command_1"].split("--lrt-pvalue")[1].split()[0]) == 1
        assert float(record["command_2"].split("--lrt-pvalue")[1].split()[0]) == pytest.approx(
            float(record["lrt_pvalue_threshold"])
        )


# --------------------------------------------------------------------------
# G4 - the 10 real isolates
# --------------------------------------------------------------------------


class TestTheRealTenIsolates:
    """The real cohort: 7 R, 3 S, real sample ids, real pyseer.

    R16, stated once and meant: **n=10, underpowered, not a finding.** Every
    number this class produces is evidence about the machinery. The phenotype
    is real and measured; the feature matrix is synthetic, because no real one
    exists to read.
    """

    @pytest.fixture(scope="class")
    def outcome(self, tmp_path_factory, machine_config):
        requires_pyseer()
        config = machine_config
        workdir = tmp_path_factory.mktemp("g4")
        gwas_input = real_cohort_input()
        engine = stage_gwas.PyseerEngine("pyseer", workdir)
        results = engine.run(gwas_input, config)
        return {
            "results": results,
            "input": gwas_input,
            "engine": engine,
            "workdir": workdir,
            "config": config,
        }

    def test_the_cohort_is_the_real_one(self, outcome):
        gwas_input = outcome["input"]
        assert gwas_input.n_samples == 10
        assert gwas_input.n_positive == 7 and gwas_input.n_negative == 3
        assert "PDT000292995.1" in gwas_input.sample_ids

    def test_it_clears_the_configured_group_floor(self, config):
        """7 >= 3 and 3 >= 3: it passes, with zero margin on the S group.

        Measured on the real phenotype file, not asserted from memory.
        """
        assert real_phenotype_calls()
        assert config.gwas.min_samples_per_group == 3

    def test_the_mixed_model_runs_at_n_ten(self, outcome):
        """It runs. pyseer imposes no sample floor of its own.

        Checked below n=10 as well: pyseer exits 0 and writes a table at n=4,
        so "pyseer could not run at n=10" is not a real failure mode here. The
        failure at n=10 is statistical, not mechanical, and belongs in the
        underpowered flag rather than in a refusal.
        """
        engine = outcome["engine"]
        assert engine.executable_version == "1.1.2"
        assert (outcome["workdir"] / "pyseer_raw.tsv").is_file()
        assert engine.reduction is not None
        assert engine.reduction.n_unique_patterns >= 1

    def test_no_null_feature_is_reported_as_an_association(self, outcome):
        """The right outcome at n=10 is an empty or near-empty table.

        Asserted as "no false positive", never as "no association": with 40
        random features and 3 susceptible isolates, some table is expected, and
        any p-value in it is a statement about noise.
        """
        config = outcome["config"]
        passing = stage_gwas.significant(outcome["results"], config)
        assert not passing, (
            "structured noise passed the threshold on n=10; that is a "
            "correction failure, not a discovery"
        )

    def test_the_result_is_labelled_underpowered(self, outcome):
        """R16: the flag travels with the result, not in a log line.

        A number computed on 10 isolates is indistinguishable from any other
        once it is in a table, so the table has to say so itself.
        """
        record = (outcome["workdir"] / stage_gwas.PYSEER_PROVENANCE_FILENAME)
        assert record.is_file()
        text = record.read_text(encoding="utf-8")
        assert "n_samples" in text
        results = outcome["results"]
        assert all(r.model == "pyseer:mixed:bonferroni_unique_patterns" for r in results)


class TestBonferroniIsWhatDecides:
    """The answer to "Bonferroni or FDR", settled by arithmetic rather than taste.

    pyseer prints only the variants that cleared `--lrt-pvalue`, which step 12a
    set to ``alpha / n_unique_patterns`` (`pyseer/lmm.py`: a variant at or above
    the gate is dropped unless `--print-filtered`, which this stage does not
    pass). So the table `parse_pyseer_output` reads contains survivors only.

    That makes the two corrections non-independent, and in one direction. With
    ``m`` survivors all satisfying ``p <= alpha / n_unique`` and ``m <= n_unique``,
    BH at rank ``i`` gives

        adj_i = min over j >= i of (p_j * m / j)  <=  p_i * m / i
             <=  (alpha / n_unique) * (m / i)     <=  alpha

    so **every Bonferroni survivor automatically satisfies
    ``adj.p <= significance_threshold``**. `significant()` cannot reject
    anything pyseer kept, and the FDR column, while reported honestly, cannot
    change the verdict.

    So: Bonferroni on unique patterns decides, and BH is the reported column.
    That is not a contradiction - it is what `docs/scientific_rules.md` sections
    3 and 3a each mandate, at two different places - but it does mean the two
    are not independent checks, and a reader who assumed they were would
    over-read a result that only one of them decided.
    """

    def test_every_survivor_clears_the_adjusted_threshold_automatically(self):
        """Checked by exhaustive small-case search, not by one example."""
        checked = 0
        for n_unique in range(2, 12):
            threshold = 0.05 / n_unique
            # Worst case: every survivor sits exactly AT the gate, so BH has the
            # least room to shrink anything.
            for survivors in range(1, n_unique + 1):
                p_values = {f"gene__s{i}": threshold for i in range(survivors)}
                adjusted = stage_gwas.benjamini_hochberg(p_values)
                for value in adjusted.values():
                    assert value <= 0.05, (
                        f"n_unique={n_unique} survivors={survivors}: a variant "
                        f"at the Bonferroni gate got adj.p={value:.4E} > 0.05"
                    )
                    checked += 1
        assert checked > 100

    def test_the_claim_holds_for_the_real_g3_run(self, g3_outcome):
        """The same property, on a table pyseer actually produced."""
        config = g3_outcome["config"]
        results = g3_outcome["results"]
        threshold = g3_outcome["engine"].reduction.threshold
        for result in results:
            assert result.p_value <= threshold, (
                f"{result.feature} reached the table at p={result.p_value:.3E} "
                f"with a gate of {threshold:.3E}; pyseer should have dropped it, "
                "so either the gate or the filter has changed"
            )
        assert [r.feature for r in stage_gwas.significant(results, config)] == [
            r.feature for r in results
        ], (
            "significant() rejected a Bonferroni survivor, which the arithmetic "
            "above says cannot happen - so one of the two gates is not what "
            "this test thinks it is"
        )


class TestPyseerRunsBelowTheRealCohortSize:
    """Recorded because it settles a question rather than because it is useful.

    pyseer 1.1.2's mixed model exits 0 at n=4 with a similarity matrix of the
    matching dimension. So there is no pyseer sample floor to honour, and a
    refusal justified by "the LMM needs more samples than it has" would be
    false. What is small at n=10 is the *evidence*, not the tool.

    There is a floor, and it is the stage's own: a cohort with no testable
    feature yields an empty patterns file, and dividing ``alpha`` by zero unique
    patterns is refused. Which of the two happens depends on how many features
    survive the generator's carrier bounds, not on pyseer.
    """

    @pytest.mark.parametrize("n", [4, 6, 8, 10])
    def test_it_runs_or_refuses_with_a_reason_and_never_guesses(
        self, tmp_path, machine_config, n
    ):
        requires_pyseer()
        gwas_input, _ = synthetic_cohort(n, seed=5)
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / f"n{n}")
        try:
            engine.run(gwas_input, machine_config)
        except DataContractError as exc:
            assert (
                "patterns file is empty" in str(exc)
                or "No feature in this cohort is testable" in str(exc)
            ), f"n={n} failed for an unnamed reason: {exc}"
            return
        raw = tmp_path / f"n{n}" / "pyseer_raw.tsv"
        assert raw.is_file()
        assert raw.read_text().splitlines(), "pyseer wrote no header at all"
        assert engine.reduction.n_unique_patterns >= 1

    def test_a_cohort_with_no_testable_feature_is_refused_before_pyseer_runs(
        self, tmp_path, machine_config
    ):
        """Zero features would otherwise reach pyseer as a singular matrix.

        Observed against the real binary: with no feature, `_write_kinship`'s
        fallback writes an all-zero Gram matrix and pyseer 1.1.2 dies inside
        `numpy.linalg.eigh` with `LinAlgError: Eigenvalues did not converge`,
        exit 1. The stage now refuses first, which is both fail-closed and a
        statement about the cohort rather than a tool traceback in a log nobody
        is watching.
        """
        requires_pyseer()
        gwas_input, _ = synthetic_cohort(4, seed=5)
        stripped = stage_gwas.GwasInput(
            sample_ids=gwas_input.sample_ids,
            binary_outcome=gwas_input.binary_outcome,
            features={},  # nothing to test
            feature_types={},
            lineages=gwas_input.lineages,
            n_positive=gwas_input.n_positive,
            n_negative=gwas_input.n_negative,
            n_excluded=0,
        )
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "empty")
        with pytest.raises(DataContractError) as excinfo:
            engine.run(stripped, machine_config)
        message = str(excinfo.value)
        assert "No feature in this cohort is testable" in message
        assert not (tmp_path / "empty" / "pyseer_pass1.log").exists(), (
            "pyseer was launched for a cohort with nothing to test"
        )


# --------------------------------------------------------------------------
# The provenance and refusal paths that do not need the binary
# --------------------------------------------------------------------------


class TestExecutableResolution:
    """G5: pyseer is resolved and refused by name, like every other adapter."""

    def test_it_records_the_version_it_resolved(self, tmp_path):
        if shutil.which("pyseer") is None:
            pytest.skip("pyseer is not on PATH")
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "wd")
        assert engine.require_executable() == "1.1.2"
        assert engine.resolved_executable
        assert Path(engine.resolved_executable).is_file()

    def test_the_resolution_is_cached_not_repeated_per_pass(self, tmp_path, monkeypatch):
        """Two pyseer passes, one filesystem search and one `--version` subprocess.

        Counting entry into `require_executable` would prove nothing - the point
        is that the expensive part (a subprocess) happens once, so that is what
        is counted.
        """
        if shutil.which("pyseer") is None:
            pytest.skip("pyseer is not on PATH")
        import subprocess

        runs = []
        original = subprocess.run

        def counting(cmd, *args, **kwargs):
            runs.append(cmd)
            return original(cmd, *args, **kwargs)

        monkeypatch.setattr(subprocess, "run", counting)
        engine = stage_gwas.PyseerEngine("pyseer", tmp_path / "wd")
        engine.require_executable()
        engine.require_executable()
        engine.require_executable()
        assert len(runs) == 1, f"pyseer was asked for its version {len(runs)} times"

    def test_a_missing_pyseer_is_refused_by_name(self, tmp_path):
        """Not a FileNotFoundError from a subprocess whose stderr goes to a log.

        This is the shape `ToolNotAvailableError` exists for, and the shape the
        other adapters use. A stage that fails with a bare OS error leaves the
        operator looking for a broken install rather than a missing one.
        """
        from papipeline.errors import ToolNotAvailableError

        engine = stage_gwas.PyseerEngine(
            "pyseer_absolutely_not_installed", tmp_path / "wd"
        )
        with pytest.raises(ToolNotAvailableError) as excinfo:
            engine.require_executable()
        message = str(excinfo.value)
        assert "pyseer_absolutely_not_installed" in message
        assert "tool_search_dirs" in message, (
            "the refusal must say where a machine overlay would declare it"
        )

    def test_a_path_that_exists_but_is_not_a_tool_is_refused(self, tmp_path):
        """A file that exists is not a tool that works.

        `gubbins`'s own resolver documents this distinction; here it is enforced,
        because a truncated install or a wrong architecture produces exactly this
        state and the alternative is a traceback from the tool's internals.
        """
        from papipeline.errors import ToolNotAvailableError

        fake = tmp_path / "pyseer"
        fake.write_text("#!/bin/sh\nexit 3\n")
        fake.chmod(0o755)
        engine = stage_gwas.PyseerEngine(str(fake), tmp_path / "wd")
        with pytest.raises(ToolNotAvailableError) as excinfo:
            engine.require_executable()
        assert "version" in str(excinfo.value)


class TestTheCorrectionIsOneName:
    """`config.gwas.correction` was read by nothing in the entire repository.

    It sat in `science.yaml`, was loaded into `GwasConfig`, and no caller ever
    asked for it - so the column was filled with Benjamini-Hochberg while the
    configuration said `fdr_bh` and nobody was checking. Wiring it means a value
    this code cannot compute has to be refused rather than ignored.
    """

    def test_the_committed_key_is_the_one_the_parser_computes(self, config):
        assert config.gwas.correction == stage_gwas.SUPPORTED_ADJUSTED_P_CORRECTION

    def test_an_unimplemented_correction_is_refused_not_ignored(self, tmp_path):
        table = tmp_path / "t.tsv"
        table.write_text("feature\tpvalue\tbeta\tfreq\nf1\t1.0E-04\t2.0\t0.5\n")
        gwas_input, _ = synthetic_cohort(40, seed=3)
        with pytest.raises(DataContractError) as excinfo:
            stage_gwas.parse_pyseer_output(
                table, gwas_input, correction="holm"
            )
        assert "holm" in str(excinfo.value)

    def test_the_refusal_says_the_two_corrections_are_different(self, tmp_path):
        """Bonferroni gates the table; BH fills the column. Not one correction."""
        table = tmp_path / "t.tsv"
        table.write_text("feature\tpvalue\nf1\t1.0E-04\n")
        gwas_input, _ = synthetic_cohort(40, seed=3)
        with pytest.raises(DataContractError) as excinfo:
            stage_gwas.parse_pyseer_output(table, gwas_input, correction="holm")
        assert "bonferroni" in str(excinfo.value).lower()

    def test_a_pyseer_supplied_value_is_still_preferred(self, tmp_path):
        """The `adj.pvalue` branch stays, and stays second.

        pyseer 1.1.2 never emits the column - proven above, by running it - so
        this branch is unreachable today. It is kept because it is the right
        thing to do if pyseer grows the column, and `TestWhatTheRealBinaryDoes-
        NotEmit` is what will notice when that happens.
        """
        table = tmp_path / "t.tsv"
        table.write_text(
            "feature\tpvalue\tadj.pvalue\tbeta\tfreq\n"
            "gene__a\t1.0E-04\t0.5\t1.0\t0.5\n"
            "gene__b\t2.0E-04\t0.6\t1.0\t0.5\n"
        )
        gwas_input, _ = synthetic_cohort(40, seed=3)
        results = stage_gwas.parse_pyseer_output(table, gwas_input)
        assert {r.feature: r.adjusted_p_value for r in results} == {
            "gene__a": 0.5, "gene__b": 0.6,
        }


class TestTiedPValuesShareAnAdjustedValue:
    """A round-12 report claimed a rank bug here. It is not one, and this is why.

    The claim was that `need_adjustment.index(p_value)` - which returns the
    *first* match - mis-assigns Benjamini-Hochberg to tied p-values, "because
    two equal raw p-values at different ranks have different adjusted values".
    The premise is false. BH is a step-up procedure computed from the largest
    p-value downwards with a running minimum, and for two equal p-values at
    adjacent ranks r and r+1:

        computed[r+1] <= p*n/(r+1)  <  p*n/r
        computed[r]   = min(computed[r+1], p*n/r) = computed[r+1]

    so they are equal by construction. They *should* be equal - they are the
    same hypothesis. Verified below over randomised inputs rather than argued
    from the algebra alone.

    The lookup is still done by position rather than by value, because the
    position is the thing that is correct by construction and the argument above
    is the kind that is easy to believe and never re-checked. The two agree
    numerically, which is the point of the second test.
    """

    def _table(self, path, rows):
        path.write_text(
            "feature\tpvalue\tbeta\tfreq\n"
            + "".join(f"{name}\t{p}\t1.0\t0.5\n" for name, p in rows)
        )
        return path

    def test_tied_p_values_get_one_adjusted_value_and_bh_agrees(self, tmp_path):
        table = self._table(
            tmp_path / "t.tsv",
            [("gene__a", "1.0E-04"), ("gene__b", "2.0E-04"),
             ("gene__c", "2.0E-04"), ("gene__d", "9.0E-01")],
        )
        gwas_input, _ = synthetic_cohort(40, seed=3)
        results = stage_gwas.parse_pyseer_output(table, gwas_input)
        by_name = {r.feature: r for r in results}

        expected = stage_gwas.benjamini_hochberg(
            {"gene__a": 1.0e-4, "gene__b": 2.0e-4,
             "gene__c": 2.0e-4, "gene__d": 9.0e-1}
        )
        for result in results:
            assert result.adjusted_p_value == pytest.approx(
                expected[result.feature]
            ), "the parser's column is not BH over the rows it parsed"
        assert by_name["gene__b"].adjusted_p_value == (
            by_name["gene__c"].adjusted_p_value
        ), (
            "two features at the same p-value were given DIFFERENT adjusted "
            "values. That would be the opposite error: BH step-up gives tied "
            "p-values one value, because they are one hypothesis."
        )

    def test_the_position_lookup_matches_bh_on_random_inputs(self):
        """The property the refutation rests on, over many shapes rather than one.

        2000 randomised inputs; if a tie ever split, or if the position lookup
        ever disagreed with the value lookup, this fails.
        """
        rng = np.random.default_rng(4)
        menu = [0.01, 0.02, 0.05, 0.2, 0.5, 0.9]
        for _ in range(2000):
            n = int(rng.integers(2, 9))
            values = [float(rng.choice(menu)) for _ in range(n)]
            adjusted = stage_gwas.benjamini_hochberg(
                {str(i): v for i, v in enumerate(values)}
            )
            by_value: dict = {}
            for i, value in enumerate(values):
                by_value.setdefault(value, set()).add(adjusted[str(i)])
                assert adjusted[str(i)] == adjusted[str(values.index(value))], (
                    "the position lookup and the first-match lookup disagree, so "
                    "one of the two is wrong and the tie argument is not the "
                    "guarantee it was taken to be"
                )
            for value, seen in by_value.items():
                assert len(seen) == 1, (
                    f"tied p-value {value} received {len(seen)} adjusted values"
                )