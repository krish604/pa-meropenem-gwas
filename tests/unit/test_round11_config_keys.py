"""The five configuration keys added in round 11, and what each one is FOR.

Each key is asserted three ways, because a key with only a DEFAULT is a comment
and a key with only a YAML value is a typo:

  * the value this repository's ``config/science.yaml`` actually carries, read
    through the loader, so the assertion tracks the file rather than a copy of
    it frozen into the test;
  * the loader's DEFAULT, asserted against a config built from an EMPTY
    section, so the fallback is pinned and cannot drift silently; and
  * that the key REACHES something. A key that is parsed and then read by nobody
    is the accepted-and-ignored failure this repository treats as a defect, so
    each of the two keys that claim to influence a tool invocation is checked
    against the argv that tool is actually handed.

The keys, and the rulings that put them here:

  1. ``analysis.recombination``          stage 8 was switched on by nothing
  2. ``variants.bcftools.mpileup.min_ireads``   ruling R6, depth-1 assemblies
  3. ``lineage.method``                  ruling R4, ST from stage 3
  4. ``cohort.subset_file``              nullable; a bounded run's member list
  5. ``phylogeny.models`` + ``asc_drop_partially_constant``   ruling R8

One of the five is still asserted in its REVERTED state, and the test says so:
``phylogeny.models`` (item 5) was added in the previous round and broke the
suite, because it declared a half of a ruling whose other half has no
implementation. Re-enabling it is a branch's job and must land in the same
commit as the code that makes it true.

``analysis.recombination`` (item 1) was reverted for the same reason and is now
asserted in its ENABLED state, because the dispatch branch that makes the key
true landed in the same commit. The test class for it is the worked example of
the rule this file is here to enforce.
"""

from __future__ import annotations

import contextlib
import os
import tempfile
from pathlib import Path

import pytest

from papipeline.adapters.gubbins import require_real_gate
from papipeline.config import (
    ACCEPTED_LINEAGE_METHODS,
    DEFAULT_COHORT_SUBSET_FILE,
    DEFAULT_LINEAGE_METHOD,
    DEFAULT_MPILEUP_MIN_IREADS,
)
from papipeline.config.loader import (
    CohortConfig,
    ConfigError,
    LineageConfig,
    PhylogenyConfig,
    load_config,
    load_machine_config,
)
from papipeline.errors import ModeNotAllowedError

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
LAPTOP = REPO / "config" / "machines" / "laptop.yaml"


@pytest.fixture(scope="module")
def config():
    return load_config(SCIENCE, machine="laptop")


@contextlib.contextmanager
def _chdir(where):
    """Run the block with a different working directory, then restore."""
    previous = os.getcwd()
    os.chdir(where)
    try:
        yield
    finally:
        os.chdir(previous)


# ---------------------------------------------------------------------------
# 1. analysis.recombination
# ---------------------------------------------------------------------------


class TestStageEightIsSwitchedOnBySomething:
    """`recombination` is enabled, and the key and the dispatch landed together.

    History, because the reason is not visible in the config alone. Round 11
    first ADDED `analysis.recombination: true`, on the reasoning that an absent
    key is not neutral: `PipelineConfig.stage_enabled` is
    `self.analysis.get(stage, False)`, so an absent key reads as False and
    `run.py` skipped stage 8 with no refusal and no warning.

    That reasoning was right about the skip and wrong about the remedy. Measured
    at that HEAD, adding the key alone accounted for **66 failing/erroring suite
    nodes**: `recombination` was ALSO in `run.UNBUILT_STAGES`, and `enabled()`
    refuses UNBUILT stages for every mode except STUB, so the key did not make
    stage 8 run -- it turned a silent skip into a hard `NotImplementedError` in
    TEST. The key was reverted; the stage then had to be BUILT.

    Both halves are now in place, in one commit: `run.derive_recombination_tables`
    dispatches the stage, `recombination` is out of `UNBUILT_STAGES`, and
    `analysis.recombination: true` is in `config/science.yaml`. These tests
    assert all three together, so the fourth combination - a key with no
    branch, or a branch with no key - cannot be committed.

    **Enabling the stage is not permission to run REAL data.** That is
    `runtime.allow_real_mode`, false on all three committed overlays, and stage
    8's REAL body refuses while it is shut. Asserted below.
    """

    def test_the_key_is_present_and_true(self, config):
        assert config.stage_enabled("recombination") is True, (
            "analysis.recombination must be true: the stage is dispatched, and "
            "`enabled()` falls back to the flag now that recombination is out of "
            "UNBUILT_STAGES. Without the key the stage is skipped - silently, "
            "which is the defect this commit closes."
        )

    def test_the_raw_key_is_really_present_not_merely_defaulted(self, config):
        """Assert the file states it, so a loader DEFAULT cannot mask its return."""
        assert config.raw["analysis"].get("recombination") is True, (
            "analysis.recombination must be stated in config/science.yaml, not "
            "supplied by some other layer; a key with only a default is a comment"
        )

    def test_an_absent_stage_key_is_false_and_not_an_error(self, config):
        """The defaulting direction is what made the omission quiet."""
        assert config.stage_enabled("no_such_stage") is False

    def test_the_stage_is_dispatched_and_not_unbuilt(self):
        """The other two halves, on the same commit as the key.

        `tests/integration/test_stage_taxonomy_is_self_verifying.py` derives
        both directions of this from `run.py`'s AST; this is the assertion that
        ties them to *this* key, so reverting any one of the three is caught
        here by name.
        """
        from papipeline.run import UNBUILT_STAGES, derive_recombination_tables

        assert "recombination" not in UNBUILT_STAGES, (
            "stage 8 is dispatched. Leaving it in UNBUILT_STAGES means "
            "`enabled()` refuses it outside STUB, so analysis.recombination: true "
            "raises NotImplementedError in TEST instead of running the stage."
        )
        assert callable(derive_recombination_tables), (
            "the dispatch branch calls run.derive_recombination_tables; without "
            "it the key enables a stage that has no branch and never runs"
        )

    def test_enabling_the_stage_does_not_enable_real_mode(self, config):
        """The two gates are different, and conflating them is the trap.

        `analysis.recombination: true` says the stage exists and is wanted.
        `runtime.allow_real_mode` says the cohort behind it may be read. Only
        the second is false here, and it is what stage 8's REAL body checks
        first.
        """
        assert bool(config.runtime.get("allow_real_mode", False)) is False
        with pytest.raises(ModeNotAllowedError) as caught:
            require_real_gate(config)
        assert "allow_real_mode" in str(caught.value)


# ---------------------------------------------------------------------------
# 2. variants.bcftools.mpileup.min_ireads
# ---------------------------------------------------------------------------


class TestMpileupMinIreads:
    """Ruling R6. `-m/--min-ireads` counts GAPPED READS, not depth.

    Verified against the pinned tool rather than remembered::

        $ bcftools mpileup | grep min-ireads
        -m, --min-ireads INT    Minimum number gapped reads for indel
                                candidates [2]

    An assembly is one consensus sequence, so there is exactly one gapped read
    per indel. At bcftools' default of 2 no assembly indel can qualify, so
    every one is dropped before it is considered and the locus reads as clean.
    """

    def test_the_science_file_states_one(self, config):
        assert config.variants_mpileup_min_ireads() == 1

    def test_the_default_is_one(self):
        assert DEFAULT_MPILEUP_MIN_IREADS == 1

    def test_an_absent_key_yields_the_default_not_an_error(self):
        """A run with no `variants` block at all is legitimate, not malformed.

        Built from a real `PipelineConfig` rather than a hand-rolled stub, so the
        accessor is exercised against the class the pipeline actually uses.
        """
        from papipeline.config.loader import (
            ConvergenceConfig,
            GwasConfig,
            PipelineConfig,
            QcConfig,
        )

        bare = PipelineConfig(
            root=REPO, raw={}, antibiotics=(), allowed_phenotypes=(),
            analysis={}, mechanisms={}, regulators={}, references={},
            antibiotic_specs={}, organism={},
            qc=QcConfig.from_dict({}),
            gwas=GwasConfig.from_dict(
                {"unique_patterns": {"method": "bonferroni", "alpha": 0.05}}
            ),
            phylogeny=PhylogenyConfig.from_dict({}),
            convergence=ConvergenceConfig.from_dict({}),
            paths={}, runtime={},
        )
        assert bare.variants_mpileup_min_ireads() == DEFAULT_MPILEUP_MIN_IREADS

    def test_the_value_reaches_the_bcftools_argv(self, config):
        """Parsed-and-ignored is a defect, so assert the command itself.

        At the pre-round-11 HEAD this was NOT PROVEN: `mpileup_command` emitted
        `-f -q -Q -a -d` and no `-m` at all, so the pinned bcftools default of 2
        applied whatever the configuration said. This test is the regression
        guard for that.
        """
        from papipeline.adapters.minimap2 import mpileup_command

        command = mpileup_command(
            "bcftools", bam=Path("in.bam"), reference=Path("ref.fa"),
            min_mapq=20, min_bq=20, annotate="FORMAT/DP,FORMAT/AD",
            max_depth=250,
            min_ireads=config.variants_mpileup_min_ireads(),
        )
        assert "-m" in command, (
            "bcftools mpileup is missing --min-ireads; the pinned default of 2 "
            "gapped reads then applies, which no assembly indel can satisfy"
        )
        assert command[command.index("-m") + 1] == "1", (
            f"-m is not threaded from configuration: {command}"
        )

    def test_the_argv_default_is_the_same_number_as_the_loader_default(self):
        """One default, stated once. Two would drift."""
        from papipeline.adapters.minimap2 import mpileup_command

        command = mpileup_command(
            "bcftools", bam=Path("in.bam"), reference=Path("ref.fa"),
            min_mapq=20, min_bq=20, annotate="x", max_depth=250,
        )
        assert command[command.index("-m") + 1] == str(DEFAULT_MPILEUP_MIN_IREADS)

    def test_min_ireads_is_distinct_from_the_depth_flag(self, config):
        """`-d` is depth. Conflating the two is the misreading R6 corrects."""
        from papipeline.adapters.minimap2 import mpileup_command

        command = mpileup_command(
            "bcftools", bam=Path("in.bam"), reference=Path("ref.fa"),
            min_mapq=20, min_bq=20, annotate="x", max_depth=250,
        )
        assert command[command.index("-d") + 1] == "250"
        assert command[command.index("-m") + 1] == "1"


# ---------------------------------------------------------------------------
# 3. lineage.method
# ---------------------------------------------------------------------------


class TestLineageMethod:
    """Ruling R4: `lineage_label` IS the MLST sequence type from stage 3.

    `tree_cut` is DEFERRED, so it is not in the accepted set: an
    accepted-values set is what stops a deferred method from being spelled into
    a config file and then believed by whoever reads the next result.
    """

    def test_the_science_file_states_st(self, config):
        assert config.lineage.method == "st"

    def test_the_default_is_st(self):
        assert LineageConfig.from_dict({}).method == DEFAULT_LINEAGE_METHOD == "st"

    def test_an_empty_section_is_not_an_error(self):
        assert LineageConfig.from_dict({}).method == "st"

    def test_tree_cut_is_not_accepted(self):
        """Deferred means unavailable, not quietly tolerated."""
        with pytest.raises(ConfigError) as excinfo:
            LineageConfig.from_dict({"method": "tree_cut"})
        assert "tree_cut" in str(excinfo.value)

    def test_the_accepted_set_holds_exactly_one_method(self):
        assert ACCEPTED_LINEAGE_METHODS == ("st",)


# ---------------------------------------------------------------------------
# 4. cohort.subset_file
# ---------------------------------------------------------------------------


class TestCohortSubsetFile:
    """Nullable, and the null case is "the whole roster", not "no samples".

    A subset names cohort members, and a member list is sample-level
    information, which `config/science.yaml` explicitly must not carry. That is
    why the key holds a PATH and why the file it names lives outside
    configuration.
    """

    def test_the_science_file_states_null(self, config):
        assert config.cohort.subset_file is None

    def test_the_default_is_none(self):
        assert DEFAULT_COHORT_SUBSET_FILE is None
        assert CohortConfig.from_dict({}).subset_file is None

    def test_an_explicit_null_is_the_same_as_absent(self):
        assert CohortConfig.from_dict({"subset_file": None}).subset_file is None

    def test_a_path_is_kept_verbatim(self):
        assert CohortConfig.from_dict(
            {"subset_file": "local/smoke_isolates.txt"}
        ).subset_file == "local/smoke_isolates.txt"

    def test_science_states_null_but_the_smoke_overlay_names_the_list(self):
        """The split: science says "whole roster", smoke says which ten.

        Asserted on the real overlays rather than on a copy, so the smoke run's
        ten cannot silently become the full cohort -- the one outcome that
        overlay exists to prevent.
        """
        assert load_config(SCIENCE, machine="laptop").cohort.subset_file is None
        assert load_config(SCIENCE, machine="bigmachine").cohort.subset_file is None

        smoke = load_config(SCIENCE, machine=REPO / "config" / "machines" / "smoke.yaml")
        assert smoke.cohort.subset_file == "local/smoke_isolates.txt"
        assert smoke.cohort.resolve_subset_file() == REPO / "local" / "smoke_isolates.txt"

    def test_an_empty_string_is_refused_rather_than_read_as_no_subset(self):
        """`null` and `""` mean different things; only one of them is valid.

        An empty subset selects zero samples, and a run that analyses nothing
        while reporting success is precisely what this key must not be able to
        express.
        """
        with pytest.raises(ConfigError) as excinfo:
            CohortConfig.from_dict({"subset_file": ""})
        assert "cohort.subset_file" in str(excinfo.value)

    def test_the_subset_list_is_never_committed(self):
        """The bounded run's member list is an input, not source.

        Asserted by SHELLING OUT to `git check-ignore`, because unlike the `data`
        rule this one is checkable: `local/` is a REAL directory in every
        worktree, so git will descend into it. The previous location under
        `data/` could not be checked at all --

            $ git check-ignore -q data/smoke_isolates.txt ; echo $?
            fatal: pathspec 'data/smoke_isolates.txt' is beyond a symbolic link
            128

        -- and an ignore rule nobody can check is an ignore rule nobody has
        checked. So this asserts exit 0 from the real command.
        """
        import subprocess

        result = subprocess.run(
            ["git", "check-ignore", "-v", "--", "local/smoke_isolates.txt"],
            cwd=REPO, capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            "git check-ignore does not match local/smoke_isolates.txt, so the "
            f"bounded run's isolate list is committable: {result.stdout}"
        )
        assert result.stdout.split("\t")[0].rstrip().endswith("local"), (
            "the list is matched by some rule other than the `local` rule, so "
            "the bare-`local` convention that also hides the directory itself is "
            f"not what is protecting it: {result.stdout!r}"
        )

    def test_the_ignore_rule_has_no_trailing_slash(self):
        """A trailing slash would match the directory but not its contents.

        Same reason the `data` and `db` rules omit it, `local/` is a REAL
        directory, so `local/` with the slash and `local` without it are both
        legal gitignore patterns that behave differently for a directory. The
        bare form is what is asserted.
        """
        rules = {
            line.strip()
            for line in (REPO / ".gitignore").read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        }
        assert "local" in rules, (
            "the `local` rule must exist and must have NO trailing slash, so "
            "local/smoke_isolates.txt is ignored as a path and not only as a "
            "member of a directory"
        )

    def test_a_relative_subset_resolves_against_the_repo_root_not_the_cwd(self):
        """The resolution, pinned. A TEMP path -- never the real list.

        This is the assertion that would fail if the loader went back to
        `Path.cwd()`: the CWD here is a TEMP directory that contains no
        `config/` and no `local/`, so a CWD-relative implementation would resolve
        this to a path that does not exist, and one that walked UP to the repo
        root would resolve it to the real `local/smoke_isolates.txt` -- a
        different file with different members. The temp directory is never
        written to by the resolution and never read from the real list.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_file = Path(tmp) / "isolates.txt"
            tmp_file.write_text("S1\nS2\n", encoding="utf-8")

            with _chdir(tmp):
                relative = os.path.relpath(tmp_file, start=REPO)
                resolved = CohortConfig.from_dict(
                    {"subset_file": relative}
                ).resolve_subset_file()

            assert resolved == tmp_file, (
                "a relative cohort.subset_file must resolve against the "
                f"repository root ({REPO}), not the working directory ({tmp}); "
                f"got {resolved}"
            )
            assert resolved.read_text(encoding="utf-8") == "S1\nS2\n", (
                f"resolution landed on {resolved}, which is not this test's own "
                "fixture -- it would mean the real member list was read"
            )

    def test_an_absolute_subset_is_left_alone(self):
        """.resolve() must not re-root a path that is already absolute."""
        with tempfile.TemporaryDirectory() as tmp:
            absolute = str(Path(tmp) / "isolates.txt")
            assert CohortConfig.from_dict(
                {"subset_file": absolute}
            ).resolve_subset_file() == Path(absolute)

    def test_no_subset_resolves_to_none_rather_than_to_the_repo_root(self):
        """The null case must not become a path pointing at the repository."""
        assert CohortConfig.from_dict({"subset_file": None}).resolve_subset_file() is None
        assert CohortConfig.from_dict({}).read_subset_file() is None

    def test_a_missing_list_refuses_instead_of_falling_back_to_the_whole_roster(self):
        """The failure mode the key exists to prevent, asserted directly.

        A list that is absent must not read as "no subset", because that
        analyses the full cohort while appearing to be the bounded run.
        """
        with tempfile.TemporaryDirectory() as tmp:
            missing = str(Path(tmp) / "not_there.txt")
            with pytest.raises(ConfigError) as excinfo:
                CohortConfig.from_dict({"subset_file": missing}).read_subset_file()
            assert "subset_file" in str(excinfo.value)

    def test_a_temp_list_reads_back_its_ids(self):
        """Reading is what a run will do; prove it on a temp file only."""
        with tempfile.TemporaryDirectory() as tmp:
            list_file = Path(tmp) / "isolates.txt"
            list_file.write_text(
                "# a comment\nS1\n\nS2\n", encoding="utf-8"
            )
            cfg = CohortConfig.from_dict({"subset_file": str(list_file)})
            assert cfg.read_subset_file() == ["S1", "S2"]

    def test_the_subset_list_is_not_tracked_today(self):
        """The real repository state, asserted for the path as git sees it now.

        `local/` is a real directory, so this one CAN be asked of git directly --
        which is the whole reason the list moved there.
        """
        import subprocess

        listed = subprocess.run(
            ["git", "ls-files", "--", "data", "local"],
            cwd=REPO, capture_output=True, text=True, check=True,
        ).stdout.strip()
        assert listed == "", (
            f"git is tracking {listed!r}; the dataset, the `db` symlink and the "
            "bounded run's isolate list are inputs, not source"
        )


# ---------------------------------------------------------------------------
# 6. the smoke overlay's one exemption from the science/machine split
# ---------------------------------------------------------------------------


class TestTheOverlayExemptionIsOneOverlayAndOneSection:
    """`load_machine_config` refuses any science-owned section in an overlay.

    That guard is real and it fired: adding `cohort:` to
    `config/machines/smoke.yaml` was rejected at `loader.py` with
    "A machine overlay may not declare sections that belong to the science
    configuration (cohort)". So the exemption added for the bounded smoke run is
    `OVERLAY_SCIENCE_EXCEPTIONS`, keyed on the overlay's own `machine:` name, and
    these tests hold it to one overlay and one section.
    """

    def test_the_exemption_names_exactly_one_overlay_and_one_section(self):
        from papipeline.config.loader import OVERLAY_SCIENCE_EXCEPTIONS

        assert OVERLAY_SCIENCE_EXCEPTIONS == {"smoke": frozenset({"cohort"})}

    @pytest.mark.parametrize("machine", ["laptop", "bigmachine"])
    def test_the_other_two_machines_are_still_refused(self, machine, tmp_path):
        """The invariant the guard exists for: two machines cannot disagree.

        Built as a real overlay file on disk and pointed at the real
        science.yaml, so the refusal comes from the loader and not from a
        hand-raised error in the test.
        """
        overlay = tmp_path / f"{machine}.yaml"
        overlay.write_text(
            f"machine: {machine}\n"
            "machine_class: development\n"
            "cohort:\n"
            "  subset_file: local/smoke_isolates.txt\n"
            "runtime:\n"
            "  threads: 4\n"
            "  memory_mb: 8192\n",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError) as excinfo:
            load_machine_config(overlay)
        assert "cohort" in str(excinfo.value)

    def test_the_smoke_exemption_does_not_open_another_section(self, tmp_path):
        """One section wide. `smoke` may not use the door to move the science."""
        from papipeline.config.loader import OVERLAY_SCIENCE_EXCEPTIONS

        overlay = tmp_path / "smoke.yaml"
        overlay.write_text(
            "machine: smoke\n"
            "machine_class: development\n"
            "phylogeny:\n"
            "  models: GTR+G\n"
            "runtime:\n"
            "  threads: 4\n"
            "  memory_mb: 8192\n",
            encoding="utf-8",
        )
        assert OVERLAY_SCIENCE_EXCEPTIONS["smoke"] == frozenset({"cohort"})
        with pytest.raises(ConfigError) as excinfo:
            load_machine_config(overlay)
        assert "phylogeny" in str(excinfo.value)

    def test_an_exempted_section_is_still_validated(self, tmp_path):
        """An exemption must not become a hole with a comment on it.

        `subset_file: ""` is refused by CohortConfig wherever it is written --
        in the science file or in the one overlay allowed to state it -- because
        an empty member list selects zero samples.
        """
        overlay = tmp_path / "smoke.yaml"
        overlay.write_text(
            "machine: smoke\n"
            "machine_class: development\n"
            "cohort:\n"
            '  subset_file: ""\n'
            "runtime:\n"
            "  threads: 4\n"
            "  memory_mb: 8192\n",
            encoding="utf-8",
        )
        with pytest.raises(ConfigError) as excinfo:
            load_machine_config(overlay)
        assert "subset_file" in str(excinfo.value)

    def test_the_real_smoke_overlay_still_declares_its_own_genome_directory(self):
        """The exemption added a capability; it must not have displaced one.

        `paths.smoke_genome_dir` is what stops a bounded run resolving genomes
        from `data/` and analysing the full cohort. If adding `cohort` had cost
        that key, the two doors would be one hole.
        """
        import yaml

        raw = yaml.safe_load(
            (REPO / "config" / "machines" / "smoke.yaml").read_text(encoding="utf-8")
        )
        assert raw["paths"]["smoke_genome_dir"] == "db/smoke_genomes"
        assert raw["runtime"]["max_samples"] == 10
        assert raw["cohort"]["subset_file"] == "local/smoke_isolates.txt"


# ---------------------------------------------------------------------------
# 5. phylogeny.models and phylogeny.asc_drop_partially_constant
# ---------------------------------------------------------------------------


class TestPhylogenyModelAndAscPreprocessing:
    """Ruling R8: +ASC, on an alignment whose constant-after-gaps columns go.

    The model is measured, not chosen -- pa-artifacts/phy2/ratios_44pairs.txt:54
    over 44 pair-alignments at equal base counts:

        GTR+G      -> GTR+F+G4      lnL -376357.853   alpha spread 1.2527
        GTR+G+ASC  -> GTR+F+ASC+G4  lnL -354053.643   alpha spread 0.3423

    R8's two halves are COUPLED, and the filter half now has an implementation:
    `adapters/iqtree.py` `drop_partially_constant_columns` runs before the
    search, so the pinned IQ-TREE no longer sees the columns it would otherwise
    reject +ASC over:

        ERROR: Invalid use of +ASC because of 20 invariant sites in the alignment

    So `models` carries `GTR+G+ASC` again, in the same commit that builds the
    filter - the coupling is the point, and one half without the other is what
    the earlier revert was correcting. These tests pin BOTH halves together.
    """

    def test_the_science_file_states_the_asc_model(self, config):
        assert config.phylogeny.models == "GTR+G+ASC", (
            "phylogeny.models must be the +ASC variant phy2 measured "
            "(pa-artifacts/phy2/ratios_44pairs.txt:54, lnL -354053.643 against "
            "-376357.853 for GTR+G), and only once the constant-after-gaps "
            "columns are actually removed -- the pinned IQ-TREE exits 2 on "
            "+ASC otherwise"
        )

    def test_the_asc_model_is_stated_because_the_filter_is_live(self, config):
        """The two halves are coupled; stating half of R8 is what broke it."""
        assert "+ASC" in config.phylogeny.models, (
            "+ASC is the measured model, and `asc_drop_partially_constant` is "
            "the filter that must run FIRST for it to be estimable"
        )
        assert config.phylogeny.asc_drop_partially_constant is True, (
            "the filter that makes +ASC estimable is off; without it the pinned "
            "IQ-TREE refuses the model by name"
        )

    def test_the_model_reaches_the_iqtree_argv_as_dash_m(self, config):
        """`models` is consumed: phylogeny.py passes it to build_command's -m."""
        from papipeline.adapters.iqtree import build_command

        command = build_command(
            "iqtree", alignment=Path("core.fa"), prefix=Path("p"),
            model=config.phylogeny.models, threads=4, seed=1,
        )
        assert command[command.index("-m") + 1] == "GTR+G+ASC"

    def test_asc_drop_partially_constant_is_stated_true(self, config):
        assert config.phylogeny.asc_drop_partially_constant is True

    def test_the_loader_default_is_false(self):
        """The conservative direction: an unconfigured caller gets the
        alignment byte-for-byte as committed, so a pinned fixture still is the
        alignment that gets analysed."""
        assert PhylogenyConfig.from_dict({}).asc_drop_partially_constant is False

    def test_it_is_not_the_same_switch_as_models(self, config):
        """Two keys, two jobs. Collapsing them would hide which one is wired.

        `models` reaches the argv. This one is a preprocessing decision the
        pipeline makes itself, so it reaches `build_tree` as an argument and
        never as an iqtree flag - asserted in
        `tests/unit/test_iqtree_adapter.py`, which is where the filter is.
        """
        from papipeline.adapters.iqtree import build_command

        command = build_command(
            "iqtree", alignment=Path("core.fa"), prefix=Path("p"),
            model=config.phylogeny.models, threads=1, seed=1,
        )
        assert "-m" in command
        # Still no `-fconst`, and still deliberately: the "constant ignoring
        # gaps/N" test is not `-fconst`, and inventing an iqtree flag for it
        # without a verified tool reference is exactly what AGENTS.md rule 1
        # forbids. The pipeline filters the alignment itself instead.
        assert "-fconst" not in command
