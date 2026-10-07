"""The real-path pyseer adapter: its flags, its gates, and its two-pass protocol.

What is under test is :mod:`papipeline.gwas_real.adapter` - the command builder
and runner for the REAL pyseer path, plus the ``gwas_real:`` settings section it
reads. Three decisions are asserted here rather than described:

* **The flags are the installed pyseer's own.** Every token the builder emits
  was checked against pyseer 1.1.2's ``--help`` and ``pyseer/__main__.py``.
  The combinations pyseer *rejects* are refused here, before a subprocess is
  launched: ``--burden`` without ``--vcf`` exits at ``__main__.py:212``,
  ``--lmm --lineage`` without ``--distances`` exits at ``:220`` ("Must also
  provide a distance matrix to report lineage effects"), ``--lmm --distances``
  without ``--lineage`` exits at ``:216``, and a phenotype whose values are not
  0/1 makes pyseer print ``Detected continuous phenotype`` at ``:244`` and fit
  a *linear* model - silently, with exit 0.

* **The threshold comes from the existing reduction.** Step 12a's
  ``stages.gwas.reduce_unique_patterns`` (which runs the vendored
  ``scripts/gwas/count_patterns.py``) is called - spied on, not read about -
  and the number it returns is what reaches pass 2's ``--lrt-pvalue``. Nothing
  here recomputes ``alpha / n``.

* **The mode gates hold.** REAL refuses while ``runtime.allow_real_mode`` is
  shut; STUB refuses because a header-only mode cannot produce a kinship matrix
  or an association table; TEST runs on the committed synthetic fixtures with an
  injected invoker and never launches a tool unless ``PAPIPELINE_TEST_PYSEER=1``
  says so explicitly (spec.md D8).
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

from papipeline.errors import DataContractError
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample

from papipeline.gwas_real import adapter
from papipeline.gwas_real.settings import AdapterSettings

REPO = Path(__file__).resolve().parents[2]

#: The committed stage-9 fixture tree (20 tips, every branch 0.05).
FIXTURE_TREE = REPO / "test_data" / "phylogeny" / "tree.nwk"

#: The committed phenotype fixture: 7 R, 10 S, plus one I, one SDD and one ND.
FIXTURE_PHENOTYPE = REPO / "test_data" / "phenotype" / "imipenem_phenotype.tsv"

#: The opt-in switch, named as a literal rather than imported from the adapter:
#: a test that borrowed the production constant could not notice it changing.
OPT_IN = "PAPIPELINE_TEST_PYSEER"

#: pyseer 1.1.2's own association-table header for an ``--lmm`` run.
PYSEER_LMM_HEADER = (
    "variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\tvariant_h2\tnotes"
)

#: The same header for an ``--lmm --lineage`` run, which is what the committed
#: ``gwas_real.lineage: true`` asks for. ``pyseer/__main__.py`` builds it in
#: order: the six fixed columns, then ``+ ['variant_h2']`` under ``--lmm``,
#: then ``+ ['lineage']`` under ``--lineage``, then ``+ ['notes']``.
PYSEER_LMM_LINEAGE_HEADER = (
    "variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tbeta-std-err\t"
    "variant_h2\tlineage\tnotes"
)

#: ``config.gwas.min_samples_per_group`` in the committed science file.
MIN_PER_GROUP = 3


# --------------------------------------------------------------------------
# The synthetic cohort: phenotype, kinship, distances, variants
# --------------------------------------------------------------------------


def binary_pairs():
    """The fixture's binary cohort: R -> 1, S -> 0, I/SDD/ND excluded.

    The exclusion is spec.md D4's - the binary R-vs-S model cannot represent an
    intermediate - and it is why the analysable cohort is 18 rather than 20.
    Whoever writes the real phenotype file does this conversion; this test
    mirrors it so the adapter is exercised over a cohort that matches the
    pipeline's own rule rather than an invented one.
    """
    pairs = []
    for line in FIXTURE_PHENOTYPE.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] == "sample_id":
            continue
        if parts[2] == "R":
            pairs.append((parts[0], 1))
        elif parts[2] == "S":
            pairs.append((parts[0], 0))
    return pairs


def manifest_over(ids):
    return SampleManifest([Sample(sample_id=s) for s in ids])


def write_phenotype(path: Path, pairs) -> Path:
    """The adapter's phenotype input: ``sample`` / ``phenotype``, values 0/1."""
    lines = ["sample\tphenotype"] + [f"{s}\t{v}" for s, v in pairs]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def feature_vectors(n: int):
    """8 features over ``n`` samples with exactly 5 distinct patterns.

    Distinct patterns, not distinct rows, are what the threshold divides by,
    so a cohort whose duplicates are visible in the data exercises the
    reduction rather than a division by the feature count.
    """
    p0 = [i % 2 for i in range(n)]
    p1 = [1 - (i % 2) for i in range(n)]
    p2 = [1 if i < max(1, n // 3) else 0 for i in range(n)]
    p3 = [1 if i >= n - max(1, n // 3) else 0 for i in range(n)]
    p4 = [1 if i % 3 == 0 else 0 for i in range(n)]
    return {
        "gene__f0": p0,
        "gene__f1": p1,
        "gene__f2": p2,
        "gene__f3": p3,
        "gene__f4": p4,
        "gene__f5": list(p0),  # same pattern as f0
        "gene__f6": list(p1),  # same pattern as f1
        "gene__f7": list(p2),  # same pattern as f2
    }


def write_rtab(path: Path, ids) -> Path:
    """A roary/piggy-style presence matrix: ``Gene`` label, samples as columns.

    The label cell must be non-empty - pyseer reads the header with
    ``str.split()`` and drops the first token, so a leading tab shifts every
    column. Same rule as ``PyseerEngine._write_rtab``.
    """
    ids = list(ids)
    features = feature_vectors(len(ids))
    lines = ["Gene\t" + "\t".join(ids)]
    for name, vector in features.items():
        lines.append(name + "\t" + "\t".join(str(v) for v in vector))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_kinship(config, ids, path: Path) -> Path:
    """The stage-10 kinship matrix over ``ids``, from the committed tree."""
    from papipeline.gwas_real import kinship

    kinship.run(
        config, manifest_over(ids), RunMode.TEST,
        tree_path=FIXTURE_TREE, out_path=path,
    )
    return path


def build_distances(config, ids, path: Path) -> Path:
    """The stage-10 *distance* matrix over ``ids`` - pyseer's ``--distances``.

    ``pyseer/input.py::load_structure`` reads it with
    ``pd.read_table(index_col=0)``, which is exactly the shape stage 10 writes,
    so the lineage half of the real path consumes stage 10's own output rather
    than a second matrix derived from the same tree.
    """
    from papipeline.stages import similarity as stage_similarity

    stage_similarity.run(
        config, manifest_over(ids), RunMode.TEST,
        tree_path=FIXTURE_TREE, out_path=path,
    )
    return path


@pytest.fixture(scope="module")
def machine_config(tmp_path_factory, pipeline_root):
    """The laptop overlay, with a writable scratch dir for the vendored helper.

    Step 12a shells out to ``count_patterns.py``, which spills ``sort`` into
    ``runtime.unique_patterns_temp_dir``; the committed overlay points that at
    ``/tmp``, and a test that wrote there would be sharing scratch with every
    other run on the machine. Same fixture shape as
    ``test_gwas_real_binary.py``.
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


class Cohort:
    """The 18-sample binary cohort, with every file the adapter reads."""

    def __init__(self, config, tmp_path: Path, pairs):
        self.pairs = list(pairs)
        self.ids = [s for s, _ in self.pairs]
        self.manifest = manifest_over(self.ids)
        self.phenotype = write_phenotype(tmp_path / "phenotype.tsv", self.pairs)
        self.variants = write_rtab(tmp_path / "features.rtab", self.ids)
        self.similarity = build_kinship(config, self.ids, tmp_path / "kinship.tsv")
        self.distances = build_distances(config, self.ids, tmp_path / "similarity.tsv")
        self.workdir = tmp_path / "pyseer_workdir"


@pytest.fixture
def cohort(machine_config, tmp_path):
    return Cohort(machine_config, tmp_path, binary_pairs())


# --------------------------------------------------------------------------
# The invoker
# --------------------------------------------------------------------------


class FakePyseer:
    """A stand-in for pyseer's stdout/stderr contract.

    Mirrors only what the adapter depends on: ``--output-patterns`` produces a
    patterns file (one pattern per tested variant, read from the ``--pres``
    matrix exactly as pyseer would derive it), and every pass prints pyseer's
    association header to the table it was given. Nothing about the *science*
    is faked - the patterns are the real presence vectors, which is what makes
    the reduction below a real reduction.
    """

    def __init__(self):
        self.calls = []

    def __call__(self, cmd, table: Path, log: Path) -> None:
        cmd = list(cmd)
        self.calls.append(cmd)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("fake pyseer: simulated stderr\n", encoding="utf-8")

        if "--output-patterns" in cmd:
            patterns = Path(cmd[cmd.index("--output-patterns") + 1])
            patterns.parent.mkdir(parents=True, exist_ok=True)
            lines = []
            if "--pres" in cmd:
                rtab = Path(cmd[cmd.index("--pres") + 1])
                for line in rtab.read_text(encoding="utf-8").splitlines()[1:]:
                    if line.strip():
                        lines.append("\t".join(line.split("\t")[1:]))
            if not lines:
                lines = ["01", "10"]
            patterns.write_text("\n".join(lines) + "\n", encoding="utf-8")

        table.parent.mkdir(parents=True, exist_ok=True)
        # The header pyseer would print for *this* command: the extra
        # ``lineage`` column exists only when ``--lineage`` travels
        # (``pyseer/__main__.py``: ``if options.lineage: header += ['lineage']``).
        header = (
            PYSEER_LMM_LINEAGE_HEADER if "--lineage" in cmd else PYSEER_LMM_HEADER
        )
        fields = ["gene__f0", "0.5", "1", "1.0E-04", "1.0", "0.1", "0.01"]
        fields += ["."] * (len(header.split("\t")) - len(fields))
        table.write_text(header + "\n" + "\t".join(fields) + "\n", encoding="utf-8")


def run_adapter(config, cohort, invoke, **kwargs):
    """``adapter.run`` over the cohort's files, with every input the committed
    config implies.

    ``distances_path`` defaults to the cohort's stage-10 distance matrix
    because ``gwas_real.lineage: true`` is the committed setting and pyseer
    exits on ``--lmm --lineage`` without ``--distances`` (__main__.py:220).
    The Cohort builds that matrix for exactly this; a caller passing
    ``distances_path=`` explicitly (the refusal tests do) still wins.
    """
    kwargs.setdefault("distances_path", cohort.distances)
    return adapter.run(
        config,
        cohort.manifest,
        RunMode.TEST,
        phenotype_path=cohort.phenotype,
        variants_path=cohort.variants,
        similarity_path=cohort.similarity,
        workdir=cohort.workdir,
        invoke=invoke,
        **kwargs,
    )


def command_of(result, index: int) -> str:
    return " ".join(result.commands[index])


def requires_pyseer():
    if os.environ.get(OPT_IN) != "1":
        pytest.skip(
            f"set {OPT_IN}=1 to invoke the real pyseer binary; spec.md D8 "
            "defines TEST mode as running no real tools, so this is opt-in"
        )
    if shutil.which("pyseer") is None:
        pytest.skip("pyseer is not on PATH; see docs/environment-arm64.md section 4")


# --------------------------------------------------------------------------
# The command builder
# --------------------------------------------------------------------------


def _build(**overrides):
    base = dict(
        phenotype=Path("phenotype.tsv"),
        variant_source="pres",
        variants=Path("features.rtab"),
        similarity=Path("kinship.tsv"),
        threads=4,
        min_af=1 / 18,
        max_af=1 - 1 / 18,
    )
    base.update(overrides)
    return adapter.build_command(**base)


class TestTheCommandBuilder:
    def test_lmm_and_similarity_are_always_on_the_command(self):
        cmd = _build()
        assert "--lmm" in cmd, "the mixed model is the whole point of this path"
        assert cmd[cmd.index("--similarity") + 1] == "kinship.tsv"

    def test_it_runs_through_the_compatibility_shim(self):
        """Not pyseer directly: 1.1.2 dies against modern scipy without it.

        The shim is the one the stage already uses; two invocations of the same
        pinned tool with two different loaders would be a difference nobody
        chose. See ``stages/gwas.py::_base_command``.
        """
        from papipeline.stages import gwas as stage_gwas

        cmd = _build()
        assert cmd[0] == sys.executable
        assert str(stage_gwas._PYSEER_COMPAT) in cmd

    @pytest.mark.parametrize("source", ["kmers", "pres", "vcf"])
    def test_the_variant_source_is_exactly_one_of_the_three(self, source):
        cmd = _build(variant_source=source)
        present = [s for s in ("--kmers", "--pres", "--vcf") if s in cmd]
        assert present == [f"--{source}"], (
            "pyseer's variant group is mutually exclusive and required "
            "(__main__.py: variant_group), so exactly one flag must travel"
        )

    def test_frequency_bounds_and_threads_are_passed_through(self):
        cmd = _build(threads=4)
        assert cmd[cmd.index("--min-af") + 1] == f"{1 / 18:.6f}"
        assert cmd[cmd.index("--max-af") + 1] == f"{1 - 1 / 18:.6f}"
        assert cmd[cmd.index("--cpu") + 1] == "4"

    def test_lineage_adds_lineage_and_its_distances_matrix(self):
        cmd = _build(lineage=True, distances=Path("similarity.tsv"))
        assert "--lineage" in cmd
        assert cmd[cmd.index("--distances") + 1] == "similarity.tsv"

    def test_the_lineage_table_is_written_into_the_workdir(self):
        """pyseer's default for ``--lineage-file`` is *the current directory*.

        Left at that default, a run started from the repository root drops
        its lineage table there - outside the run's own artefacts, over
        whatever file already has that name - so the path is always stated.
        """
        cmd = _build(
            lineage=True,
            distances=Path("similarity.tsv"),
            lineage_file=Path("workdir/lineage_effects.txt"),
        )
        assert cmd[cmd.index("--lineage-file") + 1] == "workdir/lineage_effects.txt"

    def test_a_lineage_file_without_the_model_term_is_refused(self):
        """A table pyseer only writes under ``--lineage`` is a flag that
        changes nothing when ``--lineage`` is absent."""
        with pytest.raises(DataContractError) as excinfo:
            _build(lineage=False, lineage_file=Path("wherever.txt"))
        message = str(excinfo.value)
        assert "--lineage-file" in message
        assert "under --lineage" in message

    def test_burden_groups_a_vcf_by_its_regions_file(self):
        """``--burden`` takes the regions file - it is not a bare switch.

        pyseer 1.1.2 declares it as ``--burden BURDEN`` with help "VCF regions
        to group variants by for burden testing (requires --vcf)".
        """
        cmd = _build(variant_source="vcf", variants=Path("cohort.vcf"),
                     burden=Path("regions.gff"))
        assert cmd[cmd.index("--burden") + 1] == "regions.gff"
        assert "--vcf" in cmd

    def test_covariates_travel_with_use_covariates_as_separate_tokens(self):
        cmd = _build(covariates=Path("covariates.tsv"), use_covariates=("2", "3q"))
        assert cmd[cmd.index("--covariates") + 1] == "covariates.tsv"
        at = cmd.index("--use-covariates")
        # nargs='*' in pyseer: each whitespace-separated token is one argument,
        # so a single "2 3q" string would be read as one covariate name.
        assert cmd[at + 1:at + 3] == ["2", "3q"]

    def test_the_pass1_gates_are_open_and_explicit(self):
        cmd = _build(output_patterns=Path("patterns.txt"),
                     filter_pvalue=1.0, lrt_pvalue=1.0)
        assert cmd[cmd.index("--output-patterns") + 1] == "patterns.txt"
        assert float(cmd[cmd.index("--filter-pvalue") + 1]) == 1.0
        assert float(cmd[cmd.index("--lrt-pvalue") + 1]) == 1.0

    def test_pass2_carries_the_threshold(self):
        cmd = _build(lrt_pvalue=0.01)
        assert float(cmd[cmd.index("--lrt-pvalue") + 1]) == pytest.approx(0.01)

    # -- the combinations pyseer rejects ------------------------------------

    def test_burden_without_a_vcf_is_refused(self):
        with pytest.raises(DataContractError) as excinfo:
            _build(variant_source="pres", burden=Path("regions.gff"))
        message = str(excinfo.value)
        assert "--burden" in message and "--vcf" in message

    def test_lineage_without_a_distances_matrix_is_refused(self):
        with pytest.raises(DataContractError) as excinfo:
            _build(lineage=True, distances=None)
        assert "--distances" in str(excinfo.value)

    def test_a_distances_matrix_without_lineage_is_refused(self):
        """pyseer rejects ``--lmm --distances`` unless ``--lineage`` is there.

        ``__main__.py:216``: "Must use distance matrix with fixed effects, or
        similarity matrix with random effects / Unless performing a lineage
        analysis with random effects". Passing both silently would hand the
        operator a tool traceback in place of the statement.
        """
        with pytest.raises(DataContractError) as excinfo:
            _build(lineage=False, distances=Path("similarity.tsv"))
        assert "--lineage" in str(excinfo.value)

    def test_an_unknown_variant_source_is_refused(self):
        with pytest.raises(DataContractError) as excinfo:
            _build(variant_source="vcfs")
        message = str(excinfo.value)
        for source in adapter.VARIANT_SOURCES:
            assert source in message

    def test_covariates_without_use_covariates_is_refused(self):
        """Loaded-but-unused covariates are a flag that changes nothing."""
        with pytest.raises(DataContractError) as excinfo:
            _build(covariates=Path("covariates.tsv"), use_covariates=())
        assert "--use-covariates" in str(excinfo.value)

    def test_use_covariates_without_a_covariates_file_is_refused(self):
        with pytest.raises(DataContractError) as excinfo:
            _build(covariates=None, use_covariates=("2",))
        assert "--covariates" in str(excinfo.value)

    def test_missing_threads_are_refused_naming_the_key(self):
        """Threads are a machine fact; an inferred one is a guess about hardware."""
        with pytest.raises(DataContractError) as excinfo:
            _build(threads=None)
        assert "runtime.threads" in str(excinfo.value)


# --------------------------------------------------------------------------
# The settings section
# --------------------------------------------------------------------------


class _RawConfig:
    def __init__(self, raw):
        self.raw = raw


class TestTheSettingsSection:
    def test_the_committed_science_file_declares_it(self, config):
        settings = AdapterSettings.from_config(config)
        assert settings.variant_source == "pres", (
            "the real path's default variant source is read from "
            "config gwas_real.variant_source, not chosen in code"
        )
        assert settings.lineage is True
        assert settings.burden is None
        assert settings.use_covariates == ()

    def test_an_unknown_variant_source_in_config_is_refused(self, config):
        raw = dict(config.raw)
        raw["gwas_real"] = {"variant_source": "aligned"}
        with pytest.raises(DataContractError) as excinfo:
            AdapterSettings.from_config(_RawConfig(raw))
        message = str(excinfo.value)
        for source in adapter.VARIANT_SOURCES:
            assert source in message

    def test_burden_in_config_without_a_vcf_source_is_refused(self, config):
        raw = dict(config.raw)
        raw["gwas_real"] = {"variant_source": "pres", "burden": "regions.gff"}
        with pytest.raises(DataContractError) as excinfo:
            AdapterSettings.from_config(_RawConfig(raw))
        message = str(excinfo.value)
        assert "--burden" in message and "--vcf" in message

    def test_a_section_that_is_absent_falls_back_to_the_declared_defaults(
        self, config
    ):
        raw = dict(config.raw)
        raw.pop("gwas_real", None)
        settings = AdapterSettings.from_config(_RawConfig(raw))
        assert settings.variant_source == "pres"
        assert settings.burden is None


# --------------------------------------------------------------------------
# The run: two passes over the real protocol, on synthetic fixtures
# --------------------------------------------------------------------------


class TestTheTwoPassRun:
    def test_it_runs_both_passes_and_reduces_the_patterns(self, machine_config, cohort):
        fake = FakePyseer()
        result = run_adapter(machine_config, cohort, fake)

        assert len(fake.calls) == 2, (
            "step 12a is two passes by construction: the threshold needs "
            "pyseer's own patterns file, which only pyseer can write"
        )
        pass1, pass2 = fake.calls

        # Pass 1: collect patterns with both p-value gates explicitly open.
        assert "--output-patterns" in pass1
        assert float(pass1[pass1.index("--lrt-pvalue") + 1]) == 1.0
        assert float(pass1[pass1.index("--filter-pvalue") + 1]) == 1.0

        # Pass 2: the association table, at the threshold pass 1 derived.
        assert "--output-patterns" not in pass2
        threshold = float(pass2[pass2.index("--lrt-pvalue") + 1])
        assert threshold == pytest.approx(result.reduction.threshold)

        # The reduction is real: 8 features, 5 distinct patterns.
        assert result.reduction.n_variants_tested == 8
        assert result.reduction.n_unique_patterns == 5
        assert result.reduction.threshold == pytest.approx(0.05 / 5)

        # Both passes carry the flags this path exists for.
        for cmd in (pass1, pass2):
            assert "--lmm" in cmd
            assert "--pres" in cmd
            assert cmd[cmd.index("--similarity") + 1] == str(cohort.similarity)
            assert "--lineage" in cmd
            assert cmd[cmd.index("--distances") + 1] == str(cohort.distances)
            assert cmd[cmd.index("--cpu") + 1] == str(machine_config.threads)
            # The lineage table goes into the run's workdir, never into
            # whatever directory pyseer happened to start in.
            lineage_file = Path(cmd[cmd.index("--lineage-file") + 1])
            assert lineage_file.parent == cohort.workdir
            assert lineage_file.name == "lineage_effects.txt"
            assert "--continuous" not in cmd, (
                "the outcome is binary; --continuous would force the linear "
                "model pyseer otherwise avoids for 0/1 phenotypes"
            )

        # The artefacts a reader needs afterwards.
        assert result.association_table.is_file()
        header = result.association_table.read_text().splitlines()[0]
        assert header.split("\t")[:2] == ["variant", "af"]
        assert header.split("\t")[-1] == "notes"
        assert result.patterns_file.is_file()
        assert len(result.patterns_file.read_text().splitlines()) == 8

    def test_the_threshold_comes_from_the_existing_reduction(
        self, machine_config, cohort, monkeypatch
    ):
        """Reuse, asserted by spying on the call rather than by reading prose.

        A second implementation of ``alpha / n`` in this package would be able
        to pass every other test here while disagreeing with step 12a the day
        the correction changes - so the call itself is the contract, and the
        value handed to pass 2 must be the spied return value's.
        """
        from papipeline.stages import gwas as stage_gwas

        seen = {}

        class StubReduction:
            """A double shaped like ``stage_gwas.UniquePatternReduction``.

            Only ``threshold`` is under test, but the adapter reads ``alpha``
            and ``n_unique_patterns`` too - for the log line and the
            provenance - so a stub without them fails inside production code
            rather than at the assertion this test owns.
            """

            alpha = 0.05
            n_unique_patterns = 5
            n_variants_tested = 8
            threshold = 4.2e-3

        def spy(config, patterns_file, workdir):
            seen["args"] = (patterns_file, workdir)
            return StubReduction()

        monkeypatch.setattr(stage_gwas, "reduce_unique_patterns", spy)
        fake = FakePyseer()
        result = run_adapter(machine_config, cohort, fake)

        assert seen["args"] == (result.patterns_file, cohort.workdir)
        pass2 = fake.calls[1]
        assert float(pass2[pass2.index("--lrt-pvalue") + 1]) == pytest.approx(4.2e-3)

    def test_the_provenance_names_the_tool_invocation_and_the_threshold(
        self, machine_config, cohort
    ):
        from papipeline.io.tsv import read_tsv

        result = run_adapter(machine_config, cohort, FakePyseer())
        assert result.provenance.is_file()
        record = {r["key"]: r["value"] for r in read_tsv(result.provenance)}
        # The house convention for labelling a run in an artefact is the
        # RunMode value itself (reporting.py: ``context.mode.value``; the
        # regulator report test asserts ``"REAL"``), not its lowercase name.
        assert record["mode"] == "TEST"
        assert record["variant_source"] == "pres"
        assert record["lineage"] == "true"
        assert record["n_samples"] == str(len(cohort.ids))
        assert record["kinship_source"] == str(cohort.similarity)
        assert Path(record["lineage_effects"]).parent == cohort.workdir
        assert "--similarity" in record["command_1"]
        assert "--lmm" in record["command_1"]
        assert float(record["lrt_pvalue_threshold"]) == pytest.approx(
            result.reduction.threshold
        )
        # Pass 1 gated open at 1 so every pattern reaches the file; pass 2
        # gated at exactly the threshold the reduction returned. The values
        # are floats - build_command formats them with repr(float(...)) - so
        # parsing them as ints would raise rather than assert.
        assert float(record["command_1"].split("--lrt-pvalue")[1].split()[0]) == 1.0
        assert float(record["command_2"].split("--lrt-pvalue")[1].split()[0]) == (
            pytest.approx(result.reduction.threshold)
        )

    def test_non_lineage_settings_omit_both_lineage_flags(
        self, machine_config, cohort
    ):
        settings = AdapterSettings(variant_source="pres", lineage=False)
        fake = FakePyseer()
        run_adapter(machine_config, cohort, fake, settings=settings)
        for cmd in fake.calls:
            assert "--lineage" not in cmd
            assert "--lineage-file" not in cmd
            assert "--distances" not in cmd

    def test_covariates_reach_both_passes(self, machine_config, cohort, tmp_path):
        covariates = tmp_path / "covariates.tsv"
        covariates.write_text(
            "\n".join(f"{s}\t{i * 0.5}" for i, s in enumerate(cohort.ids)) + "\n",
            encoding="utf-8",
        )
        settings = AdapterSettings(
            variant_source="pres", lineage=True,
            covariates=str(covariates), use_covariates=("2q",),
        )
        fake = FakePyseer()
        run_adapter(machine_config, cohort, fake, settings=settings)
        for cmd in fake.calls:
            assert cmd[cmd.index("--covariates") + 1] == str(covariates)
            at = cmd.index("--use-covariates")
            assert cmd[at + 1] == "2q"
            assert cmd[at + 2] in {"--distances", "--lineage"}


class TestWhatTheRunRefuses:
    def test_test_mode_without_an_invoker_is_refused_unless_opted_in(
        self, machine_config, cohort, monkeypatch
    ):
        monkeypatch.delenv(OPT_IN, raising=False)
        with pytest.raises(NotImplementedError) as excinfo:
            run_adapter(machine_config, cohort, invoke=None)
        message = str(excinfo.value)
        assert "invoke" in message
        assert OPT_IN in message, (
            "the refusal must name the opt-in, or a reader has no way to "
            "cross the D8 boundary deliberately"
        )

    def test_stub_is_refused(self, machine_config, cohort):
        with pytest.raises(NotImplementedError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.STUB,
                phenotype_path=cohort.phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
            )
        assert "STUB" in str(excinfo.value)

    def test_real_refuses_while_the_gate_is_closed(
        self, pipeline_root, cohort, monkeypatch, tmp_path
    ):
        from papipeline.config import loader as loader_module

        monkeypatch.delenv(loader_module.ALLOW_REAL_MODE_ENV, raising=False)
        config = loader_module.load_config(pipeline_root / "config" / "science.yaml")
        assert config.runtime.get("allow_real_mode") is False, (
            "the committed overlays keep REAL shut; this test is what says so"
        )
        with pytest.raises(NotImplementedError) as excinfo:
            adapter.run(
                config, cohort.manifest, RunMode.REAL,
                phenotype_path=cohort.phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=tmp_path / "real_workdir",
                invoke=FakePyseer(),
            )
        message = str(excinfo.value)
        assert "allow_real_mode" in message
        assert "runtime.allow_real_mode" in message

    def test_real_past_the_gate_needs_its_inputs_first(
        self, pipeline_root, cohort, monkeypatch, tmp_path
    ):
        """Open the gate through the real mechanism, then hit the next fact:
        the inputs are files, and an absent file is an absent input."""
        from papipeline.config import loader as loader_module

        monkeypatch.setenv(loader_module.ALLOW_REAL_MODE_ENV, "1")
        monkeypatch.setenv(loader_module.RESULTS_ROOT_ENV, str(tmp_path / "results"))
        config = loader_module.load_config(pipeline_root / "config" / "science.yaml")
        assert config.runtime.get("allow_real_mode") is True
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                config, cohort.manifest, RunMode.REAL,
                phenotype_path=tmp_path / "absent_phenotype.tsv",
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=tmp_path / "real_workdir",
                invoke=FakePyseer(),
            )
        assert "phenotype" in str(excinfo.value).lower()

    def test_threads_come_from_the_machine_overlay(
        self, pipeline_root, cohort, monkeypatch
    ):
        """An overlay-less config has no threads, and an inferred count would be
        a guess about hardware nobody described."""
        from papipeline.config import loader as loader_module

        config = loader_module.load_config(
            pipeline_root / "config" / "science.yaml", machine=None
        )
        assert config.threads is None
        with pytest.raises(DataContractError) as excinfo:
            run_adapter(config, cohort, FakePyseer())
        assert "runtime.threads" in str(excinfo.value)


class TestWhatTheInputsRefuse:
    def test_a_phenotype_sample_that_is_not_in_the_manifest_is_refused(
        self, machine_config, cohort, tmp_path
    ):
        pairs = cohort.pairs + [("TEST_PA_999", 1)]
        phenotype = write_phenotype(tmp_path / "extra.tsv", pairs)
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.TEST,
                phenotype_path=phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
            )
        message = str(excinfo.value)
        assert "TEST_PA_999" in message
        assert "manifest" in message.lower()

    def test_a_manifest_sample_missing_from_the_phenotype_is_refused(
        self, machine_config, cohort, tmp_path
    ):
        pairs = [p for p in cohort.pairs if p[0] != "TEST_PA_001"]
        phenotype = write_phenotype(tmp_path / "short.tsv", pairs)
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.TEST,
                phenotype_path=phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
            )
        assert "TEST_PA_001" in str(excinfo.value)

    def test_a_non_binary_phenotype_is_refused_before_pyseer_can_go_continuous(
        self, machine_config, cohort, tmp_path
    ):
        """pyseer switches to a linear model on any value outside 0/1 - at exit
        0, with a line of stderr nobody reads (``__main__.py:242-249``)."""
        phenotype = tmp_path / "letters.tsv"
        phenotype.write_text(
            "sample\tphenotype\n"
            + "".join(f"{s}\t{'R' if v else 'S'}\n" for s, v in cohort.pairs),
            encoding="utf-8",
        )
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.TEST,
                phenotype_path=phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
            )
        message = str(excinfo.value)
        assert "continuous" in message.lower()
        assert "0" in message and "1" in message
        assert "gwas.outcome" in message

    def test_a_group_below_the_configured_floor_is_refused(
        self, machine_config, tmp_path
    ):
        pairs = [p for p in binary_pairs() if p[1] == 0][:10]
        pairs += [p for p in binary_pairs() if p[1] == 1][:2]  # 10 S, 2 R
        ids = [s for s, _ in pairs]
        phenotype = write_phenotype(tmp_path / "small.tsv", pairs)
        similarity = build_kinship(machine_config, ids, tmp_path / "small_k.tsv")
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, manifest_over(ids), RunMode.TEST,
                phenotype_path=phenotype,
                variants_path=write_rtab(tmp_path / "small.rtab", ids),
                similarity_path=similarity,
                workdir=tmp_path / "small_workdir",
                invoke=FakePyseer(),
            )
        assert "gwas.min_samples_per_group" in str(excinfo.value)

    def test_a_similarity_matrix_over_a_different_cohort_is_refused(
        self, machine_config, cohort, tmp_path
    ):
        """pyseer intersects phenotype and similarity **silently**.

        ``pyseer/lmm.py::initialise_lmm`` does
        ``p.index.intersection(K.index)`` and then ``K.loc[p.index, p.index]``,
        so a matrix describing 20 samples against a phenotype of 18 does not
        error - it quietly analyses 18 and reports a run over "the cohort".
        The check has to happen before the tool sees either file.
        """
        full = binary_pairs() + [("TEST_PA_009", 1)]  # the excluded I isolate
        wrong = build_kinship(
            machine_config, [s for s, _ in full], tmp_path / "kinship20.tsv"
        )
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.TEST,
                phenotype_path=cohort.phenotype,
                variants_path=cohort.variants,
                similarity_path=wrong,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
            )
        message = str(excinfo.value)
        assert "TEST_PA_009" in message
        assert "pyseer" in message.lower()

    def test_lineage_settings_without_a_distances_matrix_are_refused(
        self, machine_config, cohort
    ):
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, cohort.manifest, RunMode.TEST,
                phenotype_path=cohort.phenotype,
                variants_path=cohort.variants,
                similarity_path=cohort.similarity,
                workdir=cohort.workdir,
                invoke=FakePyseer(),
                distances_path=None,
            )
        assert "--distances" in str(excinfo.value)

    def test_a_covariate_file_missing_a_cohort_sample_is_refused(
        self, machine_config, cohort, tmp_path
    ):
        covariates = tmp_path / "covariates.tsv"
        covariates.write_text(
            "\n".join(
                f"{s}\t1.0" for s in cohort.ids if s != "TEST_PA_001"
            ) + "\n",
            encoding="utf-8",
        )
        settings = AdapterSettings(
            variant_source="pres", lineage=False,
            covariates=str(covariates), use_covariates=("2q",),
        )
        with pytest.raises(DataContractError) as excinfo:
            run_adapter(machine_config, cohort, FakePyseer(), settings=settings)
        message = str(excinfo.value)
        assert "TEST_PA_001" in message
        assert "covariate" in message.lower()

    def test_an_empty_cohort_is_refused(self, machine_config, tmp_path):
        """The frequency floor is 1/n; with n=0 there is no floor to derive."""
        with pytest.raises(DataContractError) as excinfo:
            adapter.run(
                machine_config, manifest_over([]), RunMode.TEST,
                phenotype_path=tmp_path / "none.tsv",
                variants_path=tmp_path / "none.rtab",
                similarity_path=tmp_path / "none_k.tsv",
                workdir=tmp_path / "none_workdir",
                invoke=FakePyseer(),
            )
        assert "cohort" in str(excinfo.value).lower()


class TestTheLauncher:
    def test_a_failing_pyseer_is_reported_with_its_exit_code_and_log(self, tmp_path):
        """Not a bare OSError from a subprocess whose stderr went to a log.

        The stage's own `_invoke` sets this shape; the adapter's launcher has
        to say which pass failed and where pyseer's own words went.
        """
        table = tmp_path / "out.tsv"
        log = tmp_path / "out.log"
        with pytest.raises(RuntimeError) as excinfo:
            adapter._launch(
                [sys.executable, "-c", "import sys; sys.exit(3)"], table, log
            )
        message = str(excinfo.value)
        assert "exit code 3" in message
        assert str(log) in message


class TestTheShimReachesTheWorkers:
    """The compat shim must arrive in pyseer's Pool workers, not only its parent.

    pyseer fits every variant in a ``Pool`` worker when ``--cpu > 1``, and
    ``pyseer/model.py::fit_lineage_effect`` calls ``smf.Logit`` there when
    ``--lineage`` is on. ``scripts/pyseer_compat.py`` applies that alias
    inside ``main()``, which a *spawned* worker never calls - macOS's default
    start method is spawn - so without delivery through the environment the
    parent fits the null model and the worker dies with
    ``AttributeError: module 'statsmodels.formula.api' has no attribute
    'Logit'``. That is what this pair of tests pins down: the env the adapter
    builds, and the effect it has on a process it does not control.
    """

    def test_the_shim_directory_leads_pythonpath_and_existing_entries_survive(
        self, monkeypatch
    ):
        monkeypatch.setenv("PYTHONPATH", "/somewhere/else")
        env = adapter._pyseer_env()
        entries = env["PYTHONPATH"].split(os.pathsep)
        assert entries[0] == str(adapter.SHIM_ENV_DIR), (
            "the shim directory must lead, or a same-named module further "
            "along PYTHONPATH could shadow sitecustomize"
        )
        assert "/somewhere/else" in entries, (
            "PYTHONPATH is prepended to, never replaced: dropping a caller's "
            "own entries to deliver a shim trades one silent breakage for another"
        )
        # Nothing else about the environment is touched.
        assert set(env) == set(os.environ) | {"PYTHONPATH"}

    def test_a_spawned_pool_worker_applies_the_shim(self, tmp_path):
        """The defect itself, reproduced through multiprocessing.spawn.

        The probe script is deliberately *not* the compat script: a worker
        re-runs its ``__main__`` top level and must still get the shim, which
        is exactly the case ``apply_shims()``-inside-``main()`` misses.
        """
        pytest.importorskip("statsmodels")
        script = tmp_path / "worker_probe.py"
        script.write_text(
            "import multiprocessing as mp\n"
            "\n"
            "def probe():\n"
            "    import statsmodels.formula.api as smf\n"
            "    return hasattr(smf, 'Logit')\n"
            "\n"
            "if __name__ == '__main__':\n"
            "    with mp.Pool(1) as pool:\n"
            "        print(pool.apply(probe))\n",
            encoding="utf-8",
        )
        import subprocess

        proc = subprocess.run(
            [sys.executable, str(script)],
            capture_output=True,
            text=True,
            env=dict(adapter._pyseer_env()),
        )
        assert proc.stdout.strip() == "True", (
            "a spawned worker did not get scripts/pyseer_compat.py's alias; "
            "pyseer's --lineage fit would die in that worker with "
            "AttributeError on smf.Logit. stdout=%r stderr=%r"
            % (proc.stdout, proc.stderr)
        )


# --------------------------------------------------------------------------
# The real binary, opt-in, on synthetic fixtures
# --------------------------------------------------------------------------


class TestTheRealBinaryOnSyntheticFixtures:
    """R12's rule: real code paths, synthetic inputs. Opt-in, like its sibling."""

    def test_the_adapter_runs_the_real_pyseer(self, machine_config, cohort, monkeypatch):
        requires_pyseer()
        result = run_adapter(machine_config, cohort, invoke=None)
        assert result.pyseer_version, "the version was never resolved"
        assert result.association_table.is_file()
        lines = [
            line
            for line in result.association_table.read_text().splitlines()
            if line and not line.startswith("#")
        ]
        assert lines, "pyseer wrote no header; nothing was tested"
        # The committed config asks for lineage, so pyseer's header carries
        # its extra column; asserted rather than skipped past, because a
        # --lineage run that dies in a Pool worker is exactly the defect the
        # shim-delivery tests above pin down.
        assert result.settings.lineage is True
        assert lines[0].split("\t") == PYSEER_LMM_LINEAGE_HEADER.split("\t")
        assert result.reduction.n_unique_patterns >= 1
        assert result.reduction.threshold == pytest.approx(
            0.05 / result.reduction.n_unique_patterns
        )
        # Both passes really launched, and pass 2 gated on the derived value.
        assert len(result.commands) == 2
        assert "--lmm" in result.commands[0]
        assert float(
            result.commands[1][result.commands[1].index("--lrt-pvalue") + 1]
        ) == pytest.approx(result.reduction.threshold)
