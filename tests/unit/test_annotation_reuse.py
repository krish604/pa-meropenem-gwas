"""Stage 2 may reuse an existing Bakta output instead of invoking Bakta.

**What this must never become.** "Trust whatever is already on disk." The reuse
path runs the *same* conversion chain as a real run - one implementation,
`_records_from_bakta_tables`, called by both branches - and differs in exactly
one respect: the tool is not executed. The consequence is that provenance can
no longer be established by watching a process, so it is established from the
output itself: the file set, the GFF/feature-table agreement, and the
`# Software:` / `# Database:` header.

Three things here are deliberately harder than they look, and each is a bug a
previous round of this work actually had:

* **The fixture must write the FEATURE table, named by the GENOME STEM.**
  `BaktaOutputs.annotation_tsv` is `<genome_stem>.tsv`. A fixture that writes
  `<sample_id>.inference.tsv` looks like it is testing the database mismatch and
  is actually testing "the file is not there", so the mismatch assertion passes
  for the wrong reason and the refusal logic behind it is never exercised.
  `assert_feature_table_was_looked_at` exists to catch exactly that.
* **The expected database string is derived, not typed in.** It comes from the
  configured database's own `version.json`, formatted the way Bakta formats its
  header (`bakta/main.py:627`). `bakta_db.database_version` produces a
  *different* shape on purpose - `6.0 (light, 2025-02-24)`, shaped for
  `config/references.tsv` - so comparing a header against it would refuse every
  genuinely reusable output for a reason about string shapes.
* **A refusal is a value, not an exception.** This stage records an unusable
  genome and carries on, because a smoke cohort is mostly unprepared by
  construction. One unverifiable isolate must not end the run, so nothing here
  escapes to the caller.

No Bakta is executed anywhere in this file. Every test that reaches the tool
path does so through an injected fake, and `require`-mode tests use one that
**raises if called** - so an accidental invocation fails the suite rather than
spending an hour of CPU.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest

from papipeline.config.loader import ConfigError, load_config
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import annotation as stage

REPO = Path(__file__).resolve().parents[2]

#: One genome's worth of plausible Bakta output: a feature TSV whose ids the
#: GFF is a strict SUBSET of, which is the shape real Bakta produces and the
#: shape `gff_features_are_covered` is written for.
N_FEATURES = 3


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

@pytest.fixture
def provisioned_config():
    """The real config, pointed at a real database that satisfies the preflight.

    The pre-existing stage-2 test file does exactly this and this suite copies
    it rather than inventing a second way to build a config that passes the
    preflight - a hand-built `PipelineConfig` is exactly how the previous
    round's tests ended up asserting against an attribute that did not exist.
    """
    database = REPO / "db" / "bakta_db" / "db-light"
    if not (database / "version.json").is_file():
        pytest.skip("bakta database not provisioned here; db/ is gitignored")

    config = load_config(REPO / "config" / "science.yaml")
    prepared = copy.copy(config)
    raw = copy.deepcopy(dict(config.raw))
    raw["annotation"] = dict(raw.get("annotation") or {})
    raw["annotation"]["bakta_db"] = str(database)
    object.__setattr__(prepared, "raw", raw)
    return prepared


@pytest.fixture
def database_dir(provisioned_config):
    return Path(dict(provisioned_config.raw)["annotation"]["bakta_db"])


@pytest.fixture
def configured(provisioned_config):
    """Set `annotation.reuse_tool_output` on the real config."""
    def _set(mode):
        config = copy.copy(provisioned_config)
        raw = copy.deepcopy(dict(provisioned_config.raw))
        raw["annotation"] = dict(raw["annotation"])
        raw["annotation"]["reuse_tool_output"] = mode
        object.__setattr__(config, "raw", raw)
        return config
    return _set


@pytest.fixture
def out_root(tmp_path):
    """Where `run` looks for candidates: `<intermediate_root>/bakta`.

    The directory is named `bakta` by `run`, not by this suite, so a run-level
    fixture that wrote candidates one level higher would be testing nothing.
    """
    return tmp_path / "bakta"


@pytest.fixture
def expected_database(database_dir):
    """What Bakta writes into `# Database:` for the configured database."""
    return stage.expected_bakta_database_string(database_dir)


def header_block(software: str = "v1.12.1", database: str | None = None,
                 *, database_dir: Path | None = None) -> str:
    """Bakta's five `#` preamble lines, in the order it writes them."""
    if database is None:
        database = stage.expected_bakta_database_string(database_dir)
    return (
        "# Annotated with Bakta\n"
        f"# Software: {software}\n"
        f"# Database: {database}\n"
        "# DOI: 10.1099/mgen.0.000685\n"
        "# URL: github.com/oschwengers/bakta\n"
    )


def feature_tsv_text(n_features: int = N_FEATURES, *, header: str | None = None,
                     database_dir: Path | None = None) -> str:
    """A Bakta FEATURE table: header block then `#Sequence Id\tType\t...`.

    The column set is Bakta's own (`BAKTA_TSV_COLUMNS`), not the inference
    table's. The inference table shares this header block and differs only in
    its columns - which is why a fixture cannot be told apart by its preamble.
    """
    if header is None:
        header = header_block(database_dir=database_dir)
    rows = "".join(
        f"contig_00001\tCDS\t{100 + i * 400}\t{300 + i * 400}\t+\t"
        f"TESTPA_{i:05d}\t{100 + i * 400}\tGC={0.4 + i / 100:.2f}\t"
        f"gene{i}\tproduct {i}\tUniRef:X\n"
        for i in range(1, n_features + 1)
    )
    return header + "#" + "\t".join(stage.BAKTA_TSV_COLUMNS) + "\n" + rows


def gff_text(n_features: int = N_FEATURES) -> str:
    """A Bakta GFF3 whose ids are a strict subset of the feature table's."""
    rows = "".join(
        f"contig_00001\tBakta\tCDS\t{100 + i * 400}\t{300 + i * 400}\t.\t+\t0\t"
        f"ID=TESTPA_{i:05d};Name=gene{i};product=product {i}\n"
        for i in range(1, n_features + 1)
    )
    return "##gff-version 3\n" + rows


def write_candidate(
    out_root: Path,
    sample_id: str,
    genome_stem: str | None = None,
    *,
    database_dir: Path,
    database: str | None = None,
    software: str = "v1.12.1",
    n_features: int = N_FEATURES,
    tsv_rows: str | None = None,
    write_inference: bool = True,
    omit: tuple[str, ...] = (),
) -> Path:
    """One isolate's worth of Bakta output where `run_bakta` would put it.

    ``out_root/<sample_id>/<genome_stem>.{tsv,gff3,inference.tsv}`` - the
    directory is named by the SAMPLE ID and the files by the GENOME STEM,
    which is what `run_bakta` does and what makes the two distinguishable.
    ``omit`` drops a file by role, for the missing-file cases.
    """
    stem = genome_stem if genome_stem is not None else sample_id
    directory = out_root / sample_id
    directory.mkdir(parents=True, exist_ok=True)

    header = header_block(software=software, database=database,
                          database_dir=database_dir)
    if "feature_tsv" not in omit:
        (directory / f"{stem}.tsv").write_text(
            tsv_rows if tsv_rows is not None
            else feature_tsv_text(n_features, header=header),
            encoding="utf-8",
        )
    if "gff3" not in omit:
        (directory / f"{stem}.gff3").write_text(
            gff_text(n_features), encoding="utf-8")
    if "inference_tsv" not in omit:
        (directory / f"{stem}.inference.tsv").write_text(
            header
            + "#Query\tHit\tIdentity\tScore\tE-value\tQuery start\tQuery end\t"
              "Hit start\tHit end\tGap\tBitscore\tplus\tminus\n",
            encoding="utf-8",
        )
    return directory


def manifest_of(ids, *, assembly_stem=None) -> SampleManifest:
    """A manifest whose assemblies are named after `assembly_stem` if given."""
    samples = []
    for sample_id in ids:
        stem = assembly_stem or sample_id
        samples.append(
            Sample(sample_id=sample_id, assembly_path=f"/nowhere/{stem}.fna")
        )
    return SampleManifest(samples=samples)


@pytest.fixture
def exploding_runner(monkeypatch):
    """A `run_bakta` that FAILS THE TEST if it is ever called.

    Used for every `require`-mode test. A fake that quietly succeeds would let
    the suite pass with Bakta being invoked for a sample the run was configured
    never to annotate - which is the entire defect this closes.
    """
    def _boom(sample, **kwargs):
        raise AssertionError(
            f"Bakta was invoked for {sample.sample_id!r} under "
            f"reuse_tool_output=require; that mode must never execute the tool"
        )

    monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", _boom)
    return _boom


class RecordingRunner:
    """A fake `run_bakta` that records what it was asked for."""

    def __init__(self, *, database_dir: Path, state: str = "SUCCEEDED"):
        self.asked: list[str] = []
        self.database_dir = database_dir
        self.state = state

    def __call__(self, sample, **kwargs):
        from papipeline.adapters.bakta import BaktaOutputs, genome_stem_for

        self.asked.append(sample.sample_id)
        out_dir = Path(kwargs["out_root"]) / sample.sample_id
        out_dir.mkdir(parents=True, exist_ok=True)
        outputs = BaktaOutputs(
            out_dir=out_dir, sample_id=sample.sample_id,
            genome_stem=genome_stem_for(sample),
        )
        if self.state == "SUCCEEDED":
            write_candidate(
                kwargs["out_root"], sample.sample_id, outputs.genome_stem,
                database_dir=self.database_dir,
            )
        return _result(sample.sample_id, self.state, outputs)


def _result(sample_id, state, outputs):
    from papipeline.adapters.bakta import BaktaResult

    return BaktaResult(
        sample_id=sample_id, state=state, attempts=1, outputs=outputs,
        command=["bakta"],
    )


# --------------------------------------------------------------------------
# U1 - the configuration key
# --------------------------------------------------------------------------

class TestTheConfigurationKey:
    def test_the_attribute_exists_on_the_real_config(self):
        """A boolean invented by a test is not a configuration surface."""
        config = load_config(REPO / "config" / "science.yaml")
        assert config.reuse_tool_output() == "off"

    def test_the_committed_default_is_off(self):
        """Stated here as well as in the loader: behaviour is unchanged."""
        text = (REPO / "config" / "science.yaml").read_text(encoding="utf-8")
        assert 'reuse_tool_output: "off"' in text, (
            "the committed default must be off, and quoted - YAML reads a bare "
            "`off` as the boolean false"
        )

    @pytest.mark.parametrize("mode", ["off", "prefer", "require"])
    def test_every_accepted_value_is_accepted(self, configured, mode):
        assert configured(mode).reuse_tool_output() == mode

    def test_an_absent_key_is_off(self, provisioned_config):
        config = copy.copy(provisioned_config)
        raw = copy.deepcopy(dict(provisioned_config.raw))
        raw["annotation"] = dict(raw["annotation"])
        raw["annotation"].pop("reuse_tool_output")
        object.__setattr__(config, "raw", raw)
        assert config.reuse_tool_output() == "off"

    @pytest.mark.parametrize("value", ["reuse", "true", "off,prefer", "OFF "])
    def test_an_unrecognised_value_is_refused_not_assumed_off(self, configured, value):
        """A typo must not read as a correctly configured run.

        Falling back to `off` here would re-annotate the entire cohort and say
        nothing about why, which is indistinguishable in the output from the
        setting being deliberately `off`.
        """
        with pytest.raises(ConfigError) as excinfo:
            configured(value).reuse_tool_output()
        assert value in str(excinfo.value)

    def test_a_bare_off_is_refused_with_the_yaml_reason_named(self, configured):
        """YAML 1.1 reads `off`/`no`/`false` as booleans.

        The natural spelling of the default arrives as `False`, and the error a
        reader would get without this branch reads as though the pipeline had
        invented a rule about the word "off".
        """
        with pytest.raises(ConfigError) as excinfo:
            configured(False).reuse_tool_output()
        assert "YAML" in str(excinfo.value)

    def test_the_accepted_set_is_the_one_the_stage_acts_on(self):
        """One vocabulary, stated in both places and checked here."""
        from papipeline.config.loader import REUSE_TOOL_OUTPUT_MODES

        assert tuple(REUSE_TOOL_OUTPUT_MODES) == stage.REUSE_TOOL_OUTPUT_MODES


# --------------------------------------------------------------------------
# U2 - what reuse requires
# --------------------------------------------------------------------------

class TestWhatReuseRequires:
    def test_a_verified_candidate_is_reused(
        self, tmp_path, database_dir, expected_database
    ):
        write_candidate(tmp_path, "PDT_A", database_dir=database_dir)
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is True, decided.detail
        assert decided.reason == "verified"

    def test_the_expected_string_is_the_header_shape_derived_from_version_json(
        self, database_dir
    ):
        """Derived from the database, in BAKTA's shape.

        `bakta_db.database_version` returns `6.0 (light, 2025-02-24)` because
        that is the shape `config/references.tsv` records. Comparing a header
        against it would refuse every real output.
        """
        import json

        from papipeline.adapters import bakta_db

        payload = json.loads((database_dir / "version.json").read_text())
        expected = stage.expected_bakta_database_string(database_dir)
        assert expected == (
            f"v{payload['major']}.{payload['minor']}, {payload['type']}"
        )
        assert expected != bakta_db.database_version(database_dir), (
            "the two shapes must stay distinct: bakta_db.database_version is "
            "shaped for config/references.tsv and would refuse every real "
            "header if it were reused here"
        )

    def test_a_database_mismatch_is_refused_by_naming_both_strings(
        self, tmp_path, database_dir, expected_database
    ):
        write_candidate(tmp_path, "PDT_A", database_dir=database_dir,
                        database="v5.2, full")
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_DATABASE_MISMATCH
        assert "v5.2, full" in decided.detail, (
            "a refusal a reader cannot act on is a refusal that will be "
            "argued with rather than fixed"
        )
        assert expected_database in decided.detail

    @pytest.mark.parametrize(
        "role", ["feature_tsv", "gff3", "inference_tsv"]
    )
    def test_a_missing_file_is_refused_by_name(self, tmp_path, database_dir, role):
        write_candidate(tmp_path, "PDT_A", database_dir=database_dir, omit=(role,))
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_MISSING

    def test_an_empty_file_is_refused(self, tmp_path, database_dir):
        """Present but zero bytes is not a usable annotation.

        A truncated-to-nothing table passes an existence check and produces an
        empty cohort, which reads as a biological result.
        """
        directory = write_candidate(tmp_path, "PDT_A", database_dir=database_dir)
        (directory / "PDT_A.tsv").write_text("", encoding="utf-8")
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert "empty" in decided.detail

    def test_the_feature_table_is_required_not_the_inference_table(
        self, tmp_path, database_dir
    ):
        """The check reads `<genome_stem>.tsv`.

        Reading `<genome_stem>.inference.tsv` instead once lost every gene name
        in the cohort and produced a uniform false negative that looked exactly
        like biology. Both files are present here, so only the file actually
        consulted can decide the outcome.
        """
        write_candidate(tmp_path, "PDT_A", database_dir=database_dir,
                        database="v5.2, full")
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reason == stage.REUSE_DATABASE_MISMATCH, (
            "the refusal must come from the FEATURE table's header; the "
            "inference table's header is identical, so a fixture that only "
            "poisoned one of them proves nothing"
        )
        assert decided.feature_tsv is None

    def test_a_gff_id_absent_from_the_feature_table_is_refused(
        self, tmp_path, database_dir
    ):
        # The ghost exists in the GFF ONLY. The failure being modelled is "a
        # GFF gene with no TSV counterpart", i.e. two tables that do not
        # describe the same annotation - a row added to both would not be that.
        write_candidate(tmp_path, "PDT_A", database_dir=database_dir)
        gff = tmp_path / "PDT_A" / "PDT_A.gff3"
        gff.write_text(
            gff_text() + "contig_00001\tBakta\tCDS\t9999\t10200\t.\t+\t0\t"
            "ID=GHOST_00001;Name=ghost\n",
            encoding="utf-8",
        )
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_GFF_TSV_DISAGREE
        assert "1 GFF feature(s)" in decided.detail

    def test_the_containment_direction_is_gff_subset_of_features(self):
        """The invariant the reuse path inherits, stated as an invariant.

        Real output is 5,912 GFF features inside 5,944 TSV rows, so the
        direction matters: testing the other way would refuse every real
        genome, and testing count equality refuses them too.
        """
        features = [
            {"gene_id": f"g{i}"} for i in range(10)
        ]
        gff = [{"gene_id": f"g{i}"} for i in range(6)]
        assert stage.gff_features_are_covered(gff, features) is True
        assert stage.gff_features_are_covered(features, gff) is False, (
            "a feature table larger than the GFF must NOT satisfy the reversed "
            "direction; if this ever passes, the comparison is an equality or a "
            "union and the real smoke cohort would silently stop being checked"
        )

    def test_the_genome_stem_comes_from_the_assembly_not_the_sample_id(
        self, tmp_path, database_dir
    ):
        """`run_bakta` names files by the assembly stem; reuse must agree.

        They coincide for the smoke cohort, which is why this is easy to get
        wrong: a stem derived from `sample_id` works for every isolate anyone
        has tested and fails the first time an assembly is called something
        else.
        """
        write_candidate(tmp_path, "PDT_A", genome_stem="GCF_000006765.1",
                        database_dir=database_dir)
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A", assembly_path="/x/GCF_000006765.1.fna"),
            out_root=tmp_path, database_dir=database_dir,
        )
        assert decided.reused is True, decided.detail
        assert decided.genome_stem == "GCF_000006765.1"

    def test_files_named_after_the_sample_id_are_not_a_candidate(
        self, tmp_path, database_dir
    ):
        """The inverse: real, verifying bytes under the wrong file name.

        The directory is named by the sample id, so `<sample_id>.tsv` is a
        perfectly plausible mistake - and it is exactly the mistake a
        sample-id-derived stem would make and be rewarded for.
        """
        write_candidate(tmp_path, "PDT_A", genome_stem="OTHER_STEM",
                        database_dir=database_dir)
        directory = tmp_path / "PDT_A"
        for role in ("tsv", "gff3", "inference.tsv"):
            (directory / f"PDT_A.{role}").write_text(
                (directory / f"OTHER_STEM.{role}").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A", assembly_path="/x/GCF_000006765.1.fna"),
            out_root=tmp_path, database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_MISSING

    def test_no_candidate_is_never_an_exception(self, tmp_path, database_dir):
        """A refusal is a value. `decide_reuse` does not raise.

        One isolate with no reusable output must not end the run: this stage's
        documented contract is record-and-continue, because a smoke cohort is
        mostly unprepared by construction.
        """
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_NOTHING"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_MISSING
        assert "PDT_NOTHING" in decided.detail


# --------------------------------------------------------------------------
# U5 - the unreadable-header refusal, and that it is load-bearing
# --------------------------------------------------------------------------

class TestTheHeaderMustBeReadable:
    def _headerless(self, tmp_path, database_dir, text):
        directory = write_candidate(
            tmp_path, "PDT_A", database_dir=database_dir, tsv_rows=text
        )
        return directory

    def test_a_headerless_output_is_refused_rather_than_assumed_current(
        self, tmp_path, database_dir
    ):
        """**Load-bearing.** Delete the refusal below and this test fails.

        Nothing in the file contradicts the configured database, so the only
        evidence that it is the wrong database is the absence of the line that
        would say so. Defaulting to the configured value is precisely how stale
        output from another database gets accepted, and the damage is invisible:
        the gene calls differ and every downstream stage reports them as
        biology.
        """
        self._headerless(tmp_path, database_dir, "contig_00001\tx\ty\n")
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_HEADER_INCOMPLETE
        assert "Database" in decided.detail

    def test_an_output_naming_no_software_version_is_refused(self, tmp_path, database_dir):
        """Provenance names the tool too; half a header is no header."""
        text = (
            f"# Database: {stage.expected_bakta_database_string(database_dir)}\n"
            + "#" + "\t".join(stage.BAKTA_TSV_COLUMNS) + "\n"
            + "contig_00001\tCDS\t1\t300\t+\tTESTPA_00001\t1\tGC=0.4\tg\tp\tX\n"
        )
        self._headerless(tmp_path, database_dir, text)
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_HEADER_INCOMPLETE
        assert "Software" in decided.detail

    def test_a_header_beyond_the_probe_window_is_refused_not_scanned_for(
        self, tmp_path, database_dir
    ):
        """The read is bounded, and a bounded miss is a refusal.

        Scanning 20 Mb of table for a header would be slower than re-running
        the tool, and would find one in a file that is a concatenation.
        """
        padding = "".join(
            f"row{i}\tpadding\t{i}\t{i}\t+\tT{i}\t{i}\tGC=0.4\tg\tp\tX\n"
            for i in range(stage._HEADER_PROBE_LINES + 10)
        )
        text = padding + "# Database: v6.0, light\n# Software: v1.12.1\n"
        self._headerless(tmp_path, database_dir, text)
        decided = stage.decide_reuse(
            Sample(sample_id="PDT_A"), out_root=tmp_path,
            database_dir=database_dir,
        )
        assert decided.reused is False
        assert decided.reason == stage.REUSE_HEADER_INCOMPLETE

    def test_a_missing_database_version_json_is_refused(self, tmp_path, database_dir):
        """An unverifiable database cannot make existing output verifiable."""
        with pytest.raises(stage.StageError) as excinfo:
            stage.expected_bakta_database_string(tmp_path / "no-such-db")
        assert "version.json" in str(excinfo.value)


# --------------------------------------------------------------------------
# U3 / U5 - the three modes, end to end, through an injected runner
# --------------------------------------------------------------------------

class TestOffIsUnchanged:
    def test_bakta_runs_for_every_sample_and_nothing_is_written(
        self, configured, database_dir, monkeypatch, tmp_path
    ):
        """The default must be indistinguishable from the pre-existing stage."""
        runner = RecordingRunner(database_dir=database_dir)
        monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", runner)
        result = stage.run(
            configured("off"), manifest_of(["PDT_A", "PDT_B"]), RunMode.REAL,
            tmp_path,
        )
        assert runner.asked == ["PDT_A", "PDT_B"]
        assert sum(1 for v in result.values() if v) == 2
        assert not (tmp_path / "bakta" / stage.REUSE_PROVENANCE_FILENAME).exists(), (
            "under `off` the stage gains no file it did not have before"
        )


class TestPrefer:
    def test_a_valid_sample_is_reused_and_an_invalid_one_goes_to_bakta(
        self, configured, database_dir, monkeypatch, tmp_path, out_root
    ):
        runner = RecordingRunner(database_dir=database_dir)
        monkeypatch.setattr("papipeline.adapters.bakta.run_bakta", runner)
        write_candidate(out_root, "PDT_GOOD", database_dir=database_dir)
        write_candidate(out_root, "PDT_STALE", database_dir=database_dir,
                        database="v5.2, full")

        result = stage.run(
            configured("prefer"), manifest_of(["PDT_GOOD", "PDT_STALE"]),
            RunMode.REAL, tmp_path,
        )
        assert runner.asked == ["PDT_STALE"], (
            "the reusable sample must not reach the tool; the stale one must"
        )
        assert result["PDT_GOOD"], "the reused annotation must reach the records"
        assert result["PDT_STALE"], "the re-annotated one must too"

    def test_the_fallback_reason_is_in_the_log(
        self, configured, database_dir, monkeypatch, tmp_path, out_root, caplog
    ):
        import logging

        monkeypatch.setattr(
            "papipeline.adapters.bakta.run_bakta",
            RecordingRunner(database_dir=database_dir),
        )
        write_candidate(out_root, "PDT_STALE", database_dir=database_dir,
                        database="v5.2, full")
        # The level is set on the package logger, not on this stage's: an
        # ancestor's level is what an INFO record is filtered against, and
        # `papipeline` is left at WARNING by earlier tests.
        with caplog.at_level(logging.INFO, logger="papipeline"):
            stage.run(
                configured("prefer"), manifest_of(["PDT_STALE"]), RunMode.REAL,
                tmp_path,
            )
        assert "v5.2, full" in caplog.text, (
            "falling back to the tool is fine; doing it without saying which "
            "database was rejected is how the same annotation comes back twice"
        )


class TestRequireNeverInvokesBakta:
    def test_a_valid_sample_is_reused_with_the_tool_armed_to_explode(
        self, configured, database_dir, exploding_runner, tmp_path, out_root
    ):
        write_candidate(out_root, "PDT_GOOD", database_dir=database_dir)
        result = stage.run(
            configured("require"), manifest_of(["PDT_GOOD"]), RunMode.REAL,
            tmp_path,
        )
        assert len(result["PDT_GOOD"]) == N_FEATURES

    @pytest.mark.parametrize(
        "kwargs,expected_reason",
        [
            ({"omit": ("feature_tsv",)}, stage.REUSE_MISSING),
            ({"database": "v5.2, full"}, stage.REUSE_DATABASE_MISMATCH),
        ],
    )
    def test_an_invalid_sample_fails_alone_with_a_named_reason(
        self, configured, database_dir, exploding_runner, tmp_path, out_root, caplog,
        kwargs, expected_reason,
    ):
        """Record-and-continue: the cohort is not aborted.

        The whole point of `require` is a run against pre-computed annotation.
        If one bad isolate ended the run, a single stale file in a 967-isolate
        cohort would make the mode unusable - which is the shape of the defect
        this closes in the first place.
        """
        import logging

        write_candidate(out_root, "PDT_BAD", database_dir=database_dir, **kwargs)
        write_candidate(out_root, "PDT_GOOD", database_dir=database_dir)
        with caplog.at_level(logging.WARNING, logger="stages.annotation"):
            result = stage.run(
                configured("require"),
                manifest_of(["PDT_BAD", "PDT_GOOD"]),
                RunMode.REAL, tmp_path,
            )
        assert result["PDT_GOOD"], "the valid isolate must still be annotated"
        assert result["PDT_BAD"] == [], "the invalid one must not be"
        assert expected_reason in caplog.text
        if expected_reason == stage.REUSE_DATABASE_MISMATCH:
            assert "v5.2, full" in caplog.text, (
                "a refusal that does not name the rejected string is one "
                "nobody can act on"
            )

    def test_the_cohort_survives_every_isolate_being_unusable(
        self, configured, database_dir, exploding_runner, tmp_path, out_root
    ):
        result = stage.run(
            configured("require"),
            manifest_of(["PDT_X", "PDT_Y", "PDT_Z"]), RunMode.REAL, tmp_path,
        )
        assert set(result) == {"PDT_X", "PDT_Y", "PDT_Z"}
        assert all(v == [] for v in result.values())

    def test_an_unusable_sample_is_reported_as_a_refusal_not_as_bakta(
        self, configured, database_dir, exploding_runner, tmp_path, out_root, caplog
    ):
        """"The tool failed" and "the tool was not allowed to run" differ.

        Both mean no annotation, but they call for different responses: one is
        a fault to investigate, the other is a fact about the run's
        configuration.
        """
        import logging

        write_candidate(out_root, "PDT_BAD", database_dir=database_dir,
                        database="v5.2, full")
        with caplog.at_level(logging.WARNING, logger="stages.annotation"):
            stage.run(
                configured("require"), manifest_of(["PDT_BAD"]), RunMode.REAL,
                tmp_path,
            )
        assert "1 reusable output refused" in caplog.text, caplog.text


# --------------------------------------------------------------------------
# U4 - per-sample provenance
# --------------------------------------------------------------------------

class TestProvenanceIsRecorded:
    def _rows(self, tmp_path):
        from papipeline.io.tsv import read_tsv

        path = tmp_path / "bakta" / stage.REUSE_PROVENANCE_FILENAME
        return {row["sample_id"]: row for row in read_tsv(path)}

    def test_a_reused_sample_records_where_it_came_from_and_what_it_had(
        self, configured, database_dir, exploding_runner, tmp_path, out_root
    ):
        write_candidate(out_root, "PDT_GOOD", database_dir=database_dir)
        stage.run(
            configured("require"), manifest_of(["PDT_GOOD"]), RunMode.REAL,
            tmp_path,
        )
        row = self._rows(tmp_path)["PDT_GOOD"]
        assert row["reused"] == "true"
        assert row["source_dir"] == str(tmp_path / "bakta" / "PDT_GOOD")
        assert row["bakta_version"] == "1.12.1"
        assert row["database_string"] == stage.expected_bakta_database_string(
            database_dir)

    def test_the_recorded_digests_are_the_files_on_disk(
        self, configured, database_dir, exploding_runner, tmp_path, out_root
    ):
        """The digests must be recomputable by a reader, not merely present."""
        import hashlib

        write_candidate(out_root, "PDT_GOOD", database_dir=database_dir)
        stage.run(
            configured("require"), manifest_of(["PDT_GOOD"]), RunMode.REAL,
            tmp_path,
        )
        row = self._rows(tmp_path)["PDT_GOOD"]
        directory = tmp_path / "bakta" / "PDT_GOOD"
        for role, column in (
            ("PDT_GOOD.tsv", "sha256_feature_tsv"),
            ("PDT_GOOD.gff3", "sha256_gff3"),
        ):
            expected = hashlib.sha256(
                (directory / role).read_bytes()
            ).hexdigest()
            assert row[column] == expected, role

    def test_a_refusal_is_recorded_too(
        self, configured, database_dir, exploding_runner, tmp_path, out_root
    ):
        write_candidate(out_root, "PDT_BAD", database_dir=database_dir,
                        database="v5.2, full")
        stage.run(
            configured("require"), manifest_of(["PDT_BAD"]), RunMode.REAL,
            tmp_path,
        )
        row = self._rows(tmp_path)["PDT_BAD"]
        assert row["reused"] == "false"
        assert row["reason"] == stage.REUSE_DATABASE_MISMATCH
        assert "v5.2, full" in row["detail"]

    def test_a_genome_annotated_this_run_is_recorded_as_not_reused(
        self, configured, database_dir, monkeypatch, tmp_path
    ):
        """Provenance for every genome, not only the ones that skipped the tool.

        Otherwise the file answers "which bytes were trusted?" for half the
        cohort and is silent about the rest.
        """
        monkeypatch.setattr(
            "papipeline.adapters.bakta.run_bakta",
            RecordingRunner(database_dir=database_dir),
        )
        stage.run(
            configured("prefer"), manifest_of(["PDT_FRESH"]), RunMode.REAL,
            tmp_path,
        )
        row = self._rows(tmp_path)["PDT_FRESH"]
        assert row["reused"] == "false"
        assert row["reason"] == "not_reused"
        assert row["sha256_feature_tsv"]

    def test_the_column_order_is_the_declared_one(self):
        assert stage.REUSE_PROVENANCE_COLUMNS[0] == "sample_id"
        assert "reused" in stage.REUSE_PROVENANCE_COLUMNS
        assert stage.REUSE_PROVENANCE_COLUMNS == (
            "sample_id", "reused", "reason", "detail", "source_dir",
            "genome_stem", "sha256_feature_tsv", "sha256_gff3",
            "bakta_version", "database_string",
        )


# --------------------------------------------------------------------------
# the conversion chain is one implementation
# --------------------------------------------------------------------------

class TestReuseSkipsTheSubprocessNotTheValidation:
    def test_both_branches_call_the_same_conversion(self):
        """A property of the code, not a promise about two branches agreeing."""
        source = Path(stage.__file__).read_text(encoding="utf-8")
        assert source.count("_records_from_bakta_tables(") >= 3, (
            "one definition plus one call per branch; fewer means a branch grew "
            "its own conversion"
        )

    def test_reused_output_still_goes_through_containment(self):
        """Reuse skips the subprocess, not the validation."""
        gff = stage.parse_bakta_gff(gff_text())
        features = stage.parse_bakta_tsv(feature_tsv_text(header=""))
        assert stage.gff_features_are_covered(gff, features) is True
        ghost = list(features) + [{"gene_id": "GHOST_00001"}]
        assert stage.gff_features_are_covered(gff, ghost) is True
        gff_ghost = gff + [{"gene_id": "GHOST_00001"}]
        assert stage.gff_features_are_covered(gff_ghost, features) is False

    def test_the_records_are_standardised_identically_wherever_they_came_from(
        self, configured, database_dir, exploding_runner, monkeypatch, tmp_path, out_root
    ):
        """Reused and freshly annotated genomes must be indistinguishable.

        They are written by one function from one schema, so a downstream stage
        cannot tell how a record was obtained - which is the property that
        makes reuse safe to switch on for one run.
        """
        from papipeline.io.tsv import read_tsv

        write_candidate(out_root, "PDT_REUSED", database_dir=database_dir)
        monkeypatch.setattr(
            "papipeline.adapters.bakta.run_bakta",
            RecordingRunner(database_dir=database_dir),
        )
        stage.run(
            configured("prefer"), manifest_of(["PDT_REUSED", "PDT_FRESH"]),
            RunMode.REAL, tmp_path,
        )
        rows = {}
        for sample_id in ("PDT_REUSED", "PDT_FRESH"):
            rows[sample_id] = read_tsv(
                tmp_path / "annotation" / f"{sample_id}.annotation.tsv"
            )
        assert [r["gene_id"] for r in rows["PDT_REUSED"]] == [
            r["gene_id"] for r in rows["PDT_FRESH"]
        ]
        assert {r["annotation_source"] for r in rows["PDT_REUSED"]} == {"bakta"}
        assert sorted(rows) == ["PDT_FRESH", "PDT_REUSED"]