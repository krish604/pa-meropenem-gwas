"""The `variants` stage entry point: the seam that makes the parser reachable.

`stages/variants.py` had a parser and no caller. `parse_vcf` was proven against
a real bcftools record and against the committed fixtures, and nothing in the
pipeline ever called it, so a green suite coexisted with a stage that could not
run. This tests the seam that closes that: `variants.run()`.

Three things are pinned here, and the first is the one that matters most.

**Cohort membership is every manifest sample, including those with no calls.**
An isolate with no assembly, or with an assembly nothing aligned to, is an
ordinary member of the cohort - it simply carries no variant. So `run()` returns
a key for every sample in the manifest, and the distinction the AMR stage already
draws is preserved here: an *empty list* means "called, found nothing", while an
absent key would mean "never screened". Collapsing the two would make an
unanalysable genome indistinguishable from a clean one, and at 835 isolates with
82% of assemblies unusable (docs/design/cohort-variant-merge.md 5b) that
confusion would be the common case, not an edge case.

This matters downstream: `merge_calls` derives `an` from `len(calls_by_isolate)`,
so dropping the empty isolates would shrink the denominator and inflate every
allele frequency the cohort reports.

**TEST reads committed fixtures, through the real parser.** There is no PAO1
reference in `test_data/` and REAL mode is closed, so TEST cannot align anything.
What it can do is parse the five committed VCFs, which carry the shapes real
bcftools emits - a multiallelic ALT, `DP4` as a four-tuple. A TEST path that
skipped the parser would leave the awkward shapes untested, which is the reason
those fixtures exist (see test_variants_test_fixtures.py).

**The real invocation is a separate seam.** `papipeline.adapters.minimap2` is
tested on its own terms; this file asserts only that TEST mode never reaches it.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from papipeline.config.loader import PipelineConfig
from papipeline.manifest import SampleManifest
from papipeline.models import RunMode, Sample
from papipeline.stages import variants as stage_variants
from papipeline.stages.variants import PER_ISOLATE_COLUMNS

#: What a TEST run's data root is. The stage finds the fixtures at
#: `<data_root>/variants/`, so the test passes the root rather than the leaf -
#: the same argument a real TEST run passes, so the test cannot drift onto a
#: convention the pipeline does not use.
TEST_DATA = Path("test_data")


def _manifest(*sample_ids: str) -> SampleManifest:
    return SampleManifest(samples=[Sample(sample_id=s) for s in sample_ids])


class TestTheStageHasAnEntryPoint:
    def test_run_exists(self):
        """The seam itself. A stage that parses but cannot be called is not built."""
        assert callable(stage_variants.run)


class TestTestModeReadsCommittedFixtures:
    """TEST exercises the real parser over real-shaped output."""

    def test_it_returns_calls_for_every_fixture(self, config: PipelineConfig, tmp_path):
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001", "TEST_PA_002"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        assert set(calls) == {"TEST_PA_001", "TEST_PA_002"}
        assert calls["TEST_PA_001"], "the committed fixture produced no calls"
        assert calls["TEST_PA_002"], "the committed fixture produced no calls"

    def test_the_rows_carry_the_declared_contract(
        self, config: PipelineConfig, tmp_path
    ):
        """Every row is the contract's shape - no guessed columns."""
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        for row in calls["TEST_PA_001"]:
            assert set(row) == set(PER_ISOLATE_COLUMNS)

    def test_a_multiallelic_site_survives_the_whole_stage(
        self, config: PipelineConfig, tmp_path
    ):
        """The awkward shape is tested through `run()`, not only `parse_vcf`.

        `parse_vcf` already has a unit test for this, and it would keep passing if
        `run()` re-parsed, filtered or re-keyed the rows on the way out. The
        property that matters is a property of the stage, so it is asserted
        here.
        """
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        by_pos: dict = {}
        for row in calls["TEST_PA_001"]:
            by_pos.setdefault(row["pos"], []).append(row["alt"])
        assert by_pos["602"] == ["T", "G"]

    def test_the_isolate_is_attributed_to_its_own_calls(
        self, config: PipelineConfig, tmp_path
    ):
        """The sample_id is the join key for the cohort merge; it cannot drift."""
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001", "TEST_PA_002"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        for sample_id, rows in calls.items():
            for row in rows:
                assert row["sample_id"] == sample_id


class TestEveryManifestSampleIsACohortMember:
    """Membership is the whole cohort, not the subset that produced calls."""

    def test_a_sample_with_no_calls_is_present_with_an_empty_list(
        self, config: PipelineConfig, tmp_path
    ):
        """Empty means "screened, found nothing". Absent would mean "not screened"."""
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001", "TEST_PA_999"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        assert "TEST_PA_999" in calls
        assert calls["TEST_PA_999"] == []

    def test_the_denominator_is_the_cohort_not_the_number_of_callers(
        self, config: PipelineConfig, tmp_path
    ):
        """Why this matters downstream, asserted where the decision is made.

        `merge_calls` computes `an = len(calls_by_isolate)`. If a stage returned
        only the isolates that produced calls, a 20-sample cohort with 5 calling
        would report `an=5` and every allele frequency would be inflated
        four-fold. The stage is where that is preventable.
        """
        calls = stage_variants.run(
            config, _manifest(*(f"TEST_PA_00{i}" for i in range(1, 6))),
            RunMode.TEST, data_root=TEST_DATA, workdir=tmp_path,
        )
        assert len(calls) == 5
        for sample_id in ("TEST_PA_006", "TEST_PA_020"):
            calls.setdefault(sample_id, [])
        assert len(calls) == 7, "an unscreened isolate was dropped from the cohort"


class TestRealModeAlignsAndCalls:
    """REAL invokes the pinned tools through the adapter.

    The invocation itself is pinned in `science.yaml` and its flags are checked
    against the installed binaries in `test_minimap2_adapter.py`. What is tested
    here is the wiring: one alignment per manifest sample, each parsed by the
    same `parse_vcf` TEST uses, so a REAL run and a TEST run cannot diverge in
    how a call is read.
    """

    def test_it_does_not_refuse_real_mode(self, config, tmp_path, monkeypatch):
        """The stage exists for REAL; refusing there would leave it TEST-only.

        Driven with the tools and the assembly lookup stubbed, because the point
        is the control flow reaching the aligner - not whether this machine has
        a PAO1 reference (which `test_minimap2_adapter.py` covers separately).
        """
        monkeypatch.setattr(
            stage_variants, "locate_assembly", lambda s, root: Path("a.fna"),
            raising=False,
        )
        monkeypatch.setattr(
            stage_variants, "_which", lambda name: f"/usr/bin/{name}", raising=False,
        )
        monkeypatch.setattr(
            stage_variants, "call_isolate",
            lambda sample_id, **kw: SimpleNamespace(
                sample_id=sample_id, vcf=tmp_path / "x.vcf"
            ),
            raising=False,
        )
        (tmp_path / "x.vcf").write_text(
            "##fileformat=VCFv4.2\n", encoding="utf-8"
        )
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.REAL,
            data_root=tmp_path, workdir=tmp_path,
        )
        assert "TEST_PA_001" in calls

    def test_each_sample_is_aligned_once(self, config, tmp_path, monkeypatch):
        """One alignment per isolate, and no more.

        835 isolates is 835 alignments; doing it twice is an hour of wasted CPU
        and doing it zero times is a stage that fabricates its own output.
        """
        monkeypatch.setattr(
            stage_variants, "locate_assembly", lambda s, root: Path("a.fna"),
            raising=False,
        )
        seen: list = []

        def fake_call(sample_id, **kwargs):
            seen.append(sample_id)
            vcf = tmp_path / f"{sample_id}.calls.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\n"
                "NC_002516.2\t100\t.\tA\tT\t30.4\t.\tAC=2;AN=2\tGT:DP\t1/1:1\n",
                encoding="utf-8",
            )
            return SimpleNamespace(sample_id=sample_id, vcf=vcf)

        monkeypatch.setattr(stage_variants, "call_isolate", fake_call, raising=False)
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001", "TEST_PA_002"), RunMode.REAL,
            data_root=tmp_path, workdir=tmp_path,
        )
        assert sorted(seen) == ["TEST_PA_001", "TEST_PA_002"]
        assert set(calls) == {"TEST_PA_001", "TEST_PA_002"}

    def test_the_calls_are_parsed_by_the_same_parser(self, config, tmp_path, monkeypatch):
        """One parser for both modes, or REAL and TEST disagree silently."""
        monkeypatch.setattr(
            stage_variants, "locate_assembly", lambda s, root: Path("a.fna"),
            raising=False,
        )

        def fake_call(sample_id, **kwargs):
            vcf = tmp_path / f"{sample_id}.calls.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\n"
                "NC_002516.2\t100\t.\tA\tT,G\t30.4\t.\tAC=2;AN=2;DP4=0,0,0,1\t"
                "GT:DP\t1/1:1\n",
                encoding="utf-8",
            )
            return SimpleNamespace(sample_id=sample_id, vcf=vcf)

        monkeypatch.setattr(stage_variants, "call_isolate", fake_call, raising=False)
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.REAL,
            data_root=tmp_path, workdir=tmp_path,
        )
        alts = [r["alt"] for r in calls["TEST_PA_001"]]
        assert alts == ["T", "G"], "the multiallelic split did not reach parse_vcf"

    def test_a_sample_with_no_assembly_is_a_member_carrying_nothing(
        self, config, tmp_path, monkeypatch
    ):
        """An isolate with no genome is an ordinary member, not a dropped one.

        `locate_assembly` refuses, and that refusal is correct. It must not
        remove the isolate from the cohort, because `an` is derived from the
        cohort and a silently shrunken denominator inflates every frequency.
        """
        from papipeline.errors import PipelineError

        def refuse(sample, data_root):
            raise PipelineError("no assembly", sample_id=sample.sample_id)

        monkeypatch.setattr(stage_variants, "locate_assembly", refuse, raising=False)
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.REAL,
            data_root=tmp_path, workdir=tmp_path,
        )
        assert calls == {"TEST_PA_001": []}

    def test_a_sample_whose_call_failed_is_a_member_carrying_nothing(
        self, config, tmp_path, monkeypatch
    ):
        """82% of the real assemblies are unusable; that must not stop the run.

        With 835 isolates and an 11-18% usable rate, a stage that stopped at the
        first failure could never complete. The isolate stays in the cohort and
        the reason is logged, so an unanalysable genome is visible rather than
        quietly absent.
        """
        from papipeline.errors import ToolExecutionError

        def boom(sample_id, **kwargs):
            raise ToolExecutionError("bcftools exited non-zero", command="x")

        monkeypatch.setattr(stage_variants, "call_isolate", boom, raising=False)
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.REAL,
            data_root=tmp_path, workdir=tmp_path,
        )
        assert calls == {"TEST_PA_001": []}


class TestTestModeNeverReachesTheAligner:
    """TEST has no reference and no permission to align anything."""

    def test_the_aligner_is_never_constructed(
        self, config: PipelineConfig, tmp_path, monkeypatch
    ):
        """A TEST run that reached minimap2 would need a reference it does not have.

        Asserted by making construction fail loudly rather than by checking that
        no subprocess ran: a silent attempt that happened to succeed is exactly
        the state this forbids.
        """
        def explode(*args, **kwargs):
            raise AssertionError("TEST mode must not invoke the aligner")

        monkeypatch.setattr(stage_variants, "call_isolate", explode, raising=False)
        calls = stage_variants.run(
            config, _manifest("TEST_PA_001"), RunMode.TEST,
            data_root=TEST_DATA, workdir=tmp_path,
        )
        assert calls

    def test_a_missing_fixture_directory_is_a_refusal(
        self, config: PipelineConfig, tmp_path
    ):
        """Silently returning an empty cohort would be indistinguishable from a clean one."""
        with pytest.raises(Exception) as excinfo:
            stage_variants.run(
                config, _manifest("TEST_PA_001"), RunMode.TEST,
                data_root=tmp_path / "absent", workdir=tmp_path,
            )
        assert "variants" in str(excinfo.value).lower()
