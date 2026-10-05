"""The REAL Snakemake DAG resolves, and names what an operator has to provide.

**The defect this file exists for.** `rule phylogeny` used to declare
`tree.nwk` and `tree_metadata.tsv` in its `external` list -
`stage_inputs([OUT_RECOMBINATION], [METADATA, TREE_IN, TREE_META_IN])` - while
being the only rule that produces them. `stages/phylogeny.build_tree_outputs`
writes both into `CONFIG.phylogeny_dir(RESOLVED_MODE)`. So the REAL DAG asked
Snakemake to resolve as *inputs* the very files the rule was supposed to
produce, and the only way to satisfy it was to place a tree on disk by hand. A
DAG that builds while stage 9 has nothing to do is worse than one that refuses.

**Why a dry run in TEST could never have caught it.**
`tests/integration/test_snakemake_dry_run.py` dry-ran in TEST, where
`test_data/phylogeny/tree.nwk` and `tree_metadata.tsv` are committed fixtures
that exist on disk. In TEST the self-edge is invisible: the files are there, so
Snakemake is satisfied whether or not a rule declares them. The self-edge only
shows up in REAL, where the files do not exist until stage 9 writes them.

**How the REAL dry run is made.** REAL is refused at DAG-build time by
`resolve_mode`, so the gate is opened for the session with
`PIPELINE_ALLOW_REAL_MODE=1` - the override the refusal itself names. That is
not the same as setting `runtime.allow_real_mode` in an overlay, which is the
user's decision and not a test's. Nothing is executed: `--dry-run` builds the
DAG and stops.

**Why a throwaway overlay.** Two of the externals live under `data/`
(`sample_metadata.tsv`, `imipenem_phenotype.tsv`), and in a dev worktree `data/`
is a **symlink into the canonical read-only worktree** - which is also why
`config/loader.py` documents `data/phenotype/` as empty. Provisioning a real
cohort manifest there would write into another worktree's disk. So the test
redirects `paths.data_root` to a `tmp_path` through a temporary machine overlay,
the mechanism `test_snakemake_dry_run.py` already uses to probe a missing
input. The committed overlay is read, copied and patched; it is never edited.

**What the test then does about operator-provisioned externals.** It cannot
assert the REAL DAG resolves outright, because a REAL DAG legitimately needs
inputs no rule here produces. So it *provisions them itself*, iteratively: run
the dry run, read the files Snakemake named, create empty placeholders, repeat
until Snakemake stops complaining. The loop terminating is the assertion: it
means the DAG's only unsatisfied edges are the ones the operator owns, and no
edge is unsatisfiable except by hand-placing a file a rule claims to produce.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from papipeline.config.loader import (
    ALLOW_REAL_MODE_ENV,
    RESULTS_ROOT_ENV,
)
from papipeline.run import STAGE_ORDER

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = PIPELINE_ROOT / "workflow" / "Snakefile"

pytestmark = pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake is not installed; this check cannot run here",
)

#: The overlays whose DAG is checked here. `smoke` is now a KNOWN name
#: (`loader.KNOWN_MACHINES`, round 12 Phase 3 - it was reachable only as a path,
#: which the Snakefile interpolated verbatim into a shell command and so
#: resolved against the CWD), but it is deliberately NOT in this list: its
#: `data_root` is the canonical read-only worktree through the `data` symlink,
#: and `_temp_overlay` patches `paths.data_root` to redirect it. Under a smoke
#: overlay the cohort comes from `PDC_essential.tsv` and the assemblies from
#: `db/smoke_genomes/`, so redirecting `data_root` alone would not stop this
#: test provisioning into another worktree. It is covered by
#: `tests/integration/test_smoke_overlay.py`, whose fixture asserts the same
#: property.
MACHINES = ("laptop", "bigmachine")


def _temp_overlay(machine: str, tmp_path: Path) -> Path:
    """A copy of `config/machines/<machine>.yaml` with `data_root` redirected.

    `load_config` already accepts a path in place of a machine name, so this uses
    an existing seam rather than adding one. The patch is asserted to have
    applied: a silent no-op here would mean the test provisioning files into the
    canonical worktree through the `data/` symlink, which is the one thing it
    must never do.
    """
    source = (PIPELINE_ROOT / "config" / "machines" / f"{machine}.yaml").read_text(
        encoding="utf-8"
    )
    data_root = tmp_path / "data"
    (data_root / "metadata").mkdir(parents=True, exist_ok=True)
    (data_root / "phenotype").mkdir(parents=True, exist_ok=True)
    patched = source.replace("  data_root: data\n", f"  data_root: {data_root}\n")
    assert patched != source, (
        f"config/machines/{machine}.yaml no longer declares `  data_root: data`, "
        "so this test would provision the cohort manifest through the `data/` "
        "symlink into the canonical read-only worktree. Update the patch."
    )
    overlay = tmp_path / "machines" / f"{machine}.yaml"
    overlay.parent.mkdir(parents=True, exist_ok=True)
    overlay.write_text(patched, encoding="utf-8")
    return overlay


def _dry_run(overlay: Path, results_root: Path):
    return subprocess.run(
        [
            "snakemake",
            "--snakefile", str(SNAKEFILE),
            "--config", "mode=REAL", f"machine={overlay}",
            "--cores", "1",
            "--dry-run",
        ],
        cwd=PIPELINE_ROOT,
        env={
            **os.environ,
            RESULTS_ROOT_ENV: str(results_root),
            ALLOW_REAL_MODE_ENV: "1",
        },
        capture_output=True,
        text=True,
        timeout=600,
    )


def _affected_files(result: subprocess.CompletedProcess) -> list:
    """The paths Snakemake named as missing, parsed from its own message."""
    combined = result.stdout + result.stderr
    lines = combined.splitlines()
    files, collecting = [], False
    for line in lines:
        if line.strip() == "affected files:":
            collecting = True
            continue
        if not collecting:
            continue
        stripped = line.strip()
        if not stripped.startswith("/"):
            collecting = False
            continue
        files.append(stripped)
    return files


def _resolve_real_dag(machine: str, tmp_path: Path, limit: int = 25):
    """Provision the operator's externals until the REAL DAG builds.

    Returns ``(resolved, provisioned, output)``. Empty placeholders, never
    content: the point is whether the DAG's edges can be satisfied at all, not
    what the files say. A non-empty placeholder would be a fabricated result on
    disk and this is a test.
    """
    overlay = _temp_overlay(machine, tmp_path)
    results_root = tmp_path / "results"
    provisioned: list = []
    for _ in range(limit):
        result = _dry_run(overlay, results_root)
        combined = result.stdout + result.stderr
        if "MissingInputException" not in combined and "MissingOutputException" not in combined:
            return True, provisioned, combined
        missing = _affected_files(result)
        if not missing:
            return False, provisioned, combined
        for path in missing:
            target = Path(path)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("", encoding="utf-8")
            provisioned.append(path)
    return False, provisioned, ""


class TestTheRealDagResolves:
    @pytest.mark.parametrize("machine", MACHINES)
    def test_it_builds_once_the_operator_provisions_the_named_externals(
        self, machine: str, tmp_path: Path,
    ):
        resolved, provisioned, output = _resolve_real_dag(machine, tmp_path)
        assert resolved, (
            f"the REAL DAG did not resolve for machine={machine} after "
            f"{len(provisioned)} provision(s). Snakemake was still reporting "
            f"missing files:\n{output[-3000:]}"
        )

    @pytest.mark.parametrize("machine", MACHINES)
    def test_nothing_was_provisioned_outside_the_temporary_root(
        self, machine: str, tmp_path: Path,
    ):
        """The `data/` symlink guard, asserted rather than assumed.

        `data/`, `db/` and `local/` are symlinks into the canonical read-only
        worktree. A test that provisions a cohort manifest through one of them
        is writing into another agent's disk, and nothing about that failure
        would look like a failure.
        """
        _, provisioned, _ = _resolve_real_dag(machine, tmp_path)
        outside = [
            p for p in provisioned
            if not str(Path(p).resolve()).startswith(str(tmp_path.resolve()))
        ]
        assert not outside, (
            "these provisioned paths are outside tmp_path, so they reached "
            f"another worktree: {outside}"
        )

    @pytest.mark.parametrize("machine", MACHINES)
    def test_every_stage_is_scheduled_in_real(self, machine: str, tmp_path: Path):
        """Reachability, in REAL, which is where the self-edge lived.

        `test_every_stage_is_reachable_from_all` asserts this in TEST. It could
        not catch the `tree.nwk` defect precisely because TEST's fixtures make
        the self-edge resolve, so the same assertion is made here against a REAL
        DAG.
        """
        resolved, _, output = _resolve_real_dag(machine, tmp_path)
        assert resolved, f"REAL DAG did not resolve:\n{output[-3000:]}"
        scheduled = {
            m.group(1) for m in re.finditer(r"(?m)^\s*rule (\w+):", output)
        }
        missing = sorted(set(STAGE_ORDER) - scheduled)
        assert not missing, (
            f"these stages are never scheduled from `all` in REAL for "
            f"machine={machine}: {missing}"
        )

    @pytest.mark.parametrize("machine", MACHINES)
    def test_the_only_unresolved_edges_are_operator_provisioned(
        self, machine: str, tmp_path: Path,
    ):
        """What the operator has to supply, named - the deliverable of A1.

        Asserted as a set so a NEW unsatisfiable edge cannot be added without
        this failing and someone having to say why the new file is the
        operator's and not a rule's.
        """
        _, provisioned, _ = _resolve_real_dag(machine, tmp_path)
        names = {Path(p).name for p in provisioned}
        assert names == EXPECTED_REAL_EXTERNALS, (
            "the set of externals a REAL run needs has changed. Every one of "
            "these is written by a tool this pipeline does not invoke, or by "
            "the operator; a new name means either a tool gained a producer "
            "(then it should not be an external) or a rule stopped producing "
            "something it used to (then that is a bug, not an external). "
            f"unexpected={sorted(names - EXPECTED_REAL_EXTERNALS)} "
            f"missing={sorted(EXPECTED_REAL_EXTERNALS - names)}"
        )


#: Every file a REAL run needs that no rule in `workflow/Snakefile` produces.
#:
#: Each is operator-provisioned, and the reason is in the comment. This is the
#: list A1 asks to be reported, asserted so it cannot grow silently.
EXPECTED_REAL_EXTERNALS = frozenset({
    # The cohort manifest. Read by every stage, via `COHORT_INPUT`. The pipeline
    # discovers it; it does not write it. `rule synthetic_fixtures` used to
    # declare it as an OUTPUT, which meant a REAL run with no manifest had
    # twenty TEST_PA_* rows fabricated for it - see SYNTHETIC_FIXTURE_OUTPUTS in
    # the Snakefile.
    "sample_metadata.tsv",
    # Imipenem MIC/category per isolate. Read by stages 11 and 12. This is
    # clinical data (AGENTS.md: PDC_essential.tsv and the phenotype table are
    # never committed) and is supplied from the study's own records.
    "imipenem_phenotype.tsv",
    # panaroo's core gene alignment - stage 8's input. **The pipeline never runs
    # panaroo**: `papipeline/adapters/panaroo.py` imports no `subprocess`
    # (asserted in tests/unit/test_recombination_dispatch.py). This is the one
    # external with a documented command; `rule recombination` refuses by name
    # when it is absent.
    "core_gene_alignment_filtered.aln",
})

#: The five names REMOVED from :data:`EXPECTED_REAL_EXTERNALS` in round 12
#: Phase 3, with what made each one a self-edge.
#:
#: They were listed as operator-provisioned with rationales that were false:
#: "the pipeline does not invoke AMRFinderPlus" when `amr._run_by_calling_tool`
#: invokes it (`amr.py:589-591`), "the pipeline does not invoke mlst" when
#: `mlst.py:212-263` does, and so on. And the assertion could not catch the
#: falsehood, because its own provisioning loop manufactured an empty placeholder
#: for every unsatisfied edge - the very files that would have exposed the
#: self-edge this file was written to kill.
#:
#: Each is now written by the rule that was demanding it:
#:
#: * ``mlst_results.tsv``   - ``stages/mlst.py:263``
#: * ``amr_determinants.tsv`` - ``stages/amr.py:620``
#: * ``structural_variants.tsv`` - ``stages/sv.py:290``
#: * ``virulence_factors.tsv`` - ``stages/virulence.py:203``
#: * ``gwas_features.tsv``  - ``stages/gwas_features.py:751``, called by
#:   ``run.py`` inside its ``pangenome`` branch
#:
#: so they are TEST-only inputs in the Snakefile and REAL externals of nothing.
#: The set above is asserted, so a genuinely new external cannot join it without
#: someone saying why the new file is the operator's.
REMOVED_SELF_EDGES = frozenset({
    "mlst_results.tsv",
    "amr_determinants.tsv",
    "structural_variants.tsv",
    "virulence_factors.tsv",
    "gwas_features.tsv",
})


class TestTheSelfEdgesAreGone:
    """D3/D4/D5/D6, asserted on the declaration rather than on a run.

    `test_the_only_unresolved_edges_are_operator_provisioned` now fails if any of
    these COMES BACK, which is the direction that matters. These assert they are
    gone from the declaration, with the writer named, so a future edit that
    re-adds one as an external fails with a message about the self-edge rather
    than about a set comparison.
    """

    @pytest.mark.parametrize(
        "name", sorted(REMOVED_SELF_EDGES)
    )
    def test_it_is_not_in_the_external_set(self, name):
        assert name not in EXPECTED_REAL_EXTERNALS, (
            f"{name} is recorded as operator-provisioned again. It is written "
            "by the rule that demands it, so a REAL run could only satisfy the "
            "edge by hand-placing a file the rule itself produces - a DAG that "
            "builds while the stage has nothing to do."
        )

    @pytest.mark.parametrize("name", sorted(REMOVED_SELF_EDGES))
    def test_it_is_demanded_through_test_only_and_nothing_else(self, name):
        """The declaration, read off the Snakefile.

        Stated mechanically rather than by eye: the constant may appear in a
        rule's input block only through ``test_only=``, the channel whose whole
        purpose is "demanded in TEST and never in REAL". Anything else means the
        self-edge is back, and the message says which rule brought it back.
        """
        constant = {
            "mlst_results.tsv": "MLST_IN",
            "amr_determinants.tsv": "AMR_IN",
            "structural_variants.tsv": "SV_IN",
            "virulence_factors.tsv": "VF_IN",
            "gwas_features.tsv": "GWAS_IN",
        }[name]
        source = SNAKEFILE.read_text(encoding="utf-8")
        for body in re.split(r"(?m)^rule ", source)[1:]:
            if constant not in body:
                continue
            rule_name = body.split(":", 1)[0]
            input_block = body.split("output:", 1)[0]
            assert "test_only=" in input_block, (
                f"rule {rule_name!r} demands {constant} outside the "
                "`test_only` channel. It is written by that rule's own stage, so "
                "demanding it in REAL is the self-edge D3/D4/D5/D6 - the defect "
                "this file was written for."
            )


class TestTheSelfEdgeIsGone:
    """The specific defect, asserted on the declaration rather than on a run.

    `test_it_builds_once_the_operator_provisions_the_named_externals` would
    catch a regression, but only after provisioning its way to the failing rule.
    These read the Snakefile, so a re-introduced self-edge is caught in
    milliseconds and with a message that says what it is.
    """

    def _rule_body(self, name: str) -> str:
        source = SNAKEFILE.read_text(encoding="utf-8")
        match = re.search(rf"(?ms)^rule {name}:\n(.*?)(?=^rule |\Z)", source)
        assert match, f"rule {name!r} is not in the Snakefile"
        return match.group(1)

    def test_the_tree_is_an_output_of_the_rule_that_writes_it(self):
        body = self._rule_body("phylogeny")
        output = re.search(r"(?m)^    output:\n((?:        .*\n)+)", body)
        assert output, "rule phylogeny declares no output block"
        assert "TREE_IN" in output.group(1), (
            "tree.nwk is written by stages/phylogeny.build_tree_outputs and by "
            "no other rule. Declaring it as anything but an output of this rule "
            "means the DAG demands it as an input it cannot produce."
        )
        assert "TREE_META_IN" in output.group(1), (
            "tree_metadata.tsv is written by the same function and by no other "
            "rule; same reasoning as TREE_IN."
        )

    def test_the_tree_is_not_also_declared_as_an_input_by_its_producer(self):
        """A file cannot be this rule's input and its output at the same time."""
        body = self._rule_body("phylogeny")
        declared = re.search(r"(?m)^    input:\n(.*?)\n    output:", body, re.S)
        assert declared, "rule phylogeny declares no input block"
        assert "TREE_IN" not in declared.group(1), (
            "rule phylogeny declares tree.nwk as an input AND as an output. "
            "Snakemake resolves the self-edge as satisfied-by-anything and the "
            "REAL DAG stops being a check on stage 9."
        )

    def test_downstream_readers_get_an_edge_not_a_hand_placed_file(self):
        """`similarity` and `cooccurrence` read the crosswalk; they must not
        declare it external.

        They may still name `TREE_META_IN` - Snakemake links the path to the
        rule that produces it, which is the point. What would break it is
        removing `rule phylogeny`'s output declaration while these keep
        declaring the path, because then nothing produces it and the DAG stops
        on a file the operator cannot supply.
        """
        source = SNAKEFILE.read_text(encoding="utf-8")
        assert "rule phylogeny:" in source
        phylogeny = self._rule_body("phylogeny")
        assert "TREE_META_IN," in phylogeny, (
            "rule phylogeny must declare tree_metadata.tsv as an output, or the "
            "four rules that read it are demanding an operator-provisioned file"
        )

    def test_the_fixture_generator_does_not_fabricate_a_real_manifest(self):
        """`sample_metadata.tsv` is an external in REAL, not a fixture product.

        Before this, `rule synthetic_fixtures` declared `output: METADATA`
        unconditionally, so in REAL Snakemake never checked the manifest and
        would have fabricated twenty `TEST_PA_*` rows for a run whose operator
        had supplied none. The rule now declares an empty output list in REAL.
        """
        source = SNAKEFILE.read_text(encoding="utf-8")
        assert re.search(
            r"SYNTHETIC_FIXTURE_OUTPUTS\s*=\s*\(\s*\[\]\s*if RESOLVED_MODE is RunMode\.REAL",
            source,
        ), (
            "SYNTHETIC_FIXTURE_OUTPUTS must be empty in REAL. A REAL run whose "
            "operator supplied no cohort manifest would otherwise have one "
            "fabricated from the synthetic fixtures."
        )


class TestStubWritesWhatItsRuleDeclares:
    """Declaring `tree.nwk` an output obliges STUB to write it.

    Snakemake checks a rule's declared outputs after the job runs and stops on
    a missing one. So making `rule phylogeny` declare the tree as an output -
    which is what removes the self-edge - means the STUB path has to fabricate
    it, or a STUB run fails with a `MissingOutputException` naming a file the
    rule was supposed to write. This is the check that the two stayed together.
    """

    def test_a_stub_run_writes_a_tree_and_a_crosswalk(self, tmp_path: Path):
        from papipeline.config.loader import RESULTS_ROOT_ENV
        from papipeline.run import run_pipeline

        previous = os.environ.get(RESULTS_ROOT_ENV)
        os.environ[RESULTS_ROOT_ENV] = str(tmp_path)
        try:
            result = run_pipeline(mode="STUB")
        finally:
            if previous is None:
                os.environ.pop(RESULTS_ROOT_ENV, None)
            else:
                os.environ[RESULTS_ROOT_ENV] = previous

        phylogeny_dir = Path(result.outputs["tree"]).parent
        assert (phylogeny_dir / "tree.nwk").is_file(), (
            "rule phylogeny declares tree.nwk as an output, so STUB must write "
            "one or Snakemake stops on a missing output it was told to expect"
        )
        assert (phylogeny_dir / "tree_metadata.tsv").is_file()
        assert "STUB_NOT_A_TREE" in (phylogeny_dir / "tree.nwk").read_text(
            encoding="utf-8"
        ), (
            "the fabricated tree must say what it is. A file named tree.nwk "
            "holding a plausible topology is a fabricated phylogenetic result "
            "wearing a result's shape."
        )
        crosswalk = (phylogeny_dir / "tree_metadata.tsv").read_text(encoding="utf-8")
        assert crosswalk.strip().splitlines() == [
            "\t".join(("sample_id", "tree_tip_label"))
        ], (
            "a header and no rows. A crosswalk over zero samples cannot be told "
            "from a fabricated one by shape, so the only safe content is none"
        )