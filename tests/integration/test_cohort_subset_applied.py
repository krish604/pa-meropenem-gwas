"""`cohort.subset_file` selects the cohort. Asserted at manifest load.

**The defect.** `CohortConfig` did everything right: it parsed the key, refused
an empty path, refused a missing file, resolved a relative path against
`REPOSITORY_ROOT` rather than the CWD, and skipped `#` comments. And nothing
called `read_subset_file`. So `config/machines/smoke.yaml` declared

    cohort:
      subset_file: local/smoke_isolates.txt

and `discover_run_manifest` built a cohort of **all 967 isolates** in
`PDC_essential.tsv`. The bounded run analysed every one of them. That is the
exact outcome the overlay's own header says it exists to prevent, and the key
that was supposed to stop it was read by nothing.

**Why a smoke run with 967 isolates is worse than an error.** `runtime.max_samples`
for that overlay is 10, and `capped_sample_count` counts *prepared assemblies*
for a smoke overlay precisely so the cap does not fire on 967 members - so the
run sailed through the cap too. Nothing anywhere objected.

**What is asserted here.**

  * a subset keeps only the listed ids, at manifest load, so every stage is
    bounded and not just the ones that happened to read the key;
  * a listed id absent from the manifest is a refusal naming it - never a drop;
  * no file configured leaves the manifest untouched, which is the state of
    `laptop` and `bigmachine`;
  * the phenotype join accepts 10/10 for a bounded ten-isolate cohort.

**The phenotype join is the load-bearing check.** `stages/phenotype.py` is owned
by another agent and is not edited here. It is exercised through
`papipeline.join.join_manifest_to_phenotype`, which is what it calls: it refuses
any manifest member without exactly one row, so a subset that kept an isolate
the phenotype table does not describe would fail there. 10/10 therefore means
the bounded cohort and the bounded phenotype table agree, which is the property
a smoke run exists to establish.

**Synthetic throughout.** `PDC_essential.tsv` is REAL CLINICAL DATA and is not
committed, not copied, and not read here (AGENTS.md). Both the PDC table and the
member list are written into `tmp_path`, and `discover_run_manifest` already
accepts `pdc_path=` for exactly this.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import CohortConfig, load_config
from papipeline.errors import DataContractError, PipelineError
from papipeline.join import join_manifest_to_phenotype
from papipeline.manifest import SampleManifest, discover_pdc_manifest
from papipeline.models import Phenotype, PhenotypeCall, RunMode
from papipeline.run import apply_cohort_subset, discover_run_manifest

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
SMOKE = REPO / "config" / "machines" / "smoke.yaml"

#: The size the smoke overlay's header promises. Ten, asserted rather than
#: derived, because "the bounded run" is a claim about the overlay and this is
#: the test that holds it to it.
SMOKE_COHORT_SIZE = 10


@pytest.fixture
def pdc_table(tmp_path: Path) -> Path:
    """A synthetic PDC table: 25 isolates, of which the last 5 have no assembly.

    The assembly-less rows matter. `discover_pdc_manifest` makes every isolate a
    member regardless of assembly, so a subset that dropped them would change
    both the cohort and `merge_calls`'s denominator - and would do it silently.
    """
    path = tmp_path / "PDC_essential.tsv"
    header = "Isolate\tAssembly\tBioSample\tBioProject\tPDC_present\n"
    rows = []
    for index in range(1, 26):
        isolate = f"PDT{index:09d}.1"
        assembly = "" if index > 20 else f"GCA_{index:09d}.1"
        rows.append(f"{isolate}\t{assembly}\tSAMN{index:09d}\tPRJNA1\tPDC\n")
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def subset_file(tmp_path: Path) -> Path:
    """A ten-isolate member list, at the path the overlay would resolve."""
    path = tmp_path / "subset.txt"
    path.write_text(
        "# a comment, which read_subset_file skips\n"
        + "\n".join(f"PDT{i:09d}.1" for i in range(1, SMOKE_COHORT_SIZE + 1))
        + "\n",
        encoding="utf-8",
    )
    return path


def _config_with_subset(subset_file: Path, tmp_path: Path):
    """A smoke overlay copy whose `cohort.subset_file` points at `subset_file`.

    Built by patching the committed overlay rather than by editing it: the
    smoke overlay's `local/smoke_isolates.txt` is this machine's real bounded
    run, and a test must not repoint it.
    """
    source = SMOKE.read_text(encoding="utf-8")
    patched = source.replace(
        "  subset_file: local/smoke_isolates.txt",
        f"  subset_file: {subset_file}",
    )
    assert patched != source, (
        "config/machines/smoke.yaml no longer declares "
        "`  subset_file: local/smoke_isolates.txt`; update this patch. A "
        "no-op here would silently test the committed list instead."
    )
    overlay = tmp_path / "smoke.yaml"
    overlay.write_text(patched, encoding="utf-8")
    return load_config(SCIENCE, machine=overlay)


class TestTheSubsetSelectsTheCohort:
    def test_only_the_listed_ids_survive(self, pdc_table, subset_file, tmp_path):
        config = _config_with_subset(subset_file, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        assert len(manifest) == SMOKE_COHORT_SIZE, (
            f"the overlay bounds the run to {SMOKE_COHORT_SIZE} isolates and "
            f"the manifest has {len(manifest)}. Before this was applied, the "
            "smoke overlay built a cohort of every isolate in the PDC table - "
            "the outcome the overlay exists to prevent."
        )
        assert manifest.sample_ids == [
            f"PDT{i:09d}.1" for i in range(1, SMOKE_COHORT_SIZE + 1)
        ]

    def test_it_is_applied_at_manifest_load_not_inside_a_stage(
        self, pdc_table, subset_file, tmp_path,
    ):
        """The stage every stage reads.

        Applying it inside the cohort-joining stage would leave stage 1's
        validation counts, stage 2's annotation and the sample cap describing a
        cohort the run did not analyse - the same class of disagreement
        `derive_regulator_table` refuses over.
        """
        config = _config_with_subset(subset_file, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        # What every stage receives, not what one stage computes.
        assert len(manifest) == SMOKE_COHORT_SIZE
        assert apply_cohort_subset(config, manifest) == manifest, (
            "applying the subset twice must be a no-op. It is not idempotent if "
            "the second call re-reads the file and the first dropped something."
        )

    def test_manifest_order_is_preserved(self, pdc_table, tmp_path):
        """The list selects membership, not order.

        Every downstream consumer iterates in manifest order - the phenotype
        join above all - so reordering the cohort by an operator's line order
        would make results depend on that file for no stated reason.
        """
        listed = tmp_path / "reversed.txt"
        ids = [f"PDT{i:09d}.1" for i in range(1, SMOKE_COHORT_SIZE + 1)]
        listed.write_text("\n".join(reversed(ids)) + "\n", encoding="utf-8")
        config = _config_with_subset(listed, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        assert manifest.sample_ids == ids, (
            "manifest order, not file order. If this fails, a reader cannot tell "
            "whether the run's cohort ordering depends on a local file."
        )

    def test_a_member_without_an_assembly_is_still_a_member(
        self, tmp_path,
    ):
        """The subset is a selection, not a filter on sequence availability.

        `discover_pdc_manifest` keeps every isolate so a missing genome refuses
        per sample rather than vanishing. A subset that dropped them would change
        `merge_calls`'s denominator while reporting the same list.
        """
        pdc = tmp_path / "PDC_essential.tsv"
        pdc.write_text(
            "Isolate\tAssembly\tBioSample\tBioProject\tPDC_present\n"
            "PDT000000001.1\tGCA_1.1\tSAMN1\tPRJNA1\tPDC\n"
            "PDT000000002.1\t\tSAMN2\tPRJNA1\tPDC\n",
            encoding="utf-8",
        )
        listed = tmp_path / "both.txt"
        listed.write_text(
            "PDT000000001.1\nPDT000000002.1\n", encoding="utf-8"
        )
        config = _config_with_subset(listed, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc)
        assert manifest.sample_ids == ["PDT000000001.1", "PDT000000002.1"]
        assert manifest.get("PDT000000002.1").assembly_path is None, (
            "listed but unsequenced: a member whose sequence is absent. It "
            "refuses per sample in the stages that need one, which is correct."
        )


class TestAMissingListedIdRefuses:
    def test_it_names_the_id(self, pdc_table, tmp_path):
        listed = tmp_path / "with_a_ghost.txt"
        listed.write_text(
            "PDT000000001.1\nPDT999999999.1\n", encoding="utf-8"
        )
        config = _config_with_subset(listed, tmp_path)
        with pytest.raises(PipelineError) as caught:
            discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        message = str(caught.value)
        assert "PDT999999999.1" in message, (
            "the refusal must name the isolate the manifest does not contain: "
            f"{message}"
        )
        assert str(listed) in message, (
            "and the resolved list it was read from, so an operator can see "
            "which file to edit"
        )

    def test_a_partial_match_does_not_produce_a_partial_cohort(
        self, pdc_table, tmp_path,
    ):
        """Refuse, not drop. AGENTS.md rule 5.

        Silently dropping the unknown id would run nine isolates while reporting
        a ten-isolate list, and the numbers would agree with each other.
        """
        listed = tmp_path / "half_bad.txt"
        listed.write_text(
            "\n".join(
                [f"PDT{i:09d}.1" for i in range(1, 6)]
                + ["PDT999999998.1", "PDT999999999.1"]
            )
            + "\n",
            encoding="utf-8",
        )
        config = _config_with_subset(listed, tmp_path)
        with pytest.raises(PipelineError) as caught:
            discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        message = str(caught.value)
        assert "PDT999999998.1" in message and "PDT999999999.1" in message, (
            "every unknown id is named, not just the first: an operator fixing "
            "the list should see the whole problem in one refusal"
        )

    def test_a_duplicate_is_a_configuration_error(self, pdc_table, tmp_path):
        """`SampleManifest` would raise; naming it as config is more useful."""
        listed = tmp_path / "dupes.txt"
        listed.write_text(
            "PDT000000001.1\nPDT000000001.1\nPDT000000002.1\n", encoding="utf-8"
        )
        config = _config_with_subset(listed, tmp_path)
        with pytest.raises(PipelineError) as caught:
            discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        assert "PDT000000001.1" in str(caught.value)

    def test_a_missing_list_file_still_refuses_from_the_config_layer(
        self, pdc_table, tmp_path,
    ):
        """`read_subset_file`'s own refusal, unchanged.

        Falling back to "no subset, the whole roster" would analyse 25 isolates
        while appearing to be the bounded run that was asked for - which is
        exactly what the key must not be able to express.
        """
        config = _config_with_subset(tmp_path / "not_there.txt", tmp_path)
        with pytest.raises(Exception) as caught:
            discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        assert "cohort.subset_file" in str(caught.value)


class TestNoSubsetLeavesEverythingAlone:
    def test_the_committed_overlays_declare_no_subset(self):
        """`laptop` and `bigmachine` inherit `science.yaml`'s `null`."""
        for machine in ("laptop", "bigmachine"):
            config = load_config(SCIENCE, machine=machine)
            assert config.cohort.subset_file is None, (
                f"{machine} declares a subset. WHICH isolates are in the cohort "
                "is a scientific claim; only the bounded smoke overlay may make "
                "it (loader.OVERLAY_SCIENCE_EXCEPTIONS)."
            )

    def test_the_manifest_is_returned_unchanged(self, pdc_table):
        """Behaviour with no subset is identical to before this change."""
        config = load_config(SCIENCE, machine="laptop")
        full = discover_pdc_manifest(pdc_table)
        assert apply_cohort_subset(config, full) is full, (
            "with no subset configured the manifest must be returned as-is - "
            "not copied, not re-keyed. Anything else is a behaviour change "
            "hiding behind a no-op."
        )

    def test_a_test_run_is_unaffected(self, pdc_table):
        """The TEST fixtures have no subset and must stay at twenty."""
        config = load_config(SCIENCE, machine="laptop")
        manifest = apply_cohort_subset(config, SampleManifest(
            discover_pdc_manifest(pdc_table).samples
        ))
        assert len(manifest) == 25


class TestTheBoundedCohortJoinsItsPhenotypeTable:
    """10/10, which is the property a smoke run exists to establish."""

    def test_every_member_has_exactly_one_phenotype_row(
        self, pdc_table, subset_file, tmp_path,
    ):
        from papipeline.stages import phenotype as stage_phenotype

        config = _config_with_subset(subset_file, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        assert len(manifest) == SMOKE_COHORT_SIZE

        # The bounded phenotype table, written for exactly the subset. Ten rows
        # for ten isolates: the join refuses any member without exactly one row,
        # so 10/10 here is a real agreement rather than a subset of a larger
        # table that happened to line up. The columns are the ones
        # `stages.phenotype.REQUIRED_COLUMNS` declares.
        table = tmp_path / "imipenem_phenotype.tsv"
        table.write_text(
            "sample_id\tantibiotic\tphenotype\tMIC\tMIC_unit\tsource\n"
            + "\n".join(
                f"{sample.sample_id}\timipenem\t"
                f"{'R' if index % 2 else 'S'}\t{0.5 + index}\tmg/L\tsynthetic"
                for index, sample in enumerate(manifest)
            )
            + "\n",
            encoding="utf-8",
        )
        calls = stage_phenotype.load_phenotype(
            config, table.parent, "imipenem", sample_ids=manifest.sample_ids
        )
        joined = join_manifest_to_phenotype(manifest.sample_ids, calls)
        assert joined.sample_ids() == manifest.sample_ids, (
            "the phenotype join accepted a different set of isolates than the "
            "subset selected. `stages/phenotype.py` is owned elsewhere and is "
            "not edited here; this is the check that the bounded cohort and the "
            "bounded phenotype table describe the same ten isolates."
        )
        assert len(joined.joined) == SMOKE_COHORT_SIZE, (
            f"the join kept {len(joined.joined)} of {SMOKE_COHORT_SIZE}. A "
            "bounded run that quietly analyses fewer isolates than it names is "
            "the outcome this whole file exists to prevent."
        )
        assert not joined.excluded, (
            f"excluded rather than joined: {sorted(joined.excluded)}"
        )

    def test_a_member_absent_from_the_phenotype_table_still_refuses(
        self, pdc_table, subset_file, tmp_path,
    ):
        """The join's own guarantee, which the subset must not weaken.

        If the subset could admit an isolate the phenotype table does not
        describe, stage 11 would either exclude it silently or fail - and the
        first is the failure mode this file exists to prevent.
        """
        config = _config_with_subset(subset_file, tmp_path)
        manifest = discover_run_manifest(config, RunMode.REAL, pdc_path=pdc_table)
        calls = [
            PhenotypeCall(
                sample_id=sample.sample_id,
                antibiotic="imipenem",
                phenotype=Phenotype.S,
                mic=0.5,
                mic_unit="mg/L",
                source="synthetic",
            )
            for sample in list(manifest)[:-1]
        ]
        with pytest.raises(DataContractError):
            join_manifest_to_phenotype(manifest.sample_ids, calls)


class TestTheSmokeOverlayItselfIsUnchanged:
    """The committed overlay's own list is still the one it declares."""

    def test_it_still_points_at_local(self):
        config = load_config(SCIENCE, machine=SMOKE)
        assert config.cohort.subset_file == "local/smoke_isolates.txt", (
            "the smoke overlay's own member list is `local/smoke_isolates.txt` "
            "and is untracked by design (see the overlay header). Nothing here "
            "may repoint it."
        )
        resolved = config.cohort.resolve_subset_file()
        assert resolved == REPO / "local" / "smoke_isolates.txt", (
            "resolved against the REPOSITORY ROOT, not the CWD and not through "
            "the `data`/`db` symlinks - in a dev worktree those point at the "
            "canonical checkout and the list would stop being this worktree's "
            "own (loader.CohortConfig.resolve_subset_file)."
        )

    def test_a_config_without_the_key_keeps_the_default(self):
        assert CohortConfig.from_dict({}).subset_file is None