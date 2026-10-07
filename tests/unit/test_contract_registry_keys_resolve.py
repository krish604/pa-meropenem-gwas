"""`DENSE_PER_SAMPLE_STAGES` and `PER_SAMPLE_STAGES` are keyed by table, not by module.

Both sets are consulted as *names of a table's grain*:

* `contracts.stage_spec` gates its opt-in coverage check on
  `stage in DENSE_PER_SAMPLE_STAGES`, and builds the path with `table_path`,
  which resolves against `STAGE_TABLES`;
* `PER_SAMPLE_STAGES` is documented as "the coarser `has a sample_id column`
  notion", which is a property of a declared header, so its members have to be
  names that carry one.

A member that is neither a `STAGE_TABLES` key nor an `INTERNAL_TABLES` key
resolves to nothing: it is a claim about a table no registry declares, so it
can be consulted and silently mean "unknown" instead of "no". `integration`
was exactly that — the folded step that *builds* `15_master_table.tsv`, entered
under the module's name while every other folded step (`regulators`,
`structural_variants`, `mechanisms`) is entered under the file it writes.

Grain wording lives in `docs/data_contract.md`, not here: this file only pins
that every member names a table the pipeline actually declares.
"""

from __future__ import annotations

from papipeline.execution.contracts import (
    DENSE_PER_SAMPLE_STAGES,
    INTERNAL_TABLES,
    PER_SAMPLE_STAGES,
    STAGE_TABLES,
    internal_table_path,
    table_path,
)

#: Every table the pipeline declares, across both registries.
DECLARED = set(STAGE_TABLES) | set(INTERNAL_TABLES)


class TestEveryMemberNamesADeclaredTable:
    def test_dense_members_all_resolve(self):
        unknown = sorted(DENSE_PER_SAMPLE_STAGES - DECLARED)
        assert not unknown, (
            "DENSE_PER_SAMPLE_STAGES names a table no registry declares: "
            f"{unknown}. A member that cannot be resolved to STAGE_TABLES or "
            "INTERNAL_TABLES is a grain claim about a file that has no header "
            "to check, and `stage_spec` would resolve it to nothing."
        )

    def test_per_sample_members_all_resolve(self):
        unknown = sorted(PER_SAMPLE_STAGES - DECLARED)
        assert not unknown, (
            "PER_SAMPLE_STAGES names a table no registry declares: "
            f"{unknown}. The set means 'has a sample_id column', which is a "
            "property of a declared header - so the member must be a declared "
            "table name."
        )

    def test_a_dense_member_has_a_path_that_can_be_built(self, tmp_path):
        """Resolvable in the sense that matters: the path constructor agrees."""
        for stage in sorted(DENSE_PER_SAMPLE_STAGES):
            if stage in STAGE_TABLES:
                path = table_path(tmp_path, stage)
            else:
                path = internal_table_path(tmp_path, stage)
            assert path.name.endswith(".tsv"), (
                f"{stage!r} resolved to a non-TSV path: {path}"
            )


class TestTheMasterTableIsNamedByItsFile:
    """`integration` folded into reporting; its table did not follow it in.

    spec.md:351 makes the step internal to `reporting`, so nothing schedules
    "integration" as a stage — but the table it writes is declared, and it is
    the densest table in the pipeline: one row per sample, every sample, with
    `sample_id` leading (`docs/data_contract.md`: "15_master_table.tsv | one
    row per sample"). Claiming the grain under a module name left the one
    genuinely per-sample joined table out of the set that means per-sample.
    """

    def test_the_joined_table_carries_the_grain(self):
        assert "master_table" in DENSE_PER_SAMPLE_STAGES, (
            "15_master_table.tsv is one row per sample; it is the table the "
            "coverage notion was written for"
        )

    def test_the_folded_module_name_is_not_the_key(self):
        for registry_name, registry in (
            ("DENSE_PER_SAMPLE_STAGES", DENSE_PER_SAMPLE_STAGES),
            ("PER_SAMPLE_STAGES", PER_SAMPLE_STAGES),
        ):
            assert "integration" not in registry, (
                f"{registry_name} contains 'integration', which is a folded "
                "step (spec.md:351), not a declared table. Key it by the file "
                "it writes — `master_table` — like every other folded step."
            )

    def test_the_other_folded_steps_are_keyed_by_their_files(self):
        """The precedent `integration` was the exception to."""
        for folded in ("mechanisms", "regulators", "structural_variants"):
            assert folded in INTERNAL_TABLES, f"{folded} has no declared file"
            assert folded in PER_SAMPLE_STAGES, (
                f"{folded} writes a table with sample_id, so the coarse set "
                "should still name it"
            )


class TestTheCoarseSetKeepsTheStepsWhoseGrainIsNotOneRowPerSample:
    """`05_mechanisms.tsv` is one row per *mechanism call*.

    `docs/data_contract.md` states the grain; the coarse set only claims a
    `sample_id` column, which the mechanisms table has. So it stays in
    `PER_SAMPLE_STAGES` even though it must not be described as one row per
    sample — the two claims are different, and conflating them is how a
    coverage check comes to assert something false about a sample with three
    determinant mechanisms.
    """

    def test_mechanisms_is_coarse_but_not_dense(self):
        assert "mechanisms" in PER_SAMPLE_STAGES, (
            "05_mechanisms.tsv carries sample_id, so the coarse set names it"
        )

    def test_virulence_is_dense_by_the_opt_in_not_by_row_count(self):
        """`08_virulence.tsv` is one row per virulence detection.

        Kept in the dense set because the coverage check counts *distinct*
        `sample_id` values, not rows — a sample with three virulence factors
        contributes three rows and still covers itself once. The grain is
        recorded here so nobody reads dense membership as "exactly one row".
        """
        assert "virulence" in DENSE_PER_SAMPLE_STAGES
        assert "virulence" in PER_SAMPLE_STAGES
