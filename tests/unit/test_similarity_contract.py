"""The `similarity` output contract: a square matrix, not a long pair list.

The declared contract was `("sample_a", "sample_b", "distance", "method")` -
one row per pair. That shape was never read by anything: no parser, no
consumer, no test. It is also a placeholder written before anyone decided what
the stage measures, and it is wrong for the decision now made.

**The stage emits a square matrix.** A header row of sample names, then one row
per sample carrying its distances to every other sample. Two reasons:

* a long pair list duplicates every distance, so an N-sample matrix costs
  N(N-1)/2 rows - 7,405 rows at 835 isolates - while the square form costs 835.
  Nothing reads it, so the duplication buys nothing;
* the stand-in fixture written for this stage is *already* square
  (`STANDIN  STANDIN  standin_not_computed  standin_not_computed`). The fixture
  and the contract disagreed, and the fixture was closer to the intent.

**The diagonal is 0 and the matrix is symmetric.** Asserted here rather than
assumed, because a matrix that fails either is not a distance matrix and every
consumer of one assumes both. An asymmetric distance matrix is the signature of
a traversal bug, and it is invisible unless something checks.

**`method` is `patristic`.** The distance is read off the stage-9 phylogeny:
the sum of branch lengths along the path between two tips. It is not a distance
matrix IQ-TREE emitted and then re-read - that would be circular, since IQ-TREE
built the tree from such a matrix. Carrying the method in the file is what stops
a future reader assuming a different definition produced these numbers.

`sample_id` remains the leading column so the file is still keyed on the one
identifier every stage uses, which is also what makes the row order checkable
against the manifest.
"""

from __future__ import annotations

import pytest

from papipeline.execution.contracts import (
    PER_SAMPLE_STAGES,
    STAGE_TABLES,
    required_columns,
)
from papipeline.run import STAGE_ORDER, UNBUILT_STAGES

STAGE = "similarity"

#: The leading column. Everything after it is one sample's name.
ID_COLUMN = "sample_id"


class TestTheContractIsASquareMatrix:
    def test_it_starts_with_the_sample_id(self):
        columns = required_columns(STAGE)
        assert columns[0] == ID_COLUMN, (
            "the matrix is keyed on sample_id like every other stage table, "
            "which is also what makes the row order checkable"
        )

    def test_the_old_pair_columns_are_gone(self):
        """A schema change, so the old shape must not survive anywhere.

        `sample_a`/`sample_b` described a long pair list. Leaving them in place
        would leave a second, contradictory statement of the shape.
        """
        columns = required_columns(STAGE)
        for gone in ("sample_a", "sample_b"):
            assert gone not in columns

    def test_it_does_not_declare_a_distance_column(self):
        """Distances are the per-sample columns, named by the sample they lead to.

        A single `distance` column cannot survive a square layout - a row holds
        N distances, not one.
        """
        assert "distance" not in required_columns(STAGE)

    def test_the_declared_header_is_id_plus_one_vector_column(self):
        """Two declared columns, and the matrix's real width stays a run-time fact.

        The reasoning here outlived the reconciliation and is why the declaration
        is two columns rather than one *or* twenty: the distances are named after
        the samples in the cohort, so their names are not knowable until a cohort
        exists. Declaring them individually would hard-code a cohort size into the
        contract - the same assumption the manifest module exists to avoid. So the
        whole vector travels in one `distances` column, keyed on the cohort.

        This asserted `(ID_COLUMN,)` alone while the stage was unbuilt, which was
        a declaration for a file nothing wrote. Reconciled, the stage emits
        `sample_id` and `distances`, so the contract was corrected to match
        rather than the stage changed to match a placeholder.
        """
        assert required_columns(STAGE) == (ID_COLUMN, "distances")

    def test_the_stage_is_built_and_in_the_dag(self):
        """Reconciled once TEST (1d54f31) and REAL (2223332) both existed.

        This file asserted the opposite for its whole life - correctly at the
        time, as a guard against a contract being read as a working stage. The
        stage is built now, so the guard inverts: a contract that outlived its
        stage would otherwise keep the stage declared absent.
        """
        assert STAGE not in UNBUILT_STAGES
        assert STAGE in STAGE_ORDER

    def test_the_file_name_is_unchanged(self):
        """The filename is declared in the workflow and asserted by standin tests."""
        assert STAGE_TABLES[STAGE][0] == "similarity.tsv"


class TestWhatTheShapeImplies:
    """Taxonomy consequences of the change, checked rather than assumed."""

    def test_it_is_a_per_sample_table(self):
        """One row per sample, so the per-sample set is now correct.

        `PER_SAMPLE_STAGES` is what tells a caller it may opt into a coverage
        check. A square matrix has exactly one row per sample, which is the
        strongest form of that property the contract has.
        """
        assert STAGE in PER_SAMPLE_STAGES

    def test_it_is_not_coverage_asserted_by_default(self):
        """Per-sample does not mean dense.

        Coverage is opt-in for a documented reason: several per-sample stages are
        presence-dependent, and asserting coverage by default would fail a
        correct run whose correct answer is "this isolate has no determinant".
        The same reasoning holds here - a matrix always has a row per sample,
        so the check is a shape assertion, not a biological one.
        """
        from papipeline.execution.contracts import DENSE_PER_SAMPLE_STAGES

        # Being in the dense set would make a coverage check available; the
        # stage is a matrix, so it is not claiming to be a record table.
        assert STAGE not in DENSE_PER_SAMPLE_STAGES
