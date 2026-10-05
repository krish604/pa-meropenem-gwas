"""What each pangenome bucket field MEANS, pinned to a real cohort's numbers.

**Why this file exists.** `docs/environment-arm64.md:206` labels the number
5,191 as "accessory (2-9 genomes)". It is not. On the verified 10-isolate smoke
cohort the truth is:

    Core genes        (99% <= strains <= 100%)    4,828
    Soft core genes   (95% <= strains <  99%)        0
    Shell genes       (15% <= strains <  95%)    2,850
    Cloud genes       ( 0% <= strains <  15%)    2,341
    Total genes                                10,019

Read straight off panaroo 1.8.0's own `summary_statistics.txt` after a run on
the ten existing Bakta GFF3 (no inputs regenerated). 2,850 + 2,341 = 5,191, and
10,019 - 4,828 = 5,191. So 5,191 is **shell + cloud** and the "2-9 genomes"
parenthetical is wrong by the cloud's 2,341; the 2-9 figure alone is 2,850.

`n_accessory_genes` is emitted as `total - core` (`stages/pangenome.py`, in
`partition` and at the `write_outputs` metric), which is the same 5,191. **The
code was right and the prose was wrong** - and nothing tested the arithmetic, so
the two could have swapped without any test failing. These tests close that.

**Two things are deliberately NOT asserted here**, because asserting them would be
asserting a decision nobody has made:

- That `n_accessory_genes` *should* be shell-only. It is not changed; it stays
  `total - core`. Only its meaning is now pinned.
- That our core should equal panaroo's. Our core is exactly 100% of isolates;
  panaroo's is >= 99%. They coincide at n=10 and diverge at n=967. See
  `CORE_PRESENCE_THRESHOLD`.

Note this does not contradict `test_pangenome_real_caller.py`'s "no panaroo
counts are asserted here". That file is about the REAL wiring - gate order,
preflight order, cohort matching - where the cohort is a fixture. This file is
about what the bucket *labels* mean, and for that the real cohort is the only
evidence that settles it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import panaroo as adapter
from papipeline.errors import DataContractError
from papipeline.stages import pangenome as stage

# --- panaroo 1.8.0's own numbers, verbatim from its summary_statistics.txt ---
PANAROO_CORE = 4828
PANAROO_SOFT_CORE = 0
PANAROO_SHELL = 2850
PANAROO_CLOUD = 2341
PANAROO_TOTAL = 10019

#: total - core. The number `n_accessory_genes` emits.
NON_CORE = PANAROO_TOTAL - PANAROO_CORE  # == 5191
assert NON_CORE == PANAROO_SOFT_CORE + PANAROO_SHELL + PANAROO_CLOUD

COHORT = tuple(
    f"PDT{n:09d}.1"
    for n in (34122, 167133, 167135, 167136, 292995, 292998, 294804, 294805, 311294, 424983)
)
N = len(COHORT)

#: The real per-gene carrier-count histogram from that run, so the partition can
#: be exercised on the actual distribution rather than a tidy synthetic one.
#: count -> number of genes carried by exactly that many of the ten isolates.
#: 10 -> 4,828 core. 9..2 -> 2,850 shell. 1 -> 2,341 cloud.
REAL_HISTOGRAM = {10: 4828, 9: 390, 8: 223, 7: 84, 6: 138, 5: 197, 4: 325, 3: 412, 2: 1081, 1: 2341}


def _presence_from_histogram(histogram):
    """A presence map with exactly the real carrier-count distribution."""
    presence = {}
    for carriers, n_genes in histogram.items():
        for i in range(n_genes):
            presence[f"group_{carriers}_{i}"] = {COHORT[j] for j in range(carriers)}
    return presence


@pytest.fixture()
def real_presence():
    presence = _presence_from_histogram(REAL_HISTOGRAM)
    assert len(presence) == PANAROO_TOTAL
    return presence


def _summary(pangenome) -> dict:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        written = stage.write_outputs(pangenome, Path(tmp))
        lines = Path(written["pangenome_summary"]).read_text().splitlines()[1:]
    return dict(line.split("\t") for line in lines if line.strip())


class TestTheRealCohortsNumbers:
    def test_the_partition_reproduces_panaroos_core_count(self, real_presence):
        p = stage.partition(COHORT, real_presence)
        assert len(p.genes) == PANAROO_TOTAL
        assert len(p.core) == PANAROO_CORE

    def test_the_partition_reproduces_panaroos_non_core_count(self, real_presence):
        """`accessory` is total - core, so 5,191 - and 5,191 is NOT the shell."""
        p = stage.partition(COHORT, real_presence)
        assert len(p.accessory) == NON_CORE == 5191

    def test_accessory_is_not_the_shell_bucket(self, real_presence):
        """The specific trap. Accessory (5,191) > shell (2,850) by the cloud.

        This is the assertion that fails if someone "fixes" the doc by changing
        the code to shell-only, or vice versa. It pins that the two are different
        quantities that happen to be conflated in one line of prose.
        """
        p = stage.partition(COHORT, real_presence)
        shell = sum(1 for g in p.genes if 1 < len(p.presence[g]) < N)
        cloud = sum(1 for g in p.genes if len(p.presence[g]) == 1)

        assert shell == PANAROO_SHELL == 2850
        assert cloud == PANAROO_CLOUD == 2341
        assert len(p.accessory) == shell + cloud
        assert len(p.accessory) != shell
        # The "2-9 genomes" figure the doc prints next to 5,191.
        assert shell == 2850


class TestWhichEmittedFieldMeansWhat:
    """Read off the published `pangenome_summary.tsv`, not off the object."""

    def test_n_core_genes_is_exactly_100_percent_presence(self, real_presence):
        summary = _summary(stage.partition(COHORT, real_presence))
        assert summary["n_core_genes"] == "4828"

    def test_n_accessory_genes_is_total_minus_core(self, real_presence):
        summary = _summary(stage.partition(COHORT, real_presence))
        assert summary["n_accessory_genes"] == "5191"
        assert int(summary["n_accessory_genes"]) == (
            int(summary["n_genes_total"]) - int(summary["n_core_genes"])
        )

    def test_the_three_numbers_partition_the_total(self, real_presence):
        summary = _summary(stage.partition(COHORT, real_presence))
        assert summary["n_genes_total"] == "10019"
        assert int(summary["n_core_genes"]) + int(summary["n_accessory_genes"]) == int(
            summary["n_genes_total"]
        )

    def test_core_fraction_is_core_over_total_not_over_cohort(self, real_presence):
        summary = _summary(stage.partition(COHORT, real_presence))
        assert float(summary["core_fraction"]) == round(PANAROO_CORE / PANAROO_TOTAL, 6)


class TestOurCoreIsExactlyOneHundredPercentNotNinetyNine:
    def test_a_gene_in_ninety_nine_percent_is_not_core(self):
        """n=10: panaroo's Core is >=99%, which at n=10 is also 10/10.

        A 9-of-10 gene is soft-core/pan-shell for panaroo and simply non-core
        here. Both put it outside the core, so this case does not separate the
        definitions - it pins that we never quietly widened the threshold.
        """
        presence = {"in_all": set(COHORT), "in_nine": set(COHORT[:9])}
        p = stage.partition(COHORT, presence)
        assert p.core == ("in_all",)
        assert p.accessory == ("in_nine",)
        assert stage.CORE_PRESENCE_THRESHOLD == 1.0

    def test_the_threshold_is_one_and_not_a_panaroo_like_default(self):
        """A 0.95 default would have matched panaroo's soft-core line.

        Pinned so that anyone aligning this stage with panaroo has to change the
        constant deliberately, in a commit that also says which definition won.
        """
        assert stage.CORE_PRESENCE_THRESHOLD == 1.0


class TestPanarooSummaryParser:
    """The five distinct lines, on the real file's exact text."""

    REAL = (
        "Core genes\t(99% <= strains <= 100%)\t4828\n"
        "Soft core genes\t(95% <= strains < 99%)\t0\n"
        "Shell genes\t(15% <= strains < 95%)\t2850\n"
        "Cloud genes\t(0% <= strains < 15%)\t2341\n"
        "Total genes\t(0% <= strains <= 100%)\t10019\n"
    )

    def _write(self, tmp_path: Path, text: str) -> Path:
        path = tmp_path / "summary_statistics.txt"
        path.write_text(text, encoding="utf-8")
        return path

    def test_it_reads_all_five_buckets(self, tmp_path):
        b = adapter.parse_summary_statistics(self._write(tmp_path, self.REAL))
        assert b.as_dict() == {
            "core": 4828,
            "soft_core": 0,
            "shell": 2850,
            "cloud": 2341,
            "total": 10019,
        }

    def test_soft_core_is_not_read_as_core(self, tmp_path):
        """`Core genes` and `Soft core genes` overlap as prefixes.

        The counts are chosen so a prefix match or a `in`/substring match reads
        one bucket's count as the other's: 2,850 vs 1,978, far enough apart that
        a conflation cannot accidentally agree. Soft-core is 0 in the real file
        and core is 4,828, so the naive failure here is core reported as 0 - a
        working 4,828-gene core looking like a catastrophic clustering failure.

        The figures are kept self-consistent (they still sum to the total) so
        this exercises label matching rather than the add-up guard.
        """
        text = (
            "Core genes\t(99% <= strains <= 100%)\t2850\n"
            "Soft core genes\t(95% <= strains < 99%)\t1978\n"
            "Shell genes\t(15% <= strains < 95%)\t2850\n"
            "Cloud genes\t(0% <= strains < 15%)\t2341\n"
            "Total genes\t(0% <= strains <= 100%)\t10019\n"
        )
        b = adapter.parse_summary_statistics(self._write(tmp_path, text))
        assert b.core == 2850
        assert b.soft_core == 1978

    def test_a_real_summary_with_zero_soft_core_does_not_report_zero_core(self, tmp_path):
        """The n=10 case, where soft-core really is 0 and core really is 4,828."""
        b = adapter.parse_summary_statistics(self._write(tmp_path, self.REAL))
        assert b.soft_core == 0
        assert b.core == 4828

    def test_non_core_is_shell_plus_cloud(self, tmp_path):
        b = adapter.parse_summary_statistics(self._write(tmp_path, self.REAL))
        assert b.non_core == 5191
        assert b.non_core == b.shell + b.cloud + b.soft_core
        assert b.non_core != b.shell

    def test_our_core_agrees_with_panaroos_at_n_ten(self, tmp_path):
        """At n=10 the >=99% and ==100% definitions coincide, so this passes.

        The test is worth having precisely because it is the only cohort size at
        which this agreement is automatic. It is not evidence about n=967.
        """
        b = adapter.parse_summary_statistics(self._write(tmp_path, self.REAL))
        presence = _presence_from_histogram(REAL_HISTOGRAM)
        p = stage.partition(COHORT, presence)
        assert len(p.core) == b.core == 4828
        assert len(p.accessory) == b.non_core == 5191

    def test_a_missing_file_is_named(self, tmp_path):
        with pytest.raises(DataContractError, match="summary_statistics"):
            adapter.parse_summary_statistics(tmp_path / "summary_statistics.txt")

    def test_a_missing_bucket_is_a_failure_not_a_zero(self, tmp_path):
        """Rule 9: missing stays missing. A dropped line is not `soft_core=0`."""
        text = "".join(
            line + "\n"
            for line in self.REAL.splitlines()
            if not line.startswith("Shell genes")
        )
        with pytest.raises(DataContractError, match="shell"):
            adapter.parse_summary_statistics(self._write(tmp_path, text))

    def test_a_summary_that_does_not_add_up_is_refused(self, tmp_path):
        text = self.REAL.replace(
            "Cloud genes\t(0% <= strains < 15%)\t2341",
            "Cloud genes\t(0% <= strains < 15%)\t9999",
        )
        with pytest.raises(DataContractError, match="does not add up"):
            adapter.parse_summary_statistics(self._write(tmp_path, text))

    def test_an_unknown_label_does_not_silently_become_a_bucket(self, tmp_path):
        text = self.REAL.replace("Shell genes\t(15% <= strains < 95%)\t2850",
                                 "Shellish genes\t(15% <= strains < 95%)\t2850")
        with pytest.raises(DataContractError, match="shell"):
            adapter.parse_summary_statistics(self._write(tmp_path, text))
