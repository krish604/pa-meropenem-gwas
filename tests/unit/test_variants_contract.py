"""The per-isolate variant contract, derived from what bcftools actually emits.

Not designed from an expectation of bcftools - read off a real run of
`bcftools mpileup | bcftools call -mv` on `GCA_000710625.1` aligned to the pinned
PAO1 reference. The columns below are the ones that run produced, and several
would not have been guessed:

- `AC` and `AN` come from `mpileup -a FORMAT/AD` and describe the *sample*, so
  with one isolate they are almost always `2/2`. They are kept because they are
  what makes a multi-isolate merge checkable, not because they are informative
  for a single sample.
- `DP4` is a 4-tuple of ref/alt forward-reverse counts, not a single depth.
- `SGB` is Phred-scaled. Its name is a bcftools implementation detail and it is
  deliberately **not** carried into the contract.
- `MQ0F` is a fraction, not a Phred value.
- QUAL is constant at 30.4183 across every called site. That is expected for
  assembly-vs-reference calling - an assembly is a consensus sequence, so depth
  is 1 per locus - and it means **QUAL carries no per-site information here**.
  It is kept for completeness but no threshold may be derived from it.

Every isolate aligned this way produces DP=1, so these calls are *isolate versus
reference*: they conflate species-wide fixed differences with cohort-informative
polymorphisms, and only a multi-isolate merge can separate the two. That merge
is a separate stage and this contract is its input, not the GWAS's.
"""

from __future__ import annotations

import pytest

from papipeline.execution.contracts import STAGE_TABLES
from papipeline.stages.variants import PER_ISOLATE_COLUMNS, parse_vcf

# One real record from the run above, verbatim.
REAL_RECORD = (
    "NC_002516.2\t1000348\t.\tC\tT\t30.4183\t.\t"
    "DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60\t"
    "GT:PL:DP:AD\t1/1:60,3,0:1:0,1"
)
REAL_HEADER = [
    '##fileformat=VCFv4.2',
    '##contig=<ID=NC_002516.2,length=6264404>',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def _vcf(*records: str) -> str:
    return "\n".join(REAL_HEADER + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta1.sorted.bam"] + list(records))


class TestPerIsolateContractIsTheRealOne:
    def test_columns_come_from_the_real_run(self):
        """Not the PROVISIONAL placeholders."""
        assert "AC" in PER_ISOLATE_COLUMNS
        assert "AN" in PER_ISOLATE_COLUMNS
        assert "DP4" in PER_ISOLATE_COLUMNS
        assert "MQ" in PER_ISOLATE_COLUMNS

    def test_sgb_is_not_carried(self):
        """A bcftools implementation detail stays out of our contract.

        Asserted on the emitted row, not only the column tuple: the real VCF has
        SGB in INFO, so a parser that lifted every INFO field would carry it even
        with the column list unchanged. Checking only the tuple passes that.
        """
        assert "SGB" not in PER_ISOLATE_COLUMNS
        row = parse_vcf(_vcf(REAL_RECORD), sample_id="x")[0]
        assert "SGB" not in row
        assert "SGB" not in "".join(str(v) for v in row.values())

    def test_required_join_and_locus_columns_are_present(self):
        for column in ("sample_id", "chrom", "pos", "ref", "alt"):
            assert column in PER_ISOLATE_COLUMNS

    def test_the_provisional_contract_has_been_replaced(self):
        """`STAGE_TABLES` must no longer advertise the placeholder columns.

        `STAGE_TABLES[stage]` is `(filename, columns)`, so this pins the columns
        tuple directly. The old contract's `position`/`reference`/`alternate`/
        `variant_type`/`call_status` were guesses that no bcftools run produced.
        """
        _filename, columns = STAGE_TABLES["variants"]
        assert tuple(columns) == PER_ISOLATE_COLUMNS
        for placeholder in ("position", "reference", "alternate", "variant_type", "call_status"):
            assert placeholder not in columns


class TestParsingRealOutput:
    def test_a_real_record_parses(self):
        rows = parse_vcf(_vcf(REAL_RECORD), sample_id="GCA_000710625.1")
        assert len(rows) == 1
        row = rows[0]
        assert row["sample_id"] == "GCA_000710625.1"
        assert row["chrom"] == "NC_002516.2"
        assert row["pos"] == "1000348"
        assert row["ref"] == "C"
        assert row["alt"] == "T"
        assert row["AC"] == "2"
        assert row["AN"] == "2"
        assert row["MQ"] == "60"
        assert row["DP4"] == "0,0,0,1"

    def test_dp4_stays_a_four_tuple_not_a_collapsed_number(self):
        """Collapsing `0,0,0,1` to a scalar loses the ref/alt split."""
        row = parse_vcf(_vcf(REAL_RECORD), sample_id="x")[0]
        assert row["DP4"].count(",") == 3

    def test_mq0f_is_a_fraction_not_a_phred(self):
        """`MQ0F=0` means a proportion. Reading it as Phred would be nonsense."""
        row = parse_vcf(_vcf(REAL_RECORD), sample_id="x")[0]
        assert row["MQ0F"] == "0"

    def test_genotype_is_captured(self):
        row = parse_vcf(_vcf(REAL_RECORD), sample_id="x")[0]
        assert row["GT"] == "1/1"

    def test_multiallelic_sites_are_kept(self):
        """1623 C, 1440 G, 1323 T, 1099 A were all ALT alleles in the real run.

        A comma-separated ALT is two rows' worth of information in one line;
        dropping the second would silently lose half the calls.
        """
        record = REAL_RECORD.replace("\tC\tT\t", "\tC\tT,G\t")
        rows = parse_vcf(_vcf(record), sample_id="x")
        assert len(rows) == 2
        assert [r["alt"] for r in rows] == ["T", "G"]


class TestRefusals:
    def test_no_sample_id_is_refused(self):
        """The isolate is the join key; an unattributable call is not usable."""
        with pytest.raises(Exception):
            parse_vcf(_vcf(REAL_RECORD), sample_id=None)

    def test_header_only_vcf_yields_no_rows(self):
        assert parse_vcf("\n".join(REAL_HEADER), sample_id="x") == []

    def test_a_record_without_pos_is_refused(self):
        bad = REAL_RECORD.replace("\t1000348\t", "\t.\t")
        with pytest.raises(Exception):
            parse_vcf(_vcf(bad), sample_id="x")
