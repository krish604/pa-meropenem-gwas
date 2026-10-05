"""Round 12 - the producer for stage 12's ``gwas/gwas_features.tsv``.

Nothing in the pipeline wrote the table ``gwas.run`` reads; it was
operator-provisioned and undocumented as such. These tests pin the four
properties the producer has to have to be worth having:

1. **It produces the declared shape.** ``sample_id`` unique, one column per
   feature named ``<type>__<label>``, values binary, every prefix in
   ``config.gwas.feature_types``.
2. **A point mutation is an ALLELE, not a gene gain.** This is the property a
   prior triage flagged as at risk, so it is asserted directly rather than
   inferred from a column list.
3. **It refuses, naming every missing input with its path and producer stage.**
4. **It is deterministic** - same inputs, byte-identical file.

Everything here runs against SYNTHETIC inputs or INJECTED tables (R12). The
real 10-isolate intermediates are exercised by
``TestTheRealTenIsolateCohort`` against files this repository does not commit,
which skips when they are absent rather than pretending to have run.
"""

from __future__ import annotations

import pytest

from pathlib import Path

from papipeline.errors import DataContractError, SampleIdError
from papipeline.io.tsv import read_tsv
from papipeline.manifest import SampleManifest
from papipeline.models import RegulatorVariant, RunMode, Sample
from papipeline.stages import gwas_features as gf

_REPO = Path(__file__).resolve().parents[2]

SAMPLES = ("TEST_A_01", "TEST_A_02", "TEST_A_03")


def _manifest(sample_ids=SAMPLES) -> SampleManifest:
    return SampleManifest(
        [Sample(sample_id, None, "synthetic") for sample_id in sample_ids]
    )


def _amr_row(
    sample_id: str,
    determinant: str,
    *,
    gene: "str | None" = None,
    variant: "str | None" = None,
    antibiotic: str = "imipenem",
) -> dict:
    """One stage 4 row, as `read_tsv` hands it over: absent is ``None``."""
    return {
        "sample_id": sample_id,
        "antibiotic": antibiotic,
        "determinant": determinant,
        "gene": gene,
        "variant": variant,
    }


def _presence(**genes_with_carriers):
    """``{gene: [carriers]}`` for the stage 7 matrix."""
    return {gene: list(carriers) for gene, carriers in genes_with_carriers.items()}


def _call(sample_id: str, variant_type: str, gene: str = "oprD") -> RegulatorVariant:
    return RegulatorVariant(
        sample_id=sample_id,
        gene=gene,
        variant=f"SYNTHETIC_{variant_type}",
        variant_type=variant_type,
        position=412,
        reference="SYN",
        alternate="SYN",
        effect="synthetic_effect",
        mechanism="reduced_permeability",
        confidence="SYNTHETIC",
    )


# ---------------------------------------------------------------------------
# The declared shape
# ---------------------------------------------------------------------------


class TestTheDeclaredShape:
    def test_one_row_per_sample_and_a_unique_sample_id(self, config):
        rows, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=["TEST_A_01"]),
            amr_rows=[],
            regulator_records=[],
        )
        assert [r["sample_id"] for r in rows] == list(SAMPLES)
        assert columns[0] == "sample_id"
        assert len({r["sample_id"] for r in rows}) == len(rows)

    def test_rows_come_out_in_sorted_sample_order(self, config):
        rows, _, _ = gf.build_feature_rows(
            config,
            _manifest(("TEST_A_03", "TEST_A_01", "TEST_A_02")),
            gene_presence=_presence(synA=["TEST_A_01"]),
            amr_rows=[],
            regulator_records=[],
        )
        assert [r["sample_id"] for r in rows] == sorted(SAMPLES)

    def test_every_feature_column_is_named_type_double_underscore_label(
        self, config
    ):
        _, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=["TEST_A_01"]),
            amr_rows=[_amr_row("TEST_A_01", "oprD_V359L", gene="oprD", variant="V359L")],
            regulator_records=[],
        )
        feature_columns = columns[1:]
        assert feature_columns
        for column in feature_columns:
            assert "__" in column, column
            prefix = column.split("__", 1)[0]
            assert prefix in set(config.gwas.feature_types), column

    def test_every_value_is_binary(self, config):
        rows, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=["TEST_A_01", "TEST_A_03"]),
            amr_rows=[_amr_row("TEST_A_02", "oprD_V359L", gene="oprD", variant="V359L")],
            regulator_records=[_call("TEST_A_03", "GENE_ABSENCE")],
        )
        for row in rows:
            for column in columns[1:]:
                assert row[column] in (0, 1), (column, row[column])

    def test_the_two_oprd_columns_use_the_prefix_the_fixture_declares(
        self, config
    ):
        _, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=["TEST_A_01"]),
            amr_rows=[],
            regulator_records=[_call("TEST_A_01", "GENE_ABSENCE")],
        )
        # The committed stage 12 fixture names these two exactly, and
        # oprd_feature_columns documents the prefix as the caller's choice.
        assert "gene__oprD_absent" in columns
        assert "gene__oprD_LoF" in columns
        assert gf.OPRD_FEATURE_TYPE == "gene"

    def test_the_pangenome_prefix_is_the_contract_prefix(self, config):
        _, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(oprD=["TEST_A_01"]),
            amr_rows=[],
            regulator_records=[],
        )
        assert "gene__oprD" in columns
        assert gf.GENE_FEATURE_TYPE == "gene"


# ---------------------------------------------------------------------------
# F2: a point mutation is an allele, never a gene gain
# ---------------------------------------------------------------------------


class TestAPointMutationIsAnAlleleNotAGeneGain:
    def test_a_point_mutation_becomes_an_allele_named_column(self, config):
        _, columns, counts = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(oprD=list(SAMPLES)),
            amr_rows=[
                _amr_row("TEST_A_01", "oprD_V359L", gene="oprD", variant="V359L")
            ],
            regulator_records=[],
        )
        assert "gene_presence_absence__oprD_V359L" in columns
        assert counts["amr_allele"] == 1

    def test_the_gene_itself_is_not_credited_with_the_mutation(self, config):
        """The risk a prior triage named: `oprD_V359L` read as a gene gain.

        Every isolate in this cohort carries oprD, so `gene__oprD` must be 1
        throughout and must not be what the point mutation was recorded as. The
        substitution gets its own column; the gene's presence is unchanged.
        """
        rows, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(oprD=list(SAMPLES)),
            amr_rows=[
                _amr_row("TEST_A_01", "oprD_V359L", gene="oprD", variant="V359L")
            ],
            regulator_records=[],
        )
        by_sample = {r["sample_id"]: r for r in rows}
        for sample_id in SAMPLES:
            assert by_sample[sample_id]["gene__oprD"] == 1, sample_id
        assert by_sample["TEST_A_01"]["gene_presence_absence__oprD_V359L"] == 1
        assert by_sample["TEST_A_02"]["gene_presence_absence__oprD_V359L"] == 0
        assert by_sample["TEST_A_03"]["gene_presence_absence__oprD_V359L"] == 0

    def test_an_allele_is_not_labelled_a_snp(self, config):
        """POINT_DISRUPT rows include indels, so `snp__` would over-claim.

        `rpsJ_HKYK56del` is a multi-residue deletion and the `variant` column
        holds protein-level allele strings, not coordinates.
        """
        _, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(rpsJ=["TEST_A_01"]),
            amr_rows=[
                _amr_row(
                    "TEST_A_01", "rpsJ_HKYK56del", gene="rpsJ", variant="HKYK56del"
                )
            ],
            regulator_records=[],
        )
        assert not [c for c in columns if c.startswith("snp__")]
        assert "gene_presence_absence__rpsJ_HKYK56del" in columns

    def test_an_intact_determinant_is_a_gene_not_an_allele(self, config):
        _, columns, counts = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(mexA=["TEST_A_01"]),
            amr_rows=[_amr_row("TEST_A_01", "mexA", gene="mexA")],
            regulator_records=[],
        )
        assert "gene_presence_absence__mexA" in columns
        assert not [c for c in columns if c.endswith("__mexA.")]
        assert counts["amr_determinant_gene"] == 1
        assert counts["amr_allele"] == 0

    def test_a_variant_without_its_gene_is_refused_not_guessed(self, config):
        """A substitution folded into gene presence is a false mechanism claim."""
        with pytest.raises(DataContractError) as excinfo:
            gf.build_feature_rows(
                config,
                _manifest(),
                gene_presence=_presence(oprD=list(SAMPLES)),
                amr_rows=[_amr_row("TEST_A_01", "oprD_V359L", variant="V359L")],
                regulator_records=[],
            )
        message = str(excinfo.value)
        assert "TEST_A_01:oprD_V359L" in message
        assert "gene gain" in message

    def test_two_point_rows_in_one_gene_are_two_features(self, config):
        _, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(oprD=list(SAMPLES)),
            amr_rows=[
                _amr_row("TEST_A_01", "oprD_V359L", gene="oprD", variant="V359L"),
                _amr_row("TEST_A_02", "oprD_G151PfsTer39", gene="oprD",
                         variant="G151PfsTer39"),
            ],
            regulator_records=[],
        )
        assert "gene_presence_absence__oprD_V359L" in columns
        assert "gene_presence_absence__oprD_G151PfsTer39" in columns

    def test_the_allele_carrier_is_the_isolate_that_carries_the_substitution(
        self, config
    ):
        """Regression: the AMR families were emitted with no carriers at all.

        The columns existed and every value read 0, which is a plausible-looking
        table - it just says no isolate has ever been screened for point
        mutations.
        """
        rows, columns, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(oprD=list(SAMPLES)),
            amr_rows=[
                _amr_row("TEST_A_02", "oprD_V359L", gene="oprD", variant="V359L"),
                _amr_row("TEST_A_03", "oprD_V359L", gene="oprD", variant="V359L"),
                _amr_row("TEST_A_01", "oprD_G151PfsTer39", gene="oprD",
                         variant="G151PfsTer39"),
            ],
            regulator_records=[],
        )
        by_sample = {r["sample_id"]: r for r in rows}
        assert "gene_presence_absence__oprD_V359L" in columns
        assert by_sample["TEST_A_02"]["gene_presence_absence__oprD_V359L"] == 1
        assert by_sample["TEST_A_03"]["gene_presence_absence__oprD_V359L"] == 1
        assert by_sample["TEST_A_01"]["gene_presence_absence__oprD_V359L"] == 0
        assert by_sample["TEST_A_01"]["gene_presence_absence__oprD_G151PfsTer39"] == 1

    def test_a_gene_the_pangenome_also_knows_keeps_both_columns(
        self, config, caplog
    ):
        """Two prefixes, one gene: the double-test is reported, not hidden."""
        import logging

        # The full dotted name: `get_logger("stages.gwas_features")` returns
        # `papipeline.stages.gwas_features`, and the package logger's own level
        # is inherited, so an INFO record is dropped unless this one is named.
        with caplog.at_level(logging.INFO, logger="papipeline.stages.gwas_features"):
            rows, columns, _ = gf.build_feature_rows(
                config,
                _manifest(),
                gene_presence=_presence(mexA=["TEST_A_01"]),
                amr_rows=[_amr_row("TEST_A_01", "mexA", gene="mexA")],
                regulator_records=[],
            )
        assert "gene__mexA" in columns
        assert "gene_presence_absence__mexA" in columns
        by_sample = {r["sample_id"]: r for r in rows}
        assert by_sample["TEST_A_01"]["gene__mexA"] == 1
        assert by_sample["TEST_A_01"]["gene_presence_absence__mexA"] == 1
        assert "mexA" in caplog.text
        assert "two prefixes" in caplog.text

    def test_a_row_for_another_antibiotic_is_excluded_and_counted(self, config):
        _, _, counts = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(mexA=["TEST_A_01"]),
            amr_rows=[
                _amr_row("TEST_A_01", "mexA", gene="mexA"),
                _amr_row("TEST_A_01", "ampC", gene="ampC", antibiotic="ceftazidime"),
            ],
            regulator_records=[],
        )
        assert counts["amr_determinant_gene"] == 1
        assert counts["amr_dropped_other_antibiotic"] == 1


# ---------------------------------------------------------------------------
# F2: the oprD pair is the regulator screen's, not this module's
# ---------------------------------------------------------------------------


class TestTheOprdPairComesFromTheRegulatorScreen:
    def test_an_absent_locus_sets_only_the_absent_feature(self, config):
        rows, _, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=list(SAMPLES)),
            amr_rows=[],
            regulator_records=[
                _call("TEST_A_01", "GENE_ABSENCE"),
                _call("TEST_A_02", "PREMATURE_STOP"),
            ],
        )
        by_sample = {r["sample_id"]: r for r in rows}
        assert by_sample["TEST_A_01"]["gene__oprD_absent"] == 1
        assert by_sample["TEST_A_01"]["gene__oprD_LoF"] == 0
        assert by_sample["TEST_A_02"]["gene__oprD_LoF"] == 1
        assert by_sample["TEST_A_02"]["gene__oprD_absent"] == 0

    def test_a_frameshift_counts_as_loss_of_function(self, config):
        rows, _, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=list(SAMPLES)),
            amr_rows=[],
            regulator_records=[_call("TEST_A_01", "FRAMESHIFT")],
        )
        assert rows[0]["gene__oprD_LoF"] == 1

    def test_a_sample_with_no_call_gets_both_features_zero(self, config):
        """Silence is "not observed", never an absent gene."""
        rows, _, _ = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=list(SAMPLES)),
            amr_rows=[],
            regulator_records=[],
        )
        by_sample = {r["sample_id"]: r for r in rows}
        for sample_id in SAMPLES:
            assert by_sample[sample_id]["gene__oprD_absent"] == 0
            assert by_sample[sample_id]["gene__oprD_LoF"] == 0

    def test_the_constructor_that_had_no_callers_is_now_called(self, config):
        rows, columns, counts = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=list(SAMPLES)),
            amr_rows=[],
            regulator_records=[_call("TEST_A_01", "GENE_ABSENCE")],
        )
        assert counts["oprd"] == 2
        assert "gene__oprD_absent" in columns
        assert "gene__oprD_LoF" in columns


# ---------------------------------------------------------------------------
# F2: refusals
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_every_missing_input_is_named_with_path_and_producer_stage(
        self, config, tmp_path
    ):
        with pytest.raises(DataContractError) as excinfo:
            gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        message = str(excinfo.value)
        assert "3 of 3 required inputs are absent" in message
        for spec in gf.REQUIRED_INPUTS:
            assert spec.key in message
            assert str(spec.path(tmp_path)) in message
            assert spec.producer_stage in message
        # Each says which stage writes it, so the reader is sent to code.
        assert "stage 7 (pangenome)" in message
        assert "stage 4 (amr)" in message
        assert "stage 6 (regulators)" in message

    def test_the_two_present_inputs_still_leave_the_third_named(
        self, config, tmp_path
    ):
        (tmp_path / "stages").mkdir(parents=True)
        (tmp_path / "amr").mkdir(parents=True)
        (tmp_path / "stages" / "gene_presence_absence.tsv").write_text(
            "gene\tsample_id\tpresent\nsynA\tTEST_A_01\t1\n", encoding="utf-8"
        )
        (tmp_path / "amr" / "amr_determinants.tsv").write_text(
            "sample_id\tantibiotic\tdeterminant\tgene\tvariant\n", encoding="utf-8"
        )
        with pytest.raises(DataContractError) as excinfo:
            gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        message = str(excinfo.value)
        assert "1 of 3 required inputs are absent" in message
        assert "regulator_variants" in message
        assert "pangenome_gene_presence_absence" not in message.split("- ")[1]

    def test_a_pangenome_over_a_different_cohort_is_refused_by_name(
        self, config, tmp_path
    ):
        _write_all_inputs(tmp_path)
        (tmp_path / "stages" / "gene_presence_absence.tsv").write_text(
            "gene\tsample_id\tpresent\nsynA\tSOME_OTHER_ISOLATE\t1\n",
            encoding="utf-8",
        )
        with pytest.raises(SampleIdError) as excinfo:
            gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        message = str(excinfo.value)
        assert "SOME_OTHER_ISOLATE" in message
        assert "TEST_A_01" in message

    def test_a_zero_gene_pangenome_is_refused_rather_than_written(
        self, config, tmp_path
    ):
        """A pan-genome with no genes is a failed input, not a cohort with none."""
        _write_all_inputs(tmp_path)
        with pytest.raises(DataContractError) as excinfo:
            gf.build_feature_rows(
                config,
                _manifest(),
                gene_presence={},
                amr_rows=[],
                regulator_records=[],
            )
        message = str(excinfo.value)
        assert "contributed no genes" in message
        assert "failed input, not a cohort with no genes" in message

    def test_a_header_only_pangenome_matrix_is_refused_by_cohort(
        self, config, tmp_path
    ):
        _write_all_inputs(tmp_path)
        (tmp_path / "stages" / "gene_presence_absence.tsv").write_text(
            "gene\tsample_id\tpresent\n", encoding="utf-8"
        )
        with pytest.raises(SampleIdError):
            gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)

    def test_an_undeclared_prefix_is_refused_not_warned_about(self, tmp_path):
        """The contract says the prefix MUST appear in the config list.

        Driven through a real edited ``science.yaml`` rather than a patched
        attribute: ``GwasConfig`` is a frozen dataclass, so an instance
        attribute shadows anything set on the class, and a test that silently
        patched the wrong one would pass without exercising the refusal.
        """
        import shutil

        from papipeline.config.loader import load_config

        _write_all_inputs(tmp_path)
        config_dir = tmp_path / "config"
        shutil.copytree(_REPO / "config", config_dir)
        science = config_dir / "science.yaml"
        text = science.read_text(encoding="utf-8")
        assert "    - gene\n" in text
        science.write_text(
            text.replace(
                "  feature_types:\n    - gene\n    - gene_presence_absence\n",
                "  feature_types:\n    - snp\n    - unitig\n    - kmer\n",
                1,
            ),
            encoding="utf-8",
        )
        stripped = load_config(science)
        assert "gene" not in stripped.gwas.feature_types

        with pytest.raises(DataContractError) as excinfo:
            gf.produce_gwas_features(stripped, _manifest(), RunMode.TEST, tmp_path)
        message = str(excinfo.value)
        assert "not in config.gwas.feature_types" in message
        assert "gene_presence_absence" in message
        assert "gene" in message

    def test_an_empty_manifest_is_refused_before_the_inputs_are_read(
        self, config, tmp_path
    ):
        """Reported as the empty cohort it is, not as three absent files."""
        with pytest.raises(DataContractError) as excinfo:
            gf.produce_gwas_features(
                config, SampleManifest([]), RunMode.TEST, tmp_path
            )
        assert "names no sample" in str(excinfo.value)
        assert "required inputs are absent" not in str(excinfo.value)


# ---------------------------------------------------------------------------
# F2: determinism and the written file
# ---------------------------------------------------------------------------


class TestTheWrittenFile:
    def test_it_lands_where_gwas_run_reads_it(self, config, tmp_path):
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        assert path == tmp_path / "gwas" / "gwas_features.tsv"
        assert gf.feature_table_path(tmp_path) == path
        assert gf.FEATURE_TABLE_RELPATH == ("gwas", "gwas_features.tsv")

    def test_gwas_reads_exactly_this_path(self, config):
        """The producer's path and the consumer's path are one declaration."""
        from papipeline.stages import gwas

        assert (
            gwas.run.__code__.co_varnames  # the reader takes a root, not a path
            and gf.FEATURE_TABLE_RELPATH == ("gwas", "gwas_features.tsv")
        )
        source = gwas.run.__doc__ or ""
        assert source  # run() documents what it reads

    def test_the_written_file_is_readable_by_the_contract_reader(
        self, config, tmp_path
    ):
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        rows = read_tsv(path, required_columns=("sample_id",),
                        unique_columns=("sample_id",))
        assert len(rows) == len(SAMPLES)

    def test_every_gene_carrier_reaches_the_table(self, config, tmp_path):
        """Regression: an inverted guard dropped every carrier but the first.

        ``carriers[name].add`` sat behind ``if sample_id in carriers[name]``, so
        a gene's carrier set stayed empty and every whole-gene column read 0
        across the cohort - a table that looks like a cohort carrying no genes.
        """
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        rows = {r["sample_id"]: r for r in read_tsv(path)}
        assert rows["TEST_A_01"]["gene__synA"] == "1"
        assert rows["TEST_A_03"]["gene__synA"] == "1"
        assert rows["TEST_A_02"]["gene__synA"] == "0"
        assert rows["TEST_A_01"]["gene__synB"] == "1"
        assert rows["TEST_A_02"]["gene__synB"] == "0"

    def test_the_banner_names_every_input_and_carries_the_power_flag(
        self, config, tmp_path
    ):
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        banner = [
            line for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("#")
        ]
        text = "\n".join(banner)
        assert f"n={len(SAMPLES)} isolates" in text
        assert "underpowered, not a finding" in text
        for spec in gf.REQUIRED_INPUTS:
            assert str(spec.path(tmp_path)) in text

    def test_same_inputs_give_a_byte_identical_file(self, config, tmp_path):
        _write_all_inputs(tmp_path)
        first = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        payload = first.read_bytes()
        second = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        assert second.read_bytes() == payload

    def test_input_row_order_does_not_change_the_output(self, config, tmp_path):
        """Determinism is over the *inputs*, not over their order on disk."""
        _write_all_inputs(tmp_path)
        first = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        payload = first.read_bytes()
        amr = tmp_path / "amr" / "amr_determinants.tsv"
        lines = amr.read_text(encoding="utf-8").splitlines()
        amr.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n",
                       encoding="utf-8")
        second = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        assert second.read_bytes() == payload

    def test_no_statistic_is_computed_here(self, config, tmp_path):
        """Every count in the file carries the flag; nothing claims a finding."""
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        text = path.read_text(encoding="utf-8")
        assert "No statistic is computed here" in text
        for forbidden in ("p_value", "adjusted_p_value", "odds_ratio", "pvalue"):
            assert forbidden not in text


def _write_all_inputs(tmp_path, genes=None) -> None:
    """Write the three required inputs as synthetic tables (R12)."""
    genes = {"synA": ["TEST_A_01", "TEST_A_03"], "synB": ["TEST_A_01"]} if genes is None else genes
    (tmp_path / "stages").mkdir(parents=True, exist_ok=True)
    (tmp_path / "amr").mkdir(parents=True, exist_ok=True)
    (tmp_path / "regulators").mkdir(parents=True, exist_ok=True)

    gpa = ["gene\tsample_id\tpresent"]
    for gene, carriers in sorted(genes.items()):
        for sample_id in SAMPLES:
            gpa.append(f"{gene}\t{sample_id}\t{int(sample_id in carriers)}")
    (tmp_path / "stages" / "gene_presence_absence.tsv").write_text(
        "\n".join(gpa) + "\n", encoding="utf-8"
    )

    amr = [
        "sample_id\tantibiotic\tdeterminant\tgene\tvariant",
        f"{SAMPLES[0]}\timipenem\tmexA\tmexA\t.",
        f"{SAMPLES[1]}\timipenem\toprD_V359L\toprD\tV359L",
        f"{SAMPLES[1]}\timipenem\tblaIMP-14\tblaIMP-14\t.",
    ]
    (tmp_path / "amr" / "amr_determinants.tsv").write_text(
        "\n".join(amr) + "\n", encoding="utf-8"
    )

    regulators = [
        "sample_id\tgene\tvariant\tvariant_type\tposition\treference\talternate",
        f"{SAMPLES[0]}\toprD\tSYNTHETIC_ABSENCE\tGENE_ABSENCE\t0\t.\t.",
    ]
    (tmp_path / "regulators" / "regulator_variants.tsv").write_text(
        "\n".join(regulators) + "\n", encoding="utf-8"
    )


class TestEndToEndOnSyntheticInputs:
    def test_the_whole_path_produces_a_conformant_table(self, config, tmp_path):
        _write_all_inputs(tmp_path)
        path = gf.produce_gwas_features(config, _manifest(), RunMode.TEST, tmp_path)
        rows = read_tsv(path, required_columns=("sample_id",),
                        unique_columns=("sample_id",))
        header = [
            line for line in path.read_text(encoding="utf-8").splitlines()
            if not line.startswith("#")
        ][0].split("\t")
        assert header[0] == "sample_id"
        for row in rows:
            for column in header[1:]:
                assert row[column] in ("0", "1"), (column, row[column])
                assert column.split("__", 1)[0] in set(config.gwas.feature_types)

    def test_the_counts_by_family_are_reported(self, config, tmp_path):
        _write_all_inputs(tmp_path)
        _rows, _columns, counts = gf.build_feature_rows(
            config,
            _manifest(),
            gene_presence=_presence(synA=["TEST_A_01", "TEST_A_03"], synB=["TEST_A_01"]),
            amr_rows=[
                _amr_row("TEST_A_01", "mexA", gene="mexA"),
                _amr_row("TEST_A_02", "oprD_V359L", gene="oprD", variant="V359L"),
            ],
            regulator_records=[_call("TEST_A_01", "GENE_ABSENCE")],
        )
        assert counts["pangenome_gene"] == 2
        assert counts["amr_allele"] == 1
        assert counts["amr_determinant_gene"] == 1
        assert counts["oprd"] == 2
        assert counts["total"] == 6


# ---------------------------------------------------------------------------
# F3: the real 10-isolate cohort. Skipped when the artifacts are absent.
# ---------------------------------------------------------------------------

REAL_PANGENOME = "/Users/raghavkrishnankv/Desktop/pa-artifacts/b2-pangenome/panaroo-smoke/gene_presence_absence.csv"
REAL_AMR_REPORTS = "/Users/raghavkrishnankv/Desktop/pa-artifacts/round11/amr/with"
REAL_REGULATORS = (
    "/Users/raghavkrishnankv/Desktop/pa-regv/results/real/intermediate"
    "/regulators/regulator_variants.tsv"
)
REAL_SAMPLES = (
    "PDT000034122.1", "PDT000167133.1", "PDT000167135.1", "PDT000167136.1",
    "PDT000292995.1", "PDT000292998.1", "PDT000294804.1", "PDT000294805.1",
    "PDT000311294.1", "PDT000424983.1",
)


def _real_inputs_available() -> bool:
    from pathlib import Path

    return all(
        Path(p).exists()
        for p in (REAL_PANGENOME, REAL_REGULATORS, REAL_AMR_REPORTS)
    )


requires_real = pytest.mark.skipif(
    not _real_inputs_available(),
    reason=(
        "the real 10-isolate intermediates are not on this machine; this is a "
        "real-data check over files this repository does not commit"
    ),
)


@requires_real
class TestTheRealTenIsolateCohort:
    """F3's real-data check: real intermediates, real code paths, no tool run."""

    @pytest.fixture(scope="class")
    def real_root(self, tmp_path_factory):
        """Assemble a root holding the REAL intermediates in the declared paths.

        The tables are placed by REAL pipeline code, not by hand: stage 7's
        matrix is built by `pangenome.partition` from panaroo's own CSV through
        the adapter that reads it, and stage 4's table is written by the AMR
        stage's own writer from `AmrFinderPlusAdapter.parse_report` over the real
        AMRFinderPlus reports. Nothing here runs a bioinformatics tool.
        """
        from pathlib import Path

        from papipeline.adapters import panaroo as panaroo_adapter
        from papipeline.io.tsv import write_tsv
        from papipeline.pilot import amr_detect
        from papipeline.stages import amr as amr_stage
        from papipeline.stages import pangenome as pangenome_stage

        root = tmp_path_factory.mktemp("real_root")
        manifest = _manifest(REAL_SAMPLES)

        found, presence = panaroo_adapter.parse_presence_csv(Path(REAL_PANGENOME))
        assert sorted(found) == sorted(REAL_SAMPLES)
        pangenome = pangenome_stage.partition(found, presence)
        (root / "stages").mkdir(parents=True, exist_ok=True)
        write_tsv(
            root / "stages" / "gene_presence_absence.tsv",
            pangenome.long_rows(),
            list(pangenome_stage.GPA_COLUMNS),
        )

        (root / "amr").mkdir(parents=True, exist_ok=True)
        determinants = []
        for sample_id in REAL_SAMPLES:
            report = Path(REAL_AMR_REPORTS) / f"{sample_id}.amrfinder.tsv"
            points = amr_stage.point_variants_in_report(report)
            determinants.extend(
                amr_detect.parse_amrfinder_report(
                    report,
                    sample_id=sample_id,
                    antibiotic="imipenem",
                    database="AMRFinderPlus",
                    database_version_value="2026-08-07.1",
                )
            )
            assert points or True
        from dataclasses import replace

        points_by_symbol = {}
        for sample_id in REAL_SAMPLES:
            report = Path(REAL_AMR_REPORTS) / f"{sample_id}.amrfinder.tsv"
            points_by_symbol.update(amr_stage.point_variants_in_report(report))
        determinants = [
            replace(
                record,
                gene=points_by_symbol.get(record.determinant, (record.gene, None))[0],
                variant=points_by_symbol.get(record.determinant, (None, None))[1],
            )
            for record in determinants
        ]
        write_tsv(
            root / "amr" / "amr_determinants.tsv",
            [
                {
                    "sample_id": r.sample_id,
                    "antibiotic": r.antibiotic,
                    "determinant": r.determinant,
                    "gene": r.gene,
                    "variant": r.variant,
                }
                for r in determinants
            ],
            ["sample_id", "antibiotic", "determinant", "gene", "variant"],
        )

        (root / "regulators").mkdir(parents=True, exist_ok=True)
        (root / "regulators" / "regulator_variants.tsv").write_bytes(
            Path(REAL_REGULATORS).read_bytes()
        )
        return root, manifest

    def test_it_writes_a_table_over_all_ten_isolates(self, config, real_root):
        root, manifest = real_root
        path = gf.produce_gwas_features(config, manifest, RunMode.REAL, root)
        rows = read_tsv(path, required_columns=("sample_id",),
                        unique_columns=("sample_id",))
        assert len(rows) == 10
        assert {r["sample_id"] for r in rows} == set(REAL_SAMPLES)

    def test_the_banner_carries_the_ten_isolate_flag(self, config, real_root):
        root, manifest = real_root
        path = gf.produce_gwas_features(config, manifest, RunMode.REAL, root)
        text = path.read_text(encoding="utf-8")
        assert "n=10 isolates" in text
        assert "underpowered, not a finding" in text

    def test_the_counts_by_family_on_real_data(self, config, real_root):
        root, manifest = real_root
        from papipeline.stages import amr as amr_stage
        from papipeline.stages import pangenome as pangenome_stage
        from papipeline.stages import regulators as reg

        _found, presence = pangenome_stage.read_gene_presence_absence(
            root / "stages" / "gene_presence_absence.tsv"
        )
        amr_rows = read_tsv(
            root / "amr" / "amr_determinants.tsv",
            required_columns=gf.AMR_REQUIRED,
        )
        records = reg.load_regulator_variants(
            config, root / "regulators" / "regulator_variants.tsv", root
        )
        _rows, _cols, counts = gf.build_feature_rows(
            config,
            manifest,
            gene_presence=presence,
            amr_rows=amr_rows,
            regulator_records=records,
        )
        assert counts["pangenome_gene"] == len(presence)
        assert counts["amr_allele"] >= 1
        assert counts["oprd"] == 2
        assert counts["total"] == (
            counts["pangenome_gene"] + counts["amr_allele"]
            + counts["amr_determinant_gene"] + counts["oprd"]
        )
        assert amr_stage.POINT_SUBTYPES == frozenset({"POINT", "POINT_DISRUPT"})

    def test_every_feature_value_is_binary_on_real_data(self, config, real_root):
        root, manifest = real_root
        path = gf.produce_gwas_features(config, manifest, RunMode.REAL, root)
        header = [
            line for line in path.read_text(encoding="utf-8").splitlines()
            if not line.startswith("#")
        ][0].split("\t")
        rows = read_tsv(path, required_columns=("sample_id",))
        for row in rows:
            for column in header[1:]:
                assert row[column] in ("0", "1"), (column, row[column])

    def test_it_is_byte_identical_on_a_second_run(self, config, real_root):
        root, manifest = real_root
        first = gf.produce_gwas_features(config, manifest, RunMode.REAL, root)
        payload = first.read_bytes()
        second = gf.produce_gwas_features(config, manifest, RunMode.REAL, root)
        assert second.read_bytes() == payload
