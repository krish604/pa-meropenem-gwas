"""Meropenem, as configuration: enabled, declared, and reachable.

``docs/architecture.md`` ("Adding an antibiotic") says no code change is
required to add an antibiotic - three configuration edits:

1. ``config/science.yaml`` names it in ``antibiotics:`` (the file the
   architecture doc calls ``config/config.yaml``; this repository's science
   configuration is ``config/science.yaml`` and it is the file that carries
   the ``antibiotics:`` key);
2. ``config/antibiotics.tsv`` carries a row with its class and phenotype
   standard;
3. ``config/mechanisms.tsv`` extends the ``antibiotic`` column of the loci
   that apply, as a comma-separated list or ``all``.

These tests pin all three for meropenem, plus the phenotype filename
convention ``<antibiotic>_phenotype.tsv`` that stage 11 and the cohort gate
resolve through.

Synthetic only: everything here reads the committed configuration tables and
``tmp_path``. Nothing under ``data/``, ``db/`` or ``PDC_essential.tsv`` is
read, and no real-mode run happens.
"""

from __future__ import annotations

from papipeline.stages.phenotype import (
    PHENOTYPE_FILENAME_TEMPLATE,
    phenotype_path,
)


class TestMeropenemIsEnabled:
    """The science configuration names meropenem in its ``antibiotics:`` list."""

    def test_meropenem_is_in_the_antibiotics_list(self, config):
        assert "meropenem" in config.antibiotics

    def test_imipenem_stays_the_first_position(self, config):
        # `config.antibiotics[0]` is read as "the project's antibiotic" by
        # several stages and scripts as a fallback when no explicit choice is
        # configured. Enabling meropenem must not move that position.
        assert config.antibiotics[0] == "imipenem"
        assert config.antibiotics == ("imipenem", "meropenem")

    def test_require_antibiotic_accepts_meropenem(self, config):
        assert config.require_antibiotic("meropenem") == "meropenem"

    def test_an_unconfigured_antibiotic_is_still_refused(self, config):
        from papipeline.errors import AntibioticError
        import pytest

        with pytest.raises(AntibioticError, match="not configured"):
            config.require_antibiotic("ceftazidime")


class TestTheAntibioticsTable:
    """``config/antibiotics.tsv`` declares meropenem with its class and standard."""

    def test_meropenem_has_a_row(self, config):
        assert "meropenem" in config.antibiotic_specs

    def test_meropenem_is_a_carbapenem_with_a_named_standard(self, config):
        spec = config.antibiotic_specs["meropenem"]
        assert spec.antibiotic_class == "carbapenem"
        assert spec.phenotype_source_standard

    def test_both_carbapenems_share_the_same_scope(self, config):
        imi = config.antibiotic_specs["imipenem"]
        mero = config.antibiotic_specs["meropenem"]
        assert mero.mechanism_scope == imi.mechanism_scope


class TestTheMechanismsTable:
    """The loci that apply to carbapenems list meropenem in the antibiotic column."""

    #: A sample across every mechanism family the table carries, including the
    #: co-occurrence context rows whose notes explain they are in the table
    #: precisely so the analysis is not biased towards one drug class.
    CARBAPENEM_LOCI = (
        "oprD",  # reduced permeability
        "mexZ",  # efflux regulation
        "ampD",  # AmpC regulation
        "blaVIM-1",  # acquired carbapenemase
        "blaTEM-1",  # context row: not a carbapenemase, still curated context
        "aac(6')-Ib-cr",  # context row: non-beta-lactam, co-occurrence anchor
    )

    def test_every_carbapenem_locus_is_relevant_to_meropenem(self, config):
        for gene in self.CARBAPENEM_LOCI:
            spec = config.mechanism_for_gene(gene)
            assert spec.relevant_to("meropenem"), gene
            assert spec.relevant_to("imipenem"), gene

    def test_the_comma_list_does_not_declare_unlisted_antibiotics(self, config):
        oprd = config.mechanism_for_gene("oprD")
        assert not oprd.relevant_to("ceftazidime")


class TestThePhenotypeFilenameConvention:
    """The phenotype filename follows ``<antibiotic>_phenotype.tsv``."""

    def test_meropenem_resolves_to_its_own_file(self, tmp_path):
        path = phenotype_path(tmp_path, "meropenem")
        assert path.name == "meropenem_phenotype.tsv"
        assert path.name == PHENOTYPE_FILENAME_TEMPLATE.format(
            antibiotic="meropenem"
        )
        assert path.parent == tmp_path
