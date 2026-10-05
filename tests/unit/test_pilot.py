"""Tests for the pilot-100 modules.

These cover the selection rule, the PDC field parsers, the mechanism
categorisation and the analysis tables. The selection tests are the important
ones: they assert that the cohort is defined by the filesystem and that
PDC_essential.tsv row order cannot influence it.
"""

from __future__ import annotations

import pytest

from papipeline.errors import PipelineError
from papipeline.models import Phenotype
from papipeline.pilot import analysis as A
from papipeline.pilot import cohort as C
from papipeline.pilot import mechanisms as M
from papipeline.pilot import pdc_fields as P


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------


class TestGCASortKey:
    def test_numeric_not_lexicographic(self):
        """GCA_9 must sort before GCA_10, which plain string sort gets wrong."""
        accessions = ["GCA_10.1", "GCA_9.1", "GCA_2.5"]
        assert sorted(accessions, key=C.gca_sort_key) == [
            "GCA_2.5",
            "GCA_9.1",
            "GCA_10.1",
        ]

    def test_version_is_numeric(self):
        assert C.gca_sort_key("GCA_1.10") > C.gca_sort_key("GCA_1.9")

    def test_rejects_non_gca(self):
        with pytest.raises(PipelineError, match="Not a GCA"):
            C.gca_sort_key("NOT_AN_ACCESSION")


class TestSuffixDetection:
    @pytest.mark.parametrize(
        "name",
        [
            "GCA_1.1.fna",
            "GCA_1.1.fasta",
            "GCA_1.1.fa",
            "GCA_1.1.fna.gz",
            "GCA_1.1.fasta.gz",
            "GCA_1.1.fa.gz",
            "GCA_1.1.FNA",
        ],
    )
    def test_recognised(self, name):
        assert C._suffix_of(name) is not None

    @pytest.mark.parametrize("name", ["GCA_1.1.txt", "GCA_1.1.json", "GCA_1.1.fna.bz2", "readme"])
    def test_not_recognised(self, name):
        assert C._suffix_of(name) is None


class TestDiscoverAssemblies:
    def _make(self, tmp_path, names):
        for name in names:
            path = tmp_path / "data" / name.split("_")[0] / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(">c\nACGT\n")
        return tmp_path / "data"

    def test_finds_assemblies(self, tmp_path):
        root = self._make(
            tmp_path, ["GCA_000710625.1_x.fna", "GCA_000937465.2_y.fna"]
        )
        found = C.discover_assemblies(root)
        assert set(found) == {"GCA_000710625.1", "GCA_000937465.2"}

    def test_ignores_non_assembly_files(self, tmp_path):
        root = self._make(tmp_path, ["GCA_000710625.1_x.fna", "GCA_000937465.2_notes.txt"])
        found = C.discover_assemblies(root)
        assert set(found) == {"GCA_000710625.1"}

    def test_ignores_assembly_extension_without_gca(self, tmp_path):
        root = tmp_path / "data"
        (root / "random").mkdir(parents=True)
        (root / "random" / "contig.fna").write_text(">c\nACGT\n")
        assert C.discover_assemblies(root) == {}

    def test_gzipped_assemblies_are_found(self, tmp_path):
        root = self._make(tmp_path, ["GCA_000710625.1_x.fna.gz"])
        assert "GCA_000710625.1" in C.discover_assemblies(root)

    def test_missing_directory_raises(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            C.discover_assemblies(tmp_path / "absent")

    def test_selects_in_sorted_order_not_directory_order(self, tmp_path):
        names = [
            "GCA_000937465.2_b.fna",
            "GCA_000710625.1_a.fna",
            "GCA_001874795.1_c.fna",
        ]
        root = self._make(tmp_path, names)
        found = C.discover_assemblies(root)
        selected = C.select_first_n(found, n=2)
        assert [s.assembly for s in selected] == [
            "GCA_000710625.1",
            "GCA_000937465.2",
        ]
        assert [s.pilot_id for s in selected] == ["PILOT100_001", "PILOT100_002"]

    def test_refuses_to_substitute_when_short(self, tmp_path):
        root = self._make(tmp_path, ["GCA_000710625.1_a.fna"])
        found = C.discover_assemblies(root)
        with pytest.raises(PipelineError, match="refusing to substitute"):
            C.select_first_n(found, n=100)

    def test_pilot_ids_are_sequential(self, tmp_path):
        names = [f"GCA_{i:09d}.1_x.fna" for i in range(1, 6)]
        root = self._make(tmp_path, names)
        selected = C.select_first_n(C.discover_assemblies(root), n=5)
        assert [s.pilot_id for s in selected] == [f"PILOT100_00{i}" for i in range(1, 6)]


class TestPDCIndex:
    def test_indexes_on_assembly_not_row_order(self, write_tsv_file, tmp_path):
        path = write_tsv_file(
            "PDC_essential.tsv",
            "Assembly\tIsolate\nGCA_000937465.2\tB\nGCA_000710625.1\tA\n",
        )
        index = C.load_pdc_index(path)
        assert index["GCA_000710625.1"]["Isolate"] == "A"
        assert index["GCA_000937465.2"]["Isolate"] == "B"

    def test_empty_file_raises(self, write_tsv_file):
        path = write_tsv_file("PDC_essential.tsv", "")
        with pytest.raises(PipelineError, match="empty"):
            C.load_pdc_index(path)

    def test_missing_assembly_column_raises(self, write_tsv_file):
        path = write_tsv_file("PDC_essential.tsv", "Isolate\nA\n")
        with pytest.raises(PipelineError, match="no 'Assembly' column"):
            C.load_pdc_index(path)

    def test_duplicate_accession_raises(self, write_tsv_file):
        path = write_tsv_file(
            "PDC_essential.tsv",
            "Assembly\tIsolate\nGCA_1.1\tA\nGCA_1.1\tB\n",
        )
        with pytest.raises(PipelineError, match="duplicate Assembly"):
            C.load_pdc_index(path)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            C.load_pdc_index(tmp_path / "absent.tsv")


class TestJoinMetadata:
    def _selected(self):
        return [
            C.SelectedAssembly("PILOT100_001", "GCA_1.1", C.Path("a.fna"), "a.fna"),
            C.SelectedAssembly("PILOT100_002", "GCA_2.2", C.Path("b.fna"), "b.fna"),
        ]

    def test_missing_record_is_retained_and_marked(self):
        rows, counts = C.join_metadata(self._selected(), {"GCA_1.1": {"Isolate": "A"}})
        assert counts == {"selected": 2, "matched": 1, "missing": 1}
        assert len(rows) == 2
        assert rows[1]["PDC_metadata_match"] == "MISSING"

    def test_missing_record_invents_nothing(self):
        rows, _ = C.join_metadata(self._selected(), {})
        for field in ("Isolate", "BioSample", "AMR_genotypes", "AST_phenotypes_raw"):
            assert rows[0][field] == "MISSING"

    def test_raw_ast_is_retained_for_per_antibiotic_extraction(self):
        rows, _ = C.join_metadata(
            self._selected(),
            {"GCA_1.1": {"AST phenotypes": "imipenem=R,meropenem=S"}},
        )
        assert rows[0]["AST_phenotypes_raw"] == "imipenem=R,meropenem=S"

    def test_matched_record_marked(self):
        rows, _ = C.join_metadata(self._selected(), {"GCA_1.1": {"Isolate": "A"}})
        assert rows[0]["PDC_metadata_match"] == "MATCHED"


class TestSelectionProof:
    def test_detects_that_pdc_row_order_differs(self, write_tsv_file, tmp_path):
        pdc = write_tsv_file(
            "PDC_essential.tsv",
            "Assembly\tIsolate\n"
            "GCA_000937465.2\tB\nGCA_000710625.1\tA\nGCA_001874795.1\tC\n",
        )
        selected = [
            C.SelectedAssembly("PILOT100_001", "GCA_000710625.1", C.Path("a"), "a"),
        ]
        proof = C.prove_selection_is_filesystem_based(pdc, selected)
        assert proof["identical_sets"] is False
        assert proof["overlap_n"] == 1
        assert proof["only_in_pdc_roworder"] == 2


# --------------------------------------------------------------------------
# PDC field parsing
# --------------------------------------------------------------------------


class TestParseAstField:
    def test_extracts_categories(self):
        result = P.parse_ast_field("imipenem=R,meropenem=S,ceftazidime=I")
        assert result.get("imipenem") is Phenotype.R
        assert result.get("meropenem") is Phenotype.S
        assert result.get("ceftazidime") is Phenotype.I

    def test_case_insensitive_antibiotic(self):
        assert P.parse_ast_field("Imipenem=R").get("imipenem") is Phenotype.R

    def test_all_allowed_categories(self):
        text = ",".join(f"a{i}={v}" for i, v in enumerate(["R", "I", "S", "SDD", "ND"]))
        result = P.parse_ast_field(text)
        assert set(result.values.values()) == {"R", "I", "S", "SDD", "ND"}

    def test_absent_antibiotic_is_none_not_defaulted(self):
        assert P.parse_ast_field("meropenem=R").get("imipenem") is None

    def test_no_mic_is_ever_produced(self):
        """Rule 6: a category never acquires a quantitative value."""
        result = P.parse_ast_field("imipenem=R")
        assert not hasattr(result, "mic")

    def test_duplicate_conflicting_key_is_ambiguous(self):
        """A conflicting duplicate is never resolved by taking first or last."""
        result = P.parse_ast_field("ciprofloxacin=I,ciprofloxacin=S")
        assert "ciprofloxacin" in result.ambiguous
        assert result.get("ciprofloxacin") is None
        assert result.availability("ciprofloxacin") == "AMBIGUOUS_DUPLICATE_KEY"

    def test_duplicate_identical_key_is_accepted(self):
        result = P.parse_ast_field("imipenem=R,imipenem=R")
        assert result.get("imipenem") is Phenotype.R

    def test_invalid_category_is_rejected_not_coerced(self):
        result = P.parse_ast_field("imipenem=RESISTANT")
        assert "imipenem" in result.invalid
        assert result.get("imipenem") is None
        assert result.availability("imipenem") == "INVALID_CATEGORY"

    def test_unparsable_token_counted(self):
        result = P.parse_ast_field("imipenem=R,garbage,,imipenem")
        assert result.unparsed >= 1

    def test_empty_input(self):
        result = P.parse_ast_field("")
        assert result.values == {}


class TestParseGenotypeField:
    def test_parses_gene_state_pairs(self):
        result = P.parse_genotype_field("blaPDC-3=COMPLETE,mexA=COMPLETE")
        assert set(result.genes) == {"blaPDC-3", "mexA"}
        assert len(result.calls) == 2

    def test_all_observed_states_accepted(self):
        text = ",".join(
            f"g{i}={s}"
            for i, s in enumerate(
                ["COMPLETE", "PARTIAL", "PARTIAL_END_OF_CONTIG", "POINT", "HMM", "MISTRANSLATION"]
            )
        )
        result = P.parse_genotype_field(text)
        assert len(result.calls) == 6
        assert result.unknown_states == set()

    def test_unknown_state_recorded_not_dropped(self):
        result = P.parse_genotype_field("gA=WEIRD_STATE")
        assert len(result.calls) == 1
        assert "WEIRD_STATE" in result.unknown_states

    def test_variant_detection(self):
        result = P.parse_genotype_field("oprD_V359L=POINT,blaOXA-486=COMPLETE")
        variants = {c.gene for c in result.variant_calls()}
        assert "oprD_V359L" in variants
        assert "blaOXA-486" not in variants

    def test_base_gene_of_variant(self):
        result = P.parse_genotype_field("oprD_V359L=POINT")
        assert result.calls[0].base_gene == "oprD"

    def test_numbered_allele_is_not_a_variant(self):
        result = P.parse_genotype_field("blaOXA-486=COMPLETE")
        assert result.calls[0].base_gene == "blaOXA-486"
        assert not result.calls[0].is_variant

    def test_frameshift_detected_as_variant(self):
        result = P.parse_genotype_field("mexZ_Q69Ter=POINT")
        assert result.calls[0].is_variant

    def test_unparsable_token_counted(self):
        result = P.parse_genotype_field("blaPDC-3=COMPLETE,garbage")
        assert result.unparsed == 1

    def test_empty_input(self):
        assert P.parse_genotype_field("").calls == []


# --------------------------------------------------------------------------
# mechanism categorisation
# --------------------------------------------------------------------------


class TestFunctionalCategory:
    @pytest.fixture
    def families(self):
        return M.gene_family_table()

    def test_pdc_by_name_prefix(self, families):
        assert M.functional_category("blaPDC-3", {"class": "BETA-LACTAM"}, families) == M.CATEGORY_PDC

    def test_oprd_by_name(self, families):
        assert M.functional_category("oprD", {"class": "BETA-LACTAM"}, families) == M.CATEGORY_OPRD

    def test_efflux_from_database_subclass(self, families):
        assert (
            M.functional_category("mexA", {"class": "EFFLUX", "subclass": "EFFLUX"}, families)
            == M.CATEGORY_EFFLUX
        )

    def test_mbl_from_family_table(self, families):
        """AMRFinderPlus does not label MBL, so the table supplies it."""
        annotation = {"class": "BETA-LACTAM", "subclass": "CARBAPENEM"}
        assert M.functional_category("blaVIM-6", annotation, families) == M.CATEGORY_MBL
        assert M.functional_category("blaNDM-1", annotation, families) == M.CATEGORY_MBL

    def test_esbl_from_family_table(self, families):
        annotation = {"class": "BETA-LACTAM", "subclass": "CEPHALOSPORIN"}
        assert M.functional_category("blaCTX-M-15", annotation, families) == M.CATEGORY_ESBL
        assert M.functional_category("blaOXA-486", annotation, families) == M.CATEGORY_ESBL

    def test_class_a_carbapenemase_is_not_called_an_mbl(self, families):
        """KPC is a serine carbapenemase. Labelling it MBL would be wrong."""
        annotation = {"class": "BETA-LACTAM", "subclass": "CARBAPENEM"}
        assert (
            M.functional_category("blaKPC-2", annotation, families)
            == M.CATEGORY_OTHER_BETALACTAMASE
        )

    def test_other_beta_lactamase_falls_through(self, families):
        annotation = {"class": "BETA-LACTAM", "subclass": "BETA-LACTAM"}
        assert (
            M.functional_category("blaSOME-1", annotation, families)
            == M.CATEGORY_OTHER_BETALACTAMASE
        )

    def test_other_amr(self, families):
        assert (
            M.functional_category("tet(G)", {"class": "TETRACYCLINE"}, families)
            == M.CATEGORY_OTHER_AMR
        )

    def test_no_annotation_is_unmapped_not_guessed(self, families):
        assert M.functional_category("mysteryGene", None, families) == M.CATEGORY_UNMAPPED

    def test_works_without_a_family_table(self):
        """A missing table must not crash; it just cannot split ESBL/MBL."""
        annotation = {"class": "BETA-LACTAM", "subclass": "CARBAPENEM"}
        assert M.functional_category("blaVIM-6", annotation, {}) in (
            M.CATEGORY_OTHER_BETALACTAMASE,
            M.CATEGORY_OTHER_AMR,
        )


class TestFamilyForGene:
    def test_strips_allele_number(self):
        assert M.family_for_gene("blaVIM-6", {"blaVIM": "MBL"}) == "MBL"

    def test_prefix_match(self):
        assert M.family_for_gene("blaCTX-M-15", {"blaCTX": "ESBL"}) == "ESBL"

    def test_longest_prefix_wins(self):
        families = {"blaOXA": "ESBL", "blaOXAES": "MBL"}
        assert M.family_for_gene("blaOXAES-4", families) == "MBL"

    def test_allele_number_is_stripped_before_matching(self):
        """`blaOXA-486` is a distinct enzyme name, but the family is blaOXA."""
        assert M.family_for_gene("blaOXA-486", {"blaOXA": "ESBL"}) == "ESBL"

    def test_no_match(self):
        assert M.family_for_gene("mexA", {"blaVIM": "MBL"}) is None

    def test_empty_inputs(self):
        assert M.family_for_gene("", {"a": "b"}) is None
        assert M.family_for_gene("mexA", {}) is None


class TestOprdRules:
    def test_calls_are_case_insensitive(self):
        """Regression: a case-sensitive match found nothing at all."""
        assert M.oprd_calls(["oprD_V359L", "mexA"]) == ["oprD_V359L"]

    def test_terminates_are_disruptive(self):
        assert "oprD_W339Ter" in M.oprd_disrupting_genes(["oprD_W339Ter"])

    def test_frameshifts_are_disruptive(self):
        assert "oprD_Y237TerfsTer0" in M.oprd_disrupting_genes(["oprD_Y237TerfsTer0"])

    def test_missense_is_not_treated_as_disruption(self):
        """A missense change is reported, not assumed disruptive."""
        assert M.oprd_disrupting_genes(["oprD_V359L"]) == set()

    def test_locus_reported_for_missense_only(self):
        assert M.oprd_locus_present(["oprD_V359L"]) is True

    def test_locus_absent_when_no_oprd_call(self):
        assert M.oprd_locus_present(["mexA", "blaPDC-3"]) is False

    def test_other_genes_never_counted_as_oprd(self):
        assert M.oprd_disrupting_genes(["mexR_I24AfsTer94"]) == set()


class TestAssign:
    @pytest.fixture
    def config(self):
        from papipeline.config.loader import load_config
        from pathlib import Path as P

        return load_config(P(__file__).resolve().parents[2] / "config" / "science.yaml")

    def test_oprd_gene_alone_is_locus_intact(self, config):
        """Rule 1: an oprD gene is an intact locus, not reduced permeability."""
        assignment = M.assign("oprD", {"class": "BETA-LACTAM"}, config)
        assert assignment.pa_mechanism == "locus_intact"
        assert "no_disruption_reported" in assignment.evidence

    def test_oprd_with_disruption_is_reduced_permeability(self, config):
        assignment = M.assign(
            "oprD", {"class": "BETA-LACTAM"}, config, oprd_disrupted=True
        )
        assert assignment.pa_mechanism == "reduced_permeability"

    def test_known_gene_gets_knowledge_table_mechanism(self, config):
        assignment = M.assign("mexR", {"class": "EFFLUX"}, config)
        assert assignment.pa_mechanism == "efflux"

    def test_efflux_component_is_functional_but_pa_unmapped(self, config):
        """mexA is in no knowledge-table row, so no PA mechanism is asserted.

        The functional layer still calls it Efflux, from the database's own
        subclass. The two layers are independent by design.
        """
        assignment = M.assign(
            "mexA", {"class": "EFFLUX", "subclass": "EFFLUX"}, config
        )
        assert assignment.functional_category == M.CATEGORY_EFFLUX
        assert assignment.pa_mechanism is None

    def test_unknown_gene_gets_no_pa_mechanism(self, config):
        """An unlisted gene is UNMAPPED for the P. aeruginosa layer."""
        assignment = M.assign("someNovelGene", {"class": "NEW"}, config)
        assert assignment.pa_mechanism is None

    def test_evidence_is_never_causal(self, config):
        assignment = M.assign("blaVIM-6", {"class": "BETA-LACTAM"}, config)
        assert "caus" not in assignment.evidence.lower()

    def test_claim_level_never_exceeds_detected(self, config):
        assignment = M.assign("mexA", {"class": "EFFLUX"}, config)
        assert assignment.pa_evidence_level.value == "DETECTED"


# --------------------------------------------------------------------------
# assembly integrity (added after data/ proved to be corrupt)
# --------------------------------------------------------------------------


class TestAssessAssembly:
    """These conditions were all found in the real data/ directory."""

    def test_missing_file_excluded(self, tmp_path):
        verdict = C.assess_assembly(tmp_path / "absent.fna")
        assert verdict.included is False
        assert verdict.reason == "file_not_found"

    def test_zero_byte_excluded(self, tmp_path):
        path = tmp_path / "empty.fna"
        path.write_bytes(b"")
        verdict = C.assess_assembly(path)
        assert verdict.included is False
        assert verdict.reason == "zero_byte_file"

    def test_corrupt_binary_excluded(self, tmp_path):
        """The observed signature: valid FASTA, then binary garbage."""
        path = tmp_path / "corrupt.fna"
        path.write_bytes(b">c\n" + b"ACGT" * 32 + b"\n" + b"\xbf\x0e-\xbf\x8e")
        verdict = C.assess_assembly(path)
        assert verdict.included is False
        assert verdict.reason.startswith("corrupt_binary_data_at_offset_")

    def test_undersized_excluded(self, tmp_path):
        path = tmp_path / "small.fna"
        path.write_bytes(b">c\n" + b"ACGT" * 100)
        verdict = C.assess_assembly(path, min_size=4_000_000)
        assert verdict.included is False
        assert verdict.reason.startswith("undersized_")

    def test_non_fasta_excluded(self, tmp_path):
        path = tmp_path / "bin.fna"
        path.write_bytes(b"not a fasta file at all, just text\n" * 10)
        verdict = C.assess_assembly(path, min_size=None)
        assert verdict.included is False
        assert verdict.reason == "not_fasta_no_header"

    def test_oversized_excluded(self, tmp_path):
        path = tmp_path / "big.fna"
        path.write_bytes(b">c\n" + b"ACGT" * 10)
        verdict = C.assess_assembly(path, min_size=None, max_size=10)
        assert verdict.included is False
        assert verdict.reason.startswith("oversized_")

    def test_good_assembly_included(self, tmp_path):
        from papipeline.io.tsv import write_fasta

        path = tmp_path / "good.fna"
        write_fasta(path, [("c1", "ACGT" * 2000)])
        verdict = C.assess_assembly(path, min_size=1000, max_size=10_000_000)
        assert verdict.included is True
        assert verdict.reason == "ok"


class TestAssessCohort:
    def _cohort(self, tmp_path, specs):
        items = []
        for index, (name, payload) in enumerate(specs, start=1):
            path = tmp_path / name
            if payload is not None:
                path.write_bytes(payload)
            items.append(
                C.SelectedAssembly(f"PILOT100_{index:03d}", name, path, name)
            )
        return items

    def test_selection_is_unchanged_by_exclusion(self, tmp_path):
        """Excluding unusable files must not change the selection."""
        items = self._cohort(
            tmp_path,
            [
                ("GCA_1.1.fna", b">c\n" + b"ACGT" * 2000),
                ("GCA_2.2.fna", b""),
                ("GCA_3.3.fna", b">c\n" + b"\xbf" * 40),
            ],
        )
        verdicts, analysable, excluded = C.assess_cohort(
            items, min_size=1000, max_size=10_000_000
        )
        assert len(verdicts) == 3
        assert len(analysable) == 1
        assert len(excluded) == 2
        assert [i.pilot_id for i in analysable] == ["PILOT100_001"]
        assert [i.assembly for i in excluded] == ["GCA_2.2.fna", "GCA_3.3.fna"]

    def test_no_substitution_occurs(self, tmp_path):
        """Excluded assemblies are dropped, never replaced by others."""
        items = self._cohort(
            tmp_path,
            [
                ("GCA_1.1.fna", b">c\n" + b"ACGT" * 2000),
                ("GCA_2.2.fna", b""),
                ("GCA_3.3.fna", b">c\n" + b"ACGT" * 2000),
            ],
        )
        _v, analysable, _e = C.assess_cohort(
            items, min_size=1000, max_size=10_000_000
        )
        assert len(analysable) == 2
        assert all(a.assembly != "GCA_2.2.fna" for a in analysable)

    def test_breakdown_groups_by_reason_type(self, tmp_path):
        items = self._cohort(
            tmp_path,
            [
                ("GCA_1.1.fna", b">c\n" + b"ACGT" * 2000),
                ("GCA_2.2.fna", b""),
                ("GCA_3.3.fna", b">c\n" + b"ACGT" * 10),
            ],
        )
        verdicts, _a, _e = C.assess_cohort(
            items, min_size=1000, max_size=10_000_000
        )
        breakdown = C.exclusion_breakdown(verdicts)
        assert breakdown == {"zero_byte_file": 1, "undersized": 1}

    def test_all_excluded_returns_empty_analysable(self, tmp_path):
        items = self._cohort(tmp_path, [("GCA_1.1.fna", b""), ("GCA_2.2.fna", b"")])
        _v, analysable, excluded = C.assess_cohort(items)
        assert analysable == []
        assert len(excluded) == 2


# --------------------------------------------------------------------------
# frequency counting
# --------------------------------------------------------------------------


#: Annotation the AMRFinderPlus database supplies for the efflux pump
#: components. Categories come from the database's own fields, so a test
#: that wants an Efflux category must supply them.
EFFLUX_ANNOTATION = {"class": "EFFLUX", "subclass": "EFFLUX", "type": "AMR"}


def _member(pid, genes, phenotype="R", annotations=None):
    from papipeline.pilot.pdc_fields import AstField, GenotypeField

    ast = AstField()
    ast.values["imipenem"] = phenotype
    annotations = annotations or {}
    return A.CohortMember(
        pilot_id=pid,
        assembly=pid.replace("_", "."),
        genome_path=f"/tmp/{pid}.fna",
        isolate="X",
        biosample="S",
        bioproject="P",
        pdc_match="MATCHED",
        pdc_present="TRUE",
        ast=ast,
        genotypes=GenotypeField(),
        genes=set(genes),
        gene_annotations={g: dict(annotations.get(g, {})) for g in genes},
    )


class TestGeneFrequencyCounts:
    """Regression: carriers must be counted by isolate, not by category.

    The bug stored the phenotype *label* in the carrier set, so every gene
    was capped at the number of distinct categories (at most 3) and a gene
    present in all 55 isolates was reported as present in 3.
    """

    def _config(self):
        from pathlib import Path as P

        from papipeline.config.loader import load_config

        return load_config(P(__file__).resolve().parents[2] / "config" / "science.yaml")

    def test_counts_isolates_not_categories(self):
        members = [_member(f"P{i:03d}", ["mexA"]) for i in range(10)]
        rows = A.gene_frequency_rows(members, "imipenem", self._config())
        row = [r for r in rows if r["Gene"] == "mexA"][0]
        assert row["n_isolates_detected"] == 10

    def test_gene_in_all_three_categories_counts_all_carriers(self):
        members = (
            [_member("P001", ["mexA"], "R")]
            + [_member("P002", ["mexA"], "S")]
            + [_member("P003", ["mexA"], "I")]
        )
        rows = A.gene_frequency_rows(members, "imipenem", self._config())
        row = [r for r in rows if r["Gene"] == "mexA"][0]
        assert row["n_isolates_detected"] == 3
        assert row["n_R"] == 1 and row["n_S"] == 1 and row["n_I"] == 1

    def test_frequency_is_carriers_over_cohort(self):
        members = [_member(f"P{i:03d}", ["mexA"] if i < 5 else []) for i in range(10)]
        rows = A.gene_frequency_rows(members, "imipenem", self._config())
        row = [r for r in rows if r["Gene"] == "mexA"][0]
        assert row["n_isolates_in_cohort"] == 10
        assert row["frequency"] == pytest.approx(0.5)

    def test_per_category_counts_sum_to_total(self):
        members = (
            [_member("P001", ["mexA"], "R")]
            + [_member("P002", ["mexA"], "S")]
            + [_member("P003", [], "I")]
        )
        rows = A.gene_frequency_rows(members, "imipenem", self._config())
        row = [r for r in rows if r["Gene"] == "mexA"][0]
        assert row["n_R"] + row["n_S"] + row["n_I"] == row["n_isolates_detected"]


class TestMechanismCounts:
    def _config(self):
        from pathlib import Path as P

        from papipeline.config.loader import load_config

        return load_config(P(__file__).resolve().parents[2] / "config" / "science.yaml")

    def test_mechanism_carriers_are_isolates(self):
        members = [
            _member(
                f"P{i:03d}",
                ["mexA", "mexE"],
                annotations={"mexA": EFFLUX_ANNOTATION, "mexE": EFFLUX_ANNOTATION},
            )
            for i in range(8)
        ]
        fams = M.gene_family_table()
        rows = A.mechanism_summary_rows(members, "imipenem", fams)
        row = [r for r in rows if r["Mechanism"] == M.CATEGORY_EFFLUX][0]
        assert row["n_isolates"] == 8

    def test_per_category_mechanism_counts_are_isolates(self):
        members = (
            [
                _member(f"PR{i:03d}", ["mexA"], "R", annotations={"mexA": EFFLUX_ANNOTATION})
                for i in range(6)
            ]
            + [
                _member(f"PS{i:03d}", ["mexA"], "S", annotations={"mexA": EFFLUX_ANNOTATION})
                for i in range(2)
            ]
        )
        fams = M.gene_family_table()
        rows = A.mechanism_summary_rows(members, "imipenem", fams)
        row = [r for r in rows if r["Mechanism"] == M.CATEGORY_EFFLUX][0]
        assert row["n_R"] == 6
        assert row["n_S"] == 2
        assert row["n_R"] + row["n_S"] == row["n_isolates"]


class TestGenePairs:
    def test_pair_counted_once_per_isolate(self):
        members = [_member("P001", ["gA", "gB"]), _member("P002", ["gA", "gB"])]
        rows = A.gene_pair_rows(members, "imipenem")
        assert len(rows) == 1
        assert rows[0]["Gene_A"] == "gA" and rows[0]["Gene_B"] == "gB"
        assert rows[0]["n_isolates_with_both"] == 2

    def test_no_self_pairs(self):
        members = [_member("P001", ["gA", "gA"])]
        assert A.gene_pair_rows(members, "imipenem") == []

    def test_pairs_are_association_only(self):
        rows = A.gene_pair_rows([_member("P001", ["gA", "gB"])], "imipenem")
        assert rows[0]["interpretation_limit"] == (
            "co_occurrence_association_only_not_causal"
        )
