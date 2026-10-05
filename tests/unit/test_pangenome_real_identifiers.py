"""panaroo's REAL gene identifiers are ortholog group ids, not gene names.

**The gap this pins down.** `gene` in the GPA table means different things in
TEST/STUB and in REAL:

  - TEST/STUB - a **literal annotation gene name**, built by
    `stages.pangenome.build_from_annotations` from the annotation table.
  - REAL - a **panaroo ortholog group id**, of the form `group_<n>`, with
    descriptive variants (`group_<n>_<product>`) where the group is not unique.

This is a known, **accepted** TEST-fidelity gap, not a bug and not an oversight.
TEST fixtures use readable names because people read them; REAL uses group ids
because that is what panaroo clusters on. Inventing a name mapping would
fabricate a correspondence the data does not contain.

**The practical consequence, and why this file exists.** A TEST fixture can never
reproduce a REAL gene count or a REAL cluster boundary - only the code path is
comparable, never the numbers. That is only safe if REAL identifiers are pinned
to a shape, because otherwise a future change that quietly substitutes gene names
for group ids in the REAL path would pass every TEST-mode test while producing a
pangenome whose "gene" column means something else entirely - and would feed
pyseer a feature matrix keyed on unstable, non-orthologous labels.

So these tests assert the **shape** of REAL identifiers. They deliberately do not
assert counts: counts belong to a cohort, and a fixture asserting 4,828 core would
be asserting a number that only means something for one real dataset.

Identifiers below were taken from the verified real run (Bakta v1.12.1, DB v6.0
light, 10-isolate smoke cohort) - headers observed in `pan_genome_reference.fa`
were `group_1375`, `group_5706`, `group_5289`, and the `Gene` column of
`gene_presence_absence.csv` carried the same `group_<n>` form.
"""

from __future__ import annotations

import re

from papipeline.adapters import panaroo as adapter
from papipeline.stages import pangenome as stage

#: panaroo's ortholog group id. The numeric part is what makes it an id rather
#: than a name, so the suffix is asserted as digits.
GROUP_ID = re.compile(r"^group_\d+$")

#: Observed for groups whose members share a name but are not one ortholog set.
NON_UNIQUE_GROUP_ID = re.compile(r"^group_\d+_.+$")

REAL_HEADER = (
    "Gene,Non-unique Gene name,Annotation,"
    "PDT000034122.1,PDT000167133.1,PDT000167135.1"
)


class TestRealIdentifiersAreGroupIds:
    def test_a_bare_group_id_matches(self):
        assert GROUP_ID.match("group_1375")
        assert GROUP_ID.match("group_1")
        assert GROUP_ID.match("group_4828")

    def test_a_non_unique_group_id_is_still_recognisable(self):
        """panaroo appends a product name when a group is not unique.

        These must still be treated as the same ortholog group, and must not be
        mistaken for a gene name. A strict `group_<digits>` match alone would
        reject them.
        """
        assert NON_UNIQUE_GROUP_ID.match("group_1375_C4-dicarboxylate transporter")

    def test_a_gene_name_is_not_a_valid_real_identifier(self):
        """The distinction the whole file exists to protect.

        If a bare annotation gene name satisfies the identifier check, then a
        regression that swapped the REAL source for the TEST source would go
        unnoticed - and the two are indistinguishable downstream.
        """
        for name in ("abcA", "OprD", "PA2020", "ftsZ", "recA"):
            assert not GROUP_ID.match(name)
            assert not NON_UNIQUE_GROUP_ID.match(name)


class TestParserPreservesRealIdentifiers:
    def test_group_ids_survive_parsing_unchanged(self, tmp_path):
        """Identifiers must not be normalised, split or renamed on the way in.

        Anything that rewrites `group_1375` into something else breaks the link
        between the GPA table and panaroo's own `.Rtab`, which is what pyseer is
        eventually handed.
        """
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(
            "\n".join(
                [
                    REAL_HEADER,
                    "group_1375,,,1,1,0",
                    "group_5706_ABC transporter,,,1,0,0",
                    "group_5289,,,0,0,0",
                ]
            )
            + "\n"
        )
        samples, presence = adapter.parse_presence_csv(path)

        assert set(presence) == {
            "group_1375",
            "group_5706_ABC transporter",
            "group_5289",
        }
        for gene in presence:
            assert GROUP_ID.match(gene) or NON_UNIQUE_GROUP_ID.match(gene), (
                f"{gene!r} came out of the REAL parse path but is not a panaroo "
                f"group id; something is rewriting identifiers"
            )

    def test_a_gene_absent_from_every_isolate_is_still_present_as_a_key(self, tmp_path):
        """`group_5289` is 0/0/0 above.

        It must remain a key with an empty carrier set. Dropping it would make
        "absent from every isolate" indistinguishable from "never looked for",
        and the first is a real biological observation.
        """
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(
            "\n".join([REAL_HEADER, "group_5289,,,0,0,0"]) + "\n"
        )
        _, presence = adapter.parse_presence_csv(path)
        assert "group_5289" in presence
        assert presence["group_5289"] == set()

    def test_the_partition_tolerates_group_ids_verbatim(self, tmp_path):
        """Full REAL-shaped path: parse, partition, and the declared columns.

        Exercises the whole chain with realistic group ids so the fixture and the
        tool agree on shape, and asserts the cohort size comes through - which is
        the number that decides core versus accessory.
        """
        path = tmp_path / adapter.PRESENCE_CSV
        path.write_text(
            "\n".join(
                [
                    REAL_HEADER,
                    "group_1375,,,1,1,1",
                    "group_5706,,,1,1,0",
                    "group_5289,,,0,0,0",
                ]
            )
            + "\n"
        )
        cohort = [
            "PDT000034122.1",
            "PDT000167133.1",
            "PDT000167135.1",
        ]
        samples, presence = adapter.build_pangenome_from_panaroo(path.parent, cohort)
        result = stage.partition(samples, presence)

        assert samples == tuple(cohort)
        assert result.core == ("group_1375",), "one group is in all three isolates"
        assert set(result.accessory) == {"group_5706", "group_5289"}
        assert result.n_samples == 3

    def test_the_round_trip_preserves_group_ids_as_keys(self, tmp_path):
        """REAL identifiers must survive write -> read unchanged.

        A normalisation in either direction would leave the GPA table keyed on
        something other than what panaroo emitted, and the mismatch would only
        show up as a silently reduced feature count at GWAS time.
        """
        samples = ("PDT000034122.1", "PDT000167133.1", "PDT000167135.1")
        presence = {
            "group_1375": set(samples),
            "group_5706": set(samples[:2]),
        }
        written = stage.write_outputs(stage.partition(samples, presence), tmp_path)
        back_samples, back_presence = stage.read_gene_presence_absence(
            written["gene_presence_absence"]
        )
        assert back_samples == samples
        assert back_presence == presence
        for gene in back_presence:
            assert GROUP_ID.match(gene), f"{gene!r} is not a panaroo group id"
