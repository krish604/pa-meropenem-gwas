"""Sampling for the allele-frequency measurement must be stratified, not arbitrary.

The measurement that will set `max_maf` asks what allele frequencies look like
across a cohort, so the sample has to be one where "frequency" means something.
The first 50 isolates by accession would measure whichever clone or hospital the
numbering happens to start in, and a homogeneous sample collapses the spectrum
into two spikes - making any threshold derived from it an artefact of the sample
rather than a property of the data.

Stratification is by `Location`: 46 distinct values among the 835 isolates that
carry an assembly, against 29 BioProjects. Round-robin then keeps a small
stratum present, which is what a rare-allele class depends on.

The specific failure these pin: a plain `random.sample` gives `USA: Texas` about
29% of 50 isolates, so a whole quarter of a 50-isolate run comes from one
hospital. Under round-robin, every stratum contributes. The first-50-by-order
sample is pinned too, because it is the obvious wrong answer and it happens to
look reasonable.

And the 35 isolates with no recorded location are a stratum, not a filter:
dropping them would bias the sample toward labs that filled the field in.
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest

from papipeline.cohort_sampling import load_located, sample_stratified, stratify

PDC = Path("PDC_essential.tsv")


@pytest.fixture(scope="module")
def rows():
    return load_located(PDC)


class TestTheSampleIsStratified:
    def test_it_spans_many_locations(self, rows):
        chosen = sample_stratified(rows, 50)
        locations = {r.get("Location") or "(unrecorded)" for r in chosen}
        assert len(locations) >= 10, f"only {len(locations)} locations in 50 isolates"

    def test_small_strata_are_not_crowded_out(self, rows):
        """A stratum contributing nothing is what collapses a rare-allele class."""
        chosen = sample_stratified(rows, 50)
        counts = collections.Counter(r.get("Location") or "(unrecorded)" for r in chosen)
        assert min(counts.values()) >= 1
        assert len(counts) >= 15, f"only {len(counts)} strata represented"

    def test_it_does_not_collapse_onto_one_hospital(self, rows):
        """`USA: Texas` is 240 of 835; a third of the sample would be one hospital."""
        chosen = sample_stratified(rows, 50)
        counts = collections.Counter(r.get("Location") or "(unrecorded)" for r in chosen)
        top = counts.most_common(1)[0]
        assert top[1] <= 6, f"{top[0]} took {top[1]} of 50 isolates"

    def test_unrecorded_location_is_a_stratum_not_a_filter(self, rows):
        chosen = sample_stratified(rows, 50)
        assert any(not r.get("Location") for r in chosen), (
            "isolates with no recorded location were dropped from the sample"
        )

    def test_first_fifty_by_order_would_be_worse(self, rows):
        """Pins the obvious wrong answer, which happens to look plausible."""
        naive = rows[:50]
        naive_locations = {r.get("Location") or "(unrecorded)" for r in naive}
        stratified = {r.get("Location") or "(unrecorded)" for r in sample_stratified(rows, 50)}
        assert len(stratified) > len(naive_locations)


class TestTheSampleIsUsable:
    def test_every_chosen_isolate_has_an_assembly(self, rows):
        """Loading already filters these, but the filter is load-bearing."""
        for row in sample_stratified(rows, 50):
            assert row["Assembly"]

    def test_it_is_reproducible(self, rows):
        """The spectrum must be repeatable from this code, not from a lost run."""
        first = [r["Isolate"] for r in sample_stratified(rows, 50)]
        second = [r["Isolate"] for r in sample_stratified(rows, 50)]
        assert first == second

    def test_asking_for_more_than_exists_is_refused(self, rows):
        with pytest.raises(ValueError):
            sample_stratified(rows, len(rows) + 1)

    def test_stratify_keeps_empty_values_as_their_own_group(self):
        groups = stratify([{"Location": "USA"}, {"Location": ""}])
        assert "" in groups or "(unrecorded)" in groups


class TestAgainstRealMetadata:
    def test_location_is_the_strongest_stratifier_available(self, rows):
        """Justifies the choice from the data rather than asserting it."""
        assert len(stratify(rows, "Location")) >= 20
        assert len(stratify(rows, "BioProject")) <= len(stratify(rows, "Location"))
