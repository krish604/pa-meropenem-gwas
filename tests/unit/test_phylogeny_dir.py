"""Stage 10's working directory is an *output* location, not an assembly root.

`run.py` derived it as `data_root / "phylogeny"`. That conflated two different
things, and a sweep of every stage for the assembly-root defect turned it up:

* under TEST, `data_root` is `test_data`, so it happened to land on the
  committed fixture directory and looked correct;
* under REAL it wrote into the **read-only** `data/` tree;
* and once `data_root` became `assembly_root` (83d6aba, efb2311) to bound stages
  2/3/4/6, it wrote into `db/smoke_genomes/phylogeny` - inside the curated
  assembly directory.

So the path was wrong in both modes, and my own fix moved it somewhere new.
`assembly_root` answers "where do this run's assemblies live", which is the wrong
question for a directory the stage writes trees into.

The fix is one accessor that answers the right question per mode: the committed
fixture under TEST, an intermediate output directory under REAL. TEST is
unchanged, so the fixture contract does not move.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def laptop():
    return load_config(REPO / "config" / "science.yaml", machine="laptop")


@pytest.fixture(scope="module")
def smoke():
    return load_config(
        REPO / "config" / "science.yaml",
        machine=REPO / "config" / "machines" / "smoke.yaml",
    )


class TestTestModeIsUnchanged:
    def test_it_is_still_the_committed_fixture_directory(self, laptop):
        """The fixture is the TEST contract; moving it would break every
        existing consumer for no benefit."""
        assert laptop.phylogeny_dir(RunMode.TEST) == (
            REPO / "test_data" / "phylogeny"
        )


class TestRealModeWritesSomewhereWritable:
    def test_not_into_the_read_only_data_tree(self, laptop):
        directory = laptop.phylogeny_dir(RunMode.REAL)
        assert REPO / "data" not in directory.parents, (
            f"{directory} is inside data/, which is read-only; a REAL stage "
            "writing there would fail on a clean checkout"
        )

    def test_not_into_the_curated_assemblies(self, smoke):
        """The bug this sweep found.

        `data_root` became `assembly_root` so stages 2/3/4/6 would stop reading
        602 unauthorised assemblies. Stage 10 then inherited
        `db/smoke_genomes/phylogeny` - a directory of curated FASTA files, which
        is not where a tree belongs.
        """
        directory = smoke.phylogeny_dir(RunMode.REAL)
        assert smoke.machine.smoke_genome_dir() not in directory.parents, (
            f"{directory} is inside the curated genome directory; the stage "
            "writes trees and crosswalks, not sequences"
        )

    def test_it_is_under_the_intermediate_root(self, laptop):
        directory = laptop.phylogeny_dir(RunMode.REAL)
        assert directory.is_relative_to(laptop.intermediate_root(RunMode.REAL))


class TestTheOrchestratorUsesIt:
    def test_run_py_asks_for_the_directory_rather_than_deriving_one(self):
        """Pinned because the defect was a one-line derivation, and a one-line
        derivation is easy to reintroduce."""
        source = (REPO / "papipeline" / "run.py").read_text(encoding="utf-8")
        assert 'phylo_dir = data_root / "phylogeny"' not in source, (
            "run.py still derives the phylogeny directory from the assembly "
            "root, which writes stage output into a directory of assemblies"
        )
        assert "config.phylogeny_dir(resolved)" in source
