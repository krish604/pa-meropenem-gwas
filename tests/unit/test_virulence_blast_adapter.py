"""Raw blastn virulence screening against VFDB.

Three of these tests exist because the corresponding bug was *found the hard
way* against the real 4,592-sequence database, and in every case the wrong
behaviour looked like a valid scientific result rather than like a failure.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from papipeline.adapters import virulence as vf
from papipeline.errors import ToolExecutionError

#: Real VFDB header shapes, taken from `db/vfdb/sequences`. The first is the
#: common `[factor (VF) - category (VFC)]`; the second carries a hyphen *inside*
#: the factor name, which is what broke a naive split.
HEADER_PILUS = (
    ">vfdb~~~fimV~~~NP_251805~~~ (fimV) a polar peptidoglycan-binding protein "
    "involved in type IV pilus assembly "
    "[Type IV pili (VF0082) - Adherence (VFC0001)] "
    "[Pseudomonas aeruginosa PAO1]"
)
HEADER_HSI = (
    ">vfdb~~~PA1663~~~NP_250354~~~ (PA1663) transcriptional regulator "
    "[HSI-2 (VF0943) - Effector delivery system (VFC0086)] "
    "[Pseudomonas aeruginosa PAO1]"
)


def _fasta(tmp_path: Path, *headers: str) -> Path:
    path = tmp_path / "sequences"
    lines = []
    for header in headers:
        lines.append(header)
        lines.append("ATGTGTGCGCTGGATCGTAGAGAAAGGCCACTTAACAGTCAATCTGTAA")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestHeaderParsing:
    def test_it_reads_factor_gene_and_category(self, tmp_path):
        entries = vf.parse_vfdb_headers(_fasta(tmp_path, HEADER_PILUS))
        entry = entries["vfdb~~~fimV~~~NP_251805~~~"]
        assert entry.factor == "Type IV pili"
        assert entry.category == "Adherence"
        assert entry.gene == "fimV"
        assert entry.vf_id == "VF0082"
        assert entry.category_id == "VFC0001"

    def test_a_hyphen_inside_the_factor_name_does_not_split_it(self, tmp_path):
        """`HSI-2` is one factor name, not two fields.

        Splitting on `\\s*-\\s*` matched the hyphen inside `HSI-2` - `\\s*`
        admits zero whitespace - and produced gene `HSI` with category
        `2 (VF0943) - Effector delivery system (VFC0086)`. That is silently
        wrong: it invents a category in a scientific output and nothing fails.
        """
        entries = vf.parse_vfdb_headers(_fasta(tmp_path, HEADER_HSI))
        entry = entries["vfdb~~~PA1663~~~NP_250354~~~"]
        assert entry.factor == "HSI-2"
        assert entry.category == "Effector delivery system"
        assert entry.gene == "PA1663"

    def test_gene_comes_from_the_accession_not_the_prose(self, tmp_path):
        """The second prose field is `(locus)` for some entries, a product word
        for others. Reading it yielded `unknown` / `transcriptional` as a gene."""
        header = (
            ">vfdb~~~AAA92657~~~AAA92657~~~ (AAA92657) unknown protein "
            "[TraJ (VF0241) - Invasion (VFC0083)] [Escherichia coli]"
        )
        entry = vf.parse_vfdb_headers(_fasta(tmp_path, header))[
            "vfdb~~~AAA92657~~~AAA92657~~~"
        ]
        assert entry.gene == "AAA92657"

    def test_an_unannotated_header_is_skipped_and_counted(self, tmp_path):
        """A header with no `(VFnnnn)` cannot be named, so it is not screened."""
        entries = vf.parse_vfdb_headers(
            _fasta(tmp_path, HEADER_PILUS, ">vfdb~~~x~~~~~~ (x) unnamed [??]")
        )
        assert set(entries) == {"vfdb~~~fimV~~~NP_251805~~~"}

    def test_a_file_with_no_parsable_headers_is_an_error(self, tmp_path):
        """Rather than screening nothing and reporting zero carriage."""
        with pytest.raises(vf.VirulencePreflightError):
            vf.parse_vfdb_headers(_fasta(tmp_path, ">vfdb~~~a~~~~~~ (a) no block"))


class TestBlastCommand:
    def test_it_requests_subject_length_not_query_length(self, tmp_path):
        """`qlen` is the *genome*, so coverage computed from it is ~0.1%.

        blastn then exits 0 with zero findings, which reads as "this isolate
        carries no virulence factors" rather than as a broken screen. Every hit
        was silently rejected by a 90% coverage threshold.
        """
        command = vf.blastn_command(
            program="blastn",
            genome=tmp_path / "g.fna",
            database=tmp_path / "db",
            threads=4,
            evalue=1e-5,
        )
        outfmt = command[command.index("-outfmt") + 1]
        assert "slen" in outfmt
        assert "qlen" not in outfmt

    def test_threads_come_from_the_caller_not_a_literal(self, tmp_path):
        command = vf.blastn_command(
            program="blastn",
            genome=tmp_path / "g.fna",
            database=tmp_path / "db",
            threads=7,
            evalue=1e-5,
        )
        assert command[command.index("-num_threads") + 1] == "7"


class TestHitSelection:
    def _hit(self, accession, identity, length, subject_length):
        return vf.BlastHit(
            accession=accession,
            identity_pct=identity,
            alignment_length=length,
            subject_length=subject_length,
            evalue=0.0,
        )

    def _entries(self, *accessions, factor="Type IV pili", category="Adherence"):
        return {
            a: vf.VfdbEntry(
                accession=a, gene=a.split("~~~")[1], factor=factor,
                category=category, vf_id="VF0082", category_id="VFC0001",
                description="",
            )
            for a in accessions
        }

    def test_several_sequences_under_one_factor_yield_one_hit(self):
        """VFDB bundles fimV, xcpQ, algG under `Type IV pili (VF0082)`.

        Counting sequences instead of factors inflates the count by ~10x. An
        earlier version truncated the whole list to `max_hits_per_factor` and
        reported exactly *one* factor for every genome.
        """
        entries = self._entries("vfdb~~~fimV~~~", "vfdb~~~xcpQ~~~", "vfdb~~~pilA~~~")
        hits = [
            self._hit("vfdb~~~fimV~~~", 99.0, 1000, 1000),
            self._hit("vfdb~~~xcpQ~~~", 100.0, 900, 900),
            self._hit("vfdb~~~pilA~~~", 97.0, 800, 800),
        ]
        chosen = vf.select_hits(
            hits, entries, min_identity_pct=90.0, min_coverage_pct=90.0
        )
        assert len(chosen) == 1
        assert chosen[0].accession == "vfdb~~~xcpQ~~~", "best identity should win"

    def test_distinct_factors_are_all_retained(self):
        entries = {
            **self._entries("vfdb~~~a~~~", factor="Adherence"),
            **self._entries("vfdb~~~b~~~", factor="Biofilm"),
        }
        chosen = vf.select_hits(
            [self._hit("vfdb~~~a~~~", 99.0, 100, 100),
             self._hit("vfdb~~~b~~~", 99.0, 100, 100)],
            entries, min_identity_pct=90.0, min_coverage_pct=90.0,
        )
        assert len(chosen) == 2

    def test_both_thresholds_are_enforced(self):
        entries = self._entries("vfdb~~~a~~~")
        low_identity = [self._hit("vfdb~~~a~~~", 80.0, 100, 100)]
        low_coverage = [self._hit("vfdb~~~a~~~", 99.0, 50, 1000)]
        assert vf.select_hits(
            low_identity, entries, min_identity_pct=90.0, min_coverage_pct=90.0
        ) == []
        assert vf.select_hits(
            low_coverage, entries, min_identity_pct=90.0, min_coverage_pct=90.0
        ) == []

    def test_coverage_is_percent_of_the_factor(self):
        hit = self._hit("vfdb~~~a~~~", 99.0, 500, 1000)
        assert hit.coverage_pct == 50.0


class TestTableParsing:
    def test_a_malformed_row_is_dropped_not_fatal(self):
        """blast emitting an unexpected column must not lose the isolate."""
        table = (
            "qf\tvfdb~~~a~~~\ttitle\t99.0\t100\t100\t0.0\n"
            "garbage\n"
            "qf\tvfdb~~~b~~~\ttitle\t97.0\t90\t100\t0.0\n"
        )
        hits = vf.parse_blast_table(table)
        assert [h.accession for h in hits] == ["vfdb~~~a~~~", "vfdb~~~b~~~"]


class TestScreenIsolate:
    def _result(self, returncode=0, stdout="", stderr=""):
        from papipeline.adapters.external import CommandResult

        return CommandResult(
            command=["blastn"], returncode=returncode, stdout=stdout,
            stderr=stderr, duration_s=0.1,
        )

    def _runner(self, result):
        def runner(command):
            return result
        return runner

    def test_a_nonzero_exit_is_raised_not_returned(self, tmp_path):
        entries = {}
        with pytest.raises(ToolExecutionError):
            vf.screen_isolate(
                genome=tmp_path / "g.fna", entries=entries,
                database=tmp_path / "db", program="blastn", threads=1,
                evalue=1e-5, min_identity_pct=90.0, min_coverage_pct=90.0,
                database_name="VFDB", database_version="sha256:abc",
                runner=self._runner(self._result(returncode=1, stderr="boom")),
            )

    def test_a_subject_missing_from_the_header_set_is_reported(self, tmp_path):
        """An index out of step with the sequences file under-counts silently."""
        entries = {}
        table = "qf\tvfdb~~~ghost~~~~\ttitle\t99.0\t100\t100\t0.0\n"
        with pytest.raises(vf.VirulencePreflightError):
            vf.screen_isolate(
                genome=tmp_path / "g.fna", entries=entries,
                database=tmp_path / "db", program="blastn", threads=1,
                evalue=1e-5, min_identity_pct=90.0, min_coverage_pct=90.0,
                database_name="VFDB", database_version="sha256:abc",
                runner=self._runner(self._result(stdout=table)),
            )

    def test_a_carriage_row_carries_the_annotation(self, tmp_path):
        entries = {
            "vfdb~~~fimV~~~": vf.VfdbEntry(
                accession="vfdb~~~fimV~~~", gene="fimV", factor="Type IV pili",
                category="Adherence", vf_id="VF0082", category_id="VFC0001",
                description="",
            )
        }
        table = "qf\tvfdb~~~fimV~~~\ttitle\t99.9\t2760\t2760\t0.0\n"
        rows = vf.screen_isolate(
            genome=tmp_path / "g.fna", entries=entries,
            database=tmp_path / "db", program="blastn", threads=1,
            evalue=1e-5, min_identity_pct=90.0, min_coverage_pct=90.0,
            database_name="VFDB", database_version="sha256:abc",
            runner=self._runner(self._result(stdout=table)),
        )
        assert len(rows) == 1
        assert rows[0]["virulence_factor"] == "Type IV pili"
        assert rows[0]["category"] == "Adherence"
        assert rows[0]["identity_pct"] == 99.9
        assert rows[0]["coverage_pct"] == 100.0


class TestDatabaseHandling:
    def test_the_version_is_the_content_digest(self, tmp_path):
        """Not a hard-coded string: VFDB is re-released, and a stale pin
        describes bytes that are no longer on disk."""
        path = tmp_path / "sequences"
        path.write_text("ACGT\n", encoding="utf-8")
        version = vf.database_version_on_disk(path)
        assert version.startswith("sha256:")
        path.write_text("ACGA\n", encoding="utf-8")
        assert vf.database_version_on_disk(path) != version

    def test_an_existing_index_is_reused(self, tmp_path):
        """Rebuilding would rewrite a database the operator pinned."""
        sequences = tmp_path / "sequences"
        sequences.write_text(">a\nACGT\n", encoding="utf-8")
        vf.blast_db_path(sequences).write_text("index", encoding="utf-8")

        calls = []

        def runner(command):
            calls.append(command)
            raise AssertionError("makeblastdb must not run when an index exists")

        assert vf.ensure_blast_db(sequences, runner=runner).name == "sequences"
        assert calls == []

    def test_the_returned_database_name_is_the_stem_not_the_index(self, tmp_path):
        """`blastn -db x.nin` fails with 'No alias or index file found'."""
        sequences = tmp_path / "sequences"
        sequences.write_text(">a\nACGT\n", encoding="utf-8")
        vf.blast_db_path(sequences).write_text("index", encoding="utf-8")
        assert str(vf.ensure_blast_db(sequences)).endswith("sequences")