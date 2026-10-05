"""The tblastn adapter: what it writes, what it refuses, and what it records.

**Real code, synthetic inputs, injected runners.** Every test here except the
last class drives `papipeline.adapters.tblastn.run_tblastn` itself - the command
it builds, the refusals it raises and the file it writes - with a fake runner
standing in for blast, per R12. The final class runs the installed blast 2.17
and is skipped where it is absent.

**Why the table is validated before it is written.** `run_tblastn` parses its own
output with `oprd_locus.parse_tblastn_table` and refuses to put the file in place
unless that parser accepts it, because `structural_call_from_paths` re-parses the
same file. A malformed table written here would surface three stages later as a
`PipelineError` from a decision function, naming no tool. The tests below pin
both halves: it raises, AND it leaves nothing at `out_tsv` for a later stage to
misread.

**The query is searched as given.** `run_tblastn` does not pad, trim or
re-interpret the query FASTA. That is pinned directly, because the round-12 P2
comparison turned on exactly one character of that file and an adapter that
quietly normalised it would have hidden the cause.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.adapters import oprd_locus as oprd
from papipeline.adapters import tblastn as tb
from papipeline.adapters.external import CommandResult
from papipeline.errors import (
    PipelineError,
    ToolExecutionError,
    ToolNotAvailableError,
)

#: One well-formed HSP row in the declared outfmt: 14 tab-separated columns,
#: the last two the aligned sequences that place the indels.
GOOD_ROW = (
    "oprD\tcontig_1\t99.500\t300\t300\t5000\t1\t300\t100\t399\t0.0\t600.0"
    "\tMKVMKWSAIA\tMKVMKWSAIA"
)
SECOND_ROW = (
    "oprD\tcontig_1\t42.000\t120\t300\t5000\t10\t129\t900\t1019\t1.0e-20\t90.0"
    "\tMKVMKW\tMKVMKW"
)


def _result(command, returncode=0, stdout="", stderr=""):
    return CommandResult(
        command=list(command), returncode=returncode, stdout=stdout,
        stderr=stderr, duration_s=0.01,
    )


def _fake_runner(table=GOOD_ROW + "\n" + SECOND_ROW + "\n", *,
                 index_returncode=0, write_index=True,
                 search_returncode=0, stderr=""):
    """A runner standing in for both tools, keyed on argv[0]'s basename.

    ``makeblastdb`` is emulated as faithfully as the adapter needs: it creates the
    ``.nin`` the adapter then checks for, because refusing on a missing index is
    real behaviour that a fake that skipped it would never exercise.
    """
    calls = []

    def runner(command):
        command = list(command)
        calls.append(command)
        tool = Path(command[0]).name
        if tool == "makeblastdb":
            if write_index and index_returncode == 0:
                prefix = command[command.index("-out") + 1]
                Path(f"{prefix}.nin").write_text("index", encoding="utf-8")
            return _result(command, returncode=index_returncode,
                           stdout="Building a new DB\n")
        return _result(command, returncode=search_returncode, stdout=table,
                       stderr=stderr)

    runner.calls = calls
    return runner


@pytest.fixture
def inputs(tmp_path):
    query = tmp_path / "reference_protein.faa"
    query.write_text(">PA0958 oprD\n" + "MKV" * 40 + "\n", encoding="utf-8")
    assembly = tmp_path / "genome.fna"
    assembly.write_text(">contig_1\n" + "ACGT" * 200 + "\n", encoding="utf-8")
    return query, assembly, tmp_path / "work" / "tblastn.tsv"


class TestTheTableItWrites:
    def test_the_table_is_written_where_it_was_asked(self, inputs):
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        assert out_tsv.is_file()
        assert info.out_tsv == out_tsv
        assert out_tsv.read_text(encoding="utf-8") == (
            GOOD_ROW + "\n" + SECOND_ROW + "\n"
        )

    def test_the_table_is_exactly_what_the_consumer_parses(self, inputs):
        """`out_tsv` is only correct if `structural_call_from_paths` reads it.

        The consumer re-parses the file with `parse_tblastn_table`, so the test
        is that same parser, on the bytes this adapter wrote - not a count of
        lines.
        """
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        hits = oprd.parse_tblastn_table(out_tsv.read_text(encoding="utf-8"),
                                        query_codon_length=120)
        assert len(hits) == info.n_hits == 2
        assert [h.subject_id for h in hits] == ["contig_1", "contig_1"]
        assert hits[0].query_aligned == "MKVMKWSAIA"

    def test_the_row_count_is_counted_and_recorded(self, inputs):
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        assert info.n_hits == 2

    def test_an_empty_table_is_a_valid_result_and_is_written(self, inputs):
        """Zero hits is the one shape that may later be read as `absent`.

        blast exits 0 and prints nothing when it finds nothing. If this adapter
        treated empty output as a failure, an isolate with no oprD locus would
        be reported as a tool error; if it refused to write the file, the
        recorded empty search could never reach the `absent` branch.
        """
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner(table=""))
        assert out_tsv.is_file()
        assert out_tsv.read_text(encoding="utf-8") == ""
        assert info.n_hits == 0
        assert info.search.n_hits == 0
        assert info.search.is_recorded

    def test_no_partial_file_is_left_behind(self, inputs):
        query, assembly, out_tsv = inputs
        tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(), threads=1,
                       runner=_fake_runner())
        assert not list(out_tsv.parent.glob(".*partial"))

    def test_the_query_is_searched_exactly_as_given(self, inputs):
        """No padding, no trimming: `qlen` is the query's own residue count.

        This is the round-12 P2 finding made permanent. `oprd_locus.reference_protein`
        returns 443 residues for PA0958 while `measure_locus_structure` compares
        the alignment's high query coordinate against `len(reference_cds)//3`
        = 444, a count that INCLUDES the terminator - so a query without a
        trailing `*` can never satisfy `reaches_reference_terminator`. An
        adapter that normalised the query would have hidden that; this pins that
        it does not.
        """
        query, assembly, out_tsv = inputs
        residues = 7  # "MKV" * 2 + one character
        query.write_text(">x\nMKVMKVM\n", encoding="utf-8")
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner(
                                  table="x\tc\t99.0\t7\t7\t5000\t1\t7\t1\t21"
                                         "\t0.0\t50.0\tMKVMKVM\tMKVMKVM\n"))
        first = oprd.parse_tblastn_table(out_tsv.read_text(encoding="utf-8"),
                                        query_codon_length=residues)[0]
        assert first.query_length == residues
        assert query.read_text(encoding="utf-8") == ">x\nMKVMKVM\n"

    def test_the_database_is_built_beside_the_table_not_beside_the_assembly(
        self, inputs
    ):
        """`-db` cannot read a FASTA, so the adapter formats one. Beside the
        table, because the assembly is an input under `data/`."""
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        assert info.database == str(tb.database_prefix(out_tsv))
        assert Path(info.database).parent.parent == out_tsv.parent
        assert not list(assembly.parent.glob(f"{assembly.name}.n*"))
        index_command = info.makeblastdb_command
        assert "-dbtype" in index_command
        assert index_command[index_command.index("-dbtype") + 1] == "nucl"
        assert index_command[index_command.index("-out") + 1] == info.database


class TestTheCommand:
    def test_the_search_parameters_are_the_recorded_ones(self, inputs):
        """Including the two that are wrong by default.

        `-db_gencode 11` because blast's default is the standard code, which
        invents stops on a bacterial locus; `-max_target_seqs 5000` because 0 is
        rejected by blast 2.17 as legacy syntax.
        """
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        command = info.command
        assert command[command.index("-db_gencode") + 1] == "11"
        assert command[command.index("-max_target_seqs") + 1] == "5000"
        assert float(command[command.index("-evalue") + 1]) == 1.0e-3
        assert command[command.index("-db") + 1] == info.database
        assert command[command.index("-query") + 1] == str(query)

    def test_the_outfmt_is_one_argv_element_with_the_declared_fields(self, inputs):
        """The recorded `.cmd` files render it as `6\\ qseqid\\ sseqid\\ ...`.

        That is `printf '%q'` escaping spaces inside one argv element, not a
        different outfmt. This asserts the element equals the module's declared
        field list, which is what the stored tables actually contain.
        """
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        element = info.command[info.command.index("-outfmt") + 1]
        assert element == oprd.TBLASTN_OUTFMT
        assert "\\" not in element
        assert element.startswith("6 qseqid ")
        assert element.endswith(" qseq sseq")
        assert len(element.split()) - 1 == len(oprd.TBLASTN_OUTFMT_FIELDS)

    def test_threads_reach_the_command_and_the_record(self, inputs):
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=7,
                              runner=_fake_runner())
        assert info.command[info.command.index("-num_threads") + 1] == "7"
        assert ("num_threads", "7") in info.parameters

    def test_the_search_is_recorded_so_a_negative_is_auditable(self, inputs):
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=2,
                              runner=_fake_runner())
        search = info.search
        assert search.is_recorded
        assert search.program == "tblastn"
        assert search.query == str(query)
        assert search.database == info.database
        assert search.n_hits == info.n_hits
        assert "tblastn" in search.command
        assert dict(search.parameters)["db_gencode"] == "11"

    def test_the_version_and_the_digest_are_recorded(self, inputs):
        """Provenance on every result, not in a log nobody reads."""
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        import hashlib
        assert info.table_sha256 == hashlib.sha256(
            out_tsv.read_bytes()
        ).hexdigest()
        assert isinstance(info.tool_version, str) and info.tool_version
        assert isinstance(info.makeblastdb_version, str)
        row = info.as_row()
        assert row["tool"] == "tblastn"
        assert row["n_hits"] == 2
        assert row["outfmt"] == oprd.TBLASTN_OUTFMT
        assert row["parameters"].startswith("outfmt=6 qseqid")

    def test_nothing_about_the_machine_is_hard_coded(self, inputs):
        """No absolute path is baked in: the query and assembly come in as
        arguments and the executable is resolved, never spelled out."""
        query, assembly, out_tsv = inputs
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=1,
                              runner=_fake_runner())
        assert str(query) in info.command_line
        assert str(out_tsv.parent) in str(
            info.makeblastdb_command[info.makeblastdb_command.index("-out") + 1]
        )
        assert not any(
            part.startswith("/Users/") for part in info.command
        )


class TestTheRefusals:
    def test_a_missing_tblastn_is_refused_by_name(self, inputs, monkeypatch):
        """The refusal names the tool, because two tools are involved and the
        fix for each is different."""
        query, assembly, out_tsv = inputs
        monkeypatch.setattr(tb.shutil, "which", lambda name: None)
        with pytest.raises(ToolNotAvailableError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=_fake_runner())
        assert "tblastn" in str(caught.value)
        assert not out_tsv.exists()

    def test_a_missing_makeblastdb_is_refused_by_name(self, inputs, monkeypatch):
        """`-db` cannot read a FASTA, so a missing formatter is as fatal as a
        missing searcher and must not be reported as a search problem."""
        query, assembly, out_tsv = inputs
        monkeypatch.setattr(
            tb.shutil, "which",
            lambda name: f"/opt/blast/bin/{name}" if name == "tblastn" else None,
        )
        with pytest.raises(ToolNotAvailableError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=_fake_runner())
        assert "makeblastdb" in str(caught.value)
        assert caught.value.context["tool"] == "makeblastdb"

    def test_the_search_directories_are_searched_and_named(self, inputs,
                                                           monkeypatch):
        query, assembly, out_tsv = inputs
        monkeypatch.setattr(tb.shutil, "which", lambda name: None)
        found = tmp_path_tools = inputs[2].parent / "tools"
        found.mkdir(parents=True)
        for name in ("tblastn", "makeblastdb"):
            tool = found / name
            tool.write_text("#!/bin/sh\n", encoding="utf-8")
            tool.chmod(0o755)
        assert tb.resolve_tool("tblastn", [found]) == found / "tblastn"
        info = tb.run_tblastn(query, inputs[1], out_tsv,
                              tool_search_dirs=[found], threads=1,
                              runner=_fake_runner())
        assert info.executable == found / "tblastn"
        assert info.command[0] == str(found / "tblastn")
        assert tmp_path_tools.is_dir()

    def test_a_directory_that_does_not_exist_is_not_a_dead_end(self,
                                                               monkeypatch):
        monkeypatch.setattr(tb.shutil, "which", lambda name: None)
        with pytest.raises(ToolNotAvailableError) as caught:
            tb.resolve_tool("tblastn", [Path("/nonexistent/tools")])
        assert "tblastn" in str(caught.value)

    def test_a_missing_query_is_refused_before_any_tool_runs(self, inputs,
                                                             monkeypatch):
        query, assembly, out_tsv = inputs
        runner = _fake_runner()
        query.unlink()
        with pytest.raises(PipelineError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=runner)
        assert str(query) in str(caught.value)
        assert runner.calls == [], "no tool may run without its input"

    def test_a_missing_assembly_is_refused_before_any_tool_runs(self, inputs):
        query, assembly, out_tsv = inputs
        runner = _fake_runner()
        assembly.unlink()
        with pytest.raises(PipelineError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=runner)
        assert str(assembly) in str(caught.value)
        assert runner.calls == []

    def test_zero_threads_is_refused_rather_than_passed_on(self, inputs):
        """blast rejects `-num_threads 0`; passing it on would either fail the
        tool or run a search that ignored the allocation."""
        query, assembly, out_tsv = inputs
        for value in (0, -1):
            with pytest.raises(PipelineError) as caught:
                tb.run_tblastn(query, assembly, out_tsv,
                               tool_search_dirs=(), threads=value,
                               runner=_fake_runner())
            assert "threads" in str(caught.value)

    def test_a_malformed_table_is_refused_and_nothing_is_written(self, inputs):
        """A row the consumer's parser would refuse must not reach disk.

        `parse_tblastn_table` accepts 14 columns or 12 without the aligned
        sequences, and refuses anything else precisely because a different
        `-outfmt` would place the indels - and so the reading frame - silently
        wrong.
        """
        query, assembly, out_tsv = inputs
        bad = "oprD\tcontig_1\t99.5\n"
        with pytest.raises(PipelineError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=_fake_runner(table=bad))
        assert "3 fields" in str(caught.value)
        assert not out_tsv.exists(), "a malformed table must not be written"
        assert not list(out_tsv.parent.glob(".*partial"))

    def test_a_non_numeric_field_is_refused(self, inputs):
        query, assembly, out_tsv = inputs
        bad = ("oprD\tcontig_1\tNOT_A_NUMBER\t300\t300\t5000\t1\t300\t100\t399"
               "\t0.0\t600.0\tMKVM\tMKVM\n")
        with pytest.raises(PipelineError):
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=_fake_runner(table=bad))
        assert not out_tsv.exists()

    def test_a_failed_search_is_raised_with_its_stderr(self, inputs):
        query, assembly, out_tsv = inputs
        with pytest.raises(ToolExecutionError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1,
                           runner=_fake_runner(search_returncode=1,
                                              stderr="BLAST Database error"))
        assert "BLAST Database error" in str(caught.value)
        assert not out_tsv.exists()

    def test_a_failed_makeblastdb_is_raised_before_the_search(self, inputs):
        query, assembly, out_tsv = inputs
        runner = _fake_runner(index_returncode=2)
        with pytest.raises(ToolExecutionError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=runner)
        assert "makeblastdb" in str(caught.value)
        assert len(runner.calls) == 1, "tblastn must not run without a database"

    def test_a_silent_makeblastdb_that_wrote_no_index_is_refused(self, inputs):
        """Exit 0 is not success when the index is absent: tblastn would then
        search nothing, which reads as an isolate with no oprD."""
        query, assembly, out_tsv = inputs
        runner = _fake_runner(write_index=False)
        with pytest.raises(ToolExecutionError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=runner)
        assert ".nin" in str(caught.value)
        assert len(runner.calls) == 1

    def test_a_runner_returning_the_wrong_type_is_refused(self, inputs):
        """A runner that returned something else would be read as a tool
        failure with nothing raised, or worse, as success."""
        query, assembly, out_tsv = inputs

        def runner(command):
            return (0, "", "")

        with pytest.raises(PipelineError) as caught:
            tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                           threads=1, runner=runner)
        assert "CommandResult" in str(caught.value)

    def test_the_injected_runner_receives_the_keyword_protocol(self, inputs):
        """`runner(command=...)`, the house protocol `virulence._run` uses, so a
        caller already faking blastn needs no second fake."""
        query, assembly, out_tsv = inputs
        seen = []

        def runner(command):
            seen.append(command)
            if Path(command[0]).name == "makeblastdb":
                prefix = command[command.index("-out") + 1]
                Path(f"{prefix}.nin").write_text("index", encoding="utf-8")
                return _result(command, stdout="ok")
            return _result(command, stdout=GOOD_ROW + "\n")

        tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(),
                       threads=1, runner=runner)
        assert len(seen) == 2
        assert [Path(c[0]).name for c in seen] == ["makeblastdb", "tblastn"]


HAS_BLAST = bool(shutil.which("tblastn") and shutil.which("makeblastdb"))


def _synthetic_locus(length_codons: int = 100):
    """A stop-free CDS with a terminal TAA, and its protein. Generated, not
    committed, so the real-tool test needs no fixture on disk."""
    import random

    rng = random.Random(20261004)
    stops = {"TAA", "TAG", "TGA"}
    codons = ["ATG"]
    while len(codons) < length_codons - 1:
        codon = "".join(rng.choice("ACGT") for _ in range(3))
        if codon not in stops:
            codons.append(codon)
    codons.append("TAA")
    cds = "".join(codons)
    table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L",
        "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I",
        "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V",
        "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T",
        "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
        "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D",
        "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
        "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G",
        "GGC": "G", "GGA": "G", "GGG": "G",
    }
    protein = "".join(table[cds[i:i + 3]] for i in range(0, len(cds) - 3, 3))
    return cds, protein


#: Where the synthetic locus is laid into its contig, and what the contig is
#: called. Both the query id and the subject id are the production ones.
_LOCUS_CONTIG = "contig_1"
_LOCUS_START = 5000
_CONTIG_LENGTH = 20000


def _exact_row(residues: str, contig: str, subject_start: int) -> str:
    """One ``TBLASTN_OUTFMT`` row aligning the whole query to a perfect match.

    ``length`` and ``qlen`` are the query's own residue count, which is what
    blast writes and what makes this row the sharpest possible probe of the
    invariant: its high query coordinate is the query's length.
    """
    return "\t".join((
        oprd.TBLASTN_QUERY_ID, contig, "100.000", str(len(residues)),
        str(len(residues)), str(_CONTIG_LENGTH), "1", str(len(residues)),
        str(subject_start), str(subject_start + 3 * len(residues) - 1),
        "0.0", f"{2.0 * len(residues):.1f}", residues, residues,
    ))


def _assembly_with(cds: str, contig: str = _LOCUS_CONTIG,
                   start: int = _LOCUS_START) -> dict:
    seq = list("A" * _CONTIG_LENGTH)
    for offset, base in enumerate(cds):
        seq[start - 1 + offset] = base
    return {contig: "".join(seq)}


def _call_over(row: str, assembly: dict, reference_cds: str):
    """One production structural call from one table row and one assembly."""
    hits = oprd.parse_tblastn_table(
        row, query_codon_length=oprd.reference_codon_count(reference_cds))
    return oprd.structural_call(
        "SYNTH", hits, reference_cds=reference_cds, assembly=assembly,
        search=oprd.LocusSearch(
            program="tblastn", query=oprd.TBLASTN_QUERY_ID, database="test_index",
            parameters=(("evalue", "1e-3"), ("max_target_seqs", "5000"),
                        ("db_gencode", "11")),
            command="tblastn -query PA0958_oprD.faa -db index", n_hits=1,
        ))


class TestTheProductionQueryIsTheReferenceItIsMeasuredAgainst:
    """**THE INVARIANT.** The query is exactly as long as the codons it is
    measured against.

    `measure_locus_structure` decides `stop_is_reference_terminator` by asking
    whether the alignment's high query coordinate reached
    `reference_codon_count(reference_cds)` - 444 for PAO1 PA0958, the stop
    codon included. Since `qend <= qlen` always, a query one residue shorter
    can never reach that codon, so the test reads False for **every** isolate
    and an intact oprD is unreachable. Measured on the ten smoke isolates with
    the 443-residue `reference_protein` string as the query: 0 intact / 10
    disrupted, against the round-9 oracle's 2 intact / 8 disrupted.

    So the query is built in ONE place, `oprd_locus.write_tblastn_query`, from
    the same reference CDS the measurement is made against, and its length is
    compared here against the number the measurement actually published
    (`evidence.reference_aa_length`) rather than against a second formula. One
    assertion, and it is the one that failed in production.
    """

    def test_the_query_length_equals_the_codon_count_the_measurer_published(
        self, tmp_path
    ):
        cds, _protein = _synthetic_locus(length_codons=100)
        query = oprd.write_tblastn_query(tmp_path / "w" / "oprD.faa", cds)
        residues = query.read_text(encoding="utf-8").splitlines()[1]

        call = _call_over(_exact_row(residues, _LOCUS_CONTIG, _LOCUS_START),
                          _assembly_with(cds), cds)

        # `reference_aa_length` is what `measure_locus_structure` compares
        # `ref_high` against; the query is what blast reports `qlen` for. If
        # these two ever differ the terminator test is unsatisfiable.
        assert len(residues) == call.evidence.reference_aa_length
        assert len(residues) == oprd.reference_codon_count(cds)

    def test_the_reference_terminator_is_reachable_and_is_reached(self, tmp_path):
        """The invariant with its consequence: the test can now fire.

        `ref_high` is the query's own last position, so it equals
        `reference_codon_count` only because the query is that long. Everything
        downstream - an intact oprD, an in-frame deletion that keeps the native
        stop - was unreachable while it was one short.
        """
        cds, _protein = _synthetic_locus(length_codons=100)
        query = oprd.write_tblastn_query(tmp_path / "oprD.faa", cds)
        residues = query.read_text(encoding="utf-8").splitlines()[1]

        call = _call_over(_exact_row(residues, _LOCUS_CONTIG, _LOCUS_START),
                          _assembly_with(cds), cds)

        assert call.evidence.stop_is_reference_terminator is True
        assert call.evidence.lesion_type == oprd.LesionType.NONE
        assert call.verdict == oprd.StructuralVerdict.INTACT, call.reason

    def test_the_query_carries_the_reference_terminator(self, tmp_path):
        cds, _protein = _synthetic_locus(length_codons=100)
        body = oprd.tblastn_query_protein(cds)
        assert body.endswith("*")
        assert cds[-3:] == "TAA"
        assert len(body) == oprd.reference_codon_count(cds)

    def test_the_writer_uses_the_recorded_query_id_and_creates_its_parents(
        self, tmp_path
    ):
        """blast truncates a FASTA header at the first whitespace.

        A header of `>PA0958 oprD` - which is what the blastp resolver writes -
        would reach the table as `qseqid=PA0958` and rename the query the
        stored artifacts were produced with.
        """
        cds, _protein = _synthetic_locus(length_codons=100)
        query = oprd.write_tblastn_query(tmp_path / "deep" / "q.faa", cds)
        lines = query.read_text(encoding="utf-8").splitlines()
        assert lines[0] == ">" + oprd.TBLASTN_QUERY_ID
        assert " " not in lines[0][1:]
        assert len(lines) == 2

    def test_a_reference_without_a_terminator_is_refused(self):
        """A query missing its stop is not a smaller query; it is a query that
        can never satisfy the terminator test, so it is refused at the writer."""
        cds, _protein = _synthetic_locus(length_codons=100)
        with pytest.raises(PipelineError) as excinfo:
            oprd.tblastn_query_protein(cds[:-3])
        assert "terminator" in str(excinfo.value)

    def test_a_reference_that_is_not_a_whole_number_of_codons_is_refused(self):
        cds, _protein = _synthetic_locus(length_codons=100)
        with pytest.raises(PipelineError) as excinfo:
            oprd.tblastn_query_protein(cds + "A")
        assert "frame" in str(excinfo.value).lower()

    def test_the_codon_count_is_the_one_the_pinned_reference_gives(self):
        """One definition, used by both callers. 444 for PAO1 PA0958."""
        import inspect

        assert inspect.getsource(oprd.reference_codon_count).count("// 3") == 1
        assert "reference_codon_count(reference_cds)" in inspect.getsource(
            oprd.measure_locus_structure)
        assert "reference_codon_count(reference_cds)" in inspect.getsource(
            oprd.structural_call_from_paths)


@pytest.mark.skipif(
    not HAS_BLAST, reason="tblastn/makeblastdb are not installed on this machine"
)
class TestAgainstTheRealTool:
    """The flags are checked by blast itself, not by my memory of it.

    A test asserting the command *contains* `-db_gencode 11` proves only that the
    code contains it. These run the installed blast 2.17 and let it reject the
    flags if they are wrong.
    """

    def _isolate(self, tmp_path, terminator=False):
        """A synthetic isolate carrying a stop-free oprD CDS.

        ``terminator`` appends a literal ``*`` to the query. The production query
        is now written by `oprd_locus.write_tblastn_query`, which is this branch
        by construction; the flag stays so the branch that omits the terminator
        remains exercised end to end.
        """
        cds, protein = _synthetic_locus()
        query = tmp_path / "reference_protein.faa"
        body = protein + "*" if terminator else protein
        query.write_text(f">oprD_test\n{body}\n", encoding="utf-8")
        flank = "".join("ACGT" for _ in range(30))
        assembly = tmp_path / "genome.fna"
        assembly.write_text(
            f">contig_1\n{flank}{cds}{'TTTTTT' * 10}\n", encoding="utf-8"
        )
        ref_cds = tmp_path / "reference_cds.fna"
        ref_cds.write_text(f">oprD_test\n{cds}\n", encoding="utf-8")
        return query, assembly, ref_cds

    def test_blast_reports_the_query_length_the_measurer_counts(self, tmp_path):
        """**The invariant, checked by blast rather than by arithmetic.**

        `qlen` is blast's own count of the query's residues and
        `reference_aa_length` is what `measure_locus_structure` compares the
        alignment's high query coordinate against. Written this way, the two
        come from independent code - one from the tool that searched, one from
        the code that decided - so this is the assertion that would have caught
        the production defect, and it is asserted on the query
        `write_tblastn_query` actually wrote.
        """
        cds, _protein = _synthetic_locus()
        query = oprd.write_tblastn_query(tmp_path / "oprD.faa", cds)
        assembly, ref_cds = self._isolate(tmp_path, terminator=True)[1:]
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=2)
        call = oprd.structural_call_from_paths(
            "SYNTH", reference_cds_path=ref_cds, tblastn_table_path=out_tsv,
            assembly_path=assembly, search=info.search)
        hits = oprd.parse_tblastn_table(
            out_tsv.read_text(encoding="utf-8"),
            query_codon_length=oprd.reference_codon_count(cds))
        assert hits
        assert {h.query_length for h in hits} == \
            {call.evidence.reference_aa_length} == {len(cds) // 3}
        assert {h.query_id for h in hits} == {oprd.TBLASTN_QUERY_ID}
        assert max(h.query_high for h in hits) == \
            call.evidence.reference_aa_length
        assert call.verdict == oprd.StructuralVerdict.INTACT, call.reason

    def test_it_produces_a_table_the_consumer_reads(self, tmp_path):
        query, assembly, ref_cds = self._isolate(tmp_path)
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=2)
        assert info.n_hits >= 1, info.out_tsv
        hits = oprd.parse_tblastn_table(out_tsv.read_text(encoding="utf-8"),
                                        query_codon_length=len(
                                            ref_cds.read_text()) // 3)
        assert hits and hits[0].query_id == "oprD_test"
        assert info.tool_version.startswith("tblastn:")
        assert info.makeblastdb_version.startswith("makeblastdb:")

    def test_the_declared_outfmt_is_what_blast_actually_emits(self, tmp_path):
        """Fourteen columns, every one of them parseable as declared."""
        query, assembly, _ref = self._isolate(tmp_path)
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        tb.run_tblastn(query, assembly, out_tsv, tool_search_dirs=(), threads=2)
        for line in out_tsv.read_text(encoding="utf-8").splitlines():
            if line.strip():
                assert len(line.split("\t")) == len(oprd.TBLASTN_OUTFMT_FIELDS)

    def test_an_intact_locus_reads_as_intact_end_to_end(self, tmp_path):
        """Query, search, table, decision - the whole chain, on real blast.

        The query carries its terminator, because `reference_codons` counts the
        CDS codons including the stop and so a query one residue short can never
        reach the reference's last codon. That length is now the writer's job -
        `oprd_locus.write_tblastn_query` - rather than something each caller has
        to remember; see `TestTheProductionQueryIsTheReferenceItIsMeasuredAgainst`.
        """
        query, assembly, ref_cds = self._isolate(tmp_path, terminator=True)
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=2)
        call = oprd.structural_call_from_paths(
            "SYNTH", reference_cds_path=ref_cds, tblastn_table_path=out_tsv,
            assembly_path=assembly, search=info.search,
        )
        assert call.evidence.stop_is_reference_terminator is True
        assert call.verdict == oprd.StructuralVerdict.INTACT, call.reason

    def test_the_same_isolate_without_the_terminator_in_the_query_is_not_intact(
        self, tmp_path
    ):
        """The round-12 P2 finding, reproduced on synthetic input and real blast.

        Same adapter, same assembly, same reference CDS - one character removed
        from the query - and the verdict moves off `intact`. The production query
        is now built by `oprd_locus.write_tblastn_query`, so a run cannot reach
        this state; what is pinned here is the half that must not move: the
        adapter searches the query it is handed and does not pad it, so a caller
        that writes the wrong query gets a wrong verdict rather than a
        convenient normalisation.
        """
        query, assembly, ref_cds = self._isolate(tmp_path, terminator=False)
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        info = tb.run_tblastn(query, assembly, out_tsv,
                              tool_search_dirs=(), threads=2)
        call = oprd.structural_call_from_paths(
            "SYNTH", reference_cds_path=ref_cds, tblastn_table_path=out_tsv,
            assembly_path=assembly, search=info.search,
        )
        assert call.evidence.stop_is_reference_terminator is False
        assert call.verdict != oprd.StructuralVerdict.INTACT, (
            "If the adapter now pads a short query, the cause of the round-12 "
            "defect is being hidden behind a normalisation instead of fixed."
        )

    def test_a_locus_absent_from_the_assembly_reads_as_absent(self, tmp_path):
        """The negative path, end to end: an unrelated contig, a recorded empty
        search, and the one verdict that may be `absent`."""
        query, _assembly, ref_cds = self._isolate(tmp_path, terminator=True)
        other = tmp_path / "other.fna"
        other.write_text(">contig_9\n" + "ACGTTGCA" * 60 + "\n",
                         encoding="utf-8")
        out_tsv = tmp_path / "work" / "tblastn.tsv"
        info = tb.run_tblastn(query, other, out_tsv, tool_search_dirs=(),
                              threads=2)
        assert info.n_hits == 0, "an unrelated contig must yield no HSP"
        call = oprd.structural_call_from_paths(
            "SYNTH", reference_cds_path=ref_cds, tblastn_table_path=out_tsv,
            assembly_path=other, search=info.search,
        )
        assert call.verdict == oprd.StructuralVerdict.ABSENT, call.reason