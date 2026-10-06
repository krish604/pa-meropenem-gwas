"""W5: the W1 call sites, driven through `run_pipeline` - and their arity proven.

**The gap this closes.** Round 12's W1 wired new real paths into
``papipeline/run.py`` - the oprD structural screen, the GWAS feature producer, the
GWAS engine, locus coverage, the analysed count, the report context. Each was
covered by a test that called the FUNCTION directly. None drove the CALL SITE, so
for every one of them the arity and the argument names were unproven: a call site
passing four of five required arguments, or a keyword the callee does not declare,
raises ``TypeError`` only when the orchestrator reaches that line, which is after
a run has paid for every earlier stage. ``tests/unit/test_seam_stage_call_arity.py``
covers the ``stage_*.run`` entry points; this file covers the OTHER call sites W1
added, which is where the untested ones are.

Two halves, and they answer different questions:

* **PartOne** proves arity and argument NAMES by binding the arguments ``run.py``
  actually writes against the live signatures with ``inspect.Signature.bind``.
  This covers every W1 call site, costs no cohort and no tools, and cannot rot
  into passing for the wrong reason: a mismatch is a ``TypeError`` from ``bind``.
* **PartTwo** proves the same thing at RUNTIME for the call sites that can be
  reached on this machine, by driving the real orchestrator over a synthetic REAL
  cohort with an INJECTED runner (R12).

**What PartTwo injects, and the line it does not cross.** The subject is
``run.py``'s call sites and the code they call: ``write_tblastn_query``,
``adapters.tblastn.run_tblastn``, ``structural_call_from_paths``,
``write_oprd_structural_table``, ``gwas_features.produce_gwas_features``. Those
are production code over real files on real disk, and the orchestrator decides
when they run. What is injected is everything a synthetic cohort cannot supply:
the command runner (required - ``run_pipeline`` takes it as
``oprd_structural_runner`` and ``None`` would execute the real blast binaries),
and the stage functions whose external tools cannot produce a result for a
synthetic genome. Every injection is a named fixture with its reason, and each
test that depends on one says which.

**Not run, and why.** No real sample is analysed. Every isolate, assembly,
Bakta table and phenotype here is synthetic, written by this file.

**The ordering defect this turned up - NOW FIXED, and covered elsewhere.**
``derive_locus_coverage`` was called BEFORE ``stage_variants.run`` in stage 6's
body, so it read BAMs the same run had not yet written and a clean intermediate
root could never get past it. That was the stage 6 failure a REAL 10-isolate run
hit (``n_missing=10``). FIX6-PROBE reordered the body and
``tests/integration/test_stage6_coverage_after_alignment.py`` covers it, by
driving the real measurement off real indexed BAMs rather than injecting it -
which is only possible now that the caller runs first. This file therefore still
INJECTS the measurement: its subject is the eleven call sites below, and on a
synthetic genome the real measurement would be 0.0 everywhere and the honest
threshold would then suppress every regulator record.
"""

from __future__ import annotations

import ast
import copy
import dataclasses
import inspect
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
SCIENCE = REPO / "config" / "science.yaml"
SMOKE = REPO / "config" / "machines" / "smoke.yaml"

#: The W1 call sites PartOne proves. Every name here is reachable in run.py by
#: attribute or bare name; each is asserted to be FOUND, so a name that stops
#: being called fails this file rather than silently shrinking the set.
W1_CALL_SITES: Tuple[str, ...] = (
    "resolve_oprd_structure",
    "write_tblastn_query",
    "run_tblastn",
    "structural_call_from_paths",
    "write_oprd_structural_table",
    "reference_protein",
    "reference_cds_nucleotides",
    "resolve_isolate",
    "produce_gwas_features",
    "PyseerEngine",
    "_build_gwas_engine",
    "load_run_tables",
)

#: Sentinel substituted for every argument. `bind` checks arity and names, never
#: values, so one object per argument is a faithful stand-in - and using a
#: sentinel rather than a real value keeps the test from depending on what run.py
#: happens to pass today.
_SENTINEL = object()


def _resolve(node: ast.AST, module: Any) -> Tuple[Any, str]:
    """The object an `ast` name/attribute chain denotes, plus its last label."""
    if isinstance(node, ast.Name):
        return getattr(module, node.id, None), node.id
    if isinstance(node, ast.Attribute):
        base, _ = _resolve(node.value, module)
        if base is None:
            return None, ""
        return getattr(base, node.attr, None), node.attr
    return None, ""


def _call_sites(module: Any, wanted: Tuple[str, ...]):
    """Every call in run.py whose callee label is in `wanted`, with its line."""
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    found: Dict[str, List[ast.Call]] = {name: [] for name in wanted}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        target, label = _resolve(node.func, module)
        if label in found and target is not None:
            found[label].append(node)
    return found


@pytest.fixture(scope="module")
def w1_call_sites():
    import papipeline.run as run

    return _call_sites(run, W1_CALL_SITES)


class TestPartOneArityAndArgumentNames:
    """What `Signature.bind` proves, and why it is not a substitute for PartTwo."""

    def test_every_named_call_site_is_actually_called(self, w1_call_sites):
        """The non-vacuity guard. A name that is no longer called fails here."""
        never = sorted(name for name, nodes in w1_call_sites.items() if not nodes)
        assert not never, (
            f"run.py no longer calls {never}. Either the wiring was removed - in "
            "which case this file is asserting arity for calls that do not "
            "happen - or the callee was renamed, which is the defect this file "
            "exists to catch."
        )

    def test_each_call_site_binds_against_the_live_signature(self, w1_call_sites):
        import papipeline.run as run

        problems: List[str] = []
        checked = 0
        for name, nodes in sorted(w1_call_sites.items()):
            target = _resolve(
                # The bare-name callee is the function itself; the attribute one
                # is reached through the module that owns it. Both are found by
                # re-walking, so this never guesses.
                next(
                    node.func for node in nodes
                ),
                run,
            )[0]
            for node in nodes:
                checked += 1
                if any(isinstance(a, ast.Starred) for a in node.args) or any(
                    keyword.arg is None for keyword in node.keywords
                ):
                    problems.append(
                        f"{name} at run.py:{node.lineno} hides arguments behind a "
                        "star, so its arity cannot be proven by binding"
                    )
                    continue
                try:
                    inspect.signature(target).bind(
                        *([_SENTINEL] * len(node.args)),
                        **{kw.arg: _SENTINEL for kw in node.keywords if kw.arg},
                    )
                except TypeError as exc:
                    problems.append(
                        f"{name} at run.py:{node.lineno} does not bind against "
                        f"{inspect.signature(target)}: {exc}"
                    )
        assert not problems, (
            f"{len(problems)} of {checked} W1 call site(s) do not bind against "
            "the signature they call:\n  " + "\n  ".join(problems)
        )

    def test_the_binding_actually_checked_something(self, w1_call_sites):
        """A count, so "it passed" cannot mean "it inspected nothing"."""
        total = sum(len(nodes) for nodes in w1_call_sites.values())
        assert total >= len(W1_CALL_SITES), (
            f"only {total} call site(s) were inspected for {len(W1_CALL_SITES)} "
            "named targets; some target resolves to a name that is called more "
            "than once and another resolves to nothing, which the binding test "
            "would not notice on its own"
        )

    def test_a_deliberately_wrong_arity_is_detected(self):
        """The guard on the guard: `bind` must reject a call this test can make.

        A test whose own checking mechanism silently stopped checking would
        report green forever. So one wrong call is bound here and the refusal is
        required - if `bind` ever became a no-op, this fails.
        """
        def callee(required_a, required_b, *, required_kw=None):
            return None

        with pytest.raises(TypeError):
            inspect.signature(callee).bind(_SENTINEL)
        with pytest.raises(TypeError):
            inspect.signature(callee).bind(
                _SENTINEL, _SENTINEL, not_a_declared_keyword=_SENTINEL
            )

    def test_the_oprd_structural_screen_is_reached_with_its_runner(
        self, w1_call_sites
    ):
        """The seam's own parameter, proven to be threaded to the call.

        `run_pipeline` takes `oprd_structural_runner` and hands it to
        `resolve_oprd_structure`, which hands it to `adapters.tblastn.run_tblastn`.
        Every link in that chain is a keyword argument at a call site, and a
        dropped keyword silently runs the real binaries instead.
        """
        node = w1_call_sites["run_tblastn"][0]
        keywords = {kw.arg for kw in node.keywords}
        assert "runner" in keywords, (
            "run.py calls run_tblastn without passing `runner`, so an injected "
            f"runner would not reach the tool call. Keywords: {sorted(keywords)}"
        )
        outer = w1_call_sites["resolve_oprd_structure"][0]
        outer_keywords = {kw.arg for kw in outer.keywords}
        assert {"genomes_dir", "intermediate_root", "runner"} <= outer_keywords, (
            "resolve_oprd_structure is called without the inputs it needs; "
            f"keywords: {sorted(outer_keywords)}"
        )


# ── the synthetic REAL cohort ───────────────────────────────────────────────

#: Two isolates. The cohort size is irrelevant to every claim here; two is the
#: smallest that still exercises a per-sample loop.
ISOLATES: Tuple[str, ...] = ("GCA_000000001.1", "GCA_000000002.1")

#: A stop-free, deterministic stretch of nucleotides. Shape is what the aligner
#: and the parsers care about; nothing here is a biological claim.
_ASSEMBLY_SEQUENCE = "ATGGCTAGC" + "GCTGATCGATCGATCGTAGCTAGCTAGC" * 70

_DATABASE_STRING = "v6.0, light"

_FEATURE_TSV = (
    "# Annotated with Bakta\n# Software: v1.12.1\n"
    f"# Database: {_DATABASE_STRING}\n"
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
    "c1\tcds\t100\t400\t+\tTESTPA_00001\toprD\tOprD family porin\t-\n"
)
_GFF3 = (
    "##gff-version 3\n"
    "c1\tBakta\tCDS\t100\t400\t.\t+\t0\tID=TESTPA_00001;Name=TESTPA_00001\n"
)
_INFERENCE_TSV = (
    "# Software: v1.12.1\n"
    f"# Database: {_DATABASE_STRING}\n"
    "Sequence Id\tLocus Tag\tScore\tEvalue\nc1\tTESTPA_00001\t45.2\t1e-40\n"
)


@pytest.fixture()
def synthetic_real_run(tmp_path):
    """A REAL-mode pipeline whose every input is synthetic and whose roots are tmp.

    Built as a patched copy of the COMMITTED smoke overlay rather than a
    hand-written config, so the mode under test is the one the repository ships -
    including `reuse_tool_output: require`, which is what keeps Bakta unreachable:
    curated Bakta output is STAGED, so stage 2 imports annotation and never
    invokes the tool. That is the whole reason this harness uses the smoke shape,
    and it is asserted rather than assumed.

    Everything is redirected: the genome directory, the PDC table (a synthetic one
    - the real clinical file is never read or written), the subset list, the
    phenotype directory, the Bakta database and every results root.
    """
    from papipeline.config.loader import load_config
    from papipeline.models import RunMode

    genomes = tmp_path / "smoke_genomes"
    genomes.mkdir(parents=True)
    for isolate in ISOLATES:
        (genomes / f"{isolate}.fna").write_text(
            f">c1\n{_ASSEMBLY_SEQUENCE}\n", encoding="utf-8"
        )

    pdc = tmp_path / "PDC_synthetic.tsv"
    pdc.write_text(
        "Isolate\tAssembly\tBioSample\n"
        + "".join(
            f"{isolate}\tGCA_{index:09d}.1\tSAMN{index:08d}\n"
            for index, isolate in enumerate(ISOLATES)
        ),
        encoding="utf-8",
    )
    subset = tmp_path / "isolates.txt"
    subset.write_text("\n".join(ISOLATES) + "\n", encoding="utf-8")

    database = tmp_path / "bakta-db" / "db-light"
    database.mkdir(parents=True)
    (database / "bakta.db").write_text("stub\n", encoding="utf-8")
    (database / "version.json").write_text(
        '{"major": 6, "minor": 0, "type": "light", "date": "2025-02-24"}',
        encoding="utf-8",
    )

    results = tmp_path / "results"
    phenotype = tmp_path / "phenotypes"
    phenotype.mkdir()
    (phenotype / "phenotypes.tsv").write_text(
        "sample_id\tmeropenem\tMIC\n"
        + "".join(f"{isolate}\tR\t8\n" for isolate in ISOLATES),
        encoding="utf-8",
    )

    overlay = yaml.safe_load(SMOKE.read_text(encoding="utf-8"))
    assert overlay["annotation"]["reuse_tool_output"] == "require", (
        "the harness relies on `require` to keep Bakta unreachable; the committed "
        f"smoke overlay says {overlay['annotation']['reuse_tool_output']!r}"
    )
    overlay["paths"].update({
        "smoke_genome_dir": str(genomes),
        "pdc_table": str(pdc),
        "smoke_phenotype_dir": str(phenotype),
        "data_root": str(tmp_path / "data"),
        "results_root": str(results),
        "reports_root": str(results / "reports"),
        "intermediate_root": str(results / "intermediate"),
        "status_dir": str(tmp_path / "status"),
    })
    overlay["cohort"]["subset_file"] = str(subset)
    overlay["runtime"]["allow_real_mode"] = True
    # The observatory writes to a path under $HOME by default; a test must not
    # publish into the user's real store.
    overlay["runtime"]["observatory"] = {"enabled": False}
    overlay_path = tmp_path / "smoke-w5.yaml"
    overlay_path.write_text(yaml.safe_dump(overlay), encoding="utf-8")

    config = load_config(SCIENCE, machine=overlay_path)
    raw = copy.deepcopy(dict(config.raw or {}))
    raw.setdefault("annotation", {})
    raw["annotation"]["bakta_db"] = str(database)
    config = dataclasses.replace(config, raw=raw)

    intermediate = Path(config.intermediate_root(RunMode.REAL))

    for isolate in ISOLATES:
        curated = intermediate / "bakta" / isolate
        curated.mkdir(parents=True, exist_ok=True)
        (curated / f"{isolate}.tsv").write_text(_FEATURE_TSV, encoding="utf-8")
        (curated / f"{isolate}.gff3").write_text(_GFF3, encoding="utf-8")
        (curated / f"{isolate}.inference.tsv").write_text(
            _INFERENCE_TSV, encoding="utf-8"
        )
        # `.faa` is not required by `decide_reuse`, but
        # `run._oprd_locus_inputs` names it as one of the four declared inputs of
        # the structural screen, so a curated set without it would refuse there.
        (curated / f"{isolate}.faa").write_text(
            ">TESTPA_00001_oprD\nMKKIAVTQ\n", encoding="utf-8"
        )

    return {
        "config": config,
        "overlay_path": overlay_path,
        "intermediate": intermediate,
        "genomes": genomes,
        "results": results,
    }


@pytest.fixture()
def injected_locus_coverage(monkeypatch, synthetic_real_run):
    """Coverage 1.0 for every (isolate, locus), injected at the measurement.

    `derive_locus_coverage` reads coverage off BAMs stage 6 wrote, and this
    harness cannot produce one: a synthetic genome aligned against the pinned
    reference yields no alignment, so every locus would measure 0.0. That is not
    a neutral value here - with `min_locus_coverage` now wired at 0.9, a coverage
    of 0.0 correctly suppresses every regulator record, and stage 6 then refuses
    with its own honest "ran and found nothing". **The W1 wiring is what makes
    this harness hard**, which is worth stating plainly rather than working
    around silently.

    So the MEASUREMENT is injected and everything downstream of it is real:
    `stage_regulators.run` receives a non-`None` `locus_coverage` together with
    `min_locus_coverage`, which is the exact combination the wiring exists to
    deliver and the one BACKLOG recorded as unreachable. The suppression
    threshold itself is exercised for real by
    `test_the_min_locus_coverage_threshold_suppresses_a_poorly_covered_locus`.

    Coverage is derived from the config's own regulator list rather than typed
    in, so a locus added to `science.yaml` is covered automatically.
    """
    import papipeline.run as run

    genes = [
        gene for gene, spec in dict(synthetic_real_run["config"].regulators).items()
        if spec.screenable
    ]
    assert genes, "no screenable regulator loci in the config; nothing to cover"

    def derive(config, manifest, *, variants_workdir):
        # The signature is the live one: `derive_locus_coverage` takes no `mode`.
        # `variants_workdir` is the keyword run.py:1222 passes; asserted present so
        # a change to that call site fails here rather than silently dropping the
        # measurement.
        assert isinstance(variants_workdir, Path)
        return {
            (sample_id, gene): 1.0
            for sample_id in manifest.sample_ids
            for gene in genes
        }

    monkeypatch.setattr(run, "derive_locus_coverage", derive)
    return genes


@pytest.fixture()
def fake_blast_runner():
    """An injected command runner that writes what the blast tools would write.

    `run_pipeline`'s `oprd_structural_runner` parameter exists precisely so this
    can exist: `None` runs the real `makeblastdb` and `tblastn`. Two obligations
    are real rather than convenient, and both are properties of the production
    code that a sloppier fake would miss:

    * `run_tblastn` checks that `makeblastdb` actually WROTE `<prefix>.nin` and
      raises `ToolExecutionError` if it did not, so the fake writes it;
    * `run_tblastn` parses the captured stdout with `parse_tblastn_table` and
      refuses a row it would later reject, so the fake returns a table the real
      parser accepts.

    An empty table is a VALID tblastn result meaning zero hits, which is what
    makes the zero-hit path (`absent`/`not_assessed`) reachable without
    fabricating an alignment.
    """
    from papipeline.adapters.external import CommandResult

    calls: List[List[str]] = []

    def runner(command):
        calls.append(list(command))
        program = Path(command[0]).name
        if program == "makeblastdb":
            prefix = command[command.index("-out") + 1]
            Path(f"{prefix}.nin").write_text("synthetic index\n", encoding="utf-8")
        return CommandResult(
            command=list(command), returncode=0, stdout="", stderr="",
            duration_s=0.0,
        )

    runner.calls = calls
    return runner


@pytest.fixture()
def injected_amr(monkeypatch):
    """AMR's declared table with one row, because a synthetic genome has none.

    Stage 4 runs AMRFinderPlus for real. Against a synthetic genome it correctly
    finds nothing, and its output contract then refuses a header-only table - so
    the run stops at stage 4 and never reaches any W1 call site. This writes the
    contract stage 4 itself declares, which is the minimum a later stage needs
    and is honest about what it is: a stand-in for a tool result, not for the
    call site under test.
    """
    import papipeline.stages.amr as stage_amr
    from papipeline.execution.contracts import STAGE_TABLES
    from papipeline.io.tsv import write_tsv

    columns = list(STAGE_TABLES["amr"][1])
    row = {column: "-" for column in columns}
    row["sample_id"] = ISOLATES[0]
    row["gene"] = "blaOXA-1"
    row["source"] = "amrfinderplus"

    def run(config, manifest, mode, tool_output_root, *args, **kwargs):
        root = Path(config.intermediate_root(mode))
        write_tsv(root / "stages" / "04_amr.tsv", [row], columns)
        return {}

    monkeypatch.setattr(stage_amr, "run", run)
    return row

# ── PartTwo: the call sites, driven through the orchestrator ────────────────

class TestPartTwoTheOprDStructuralScreenIsDrivenByRunPipeline:
    """`only=["variants"]` reaches `resolve_oprd_structure` and every call below it.

    This is the W1 call site with its own injection parameter, so it is the one
    that can be driven end to end on this machine. Everything asserted here is
    production code having run: the orchestrator chose the stage, wrote the
    query, handed it to the adapter, read the verdict back off disk, and wrote
    the table it declares.
    """

    @staticmethod
    def _drive(env, runner):
        """Run the orchestrator, and return ``(result_or_None, stage_error_or_None)``.

        The drive COMPLETES. A synthetic genome has no variants against the
        pinned reference, so `injected_variants` supplies one inside a screened
        regulator locus and `injected_locus_coverage` supplies full coverage -
        without both, stage 6 writes a header-only `variants.tsv` and its own
        output contract refuses it, one line before the call sites under test.
        With them, stage 6 completes and the whole variants branch runs for real.

        The `(result, error)` shape exists for one reason:
        `test_the_drive_completes_and_reaches_every_call_site` is the arity
        regression guard, and it needs the exception rather than a bare raise in
        order to say what KIND of failure it was. A `TypeError` here is the
        defect this file was written to find, and one was found and fixed.
        """
        from papipeline.errors import StageError
        from papipeline.run import run_pipeline

        try:
            result = run_pipeline(
                config_path=SCIENCE,
                config=env["config"],
                mode="REAL",
                only=["variants"],
                write_html=False,
                machine=env["overlay_path"],
                oprd_structural_runner=runner,
            )
        except StageError as exc:
            return None, exc
        return result, None

    def test_the_drive_completes_and_reaches_every_call_site(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """The arity regression guard, and the reason the rest of this class works.

        A `TypeError` from any call site is the defect this file exists to find,
        and one was found and fixed: `derive_regulator_table` was called with a
        `min_locus_coverage` keyword it did not declare, which raised on every
        REAL run of stage 6. A drive that raised a `TypeError` would produce all
        the artefacts below it never reached, so those assertions would not fire;
        requiring the failure to be the contract refusal keeps them meaningful.
        """
        result, error = self._drive(synthetic_real_run, fake_blast_runner)
        assert error is None, (
            f"the drive did not complete: {error}\n\nIf this is a TypeError, a "
            "run.py call site is missing or misspelling an argument - the exact "
            "defect this file exists to catch, and one that WAS found here: "
            "`derive_regulator_table` was called with a `min_locus_coverage` "
            "keyword it did not declare, raising on every REAL run of stage 6. "
            "If it is a DataContractError, a stage wrote a table its own "
            "contract refused, which means an injected table is no longer "
            "shaped the way the stage expects."
        )
        assert result is not None
        assert result.outputs.get("oprd_structural_calls") is not None, (
            "the run completed but recorded no `oprd_structural_calls`, so "
            "`write_oprd_structural_table` was not reached through the "
            f"orchestrator. Stage status: {result.stage_status}"
        )
        assert result.stage_status.get("variants") == "completed", (
            "stage 6 did not complete, so the assertions below are about a run "
            f"that stopped early: {result.stage_status}"
        )

    def test_write_tblastn_query_wrote_a_query_carrying_the_reference_stop(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """`write_tblastn_query`, reached through the orchestrator.

        The property asserted is the one that was broken in production: the query
        carries the reference's terminal stop, so its residue count equals
        `reference_codon_count` rather than the protein's 443. A query without
        the terminator is what made `stop_is_reference_terminator`
        unsatisfiable and every isolate read `disrupted`.
        """
        from papipeline.adapters.oprd_locus import reference_codon_count

        self._drive(synthetic_real_run, fake_blast_runner)

        # Precisely the structural screen's query, at the path run.py:2746
        # builds: `oprd_locus/<sample_id>/tblastn/reference_protein.faa`.
        #
        # NOT a bare rglob, because a second producer writes a file of the SAME
        # NAME elsewhere in the same tree - `resolve_isolate` writes the blastp
        # query at `oprd_locus/<sample_id>/reference_protein.faa`, and that one
        # carries 443 residues with no terminator. A glob wide enough to catch
        # both finds the wrong one first and reports the terminator defect that
        # OPRD-FIX fixed. The two files coexisting is a real finding in its own
        # right (see BACKLOG round 12) but it is not what this test is about.
        work_root = synthetic_real_run["intermediate"] / "oprd_locus"
        queries = sorted(work_root.glob("*/tblastn/reference_protein.faa"))
        assert len(queries) == len(ISOLATES), (
            f"expected one structural query per isolate at "
            f"`oprd_locus/*/tblastn/reference_protein.faa`, found {len(queries)}: "
            f"{[str(q) for q in queries]}. `write_tblastn_query` was not reached "
            "through run_pipeline."
        )
        expected = reference_codon_count(_reference_cds(synthetic_real_run["config"]))
        for query in queries:
            residues = sum(
                len(line.strip())
                for line in query.read_text(encoding="utf-8").splitlines()
                if line.strip() and not line.startswith(">")
            )
            assert residues == expected, (
                f"{query.name} carries {residues} residues but the reference has "
                f"{expected} codons; a query without the terminal stop is what "
                "made the terminator invariant unsatisfiable in production"
            )

    def test_run_tblastn_reached_the_injected_runner_once_per_tool_per_isolate(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """`run_tblastn`, reached through the runner the orchestrator was given.

        This is also arity evidence at runtime: the runner is called with the
        house protocol `runner(command=...)`, once per tool per isolate, and the
        thread count the RUN was configured with is the one passed to tblastn -
        the last of which is a wiring claim no direct call to the adapter could
        make.
        """
        self._drive(synthetic_real_run, fake_blast_runner)

        programs = [Path(call[0]).name for call in fake_blast_runner.calls]
        assert programs.count("makeblastdb") == len(ISOLATES), (
            f"makeblastdb was called {programs.count('makeblastdb')} times for "
            f"{len(ISOLATES)} isolates: {programs}"
        )
        assert programs.count("tblastn") == len(ISOLATES), (
            f"tblastn was called {programs.count('tblastn')} times for "
            f"{len(ISOLATES)} isolates: {programs}"
        )
        configured = int(
            synthetic_real_run["config"].runtime.get("threads", 1) or 1
        )
        for call in fake_blast_runner.calls:
            if Path(call[0]).name != "tblastn":
                continue
            assert any(
                arg.endswith("reference_protein.faa") for arg in call
            ), f"tblastn was not given the query run.py wrote: {call}"
            assert call[call.index("-num_threads") + 1] == str(configured), (
                f"tblastn was given -num_threads "
                f"{call[call.index('-num_threads') + 1]}, which is not the count "
                f"the run was configured with ({configured})"
            )

    def test_the_injected_runner_was_used_so_no_real_blast_binary_ran(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """The guard on the guard: `None` would run the real tools.

        `run_pipeline(oprd_structural_runner=None)` runs `makeblastdb` and
        `tblastn` for real. This asserts the parameter was honoured rather than
        accepted-and-ignored, which is the failure mode an injected runner has -
        one that would make every other test in this class pass while executing
        real binaries.
        """
        self._drive(synthetic_real_run, fake_blast_runner)
        assert fake_blast_runner.calls, (
            "the injected runner was never called, so the parameter did not reach "
            "the tool call and the real blast binaries would have run instead"
        )

    def test_structural_call_from_paths_was_given_the_paths_it_declares(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """The reference CDS and the search table exist where the callee reads them.

        `structural_call_from_paths` is handed PATHS, not parsed rows, so the
        only way it can have been reached correctly is if the caller wrote the
        reference CDS and the adapter wrote the table first. Both are asserted by
        their presence on disk, which is what the callee opens.
        """
        from papipeline.run import OPRD_REFERENCE_CDS_NAME

        work_root = synthetic_real_run["intermediate"] / "oprd_locus"
        self._drive(synthetic_real_run, fake_blast_runner)

        reference = work_root / OPRD_REFERENCE_CDS_NAME
        assert reference.is_file(), (
            f"{OPRD_REFERENCE_CDS_NAME} was not written; "
            "structural_call_from_paths is given its path and cannot proceed"
        )
        tables = sorted(work_root.rglob("tblastn.tsv"))
        assert len(tables) == len(ISOLATES), (
            f"expected one tblastn table per isolate, found {len(tables)}: "
            f"{[str(t) for t in tables]}"
        )
        # A zero-hit search is a valid, non-empty ANSWER, and it is the only
        # shape that may later be read as `absent` - so the table must exist and
        # must be empty, and a table the adapter did not write would be neither.
        for table in tables:
            assert table.stat().st_size == 0, (
                f"{table} is non-empty; the injected runner returns a zero-hit "
                "table, so anything here means the table was written by something "
                "other than the adapter this test drives"
            )

    def test_write_oprd_structural_table_wrote_a_verdict_row_per_isolate(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage
    ):
        """`write_oprd_structural_table`, reached through the orchestrator.

        Asserted on the FILE rather than on `result.outputs`, because the run
        refuses afterwards (see `_drive`) and never returns a result object to
        inspect. The table is what the call site produced either way, and the
        path is the one the orchestrator builds.
        """
        from papipeline.run import OPRD_STRUCTURAL_TABLE_NAME

        self._drive(synthetic_real_run, fake_blast_runner)

        table = (
            synthetic_real_run["intermediate"] / "oprd_locus"
            / OPRD_STRUCTURAL_TABLE_NAME
        )
        assert table.is_file(), (
            f"{table} was not written, so `write_oprd_structural_table` was not "
            "reached through run_pipeline"
        )
        lines = table.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        # `structural_verdict`, not `verdict`: the column name is taken from
        # OPRD_STRUCTURAL_COLUMNS rather than typed in, so a rename of the
        # contract fails here instead of passing on a name that no longer exists.
        from papipeline.run import OPRD_STRUCTURAL_COLUMNS

        assert header == list(OPRD_STRUCTURAL_COLUMNS), (
            "the structural table's header does not match the columns "
            f"OPRD_STRUCTURAL_COLUMNS declares: {header}"
        )
        verdict_column = header.index("structural_verdict")
        rows = [line.split("\t") for line in lines[1:] if line.strip()]
        assert sorted(row[0] for row in rows) == sorted(ISOLATES), (
            f"the structural table carries {[r[0] for r in rows]} for a cohort "
            f"of {sorted(ISOLATES)}; one row per isolate is the whole contract"
        )
        # And each row carries a decision, not a blank. The injected runner
        # returns a zero-hit search, so the verdict is the not-assessed/absent
        # family - what matters is that it is a NAMED value rather than an
        # empty cell, because an empty cell is what "measured nothing" looks
        # like and is indistinguishable from an unwritten row.
        verdicts = [row[verdict_column] for row in rows]
        assert all(verdicts), (
            f"every structural verdict is empty: {verdicts}. An empty verdict "
            "column means the screen ran and decided nothing."
        )

    def test_bakta_was_never_invoked_on_the_real_path_this_drives(
        self, synthetic_real_run, fake_blast_runner, injected_amr,
        injected_variants, injected_locus_coverage, monkeypatch
    ):
        """`require` holds on the REAL path this harness drives.

        The curated Bakta output is staged, so stage 2 imports annotation. If a
        change made it fall through to the tool this raises - which is the whole
        reason the harness is built on the smoke overlay's `require` rather than
        on any other configuration.
        """
        import papipeline.adapters.bakta as bakta_adapter

        def forbidden(*args, **kwargs):
            raise AssertionError(
                "Bakta was invoked under reuse_tool_output=require; curated "
                "output is staged, so the tool must not be reached"
            )

        monkeypatch.setattr(bakta_adapter, "run_bakta", forbidden)
        self._drive(synthetic_real_run, fake_blast_runner)
        # No assertion: `forbidden` raising IS the failure.


def _reference_cds(config) -> str:
    """The pinned reference oprD CDS, through production code."""
    from papipeline.adapters.oprd_locus import reference_cds_nucleotides
    from papipeline.run import OPRD_LOCUS_TAG

    return reference_cds_nucleotides(
        config.reference_gff(), config.reference_fasta(), OPRD_LOCUS_TAG
    )


@pytest.fixture()
def injected_variants(monkeypatch, synthetic_real_run):
    """One `variants.tsv` row that lands INSIDE a screened regulator locus.

    Stage 6 calls the variants caller against the PINNED reference, and a
    synthetic sequence unrelated to it yields zero calls - so the stage writes a
    header-only table, which stops the run twice: `derive_regulator_table` reads
    it at `papipeline/run.py:2156` and refuses an empty one, and stage 6's own
    output contract refuses it too. Both refusals sit BEFORE
    `resolve_oprd_structure` at :1301, so without this the oprD screen is
    unreachable on synthetic input at all.

    The row's position is not typed in: it is read out of the pinned reference's
    own GFF through the production interval loader, and the reference allele is
    read out of the pinned FASTA with `pysam`. A synthetic variant placed
    anywhere else would simply fall outside the screened loci, and the regulator
    screen would then refuse with its own honest "ran and found nothing" message
    - which is correct behaviour and would stop the run one line later.

    So this replaces the TABLE the caller produced, not the reading of it: the
    screen, `derive_regulator_table` and `resolve_oprd_structure` all still run
    for real afterwards.
    """
    import pysam

    import papipeline.stages.variants as stage_variants
    from papipeline.adapters import gff as gff_adapter
    from papipeline.io.tsv import write_tsv

    config = synthetic_real_run["config"]
    specs = {
        gene: spec for gene, spec in dict(config.regulators).items()
        if spec.screenable
    }
    intervals = gff_adapter.load_gene_intervals(
        config.reference_gff(), [spec.locus_tag for spec in specs.values()]
    )
    located = [
        (gene, intervals[spec.locus_tag])
        for gene, spec in sorted(specs.items())
        if spec.locus_tag in intervals
    ]
    assert located, (
        "no regulator locus was found in the pinned reference GFF, so this "
        "harness cannot place a variant inside a screened locus"
    )
    gene, interval = located[0]
    position = interval.start + 1  # 1-based, strictly inside the CDS
    with pysam.FastaFile(str(config.reference_fasta())) as fasta:
        reference_base = fasta.fetch(interval.contig, position - 1, position).upper()

    # Column names come from the stage's own PER_ISOLATE_COLUMNS rather than
    # being typed in: they are `chrom`/`pos`, not `contig`/`position`, and a
    # row carrying the wrong key would be silently all "-" and land outside every
    # screened locus - which is a silent wrong answer, not a failure.
    columns = list(stage_variants.PER_ISOLATE_COLUMNS)
    row = {column: "-" for column in columns}
    row.update({
        "sample_id": ISOLATES[0],
        "chrom": interval.contig,
        "pos": str(position),
        "ref": reference_base,
        "alt": "ACGT"[("ACGT".index(reference_base) + 1) % 4],
    })
    assert set(row) == set(columns), (row, columns)

    def run(config, manifest, mode, *, data_root, workdir, statuses):
        # `data_root`/`workdir`/`statuses` are the keywords run.py actually
        # passes (see run.py:1222); naming them keeps the double's signature
        # honest about what the call site supplies, so a change there fails here
        # too. `statuses` is an out-param: the real stage fills it and stage 6a
        # reads it back, so the double fills it as well rather than leaving the
        # second half of that handshake untested.
        assert isinstance(data_root, Path) and isinstance(workdir, Path)
        assert isinstance(statuses, dict)
        # The RETURN value is what matters: run.py:1234 rewrites `variants.tsv`
        # from it via `write_tsv`, so a double that only wrote the file would be
        # overwritten with a header-only table one line later.
        for sample_id in ISOLATES:
            statuses[sample_id] = stage_variants.STATUS_CALLED
        return {ISOLATES[0]: [row], ISOLATES[1]: []}

    monkeypatch.setattr(stage_variants, "run", run)
    return row
