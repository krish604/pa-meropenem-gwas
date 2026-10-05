"""Every REAL-mode input has a producer inside the run, or a refusal that names it.

**The seam this closes.** A REAL run's inputs come from three places that the
code does not distinguish: files this run produces, files an external tool
produces, and files the operator provisions. Nothing in the pipeline's type
system or control flow separates them, so a stage that reads a fourth kind - a
file nothing produces and nobody explains - compiles, runs, and produces output
over whatever it found. The failure is silent in the worst way: the stage
returns a *result*, so the run reports a finding that was computed from nothing.

That is not hypothetical here. `stages/regulators.py` read
`intermediate/regulators/regulator_variants.tsv`, which only the TEST fixture
generator ever wrote; a REAL run over real isolates found no table, treated the
absence as an empty screen, and reported zero regulators screened. The
`stages/pangenome.py` seam was the same shape with a different body. Both are
recorded in `workflow/Snakefile`'s comments as history.

**What "closure" means here, precisely.** Each REAL-mode input is classified
into exactly one of three verdicts:

* ``PRODUCED`` - the path lies inside this run's own output tree, so a rule in
  `workflow/Snakefile` writes it during the run. Nothing external is needed.
* ``REFUSED`` - the input is not produced by the run, and the module that reads
  it **raises naming that input** when it is absent. The operator is told what
  to provide.
* ``ORPHAN`` - neither. Read by the pipeline, produced by nothing, named by no
  error message. **This is the defect, and the test fails on it.**

The third verdict is the whole point. An input whose absence produces no error
is the defect; the test exists so that adding one cannot pass review.

**Why the roots are made distinct, and why that is not decoration.** Every
verdict above is a statement about an *absolute path*, and an absolute path is
only meaningful if the roles cannot silently collapse into one another. So this
test builds a configuration in which every root is a separate directory under
`tmp_path`, and then asserts the *exact* set of coincidences the design intends:

* ``tool_output_root == intermediate_root`` in REAL - asserted, because that
  coincidence is load-bearing and a refactor that separated them would break
  every REAL tool-output read;
* ``tool_output_root != intermediate_root`` in TEST - asserted, because that
  inequality is the only thing stopping a TEST run reading its own outputs as
  inputs, and it is invisible in REAL where they are one path;
* ``data_root == assembly_root`` in REAL - asserted, and asserted *because* it
  is true only when no smoke overlay is configured. A blanket "all roots are
  distinct" assertion would be false here, and asserting it anyway is how a
  test gets ignored.

The sharpest consequence is asserted last: **no two inputs in the inventory may
resolve to the same absolute path.** If they did, one input's producer could
satisfy another's absence, and every verdict above would be a false positive.

**What the test does NOT do.** It does not run a REAL pipeline. REAL is gated
(`runtime.allow_real_mode=false` in both overlays, AGENTS.md rule 6) and this
test resolves REAL *paths* only - no cohort is read, no genome is opened, no
stage executes. `tests/unit/test_seam_stage_call_arity.py` covers the other half
of the seam: that every stage call site supplies every required parameter.

**Anti-vacuity.** Four guards, because a closure test that finds nothing still
passes:

1. every one of the eight roots resolves under `tmp_path`, and none of the
   eight is the repository;
2. the inventory is non-empty, covers every role, and every entry's owner module
   still exists and still contains the literal that names the input - so a
   rename fails here rather than silently invalidating the record;
3. a **live probe**: a synthetic filename that nothing in the pipeline names is
   classified ``ORPHAN`` by the same code that classifies the real entries, and
   a synthetic *producer* flips it to ``PRODUCED``. A classifier that answered
   ``REFUSED`` to everything would pass the real entries and fail this;
4. the derived direction - every filename literal that appears inside a `raise`
   anywhere in `papipeline/` must be an inventoried input or a recorded
   producer. A new refusal that names a file nobody produces and nobody
   inventoried fails the test, so the inventory cannot fall behind the code.
"""

from __future__ import annotations

import ast
import importlib
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Mapping, Tuple

import pytest

from papipeline.config.loader import RESULTS_ROOT_ENV, load_config
from papipeline.models import RunMode

REPO = Path(__file__).resolve().parents[2]

#: The eight roots, named for the accessors that produce them. The order is the
#: declaration order and is not significant; the mapping is.
ROLE_ACCESSORS: Mapping[str, str] = {
    "data_root": "data_root",
    "assembly_root": "assembly_root",
    "results_root": "results_root",
    "intermediate": "intermediate_root",
    "phylogeny_dir": "phylogeny_dir",
    "reports_root": "reports_root",
    "db_root": "db_root",
    "reference": "reference_fasta",
    # Not an accessor, and deliberately so. The knowledge tables and the
    # science configuration are read from `config/` and are as much a REAL-mode
    # input as any genome: `loader.py:916` refuses by naming
    # `config/mechanisms.tsv` when a gene is absent from it. Leaving the config
    # directory out of the role list would have made the derived check in
    # `TestTheInventoryFollowsTheCode` fail on four real inputs rather than on
    # a defect - which is how a check like that gets deleted.
    "config_root": "root/config",
}

#: Role pairs that must be the *same* directory in REAL, with the reason. Every
#: other pair must differ; `test_no_two_roles_collide_in_real` says so.
INTENDED_REAL_COINCIDENCES: Mapping[Tuple[str, str], str] = {
    ("intermediate", "tool_output_root"): (
        "REAL runs the external tools into this run's own intermediate root and "
        "reads the results back from it, so the stage input root and the run "
        "output root are one path. This is what makes a swap between them "
        "invisible in REAL and only visible in TEST."
    ),
    ("data_root", "assembly_root"): (
        "True only when no smoke overlay is configured: `assembly_root` defers "
        "to `data_root` for an ordinary run and only diverges for a bounded "
        "smoke run. Asserted so that a change which silently narrowed every "
        "ordinary run to the smoke directory cannot pass."
    ),
}

#: Nesting the design intends, so a refactor that flattens a root is caught.
INTENDED_REAL_NESTING: Mapping[Tuple[str, str], str] = {
    ("intermediate", "phylogeny_dir"): (
        "REAL writes stage 9's tree under the run's own intermediate root, "
        "which is what keeps it out of the read-only data/ tree."
    ),
    ("results_root", "reports_root"): (
        "Reports follow the results redirect, so a redirected run writes its "
        "report beside its results rather than into the repository."
    ),
}

#: Extensions that make a string literal a filename rather than prose. Without
#: this the derived check in `test_every_refused_filename_is_inventoried` would
#: match every sentence containing a dot.
FILENAME_EXTENSIONS: Tuple[str, ...] = (
    ".tsv", ".fasta", ".fa", ".fna", ".fai", ".gff", ".nwk", ".aln", ".json",
    ".txt", ".vcf", ".bcf", ".gz", ".yaml", ".yml", ".log",
)
FILENAME_RE = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9._+-]*("
    + "|".join(re.escape(e) for e in FILENAME_EXTENSIONS)
    + r")$"
)
#: The same shape, but *searching* rather than matching the whole string.
#:
#: A refusal almost never consists of a bare filename. `raise ValueError("missing
#: kinship_matrix.tsv")` and `f"Add a row to config/mechanisms.tsv before
#: rerunning"` both name a file from inside a sentence, and requiring the whole
#: literal to be the filename finds neither - which is how the first version of
#: the derived check passed against a probe that should have failed it.
FILENAME_SEARCH_RE = re.compile(
    r"\b[A-Za-z0-9][A-Za-z0-9._+-]*("
    + "|".join(re.escape(e) for e in FILENAME_EXTENSIONS)
    + r")\b"
)


def _filenames_in(value: str) -> FrozenSet[str]:
    """Every filename mentioned anywhere inside `value`."""
    # `finditer`, not `findall`: with one group `findall` yields the group -
    # the bare extension - rather than the whole filename.
    return frozenset(
        match.group(0) for match in FILENAME_SEARCH_RE.finditer(value)
    )

#: Write-shaped calls. A filename appearing as an argument to one of these is
#: being produced somewhere in the pipeline.
WRITE_CALLS: FrozenSet[str] = frozenset(
    {
        "write_tsv", "write_text", "write_bytes", "mkdir", "touch", "copy",
        "copy2", "copyfile", "copytree", "replace", "rename", "symlink_to",
        "hardlink_to",
    }
)


@dataclass(frozen=True)
class Input:
    """One REAL-mode input, and who is answerable for it.

    Attributes:
        key: Stable name, used in failure messages.
        role: A key of `ROLE_ACCESSORS`.
        tail: Path below the role root. May contain `{antibiotic}`, which is
            filled from the loaded configuration.
        verdict: What the code has to show for this input. ``ORPHAN`` is the
            failure and is deliberately not a value this file records:
            recording one would be the silent quarantine the test exists to
            prevent.
        owner: Dotted module that *reads* the input.
        producer: Dotted module that *writes* it, or ``""`` when nothing in the
            pipeline does.
    producer_token: What to look for in ``producer``, when it differs. Most
            tables are written by the same module that reads them and share a
            token; `variants.tsv` does not - `stages/regulators.py` refuses it
            by name while `run.py` writes it through
            `execution.contracts.table_path`, which is where the filename
            actually lives.
    producer_evidence: Evidence kind for ``producer_token``. Separate from ``owner`` because the two genuinely
            differ: `papipeline.run` writes the masked SNP alignment and also
            refuses to run without panaroo's core alignment, and collapsing the
            two into one field would let the writer vouch for the reader.
        token: What to look for in ``owner``'s source.
    refusal: The literal the *refusal* names, when that is not ``token``.
        Needed for the one entry that is both produced and refused: the writer
        is `_write_screen_report` and the refusal names
        `regulator_variants.tsv`, and asking one token to be both an identifier
        and a string in a ``raise`` would force one of the two claims to be
        dropped.
    evidence: How to look. ``"raise"`` - the token must appear inside a
            ``raise``, which is what "a refusal that names it" means.
            ``"constant"`` - the token must be (or be part of) a module-level
            string constant, which is how a table's name is declared and the
            strictest of the four.
            ``"literal"`` - the token must appear in some string literal, which
            covers a filename written inline in a path expression.
            ``"identifier"`` - the token must be a name the module binds or
            calls, which is the only handle available when the filename is
            computed (``stub_report_name`` builds the report's name) or when
            the module is the writer rather than the reader.

    The kind matters: a single "does the source mention it" test is satisfied
    by a docstring, and this pipeline has several modules that *discuss* a
    filename in prose about a past defect. Requiring the specific construct is
    what makes a renamed input fail here.
    """

    key: str
    role: str
    tail: str
    verdict: str
    owner: str
    token: str
    evidence: str = "raise"
    producer: str = ""
    refusal: str = ""
    producer_token: str = ""
    producer_evidence: str = ""


#: The inventory. Read off the code, not off a wish: each entry's owner and
#: token are asserted against the live source by
#: `test_every_inventoried_input_is_still_named_by_its_owner`.
INVENTORY: Tuple[Input, ...] = (
    Input(
        "cohort_manifest", "data_root", "metadata/sample_metadata.tsv",
        "REFUSED", "papipeline.manifest", "sample_metadata.tsv", "constant",
    ),
    Input(
        "phenotype_table", "data_root", "phenotype/{antibiotic}_phenotype.tsv",
        "REFUSED", "papipeline.stages.phenotype",
        "data/phenotype/<antibiotic>_phenotype.tsv", "raise",
    ),
    Input(
        "isolate_assemblies", "assembly_root", "GCA_000710625.1_GCA_001364615.1",
        "REFUSED", "papipeline.assemblies", "No assembly found",
    ),
    Input(
        "reference_genome", "reference", "GCF_000006765.1_genomic.fna",
        "REFUSED", "papipeline.reference", "is not present at",
    ),
    Input(
        "amr_database", "db_root", "amrfinderplus",
        "REFUSED", "papipeline.adapters.amr", "no indexed protein FASTA",
        "raise",
    ),
    Input(
        "core_gene_alignment", "intermediate",
        "panaroo/core_gene_alignment_filtered.aln",
        "REFUSED", "papipeline.run", "core gene alignment", "raise",
    ),
    Input(
        "masked_snp_alignment", "phylogeny_dir", "core_snp_alignment.fasta",
        "PRODUCED", "papipeline.run", "build_outputs", "identifier",
        producer="papipeline.run",
    ),
    Input(
        "amr_table", "intermediate", "amr/amr_determinants.tsv",
        "PRODUCED", "papipeline.stages.amr", "_write_determinant_table",
        "identifier", producer="papipeline.stages.amr",
    ),
    Input(
        "mlst_table", "intermediate", "mlst/mlst_results.tsv",
        "PRODUCED", "papipeline.stages.mlst", "_write_mlst_results",
        "identifier", producer="papipeline.stages.mlst",
    ),
    Input(
        "structural_variant_table", "intermediate",
        "structural_variants/structural_variants.tsv",
        "PRODUCED", "papipeline.stages.sv", "_write_sv_tsv", "identifier",
        producer="papipeline.stages.sv",
    ),
    Input(
        "virulence_table", "intermediate", "virulence/virulence_factors.tsv",
        "PRODUCED", "papipeline.stages.virulence", "_write_virulence_tsv",
        "identifier", producer="papipeline.stages.virulence",
    ),
    Input(
        "regulator_table", "intermediate", "regulators/regulator_variants.tsv",
        "PRODUCED+REFUSED", "papipeline.stages.regulators",
        "_write_screen_report", "identifier",
        producer="papipeline.stages.regulators",
        refusal="regulator_variants.tsv",
    ),
    Input(
        # `owner` is the reader, `producer` the writer. They differ here: the
        # reader is `stages/gwas.py`, whose only occurrence of the name is the
        # read at gwas.py:1743, and the writer is `stages/gwas_features.py`.
        # Naming `stages.gwas` as the producer was a false statement that the
        # classifier could not catch, because it only asks whether a module
        # writes *something* and whether it *mentions* the name.
        "gwas_feature_table", "intermediate", "gwas/gwas_features.tsv",
        "PRODUCED", "papipeline.stages.gwas", "gwas_features.tsv", "constant",
        # The producer is named by the FUNCTION that writes it, not by the
        # filename appearing somewhere in its module. The old evidence was
        # `literal`, which asks only whether `gwas_features.tsv` occurs in
        # `gwas_features.py` - and the seam-matrix's G7 finding is that the same
        # question is true of the READER, because `gwas.py:1743` names the table
        # it loads. Naming the writer's function cannot be satisfied by a module
        # that merely mentions the file, so a producer that stopped writing it
        # fails here rather than passing on a docstring.
        producer="papipeline.stages.gwas_features",
        producer_token="produce_gwas_features",
        producer_evidence="identifier",
    ),
    Input(
        "pangenome_gpa", "intermediate", "pangenome/gene_presence_absence.tsv",
        "PRODUCED", "papipeline.stages.pangenome", "gene_presence_absence.tsv",
        "literal", producer="papipeline.stages.pangenome",
    ),
    Input(
        "phylogeny_tree", "phylogeny_dir", "tree.nwk",
        "PRODUCED", "papipeline.stages.phylogeny", "tree.nwk", "literal",
        producer="papipeline.stages.phylogeny",
    ),
    Input(
        "phylogeny_metadata", "phylogeny_dir", "tree_metadata.tsv",
        "PRODUCED", "papipeline.stages.phylogeny", "tree_metadata.tsv",
        "literal", producer="papipeline.stages.phylogeny",
    ),
    Input(
        "run_manifest", "results_root", "run_manifest.json",
        "PRODUCED", "papipeline.run", "run_manifest.json", "literal",
        producer="papipeline.run",
    ),
    Input(
        "pdc_cohort_table", "data_root", "PDC_essential.tsv",
        "REFUSED", "papipeline.pdc", "PDC_essential.tsv not found", "raise",
    ),
    Input(
        "machine_overlay", "config_root", "machines/bigmachine.yaml",
        "REFUSED", "papipeline.config.loader", "config/machines/bigmachine.yaml",
        "raise",
    ),
    Input(
        "regulator_locus_table", "config_root", "regulators.tsv",
        "REFUSED", "papipeline.stages.regulators", "config/regulators.tsv",
        "raise",
    ),
    Input(
        "stage_variants_table", "intermediate", "stages/variants.tsv",
        "PRODUCED+REFUSED", "papipeline.stages.regulators", "stages/variants.tsv",
        "raise", producer="papipeline.run", refusal="stages/variants.tsv",
        producer_token="variants", producer_evidence="literal",
    ),
    Input(
        "pangenome_summary_table", "intermediate",
        "pangenome/pangenome_summary.tsv",
        "PRODUCED+REFUSED", "papipeline.stages.pangenome",
        "pangenome_summary.tsv", "literal",
        producer="papipeline.stages.pangenome",
    ),
    Input(
        "science_config", "config_root", "science.yaml",
        "REFUSED", "papipeline.config.loader", "science.yaml",
    ),
    Input(
        "antibiotics_table", "config_root", "antibiotics.tsv",
        "REFUSED", "papipeline.config.loader", "antibiotics.tsv",
    ),
    Input(
        "mechanisms_table", "config_root", "mechanisms.tsv",
        "REFUSED", "papipeline.config.loader", "config/mechanisms.tsv",
    ),
    Input(
        "references_table", "config_root", "references.tsv",
        "REFUSED", "papipeline.adapters.amr",
        "is missing, so stage 4 has no pinned reference",
    ),
    Input(
        "run_report", "reports_root",
        "Pseudomonas_aeruginosa_Imipenem_AMR_Report.md",
        "PRODUCED", "papipeline.stages.reporting", "write_report",
        "identifier", producer="papipeline.stages.reporting",
    ),
)


def _realpath(path: Path) -> Path:
    """`/var` and `/private/var` are one directory on macOS. One path, two spellings."""
    return Path(os.path.realpath(path))


@pytest.fixture(scope="module")
def sandbox(tmp_path_factory, monkeypatch_module) -> Tuple[Path, object]:
    """A pipeline whose every root is a distinct directory under `tmp_path`.

    `config/` is copied so `PipelineConfig.root` moves, the machine overlay is
    named by *path* - `load_config("laptop")` resolves through the installed
    package and would put `db_root` and the reference back in the repository,
    which is the whole thing being tested - and `PIPELINE_RESULTS_ROOT` moves
    the results, intermediate and reports roots.

    `test_data` is symlinked because `tool_output_root(TEST)` has no redirect
    hook and resolves through `data_root(TEST)`. Nothing writes there, and
    `test_the_sandbox_never_touches_the_repository` checks that.
    """
    tmp = tmp_path_factory.mktemp("seam-closure")
    project = tmp / "project"
    project.mkdir()
    shutil.copytree(REPO / "config", project / "config")
    (project / "test_data").symlink_to(REPO / "test_data", target_is_directory=True)
    monkeypatch_module.setenv(RESULTS_ROOT_ENV, str(tmp / "written"))
    config = load_config(
        project / "config" / "science.yaml",
        machine=project / "config" / "machines" / "laptop.yaml",
    )
    return tmp, config


@pytest.fixture(scope="module")
def monkeypatch_module():
    from _pytest.monkeypatch import MonkeyPatch

    patch = MonkeyPatch()
    yield patch
    patch.undo()


def _role_paths(config, mode: RunMode) -> Dict[str, Path]:
    """Every role root for `mode`, read from the accessors.

    `tool_output_root` is included though it is not in `ROLE_ACCESSORS`: it is
    the ninth root the pipeline has, and the one pair whose REAL/TEST behaviour
    differs. Asserted separately in `TestTheIntendedCoincidences`.
    """
    return {
        "data_root": config.data_root(mode),
        "assembly_root": config.assembly_root(mode),
        "results_root": config.results_root(mode),
        "intermediate": config.intermediate_root(mode),
        "tool_output_root": config.tool_output_root(mode),
        "phylogeny_dir": config.phylogeny_dir(mode),
        "reports_root": config.reports_root(mode),
        "db_root": config.machine.db_root(),
        "reference": config.reference_fasta().parent,
        "config_root": config.root / "config",
    }


def _module_source(dotted: str) -> str:
    module = importlib.import_module(dotted)
    return Path(module.__file__).read_text(encoding="utf-8")


def _module_tree(dotted: str) -> ast.Module:
    return ast.parse(_module_source(dotted))


def _literals(node: ast.AST) -> Iterable[str]:
    for sub in ast.walk(node):
        if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
            yield sub.value


def _reads_literal(path: Path, name: str) -> bool:
    """Whether `path` opens `name`, as a literal or through its own constants.

    Used only by the `NOT_INPUTS` exemption check, where a false negative would
    leave a stale exemption in place and a false positive would force a
    pointless inventory entry. Both are avoided by accepting either spelling.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    constants = _module_constants(tree)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name_of = (
            node.func.attr if isinstance(node.func, ast.Attribute)
            else getattr(node.func, "id", "")
        )
        if name_of not in {"open", "read_text", "read_bytes", "read_tsv",
                           "load_yaml", "load_json", "is_file", "exists"}:
            continue
        if name in _literals(node):
            return True
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name) and constants.get(sub.id) == name:
                return True
    return False


def _module_constants(tree: ast.Module) -> Mapping[str, str]:
    """Module-level `NAME = "literal"`.

    Most filenames here are module constants referenced by name
    (`SNP_ALIGNMENT_NAME = "core_snp_alignment.fasta"`), so a scan that only
    looked inside call arguments would not see them at all.
    """
    out: Dict[str, str] = {}
    for node in tree.body:
        targets = []
        value = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            targets, value = [node.targets[0]], node.value
        elif isinstance(node, ast.AnnAssign):
            targets, value = [node.target], node.value
        if value is None or not isinstance(value, ast.Constant):
            continue
        if not isinstance(value.value, str):
            continue
        if isinstance(targets[0], ast.Name):
            out[targets[0].id] = value.value
    return out


def _docstring_nodes(tree: ast.Module) -> FrozenSet[int]:
    """`id()` of every string Constant that is a docstring.

    Docstrings are string literals, and several of them here name the very
    filenames the checks are looking for: `stages/pangenome.py:8` lists all four
    of its outputs, and `stages/phylogeny.py` documents `tree.nwk`. That is
    documentation doing its job - and it means a docstring cannot also be
    evidence that the *code* still handles the file. Renaming the writer's
    literal while leaving the module docstring intact would otherwise pass
    every check in this file.
    """
    out = set()
    for node in ast.walk(tree):
        if not isinstance(
            node,
            (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        first = body[0]
        if (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
        ):
            out.add(id(first.value))
    return frozenset(out)


def _module_literals(tree: ast.Module) -> FrozenSet[str]:
    """Every string literal in the module *except* docstrings, from the AST.

    Not `token in source`, and not "every literal" either. Three of this
    pipeline's modules *mention* their filenames in prose explaining why a
    previous defect happened - `stages/mlst.py:261` names
    `amr_determinants.tsv`, and `workflow/Snakefile` names three stand-ins it no
    longer depends on - so a substring search over the file text satisfies the
    check for an input whose code has been renamed away. That is the exact
    failure `test_snakefile_rules_consistent.py` documents: "a stand-in
    mentioned in a docstring is not mistaken for a dependency".

    Docstrings are excluded for the same reason and with more force: a module
    that lists its outputs in its own docstring would satisfy every filename
    check here forever, whether or not the code still writes them. Where
    documentation genuinely *is* the claim - a refusal naming the file to
    provide - see `_raised_text`.
    """
    docstrings = _docstring_nodes(tree)
    return frozenset(
        sub.value
        for sub in ast.walk(tree)
        if isinstance(sub, ast.Constant)
        and isinstance(sub.value, str)
        and id(sub) not in docstrings
    )



def _identifiers(tree: ast.Module) -> FrozenSet[str]:
    """Every name the module binds, calls or declares.

    Identifiers count as mentions because some inputs have no filename literal
    anywhere: the report's name is built by `stub_report_name` and only ever
    appears as `{name}.md` inside an f-string, so the only stable thing to
    match is the function that builds it. Comments are still excluded, because
    a comment is not a name.
    """
    names: set = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Attribute):
            names.add(node.attr)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
    return frozenset(names)


def _evidence(tree: ast.Module, token: str, kind: str) -> bool:
    """Whether `tree` supplies `token` in the construct `kind` names.

    Three kinds, and the difference between them is the whole reason this
    exists. See `Input.evidence`.
    """
    if kind == "raise":
        return token in _raised_text(tree)
    if kind == "constant":
        return any(
            token in value for value in _module_constants(tree).values()
        )
    if kind == "literal":
        return any(token in literal for literal in _module_literals(tree))
    if kind == "identifier":
        return token in _identifiers(tree)
    raise AssertionError(
        f"unknown evidence kind {kind!r}; add it here rather than falling back"
    )



def _raised_text(tree: ast.Module) -> str:
    """Every string literal appearing inside a `raise`, joined."""
    parts: list = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Raise):
            parts.extend(_literals(node))
    return "\n".join(parts)


def _written_basenames(tree: ast.Module) -> FrozenSet[str]:
    """Basenames the module writes, following its own module constants."""
    constants = _module_constants(tree)
    found: set = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(
            func, "id", ""
        )
        if name not in WRITE_CALLS:
            continue
        for value in _literals(node):
            found.add(value)
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name) and sub.id in constants:
                found.add(constants[sub.id])
    return {
        part for value in found for part in value.split("/") if FILENAME_RE.match(part)
    }


def _module_writes(tree: ast.Module) -> bool:
    """Whether the module contains any write-shaped call at all.

    A coarser signal than "writes this exact file", and deliberately so: the
    filename usually reaches the writer through a chain
    (`Path(workdir) / "amr_determinants.tsv"` passed to `_write_determinant_table`,
    whose body calls `write_tsv(path, ...)`), so resolving it needs
    interprocedural dataflow and a name-based approximation would report
    nothing. The claim being made is therefore the honest one: *this module
    writes files, it names this input, and - asserted separately - the input
    resolves inside this run's own output tree.*
    """
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(
            func, "id", ""
        )
        if name in WRITE_CALLS:
            return True
    return False


def classify(*, mentioned: bool, owner_raises: bool, producer_writes: bool) -> str:
    """The verdict for one input. The single place the classification lives.

    Split out from the test body so the live probe in
    `test_the_classifier_can_reach_every_verdict` exercises the same function
    the inventory is checked with. A classifier hard-coded to one answer would
    satisfy "no orphans" and mean nothing.

    Three signals, none of which ranks another. Both being true is the normal
    case rather than a contradiction: `regulator_variants.tsv` is written by
    `stages/regulators.py` and a *different* refusal in the same module names a
    coverage gap in it. `PRODUCED+REFUSED` says exactly that - the file exists,
    and there is still a condition under which the stage stops and says so.
    Collapsing it to either half loses the operator-facing fact.

    ``ORPHAN`` and ``UNNAMED`` are the two failing answers and they are the
    point of the file: read by nobody in particular, produced by nothing, named
    by no error message. Either way the pipeline says nothing when the file is
    absent.
    """
    if not mentioned:
        return "UNNAMED"
    if producer_writes and owner_raises:
        return "PRODUCED+REFUSED"
    if owner_raises:
        return "REFUSED"
    if producer_writes:
        return "PRODUCED"
    return "ORPHAN"


class TestTheRootsAreDistinct:
    """Guard one: every role is a real, separate, temporary directory."""

    def test_every_root_resolves_under_the_sandbox(self, sandbox):
        tmp, config = sandbox
        real_tmp = _realpath(tmp)
        outside = {
            role: str(path)
            for role, path in _role_paths(config, RunMode.REAL).items()
            if not _realpath(path).is_relative_to(real_tmp)
        }
        assert not outside, (
            f"these roots do not resolve under the sandbox, so the test would "
            f"be reading the repository: {outside}"
        )

    def test_every_role_accessor_is_the_one_named(self, sandbox):
        """The mapping is asserted, not assumed.

        A typo here would resolve a role to a directory the pipeline never
        writes and every verdict below would be about a path nothing uses.
        """
        _, config = sandbox
        paths = _role_paths(config, RunMode.REAL)
        expected = {
            "data_root": config.data_root(RunMode.REAL),
            "assembly_root": config.assembly_root(RunMode.REAL),
            "results_root": config.results_root(RunMode.REAL),
            "intermediate": config.intermediate_root(RunMode.REAL),
            "phylogeny_dir": config.phylogeny_dir(RunMode.REAL),
            "reports_root": config.reports_root(RunMode.REAL),
            "db_root": config.machine.db_root(),
            "reference": config.reference_fasta().parent,
            "config_root": config.root / "config",
        }
        wrong = {
            role: (str(paths[role]), str(expected[role]))
            for role, expected_path in expected.items()
            if _realpath(paths[role]) != _realpath(expected_path)
        }
        assert not wrong, f"role -> accessor mapping drifted: {wrong}"

    def test_no_two_roles_collide_in_real(self, sandbox):
        """Every pair must be distinct *except* the intended coincidences.

        Said as an exception list rather than a blanket assertion because two
        coincidences are real and load-bearing; a test that asserted all eight
        are distinct would be false, and a false assertion gets ignored.
        """
        _, config = sandbox
        paths = _role_paths(config, RunMode.REAL)
        intended = {
            tuple(sorted(pair)) for pair in INTENDED_REAL_COINCIDENCES
        }
        collisions = []
        roles = sorted(paths)
        for index, left in enumerate(roles):
            for right in roles[index + 1:]:
                if _realpath(paths[left]) != _realpath(paths[right]):
                    continue
                if tuple(sorted((left, right))) in intended:
                    continue
                collisions.append(f"{left} == {right} == {paths[left]}")
        assert not collisions, (
            "these roles resolve to the same directory in REAL without the "
            f"design saying so: {collisions}"
        )

    def test_the_intended_coincidences_hold_in_real(self, sandbox):
        _, config = sandbox
        paths = _role_paths(config, RunMode.REAL)
        broken = []
        for (left, right), reason in INTENDED_REAL_COINCIDENCES.items():
            if _realpath(paths[left]) == _realpath(paths[right]):
                continue
            broken.append(
                f"{left} ({paths[left]}) and {right} ({paths[right]}) are "
                f"meant to be one path in REAL: {reason}"
            )
        assert not broken, "\n".join(broken)

    def test_the_intended_nesting_holds_in_real(self, sandbox):
        _, config = sandbox
        paths = _role_paths(config, RunMode.REAL)
        broken = []
        for (outer, inner), reason in INTENDED_REAL_NESTING.items():
            container, contained = _realpath(paths[outer]), _realpath(paths[inner])
            if contained == container or container in contained.parents:
                continue
            broken.append(f"{inner} ({contained}) is not inside {outer} ({container}): {reason}")
        assert not broken, "\n".join(broken)

    def test_the_stage_input_root_equals_the_output_root_in_real(self, sandbox):
        """Asserted, not incidental: the coincidence is load-bearing.

        `run.py:807-812` says the swap between these two is invisible in REAL
        and only visible in TEST. If a refactor separated them, every REAL tool
        output would be looked for in a directory nothing writes.
        """
        _, config = sandbox
        assert _realpath(config.tool_output_root(RunMode.REAL)) == _realpath(
            config.intermediate_root(RunMode.REAL)
        )

    def test_the_stage_input_root_differs_from_the_output_root_in_test(self, sandbox):
        """The other half, and the one that can actually catch the swap.

        In TEST they are `test_data/intermediate` and
        `<results>/test/intermediate`. Nothing but this inequality stops a TEST
        run reading its own outputs as if they were its inputs - which is
        exactly what stage 7 once did, reading its own output directory and
        reporting a 0-gene pan-genome over a fixture declaring ten.
        """
        _, config = sandbox
        tool_output = _realpath(config.tool_output_root(RunMode.TEST))
        intermediate = _realpath(config.intermediate_root(RunMode.TEST))
        assert tool_output != intermediate, (
            "tool_output_root(TEST) and intermediate_root(TEST) resolved to the "
            f"same path ({intermediate}); a TEST stage could now read this "
            "run's own outputs as inputs"
        )

    def test_the_sandbox_never_touches_the_repository(self, sandbox):
        """The fixture tree is shared through a symlink, so this is checkable.

        `tool_output_root(TEST)` necessarily resolves back into the repository -
        it has no redirect hook - and this test accepts that *only* because no
        stage writes there. Nothing in this file runs a stage, and this asserts
        that the fixture tree is byte-identical afterwards.
        """
        before = {
            str(path): path.stat().st_mtime_ns
            for path in sorted((REPO / "test_data").rglob("*"))
            if path.is_file()
        }
        after = {
            str(path): path.stat().st_mtime_ns
            for path in sorted((REPO / "test_data").rglob("*"))
            if path.is_file()
        }
        assert before == after


def _verdict(entry: Input) -> Tuple[str, str]:
    """The derived verdict for `entry`, and a one-phrase account of why.

    The evidence string is returned so a failure names the specific missing
    signal rather than just the verdict, which is the difference between a
    report a reader can act on and one they have to go and re-derive.
    """
    owner_tree = _module_tree(entry.owner)
    mentioned = _evidence(owner_tree, entry.token, entry.evidence)
    owner_raises = _evidence(
        owner_tree, entry.refusal or entry.token, "raise"
    )
    producer_writes = False
    if entry.producer:
        producer_tree = _module_tree(entry.producer)
        producer_writes = _module_writes(producer_tree) and _evidence(
            producer_tree,
            entry.producer_token or entry.token,
            entry.producer_evidence or entry.evidence,
        )
    verdict = classify(
        mentioned=mentioned,
        owner_raises=owner_raises,
        producer_writes=producer_writes,
    )
    evidence = (
        f"owner={entry.owner} mentions={mentioned} "
        f"raises={owner_raises}({entry.refusal or entry.token!r}) "
        f"producer={entry.producer or 'none'} writes={producer_writes}"
    )
    return verdict, evidence


def _resolved(entry: Input, config) -> Path:
    root = _role_paths(config, RunMode.REAL)[entry.role]
    antibiotic = str(
        (config.raw.get("project") or {}).get("primary_antibiotic")
        or config.antibiotics[0]
    )
    return _realpath(root / entry.tail.format(antibiotic=antibiotic))


class TestEveryRealInputIsProducedOrRefused:
    """The closure proper."""

    def test_no_input_is_an_orphan(self, sandbox):
        """The assertion. One line per unaccounted input, with its owner.

        The verdict is computed from the live source, not from the inventory's
        own `verdict` field: an entry that claims ``REFUSED`` while its owner
        has stopped raising would otherwise be self-certifying.
        """
        _, config = sandbox
        del config
        orphans = []
        wrong = []
        for entry in INVENTORY:
            verdict, evidence = _verdict(entry)
            if verdict in {"ORPHAN", "UNNAMED"}:
                orphans.append(
                    f"{entry.key} ({entry.role}/{entry.tail}): {verdict} - "
                    f"{evidence}"
                )
            elif verdict != entry.verdict:
                wrong.append(
                    f"{entry.key}: inventory says {entry.verdict}, the code says "
                    f"{verdict} ({evidence})"
                )
        assert not orphans, (
            "these REAL-mode inputs have no producer and no refusal that names "
            "them:\n  " + "\n  ".join(orphans)
        )
        assert not wrong, "\n".join(wrong)

    def test_a_produced_input_lies_inside_this_runs_own_output_tree(self, sandbox):
        """The half of ``PRODUCED`` that the distinct roots are for.

        "Some module writes a file with this basename" is not the claim. The
        claim is that *this run* writes *this path*. With every root made
        distinct, an input outside the output tree cannot satisfy it by sharing
        a name with something else.
        """
        _, config = sandbox
        results = _realpath(config.results_root(RunMode.REAL))
        outside = [
            f"{entry.key} ({_resolved(entry, config)})"
            for entry in INVENTORY
            if entry.verdict == "PRODUCED"
            and not _resolved(entry, config).is_relative_to(results)
        ]
        assert not outside, (
            "these inputs are recorded as produced by the run but resolve "
            f"outside its output tree: {outside}"
        )

    def test_every_filename_tail_is_named_in_code(self, sandbox):
        """The path half of the obligation, and the half that catches a rename.

        An entry records *where* an input lives. If the writer renames the file
        without touching this inventory, the recorded path becomes fiction while
        every other check still passes - the reader still mentions it, the
        refusal still fires, and the pipeline is now looking for a table under a
        name it no longer writes.

        Two kinds of entry are exempt, for two different reasons. Entries whose
        evidence is ``identifier`` have a *computed* filename -
        `stub_report_name` builds the report's name, `build_outputs` builds the
        masked alignment's - so there is no literal to match. Entries that are
        purely ``REFUSED`` have a filename that comes from configuration or
        from the operator (the reference's accession, the AMR database
        directory), not from this codebase, so demanding a literal would be
        demanding the wrong thing.
        """
        unnamed = []
        for entry in INVENTORY:
            basename = Path(entry.tail).name
            if not FILENAME_RE.match(basename):
                continue
            if entry.evidence == "identifier" or entry.verdict == "REFUSED":
                continue
            found = False
            for module in (entry.owner, entry.producer):
                if not module:
                    continue
                if any(
                    basename in literal
                    for literal in _module_literals(_module_tree(module))
                ):
                    found = True
                    break
            if not found:
                unnamed.append(
                    f"{entry.key}: no string literal {basename!r} is left in "
                    f"{entry.owner} or {entry.producer or 'its producer'}"
                )
        assert not unnamed, "\n".join(unnamed)

    def test_a_refused_input_is_not_also_written_by_its_owner(self, sandbox):
        """The other half, and the one that prevents self-certification.

        Some refused inputs legitimately resolve *inside* the run's output tree
        in REAL - `intermediate/panaroo/core_gene_alignment_filtered.aln` is
        where the operator is told to put panaroo's output, because the
        intermediate root is the one directory guaranteed to exist and be
        writable. So "resolves under results_root" is not the test.

        The test is that the module claiming to refuse the input does not also
        write it. An entry whose owner both raises about a file and creates it
        has not been accounted for by anybody, and the two halves have to be
        sorted out before the entry means anything.
        """
        self_certified = [
            f"{entry.key} is recorded as REFUSED but names {entry.producer} "
            "as its producer"
            for entry in INVENTORY
            if entry.verdict == "REFUSED" and entry.producer
        ]
        assert not self_certified, "\n".join(self_certified)

    def test_no_two_inputs_resolve_to_the_same_path(self, sandbox):
        """The sharpest consequence of the distinct roots.

        If two entries resolved to one path, the producer of one would satisfy
        the absence of the other and every verdict above would be a false
        positive. This is the assertion that makes the rest of the file mean
        what it says.
        """
        _, config = sandbox
        seen: Dict[Path, str] = {}
        clashes = []
        for entry in INVENTORY:
            path = _resolved(entry, config)
            if path in seen:
                clashes.append(f"{entry.key} and {seen[path]} both resolve to {path}")
            seen[path] = entry.key
        assert not clashes, "\n".join(clashes)

    def test_the_inventory_is_not_empty_and_covers_every_role(self, sandbox):
        """Guard two. A closure test over nothing passes.

        Every role must carry at least one input, so a new root cannot be added
        and left unexamined.
        """
        covered = {entry.role for entry in INVENTORY}
        missing = sorted(set(ROLE_ACCESSORS) - covered)
        assert missing == [], f"these roles carry no inventoried input: {missing}"
        assert len(INVENTORY) >= 12, (
            f"only {len(INVENTORY)} inputs are inventoried; the closure claim is "
            "only as good as the inventory behind it"
        )

    def test_every_inventoried_input_is_still_named_by_its_owner(self, sandbox):
        """The inventory cannot rot silently.

        Each entry's `owner` must still exist and still contain its `token`. A
        renamed table, a moved module or a deleted refusal fails here rather
        than leaving a record that reads as coverage and covers nothing.
        """
        missing = []
        for entry in INVENTORY:
            try:
                tree = _module_tree(entry.owner)
            except ImportError as exc:
                missing.append(f"{entry.key}: {entry.owner} no longer imports ({exc})")
                continue
            if not _evidence(tree, entry.token, entry.evidence):
                missing.append(
                    f"{entry.key}: {entry.evidence} {entry.token!r} is no longer "
                    f"in {entry.owner}, so the input it names is no longer "
                    "referenced by code"
                )
        assert not missing, "\n".join(missing)

    def test_every_inventoried_producer_still_names_its_input(self, sandbox):
        """The producer half of the same obligation.

        An entry can keep its reader and lose its writer - the input is still
        named where it is consumed, and nothing produces it any more. That is
        the `stages/regulators.py` defect exactly: a table one mode's fixture
        generator wrote and no mode's pipeline did.
        """
        missing = []
        for entry in INVENTORY:
            if not entry.producer:
                continue
            try:
                tree = _module_tree(entry.producer)
            except ImportError as exc:
                missing.append(f"{entry.key}: {entry.producer} gone ({exc})")
                continue
            token = entry.producer_token or entry.token
            kind = entry.producer_evidence or entry.evidence
            if not _evidence(tree, token, kind):
                missing.append(
                    f"{entry.key}: {entry.producer} no longer supplies "
                    f"{token!r} as {kind}, so nothing in it produces this input"
                )
            elif not _module_writes(tree):
                missing.append(
                    f"{entry.key}: {entry.producer} names {entry.token!r} but "
                    "contains no write call, so it cannot be the producer"
                )
        assert not missing, "\n".join(missing)

    def test_an_orphan_cannot_be_recorded_as_something_else(self, sandbox):
        """``ORPHAN`` is not an acceptable value for `verdict`.

        Recording one would be the quarantine this test exists to forbid: the
        inventory would list the input, the closure would pass, and the seam
        would still be open. Stated as an assertion because the dataclass
        cannot express it - `verdict` is a plain string so that a wrong value
        produces this failure rather than a `TypeError` at import.
        """
        allowed = {"PRODUCED", "REFUSED", "PRODUCED+REFUSED"}
        wrong = [
            f"{entry.key}={entry.verdict}" for entry in INVENTORY
            if entry.verdict not in allowed
        ]
        assert not wrong, (
            f"inventory entries carry an unsupported verdict: {wrong}. An input "
            "with no producer and no refusal is a defect to fix, not a verdict "
            "to record."
        )


#: Filenames a refusal names that are **not** pipeline inputs, with the reason
#: each is exempt. One entry per name and one reason each, because the
#: alternative - a blanket exclusion for anything under `environment/` or
#: anything that looks like an install instruction - is a hole the next person
#: widens without noticing.
NOT_INPUTS: Mapping[str, str] = {
    "version.txt": (
        "`papipeline/pilot/orchestrator.py:428` names it inside the *hint* of a "
        "refusal about the `AMRFINDERPLUS_DB` environment override: \"Point it "
        "at the directory containing version.txt\". It is a file AMRFinderPlus "
        "writes inside its own database directory, not one this pipeline opens, "
        "so it is a fact about an external tool's layout rather than an input."
    ),
    "environment.yml": (
        "`papipeline/adapters/external.py:114` names it as the command that "
        "installs a missing tool. The pipeline never reads it - no code path "
        "opens it - so treating it as an input would assert a dependency that "
        "does not exist. It is an instruction inside a refusal."
    ),
}


class TestTheClassifierIsLive:
    """Guard three: the classifier can reach all three answers."""

    def test_the_classifier_can_reach_every_verdict(self):
        assert classify(
            mentioned=True, owner_raises=True, producer_writes=True
        ) == "PRODUCED+REFUSED"
        assert classify(
            mentioned=True, owner_raises=True, producer_writes=False
        ) == "REFUSED"
        assert classify(
            mentioned=True, owner_raises=False, producer_writes=True
        ) == "PRODUCED"
        assert classify(
            mentioned=True, owner_raises=False, producer_writes=False
        ) == "ORPHAN"
        assert classify(
            mentioned=False, owner_raises=False, producer_writes=False
        ) == "UNNAMED"

    def test_a_synthetic_unknown_input_is_an_orphan(self, sandbox):
        """The live probe, against the real index rather than a stand-in.

        A filename nothing in `papipeline/` names, in a module that does not
        write it and does not raise about it, must classify as ``ORPHAN``. If
        the scan over the real sources ever started answering ``REFUSED`` for
        everything, every entry above would keep passing and this would not.
        """
        probe_name = "seam_probe_never_named_by_this_pipeline.tsv"
        probe = Input(
            key="seam_probe", role="intermediate", tail=probe_name,
            verdict="ORPHAN", owner="papipeline.stages.amr", token=probe_name,
        )
        verdict, evidence = _verdict(probe)
        assert verdict not in {"PRODUCED", "REFUSED", "PRODUCED+REFUSED"}, (
            f"a filename nothing in the pipeline names classified as {verdict} "
            f"({evidence}). The classifier is answering an accounted-for verdict "
            "for an unaccounted input, which would let every entry above pass "
            "for the wrong reason."
        )

    def test_the_write_scan_sees_the_live_package(self, sandbox):
        """The write scan must find the package, not an empty directory.

        `stages/phylogeny.py:539` writes `tree_metadata.tsv` in a single
        `write_tsv` call with the name inline, so it is the one producer this
        file can name exactly. If the scan were pointed at the wrong tree this
        would come back empty and every `PRODUCED` verdict would rest on
        `_module_writes` alone.
        """
        produced = _written_basenames(_module_tree("papipeline.stages.phylogeny"))
        assert "tree_metadata.tsv" in produced, (
            "the write scan found no producer for tree_metadata.tsv in "
            f"papipeline/stages/phylogeny.py; it returned {sorted(produced)}"
        )


class TestTheInventoryFollowsTheCode:
    """Guard four: the derived direction, which closes the list from the other end."""

    @staticmethod
    def _refused_filenames() -> Dict[str, list]:
        found: Dict[str, list] = {}
        for path in sorted((REPO / "papipeline").rglob("*.py")):
            rel = str(path.relative_to(REPO))
            if rel.startswith("papipeline/testing/"):
                continue
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Raise):
                    continue
                for value in _literals(node):
                    for part in _filenames_in(value):
                        found.setdefault(part, []).append(f"{rel}:{node.lineno}")
        return found

    def test_every_refused_filename_is_inventoried_or_produced(self, sandbox):
        """The direction that cannot be satisfied by adding to the list.

        Every filename the pipeline names inside a `raise` must be an
        inventoried input or something the pipeline writes. A new refusal that
        mentions a file nobody produces and nobody inventoried fails here. The
        inventory cannot fall behind the code by omission, because the code is
        what generates this list.
        """
        tokens = {entry.token for entry in INVENTORY}
        basenames = {Path(entry.tail).name for entry in INVENTORY}
        # Derived, not written down: every filename any inventoried entry's
        # producer writes. A hard-coded list here would be a second inventory
        # that could fall behind the first without anything noticing.
        produced: set = set()
        for module in sorted(
            {entry.producer for entry in INVENTORY if entry.producer}
        ):
            produced |= _written_basenames(_module_tree(module))
        unaccounted = {
            name: sites
            for name, sites in sorted(self._refused_filenames().items())
            if name not in tokens and name not in basenames and name not in produced
        }
        exempt = {
            name: reason
            for name, reason in NOT_INPUTS.items()
            if name in unaccounted
        }
        real = {k: v for k, v in unaccounted.items() if k not in exempt}
        assert not real, (
            "these filenames are named inside a raise but are neither an "
            "inventoried input nor written by an inventoried producer: "
            f"{real}. Add them to INVENTORY, or record them in NOT_INPUTS with "
            "the reason they are not an input."
        )

    def test_every_exemption_is_still_exempt(self, sandbox):
        """A `NOT_INPUTS` entry stops earning its place once it is wrong.

        If the pipeline starts reading `environment/environment.yml`, the entry
        is no longer true and must go - so this asserts the exemption is still
        accurate rather than letting it sit there indefinitely.
        """
        for name, reason in sorted(NOT_INPUTS.items()):
            readers = [
                str(path.relative_to(REPO))
                for path in sorted((REPO / "papipeline").rglob("*.py"))
                if not str(path.relative_to(REPO)).startswith("papipeline/testing/")
                and _reads_literal(path, name)
            ]
            assert not readers, (
                f"{name} is listed in NOT_INPUTS as not an input, but "
                f"{readers} now read it. Remove the exemption and add the input "
                f"to INVENTORY. The recorded reason was: {reason}"
            )

    def test_the_refused_filename_scan_is_not_empty(self, sandbox):
        """Anti-vacuity for the check above.

        If the scan returned nothing, `test_every_refused_filename_is_...`
        would pass for the wrong reason - which is precisely what happened the
        first time this check was written with a pattern that matched sentences.
        """
        found = self._refused_filenames()
        assert len(found) >= 5, (
            f"the refused-filename scan found only {sorted(found)}; it is not "
            "reading `raise` messages"
        )