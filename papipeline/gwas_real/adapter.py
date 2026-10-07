"""The real-path pyseer adapter: flags, gates, and the two-pass protocol.

This is the REAL counterpart to ``stages/gwas.py``'s ``PyseerEngine``, and it
is file-based on purpose: the real path is given a phenotype file, a variant
file (``--kmers``/``--pres``/``--vcf``), the stage-10 kinship matrix and, when
lineage enters the model, stage 10's distance matrix - not the feature table
and Gram-matrix kinship the TEST engine writes for itself.

Every flag below was checked against the installed pyseer 1.1.2
(``pyseer --help`` and ``pyseer/__main__.py``); none is invented:

* ``--lmm`` with ``--similarity`` is the mixed model, and pyseer *requires*
  that pairing (``__main__.py:216``: "Must use distance matrix with fixed
  effects, or similarity matrix with random effects").
* ``--lineage`` additionally requires ``--distances`` (``__main__.py:220``:
  "Must also provide a distance matrix to report lineage effects"), which is
  why lineage here consumes stage 10's ``similarity.tsv`` - and
  ``pyseer/input.py::load_structure`` reads it with
  ``pd.read_table(index_col=0)``, exactly the shape stage 10 writes.
* ``--burden`` is a *value* (the VCF regions file) and refuses to run without
  ``--vcf`` (``__main__.py:212``).
* ``--covariates`` needs ``--use-covariates``, whose arguments are separate
  argv tokens (``nargs='*'``).

Two refusals here are about pyseer's silence rather than its errors. pyseer
**silently intersects** the phenotype against the similarity matrix
(``pyseer/lmm.py::initialise_lmm``: ``p.index.intersection(K.index)``) - so a
cohort mismatch would analyse a different set of samples and still exit 0 - and
it **silently switches to a continuous linear model** the moment any phenotype
value is not 0/1 (``__main__.py:244``: "Detected continuous phenotype"). Both
are checked before a subprocess exists.

The significance threshold is step 12a's, reused rather than reimplemented:
``stages.gwas.reduce_unique_patterns`` runs the vendored
``scripts/gwas/count_patterns.py`` over pyseer's own ``--output-patterns``
file, and its ``alpha / n_unique_patterns`` is what pass 2 gates on. This
module computes no threshold of its own.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, List, Mapping, Optional, Sequence, Tuple

from ..errors import DataContractError
from ..execution.contracts import STAGE_TABLES
from ..io.tsv import write_tsv
from ..logging_utils import get_logger
from ..manifest import SampleManifest
from ..models import RunMode
from ..stages import gwas as stage_gwas

from .settings import AdapterSettings, VARIANT_SOURCES

LOGGER = get_logger("gwas_real.adapter")

#: Stage 10's distance matrix, under the name ``execution.contracts``
#: declares for it (``table_path``: "The single source of that location").
#: Read from the contract rather than retyped here: a refusal message must
#: name the file the pipeline actually writes, and a literal in this module
#: would be a second copy of that name, free to drift from the first.
STAGE10_DISTANCES_FILENAME = STAGE_TABLES["similarity"][0]

#: pyseer's phenotype column, which the phenotype file this path reads declares.
PHENOTYPE_COLUMN = "phenotype"

#: Opt-in switch for launching the real binary in TEST mode. Named here and in
#: the tests as the same literal string; spec.md D8 makes the default "no real
#: tools in TEST", and an env var is the explicit crossing of that line.
OPT_IN_ENV = "PAPIPELINE_TEST_PYSEER"

#: What the adapter writes beside the association table.
PROVENANCE_FILENAME = "pyseer_real_provenance.tsv"

#: pyseer's own name for the lineage-effects table (``--lineage-file``'s
#: documented default), written into the *workdir* rather than the directory
#: the process started in.
LINEAGE_FILENAME = "lineage_effects.txt"

#: How pyseer is invoked: the interpreter runs the compatibility shim, for the
#: reason ``stages/gwas.py::_base_command`` documents at length - the pinned
#: 1.1.2 calls scipy/statsmodels names that no longer exist, and the shim
#: rebinds them to the identical functions without changing any fitted value.
#: One shim, two call sites, so the two cannot diverge.
_PYSEER_COMPAT = stage_gwas._PYSEER_COMPAT

#: The directory whose ``sitecustomize.py`` applies that same shim in every
#: interpreter started for a pass - including the ``multiprocessing`` workers
#: pyseer spawns when ``--cpu > 1``. On macOS the start method is ``spawn``,
#: so a worker re-imports the compat script's *top level* (which defines
#: ``apply_shims`` without calling it) and would otherwise reach
#: ``pyseer/model.py::fit_lineage_effect``'s ``smf.Logit`` unarmed. See
#: ``papipeline/gwas_real/_pyseer_shim/sitecustomize.py`` for the full
#: chain; the logic lives in the compat script and is only *delivered* here.
SHIM_ENV_DIR = Path(__file__).resolve().parent / "_pyseer_shim"

#: The invocation signature the injected TEST invoker must satisfy:
#: ``(argv, table_path, log_path) -> None``. Mirrors what pyseer itself does -
#: association table on stdout, progress on stderr.
Invoker = Callable[[Sequence[str], Path, Path], None]


# --------------------------------------------------------------------------
# The command
# --------------------------------------------------------------------------


def build_command(
    *,
    phenotype: Path,
    variant_source: str,
    variants: Path,
    similarity: Path,
    threads: Optional[int],
    min_af: float,
    max_af: float,
    covariates: Optional[Path] = None,
    use_covariates: Sequence[str] = (),
    distances: Optional[Path] = None,
    lineage: bool = False,
    lineage_file: Optional[Path] = None,
    burden: Optional[Path] = None,
    output_patterns: Optional[Path] = None,
    filter_pvalue: Optional[float] = None,
    lrt_pvalue: Optional[float] = None,
) -> List[str]:
    """The pyseer argv for one pass, or a refusal for a combination pyseer rejects.

    Refuses rather than delegating each impossible combination to the tool
    because every one of them exits in ``pyseer/__main__.py`` *before* any
    analysis - an argparse-shaped failure recorded in a log, instead of a
    statement about the configuration that caused it.

    Args:
        threads: From ``config.threads`` (the machine overlay). ``None`` is
            refused rather than defaulted: the process count is a fact about
            hardware, and pyseer's own default of 1 would be a silent guess
            about this machine.
        min_af / max_af: Frequency bounds, derived by the caller from the
            cohort size (``1/n`` and ``1 - 1/n``), same rule as the stage.
        lrt_pvalue / filter_pvalue: Per-pass gates. Step 12a passes both open
            (1) so every tested variant reaches the patterns file; step 12b
            passes the derived threshold.
    """
    source = str(variant_source).strip().lower()
    if source not in VARIANT_SOURCES:
        raise DataContractError(
            f"Unknown variant source {variant_source!r}; pyseer accepts exactly "
            f"one of {', '.join(VARIANT_SOURCES)} (its variant group is "
            "mutually exclusive and required).",
            variant_source=str(variant_source),
            supported=",".join(VARIANT_SOURCES),
        )
    if burden is not None and source != "vcf":
        raise DataContractError(
            "pyseer refuses `--burden` without `--vcf` (__main__.py: "
            "'burden testing (requires --vcf)'), so a burden run needs "
            "variant_source=vcf; grouping variants by region is only defined "
            "over a VCF's sites.",
            variant_source=source,
            burden=str(burden),
        )
    if lineage and distances is None:
        raise DataContractError(
            "pyseer refuses `--lineage` without `--distances` (__main__.py: "
            "'Must also provide a distance matrix to report lineage effects'), "
            f"and the matrix stage 10 writes ({STAGE10_DISTANCES_FILENAME}) is "
            "the one this path uses. Pass distances_path=, or turn lineage off.",
            key="gwas_real.lineage",
        )
    if distances is not None and not lineage:
        raise DataContractError(
            "pyseer refuses `--lmm --distances` unless `--lineage` is present "
            "(__main__.py: 'Must use distance matrix with fixed effects, or "
            "similarity matrix with random effects / Unless performing a "
            "lineage analysis with random effects'), so a distances matrix "
            "travels only with lineage=True.",
            distances=str(distances),
        )
    if lineage_file is not None and not lineage:
        raise DataContractError(
            "--lineage-file names where pyseer writes the lineage-effects "
            "table, which it only writes under --lineage; without the model "
            "term the flag would change nothing.",
            lineage_file=str(lineage_file),
        )
    if threads is None:
        raise DataContractError(
            "runtime.threads is not set, so pyseer has no process count to be "
            "given (--cpu). It is a machine fact and belongs in the machine "
            "overlay; inferring one here would be a guess about hardware "
            "nobody on this machine described.",
            key="runtime.threads",
        )
    if covariates is None and use_covariates:
        raise DataContractError(
            "--use-covariates names columns of a file that was not given; pass "
            "--covariates (covariates= here) as well, or drop the column list. "
            "pyseer would ignore the list entirely ('Default: load covariates "
            "but don't use them'), and a flag that changes nothing is worse "
            "than an error.",
            use_covariates=" ".join(str(t) for t in use_covariates),
        )
    if covariates is not None and not use_covariates:
        raise DataContractError(
            "--covariates loads a file pyseer then uses no column of unless "
            "--use-covariates names one (__main__.py: 'Default: load "
            "covariates but don't use them'); either both or neither.",
            covariates=str(covariates),
        )

    command: List[str] = [
        sys.executable,
        str(_PYSEER_COMPAT),
        "--phenotypes", str(phenotype),
        "--phenotype-column", PHENOTYPE_COLUMN,
        f"--{source}", str(variants),
        "--similarity", str(similarity),
        "--lmm",
        # Bounds from the cohort size, never from a configured constant: a
        # feature carried by fewer isolates cannot be tested, and a fixed floor
        # would drift with n. Mirrors PyseerEngine._base_command.
        "--min-af", f"{float(min_af):.6f}",
        "--max-af", f"{float(max_af):.6f}",
        "--cpu", str(max(1, int(threads))),
    ]
    if covariates is not None:
        command += [
            "--covariates", str(covariates),
            "--use-covariates", *[str(t) for t in use_covariates],
        ]
    if distances is not None:
        command += ["--distances", str(distances)]
    if lineage:
        command += ["--lineage"]
        # Stated, never left at pyseer's default: the default writes the
        # table into whatever directory the process happened to start in,
        # which for a test or a snakemake worker is the repository root, not
        # this run's workdir.
        if lineage_file is not None:
            command += ["--lineage-file", str(lineage_file)]
    if burden is not None:
        command += ["--burden", str(burden)]
    if output_patterns is not None:
        command += ["--output-patterns", str(output_patterns)]
    if filter_pvalue is not None:
        command += ["--filter-pvalue", repr(float(filter_pvalue))]
    if lrt_pvalue is not None:
        command += ["--lrt-pvalue", repr(float(lrt_pvalue))]
    return command


# --------------------------------------------------------------------------
# Input reading and the checks pyseer will not make
# --------------------------------------------------------------------------


def _read_phenotype(path: Path) -> List[Tuple[str, int]]:
    """The binary phenotype, in file order: ``[(sample, 0|1), ...]``.

    Refuses three things pyseer would accept or silently reinterpret: a missing
    file, a value outside {0, 1} (which flips pyseer to a *continuous* model at
    exit 0 - ``__main__.py:242-249``), and a duplicated sample id.
    """
    path = Path(path)
    if not path.is_file():
        raise DataContractError(
            f"Binary phenotype file not found at {path}. The real pyseer path "
            "reads a two-column TSV (`sample`, `phenotype`) whose values are "
            "0/1 for the analysed cohort - produced from the stage-11 calls by "
            "config.gwas.outcome (positive/negative/excluded), with "
            "non-binary calls excluded rather than folded in.",
            path=str(path),
        )
    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not lines:
        raise DataContractError(
            f"Phenotype file at {path} holds no rows, so there is no outcome "
            "to associate.",
            path=str(path),
        )
    header = lines[0].split("\t")
    if "sample" not in header or PHENOTYPE_COLUMN not in header:
        raise DataContractError(
            f"Phenotype file at {path} must have the columns `sample` and "
            f"`{PHENOTYPE_COLUMN}`; found {', '.join(header)}.",
            path=str(path),
            columns=",".join(header),
        )
    sample_col = header.index("sample")
    value_col = header.index(PHENOTYPE_COLUMN)

    pairs: List[Tuple[str, int]] = []
    seen = set()
    for lineno, line in enumerate(lines[1:], start=2):
        parts = line.split("\t")
        if len(parts) <= max(sample_col, value_col):
            raise DataContractError(
                f"Phenotype row {lineno} of {path} has {len(parts)} field(s); "
                "the file is truncated.",
                path=str(path),
                line=lineno,
            )
        sample = parts[sample_col].strip()
        if sample in seen:
            raise DataContractError(
                f"Duplicate sample {sample!r} in the phenotype file at {path}; "
                "one sample, one outcome.",
                path=str(path),
                sample_id=sample,
            )
        seen.add(sample)
        raw = parts[value_col].strip()
        try:
            value = float(raw)
        except ValueError:
            value = None
        if value not in (0.0, 1.0):
            raise DataContractError(
                f"Sample {sample!r} has phenotype {raw!r} in {path}. pyseer's "
                "model here is binary, and pyseer silently switches to a "
                "continuous linear model when any value is not 0 or 1 (it "
                "prints 'Detected continuous phenotype' and exits 0). The "
                "outcome levels are config.gwas.outcome (positive / negative / "
                "excluded); I, SDD and ND are excluded from this file rather "
                "than mapped onto 0 or 1.",
                path=str(path),
                sample_id=sample,
                value=raw,
                key="gwas.outcome",
            )
        pairs.append((sample, int(value)))
    return pairs


def _read_matrix_ids(path: Path, what: str) -> List[str]:
    """Sample ids from a square matrix's header (``sample_id`` + one column each)."""
    path = Path(path)
    if not path.is_file():
        raise DataContractError(
            f"{what} matrix not found at {path}. The real pyseer path reads it "
            "as an input; a run cannot produce one it was never given.",
            path=str(path),
            matrix=what,
        )
    first = path.read_text(encoding="utf-8").splitlines()
    header = next((line for line in first if line.strip()), "")
    columns = header.split("\t")
    if len(columns) < 2:
        raise DataContractError(
            f"{what} matrix at {path} has no sample columns, so it is not the "
            "square matrix this path expects.",
            path=str(path),
            matrix=what,
        )
    return [c.strip() for c in columns[1:]]


def _read_covariate_ids(path: Path) -> List[str]:
    """First column of pyseer's headerless covariates file."""
    path = Path(path)
    if not path.is_file():
        raise DataContractError(
            f"Covariates file not found at {path}.", path=str(path)
        )
    ids = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            ids.append(line.split("\t")[0].strip())
    if not ids:
        raise DataContractError(
            f"Covariates file at {path} holds no rows.", path=str(path)
        )
    return ids


def _require_same_cohort(
    expected: Sequence[str], found: Sequence[str], what: str, *, extra_hint: str
) -> None:
    """Exact set equality between the cohort and a matrix's samples.

    pyseer will not make this check: ``pyseer/lmm.py::initialise_lmm`` does
    ``p.index.intersection(K.index)`` and reindexes the kernel to whatever
    survives, so a mismatch analyses a *different* cohort and still exits 0.
    The sets - not the orders - are what pyseer's reindex makes equivalent, so
    the sets are what is compared.
    """
    expected_set, found_set = set(expected), set(found)
    missing = [s for s in expected if s not in found_set]
    extra = [s for s in found if s not in expected_set]
    if not missing and not extra:
        return
    parts = []
    if missing:
        parts.append(
            f"{len(missing)} cohort sample(s) absent from the {what}: "
            + ", ".join(missing[:10])
            + (" ..." if len(missing) > 10 else "")
        )
    if extra:
        parts.append(
            f"{len(extra)} {what} sample(s) not in the cohort: "
            + ", ".join(extra[:10])
            + (" ..." if len(extra) > 10 else "")
        )
    raise DataContractError(
        "; ".join(parts)
        + f". {extra_hint}",
        missing=",".join(missing[:10]),
        extra=",".join(extra[:10]),
        matrix=what,
        n_cohort=len(expected_set),
        n_matrix=len(found_set),
    )


# --------------------------------------------------------------------------
# The run
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class RealGwasRun:
    """What one adapter run produced and with what arguments."""

    mode: RunMode
    settings: AdapterSettings
    commands: Tuple[Tuple[str, ...], ...]
    reduction: Any  # stage_gwas.UniquePatternReduction
    patterns_file: Path
    pass1_table: Path
    association_table: Path
    provenance: Path
    pyseer_version: Optional[str]

    @property
    def threshold(self) -> float:
        """``alpha / n_unique_patterns`` - step 12a's number, not this module's."""
        return self.reduction.threshold


def _pyseer_env() -> Mapping[str, str]:
    """The environment a pyseer pass runs in: this process's, plus the shim.

    ``PYTHONPATH`` is *prepended* to rather than replaced, because a caller
    (a test, a wrapper script) may already be relying on its own entries, and
    shadowing those to deliver a shim would trade one silent breakage for
    another. Nothing else is touched: threads and paths are arguments, not
    environment.
    """
    env = dict(os.environ)
    entry = str(SHIM_ENV_DIR)
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = entry if not existing else entry + os.pathsep + existing
    return env


def _launch(cmd: Sequence[str], table: Path, log: Path) -> None:
    """Run pyseer, keeping its stdout table and stderr log apart.

    The separation is the stage's rule and not a style choice: pyseer writes
    the association table to stdout and progress to stderr, and merging them
    is how a run ends up reporting zero associations while logging thousands
    of lines. A non-zero exit names the pass and the log - a bare
    ``FileNotFoundError`` from a redirected subprocess would send an operator
    looking for a broken install rather than at pyseer's own words.

    The environment carries the shim into the process *and its workers*; see
    :func:`_pyseer_env`.
    """
    import subprocess

    table.parent.mkdir(parents=True, exist_ok=True)
    if table.exists():
        table.unlink()
    LOGGER.info("pyseer: %s", " ".join(str(c) for c in cmd))
    with table.open("w") as out, log.open("w") as err:
        proc = subprocess.run(
            [str(c) for c in cmd], stdout=out, stderr=err, env=_pyseer_env()
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"pyseer failed with exit code {proc.returncode}; see {log}"
        )


def _require_pyseer(config: Any, workdir: Path) -> Tuple[str, str]:
    """Resolve pyseer by name and record the version it reports about itself.

    Reuses the stage's resolver rather than growing a second one: it refuses
    by name in ``ToolNotAvailableError``, says where a machine overlay would
    declare the tool, and caches the version subprocess so the two passes pay
    for one. The command still runs the compatibility shim - this is the
    availability and provenance check, not the loader (see module docstring).
    """
    candidates = config.tool_candidates("pyseer")
    executable = str(candidates[0]) if candidates else "pyseer"
    engine = stage_gwas.PyseerEngine(executable, Path(workdir))
    version = engine.require_executable()
    return str(engine.resolved_executable or executable), version


def run(
    config: Any,
    manifest: SampleManifest,
    mode: RunMode,
    *,
    phenotype_path: Path,
    variants_path: Path,
    similarity_path: Path,
    workdir: Path,
    distances_path: Optional[Path] = None,
    covariates_path: Optional[Path] = None,
    settings: Optional[AdapterSettings] = None,
    invoke: Optional[Invoker] = None,
) -> RealGwasRun:
    """Run pyseer's mixed model over the real-path inputs, in two passes.

    Args:
        config: Loaded configuration. ``gwas_real:`` supplies the settings,
            ``config.threads`` the process count, ``config.gwas`` the group
            floor and step 12a's alpha.
        manifest: The cohort. Must match the phenotype file exactly.
        mode: ``TEST`` needs an injected ``invoke`` (or the explicit
            ``PAPIPELINE_TEST_PYSEER=1`` opt-in, which launches the real
            binary against these synthetic fixtures); ``REAL`` is gated on
            ``runtime.allow_real_mode``; ``STUB`` is refused.
        phenotype_path: Binary 0/1 phenotype, ``sample``/``phenotype``.
        variants_path: The variant file named by ``settings.variant_source``.
        similarity_path: The stage-10 kinship matrix (``--similarity``).
        workdir: Where inputs, logs, the patterns file and the table go.
        distances_path: Stage-10 ``similarity.tsv``, required when lineage is on.
        covariates_path: Overrides ``settings.covariates`` when given.
        settings: Overrides the configured ``gwas_real:`` section.
        invoke: TEST-mode stand-in for pyseer's stdout/stderr contract.

    Returns:
        The two commands, the unique-pattern reduction the threshold came
        from, and the paths of every artefact written.

    Raises:
        NotImplementedError: A mode gate (STUB, REAL closed, TEST without an
            invoker or the opt-in).
        DataContractError: A configuration or input pyseer would either reject
            or - worse - silently reinterpret.
        RuntimeError: pyseer itself exited non-zero.
    """
    if mode is RunMode.STUB:
        raise NotImplementedError(
            "STUB mode produces contract headers and zero rows; it does not "
            "run pyseer, and a header-only association table would read as "
            "tools ran. Use TEST on the committed synthetic fixtures, or a "
            "REAL run with runtime.allow_real_mode open. "
            f"(mode={getattr(mode, 'value', mode)})"
        )
    if mode is RunMode.REAL:
        if not bool(config.runtime.get("allow_real_mode", False)):
            raise NotImplementedError(
                "REAL-mode pyseer is gated: runtime.allow_real_mode is false. "
                "Set it in the machine overlay, or open it for one session with "
                "PIPELINE_ALLOW_REAL_MODE=1. The committed overlays keep it shut "
                "so a REAL run cannot begin by accident. "
                f"(mode={getattr(mode, 'value', mode)})"
            )
    elif invoke is None and os.environ.get(OPT_IN_ENV) != "1":
        raise NotImplementedError(
            "TEST-mode pyseer needs an injected `invoke`: spec.md D8 defines "
            "TEST as running no real tools, and without an invoker this run "
            f"would launch the real binary. To do exactly that against these "
            f"synthetic fixtures instead, set {OPT_IN_ENV}=1. "
            f"(mode={getattr(mode, 'value', mode)})"
        )

    if settings is None:
        settings = AdapterSettings.from_config(config)
    else:
        settings.validate()

    pairs = _read_phenotype(Path(phenotype_path))
    sample_ids = [str(s) for s in manifest.sample_ids]
    if not sample_ids:
        raise DataContractError(
            "The cohort is empty, so there is no frequency floor (1/n) to "
            "derive and no outcome to associate.",
            n_samples=0,
        )
    _require_same_cohort(
        sample_ids,
        [s for s, _ in pairs],
        "phenotype file",
        extra_hint=(
            "The manifest decides cohort membership and the phenotype file "
            "states each member's outcome; a sample in one and not the other "
            "is a mismatch to stop on, never a silent drop."
        ),
    )

    n_positive = sum(v for _, v in pairs)
    n_negative = len(pairs) - n_positive
    floor = config.gwas.min_samples_per_group
    if n_positive < floor or n_negative < floor:
        raise DataContractError(
            f"Too few samples in one outcome group to run GWAS "
            f"({n_positive} positive, {n_negative} negative, minimum {floor} "
            "per group).",
            n_positive=n_positive,
            n_negative=n_negative,
            minimum=floor,
            minimum_key="gwas.min_samples_per_group",
        )

    found = _read_matrix_ids(Path(similarity_path), "similarity (kinship)")
    _require_same_cohort(
        sample_ids,
        found,
        "similarity matrix",
        extra_hint=(
            "pyseer intersects the phenotype against the similarity matrix "
            "silently (pyseer/lmm.py: p.index.intersection(K.index)), so a "
            "mismatch would analyse a different cohort and still exit 0. The "
            "matrix must describe exactly these samples."
        ),
    )

    covariates: Optional[Path]
    if covariates_path is not None:
        covariates = Path(covariates_path)
    elif settings.covariates:
        covariates = Path(settings.covariates)
    else:
        covariates = None

    distances: Optional[Path] = None
    if settings.lineage:
        if distances_path is None:
            raise DataContractError(
                "gwas_real.lineage is on, so a distances matrix is required: "
                "pyseer refuses `--lineage` without `--distances` "
                "(__main__.py: 'Must also provide a distance matrix to report "
                f"lineage effects'). Stage 10's {STAGE10_DISTANCES_FILENAME} "
                "is the matrix this path uses; pass distances_path=, or turn "
                "lineage off.",
                key="gwas_real.lineage",
            )
        distances = Path(distances_path)
        found = _read_matrix_ids(distances, "distances")
        _require_same_cohort(
            sample_ids,
            found,
            "distances matrix",
            extra_hint=(
                "pyseer exits on a lineage/similarity set mismatch "
                "(__main__.py: 'Lineage file and similarity matrix contain "
                "different sets of samples'), so the check is made here first."
            ),
        )

    if covariates is not None:
        cov_ids = _read_covariate_ids(covariates)
        missing = [s for s in sample_ids if s not in set(cov_ids)]
        if missing:
            raise DataContractError(
                f"{len(missing)} cohort sample(s) are absent from the covariates "
                f"file at {covariates}: "
                + ", ".join(missing[:10])
                + (" ..." if len(missing) > 10 else "")
                + ". pyseer requires every sample with a phenotype to be "
                "present in the covariate file and exits 1 otherwise; the "
                "check is made here so the message names the samples.",
                missing=",".join(missing[:10]),
                path=str(covariates),
            )

    variants = Path(variants_path)
    if not variants.is_file():
        raise DataContractError(
            f"Variant file not found at {variants}. The real pyseer path reads "
            f"a --{settings.variant_source} file as an input; it does not "
            "produce one.",
            path=str(variants),
            variant_source=settings.variant_source,
        )

    threads = config.threads
    if threads is None:
        raise DataContractError(
            "runtime.threads is not set, so pyseer has no process count to be "
            "given (--cpu). It is a machine fact and belongs in the machine "
            "overlay; inferring one here would be a guess about hardware "
            "nobody on this machine described.",
            key="runtime.threads",
        )

    n = len(sample_ids)
    min_af = 1.0 / n
    max_af = 1.0 - min_af
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    patterns_file = workdir / stage_gwas.PATTERNS_FILENAME

    version: Optional[str] = None
    if invoke is None:
        # Resolved before pass 1, exactly as the stage does it: a missing tool
        # is refused by name rather than as a traceback inside a log.
        _resolved, version = _require_pyseer(config, workdir)

    common = dict(
        phenotype=Path(phenotype_path),
        variant_source=settings.variant_source,
        variants=variants,
        similarity=Path(similarity_path),
        threads=threads,
        min_af=min_af,
        max_af=max_af,
        covariates=covariates,
        use_covariates=settings.use_covariates,
        distances=distances,
        lineage=settings.lineage,
        lineage_file=(
            workdir / LINEAGE_FILENAME if settings.lineage else None
        ),
        burden=(Path(settings.burden) if settings.burden else None),
    )

    # --- pass 1: record a pattern for every eligible variant ---------------
    if patterns_file.exists():
        patterns_file.unlink()
    pass1_table = workdir / "pyseer_pass1.tsv"
    LOGGER.info(
        "Step 12a (real path): first pyseer pass, to record a pattern for "
        "every eligible variant. Its association table is not read; only the "
        "patterns file is used."
    )
    pass1 = build_command(
        **common,
        output_patterns=patterns_file,
        # Both gates stated open rather than left at pyseer's defaults, so a
        # future default change cannot alter how many patterns the file holds.
        # They do not bound the file either way - pyseer writes a pattern for
        # every non-prefilter variant before evaluating its own p-value gate -
        # but stating them keeps that fact checkable.
        filter_pvalue=1.0,
        lrt_pvalue=1.0,
    )
    runner = invoke if invoke is not None else _launch
    runner(pass1, pass1_table, workdir / "pyseer_pass1.log")

    # Step 12a, reused: the vendored count_patterns.py over pyseer's own file.
    reduction = stage_gwas.reduce_unique_patterns(config, patterns_file, workdir)

    # --- pass 2: associate at the threshold 12a derived ---------------------
    association_table = workdir / "pyseer_raw.tsv"
    LOGGER.info(
        "Step 12b (real path): second pyseer pass at --lrt-pvalue %.6E "
        "(alpha %s / %d unique patterns)",
        reduction.threshold, reduction.alpha, reduction.n_unique_patterns,
    )
    pass2 = build_command(
        **common,
        lrt_pvalue=reduction.threshold,
    )
    runner(pass2, association_table, workdir / "pyseer.log")

    provenance = _write_provenance(
        workdir=workdir,
        mode=mode,
        settings=settings,
        commands=(pass1, pass2),
        reduction=reduction,
        n_samples=n,
        similarity=Path(similarity_path),
        distances=distances,
        lineage_file=common.get("lineage_file"),
        covariates=covariates,
        pyseer_version=version,
    )
    return RealGwasRun(
        mode=mode,
        settings=settings,
        commands=(tuple(pass1), tuple(pass2)),
        reduction=reduction,
        patterns_file=patterns_file,
        pass1_table=pass1_table,
        association_table=association_table,
        provenance=provenance,
        pyseer_version=version,
    )


def _write_provenance(
    *,
    workdir: Path,
    mode: RunMode,
    settings: AdapterSettings,
    commands: Sequence[Sequence[str]],
    reduction: Any,
    n_samples: int,
    similarity: Path,
    distances: Optional[Path],
    lineage_file: Optional[Path],
    covariates: Optional[Path],
    pyseer_version: Optional[str],
) -> Path:
    """Record which invocation produced the table, and at which threshold.

    Same rationale as ``PyseerEngine.write_provenance``: every p-value in the
    table comes from a python package that has already changed once relative to
    the libraries beside it, so the version and the exact argv travel with the
    artefact. Written even when the table is empty - "no association, from
    pyseer X, with these arguments" is the statement a reader needs, and a run
    that failed leaves no file at all, which is a different fact.
    """
    path = workdir / PROVENANCE_FILENAME
    rows = [
        {"key": "mode", "value": str(getattr(mode, "value", mode))},
        {"key": "variant_source", "value": settings.variant_source},
        {"key": "lineage", "value": str(settings.lineage).lower()},
        {"key": "burden", "value": settings.burden or "none"},
        {"key": "pyseer_version", "value": pyseer_version or "unresolved"},
        {
            "key": "pyseer_invoked_as",
            "value": " ".join([sys.executable, str(_PYSEER_COMPAT)]),
        },
        {"key": "kinship_source", "value": str(similarity)},
        {"key": "distances_source", "value": str(distances) if distances else "none"},
        {
            "key": "lineage_effects",
            "value": str(lineage_file) if lineage_file else "none",
        },
        {"key": "covariates", "value": str(covariates) if covariates else "none"},
        {"key": "n_samples", "value": str(n_samples)},
        {
            "key": "n_unique_patterns",
            "value": str(reduction.n_unique_patterns),
        },
        {
            "key": "lrt_pvalue_threshold",
            "value": repr(reduction.threshold),
        },
    ]
    for index, command in enumerate(commands, start=1):
        rows.append({"key": f"command_{index}", "value": " ".join(command)})
    write_tsv(
        path,
        rows,
        ("key", "value"),
        header_comment=[
            "What produced the association table on the REAL path: the tool "
            "version, the shim it ran through, the exact argv of each pass, "
            "and the unique-pattern threshold pass 2 gated on.",
            "A result whose tool version is unknown cannot be reproduced.",
        ],
    )
    return path


__all__ = [
    "AdapterSettings",
    "OPT_IN_ENV",
    "PHENOTYPE_COLUMN",
    "PROVENANCE_FILENAME",
    "RealGwasRun",
    "VARIANT_SOURCES",
    "build_command",
    "run",
]
