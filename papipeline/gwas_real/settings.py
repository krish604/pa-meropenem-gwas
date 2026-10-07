"""The ``gwas_real:`` settings section: what the REAL pyseer path is given.

These are *scientific* choices - which variant file pyseer is pointed at,
whether lineage enters the model, whether the run is a grouped burden test -
so they live in ``config/science.yaml`` like every other scientific choice, and
are read from ``config.raw`` here rather than parsed into the loader: the
loader is not a path this work package owns, and a section read through
``raw`` is how ``stages/gwas.py::build_input`` reads ``project`` today.

Nothing in this module is environmental. Paths a caller passes explicitly win
over the configured ones; threads and memory come from the machine overlay via
``config.threads``; the phenotype, variant and kinship files are inputs, not
configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Tuple

from ..errors import DataContractError

#: pyseer's variant group, mutually exclusive and required
#: (``pyseer/__main__.py``: ``variant_group = variants.add_mutually_exclusive_group(required=True)``).
#: Exactly one of these flags reaches the command.
VARIANT_SOURCES: Tuple[str, ...] = ("kmers", "pres", "vcf")

#: The section's name in ``config/science.yaml``.
SECTION = "gwas_real"


def _tokens(value: Any) -> Tuple[str, ...]:
    """``--use-covariates`` arguments, flattened to one tuple of tokens.

    pyseer declares the flag with ``nargs='*'``, so each whitespace-separated
    token is one argument - the help text's own example is the string
    ``"2 3q"`` which a shell splits into two. Accepting both a list and a
    single string keeps a YAML value from silently becoming one token.
    """
    if value is None:
        return ()
    if isinstance(value, str):
        value = [value]
    tokens: list = []
    for entry in value:
        tokens.extend(str(entry).split())
    return tuple(tokens)


@dataclass(frozen=True)
class AdapterSettings:
    """What the real-path command is built from.

    Attributes:
        variant_source: ``kmers``, ``pres`` or ``vcf`` - the one variant flag
            pyseer accepts.
        lineage: Report lineage effects. Requires a ``--distances`` matrix;
            stage 10's ``similarity.tsv`` is the one this path uses, because
            ``pyseer/input.py::load_structure`` reads exactly the shape stage
            10 writes.
        burden: Path to the VCF regions file for grouped burden testing, or
            ``None`` for a single-variant scan. pyseer declares
            ``--burden BURDEN`` as a *value* (the regions file), not a switch,
            and refuses it without ``--vcf``.
        covariates: Path to pyseer's headerless covariates file, or ``None``.
        use_covariates: Which columns of that file enter the model, as pyseer
            encodes them (``2`` categorical, ``2q`` quantitative).
    """

    variant_source: str = "pres"
    lineage: bool = True
    burden: Optional[str] = None
    covariates: Optional[str] = None
    use_covariates: Tuple[str, ...] = field(default=())

    @classmethod
    def from_config(cls, config: Any) -> "AdapterSettings":
        """Read the ``gwas_real:`` section, falling back to the defaults above.

        An absent section is not an error - the defaults are the section's own
        committed values, stated once - but a *present* value this adapter
        cannot honour is: a variant source pyseer has no flag for, or a burden
        file without ``--vcf``, would reach the tool as an exit code in a log
        nobody is watching.
        """
        raw: Mapping[str, Any] = (getattr(config, "raw", None) or {}).get(SECTION) or {}
        settings = cls(
            variant_source=str(raw.get("variant_source", "pres")).strip().lower(),
            lineage=bool(raw.get("lineage", True)),
            burden=(str(raw["burden"]) if raw.get("burden") else None),
            covariates=(str(raw["covariates"]) if raw.get("covariates") else None),
            use_covariates=_tokens(raw.get("use_covariates")),
        )
        settings.validate()
        return settings

    def validate(self) -> None:
        """Refuse a combination pyseer would refuse, before anything launches."""
        if self.variant_source not in VARIANT_SOURCES:
            raise DataContractError(
                f"gwas_real.variant_source is {self.variant_source!r}; pyseer "
                f"accepts exactly one of {', '.join(VARIANT_SOURCES)} "
                "(its variant group is mutually exclusive and required). An "
                "unknown name would reach the tool as an argparse error in a "
                "log nobody is watching.",
                key=f"{SECTION}.variant_source",
                value=self.variant_source,
                supported=",".join(VARIANT_SOURCES),
            )
        if self.burden and self.variant_source != "vcf":
            raise DataContractError(
                f"A burden test (gwas_real.burden = {self.burden!r}) requires "
                "gwas_real.variant_source: vcf, because pyseer refuses "
                "`--burden` without `--vcf` (__main__.py: 'burden testing "
                "(requires --vcf)'). Grouping variants by region is only "
                "defined over a VCF's sites.",
                key=f"{SECTION}.burden",
                variant_source=self.variant_source,
            )
        if bool(self.covariates) != bool(self.use_covariates):
            raise DataContractError(
                "gwas_real.covariates and gwas_real.use_covariates must be "
                "given together: pyseer loads a covariates file but uses no "
                "column of it unless --use-covariates names one, so either "
                "half on its own is a flag that changes nothing.",
                key=f"{SECTION}.covariates",
                covariates=self.covariates,
                use_covariates=" ".join(self.use_covariates),
            )


__all__ = ["AdapterSettings", "SECTION", "VARIANT_SOURCES"]
