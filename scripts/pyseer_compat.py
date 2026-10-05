#!/usr/bin/env python3
"""Run pyseer 1.1.2 against this environment's library versions.

pyseer 1.1.2 was released in 2021. The pilot100 environment has scipy 1.17,
numpy 2.4, pandas 3.0 and statsmodels 0.15, and pyseer uses three names that
no longer exist:

* ``smf.Logit(p, v)`` / ``smf.OLS(p, v)`` - it calls the *formula* API with
  raw ``(endog, exog)`` arrays. The formula wrappers parse their first
  argument as a formula string, so they raise ``ValueError: not enough
  values to unpack``. The array-taking classes are in ``statsmodels.api``.
* ``scipy.arange`` - numpy's namespace, removed from scipy long ago.
  ``pyseer.fastlmm`` uses it.

Every fix below binds a name that no longer exists to the identical function
it used to refer to. None of them change a fitted value, a test statistic or
a p-value; they only let the tool start. Results are produced by pyseer's own
model code.

Usage is identical to ``pyseer``:

    python scripts/pyseer_compat.py --phenotypes p.tsv --pres x.rtab ...
"""

from __future__ import annotations

import sys


def apply_shims() -> None:
    import numpy
    import scipy
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    # statsmodels: formula API called with arrays.
    for name in ("Logit", "OLS"):
        if not hasattr(smf, name):
            setattr(smf, name, getattr(sm, name))

    # scipy: numpy aliases removed from the scipy namespace. pyseer's
    # fastlmm uses several (scipy.arange, scipy.isreal, ...). Rather than
    # enumerate them, any name scipy can no longer resolve falls back to the
    # identically named numpy attribute.
    _scipy_getattr = scipy.__getattr__

    def _getattr(name: str):
        try:
            return _scipy_getattr(name)
        except (KeyError, AttributeError):
            if hasattr(numpy, name):
                return getattr(numpy, name)
            raise

    scipy.__getattr__ = _getattr
    if not hasattr(scipy, "linalg") and hasattr(numpy, "linalg"):
        scipy.linalg = numpy.linalg


def main() -> int:
    apply_shims()
    from pyseer.__main__ import main as pyseer_main
    return pyseer_main()


if __name__ == "__main__":
    sys.exit(main())
