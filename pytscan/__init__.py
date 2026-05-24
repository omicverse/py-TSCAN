"""pytscan — pure-Python port of TSCAN (Ji & Ji 2016).

Public API (mirrors R one-to-one for the functional path, plus an
AnnData-native class):

    Class:
        TSCAN(adata)            # .preprocess().exprmclust().order()

    Functional (R-mirror):
        preprocess(data, ...)
        exprmclust(data, ...)
        TSCANorder(mclustobj, ...)
        difftest(data, order, df=3)
        orderscore(subpopulation, orders)

The internal state-machine + numerical results match TSCAN 2.0.0 on the
canonical `lpsdata` fixture under `data/manifest.yaml`'s pre-registered
parity gate (Pearson ≥ 0.99 on pseudotime; ARI ≥ 0.95 on cluster IDs).
"""

from __future__ import annotations

__version__ = "0.1.0"

from .core import TSCAN
from .exprmclust import ExprMclustResult, exprmclust
from .ordering import TSCANorder, TSCANorderResult
from .orderscore import orderscore
from .preprocess import preprocess
from .plotting import plotmclust, plotPseudotime

# Lazy imports — difftest pulls pygam which is optional
def difftest(*args, **kwargs):
    """See pytscan.difftest.difftest. Imported lazily (requires pygam)."""
    from .difftest import difftest as _difftest
    return _difftest(*args, **kwargs)


__all__ = [
    "TSCAN",
    "ExprMclustResult",
    "TSCANorder",
    "TSCANorderResult",
    "exprmclust",
    "orderscore",
    "preprocess",
    "difftest",
    "plotmclust",
    "plotPseudotime",
    "__version__",
]
