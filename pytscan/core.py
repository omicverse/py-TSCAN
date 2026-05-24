"""TSCAN — AnnData-native class wrapper around the functional API.

Provides:
    >>> tsc = TSCAN(adata)
    >>> tsc.preprocess().exprmclust().order()
    >>> adata.obs['Pseudotime'], adata.obs['State']

The functional API in `preprocess`, `exprmclust`, `ordering`, `difftest`,
`orderscore` mirrors R one-to-one for users who want the original R-style
workflow.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .exprmclust import ExprMclustResult, exprmclust as _exprmclust
from .ordering import TSCANorder as _TSCANorder
from .preprocess import preprocess as _preprocess


class TSCAN:
    """AnnData-native TSCAN trajectory wrapper."""

    def __init__(self, adata, layer: str | None = None):
        """Args:
            adata: anndata.AnnData (cells × genes; counts in .X or `layer`).
            layer: which layer to use (default .X).
        """
        try:
            import anndata as ad  # noqa: F401
        except ImportError:  # pragma: no cover
            raise ImportError("pytscan.TSCAN requires `anndata`; install with `pip install anndata`")
        self.adata = adata
        self._layer = layer
        self._procdata: pd.DataFrame | None = None
        self._exprmclust: ExprMclustResult | None = None
        self._order = None

    # ---------------- chained API ---------------------------------------

    def preprocess(self, **kwargs) -> "TSCAN":
        # AnnData is cells × genes; R convention is genes × cells.
        X = self.adata.X if self._layer is None else self.adata.layers[self._layer]
        if hasattr(X, "toarray"):
            X = X.toarray()
        gene_names = self.adata.var_names
        cell_names = self.adata.obs_names
        df = pd.DataFrame(np.asarray(X).T, index=gene_names, columns=cell_names)
        self._procdata = _preprocess(df, **kwargs)
        return self

    def exprmclust(self, **kwargs) -> "TSCAN":
        if self._procdata is None:
            raise RuntimeError("Call .preprocess() before .exprmclust()")
        self._exprmclust = _exprmclust(self._procdata, **kwargs)
        # write cluster ID back to adata.obs (R order)
        ser = pd.Series(
            self._exprmclust.clusterid,
            index=self._exprmclust.cell_names,
            name="TSCAN_clusterid",
        )
        self.adata.obs["TSCAN_clusterid"] = ser.reindex(self.adata.obs_names)
        return self

    def order(self, **kwargs) -> "TSCAN":
        if self._exprmclust is None:
            raise RuntimeError("Call .exprmclust() before .order()")
        res = _TSCANorder(self._exprmclust, orderonly=False, **kwargs)
        self._order = res
        # write pseudotime + state back to adata.obs
        if res.df is not None:
            ser_p = pd.Series(
                res.df["Pseudotime"].values,
                index=res.df["sample_name"].values,
                name="Pseudotime",
            )
            ser_s = pd.Series(
                res.df["State"].values,
                index=res.df["sample_name"].values,
                name="State",
            )
            self.adata.obs["Pseudotime"] = ser_p.reindex(self.adata.obs_names)
            self.adata.obs["State"] = ser_s.reindex(self.adata.obs_names)
        return self

    # ---------------- convenience accessors ------------------------------

    @property
    def procdata(self) -> pd.DataFrame:
        if self._procdata is None:
            raise RuntimeError("Call .preprocess() first")
        return self._procdata

    @property
    def mclustobj(self) -> ExprMclustResult:
        if self._exprmclust is None:
            raise RuntimeError("Call .exprmclust() first")
        return self._exprmclust

    @property
    def ordering(self) -> Any:
        if self._order is None:
            raise RuntimeError("Call .order() first")
        return self._order
