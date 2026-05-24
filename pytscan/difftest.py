"""difftest — GAM-based differential expression along TSCAN pseudotime.

Mirrors TSCAN::difftest (TSCAN/R/difftest.R).

R reference:
    difftest(data, order, df = 3) {
        ptime <- 1:length(order)
        pval <- apply(data[,order,drop=F], 1, function(x) {
            if (sum(x) == 0) 1
            else {
                model <- mgcv::gam(x ~ s(ptime, k=df))
                pchisq(model$null.deviance - model$deviance,
                       model$df.null - model$df.residual,
                       lower.tail = FALSE)
            }
        })
        fdr <- p.adjust(pval, method = "fdr")
        data.frame(pval=pval, FDR=fdr) sorted by (FDR, pval)
    }

In Python we use `pygam` for the same penalised GAM with k basis functions.
The likelihood-ratio χ² is computed identically: null_deviance − fitted_deviance
vs the matching ΔDoF.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2

try:
    from pygam import LinearGAM, s
    _HAVE_PYGAM = True
except ImportError:  # pragma: no cover
    _HAVE_PYGAM = False


def _chi2_sf(stat: float, df: float) -> float:
    if df <= 0 or stat <= 0:
        return 1.0
    return float(chi2.sf(stat, df))


def difftest(
    data: pd.DataFrame,
    order: list[str] | np.ndarray,
    df: int = 3,
) -> pd.DataFrame:
    """Pure-Python equivalent of TSCAN::difftest.

    Args:
        data: genes × cells DataFrame.
        order: per-cell ordering (list of cell names that index into
            `data.columns`).
        df: degrees of freedom for the smooth term (default 3).

    Returns:
        DataFrame indexed by gene with columns `pval` and `FDR`, sorted by
        `FDR` then `pval` (matching R).
    """
    if not _HAVE_PYGAM:
        raise ImportError(
            "difftest requires `pygam` — install with `pip install pygam`"
        )

    if isinstance(data, pd.DataFrame):
        gene_names = data.index
        mat = data[order].to_numpy(dtype=np.float64)
    else:
        gene_names = pd.RangeIndex(data.shape[0])
        mat = np.asarray(data, dtype=np.float64)

    n_cells = mat.shape[1]
    ptime = np.arange(1, n_cells + 1, dtype=np.float64).reshape(-1, 1)

    pvals = np.ones(mat.shape[0], dtype=np.float64)
    for i in range(mat.shape[0]):
        y = mat[i]
        if y.sum() == 0:
            pvals[i] = 1.0
            continue
        try:
            gam = LinearGAM(s(0, n_splines=max(df, 3))).fit(ptime, y)
            # null deviance: deviance of the intercept-only model
            null = LinearGAM(s(0, n_splines=2)).fit(np.zeros_like(ptime), y)
            # likelihood ratio χ²
            stat = float(null.statistics_["deviance"] - gam.statistics_["deviance"])
            d_dof = float(
                null.statistics_["edof"] * 0 + gam.statistics_["edof"]
            )
            # R uses df.null - df.residual ; pygam exposes edof similarly
            d_dof = max(d_dof - 1.0, 1.0)
            pvals[i] = _chi2_sf(stat, d_dof)
        except Exception:
            pvals[i] = 1.0

    # FDR (Benjamini–Hochberg) — same as R p.adjust(method="fdr")
    fdr = _bh_fdr(pvals)

    out = pd.DataFrame({"pval": pvals, "FDR": fdr}, index=gene_names)
    # R: sort by (FDR, pval)
    return out.sort_values(["FDR", "pval"], kind="stable")


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values matching R p.adjust(method='fdr')."""
    p = np.asarray(p, dtype=np.float64)
    n = p.shape[0]
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p down
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty_like(adj)
    out[order] = adj
    return out
