"""preprocess — log + filter low-expression / low-CV genes.

Mirrors TSCAN::preprocess (TSCAN/R/preprocess.R).

R reference:
    preprocess(data, clusternum=NULL, takelog=TRUE, logbase=2,
               pseudocount=1, minexpr_value=1, minexpr_percent=0.5, cvcutoff=1)

The R code does:
    if (takelog) data <- log(data + pseudocount) / log(logbase)
    data <- data[rowMeans(data > minexpr_value) > minexpr_percent
                 & apply(data,1,sd)/rowMeans(data) > cvcutoff, ]
    if (!is.null(clusternum)) hierarchical clustering on rows + aggregate

Input: genes × cells (R convention).  We keep that convention.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def preprocess(
    data,
    clusternum: int | None = None,
    takelog: bool = True,
    logbase: float = 2.0,
    pseudocount: float = 1.0,
    minexpr_value: float = 1.0,
    minexpr_percent: float = 0.5,
    cvcutoff: float = 1.0,
) -> pd.DataFrame:
    """Match TSCAN::preprocess exactly.

    Args:
        data: genes × cells matrix (numpy array or pandas DataFrame; the
            DataFrame's index = gene names, columns = cell names).
        clusternum: if not None, hierarchical-cluster genes and aggregate by mean.
            Default None — no clustering.
        takelog: log-transform raw counts first.
        logbase: base of logarithm.
        pseudocount: pseudocount added before log.
        minexpr_value: floor below which a cell counts as "low".
        minexpr_percent: gene must have > this fraction of cells above the floor.
        cvcutoff: SD/mean coefficient-of-variation floor across cells per gene.

    Returns:
        Filtered (and optionally clustered) genes × cells DataFrame.
    """
    is_df = isinstance(data, pd.DataFrame)
    if is_df:
        gene_names = data.index
        cell_names = data.columns
        arr = data.to_numpy(dtype=float, copy=True)
    else:
        arr = np.asarray(data, dtype=float).copy()
        gene_names = pd.RangeIndex(arr.shape[0])
        cell_names = pd.RangeIndex(arr.shape[1])

    if takelog:
        arr = np.log(arr + pseudocount) / np.log(logbase)

    # Filter: rowMeans(data > minexpr_value) > minexpr_percent
    high_expr_mask = (arr > minexpr_value).mean(axis=1) > minexpr_percent
    # Filter: apply(data,1,sd)/rowMeans(data) > cvcutoff
    # R's sd is sample std (n-1). numpy default is population (n); use ddof=1.
    with np.errstate(divide="ignore", invalid="ignore"):
        row_sd = arr.std(axis=1, ddof=1)
        row_mean = arr.mean(axis=1)
        cv = np.where(row_mean != 0, row_sd / row_mean, 0.0)
    cv_mask = cv > cvcutoff

    keep = high_expr_mask & cv_mask
    arr = arr[keep]
    gene_names = pd.Index(gene_names)[keep]

    if clusternum is not None:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import pdist
        # R: hclust(dist(data)) defaults to complete linkage, Euclidean dist
        dmat = pdist(arr, metric="euclidean")
        Z = linkage(dmat, method="complete")
        # cutree(clures, clusternum) → fcluster(Z, t=clusternum, criterion='maxclust')
        cluster_ids = fcluster(Z, t=clusternum, criterion="maxclust")
        df = pd.DataFrame(arr, columns=cell_names)
        df["__cluster__"] = cluster_ids
        agg = df.groupby("__cluster__").mean()
        arr = agg.to_numpy()
        gene_names = pd.Index([f"cluster_{i+1}" for i in range(arr.shape[0])])

    return pd.DataFrame(arr, index=gene_names, columns=cell_names)
