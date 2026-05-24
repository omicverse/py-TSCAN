"""exprmclust — PCA + adaptive PC selection + Mclust + MST on cluster centers.

Mirrors TSCAN::exprmclust (TSCAN/R/exprmclust.R).

R reference:
    exprmclust(data, clustermethod='mclust', clusternum=2:9,
               modelNames='VVV', reduce=TRUE, cluster=NULL)

    set.seed(12345)
    if (reduce) {
      sdev <- prcomp(t(data), scale=TRUE)$sdev[1:20]
      x <- 1:20
      optpoint <- which.min(sapply(2:10, function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(sdev ~ x + x2)$residuals^2)
      }))
      pcadim <- optpoint + 1
      tmpdata <- t(apply(data, 1, scale))        # row-wise z-score per gene
      tmppc   <- prcomp(t(tmpdata), scale=TRUE)
      pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
    }
    Mclust(pcareduceres, G=clusternum, modelNames=modelNames)
    clucenter[cid,] <- colMeans(pcareduceres[clusterid==cid,])
    dp <- as.matrix(dist(clucenter))
    minimum.spanning.tree(graph.adjacency(dp, weighted=TRUE))

We reproduce this end-to-end:
- piecewise-linear elbow on prcomp(t(data), scale=T)$sdev[1:20]
- two-stage standardisation: per-gene z-score then prcomp(scale=T)
- Mclust with modelNames="VVV", G=2:9, internal seed 12345
- complete graph + scipy.sparse.csgraph.minimum_spanning_tree
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree


# ----------------------------------------------------------------------------- #
# Helpers                                                                       #
# ----------------------------------------------------------------------------- #

def _prcomp_sdev(x: np.ndarray, scale: bool = True) -> np.ndarray:
    """Mirror R prcomp(...)$sdev exactly.

    R: prcomp centers, optionally scales by sample SD (ddof=1), then SVD on
    the resulting matrix. sdev[i] = singular_value[i] / sqrt(n-1).

    acceleration: §1.5 Skip U/V computation (E, exact). Only s is needed
    for the elbow detection; np.linalg.svd(..., compute_uv=False) is faster.
    """
    x = np.asarray(x, dtype=np.float64)
    n, p = x.shape
    xc = x - x.mean(axis=0, keepdims=True)
    if scale:
        sd = xc.std(axis=0, ddof=1)
        sd_safe = np.where(sd > 0, sd, 1.0)
        xc = xc / sd_safe
    s = np.linalg.svd(xc, full_matrices=False, compute_uv=False)
    return s / np.sqrt(n - 1)


def _prcomp_rotation_scores(
    x: np.ndarray,
    pcadim: int,
    scale: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (scores, rotation) matching R's prcomp.

    scores  = (x_centered_scaled) %*% rotation
    rotation = right singular vectors (columns)
    Sign convention: match R prcomp by enforcing the largest |element| of
    each column of rotation to be POSITIVE (R does the same via the
    underlying LAPACK choices on most platforms).
    """
    x = np.asarray(x, dtype=np.float64)
    n, p = x.shape
    xc = x - x.mean(axis=0, keepdims=True)
    if scale:
        sd = xc.std(axis=0, ddof=1)
        sd_safe = np.where(sd > 0, sd, 1.0)
        xc = xc / sd_safe
    u, s, vt = np.linalg.svd(xc, full_matrices=False)
    rotation = vt.T
    # Sign: R fixes the sign of each column so the first element with max abs
    # is non-negative. We mimic that to keep pcareduceres comparable column-wise.
    for j in range(rotation.shape[1]):
        max_idx = int(np.argmax(np.abs(rotation[:, j])))
        if rotation[max_idx, j] < 0:
            rotation[:, j] = -rotation[:, j]
            u[:, j] = -u[:, j]
    scores = xc @ rotation[:, :pcadim]
    return scores, rotation


def _piecewise_elbow_pcadim(sdev20: np.ndarray) -> int:
    """Replicate R's elbow detection.

    R code:
        x <- 1:20
        optpoint <- which.min(sapply(2:10, function(i) {
            x2 <- pmax(0, x - i)
            sum(lm(sdev ~ x + x2)$residuals^2)
        }))
        pcadim <- optpoint + 1

    sapply(2:10, ...) returns a length-9 vector indexed 1..9 by which.min,
    so optpoint = (i - 2 + 1) = i - 1 when i is the winning breakpoint.
    Therefore pcadim = optpoint + 1 = i.

    Fit y = a + b*x + c*max(0, x - i) for i in 2..10; pick i with min RSS;
    return that i directly. (Bug fix: earlier draft returned i+1.)
    """
    x = np.arange(1, 21, dtype=np.float64)
    sdev = np.asarray(sdev20, dtype=np.float64)
    if sdev.shape[0] != 20:
        raise ValueError(f"sdev20 must have length 20, got {sdev.shape[0]}")
    best_rss = np.inf
    best_i = 2
    for i in range(2, 11):
        x2 = np.maximum(0.0, x - i)
        X = np.column_stack([np.ones_like(x), x, x2])
        beta, *_ = np.linalg.lstsq(X, sdev, rcond=None)
        resid = sdev - X @ beta
        rss = float(resid @ resid)
        if rss < best_rss:
            best_rss = rss
            best_i = i
    return best_i


def _zscore_rows(data: np.ndarray) -> np.ndarray:
    """t(apply(data, 1, scale)) — z-score each row using R's sample SD."""
    mu = data.mean(axis=1, keepdims=True)
    sd = data.std(axis=1, ddof=1, keepdims=True)
    sd_safe = np.where(sd > 0, sd, 1.0)
    return (data - mu) / sd_safe


def _mst_on_centers(clucenter: np.ndarray):
    """Build the MST on a small G×G dense Euclidean-distance matrix.

    Returns a networkx.Graph with G nodes and G-1 edges, weighted by distance.
    Node labels are 1..G (matching R's 1-based cluster ids).
    """
    import networkx as nx
    G = clucenter.shape[0]
    # dense pairwise
    diff = clucenter[:, None, :] - clucenter[None, :, :]
    dp = np.sqrt((diff ** 2).sum(-1))
    # MST
    sp = csr_matrix(dp)
    mst = minimum_spanning_tree(sp).toarray()
    # symmetrise; the MST as returned is upper-triangular
    mst = mst + mst.T
    g = nx.Graph()
    for i in range(G):
        g.add_node(i + 1)
    for i in range(G):
        for j in range(i + 1, G):
            if mst[i, j] > 0:
                g.add_edge(i + 1, j + 1, weight=float(mst[i, j]))
    return g


# ----------------------------------------------------------------------------- #
# Public API                                                                    #
# ----------------------------------------------------------------------------- #

@dataclass
class ExprMclustResult:
    """Mirror R's exprmclust return list."""
    pcareduceres: np.ndarray             # n_cells × pcadim
    MSTtree: object                       # networkx.Graph
    clusterid: np.ndarray                 # length n_cells, ints 1..G
    clucenter: np.ndarray                 # G × pcadim
    pcadim: int
    cell_names: pd.Index

    def to_dict(self) -> dict:
        return {
            "pcareduceres": self.pcareduceres,
            "MSTtree": self.MSTtree,
            "clusterid": self.clusterid,
            "clucenter": self.clucenter,
            "pcadim": self.pcadim,
            "cell_names": list(self.cell_names),
        }


def exprmclust(
    data,
    clustermethod: str = "mclust",
    clusternum: Iterable[int] | int = range(2, 10),
    modelNames: str = "VVV",
    reduce: bool = True,
    cluster: np.ndarray | None = None,
    seed: int = 12345,
) -> ExprMclustResult:
    """Pure-Python equivalent of TSCAN::exprmclust.

    Args:
        data: genes × cells matrix (DataFrame or 2D array).
        clustermethod: 'mclust' (default) or 'kmeans'.
        clusternum: candidate G values for mclust; or fixed K for kmeans.
        modelNames: Mclust covariance parameterisation. Default 'VVV'
            (ellipsoidal, varying volume/shape/orientation).
        reduce: whether to apply PCA before clustering.
        cluster: optional pre-specified cluster labels (1-indexed ints).
        seed: matches R's `set.seed(12345)`. Default 12345.
    """
    # ----------- accept DataFrame or ndarray ---------------------------------
    if isinstance(data, pd.DataFrame):
        cell_names = data.columns
        arr = data.to_numpy(dtype=np.float64, copy=True)
    else:
        arr = np.asarray(data, dtype=np.float64)
        cell_names = pd.RangeIndex(arr.shape[1])

    rng = np.random.default_rng(seed)
    _ = rng  # currently mclust_py uses its own seed argument

    # ----------- PCA reduction ----------------------------------------------
    if reduce:
        # First pass: sdev on prcomp(t(data), scale=T)[1:20] → elbow → pcadim
        sdev = _prcomp_sdev(arr.T, scale=True)
        # R uses [1:20]; if fewer than 20 components, pad with zeros (matches
        # R behaviour for tall matrices).
        sdev20 = np.zeros(20)
        k = min(20, sdev.shape[0])
        sdev20[:k] = sdev[:k]
        pcadim = _piecewise_elbow_pcadim(sdev20)

        # Second pass: row-wise z-score then PCA on cells × genes
        tmpdata = _zscore_rows(arr)
        scores, _ = _prcomp_rotation_scores(tmpdata.T, pcadim=pcadim, scale=True)
        pcareduceres = scores  # n_cells × pcadim
    else:
        pcadim = arr.shape[0]
        pcareduceres = arr.T  # cells × genes

    # ----------- clustering --------------------------------------------------
    if clustermethod == "mclust":
        if cluster is None:
            # restrict to G > 1
            if isinstance(clusternum, range):
                Gs = [g for g in clusternum if g > 1]
            else:
                Gs = [g for g in clusternum if g > 1]
            # Import here so the package imports cleanly without mclust_py
            from mclust_py import Mclust
            res = Mclust(pcareduceres, G=Gs, model_names=[modelNames])
            clusterid = res.classification.astype(np.int64).copy()
            # mclust_py returns labels in 0..G-1 or 1..G depending on version;
            # we normalise to 1..G (R convention).
            if clusterid.min() == 0:
                clusterid = clusterid + 1
            clunum = int(clusterid.max())
        else:
            cluster = np.asarray(cluster, dtype=np.int64)
            clusterid = cluster.copy()
            clunum = int(np.unique(cluster).size)
    elif clustermethod == "kmeans":
        from sklearn.cluster import KMeans
        K = int(clusternum) if np.isscalar(clusternum) else len(list(clusternum))
        km = KMeans(n_clusters=K, random_state=seed, n_init=10).fit(pcareduceres)
        clusterid = km.labels_.astype(np.int64) + 1  # 1-indexed
        clunum = K
    else:
        raise ValueError(f"clustermethod must be 'mclust' or 'kmeans', got {clustermethod!r}")

    # ----------- cluster centres -------------------------------------------
    clucenter = np.zeros((clunum, pcareduceres.shape[1]), dtype=np.float64)
    for cid in range(1, clunum + 1):
        members = clusterid == cid
        if members.any():
            clucenter[cid - 1] = pcareduceres[members].mean(axis=0)

    # ----------- MST on centres --------------------------------------------
    mst = _mst_on_centers(clucenter)

    return ExprMclustResult(
        pcareduceres=pcareduceres,
        MSTtree=mst,
        clusterid=clusterid,
        clucenter=clucenter,
        pcadim=pcadim,
        cell_names=pd.Index(cell_names),
    )
