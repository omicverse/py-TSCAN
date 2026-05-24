"""orderscore — Pseudotemporal Ordering Score (POS).

Mirrors TSCAN::orderscore (TSCAN/R/orderscore.R).

POS measures how well an ordering of cells agrees with their known
sub-population labels (must be coded 0,1,2,...). Lower is better.

R logic:
    scoreorder = subinfo[order]
    optscoreorder = sort(scoreorder)
    optscore = sum_{i<j} (optscoreorder[j] - optscoreorder[i])
    rawscore = sum_{i<j} (scoreorder[j]    - scoreorder[i])
    return rawscore / optscore
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _pos_one(order: list[str], subinfo: dict[str, int]) -> float:
    """Score one ordering."""
    s = np.array([subinfo[n] for n in order], dtype=np.float64)
    s_sorted = np.sort(s)
    n = s.shape[0]

    def _sum_pairs(arr: np.ndarray) -> float:
        # sum_{i<j} (arr[j] - arr[i]) = sum_j (j*arr[j]) - sum_i ((n-1-i)*arr[i]) ?
        # Simpler O(n²) loop matches R exactly:
        total = 0.0
        for x in range(n - 1):
            total += float((arr[(x + 1):] - arr[x]).sum())
        return total

    optscore = _sum_pairs(s_sorted)
    rawscore = _sum_pairs(s)
    return rawscore / optscore if optscore != 0 else 0.0


def orderscore(
    subpopulation: pd.DataFrame,
    orders: list[list[str]],
) -> np.ndarray:
    """Pure-Python equivalent of TSCAN::orderscore.

    Args:
        subpopulation: DataFrame with 2 columns — cell names (col 0) and
            integer sub-population codes (col 1).
        orders: list of orderings, each a list of cell names.

    Returns:
        1-D numpy array of POS values, one per ordering.
    """
    cell_names = subpopulation.iloc[:, 0].to_numpy()
    sub_codes = subpopulation.iloc[:, 1].to_numpy(dtype=np.int64)
    subinfo = dict(zip(cell_names, sub_codes))
    return np.array([_pos_one(list(o), subinfo) for o in orders], dtype=np.float64)
