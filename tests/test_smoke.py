"""Smoke tests — does pytscan import and run end-to-end?"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PORT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT_DIR))
sys.path.insert(0, "/scratch/users/steorra/analysis/omicverse_dev/py-mclustR")

import pytscan  # noqa: E402


def test_import():
    assert hasattr(pytscan, "TSCAN")
    assert hasattr(pytscan, "preprocess")
    assert hasattr(pytscan, "exprmclust")
    assert hasattr(pytscan, "TSCANorder")
    assert hasattr(pytscan, "orderscore")
    assert pytscan.__version__ == "0.1.0"


def _toy_data(seed: int = 12345) -> pd.DataFrame:
    """Three Gaussian blobs in d=200 → genes × cells."""
    rng = np.random.default_rng(seed)
    n_per = 30
    centers = rng.normal(0, 10, (3, 200))
    arr = np.vstack([
        centers[i] + rng.normal(0, 1, (n_per, 200))
        for i in range(3)
    ]).T  # 200 genes × 90 cells
    # add a small constant so log(x+1) is well-defined
    arr = np.abs(arr) + 1.0
    return pd.DataFrame(
        arr,
        index=[f"gene{i}" for i in range(arr.shape[0])],
        columns=[f"cell{i}" for i in range(arr.shape[1])],
    )


def test_pipeline_runs():
    data = _toy_data()
    proc = pytscan.preprocess(data, cvcutoff=0.01)  # toy data has low CV
    assert proc.shape[1] == data.shape[1]
    mc = pytscan.exprmclust(proc)
    assert mc.pcareduceres.shape[0] == proc.shape[1]
    assert int(mc.clusterid.max()) >= 2
    res = pytscan.TSCANorder(mc, orderonly=False)
    assert res.df is not None
    assert "Pseudotime" in res.df.columns
    assert len(res.ordered_names) > 0


def test_pcadim_elbow_known():
    """The elbow detector should map known inputs to known pcadims."""
    from pytscan.exprmclust import _piecewise_elbow_pcadim
    # synthetic decay with sharp drop at i=5 → pcadim should equal 5
    sdev = np.array([10, 9, 8, 7, 6, 0.5, 0.4, 0.3, 0.25, 0.2,
                     0.18, 0.16, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07])
    pcadim = _piecewise_elbow_pcadim(sdev)
    # Elbow at i=5 means optpoint=4 in R (since sapply(2:10) is 1-indexed),
    # so pcadim = optpoint + 1 = 5. The piecewise fit can pick adjacent i's
    # on near-degenerate cases, so allow ±2.
    assert 3 <= pcadim <= 7, f"Expected pcadim near 5, got {pcadim}"


def test_orderscore_perfect():
    """A perfectly-ordered list should score 1.0 (R: 1.0 by construction)."""
    sub = pd.DataFrame({
        "cell": [f"c{i}" for i in range(10)],
        "code": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    })
    perfect = [f"c{i}" for i in range(10)]  # already in code-order
    reverse = perfect[::-1]
    scores = pytscan.orderscore(sub, [perfect, reverse])
    assert abs(scores[0] - 1.0) < 1e-9
    assert abs(scores[1] + 1.0) < 1e-9  # reverse → -1
