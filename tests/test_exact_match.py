"""The pre-registered parity gate against R TSCAN on canonical lpsdata fixture.

Uses two modes:
  - default mode: run pytscan end-to-end and compare to R
  - verification mode: feed R's clusterid into pytscan, isolate the
    ordering algorithm (proves bit-identity given identical clusters)
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

PORT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT_DIR))
sys.path.insert(0, "/scratch/users/steorra/analysis/omicverse_dev/py-mclustR")
sys.path.insert(0, "/scratch/users/steorra/analysis/omicverse_traj_dev/omicverse-rebuild/engine")

from parity_metrics import compute_parity, is_pass  # noqa: E402

import pytscan  # noqa: E402


@pytest.fixture(scope="session")
def manifest():
    with open(PORT_DIR / "data" / "manifest.yaml") as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def reference_output():
    """Run / cache R reference."""
    cache = PORT_DIR / "data" / "reference_output.json"
    if not cache.exists():
        subprocess.run(
            ["conda", "run", "-p", "/scratch/users/steorra/env/CMAP",
             "Rscript", str(PORT_DIR / "tests" / "r_reference_driver.R"),
             str(PORT_DIR / "data" / "fixture_lpsdata.rds"),
             str(cache)],
            check=True, cwd=PORT_DIR,
        )
    return json.loads(cache.read_text())


@pytest.fixture(scope="session")
def fixture_data():
    import pyreadr
    return next(iter(pyreadr.read_r(str(PORT_DIR / "data" / "fixture_lpsdata.rds")).values()))


def test_pcadim_exact_match(reference_output, fixture_data):
    """pcadim (deterministic): must match R exactly."""
    proc = pytscan.preprocess(fixture_data)
    mc = pytscan.exprmclust(proc)
    assert mc.pcadim == int(reference_output["pcadim"]), (
        f"pcadim {mc.pcadim} != R {reference_output['pcadim']}"
    )


def test_pcareduceres_bit_equivalent(reference_output, fixture_data):
    """PCA scores must match R up to column sign-flip (max abs err < 1e-10)."""
    proc = pytscan.preprocess(fixture_data)
    mc = pytscan.exprmclust(proc)
    r_pca = np.array(reference_output["pcareduceres"])
    p_pca = mc.pcareduceres
    assert r_pca.shape == p_pca.shape
    # sign-flip per column to align
    for j in range(r_pca.shape[1]):
        if np.corrcoef(r_pca[:, j], p_pca[:, j])[0, 1] < 0:
            p_pca[:, j] = -p_pca[:, j]
    err = float(np.abs(r_pca - p_pca).max())
    assert err < 1e-10, f"max abs err {err:.4e} ≥ 1e-10 — PCA diverged"


def test_pseudotime_parity_default_mode(manifest, reference_output, fixture_data):
    """Default mode: full Python pipeline vs R end-to-end."""
    proc = pytscan.preprocess(fixture_data)
    mc = pytscan.exprmclust(proc)
    ord_res = pytscan.TSCANorder(mc, orderonly=False)

    pt_py = pd.Series(np.nan, index=proc.columns)
    pt_py.loc[ord_res.df["sample_name"].values] = ord_res.df["Pseudotime"].values
    pt_r = pd.Series(reference_output["Pseudotime"],
                     index=reference_output["cell_names"]).reindex(proc.columns).values

    mask = np.isfinite(pt_r) & np.isfinite(pt_py.values)
    metric = compute_parity(pt_r[mask], pt_py.values[mask], "ordinal")
    assert is_pass(metric, "ordinal", manifest["parity_threshold"]), (
        f"Default-mode pseudotime Pearson {metric:.4f} < {manifest['parity_threshold']}"
    )


def test_pseudotime_parity_verification_mode(manifest, reference_output, fixture_data):
    """Verification mode: feed R's clusterid → prove ordering is bit-identical."""
    proc = pytscan.preprocess(fixture_data)
    cell_names = reference_output["cell_names"]
    r_clu_aligned = pd.Series(
        reference_output["clusterid"], index=cell_names
    ).reindex(proc.columns).values

    mc = pytscan.exprmclust(proc, cluster=r_clu_aligned)
    ord_res = pytscan.TSCANorder(mc, orderonly=False)

    pt_py = pd.Series(np.nan, index=proc.columns)
    pt_py.loc[ord_res.df["sample_name"].values] = ord_res.df["Pseudotime"].values
    pt_r = pd.Series(reference_output["Pseudotime"],
                     index=cell_names).reindex(proc.columns).values

    mask = np.isfinite(pt_r) & np.isfinite(pt_py.values)
    metric = compute_parity(pt_r[mask], pt_py.values[mask], "ordinal")
    assert metric >= 0.99, (
        f"Verification-mode pseudotime Pearson {metric:.6f} < 0.99 — "
        "ordering algorithm has a bug INDEPENDENT of the Mclust dependency"
    )
    # also assert the ordering is exactly identical
    ord_r = list(reference_output["ordered_names"])
    ord_py = list(ord_res.ordered_names)
    assert ord_r == ord_py, "Verification mode: ordered cell list diverges from R"
