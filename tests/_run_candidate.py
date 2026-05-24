"""Candidate runner — invoked under conda env omicdev.

Usage: python _run_candidate.py <fixture.rds-or-csv> <output.json>

Loads the fixture (RDS via pyreadr or CSV), runs the pytscan pipeline,
emits JSON whose keys match `data/manifest.yaml::outputs[*].location_candidate`.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Make pytscan importable from the port dir and mclust_py from py-mclustR
_PORT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT_DIR))
sys.path.insert(0, "/scratch/users/steorra/analysis/omicverse_dev/py-mclustR")

from pytscan import preprocess, exprmclust, TSCANorder  # noqa: E402


def _load_fixture(path: Path) -> pd.DataFrame:
    if path.suffix == ".rds":
        try:
            import pyreadr
        except ImportError as e:
            raise RuntimeError("Reading .rds requires pyreadr — pip install pyreadr") from e
        res = pyreadr.read_r(str(path))
        # pyreadr returns OrderedDict; first value is the data
        df = next(iter(res.values()))
        return df
    if path.suffix == ".csv":
        return pd.read_csv(path, index_col=0)
    raise ValueError(f"Unknown fixture extension: {path.suffix}")


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: python _run_candidate.py <fixture> <output.json>")
        sys.exit(2)
    fixture_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    print(f"[cand] loading fixture {fixture_path}")
    lpsdata = _load_fixture(fixture_path)
    print(f"[cand] fixture dims: {lpsdata.shape}")

    t0 = time.perf_counter()
    procdata = preprocess(lpsdata)
    t_pre = time.perf_counter() - t0
    print(f"[cand] preprocess: {t_pre:.3f}s -> {procdata.shape}")

    t0 = time.perf_counter()
    mc = exprmclust(procdata)
    t_mc = time.perf_counter() - t0
    print(f"[cand] exprmclust: {t_mc:.3f}s; G={int(mc.clusterid.max())}; "
          f"pcareduceres dims {mc.pcareduceres.shape}; pcadim={mc.pcadim}")

    t0 = time.perf_counter()
    ord_res = TSCANorder(mc, orderonly=False)
    t_ord = time.perf_counter() - t0
    print(f"[cand] TSCANorder: {t_ord:.3f}s -> ordered {len(ord_res.ordered_names)} cells")

    # ---- build a per-cell vector keyed to procdata column order ------
    all_cells = list(procdata.columns)
    pseudotime_per_cell = pd.Series(index=all_cells, dtype=float)
    pseudotime_per_cell[:] = np.nan
    if ord_res.df is not None:
        pseudotime_per_cell.loc[ord_res.df["sample_name"].values] = \
            ord_res.df["Pseudotime"].values
    state_per_cell = pd.Series(index=all_cells, dtype=float)
    state_per_cell[:] = np.nan
    if ord_res.df is not None:
        state_per_cell.loc[ord_res.df["sample_name"].values] = \
            ord_res.df["State"].values

    out = {
        "cell_names": all_cells,
        "Pseudotime": pseudotime_per_cell.tolist(),
        "State": state_per_cell.astype("Int64").tolist(),
        "ordered_names": list(ord_res.ordered_names),
        "ordered_pseudotime": list(range(1, len(ord_res.ordered_names) + 1)),
        "clusterid": [int(x) for x in mc.clusterid],
        "pcadim": int(mc.pcadim),
        "pcareduceres": mc.pcareduceres.tolist(),
        "clucenter": mc.clucenter.tolist(),
        "procdata_dim": list(procdata.shape),
        "procdata_rownames": list(procdata.index),
        "timings": {
            "preprocess": t_pre,
            "exprmclust": t_mc,
            "TSCANorder": t_ord,
        },
    }
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"[cand] wrote {output_path}")


if __name__ == "__main__":
    main()
