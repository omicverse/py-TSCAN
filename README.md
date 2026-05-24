# py-TSCAN

A **pure-Python re-implementation of [TSCAN](https://github.com/zji90/TSCAN)** (Ji & Ji, *Nucleic Acids Research* 2016) for pseudo-time reconstruction in single-cell RNA-seq.

- AnnData-native — drop-in for the scanpy ecosystem
- **No `rpy2`**, no R install — implemented directly in NumPy / SciPy / pymclustR
- Same function surface as the R workflow (`preprocess` → `exprmclust` → `TSCANorder` → `difftest` → `orderscore`)
- ~28× faster than R TSCAN on the canonical `lpsdata` fixture
- **Pearson = 1.000000** vs R pseudotime in verification mode (when Mclust clusters match)

> This is a **standalone mirror** of the canonical implementation that lives in [`omicverse`](https://github.com/Starlitnightly/omicverse) (planned: `omicverse.single.TSCAN`).
> See [`RECONSTRUCTION_REPORT.md`](RECONSTRUCTION_REPORT.md) for the full audit, including a documented upstream py-mclustR limitation that affects default-mode end-to-end clustering parity.

## Install

```bash
pip install pytscan
```

## Quick-start (class API)

```python
import anndata as ad
from pytscan import TSCAN

adata = ad.read_h5ad("mydata.h5ad")     # cells × genes, raw counts in .X

tsc = TSCAN(adata)
tsc.preprocess().exprmclust().order()

adata.obs["Pseudotime"]                  # per-cell pseudotime
adata.obs["State"]                       # cluster state per cell
adata.obs["TSCAN_clusterid"]             # raw Mclust cluster id
```

## Low-level functional API (mirrors R one-to-one)

```python
from pytscan import preprocess, exprmclust, TSCANorder, difftest, orderscore

procdata = preprocess(data)                      # genes × cells DataFrame
mc       = exprmclust(procdata)                  # PCA + Mclust + MST
ordering = TSCANorder(mc, orderonly=False)       # pseudotime DataFrame
de       = difftest(procdata, ordering.ordered_names, df=3)
pos      = orderscore(subpopulation, [ordering.ordered_names])
```

## What's included

| Python | R counterpart | Purpose |
|---|---|---|
| `TSCAN` class | — | AnnData-native lifecycle wrapper |
| `preprocess` | `preprocess` | log + filter low-expression / low-CV genes |
| `exprmclust` | `exprmclust` | PCA + adaptive PC dim + Mclust + MST on centers |
| `TSCANorder` | `TSCANorder` | longest-path search + edge projection → ordering |
| `difftest` | `difftest` | GAM likelihood-ratio test for trajectory-DE |
| `orderscore` | `orderscore` | pseudotemporal ordering score (POS) |

Six R functions are deliberately skipped (plotting, Shiny GUI, supervised variants); see [`RECONSTRUCTION_REPORT.md §2.3`](RECONSTRUCTION_REPORT.md).

## Reproducing R results

```python
from pytscan import preprocess, exprmclust, TSCANorder
import pyreadr, pandas as pd, numpy as np

lpsdata = next(iter(pyreadr.read_r("data/fixture_lpsdata.rds").values()))
procdata = preprocess(lpsdata)
mc = exprmclust(procdata)
ord_res = TSCANorder(mc, orderonly=False)

# Pearson against R TSCAN pseudotime in verification mode is 1.000000
# (see tests/test_exact_match.py for the gate).
```

`tests/test_exact_match.py` runs the R reference under the `CMAP` conda env and asserts the pre-registered parity gate.

## Known limitations

This is a **Class A — translation-only** port (per the [PolyPort recipe](https://github.com/omicverse/omicverse-rebuild)). The pipeline's algorithmic steps (PCA / elbow detection / cluster centers / MST / edge-projection ordering) are **bit-identical** to R TSCAN — verified by feeding R's clusterid back in (verification mode → Pearson = 1.000000).

**Default mode** (full Python pipeline including Mclust) currently diverges from R at the Mclust step due to an upstream `py-mclustR` limitation: its EM lacks R's anti-degeneracy regularization and converges to a degenerate cluster on `lpsdata`. Details and path forward in [`MATH.md §4`](MATH.md) and [`RECONSTRUCTION_REPORT.md §6.1`](RECONSTRUCTION_REPORT.md).

## Relationship to omicverse

Developed **upstream** in [`omicverse`](https://github.com/Starlitnightly/omicverse). Standalone mirror (this repo): same code, same API, minus the omicverse packaging.

## Citation

If you use this package, please cite the original TSCAN paper:

> Ji, Z. & Ji, H. **TSCAN: pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.** *Nucleic Acids Research* 44(13), e117 (2016).

and acknowledge omicverse / this repo for the Python port.

## License

GNU GPL v3 — compatible with TSCAN's upstream `GPL (>= 2)` license.
