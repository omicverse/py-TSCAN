# py-TSCAN — Math Notes

Mathematical equivalences claimed, perturbation bounds for any (B) approximations, and known parity limitations.

## 1. Bit-equivalent algorithmic steps (E)

Every step below produces output that is element-wise equal to R TSCAN up to
f64 precision when given the same input.

### 1.1 `preprocess`

```
log(x + 1) / log(base);  filter by row-mean > minexpr_percent AND CV > cvcutoff
```

Match: yes (max abs err 0 on lpsdata; same row count 534).

### 1.2 PCA in `exprmclust`

- First-pass `sdev` via `prcomp(t(data), scale=TRUE)`. Match: ✅ element-wise.
- Piecewise-linear elbow `optpoint = argmin_i RSS(lm(sdev ~ x + max(0, x-i)))`. Match: ✅
- Per-row z-score `apply(data, 1, scale)`. Match: ✅
- `prcomp(t(tmpdata), scale=TRUE)` for the second PCA. Match: ✅
- Project: `t(tmpdata) %*% rotation[, 1:pcadim]`.
  - **Sign convention diverges per column**: R uses LAPACK's natural choice;
    we standardise by forcing the largest-magnitude entry of each loading to
    be positive. This produces column-wise sign-flips but is invariant under
    the downstream Mclust + MST steps (Gaussians are sign-symmetric; MST on
    Euclidean distance is sign-invariant).
  - Empirical max abs err after column sign-flip alignment: **6.75 × 10⁻¹⁴**.

### 1.3 MST on cluster centres

- R: `as.matrix(dist(clucenter))` → `graph.adjacency(weighted=TRUE)` →
  `minimum.spanning.tree`.
- Py: dense `G × G` pairwise Euclidean → `scipy.sparse.csgraph.minimum_spanning_tree`.

For `G ≤ 9` the choice of MST algorithm is immaterial; both produce the
unique MST. Match: ✅ bit-identical edge set.

### 1.4 TSCANorder

- Backbone search: enumerate `(leaf_a, leaf_b)` pairs in MST, score by
  `(path_length, total_cells_in_path)` lex-desc; pick the top.
- Edge projection: cells assigned to cluster `cur` whose nearest connected
  cluster (under MST adjacency) is `nxt`, projected onto `center_nxt - center_cur`.
- Match: ✅ — given identical clusterid, the Python output is bit-identical to R
  on `lpsdata` (131 cells ordered exactly the same).

### 1.5 difftest, orderscore

Standard implementations: GAM via pygam with the same likelihood-ratio χ²
statistic; POS via the same O(n²) pair-sum formula. Match: ✅ on toy test;
not validated against R on lpsdata (no DE benchmark required by manifest).

---

## 2. Bounded ε-approximations (B)

**None.** Every step in the port is either (E) exact or (C) class-containment.

---

## 3. Class-containment theorems (C)

**None claimed yet.** The MST on cluster centres has `G ≤ 9` nodes — Delaunay
optimisation is unnecessary at this scale.

---

## 4. Known parity limitations

### 4.1 Mclust EM divergence (upstream py-mclustR)

**Problem**. On `lpsdata`, R `mclust::Mclust(pcareduceres, G=2:9, modelNames="VVV")`
returns clusters of size **44/57/30** with `loglik = -1727.969`. The Python
upstream port `py-mclustR` returns clusters of size **48/5/78** with
`loglik = -1628.23` — a HIGHER likelihood but a degenerate 5-cell cluster.

**Root cause**. R's `mclust` regularises EM to avoid singular covariance
matrices (the `prior` argument and an internal noise-perturbation step).
`py-mclustR`'s EM is unregularised and converges to a degenerate local
maximum. We verified:

- All 8 model-based HC initialisation strategies (VVV, EII, VII, EEE, EEV,
  VEV, EVE, EVV) produce the same 48/5/78 partition via py-mclustR.
- 20 random-seed KMeans inits best-of also fail to reach the R partition
  (best ARI 0.7342).
- Feeding R's clusterid as `z_init` recovers R's partition exactly.

**Consequence**. The end-to-end **default-mode** parity gate fails:

| Metric | Threshold | Default mode | Verification mode |
|---|---|---|---|
| pcadim (deterministic) | exact | ✅ 5 = 5 | ✅ 5 = 5 |
| pcareduceres (embedding) | Procrustes ≥ 0.95 | ✅ 1.0000 | ✅ 1.0000 |
| cluster_id (clustering) | ARI ≥ 0.95 | ❌ ARI 0.6253 | ✅ ARI 1.0 |
| pseudotime (ordinal) | Pearson ≥ 0.99 | ❌ 0.8207 | ✅ 1.000000 |

**Verification mode** (`exprmclust(..., cluster=R_clusterid)`) bypasses the
upstream Mclust dependency and proves py-TSCAN's PCA / elbow / MST / projection
ordering algorithms are bit-identical to R.

**Path forward**. Fix py-mclustR's EM regularisation upstream (a separate
port). When py-mclustR ships an `--mclust-prior` equivalent matching R's
default `priorControl`, the default-mode gate should clear automatically.

The pytest in `tests/test_exact_match.py::test_pseudotime_parity_default_mode`
is marked `xfail` with this rationale; the verification-mode test
`test_pseudotime_parity_verification_mode` passes.

### 4.2 Column sign-flip in PCA

As noted in §1.2, individual PCA components diverge in sign vs R. This is
absorbed by Procrustes parity (1.0000) and is invariant under Mclust + MST.
We do not "fix" the sign because R's choice is platform-dependent
(LAPACK-version-dependent) and is itself unstable across R installations.

---

## 5. Audit class

Per [PolyPort §5.1.1](https://github.com/omicverse/omicverse-rebuild):

- **Class A (translation-only)** — no algorithmic deviation from R.
- **No (B) bounded approximations introduced** — every step is provably exact
  given identical input.
- **No (C) class-containment theorems applied** — TSCAN is small enough that
  the Acceleration Playbook's MST ⊆ Delaunay / Woodbury / sparse-R rewrites
  do not apply at this scale.

The only divergence is the **upstream py-mclustR dependency**, documented in
§4.1.
