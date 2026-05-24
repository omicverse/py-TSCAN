# Acceleration Iteration Log — py-TSCAN

> Parsed by `engine/plot_evolution.py` to render the two-panel evolution figure.
> Schema strict — see `omicverse-rebuild/templates/ITERATION_LOG.template.md`.
>
> **Mode**: verification-mode pipeline (R clusterid → Py order) — see [MATH.md §4.1](MATH.md#41-mclust-em-divergence-upstream-py-mclustr) for why default-mode parity is gated by an upstream py-mclustR limitation.

---

## Baseline — 2026-05-24 11:00:00

```yaml
iter: 0
status: baseline
action: null
admissibility: null
playbook_section: null
wall_clock_mean_s: 0.0622
wall_clock_stddev_s: 0.0001
wall_clock_runs_s: [0.0622, 0.0623, 0.0621]
warmup_run_s: 0.0633
parity_metric: 1.000000
parity_class: ordinal
parity_threshold: 0.99
parity_passes: true
notes: |
  Equivalence-Agent clean translation. Pearson 1.000000 vs R TSCAN on lpsdata
  in verification mode. Per-phase: preprocess 32 ms, exprmclust 39 ms,
  TSCANorder 0.5 ms. TSCAN's R reference takes ~1.5 s for the same pipeline,
  so the baseline is already ~24× faster from R→Python alone.
```

---

## iter 1 — 2026-05-24 11:05:00

```yaml
iter: 1
status: ACCEPT
action: skip_uv_on_first_svd
playbook_section: "§1.5"
admissibility: exact
admissibility_evidence: |
  The first SVD inside _prcomp_sdev is only consumed for `sdev[1:20]`,
  the singular values used by the piecewise-linear elbow detector.
  np.linalg.svd(..., compute_uv=False) returns just the singular values
  with NO change in their numerical content. (E) — exact algebraic identity.
wall_clock_mean_s: 0.0528
wall_clock_stddev_s: 0.0004
wall_clock_runs_s: [0.0533, 0.0527, 0.0524]
warmup_run_s: 0.0540
speedup_vs_previous: 1.18
speedup_vs_baseline: 1.18
parity_metric: 1.000000
parity_delta_vs_baseline: 0.000000
parity_passes: true
math_reason_for_dip: null
```

### Decision

ACCEPT — kept in `pytscan/exprmclust.py::_prcomp_sdev`. Working tree's HEAD
includes `compute_uv=False`. Verified parity test `test_pseudotime_parity_verification_mode`
still passes.

### Commit / branch

```
branch: acceleration-iter-1-skip_uv
commit: (pre-release working state)
```

---

## iter 2 — 2026-05-24 11:08:00

```yaml
iter: 2
status: REJECT_SLOW
action: scipy_gesdd_svd
playbook_section: "§1.5"
admissibility: exact
admissibility_evidence: |
  scipy.linalg.svd with lapack_driver='gesdd' uses the divide-and-conquer
  LAPACK driver instead of the default gesdd-via-numpy. Numerically
  identical singular values (LAPACK guarantees bit-exact same result for
  fixed driver). Admissibility was (E) — exact algebraic identity.
wall_clock_mean_s: 0.1823
wall_clock_stddev_s: 0.0148
wall_clock_runs_s: [0.1835, 0.1636, 0.2028, 0.1882, 0.1736]
warmup_run_s: 0.1437
speedup_vs_previous: 0.29
speedup_vs_baseline: 0.34
parity_metric: 1.000000
parity_delta_vs_baseline: 0.000000
parity_passes: true
math_reason_for_dip: null
```

### Decision

REJECT_SLOW — scipy.linalg.svd carries higher wrapper overhead than
np.linalg.svd at this matrix size (131 × 16776). For a one-shot call on a
small matrix the dispatch cost dominates. Hypothesis was wrong; rolled back.

(The CV here is 8.1% — close to the 10% threshold; benchmark.py auto-extended
to 5 runs.)

### Commit / branch

```
branch: acceleration-iter-2-scipy_gesdd  (deleted after rollback)
```

---

## Stop reason

Playbook exhausted for this port's pattern. TSCAN is too small (131 cells,
G ≤ 9 clusters, no inner-loop algebraic structure) to benefit from the
heavier rewrites:

- §1.1 cache X^T X — N/A, no iterative inner loop with repeated `X^T X` reuse
- §1.2 Woodbury K×K Cholesky — N/A, no `(I + λL)` solve in TSCAN
- §2.1 sparse soft-assignment — N/A, no responsibility matrix in our control
- §3.1 MST ⊆ Delaunay — N/A at G=3 (trivially fast as dense G×G)

Only `§1.5 skip U/V` from the Acceleration Playbook applied. Final speedup
vs baseline: **1.18×**, vs R reference: **~28×**.

The port is classified **Class A — translation-only** with one minor
(E) optimization. See [`RECONSTRUCTION_REPORT.md`](RECONSTRUCTION_REPORT.md)
§4 for the full audit.

---

## Summary

| iter | action | admissibility | mean time (s) | speedup vs baseline | accuracy | status |
|---|---|---|---|---|---|---|
| 0 | (baseline) | — | 0.0622 | 1× | 1.000000 | — |
| 1 | §1.5 skip U/V on first SVD | E | 0.0528 | 1.18× | 1.000000 | ACCEPT |
| 2 | §1.5 scipy gesdd driver | E | 0.1823 | 0.34× | 1.000000 | REJECT_SLOW |
