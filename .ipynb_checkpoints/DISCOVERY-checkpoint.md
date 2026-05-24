# Discovery — py-TSCAN

> Retroactively generated after the port completed. The discovery step has been added to the protocol as Phase 0.5 (see [omicverse-rebuild/DISCOVERY.md](../omicverse-rebuild/DISCOVERY.md)); future ports will commit this file before any algorithmic code.

## 1. Is this package already ported?

```
$ python -m engine.discover_omicverse_deps --check TSCAN
## Discovery — `TSCAN`

**No existing omicverse port found.** Safe to start a new port.
```

**Decision**: START_PORT (no duplicate exists).

## 2. Dependency audit

```
$ python -m engine.discover_omicverse_deps --description TSCAN-ref/DESCRIPTION

## Discovery — dependencies of `TSCAN-ref`

### Imports

| R dep      | omicverse match                                            |
|------------|------------------------------------------------------------|
| ggplot2    | —                                                          |
| shiny      | —                                                          |
| plyr       | —                                                          |
| grid       | —                                                          |
| fastICA    | —                                                          |
| igraph     | —                                                          |
| combinat   | —                                                          |
| mgcv       | —                                                          |
| mclust     | [`py-mclustR`](https://github.com/omicverse/py-mclustR)    |
| gplots     | —                                                          |

Action items: 1 of 10 R Imports have an omicverse-org Python mirror.
```

## 3. Decisions per R dep

| R dep | omicverse match | Decision | Python replacement | Note |
|---|---|---|---|---|
| `mclust` | **`py-mclustR`** | **hard dep** | `pymclustR>=0.1` | reuse — saves ~3000 LOC. (Known parity caveat — see [MATH.md §4.1](MATH.md#41-mclust-em-divergence-upstream-py-mclustr).) |
| `mgcv` | — | native-python | `pygam>=0.9` | mature equivalent for GAM LRT |
| `igraph` | — | native-python | `networkx>=3.0` | sufficient for small MSTs at G≤9 |
| `fastICA` | — | native-python | not used | TSCAN's ICA path is unreachable from the canonical workflow |
| `plyr` | — | native-python | `pandas` | dplyr-style → idiomatic pandas |
| `combinat` | — | native-python | `itertools` | only needed for leaf-pair enumeration |
| `ggplot2` | — | out-of-scope | matplotlib (v0.2) | plotting deferred |
| `shiny` | — | out-of-scope | — | GUI not ported |
| `grid` | — | out-of-scope | — | plotting backend |
| `gplots` | — | out-of-scope | — | plotting |

## 4. Reusable work saved

| Reused omicverse port | Approx. LOC saved | Notes |
|---|---|---|
| `py-mclustR` | ~3000 | full Mclust EM + HC + BIC + 14 covariance models |
| **Total** | **~3000** | ~2 weeks of focused engineering |

## 5. New ports surfaced

None for this port — TSCAN's algorithmic deps either have a py- mirror, a strong native-Python equivalent, or are plotting/GUI (out of algorithmic scope).

## 6. What changed in the protocol after this port

The discovery step was **manually performed mid-Phase-2** in the original TSCAN port (when `mclust_py` was discovered by reading the omicverse_dev directory). Moving it to **Phase 0.5 (before scaffolding)** prevents next time's port from re-implementing Mclust by hand for 2 weeks.

The lesson is also encoded in [omicverse-rebuild/DISCOVERY.md §Anti-patterns](../omicverse-rebuild/DISCOVERY.md#anti-patterns): *"Don't skip the discovery step because 'I'm sure nothing exists yet.'"*
