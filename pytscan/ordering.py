"""TSCANorder — pseudotime ordering via MST longest path + cluster-edge projection.

Mirrors TSCAN::TSCANorder (TSCAN/R/TSCANorder.R, 166 lines).

Algorithm:
1. Use the MST from exprmclust. Find the longest path (in #vertices, tie-break
   by total cluster size) between two leaf vertices. That path is `MSTorder`.
2. For each consecutive pair (cur, nxt) along MSTorder:
   a. In cluster `cur`: keep cells whose closest *neighbouring cluster* in the
      MST is `nxt`. Project those cells onto (center_nxt - center_cur). Sort.
   b. In cluster `nxt`: keep cells whose closest neighbouring cluster is `cur`.
      Project onto the same vector. Sort. Concatenate after (a).
3. Concatenate all edges → ordered cell list. Pseudotime = rank.

This is **the canonical TSCAN ordering** — for the default backbone path,
all cells from clusters NOT on the backbone are simply dropped (matching R
behaviour).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd

from .exprmclust import ExprMclustResult


# ----------------------------------------------------------------------------- #
# Longest-path search on the MST                                                #
# ----------------------------------------------------------------------------- #

def _find_backbone(
    mst: nx.Graph,
    clusterid: np.ndarray,
    startcluster: int | None = None,
) -> tuple[list[int], list[tuple[int, int]]]:
    """Return (backbone_path, list_of_branch_endpoints).

    backbone = longest leaf-to-leaf shortest path (R: `optcomb`).
    Tie-breaker: total #cells along the path.

    R logic:
        allcomb = expand.grid(leaves, leaves); keep i<j
        numres  = (len(path), sum(clutable[path]))
        sort by (length desc, total cells desc) → pick first → MSTorder
        branchcomb = the rest
    """
    leaves = [n for n in mst.nodes if mst.degree(n) == 1]
    # Cluster size table — match R `table(clusterid)` 1-indexed
    clutable: dict[int, int] = {}
    for cid in clusterid:
        clutable[int(cid)] = clutable.get(int(cid), 0) + 1

    if startcluster is None:
        allcomb = [
            (i, j) for i, j in itertools.combinations(leaves, 2)
            if i < j
        ]
    else:
        allcomb = [(startcluster, j) for j in leaves if j != startcluster]

    if not allcomb:
        # single-cluster or single-leaf — degenerate
        return list(mst.nodes), []

    scored: list[tuple[int, int, list[int], tuple[int, int]]] = []
    for (a, b) in allcomb:
        path = nx.shortest_path(mst, source=a, target=b)
        length = len(path)
        total_cells = sum(clutable.get(int(c), 0) for c in path)
        scored.append((length, total_cells, list(path), (a, b)))

    # sort by (length desc, total_cells desc)
    scored.sort(key=lambda t: (t[0], t[1]), reverse=True)

    best_path = scored[0][2]
    branch_combs = [s[3] for s in scored[1:]]
    return best_path, branch_combs


# ----------------------------------------------------------------------------- #
# The R `internalorderfunc` translated                                          #
# ----------------------------------------------------------------------------- #

def _ordering_along_edges(
    internal_order: list[int],
    mst_in_out: bool,
    clucenter: np.ndarray,
    clusterid: np.ndarray,
    pcareduceres: np.ndarray,
    cell_names: pd.Index,
    adjmat: np.ndarray,
) -> list[str]:
    """Direct translation of TSCAN's `internalorderfunc`."""
    ordered: list[str] = []
    n_steps = len(internal_order) - 1

    # numpy ndarrays indexed by cluster id (1..G)
    G = clucenter.shape[0]
    cell_names_arr = np.asarray(cell_names)

    for i in range(n_steps):
        cur = internal_order[i]
        nxt = internal_order[i + 1]
        cur_center = clucenter[cur - 1]
        nxt_center = clucenter[nxt - 1]

        # ---- (a) cur cluster cells whose nearest neighbour is nxt ---------
        cur_mask = clusterid == cur
        cur_res = pcareduceres[cur_mask]
        cur_cellnames = cell_names_arr[cur_mask]

        if mst_in_out:
            # connected clusters via MST adjacency
            connect = [j + 1 for j in range(G) if adjmat[cur - 1, j] == 1]
        else:
            if i == 0:
                connect = [nxt]
            else:
                connect = [nxt, internal_order[i - 1]]

        if not connect:
            continue

        # distance from each cell to each connected cluster centre
        dists = np.column_stack([
            ((cur_res - clucenter[c - 1]) ** 2).sum(axis=1)
            for c in connect
        ])
        nearest_idx = dists.argmin(axis=1)
        # Which connect-index is `nxt`?
        nxt_in_connect = connect.index(nxt)
        edge_cell_mask = nearest_idx == nxt_in_connect
        edge_cells = cur_cellnames[edge_cell_mask]
        edge_pca = cur_res[edge_cell_mask]

        if len(edge_cells) > 0:
            difvec = nxt_center - cur_center
            pos = edge_pca @ difvec
            order_idx = np.argsort(pos, kind="stable")
            ordered.extend(edge_cells[order_idx].tolist())

        # ---- (b) nxt cluster cells whose nearest neighbour is cur ---------
        nxt_mask = clusterid == nxt
        nxt_res = pcareduceres[nxt_mask]
        nxt_cellnames = cell_names_arr[nxt_mask]

        if mst_in_out:
            connect2 = [j + 1 for j in range(G) if adjmat[nxt - 1, j] == 1]
        else:
            if i == n_steps - 1:
                connect2 = [cur]
            else:
                connect2 = [cur, internal_order[i + 2]]

        if not connect2:
            continue

        if len(connect2) == 1:
            edge_cell_mask = np.ones(nxt_res.shape[0], dtype=bool)
        else:
            dists2 = np.column_stack([
                ((nxt_res - clucenter[c - 1]) ** 2).sum(axis=1)
                for c in connect2
            ])
            nearest_idx2 = dists2.argmin(axis=1)
            cur_in_connect2 = connect2.index(cur)
            edge_cell_mask = nearest_idx2 == cur_in_connect2

        edge_cells2 = nxt_cellnames[edge_cell_mask]
        edge_pca2 = nxt_res[edge_cell_mask]
        if len(edge_cells2) > 0:
            difvec = nxt_center - cur_center
            pos2 = edge_pca2 @ difvec
            order_idx2 = np.argsort(pos2, kind="stable")
            ordered.extend(edge_cells2[order_idx2].tolist())

    return ordered


# ----------------------------------------------------------------------------- #
# Public API                                                                    #
# ----------------------------------------------------------------------------- #

@dataclass
class TSCANorderResult:
    """Either a list of cell names (orderonly=True) or a DataFrame with
    sample_name, State, Pseudotime."""
    ordered_names: list[str]
    df: pd.DataFrame | None
    MSTorder: list[int]

    def __iter__(self):
        return iter(self.ordered_names)

    def __len__(self):
        return len(self.ordered_names)


def TSCANorder(
    mclustobj: ExprMclustResult,
    startcluster: int | None = None,
    MSTorder: list[int] | None = None,
    orderonly: bool = False,
    flip: bool = False,
    listbranch: bool = False,
    divide: bool = True,
) -> TSCANorderResult | dict:
    """Pure-Python equivalent of TSCAN::TSCANorder.

    Args mirror R defaults exactly. Returns:
        - if orderonly=True: TSCANorderResult with .ordered_names only;
        - else:              TSCANorderResult with .df = DataFrame
                             (sample_name, State, Pseudotime).
        - if listbranch=True and MSTorder is None: a dict
          {'backbone <path>': res, 'branch: <path>': res, ...}.
    """
    if MSTorder is not None and len(MSTorder) == 1:
        raise ValueError("MSTorder is not a path!")

    clucenter = mclustobj.clucenter
    clusterid = mclustobj.clusterid
    pcareduceres = mclustobj.pcareduceres
    cell_names = mclustobj.cell_names
    mst = mclustobj.MSTtree

    # Build dense adjacency matrix (G×G) of the MST
    G = clucenter.shape[0]
    adjmat = np.zeros((G, G), dtype=np.int64)
    for u, v in mst.edges():
        adjmat[u - 1, v - 1] = 1
        adjmat[v - 1, u - 1] = 1

    # -------------- determine MSTorder + branchcomb ---------------------
    orderinMST: bool
    branchcomb: list[tuple[int, int]] = []
    if MSTorder is None:
        MSTorder, branchcomb = _find_backbone(mst, clusterid, startcluster)
        orderinMST = True
        if flip:
            MSTorder = list(reversed(MSTorder))
    else:
        # validate the path is in the MST
        edges_in_mst = all(
            adjmat[MSTorder[i] - 1, MSTorder[i + 1] - 1] == 1
            for i in range(len(MSTorder) - 1)
        )
        if divide:
            orderinMST = bool(edges_in_mst)
        else:
            orderinMST = False

    # -------------- run the ordering ------------------------------------
    if not orderinMST:
        ordered_names = _ordering_along_edges(
            MSTorder, False, clucenter, clusterid,
            pcareduceres, cell_names, adjmat,
        )
        result_names = ordered_names
    else:
        if branchcomb and listbranch:
            allres: dict = {}
            backbone_names = _ordering_along_edges(
                MSTorder, True, clucenter, clusterid,
                pcareduceres, cell_names, adjmat,
            )
            backbone_key = f"backbone {','.join(str(x) for x in MSTorder)}"
            allres[backbone_key] = _wrap(backbone_names, clusterid,
                                         cell_names, MSTorder, orderonly)
            for (a, b) in branchcomb:
                br_path = nx.shortest_path(mst, source=a, target=b)
                if flip:
                    br_path = list(reversed(br_path))
                br_names = _ordering_along_edges(
                    br_path, True, clucenter, clusterid,
                    pcareduceres, cell_names, adjmat,
                )
                allres[f"branch: {','.join(str(x) for x in br_path)}"] = _wrap(
                    br_names, clusterid, cell_names, br_path, orderonly,
                )
            return allres
        else:
            ordered_names = _ordering_along_edges(
                MSTorder, True, clucenter, clusterid,
                pcareduceres, cell_names, adjmat,
            )
            result_names = ordered_names

    return _wrap(result_names, clusterid, cell_names, MSTorder, orderonly)


def _wrap(
    ordered_names: list[str],
    clusterid: np.ndarray,
    cell_names: pd.Index,
    MSTorder: list[int],
    orderonly: bool,
) -> TSCANorderResult:
    if orderonly:
        return TSCANorderResult(
            ordered_names=ordered_names,
            df=None,
            MSTorder=list(MSTorder),
        )
    # Build per-cell state lookup
    state_by_cell = dict(zip(cell_names, clusterid))
    df = pd.DataFrame({
        "sample_name": ordered_names,
        "State": [int(state_by_cell[n]) for n in ordered_names],
        "Pseudotime": np.arange(1, len(ordered_names) + 1, dtype=np.int64),
    })
    return TSCANorderResult(
        ordered_names=ordered_names,
        df=df,
        MSTorder=list(MSTorder),
    )
