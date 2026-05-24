"""Visualisation — 1:1 port of TSCAN::plotmclust + a pseudotime helper."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    geom_point,
    geom_segment,
    geom_text,
    ggplot,
    labs,
    scale_color_manual,
    theme_classic,
)

from .exprmclust import ExprMclustResult


def _gg_color_hue(n: int) -> list[str]:
    """Set1 brewer for ≤9, else extend with ggplot2 default hues."""
    set1 = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
        "#FFFF33", "#A65628", "#F781BF", "#999999",
    ]
    if n <= 9:
        return set1[:n]
    return set1 + [f"#{int(255 * (i / n)):02x}80a0" for i in range(n - 9)]


def plotmclust(
    result: ExprMclustResult,
    *,
    x: int = 1,
    y: int = 2,
    show_tree: bool = True,
    show_cluster: bool = True,
    cell_point_size: float = 2.0,
):
    """1:1 port of TSCAN::plotmclust.

    Args:
        result: output of ``pytscan.exprmclust``.
        x, y: 1-indexed PCA components.
        show_tree: overlay the MST between cluster centres.
        show_cluster: colour cells by cluster id.
    """
    xi, yi = int(x) - 1, int(y) - 1
    coords = np.asarray(result.pcareduceres, dtype=np.float64)[:, [xi, yi]]
    df = pd.DataFrame(coords, columns=[f"PC_dim_{x}", f"PC_dim_{y}"])
    df["State"] = [str(c) for c in result.clusterid]
    x_col, y_col = df.columns[0], df.columns[1]

    unique = sorted(df["State"].unique(), key=int)
    palette = dict(zip(unique, _gg_color_hue(len(unique))))

    if show_cluster:
        p = (
            ggplot(df, aes(x=x_col, y=y_col, colour="State"))
            + geom_point(size=cell_point_size)
            + scale_color_manual(values=palette)
        )
    else:
        p = ggplot(df, aes(x=x_col, y=y_col)) + geom_point(size=cell_point_size)

    if show_tree:
        clu_center = np.asarray(result.clucenter, dtype=np.float64)[:, [xi, yi]]
        # MST adjacency from networkx graph
        import networkx as nx
        mst = result.MSTtree
        n_clu = clu_center.shape[0]
        edges = []
        for (u, v) in mst.edges():
            ui = int(u) - 1 if isinstance(u, (int, np.integer)) and u > 0 else int(u)
            vi = int(v) - 1 if isinstance(v, (int, np.integer)) and v > 0 else int(v)
            if 0 <= ui < n_clu and 0 <= vi < n_clu:
                edges.append({
                    "x": clu_center[ui, 0], "y": clu_center[ui, 1],
                    "xend": clu_center[vi, 0], "yend": clu_center[vi, 1],
                })
        if edges:
            edge_df = pd.DataFrame(edges)
            p = p + geom_segment(
                aes(x="x", y="y", xend="xend", yend="yend"),
                data=edge_df,
                size=1.0,
                colour="black",
            )

        # Cluster id labels
        label_df = pd.DataFrame({
            "x": clu_center[:, 0],
            "y": clu_center[:, 1],
            "id": [str(i + 1) for i in range(n_clu)],
        })
        p = p + geom_text(
            aes(x="x", y="y", label="id"),
            data=label_df,
            size=8,
            colour="black",
        )

    p = p + theme_classic() + labs(x=x_col, y=y_col, colour="State")
    return p


def plotPseudotime(
    result: ExprMclustResult,
    pseudotime,
    *,
    x: int = 1,
    y: int = 2,
    cell_point_size: float = 2.0,
):
    """Convenience: PCA scatter coloured by pseudotime."""
    from ggplot2_py import scale_color_gradientn

    xi, yi = int(x) - 1, int(y) - 1
    coords = np.asarray(result.pcareduceres, dtype=np.float64)[:, [xi, yi]]
    df = pd.DataFrame(coords, columns=[f"PC_dim_{x}", f"PC_dim_{y}"])
    df["pseudotime"] = np.asarray(pseudotime, dtype=np.float64)
    x_col, y_col = df.columns[0], df.columns[1]
    viridis = [
        "#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
        "#1F9E89", "#35B779", "#6CCE59", "#B4DE2C", "#FDE725",
    ]
    return (
        ggplot(df, aes(x=x_col, y=y_col, colour="pseudotime"))
        + geom_point(size=cell_point_size)
        + scale_color_gradientn(colours=viridis)
        + theme_classic()
        + labs(x=x_col, y=y_col, colour="pseudotime")
    )
