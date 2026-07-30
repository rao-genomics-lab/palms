"""Render the provenance DAG as a matplotlib figure (in-app visualization).

Uses networkx (already a dependency) for a layered, top-down layout — setup at
the top, artifacts in the middle, terminal plot/export leaves at the bottom —
colored by node kind, with stale nodes outlined. No graphviz needed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from xenium_viewer.utils.prov_graph import ProvGraph, SETUP, ARTIFACT, TERMINAL, NOTE

_FILL = {SETUP: "#dCE8ff", ARTIFACT: "#dcf3dc", TERMINAL: "#eeeeee",
         NOTE: "#fff6e5"}
_EDGE = {SETUP: "#4477cc", ARTIFACT: "#3a3", TERMINAL: "#999999",
         NOTE: "#cc9933"}


def _wrap(text: str, width: int = 22) -> str:
    import textwrap
    return "\n".join(textwrap.wrap(text, width)) or text


def render_dag(graph: ProvGraph, path: Optional[str | Path] = None):
    """Build (and optionally save) a matplotlib Figure of the provenance DAG.

    Returns the Figure, or None if the graph is empty.
    """
    import networkx as nx
    import matplotlib.pyplot as plt

    if len(graph) == 0:
        return None

    G = nx.DiGraph()
    for node in graph.nodes():
        G.add_node(node.id)
    for node in graph.nodes():
        for d in node.deps:
            if d in graph:
                G.add_edge(d, node.id)

    # Layered top-down layout from topological generations.
    layers = list(nx.topological_generations(G))
    pos, layer_of = {}, {}
    for li, gen in enumerate(layers):
        gen = sorted(gen)
        n = len(gen)
        for j, nid in enumerate(gen):
            pos[nid] = (j - (n - 1) / 2.0, -li)
            layer_of[nid] = li

    max_w = max((len(g) for g in layers), default=1)
    fig, ax = plt.subplots(figsize=(max(8.0, max_w * 2.4),
                                    max(4.0, len(layers) * 1.5)))

    nx.draw_networkx_edges(
        G, pos, ax=ax, arrows=True, arrowstyle="-|>", arrowsize=12,
        edge_color="#888888", width=1.2, node_size=2600,
    )

    order = [n for gen in layers for n in sorted(gen)]
    node_fill = [_FILL.get(graph.get(n).kind, "#ffffff") for n in order]
    edge_col = ["#cc6600" if graph.get(n).stale else _EDGE.get(graph.get(n).kind, "#666")
                for n in order]
    edge_w = [2.6 if graph.get(n).stale else 1.3 for n in order]
    nx.draw_networkx_nodes(
        G, pos, ax=ax, nodelist=order, node_color=node_fill, node_shape="s",
        node_size=2600, edgecolors=edge_col, linewidths=edge_w,
    )

    labels = {}
    for n in order:
        node = graph.get(n)
        lab = _wrap(node.label or n)
        if node.stale:
            lab += "\n⚠ stale"
        labels[n] = lab
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax, font_size=7)

    # Legend
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=_FILL[SETUP], edgecolor=_EDGE[SETUP], label="setup"),
        Patch(facecolor=_FILL[ARTIFACT], edgecolor=_EDGE[ARTIFACT], label="artifact"),
        Patch(facecolor=_FILL[TERMINAL], edgecolor=_EDGE[TERMINAL], label="terminal (plot/export)"),
        Patch(facecolor=_FILL[NOTE], edgecolor=_EDGE[NOTE], label="note (viewer state)"),
        Patch(facecolor="white", edgecolor="#cc6600", linewidth=2.4, label="stale"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.9)

    ax.set_title("Analysis provenance graph", fontsize=10)
    ax.margins(0.12)
    ax.axis("off")
    fig.tight_layout()

    if path is not None:
        fig.savefig(str(path), dpi=200, bbox_inches="tight")
    return fig
