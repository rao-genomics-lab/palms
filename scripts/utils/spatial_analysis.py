"""
Spatial analysis utilities for the Xenium viewer.

Provides:
  - Spatial neighbor graph construction via squidpy
  - Ligand-receptor interaction analysis
  - L-R result plotting
  - Neighborhood enrichment analysis
  - Co-occurrence analysis
"""

from __future__ import annotations

import numpy as np

# Patch for omnipath compatibility with numpy >= 2.0 (np.NAN was removed)
if not hasattr(np, 'NAN'):
    np.NAN = np.nan

import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import warnings


def compute_spatial_neighbors(
    adata_norm: sc.AnnData,
    n_neighs: int = 6,
) -> None:
    """Compute spatial neighbor graph (modifies adata_norm in-place)."""
    sq.gr.spatial_neighbors(
        adata_norm, coord_type='generic', n_neighs=n_neighs,
    )


def run_ligrec(
    adata_norm: sc.AnnData,
    cluster_key: str,
    n_perms: int = 1000,
    threshold: float = 0.01,
    seed: int = 42,
    interactions_params: dict | None = None,
) -> dict:
    """Run ligand-receptor interaction analysis.

    Parameters
    ----------
    interactions_params : dict or None
        Passed to ``sq.gr.ligrec(interactions_params=...)``.
        Keys may include ``"include"`` (tuple of InteractionDataset enums)
        and ``"resources"`` (e.g. ``"CellPhoneDB"``).

    Returns dict with keys 'means', 'pvalues' (both DataFrames),
    and 'warning' (str or None).
    """
    warning_msg = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = sq.gr.ligrec(
                adata_norm,
                cluster_key=cluster_key,
                n_perms=n_perms,
                threshold=threshold,
                seed=seed,
                use_raw=False,
                copy=True,
                transmitter_params={"categories": "ligand"},
                receiver_params={"categories": "receptor"},
                interactions_params=interactions_params or {},
            )
        means = result['means']
        pvalues = result['pvalues']

        if means.shape[0] == 0:
            warning_msg = (
                "No ligand-receptor interactions found. "
                "The 480-gene Xenium panel may not contain enough L-R pairs."
            )
    except Exception as e:
        means = pd.DataFrame()
        pvalues = pd.DataFrame()
        warning_msg = f"L-R analysis failed: {e}"

    return {
        'means': means,
        'pvalues': pvalues,
        'warning': warning_msg,
    }


def run_nhood_enrichment(
    adata_norm: sc.AnnData,
    cluster_key: str,
    n_perms: int = 1000,
    seed: int = 42,
) -> dict:
    """Run neighborhood enrichment analysis.

    Returns dict with keys 'zscore' (NxN ndarray), 'count' (NxN ndarray),
    'clusters' (list of str), and 'warning' (str or None).
    """
    warning_msg = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sq.gr.nhood_enrichment(
                adata_norm,
                cluster_key=cluster_key,
                n_perms=n_perms,
                seed=seed,
            )
        result = adata_norm.uns[f'{cluster_key}_nhood_enrichment']
        zscore = np.array(result['zscore'])
        count = np.array(result['count'])
        clusters = list(adata_norm.obs[cluster_key].cat.categories.astype(str))
    except Exception as e:
        zscore = np.array([])
        count = np.array([])
        clusters = []
        warning_msg = f"Neighborhood enrichment failed: {e}"

    return {
        'zscore': zscore,
        'count': count,
        'clusters': clusters,
        'warning': warning_msg,
    }


def make_nhood_enrichment_plot(
    result: dict,
    mode: str = 'zscore',
    cluster_filter: list[str] | None = None,
    cluster_labels: dict | None = None,
) -> plt.Figure:
    """Create a neighborhood enrichment heatmap. Returns matplotlib Figure.

    Parameters
    ----------
    result : dict
        Output from run_nhood_enrichment().
    mode : str
        'zscore' or 'count'.
    cluster_filter : list of str or None
        If provided, subset to only these cluster labels.
    cluster_labels : dict or None
        Maps cluster ID -> display label. Applied to axis tick labels.
    """
    import seaborn as sns

    matrix = result[mode].copy()
    clusters = list(result['clusters'])

    # Build DataFrame for labeling
    df = pd.DataFrame(matrix, index=clusters, columns=clusters)

    # Apply cluster labels to index/columns
    if cluster_labels:
        label_map = {str(k): v for k, v in cluster_labels.items()}
        df.rename(index=label_map, columns=label_map, inplace=True)

    # Subset if cluster filter provided
    if cluster_filter:
        keep = [c for c in cluster_filter if c in df.index]
        if keep:
            df = df.loc[keep, keep]

    if mode == 'zscore':
        cmap = 'coolwarm'
        vmax = max(abs(df.values.min()), abs(df.values.max()), 1.0)
        vmin = -vmax
        fmt = '.1f'
        title = 'Neighborhood Enrichment (z-score)'
    else:
        cmap = 'YlOrRd'
        vmin = None
        vmax = None
        fmt = '.0f'
        title = 'Neighborhood Enrichment (count)'

    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.6), max(5, len(df) * 0.5)))
    sns.heatmap(
        df, annot=True, fmt=fmt, cmap=cmap,
        vmin=vmin, vmax=vmax,
        linewidths=0.5, ax=ax,
        square=True,
    )
    ax.set_title(title)
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Cluster')
    fig.tight_layout()
    return fig


def make_ligrec_plot(
    result: dict,
    pvalue_threshold: float = 0.05,
    source_groups: list[str] | None = None,
    target_groups: list[str] | None = None,
    cluster_labels: dict | None = None,
) -> plt.Figure:
    """Create a ligand-receptor interaction dot plot. Returns matplotlib Figure."""
    # Apply cluster labels to the result's MultiIndex column levels
    if cluster_labels:
        label_map = {str(k): v for k, v in cluster_labels.items()}
        plot_result = {}
        for key in ('means', 'pvalues'):
            df = result[key]
            if df is not None and not df.empty and isinstance(df.columns, pd.MultiIndex):
                new_levels = []
                for level in df.columns.levels:
                    new_levels.append(level.map(lambda x: label_map.get(str(x), x)))
                new_cols = df.columns.set_levels(new_levels)
                df = df.copy()
                df.columns = new_cols
            plot_result[key] = df
        # Map source/target groups to original IDs for filtering, then use labelled result
        if source_groups:
            source_groups = [label_map.get(str(g), g) for g in source_groups]
        if target_groups:
            target_groups = [label_map.get(str(g), g) for g in target_groups]
    else:
        plot_result = result

    sq.pl.ligrec(
        plot_result,
        pvalue_threshold=pvalue_threshold,
        source_groups=source_groups,
        target_groups=target_groups,
        show=False,
    )
    return plt.gcf()


def run_co_occurrence(
    adata_norm: sc.AnnData,
    cluster_key: str,
    interval: int = 50,
    seed: int = 42,
) -> dict:
    """Run spatial co-occurrence analysis.

    Returns dict with keys 'occ' (n_clusters x n_clusters x n_intervals-1),
    'interval' (n_intervals,), 'clusters' (list of str), and 'warning' (str or None).
    """
    warning_msg = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sq.gr.co_occurrence(
                adata_norm,
                cluster_key=cluster_key,
                interval=interval,
            )
        result = adata_norm.uns[f'{cluster_key}_co_occurrence']
        occ = np.array(result['occ'])
        interval_arr = np.array(result['interval'])
        clusters = list(adata_norm.obs[cluster_key].cat.categories.astype(str))
    except Exception as e:
        occ = np.array([])
        interval_arr = np.array([])
        clusters = []
        warning_msg = f"Co-occurrence analysis failed: {e}"

    return {
        'occ': occ,
        'interval': interval_arr,
        'clusters': clusters,
        'warning': warning_msg,
    }


def make_co_occurrence_plot(
    result: dict,
    clusters_to_plot: list[str] | None = None,
    target_clusters: list[str] | None = None,
    cluster_colors: dict | None = None,
    cluster_labels: dict | None = None,
) -> plt.Figure:
    """Create co-occurrence line plots. Returns matplotlib Figure.

    Parameters
    ----------
    result : dict
        Output from run_co_occurrence().
    clusters_to_plot : list of str or None
        Cluster IDs to create subplots for. If None, plot all clusters.
    target_clusters : list of str or None
        Cluster IDs to draw as target lines on each subplot.
        If None, draw all clusters (current behavior).
    cluster_colors : dict or None
        Maps int cluster_id -> RGBA numpy array (from CellColorManager).
        If provided, line colors will match the napari cell coloring palette.
    """
    import seaborn as sns

    occ = result['occ']
    interval = result['interval']
    clusters = list(result['clusters'])
    distances = interval[1:]  # bin right edges

    # Determine which clusters get subplots
    if clusters_to_plot:
        plot_ids = [c for c in clusters_to_plot if c in clusters]
    else:
        plot_ids = clusters

    n_plots = len(plot_ids)
    if n_plots == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No clusters to plot", ha='center', va='center',
                transform=ax.transAxes)
        return fig

    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5 * n_cols, 4 * n_rows),
                             squeeze=False)

    if cluster_colors is not None:
        # Convert int cluster_id -> RGBA to string cluster_id -> RGB tuple
        color_map = {str(cid): tuple(rgba[:3]) for cid, rgba in cluster_colors.items()}
        # Fallback for any cluster not in the dict
        fallback = sns.color_palette("tab20", n_colors=len(clusters))
        for i, c in enumerate(clusters):
            if c not in color_map:
                color_map[c] = fallback[i]
    else:
        palette = sns.color_palette("tab20", n_colors=len(clusters))
        color_map = {c: palette[i] for i, c in enumerate(clusters)}

    # Determine which clusters to draw as target lines
    if target_clusters:
        targets = [c for c in target_clusters if c in clusters]
    else:
        targets = clusters

    # Build label map for display names
    label_map = {}
    if cluster_labels:
        label_map = {str(k): v for k, v in cluster_labels.items()}

    for idx, query_cluster in enumerate(plot_ids):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]
        qi = clusters.index(query_cluster)

        for target_cluster in targets:
            ti = clusters.index(target_cluster)
            display_target = label_map.get(target_cluster, target_cluster)
            ax.plot(distances, occ[qi, ti, :],
                    label=display_target, color=color_map[target_cluster],
                    linewidth=1.2)

        ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        display_query = label_map.get(query_cluster, query_cluster)
        ax.set_title(f"Cluster {display_query}", fontsize=10)
        ax.set_xlabel("Distance")
        ax.set_ylabel("Co-occurrence score")

    # Remove empty subplots
    for idx in range(n_plots, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].set_visible(False)

    # Single legend outside the last used subplot
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', title='Target cluster',
               fontsize=8, title_fontsize=9, bbox_to_anchor=(1.0, 0.5))

    fig.tight_layout(rect=[0, 0, 0.88, 1.0])
    return fig
