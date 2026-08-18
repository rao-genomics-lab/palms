"""
Spatial analysis utilities for the Xenium viewer.

Provides:
  - Spatial neighbor graph construction via squidpy
  - L-R and co-occurrence result plotting
  - Neighborhood enrichment analysis

The ``ligrec`` and ``co_occurrence`` runners were deleted when those tabs moved
onto the step executor: they are one squidpy call each, and the templates in
``tabs/tab_ligrec.py`` / ``tabs/tab_co_occurrence.py`` now *are* that call — the
same string the viewer executes and the notebook records.
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
    """Compute spatial neighbor graph (modifies adata_norm in-place).

    ``spatial_neighbors_knn`` replaces ``spatial_neighbors(coord_type='generic')``,
    removed in squidpy 1.9. Measured on 1.8.2, the two produce byte-identical
    ``spatial_connectivities`` and ``spatial_distances`` and the same
    ``uns['spatial_neighbors']`` — ``coord_type='generic'`` with ``n_neighs=k``
    *is* k-nearest-neighbours. Keep this in step with the ``spatial_neighbors``
    template, which is what the exported notebook replays.
    """
    sq.gr.spatial_neighbors_knn(adata_norm, n_neighs=n_neighs)



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
    annotate: bool = False,
) -> plt.Figure:
    """Create a neighborhood enrichment heatmap matching squidpy native style.

    Parameters
    ----------
    result : dict
        Output from run_nhood_enrichment().
    mode : str
        'zscore' or 'count'.
    cluster_filter : list of str or None
        If provided, subset to only these cluster labels.
    cluster_labels : dict or None
        Maps cluster ID -> display label. Applied to category bar tick labels.
    annotate : bool
        If True, overlay numeric values on each cell.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from xenium_viewer.utils.coloring import CLUSTER_PALETTE

    matrix = result[mode].copy()
    clusters = list(result['clusters'])

    # Build DataFrame for subsetting
    df = pd.DataFrame(matrix, index=clusters, columns=clusters)

    # Subset if cluster filter provided (before label mapping)
    if cluster_filter:
        keep = [c for c in cluster_filter if c in df.index]
        if keep:
            df = df.loc[keep, keep]

    display_clusters = list(df.index)
    n = len(display_clusters)

    # Build display labels
    if cluster_labels:
        label_map = {str(k): v for k, v in cluster_labels.items()}
    else:
        label_map = {}
    display_names = [label_map.get(c, c) for c in display_clusters]

    # Build cluster colors (cycle palette for >20 clusters)
    palette = CLUSTER_PALETTE
    cluster_colors = [tuple(palette[i % len(palette)][:3]) for i in range(n)]

    # Colormap / normalization
    data = df.values
    if mode == 'zscore':
        cmap = plt.get_cmap('viridis')
        vmin = np.nanmin(data)
        vmax = np.nanmax(data)
        title = 'Neighborhood Enrichment (z-score)'
    else:
        cmap = plt.get_cmap('viridis')
        vmin = np.nanmin(data)
        vmax = np.nanmax(data)
        title = 'Neighborhood Enrichment (count)'

    # Figure size matching squidpy: (2*n//3, 2*n//3) with minimum (4, 4)
    sz = max(4, 2 * n // 3)
    fig, ax = plt.subplots(figsize=(sz, sz), constrained_layout=False)

    # Main heatmap via imshow
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)

    # Optional text annotations (squidpy style: white/black based on luminance)
    if annotate:
        fmt = '.2f' if mode == 'zscore' else '.0f'
        norm_data = (data - vmin) / (vmax - vmin + 1e-12)
        cmap_array = cmap(norm_data)
        for i in range(n):
            for j in range(n):
                rgb = cmap_array[i, j, :3]
                lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
                color = 'white' if lum < 0.5 else 'black'
                ax.text(j, i, f'{data[i, j]:{fmt}}',
                        ha='center', va='center', fontsize=6, color=color)

    # Divider for category bars and colorbar
    divider = make_axes_locatable(ax)

    # --- Top category bar ---
    ax_top = divider.append_axes('top', size='3%', pad=0.0)
    cat_cmap = ListedColormap(cluster_colors)
    bounds = np.arange(n + 1)
    cat_norm = BoundaryNorm(bounds, cat_cmap.N)
    cb_top = fig.colorbar(
        plt.cm.ScalarMappable(norm=cat_norm, cmap=cat_cmap),
        cax=ax_top, orientation='horizontal',
    )
    cb_top.set_ticks([])
    ax_top.set_title(title, fontsize=10)
    ax_top.xaxis.set_ticks_position('top')

    # --- Left category bar ---
    ax_left = divider.append_axes('left', size='3%', pad=0.0)
    cb_left = fig.colorbar(
        plt.cm.ScalarMappable(norm=cat_norm, cmap=cat_cmap),
        cax=ax_left, orientation='vertical',
    )
    cb_left.set_ticks(np.arange(n) + 0.5)
    cb_left.set_ticklabels(display_names)
    ax_left.invert_yaxis()
    ax_left.yaxis.set_ticks_position('left')
    ax_left.tick_params(axis='y', length=0, labelsize=8)

    # --- Right colorbar ---
    ax_cbar = divider.append_axes('right', size='3%', pad='2%')
    cbar = fig.colorbar(im, cax=ax_cbar)
    n_ticks = 5
    tick_vals = np.linspace(vmin, vmax, n_ticks)
    cbar.set_ticks(tick_vals)
    cbar.set_ticklabels([f'{v:0.2f}' for v in tick_vals])

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
        from xenium_viewer.utils.coloring import CLUSTER_PALETTE
        for i, c in enumerate(clusters):
            if c not in color_map:
                rgba = CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
                color_map[c] = tuple(rgba[:3].tolist())
    else:
        from xenium_viewer.utils.coloring import CLUSTER_PALETTE
        color_map = {}
        for i, c in enumerate(clusters):
            rgba = CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
            color_map[c] = tuple(rgba[:3].tolist())

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
