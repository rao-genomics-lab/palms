"""
Spatial analysis utilities for the Xenium viewer.

Provides:
  - Spatial neighbor graph construction via squidpy
  - Ligand-receptor interaction analysis
  - L-R result plotting
  - Neighborhood enrichment analysis
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
) -> dict:
    """Run ligand-receptor interaction analysis.

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
    """
    import seaborn as sns

    matrix = result[mode].copy()
    clusters = list(result['clusters'])

    # Build DataFrame for labeling
    df = pd.DataFrame(matrix, index=clusters, columns=clusters)

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
) -> plt.Figure:
    """Create a ligand-receptor interaction dot plot. Returns matplotlib Figure."""
    sq.pl.ligrec(
        result,
        pvalue_threshold=pvalue_threshold,
        source_groups=source_groups,
        target_groups=target_groups,
        show=False,
    )
    return plt.gcf()
