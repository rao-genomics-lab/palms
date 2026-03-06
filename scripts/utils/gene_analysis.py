"""
Gene analysis utilities for the Xenium viewer.

Provides:
  - Log-normalization with caching
  - Clustering alignment to adata obs
  - Rank genes (Wilcoxon / t-test / logreg)
  - Dotplot of top marker genes per cluster
  - Rank genes panel plot
  - ROI-based differential expression
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from typing import Optional

# Module-level cache: id(adata) -> normalized copy
_norm_cache: dict[int, sc.AnnData] = {}


def get_normalized_adata(adata: sc.AnnData) -> sc.AnnData:
    """Return a log-normalized copy of adata (cached by object id)."""
    key = id(adata)
    if key in _norm_cache:
        return _norm_cache[key]
    copy = adata.copy()
    sc.pp.normalize_total(copy, target_sum=1e4)
    sc.pp.log1p(copy)
    sc.pp.pca(copy)
    _norm_cache[key] = copy
    return copy


def add_clustering_to_obs(
    adata_norm: sc.AnnData,
    adata_orig: sc.AnnData,
    clustering_series: pd.Series,
    key_name: str,
) -> None:
    """Align clustering_series to adata_norm.obs via cell_id and store as categorical."""
    if 'cell_id' in adata_orig.obs.columns:
        cell_ids = adata_orig.obs['cell_id'].values
        aligned = clustering_series.reindex(cell_ids)
    else:
        aligned = clustering_series.reindex(adata_orig.obs_names)
    # Handle both integer and string cluster IDs (e.g. imported clusterings)
    # Fill missing with a sentinel, then convert to string categorical
    vals = aligned.values
    try:
        vals = np.where(pd.isna(vals), -1, vals).astype(int).astype(str)
    except (ValueError, TypeError):
        vals = np.where(pd.isna(vals), '-1', vals).astype(str)
    adata_norm.obs[key_name] = pd.Categorical(vals)


def run_rank_genes(
    adata_norm: sc.AnnData,
    groupby: str,
    method: str = 'wilcoxon',
    n_genes: int = 25,
) -> pd.DataFrame:
    """Run rank_genes_groups and return full results DataFrame."""
    sc.tl.rank_genes_groups(adata_norm, groupby, method=method, n_genes=n_genes)
    return sc.get.rank_genes_groups_df(adata_norm, group=None)


def make_rank_genes_dotplot(
    adata_norm: sc.AnnData,
    groupby: str,
    n_genes: int = 5,
    cluster_labels: Optional[dict] = None,
    dendrogram: bool = True,
) -> plt.Figure:
    """Create a dotplot of top marker genes per cluster. Returns matplotlib Figure."""
    # Work on a copy if we need to rename categories
    if cluster_labels:
        ad = adata_norm.copy()
        cat = ad.obs[groupby].cat
        new_cats = [cluster_labels.get(c, cluster_labels.get(int(c), c) if c.lstrip('-').isdigit() else c)
                    for c in cat.categories]
        ad.obs[groupby] = ad.obs[groupby].cat.rename_categories(new_cats)
        # Re-run rank_genes so the group names match
        method = adata_norm.uns.get('rank_genes_groups', {}).get('params', {}).get('method', 'wilcoxon')
        sc.tl.rank_genes_groups(ad, groupby, method=method, n_genes=n_genes)
    else:
        ad = adata_norm

    if dendrogram:
        try:
            sc.tl.dendrogram(ad, groupby)
        except Exception:
            dendrogram = False

    dp = sc.pl.rank_genes_groups_dotplot(
        ad, groupby=groupby, n_genes=n_genes,
        dendrogram=dendrogram, show=False, return_fig=True,
    )
    return dp


def make_rank_genes_plot(
    adata_norm: sc.AnnData,
    n_genes: int = 25,
    cluster_labels: Optional[dict] = None,
) -> plt.Figure:
    """Create the standard rank_genes_groups panel plot. Returns matplotlib Figure."""
    if cluster_labels:
        ad = adata_norm.copy()
        groupby = ad.uns.get('rank_genes_groups', {}).get('params', {}).get('groupby')
        if groupby and groupby in ad.obs.columns:
            cat = ad.obs[groupby].cat
            new_cats = [cluster_labels.get(c, cluster_labels.get(int(c), c) if c.lstrip('-').isdigit() else c)
                        for c in cat.categories]
            ad.obs[groupby] = ad.obs[groupby].cat.rename_categories(new_cats)
            method = ad.uns.get('rank_genes_groups', {}).get('params', {}).get('method', 'wilcoxon')
            sc.tl.rank_genes_groups(ad, groupby, method=method, n_genes=n_genes)
        sc.pl.rank_genes_groups(ad, n_genes=n_genes, show=False)
    else:
        sc.pl.rank_genes_groups(adata_norm, n_genes=n_genes, show=False)
    return plt.gcf()


def compute_roi_deg(
    adata: sc.AnnData,
    centroids_yx: np.ndarray,
    roi_polygons: list,
    pixel_size: float,
    cluster_mask: Optional[np.ndarray] = None,
    method: str = 'wilcoxon',
) -> pd.DataFrame:
    """Differential expression between ROI regions.

    Parameters
    ----------
    adata : raw-count AnnData
    centroids_yx : (N, 2) pixel coords (y, x)
    roi_polygons : list of Nx2 arrays in napari (y, x) coords
    pixel_size : microns per pixel
    cluster_mask : optional bool array (n_obs,) — True = include cell
    method : 'wilcoxon' or 't-test'

    Returns
    -------
    DataFrame with columns: group, names, scores, logfoldchanges, pvals, pvals_adj
    """
    from shapely.geometry import Polygon as ShapelyPolygon
    from shapely import contains_xy

    region_labels = np.full(adata.n_obs, '', dtype=object)

    for i, poly_yx in enumerate(roi_polygons):
        poly_xy = poly_yx[:, ::-1]
        shapely_poly = ShapelyPolygon(poly_xy)
        if not shapely_poly.is_valid:
            shapely_poly = shapely_poly.buffer(0)
        inside = contains_xy(shapely_poly, centroids_yx[:, 1], centroids_yx[:, 0])
        if cluster_mask is not None:
            inside = inside & cluster_mask
        region_labels[inside] = f"Region {i + 1}"

    # Keep only cells inside at least one ROI
    mask = region_labels != ''
    if mask.sum() == 0:
        return pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'])

    subset = adata[mask].copy()
    subset.obs['roi_region'] = pd.Categorical(region_labels[mask])

    unique_regions = subset.obs['roi_region'].cat.categories
    if len(unique_regions) < 2:
        return pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'])

    # Normalize the subset
    sc.pp.normalize_total(subset, target_sum=1e4)
    sc.pp.log1p(subset)

    sc.tl.rank_genes_groups(subset, 'roi_region', method=method, reference='rest')
    return sc.get.rank_genes_groups_df(subset, group=None)
