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

import itertools
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from typing import Optional, Callable

# Module-level cache: id(adata) -> normalized copy
_norm_cache: dict[int, sc.AnnData] = {}


def _build_label_mapping(categories, cluster_labels: dict) -> dict[str, str]:
    """Build str(cluster_id) -> label mapping, handling int or str keys."""
    mapping = {}
    for c in categories:
        label = cluster_labels.get(
            c, cluster_labels.get(int(c) if str(c).lstrip('-').isdigit() else c, str(c))
        )
        mapping[str(c)] = str(label)
    return mapping


def _relabel_axes(fig: plt.Figure, mapping: dict[str, str]) -> None:
    """Replace integer cluster IDs in axis tick labels with user-assigned labels."""
    for ax in fig.get_axes():
        for get_fn, set_fn in (
            (ax.get_xticklabels, ax.set_xticklabels),
            (ax.get_yticklabels, ax.set_yticklabels),
        ):
            ticks = get_fn()
            if not ticks:
                continue
            texts = [t.get_text() for t in ticks]
            new_texts = [mapping.get(t, t) for t in texts]
            if new_texts != texts:
                set_fn(new_texts,
                       rotation=ticks[0].get_rotation(),
                       fontsize=ticks[0].get_fontsize())


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
    label_mapping = _build_label_mapping(
        adata_norm.obs[groupby].cat.categories, cluster_labels or {}
    )

    if dendrogram:
        try:
            sc.tl.dendrogram(adata_norm, groupby)
        except Exception:
            dendrogram = False

    dp = sc.pl.rank_genes_groups_dotplot(
        adata_norm, groupby=groupby, n_genes=n_genes,
        dendrogram=dendrogram, show=False, return_fig=True,
    )
    if label_mapping:
        if isinstance(dp, plt.Figure):
            fig = dp
        elif hasattr(dp, 'fig') and dp.fig is not None:
            fig = dp.fig
        elif hasattr(dp, 'figure') and dp.figure is not None:
            fig = dp.figure
        else:
            fig = None
        if fig is not None:
            _relabel_axes(fig, label_mapping)
    return dp


def make_rank_genes_plot(
    adata_norm: sc.AnnData,
    n_genes: int = 25,
    cluster_labels: Optional[dict] = None,
) -> plt.Figure:
    """Create the standard rank_genes_groups panel plot. Returns matplotlib Figure."""
    groupby = adata_norm.uns.get('rank_genes_groups', {}).get('params', {}).get('groupby')
    sc.pl.rank_genes_groups(adata_norm, n_genes=n_genes, show=False)
    fig = plt.gcf()
    if cluster_labels and groupby and groupby in adata_norm.obs.columns:
        label_mapping = _build_label_mapping(
            adata_norm.obs[groupby].cat.categories, cluster_labels
        )
        _relabel_axes(fig, label_mapping)
    return fig


def run_celltypist_annotation(adata, model_name):
    """Run CellTypist annotation and return per-cell predictions with confidence.

    Parameters
    ----------
    adata : AnnData
        Log-normalized AnnData (from get_normalized_adata()).
    model_name : str
        Name of a CellTypist model (e.g. 'Immune_All_Low.pkl').

    Returns
    -------
    (predictions, confidence) : tuple[pd.Series, pd.Series]
        predictions: cell type strings, indexed by obs_names.
        confidence: max probability per cell (0–1), indexed by obs_names.
    """
    import celltypist
    from celltypist import models
    model = models.Model.load(model_name)
    result = celltypist.annotate(adata, model=model, majority_voting=False)
    predictions = result.predicted_labels['predicted_labels']
    confidence = result.probability_matrix.max(axis=1)
    return predictions, confidence


def load_reference_h5ad(path: str | Path) -> sc.AnnData:
    """Load a reference scRNA-seq h5ad file."""
    path = Path(path).expanduser()
    return sc.read_h5ad(path)


def get_annotation_columns(ref_adata: sc.AnnData, max_cardinality: int = 200) -> list[str]:
    """Return obs columns suitable as annotation labels (categorical/object/string, 2–200 unique)."""
    cols = []
    for col in ref_adata.obs.columns:
        dtype = ref_adata.obs[col].dtype
        if isinstance(dtype, pd.CategoricalDtype) or dtype == object or dtype.name == 'string':
            n_unique = ref_adata.obs[col].nunique()
            if 2 <= n_unique <= max_cardinality:
                cols.append(col)
    return cols


def run_label_transfer(
    xenium_adata: sc.AnnData,
    ref_adata: sc.AnnData,
    annotation_col: str,
) -> tuple[pd.Series, int]:
    """Transfer labels from reference to Xenium data via sc.tl.ingest().

    Parameters
    ----------
    xenium_adata : AnnData
        Raw-count Xenium AnnData.
    ref_adata : AnnData
        Reference scRNA-seq AnnData with annotation_col in obs.
    annotation_col : str
        Column in ref_adata.obs containing cell type labels.

    Returns
    -------
    (predictions, n_common_genes) : tuple[pd.Series, int]
        predictions: transferred labels indexed by xenium obs_names.
        n_common_genes: number of common genes used.
    """
    # Find common genes
    common_genes = list(set(xenium_adata.var_names) & set(ref_adata.var_names))
    n_common = len(common_genes)
    if n_common < 10:
        raise ValueError(
            f"Only {n_common} common genes between Xenium and reference — "
            f"need at least 10 for label transfer."
        )
    common_genes = sorted(common_genes)

    # Subset both to common genes
    ref_sub = ref_adata[:, common_genes].copy()
    xen_sub = xenium_adata[:, common_genes].copy()

    # Preprocess reference: normalize → log1p → HVG → PCA → neighbors → UMAP
    sc.pp.normalize_total(ref_sub, target_sum=1e4)
    sc.pp.log1p(ref_sub)
    sc.pp.highly_variable_genes(ref_sub)
    sc.pp.pca(ref_sub)
    sc.pp.neighbors(ref_sub)
    sc.tl.umap(ref_sub)

    # Preprocess Xenium subset: normalize → log1p
    sc.pp.normalize_total(xen_sub, target_sum=1e4)
    sc.pp.log1p(xen_sub)

    # Label transfer via ingest
    sc.tl.ingest(xen_sub, ref_sub, obs=annotation_col)

    predictions = xen_sub.obs[annotation_col]
    predictions.index = xenium_adata.obs_names
    return predictions, n_common


def build_llm_annotation_prompt(rank_df: pd.DataFrame, n_genes: int = 10) -> str:
    """Build an LLM prompt from rank genes results asking for cell type annotation."""
    cluster_lines = []
    for group, grp_df in rank_df.groupby("group"):
        top_genes = grp_df.head(n_genes)["names"].tolist()
        cluster_lines.append(f"Cluster {group}: {', '.join(top_genes)}")

    cluster_text = "\n".join(cluster_lines)
    return (
        "You are a bioinformatics expert. I have single-cell RNA-seq data clustered "
        f"into {len(cluster_lines)} groups. Below are the top {n_genes} differentially "
        "expressed marker genes for each cluster (ranked by statistical significance).\n\n"
        "Identify the most likely cell type for each cluster based on these marker genes.\n\n"
        "Return ONLY a valid JSON object mapping cluster ID (as string) to cell type name. "
        'Example: {"0": "CD8+ T cells", "1": "Fibroblasts", ...}\n\n'
        f"{cluster_text}"
    )


def parse_llm_annotation_response(stdout: str) -> dict:
    """Extract a JSON dict from LLM stdout, handling markdown fences and extra text."""
    # Try to find JSON in markdown code blocks first
    fence_match = re.search(r"```(?:json)?\s*\n?({.*?})\s*\n?```", stdout, re.DOTALL)
    if fence_match:
        return json.loads(fence_match.group(1))

    # Try to find a JSON object directly
    brace_match = re.search(r"\{[^{}]*\}", stdout, re.DOTALL)
    if brace_match:
        return json.loads(brace_match.group(0))

    raise ValueError(f"Could not parse JSON from LLM response:\n{stdout[:500]}")


def run_llm_annotation(rank_df: pd.DataFrame, cli: str, n_genes: int = 10) -> dict:
    """Run LLM-based cluster annotation via a local CLI tool.

    Parameters
    ----------
    rank_df : DataFrame
        Output of sc.get.rank_genes_groups_df (columns: group, names, scores, ...).
    cli : str
        One of 'claude', 'gemini', 'codex'.
    n_genes : int
        Number of top genes per cluster to include in the prompt.

    Returns
    -------
    dict mapping cluster ID (str) to cell type name (str).
    """
    prompt = build_llm_annotation_prompt(rank_df, n_genes=n_genes)

    if cli == "claude":
        cmd = ["claude", "-p", prompt, "--model", "haiku", "--allowedTools", "", "--max-turns", "1"]
    elif cli == "gemini":
        cmd = ["gemini", "-p", prompt]
    elif cli == "codex":
        cmd = ["codex", "-q", prompt, "--full-auto"]
    else:
        raise ValueError(f"Unknown CLI: {cli}")

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(f"{cli} failed (exit {result.returncode}): {result.stderr[:500]}")

    return parse_llm_annotation_response(result.stdout)


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


def compute_arms_tile_deg(
    adata: sc.AnnData,
    centroids_yx: np.ndarray,
    tile_polygons_yx: list,
    cluster_ids: np.ndarray,
    min_cells_per_cluster: int = 10,
    method: str = 'wilcoxon',
    cluster_mask: Optional[np.ndarray] = None,
) -> tuple:
    """Differential expression between ARMS tile clusters.

    Parameters
    ----------
    adata : raw-count AnnData
    centroids_yx : (N, 2) pixel coords (y, x)
    tile_polygons_yx : list of Nx2 arrays in (y, x) pixel coords (already transformed)
    cluster_ids : array of cluster IDs, one per tile polygon
    min_cells_per_cluster : minimum cells to keep a cluster
    method : 'wilcoxon' or 't-test'

    Returns
    -------
    (deg_df, summary_dict, adata_norm) where summary_dict = {cluster_id: n_cells}
    and adata_norm is the normalized subset AnnData (for pairwise volcano plots).
    When no valid clusters exist, adata_norm is None.
    """
    from shapely.geometry import Polygon as ShapelyPolygon
    from shapely import contains_xy

    cell_cluster = np.full(adata.n_obs, -1, dtype=int)

    for i, poly_yx in enumerate(tile_polygons_yx):
        poly_xy = poly_yx[:, ::-1]
        shapely_poly = ShapelyPolygon(poly_xy)
        if not shapely_poly.is_valid:
            shapely_poly = shapely_poly.buffer(0)
        inside = contains_xy(shapely_poly, centroids_yx[:, 1], centroids_yx[:, 0])
        if cluster_mask is not None:
            inside = inside & cluster_mask
        cell_cluster[inside] = cluster_ids[i]

    # Drop cells not in any tile
    mask = cell_cluster >= 0
    if mask.sum() == 0:
        empty = pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'])
        return (empty, {}, None)

    subset = adata[mask].copy()
    subset.obs['arms_cluster'] = pd.Categorical(cell_cluster[mask].astype(str))

    # Count cells per cluster and drop small ones
    cluster_counts = subset.obs['arms_cluster'].value_counts()
    summary = {int(k): int(v) for k, v in cluster_counts.items() if k != '-1'}
    keep_clusters = [c for c, n in cluster_counts.items() if n >= min_cells_per_cluster and c != '-1']

    if len(keep_clusters) < 2:
        empty = pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'])
        return (empty, summary, None)

    subset = subset[subset.obs['arms_cluster'].isin(keep_clusters)].copy()
    subset.obs['arms_cluster'] = pd.Categorical(subset.obs['arms_cluster'].values)

    sc.pp.normalize_total(subset, target_sum=1e4)
    sc.pp.log1p(subset)

    sc.tl.rank_genes_groups(subset, 'arms_cluster', method=method, reference='rest')
    deg_df = sc.get.rank_genes_groups_df(subset, group=None)
    return (deg_df, summary, subset)


def run_pairwise_deg(
    adata_norm: sc.AnnData,
    groupby: str,
    group_a: str,
    group_b: str,
    method: str = 'wilcoxon',
) -> pd.DataFrame:
    """Run DEG for group_a vs group_b and return results DataFrame."""
    sc.tl.rank_genes_groups(
        adata_norm, groupby,
        groups=[str(group_a)], reference=str(group_b),
        method=method, n_genes=adata_norm.n_vars,
    )
    return sc.get.rank_genes_groups_df(adata_norm, group=str(group_a))


def make_volcano_plot(
    df: pd.DataFrame,
    group_a: str,
    group_b: str,
    lfc_thresh: float = 0.5,
    pval_thresh: float = 0.05,
    n_label: int = 10,
) -> plt.Figure:
    """Create a volcano plot from DEG results. Returns matplotlib Figure."""
    lfc = df['logfoldchanges'].values.astype(float)
    padj = df['pvals_adj'].values.astype(float)
    names = df['names'].values

    # Clip tiny p-values to avoid -log10(0)
    padj_clipped = np.clip(padj, 1e-300, 1.0)
    neg_log10 = -np.log10(padj_clipped)

    # Classify genes
    sig = padj < pval_thresh
    up = sig & (lfc > lfc_thresh)
    down = sig & (lfc < -lfc_thresh)
    ns = ~(up | down)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(lfc[ns], neg_log10[ns], s=4, alpha=0.5, c='#aaaaaa', edgecolors='none', label='NS')
    ax.scatter(lfc[up], neg_log10[up], s=4, alpha=0.5, c='#d62728', edgecolors='none', label='Up')
    ax.scatter(lfc[down], neg_log10[down], s=4, alpha=0.5, c='#1f77b4', edgecolors='none', label='Down')

    # Threshold lines
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='gray', alpha=0.5)
    ax.axvline(lfc_thresh, linestyle='--', color='gray', alpha=0.5)
    ax.axvline(-lfc_thresh, linestyle='--', color='gray', alpha=0.5)

    # Label top significant genes
    sig_mask = up | down
    if sig_mask.any():
        sig_idx = np.where(sig_mask)[0]
        # Sort by padj (smallest first)
        order = sig_idx[np.argsort(padj[sig_idx])]
        for idx in order[:n_label]:
            ax.annotate(
                names[idx],
                (lfc[idx], neg_log10[idx]),
                fontsize=7, textcoords='offset points', xytext=(4, 4),
            )

    ax.set_xlabel('log2 fold change')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(f'{group_a} vs {group_b}')
    ax.legend(markerscale=3, framealpha=0.8)
    fig.tight_layout()
    return fig


def generate_all_volcano_plots(
    adata_norm: sc.AnnData,
    groupby: str,
    method: str = 'wilcoxon',
    output_dir: str | Path = '.',
    lfc_thresh: float = 0.5,
    pval_thresh: float = 0.05,
    n_label: int = 10,
    progress_callback: Optional[Callable] = None,
) -> int:
    """Generate volcano PNGs for all pairwise cluster comparisons.

    Returns the number of plots generated.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    groups = sorted(
        [g for g in adata_norm.obs[groupby].cat.categories if str(g) != '-1'],
        key=lambda x: (int(x) if str(x).lstrip('-').isdigit() else 0, str(x)),
    )
    pairs = list(itertools.combinations(groups, 2))
    total = len(pairs)

    # Use non-interactive backend for thread safety
    backend = matplotlib.get_backend()

    for i, (a, b) in enumerate(pairs):
        if progress_callback:
            progress_callback(i, total, str(a), str(b))

        df = run_pairwise_deg(adata_norm, groupby, str(a), str(b), method=method)
        fig = make_volcano_plot(df, str(a), str(b), lfc_thresh, pval_thresh, n_label)
        fig.savefig(output_dir / f'volcano_{a}_vs_{b}.png', dpi=300)
        plt.close(fig)

    return total
