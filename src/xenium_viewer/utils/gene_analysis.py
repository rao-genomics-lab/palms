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


# ── where a ranking lives ────────────────────────────────────────────────────
# scanpy overwrites ``uns['rank_genes_groups']`` in place, so a session that
# ranked two clusterings kept only the last, and no ``sc.pl.rank_genes_groups*``
# call could reach a specific one — they all take ``key=``. The step template
# passes ``key_added=``; everything that reads a ranking back derives the same
# name from the clustering it belongs to, through these two functions.

RANK_GENES_PREFIX = "rank_genes"
#: What scanpy writes with no ``key_added``. Still read, never written.
LEGACY_RANK_KEY = "rank_genes_groups"


def rank_genes_key(groupby: str) -> str:
    """The ``uns`` slot holding the ranking for clustering *groupby*."""
    return f"{RANK_GENES_PREFIX}_{groupby}"


def resolve_rank_key(adata, groupby: Optional[str]) -> str:
    """The keyed slot if *adata* has one, else scanpy's default.

    The fallback is what keeps caches written before the keying readable: they
    hold a single ``uns['rank_genes_groups']`` and no keyed slot at all.
    """
    if groupby:
        keyed = rank_genes_key(groupby)
        if adata is not None and keyed in adata.uns:
            return keyed
    return LEGACY_RANK_KEY


def run_rank_genes(
    adata_norm: sc.AnnData,
    groupby: str,
    method: str = 'wilcoxon',
    n_genes: int = 25,
) -> pd.DataFrame:
    """Run rank_genes_groups and return full results DataFrame."""
    key = rank_genes_key(groupby)
    sc.tl.rank_genes_groups(adata_norm, groupby, method=method, n_genes=n_genes,
                            key_added=key)
    return sc.get.rank_genes_groups_df(adata_norm, group=None, key=key)


def make_rank_genes_dotplot(
    adata_norm: sc.AnnData,
    groupby: str,
    n_genes: int = 5,
    cluster_labels: Optional[dict] = None,
    dendrogram: bool = True,
    key: Optional[str] = None,
) -> plt.Figure:
    """Create a dotplot of top marker genes per cluster. Returns matplotlib Figure.

    *key* is the ``uns`` slot the ranking was written to; it defaults through
    :func:`resolve_rank_key`, so an ``adata_norm`` restored from a cache written
    before the keying still plots.
    """
    key = key or resolve_rank_key(adata_norm, groupby)
    label_mapping = _build_label_mapping(
        adata_norm.obs[groupby].cat.categories, cluster_labels or {}
    )

    if dendrogram:
        try:
            sc.tl.dendrogram(adata_norm, groupby)
        except Exception:
            dendrogram = False

    dp = sc.pl.rank_genes_groups_dotplot(
        adata_norm, groupby=groupby, n_genes=n_genes, key=key,
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
    key: Optional[str] = None,
) -> plt.Figure:
    """Create the standard rank_genes_groups panel plot. Returns matplotlib Figure.

    *key* names the ranking to draw; the groupby is read back out of scanpy's
    own ``params`` for that slot rather than being passed in twice.
    """
    key = key or LEGACY_RANK_KEY
    groupby = adata_norm.uns.get(key, {}).get('params', {}).get('groupby')
    sc.pl.rank_genes_groups(adata_norm, n_genes=n_genes, key=key, show=False)
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

    Retained as the reference implementation the ``roi_deg`` step is checked
    against (``tests/test_spatial_roi_steps.py``). The viewer no longer calls
    it: the ROI DEG tab runs a template of the same shapely + scanpy code, so
    that what it executes is what the notebook records.

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
        return pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']), None

    subset = adata[mask].copy()
    subset.obs['roi_region'] = pd.Categorical(region_labels[mask])

    unique_regions = subset.obs['roi_region'].cat.categories
    if len(unique_regions) < 2:
        return pd.DataFrame(columns=['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']), None

    # Normalize the subset
    sc.pp.normalize_total(subset, target_sum=1e4)
    sc.pp.log1p(subset)

    sc.tl.rank_genes_groups(subset, 'roi_region', method=method, reference='rest', key_added=method)
    deg_df = sc.get.rank_genes_groups_df(subset, group=None, key=method)
    return deg_df, subset


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

    sc.tl.rank_genes_groups(subset, 'arms_cluster', method=method, reference='rest', key_added=method)
    deg_df = sc.get.rank_genes_groups_df(subset, group=None, key=method)
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
    n_label: int = 20,
) -> plt.Figure:
    """Create an EnhancedVolcano-style volcano plot from DEG results.

    Genes passing both thresholds are coloured (red = up, blue = down);
    all others are grey.  The x-axis is symmetric around zero.  Points
    outside the auto-computed display range are shown as directional
    triangle markers pinned to the axis edges.
    """
    try:
        from adjustText import adjust_text
        _has_adjusttext = True
    except ImportError:
        _has_adjusttext = False

    import seaborn as sns

    lfc   = df['logfoldchanges'].values.astype(float)
    padj  = df['pvals_adj'].values.astype(float)
    names = df['names'].values

    padj_clipped = np.clip(padj, 1e-300, 1.0)
    neg_log10    = -np.log10(padj_clipped)

    # 3-category: only genes passing BOTH thresholds get colour
    sig_p   = padj < pval_thresh
    up_mask = sig_p & (lfc >  lfc_thresh)
    dn_mask = sig_p & (lfc < -lfc_thresh)
    ns_mask = ~(up_mask | dn_mask)

    # Auto-compute symmetric display limits
    finite_lfc = lfc[np.isfinite(lfc)]
    finite_y   = neg_log10[np.isfinite(neg_log10)]
    xlim_val = max(
        float(np.nanpercentile(np.abs(finite_lfc), 99)) * 1.1 if len(finite_lfc) else lfc_thresh * 3,
        lfc_thresh * 2.5,
    )
    ylim_val = max(
        float(np.nanpercentile(finite_y, 99)) * 1.1 if len(finite_y) else 5.0,
        -np.log10(pval_thresh) * 2.5,
    )

    in_range = (
        (np.abs(lfc) <= xlim_val) & (neg_log10 <= ylim_val)
        & np.isfinite(lfc) & np.isfinite(neg_log10)
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # ── In-range scatter ──────────────────────────────────────────────────
    for mask, color, label in [
        (ns_mask & in_range, '#B3B3B3', f'NS (n={ns_mask.sum()})'),
        (up_mask & in_range, '#DC0000', f'Up-regulated (n={up_mask.sum()})'),
        (dn_mask & in_range, '#4DBBD5', f'Down-regulated (n={dn_mask.sum()})'),
    ]:
        if mask.any():
            ax.scatter(lfc[mask], neg_log10[mask], s=15, alpha=0.7, c=color,
                       edgecolors='none', label=label, zorder=2, rasterized=True)

    # ── Outlier triangle markers pinned to axis edges ─────────────────────
    _edge_x = xlim_val * 0.97
    _edge_y = ylim_val * 0.97
    for mask, color in [
        (ns_mask & ~in_range, '#B3B3B3'),
        (up_mask & ~in_range, '#DC0000'),
        (dn_mask & ~in_range, '#4DBBD5'),
    ]:
        if not mask.any():
            continue
        lfc_m   = lfc[mask]
        neg_m   = np.clip(neg_log10[mask], 0, _edge_y)
        right   = lfc_m >  xlim_val
        left    = lfc_m < -xlim_val
        top_only = (neg_log10[mask] > ylim_val) & ~right & ~left
        if right.any():
            ax.scatter(np.full(right.sum(), _edge_x), neg_m[right],
                       marker='>', s=30, c=color, edgecolors='none', alpha=0.9, zorder=3)
        if left.any():
            ax.scatter(np.full(left.sum(), -_edge_x), neg_m[left],
                       marker='<', s=30, c=color, edgecolors='none', alpha=0.9, zorder=3)
        if top_only.any():
            ox = np.clip(lfc_m[top_only], -_edge_x, _edge_x)
            ax.scatter(ox, np.full(top_only.sum(), _edge_y),
                       marker='^', s=30, c=color, edgecolors='none', alpha=0.9, zorder=3)

    # ── Threshold dashed lines ────────────────────────────────────────────
    ax.axhline(-np.log10(pval_thresh), linestyle='--', color='#333333', linewidth=0.8, alpha=0.8)
    ax.axvline( lfc_thresh,            linestyle='--', color='#333333', linewidth=0.8, alpha=0.8)
    ax.axvline(-lfc_thresh,            linestyle='--', color='#333333', linewidth=0.8, alpha=0.8)

    # ── Gene labels (within display range, sorted by padj) ───────────────
    sig_in_range = (up_mask | dn_mask) & in_range
    texts = []
    if sig_in_range.any():
        idx = np.where(sig_in_range)[0]
        for i in idx[np.argsort(padj[idx])][:n_label]:
            texts.append(ax.text(lfc[i], neg_log10[i], names[i],
                                 fontsize=8, ha='center', va='bottom'))
    if texts and _has_adjusttext:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#777777', lw=0.6),
                    expand=(1.2, 1.5))

    # ── Axes, title, theme ────────────────────────────────────────────────
    ax.set_xlim(-xlim_val, xlim_val)
    ax.set_ylim(0, ylim_val)
    ax.set_xlabel(r'$\log_2$ fold change', fontsize=12)
    ax.set_ylabel(r'$-\log_{10}$(adjusted p-value)', fontsize=12)
    fig.suptitle(f'{group_a}  vs  {group_b}', fontsize=14, fontweight='bold')
    ax.set_title(
        f'$p_{{adj}}$ cutoff: {pval_thresh}   |   $|log_2FC|$ cutoff: {lfc_thresh}',
        fontsize=9, color='#555555', pad=6,
    )
    sns.despine(ax=ax, top=True, right=True)
    ax.tick_params(labelsize=10)
    ax.legend(loc='upper right', markerscale=2, fontsize=9,
              framealpha=0.8, edgecolor='#cccccc')
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
    formats=None,
) -> int:
    """Generate volcano plots for all pairwise cluster comparisons.

    Returns the number of comparisons plotted (not the number of files: each
    one is written in every format in *formats*, which defaults to the same
    PNG + PDF pair every other figure in the viewer gets).
    """
    from xenium_viewer.utils.plot_output import DEFAULT_FORMATS, save_figure
    formats = list(formats or DEFAULT_FORMATS)
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
        save_figure(fig, [output_dir / f'volcano_{a}_vs_{b}.{ext}'
                          for ext in formats])
        plt.close(fig)

    return total
