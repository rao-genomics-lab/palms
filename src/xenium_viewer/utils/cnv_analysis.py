"""
CNV inference utilities for the Xenium viewer.

Thin orchestration around the ``insitucnv``/``infercnvpy`` packages
(https://github.com/Moldia/InSituCNV), which implement the InSituCNV method
for inferring copy-number variation from image-based spatial transcriptomics
data. This module has no Qt/napari imports and never imports insitucnv or
infercnvpy at module load time — both are optional dependencies
(``pip install -e ".[cnv]"``).

Default parameters (``smoothing_neighbors``, ``window_size``, ``step``,
``lfc_clip``) match InSituCNV's own reference notebook
(``notebooks/run_insitucnv.ipynb``), not a Xenium-specific guess.
``window_size``/``step`` are gene counts, computed independently per
chromosome: infercnvpy slides a window of ``window_size`` genes along
each chromosome's genes (ordered by genomic position), stepping by
``step``. A chromosome with fewer genes than ``window_size`` doesn't get
dropped — infercnvpy falls back to a single window averaging all of that
chromosome's available genes. So a larger window mainly trades away
sub-chromosomal resolution (most chromosomes on a small panel collapse to
one whole-chromosome average) for a less noisy per-window estimate, which
suits CNV signal that's dominated by whole-chromosome/arm-level events.
This module reports how many genes/windows were actually used so callers
can judge result quality and retune if needed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from xenium_viewer.utils.gene_analysis import get_normalized_adata, add_clustering_to_obs
from xenium_viewer.utils.adata_persistence import _convert_adata_arrow_strings

CNV_REFERENCE_OBS_KEY = "cnv_reference"


def _patch_matplotlib_cm_compat() -> None:
    """Compatibility shim: insitucnv 0.1.0 calls the removed
    ``matplotlib.cm.get_cmap()`` API (deprecated since matplotlib 3.7,
    removed in 3.9+) when building its own cluster color palette. Alias it
    back to ``matplotlib.pyplot.get_cmap`` if missing, matching this repo's
    other third-party compatibility patches (see CLAUDE.md)."""
    import matplotlib.cm as _cm
    if not hasattr(_cm, "get_cmap"):
        import matplotlib.pyplot as _plt
        _cm.get_cmap = _plt.get_cmap


def run_cnv_pipeline(
    adata,
    reference_series: pd.Series,
    reference_categories: list[str],
    reference_clustering_name: str = "",
    n_neighbors: int = 15,
    smoothing_neighbors: int = 20,
    window_size: int = 60,
    step: int = 10,
    lfc_clip: float = 4.0,
    resolution: float = 0.2,
    analyze_categories: list[str] | None = None,
    backend: str = "infercnv",
    copykat_output_dir: str | None = None,
) -> dict:
    """Run the InSituCNV pipeline on ``adata`` (raw counts expected in ``.X``).

    Parameters
    ----------
    adata : AnnData
        Raw-count AnnData (e.g. ctx.adata).
    reference_series : pd.Series
        Clustering/annotation assignment indexed by cell_id (e.g. a value
        from ``ctx.clusterings``), identifying the "normal" reference
        population via ``reference_categories``.
    reference_categories : list[str]
        Category values within ``reference_series`` to treat as the normal
        reference for inferCNV.
    reference_clustering_name : str
        Human-readable name of the source clustering (e.g. the key the user
        picked from ``ctx.clusterings``), stored in the result purely for
        display/reproducibility — not used internally (the working copy
        always uses a fixed internal obs column, ``CNV_REFERENCE_OBS_KEY``).
    n_neighbors : int
        Neighbors for the expression PCA graph used for smoothing.
    smoothing_neighbors : int
        Neighbors used by InSituCNV's graph-smoothing step.
    window_size, step : int
        infercnvpy sliding-window parameters, in number of genes (matches
        InSituCNV's reference notebook defaults — see module docstring for
        how this behaves on chromosomes with few genes).
    lfc_clip : float
        infercnvpy log-fold-change clipping value.
    resolution : float
        Leiden clustering resolution for CNV subclones. May need tuning
        per dataset — InSituCNV's own notebook evaluates several
        resolutions and picks one after reviewing the results, rather than
        recommending a single universal default.
    analyze_categories : list[str] or None
        Category values within ``reference_series`` naming the cell types to
        analyze. When given, the analysis is restricted to these cells plus
        the reference population (``reference_categories``) — every other cell
        is dropped before inference so the CNV profile, score, clustering, and
        heatmap only cover the selected cells. ``None`` (or empty) analyzes all
        cells (the previous behavior).
    backend : str
        CNV-calling engine: ``"infercnv"`` (default, infercnvpy) or ``"copykat"``
        (the insituCNV-copykat fork's R-based CopyKAT, run on the smoothed ``"M"``
        layer). Both write ``adata.obsm["X_cnv"]`` so the rest of the pipeline is
        shared; cluster columns are namespaced (``cnv_leiden_res*`` vs
        ``copykat_leiden_res*``).
    copykat_output_dir : str or None
        Directory CopyKAT writes its R-side side-effect files into (ignored for
        inferCNV). Defaults to the current working directory.

    Returns
    -------
    dict with keys: adata_cnv, cluster_key, cluster_series, cnv_score,
    n_genes_total, n_genes_mapped, n_windows, n_cells, reference_obs_key,
    reference_clustering_name, reference_categories, analyze_categories, params.
    """
    if backend not in ("infercnv", "copykat"):
        raise ValueError(f"Unknown CNV backend '{backend}'. Choose 'infercnv' or 'copykat'.")

    try:
        from insitucnv.tl import (
            prepare_cnv_input,
            run_infercnv,
            compute_cnv_neighbors,
            cluster_cnv_resolutions,
        )
    except ImportError as e:
        raise ImportError(
            "insitucnv is not installed. Run: pip install -e '.[cnv]'"
        ) from e

    _patch_matplotlib_cm_compat()

    if not reference_categories:
        raise ValueError("At least one reference category must be selected.")

    # Optionally restrict the whole analysis to the selected cell types plus the
    # reference population (inferCNV needs the reference cells as its baseline).
    if analyze_categories:
        include = {str(c) for c in analyze_categories} | {str(c) for c in reference_categories}
        idx = adata.obs["cell_id"].values if "cell_id" in adata.obs.columns else adata.obs_names
        mask = reference_series.reindex(idx).astype(str).isin(include).to_numpy()
        if not mask.any():
            raise ValueError("No cells match the selected cell types + reference population.")
        adata = adata[mask].copy()

    n_genes_total = adata.n_vars

    adata_work = get_normalized_adata(adata).copy()
    sc.pp.neighbors(adata_work, n_neighbors=n_neighbors)
    adata_work.layers["raw_counts"] = adata.X.copy()
    add_clustering_to_obs(adata_work, adata, reference_series, CNV_REFERENCE_OBS_KEY)

    adata_work = prepare_cnv_input(
        adata_work,
        raw_layer="raw_counts",
        smoothing_neighbors=smoothing_neighbors,
        add_gene_positions=True,
        drop_unmapped_genes=True,
        copy=False,
    )
    n_genes_mapped = adata_work.n_vars

    # infercnvpy does numpy-style fancy indexing on var_names/var columns,
    # which breaks on pandas 3.0's PyArrow-backed string dtype (the same
    # class of issue _convert_adata_arrow_strings already patches for the
    # zarr writer — see loader.py's _convert_arrow_strings).
    _convert_adata_arrow_strings(adata_work)

    # ── CNV calling: inferCNV or CopyKAT ────────────────────────────────
    # Both engines write adata.obsm["X_cnv"] and adata.uns["cnv"]["chr_pos"],
    # so everything downstream (neighbors, clustering, heatmap) is shared.
    if backend == "copykat":
        try:
            from insitucnv.tl import run_copykat
        except ImportError as e:
            raise RuntimeError(
                "CopyKAT backend requires the insituCNV-copykat fork with rpy2.\n"
                "Install: pip install 'git+https://github.com/sraorao/insituCNV-copykat.git' rpy2"
            ) from e
        # The copykat R package is GitHub-only (not in environment.yml); make sure
        # it's installed before the first CopyKAT run (idempotent, one-time fetch).
        from xenium_viewer.install_copykat import ensure_copykat_installed
        ensure_copykat_installed()
        try:
            # Fork defaults: input_layer='M' (spatially-smoothed), genome hg20,
            # win_size 25, ks_cut 0.1, ngene_chr 5, min_gene_per_cell 5, etc.
            run_copykat(
                adata_work,
                reference_key=CNV_REFERENCE_OBS_KEY,
                reference_categories=[str(c) for c in reference_categories],
                input_layer="M",
                output_dir=copykat_output_dir,
                copy=False,
            )
        except ImportError as e:
            raise RuntimeError(
                "CopyKAT backend requires rpy2 + R + the copykat R package.\n"
                "Install the r-* deps and: R -e "
                "'remotes::install_github(\"navinlabcode/copykat\", dependencies=FALSE)'"
            ) from e
        key_prefix = "copykat_leiden_res"
    else:
        run_infercnv(
            adata_work,
            reference_key=CNV_REFERENCE_OBS_KEY,
            reference_categories=[str(c) for c in reference_categories],
            window_size=window_size,
            step=step,
            lfc_clip=lfc_clip,
            calculate_gene_values=True,
            copy=False,
        )
        key_prefix = "cnv_leiden_res"
    n_windows = adata_work.obsm["X_cnv"].shape[1]

    compute_cnv_neighbors(adata_work, copy=False)
    cluster_keys = cluster_cnv_resolutions(
        adata_work, [resolution], key_prefix=key_prefix, dendrogram=False, copy=False
    )
    cluster_key = cluster_keys[0]

    cluster_series = adata_work.obs[cluster_key].copy()
    if "cell_id" in adata_work.obs.columns:
        cluster_series.index = adata_work.obs["cell_id"].values
    cluster_series.name = cluster_key

    X_cnv = adata_work.obsm["X_cnv"]
    X_cnv_dense = X_cnv.toarray() if hasattr(X_cnv, "toarray") else np.asarray(X_cnv)
    cnv_score = pd.Series(
        np.abs(X_cnv_dense).mean(axis=1),
        index=(adata_work.obs["cell_id"].values if "cell_id" in adata_work.obs.columns else adata_work.obs_names),
        name="cnv_score",
    )

    return {
        "adata_cnv": adata_work,
        "cluster_key": cluster_key,
        "cluster_series": cluster_series,
        "cnv_score": cnv_score,
        "n_genes_total": int(n_genes_total),
        "n_genes_mapped": int(n_genes_mapped),
        "n_windows": int(n_windows),
        "n_cells": int(adata_work.n_obs),
        "backend": backend,
        "reference_obs_key": CNV_REFERENCE_OBS_KEY,
        "reference_clustering_name": reference_clustering_name,
        "reference_categories": [str(c) for c in reference_categories],
        "analyze_categories": [str(c) for c in analyze_categories] if analyze_categories else [],
        "params": {
            "n_neighbors": n_neighbors,
            "smoothing_neighbors": smoothing_neighbors,
            "window_size": window_size,
            "step": step,
            "lfc_clip": lfc_clip,
            "resolution": resolution,
        },
    }


def make_cnv_heatmap(adata_cnv, groupby: str):
    """Build an infercnvpy chromosome heatmap figure for ``groupby``.

    Returns a matplotlib Figure (caller is responsible for
    ``ctx.auto_save_plot(fig, ...)``), matching the
    make_rank_genes_dotplot-style "return a Figure" convention.
    """
    import infercnvpy as cnv
    import matplotlib.pyplot as plt

    # Heatmap settings match the insituCNV-copykat fork's plot_chromosome_heatmap
    # (dendrogram + a fixed ±0.4 CNV colour range), used for both backends.
    cnv.pl.chromosome_heatmap(
        adata_cnv, groupby=groupby, dendrogram=True, vmin=-0.4, vmax=0.4, show=False
    )
    return plt.gcf()


def subsample_indices(
    reference_series: pd.Series,
    reference_ids: list[str],
    analyze_ids: list[str] | None,
    obs_index,
    max_cells: int,
    seed: int = 0,
) -> np.ndarray:
    """Boolean mask (aligned to ``obs_index``) selecting up to ``max_cells`` cells.

    Used for CopyKAT, which is too slow to run on a full Xenium sample. All
    reference (baseline) cells are kept first — CopyKAT needs them to seed its
    diploid baseline — then the remaining slots are filled with a seeded random
    draw of the analyzed (non-reference) cells. If reference cells alone exceed
    ``max_cells`` they are themselves subsampled. Returns all-True when the cell
    count is already <= ``max_cells``.
    """
    aligned = reference_series.reindex(obs_index).astype("object")
    values = np.array([str(v) for v in aligned.to_numpy()], dtype=object)
    n = len(values)
    keep = np.zeros(n, dtype=bool)
    if n <= max_cells:
        keep[:] = True
        return keep

    ref_set = {str(r) for r in reference_ids}
    if analyze_ids:
        analyze_set = {str(a) for a in analyze_ids} | ref_set
        in_scope = np.array([v in analyze_set for v in values], dtype=bool)
    else:
        in_scope = np.ones(n, dtype=bool)

    is_ref = np.array([v in ref_set for v in values], dtype=bool) & in_scope
    is_other = in_scope & ~is_ref
    rng = np.random.RandomState(seed)

    ref_idx = np.flatnonzero(is_ref)
    if ref_idx.size >= max_cells:
        chosen = rng.choice(ref_idx, size=max_cells, replace=False)
        keep[chosen] = True
        return keep

    keep[ref_idx] = True
    remaining = max_cells - ref_idx.size
    other_idx = np.flatnonzero(is_other)
    if other_idx.size > remaining:
        other_idx = rng.choice(other_idx, size=remaining, replace=False)
    keep[other_idx] = True
    return keep
