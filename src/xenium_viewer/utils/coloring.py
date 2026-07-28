"""
Cell coloring utilities for the Xenium viewer.

Strategy: color the cell_labels raster mask using napari's DirectLabelColormap.
Each pixel value k (1..N_cells) corresponds to a cell; we map it to an RGBA
color derived from gene expression or cluster assignment.

Key class: CellColorManager
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from functools import lru_cache
from typing import Optional

# Matplotlib colormaps for gene expression
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 256 maximally-distinct categorical colors for clusters (from colorcet glasbey_dark)
import colorcet as _cc
CLUSTER_PALETTE = np.array(
    [[int(h[i:i+2], 16) / 255.0 for i in (1, 3, 5)] + [1.0]
     for h in _cc.glasbey_dark],
    dtype=np.float32,
)

# 10 maximally-distinct RGBA colors for multi-gene transcript overlay
TRANSCRIPT_PALETTE = np.array([
    [1.00, 1.00, 0.00, 1.0],  # yellow
    [0.00, 1.00, 1.00, 1.0],  # cyan
    [1.00, 0.00, 1.00, 1.0],  # magenta
    [1.00, 0.60, 0.00, 1.0],  # orange
    [0.00, 1.00, 0.00, 1.0],  # green
    [0.30, 0.70, 1.00, 1.0],  # sky blue
    [1.00, 0.00, 0.00, 1.0],  # red
    [0.60, 0.30, 0.90, 1.0],  # violet
    [1.00, 0.50, 0.70, 1.0],  # pink
    [0.60, 0.40, 0.20, 1.0],  # brown
], dtype=np.float32)

AVAILABLE_COLORMAPS = ["viridis", "magma", "plasma", "RdBu_r", "YlOrRd"]


# ── Categorical palettes for patch overlays ─────────────────────────────────
def _mpl_palette(name: str, n: int) -> np.ndarray:
    cmap = plt.get_cmap(name)
    return np.array(
        [list(cmap(i)[:3]) + [1.0] for i in range(n)],
        dtype=np.float32,
    )


TAB10_PALETTE = _mpl_palette("tab10", 10)
TAB20_PALETTE = _mpl_palette("tab20", 20)
SET1_PALETTE = _mpl_palette("Set1", 9)
SET3_PALETTE = _mpl_palette("Set3", 12)

# ARMS palette: RColorBrewer Set1(8) + Set2(8) + Dark2(8) — matches the
# cluster colours from the ARMS R package (create_cluster_colors()).
# 1-based indexing: cluster 1 → position 0, etc.
_ARMS_HEX = [
    # Set1 (8)
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
    # Set2 (8)
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
    # Dark2 (8)
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666",
]
ARMS_PALETTE = np.array(
    [[int(h[i:i+2], 16) / 255.0 for i in (1, 3, 5)] + [1.0]
     for h in _ARMS_HEX],
    dtype=np.float32,
)

PATCH_PALETTES = {
    "tab10": TAB10_PALETTE,
    "tab20": TAB20_PALETTE,
    "glasbey_dark": CLUSTER_PALETTE,
    "Set1": SET1_PALETTE,
    "Set3": SET3_PALETTE,
    "ARMS (Set1+Set2+Dark2)": ARMS_PALETTE,
}


class CellColorManager:
    """
    Manages RGBA color arrays for the cell_labels layer.

    The color array has shape (max_label + 1, 4) where index 0 = background
    (transparent). Index k maps to the RGBA color of cell with label value k.

    Parameters
    ----------
    adata : anndata.AnnData
        The cells AnnData table from sdata["table"].
    label_to_obs : np.ndarray
        Mapping from label integer value to obs row index.
        Shape: (max_label + 1,). Value -1 means no cell.
    """

    def __init__(self, adata, label_to_obs: np.ndarray):
        self.adata = adata
        self.label_to_obs = label_to_obs
        self._gene_cache: dict[tuple, np.ndarray] = {}
        self._cluster_cache: dict[str, np.ndarray] = {}
        self._continuous_cache: dict[str, np.ndarray] = {}
        self._max_label = len(label_to_obs) - 1

    def invalidate_cluster_cache(self, name: str | None = None) -> None:
        """Drop cached cluster colors for *name* (or all of them).

        ``get_cluster_colors`` caches on the series' ``name``, which is the
        clustering key. Any producer that *replaces* the series behind an
        existing key — re-running Leiden at the same resolution, re-importing a
        file, a new CNV run — must call this, or the raster keeps the previous
        run's colors while the legend and cluster filter are rebuilt from the
        new assignment. The mismatch reads as a clustering that was only
        partially overwritten.
        """
        if name is None:
            self._cluster_cache.clear()
        else:
            self._cluster_cache.pop(name, None)

    def _empty_color_array(self) -> np.ndarray:
        """Return transparent (alpha=0) RGBA array for all labels."""
        arr = np.zeros((self._max_label + 1, 4), dtype=np.float32)
        return arr

    def get_gene_colors(
        self,
        gene_name: str,
        colormap: str = "viridis",
        min_alpha: float = 0.0,
        zero_alpha: float = 0.0,
    ) -> np.ndarray:
        """
        Build RGBA color array for a gene's expression.

        Parameters
        ----------
        gene_name : str
            Must be in adata.var_names.
        colormap : str
            Matplotlib colormap name.
        min_alpha : float
            Alpha for cells with non-zero expression (0–1).
        zero_alpha : float
            Alpha for cells with zero expression (0 = transparent).

        Returns
        -------
        np.ndarray, shape (max_label + 1, 4), dtype float32
        """
        cache_key = (gene_name, colormap)
        if cache_key in self._gene_cache:
            return self._gene_cache[cache_key]

        if gene_name not in self.adata.var_names:
            raise ValueError(f"Gene '{gene_name}' not found in AnnData")

        # Extract expression for this gene (sparse → dense)
        gene_idx = self.adata.var_names.get_loc(gene_name)
        X = self.adata.X
        if hasattr(X, "toarray"):
            expr = np.asarray(X[:, gene_idx].toarray()).ravel().astype(np.float32)
        else:
            expr = np.asarray(X[:, gene_idx]).ravel().astype(np.float32)

        # Normalise non-zero cells to [0, 1] using min/max of non-zero values
        # This spreads the full colormap range across expressing cells
        nonzero = expr > 0
        expr_norm = np.zeros_like(expr)
        if nonzero.any():
            vmin = expr[nonzero].min()
            vmax = expr[nonzero].max()
            if vmax > vmin:
                expr_norm[nonzero] = (expr[nonzero] - vmin) / (vmax - vmin)
            else:
                expr_norm[nonzero] = 1.0  # all same value → top of colormap

        # Map via colormap
        cmap = plt.get_cmap(colormap)
        rgba_obs = cmap(expr_norm).astype(np.float32)  # shape (N_cells, 4)

        # Set alpha: transparent for zero-expression cells
        zero_mask = expr == 0
        rgba_obs[zero_mask, 3] = zero_alpha
        rgba_obs[~zero_mask, 3] = 1.0

        # Build label-indexed array
        color_arr = self._empty_color_array()
        valid_mask = self.label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = self.label_to_obs[valid_labels]
        color_arr[valid_labels] = rgba_obs[obs_indices]

        self._gene_cache[cache_key] = color_arr
        return color_arr

    def get_continuous_colors(
        self,
        values: np.ndarray,
        colormap: str = "viridis",
        cache_key: Optional[str] = None,
    ) -> np.ndarray:
        """
        Build RGBA color array for an arbitrary continuous per-cell score
        (e.g. a CNV burden score) that isn't a gene in ``adata.var_names``.

        Parameters
        ----------
        values : np.ndarray
            Length adata.n_obs, in adata.obs row order (not gene expression).
        colormap : str
            Matplotlib colormap name.
        cache_key : str, optional
            If given, cache the result under this key.

        Returns
        -------
        np.ndarray, shape (max_label + 1, 4), dtype float32
        """
        if cache_key is not None and cache_key in self._continuous_cache:
            return self._continuous_cache[cache_key]

        expr = np.asarray(values, dtype=np.float32).ravel()
        if expr.shape[0] != self.adata.n_obs:
            raise ValueError(
                f"values length {expr.shape[0]} does not match adata.n_obs {self.adata.n_obs}"
            )

        # Normalise over all finite values to [0, 1] — unlike gene expression,
        # a score of 0 is meaningful here, so we don't restrict to nonzero.
        finite = np.isfinite(expr)
        expr_norm = np.zeros_like(expr)
        if finite.any():
            vmin = expr[finite].min()
            vmax = expr[finite].max()
            if vmax > vmin:
                expr_norm[finite] = (expr[finite] - vmin) / (vmax - vmin)
            else:
                expr_norm[finite] = 1.0

        # Map via colormap
        cmap = plt.get_cmap(colormap)
        rgba_obs = cmap(expr_norm).astype(np.float32)

        # Alpha: transparent for non-finite (missing) values
        rgba_obs[:, 3] = 1.0
        rgba_obs[~finite, 3] = 0.0

        # Build label-indexed array
        color_arr = self._empty_color_array()
        valid_mask = self.label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = self.label_to_obs[valid_labels]
        color_arr[valid_labels] = rgba_obs[obs_indices]

        if cache_key is not None:
            self._continuous_cache[cache_key] = color_arr
        return color_arr

    def get_gene_colors_filtered(
        self,
        gene_name: str,
        cluster_series: pd.Series,
        cluster_id: int,
        colormap: str = "viridis",
    ) -> np.ndarray:
        """
        Build RGBA color array for gene expression, filtered to a single cluster.

        Cells NOT in the selected cluster get alpha=0 (transparent).

        Parameters
        ----------
        gene_name : str
        cluster_series : pd.Series
            Cluster assignments indexed by cell barcode.
        cluster_id : int
            Cluster ID to keep visible.
        colormap : str

        Returns
        -------
        np.ndarray, shape (max_label + 1, 4), dtype float32
        """
        # Align cluster assignments to adata obs
        if 'cell_id' in self.adata.obs.columns:
            cell_ids = self.adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids)
        else:
            clusters_aligned = cluster_series.reindex(self.adata.obs_names)
        # Handle both integer and string cluster IDs
        try:
            filled = clusters_aligned.fillna(-1)
            cluster_values = filled.values.astype(np.int32)
            in_cluster = cluster_values == cluster_id
        except (ValueError, TypeError):
            cluster_values = clusters_aligned.values.astype(str)
            in_cluster = cluster_values == str(cluster_id)

        # Extract expression for this gene (sparse → dense)
        gene_idx = self.adata.var_names.get_loc(gene_name)
        X = self.adata.X
        if hasattr(X, "toarray"):
            expr = np.asarray(X[:, gene_idx].toarray()).ravel().astype(np.float32)
        else:
            expr = np.asarray(X[:, gene_idx]).ravel().astype(np.float32)

        # Normalise using min/max of non-zero expression WITHIN the cluster
        nonzero_in_cluster = in_cluster & (expr > 0)
        expr_norm = np.zeros_like(expr)
        if nonzero_in_cluster.any():
            vmin = expr[nonzero_in_cluster].min()
            vmax = expr[nonzero_in_cluster].max()
            if vmax > vmin:
                expr_norm[nonzero_in_cluster] = (expr[nonzero_in_cluster] - vmin) / (vmax - vmin)
            else:
                expr_norm[nonzero_in_cluster] = 1.0

        # Map via colormap
        cmap = plt.get_cmap(colormap)
        rgba_obs = cmap(expr_norm).astype(np.float32)

        # Alpha: transparent for zero-expression and out-of-cluster cells
        rgba_obs[:, 3] = 0.0
        rgba_obs[nonzero_in_cluster, 3] = 1.0

        # Build label-indexed array
        color_arr = self._empty_color_array()
        valid_mask = self.label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = self.label_to_obs[valid_labels]
        color_arr[valid_labels] = rgba_obs[obs_indices]

        return color_arr

    def get_cluster_colors(
        self,
        cluster_series: pd.Series,
        palette: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, dict]:
        """
        Build RGBA color array for a cluster assignment.

        Parameters
        ----------
        cluster_series : pd.Series
            Cluster assignments indexed by cell barcode.
            Values should be integer cluster IDs (1-based is fine).
        palette : np.ndarray, optional
            Shape (K, 4) RGBA palette. Defaults to CLUSTER_PALETTE.

        Returns
        -------
        color_arr : np.ndarray, shape (max_label + 1, 4)
        cluster_to_color : dict mapping cluster_id -> RGBA tuple
        """
        cache_key = cluster_series.name or "unnamed"
        if cache_key in self._cluster_cache:
            return self._cluster_cache[cache_key]

        if palette is None:
            palette = CLUSTER_PALETTE

        # Align cluster_series with adata obs via cell_id column.
        # Cluster CSVs are indexed by cell barcode (e.g. 'aaaagflk-1'),
        # but adata.obs_names are integer indices ('0', '1', ...).
        # Use adata.obs['cell_id'] as the join key.
        if 'cell_id' in self.adata.obs.columns:
            cell_ids = self.adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids)
        else:
            clusters_aligned = cluster_series.reindex(self.adata.obs_names)

        # Handle both integer and string cluster IDs
        try:
            filled = clusters_aligned.fillna(-1)
            cluster_values = filled.values.astype(np.int32)
        except (ValueError, TypeError):
            import pandas as pd
            # String-valued clusters: factorize to integers
            raw = clusters_aligned.values
            codes, uniques = pd.factorize(raw)
            cluster_values = codes.astype(np.int32)  # -1 for NaN

        unique_clusters = sorted(set(cluster_values[cluster_values >= 0]))
        cluster_to_color = {}
        for i, c in enumerate(unique_clusters):
            cluster_to_color[c] = palette[i % len(palette)]

        # Build RGBA array per obs
        rgba_obs = np.zeros((len(self.adata), 4), dtype=np.float32)
        for obs_idx, c in enumerate(cluster_values):
            if c >= 0:
                rgba_obs[obs_idx] = cluster_to_color[c]

        # Build label-indexed array
        color_arr = self._empty_color_array()
        valid_mask = self.label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = self.label_to_obs[valid_labels]
        color_arr[valid_labels] = rgba_obs[obs_indices]

        self._cluster_cache[cache_key] = (color_arr, cluster_to_color)
        return color_arr, cluster_to_color

    def build_direct_label_colormap(self, color_arr: np.ndarray):
        """
        Convert an RGBA array to a napari DirectLabelColormap.

        Parameters
        ----------
        color_arr : np.ndarray, shape (N+1, 4), dtype float32

        Returns
        -------
        napari.utils.colormaps.label_colormap.DirectLabelColormap
        """
        from napari.utils.colormaps import DirectLabelColormap

        # Build color_dict from nonzero-alpha entries only (performance)
        nonzero = np.where(color_arr[:, 3] > 0)[0]
        labels_py = nonzero.tolist()                    # bulk C-level int conversion
        colors_py = color_arr[nonzero].tolist()         # bulk C-level float conversion
        color_dict = {label: tuple(color) for label, color in zip(labels_py, colors_py)}
        # Background and unlabelled cells → transparent
        color_dict[None] = (0.0, 0.0, 0.0, 0.0)

        return DirectLabelColormap(color_dict=color_dict)

    def apply_to_labels_layer(self, labels_layer, color_arr: np.ndarray):
        """
        Apply an RGBA color array to a napari Labels layer.

        Must be called from the main Qt thread.

        Parameters
        ----------
        labels_layer : napari.layers.Labels
        color_arr : np.ndarray, shape (max_label + 1, 4)
        """
        colormap = self.build_direct_label_colormap(color_arr)
        labels_layer.colormap = colormap
        labels_layer.refresh()
