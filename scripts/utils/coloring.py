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

# Categorical palettes for clusters (max 20 clusters)
CLUSTER_PALETTE = np.array([
    [0.12, 0.47, 0.71, 1.0],  # blue
    [1.00, 0.50, 0.05, 1.0],  # orange
    [0.17, 0.63, 0.17, 1.0],  # green
    [0.84, 0.15, 0.16, 1.0],  # red
    [0.58, 0.40, 0.74, 1.0],  # purple
    [0.55, 0.34, 0.29, 1.0],  # brown
    [0.89, 0.47, 0.76, 1.0],  # pink
    [0.50, 0.50, 0.50, 1.0],  # gray
    [0.74, 0.74, 0.13, 1.0],  # olive
    [0.09, 0.75, 0.81, 1.0],  # cyan
    [0.65, 0.81, 0.89, 1.0],  # light blue
    [0.98, 0.60, 0.60, 1.0],  # light red
    [0.70, 0.87, 0.54, 1.0],  # light green
    [0.99, 0.75, 0.44, 1.0],  # light orange
    [0.79, 0.70, 0.84, 1.0],  # light purple
    [0.70, 0.49, 0.38, 1.0],  # light brown
    [0.99, 0.71, 0.94, 1.0],  # light pink
    [0.90, 0.90, 0.90, 1.0],  # light gray
    [0.90, 0.90, 0.60, 1.0],  # light olive
    [0.73, 0.93, 0.96, 1.0],  # light cyan
], dtype=np.float32)

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
        self._max_label = len(label_to_obs) - 1

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

        # Normalise to [0, 1]
        vmax = expr.max()
        if vmax > 0:
            expr_norm = expr / vmax
        else:
            expr_norm = expr.copy()

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
        # Get the full gene coloring (uses cache)
        color_arr = self.get_gene_colors(gene_name, colormap=colormap).copy()

        # Align cluster assignments to adata obs (same pattern as get_cluster_colors)
        if 'cell_id' in self.adata.obs.columns:
            cell_ids = self.adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids, fill_value=-1)
        else:
            clusters_aligned = cluster_series.reindex(self.adata.obs_names, fill_value=-1)
        cluster_values = clusters_aligned.values.astype(np.int32)

        # Build mask of obs indices NOT in the selected cluster
        not_in_cluster = cluster_values != cluster_id

        # Set alpha=0 for labels whose obs is not in the selected cluster
        valid_mask = self.label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = self.label_to_obs[valid_labels]
        labels_to_clear = valid_labels[not_in_cluster[obs_indices]]
        color_arr[labels_to_clear, 3] = 0.0

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
            clusters_aligned = cluster_series.reindex(cell_ids, fill_value=-1)
        else:
            clusters_aligned = cluster_series.reindex(self.adata.obs_names, fill_value=-1)
        cluster_values = clusters_aligned.values.astype(np.int32)

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
        color_dict = {
            int(label): tuple(color_arr[label].tolist())
            for label in nonzero
        }
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
