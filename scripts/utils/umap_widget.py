"""
Linked matplotlib UMAP scatter window for the Xenium viewer.

Opens as a separate matplotlib figure (not embedded in napari), showing all
318K cells as a scatter plot. Colors sync with the napari Labels layer.

Usage:
    from utils.umap_widget import UMAPWindow
    umap_win = UMAPWindow(umap_df)
    umap_win.show()
    umap_win.color_by_gene("MMSET", color_arr, label_to_obs)
"""

from __future__ import annotations

from typing import Optional, Callable

import numpy as np
import pandas as pd
# matplotlib.pyplot imported lazily in _build_figure() so napari's Qt
# is already running before any matplotlib Qt canvas is created.
from matplotlib.figure import Figure
from matplotlib.axes import Axes


class UMAPWindow:
    """
    Interactive UMAP scatter plot linked to the napari Labels layer.

    Parameters
    ----------
    umap_df : pd.DataFrame
        UMAP coordinates with columns ['UMAP_1', 'UMAP_2'], indexed by cell barcode.
    adata_obs_names : pd.Index
        Cell barcodes from sdata["table"].obs_names (may have 91 more than UMAP).
    on_cluster_selected : Callable, optional
        Callback when user clicks a cluster point in the UMAP.
        Called with (obs_indices: np.ndarray).
    """

    def __init__(
        self,
        umap_df: pd.DataFrame,
        adata_obs_names: pd.Index,
        on_cluster_selected: Optional[Callable] = None,
    ):
        self.umap_df = umap_df
        self.adata_obs_names = adata_obs_names
        self.on_cluster_selected = on_cluster_selected

        # Align UMAP to adata obs (reindex handles the 91 missing cells)
        self._aligned_umap = umap_df.reindex(adata_obs_names)
        self._xy = self._aligned_umap[["UMAP_1", "UMAP_2"]].values.astype(np.float32)
        self._valid = ~np.isnan(self._xy[:, 0])  # mask of cells in UMAP

        self._fig: Optional[Figure] = None
        self._ax: Optional[Axes] = None
        self._scatter = None
        self._current_colors: Optional[np.ndarray] = None

    def show(self):
        """Open or bring the UMAP window to focus."""
        import matplotlib.pyplot as plt
        if self._fig is None or not plt.fignum_exists(self._fig.number):
            self._build_figure()
        self._fig.canvas.manager.window.activateWindow()
        self._fig.canvas.manager.window.raise_()
        self._fig.canvas.draw_idle()

    def _build_figure(self):
        import matplotlib.pyplot as plt
        self._fig, self._ax = plt.subplots(figsize=(8, 6))
        self._fig.canvas.manager.set_window_title("Xenium Viewer — UMAP")
        self._ax.set_xlabel("UMAP 1")
        self._ax.set_ylabel("UMAP 2")
        self._ax.set_title("UMAP (gene expression)")
        self._ax.set_aspect("equal")

        # Default: all cells gray
        n_valid = self._valid.sum()
        default_colors = np.full((n_valid, 4), [0.7, 0.7, 0.7, 0.5], dtype=np.float32)

        xy_valid = self._xy[self._valid]
        self._scatter = self._ax.scatter(
            xy_valid[:, 0],
            xy_valid[:, 1],
            c=default_colors,
            s=1,
            rasterized=True,
            linewidths=0,
        )
        self._fig.tight_layout()
        self._fig.canvas.draw_idle()
        self._fig.show()

    def color_by_gene(
        self,
        gene_name: str,
        color_arr: np.ndarray,
        label_to_obs: np.ndarray,
    ):
        """
        Color UMAP points using the same RGBA array as the Labels layer.

        Parameters
        ----------
        gene_name : str
        color_arr : np.ndarray, shape (max_label + 1, 4)
            From CellColorManager.get_gene_colors()
        label_to_obs : np.ndarray
            label → obs index mapping
        """
        rgba_obs = self._labels_color_arr_to_obs_colors(color_arr, label_to_obs)
        self._update_scatter_colors(rgba_obs)
        self._ax.set_title(f"UMAP — {gene_name}")
        self._fig.canvas.draw_idle()

    def color_by_cluster(
        self,
        clustering_name: str,
        color_arr: np.ndarray,
        label_to_obs: np.ndarray,
    ):
        """Color UMAP points by cluster assignment."""
        rgba_obs = self._labels_color_arr_to_obs_colors(color_arr, label_to_obs)
        self._update_scatter_colors(rgba_obs)
        self._ax.set_title(f"UMAP — {clustering_name}")
        self._fig.canvas.draw_idle()

    def _labels_color_arr_to_obs_colors(
        self,
        color_arr: np.ndarray,
        label_to_obs: np.ndarray,
    ) -> np.ndarray:
        """
        Convert a label-indexed color array to an obs-indexed color array.

        color_arr shape: (max_label + 1, 4)
        Returns rgba_obs shape: (N_cells, 4)
        """
        n_obs = len(self.adata_obs_names)
        rgba_obs = np.zeros((n_obs, 4), dtype=np.float32)
        rgba_obs[:, 3] = 0.8  # default alpha

        # Invert the label_to_obs mapping: obs_idx -> label
        # label_to_obs[label] = obs_idx, so we need obs→label
        max_label = len(label_to_obs) - 1
        obs_to_label = np.full(n_obs, -1, dtype=np.int32)
        for label_val in range(1, max_label + 1):
            obs_idx = label_to_obs[label_val]
            if 0 <= obs_idx < n_obs:
                obs_to_label[obs_idx] = label_val

        valid_obs = obs_to_label >= 0
        rgba_obs[valid_obs] = color_arr[obs_to_label[valid_obs]]

        return rgba_obs

    def _update_scatter_colors(self, rgba_obs: np.ndarray):
        """Update scatter plot colors (only for cells present in UMAP)."""
        import matplotlib.pyplot as plt
        if self._fig is None or not plt.fignum_exists(self._fig.number):
            self._build_figure()

        rgba_valid = rgba_obs[self._valid]
        self._scatter.set_facecolor(rgba_valid)
        self._current_colors = rgba_obs

    def highlight_obs_indices(self, obs_indices: np.ndarray, dim_others: bool = True):
        """
        Highlight specific cells in the UMAP (e.g. from a napari selection).

        Parameters
        ----------
        obs_indices : np.ndarray
            Row indices into adata.obs to highlight.
        dim_others : bool
            If True, reduce alpha of non-selected cells.
        """
        if self._scatter is None or self._current_colors is None:
            return

        colors = self._current_colors.copy()
        if dim_others:
            colors[:, 3] *= 0.1  # dim all
            colors[obs_indices, 3] = 1.0  # restore selected

        self._update_scatter_colors(colors)
        self._fig.canvas.draw_idle()

    def reset_highlights(self):
        """Restore all cells to full opacity."""
        if self._current_colors is not None:
            self._update_scatter_colors(self._current_colors)
            self._fig.canvas.draw_idle()
