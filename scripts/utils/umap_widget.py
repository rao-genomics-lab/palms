"""
Linked UMAP scatter plot embedded as a napari dock widget.

Uses matplotlib's FigureCanvasQTAgg to render directly inside napari's
Qt window — no separate window needed, no ICE/X11 conflicts.

Usage:
    from utils.umap_widget import UMAPWidget
    widget = UMAPWidget(umap_df, adata_obs_names)
    viewer.window.add_dock_widget(widget, name="UMAP", area="bottom")
    widget.color_by_gene("MMSET", color_arr, label_to_obs)
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from qtpy.QtWidgets import QWidget, QVBoxLayout


class UMAPWidget(QWidget):
    """
    UMAP scatter plot as a Qt widget (embeddable in napari).

    Parameters
    ----------
    umap_df : pd.DataFrame
        UMAP coordinates with columns ['UMAP_1', 'UMAP_2'], indexed by cell barcode.
    adata_obs_names : pd.Index
        Cell barcodes from sdata["table"].obs_names (may have 91 more than UMAP).
    """

    def __init__(
        self,
        umap_df: pd.DataFrame,
        adata_obs_names: pd.Index,
        parent=None,
    ):
        super().__init__(parent)
        self.adata_obs_names = adata_obs_names

        # Align UMAP to adata obs (reindex handles the 91 missing cells)
        aligned = umap_df.reindex(adata_obs_names)
        self._xy = aligned[["UMAP_1", "UMAP_2"]].values.astype(np.float32)
        self._valid = ~np.isnan(self._xy[:, 0])  # mask of cells in UMAP

        self._current_colors: Optional[np.ndarray] = None

        # Build the embedded matplotlib figure
        self._fig = Figure(figsize=(6, 5), dpi=100)
        self._canvas = FigureCanvasQTAgg(self._fig)
        self._ax = self._fig.add_subplot(111)

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._canvas)
        self.setLayout(layout)

        self._build_scatter()

    def _build_scatter(self):
        """Draw the initial gray scatter."""
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
        self._ax.set_xlabel("UMAP 1")
        self._ax.set_ylabel("UMAP 2")
        self._ax.set_title("UMAP (gene expression)")
        self._ax.set_aspect("equal")
        self._fig.tight_layout()
        self._canvas.draw_idle()

    # ── Public coloring API ──────────────────────────────────────────────────

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
        label_to_obs : np.ndarray
        """
        rgba_obs = self._labels_color_arr_to_obs_colors(color_arr, label_to_obs)
        self._update_scatter_colors(rgba_obs)
        self._ax.set_title(f"UMAP — {gene_name}")
        self._canvas.draw_idle()

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
        self._canvas.draw_idle()

    # ── Selection highlighting ───────────────────────────────────────────────

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
        self._canvas.draw_idle()

    def reset_highlights(self):
        """Restore all cells to full opacity."""
        if self._current_colors is not None:
            self._update_scatter_colors(self._current_colors)
            self._canvas.draw_idle()

    # ── Internal helpers ─────────────────────────────────────────────────────

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
        rgba_valid = rgba_obs[self._valid]
        self._scatter.set_facecolor(rgba_valid)
        self._current_colors = rgba_obs
