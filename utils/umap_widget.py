"""
Interactive UMAP scatter plot as a second napari viewer window.

Uses a napari Points layer for the UMAP, giving native pan/zoom and
hover-to-inspect. When cells are colored by cluster, hovering over a
point in the UMAP shows the cluster ID in the status bar.

Usage:
    from utils.umap_widget import UMAPViewer
    umap_viewer = UMAPViewer(umap_df, cell_ids)
    umap_viewer.color_by_cluster("graphclust", color_arr, label_to_obs,
                                  cluster_ids_per_obs)
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import napari


class UMAPViewer:
    """
    UMAP scatter plot in a dedicated napari Viewer window.

    Points are displayed as a napari Points layer. When cluster coloring
    is applied, hovering over a point shows cluster ID + cell ID in the
    status bar.

    Parameters
    ----------
    umap_df : pd.DataFrame
        UMAP coordinates with columns ['UMAP_1', 'UMAP_2'], indexed by cell barcode.
    cell_ids : np.ndarray
        Cell IDs from sdata["table"].obs['cell_id'].
    """

    def __init__(
        self,
        umap_df: pd.DataFrame,
        cell_ids: np.ndarray,
    ):
        self.n_cells = len(cell_ids)
        self._cell_ids = cell_ids

        # Align UMAP to adata obs via cell_id (reindex handles the 91 missing cells)
        aligned = umap_df.reindex(cell_ids)
        self._xy = aligned[["UMAP_1", "UMAP_2"]].values.astype(np.float32)
        self._valid = ~np.isnan(self._xy[:, 0])  # mask of cells in UMAP
        self._valid_indices = np.where(self._valid)[0]

        # Coordinates for napari: (y, x) ordering, only valid cells
        xy_valid = self._xy[self._valid]
        self._points_data = np.column_stack([xy_valid[:, 1], xy_valid[:, 0]])

        # Current state
        self._current_colors: Optional[np.ndarray] = None
        self._cluster_ids: Optional[np.ndarray] = None  # per valid-point cluster IDs
        self._clustering_name: Optional[str] = None

        # Create the UMAP viewer window (deferred — created on first show())
        self._viewer: Optional[napari.Viewer] = None
        self._points_layer = None
        self._pending_title: Optional[str] = None

    def _viewer_is_alive(self) -> bool:
        """Check if the UMAP viewer window still exists."""
        if self._viewer is None:
            return False
        try:
            # Access the Qt window — raises RuntimeError if deleted
            self._viewer.window._qt_window.isVisible()
            return True
        except (RuntimeError, AttributeError):
            # Qt C++ object has been deleted
            self._viewer = None
            self._points_layer = None
            return False

    def _ensure_viewer(self):
        """Create (or recreate) the napari UMAP viewer."""
        if self._viewer_is_alive():
            return

        self._viewer = napari.Viewer(
            title="UMAP",
            ndisplay=2,
            show=True,
        )

        # Initial gray points
        n_valid = self._valid.sum()
        default_colors = np.full((n_valid, 4), [0.7, 0.7, 0.7, 0.5], dtype=np.float32)

        self._points_layer = self._viewer.add_points(
            self._points_data,
            name="UMAP cells",
            size=0.15,
            face_color=default_colors,
            border_color="transparent",
            opacity=1.0,
        )

        # Hide axes ticks/grid for a cleaner look
        self._viewer.axes.visible = False

        # Reset view to fit all points
        self._viewer.reset_view()

        # Register hover callback — use cursor.events.position so our status
        # text overwrites napari's default layer status updates
        self._viewer.cursor.events.position.connect(
            lambda event: self._handle_hover()
        )

        # Reapply colors if we had them before the window was closed
        if self._current_colors is not None:
            rgba_valid = self._current_colors[self._valid]
            self._points_layer.face_color = rgba_valid

        # Apply pending title if set
        if self._pending_title is not None:
            self._viewer.title = self._pending_title

    def show(self):
        """Open (or bring to front) the UMAP viewer window."""
        self._ensure_viewer()
        self._viewer.window.activate()

    def set_point_size(self, size: float):
        """Set the UMAP point size."""
        if self._points_layer is not None:
            self._points_layer.size = size

    # ── Public coloring API ──────────────────────────────────────────────────

    def color_by_gene(
        self,
        gene_name: str,
        color_arr: np.ndarray,
        label_to_obs: np.ndarray,
    ):
        """
        Color UMAP points using the same RGBA array as the Labels layer.

        Does NOT auto-open the UMAP window. Colors are stored and applied
        when the user manually opens the window via show().

        Parameters
        ----------
        gene_name : str
        color_arr : np.ndarray, shape (max_label + 1, 4)
        label_to_obs : np.ndarray
        """
        rgba_obs = self._labels_color_arr_to_obs_colors(color_arr, label_to_obs)
        self._current_colors = rgba_obs

        # Clear cluster hover info
        self._cluster_ids = None
        self._clustering_name = None

        self._pending_title = f"UMAP — {gene_name}"

        # Update viewer only if already open
        if self._viewer_is_alive():
            rgba_valid = rgba_obs[self._valid]
            self._points_layer.face_color = rgba_valid
            self._viewer.title = self._pending_title

    def color_by_cluster(
        self,
        clustering_name: str,
        color_arr: np.ndarray,
        label_to_obs: np.ndarray,
        cluster_ids_per_obs: Optional[np.ndarray] = None,
    ):
        """
        Color UMAP points by cluster assignment.

        Does NOT auto-open the UMAP window. Colors are stored and applied
        when the user manually opens the window via show().

        Parameters
        ----------
        clustering_name : str
        color_arr : np.ndarray, shape (max_label + 1, 4)
        label_to_obs : np.ndarray
        cluster_ids_per_obs : np.ndarray, optional
            Shape (n_obs,), cluster ID for each obs row. If provided,
            enables hover-to-see-cluster-ID.
        """
        rgba_obs = self._labels_color_arr_to_obs_colors(color_arr, label_to_obs)
        self._current_colors = rgba_obs

        # Store cluster IDs for hover
        if cluster_ids_per_obs is not None and len(cluster_ids_per_obs) == self.n_cells:
            self._cluster_ids = cluster_ids_per_obs[self._valid]
        else:
            self._cluster_ids = None
        self._clustering_name = clustering_name

        self._pending_title = f"UMAP — {clustering_name}"

        # Update viewer only if already open
        if self._viewer_is_alive():
            rgba_valid = rgba_obs[self._valid]
            self._points_layer.face_color = rgba_valid
            self._viewer.title = self._pending_title

    # ── Hover handler ────────────────────────────────────────────────────────

    def _handle_hover(self):
        """Show cluster ID and cell ID in the UMAP viewer status bar on hover."""
        from qtpy.QtCore import QTimer

        if self._points_layer is None or self._viewer is None:
            return

        # Get the value under the cursor (point index or None)
        val = self._points_layer.get_value(
            self._viewer.cursor.position,
            view_direction=None,
            dims_displayed=list(range(self._points_layer.ndim)),
            world=True,
        )

        # Points layer returns (point_index, None) or (None, None)
        if val is None:
            return
        point_idx = val[0] if isinstance(val, tuple) else val
        if point_idx is None:
            return

        # Map point index → obs index → cell_id
        obs_idx = self._valid_indices[point_idx]
        cell_id = self._cell_ids[obs_idx]

        parts = [f"Cell: {cell_id}"]
        if self._cluster_ids is not None and self._clustering_name:
            cid = self._cluster_ids[point_idx]
            parts.append(f"{self._clustering_name}: {cid}")

        text = " | ".join(parts)
        viewer = self._viewer
        QTimer.singleShot(0, lambda: setattr(viewer, 'status', text))

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
        if self._points_layer is None or self._current_colors is None:
            return

        colors = self._current_colors.copy()
        if dim_others:
            colors[:, 3] *= 0.1  # dim all
            colors[obs_indices, 3] = 1.0  # restore selected

        rgba_valid = colors[self._valid]
        self._points_layer.face_color = rgba_valid

    def reset_highlights(self):
        """Restore all cells to full opacity."""
        if self._current_colors is not None and self._points_layer is not None:
            rgba_valid = self._current_colors[self._valid]
            self._points_layer.face_color = rgba_valid

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
        n_obs = self.n_cells
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
