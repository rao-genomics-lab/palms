"""Tab: Custom Cell Segmentation.

Lets the user swap the native Xenium cell segmentation for a custom one
produced by the two-stage preprocessing pipeline:

  Stage 1: Rscript scripts/extract_seurat_segmentation.R <rds> <stage1_dir>
  Stage 2: python  scripts/build_custom_segmentation.py  <xenium_dir> <stage1_dir>

The user selects the resulting `custom_segmentation.h5ad` file; the companion
`custom_labels.zarr` is found automatically from the metadata JSON in the same
directory.  All downstream analysis (coloring, ROI, clustering, spatial) then
operates on the custom cells.

Once loaded, the custom segmentation is cached inside the SpatialData zarr store
(sdata.labels["custom_cell_labels"] + sdata.tables["custom_table"]) so it can
be restored automatically on the next launch without selecting the h5ad again.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import PushButton
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QFileDialog, QScrollArea,
)
from napari.qt.threading import thread_worker
from palms.tabs._helpers import StatusProxy

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    status = StatusProxy(ctx.viewer)

    # ── Status display ────────────────────────────────────────────────────────
    active_label = QLabel("Active segmentation: Xenium (native)")
    active_label.setWordWrap(True)

    info_label = QLabel(
        "To use custom segmentation:\n"
        "  1. Run extract_seurat_segmentation.R\n"
        "  2. Run build_custom_segmentation.py\n"
        "  3. Click 'Load Custom Segmentation...'"
    )
    info_label.setWordWrap(True)

    # ── Buttons ───────────────────────────────────────────────────────────────
    load_btn = PushButton(label="Load Custom Segmentation...", enabled=True)
    revert_btn = PushButton(label="Revert to Xenium Segmentation", enabled=False)
    update_sdata_btn = PushButton(label="Update SpatialData on disk", enabled=True)

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def _on_load():
        from palms.utils.adata_persistence import (
            load_custom_seg_from_sdata,
            CUSTOM_LABELS_KEY,
            CUSTOM_TABLE_KEY,
        )

        # Pre-check: if custom seg is already cached in sdata, offer to load it
        if (
            ctx.sdata is not None
            and CUSTOM_LABELS_KEY in ctx.sdata.labels
            and CUSTOM_TABLE_KEY in ctx.sdata.tables
        ):
            from qtpy.QtWidgets import QMessageBox
            reply = QMessageBox.question(
                None,
                "Cached segmentation available",
                "A custom segmentation is cached in SpatialData.\n"
                "Load from cache (fast) or select a new h5ad file?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes,
            )
            if reply == QMessageBox.StandardButton.Yes:
                adata, scales = load_custom_seg_from_sdata(ctx.sdata)
                if adata is not None:
                    try:
                        _apply_custom_segmentation(ctx, adata, scales)
                        active_label.setText(
                            f"Active segmentation: Custom (from cache)\n"
                            f"  {adata.n_obs:,} cells × {adata.n_vars} genes"
                        )
                        revert_btn.enabled = True
                        status.value = (
                            f"Custom segmentation loaded from cache: {adata.n_obs:,} cells"
                        )
                    except Exception as exc:
                        status.value = f"ERROR loading from cache: {exc}"
                        import traceback
                        traceback.print_exc()
                    return
                # Cache load failed — fall through to file dialog

        h5ad_path_str, _ = QFileDialog.getOpenFileName(
            None,
            "Select custom_segmentation.h5ad",
            str(ctx.data_path) if ctx.data_path else "",
            "AnnData files (*.h5ad)",
        )
        if not h5ad_path_str:
            return

        h5ad_path = Path(h5ad_path_str)

        # Locate zarr via JSON metadata in same directory
        json_path = h5ad_path.parent / "custom_segmentation.json"
        if json_path.exists():
            meta = json.loads(json_path.read_text())
            zarr_name = meta.get("label_raster", "custom_labels.zarr")
            zarr_path = h5ad_path.parent / zarr_name
        else:
            zarr_path = h5ad_path.parent / "custom_labels.zarr"

        if not zarr_path.exists():
            status.value = f"ERROR: {zarr_path} not found alongside h5ad"
            return

        # Validate h5ad
        try:
            import anndata
            _adata = anndata.read_h5ad(str(h5ad_path), backed="r")
            if "cell_id" not in _adata.obs.columns:
                status.value = "ERROR: h5ad must have obs['cell_id'] (integer label)"
                _adata.file.close()
                return
            if "spatial" not in _adata.obsm:
                status.value = "ERROR: h5ad must have obsm['spatial'] (xy µm)"
                _adata.file.close()
                return
            _adata.file.close()
        except Exception as exc:
            status.value = f"ERROR loading h5ad: {exc}"
            return

        load_btn.enabled = False
        status.value = "Loading custom segmentation..."

        @thread_worker
        def _run():
            import anndata
            import zarr
            new_adata = anndata.read_h5ad(str(h5ad_path))
            z = zarr.open(str(zarr_path), mode="r")
            scale_keys = sorted(
                [k for k in z.keys() if k.isdigit()],
                key=lambda s: int(s),
            )
            if not scale_keys:
                raise RuntimeError(f"No numeric scale keys found in {zarr_path}")
            scales = [np.asarray(z[k]) for k in scale_keys]
            return new_adata, scales

        def _on_done(result):
            new_adata, scales = result
            try:
                _apply_custom_segmentation(ctx, new_adata, scales)
                active_label.setText(
                    f"Active segmentation: Custom ({h5ad_path.name})\n"
                    f"  {new_adata.n_obs:,} cells × {new_adata.n_vars} genes"
                )
                revert_btn.enabled = True
                status.value = "Caching custom segmentation to SpatialData..."

                # Asynchronously cache to sdata
                from palms.utils.adata_persistence import save_custom_seg_to_sdata

                @thread_worker
                def _save():
                    # A different segmentation than whatever is cached, and
                    # _apply_custom_segmentation already removed the old layer.
                    save_custom_seg_to_sdata(
                        ctx, new_adata, scales, replace_backed=True,
                    )

                def _save_done(_):
                    status.value = (
                        f"Custom segmentation loaded and cached: {new_adata.n_obs:,} cells"
                    )

                def _save_error(exc):
                    # napari's `errored` emits the exception itself, not an
                    # exc_info triple — indexing it raised TypeError and masked
                    # whatever actually went wrong.
                    status.value = f"Custom segmentation loaded (cache failed: {exc})"
                    from palms.utils.reporting import report_write_failure
                    report_write_failure(exc, "custom segmentation")

                save_worker = _save()
                save_worker.returned.connect(_save_done)
                save_worker.errored.connect(_save_error)
                ctx.state["_active_save_worker"] = save_worker  # prevent GC
                save_worker.start()

            except Exception as exc:
                status.value = f"ERROR swapping segmentation: {exc}"
                import traceback
                traceback.print_exc()
            finally:
                load_btn.enabled = True

        worker = _run()
        worker.returned.connect(_on_done)
        worker.start()

    def _on_revert():
        revert_btn.enabled = False
        load_btn.enabled = False
        status.value = "Reverting to native Xenium segmentation..."
        try:
            _revert_xenium_segmentation(ctx)
            active_label.setText("Active segmentation: Xenium (native)")
            status.value = "Reverted to Xenium native segmentation"
        except Exception as exc:
            status.value = f"ERROR reverting: {exc}"
            import traceback
            traceback.print_exc()
        finally:
            load_btn.enabled = True

    def _on_update_sdata():
        update_sdata_btn.enabled = False
        status.value = "Saving to SpatialData..."

        @thread_worker
        def _run():
            if ctx.segmentation_source == "custom":
                from palms.utils.adata_persistence import save_custom_seg_to_sdata
                scales = list(ctx.cell_labels_layer.data)
                save_custom_seg_to_sdata(ctx, ctx.adata, scales)
            else:
                from palms.utils.adata_persistence import _persist_table
                _persist_table(ctx)

        def _done(_):
            update_sdata_btn.enabled = True
            status.value = "SpatialData updated on disk."

        def _error(exc):
            update_sdata_btn.enabled = True
            status.value = f"ERROR updating SpatialData: {exc}"

        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_error)
        w.start()

    # ── Connect ───────────────────────────────────────────────────────────────
    load_btn.changed.connect(_on_load)
    revert_btn.changed.connect(_on_revert)
    update_sdata_btn.changed.connect(_on_update_sdata)

    # ── Session restore ────────────────────────────────────────────────────────
    def _restore_session(session):
        if session.get("segmentation_source") != "custom":
            return
        from palms.utils.adata_persistence import load_custom_seg_from_sdata
        adata, scales = load_custom_seg_from_sdata(ctx.sdata)
        if adata is None:
            active_label.setText(
                "Active segmentation: Custom (cache not found — please reload)"
            )
            return
        try:
            _apply_custom_segmentation(ctx, adata, scales)
            active_label.setText(
                f"Active segmentation: Custom (restored from cache)\n"
                f"  {adata.n_obs:,} cells × {adata.n_vars} genes"
            )
            revert_btn.enabled = True
        except Exception as exc:
            active_label.setText(f"Active segmentation: Restore failed — {exc}")
            import traceback
            traceback.print_exc()

    # ── Layout ────────────────────────────────────────────────────────────────
    container = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)
    layout.addWidget(active_label)
    layout.addWidget(info_label)
    layout.addWidget(load_btn.native)
    layout.addWidget(revert_btn.native)
    layout.addWidget(update_sdata_btn.native)
    layout.addStretch()
    container.setLayout(layout)

    scroll = QScrollArea()
    scroll.setWidget(container)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.Shape.NoFrame)

    return scroll, {"restore_session": _restore_session}


# ── Segmentation swap helpers ─────────────────────────────────────────────────

def _build_label_to_obs(adata) -> np.ndarray:
    """Build label_value → obs row index mapping from adata.obs['cell_id']."""
    if "cell_id" in adata.obs.columns:
        label_values = adata.obs["cell_id"].values.astype(np.int32)
    else:
        # Assume sequential 1..N
        n = len(adata.obs)
        lto = np.arange(-1, n, dtype=np.int32)
        return lto

    max_label = int(label_values.max())
    lto = np.full(max_label + 1, -1, dtype=np.int32)
    for obs_idx, lv in enumerate(label_values):
        if lv > 0:
            lto[lv] = obs_idx
    return lto


def _apply_custom_segmentation(ctx: ViewerContext, new_adata, scales: list):
    """Replace the cell labels layer and all derived state with custom data.

    scales: list of arrays (highest→lowest resolution), either numpy arrays
    (from source zarr) or dask arrays (from sdata DataTree on restore).
    """
    from palms.utils.coloring import CellColorManager

    # ── Replace labels layer in napari viewer ─────────────────────────────────
    old_layer = ctx.cell_labels_layer
    if old_layer is not None and old_layer in ctx.viewer.layers:
        ctx.viewer.layers.remove(old_layer)

    new_layer = ctx.viewer.add_labels(
        scales, name="cell_labels_custom", opacity=0.8
    )
    new_layer.contour = 0

    # ── Build new label_to_obs and CellColorManager ───────────────────────────
    new_l2o = _build_label_to_obs(new_adata)
    new_cm = CellColorManager(new_adata, new_l2o)

    # ── Compute new centroids ─────────────────────────────────────────────────
    centroids_um = np.asarray(new_adata.obsm["spatial"], dtype=np.float64)
    centroids_px = centroids_um / ctx.pixel_size
    new_centroids_yx = centroids_px[:, ::-1]

    # ── Update ctx ────────────────────────────────────────────────────────────
    ctx.cell_labels_layer = new_layer
    ctx.adata = new_adata
    ctx.label_to_obs = new_l2o
    ctx.color_manager = new_cm
    ctx.centroids_yx = new_centroids_yx
    ctx.segmentation_source = "custom"
    ctx.state["segmentation_source"] = "custom"

    # ── Reload clusterings from new adata ─────────────────────────────────────
    from palms.utils.adata_persistence import load_custom_clusterings_from_adata
    new_clusterings = load_custom_clusterings_from_adata(new_adata)

    # Also expose any cluster-like obs columns (seurat_clusters, cell_type, etc.)
    # Re-index by cell_id so get_cluster_colors can align them via reindex(cell_ids).
    import pandas as pd
    cell_ids_for_idx = (
        new_adata.obs["cell_id"].values if "cell_id" in new_adata.obs.columns else None
    )
    for col in new_adata.obs.columns:
        if col in ("cell_id", "region"):
            continue
        s = new_adata.obs[col]
        if s.dtype.kind in ("O", "U") or isinstance(s.dtype, pd.CategoricalDtype):
            if col not in new_clusterings:
                if cell_ids_for_idx is not None:
                    new_clusterings[col] = pd.Series(
                        s.values, index=cell_ids_for_idx, name=col
                    )
                else:
                    new_clusterings[col] = s

    ctx.clusterings.clear()
    ctx.clusterings.update(new_clusterings)
    ctx.clustering_names = list(new_clusterings.keys())

    # ── Refresh clustering UI across all tabs ─────────────────────────────────
    ctx.refresh_clustering_choices()

    # ── Clear stale analysis caches ───────────────────────────────────────────
    for key in ("nhood_result", "nhood_fig", "co_result", "co_fig",
                "rank_genes_df", "rank_genes_adata_norm", "rank_genes_groupby",
                "ligrec_result", "annot_dist_distances"):
        ctx.state.pop(key, None)

    # ── Reset label_to_cluster (cell coloring) ───────────────────────────────
    ctx.state["label_to_cluster"] = None
    ctx.state["active_clustering_name"] = None

    # ── Update gene names ─────────────────────────────────────────────────────
    ctx.gene_names = list(new_adata.var_names)


def _extract_dt_scales(dt) -> list:
    """Extract ordered dask arrays from a spatialdata DataTree (multiscale labels)."""
    import re

    def _sort_key(name):
        nums = re.findall(r'\d+', name)
        return int(nums[0]) if nums else 0

    scales = []
    for name in sorted(dt.children.keys(), key=_sort_key):
        child = dt.children[name]
        ds = getattr(child, 'ds', None)
        if ds is None:
            continue
        if 'image' in ds:
            scales.append(ds['image'].data)
        elif ds.data_vars:
            first = next(iter(ds.data_vars))
            scales.append(ds[first].data)
    return scales


def _revert_xenium_segmentation(ctx: ViewerContext):
    """Restore native Xenium cell labels layer and adata."""
    from palms.utils.coloring import CellColorManager

    if ctx.sdata is None:
        raise RuntimeError("sdata not available — cannot revert")

    # ── Reload native adata and label_to_obs ─────────────────────────────────
    loader_mod = _import_loader()
    new_adata = ctx.sdata["table"]
    new_l2o = loader_mod.get_label_to_obs_mapping(ctx.sdata)

    # ── Reload native label scales from sdata DataTree ───────────────────────
    scales = _extract_dt_scales(ctx.sdata.labels["cell_labels"])

    # ── Replace labels layer ──────────────────────────────────────────────────
    old_layer = ctx.cell_labels_layer
    if old_layer is not None and old_layer in ctx.viewer.layers:
        ctx.viewer.layers.remove(old_layer)

    new_layer = ctx.viewer.add_labels(scales, name="cell_labels", opacity=0.8)
    new_layer.contour = 0

    # ── Build CellColorManager ────────────────────────────────────────────────
    new_cm = CellColorManager(new_adata, new_l2o)

    # ── Recompute centroids ────────────────────────────────────────────────────
    centroids_um = np.asarray(new_adata.obsm["spatial"], dtype=np.float64)
    centroids_px = centroids_um / ctx.pixel_size
    new_centroids_yx = centroids_px[:, ::-1]

    # ── Update ctx ────────────────────────────────────────────────────────────
    ctx.cell_labels_layer = new_layer
    ctx.adata = new_adata
    ctx.label_to_obs = new_l2o
    ctx.color_manager = new_cm
    ctx.centroids_yx = new_centroids_yx
    ctx.segmentation_source = "xenium"
    ctx.state["segmentation_source"] = "xenium"

    # ── Reload original clusterings ───────────────────────────────────────────
    from palms.utils.adata_persistence import load_custom_clusterings_from_adata
    clusterings = loader_mod.load_clusterings(ctx.data_path)
    custom_from_adata = load_custom_clusterings_from_adata(new_adata)
    if custom_from_adata:
        clusterings.update(custom_from_adata)

    ctx.clusterings.clear()
    ctx.clusterings.update(clusterings)
    ctx.clustering_names = list(clusterings.keys())

    ctx.refresh_clustering_choices()

    # ── Clear stale caches ────────────────────────────────────────────────────
    for key in ("nhood_result", "nhood_fig", "co_result", "co_fig",
                "rank_genes_df", "rank_genes_adata_norm", "rank_genes_groupby",
                "ligrec_result", "annot_dist_distances"):
        ctx.state.pop(key, None)

    ctx.state["label_to_cluster"] = None
    ctx.state["active_clustering_name"] = None
    ctx.gene_names = list(new_adata.var_names)


def _import_loader():
    """Import 01_load_sdata.py (numeric prefix — importlib required)."""
    import importlib.util, sys
    scripts_dir = Path(__file__).parent.parent
    if "load_sdata" not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            "load_sdata", scripts_dir / "01_load_sdata.py"
        )
        mod = importlib.util.module_from_spec(spec)
        sys.modules["load_sdata"] = mod
        spec.loader.exec_module(mod)
    return sys.modules["load_sdata"]
