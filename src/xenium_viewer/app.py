"""
Xenium Viewer — Linux equivalent of Xenium Explorer.

Usage:
    conda activate xenium_viewer
    python scripts/02_xenium_viewer.py [/path/to/xenium/output]

If no path is given, a file dialog opens to select the dataset directory.

Opens a napari window with:
  - 4-channel morphology_focus image (DAPI, markers, 18S, SMA/Vim)
  - Cell and nucleus labels (raster masks) for fast coloring
  - Transcript points layer (populated on demand, supports up to 10 genes)
  - A docked control panel (gene/cluster selection, colormap, QV threshold)
  - A linked matplotlib UMAP window (deferred until first coloring)

Performance notes:
  - Morphology TIFFs are rendered via a 5-level software pyramid (no internal pyramid)
  - Cell boundaries (318K shapes) are skipped; use contour=2 on labels layer
  - Transcript loading uses feather cache (run 00_preprocess_transcripts.py first)
  - Label coloring uses DirectLabelColormap for O(nonzero) construction
  - Second launch uses zarr cache for ~60-70% faster startup
"""

import gc
import os
import sys
import json
import time
import argparse
import warnings
from datetime import datetime
from pathlib import Path

# ─── Prevent ICE/X11 EPIPE crash on Linux ────────────────────────────────
# libICE's default IO error handler calls exit() when it encounters a
# broken pipe on the X11 session manager socket. Override it to a no-op
# so the application survives the (harmless) session manager disconnect.
os.environ['SESSION_MANAGER'] = ''  # disable SM connection entirely

import ctypes
import ctypes.util
_ice_lib_path = ctypes.util.find_library('ICE')
if _ice_lib_path:
    try:
        _libICE = ctypes.CDLL(_ice_lib_path)
        _ICE_IO_ERROR_HANDLER = ctypes.CFUNCTYPE(None, ctypes.c_void_p)
        def _ice_io_error_noop(conn):
            pass  # swallow the error instead of calling exit()
        _ice_handler = _ICE_IO_ERROR_HANDLER(_ice_io_error_noop)
        _libICE.IceSetIOErrorHandler(_ice_handler)
    except (OSError, AttributeError):
        pass  # libICE not available — not a problem
# ──────────────────────────────────────────────────────────────────────────

import numpy as np
import napari

# Suppress non-critical warnings from spatialdata stack
warnings.filterwarnings("ignore", category=UserWarning, module="spatialdata")
warnings.filterwarnings("ignore", category=FutureWarning)

from xenium_viewer.utils.coloring import CellColorManager
from xenium_viewer.utils.transcript_index import TranscriptLoader
from xenium_viewer.utils.umap_widget import UMAPViewer
from xenium_viewer.utils.viewer_context import ViewerContext
from xenium_viewer.tabs._helpers import create_shared_helpers, create_preferences_menu, create_file_menu, create_view_menu
from xenium_viewer.utils.prov_graph import ProvGraph
from xenium_viewer import loader as _loader_mod
from xenium_viewer import preprocess as _preprocess_mod

# ─── Channel metadata ───────────────────────────────────────────────────────
CHANNEL_NAMES = [
    "DAPI",
    "ATP1A1-CD45-E-Cadherin",
    "18S",
    "AlphaSMA-Vimentin",
]
# Default contrast limits per channel (can be tuned)
CHANNEL_CONTRAST = [
    (0, 5000),
    (0, 3000),
    (0, 4000),
    (0, 3000),
]
CHANNEL_COLORMAPS = ["blue", "green", "red", "magenta"]


def _parse_args():
    """Parse CLI arguments and return the data directory path."""
    parser = argparse.ArgumentParser(
        description="Xenium Linux Viewer — napari-based spatial transcriptomics viewer"
    )
    parser.add_argument(
        "data_dir", nargs="?", default=None,
        help="Path to Xenium output directory (opens file dialog if omitted)",
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Skip zarr cache, load from raw Xenium output",
    )
    args = parser.parse_args()

    if args.data_dir:
        data_path = Path(args.data_dir)
    else:
        # Show a file dialog to pick the directory
        from qtpy.QtWidgets import QApplication, QFileDialog
        app = QApplication.instance() or QApplication(sys.argv)
        data_path_str = QFileDialog.getExistingDirectory(
            None, "Select Xenium Output Directory"
        )
        if not data_path_str:
            print("No directory selected. Starting with empty viewer.")
            return None, args.no_cache
        data_path = Path(data_path_str)

    # Validate only when a path was given
    experiment_file = data_path / "experiment.xenium"
    if not experiment_file.exists():
        print(f"Error: {experiment_file} not found. Is this a Xenium output directory?")
        sys.exit(1)

    return data_path, args.no_cache


def _read_pixel_size(data_path: Path) -> float:
    """Read pixel_size from experiment.xenium."""
    experiment_file = data_path / "experiment.xenium"
    with open(experiment_file) as f:
        meta = json.load(f)
    return float(meta["pixel_size"])


def _find_xenium_datasets(folder: Path) -> list:
    """Return Xenium dataset paths at or immediately under *folder*."""
    if (folder / "experiment.xenium").exists():
        return [folder]
    return sorted(
        p for p in folder.iterdir()
        if p.is_dir() and (p / "experiment.xenium").exists()
    )


from xenium_viewer.tabs._helpers import make_progress_dialog as _make_progress_dialog


from qtpy.QtCore import QThread, Signal as QtSignal


class PreprocessWorker(QThread):
    progress = QtSignal(int, str)   # (percent 0-100, message)
    finished = QtSignal(bool, str)  # (success, message)

    def __init__(self, datasets: list):
        super().__init__()
        self.datasets = datasets

    def run(self):
        preprocess = _preprocess_mod.preprocess
        loader_mod = _loader_mod

        n = len(self.datasets)
        total_steps = n * 2
        for i, ds in enumerate(self.datasets):
            # Step A: zarr cache (skip if already fresh)
            cache_path = ds / "sdata_cached.zarr"
            experiment_path = ds / "experiment.xenium"
            zarr_is_fresh = False
            if cache_path.exists():
                if not experiment_path.exists():
                    zarr_is_fresh = True          # no experiment.xenium to compare against
                elif cache_path.stat().st_mtime >= experiment_path.stat().st_mtime:
                    zarr_is_fresh = True

            pct = int(i * 2 * 100 / total_steps)
            if zarr_is_fresh:
                self.progress.emit(pct, f"[{i+1}/{n}] Zarr cache already up to date, skipping: {ds.name}")
            else:
                self.progress.emit(pct, f"[{i+1}/{n}] Creating zarr cache: {ds.name}")
                try:
                    loader_mod.load_sdata(ds, use_cache=True)
                except Exception as exc:
                    self.finished.emit(False, f"Zarr creation failed for {ds.name}:\n{exc}")
                    return
            # Step B: transcript feathers
            pct = int((i * 2 + 1) * 100 / total_steps)
            self.progress.emit(pct, f"[{i+1}/{n}] Preprocessing transcripts: {ds.name}")
            try:
                preprocess(
                    parquet_path=ds / "transcripts.parquet",
                    cache_dir=ds / "transcript_cache",
                )
            except Exception as exc:
                self.finished.emit(False, f"Transcript preprocessing failed for {ds.name}:\n{exc}")
                return
        self.finished.emit(True, f"Done. {n} dataset(s) preprocessed.")


def _extract_dt_scales(dt):
    """
    Extract an ordered list of dask arrays from a spatialdata DataTree.

    spatialdata 0.7.x stores multiscale data as a DataTree with children
    named 'scale0', 'scale1', ..., each holding a Dataset with a single
    variable called 'image'.

    Returns a list of dask arrays sorted from highest to lowest resolution,
    suitable for passing to napari's multiscale image/labels API.
    """
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
        # The data variable is always named 'image' in spatialdata 0.7.x
        if 'image' in ds:
            scales.append(ds['image'].data)
        elif ds.data_vars:
            first = next(iter(ds.data_vars))
            scales.append(ds[first].data)

    return scales


def _add_layers_manually(viewer, sdata):
    """Add images and labels from sdata DataTree to the napari viewer."""
    # ── Morphology image (multiscale pyramid from sdata) ─────────────────────
    if "morphology_focus" in sdata.images:
        print("  Adding morphology_focus (multiscale)...")
        scales = _extract_dt_scales(sdata.images["morphology_focus"])
        if scales:
            # Detect number of channels from the highest-res level
            n_ch = scales[0].shape[0]  # CYX layout
            if n_ch == len(CHANNEL_COLORMAPS):
                viewer.add_image(
                    scales,
                    name="morphology_focus",
                    channel_axis=0,
                    colormap=CHANNEL_COLORMAPS,
                    contrast_limits=CHANNEL_CONTRAST,
                    visible=True,
                )
            else:
                # Different channel count — use defaults
                print(f"  morphology_focus has {n_ch} channel(s) (expected {len(CHANNEL_COLORMAPS)}), using defaults")
                viewer.add_image(
                    scales,
                    name="morphology_focus",
                    channel_axis=0 if n_ch > 1 else None,
                    visible=True,
                )
        else:
            print("  Warning: could not extract morphology_focus scales")

    # ── Labels (multiscale pyramid from sdata) ────────────────────────────
    for key in ["cell_labels", "nucleus_labels"]:
        if key in sdata.labels:
            print(f"  Adding {key} (multiscale)...")
            scales = _extract_dt_scales(sdata.labels[key])
            if scales:
                viewer.add_labels(scales, name=key)
            else:
                print(f"  Warning: could not extract {key} scales")


def _build_control_panel(ctx: ViewerContext):
    """Build the tabbed control panel using per-tab modules.

    Returns (tab_widget, state, he_state, restore_session).
    """
    from qtpy.QtWidgets import QTabWidget
    from qtpy.QtCore import QTimer

    state = ctx.state
    he_state = ctx.he_state

    # ── Import tab builders ──────────────────────────────────────────────
    from xenium_viewer.tabs.tab_clustering import build_tab as build_clustering_tab
    from xenium_viewer.tabs.tab_cell_coloring import build_tab as build_cell_coloring_tab
    from xenium_viewer.tabs.tab_transcripts import build_tab as build_transcripts_tab
    from xenium_viewer.tabs.tab_umap import build_tab as build_umap_tab
    from xenium_viewer.tabs.tab_roi import build_tab as build_roi_tab
    from xenium_viewer.tabs.tab_he_registration import build_tab as build_he_registration_tab
    from xenium_viewer.tabs.tab_gene_analysis import build_tab as build_gene_analysis_tab
    from xenium_viewer.tabs.tab_marker_genes import build_tab as build_marker_genes_tab
    from xenium_viewer.tabs.tab_ligrec import build_tab as build_ligrec_tab
    from xenium_viewer.tabs.tab_nhood import build_tab as build_nhood_tab
    from xenium_viewer.tabs.tab_co_occurrence import build_tab as build_co_occurrence_tab
    from xenium_viewer.tabs.tab_arms import build_tab as build_arms_tab
    from xenium_viewer.tabs.tab_gene_correlation import build_tab as build_gene_correlation_tab
    from xenium_viewer.tabs.tab_cnv import build_tab as build_cnv_tab
    from xenium_viewer.tabs.tab_novae import build_tab as build_novae_tab
    from xenium_viewer.tabs.tab_notebook import build_tab as build_notebook_tab
    from xenium_viewer.tabs.tab_annotations import build_tab as build_annotations_tab
    from xenium_viewer.tabs.tab_annot_nhood import build_tab as build_annot_nhood_tab
    from xenium_viewer.tabs.tab_annot_distance import build_tab as build_annot_distance_tab
    from xenium_viewer.tabs.tab_segmentation import build_tab as build_segmentation_tab
    from xenium_viewer.tabs.tab_external_images import build_tab as build_external_images_tab
    from xenium_viewer.tabs.tab_patch_overlays import build_tab as build_patch_overlays_tab
    from xenium_viewer.tabs.tab_crop_dataset import build_tab as build_crop_dataset_tab

    # ── Build Cell Coloring first (creates cross-tab widgets) ────────────
    coloring_widget, coloring_exports = build_cell_coloring_tab(ctx)
    # ctx.clustering_widget etc. are now set

    # ── Build analysis tabs (they register their clustering widgets on ctx) ─
    ga_widget, ga_exports = build_gene_analysis_tab(ctx)
    mg_widget, mg_exports = build_marker_genes_tab(ctx)
    lr_widget, lr_exports = build_ligrec_tab(ctx)
    nhood_widget, nhood_exports = build_nhood_tab(ctx)
    co_widget, co_exports = build_co_occurrence_tab(ctx)
    cnv_widget, cnv_exports = build_cnv_tab(ctx)
    novae_widget, novae_exports = build_novae_tab(ctx)

    # ── Create shared helpers (needs all widgets registered) ─────────────
    create_shared_helpers(ctx)

    # ── Seed the code preamble so data_path/imports are always cell #1 ────
    # Emitted up front (when recording is on) so every recorded step has a
    # correct, self-contained preamble to depend on, even if the first action
    # a user takes doesn't chain through normalize/clustering.
    if ctx.state.get("record_code"):
        ctx.record_preamble()

    # ── Populate initial cluster checkboxes ──────────────────────────────
    ctx.repopulate_cluster_checkboxes()

    # ── Create Preferences menu ──────────────────────────────────────────
    create_preferences_menu(ctx)

    # ── Build remaining tabs ─────────────────────────────────────────────
    clustering_widget, clustering_exports = build_clustering_tab(ctx)
    transcripts_widget, transcripts_exports = build_transcripts_tab(ctx)
    umap_widget, umap_exports = build_umap_tab(ctx)
    roi_widget, roi_exports = build_roi_tab(ctx)
    he_widget, he_exports = build_he_registration_tab(ctx)
    arms_widget, arms_exports = build_arms_tab(ctx)
    corr_widget, corr_exports = build_gene_correlation_tab(ctx)
    notebook_widget, notebook_exports = build_notebook_tab(ctx)
    annot_widget, annot_exports = build_annotations_tab(ctx)
    annot_nhood_widget, annot_nhood_exports = build_annot_nhood_tab(ctx)
    annot_dist_widget, annot_dist_exports = build_annot_distance_tab(ctx)
    seg_widget, seg_exports = build_segmentation_tab(ctx)
    ext_img_widget, ext_img_exports = build_external_images_tab(ctx)
    patch_widget, patch_exports = build_patch_overlays_tab(ctx)
    crop_widget, crop_exports = build_crop_dataset_tab(ctx)

    # ── Mouse hover: show cluster ID in status bar ───────────────────────
    if ctx.cell_labels_layer is not None:
        _cursor_timer = QTimer()
        _cursor_timer.setSingleShot(True)
        _cursor_timer.setInterval(80)

        def _do_cursor_lookup():
            if ctx.cell_labels_layer is None:
                return
            lut = state["label_to_cluster"]
            if lut is None:
                return
            label_val = ctx.cell_labels_layer.get_value(
                ctx.viewer.cursor.position,
                view_direction=None,
                dims_displayed=list(range(ctx.cell_labels_layer.ndim)),
                world=True,
            )
            if isinstance(label_val, tuple):
                label_val = label_val[1]
            if label_val is not None and 0 < int(label_val) < len(lut):
                cid = lut[int(label_val)]
                name = state["active_clustering_name"] or "cluster"
                if cid >= 0:
                    raw_map = state.get('_cluster_id_to_raw')
                    raw_cid = raw_map[cid] if raw_map and cid in raw_map else cid
                    labels = ctx.get_active_labels()
                    display_cid = labels.get(raw_cid, labels.get(str(raw_cid), raw_cid))
                    if str(display_cid) != str(raw_cid):
                        text = f"Cell {int(label_val)} \u2014 {name}: {raw_cid} ({display_cid})"
                    else:
                        text = f"Cell {int(label_val)} \u2014 {name}: {display_cid}"
                else:
                    text = f"Cell {int(label_val)} \u2014 {name}: unassigned"
                QTimer.singleShot(0, lambda t=text: setattr(ctx.viewer, 'status', t))

        _cursor_timer.timeout.connect(_do_cursor_lookup)

        def _on_cursor_move(event):
            if ctx.cell_labels_layer is None:
                return
            _cursor_timer.start()  # restart resets the 80ms window

        ctx.viewer.cursor.events.position.connect(_on_cursor_move)

    # ── Assemble tabbed control panel (grouped) ──────────────────────────
    cells_tabs = QTabWidget()
    cells_tabs.addTab(clustering_widget,  "Clustering")
    cells_tabs.addTab(coloring_widget,    "Coloring")
    cells_tabs.addTab(transcripts_widget, "Transcripts")
    cells_tabs.addTab(umap_widget,        "UMAP")

    genes_tabs = QTabWidget()
    genes_tabs.addTab(ga_widget,   "Rank Genes")
    genes_tabs.addTab(mg_widget,   "Markers")
    genes_tabs.addTab(corr_widget, "Correlation")
    genes_tabs.addTab(cnv_widget,  "CNV")

    spatial_tabs = QTabWidget()
    spatial_tabs.addTab(roi_widget,        "ROI DEG")
    spatial_tabs.addTab(lr_widget,         "Lig-Rec")
    spatial_tabs.addTab(nhood_widget,      "Nhood Enrich")
    spatial_tabs.addTab(co_widget,         "Co-occur")
    spatial_tabs.addTab(novae_widget,      "Domains")
    spatial_tabs.addTab(annot_nhood_widget,"Annot Nhood")
    spatial_tabs.addTab(annot_dist_widget, "Annot Dist")

    images_tabs = QTabWidget()
    images_tabs.addTab(he_widget,      "H&E")
    images_tabs.addTab(arms_widget,    "ARMS")
    images_tabs.addTab(ext_img_widget, "Ext Images")
    images_tabs.addTab(patch_widget,   "Patches")

    tools_tabs = QTabWidget()
    tools_tabs.addTab(annot_widget,    "Annotations")
    tools_tabs.addTab(seg_widget,      "Segmentation")
    tools_tabs.addTab(crop_widget,     "Crop Dataset")
    tools_tabs.addTab(notebook_widget, "Notebook")

    for _group in (cells_tabs, genes_tabs, spatial_tabs, images_tabs, tools_tabs):
        _group.setTabPosition(QTabWidget.South)

    tab_widget = QTabWidget()
    tab_widget.addTab(cells_tabs,   "Cells")
    tab_widget.addTab(genes_tabs,   "Genes")
    tab_widget.addTab(spatial_tabs, "Spatial")
    tab_widget.addTab(images_tabs,  "Images")
    tab_widget.addTab(tools_tabs,   "Tools")

    # ── Compose session restore from per-tab restorers ───────────────────
    all_exports = [
        clustering_exports, coloring_exports, transcripts_exports,
        umap_exports, roi_exports, he_exports, ga_exports, mg_exports,
        lr_exports, nhood_exports, co_exports, cnv_exports, novae_exports, arms_exports, corr_exports,
        notebook_exports, annot_exports, annot_nhood_exports, annot_dist_exports,
        seg_exports, ext_img_exports, patch_exports, crop_exports,
    ]

    def restore_session(session):
        for exports in all_exports:
            restorer = exports.get("restore_session")
            if restorer:
                restorer(session)

    return tab_widget, state, he_state, restore_session


def _load_dataset(data_path: Path, no_cache: bool) -> dict:
    """Load all data for a Xenium dataset. No Qt or napari calls.

    Returns a dict with: pixel_size, sdata, adata, umap_df, clusterings,
    label_to_obs, gene_names, clustering_names, color_manager, transcript_loader.
    """
    pixel_size = _read_pixel_size(data_path)
    print(f"Pixel size: {pixel_size} um/px")

    # Repair zarr store if consolidated metadata is out of sync with disk.
    # A brief bug wrote adata_norm to sdata.tables without spatialdata_attrs,
    # which corrupted the consolidated_metadata in zarr.json even after the
    # directory was deleted, making sdata["table"] inaccessible on reload.
    if not no_cache:
        zarr_cache = data_path / "sdata_cached.zarr"
        if zarr_cache.exists():
            import shutil, json, zarr as _zarr
            # Remove any stray adata_norm directory
            bad_adata_norm = zarr_cache / "tables" / "adata_norm"
            need_consolidate = bad_adata_norm.exists()
            if need_consolidate:
                shutil.rmtree(bad_adata_norm, ignore_errors=True)
                print("  Removed invalid tables/adata_norm from zarr store")

            # Also detect stale consolidated metadata: tables/table/ exists on
            # disk but is absent from the inline zarr.json consolidated_metadata.
            if not need_consolidate:
                table_dir = zarr_cache / "tables" / "table"
                zarr_json = zarr_cache / "zarr.json"
                if table_dir.exists() and zarr_json.exists():
                    try:
                        root = json.loads(zarr_json.read_text())
                        tables_cm = (root.get("consolidated_metadata", {})
                                         .get("metadata", {})
                                         .get("tables", {})
                                         .get("consolidated_metadata", {})
                                         .get("metadata", {}))
                        if "table" not in tables_cm:
                            need_consolidate = True
                    except Exception:
                        pass

            if need_consolidate:
                try:
                    _zarr.consolidate_metadata(str(zarr_cache))
                    print("  Rebuilt zarr consolidated metadata")
                except Exception as e:
                    print(f"  Warning: could not consolidate zarr metadata: {e}")

    loader_mod = _loader_mod

    t0 = time.perf_counter()
    print("Loading SpatialData...")
    sdata = loader_mod.load_sdata(data_path, use_cache=not no_cache)
    print(f"  SpatialData loaded in {time.perf_counter() - t0:.1f}s")

    t0 = time.perf_counter()
    print("Loading UMAP...")
    umap_df = loader_mod.load_umap(data_path)
    print(f"  UMAP loaded in {time.perf_counter() - t0:.1f}s")

    t0 = time.perf_counter()
    print("Loading cluster assignments...")
    clusterings = loader_mod.load_clusterings(data_path)
    print(f"  Clusterings loaded in {time.perf_counter() - t0:.1f}s")

    print("Building label->obs mapping...")
    label_to_obs = loader_mod.get_label_to_obs_mapping(sdata)

    adata = sdata["table"]

    # Store UMAP coordinates in adata.obsm for portability
    if umap_df is not None and not umap_df.empty:
        umap_cols = [c for c in umap_df.columns if c.startswith("UMAP")]
        if len(umap_cols) >= 2:
            # umap_df is indexed by cell barcode; adata.obs.index may differ
            import pandas as _pd
            if "cell_id" in adata.obs.columns:
                cell_id_to_idx = _pd.Series(adata.obs.index, index=adata.obs["cell_id"].values)
                umap_reindexed = umap_df[umap_cols[:2]].rename(index=cell_id_to_idx).reindex(adata.obs.index)
            else:
                umap_reindexed = umap_df[umap_cols[:2]].reindex(adata.obs.index)
            adata.obsm["X_umap"] = umap_reindexed.values.astype(np.float32)
    elif "X_umap" in adata.obsm:
        # No raw analysis/umap/.../projection.csv (e.g. a Crop Dataset export
        # has no 'analysis/' folder at all), but the table already carries
        # UMAP coordinates embedded from a previous load of the source
        # dataset — rebuild umap_df from that instead of leaving the UMAP
        # tab empty.
        import pandas as _pd
        index = adata.obs["cell_id"].values if "cell_id" in adata.obs.columns else adata.obs_names
        umap_df = _pd.DataFrame(
            adata.obsm["X_umap"][:, :2], columns=["UMAP_1", "UMAP_2"], index=index,
        )
        print(f"  Reconstructed UMAP from adata.obsm['X_umap']: {len(umap_df)} cells")

    if umap_df is None:
        import pandas as _pd
        umap_df = _pd.DataFrame(columns=["UMAP_1", "UMAP_2"])

    # Load custom clusterings previously saved into adata.obs
    from xenium_viewer.utils.adata_persistence import load_custom_clusterings_from_adata
    custom_from_adata = load_custom_clusterings_from_adata(adata)
    if custom_from_adata:
        clusterings.update(custom_from_adata)
        print(f"  Loaded {len(custom_from_adata)} custom clustering(s) from adata.obs")

    gene_names = list(adata.var_names)
    clustering_names = list(clusterings.keys())
    print(f"Genes: {len(gene_names)}, Clusterings: {len(clustering_names)}")

    color_manager = CellColorManager(adata, label_to_obs)
    transcript_loader = TranscriptLoader(
        cache_dir=data_path / "transcript_cache",
        parquet_path=data_path / "transcripts.parquet",
        pixel_size=pixel_size,
    )

    return {
        "pixel_size": pixel_size,
        "sdata": sdata,
        "adata": adata,
        "umap_df": umap_df,
        "clusterings": clusterings,
        "label_to_obs": label_to_obs,
        "gene_names": gene_names,
        "clustering_names": clustering_names,
        "color_manager": color_manager,
        "transcript_loader": transcript_loader,
    }


def _populate_viewer(viewer, data: dict) -> dict:
    """Add layers to the napari viewer from loaded data.

    Returns a dict with: cell_labels_layer, transcript_layer, roi_layer,
    morph_thumb, morph_full_shape_yx, centroids_yx.
    """
    sdata = data["sdata"]
    adata = data["adata"]
    pixel_size = data["pixel_size"]

    t0 = time.perf_counter()
    _add_layers_manually(viewer, sdata)
    print(f"  Layers added in {time.perf_counter() - t0:.1f}s")

    # ── Extract lowest-res morphology thumbnail for coarse tissue alignment ───
    morph_thumb = None
    morph_full_shape_yx = None
    if "morphology_focus" in sdata.images:
        morph_scales = _extract_dt_scales(sdata.images["morphology_focus"])
        if morph_scales:
            morph_thumb = morph_scales[-1].compute()  # (C, Y, X), uint16
            morph_full = morph_scales[0]
            morph_full_shape_yx = (morph_full.shape[-2], morph_full.shape[-1])
            print(f"  Morphology thumbnail for coarse align: {morph_thumb.shape}")

    # ── Get label layers ─────────────────────────────────────────────────────
    cell_labels_layer = None
    nucleus_labels_layer = None
    for layer in viewer.layers:
        name = layer.name.lower()
        if "cell_label" in name and "nucleus" not in name:
            cell_labels_layer = layer
        elif "nucleus_label" in name:
            nucleus_labels_layer = layer

    if nucleus_labels_layer is not None:
        nucleus_labels_layer.contour = 2
        nucleus_labels_layer.opacity = 0.5

    if cell_labels_layer is not None:
        cell_labels_layer.contour = 0
        cell_labels_layer.opacity = 0.8

    # ── Hide all layers except cell_labels at startup ────────────────────────
    for layer in viewer.layers:
        layer.visible = (layer is cell_labels_layer)

    # ── Precompute cell centroids in pixel coordinates ──────────────────────
    centroids_um = adata.obsm['spatial']  # shape (N, 2), columns = (x, y)
    centroids_px = centroids_um / pixel_size
    centroids_yx = centroids_px[:, ::-1]  # shape (N, 2), columns = (y, x) for napari

    # ── ROI shapes layer ─────────────────────────────────────────────────────
    roi_layer = viewer.add_shapes(
        data=[], name="ROIs", shape_type="polygon",
        edge_color="white", face_color=[1, 1, 1, 0.1], edge_width=2,
    )

    # ── Annotation shapes layer ───────────────────────────────────────────────
    annotation_layer = viewer.add_shapes(
        data=[], name="Annotations",
        properties={"annotation_type": []},
        text={
            "string": "{annotation_type}",
            "size": 10,
            "color": "white",
            "anchor": "upper_left",
        },
        shape_type="polygon",
        edge_color="yellow", face_color=[1, 1, 0, 0.08], edge_width=2,
    )

    # ── Crop Dataset shapes layer ─────────────────────────────────────────────
    crop_layer = viewer.add_shapes(
        data=[], name="Crop Regions", shape_type="polygon",
        edge_color="orange", face_color=[1, 0.6, 0, 0.08], edge_width=2,
    )

    # ── Transcript points layer ───────────────────────────────────────────────
    transcript_layer = viewer.add_points(
        np.empty((0, 2), dtype=np.float32),
        name="transcripts",
        size=4,
        face_color="yellow",
        border_color="transparent",
        opacity=0.7,
        visible=False,
    )

    # ── Transcript density heatmap layer ─────────────────────────────────────
    transcript_bins_layer = viewer.add_image(
        np.zeros((1, 1), dtype=np.float32),
        name="transcript_density",
        colormap="hot",
        opacity=0.7,
        visible=False,
    )

    return {
        "cell_labels_layer": cell_labels_layer,
        "transcript_layer": transcript_layer,
        "transcript_bins_layer": transcript_bins_layer,
        "roi_layer": roi_layer,
        "annotation_layer": annotation_layer,
        "crop_layer": crop_layer,
        "morph_thumb": morph_thumb,
        "morph_full_shape_yx": morph_full_shape_yx,
        "centroids_yx": centroids_yx,
    }


def _snapshot_layers(ctx: ViewerContext) -> dict:
    """Capture napari layer data while Qt objects are still alive.

    Returns a snapshot dict suitable for passing to save_session().
    """
    snapshot = {}

    # ROI polygons
    snapshot["roi_data"] = [
        np.asarray(p, dtype=np.float64) for p in ctx.roi_layer.data
    ] if ctx.roi_layer is not None else []

    # H&E landmark layer data
    he_state = ctx.he_state
    xen_lm = he_state.get("xenium_lm_layer")
    he_lm = he_state.get("he_lm_layer")
    snapshot["xenium_landmarks"] = (
        np.asarray(xen_lm.data, dtype=np.float64)
        if xen_lm is not None and len(xen_lm.data) > 0 else None
    )
    snapshot["he_landmarks"] = (
        np.asarray(he_lm.data, dtype=np.float64)
        if he_lm is not None and len(he_lm.data) > 0 else None
    )

    # ARMS overlay landmark data
    arms_state = ctx.arms_state
    arms_xen_lm = arms_state.get("xenium_lm_layer")
    arms_he_lm = arms_state.get("he_lm_layer")
    snapshot["arms_xenium_landmarks"] = (
        np.asarray(arms_xen_lm.data, dtype=np.float64)
        if arms_xen_lm is not None and len(arms_xen_lm.data) > 0 else None
    )
    snapshot["arms_he_landmarks"] = (
        np.asarray(arms_he_lm.data, dtype=np.float64)
        if arms_he_lm is not None and len(arms_he_lm.data) > 0 else None
    )
    # ARMS state (exclude non-serializable napari layers)
    snapshot["arms_state"] = {
        k: v for k, v in arms_state.items()
        if k not in ("he_layer", "he_tif", "xenium_lm_layer", "he_lm_layer", "shapes_layer")
    }

    # External images — UI-only residuals (pixels + affine live in sdata.images)
    ext_ui = []
    for entry in (ctx.external_images_state or []):
        affine_matrix = None
        lyr = entry.get("layer_ref")
        if lyr is not None:
            try:
                m = np.asarray(lyr.affine.affine_matrix, dtype=np.float64)
                if not np.allclose(m, np.eye(m.shape[0]), atol=1e-6):
                    affine_matrix = m.tolist()
            except Exception:
                pass
        # Serialize channel states (strip non-JSON-safe types)
        ch_states = None
        raw_ch = entry.get("channel_states")
        if raw_ch:
            ch_states = [
                {
                    "visible": bool(cs.get("visible", True)),
                    "color": [float(c) for c in cs["color"][:3]],
                    "clim": [float(cs["clim"][0]), float(cs["clim"][1])],
                    "data_min": float(cs.get("data_min", 0)),
                    "data_max": float(cs.get("data_max", 65535)),
                }
                for cs in raw_ch
            ]
        ext_ui.append({
            "element_name": entry.get("element_name"),
            "path": entry.get("path"),
            "affine_source_name": entry.get("affine_source_name"),
            "affine_matrix": affine_matrix,
            "opacity": float(entry.get("opacity", 1.0)),
            "channel_states": ch_states,
            "flip_v": bool(entry.get("flip_v", False)),
            "flip_h": bool(entry.get("flip_h", False)),
        })
    snapshot["external_images_ui"] = ext_ui

    # Patch overlays — UI-only residuals (geometry + clusters live in sdata.shapes)
    patch_ui = []
    for entry in (ctx.patch_overlays_state or []):
        # Capture current affine matrix so it can be restored immediately,
        # even before the source layer (e.g. H&E) finishes loading.
        affine_matrix = None
        lyr = entry.get("shapes_layer")
        if lyr is not None:
            try:
                m = np.asarray(lyr.affine.affine_matrix, dtype=np.float64)
                if not np.allclose(m, np.eye(m.shape[0]), atol=1e-6):
                    affine_matrix = m.tolist()
            except Exception:
                pass
        patch_ui.append({
            "element_name": entry.get("element_name"),
            "source_path": entry.get("source_path"),
            "source_kind": entry.get("source_kind"),
            "active_cluster_column": entry.get("active_cluster_column"),
            "palette_name": entry.get("palette_name"),
            "patch_size_px": int(entry.get("patch_size_px", 0)),
            "affine_source_name": entry.get("affine_source_name"),
            "affine_matrix": affine_matrix,
            "outline_only": bool(entry.get("outline_only", False)),
            "edge_width": int(entry.get("edge_width", 2)),
            "opacity": float(entry.get("opacity", 0.8)),
            "confidence_threshold": float(entry.get("confidence_threshold", 0.0)),
            "hidden_cluster_ids": sorted(
                int(cid) for cid in (entry.get("hidden_cluster_ids") or set())
            ),
        })
    snapshot["patch_overlays_ui"] = patch_ui

    return snapshot


def _make_initial_state(gene_names: list, clustering_names: list) -> dict:
    """Return a fresh viewer state dict."""
    return {
        "current_gene": gene_names[0] if gene_names else None,
        "current_clustering": clustering_names[0] if clustering_names else None,
        "current_colormap": "viridis",
        "show_transcripts": False,
        "min_qv": 20,
        "color_mode": "Gene Expression",
        "transcript_genes": [],
        "filter_by_cluster": False,
        "label_to_cluster": None,
        "active_clustering_name": None,
        "nhood_result": None,
        "nhood_fig": None,
        "co_result": None,
        "co_fig": None,
        "plot_format": "svg",
        "plot_font_size": 10,
        "record_code": True,
        "code_journal": [],
        "code_journal_tags": set(),
        "prov_graph": ProvGraph(),
        "_legacy_counter": 0,
        # Stable filename so the recorded code accumulates into one file across
        # sessions (was a per-launch code_<timestamp>.py). The .ipynb sidecar is
        # written by the notebook export.
        "code_file": "analysis.py",
        "custom_clusterings": {},
    }


def _make_initial_he_state() -> dict:
    """Return a fresh H&E registration state dict."""
    return {
        "he_layer": None,
        "he_tif": None,
        "he_filename": None,
        "he_path": None,
        "he_shape_yx": None,
        "xenium_lm_layer": None,
        "he_lm_layer": None,
        "affine_3x3": None,
        "coarse_affine": None,
        "flip_v": False,
        "flip_h": False,
    }


def _make_initial_arms_state() -> dict:
    """Return a fresh ARMS overlay state dict."""
    return {
        "he_layer": None,
        "he_tif": None,
        "he_filename": None,
        "he_path": None,
        "he_shape_yx": None,
        "xenium_lm_layer": None,
        "he_lm_layer": None,
        "affine_3x3": None,
        "flip_v": False,
        "flip_h": False,
        "shapes_layer": None,
        "tile_names": None,
        "cluster_ids": None,
        "geojson_path": None,
        "csv_path": None,
        "cluster_checkboxes": {},
    }


def _do_full_init(viewer, data_path: Path, no_cache: bool, _app: dict) -> ViewerContext:
    """Load a dataset and build the full UI. Returns a fresh ViewerContext."""
    from xenium_viewer.utils.session import load_session

    data = _load_dataset(data_path, no_cache)
    adata = data["adata"]

    # Migrate old viewer_session data into adata (one-time, idempotent)
    zarr_path = data_path / "sdata_cached.zarr"
    if not no_cache and zarr_path.exists():
        from xenium_viewer.utils.adata_persistence import migrate_old_session_to_adata
        migrate_old_session_to_adata(zarr_path, data["sdata"], adata)

    cell_ids = adata.obs['cell_id'].values
    umap_viewer = UMAPViewer(data["umap_df"], cell_ids)
    layers = _populate_viewer(viewer, data)

    ctx = ViewerContext(
        viewer=viewer,
        adata=adata,
        sdata=data["sdata"],
        clusterings=data["clusterings"],
        color_manager=data["color_manager"],
        transcript_loader=data["transcript_loader"],
        umap_viewer=umap_viewer,
        label_to_obs=data["label_to_obs"],
        centroids_yx=layers["centroids_yx"],
        pixel_size=data["pixel_size"],
        data_path=data_path,
        no_cache=no_cache,
        gene_names=data["gene_names"],
        clustering_names=data["clustering_names"],
        cell_labels_layer=layers["cell_labels_layer"],
        transcript_layer=layers["transcript_layer"],
        transcript_bins_layer=layers["transcript_bins_layer"],
        roi_layer=layers["roi_layer"],
        annotation_layer=layers["annotation_layer"],
        crop_layer=layers["crop_layer"],
        morph_thumb=layers["morph_thumb"],
        morph_full_shape_yx=layers["morph_full_shape_yx"],
    )
    ctx.state = _make_initial_state(data["gene_names"], data["clustering_names"])
    ctx.he_state = _make_initial_he_state()
    ctx.arms_state = _make_initial_arms_state()

    # Remove old dock widget if present
    if _app["dock_widget"] is not None:
        try:
            _old_dw = _app["dock_widget"]
            if hasattr(_old_dw, "visibilityChanged"):
                try:
                    _old_dw.visibilityChanged.disconnect()
                except Exception:
                    pass
            viewer.window.remove_dock_widget(_old_dw)
        except Exception:
            pass
        _app["dock_widget"] = None

    panel, _state, _he_state, restore_fn = _build_control_panel(ctx)
    _app["dock_widget"] = viewer.window.add_dock_widget(
        panel, name="Xenium Controls", area="right"
    )
    _app["restore_fn"] = restore_fn

    # Sync the View menu checkbox when the dock is closed via its own close button
    _dw = _app["dock_widget"]
    if hasattr(_dw, "visibilityChanged"):
        def _on_dock_visibility(visible):
            action = _app.get("controls_action")
            if action is not None and action.isChecked() != visible:
                action.setChecked(visible)
        _dw.visibilityChanged.connect(_on_dock_visibility)

    # Load analysis results from adata.uns (nhood, co-occ, ligrec, rank genes)
    # Must run BEFORE restore_fn so tab _restore_session handlers find
    # the results in ctx.state and can enable their Show/Export buttons.
    from xenium_viewer.utils.adata_persistence import (
        load_analysis_results_from_adata,
        load_rank_genes_from_adata,
        load_cnv_results_from_adata,
        load_rois_from_sdata,
        load_landmarks_from_sdata,
        load_arms_tiles_from_sdata,
        migrate_landmarks_to_sdata,
    )
    analysis = load_analysis_results_from_adata(adata)
    for key in ("nhood_result", "co_result", "ligrec_result"):
        if analysis.get(key) is not None and ctx.state.get(key) is None:
            ctx.state[key] = analysis[key]

    rg_df, rg_adata_norm, rg_groupby = load_rank_genes_from_adata(adata, data["sdata"])
    if rg_df is not None and ctx.state.get('rank_genes_df') is None:
        ctx.state['rank_genes_df'] = rg_df
        ctx.state['rank_genes_adata_norm'] = rg_adata_norm
        ctx.state['rank_genes_groupby'] = rg_groupby

    cnv_results = load_cnv_results_from_adata(adata, data["sdata"])
    if cnv_results is not None and not ctx.state.get('cnv_results'):
        ctx.state['cnv_results'] = cnv_results

    # Load ROIs from sdata.shapes['rois'] (new format); zarr arrays in
    # load_session() serve as fallback for datasets not yet saved with new code.
    sdata_rois = load_rois_from_sdata(data["sdata"]) if not no_cache else []

    # Load DEG results from sdata sidecar parquets
    from xenium_viewer.utils.adata_persistence import (
        load_roi_deg_from_sdata, load_arms_tile_deg_from_sdata,
        load_cluster_labels_from_sdata,
    )
    sdata_roi_deg = load_roi_deg_from_sdata(data["sdata"]) if not no_cache else None
    sdata_arms_tile_deg = load_arms_tile_deg_from_sdata(data["sdata"]) if not no_cache else None
    sdata_cluster_labels = load_cluster_labels_from_sdata(data["sdata"]) if not no_cache else {}

    # Restore session if available
    if not no_cache and zarr_path.exists():
        session = load_session(zarr_path)
        if session is not None:
            # One-time migration: copy zarr landmark arrays and GeoJSON/CSV tile
            # data into sdata.shapes so they are portable without external files.
            migrate_landmarks_to_sdata(zarr_path, data["sdata"], session)
            # Load landmarks/tiles from sdata (captures newly migrated data too)
            _sdata = data["sdata"]
            he_xen_lm   = load_landmarks_from_sdata(_sdata, 'he_xenium_landmarks')
            he_he_lm    = load_landmarks_from_sdata(_sdata, 'he_he_landmarks')
            arms_xen_lm = load_landmarks_from_sdata(_sdata, 'arms_xenium_landmarks')
            arms_he_lm  = load_landmarks_from_sdata(_sdata, 'arms_he_landmarks')
            arms_tiles  = load_arms_tiles_from_sdata(_sdata)
            # sdata values take priority over zarr array fallbacks
            if sdata_rois:
                session['rois'] = sdata_rois
            if he_xen_lm is not None:
                session['xenium_landmarks'] = he_xen_lm
            if he_he_lm is not None:
                session['he_landmarks'] = he_he_lm
            if arms_xen_lm is not None:
                session['arms_xenium_landmarks'] = arms_xen_lm
            if arms_he_lm is not None:
                session['arms_he_landmarks'] = arms_he_lm
            if arms_tiles[0]:
                session['arms_tiles_sdata'] = arms_tiles
            # Inject DEG results from sdata (override None values from load_session)
            if sdata_roi_deg is not None and session.get('roi_deg_df') is None:
                session['roi_deg_df'] = sdata_roi_deg
            if sdata_arms_tile_deg is not None and session.get('arms_tile_deg_df') is None:
                session['arms_tile_deg_df'] = sdata_arms_tile_deg
            # Migrate cluster labels from session attrs → sdata obs columns (one-time)
            if not sdata_cluster_labels and session.get('cluster_labels'):
                from xenium_viewer.utils.adata_persistence import save_cluster_labels_to_sdata
                for _ck, _ld in session['cluster_labels'].items():
                    if isinstance(_ld, dict):
                        save_cluster_labels_to_sdata(ctx, _ck, _ld)
                sdata_cluster_labels = load_cluster_labels_from_sdata(data["sdata"])
                if sdata_cluster_labels:
                    print(f"  Migrated cluster labels for: {list(sdata_cluster_labels)}")
            # Inject cluster labels from sdata (sdata obs columns are authoritative)
            if sdata_cluster_labels:
                merged = dict(session.get('cluster_labels') or {})
                merged.update(sdata_cluster_labels)  # sdata wins on conflict
                session['cluster_labels'] = merged
            # Restore the reproducible-code provenance graph so the analysis
            # notebook accumulates across sessions (re-derive the flat journal
            # and rewrite analysis.py from the restored graph).
            if session.get("prov_graph"):
                from xenium_viewer.utils.prov_graph import ProvGraph, graph_to_cells
                _g = ProvGraph.from_list(session["prov_graph"])
                ctx.state["prov_graph"] = _g
                # Re-emit the preamble for THIS launch's data_path (upsert; a
                # no-op when the path is unchanged).
                if ctx.state.get("record_code"):
                    ctx.record_preamble()
                ctx.state["code_journal"] = [
                    c.source for c in graph_to_cells(_g) if c.cell_type == "code"
                ]
                try:
                    _cp = ctx.data_path / ctx.state.get("code_file", "analysis.py")
                    with open(_cp, "w") as _f:
                        _f.write("\n".join(ctx.state["code_journal"]) + "\n")
                except OSError:
                    pass
                _sync = ctx.state.get("_notebook_sync_fn")
                if _sync:
                    _sync()
                print(f"Restored code provenance graph: {len(_g)} node(s)")
            print("Restoring session from zarr cache...")
            restore_fn(session)
            print("Session restored.")
    elif sdata_rois or sdata_roi_deg is not None or sdata_arms_tile_deg is not None or sdata_cluster_labels:
        partial_session = {'rois': sdata_rois}
        if sdata_roi_deg is not None:
            partial_session['roi_deg_df'] = sdata_roi_deg
        if sdata_arms_tile_deg is not None:
            partial_session['arms_tile_deg_df'] = sdata_arms_tile_deg
        if sdata_cluster_labels:
            partial_session['cluster_labels'] = sdata_cluster_labels
        restore_fn(partial_session)

    viewer.title = f"Xenium Viewer — {data_path.name}"

    # ── Minimap overlay ───────────────────────────────────────────────────────
    if ctx.morph_thumb is not None and ctx.morph_full_shape_yx is not None:
        try:
            from xenium_viewer.utils.minimap_widget import MinimapWidget
            canvas_native = viewer.window._qt_viewer.canvas.native
            minimap = MinimapWidget(
                ctx.viewer, ctx.morph_thumb, ctx.morph_full_shape_yx, canvas_native
            )
            minimap.show()
            _app["minimap"] = minimap
            act = _app.get("minimap_action")
            if act is not None:
                act.setEnabled(True)
                act.setChecked(True)
        except Exception as exc:
            print(f"  Warning: minimap could not be created: {exc}")
            act = _app.get("minimap_action")
            if act is not None:
                act.setEnabled(False)
                act.setChecked(False)

    return ctx


def _push_to_console(viewer, ctx):
    """Inject key variables into napari's IPython console."""
    try:
        console = viewer.window._qt_viewer.console
        console.push({
            'viewer': viewer,
            'ctx': ctx,
            'adata': ctx.adata,
            'sdata': ctx.sdata,
            'clusterings': ctx.clusterings,
            'color_manager': ctx.color_manager,
            'gene_names': ctx.gene_names,
            'data_path': ctx.data_path,
        })
        console.execute(
            'print("\\nXenium Viewer variables available:\\n"'
            '"  adata, sdata, viewer, ctx, clusterings,\\n"'
            '"  color_manager, gene_names, data_path\\n"'
            '"\\nTip: recorded code is in ctx.data_path / ctx.state[\'code_file\']")',
            hidden=True,
        )
    except Exception as exc:
        print(f"  Warning: could not push variables to console: {exc}")


def run_viewer(data_path=None, no_cache: bool = False):
    print("=" * 60)
    print("Xenium Linux Viewer")
    if data_path:
        print(f"Dataset: {data_path}")
    else:
        print("No dataset — starting with empty viewer.")
    print("=" * 60)

    t_start = time.perf_counter()

    # ── Napari viewer (must be created before any QWidgets) ──────────────────
    viewer = napari.Viewer(title="Xenium Viewer")

    # ── Improve console accessibility ─────────────────────────────────────
    # Cap layer controls/list height so the console button stays visible,
    # and start with a larger window to allow console resizing.
    viewer.window.resize(1400, 900)
    qt_viewer = viewer.window._qt_viewer
    qt_viewer.dockLayerControls.setMaximumHeight(200)
    qt_viewer.dockLayerList.setMaximumHeight(200)

    # Force console to 300px when opened (napari doesn't allocate space by default)
    from qtpy.QtCore import Qt as _Qt, QTimer as _QTimer

    def _on_console_visibility(visible):
        if visible:
            _QTimer.singleShot(0, lambda: viewer.window._qt_window.resizeDocks(
                [qt_viewer.dockConsole], [300], _Qt.Vertical
            ))

    qt_viewer.dockConsole.visibilityChanged.connect(_on_console_visibility)

    # ── App-level mutable container ──────────────────────────────────────────
    _app = {
        "dock_widget": None,
        "restore_fn": None,
        "snapshot": {},
        "reload_in_progress": False,
    }

    ctx = None  # set after first dataset load

    # ── Detached CopyKAT background jobs: close/switch guard helpers ──────────
    def _running_bg_jobs():
        """Jobs whose done-file hasn't appeared yet (still running / detached)."""
        if ctx is None:
            return []
        jobs = ctx.state.get("_cnv_bg_jobs", []) if isinstance(ctx.state, dict) else []
        out = []
        for j in jobs:
            try:
                if not j["done_file"].exists():
                    out.append(j)
            except Exception:
                pass
        return out

    def _ask_close_bg(jobs):
        """Return 'stop' | 'continue' | 'cancel' for running CopyKAT jobs."""
        from qtpy.QtWidgets import QMessageBox
        box = QMessageBox(viewer.window._qt_window)
        box.setWindowTitle("CopyKAT analysis running")
        box.setText(
            f"{len(jobs)} CopyKAT background analysis(es) are still running.\n\n"
            "Stop them now, or let them continue in the background after the app closes "
            "(results will be picked up next time you open this dataset)?"
        )
        stop_btn = box.addButton("Stop analysis", QMessageBox.DestructiveRole)
        cont_btn = box.addButton("Continue in background", QMessageBox.AcceptRole)
        cancel_btn = box.addButton("Cancel", QMessageBox.RejectRole)
        box.setDefaultButton(cont_btn)
        box.exec_()
        clicked = box.clickedButton()
        if clicked is stop_btn:
            return "stop"
        if clicked is cancel_btn:
            return "cancel"
        return "continue"

    def _kill_bg_jobs(jobs):
        import os as _os
        import signal as _signal
        for j in jobs:
            try:
                _os.killpg(j["pgid"], _signal.SIGTERM)
            except Exception:
                try:
                    j["proc"].terminate()
                except Exception:
                    pass

    # ── Open Dataset callback ────────────────────────────────────────────────
    def _on_open_dataset():
        nonlocal ctx
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        if _app["reload_in_progress"]:
            return

        running = _running_bg_jobs()
        if running:
            choice = _ask_close_bg(running)
            if choice == "cancel":
                return
            if choice == "stop":
                _kill_bg_jobs(running)
        _app["reload_in_progress"] = True

        try:
            new_path_str = QFileDialog.getExistingDirectory(
                None, "Select Xenium Output Directory"
            )
            if not new_path_str:
                return

            new_path = Path(new_path_str)

            if not (new_path / "experiment.xenium").exists():
                QMessageBox.warning(
                    None,
                    "Invalid Directory",
                    f"No experiment.xenium found in:\n{new_path}\n\n"
                    "Please select a valid Xenium output directory.",
                )
                return

            print(f"\nOpening dataset: {new_path}")

            if ctx is not None:
                # 1. Snapshot current layer data while layers still alive
                _app["snapshot"] = _snapshot_layers(ctx)

                # 2. Save session for current dataset
                zarr_path_old = ctx.data_path / "sdata_cached.zarr"
                if not ctx.no_cache and zarr_path_old.exists():
                    from xenium_viewer.utils.adata_persistence import _persist_table
                    _persist_table(ctx)
                    ctx.state["segmentation_source"] = ctx.segmentation_source
                    from xenium_viewer.utils.session import save_session
                    save_session(zarr_path_old, ctx.state, ctx.he_state, _app["snapshot"])
                    try:
                        from xenium_viewer.utils.notebook_export import write_graph_notebook
                        _g = ctx.state.get("prov_graph")
                        if _g is not None and len(_g):
                            write_graph_notebook(_g, ctx.data_path / "analysis_notebook.ipynb")
                    except Exception:
                        pass
                    print("Session saved for previous dataset.")

                # 3. Close UMAP second window if open
                if ctx.umap_viewer is not None:
                    try:
                        ctx.umap_viewer.close()
                    except Exception:
                        pass
                    ctx.umap_viewer = None

                # 4. Increment generation counter (stale-worker guard)
                ctx.dataset_generation += 1

                # 5. Release old dataset objects before loading new one
                # Close any open external-image TiffFile handles and disconnect
                # mirrored affine callbacks so stale refs don't linger.
                for _entry in (ctx.external_images_state or []):
                    try:
                        cb = _entry.get("affine_disconnect")
                        if cb is not None:
                            cb()
                    except Exception:
                        pass
                    try:
                        tif = _entry.get("tif")
                        if tif is not None:
                            tif.close()
                    except Exception:
                        pass
                ctx.external_images_state = []
                for _entry in (ctx.patch_overlays_state or []):
                    try:
                        cb = _entry.get("affine_disconnect")
                        if cb is not None:
                            cb()
                    except Exception:
                        pass
                ctx.patch_overlays_state = []

                ctx.sdata = None
                ctx.adata = None
                ctx.clusterings = None
                ctx.color_manager = None
                ctx.transcript_loader = None
                ctx.label_to_obs = None
                ctx.gene_names = None
                ctx.clustering_names = None
                ctx.centroids_yx = None
                ctx.cell_labels_layer = None
                ctx.transcript_layer = None
                ctx.roi_layer = None
                ctx.annotation_layer = None
                ctx.crop_layer = None
                ctx.morph_thumb = None
                gc.collect()

            # Clean up old minimap before creating a new one
            old_minimap = _app.get("minimap")
            if old_minimap is not None:
                try:
                    old_minimap.hide()
                    old_minimap.deleteLater()
                except Exception:
                    pass
                _app["minimap"] = None
            act = _app.get("minimap_action")
            if act is not None:
                act.setEnabled(False)
                act.setChecked(False)

            # 6. Clear all layers
            viewer.layers.clear()

            # 7. Full init (loads data, builds ctx, control panel, restores session)
            try:
                ctx = _do_full_init(viewer, new_path, no_cache, _app)
                _push_to_console(viewer, ctx)
            except Exception as exc:
                QMessageBox.critical(
                    None,
                    "Dataset Load Error",
                    f"Failed to load dataset:\n{new_path}\n\nError:\n{exc}",
                )
                return

            print(f"Dataset opened: {new_path.name}")

        finally:
            _app["reload_in_progress"] = False

    # ── Preprocess Dataset callback ──────────────────────────────────────────
    def _on_preprocess_dataset():
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        folder_str = QFileDialog.getExistingDirectory(
            None, "Select Xenium Dataset or Parent Folder"
        )
        if not folder_str:
            return
        folder = Path(folder_str)

        datasets = _find_xenium_datasets(folder)
        if not datasets:
            QMessageBox.warning(
                None, "No Datasets Found",
                f"No Xenium datasets (experiment.xenium) found in:\n{folder}",
            )
            return

        session_line = (
            "\nThe current viewer session will be cleared to free memory."
            if ctx is not None else ""
        )
        if len(datasets) == 1:
            msg = (
                f"Preprocess dataset:\n  {datasets[0].name}\n\n"
                f"This will create the zarr cache and per-gene transcript feathers.{session_line}"
            )
        else:
            names = "\n".join(f"  \u2022 {d.name}" for d in datasets)
            msg = (
                f"Preprocess {len(datasets)} datasets:\n{names}\n\n"
                f"This will create zarr caches and per-gene transcript feathers.{session_line}"
            )
        reply = QMessageBox.question(
            None, "Preprocess Dataset", msg,
            QMessageBox.Ok | QMessageBox.Cancel,
        )
        if reply != QMessageBox.Ok:
            return

        # Clear viewer to free memory before heavy I/O
        viewer.layers.clear()
        if ctx is not None:
            ctx.sdata = None;  ctx.adata = None;  ctx.clusterings = None
            ctx.color_manager = None;  ctx.transcript_loader = None
            ctx.label_to_obs = None;  ctx.gene_names = None
            ctx.clustering_names = None;  ctx.centroids_yx = None
            ctx.cell_labels_layer = None;  ctx.transcript_layer = None
            ctx.roi_layer = None;  ctx.annotation_layer = None;  ctx.crop_layer = None
            ctx.morph_thumb = None
        gc.collect()

        dlg, bar, lbl = _make_progress_dialog("Preprocessing Datasets")
        worker = PreprocessWorker(datasets)

        def _on_progress(pct, msg):
            bar.setValue(pct)
            lbl.setText(msg)

        def _on_finished(success, msg):
            dlg.accept()
            if success:
                QMessageBox.information(
                    None, "Preprocessing Complete",
                    msg + "\n\nUse File \u2192 Open Dataset to load a dataset.",
                )
            else:
                QMessageBox.critical(None, "Preprocessing Failed", msg)

        worker.progress.connect(_on_progress)
        worker.finished.connect(_on_finished)
        worker.start()
        dlg.exec_()  # blocks Qt event loop until dlg.accept() is called

    # ── File / View menus — added once to napari's native menus ─────────────
    create_file_menu(viewer, _on_open_dataset, _on_preprocess_dataset)
    create_view_menu(viewer, _app)

    # ── Warn on close if a detached CopyKAT job is still running ─────────────
    # aboutToQuit fires too late to veto, so filter the window's Close event.
    from qtpy.QtCore import QObject as _QObject, QEvent as _QEvent

    class _CloseGuard(_QObject):
        def eventFilter(self, obj, event):
            if event.type() == _QEvent.Close:
                running = _running_bg_jobs()
                if running:
                    choice = _ask_close_bg(running)
                    if choice == "cancel":
                        event.ignore()
                        return True  # veto the close
                    if choice == "stop":
                        _kill_bg_jobs(running)
                    # 'continue' leaves the detached process running
            return False

    _close_guard = _CloseGuard(viewer.window._qt_window)
    viewer.window._qt_window.installEventFilter(_close_guard)
    _app["_close_guard"] = _close_guard  # keep a reference alive

    # ── Snapshot layer data before Qt teardown, then save on exit ────────────
    if not no_cache:
        def _on_viewer_closing(_event=None):
            if ctx is not None:
                _app["snapshot"] = _snapshot_layers(ctx)

        from qtpy.QtWidgets import QApplication
        QApplication.instance().aboutToQuit.connect(_on_viewer_closing)

    if data_path is not None:
        ctx = _do_full_init(viewer, data_path, no_cache, _app)
        _push_to_console(viewer, ctx)
        total_time = time.perf_counter() - t_start
        print(f"\nViewer ready in {total_time:.1f}s. Close the napari window to exit.")
    else:
        print("Viewer ready (no dataset loaded). Use File menu to open or preprocess.")

    napari.run()

    # ── Save session state on exit ────────────────────────────────────────
    if ctx is not None and not no_cache:
        final_zarr_path = ctx.data_path / "sdata_cached.zarr"
        if final_zarr_path.exists():
            from xenium_viewer.utils.adata_persistence import _persist_table, save_rois_to_sdata
            _persist_table(ctx)
            roi_data = _app["snapshot"].get("roi_data", [])
            save_rois_to_sdata(ctx, roi_data)
            ctx.state["segmentation_source"] = ctx.segmentation_source
            from xenium_viewer.utils.session import save_session
            save_session(final_zarr_path, ctx.state, ctx.he_state, _app["snapshot"])
            try:
                from xenium_viewer.utils.notebook_export import write_graph_notebook
                _g = ctx.state.get("prov_graph")
                if _g is not None and len(_g):
                    write_graph_notebook(_g, ctx.data_path / "analysis_notebook.ipynb")
            except Exception:
                pass
            print("Session saved to zarr cache.")


def main():
    """Console-script entry point — parses argv and launches the viewer."""
    data_path, no_cache = _parse_args()
    run_viewer(data_path, no_cache=no_cache)


if __name__ == "__main__":
    main()
