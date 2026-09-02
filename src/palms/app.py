"""
PALMS — Provenance-Aware Linking of Multimodal Spatial-omics.

Registers spatial transcriptomics, histology and genomic overlays into one
coordinate space, and records every action as replayable code. Reads Xenium 3.x
output; an open alternative to Xenium Explorer, which has no Linux build.

Usage:
    conda activate palms
    palms [/path/to/xenium/output]

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

import ctypes
import ctypes.util

# Safe to import at module scope: it pulls in nothing, and defers its own
# napari_mcp import until --mcp is actually used.
from palms import dev_mcp

_IS_LINUX = sys.platform.startswith('linux')   # WSL reports 'linux' too

# ─── Prevent ICE/X11 EPIPE crash on Linux ────────────────────────────────
# libICE's default IO error handler calls exit() when it encounters a
# broken pipe on the X11 session manager socket. Override it to a no-op
# so the application survives the (harmless) session manager disconnect.
#
# Guarded to Linux: macOS has no session manager, and clearing SESSION_MANAGER
# there would be an unexplained edit to the user's environment.
if _IS_LINUX:
    os.environ['SESSION_MANAGER'] = ''  # disable SM connection entirely

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


# ─── Warn before the GLX double-load aborts the process ───────────────────
# Must run BEFORE `import napari`: that import is what triggers the abort, so a
# check placed after it never executes. The logic lives in utils.gl_check so it
# can be tested without importing napari. It reports only — see that module for
# why preloading conda's libGLX.so.0 does not repair the collision.
from palms.utils.gl_check import warn_if_libglx_will_collide

warn_if_libglx_will_collide()
# ──────────────────────────────────────────────────────────────────────────

import numpy as np
import napari

# Suppress non-critical warnings from spatialdata stack
warnings.filterwarnings("ignore", category=UserWarning, module="spatialdata")
warnings.filterwarnings("ignore", category=FutureWarning)

from palms.utils.coloring import CellColorManager
from palms.utils.transcript_index import TranscriptLoader
from palms.utils.umap_widget import UMAPViewer
from palms.utils.viewer_context import ViewerContext
from palms.utils.units import layer_affine_px
from palms.tabs._helpers import (
    create_shared_helpers, create_preferences_menu, create_file_menu,
    create_view_menu, ensure_plots_dock as _ensure_plots_dock,
    reveal_plots_dock as _reveal_plots_dock,
)
from palms.utils.plot_output import DEFAULT_FORMATS
from palms.utils.prov_graph import ProvGraph
from palms import loader as _loader_mod
from palms import preprocess as _preprocess_mod

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
    parser.add_argument(
        "--no-user-templates", action="store_true",
        help="Ignore customised analysis templates and run the shipped ones",
    )
    parser.add_argument(
        "--mcp", nargs="?", type=int, const=dev_mcp.DEFAULT_PORT, default=None,
        metavar="PORT",
        help="Dev only: expose this viewer to an AI assistant over an "
             "unauthenticated localhost MCP bridge that can execute arbitrary "
             "Python in this process. Requires the 'mcp' extra.",
    )
    args = parser.parse_args()

    # Applied before anything resolves a template. This is the first thing to
    # try when a result is in doubt, so it must not require finding and moving
    # files — "run it again with --no-user-templates" has to be a one-liner.
    if args.no_user_templates:
        from palms.utils.step_templates import set_overrides_enabled
        set_overrides_enabled(False)

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
            return None, args.no_cache, args.mcp
        data_path = Path(data_path_str)

    # Validate only when a path was given
    experiment_file = data_path / "experiment.xenium"
    if not experiment_file.exists():
        print(f"Error: {experiment_file} not found. Is this a Xenium output directory?")
        sys.exit(1)

    return data_path, args.no_cache, args.mcp


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


from palms.tabs._helpers import make_progress_dialog as _make_progress_dialog


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


def _warn_if_pyramid_is_not_stored(sdata):
    """Say so before building layers from a pyramid that is still a computation.

    Normally every level of an image element is an array on disk, and napari
    reads a few chunks of the smallest one. When the dataset did not come from a
    cache (`--no-cache`, or a build whose write failed) only `scale0` holds data;
    the rest is a chained `coarsen().mean()`. napari draws the *smallest* level
    first, so adding the layer walks that whole chain, and the chain is not
    streamable — each level is rechunked either side of the coarsen, which is a
    many-to-many gather, so a full level of float64 has to be live at once. On a
    full slide that is ~24 GB for `scale1` alone, and capping dask's concurrency
    does not change it (measured: 23.1 GB uncapped, 25.2 GB at 4 workers, both
    dead at a 40 GB cap).

    There is no cheap fix here — the pyramid genuinely has to be computed — so
    this warns instead of pretending otherwise. Without a cache to read, a
    full-slide dataset needs that memory, and the failure mode is a killed
    session rather than an error, which is worth a sentence of warning.
    """
    from palms.utils.raster_io import level_is_computed

    # Every image element, not just morphology_focus: he_image, arms_he_image
    # and the ext_* images reach napari the same way and cost the same when
    # their levels are a computation. Scoping this to one element is why an H&E
    # could kill a session without the viewer having said anything first.
    for name in list(sdata.images or {}):
        element = sdata.images.get(name)
        if element is None or not hasattr(element, "children"):
            continue
        scales = _extract_dt_scales(element)
        if len(scales) < 2 or not level_is_computed(scales[-1]):
            continue
        msg = (f"{name} has no stored pyramid — its lower levels will be "
               "computed from the full-resolution image while the layers are built. "
               "On a full slide this needs tens of GB; if the session dies without "
               "a traceback, check `journalctl -u systemd-oomd`.")
        print(f"  Warning: {msg}")
        try:
            from palms.utils import reporting
            reporting.get_logger().warning(msg)
        except Exception:
            pass


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
    from palms.tabs.tab_clustering import build_tab as build_clustering_tab
    from palms.tabs.tab_cell_coloring import build_tab as build_cell_coloring_tab
    from palms.tabs.tab_transcripts import build_tab as build_transcripts_tab
    from palms.tabs.tab_umap import build_tab as build_umap_tab
    from palms.tabs.tab_roi import build_tab as build_roi_tab
    from palms.tabs.tab_he_registration import build_tab as build_he_registration_tab
    from palms.tabs.tab_gene_analysis import build_tab as build_gene_analysis_tab
    from palms.tabs.tab_marker_genes import build_tab as build_marker_genes_tab
    from palms.tabs.tab_ligrec import build_tab as build_ligrec_tab
    from palms.tabs.tab_nhood import build_tab as build_nhood_tab
    from palms.tabs.tab_co_occurrence import build_tab as build_co_occurrence_tab
    from palms.tabs.tab_arms import build_tab as build_arms_tab
    from palms.tabs.tab_gene_correlation import build_tab as build_gene_correlation_tab
    from palms.tabs.tab_cnv import build_tab as build_cnv_tab
    from palms.tabs.tab_novae import build_tab as build_novae_tab
    from palms.tabs.tab_notebook import build_tab as build_notebook_tab
    from palms.tabs.tab_annotations import build_tab as build_annotations_tab
    from palms.tabs.tab_annot_nhood import build_tab as build_annot_nhood_tab
    from palms.tabs.tab_annot_distance import build_tab as build_annot_distance_tab
    from palms.tabs.tab_segmentation import build_tab as build_segmentation_tab
    from palms.tabs.tab_external_images import build_tab as build_external_images_tab
    from palms.tabs.tab_patch_overlays import build_tab as build_patch_overlays_tab
    from palms.tabs.tab_crop_dataset import build_tab as build_crop_dataset_tab
    from palms.tabs.tab_dataset import build_tab as build_dataset_tab
    from palms.tabs.tab_cache import build_tab as build_cache_tab
    from palms.tabs.tab_templates import build_tab as build_templates_tab

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
    dataset_widget, dataset_exports = build_dataset_tab(ctx)
    cache_widget, cache_exports = build_cache_tab(ctx)
    templates_widget, templates_exports = build_templates_tab(ctx)

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
    tools_tabs.addTab(dataset_widget,  "Dataset")
    tools_tabs.addTab(cache_widget,    "Cache")
    tools_tabs.addTab(templates_widget, "Templates")

    for _group in (cells_tabs, genes_tabs, spatial_tabs, images_tabs, tools_tabs):
        _group.setTabPosition(QTabWidget.TabPosition.South)

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
        dataset_exports, cache_exports, templates_exports,
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
    # Start the per-dataset log before anything can fail. Write failures used to
    # go only to stdout, which a GUI user never reads — so when the cache was
    # being corrupted, the warnings that would have explained it were lost.
    from palms.utils.reporting import reset_failures, setup_logging
    reset_failures()
    log_file = setup_logging(data_path)
    if log_file is not None:
        print(f"Logging to {log_file}")

    pixel_size = _read_pixel_size(data_path)
    print(f"Pixel size: {pixel_size} um/px")

    # Finish any write interrupted by a crash, and fix bookkeeping that no
    # longer matches disk, before anything tries to open the store. This
    # generalises an earlier fix that only knew about tables/table and a stray
    # tables/adata_norm; utils/cache_repair diffs the whole consolidated
    # metadata against the filesystem, and repairs without discarding.
    if not no_cache:
        zarr_cache = data_path / "sdata_cached.zarr"
        if zarr_cache.exists():
            from palms.utils import cache_repair
            report = cache_repair.verify(zarr_cache)
            if not report.ok:
                result = cache_repair.repair(zarr_cache, report)
                for action in result.actions:
                    print(f"  Cache: {action}")
                for failure in result.failures:
                    print(f"  Cache warning: {failure}")

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
    from palms.utils.adata_persistence import load_custom_clusterings_from_adata
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


def _install_unit_scaling(viewer, pixel_size: float) -> None:
    """Put this viewer's world into micrometres and turn the scale bar on.

    napari 0.8 has no ``scale_bar.unit``; the unit lives on the layers
    (``layer.units``, pint-backed) and propagates to ``viewer.dims.units``, while
    the magnitude has to live in ``layer.scale``. A scaled unit string is not an
    option — ``"0.2125 um"`` is accepted and silently loses its magnitude.
    """
    from palms.utils import units as _units

    def _stamp(layer):
        if _units.apply_to_layer(layer, pixel_size):
            return
        # The one case napari's own warning would have been right about. Reported
        # here instead, because this names the layer and that one does not — and
        # because the window that suppresses it closes on the next line.
        from palms.utils import reporting
        reporting.report_layer_scaling_failure(str(getattr(layer, "name", "?")))

    # napari evaluates unit consistency mid-insertion, before the new layer has
    # been stamped, and warns about a state that is over by the next draw. See
    # utils/units.py for the trace and why connection order cannot fix it.
    def _on_inserting(_event=None):
        _units.quiet_insertion().__enter__()

    def _on_inserted(event):
        try:
            _stamp(event.value)
        finally:
            _units.quiet_insertion().__exit__(None, None, None)

    for layer in viewer.layers:
        _stamp(layer)
    viewer.layers.events.inserting.connect(_on_inserting)
    viewer.layers.events.inserted.connect(_on_inserted)

    viewer.scale_bar.visible = True
    viewer.scale_bar.colored = False


def _populate_viewer(viewer, data: dict) -> dict:
    """Add layers to the napari viewer from loaded data.

    Returns a dict with: cell_labels_layer, transcript_layer, roi_layer,
    morph_thumb, morph_full_shape_yx, centroids_yx.
    """
    sdata = data["sdata"]
    adata = data["adata"]
    pixel_size = data["pixel_size"]

    _warn_if_pyramid_is_not_stored(sdata)

    # ── Physical units on the canvas ─────────────────────────────────────────
    # Installed before the first layer exists, and left connected, so *every*
    # layer any tab adds later is in micrometres too. There are more than twenty
    # add_image/add_labels/add_points/add_shapes sites across eight modules; one
    # of them missed would put that layer in pixels while the rest are in µm,
    # which misplaces it rather than merely mislabelling it. See utils/units.py.
    _install_unit_scaling(viewer, pixel_size)

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

    # The preview draws into its *own* layer rather than borrowing the one
    # above. Two layers make the two states mutually exclusive by construction:
    # renaming a single layer to say "preview" races the worker that fills it,
    # and a rename that lands late puts the wrong caption over the right data.
    # The name is the label — a napari layer name is user-visible, and this one
    # has to say what it is even when the Transcripts tab is not open.
    transcript_bins_preview_layer = viewer.add_image(
        np.zeros((1, 1), dtype=np.float32),
        name="transcript_density (PREVIEW - not recorded)",
        colormap="hot",
        opacity=0.7,
        visible=False,
    )

    return {
        "cell_labels_layer": cell_labels_layer,
        "transcript_layer": transcript_layer,
        "transcript_bins_layer": transcript_bins_layer,
        "transcript_bins_preview_layer": transcript_bins_preview_layer,
        "roi_layer": roi_layer,
        "annotation_layer": annotation_layer,
        "crop_layer": crop_layer,
        "morph_thumb": morph_thumb,
        "morph_full_shape_yx": morph_full_shape_yx,
        "centroids_yx": centroids_yx,
    }


def _load_prov_graph_items(data_path, session: dict) -> list:
    """The serialized provenance graph, sidecar first, session attr second.

    Two writers, two cadences: ``_helpers._save_prov_graph`` rewrites the
    sidecar on every recorded step, while ``save_session`` writes the attr only
    on a dataset switch or at exit. The sidecar is therefore never *behind* the
    attr and is usually ahead of it, so it wins whenever it parses. A dataset
    that predates the sidecar, or one copied without ``viewer_cache/``, still
    restores from the attr.
    """
    from palms.tabs._helpers import PROV_GRAPH_SIDECAR
    from palms.utils.adata_persistence import sidecar_dir

    attr_items = session.get("prov_graph") or []
    if data_path is None:
        return list(attr_items)
    sidecar = sidecar_dir(data_path) / PROV_GRAPH_SIDECAR
    try:
        items = json.loads(sidecar.read_text())
    except FileNotFoundError:
        return list(attr_items)
    except (OSError, ValueError) as e:
        print(f"  Provenance sidecar unreadable ({e}); using the session attr.")
        return list(attr_items)
    if not isinstance(items, list) or not items:
        return list(attr_items)
    if len(items) < len(attr_items):
        # The sidecar is written on every recorded step and the attr only at
        # exit, so the sidecar is never legitimately *smaller*. If it is,
        # something wrote a partial graph (see prov_graph_restored) and the attr
        # is the better record — losing an analysis to a stale file is the one
        # outcome worth being conservative about.
        #
        # The Notebook tab's "Drop Stale Nodes" is the one deliberate shrink, and
        # it does not arrive here as a mismatch: tab_notebook._persist_pruned_graph
        # writes the pruned graph to the session attr *and* the sidecar in the
        # same action, so both copies are the same size and this branch never
        # sees it. A lone shrink still means what it always meant — a partial
        # write — and is still refused.
        print(f"  Provenance sidecar has {len(items)} node(s) but the session "
              f"attr has {len(attr_items)}; using the session attr.")
        return list(attr_items)
    if len(items) != len(attr_items):
        print(f"  Provenance graph: {len(items)} node(s) from {sidecar.name} "
              f"(session attr had {len(attr_items)})")
    return items


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
                # Stored in pixels, like every other affine on disk — the layer's
                # own copy is in µm. See utils/units.py.
                m = layer_affine_px(lyr, ctx.pixel_size)
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
                m = layer_affine_px(lyr, ctx.pixel_size)          # pixels, see above
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
        # Every plot is written in each of these; see utils/plot_output.py.
        "plot_formats": list(DEFAULT_FORMATS),
        "plot_font_size": 10,
        "umap_genes": [],
        # Global CPU-core budget for parallel analyses (currently CopyKAT's
        # parallelDist). Default to half the machine's cores, leaving headroom.
        "n_cores": max(1, (os.cpu_count() or 2) // 2),
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
    from palms.utils.session import load_session

    data = _load_dataset(data_path, no_cache)
    adata = data["adata"]

    # Migrate old viewer_session data into adata (one-time, idempotent)
    zarr_path = data_path / "sdata_cached.zarr"
    if not no_cache and zarr_path.exists():
        from palms.utils.adata_persistence import migrate_old_session_to_adata
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
        transcript_bins_preview_layer=layers["transcript_bins_preview_layer"],
        roi_layer=layers["roi_layer"],
        annotation_layer=layers["annotation_layer"],
        crop_layer=layers["crop_layer"],
        morph_thumb=layers["morph_thumb"],
        morph_full_shape_yx=layers["morph_full_shape_yx"],
    )
    ctx.state = _make_initial_state(data["gene_names"], data["clustering_names"])
    ctx.he_state = _make_initial_he_state()
    ctx.arms_state = _make_initial_arms_state()
    # Bound here rather than once at startup because every reload builds a new
    # ViewerContext; the callable itself lives in _app and outlives them all.
    ctx.reload_dataset = _app.get("reload_current_dataset")

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

    # Remove the old Plots dock too — its cards hold figures from the dataset
    # being replaced, and its gallery would otherwise mix two datasets' plots.
    if _app.get("plots_dock") is not None:
        try:
            viewer.window.remove_dock_widget(_app["plots_dock"])
        except Exception:
            pass
        _app["plots_dock"] = None
        _app["plots_panel"] = None

    panel, _state, _he_state, restore_fn = _build_control_panel(ctx)
    _app["dock_widget"] = viewer.window.add_dock_widget(
        panel, name="Controls", area="right"
    )
    _app["restore_fn"] = restore_fn

    # The Plots dock: every figure any tab produces lands here. Built after the
    # control panel so it docks below it, and hidden until the first plot —
    # ``ctx.show_plot`` reveals it.
    from palms.utils.plots_panel import PlotsPanel
    plots_panel = PlotsPanel()
    _app["plots_panel"] = plots_panel
    ctx.plots_panel = plots_panel
    # ensure_plots_dock owns the dock's lifecycle — it also re-creates the dock
    # if napari destroyed it, which its title-bar close button does.
    plots_dock = _ensure_plots_dock(viewer, _app)
    if plots_dock is not None:
        plots_dock.setVisible(False)
    ctx.reveal_plots_dock = lambda: _reveal_plots_dock(viewer, _app)

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
    from palms.utils.adata_persistence import (
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
    from palms.utils.adata_persistence import (
        load_roi_deg_from_sdata, load_arms_tile_deg_from_sdata,
        load_cluster_labels_from_sdata,
    )
    sdata_roi_deg = load_roi_deg_from_sdata(data["sdata"]) if not no_cache else None
    sdata_arms_tile_deg = load_arms_tile_deg_from_sdata(data["sdata"]) if not no_cache else None
    sdata_cluster_labels = load_cluster_labels_from_sdata(data["sdata"]) if not no_cache else {}

    # Restore. A store with no `viewer_session` is still a store worth restoring
    # from: a crop export, a cache built by `palms-build-cache`, a session node
    # deleted in Tools -> Dataset, a recovered cache. All of those hold the
    # elements — the registered H&E, the ROIs, the landmarks, the cluster label
    # columns — and used to open with none of it shown, because the whole restore
    # hung off `load_session()` returning something. It is the elements that are
    # authoritative; the session only adds what has no element to live in.
    stored = load_session(zarr_path) if (not no_cache and zarr_path.exists()) else None
    have_session = stored is not None
    session = stored or {}

    # Which segmentation the recorded preamble binds ``adata`` from. Seeded here,
    # from the session, because the preamble is re-emitted below *before*
    # ``restore_fn`` reaches Tools → Segmentation's own restore handler: emitting
    # the Xenium form first and letting that handler correct it would upsert a
    # changed preamble on every launch, flagging the entire notebook stale for
    # nothing — the same defect a manual dataset rename used to cause.
    ctx.state.setdefault(
        "segmentation_source", session.get("segmentation_source", "xenium")
    )

    if have_session:
        # One-time migration of legacy zarr landmark arrays and GeoJSON/CSV tile
        # data into sdata.shapes. Gated because it reads the group it migrates
        # from — it self-guards on the group's absence, so this is documentation.
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
        from palms.utils.adata_persistence import save_cluster_labels_to_sdata
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
    # The sidecar is rewritten on every recorded step; the session attr
    # only by save_session (dataset switch / exit), so it is behind
    # whenever the last run ended in a crash, a kill, or a still-open
    # viewer. Prefer the sidecar, fall back to the attr.
    _prov_items = _load_prov_graph_items(ctx.data_path, session)
    if _prov_items:
        from palms.utils.prov_graph import ProvGraph, graph_to_cells
        _g = ProvGraph.from_list(_prov_items)
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

    print("Restoring session from zarr cache..." if have_session
          else "No stored session; restoring from the store's own elements...")
    restore_fn(session)
    print("Session restored.")

    # Only now may the graph be written back. Everything above this line runs
    # with whatever the tabs seeded — at minimum the preamble emitted during
    # construction — and persisting *that* would overwrite the session's real
    # graph with a one-node stub, which the next launch would then prefer.
    ctx.state["prov_graph_restored"] = True
    ctx.save_prov_graph()

    viewer.title = f"PALMS — {data_path.name}"

    # ── Minimap overlay ───────────────────────────────────────────────────────
    if ctx.morph_thumb is not None and ctx.morph_full_shape_yx is not None:
        try:
            from palms.utils.minimap_widget import MinimapWidget
            canvas_native = viewer.window._qt_viewer.canvas.native
            minimap = MinimapWidget(
                ctx.viewer, ctx.morph_thumb, ctx.morph_full_shape_yx, canvas_native,
                pixel_size=ctx.pixel_size,
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
    # The one place a *fresh* ctx is published to an interactive surface — this
    # runs on first load and again on every dataset switch — so it is also
    # where the --mcp bridge's namespace gets its handle. See dev_mcp.
    dev_mcp.publish_context(viewer, ctx)
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
            'print("\\nPALMS variables available:\\n"'
            '"  adata, sdata, viewer, ctx, clusterings,\\n"'
            '"  color_manager, gene_names, data_path\\n"'
            '"\\nTip: recorded code is in ctx.data_path / ctx.state[\'code_file\']")',
            hidden=True,
        )
    except Exception as exc:
        print(f"  Warning: could not push variables to console: {exc}")


def run_viewer(data_path=None, no_cache: bool = False, mcp_port: int | None = None):
    print("=" * 60)
    print("Xenium Linux Viewer")
    if data_path:
        print(f"Dataset: {data_path}")
    else:
        print("No dataset — starting with empty viewer.")
    print("=" * 60)

    t_start = time.perf_counter()

    # ── Napari viewer (must be created before any QWidgets) ──────────────────
    viewer = napari.Viewer(title="PALMS")

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
        "plots_dock": None,
        "plots_panel": None,
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
        stop_btn = box.addButton("Stop analysis", QMessageBox.ButtonRole.DestructiveRole)
        cont_btn = box.addButton("Continue in background", QMessageBox.ButtonRole.AcceptRole)
        cancel_btn = box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        box.setDefaultButton(cont_btn)
        box.exec()
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

    # ── Open / reload dataset ────────────────────────────────────────────────
    def _load_dataset_into_viewer(new_path: Path) -> bool:
        """Tear down the current dataset and load *new_path* in its place.

        Passing the dataset already open reloads it, which is how the Cache tab
        makes recovered elements visible: they were written straight into the
        zarr, so the in-memory sdata, the layers and the tab widgets know
        nothing about them until everything is rebuilt from disk.

        Returns True if a dataset was loaded.
        """
        nonlocal ctx
        from qtpy.QtWidgets import QMessageBox

        if _app["reload_in_progress"]:
            return False

        running = _running_bg_jobs()
        if running:
            choice = _ask_close_bg(running)
            if choice == "cancel":
                return False
            if choice == "stop":
                _kill_bg_jobs(running)
        _app["reload_in_progress"] = True

        try:
            print(f"\nOpening dataset: {new_path}")

            if ctx is not None:
                # 1. Snapshot current layer data while layers still alive
                _app["snapshot"] = _snapshot_layers(ctx)

                # 2. Save session for current dataset
                zarr_path_old = ctx.data_path / "sdata_cached.zarr"
                if not ctx.no_cache and zarr_path_old.exists():
                    from palms.utils.adata_persistence import _persist_table
                    _persist_table(ctx)
                    ctx.state["segmentation_source"] = ctx.segmentation_source
                    from palms.utils.session import save_session
                    save_session(zarr_path_old, ctx.state, ctx.he_state, _app["snapshot"])
                    try:
                        from palms.utils.notebook_export import write_graph_notebook
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
                return False

            print(f"Dataset opened: {new_path.name}")
            return True

        finally:
            _app["reload_in_progress"] = False

    def _on_open_dataset():
        from qtpy.QtWidgets import QFileDialog, QMessageBox

        if _app["reload_in_progress"]:
            return
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
        _load_dataset_into_viewer(new_path)

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
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
        )
        if reply != QMessageBox.StandardButton.Ok:
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
        dlg.exec()  # blocks Qt event loop until dlg.accept() is called

    # ── File / View menus — added once to napari's native menus ─────────────
    create_file_menu(viewer, _on_open_dataset, _on_preprocess_dataset)
    create_view_menu(viewer, _app)

    # Let tabs reload the open dataset. The Cache tab needs it: recovered
    # elements are written straight into the zarr, so nothing in memory — the
    # sdata, the layers, the clustering combo boxes — knows they exist until
    # everything is rebuilt from disk. Rebound on every reload because
    # _do_full_init returns a fresh ViewerContext.
    def _reload_current_dataset() -> bool:
        current = ctx.data_path if ctx is not None else None
        if current is None:
            return False
        return _load_dataset_into_viewer(Path(current))

    # Registered before the first _do_full_init below, so every ViewerContext
    # ever built picks it up from _app.
    _app["reload_current_dataset"] = _reload_current_dataset

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

    # After the load, not before: every bridge call marshals onto the Qt main
    # thread, so a bridge started ahead of _do_full_init would accept requests
    # it could not service until the dataset had finished loading anyway.
    _mcp_server = dev_mcp.start_bridge(viewer, mcp_port) if mcp_port else None

    napari.run()

    # ── Save session state on exit ────────────────────────────────────────
    if ctx is not None and not no_cache:
        final_zarr_path = ctx.data_path / "sdata_cached.zarr"
        if final_zarr_path.exists():
            from palms.utils.adata_persistence import _persist_table, save_rois_to_sdata
            _persist_table(ctx)
            roi_data = _app["snapshot"].get("roi_data", [])
            save_rois_to_sdata(ctx, roi_data)
            ctx.state["segmentation_source"] = ctx.segmentation_source
            from palms.utils.session import save_session
            save_session(final_zarr_path, ctx.state, ctx.he_state, _app["snapshot"])
            try:
                from palms.utils.notebook_export import write_graph_notebook
                _g = ctx.state.get("prov_graph")
                if _g is not None and len(_g):
                    write_graph_notebook(_g, ctx.data_path / "analysis_notebook.ipynb")
            except Exception:
                pass
            print("Session saved to zarr cache.")


def main():
    """Console-script entry point — parses argv and launches the viewer."""
    data_path, no_cache, mcp_port = _parse_args()
    from palms.loader import CacheLoadAborted, NoRawSourceError
    try:
        run_viewer(data_path, no_cache=no_cache, mcp_port=mcp_port)
    except CacheLoadAborted as e:
        # The user declined to rebuild, or we had no way to ask. Exiting
        # quietly is the point: the cache is left exactly as it was.
        print(f"\n{e}")
        raise SystemExit(1) from None
    except NoRawSourceError as e:
        # The message names the dataset, why it cannot be rebuilt and what to
        # try instead; a traceback on top of that adds nothing.
        print(f"\n{e}")
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
