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

import os
import sys
import json
import time
import argparse
import warnings
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

# ─── Path setup ────────────────────────────────────────────────────────────
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.coloring import CellColorManager
from utils.transcript_index import TranscriptLoader
from utils.umap_widget import UMAPViewer
from utils.viewer_context import ViewerContext
from tabs._helpers import create_shared_helpers, create_preferences_menu

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
            print("No directory selected. Exiting.")
            sys.exit(0)
        data_path = Path(data_path_str)

    # Validate
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


def _import_loader():
    """Import 01_load_sdata as a module regardless of the numeric prefix."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "load_sdata", SCRIPTS_DIR / "01_load_sdata.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
    from tabs.tab_clustering import build_tab as build_clustering_tab
    from tabs.tab_cell_coloring import build_tab as build_cell_coloring_tab
    from tabs.tab_transcripts import build_tab as build_transcripts_tab
    from tabs.tab_umap import build_tab as build_umap_tab
    from tabs.tab_roi import build_tab as build_roi_tab
    from tabs.tab_he_registration import build_tab as build_he_registration_tab
    from tabs.tab_gene_analysis import build_tab as build_gene_analysis_tab
    from tabs.tab_ligrec import build_tab as build_ligrec_tab
    from tabs.tab_nhood import build_tab as build_nhood_tab
    from tabs.tab_co_occurrence import build_tab as build_co_occurrence_tab
    from tabs.tab_arms import build_tab as build_arms_tab

    # ── Build Cell Coloring first (creates cross-tab widgets) ────────────
    coloring_widget, coloring_exports = build_cell_coloring_tab(ctx)
    # ctx.clustering_widget etc. are now set

    # ── Build analysis tabs (they register their clustering widgets on ctx) ─
    ga_widget, ga_exports = build_gene_analysis_tab(ctx)
    lr_widget, lr_exports = build_ligrec_tab(ctx)
    nhood_widget, nhood_exports = build_nhood_tab(ctx)
    co_widget, co_exports = build_co_occurrence_tab(ctx)

    # ── Create shared helpers (needs all widgets registered) ─────────────
    create_shared_helpers(ctx)

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

    # ── Mouse hover: show cluster ID in status bar ───────────────────────
    if ctx.cell_labels_layer is not None:
        def _on_cursor_move(event):
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
                    text = f"Cell {int(label_val)} \u2014 {name}: {display_cid}"
                else:
                    text = f"Cell {int(label_val)} \u2014 {name}: unassigned"
                QTimer.singleShot(0, lambda t=text: setattr(ctx.viewer, 'status', t))

        ctx.viewer.cursor.events.position.connect(_on_cursor_move)

    # ── Assemble tabbed control panel ────────────────────────────────────
    tab_widget = QTabWidget()
    tab_widget.addTab(clustering_widget, "Clustering")
    tab_widget.addTab(coloring_widget, "Cell Coloring")
    tab_widget.addTab(transcripts_widget, "Transcripts")
    tab_widget.addTab(umap_widget, "UMAP")
    tab_widget.addTab(roi_widget, "ROI Analysis")
    tab_widget.addTab(he_widget, "H&E Registration")
    tab_widget.addTab(ga_widget, "Gene Analysis")
    tab_widget.addTab(lr_widget, "Ligand-Receptor")
    tab_widget.addTab(nhood_widget, "Nhood Enrichment")
    tab_widget.addTab(co_widget, "Co-occurrence")
    tab_widget.addTab(arms_widget, "ARMS Overlay")

    # ── Compose session restore from per-tab restorers ───────────────────
    all_exports = [
        clustering_exports, coloring_exports, transcripts_exports,
        umap_exports, roi_exports, he_exports, ga_exports,
        lr_exports, nhood_exports, co_exports, arms_exports,
    ]

    def restore_session(session):
        for exports in all_exports:
            restorer = exports.get("restore_session")
            if restorer:
                restorer(session)

    return tab_widget, state, he_state, restore_session


def main(data_path: Path, no_cache: bool = False):
    print("=" * 60)
    print("Xenium Linux Viewer")
    print(f"Dataset: {data_path}")
    print("=" * 60)

    t_start = time.perf_counter()

    # ── Read pixel size ──────────────────────────────────────────────────────
    pixel_size = _read_pixel_size(data_path)
    print(f"Pixel size: {pixel_size} um/px")

    # ── Load data ────────────────────────────────────────────────────────────
    loader_mod = _import_loader()

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
    gene_names = list(adata.var_names)
    clustering_names = list(clusterings.keys())

    print(f"Genes: {len(gene_names)}, Clusterings: {len(clustering_names)}")

    # ── Managers ─────────────────────────────────────────────────────────────
    color_manager = CellColorManager(adata, label_to_obs)
    transcript_loader = TranscriptLoader(
        cache_dir=data_path / "transcript_cache",
        parquet_path=data_path / "transcripts.parquet",
        pixel_size=pixel_size,
    )

    # ── Napari viewer (must be created before any QWidgets) ──────────────────
    t0 = time.perf_counter()
    print("Opening napari...")
    viewer = napari.Viewer(title=f"Xenium Viewer — {data_path.name}")

    # ── UMAP viewer (separate napari window, deferred until first coloring) ──
    cell_ids = adata.obs['cell_id'].values
    umap_viewer = UMAPViewer(umap_df, cell_ids)
    print(f"  Viewer + UMAP viewer created in {time.perf_counter() - t0:.1f}s")

    # ── Add layers from sdata ─────────────────────────────────────────────────
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

    # ── Build ViewerContext ──────────────────────────────────────────────────
    ctx = ViewerContext(
        viewer=viewer,
        adata=adata,
        sdata=sdata,
        clusterings=clusterings,
        color_manager=color_manager,
        transcript_loader=transcript_loader,
        umap_viewer=umap_viewer,
        label_to_obs=label_to_obs,
        centroids_yx=centroids_yx,
        pixel_size=pixel_size,
        data_path=data_path,
        no_cache=no_cache,
        gene_names=gene_names,
        clustering_names=clustering_names,
        cell_labels_layer=cell_labels_layer,
        transcript_layer=transcript_layer,
        roi_layer=roi_layer,
        morph_thumb=morph_thumb,
        morph_full_shape_yx=morph_full_shape_yx,
    )

    # Initialize mutable state dicts
    ctx.state = {
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
        "plot_format": "png",
        "plot_font_size": 10,
        "record_code": True,
        "code_journal": [],
        "code_journal_tags": set(),
        "custom_clusterings": {},
    }

    ctx.he_state = {
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

    ctx.arms_state = {
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

    # ── Control panel ────────────────────────────────────────────────────────
    panel, _state, _he_state, restore_session = _build_control_panel(ctx)
    viewer.window.add_dock_widget(panel, name="Xenium Controls", area="right")

    # ── Restore saved session ────────────────────────────────────────────────
    zarr_path = data_path / "sdata_cached.zarr"
    if not no_cache and zarr_path.exists():
        from utils.session import load_session
        session = load_session(zarr_path)
        if session is not None:
            print("Restoring session from zarr cache...")
            restore_session(session)
            print("Session restored.")

    # ── Snapshot layer data before Qt teardown, then save on exit ────────────
    _snapshot = {}  # filled by the about_to_close handler
    _arms_state = ctx.arms_state

    if not no_cache:
        def _on_viewer_closing(_event=None):
            """Capture napari layer data while Qt objects are still alive."""
            # ROI polygons
            _snapshot["roi_data"] = [
                np.asarray(p, dtype=np.float64) for p in roi_layer.data
            ] if roi_layer is not None else []

            # Landmark layer data
            xen_lm = _he_state.get("xenium_lm_layer")
            he_lm = _he_state.get("he_lm_layer")
            _snapshot["xenium_landmarks"] = (
                np.asarray(xen_lm.data, dtype=np.float64) if xen_lm is not None and len(xen_lm.data) > 0 else None
            )
            _snapshot["he_landmarks"] = (
                np.asarray(he_lm.data, dtype=np.float64) if he_lm is not None and len(he_lm.data) > 0 else None
            )

            # ARMS overlay landmark data
            try:
                arms_xen_lm = _arms_state.get("xenium_lm_layer")
                arms_he_lm = _arms_state.get("he_lm_layer")
                _snapshot["arms_xenium_landmarks"] = (
                    np.asarray(arms_xen_lm.data, dtype=np.float64)
                    if arms_xen_lm is not None and len(arms_xen_lm.data) > 0 else None
                )
                _snapshot["arms_he_landmarks"] = (
                    np.asarray(arms_he_lm.data, dtype=np.float64)
                    if arms_he_lm is not None and len(arms_he_lm.data) > 0 else None
                )
                # ARMS state (exclude non-serializable napari layers)
                _snapshot["arms_state"] = {
                    k: v for k, v in _arms_state.items()
                    if k not in ("he_layer", "he_tif", "xenium_lm_layer", "he_lm_layer", "shapes_layer")
                }
            except NameError:
                _snapshot["arms_xenium_landmarks"] = None
                _snapshot["arms_he_landmarks"] = None
                _snapshot["arms_state"] = {}

        from qtpy.QtWidgets import QApplication
        QApplication.instance().aboutToQuit.connect(_on_viewer_closing)

    total_time = time.perf_counter() - t_start
    print(f"\nViewer ready in {total_time:.1f}s. Close the napari window to exit.")
    napari.run()

    # ── Save session state on exit ────────────────────────────────────────
    if not no_cache and zarr_path.exists():
        from utils.session import save_session
        save_session(zarr_path, _state, _he_state, _snapshot)
        print("Session saved to zarr cache.")


if __name__ == "__main__":
    data_path, no_cache = _parse_args()
    main(data_path, no_cache=no_cache)
