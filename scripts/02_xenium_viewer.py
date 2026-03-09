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
from napari.qt.threading import thread_worker

# Suppress non-critical warnings from spatialdata stack
warnings.filterwarnings("ignore", category=UserWarning, module="spatialdata")
warnings.filterwarnings("ignore", category=FutureWarning)

# ─── Path setup ────────────────────────────────────────────────────────────
SCRIPTS_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.coloring import (
    CellColorManager, AVAILABLE_COLORMAPS, CLUSTER_PALETTE, TRANSCRIPT_PALETTE,
)
from utils.transcript_index import TranscriptLoader
from utils.umap_widget import UMAPViewer
from utils.registration import (
    load_he_pyramid, compute_landmark_affine, save_landmarks, load_landmarks,
    extract_tissue_mask_fluorescence, extract_tissue_mask_he, compute_coarse_affine,
)
from utils.gene_analysis import (
    get_normalized_adata, add_clustering_to_obs, run_rank_genes,
    make_rank_genes_dotplot, make_rank_genes_plot, compute_roi_deg,
)
from utils.spatial_analysis import (
    compute_spatial_neighbors, run_ligrec, make_ligrec_plot,
    run_nhood_enrichment, make_nhood_enrichment_plot,
    run_co_occurrence, make_co_occurrence_plot,
)

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

    # ── Control panel ────────────────────────────────────────────────────────
    panel, _state, _he_state, restore_session = _build_control_panel(
        viewer=viewer,
        gene_names=gene_names,
        clustering_names=clustering_names,
        clusterings=clusterings,
        color_manager=color_manager,
        transcript_loader=transcript_loader,
        cell_labels_layer=cell_labels_layer,
        transcript_layer=transcript_layer,
        umap_viewer=umap_viewer,
        label_to_obs=label_to_obs,
        roi_layer=roi_layer,
        centroids_yx=centroids_yx,
        pixel_size=pixel_size,
        data_path=data_path,
        morph_thumb=morph_thumb,
        morph_full_shape_yx=morph_full_shape_yx,
        adata=adata,
        sdata=sdata,
        no_cache=no_cache,
    )
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

        from qtpy.QtWidgets import QApplication
        QApplication.instance().aboutToQuit.connect(_on_viewer_closing)

    total_time = time.perf_counter() - t_start
    print(f"\nViewer ready in {total_time:.1f}s. Close the napari window to exit.")
    napari.run()

    # ── Auto-save reproducible code ──────────────────────────────────────────
    if _state.get("record_code") and _state["code_journal"]:
        code_path = data_path / "code.py"
        with open(code_path, 'w') as f:
            f.write("\n".join(_state["code_journal"]) + "\n")
        print(f"Reproducible code saved to {code_path}")

    # ── Save session state on exit ────────────────────────────────────────────
    if not no_cache and zarr_path.exists():
        from utils.session import save_session
        save_session(zarr_path, _state, _he_state, _snapshot)
        print("Session saved to zarr cache.")


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

    # ── Labels (multiscale pyramid from sdata) ────────────────────────────────
    for key in ["cell_labels", "nucleus_labels"]:
        if key in sdata.labels:
            print(f"  Adding {key} (multiscale)...")
            scales = _extract_dt_scales(sdata.labels[key])
            if scales:
                viewer.add_labels(scales, name=key)
            else:
                print(f"  Warning: could not extract {key} scales")


def _build_control_panel(
    viewer,
    gene_names: list,
    clustering_names: list,
    clusterings: dict,
    color_manager: CellColorManager,
    transcript_loader: TranscriptLoader,
    cell_labels_layer,
    transcript_layer,
    umap_viewer: UMAPViewer,
    label_to_obs: np.ndarray,
    roi_layer=None,
    centroids_yx: np.ndarray = None,
    pixel_size: float = 1.0,
    data_path: Path = None,
    morph_thumb=None,
    morph_full_shape_yx=None,
    adata=None,
    sdata=None,
    no_cache: bool = False,
):
    """Build and return a magicgui Container docked widget."""
    from magicgui.widgets import (
        Container, ComboBox, CheckBox, PushButton, Label,
        Slider, RadioButtons, FloatSpinBox,
    )
    from qtpy.QtWidgets import (
        QListWidget, QHBoxLayout, QWidget, QVBoxLayout, QLabel, QTabWidget,
        QCheckBox, QScrollArea, QGridLayout, QGroupBox,
    )

    # ── State ────────────────────────────────────────────────────────────────
    _state = {
        "current_gene": gene_names[0] if gene_names else None,
        "current_clustering": clustering_names[0] if clustering_names else None,
        "current_colormap": "viridis",
        "show_transcripts": False,
        "min_qv": 20,
        "color_mode": "Gene Expression",
        "transcript_genes": [],  # list of active transcript gene names
        "filter_by_cluster": False,
        "label_to_cluster": None,  # np.ndarray: label_value -> cluster_id (or -1)
        "active_clustering_name": None,  # name shown on hover
        "nhood_result": None,  # dict from run_nhood_enrichment()
        "nhood_fig": None,  # matplotlib Figure
        "co_result": None,  # dict from run_co_occurrence()
        "co_fig": None,  # matplotlib Figure
        "plot_format": "png",  # "png" or "svg" — set via Preferences menu
        "plot_font_size": 10,  # matplotlib font.size — set via Preferences menu
        "record_code": True,        # toggle via Preferences
        "code_journal": [],         # list of code block strings
        "code_journal_tags": set(), # dedup tags for preamble sections
        "custom_clusterings": {},   # custom clusterings to persist (Leiden, imported)
    }

    # ── Preferences menu (plot format) ────────────────────────────────────────
    from qtpy.QtWidgets import QActionGroup, QMenu
    from qtpy.QtGui import QAction

    menu_bar = viewer.window._qt_window.menuBar()
    prefs_menu = QMenu("Preferences", menu_bar)
    menu_bar.addMenu(prefs_menu)

    format_menu = prefs_menu.addMenu("Plot format")
    format_group = QActionGroup(format_menu)
    format_group.setExclusive(True)

    png_action = QAction("PNG", format_group, checkable=True, checked=True)
    svg_action = QAction("SVG", format_group, checkable=True)
    format_menu.addAction(png_action)
    format_menu.addAction(svg_action)

    def _on_format_changed(action):
        _state["plot_format"] = action.text().lower()

    format_group.triggered.connect(_on_format_changed)

    fontsize_menu = prefs_menu.addMenu("Plot font size")
    fontsize_group = QActionGroup(fontsize_menu)
    fontsize_group.setExclusive(True)

    for sz in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20):
        act = QAction(str(sz), fontsize_group, checkable=True,
                      checked=(sz == 10))
        fontsize_menu.addAction(act)

    def _on_fontsize_changed(action):
        _state["plot_font_size"] = int(action.text())

    fontsize_group.triggered.connect(_on_fontsize_changed)

    # --- Record code checkbox ---
    record_action = QAction("Record reproducible code", prefs_menu, checkable=True, checked=True)
    prefs_menu.addAction(record_action)

    def _on_record_toggled(checked):
        _state["record_code"] = checked
        if checked:
            _state["code_journal"].clear()
            _state["code_journal_tags"].clear()

    record_action.toggled.connect(_on_record_toggled)

    # --- Save code action ---
    save_code_action = QAction("Save recorded code...", prefs_menu)
    prefs_menu.addAction(save_code_action)

    def _on_save_code():
        if not _state["code_journal"]:
            return
        from qtpy.QtWidgets import QFileDialog as _QFD
        path, _ = _QFD.getSaveFileName(
            None, "Save Reproducible Code", "analysis.py", "Python Files (*.py)",
        )
        if path:
            with open(path, 'w') as f:
                f.write("\n".join(_state["code_journal"]) + "\n")

    save_code_action.triggered.connect(_on_save_code)

    def _apply_plot_font_size():
        import matplotlib.pyplot as _plt
        _plt.rcParams['font.size'] = _state.get("plot_font_size", 10)

    # ── Reproducible code journal helpers ──────────────────────────────────
    def _record_code(code: str, tag: str = None):
        """Append a code block to the journal. If tag is given, skip if already emitted."""
        if not _state.get("record_code"):
            return
        if tag:
            if tag in _state["code_journal_tags"]:
                return
            _state["code_journal_tags"].add(tag)
        _state["code_journal"].append(code)

    def _record_preamble():
        """Emit imports + data loading (once)."""
        _record_code(
            "import scanpy as sc\n"
            "import squidpy as sq\n"
            "import matplotlib.pyplot as plt\n"
            "import pandas as pd\n"
            "import numpy as np\n"
            f"\nplt.rcParams['font.size'] = {_state.get('plot_font_size', 10)}\n"
            f"\n# Load data\n"
            "from spatialdata_io import xenium\n"
            f"sdata = xenium(\"{data_path}\")\n"
            "adata = sdata[\"table\"].copy()",
            tag="preamble"
        )

    def _record_normalize():
        """Emit normalization (once)."""
        _record_preamble()
        _record_code(
            "\n# Normalize, log-transform, PCA\n"
            "sc.pp.normalize_total(adata)\n"
            "sc.pp.log1p(adata)\n"
            "sc.pp.pca(adata)",
            tag="normalize"
        )

    def _record_clustering(key):
        """Emit clustering assignment (once per key)."""
        _record_normalize()
        dir_name = f"gene_expression_{key}"
        csv_path = os.path.join(data_path, "analysis", "clustering", dir_name, "clusters.csv")
        if not os.path.exists(csv_path):
            csv_path = os.path.join(data_path, "analysis", "clustering", key, "clusters.csv")
        _record_code(
            f"\n# Add clustering: {key}\n"
            f"clust_df = pd.read_csv(\"{csv_path}\", index_col=0)\n"
            f"adata.obs[\"{key}\"] = pd.Categorical("
            f"clust_df.reindex(adata.obs_names).iloc[:, 0].astype(str).values)",
            tag=f"clustering_{key}"
        )

    def _record_spatial_neighbors(n_neighs):
        """Emit spatial neighbor computation (once per n_neighs)."""
        _record_code(
            f"\n# Compute spatial neighbors (k={n_neighs})\n"
            "adata.obsm['spatial'] = adata.obsm.get('spatial', "
            "np.column_stack([adata.obs['x_centroid'], adata.obs['y_centroid']]))\n"
            f"sq.gr.spatial_neighbors(adata, n_neighs={n_neighs}, coord_type=\"generic\")",
            tag=f"spatial_neighbors_{n_neighs}"
        )

    # ── Cell Coloring widgets ────────────────────────────────────────────────
    mode_widget = RadioButtons(
        label="Color cells by",
        choices=["Gene Expression", "Cluster"],
        value="Gene Expression",
    )

    gene_widget = ComboBox(
        label="Gene",
        choices=gene_names,
        value=gene_names[0] if gene_names else None,
    )

    colormap_widget = ComboBox(
        label="Colormap",
        choices=AVAILABLE_COLORMAPS,
        value="viridis",
    )

    clustering_widget = ComboBox(
        label="Clustering",
        choices=clustering_names,
        value=clustering_names[0] if clustering_names else None,
        enabled=False,
    )

    # ── Cluster filter widgets ───────────────────────────────────────────────
    filter_check = CheckBox(label="Filter by cluster", value=False, enabled=True)

    # Scrollable checkbox grid for multi-cluster selection
    cluster_filter_container = QWidget()
    cluster_filter_grid = QGridLayout()
    cluster_filter_grid.setContentsMargins(0, 0, 0, 0)
    cluster_filter_container.setLayout(cluster_filter_grid)

    cluster_scroll = QScrollArea()
    cluster_scroll.setWidget(cluster_filter_container)
    cluster_scroll.setWidgetResizable(True)
    cluster_scroll.setMaximumHeight(150)
    cluster_scroll.setEnabled(False)

    select_all_btn = PushButton(label="Select All", enabled=False)
    deselect_all_btn = PushButton(label="Deselect All", enabled=False)

    _state["cluster_checkboxes"] = {}  # int -> QCheckBox

    # ── White background toggle ─────────────────────────────────────────────
    bg_white_check = CheckBox(label="White background", value=False)

    apply_color_button = PushButton(label="Apply Cell Coloring", enabled=True)

    # ── Transcript Overlay widgets (multi-gene) ──────────────────────────────
    transcript_gene_widget = ComboBox(
        label="Transcript gene",
        choices=gene_names,
        value=gene_names[0] if gene_names else None,
    )

    add_gene_button = PushButton(label="Add Gene", enabled=True)
    remove_gene_button = PushButton(label="Remove Selected", enabled=True)
    clear_genes_button = PushButton(label="Clear All", enabled=True)

    # Build the gene list widget using Qt directly
    gene_list_qt = QListWidget()
    gene_list_qt.setMaximumHeight(150)

    # Legend label for gene↔color mapping
    legend_label_qt = QLabel("")
    legend_label_qt.setWordWrap(True)

    transcript_check = CheckBox(label="Show transcripts", value=False)
    qv_slider = Slider(label="Min QV", min=0, max=40, value=20)
    apply_transcripts_button = PushButton(label="Apply Transcripts", enabled=True)

    # ── UMAP widgets ────────────────────────────────────────────────────────
    show_umap_button = PushButton(label="Show UMAP Window", enabled=True)
    umap_size_slider = Slider(label="UMAP pt size", min=1, max=50, value=15)

    def on_show_umap():
        umap_viewer.show()

    def on_umap_size_change(value):
        umap_viewer.set_point_size(value / 100.0)

    show_umap_button.clicked.connect(on_show_umap)
    umap_size_slider.changed.connect(on_umap_size_change)

    # ── Status — all messages go to napari status bar ──────────────────────
    def _set_status(msg: str):
        viewer.status = msg

    # Keep a dummy object so `.value = ...` assignments still work everywhere
    class _StatusProxy:
        @property
        def value(self):
            return viewer.status
        @value.setter
        def value(self, msg):
            _set_status(msg)

    status_label = _StatusProxy()
    ga_status = _StatusProxy()
    roi_deg_status = _StatusProxy()
    lr_status = _StatusProxy()
    he_status_label = _StatusProxy()
    reg_status_label = _StatusProxy()

    # ── Helper: update legend text ────────────────────────────────────────────
    def _update_legend():
        genes = _state["transcript_genes"]
        if not genes:
            legend_label_qt.setText("")
            return
        parts = []
        color_names = [
            "Yellow", "Cyan", "Magenta", "Orange", "Green",
            "Sky Blue", "Red", "Violet", "Pink", "Brown",
        ]
        for i, g in enumerate(genes):
            parts.append(f"{color_names[i % len(color_names)]}: {g}")
        legend_label_qt.setText(" | ".join(parts))

    # ── Helper: repopulate cluster checkboxes ──────────────────────────────
    def _repopulate_cluster_checkboxes():
        # Clear existing checkboxes
        while cluster_filter_grid.count():
            item = cluster_filter_grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        _state["cluster_checkboxes"].clear()

        key = clustering_widget.value
        if not key or key not in clusterings:
            return
        # Handle both integer and string cluster IDs
        raw_ids = clusterings[key].dropna().unique().tolist()
        try:
            ids = sorted([int(x) for x in raw_ids])
        except (ValueError, TypeError):
            ids = sorted(raw_ids, key=lambda x: str(x))
        try:
            labels = _get_labels_for(key)
        except NameError:
            labels = {}
        cols = 3
        for i, cid in enumerate(ids):
            display = str(labels.get(cid, labels.get(str(cid), cid)))
            cb = QCheckBox(display)
            cb.setChecked(True)
            cb.setEnabled(filter_check.value)
            cluster_filter_grid.addWidget(cb, i // cols, i % cols)
            _state["cluster_checkboxes"][cid] = cb

    def _get_selected_cluster_ids():
        """Return set of cluster IDs whose checkboxes are checked."""
        return {cid for cid, cb in _state["cluster_checkboxes"].items() if cb.isChecked()}

    def _make_cluster_mask(aligned_values, selected_ids):
        """Build a boolean mask: True where aligned_values is in selected_ids.

        Handles both int and string cluster IDs.
        """
        sel = {str(s) for s in selected_ids}
        return np.array([str(v) in sel for v in aligned_values], dtype=bool)

    def _on_select_all():
        for cb in _state["cluster_checkboxes"].values():
            cb.setChecked(True)

    def _on_deselect_all():
        for cb in _state["cluster_checkboxes"].values():
            cb.setChecked(False)

    select_all_btn.clicked.connect(_on_select_all)
    deselect_all_btn.clicked.connect(_on_deselect_all)

    # Populate checkboxes for initial clustering
    _repopulate_cluster_checkboxes()

    # ── Leiden Clustering widgets ──────────────────────────────────────────
    from qtpy.QtWidgets import QTextEdit as _QTextEdit

    leiden_n_neighbors = Slider(label="n_neighbors", min=5, max=50, value=15)
    leiden_n_pcs = Slider(label="n_pcs", min=10, max=50, value=40)
    leiden_resolution = FloatSpinBox(label="resolution", min=0.1, max=5.0, step=0.1, value=1.0)
    leiden_hvg_check = CheckBox(label="Use HVGs only", value=False)
    leiden_n_hvgs = Slider(label="n_top_genes", min=500, max=4000, value=2000, enabled=False)
    leiden_scale_check = CheckBox(label="Scale (max_value=10)", value=False)

    def _on_hvg_toggle(val):
        leiden_n_hvgs.enabled = val
    leiden_hvg_check.changed.connect(_on_hvg_toggle)

    leiden_run_button = PushButton(label="Run Leiden Clustering", enabled=True)
    leiden_import_button = PushButton(label="Import Clustering...", enabled=True)
    leiden_export_button = PushButton(label="Export Clustering...", enabled=True)
    leiden_edit_labels_button = PushButton(label="Edit Cluster Labels...", enabled=True)

    leiden_status_text = _QTextEdit()
    leiden_status_text.setReadOnly(True)
    leiden_status_text.setFontFamily("monospace")
    leiden_status_text.setMaximumHeight(150)

    leiden_status = _StatusProxy()

    def _refresh_clustering_choices():
        """Update all clustering ComboBoxes with current clustering_names."""
        names = list(clusterings.keys())
        for combo in [clustering_widget, ga_clustering_widget, lr_clustering_widget,
                      ne_clustering_widget, co_clustering_widget]:
            old_val = combo.value
            combo.choices = names
            if old_val in names:
                combo.value = old_val

    def on_run_leiden():
        n_neighbors = leiden_n_neighbors.value
        n_pcs = leiden_n_pcs.value
        resolution = leiden_resolution.value
        use_hvg = leiden_hvg_check.value
        do_scale = leiden_scale_check.value
        n_hvgs = leiden_n_hvgs.value
        leiden_run_button.enabled = False
        leiden_status.value = "Running Leiden clustering..."

        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_leiden_ready})
        def _run():
            import scanpy as sc
            if use_hvg or do_scale:
                adata_work = _adata.copy()
                sc.pp.normalize_total(adata_work, target_sum=1e4)
                sc.pp.log1p(adata_work)
                if use_hvg:
                    sc.pp.highly_variable_genes(adata_work, n_top_genes=n_hvgs, flavor="seurat")
                    adata_work = adata_work[:, adata_work.var.highly_variable].copy()
                if do_scale:
                    sc.pp.scale(adata_work, max_value=10)
                sc.pp.pca(adata_work)
            else:
                adata_work = get_normalized_adata(_adata)
            sc.pp.neighbors(adata_work, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.leiden(adata_work, resolution=resolution, key_added="leiden")
            import pandas as pd
            cell_ids = _adata.obs['cell_id'].values if 'cell_id' in _adata.obs.columns else _adata.obs_names
            series = pd.Series(
                adata_work.obs['leiden'].astype(int).values,
                index=cell_ids,
                name="leiden",
            )
            n_clusters = series.nunique()
            return series, n_clusters, resolution, n_neighbors, n_pcs, use_hvg, do_scale, n_hvgs
        _run()

    def _on_leiden_ready(result):
        series, n_clusters, resolution, n_neighbors, n_pcs, use_hvg, do_scale, n_hvgs = result
        key = f"leiden_r{resolution}"
        clusterings[key] = series
        _state["custom_clusterings"][key] = series
        _refresh_clustering_choices()

        leiden_status_text.setPlainText(
            f"Leiden clustering complete\n"
            f"  Key: {key}\n"
            f"  Clusters: {n_clusters}\n"
            f"  n_neighbors: {n_neighbors}\n"
            f"  n_pcs: {n_pcs}\n"
            f"  resolution: {resolution}\n"
            f"  HVGs: {n_hvgs if use_hvg else 'all genes'}\n"
            f"  Scaled: {'yes (max=10)' if do_scale else 'no'}"
        )
        leiden_status.value = f"Leiden done: {n_clusters} clusters ({key})"
        leiden_run_button.enabled = True

        if use_hvg or do_scale:
            code_lines = [
                "\n# Custom preprocessing for Leiden",
                "adata_leiden = adata.copy()",
                "sc.pp.normalize_total(adata_leiden, target_sum=1e4)",
                "sc.pp.log1p(adata_leiden)",
            ]
            if use_hvg:
                code_lines.append(f"sc.pp.highly_variable_genes(adata_leiden, n_top_genes={n_hvgs}, flavor='seurat')")
                code_lines.append("adata_leiden = adata_leiden[:, adata_leiden.var.highly_variable].copy()")
            if do_scale:
                code_lines.append("sc.pp.scale(adata_leiden, max_value=10)")
            code_lines.append("sc.pp.pca(adata_leiden)")
            code_lines.append(f"sc.pp.neighbors(adata_leiden, n_neighbors={n_neighbors}, n_pcs={n_pcs})")
            code_lines.append(f'sc.tl.leiden(adata_leiden, resolution={resolution}, key_added="{key}")')
            _record_code("\n".join(code_lines), tag=f"leiden_{key}")
        else:
            _record_normalize()
            _record_code(
                f"\n# Leiden clustering (n_neighbors={n_neighbors}, n_pcs={n_pcs}, resolution={resolution})\n"
                f"sc.pp.neighbors(adata, n_neighbors={n_neighbors}, n_pcs={n_pcs})\n"
                f"sc.tl.leiden(adata, resolution={resolution}, key_added=\"{key}\")",
                tag=f"leiden_{key}"
            )

    leiden_run_button.clicked.connect(on_run_leiden)

    # ── Import / Export clustering callbacks ─────────────────────────────────
    def _on_import_clustering():
        import pandas as pd
        path, _ = QFileDialog.getOpenFileName(
            None, "Import Clustering", "",
            "CSV/TSV Files (*.csv *.tsv *.txt);;All Files (*)",
        )
        if not path:
            return
        df = pd.read_csv(path, sep=None, engine='python')
        if 'cell_id' in df.columns and 'group' in df.columns:
            series = pd.Series(df['group'].values, index=df['cell_id'].values)
        else:
            series = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0].values)
        name = Path(path).stem
        clusterings[name] = series
        _state["custom_clusterings"][name] = series
        _refresh_clustering_choices()
        clustering_widget.value = name
        leiden_status.value = f"Imported '{name}' ({series.nunique()} groups, {len(series)} cells)"

    def _on_export_clustering():
        import pandas as pd
        clustering_key = clustering_widget.value
        if not clustering_key or clustering_key not in clusterings:
            leiden_status.value = "No clustering selected"
            return
        series = clusterings[clustering_key]
        labels = _get_active_labels()
        if labels:
            mapped = series.map(lambda x: labels.get(x, labels.get(str(x),
                                labels.get(int(x) if str(x).lstrip('-').isdigit() else x, x))))
        else:
            mapped = series
        df = pd.DataFrame({'cell_id': mapped.index, 'group': mapped.values})
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Clustering", f"{clustering_key}.csv",
            "CSV Files (*.csv);;TSV Files (*.tsv)",
        )
        if not path:
            return
        sep = '\t' if path.endswith('.tsv') else ','
        df.to_csv(path, index=False, sep=sep)
        leiden_status.value = f"Exported {len(df)} cells to {path}"

    leiden_import_button.clicked.connect(_on_import_clustering)
    leiden_export_button.clicked.connect(_on_export_clustering)

    # ── Cell Coloring callbacks ─────────────────────────────────────────────
    def on_mode_change(value):
        _state["color_mode"] = value
        is_gene = (value == "Gene Expression")
        gene_widget.enabled = is_gene
        colormap_widget.enabled = is_gene
        clustering_widget.enabled = (value == "Cluster") or filter_check.value
        _set_cluster_filter_enabled(filter_check.value)

    def _set_cluster_filter_enabled(enabled):
        cluster_scroll.setEnabled(enabled)
        select_all_btn.enabled = enabled
        deselect_all_btn.enabled = enabled
        for cb in _state["cluster_checkboxes"].values():
            cb.setEnabled(enabled)

    def on_filter_change(value):
        _state["filter_by_cluster"] = value
        clustering_widget.enabled = (_state["color_mode"] == "Cluster") or value
        _set_cluster_filter_enabled(value)

    def on_clustering_change(value):
        _repopulate_cluster_checkboxes()
        # Sync all analysis tabs to the same clustering
        for combo in [ga_clustering_widget, lr_clustering_widget,
                      ne_clustering_widget, co_clustering_widget]:
            if value in [c for c in combo.choices]:
                combo.value = value

    def on_apply_color():
        if cell_labels_layer is None:
            status_label.value = "No cell_labels layer found"
            return

        mode = _state["color_mode"]
        status_label.value = "Computing cell colors..."
        apply_color_button.enabled = False

        if mode == "Gene Expression":
            gene = gene_widget.value
            cmap = colormap_widget.value
            _state["current_gene"] = gene
            _state["current_colormap"] = cmap

            use_filter = filter_check.value
            selected_ids = _get_selected_cluster_ids() if use_filter else None
            clustering_key = clustering_widget.value if use_filter else None

            @thread_worker(connect={"returned": _on_gene_colors_ready})
            def compute_gene():
                color_arr = color_manager.get_gene_colors(gene, colormap=cmap)
                return gene, color_arr, selected_ids, clustering_key
            compute_gene()

        else:  # Cluster
            clustering_key = clustering_widget.value
            _state["current_clustering"] = clustering_key
            cluster_series = clusterings[clustering_key]
            cluster_series.name = clustering_key

            use_filter = filter_check.value
            selected_ids = _get_selected_cluster_ids() if use_filter else None

            @thread_worker(connect={"returned": _on_cluster_colors_ready})
            def compute_cluster():
                return clustering_key, color_manager.get_cluster_colors(cluster_series), selected_ids

            compute_cluster()

    def _on_gene_colors_ready(result):
        gene, color_arr, selected_ids, clustering_key = result
        # If filtering by clusters, zero out non-selected cells
        if selected_ids is not None and clustering_key:
            _, label_to_cluster_arr = _get_cluster_ids_per_obs(clustering_key)
            int_ids = _translate_selected_ids_to_int(selected_ids)
            color_arr = color_arr.copy()
            mask_out = ~np.isin(label_to_cluster_arr, int_ids)
            valid_range = min(len(mask_out), len(color_arr))
            color_arr[:valid_range][mask_out[:valid_range]] = 0
            filter_desc = f" (clusters: {sorted(selected_ids)})"
        else:
            filter_desc = ""
        color_manager.apply_to_labels_layer(cell_labels_layer, color_arr)
        umap_viewer.color_by_gene(gene, color_arr, label_to_obs)
        # Clear cluster hover lookup (no longer showing clusters)
        _state["label_to_cluster"] = None
        _state["active_clustering_name"] = None
        status_label.value = f"Cells colored by gene: {gene}{filter_desc}"
        apply_color_button.enabled = True

    def _get_cluster_ids_per_obs(clustering_key):
        """Return per-obs cluster IDs aligned to adata, plus label lookup.

        For string-valued clusterings (e.g. imported), cluster IDs are encoded
        as integers via factorization. The mappings are stored in _state for
        hover display and filter translation.

        _state['_cluster_id_to_raw'] : dict[int, raw_id] or None
            Maps factorized int → original cluster value (for hover display).
        _state['_cluster_raw_to_id'] : dict[raw_id, int] or None
            Maps original cluster value → factorized int (for filter matching).
        """
        cluster_series = clusterings[clustering_key]
        adata = color_manager.adata
        if 'cell_id' in adata.obs.columns:
            cell_ids = adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids)
        else:
            clusters_aligned = cluster_series.reindex(adata.obs_names)

        # Try direct int conversion; if it fails, factorize string values
        import pandas as pd
        try:
            filled = clusters_aligned.fillna(-1)
            cluster_values = filled.values.astype(np.int32)
            _state['_cluster_id_to_raw'] = None
            _state['_cluster_raw_to_id'] = None
        except (ValueError, TypeError):
            codes, uniques = pd.factorize(clusters_aligned.values)
            cluster_values = codes.astype(np.int32)  # -1 for NaN stays -1
            _state['_cluster_id_to_raw'] = {int(i): u for i, u in enumerate(uniques)}
            _state['_cluster_raw_to_id'] = {u: int(i) for i, u in enumerate(uniques)}

        # Also build label -> cluster lookup for spatial hover
        max_label = len(label_to_obs) - 1
        label_to_cluster = np.full(max_label + 1, -1, dtype=np.int32)
        valid_mask = label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = label_to_obs[valid_labels]
        label_to_cluster[valid_labels] = cluster_values[obs_indices]
        return cluster_values, label_to_cluster

    def _translate_selected_ids_to_int(selected_ids):
        """Convert selected_ids (raw cluster IDs) to factorized ints if needed.

        When the clustering has string IDs, _state['_cluster_raw_to_id'] maps
        the raw values to the factorized integers used in label_to_cluster.
        """
        raw_to_id = _state.get('_cluster_raw_to_id')
        if raw_to_id is None:
            return list(selected_ids)
        return [raw_to_id[sid] for sid in selected_ids if sid in raw_to_id]

    def _on_cluster_colors_ready(result):
        clustering_key, (color_arr, cluster_to_color), selected_ids = result
        # Store cluster colors for use in analysis plots (e.g. co-occurrence)
        _state["cluster_to_color"] = cluster_to_color
        # Get per-obs cluster IDs for both UMAP hover and spatial hover
        cluster_ids_per_obs, label_to_cluster = _get_cluster_ids_per_obs(clustering_key)

        # If filtering by selected clusters, zero out all other cells
        if selected_ids is not None:
            int_ids = _translate_selected_ids_to_int(selected_ids)
            color_arr = color_arr.copy()
            mask_out = ~np.isin(label_to_cluster, int_ids)
            valid_range = min(len(mask_out), len(color_arr))
            color_arr[:valid_range][mask_out[:valid_range]] = 0

        color_manager.apply_to_labels_layer(cell_labels_layer, color_arr)
        umap_viewer.color_by_cluster(
            clustering_key, color_arr, label_to_obs,
            cluster_ids_per_obs=cluster_ids_per_obs,
        )
        _state["label_to_cluster"] = label_to_cluster
        _state["active_clustering_name"] = clustering_key
        filter_desc = f" (clusters: {sorted(selected_ids)})" if selected_ids is not None else ""
        status_label.value = f"Cells colored by cluster: {clustering_key}{filter_desc}"
        apply_color_button.enabled = True

    # ── Transcript Overlay callbacks (multi-gene) ────────────────────────────
    def on_add_gene():
        gene = transcript_gene_widget.value
        genes = _state["transcript_genes"]
        if gene in genes:
            status_label.value = f"'{gene}' already in list"
            return
        if len(genes) >= 10:
            status_label.value = "Maximum 10 genes reached"
            return
        genes.append(gene)
        gene_list_qt.addItem(gene)
        _update_legend()

    def on_remove_gene():
        selected = gene_list_qt.currentRow()
        if selected >= 0:
            gene_list_qt.takeItem(selected)
            _state["transcript_genes"].pop(selected)
            _update_legend()

    def on_clear_genes():
        _state["transcript_genes"].clear()
        gene_list_qt.clear()
        _update_legend()

    def on_apply_transcripts():
        genes = _state["transcript_genes"]
        if transcript_check.value and genes:
            status_label.value = f"Loading transcripts for {len(genes)} gene(s)..."
            apply_transcripts_button.enabled = False

            @thread_worker(connect={"returned": _on_transcripts_ready})
            def fetch():
                return transcript_loader.get_multi_gene_points(genes)

            fetch()
        elif transcript_check.value and not genes:
            status_label.value = "No genes in list — add genes first"
        else:
            transcript_layer.visible = False
            status_label.value = "Transcripts hidden"

    def _on_transcripts_ready(result):
        points, colors = result
        transcript_layer.data = points
        transcript_layer.face_color = colors
        transcript_layer.visible = True
        genes = _state["transcript_genes"]
        status_label.value = (
            f"Transcripts: {', '.join(genes)} ({len(points):,} spots total)"
        )
        apply_transcripts_button.enabled = True

    # ── Mouse hover: show cluster ID in status bar ─────────────────────────
    # napari's own layer status updates overwrite viewer.status on every
    # cursor move. We defer our update with QTimer.singleShot(0) so it
    # runs on the next event-loop tick, after napari's updates finish.
    from qtpy.QtCore import QTimer

    if cell_labels_layer is not None:
        def _on_cursor_move(event):
            lut = _state["label_to_cluster"]
            if lut is None:
                return  # not in cluster mode
            label_val = cell_labels_layer.get_value(
                viewer.cursor.position,
                view_direction=None,
                dims_displayed=list(range(cell_labels_layer.ndim)),
                world=True,
            )
            # multiscale labels return (data_level, label_value)
            if isinstance(label_val, tuple):
                label_val = label_val[1]
            if label_val is not None and 0 < int(label_val) < len(lut):
                cid = lut[int(label_val)]
                name = _state["active_clustering_name"] or "cluster"
                if cid >= 0:
                    # Map factorized int back to raw cluster name if needed
                    raw_map = _state.get('_cluster_id_to_raw')
                    raw_cid = raw_map[cid] if raw_map and cid in raw_map else cid
                    labels = _get_active_labels()
                    display_cid = labels.get(raw_cid, labels.get(str(raw_cid), raw_cid))
                    text = f"Cell {int(label_val)} \u2014 {name}: {display_cid}"
                else:
                    text = f"Cell {int(label_val)} \u2014 {name}: unassigned"
                QTimer.singleShot(0, lambda t=text: setattr(viewer, 'status', t))

        viewer.cursor.events.position.connect(_on_cursor_move)

    # ── ROI Analysis widgets ──────────────────────────────────────────────────
    from qtpy.QtWidgets import QTextEdit, QFileDialog

    def _get_plot_save_path(title: str, default_stem: str) -> str | None:
        """Open a save dialog using the preferred plot format (png/svg)."""
        fmt = _state.get("plot_format", "png")
        filter_str = f"{fmt.upper()} Files (*.{fmt});;All Files (*)"
        path, _ = QFileDialog.getSaveFileName(
            None, title, f"{default_stem}.{fmt}", filter_str,
        )
        return path if path else None

    roi_calc_button = PushButton(label="Calculate Expression", enabled=True)
    roi_export_button = PushButton(label="Export CSV", enabled=False)
    roi_text = QTextEdit()
    roi_text.setReadOnly(True)
    roi_text.setFontFamily("monospace")
    roi_text.setMaximumHeight(300)

    def on_calculate_roi():
        from shapely.geometry import Polygon as ShapelyPolygon
        from shapely import contains_xy

        gene = gene_widget.value
        if gene is None:
            roi_text.setPlainText("No gene selected.")
            return

        polygons = roi_layer.data if roi_layer is not None else []
        if len(polygons) == 0:
            roi_text.setPlainText("No ROI polygons drawn.\nUse the Shapes layer to draw polygons.")
            return

        # Get expression for current gene
        adata = color_manager.adata
        gene_idx = adata.var_names.get_loc(gene)
        X = adata.X
        if hasattr(X, "toarray"):
            expr = np.asarray(X[:, gene_idx].toarray()).ravel().astype(np.float32)
        else:
            expr = np.asarray(X[:, gene_idx]).ravel().astype(np.float32)

        # Build cluster mask if filter is active
        use_filter = filter_check.value
        cluster_mask = None  # None = no filtering, else bool array (n_obs,)
        filter_desc = ""
        if use_filter:
            clustering_key = clustering_widget.value
            selected_ids = _get_selected_cluster_ids()
            cluster_series = clusterings[clustering_key]
            if 'cell_id' in adata.obs.columns:
                cell_ids_arr = adata.obs['cell_id'].values
                clusters_aligned = cluster_series.reindex(cell_ids_arr)
            else:
                clusters_aligned = cluster_series.reindex(adata.obs_names)
            cluster_mask = _make_cluster_mask(clusters_aligned.values, selected_ids)
            filter_desc = f" ({clustering_key} clusters: {sorted(selected_ids)})"

        from scipy import stats
        from itertools import combinations

        lines = [f"Gene: {gene}{filter_desc}", ""]
        roi_results = []  # list of dicts for export
        region_exprs = []  # list of (region_idx, expression_array) for t-tests

        for i, poly_yx in enumerate(polygons):
            # poly_yx is Nx2 in napari (y, x) coords; shapely needs (x, y)
            poly_xy = poly_yx[:, ::-1]
            shapely_poly = ShapelyPolygon(poly_xy)
            if not shapely_poly.is_valid:
                shapely_poly = shapely_poly.buffer(0)

            # Check which centroids are inside (intersect with cluster filter if active)
            inside = contains_xy(shapely_poly, centroids_yx[:, 1], centroids_yx[:, 0])
            if cluster_mask is not None:
                inside = inside & cluster_mask
            inside_idx = np.where(inside)[0]
            n_cells = len(inside_idx)

            if n_cells == 0:
                lines.append(f"Region {i+1}: 0 cells")
                region_exprs.append((i + 1, np.array([], dtype=np.float32)))
            else:
                region_expr = expr[inside_idx]
                lines.append(
                    f"Region {i+1}: {n_cells} cells, "
                    f"mean={region_expr.mean():.2f}, "
                    f"median={np.median(region_expr):.2f}, "
                    f"std={region_expr.std():.2f}, "
                    f"min={region_expr.min():.0f}, "
                    f"max={region_expr.max():.0f}"
                )
                region_exprs.append((i + 1, region_expr))

            # Store for export
            for idx in inside_idx:
                # centroids_yx is (y, x) in pixels; convert back to microns (x, y)
                x_um = centroids_yx[idx, 1] * pixel_size
                y_um = centroids_yx[idx, 0] * pixel_size
                cell_id = adata.obs['cell_id'].values[idx] if 'cell_id' in adata.obs.columns else str(idx)
                roi_results.append({
                    "region_id": i + 1,
                    "cell_id": cell_id,
                    "x_centroid_um": x_um,
                    "y_centroid_um": y_um,
                    "expression": expr[idx],
                })

        # ── Significance testing (Welch's t-test, pairwise) ──────────────
        # Only between regions with >= 2 cells
        testable = [(r, e) for r, e in region_exprs if len(e) >= 2]
        pairs = list(combinations(testable, 2))
        if pairs:
            lines.append("")
            lines.append("── Pairwise Welch's t-tests ──")
            raw_pvals = []
            pair_labels = []
            for (r1, e1), (r2, e2) in pairs:
                t_stat, p_val = stats.ttest_ind(e1, e2, equal_var=False)
                raw_pvals.append(p_val)
                pair_labels.append((r1, r2, t_stat, p_val))

            # Benjamini-Hochberg correction if >1 comparison
            n_tests = len(raw_pvals)
            if n_tests > 1:
                sorted_idx = np.argsort(raw_pvals)
                adjusted = np.empty(n_tests, dtype=np.float64)
                for rank_pos, orig_idx in enumerate(sorted_idx):
                    adjusted[orig_idx] = raw_pvals[orig_idx] * n_tests / (rank_pos + 1)
                # Enforce monotonicity (step-up) and cap at 1.0
                adjusted_sorted = adjusted[sorted_idx]
                for j in range(n_tests - 2, -1, -1):
                    adjusted_sorted[j] = min(adjusted_sorted[j], adjusted_sorted[j + 1])
                adjusted[sorted_idx] = adjusted_sorted
                adjusted = np.minimum(adjusted, 1.0)

                for k, (r1, r2, t_stat, p_raw) in enumerate(pair_labels):
                    p_adj = adjusted[k]
                    sig = " *" if p_adj < 0.05 else ""
                    lines.append(
                        f"  Region {r1} vs {r2}: t={t_stat:.3f}, "
                        f"p={p_raw:.2e}, p_adj(BH)={p_adj:.2e}{sig}"
                    )
                lines.append(f"  ({n_tests} comparisons, Benjamini-Hochberg correction)")
            else:
                r1, r2, t_stat, p_val = pair_labels[0]
                sig = " *" if p_val < 0.05 else ""
                lines.append(
                    f"  Region {r1} vs {r2}: t={t_stat:.3f}, p={p_val:.2e}{sig}"
                )

        roi_text.setPlainText("\n".join(lines))
        _state["roi_results"] = roi_results
        _state["roi_gene"] = gene
        roi_export_button.enabled = len(roi_results) > 0

    def on_export_csv():
        results = _state.get("roi_results", [])
        if not results:
            return
        import csv
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ROI Data", f"roi_{_state.get('roi_gene', 'gene')}.csv",
            "CSV Files (*.csv)",
        )
        if not path:
            return
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["region_id", "cell_id", "x_centroid_um", "y_centroid_um", "expression"])
            writer.writeheader()
            writer.writerows(results)
        status_label.value = f"Exported {len(results)} cells to {path}"

    roi_calc_button.clicked.connect(on_calculate_roi)
    roi_export_button.clicked.connect(on_export_csv)

    # ── ROI DEG widgets (extend ROI Analysis tab) ────────────────────────────
    from qtpy.QtWidgets import QLabel as QtLabel

    roi_deg_method_widget = ComboBox(
        label="DEG Method", choices=["wilcoxon", "t-test"], value="wilcoxon",
    )
    roi_deg_filter_check = CheckBox(label="Filter by cluster", value=False)
    roi_deg_button = PushButton(label="Run ROI DEG", enabled=True)
    roi_deg_text = QTextEdit()
    roi_deg_text.setReadOnly(True)
    roi_deg_text.setFontFamily("monospace")
    roi_deg_text.setMaximumHeight(250)
    roi_deg_export_button = PushButton(label="Export DEG CSV...", enabled=False)
    # roi_deg_status assigned above as _StatusProxy

    def on_roi_deg():
        polygons = roi_layer.data if roi_layer is not None else []
        if len(polygons) < 2:
            roi_deg_status.value = "Need at least 2 ROI polygons drawn"
            return

        roi_deg_status.value = "Running differential expression..."
        roi_deg_button.enabled = False

        use_filter = roi_deg_filter_check.value
        cluster_mask = None
        if use_filter:
            clustering_key = clustering_widget.value
            selected_ids = _get_selected_cluster_ids()
            cluster_series = clusterings[clustering_key]
            _adata = adata if adata is not None else color_manager.adata
            if 'cell_id' in _adata.obs.columns:
                cell_ids_arr = _adata.obs['cell_id'].values
                clusters_aligned = cluster_series.reindex(cell_ids_arr)
            else:
                clusters_aligned = cluster_series.reindex(_adata.obs_names)
            cluster_mask = _make_cluster_mask(clusters_aligned.values, selected_ids)

        method = roi_deg_method_widget.value
        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_roi_deg_ready})
        def _run():
            return compute_roi_deg(
                _adata, centroids_yx, polygons, pixel_size,
                cluster_mask=cluster_mask, method=method,
            )
        _run()

    def _on_roi_deg_ready(df):
        _state["roi_deg_df"] = df
        roi_deg_button.enabled = True
        if df.empty:
            roi_deg_text.setPlainText("No significant results or insufficient cells in ROIs.")
            roi_deg_status.value = "DEG: no results"
            roi_deg_export_button.enabled = False
            return
        # Show top 50 rows
        preview = df.head(50).to_string(index=False)
        roi_deg_text.setPlainText(preview)
        roi_deg_status.value = f"DEG complete: {len(df)} gene-group results"
        roi_deg_export_button.enabled = True

    def on_export_roi_deg():
        df = _state.get("roi_deg_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ROI DEG Results", "roi_deg_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        roi_deg_status.value = f"Exported {len(df)} rows to {path}"

    roi_deg_button.clicked.connect(on_roi_deg)
    roi_deg_export_button.clicked.connect(on_export_roi_deg)

    # ── Gene Analysis widgets (Tab 6) ─────────────────────────────────────────
    ga_clustering_widget = ComboBox(
        label="Clustering", choices=clustering_names,
        value=clustering_names[0] if clustering_names else None,
    )
    ga_method_widget = ComboBox(
        label="Method", choices=["wilcoxon", "t-test", "logreg"], value="wilcoxon",
    )
    ga_n_genes_slider = Slider(label="Top N genes", min=5, max=50, value=25)
    ga_run_button = PushButton(label="Run Rank Genes", enabled=True)

    ga_dotplot_n_slider = Slider(label="Genes per cluster", min=3, max=20, value=5)
    ga_dendro_check = CheckBox(label="Dendrogram", value=True)
    ga_dotplot_button = PushButton(label="Show Dotplot", enabled=False)
    ga_edit_labels_button = PushButton(label="Edit Cluster Labels...", enabled=False)
    ga_save_dotplot_button = PushButton(label="Save Dotplot as PNG...", enabled=False)

    ga_rank_plot_button = PushButton(label="Show Rank Genes Plot", enabled=False)

    ga_results_text = QTextEdit()
    ga_results_text.setReadOnly(True)
    ga_results_text.setFontFamily("monospace")
    ga_results_text.setMaximumHeight(300)
    ga_export_button = PushButton(label="Export Full Results CSV...", enabled=False)
    # ga_status assigned above as _StatusProxy

    def on_run_rank_genes():
        ga_status.value = "Running rank genes (normalizing + computing)..."
        ga_run_button.enabled = False

        clustering_key = ga_clustering_widget.value
        method = ga_method_widget.value
        n_genes = ga_n_genes_slider.value
        _state["_rg_method"] = method
        _state["_rg_n_genes"] = n_genes

        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_rank_genes_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, clusterings[clustering_key], clustering_key)
            df = run_rank_genes(adata_norm, clustering_key, method=method, n_genes=n_genes)
            return df, adata_norm, clustering_key
        _run()

    def _on_rank_genes_ready(result):
        df, adata_norm, clustering_key = result
        _state["rank_genes_df"] = df
        _state["rank_genes_adata_norm"] = adata_norm
        _state["rank_genes_groupby"] = clustering_key
        ga_run_button.enabled = True
        ga_dotplot_button.enabled = True
        ga_rank_plot_button.enabled = True
        ga_edit_labels_button.enabled = True
        ga_export_button.enabled = True
        # Show top 50 rows
        preview = df.head(50).to_string(index=False)
        ga_results_text.setPlainText(preview)
        ga_status.value = f"Rank genes done: {len(df)} results ({clustering_key}, {ga_method_widget.value})"

        # Record reproducible code
        _rg_method = _state.get("_rg_method", "wilcoxon")
        _rg_n = _state.get("_rg_n_genes", 25)
        _record_clustering(clustering_key)
        _record_code(
            f"\n# Rank genes: method={_rg_method}, groupby={clustering_key}, n_genes={_rg_n}\n"
            f"sc.tl.rank_genes_groups(adata, groupby=\"{clustering_key}\", "
            f"method=\"{_rg_method}\", n_genes={_rg_n})\n"
            f"rank_df = sc.get.rank_genes_groups_df(adata, group=None)"
        )

    def on_show_dotplot():
        adata_norm = _state.get("rank_genes_adata_norm")
        groupby = _state.get("rank_genes_groupby")
        if adata_norm is None or groupby is None:
            ga_status.value = "Run rank genes first"
            return
        ga_status.value = "Generating dotplot..."
        ga_dotplot_button.enabled = False

        n_genes = ga_dotplot_n_slider.value
        dendro = ga_dendro_check.value
        labels = _get_labels_for(groupby)

        @thread_worker(connect={"returned": _on_dotplot_ready})
        def _run():
            _apply_plot_font_size()
            fig = make_rank_genes_dotplot(
                adata_norm, groupby, n_genes=n_genes,
                cluster_labels=labels, dendrogram=dendro,
            )
            return fig
        _run()

    def _on_dotplot_ready(fig):
        _state["dotplot_fig"] = fig
        ga_dotplot_button.enabled = True
        ga_save_dotplot_button.enabled = True
        import matplotlib.pyplot as _plt
        _plt.show(block=False)
        ga_status.value = "Dotplot displayed"

        # Record reproducible code
        _dp_n = ga_dotplot_n_slider.value
        _dp_dendro = ga_dendro_check.value
        _dp_groupby = _state.get("rank_genes_groupby", "")
        _record_code(
            f"\n# Dotplot (n_genes={_dp_n}, dendrogram={_dp_dendro})\n"
            + (f"sc.tl.dendrogram(adata, groupby=\"{_dp_groupby}\")\n" if _dp_dendro else "")
            + f"sc.pl.rank_genes_groups_dotplot(adata, n_genes={_dp_n}, "
            f"dendrogram={_dp_dendro})\nplt.show()"
        )

    # ── Per-clustering label storage and helpers ───────────────────────────
    # _state["cluster_labels"] is a nested dict: {clustering_name: {cluster_id: label}}
    if "cluster_labels" not in _state or not isinstance(_state.get("cluster_labels"), dict):
        _state["cluster_labels"] = {}

    def _get_active_labels():
        """Get cluster_labels dict for the currently active clustering."""
        key = _state.get("active_clustering_name") or clustering_widget.value
        if key:
            all_labels = _state.get("cluster_labels", {})
            if isinstance(all_labels, dict):
                return all_labels.get(key, {})
        return {}

    def _get_labels_for(clustering_key):
        """Get cluster_labels dict for a specific clustering."""
        all_labels = _state.get("cluster_labels", {})
        if isinstance(all_labels, dict):
            return all_labels.get(clustering_key, {})
        return {}

    def _build_label_editor_dialog(clustering_key):
        """Open a multi-column label editor dialog for the given clustering.

        Returns True if labels were updated, False if cancelled.
        """
        from qtpy.QtWidgets import QDialog, QGridLayout, QLineEdit, QDialogButtonBox

        if not clustering_key or clustering_key not in clusterings:
            return False
        cluster_series = clusterings[clustering_key]
        # Handle both integer and string cluster IDs
        ids = sorted(cluster_series.dropna().unique().tolist(), key=lambda x: (str(x),))
        existing = _get_labels_for(clustering_key)

        dialog = QDialog()
        dialog.setWindowTitle(f"Edit Cluster Labels \u2014 {clustering_key}")

        # Multi-column layout: ~10 rows per column, up to 3 columns
        n_cols = min(3, max(1, (len(ids) + 9) // 10))
        n_per_col = (len(ids) + n_cols - 1) // n_cols

        outer_layout = QVBoxLayout()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_content = QWidget()
        grid = QGridLayout()

        edits = {}
        for i, cid in enumerate(ids):
            col_idx = i // n_per_col
            row_idx = i % n_per_col
            grid.addWidget(QtLabel(f"{cid}:"), row_idx, col_idx * 2)
            edit = QLineEdit(str(existing.get(cid, existing.get(str(cid), cid))))
            grid.addWidget(edit, row_idx, col_idx * 2 + 1)
            edits[cid] = edit

        scroll_content.setLayout(grid)
        scroll.setWidget(scroll_content)
        outer_layout.addWidget(scroll)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        outer_layout.addWidget(buttons)
        dialog.setLayout(outer_layout)
        dialog.resize(min(800, 250 * n_cols), min(600, 30 * n_per_col + 60))

        if dialog.exec_() == QDialog.Accepted:
            new_labels = {cid: e.text() for cid, e in edits.items()}
            if "cluster_labels" not in _state or not isinstance(_state["cluster_labels"], dict):
                _state["cluster_labels"] = {}
            _state["cluster_labels"][clustering_key] = new_labels
            return True
        return False

    def _open_label_editor():
        clustering_key = ga_clustering_widget.value
        if _build_label_editor_dialog(clustering_key):
            ga_status.value = f"Labels updated for {clustering_key}"

    def _open_clustering_label_editor():
        clustering_key = clustering_widget.value
        if not clustering_key or clustering_key not in clusterings:
            leiden_status.value = "No clustering selected"
            return
        if _build_label_editor_dialog(clustering_key):
            leiden_status.value = f"Labels updated for {clustering_key}"
            _repopulate_cluster_checkboxes()

    leiden_edit_labels_button.clicked.connect(_open_clustering_label_editor)

    def on_save_dotplot():
        fig = _state.get("dotplot_fig")
        if fig is None:
            return
        path = _get_plot_save_path("Save Dotplot", "dotplot")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        ga_status.value = f"Dotplot saved to {path}"

    def on_show_rank_plot():
        adata_norm = _state.get("rank_genes_adata_norm")
        if adata_norm is None:
            ga_status.value = "Run rank genes first"
            return
        import matplotlib.pyplot as _plt
        _apply_plot_font_size()
        _rp_n = ga_n_genes_slider.value
        groupby = _state.get("rank_genes_groupby", "")
        labels = _get_labels_for(groupby)
        fig = make_rank_genes_plot(adata_norm, n_genes=_rp_n, cluster_labels=labels)
        _plt.show(block=False)
        ga_status.value = "Rank genes plot displayed"

        # Record reproducible code
        _record_code(
            f"\n# Rank genes panel plot (n_genes={_rp_n})\n"
            f"sc.pl.rank_genes_groups(adata, n_genes={_rp_n})\nplt.show()"
        )

    def on_export_rank_genes():
        df = _state.get("rank_genes_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Rank Genes Results", "rank_genes_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        ga_status.value = f"Exported {len(df)} rows to {path}"

    ga_run_button.clicked.connect(on_run_rank_genes)
    ga_dotplot_button.clicked.connect(on_show_dotplot)
    ga_edit_labels_button.clicked.connect(_open_label_editor)
    ga_save_dotplot_button.clicked.connect(on_save_dotplot)
    ga_rank_plot_button.clicked.connect(on_show_rank_plot)
    ga_export_button.clicked.connect(on_export_rank_genes)

    # ── Ligand-Receptor widgets (Tab 7) ───────────────────────────────────────
    lr_clustering_widget = ComboBox(
        label="Clustering", choices=clustering_names,
        value=clustering_names[0] if clustering_names else None,
    )
    lr_perms_slider = Slider(label="Permutations", min=100, max=1000, value=1000)
    lr_neighs_slider = Slider(label="N neighbors", min=3, max=20, value=6)
    lr_run_button = PushButton(label="Run L-R Analysis", enabled=True)

    # ── Interaction database filter checkboxes ─────────────────────────────
    lr_ds_group = QGroupBox("Interaction datasets")
    lr_ds_layout = QGridLayout()
    lr_ds_omnipath = QCheckBox("OmniPath"); lr_ds_omnipath.setChecked(True)
    lr_ds_ligrecextra = QCheckBox("LigRecExtra"); lr_ds_ligrecextra.setChecked(True)
    lr_ds_pathwayextra = QCheckBox("PathwayExtra"); lr_ds_pathwayextra.setChecked(True)
    lr_ds_kinaseextra = QCheckBox("KinaseExtra"); lr_ds_kinaseextra.setChecked(True)
    lr_ds_layout.addWidget(lr_ds_omnipath, 0, 0)
    lr_ds_layout.addWidget(lr_ds_ligrecextra, 0, 1)
    lr_ds_layout.addWidget(lr_ds_pathwayextra, 1, 0)
    lr_ds_layout.addWidget(lr_ds_kinaseextra, 1, 1)
    lr_ds_group.setLayout(lr_ds_layout)
    lr_cpdb_only = QCheckBox("CellPhoneDB only")
    lr_cpdb_only.setChecked(False)

    # lr_status assigned above as _StatusProxy

    lr_results_text = QTextEdit()
    lr_results_text.setReadOnly(True)
    lr_results_text.setFontFamily("monospace")
    lr_results_text.setMaximumHeight(250)

    lr_pval_widget = ComboBox(
        label="P-value threshold", choices=["0.001", "0.005", "0.01", "0.05"],
        value="0.05",
    )
    lr_plot_button = PushButton(label="Show L-R Plot", enabled=False)
    lr_save_plot_button = PushButton(label="Save L-R Plot as PNG...", enabled=False)
    lr_export_means_button = PushButton(label="Export Means CSV...", enabled=False)
    lr_export_pvals_button = PushButton(label="Export P-values CSV...", enabled=False)

    def on_run_ligrec():
        lr_status.value = "Running L-R analysis... (this may take several minutes)"
        lr_run_button.enabled = False

        clustering_key = lr_clustering_widget.value
        n_perms = lr_perms_slider.value
        n_neighs = lr_neighs_slider.value
        _state["_lr_params"] = {
            "clustering_key": clustering_key,
            "n_perms": n_perms,
            "n_neighs": n_neighs,
        }
        _adata = adata if adata is not None else color_manager.adata

        # Build interactions_params from dataset checkboxes
        from omnipath.constants import InteractionDataset
        include = []
        if lr_ds_omnipath.isChecked():
            include.append(InteractionDataset.OMNIPATH)
        if lr_ds_ligrecextra.isChecked():
            include.append(InteractionDataset.LIGREC_EXTRA)
        if lr_ds_pathwayextra.isChecked():
            include.append(InteractionDataset.PATHWAY_EXTRA)
        if lr_ds_kinaseextra.isChecked():
            include.append(InteractionDataset.KINASE_EXTRA)
        interactions_params = {"include": tuple(include)} if include else {}
        if lr_cpdb_only.isChecked():
            interactions_params["resources"] = "CellPhoneDB"

        @thread_worker(connect={"returned": _on_ligrec_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, clusterings[clustering_key], clustering_key)
            # Set spatial coordinates for squidpy
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            compute_spatial_neighbors(adata_norm, n_neighs=n_neighs)
            result = run_ligrec(adata_norm, clustering_key, n_perms=n_perms,
                                interactions_params=interactions_params)
            return result
        _run()

    def _on_ligrec_ready(result):
        _state["ligrec_result"] = result
        lr_run_button.enabled = True

        warning = result.get('warning')
        means = result['means']
        pvalues = result['pvalues']

        if warning:
            lr_status.value = f"L-R: {warning}"
            lr_plot_button.enabled = False
            lr_save_plot_button.enabled = False
            lr_export_means_button.enabled = False
            lr_export_pvals_button.enabled = False
            return

        # Summary
        n_interactions = means.shape[0]
        pval_thresh = float(lr_pval_widget.value)
        n_sig = (pvalues < pval_thresh).sum().sum() if not pvalues.empty else 0

        lines = [
            f"L-R interactions found: {n_interactions}",
            f"Significant (p < {pval_thresh}): {n_sig}",
            "",
        ]
        if not means.empty:
            lines.append("Top interactions by mean expression:")
            # Show top 20 by max mean across cluster pairs
            top_means = means.max(axis=1).sort_values(ascending=False).head(20)
            for idx, val in top_means.items():
                lines.append(f"  {idx}: {val:.4f}")

        lr_results_text.setPlainText("\n".join(lines))
        lr_status.value = f"L-R done: {n_interactions} interactions, {n_sig} significant"
        lr_plot_button.enabled = n_interactions > 0
        lr_save_plot_button.enabled = False
        lr_export_means_button.enabled = not means.empty
        lr_export_pvals_button.enabled = not pvalues.empty

        # Record reproducible code
        _lr_p = _state.get("_lr_params", {})
        _lr_ck = _lr_p.get("clustering_key", "")
        _lr_np = _lr_p.get("n_perms", 1000)
        _lr_nn = _lr_p.get("n_neighs", 6)
        _record_clustering(_lr_ck)
        _record_spatial_neighbors(_lr_nn)
        _record_code(
            f"\n# Ligand-receptor analysis (n_perms={_lr_np})\n"
            f"sq.gr.ligrec(\n"
            f"    adata, cluster_key=\"{_lr_ck}\", n_perms={_lr_np},\n"
            f"    threshold=0.01, seed=42,\n"
            f"    transmitter_params={{\"categories\": \"ligand\"}},\n"
            f"    receiver_params={{\"categories\": \"receptor\"}},\n"
            f")"
        )

    def _get_cluster_filter():
        """Return list of cluster ID strings based on cell coloring filter, or None."""
        if not _state.get("filter_by_cluster"):
            return None
        selected = _get_selected_cluster_ids()
        if not selected:
            return None
        return sorted(str(cid) for cid in selected)

    def on_show_lr_plot():
        result = _state.get("ligrec_result")
        if result is None:
            return
        pval_thresh = float(lr_pval_widget.value)
        groups = _get_cluster_filter()
        lr_ck = _state.get("_lr_params", {}).get("clustering_key", "")
        labels = _get_labels_for(lr_ck)
        import matplotlib.pyplot as _plt
        _apply_plot_font_size()
        try:
            fig = make_ligrec_plot(
                result, pvalue_threshold=pval_thresh,
                source_groups=groups, target_groups=groups,
                cluster_labels=labels,
            )
            _state["ligrec_fig"] = fig
            _plt.show(block=False)
            lr_save_plot_button.enabled = True
            if groups:
                lr_status.value = f"L-R plot displayed (clusters: {', '.join(groups)})"
            else:
                lr_status.value = "L-R plot displayed"

            # Record reproducible code
            _lr_ck = _state.get("_lr_params", {}).get("clustering_key", "")
            _record_code(
                f"\n# L-R dotplot (pvalue_threshold={pval_thresh})\n"
                f"sq.pl.ligrec(adata, cluster_key=\"{_lr_ck}\", "
                f"pvalue_threshold={pval_thresh}"
                + (f", source_groups={groups}, target_groups={groups}" if groups else "")
                + ")\nplt.show()"
            )
        except Exception as e:
            lr_status.value = f"Plot error: {e}"

    def on_save_lr_plot():
        fig = _state.get("ligrec_fig")
        if fig is None:
            return
        path = _get_plot_save_path("Save L-R Plot", "ligrec_plot")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        lr_status.value = f"L-R plot saved to {path}"

    def on_export_lr_means():
        result = _state.get("ligrec_result")
        if result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export L-R Means", "ligrec_means.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        result['means'].to_csv(path)
        lr_status.value = f"Means exported to {path}"

    def on_export_lr_pvals():
        result = _state.get("ligrec_result")
        if result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export L-R P-values", "ligrec_pvalues.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        result['pvalues'].to_csv(path)
        lr_status.value = f"P-values exported to {path}"

    lr_run_button.clicked.connect(on_run_ligrec)
    lr_plot_button.clicked.connect(on_show_lr_plot)
    lr_save_plot_button.clicked.connect(on_save_lr_plot)
    lr_export_means_button.clicked.connect(on_export_lr_means)
    lr_export_pvals_button.clicked.connect(on_export_lr_pvals)

    # ── Neighborhood Enrichment widgets (Tab 8) ──────────────────────────────
    ne_clustering_widget = ComboBox(
        label="Clustering", choices=clustering_names,
        value=clustering_names[0] if clustering_names else None,
    )
    ne_perms_slider = Slider(label="Permutations", min=100, max=1000, value=1000)
    ne_neighs_slider = Slider(label="N neighbors", min=3, max=20, value=6)
    ne_run_button = PushButton(label="Run Nhood Enrichment", enabled=True)
    ne_mode_widget = ComboBox(
        label="Display mode", choices=["zscore", "count"], value="zscore",
    )
    ne_results_text = QTextEdit()
    ne_results_text.setReadOnly(True)
    ne_results_text.setFontFamily("monospace")
    ne_results_text.setMaximumHeight(250)
    ne_plot_button = PushButton(label="Show Heatmap", enabled=False)
    ne_save_plot_button = PushButton(label="Save Plot as PNG...", enabled=False)
    ne_export_button = PushButton(label="Export Z-scores CSV...", enabled=False)
    ne_status = _StatusProxy()

    def on_run_nhood():
        ne_status.value = "Running neighborhood enrichment... (this may take a minute)"
        ne_run_button.enabled = False

        clustering_key = ne_clustering_widget.value
        n_perms = ne_perms_slider.value
        n_neighs = ne_neighs_slider.value
        _state["_ne_params"] = {
            "n_perms": n_perms,
            "n_neighs": n_neighs,
        }
        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_nhood_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, clusterings[clustering_key], clustering_key)
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            compute_spatial_neighbors(adata_norm, n_neighs=n_neighs)
            result = run_nhood_enrichment(adata_norm, clustering_key, n_perms=n_perms)
            result['_adata_norm'] = adata_norm
            result['_cluster_key'] = clustering_key
            return result
        _run()

    def _on_nhood_ready(result):
        _state["nhood_result"] = result
        ne_run_button.enabled = True

        warning = result.get('warning')
        if warning:
            ne_status.value = f"Nhood: {warning}"
            ne_results_text.setPlainText(warning)
            ne_plot_button.enabled = False
            ne_save_plot_button.enabled = False
            ne_export_button.enabled = False
            return

        zscore = result['zscore']
        clusters = result['clusters']
        n = len(clusters)

        # Summary: top enriched and depleted pairs
        lines = [
            f"Neighborhood enrichment: {n}x{n} matrix ({n} clusters)",
            f"Clusters: {', '.join(clusters)}",
            "",
        ]

        if zscore.size > 0:
            # Find top enriched pairs (off-diagonal)
            pairs = []
            for i in range(n):
                for j in range(n):
                    if i != j:
                        pairs.append((clusters[i], clusters[j], zscore[i, j]))
            pairs.sort(key=lambda x: x[2], reverse=True)

            lines.append("Top 10 enriched pairs (z-score):")
            for c1, c2, z in pairs[:10]:
                lines.append(f"  {c1} <-> {c2}: {z:.2f}")

            lines.append("")
            lines.append("Top 10 depleted pairs (z-score):")
            for c1, c2, z in pairs[-10:]:
                lines.append(f"  {c1} <-> {c2}: {z:.2f}")

        ne_results_text.setPlainText("\n".join(lines))
        ne_status.value = f"Nhood enrichment done: {n} clusters"
        ne_plot_button.enabled = True
        ne_export_button.enabled = True

        # Record reproducible code
        _ne_ck = result.get('_cluster_key', '')
        _ne_p = _state.get("_ne_params", {})
        _ne_np = _ne_p.get("n_perms", 1000)
        _ne_nn = _ne_p.get("n_neighs", 6)
        _record_clustering(_ne_ck)
        _record_spatial_neighbors(_ne_nn)
        _record_code(
            f"\n# Neighborhood enrichment (n_perms={_ne_np})\n"
            f"sq.gr.nhood_enrichment(adata, cluster_key=\"{_ne_ck}\", "
            f"n_perms={_ne_np}, seed=42)"
        )

    def on_show_nhood_plot():
        result = _state.get("nhood_result")
        if result is None:
            return
        groups = _get_cluster_filter()
        ne_ck = result.get('_cluster_key', ne_clustering_widget.value)
        labels = _get_labels_for(ne_ck)
        import matplotlib.pyplot as _plt
        _apply_plot_font_size()
        try:
            if groups or labels:
                # Custom plot with cluster filter and/or labels
                fig = make_nhood_enrichment_plot(
                    result, mode=ne_mode_widget.value,
                    cluster_filter=groups, cluster_labels=labels,
                )
            else:
                # Try squidpy native when no filter/labels active
                adata_norm = result.get('_adata_norm')
                cluster_key = result.get('_cluster_key')
                if adata_norm is not None and cluster_key is not None:
                    import squidpy as _sq
                    _sq.pl.nhood_enrichment(
                        adata_norm, cluster_key=cluster_key,
                        mode=ne_mode_widget.value,
                    )
                    fig = _plt.gcf()
                else:
                    # Session restore fallback
                    fig = make_nhood_enrichment_plot(
                        result, mode=ne_mode_widget.value,
                    )
            _state["nhood_fig"] = fig
            _plt.show(block=False)
            ne_save_plot_button.enabled = True
            if groups:
                ne_status.value = f"Heatmap displayed (clusters: {', '.join(groups)})"
            else:
                ne_status.value = "Heatmap displayed"

            # Record reproducible code
            _ne_mode = ne_mode_widget.value
            _ne_ck = result.get('_cluster_key', '')
            _record_code(
                f"\n# Nhood enrichment heatmap (mode={_ne_mode})\n"
                f"sq.pl.nhood_enrichment(adata, cluster_key=\"{_ne_ck}\", "
                f"mode=\"{_ne_mode}\")\nplt.show()"
            )
        except Exception as e:
            ne_status.value = f"Plot error: {e}"

    def on_save_nhood_plot():
        fig = _state.get("nhood_fig")
        if fig is None:
            return
        path = _get_plot_save_path("Save Nhood Enrichment Plot", "nhood_enrichment")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        ne_status.value = f"Plot saved to {path}"

    def on_export_nhood():
        result = _state.get("nhood_result")
        if result is None:
            return
        import pandas as _pd
        zscore = result['zscore']
        clusters = result['clusters']
        df = _pd.DataFrame(zscore, index=clusters, columns=clusters)
        # Apply cluster filter if active
        groups = _get_cluster_filter()
        if groups:
            keep = [c for c in groups if c in df.index]
            if keep:
                df = df.loc[keep, keep]
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Z-scores", "nhood_zscore.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path)
        ne_status.value = f"Z-scores exported to {path}"

    ne_run_button.clicked.connect(on_run_nhood)
    ne_plot_button.clicked.connect(on_show_nhood_plot)
    ne_save_plot_button.clicked.connect(on_save_nhood_plot)
    ne_export_button.clicked.connect(on_export_nhood)

    # ── Co-occurrence widgets (Tab 9) ────────────────────────────────────────
    co_clustering_widget = ComboBox(
        label="Clustering", choices=clustering_names,
        value=clustering_names[0] if clustering_names else None,
    )
    co_interval_slider = Slider(label="Distance bins", min=10, max=100, value=50)
    co_run_button = PushButton(label="Run Co-occurrence", enabled=True)
    co_results_text = QTextEdit()
    co_results_text.setReadOnly(True)
    co_results_text.setFontFamily("monospace")
    co_results_text.setMaximumHeight(250)
    co_plot_button = PushButton(label="Show Co-occurrence Plot", enabled=False)
    co_save_plot_button = PushButton(label="Save Plot as PNG...", enabled=False)
    co_export_button = PushButton(label="Export CSV...", enabled=False)
    co_filter_targets = CheckBox(label="Filter targets", value=False, enabled=True)
    co_status = _StatusProxy()

    def on_run_co_occurrence():
        co_status.value = "Running co-occurrence analysis... (this may take a minute)"
        co_run_button.enabled = False

        clustering_key = co_clustering_widget.value
        interval = co_interval_slider.value
        _state["_co_params"] = {"interval": interval}
        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_co_occurrence_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, clusterings[clustering_key], clustering_key)
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            result = run_co_occurrence(adata_norm, clustering_key, interval=interval)
            result['_adata_norm'] = adata_norm
            result['_cluster_key'] = clustering_key
            return result
        _run()

    def _on_co_occurrence_ready(result):
        _state["co_result"] = result
        co_run_button.enabled = True

        warning = result.get('warning')
        if warning:
            co_status.value = f"Co-occurrence: {warning}"
            co_results_text.setPlainText(warning)
            co_plot_button.enabled = False
            co_save_plot_button.enabled = False
            co_export_button.enabled = False
            return

        occ = result['occ']
        interval_arr = result['interval']
        clusters = result['clusters']
        n = len(clusters)

        lines = [
            f"Co-occurrence analysis: {n} clusters",
            f"Clusters: {', '.join(clusters)}",
            f"Distance range: {interval_arr[0]:.1f} – {interval_arr[-1]:.1f}",
            f"Number of distance bins: {len(interval_arr) - 1}",
            "",
            "Use 'Show Co-occurrence Plot' to visualize.",
            "Filter clusters via Cell Coloring tab to plot a subset.",
        ]

        co_results_text.setPlainText("\n".join(lines))
        co_status.value = f"Co-occurrence done: {n} clusters"
        co_plot_button.enabled = True
        co_export_button.enabled = True

        # Record reproducible code
        _co_ck = result.get('_cluster_key', '')
        _co_iv = _state.get("_co_params", {}).get("interval", 50)
        _record_clustering(_co_ck)
        _record_code(
            f"\n# Co-occurrence (interval={_co_iv})\n"
            f"sq.gr.co_occurrence(adata, cluster_key=\"{_co_ck}\", "
            f"interval={_co_iv})"
        )

    def on_show_co_plot():
        result = _state.get("co_result")
        if result is None:
            return
        groups = _get_cluster_filter()
        filter_targets = co_filter_targets.value and groups
        cc = _state.get("cluster_to_color")
        co_ck = result.get('_cluster_key', co_clustering_widget.value)
        labels = _get_labels_for(co_ck)
        import matplotlib.pyplot as _plt
        _apply_plot_font_size()
        try:
            if filter_targets or cc is not None or labels:
                # Use custom plot when we have cluster colors, filtered targets, or labels
                fig = make_co_occurrence_plot(
                    result,
                    clusters_to_plot=groups if filter_targets else groups,
                    target_clusters=groups if filter_targets else None,
                    cluster_colors=cc,
                    cluster_labels=labels,
                )
            else:
                # Try squidpy native (query filtered, all targets)
                adata_norm = result.get('_adata_norm')
                cluster_key = result.get('_cluster_key')
                if adata_norm is not None and cluster_key is not None:
                    import squidpy as _sq
                    _sq.pl.co_occurrence(
                        adata_norm, cluster_key=cluster_key, clusters=groups,
                    )
                    fig = _plt.gcf()
                else:
                    # Session restore fallback
                    fig = make_co_occurrence_plot(result, clusters_to_plot=groups, cluster_colors=cc)
            _state["co_fig"] = fig
            _plt.show(block=False)
            co_save_plot_button.enabled = True
            if groups:
                co_status.value = f"Co-occurrence plot (clusters: {', '.join(groups)})"
            else:
                co_status.value = "Co-occurrence plot displayed"

            # Record reproducible code
            _co_ck = result.get('_cluster_key', '')
            _record_code(
                f"\n# Co-occurrence plot\n"
                f"sq.pl.co_occurrence(adata, cluster_key=\"{_co_ck}\""
                + (f", clusters={groups}" if groups else "")
                + ")\nplt.show()"
            )
        except Exception as e:
            co_status.value = f"Plot error: {e}"

    def on_save_co_plot():
        fig = _state.get("co_fig")
        if fig is None:
            return
        path = _get_plot_save_path("Save Co-occurrence Plot", "co_occurrence")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        co_status.value = f"Plot saved to {path}"

    def on_export_co():
        result = _state.get("co_result")
        if result is None:
            return
        import pandas as _pd

        occ = result['occ']
        interval_arr = result['interval']
        clusters = result['clusters']
        distances = interval_arr[1:]

        groups = _get_cluster_filter()

        rows = []
        for i, src in enumerate(clusters):
            if groups and src not in groups:
                continue
            for j, tgt in enumerate(clusters):
                for k, d in enumerate(distances):
                    rows.append({
                        'source_cluster': src,
                        'target_cluster': tgt,
                        'distance': d,
                        'co_occurrence': occ[i, j, k],
                    })
        df = _pd.DataFrame(rows)
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Co-occurrence", "co_occurrence.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        co_status.value = f"Co-occurrence exported to {path}"

    co_run_button.clicked.connect(on_run_co_occurrence)
    co_plot_button.clicked.connect(on_show_co_plot)
    co_save_plot_button.clicked.connect(on_save_co_plot)
    co_export_button.clicked.connect(on_export_co)

    # ── White background callback ────────────────────────────────────────────
    def on_bg_change(value):
        viewer.window._qt_viewer.canvas.bgcolor = (1, 1, 1, 1) if value else (0, 0, 0, 1)

    bg_white_check.changed.connect(on_bg_change)

    # ── Wire events ──────────────────────────────────────────────────────────
    mode_widget.changed.connect(on_mode_change)
    filter_check.changed.connect(on_filter_change)
    clustering_widget.changed.connect(on_clustering_change)
    apply_color_button.clicked.connect(on_apply_color)
    add_gene_button.clicked.connect(on_add_gene)
    remove_gene_button.clicked.connect(on_remove_gene)
    clear_genes_button.clicked.connect(on_clear_genes)
    apply_transcripts_button.clicked.connect(on_apply_transcripts)

    # ── Helper: build a QWidget from a list of magicgui widgets ─────────────
    def _make_tab(*widgets_and_natives):
        """Pack magicgui widgets and raw QWidgets into a single QWidget."""
        tab = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        for w in widgets_and_natives:
            if hasattr(w, "native"):
                layout.addWidget(w.native)
            else:
                layout.addWidget(w)  # already a QWidget
        layout.addStretch()
        tab.setLayout(layout)
        return tab

    # ── Transcript gene list buttons row ─────────────────────────────────────
    btn_row = QWidget()
    btn_layout = QHBoxLayout()
    btn_layout.setContentsMargins(0, 0, 0, 0)
    btn_layout.addWidget(add_gene_button.native)
    btn_layout.addWidget(remove_gene_button.native)
    btn_layout.addWidget(clear_genes_button.native)
    btn_row.setLayout(btn_layout)

    # ── Assemble tabbed control panel ─────────────────────────────────────────
    tab_widget = QTabWidget()

    # Import/Export buttons row
    leiden_io_row = QWidget()
    leiden_io_layout = QHBoxLayout()
    leiden_io_layout.setContentsMargins(0, 0, 0, 0)
    leiden_io_layout.addWidget(leiden_import_button.native)
    leiden_io_layout.addWidget(leiden_export_button.native)
    leiden_io_row.setLayout(leiden_io_layout)

    # Tab 0: Clustering (Leiden)
    tab_widget.addTab(
        _make_tab(
            leiden_n_neighbors,
            leiden_n_pcs,
            leiden_resolution,
            leiden_hvg_check,
            leiden_n_hvgs,
            leiden_scale_check,
            leiden_run_button,
            leiden_status_text,
            leiden_io_row,
            leiden_edit_labels_button,
        ),
        "Clustering",
    )

    # Select All / Deselect All buttons row
    cluster_btn_row = QWidget()
    cluster_btn_layout = QHBoxLayout()
    cluster_btn_layout.setContentsMargins(0, 0, 0, 0)
    cluster_btn_layout.addWidget(select_all_btn.native)
    cluster_btn_layout.addWidget(deselect_all_btn.native)
    cluster_btn_row.setLayout(cluster_btn_layout)

    # Tab 1: Cell Coloring
    tab_widget.addTab(
        _make_tab(
            bg_white_check,
            mode_widget,
            gene_widget,
            colormap_widget,
            clustering_widget,
            filter_check,
            cluster_btn_row,
            cluster_scroll,
            apply_color_button,
        ),
        "Cell Coloring",
    )

    # Tab 2: Transcripts
    tab_widget.addTab(
        _make_tab(
            transcript_gene_widget,
            transcript_check,
            qv_slider,
            apply_transcripts_button,
            gene_list_qt,
            btn_row,
            legend_label_qt,
        ),
        "Transcripts",
    )

    # Tab 3: UMAP
    tab_widget.addTab(
        _make_tab(
            show_umap_button,
            umap_size_slider,
        ),
        "UMAP",
    )

    # Tab 4: ROI Analysis (with DEG section)
    roi_deg_header = QtLabel("── Differential Expression Between ROIs ──")
    roi_deg_header.setStyleSheet("font-weight: bold; margin-top: 10px;")
    tab_widget.addTab(
        _make_tab(
            roi_calc_button,
            roi_export_button,
            roi_deg_header,
            roi_deg_method_widget,
            roi_deg_filter_check,
            roi_deg_button,
            roi_deg_export_button,
        ),
        "ROI Analysis",
    )

    # ── Tab 5: H&E Registration ──────────────────────────────────────────────
    _he_state = {
        "he_layer": None,
        "he_tif": None,       # keep TiffFile alive for zarr store
        "he_filename": None,
        "he_path": None,      # full path to H&E file for session persistence
        "he_shape_yx": None,  # (H, W) of full-res H&E for flip computation
        "xenium_lm_layer": None,
        "he_lm_layer": None,
        "affine_3x3": None,
        "coarse_affine": None,  # coarse tissue-outline alignment
        "flip_v": False,
        "flip_h": False,
    }

    he_load_button = PushButton(label="Load H&E Image...", enabled=True)
    he_flip_v = CheckBox(label="Flip vertically", value=False)
    he_flip_h = CheckBox(label="Flip horizontally", value=False)
    # he_status_label assigned above as _StatusProxy
    he_opacity_slider = Slider(label="H&E opacity", min=0, max=100, value=70)
    he_opacity_slider.enabled = False

    coarse_align_button = PushButton(label="Coarse Align", enabled=False)

    add_xenium_lm_button = PushButton(label="Add Xenium Landmark", enabled=False)
    add_he_lm_button = PushButton(label="Add H&E Landmark", enabled=False)
    clear_lm_button = PushButton(label="Clear All", enabled=False)

    register_button = PushButton(label="Compute Registration", enabled=False)
    # reg_status_label assigned above as _StatusProxy

    reg_residuals_qt = QTextEdit()
    reg_residuals_qt.setReadOnly(True)
    reg_residuals_qt.setFontFamily("monospace")
    reg_residuals_qt.setMaximumHeight(150)

    save_lm_button = PushButton(label="Save Landmarks...", enabled=False)
    load_lm_button = PushButton(label="Load Landmarks...", enabled=True)
    save_affine_button = PushButton(label="Save Affine...", enabled=False)

    def _build_flip_affine():
        """Build a 3x3 affine that flips the H&E image around its center.

        Uses the full-resolution shape stored in _he_state to compute
        the center, then mirrors along the Y and/or X axis as requested
        by the flip checkboxes. Returns identity if no flips are active
        or no H&E is loaded.
        """
        shape = _he_state.get("he_shape_yx")
        if shape is None:
            return np.eye(3)
        h, w = shape
        M = np.eye(3)
        if he_flip_v.value:
            # Reflect Y around center: y' = h - 1 - y
            M = np.array([[  -1, 0, h - 1],
                          [   0, 1,     0],
                          [   0, 0,     1]], dtype=np.float64) @ M
        if he_flip_h.value:
            # Reflect X around center: x' = w - 1 - x
            M = np.array([[ 1,  0,     0],
                          [ 0, -1, w - 1],
                          [ 0,  0,     1]], dtype=np.float64) @ M
        return M

    def _apply_he_affine():
        """Compose flip * registration and apply to H&E + H&E landmark layers.

        Priority: fine (landmark) > coarse (tissue outline) > identity.
        Fine replaces coarse entirely (not composed) because the landmark-based
        similarity transform already captures scale+rotation+translation.
        """
        flip = _build_flip_affine()
        fine = _he_state["affine_3x3"]
        coarse = _he_state["coarse_affine"]

        if fine is not None:
            combined = fine @ flip
        elif coarse is not None:
            combined = coarse @ flip
        else:
            combined = flip

        if _he_state["he_layer"] is not None:
            _he_state["he_layer"].affine = combined
        if _he_state["he_lm_layer"] is not None:
            _he_state["he_lm_layer"].affine = combined

    def on_flip_changed(_value=None):
        _apply_he_affine()
        _he_state["flip_v"] = he_flip_v.value
        _he_state["flip_h"] = he_flip_h.value
        _save_he_affine_to_sdata()
        flips = []
        if he_flip_v.value:
            flips.append("V")
        if he_flip_h.value:
            flips.append("H")
        if flips:
            he_status_label.value = f"Flip applied: {'+'.join(flips)}"

    he_flip_v.changed.connect(on_flip_changed)
    he_flip_h.changed.connect(on_flip_changed)

    def _check_landmark_count(*_args):
        """Enable 'Compute Registration' when both layers have >= 3 points."""
        xen = _he_state["xenium_lm_layer"]
        he = _he_state["he_lm_layer"]
        if xen is not None and he is not None:
            n = min(len(xen.data), len(he.data))
            register_button.enabled = n >= 3
            save_lm_button.enabled = n >= 1

    def _create_landmark_layers():
        """Create the two Points layers for landmark placement."""
        if _he_state["xenium_lm_layer"] is not None:
            return  # already created
        xen_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="Xenium Landmarks",
            size=30,
            face_color="cyan",
            symbol="cross",
            border_color="cyan",
            border_width=0.1,
            border_width_is_relative=True,
            opacity=1.0,
        )
        he_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="H&E Landmarks",
            size=30,
            face_color="red",
            symbol="cross",
            border_color="red",
            border_width=0.1,
            border_width_is_relative=True,
            opacity=1.0,
        )
        xen_lm.events.data.connect(_check_landmark_count)
        he_lm.events.data.connect(_check_landmark_count)
        _he_state["xenium_lm_layer"] = xen_lm
        _he_state["he_lm_layer"] = he_lm
        add_xenium_lm_button.enabled = True
        add_he_lm_button.enabled = True
        clear_lm_button.enabled = True

    def on_load_he():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load H&E Image", default_dir,
            "Image Files (*.ome.tif *.tif *.tiff *.svs);;All Files (*)",
        )
        if not path:
            return
        he_status_label.value = "Loading H&E..."
        he_load_button.enabled = False

        @thread_worker(connect={"returned": _on_he_loaded})
        def load_task():
            return load_he_pyramid(path), path

        load_task()

    def _save_he_to_sdata(pyramid, he_filename):
        """Persist H&E image to sdata zarr cache as images/he_image."""
        if sdata is None or no_cache:
            return
        try:
            from spatialdata.models import Image2DModel

            # Take the full-res level, convert (Y, X, C) → (C, Y, X)
            base = np.asarray(pyramid[0])
            if base.ndim == 3 and base.shape[-1] in (3, 4):
                base_cyx = np.transpose(base, (2, 0, 1))
            else:
                base_cyx = base

            # Build multiscale spatialdata element with pyramid
            parsed = Image2DModel.parse(
                base_cyx.astype(np.uint8),
                dims=("c", "y", "x"),
                scale_factors=[2, 2, 2, 2],
                chunks=(3, 1024, 1024),
            )

            # Remove old he_image if present
            if "he_image" in sdata.images:
                del sdata.images["he_image"]

            sdata.images["he_image"] = parsed
            sdata.write_element("he_image", overwrite=True)

            # Store metadata in zarr attrs
            zarr_path = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_path), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            store["viewer_session"].attrs["he_filename"] = he_filename
            store["viewer_session"].attrs["he_shape_yx"] = list(base.shape[:2])

            print(f"  H&E image saved to sdata zarr cache ({base_cyx.shape})")
        except Exception as e:
            print(f"  Warning: could not save H&E to sdata: {e}")

    def _save_he_affine_to_sdata():
        """Update the affine transformation on he_image in sdata zarr cache."""
        if sdata is None or no_cache or "he_image" not in sdata.images:
            return
        try:
            from spatialdata.transformations import Affine as SdAffine, set_transformation

            # Build the combined affine (flip + registration) in (y, x) convention
            flip = _build_flip_affine()
            fine = _he_state["affine_3x3"]
            coarse = _he_state["coarse_affine"]
            if fine is not None:
                combined = fine @ flip
            elif coarse is not None:
                combined = coarse @ flip
            else:
                combined = flip

            sd_affine = SdAffine(combined, input_axes=("y", "x"), output_axes=("y", "x"))
            set_transformation(sdata.images["he_image"], sd_affine, "global")
            sdata.write_transformations("he_image")

            # Also persist flip + affine state to viewer_session attrs
            zarr_path = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_path), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            sess = store["viewer_session"]
            sess.attrs["flip_v"] = bool(_he_state.get("flip_v", False))
            sess.attrs["flip_h"] = bool(_he_state.get("flip_h", False))
            if fine is not None:
                sess.attrs["affine_3x3"] = fine.tolist()
            if coarse is not None:
                sess.attrs["coarse_affine"] = coarse.tolist()
        except Exception as e:
            print(f"  Warning: could not save H&E affine: {e}")

    def _on_he_loaded(result):
        (pyramid, tif), path = result
        # Remove old H&E layer if present
        if _he_state["he_layer"] is not None:
            try:
                viewer.layers.remove(_he_state["he_layer"])
            except ValueError:
                pass

        _he_state["he_tif"] = tif
        _he_state["he_filename"] = Path(path).name
        _he_state["he_path"] = str(path)
        # Store full-res spatial shape (Y, X) for flip affine computation
        base = pyramid[0]
        _he_state["he_shape_yx"] = (base.shape[0], base.shape[1])

        he_layer = viewer.add_image(
            pyramid,
            name=f"H&E ({Path(path).name})",
            rgb=True,
            blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        _he_state["he_layer"] = he_layer
        _he_state["affine_3x3"] = None
        _he_state["coarse_affine"] = None

        # Apply any pre-selected flip
        _apply_he_affine()

        # Create landmark layers if not yet present
        _create_landmark_layers()

        # Persist H&E image to sdata zarr cache
        _save_he_to_sdata(pyramid, Path(path).name)

        he_opacity_slider.enabled = True
        he_load_button.enabled = True
        coarse_align_button.enabled = morph_thumb is not None
        shape_str = "x".join(str(s) for s in pyramid[0].shape)
        he_status_label.value = f"H&E loaded: {Path(path).name} ({shape_str}, {len(pyramid)} levels)"

    def on_he_opacity(value):
        if _he_state["he_layer"] is not None:
            _he_state["he_layer"].opacity = value / 100.0

    def on_coarse_align():
        if _he_state["he_layer"] is None:
            reg_status_label.value = "Load H&E image first"
            return
        if morph_thumb is None:
            reg_status_label.value = "No morphology data available"
            return

        reg_status_label.value = "Computing coarse alignment..."
        coarse_align_button.enabled = False

        @thread_worker(connect={"returned": _on_coarse_done})
        def _compute_coarse():
            he_pyramid = _he_state["he_layer"].data
            he_low = np.asarray(he_pyramid[-1])  # (Y, X, 3) RGB, lowest res

            # Extract tissue masks
            morph_mask = extract_tissue_mask_fluorescence(morph_thumb)
            he_mask = extract_tissue_mask_he(he_low)

            # Compute downsample factors
            target_ds = morph_full_shape_yx[0] / morph_thumb.shape[1]  # full_Y / thumb_Y
            he_full_shape = _he_state["he_shape_yx"]
            source_ds = he_full_shape[0] / he_low.shape[0]  # full_Y / thumb_Y

            coarse_affine_yx = compute_coarse_affine(
                target_mask=morph_mask,
                source_mask=he_mask,
                target_downsample=target_ds,
                source_downsample=source_ds,
            )
            return coarse_affine_yx

        _compute_coarse()

    def _on_coarse_done(coarse_affine):
        _he_state["coarse_affine"] = coarse_affine
        _he_state["affine_3x3"] = None  # clear any previous fine registration
        _apply_he_affine()
        _save_he_affine_to_sdata()
        coarse_align_button.enabled = True
        scale = np.sqrt(coarse_affine[0, 0]**2 + coarse_affine[0, 1]**2)
        reg_status_label.value = f"Coarse aligned (scale={scale:.4f}). Place landmarks to refine."
        reg_residuals_qt.setPlainText(
            f"Coarse tissue-outline alignment applied.\n"
            f"Scale: {scale:.4f}\n"
            f"Place >= 3 matching landmarks, then click 'Compute Registration'\n"
            f"to refine alignment."
        )

    def on_add_xenium_lm():
        lm = _he_state["xenium_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.value = "Click on a feature in the Xenium image"

    def on_add_he_lm():
        lm = _he_state["he_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.value = "Click on the same feature in the H&E image"

    def on_clear_lm():
        for key in ("xenium_lm_layer", "he_lm_layer"):
            lm = _he_state[key]
            if lm is not None:
                lm.selected_data = set()
                lm.data = np.empty((0, 2), dtype=np.float64)
        _he_state["affine_3x3"] = None
        _he_state["coarse_affine"] = None
        # Reset to flip-only affine (or identity if no flips)
        _apply_he_affine()
        reg_residuals_qt.clear()
        reg_status_label.value = "Landmarks cleared"
        register_button.enabled = False
        save_lm_button.enabled = False
        save_affine_button.enabled = False

    def on_register():
        xen_pts = _he_state["xenium_lm_layer"].data
        he_pts = _he_state["he_lm_layer"].data
        n = min(len(xen_pts), len(he_pts))
        if n < 3:
            reg_status_label.value = "Need at least 3 paired landmarks"
            return

        xen_pts = np.asarray(xen_pts[:n], dtype=np.float64)
        he_pts = np.asarray(he_pts[:n], dtype=np.float64)

        affine, residuals = compute_landmark_affine(xen_pts, he_pts)
        _he_state["affine_3x3"] = affine

        # Apply registration + flip to H&E image and landmark layers
        _apply_he_affine()

        # Display residuals
        lines = [f"Registration: {n} landmarks, similarity transform"]
        lines.append(f"Mean residual: {residuals.mean():.1f} px ({residuals.mean() * pixel_size:.1f} um)")
        lines.append(f"Max  residual: {residuals.max():.1f} px ({residuals.max() * pixel_size:.1f} um)")
        lines.append("")
        for i, r in enumerate(residuals):
            lines.append(f"  Landmark {i+1}: {r:.1f} px ({r * pixel_size:.1f} um)")
        # Show scale factor from the affine
        scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
        lines.append(f"\nScale factor: {scale:.4f}")
        reg_residuals_qt.setPlainText("\n".join(lines))

        reg_status_label.value = f"Registered ({n} landmarks, mean residual {residuals.mean():.1f} px)"
        _save_he_affine_to_sdata()
        save_affine_button.enabled = True

    def on_save_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Landmarks", default_dir + "/landmarks.json",
            "JSON Files (*.json)",
        )
        if not path:
            return
        xen_pts = np.asarray(_he_state["xenium_lm_layer"].data, dtype=np.float64)
        he_pts = np.asarray(_he_state["he_lm_layer"].data, dtype=np.float64)
        save_landmarks(
            path, xen_pts, he_pts,
            affine=_he_state["affine_3x3"],
            he_filename=_he_state["he_filename"],
        )
        reg_status_label.value = f"Landmarks saved to {Path(path).name}"

    def on_load_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load Landmarks", default_dir,
            "JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        data = load_landmarks(path)
        # Ensure landmark layers exist
        _create_landmark_layers()
        _he_state["xenium_lm_layer"].data = data["xenium_landmarks_yx"]
        _he_state["he_lm_layer"].data = data["he_landmarks_yx"]
        if "affine_3x3_yx" in data:
            affine = data["affine_3x3_yx"]
            _he_state["affine_3x3"] = affine
            _apply_he_affine()
            save_affine_button.enabled = True
            # Show scale factor
            scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
            reg_residuals_qt.setPlainText(f"Loaded affine (scale={scale:.4f})")
        if "he_filename" in data:
            _he_state["he_filename"] = data["he_filename"]
        n = min(len(data["xenium_landmarks_yx"]), len(data["he_landmarks_yx"]))
        reg_status_label.value = f"Loaded {n} landmarks from {Path(path).name}"

    def on_save_affine():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Affine", default_dir + "/he_affine.json",
            "JSON Files (*.json)",
        )
        if not path:
            return
        affine = _he_state["affine_3x3"]
        if affine is None:
            return
        with open(path, "w") as f:
            json.dump({"affine_3x3_yx": affine.tolist()}, f, indent=2)
        reg_status_label.value = f"Affine saved to {Path(path).name}"

    # Wire H&E Registration events
    he_load_button.clicked.connect(on_load_he)
    he_opacity_slider.changed.connect(on_he_opacity)
    coarse_align_button.clicked.connect(on_coarse_align)
    add_xenium_lm_button.clicked.connect(on_add_xenium_lm)
    add_he_lm_button.clicked.connect(on_add_he_lm)
    clear_lm_button.clicked.connect(on_clear_lm)
    register_button.clicked.connect(on_register)
    save_lm_button.clicked.connect(on_save_landmarks)
    load_lm_button.clicked.connect(on_load_landmarks)
    save_affine_button.clicked.connect(on_save_affine)

    # Landmark buttons row
    lm_btn_row = QWidget()
    lm_btn_layout = QHBoxLayout()
    lm_btn_layout.setContentsMargins(0, 0, 0, 0)
    lm_btn_layout.addWidget(add_xenium_lm_button.native)
    lm_btn_layout.addWidget(add_he_lm_button.native)
    lm_btn_layout.addWidget(clear_lm_button.native)
    lm_btn_row.setLayout(lm_btn_layout)

    # Save/load buttons row
    io_btn_row = QWidget()
    io_btn_layout = QHBoxLayout()
    io_btn_layout.setContentsMargins(0, 0, 0, 0)
    io_btn_layout.addWidget(save_lm_button.native)
    io_btn_layout.addWidget(load_lm_button.native)
    io_btn_layout.addWidget(save_affine_button.native)
    io_btn_row.setLayout(io_btn_layout)

    # Flip checkboxes row
    flip_row = QWidget()
    flip_layout = QHBoxLayout()
    flip_layout.setContentsMargins(0, 0, 0, 0)
    flip_layout.addWidget(he_flip_v.native)
    flip_layout.addWidget(he_flip_h.native)
    flip_row.setLayout(flip_layout)

    # Tab 5: H&E Registration
    tab_widget.addTab(
        _make_tab(
            he_load_button,
            flip_row,
            he_opacity_slider,
            coarse_align_button,
            lm_btn_row,
            register_button,
            io_btn_row,
        ),
        "H&E Registration",
    )

    # Tab 6: Gene Analysis
    ga_dotplot_btn_row = QWidget()
    ga_dotplot_btn_layout = QHBoxLayout()
    ga_dotplot_btn_layout.setContentsMargins(0, 0, 0, 0)
    ga_dotplot_btn_layout.addWidget(ga_dotplot_button.native)
    ga_dotplot_btn_layout.addWidget(ga_edit_labels_button.native)
    ga_dotplot_btn_layout.addWidget(ga_save_dotplot_button.native)
    ga_dotplot_btn_row.setLayout(ga_dotplot_btn_layout)

    tab_widget.addTab(
        _make_tab(
            ga_clustering_widget,
            ga_method_widget,
            ga_n_genes_slider,
            ga_run_button,
            ga_dotplot_n_slider,
            ga_dendro_check,
            ga_dotplot_btn_row,
            ga_rank_plot_button,
            ga_export_button,
        ),
        "Gene Analysis",
    )

    # Tab 7: Ligand-Receptor
    lr_export_btn_row = QWidget()
    lr_export_btn_layout = QHBoxLayout()
    lr_export_btn_layout.setContentsMargins(0, 0, 0, 0)
    lr_export_btn_layout.addWidget(lr_export_means_button.native)
    lr_export_btn_layout.addWidget(lr_export_pvals_button.native)
    lr_export_btn_row.setLayout(lr_export_btn_layout)

    lr_plot_btn_row = QWidget()
    lr_plot_btn_layout = QHBoxLayout()
    lr_plot_btn_layout.setContentsMargins(0, 0, 0, 0)
    lr_plot_btn_layout.addWidget(lr_plot_button.native)
    lr_plot_btn_layout.addWidget(lr_save_plot_button.native)
    lr_plot_btn_row.setLayout(lr_plot_btn_layout)

    tab_widget.addTab(
        _make_tab(
            lr_clustering_widget,
            lr_perms_slider,
            lr_neighs_slider,
            lr_ds_group,
            lr_cpdb_only,
            lr_run_button,
            lr_pval_widget,
            lr_plot_btn_row,
            lr_export_btn_row,
        ),
        "Ligand-Receptor",
    )

    # Tab 8: Neighborhood Enrichment
    ne_plot_btn_row = QWidget()
    ne_plot_btn_layout = QHBoxLayout()
    ne_plot_btn_layout.setContentsMargins(0, 0, 0, 0)
    ne_plot_btn_layout.addWidget(ne_plot_button.native)
    ne_plot_btn_layout.addWidget(ne_save_plot_button.native)
    ne_plot_btn_row.setLayout(ne_plot_btn_layout)

    tab_widget.addTab(
        _make_tab(
            ne_clustering_widget,
            ne_perms_slider,
            ne_neighs_slider,
            ne_run_button,
            ne_mode_widget,
            ne_results_text,
            ne_plot_btn_row,
            ne_export_button,
        ),
        "Nhood Enrichment",
    )

    # Tab 9: Co-occurrence
    co_plot_btn_row = QWidget()
    co_plot_btn_layout = QHBoxLayout()
    co_plot_btn_layout.setContentsMargins(0, 0, 0, 0)
    co_plot_btn_layout.addWidget(co_plot_button.native)
    co_plot_btn_layout.addWidget(co_save_plot_button.native)
    co_plot_btn_row.setLayout(co_plot_btn_layout)

    tab_widget.addTab(
        _make_tab(
            co_clustering_widget,
            co_interval_slider,
            co_run_button,
            co_results_text,
            co_filter_targets,
            co_plot_btn_row,
            co_export_button,
        ),
        "Co-occurrence",
    )

    # ── Session restore function ─────────────────────────────────────────────
    def restore_session(session):
        """Apply loaded session data to the viewer state."""
        # ROIs
        rois = session.get("rois", [])
        if rois and roi_layer is not None:
            roi_layer.data = rois
            print(f"  Restored {len(rois)} ROI polygons")

        # Cluster labels (per-clustering nested dict)
        cl = session.get("cluster_labels")
        if cl and isinstance(cl, dict):
            _state["cluster_labels"] = cl
            n_clusterings = len(cl)
            n_labels = sum(len(v) for v in cl.values() if isinstance(v, dict))
            print(f"  Restored cluster labels: {n_labels} labels across {n_clusterings} clustering(s)")

        # Custom clusterings (Leiden, imported)
        cc = session.get("custom_clusterings", {})
        if cc:
            for name, series in cc.items():
                clusterings[name] = series
                _state["custom_clusterings"][name] = series
            _refresh_clustering_choices()
            print(f"  Restored {len(cc)} custom clustering(s): {', '.join(cc.keys())}")

        # Analysis results: rank genes
        rg = session.get("rank_genes_df")
        if rg is not None:
            _state["rank_genes_df"] = rg
            ga_export_button.enabled = True
            print(f"  Restored rank genes ({len(rg)} rows)")

        # Analysis results: ROI DEG
        rd = session.get("roi_deg_df")
        if rd is not None:
            _state["roi_deg_df"] = rd
            roi_deg_export_button.enabled = True
            print(f"  Restored ROI DEG ({len(rd)} rows)")

        # Analysis results: ligrec
        lm = session.get("ligrec_means")
        lp = session.get("ligrec_pvalues")
        if lm is not None or lp is not None:
            _state["ligrec_result"] = {
                "means": lm if lm is not None else __import__("pandas").DataFrame(),
                "pvalues": lp if lp is not None else __import__("pandas").DataFrame(),
                "warning": None,
            }
            lr_export_means_button.enabled = lm is not None
            lr_export_pvals_button.enabled = lp is not None
            lr_plot_button.enabled = lm is not None and not lm.empty
            print(f"  Restored L-R results")

        # Analysis results: nhood enrichment
        nh = session.get("nhood_result")
        if nh is not None:
            _state["nhood_result"] = nh
            ne_plot_button.enabled = True
            ne_export_button.enabled = True
            n = len(nh.get('clusters', []))
            print(f"  Restored nhood enrichment ({n} clusters)")

        # Analysis results: co-occurrence
        co = session.get("co_result")
        if co is not None:
            _state["co_result"] = co
            co_plot_button.enabled = True
            co_export_button.enabled = True
            n = len(co.get('clusters', []))
            print(f"  Restored co-occurrence ({n} clusters)")

        # H&E restore — load from sdata zarr cache
        if sdata is not None and "he_image" in sdata.images:
            he_flip_v.value = session.get("flip_v", False)
            he_flip_h.value = session.get("flip_h", False)

            _session_he_data = {
                "affine_3x3": session.get("affine_3x3"),
                "coarse_affine": session.get("coarse_affine"),
                "xenium_landmarks": session.get("xenium_landmarks"),
                "he_landmarks": session.get("he_landmarks"),
                "he_filename": session.get("he_filename", "H&E"),
                "he_shape_yx": session.get("he_shape_yx"),
            }

            he_status_label.value = "Restoring H&E from cache..."

            @thread_worker(connect={"returned": lambda result: _on_he_restored_from_sdata(result, _session_he_data)})
            def _load_he_from_sdata():
                he_dt = sdata.images["he_image"]
                pyramid = _extract_dt_scales(he_dt)
                # Convert from CYX to YXC (RGB) for napari
                pyramid_rgb = []
                for arr in pyramid:
                    computed = arr.compute() if hasattr(arr, 'compute') else np.asarray(arr)
                    if computed.ndim == 3 and computed.shape[0] in (3, 4):
                        computed = np.transpose(computed, (1, 2, 0))
                    pyramid_rgb.append(computed)
                return pyramid_rgb

            _load_he_from_sdata()
        elif session.get("he_filename"):
            print(f"  Warning: H&E image not found in sdata cache, skipping H&E restore")

    def _on_he_restored_from_sdata(pyramid_rgb, session_he_data):
        """Callback after H&E loads from sdata cache — apply saved affine."""
        if _he_state["he_layer"] is not None:
            try:
                viewer.layers.remove(_he_state["he_layer"])
            except ValueError:
                pass

        he_filename = session_he_data.get("he_filename", "H&E")
        _he_state["he_tif"] = None
        _he_state["he_filename"] = he_filename
        _he_state["he_path"] = None  # loaded from cache, no file path
        base = pyramid_rgb[0]
        _he_state["he_shape_yx"] = session_he_data.get("he_shape_yx") or (base.shape[0], base.shape[1])

        he_layer = viewer.add_image(
            pyramid_rgb,
            name=f"H&E ({he_filename})",
            rgb=True,
            blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        _he_state["he_layer"] = he_layer

        # Restore saved affines
        _he_state["affine_3x3"] = session_he_data.get("affine_3x3")
        _he_state["coarse_affine"] = session_he_data.get("coarse_affine")
        _apply_he_affine()

        # Create landmark layers and populate with saved data
        _create_landmark_layers()
        xen_lm = session_he_data.get("xenium_landmarks")
        he_lm = session_he_data.get("he_landmarks")
        if xen_lm is not None and _he_state["xenium_lm_layer"] is not None:
            _he_state["xenium_lm_layer"].data = xen_lm
        if he_lm is not None and _he_state["he_lm_layer"] is not None:
            _he_state["he_lm_layer"].data = he_lm

        he_opacity_slider.enabled = True
        he_load_button.enabled = True
        coarse_align_button.enabled = morph_thumb is not None
        save_affine_button.enabled = _he_state["affine_3x3"] is not None

        has_affine = _he_state["affine_3x3"] is not None or _he_state["coarse_affine"] is not None
        he_status_label.value = f"H&E restored: {he_filename}" + (" (with registration)" if has_affine else "")
        print(f"  Restored H&E from cache: {he_filename}" + (" with registration" if has_affine else ""))

    return tab_widget, _state, _he_state, restore_session


if __name__ == "__main__":
    data_path, no_cache = _parse_args()
    main(data_path, no_cache=no_cache)
