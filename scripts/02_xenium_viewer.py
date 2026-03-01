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
    panel = _build_control_panel(
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
    )
    viewer.window.add_dock_widget(panel, name="Xenium Controls", area="right")

    total_time = time.perf_counter() - t_start
    print(f"\nViewer ready in {total_time:.1f}s. Close the napari window to exit.")
    napari.run()


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
            viewer.add_image(
                scales,
                name="morphology_focus",
                channel_axis=0,
                colormap=CHANNEL_COLORMAPS,
                contrast_limits=CHANNEL_CONTRAST,
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
):
    """Build and return a magicgui Container docked widget."""
    from magicgui.widgets import (
        Container, ComboBox, CheckBox, PushButton, Label,
        Slider, RadioButtons,
    )
    from qtpy.QtWidgets import (
        QListWidget, QHBoxLayout, QWidget, QVBoxLayout, QLabel, QTabWidget,
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
    }

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

    # Get initial cluster IDs for the first clustering
    _initial_cluster_ids = []
    if clustering_names:
        _initial_cluster_ids = sorted(
            clusterings[clustering_names[0]].dropna().unique().astype(int).tolist()
        )
    cluster_id_widget = ComboBox(
        label="Cluster ID",
        choices=[str(c) for c in _initial_cluster_ids],
        value=str(_initial_cluster_ids[0]) if _initial_cluster_ids else None,
        enabled=False,
    )

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

    # ── Status ────────────────────────────────────────────────────────────────
    status_label = Label(value="Ready")

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

    # ── Helper: repopulate cluster ID choices ─────────────────────────────────
    def _repopulate_cluster_ids():
        key = clustering_widget.value
        if key and key in clusterings:
            ids = sorted(clusterings[key].dropna().unique().astype(int).tolist())
            cluster_id_widget.choices = [str(c) for c in ids]
            if ids:
                cluster_id_widget.value = str(ids[0])

    # ── Cell Coloring callbacks ─────────────────────────────────────────────
    def on_mode_change(value):
        _state["color_mode"] = value
        is_gene = (value == "Gene Expression")
        gene_widget.enabled = is_gene
        colormap_widget.enabled = is_gene
        filter_check.enabled = is_gene
        clustering_widget.enabled = (value == "Cluster") or (is_gene and filter_check.value)
        cluster_id_widget.enabled = is_gene and filter_check.value

    def on_filter_change(value):
        _state["filter_by_cluster"] = value
        is_gene = (_state["color_mode"] == "Gene Expression")
        clustering_widget.enabled = (not is_gene) or value
        cluster_id_widget.enabled = is_gene and value

    def on_clustering_change(value):
        _repopulate_cluster_ids()

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
            if use_filter:
                clustering_key = clustering_widget.value
                cluster_id = int(cluster_id_widget.value)
                cluster_series = clusterings[clustering_key]
                cluster_series.name = clustering_key

                @thread_worker(connect={"returned": _on_gene_colors_ready})
                def compute_gene_filtered():
                    return gene, color_manager.get_gene_colors_filtered(
                        gene, cluster_series, cluster_id, colormap=cmap
                    )
                compute_gene_filtered()
            else:
                @thread_worker(connect={"returned": _on_gene_colors_ready})
                def compute_gene():
                    return gene, color_manager.get_gene_colors(gene, colormap=cmap)
                compute_gene()

        else:  # Cluster
            clustering_key = clustering_widget.value
            _state["current_clustering"] = clustering_key
            cluster_series = clusterings[clustering_key]
            cluster_series.name = clustering_key

            @thread_worker(connect={"returned": _on_cluster_colors_ready})
            def compute_cluster():
                return clustering_key, color_manager.get_cluster_colors(cluster_series)

            compute_cluster()

    def _on_gene_colors_ready(result):
        gene, color_arr = result
        color_manager.apply_to_labels_layer(cell_labels_layer, color_arr)
        umap_viewer.color_by_gene(gene, color_arr, label_to_obs)
        # Clear cluster hover lookup (no longer showing clusters)
        _state["label_to_cluster"] = None
        _state["active_clustering_name"] = None
        status_label.value = f"Cells colored by gene: {gene}"
        apply_color_button.enabled = True

    def _get_cluster_ids_per_obs(clustering_key):
        """Return per-obs cluster IDs aligned to adata, plus label lookup."""
        cluster_series = clusterings[clustering_key]
        adata = color_manager.adata
        if 'cell_id' in adata.obs.columns:
            cell_ids = adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids, fill_value=-1)
        else:
            clusters_aligned = cluster_series.reindex(adata.obs_names, fill_value=-1)
        cluster_values = clusters_aligned.values.astype(np.int32)

        # Also build label -> cluster lookup for spatial hover
        max_label = len(label_to_obs) - 1
        label_to_cluster = np.full(max_label + 1, -1, dtype=np.int32)
        valid_mask = label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = label_to_obs[valid_labels]
        label_to_cluster[valid_labels] = cluster_values[obs_indices]
        return cluster_values, label_to_cluster

    def _on_cluster_colors_ready(result):
        clustering_key, (color_arr, cluster_to_color) = result
        color_manager.apply_to_labels_layer(cell_labels_layer, color_arr)
        # Get per-obs cluster IDs for both UMAP hover and spatial hover
        cluster_ids_per_obs, label_to_cluster = _get_cluster_ids_per_obs(clustering_key)
        umap_viewer.color_by_cluster(
            clustering_key, color_arr, label_to_obs,
            cluster_ids_per_obs=cluster_ids_per_obs,
        )
        _state["label_to_cluster"] = label_to_cluster
        _state["active_clustering_name"] = clustering_key
        status_label.value = f"Cells colored by cluster: {clustering_key}"
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
                    text = f"Cell {int(label_val)} — {name}: {cid}"
                else:
                    text = f"Cell {int(label_val)} — {name}: unassigned"
                QTimer.singleShot(0, lambda t=text: setattr(viewer, 'status', t))

        viewer.cursor.events.position.connect(_on_cursor_move)

    # ── ROI Analysis widgets ──────────────────────────────────────────────────
    from qtpy.QtWidgets import QTextEdit, QFileDialog
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

        lines = [f"Gene: {gene}", ""]
        roi_results = []  # list of dicts for export

        for i, poly_yx in enumerate(polygons):
            # poly_yx is Nx2 in napari (y, x) coords; shapely needs (x, y)
            poly_xy = poly_yx[:, ::-1]
            shapely_poly = ShapelyPolygon(poly_xy)
            if not shapely_poly.is_valid:
                shapely_poly = shapely_poly.buffer(0)

            # Check which centroids are inside
            inside = contains_xy(shapely_poly, centroids_yx[:, 1], centroids_yx[:, 0])
            inside_idx = np.where(inside)[0]
            n_cells = len(inside_idx)

            if n_cells == 0:
                lines.append(f"Region {i+1}: 0 cells")
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

    # Tab 1: Cell Coloring
    tab_widget.addTab(
        _make_tab(
            mode_widget,
            gene_widget,
            colormap_widget,
            clustering_widget,
            filter_check,
            cluster_id_widget,
            apply_color_button,
            status_label,
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

    # Tab 4: ROI Analysis
    tab_widget.addTab(
        _make_tab(
            roi_calc_button,
            roi_text,
            roi_export_button,
        ),
        "ROI Analysis",
    )

    return tab_widget


if __name__ == "__main__":
    data_path, no_cache = _parse_args()
    main(data_path, no_cache=no_cache)
