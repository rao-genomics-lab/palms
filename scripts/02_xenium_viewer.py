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
        data_path=data_path,
        morph_thumb=morph_thumb,
        morph_full_shape_yx=morph_full_shape_yx,
        adata=adata,
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
    data_path: Path = None,
    morph_thumb=None,
    morph_full_shape_yx=None,
    adata=None,
):
    """Build and return a magicgui Container docked widget."""
    from magicgui.widgets import (
        Container, ComboBox, CheckBox, PushButton, Label,
        Slider, RadioButtons,
    )
    from qtpy.QtWidgets import (
        QListWidget, QHBoxLayout, QWidget, QVBoxLayout, QLabel, QTabWidget,
        QCheckBox, QScrollArea, QGridLayout,
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
        ids = sorted(clusterings[key].dropna().unique().astype(int).tolist())
        cols = 3
        for i, cid in enumerate(ids):
            cb = QCheckBox(str(cid))
            cb.setChecked(True)
            cb.setEnabled(filter_check.value)
            cluster_filter_grid.addWidget(cb, i // cols, i % cols)
            _state["cluster_checkboxes"][cid] = cb

    def _get_selected_cluster_ids():
        """Return set of cluster IDs whose checkboxes are checked."""
        return {cid for cid, cb in _state["cluster_checkboxes"].items() if cb.isChecked()}

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
            color_arr = color_arr.copy()
            mask_out = ~np.isin(label_to_cluster_arr, list(selected_ids))
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
        clustering_key, (color_arr, cluster_to_color), selected_ids = result
        # Get per-obs cluster IDs for both UMAP hover and spatial hover
        cluster_ids_per_obs, label_to_cluster = _get_cluster_ids_per_obs(clustering_key)

        # If filtering by selected clusters, zero out all other cells
        if selected_ids is not None:
            color_arr = color_arr.copy()
            mask_out = ~np.isin(label_to_cluster, list(selected_ids))
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
                clusters_aligned = cluster_series.reindex(cell_ids_arr, fill_value=-1)
            else:
                clusters_aligned = cluster_series.reindex(adata.obs_names, fill_value=-1)
            cluster_mask = np.isin(clusters_aligned.values.astype(np.int32), list(selected_ids))
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
                clusters_aligned = cluster_series.reindex(cell_ids_arr, fill_value=-1)
            else:
                clusters_aligned = cluster_series.reindex(_adata.obs_names, fill_value=-1)
            cluster_mask = np.isin(clusters_aligned.values.astype(np.int32), list(selected_ids))

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
        labels = _state.get("cluster_labels")

        @thread_worker(connect={"returned": _on_dotplot_ready})
        def _run():
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

    def _open_label_editor():
        from qtpy.QtWidgets import QDialog, QGridLayout, QLineEdit, QDialogButtonBox

        clustering_key = ga_clustering_widget.value
        cluster_series = clusterings[clustering_key]
        ids = sorted(cluster_series.dropna().unique().astype(int).tolist())
        existing = _state.get("cluster_labels", {})

        dialog = QDialog()
        dialog.setWindowTitle("Edit Cluster Labels")
        grid = QGridLayout()
        edits = {}
        for i, cid in enumerate(ids):
            grid.addWidget(QtLabel(f"Cluster {cid}:"), i, 0)
            edit = QLineEdit(existing.get(cid, str(cid)))
            grid.addWidget(edit, i, 1)
            edits[cid] = edit
        from qtpy.QtWidgets import QDialogButtonBox as QDBBox
        buttons = QDBBox(QDBBox.Ok | QDBBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        grid.addWidget(buttons, len(ids), 0, 1, 2)
        dialog.setLayout(grid)

        if dialog.exec_() == QDialog.Accepted:
            _state["cluster_labels"] = {cid: e.text() for cid, e in edits.items()}
            ga_status.value = f"Labels updated for {len(ids)} clusters"

    def on_save_dotplot():
        fig = _state.get("dotplot_fig")
        if fig is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Dotplot", "dotplot.png", "PNG Files (*.png);;All Files (*)",
        )
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
        fig = make_rank_genes_plot(adata_norm, n_genes=ga_n_genes_slider.value)
        _plt.show(block=False)
        ga_status.value = "Rank genes plot displayed"

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
        _adata = adata if adata is not None else color_manager.adata

        @thread_worker(connect={"returned": _on_ligrec_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, clusterings[clustering_key], clustering_key)
            # Set spatial coordinates for squidpy
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            compute_spatial_neighbors(adata_norm, n_neighs=n_neighs)
            result = run_ligrec(adata_norm, clustering_key, n_perms=n_perms)
            return result
        _run()

    def _on_ligrec_ready(result):
        _state["ligrec_result"] = result
        lr_run_button.enabled = True

        warning = result.get('warning')
        means = result['means']
        pvalues = result['pvalues']

        if warning:
            lr_results_text.setPlainText(warning)
            lr_status.value = "L-R analysis: warning (see results)"
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

    def on_show_lr_plot():
        result = _state.get("ligrec_result")
        if result is None:
            return
        pval_thresh = float(lr_pval_widget.value)
        import matplotlib.pyplot as _plt
        try:
            fig = make_ligrec_plot(result, pvalue_threshold=pval_thresh)
            _state["ligrec_fig"] = fig
            _plt.show(block=False)
            lr_save_plot_button.enabled = True
            lr_status.value = "L-R plot displayed"
        except Exception as e:
            lr_status.value = f"Plot error: {e}"

    def on_save_lr_plot():
        fig = _state.get("ligrec_fig")
        if fig is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Save L-R Plot", "ligrec_plot.png", "PNG Files (*.png);;All Files (*)",
        )
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
        "he_shape_yx": None,  # (H, W) of full-res H&E for flip computation
        "xenium_lm_layer": None,
        "he_lm_layer": None,
        "affine_3x3": None,
        "coarse_affine": None,  # coarse tissue-outline alignment
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
            lr_run_button,
            lr_pval_widget,
            lr_plot_btn_row,
            lr_export_btn_row,
        ),
        "Ligand-Receptor",
    )

    return tab_widget


if __name__ == "__main__":
    data_path, no_cache = _parse_args()
    main(data_path, no_cache=no_cache)
