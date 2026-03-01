"""
Xenium Viewer — Linux equivalent of Xenium Explorer.

Usage:
    conda activate xenium_viewer
    python scripts/02_xenium_viewer.py

Opens a napari window with:
  - 4-channel morphology_focus image (DAPI, markers, 18S, SMA/Vim)
  - Cell and nucleus labels (raster masks) for fast coloring
  - Transcript points layer (populated on demand per gene)
  - A docked control panel (gene/cluster selection, colormap, QV threshold)
  - A linked matplotlib UMAP window

Performance notes:
  - Morphology TIFFs are rendered via a 5-level software pyramid (no internal pyramid)
  - Cell boundaries (318K shapes) are kept hidden; use contour=2 on labels layer
  - Transcript loading uses feather cache (run 00_preprocess_transcripts.py first)
  - Label coloring uses DirectLabelColormap for O(nonzero) construction
"""

import os
import sys
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
DATA_DIR = SCRIPTS_DIR.parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.coloring import CellColorManager, AVAILABLE_COLORMAPS, CLUSTER_PALETTE
from utils.transcript_index import TranscriptLoader
from utils.umap_widget import UMAPWidget

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



def _import_loader():
    """Import 01_load_sdata as a module regardless of the numeric prefix."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "load_sdata", SCRIPTS_DIR / "01_load_sdata.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    print("=" * 60)
    print("Xenium Linux Viewer")
    print("=" * 60)

    # ── Load data ────────────────────────────────────────────────────────────
    loader_mod = _import_loader()
    print("Loading SpatialData...")
    sdata = loader_mod.load_sdata()

    print("Loading UMAP...")
    umap_df = loader_mod.load_umap()

    print("Loading cluster assignments...")
    clusterings = loader_mod.load_clusterings()

    print("Building label→obs mapping...")
    label_to_obs = loader_mod.get_label_to_obs_mapping(sdata)

    adata = sdata["table"]
    gene_names = list(adata.var_names)
    clustering_names = list(clusterings.keys())

    print(f"Genes: {len(gene_names)}, Clusterings: {len(clustering_names)}")

    # ── Managers ─────────────────────────────────────────────────────────────
    color_manager = CellColorManager(adata, label_to_obs)
    transcript_loader = TranscriptLoader(cache_dir=SCRIPTS_DIR / "transcript_cache")

    # ── Napari viewer (must be created before any QWidgets) ──────────────────
    print("Opening napari...")
    viewer = napari.Viewer(title="Xenium Linux Viewer")

    # ── UMAP widget (needs QApplication from napari.Viewer) ──────────────────
    umap_widget = UMAPWidget(umap_df, adata.obs_names)

    # ── Add layers from sdata ─────────────────────────────────────────────────
    _add_layers_manually(viewer, sdata)

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
        umap_widget=umap_widget,
        label_to_obs=label_to_obs,
    )
    viewer.window.add_dock_widget(panel, name="Xenium Controls", area="right")

    # ── UMAP dock widget ───────────────────────────────────────────────────────
    viewer.window.add_dock_widget(umap_widget, name="UMAP", area="bottom")

    print("\nViewer ready. Close the napari window to exit.")
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
    umap_widget: UMAPWidget,
    label_to_obs: np.ndarray,
):
    """Build and return a magicgui Container docked widget."""
    from magicgui import magicgui
    from magicgui.widgets import (
        Container, ComboBox, CheckBox, PushButton, Label,
        Slider, RadioButtons, SpinBox,
    )

    # ── State ────────────────────────────────────────────────────────────────
    _state = {
        "current_gene": gene_names[0] if gene_names else None,
        "current_clustering": clustering_names[0] if clustering_names else None,
        "current_colormap": "viridis",
        "show_transcripts": False,
        "min_qv": 20,
        "color_mode": "Gene Expression",
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

    apply_color_button = PushButton(label="Apply Cell Coloring", enabled=True)

    # ── Transcript Overlay widgets ────────────────────────────────────────────
    transcript_gene_widget = ComboBox(
        label="Transcript gene",
        choices=gene_names,
        value=gene_names[0] if gene_names else None,
    )

    transcript_check = CheckBox(label="Show transcripts", value=False)

    qv_slider = Slider(label="Min QV", min=0, max=40, value=20)

    apply_transcripts_button = PushButton(label="Apply Transcripts", enabled=True)

    # ── Status ────────────────────────────────────────────────────────────────
    status_label = Label(value="Ready")

    # ── Cell Coloring callbacks ─────────────────────────────────────────────
    def on_mode_change(value):
        _state["color_mode"] = value
        gene_widget.enabled = (value == "Gene Expression")
        colormap_widget.enabled = (value == "Gene Expression")
        clustering_widget.enabled = (value == "Cluster")

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
        umap_widget.color_by_gene(gene, color_arr, label_to_obs)
        status_label.value = f"Cells colored by gene: {gene}"
        apply_color_button.enabled = True

    def _on_cluster_colors_ready(result):
        clustering_key, (color_arr, cluster_to_color) = result
        color_manager.apply_to_labels_layer(cell_labels_layer, color_arr)
        umap_widget.color_by_cluster(clustering_key, color_arr, label_to_obs)
        status_label.value = f"Cells colored by cluster: {clustering_key}"
        apply_color_button.enabled = True

    # ── Transcript Overlay callbacks ──────────────────────────────────────────
    def on_apply_transcripts():
        if transcript_check.value:
            gene = transcript_gene_widget.value
            status_label.value = f"Loading transcripts for {gene}..."
            apply_transcripts_button.enabled = False

            @thread_worker(connect={"returned": _on_transcripts_ready})
            def fetch():
                return gene, transcript_loader.get_points_array(gene)

            fetch()
        else:
            transcript_layer.visible = False
            status_label.value = "Transcripts hidden"

    def _on_transcripts_ready(result):
        gene, points = result
        transcript_layer.data = points
        transcript_layer.visible = True
        status_label.value = f"Transcripts: {gene} ({len(points):,} spots)"
        apply_transcripts_button.enabled = True

    # ── Wire events ──────────────────────────────────────────────────────────
    mode_widget.changed.connect(on_mode_change)
    apply_color_button.clicked.connect(on_apply_color)
    apply_transcripts_button.clicked.connect(on_apply_transcripts)

    # ── Assemble container ───────────────────────────────────────────────────
    container = Container(
        widgets=[
            Label(value="─── Cell Coloring ───"),
            mode_widget,
            gene_widget,
            colormap_widget,
            clustering_widget,
            apply_color_button,
            Label(value=""),
            Label(value="─── Transcript Overlay ───"),
            transcript_gene_widget,
            transcript_check,
            qv_slider,
            apply_transcripts_button,
            Label(value=""),
            status_label,
        ]
    )
    return container


if __name__ == "__main__":
    main()
