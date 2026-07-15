"""ViewerContext — shared state dataclass passed to every tab module."""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class ViewerContext:
    """Holds all objects shared across tabs.

    Core objects (set once at startup):
    """
    # ── Core objects ──────────────────────────────────────────────────────────
    viewer: Any = None
    adata: Any = None
    sdata: Any = None
    clusterings: dict = field(default_factory=dict)
    color_manager: Any = None
    transcript_loader: Any = None
    umap_viewer: Any = None
    label_to_obs: np.ndarray | None = None
    centroids_yx: np.ndarray | None = None
    pixel_size: float = 1.0
    data_path: Path | None = None
    no_cache: bool = False

    # ── Gene / clustering name lists ─────────────────────────────────────────
    gene_names: list = field(default_factory=list)
    clustering_names: list = field(default_factory=list)

    # ── Napari layers ────────────────────────────────────────────────────────
    cell_labels_layer: Any = None
    transcript_layer: Any = None
    transcript_bins_layer: Any = None   # Image layer for transcript density heatmap
    roi_layer: Any = None
    annotation_layer: Any = None        # Named tissue annotation shapes
    crop_layer: Any = None              # Crop Dataset polygons (not session-persisted)

    # ── Morphology data (for coarse align) ──────────────────────────────────
    morph_thumb: Any = None
    morph_full_shape_yx: tuple | None = None

    # ── Mutable state dicts ──────────────────────────────────────────────────
    state: dict = field(default_factory=dict)
    he_state: dict = field(default_factory=dict)
    arms_state: dict = field(default_factory=dict)
    external_images_state: list = field(default_factory=list)
    patch_overlays_state: list = field(default_factory=list)

    # ── Cross-tab widget references (set by tabs that create them) ───────────
    clustering_widget: Any = None    # ComboBox — created by cell coloring tab
    gene_widget: Any = None          # ComboBox — created by cell coloring tab
    filter_check: Any = None         # CheckBox — created by cell coloring tab
    colormap_widget: Any = None      # ComboBox — created by cell coloring tab
    # Per-analysis-tab clustering widgets (set during tab construction)
    ga_clustering_widget: Any = None
    lr_clustering_widget: Any = None
    ne_clustering_widget: Any = None
    co_clustering_widget: Any = None
    annot_nhood_clustering_widget: Any = None
    annot_dist_clustering_widget: Any = None
    mg_clustering_widget: Any = None      # marker genes tab
    cnv_clustering_widget: Any = None     # CNV tab — reference-population picker

    # ── Shared helper callables (attached by create_shared_helpers) ──────────
    record_code: Any = None
    record_preamble: Any = None
    record_normalize: Any = None
    record_clustering: Any = None
    record_spatial_neighbors: Any = None
    refresh_clustering_choices: Any = None
    auto_save_plot: Any = None
    repopulate_cluster_checkboxes: Any = None
    get_selected_cluster_ids: Any = None
    make_cluster_mask: Any = None
    get_active_labels: Any = None
    get_labels_for: Any = None
    build_label_editor_dialog: Any = None
    apply_plot_font_size: Any = None
    get_cluster_ids_per_obs: Any = None
    translate_selected_ids_to_int: Any = None
    get_cluster_filter: Any = None
    set_status: Any = None

    # ── Cluster filter UI widgets (created by cell coloring tab) ─────────────
    cluster_scroll: Any = None
    cluster_filter_grid: Any = None
    select_all_btn: Any = None
    deselect_all_btn: Any = None

    # ── Dataset generation counter (incremented on dataset reload) ────────────
    dataset_generation: int = 0

    # ── Segmentation source ───────────────────────────────────────────────────
    segmentation_source: str = "xenium"   # "xenium" | "custom"
