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
    plots_panel: Any = None          # the Plots dock gallery; set by app.py
    label_to_obs: np.ndarray | None = None
    centroids_yx: np.ndarray | None = None
    # The unfiltered table and its label map. ``adata``/``label_to_obs``
    # above are what the analysis is *about*, which Tools -> QC may narrow
    # to a subset; these stay pointed at every cell, and are what results
    # are persisted into. Equal to ``adata``/``label_to_obs`` while no
    # filter is in force. Set wherever those two are.
    full_adata: Any = None
    full_label_to_obs: np.ndarray | None = None
    umap_df: Any = None              # UMAP frame, kept so a rebind can re-index it
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
    # Display-only twin of the layer above, filled by the Transcripts tab's
    # preview. Never written by the recorded path; see tab_transcripts.py.
    transcript_bins_preview_layer: Any = None
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
    # Per-tab gene pickers. Bound to var_names at build time, so every one
    # of them needs re-populating when a gene filter shrinks the panel --
    # see ``refresh_gene_choices``.
    corr_gene_a_widget: Any = None
    corr_gene_b_widget: Any = None
    transcript_gene_widget: Any = None
    transcript_density_gene_widget: Any = None
    umap_gene_widget: Any = None
    cnv_clustering_widget: Any = None     # CNV tab — reference-population picker

    # ── Step execution (attached by create_shared_helpers) ───────────────────
    # ``executor`` is the StepExecutor whose namespace mirrors the exported
    # notebook's globals; ``run_step`` runs a Step through it and performs the
    # recorder side-effects. Prefer ``run_step`` over ``record_node`` in new
    # code — it is what makes the executed and recorded source identical.
    executor: Any = None
    run_step: Any = None
    # ``preview_step`` is the display-only sibling: it executes a Step's source
    # for a picture the user is only looking at, into a *copy* of the namespace,
    # and records nothing. It is not an alternative to ``run_step`` — it exists
    # so display code can reuse a template's text instead of keeping a second
    # implementation of the same computation. Its call sites are an allow-list
    # in ``tests/test_preview_never_records.py``; do not add one lightly.
    preview_step: Any = None

    # Reloads the open dataset from disk, rebuilding layers, managers and every
    # tab widget. Attached by app.py. Needed when something changes the zarr
    # store behind the viewer's back — the Cache tab's recovery writes elements
    # straight into it, so nothing in memory knows they exist until this runs.
    reload_dataset: Any = None

    # ── Shared helper callables (attached by create_shared_helpers) ──────────
    record_code: Any = None
    record_node: Any = None
    record_preamble: Any = None
    record_environment: Any = None  # versions + seeds; recorded with the preamble
    ensure_normalized: Any = None   # runs the "normalize" step, returns adata_norm
    ensure_spatial_neighbors: Any = None  # runs the "spatial_neighbors" step on adata_norm
    record_clustering: Any = None
    refresh_clustering_choices: Any = None
    refresh_gene_choices: Any = None      # re-populate every gene ComboBox
    # QC filtering (Tools -> QC). ``cell_root`` names the node a step that
    # reads cells must depend on -- "qc_filter" when a filter is in force,
    # "preamble" otherwise -- and is the one definition of it.
    cell_root: Any = None
    ensure_qc_filter: Any = None     # re-apply the stored cutoffs; idempotent
    apply_qc_filter: Any = None      # run + record a filter from a Preview
    clear_qc_filter: Any = None      # revert to every cell, and drop the node
    # The one route a figure takes: saves it under <data_path>/plots/ in every
    # configured format and shows it in the Plots dock. Replaced auto_save_plot,
    # which only ever wrote a file, in one format, from six of eighteen sites.
    show_plot: Any = None
    plot_paths: Any = None            # where a stem will be written
    # Shows the Plots dock, re-creating it if napari destroyed it (its
    # title-bar close button does). Bound by app.py.
    reveal_plots_dock: Any = None
    recorded_plot_paths: Any = None   # the same, relative to data_path
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
