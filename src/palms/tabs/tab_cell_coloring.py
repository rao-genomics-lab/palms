"""Tab 1: Cell Coloring — gene/cluster color mode, filter checkboxes."""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, RadioButtons, Slider
from qtpy.QtWidgets import (
    QWidget, QHBoxLayout, QScrollArea, QGridLayout,
)
from napari.qt.threading import thread_worker
from palms.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from palms.utils.prov_graph import NOTE, TERMINAL

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext

from palms.utils.coloring import AVAILABLE_COLORMAPS


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    mode_widget = RadioButtons(
        label="Color cells by",
        choices=["Gene Expression", "Cluster"],
        value="Gene Expression",
    )

    gene_widget = ComboBox(
        label="Gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )

    colormap_widget = ComboBox(
        label="Colormap",
        choices=AVAILABLE_COLORMAPS,
        value="viridis",
    )

    clustering_widget = ComboBox(
        label="Clustering",
        choices=ctx.clustering_names,
        enabled=False,
        **combo_value_kwargs(ctx.clustering_names),
    )

    filter_check = CheckBox(label="Filter by cluster", value=False, enabled=True)

    min_cells_widget = Slider(
        label="Min cluster size", min=100, max=10000, value=500, step=100
    )
    filter_small_btn = PushButton(label="Filter Small Clusters", enabled=True)

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

    state["cluster_checkboxes"] = {}

    bg_white_check = CheckBox(label="White background", value=False)
    apply_color_button = PushButton(label="Apply Cell Coloring", enabled=True)
    edit_labels_button = PushButton(label="Edit Cluster Labels...", enabled=True)

    status_label = StatusProxy(ctx.viewer)

    # Register widgets on ctx for cross-tab access
    ctx.clustering_widget = clustering_widget
    ctx.gene_widget = gene_widget
    ctx.filter_check = filter_check
    ctx.colormap_widget = colormap_widget
    ctx.cluster_scroll = cluster_scroll
    ctx.cluster_filter_grid = cluster_filter_grid
    ctx.select_all_btn = select_all_btn
    ctx.deselect_all_btn = deselect_all_btn

    # ── Select All / Deselect All ────────────────────────────────────────
    def _on_select_all():
        for cb in state["cluster_checkboxes"].values():
            cb.setChecked(True)

    def _on_deselect_all():
        for cb in state["cluster_checkboxes"].values():
            cb.setChecked(False)

    select_all_btn.clicked.connect(_on_select_all)
    deselect_all_btn.clicked.connect(_on_deselect_all)

    # ── Callbacks ────────────────────────────────────────────────────────
    def _set_cluster_filter_enabled(enabled):
        cluster_scroll.setEnabled(enabled)
        select_all_btn.enabled = enabled
        deselect_all_btn.enabled = enabled
        for cb in state["cluster_checkboxes"].values():
            cb.setEnabled(enabled)

    def on_mode_change(value):
        state["color_mode"] = value
        is_gene = (value == "Gene Expression")
        gene_widget.enabled = is_gene
        colormap_widget.enabled = is_gene
        clustering_widget.enabled = (value == "Cluster") or filter_check.value
        _set_cluster_filter_enabled(filter_check.value)

    def on_filter_change(value):
        state["filter_by_cluster"] = value
        clustering_widget.enabled = (state["color_mode"] == "Cluster") or value
        _set_cluster_filter_enabled(value)

    def on_clustering_change(value):
        ctx.repopulate_cluster_checkboxes()
        for combo in [ctx.ga_clustering_widget, ctx.lr_clustering_widget,
                      ctx.ne_clustering_widget, ctx.co_clustering_widget]:
            if combo is not None and value in [c for c in combo.choices]:
                combo.value = value

    def _on_gene_colors_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        gene, color_arr, selected_ids, clustering_key = result
        if selected_ids is not None and clustering_key:
            _, label_to_cluster_arr = ctx.get_cluster_ids_per_obs(clustering_key)
            int_ids = ctx.translate_selected_ids_to_int(selected_ids)
            color_arr = color_arr.copy()
            mask_out = ~np.isin(label_to_cluster_arr, int_ids)
            valid_range = min(len(mask_out), len(color_arr))
            color_arr[:valid_range][mask_out[:valid_range]] = 0
            filter_desc = f" (clusters: {sorted(selected_ids)})"
        else:
            filter_desc = ""
        ctx.color_manager.apply_to_labels_layer(ctx.cell_labels_layer, color_arr)
        ctx.umap_viewer.color_by_gene(gene, color_arr, ctx.label_to_obs)
        state["label_to_cluster"] = None
        state["active_clustering_name"] = None
        status_label.value = f"Cells colored by gene: {gene}{filter_desc}"
        ctx.record_node(
            "plot:spatial_gene",
            f"\n# Spatial plot colored by gene expression"
            + (f"  (viewer filter: {filter_desc})" if filter_desc else "") + "\n"
            f"sc.pl.embedding(adata, basis=\"spatial\", color=\"{gene}\", "
            f"cmap=\"{state['current_colormap']}\")",
            deps=["preamble"],
            kind=TERMINAL,
            label=f"Spatial plot: {gene}",
        )
        apply_color_button.enabled = True

    def _on_cluster_colors_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        (clustering_key, colormap, color_arr, cluster_to_color,
         label_to_cluster, cluster_ids_per_obs, selected_ids) = result

        state["cluster_to_color"] = cluster_to_color
        state["label_to_cluster"] = label_to_cluster
        state["active_clustering_name"] = clustering_key

        ctx.cell_labels_layer.colormap = colormap
        ctx.cell_labels_layer.refresh()

        ctx.umap_viewer.color_by_cluster(
            clustering_key, color_arr, ctx.label_to_obs,
            cluster_ids_per_obs=cluster_ids_per_obs,
        )
        filter_desc = f" (clusters: {sorted(selected_ids)})" if selected_ids is not None else ""
        status_label.value = f"Cells colored by cluster: {clustering_key}{filter_desc}"
        ctx.record_clustering(clustering_key)
        ctx.record_node(
            "plot:spatial_cluster",
            f"\n# Spatial plot colored by cluster"
            + (f"  (viewer filter: clusters {sorted(selected_ids)})" if selected_ids else "")
            + "\n"
            f"sc.pl.embedding(adata, basis=\"spatial\", color=\"{clustering_key}\")",
            deps=[f"clustering:{clustering_key}"],
            kind=TERMINAL,
            label=f"Spatial plot: {clustering_key}",
        )
        apply_color_button.enabled = True

    def on_apply_color():
        if ctx.cell_labels_layer is None:
            status_label.value = "No cell_labels layer found"
            return

        mode = state["color_mode"]
        status_label.value = "Computing cell colors..."
        apply_color_button.enabled = False
        gen = ctx.dataset_generation

        if mode == "Gene Expression":
            gene = gene_widget.value
            cmap = colormap_widget.value
            state["current_gene"] = gene
            state["current_colormap"] = cmap

            use_filter = filter_check.value
            selected_ids = ctx.get_selected_cluster_ids() if use_filter else None
            c_key = clustering_widget.value if use_filter else None

            @thread_worker
            def compute_gene():
                color_arr = ctx.color_manager.get_gene_colors(gene, colormap=cmap)
                return gene, color_arr, selected_ids, c_key

            worker = compute_gene()
            worker.returned.connect(lambda result: _on_gene_colors_ready(result, gen))
            worker.start()

        else:  # Cluster
            clustering_key = clustering_widget.value
            state["current_clustering"] = clustering_key
            cluster_series = ctx.clusterings[clustering_key]
            cluster_series.name = clustering_key

            use_filter = filter_check.value
            selected_ids = ctx.get_selected_cluster_ids() if use_filter else None

            @thread_worker
            def compute_cluster():
                color_arr, cluster_to_color = ctx.color_manager.get_cluster_colors(cluster_series)
                cluster_ids_per_obs, label_to_cluster = ctx.get_cluster_ids_per_obs(clustering_key)

                if selected_ids is not None:
                    int_ids = ctx.translate_selected_ids_to_int(selected_ids)
                    color_arr = color_arr.copy()
                    mask_out = ~np.isin(label_to_cluster, int_ids)
                    valid_range = min(len(mask_out), len(color_arr))
                    color_arr[:valid_range][mask_out[:valid_range]] = 0

                colormap = ctx.color_manager.build_direct_label_colormap(color_arr)
                return (clustering_key, colormap, color_arr, cluster_to_color,
                        label_to_cluster, cluster_ids_per_obs, selected_ids)

            worker = compute_cluster()
            worker.returned.connect(lambda result: _on_cluster_colors_ready(result, gen))
            worker.start()

    def on_bg_change(value):
        ctx.viewer.window._qt_viewer.canvas.bgcolor = (1, 1, 1, 1) if value else (0, 0, 0, 1)
        ctx.record_node(
            "viewer:background",
            f"\n# Viewer background set to {'white' if value else 'black'} (display only)",
            deps=["preamble"],
            kind=NOTE,
            label="Viewer background",
        )

    def _on_filter_small_clusters():
        clustering_key = clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            status_label.value = "No clustering selected"
            return
        threshold = min_cells_widget.value
        series = ctx.clusterings[clustering_key]
        counts = series.value_counts()
        # Enable filter checkbox so the filter is actually applied downstream
        filter_check.value = True
        # Repopulate checkboxes if not yet populated
        if not state["cluster_checkboxes"]:
            ctx.repopulate_cluster_checkboxes()
        # Uncheck clusters smaller than threshold
        n_excluded = 0
        for cid, cb in state["cluster_checkboxes"].items():
            cell_count = counts.get(cid, counts.get(str(cid), 0))
            include = int(cell_count) >= threshold
            cb.setChecked(include)
            if not include:
                n_excluded += 1
        status_label.value = (
            f"Size filter: {n_excluded} cluster(s) with < {threshold} cells excluded"
        )
        ctx.record_node(
            "viewer:size_filter",
            f"\n# Cluster size filter (viewer): min_cells={threshold}, "
            f"{n_excluded} cluster(s) excluded from the display",
            deps=["preamble"],
            kind=NOTE,
            label="Cluster size filter",
        )

    def _on_edit_labels():
        clustering_key = clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            status_label.value = "No clustering selected"
            return
        if ctx.build_label_editor_dialog(clustering_key):
            status_label.value = f"Labels updated for {clustering_key}"
            ctx.repopulate_cluster_checkboxes()

    # ── Wire events ──────────────────────────────────────────────────────
    mode_widget.changed.connect(on_mode_change)
    filter_check.changed.connect(on_filter_change)
    clustering_widget.changed.connect(on_clustering_change)
    apply_color_button.clicked.connect(on_apply_color)
    bg_white_check.changed.connect(on_bg_change)
    edit_labels_button.clicked.connect(_on_edit_labels)
    filter_small_btn.clicked.connect(_on_filter_small_clusters)

    # Select All / Deselect All buttons row
    cluster_btn_row = QWidget()
    cluster_btn_layout = QHBoxLayout()
    cluster_btn_layout.setContentsMargins(0, 0, 0, 0)
    cluster_btn_layout.addWidget(select_all_btn.native)
    cluster_btn_layout.addWidget(deselect_all_btn.native)
    cluster_btn_row.setLayout(cluster_btn_layout)

    widget = make_tab(
        bg_white_check,
        mode_widget,
        gene_widget,
        colormap_widget,
        clustering_widget,
        filter_check,
        min_cells_widget,
        filter_small_btn,
        cluster_btn_row,
        cluster_scroll,
        apply_color_button,
        edit_labels_button,
    )

    return widget, {
        "clustering_widget": clustering_widget,
        "gene_widget": gene_widget,
        "filter_check": filter_check,
    }
