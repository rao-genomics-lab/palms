"""Tab: CNV inference (InSituCNV / infercnvpy)."""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, SpinBox, FloatSpinBox
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QLabel, QScrollArea, QGridLayout, QCheckBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar, combo_value_kwargs

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    cnv_clustering_widget = ComboBox(
        label="Reference clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.cnv_clustering_widget = cnv_clustering_widget

    cnv_ref_label = QLabel("Select reference (\"normal\") clusters:")
    cnv_select_all_btn = PushButton(label="Select All")
    cnv_deselect_all_btn = PushButton(label="Deselect All")

    cnv_cluster_container = QWidget()
    cnv_cluster_grid = QGridLayout()
    cnv_cluster_grid.setContentsMargins(0, 0, 0, 0)
    cnv_cluster_container.setLayout(cnv_cluster_grid)

    cnv_cluster_scroll = QScrollArea()
    cnv_cluster_scroll.setWidget(cnv_cluster_container)
    cnv_cluster_scroll.setWidgetResizable(True)
    cnv_cluster_scroll.setMaximumHeight(150)

    state["cnv_reference_checkboxes"] = {}

    def _cnv_select_all():
        for cb in state["cnv_reference_checkboxes"].values():
            cb.setChecked(True)

    def _cnv_deselect_all():
        for cb in state["cnv_reference_checkboxes"].values():
            cb.setChecked(False)

    cnv_select_all_btn.clicked.connect(_cnv_select_all)
    cnv_deselect_all_btn.clicked.connect(_cnv_deselect_all)

    def _repopulate_cnv_reference_checkboxes():
        grid = cnv_cluster_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["cnv_reference_checkboxes"].clear()

        key = cnv_clustering_widget.value
        if not key or key not in ctx.clusterings:
            return
        raw_ids = ctx.clusterings[key].dropna().unique().tolist()
        try:
            ids = sorted([int(x) for x in raw_ids])
        except (ValueError, TypeError):
            ids = sorted(raw_ids, key=lambda x: str(x))
        labels = ctx.get_labels_for(key) if ctx.get_labels_for else {}
        cols = 3
        for i, cid in enumerate(ids):
            display = str(labels.get(cid, labels.get(str(cid), cid)))
            cb = QCheckBox(display)
            cb.setChecked(False)  # reference population must be explicitly chosen
            grid.addWidget(cb, i // cols, i % cols)
            state["cnv_reference_checkboxes"][cid] = cb

    cnv_clustering_widget.changed.connect(lambda _: _repopulate_cnv_reference_checkboxes())
    _repopulate_cnv_reference_checkboxes()

    def _get_cnv_reference_ids():
        """Return list of checked reference cluster IDs (str)."""
        cbs = state["cnv_reference_checkboxes"]
        return [str(cid) for cid, cb in cbs.items() if cb.isChecked()]

    # ── Parameters ──────────────────────────────────────────────────────
    cnv_n_neighbors = SpinBox(label="Neighbors (expression graph)", min=5, max=100, value=15)
    cnv_smoothing_neighbors = SpinBox(label="Smoothing neighbors", min=5, max=200, value=30)
    cnv_window_size = SpinBox(label="Window size (genes)", min=2, max=200, value=10)
    cnv_step = SpinBox(label="Window step", min=1, max=50, value=2)
    cnv_resolution = FloatSpinBox(label="CNV cluster resolution", min=0.05, max=2.0, step=0.05, value=0.2)

    run_button = PushButton(label="Run CNV Inference", enabled=True)
    heatmap_button = PushButton(label="Show Chromosome Heatmap", enabled=False)
    score_color_button = PushButton(label="Color Cells by CNV Score", enabled=False)

    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(220)

    cnv_status = StatusProxy(ctx.viewer)
    cnv_progress = make_progress_bar()

    def _align_score_to_obs_order(score_series):
        adata = ctx.adata
        if 'cell_id' in adata.obs.columns:
            cell_ids = adata.obs['cell_id'].values
            aligned = score_series.reindex(cell_ids)
        else:
            aligned = score_series.reindex(adata.obs_names)
        return aligned.to_numpy(dtype=float)

    def on_run_cnv():
        if ctx.adata is None:
            cnv_status.value = "No AnnData loaded"
            return

        reference_key = cnv_clustering_widget.value
        reference_ids = _get_cnv_reference_ids()
        if not reference_key or reference_key not in ctx.clusterings:
            cnv_status.value = "Select a reference clustering first"
            return
        if not reference_ids:
            cnv_status.value = "Select at least one reference (\"normal\") cluster"
            return

        run_button.enabled = False
        results_text.setPlainText("Running CNV inference...")

        _adata = ctx.adata
        reference_series = ctx.clusterings[reference_key]
        n_neighbors = cnv_n_neighbors.value
        smoothing_neighbors = cnv_smoothing_neighbors.value
        window_size = cnv_window_size.value
        step = cnv_step.value
        resolution = cnv_resolution.value
        gen = ctx.dataset_generation

        @thread_worker
        def _run():
            from xenium_viewer.utils.cnv_analysis import run_cnv_pipeline
            return run_cnv_pipeline(
                _adata,
                reference_series,
                reference_ids,
                reference_clustering_name=reference_key,
                n_neighbors=n_neighbors,
                smoothing_neighbors=smoothing_neighbors,
                window_size=window_size,
                step=step,
                resolution=resolution,
            )

        worker = _run()
        _timer, _ = attach_spinner(
            worker,
            lambda m: setattr(cnv_status, "value", m),
            "Running CNV inference...",
            progress_bar=cnv_progress,
        )
        state["_cnv_spinner_timer"] = _timer  # keep reference to prevent GC

        worker.returned.connect(lambda result: _on_cnv_ready(result, gen))
        worker.errored.connect(_on_cnv_error)
        worker.start()

    def _on_cnv_ready(result, gen):
        run_button.enabled = True
        if ctx.dataset_generation != gen:
            return  # dataset reloaded while worker ran

        state["cnv_result"] = result
        key = result["cluster_key"]
        series = result["cluster_series"]

        # Store in clusterings (invalidate stale color cache for this key first)
        ctx.color_manager._cluster_cache.pop(key, None)
        ctx.clusterings[key] = series
        if "custom_clusterings" not in state:
            state["custom_clusterings"] = {}
        state["custom_clusterings"][key] = series

        ctx.refresh_clustering_choices()

        # Auto-apply cluster coloring in background thread
        @thread_worker
        def _apply_colors():
            color_arr, cluster_to_color = ctx.color_manager.get_cluster_colors(series)
            cluster_ids_per_obs, label_to_cluster = ctx.get_cluster_ids_per_obs(key)
            colormap = ctx.color_manager.build_direct_label_colormap(color_arr)
            return colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs

        def _on_colors_ready(color_result):
            colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs = color_result
            state["cluster_to_color"] = cluster_to_color
            state["label_to_cluster"] = label_to_cluster
            state["active_clustering_name"] = key
            if ctx.cell_labels_layer is not None:
                ctx.cell_labels_layer.colormap = colormap
                ctx.cell_labels_layer.refresh()
            ctx.umap_viewer.color_by_cluster(
                key, color_arr, ctx.label_to_obs,
                cluster_ids_per_obs=cluster_ids_per_obs,
            )
            cnv_status.value = f"CNV inference applied: {series.nunique()} CNV clusters"

        color_worker = _apply_colors()
        color_worker.returned.connect(_on_colors_ready)
        color_worker.start()

        # Persist
        from xenium_viewer.utils.adata_persistence import save_clustering_to_adata, save_cnv_results_to_adata
        save_clustering_to_adata(ctx, key, series)
        save_cnv_results_to_adata(ctx, result)

        # Record code
        reference_clustering_name = result["reference_clustering_name"]
        ctx.record_clustering(reference_clustering_name)
        params = result["params"]
        ref_obs_key = result["reference_obs_key"]
        ref_repr = repr(result["reference_categories"])
        ctx.record_code(
            f"\n# CNV inference (InSituCNV / infercnvpy)\n"
            f"from insitucnv.tl import prepare_cnv_input, run_infercnv, compute_cnv_neighbors, cluster_cnv_resolutions\n"
            f"adata.obs['{ref_obs_key}'] = adata.obs['{reference_clustering_name}']  # reference clustering\n"
            f"adata.layers['raw_counts'] = adata.X.copy()\n"
            f"sc.pp.normalize_total(adata); sc.pp.log1p(adata); sc.pp.pca(adata)\n"
            f"sc.pp.neighbors(adata, n_neighbors={params['n_neighbors']})\n"
            f"adata = prepare_cnv_input(adata, raw_layer='raw_counts', "
            f"smoothing_neighbors={params['smoothing_neighbors']}, add_gene_positions=True, drop_unmapped_genes=True, copy=False)\n"
            f"run_infercnv(adata, reference_key='{ref_obs_key}', reference_categories={ref_repr}, "
            f"window_size={params['window_size']}, step={params['step']}, calculate_gene_values=True, copy=False)\n"
            f"compute_cnv_neighbors(adata, copy=False)\n"
            f"cluster_cnv_resolutions(adata, [{params['resolution']}], copy=False)\n"
            f"# result: adata.obs['{key}']",
            tag="cnv_inference",
        )

        heatmap_button.enabled = True
        score_color_button.enabled = True

        results_text.setPlainText(
            f"CNV inference complete\n"
            f"  Reference clustering: {reference_clustering_name}\n"
            f"  Reference clusters: {', '.join(result['reference_categories'])}\n"
            f"  Genes: {result['n_genes_mapped']} / {result['n_genes_total']} mapped to genome\n"
            f"  CNV windows: {result['n_windows']}\n"
            f"  CNV clusters found: {series.nunique()}\n"
            f"  CNV score range: {result['cnv_score'].min():.3f} - {result['cnv_score'].max():.3f}\n"
            f"\nResult stored as: {key}\n"
            f"Use Cell Coloring tab to re-apply or switch colorings, or the buttons\n"
            f"below to view the chromosome heatmap / color by CNV score."
        )

    def _on_cnv_error(exc):
        run_button.enabled = True
        cnv_status.value = f"CNV inference error: {exc}"
        results_text.setPlainText(f"Error running CNV inference:\n{exc}")

    def on_show_heatmap():
        result = state.get("cnv_result")
        if result is None:
            return
        adata_cnv = result.get("adata_cnv")
        if adata_cnv is None:
            cnv_status.value = "No cached CNV profile available — rerun CNV inference to view the heatmap"
            return
        cnv_status.value = "Building chromosome heatmap..."
        heatmap_button.enabled = False

        cluster_key = result["cluster_key"]

        @thread_worker
        def _build():
            from xenium_viewer.utils.cnv_analysis import make_cnv_heatmap
            ctx.apply_plot_font_size()
            return make_cnv_heatmap(adata_cnv, cluster_key)

        def _on_ready(fig):
            heatmap_button.enabled = True
            import matplotlib.pyplot as _plt
            state["cnv_heatmap_fig"] = fig
            _plt.show(block=False)
            path = ctx.auto_save_plot(fig, "cnv_heatmap")
            cnv_status.value = f"CNV chromosome heatmap displayed — saved to {path}"
            ctx.record_code(
                f"\n# CNV chromosome heatmap\n"
                f"import infercnvpy as cnv\n"
                f"cnv.pl.chromosome_heatmap(adata, groupby='{cluster_key}')"
            )

        def _on_error(exc):
            heatmap_button.enabled = True
            cnv_status.value = f"Heatmap error: {exc}"

        worker = _build()
        worker.returned.connect(_on_ready)
        worker.errored.connect(_on_error)
        worker.start()

    def on_color_by_score():
        result = state.get("cnv_result")
        if result is None:
            return
        score_series = result.get("cnv_score")
        if score_series is None:
            cnv_status.value = "No CNV score available — rerun CNV inference"
            return
        cnv_status.value = "Coloring cells by CNV score..."
        score_color_button.enabled = False

        @thread_worker
        def _build():
            values = _align_score_to_obs_order(score_series)
            color_arr = ctx.color_manager.get_continuous_colors(
                values, colormap="viridis", cache_key="cnv_score",
            )
            colormap = ctx.color_manager.build_direct_label_colormap(color_arr)
            return colormap

        def _on_ready(colormap):
            score_color_button.enabled = True
            if ctx.cell_labels_layer is not None:
                ctx.cell_labels_layer.colormap = colormap
                ctx.cell_labels_layer.refresh()
            cnv_status.value = "Cells colored by CNV score (viridis)"

        def _on_error(exc):
            score_color_button.enabled = True
            cnv_status.value = f"CNV score coloring error: {exc}"

        worker = _build()
        worker.returned.connect(_on_ready)
        worker.errored.connect(_on_error)
        worker.start()

    run_button.clicked.connect(on_run_cnv)
    heatmap_button.clicked.connect(on_show_heatmap)
    score_color_button.clicked.connect(on_color_by_score)

    cnv_sel_btn_row = QWidget()
    cnv_sel_btn_layout = QHBoxLayout()
    cnv_sel_btn_layout.setContentsMargins(0, 0, 0, 0)
    cnv_sel_btn_layout.addWidget(cnv_select_all_btn.native)
    cnv_sel_btn_layout.addWidget(cnv_deselect_all_btn.native)
    cnv_sel_btn_row.setLayout(cnv_sel_btn_layout)

    widget = make_tab(
        cnv_clustering_widget,
        cnv_ref_label,
        cnv_sel_btn_row,
        cnv_cluster_scroll,
        cnv_n_neighbors,
        cnv_smoothing_neighbors,
        cnv_window_size,
        cnv_step,
        cnv_resolution,
        run_button,
        cnv_progress,
        results_text,
        heatmap_button,
        score_color_button,
    )

    def _restore_session(session):
        result = state.get("cnv_result")
        if result is None:
            return
        heatmap_button.enabled = True
        score_color_button.enabled = True
        results_text.setPlainText(
            f"CNV inference (restored from previous session)\n"
            f"  Reference clustering: {result.get('reference_clustering_name')}\n"
            f"  Reference clusters: {', '.join(result.get('reference_categories', []))}\n"
            f"  Genes: {result.get('n_genes_mapped')} / {result.get('n_genes_total')} mapped to genome\n"
            f"  CNV windows: {result.get('n_windows')}\n"
            f"  CNV cluster key: {result.get('cluster_key')}\n"
            f"\nUse Cell Coloring tab to re-apply the CNV clustering, or the buttons\n"
            f"below to view the chromosome heatmap / color by CNV score."
        )
        print(f"  Restored CNV results (cluster key: {result.get('cluster_key')})")

    return widget, {"cnv_clustering_widget": cnv_clustering_widget, "restore_session": _restore_session}
