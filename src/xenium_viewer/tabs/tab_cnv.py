"""Tab: CNV inference (InSituCNV / infercnvpy)."""

from __future__ import annotations
import re
from typing import TYPE_CHECKING

import pandas as pd
from magicgui.widgets import ComboBox, PushButton, SpinBox, FloatSpinBox
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QLabel, QScrollArea, QGridLayout, QCheckBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL

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

    # ── Cell types to analyze (limit CNV to a subset) ───────────────────
    cnv_analyze_label = QLabel("Cell types to analyze (CNV subclones):")
    cnv_analyze_hint = QLabel(
        "Only these cell types plus the reference are included in the analysis;\n"
        "leave all checked to analyze the whole tissue."
    )
    cnv_analyze_select_all_btn = PushButton(label="Select All")
    cnv_analyze_deselect_all_btn = PushButton(label="Deselect All")

    cnv_analyze_container = QWidget()
    cnv_analyze_grid = QGridLayout()
    cnv_analyze_grid.setContentsMargins(0, 0, 0, 0)
    cnv_analyze_container.setLayout(cnv_analyze_grid)

    cnv_analyze_scroll = QScrollArea()
    cnv_analyze_scroll.setWidget(cnv_analyze_container)
    cnv_analyze_scroll.setWidgetResizable(True)
    cnv_analyze_scroll.setMaximumHeight(150)

    state["cnv_analyze_checkboxes"] = {}

    def _cnv_analyze_select_all():
        for cb in state["cnv_analyze_checkboxes"].values():
            cb.setChecked(True)

    def _cnv_analyze_deselect_all():
        for cb in state["cnv_analyze_checkboxes"].values():
            cb.setChecked(False)

    cnv_analyze_select_all_btn.clicked.connect(_cnv_analyze_select_all)
    cnv_analyze_deselect_all_btn.clicked.connect(_cnv_analyze_deselect_all)

    def _repopulate_cnv_analyze_checkboxes():
        grid = cnv_analyze_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["cnv_analyze_checkboxes"].clear()

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
            cb.setChecked(True)  # analyze all cell types by default
            grid.addWidget(cb, i // cols, i % cols)
            state["cnv_analyze_checkboxes"][cid] = cb

    _repopulate_cnv_analyze_checkboxes()

    def _get_cnv_analyze_ids():
        """Return list of checked analyze cluster IDs (str)."""
        cbs = state["cnv_analyze_checkboxes"]
        return [str(cid) for cid, cb in cbs.items() if cb.isChecked()]

    def _all_cluster_ids(key):
        """All (str) cluster IDs present in clustering ``key``."""
        if not key or key not in ctx.clusterings:
            return set()
        return {str(x) for x in ctx.clusterings[key].dropna().unique().tolist()}

    # The reference grid is already rebuilt on dropdown change (see above);
    # rebuild the analyze grid on the same signal.
    cnv_clustering_widget.changed.connect(lambda _: _repopulate_cnv_analyze_checkboxes())

    # ── Parameters ──────────────────────────────────────────────────────
    # Defaults match InSituCNV's own reference notebook (run_insitucnv.ipynb).
    cnv_n_neighbors = SpinBox(label="Neighbors (expression graph)", min=5, max=100, value=15)
    cnv_smoothing_neighbors = SpinBox(label="Smoothing neighbors", min=5, max=200, value=20)
    cnv_window_size = SpinBox(label="Window size (genes)", min=2, max=200, value=60)
    cnv_step = SpinBox(label="Window step", min=1, max=50, value=10)
    cnv_resolution = FloatSpinBox(label="CNV cluster resolution", min=0.05, max=2.0, step=0.05, value=0.2)
    cnv_resolution.tooltip = (
        "InSituCNV's own notebook evaluates several resolutions (e.g. 0.1, 0.2, 0.3) "
        "and picks one per dataset after reviewing the results — this default may not "
        "be right for your data."
    )
    cnv_resolution_hint = QLabel(
        "Default may need tuning per dataset — check the chromosome heatmap and\n"
        "cluster count after running, and re-run with a different value if clusters\n"
        "look too coarse or too fragmented."
    )

    run_button = PushButton(label="Run CNV Inference", enabled=True)
    heatmap_res_widget = ComboBox(label="Heatmap resolution", choices=[], **combo_value_kwargs([]))
    heatmap_button = PushButton(label="Save Chromosome Heatmap (PDF/PNG)", enabled=False)
    score_color_button = PushButton(label="Color Cells by CNV Score", enabled=False)

    def _key_to_res_label(key: str) -> str:
        """'cnv_leiden_res0.2' -> 'res 0.2' (falls back to the raw key)."""
        m = re.search(r"res([0-9.]+)$", str(key))
        return f"res {m.group(1)}" if m else str(key)

    def _set_heatmap_choices(keys, select=None):
        """Populate the heatmap-resolution ComboBox from accumulated cluster keys.

        When a CNV profile is loaded, only keys actually present as columns in it
        are offered (so a stale/partial restore can't select a missing column).
        Choices are (label, value) tuples: friendly 'res 0.2' shown, raw key used.
        """
        keys = [k for k in (keys or []) if k]
        result = state.get("cnv_result")
        adata_cnv = result.get("adata_cnv") if result else None
        if adata_cnv is not None:
            cols = set(adata_cnv.obs.columns)
            keys = [k for k in keys if k in cols]
        choices = [(_key_to_res_label(k), k) for k in keys]
        heatmap_res_widget.choices = choices
        if select in keys:
            heatmap_res_widget.value = select
        elif keys:
            heatmap_res_widget.value = keys[-1]

    def _cnv_signature(result):
        """Core CNV parameters the profile depends on (resolution excluded).

        Must stay term-for-term identical to the signature recomputed in
        load_cnv_results_from_adata, including the trailing analyzed-type set.
        """
        p = result["params"]
        return (
            result["reference_clustering_name"],
            tuple(result["reference_categories"]),
            p.get("n_neighbors"),
            p.get("smoothing_neighbors"),
            p.get("window_size"),
            p.get("step"),
            p.get("lfc_clip"),
            tuple(sorted(result.get("analyze_categories") or [])),
        )

    def _res_from_key(key):
        m = re.search(r"res([0-9.]+)$", str(key))
        try:
            return float(m.group(1)) if m else None
        except (TypeError, ValueError):
            return None

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
        analyze_ids = _get_cnv_analyze_ids()
        if not analyze_ids:
            cnv_status.value = "Select at least one cell type to analyze"
            return
        # All cell types selected ⇒ no restriction (keeps default behavior / signature).
        analyze_categories = (
            None if set(analyze_ids) >= _all_cluster_ids(reference_key) else list(analyze_ids)
        )

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
                analyze_categories=analyze_categories,
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

        key = result["cluster_key"]
        series = result["cluster_series"]

        # Accumulate resolutions across runs that share the same core CNV params.
        # The CNV profile (X_cnv, gene positions) is identical for a fixed set of
        # core params, so we keep ONE retained adata_cnv and grow it by one obs
        # column per resolution. If core params change, that's a different profile
        # — reset the accumulated set (its heatmaps need the earlier profile).
        sig = _cnv_signature(result)
        prev = state.get("cnv_result")
        param_change_note = ""
        if (prev is not None and prev.get("signature") == sig
                and prev.get("adata_cnv") is not None):
            shared = prev["adata_cnv"]
            shared.obs[key] = (
                pd.Series(result["adata_cnv"].obs[key].values,
                          index=result["adata_cnv"].obs_names)
                .reindex(shared.obs_names).values
            )
            result["adata_cnv"] = shared
            result["cluster_keys"] = list(dict.fromkeys(prev.get("cluster_keys", []) + [key]))
        else:
            result["cluster_keys"] = [key]
            if prev is not None and prev.get("signature") is not None and prev.get("signature") != sig:
                param_change_note = (
                    " (CNV parameters changed — previous resolutions cleared; "
                    "their heatmaps need the earlier profile, re-run to restore them)"
                )
        result["signature"] = sig
        result["resolutions"] = sorted(
            {r for r in (_res_from_key(k) for k in result["cluster_keys"]) if r is not None}
        )
        state["cnv_result"] = result

        # Store the new clustering (invalidate stale color cache for this key first)
        ctx.color_manager._cluster_cache.pop(key, None)
        ctx.clusterings[key] = series
        if "custom_clusterings" not in state:
            state["custom_clusterings"] = {}
        state["custom_clusterings"][key] = series

        ctx.refresh_clustering_choices()
        _set_heatmap_choices(result["cluster_keys"], select=key)

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
        analyze_cats = result.get("analyze_categories") or []
        if analyze_cats:
            include_repr = repr(sorted(set(analyze_cats) | set(result["reference_categories"])))
            subset_line = (
                f"adata = adata[adata.obs['{reference_clustering_name}'].astype(str).isin({include_repr})].copy()"
                f"  # limit CNV to selected cell types + reference\n"
            )
        else:
            subset_line = ""
        ctx.record_node(
            "cnv",
            f"\n# CNV inference (InSituCNV / infercnvpy)\n"
            f"from insitucnv.tl import prepare_cnv_input, run_infercnv, compute_cnv_neighbors, cluster_cnv_resolutions\n"
            f"{subset_line}"
            f"adata.obs['{ref_obs_key}'] = adata.obs['{reference_clustering_name}']  # reference clustering\n"
            f"adata.layers['raw_counts'] = sdata['table'].X.copy()  # raw counts (pre-normalization)\n"
            f"sc.pp.normalize_total(adata); sc.pp.log1p(adata); sc.pp.pca(adata)\n"
            f"sc.pp.neighbors(adata, n_neighbors={params['n_neighbors']})\n"
            f"adata = prepare_cnv_input(adata, raw_layer='raw_counts', "
            f"smoothing_neighbors={params['smoothing_neighbors']}, add_gene_positions=True, drop_unmapped_genes=True, copy=False)\n"
            f"run_infercnv(adata, reference_key='{ref_obs_key}', reference_categories={ref_repr}, "
            f"window_size={params['window_size']}, step={params['step']}, calculate_gene_values=True, copy=False)\n"
            f"compute_cnv_neighbors(adata, copy=False)\n"
            f"cluster_cnv_resolutions(adata, {result['resolutions']!r}, copy=False)\n"
            f"# result: adata.obs[{result['cluster_keys']!r}]",
            deps=[f"clustering:{reference_clustering_name}"],
            label="CNV inference",
        )

        heatmap_button.enabled = True
        score_color_button.enabled = True

        res_line = ", ".join(str(r) for r in result["resolutions"])
        if analyze_cats:
            labels = ctx.get_labels_for(reference_clustering_name) if ctx.get_labels_for else {}
            analyze_line = ", ".join(
                str(labels.get(c, labels.get(str(c), c))) for c in analyze_cats
            )
            cells_line = f"  Cells analyzed: {result.get('n_cells')} (cell types: {analyze_line})\n"
        else:
            cells_line = f"  Cells analyzed: {result.get('n_cells')} (all cell types)\n"
        results_text.setPlainText(
            f"CNV inference complete{param_change_note}\n"
            f"  Reference clustering: {reference_clustering_name}\n"
            f"  Reference clusters: {', '.join(result['reference_categories'])}\n"
            f"{cells_line}"
            f"  Genes: {result['n_genes_mapped']} / {result['n_genes_total']} mapped to genome\n"
            f"  CNV windows: {result['n_windows']}\n"
            f"  CNV clusters found (res {params['resolution']}): {series.nunique()}\n"
            f"  CNV score range: {result['cnv_score'].min():.3f} - {result['cnv_score'].max():.3f}\n"
            f"\nResult stored as: {key}\n"
            f"  Resolutions available for heatmap: {res_line}\n"
            f"Use Cell Coloring tab to re-apply or switch colorings; pick a resolution\n"
            f"below to save its chromosome heatmap, or color cells by CNV score."
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

        cluster_key = heatmap_res_widget.value or result["cluster_key"]
        if cluster_key not in adata_cnv.obs.columns:
            cnv_status.value = (
                f"Resolution '{_key_to_res_label(cluster_key)}' is not in the current CNV "
                f"profile — re-run CNV inference at that resolution to save its heatmap"
            )
            return

        cnv_status.value = f"Building chromosome heatmap ({_key_to_res_label(cluster_key)})..."
        heatmap_button.enabled = False
        safe = re.sub(r"[^0-9A-Za-z._-]", "_", cluster_key)

        @thread_worker
        def _build():
            import os
            from xenium_viewer.utils.cnv_analysis import make_cnv_heatmap
            ctx.apply_plot_font_size()
            fig = make_cnv_heatmap(adata_cnv, cluster_key)
            plots_dir = os.path.join(ctx.data_path, "plots")
            os.makedirs(plots_dir, exist_ok=True)
            png_path = os.path.join(plots_dir, f"cnv_heatmap_{safe}.png")
            pdf_path = os.path.join(plots_dir, f"cnv_heatmap_{safe}.pdf")
            fig.savefig(png_path, dpi=300, bbox_inches="tight")
            fig.savefig(pdf_path, bbox_inches="tight")
            import matplotlib.pyplot as _plt
            _plt.close(fig)
            return png_path, pdf_path

        def _on_ready(paths):
            heatmap_button.enabled = True
            png_path, pdf_path = paths
            cnv_status.value = f"CNV chromosome heatmap saved to {png_path} and {pdf_path}"
            ctx.record_node(
                f"plot:cnv_heatmap:{cluster_key}",
                f"\n# CNV chromosome heatmap ({_key_to_res_label(cluster_key)})\n"
                f"import infercnvpy as cnv\n"
                f"cnv.pl.chromosome_heatmap(adata, groupby='{cluster_key}', show=False)\n"
                f"plt.savefig('cnv_heatmap_{safe}.png', dpi=300, bbox_inches='tight')\n"
                f"plt.savefig('cnv_heatmap_{safe}.pdf', bbox_inches='tight')",
                deps=["cnv"],
                kind=TERMINAL,
                label=f"CNV chromosome heatmap ({_key_to_res_label(cluster_key)})",
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

    cnv_analyze_btn_row = QWidget()
    cnv_analyze_btn_layout = QHBoxLayout()
    cnv_analyze_btn_layout.setContentsMargins(0, 0, 0, 0)
    cnv_analyze_btn_layout.addWidget(cnv_analyze_select_all_btn.native)
    cnv_analyze_btn_layout.addWidget(cnv_analyze_deselect_all_btn.native)
    cnv_analyze_btn_row.setLayout(cnv_analyze_btn_layout)

    widget = make_tab(
        cnv_clustering_widget,
        cnv_ref_label,
        cnv_sel_btn_row,
        cnv_cluster_scroll,
        cnv_analyze_label,
        cnv_analyze_hint,
        cnv_analyze_btn_row,
        cnv_analyze_scroll,
        cnv_n_neighbors,
        cnv_smoothing_neighbors,
        cnv_window_size,
        cnv_step,
        cnv_resolution,
        cnv_resolution_hint,
        run_button,
        cnv_progress,
        results_text,
        heatmap_res_widget,
        heatmap_button,
        score_color_button,
    )

    def _restore_session(session):
        result = state.get("cnv_result")
        if result is None:
            return
        heatmap_button.enabled = True
        score_color_button.enabled = True
        cluster_keys = result.get("cluster_keys") or [result.get("cluster_key")]
        _set_heatmap_choices(cluster_keys, select=result.get("cluster_key"))
        res_line = ", ".join(
            str(r) for r in sorted(
                {r for r in (_res_from_key(k) for k in cluster_keys) if r is not None}
            )
        )
        analyze_cats = result.get("analyze_categories") or []
        analyze_line = (
            f"  Analyzed cell types: {', '.join(analyze_cats)}\n" if analyze_cats
            else "  Analyzed cell types: all\n"
        )
        results_text.setPlainText(
            f"CNV inference (restored from previous session)\n"
            f"  Reference clustering: {result.get('reference_clustering_name')}\n"
            f"  Reference clusters: {', '.join(result.get('reference_categories', []))}\n"
            f"{analyze_line}"
            f"  Genes: {result.get('n_genes_mapped')} / {result.get('n_genes_total')} mapped to genome\n"
            f"  CNV windows: {result.get('n_windows')}\n"
            f"  Resolutions available for heatmap: {res_line or '(none)'}\n"
            f"\nUse Cell Coloring tab to re-apply the CNV clustering; pick a resolution\n"
            f"below to save its chromosome heatmap, or color cells by CNV score."
        )
        print(f"  Restored CNV results (cluster keys: {cluster_keys})")

    return widget, {"cnv_clustering_widget": cnv_clustering_widget, "restore_session": _restore_session}
