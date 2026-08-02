"""Tab 9: Co-occurrence analysis."""

from __future__ import annotations
from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QFileDialog,
    QLabel, QScrollArea, QGridLayout, QCheckBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_tqdm_progress, qt_tqdm_context, make_progress_bar, combo_value_kwargs
from xenium_viewer.utils.prov_graph import ARTIFACT, TERMINAL
from xenium_viewer.utils.steps import Step, StepError, coerce
from xenium_viewer.utils.step_templates import (
    Preview, builtin_spec, builtin_text, step_template as _resolved,
)

TEMPLATE_ID = "spatial.cooccur"

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# Executed and recorded from one string. The old recorded cell ran on `adata`
# (raw counts, and it rebuilt obsm['spatial'] by hand from x_centroid/y_centroid)
# while the viewer ran on the normalised copy.
_COOCCUR_TEMPLATE = builtin_text("spatial.cooccur")


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    co_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.co_clustering_widget = co_clustering_widget

    co_interval_slider = Slider(label="Distance bins", min=10, max=100, value=50)
    co_run_button = PushButton(label="Run Co-occurrence", enabled=True)
    co_results_text = QTextEdit()
    co_results_text.setReadOnly(True)
    co_results_text.setFontFamily("monospace")
    co_results_text.setMaximumHeight(250)
    co_plot_button = PushButton(label="Show Co-occurrence Plot", enabled=False)
    co_export_button = PushButton(label="Export CSV...", enabled=False)
    co_filter_targets = CheckBox(label="Filter targets", value=False, enabled=True)

    # ── Cluster selector ─────────────────────────────────────────────────
    co_cluster_label = QLabel("Select clusters:")
    co_select_all_btn = PushButton(label="Select All")
    co_deselect_all_btn = PushButton(label="Deselect All")

    co_cluster_container = QWidget()
    co_cluster_grid = QGridLayout()
    co_cluster_grid.setContentsMargins(0, 0, 0, 0)
    co_cluster_container.setLayout(co_cluster_grid)

    co_cluster_scroll = QScrollArea()
    co_cluster_scroll.setWidget(co_cluster_container)
    co_cluster_scroll.setWidgetResizable(True)
    co_cluster_scroll.setMaximumHeight(150)

    state["co_cluster_checkboxes"] = {}

    def _co_select_all():
        for cb in state["co_cluster_checkboxes"].values():
            cb.setChecked(True)

    def _co_deselect_all():
        for cb in state["co_cluster_checkboxes"].values():
            cb.setChecked(False)

    co_select_all_btn.clicked.connect(_co_select_all)
    co_deselect_all_btn.clicked.connect(_co_deselect_all)

    def _repopulate_co_cluster_checkboxes():
        grid = co_cluster_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["co_cluster_checkboxes"].clear()

        key = co_clustering_widget.value
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
            cb.setChecked(True)
            grid.addWidget(cb, i // cols, i % cols)
            state["co_cluster_checkboxes"][cid] = cb

    co_clustering_widget.changed.connect(lambda _: _repopulate_co_cluster_checkboxes())
    _repopulate_co_cluster_checkboxes()

    def _get_co_selected_clusters():
        """Return sorted list of checked cluster IDs (str), or None if all checked."""
        cbs = state["co_cluster_checkboxes"]
        if not cbs:
            return None
        selected = [str(cid) for cid, cb in cbs.items() if cb.isChecked()]
        if len(selected) == len(cbs):
            return None
        return sorted(selected, key=lambda x: (int(x) if x.isdigit() else x))

    co_status = StatusProxy(ctx.viewer)
    co_progress = make_progress_bar()

    from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
    from xenium_viewer.utils.spatial_analysis import make_co_occurrence_plot

    def _cooccur_preview() -> Preview:
        """What "Run Co-occurrence" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane.
        """
        return Preview(
            list(builtin_spec(TEMPLATE_ID).blocks),
            {
                "cluster_key": co_clustering_widget.value,
                "interval": coerce(co_interval_slider.value),
            },
        )

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _cooccur_preview

    def on_run_co_occurrence():
        co_status.value = "Running co-occurrence analysis..."
        co_run_button.enabled = False

        blocks, params, _ = _cooccur_preview()
        clustering_key = params["cluster_key"]
        interval = params["interval"]
        state["_co_params"] = {"interval": interval}
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        # The clustering must exist as a node, and in adata.obs, before a step
        # can declare it as a dependency and read it. Both on the GUI thread.
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata, ctx.clusterings[clustering_key], clustering_key)

        step = Step(
            id=f"cooccur:{clustering_key}",
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=[f"clustering:{clustering_key}"],
            kind=ARTIFACT,
            label=f"Co-occurrence: {clustering_key}",
            outputs=["adata_norm"],
        )

        _progress = [None]  # filled after worker is created

        @thread_worker
        def _run():
            ctx.ensure_normalized()
            try:
                with qt_tqdm_context(_progress[0], "Co-occurrence: "):
                    adata_norm = ctx.run_step(step)["adata_norm"]
            except StepError as e:
                return {'occ': np.array([]), 'interval': np.array([]),
                        'clusters': [], 'warning': str(e)}
            uns = adata_norm.uns[f'{clustering_key}_co_occurrence']
            return {
                'occ': np.array(uns['occ']),
                'interval': np.array(uns['interval']),
                'clusters': list(adata_norm.obs[clustering_key].cat.categories.astype(str)),
                'warning': None,
                '_adata_norm': adata_norm,
                '_cluster_key': clustering_key,
            }

        worker = _run()
        _progress[0], state['_progress_timer'] = attach_tqdm_progress(
            worker,
            lambda m: setattr(co_status, 'value', m),
            "Co-occurrence: ",
            progress_bar=co_progress,
        )
        worker.returned.connect(_on_co_occurrence_ready)
        worker.start()

    def _on_co_occurrence_ready(result):
        state["co_result"] = result
        co_run_button.enabled = True

        # Persist to adata.uns immediately
        from xenium_viewer.utils.adata_persistence import save_co_occurrence_to_adata
        save_co_occurrence_to_adata(ctx, result)

        warning = result.get('warning')
        if warning:
            co_status.value = f"Co-occurrence: {warning}"
            co_results_text.setPlainText(warning)
            co_plot_button.enabled = False
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
            "Select clusters above to choose subplots. Use 'Filter targets' + Cell Coloring to restrict lines.",
        ]

        co_results_text.setPlainText("\n".join(lines))
        co_status.value = f"Co-occurrence done: {n} clusters"
        co_plot_button.enabled = True
        co_export_button.enabled = True

        # Recording happened inside ctx.run_step(), which recorded the source it ran.

    def on_show_co_plot():
        result = state.get("co_result")
        if result is None:
            return
        subplot_clusters = _get_co_selected_clusters()       # local checkboxes → subplots
        target_filter = ctx.get_cluster_filter() if co_filter_targets.value else None  # Cell Coloring → lines
        cc = state.get("cluster_to_color")
        co_ck = result.get('_cluster_key', co_clustering_widget.value)
        labels = ctx.get_labels_for(co_ck)
        import matplotlib.pyplot as _plt
        ctx.apply_plot_font_size()
        try:
            fig = make_co_occurrence_plot(
                result,
                clusters_to_plot=subplot_clusters,
                target_clusters=target_filter,
                cluster_colors=cc,
                cluster_labels=labels,
            )
            state["co_fig"] = fig
            _plt.show(block=False)
            path = ctx.auto_save_plot(fig, "co_occurrence")
            if subplot_clusters:
                co_status.value = f"Co-occurrence plot (subplots: {', '.join(subplot_clusters)}) — saved to {path}"
            else:
                co_status.value = f"Co-occurrence plot displayed — saved to {path}"

            _co_ck = result.get('_cluster_key', '')
            _co_fmt = ctx.state.get("plot_format", "svg")
            ctx.record_node(
                f"plot:cooccur:{_co_ck}",
                f"\n# Co-occurrence plot\n"
                f"sq.pl.co_occurrence(adata_norm, cluster_key=\"{_co_ck}\""
                + (f", clusters={subplot_clusters}" if subplot_clusters else "")
                + ")\n"
                + (f"# filter_targets={target_filter}\n" if target_filter else "")
                + f"fig = plt.gcf()\n"
                + f"fig.savefig(\"co_occurrence.{_co_fmt}\", dpi=300, bbox_inches='tight')",
                deps=[f"cooccur:{_co_ck}"],
                kind=TERMINAL,
                label="Co-occurrence plot",
            )
        except Exception as e:
            co_status.value = f"Plot error: {e}"

    def on_export_co():
        result = state.get("co_result")
        if result is None:
            return
        import pandas as _pd

        occ = result['occ']
        interval_arr = result['interval']
        clusters = result['clusters']
        distances = interval_arr[1:]

        groups = _get_co_selected_clusters()

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
        _co_ck = result.get('_cluster_key', '')
        _fname = os.path.basename(path)
        ctx.record_node(
            f"export:cooccur:{_co_ck}",
            f"\n# Export co-occurrence data\n"
            f"_res = adata_norm.uns[\"{_co_ck}_co_occurrence\"]\n"
            f"_occ, _iv = _res[\"occ\"], _res[\"interval\"]\n"
            f"_cats = list(adata_norm.obs[\"{_co_ck}\"].cat.categories)\n"
            f"_rows = [\n"
            f"    {{\"source_cluster\": s, \"target_cluster\": t, "
            f"\"distance\": _iv[k + 1], \"co_occurrence\": _occ[i, j, k]}}\n"
            f"    for i, s in enumerate(_cats) for j, t in enumerate(_cats) "
            f"for k in range(len(_iv) - 1)\n"
            f"]\n"
            f"pd.DataFrame(_rows).to_csv(\"{_fname}\", index=False)",
            deps=[f"cooccur:{_co_ck}"],
            kind=TERMINAL,
            label="Export co-occurrence data",
        )

    co_run_button.clicked.connect(on_run_co_occurrence)
    co_plot_button.clicked.connect(on_show_co_plot)
    co_export_button.clicked.connect(on_export_co)

    co_sel_btn_row = QWidget()
    co_sel_btn_layout = QHBoxLayout()
    co_sel_btn_layout.setContentsMargins(0, 0, 0, 0)
    co_sel_btn_layout.addWidget(co_select_all_btn.native)
    co_sel_btn_layout.addWidget(co_deselect_all_btn.native)
    co_sel_btn_row.setLayout(co_sel_btn_layout)

    widget = make_tab(
        co_clustering_widget,
        co_cluster_label,
        co_sel_btn_row,
        co_cluster_scroll,
        co_interval_slider,
        co_run_button,
        co_progress,
        co_results_text,
        co_filter_targets,
        co_plot_button,
        co_export_button,
    )

    def _restore_session(session):
        # Check state first (populated from adata.uns at startup), fall back to session
        co = state.get("co_result") or session.get("co_result")
        if co is not None:
            state["co_result"] = co
            co_plot_button.enabled = True
            co_export_button.enabled = True
            n = len(co.get('clusters', []))
            print(f"  Restored co-occurrence ({n} clusters)")

    return widget, {"co_clustering_widget": co_clustering_widget, "restore_session": _restore_session}
