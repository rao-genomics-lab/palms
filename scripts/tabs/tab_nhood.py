"""Tab 8: Neighborhood Enrichment."""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    ne_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        value=ctx.clustering_names[0] if ctx.clustering_names else None,
    )
    # Register on ctx for cross-tab refresh
    ctx.ne_clustering_widget = ne_clustering_widget

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
    ne_status = StatusProxy(ctx.viewer)

    from utils.gene_analysis import get_normalized_adata, add_clustering_to_obs
    from utils.spatial_analysis import (
        compute_spatial_neighbors, run_nhood_enrichment, make_nhood_enrichment_plot,
    )

    def on_run_nhood():
        ne_status.value = "Running neighborhood enrichment... (this may take a minute)"
        ne_run_button.enabled = False

        clustering_key = ne_clustering_widget.value
        n_perms = ne_perms_slider.value
        n_neighs = ne_neighs_slider.value
        state["_ne_params"] = {"n_perms": n_perms, "n_neighs": n_neighs}
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        @thread_worker(connect={"returned": _on_nhood_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, ctx.clusterings[clustering_key], clustering_key)
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            compute_spatial_neighbors(adata_norm, n_neighs=n_neighs)
            result = run_nhood_enrichment(adata_norm, clustering_key, n_perms=n_perms)
            result['_adata_norm'] = adata_norm
            result['_cluster_key'] = clustering_key
            return result
        _run()

    def _on_nhood_ready(result):
        state["nhood_result"] = result
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

        lines = [
            f"Neighborhood enrichment: {n}x{n} matrix ({n} clusters)",
            f"Clusters: {', '.join(clusters)}",
            "",
        ]

        if zscore.size > 0:
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

        _ne_ck = result.get('_cluster_key', '')
        _ne_p = state.get("_ne_params", {})
        _ne_np = _ne_p.get("n_perms", 1000)
        _ne_nn = _ne_p.get("n_neighs", 6)
        ctx.record_clustering(_ne_ck)
        ctx.record_spatial_neighbors(_ne_nn)
        ctx.record_code(
            f"\n# Neighborhood enrichment (n_perms={_ne_np})\n"
            f"sq.gr.nhood_enrichment(adata, cluster_key=\"{_ne_ck}\", "
            f"n_perms={_ne_np}, seed=42)"
        )

    def on_show_nhood_plot():
        result = state.get("nhood_result")
        if result is None:
            return
        groups = ctx.get_cluster_filter()
        ne_ck = result.get('_cluster_key', ne_clustering_widget.value)
        labels = ctx.get_labels_for(ne_ck)
        import matplotlib.pyplot as _plt
        ctx.apply_plot_font_size()
        try:
            if groups or labels:
                fig = make_nhood_enrichment_plot(
                    result, mode=ne_mode_widget.value,
                    cluster_filter=groups, cluster_labels=labels,
                )
            else:
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
                    fig = make_nhood_enrichment_plot(
                        result, mode=ne_mode_widget.value,
                    )
            state["nhood_fig"] = fig
            _plt.show(block=False)
            ne_save_plot_button.enabled = True
            if groups:
                ne_status.value = f"Heatmap displayed (clusters: {', '.join(groups)})"
            else:
                ne_status.value = "Heatmap displayed"

            _ne_mode = ne_mode_widget.value
            _ne_ck = result.get('_cluster_key', '')
            ctx.record_code(
                f"\n# Nhood enrichment heatmap (mode={_ne_mode})\n"
                f"sq.pl.nhood_enrichment(adata, cluster_key=\"{_ne_ck}\", "
                f"mode=\"{_ne_mode}\")\nplt.show()"
            )
        except Exception as e:
            ne_status.value = f"Plot error: {e}"

    def on_save_nhood_plot():
        fig = state.get("nhood_fig")
        if fig is None:
            return
        path = ctx.get_plot_save_path("Save Nhood Enrichment Plot", "nhood_enrichment")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        ne_status.value = f"Plot saved to {path}"
        ctx.record_code(f"\n# Save nhood enrichment plot\n# nhood_enrichment -> \"{path}\"")

    def on_export_nhood():
        result = state.get("nhood_result")
        if result is None:
            return
        import pandas as _pd
        zscore = result['zscore']
        clusters = result['clusters']
        df = _pd.DataFrame(zscore, index=clusters, columns=clusters)
        groups = ctx.get_cluster_filter()
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
        ctx.record_code(f"\n# Export nhood z-scores\n# nhood_zscore.csv -> \"{path}\"")

    ne_run_button.clicked.connect(on_run_nhood)
    ne_plot_button.clicked.connect(on_show_nhood_plot)
    ne_save_plot_button.clicked.connect(on_save_nhood_plot)
    ne_export_button.clicked.connect(on_export_nhood)

    ne_plot_btn_row = QWidget()
    ne_plot_btn_layout = QHBoxLayout()
    ne_plot_btn_layout.setContentsMargins(0, 0, 0, 0)
    ne_plot_btn_layout.addWidget(ne_plot_button.native)
    ne_plot_btn_layout.addWidget(ne_save_plot_button.native)
    ne_plot_btn_row.setLayout(ne_plot_btn_layout)

    widget = make_tab(
        ne_clustering_widget,
        ne_perms_slider,
        ne_neighs_slider,
        ne_run_button,
        ne_mode_widget,
        ne_results_text,
        ne_plot_btn_row,
        ne_export_button,
    )

    def _restore_session(session):
        nh = session.get("nhood_result")
        if nh is not None:
            state["nhood_result"] = nh
            ne_plot_button.enabled = True
            ne_export_button.enabled = True
            n = len(nh.get('clusters', []))
            print(f"  Restored nhood enrichment ({n} clusters)")

    return widget, {"ne_clustering_widget": ne_clustering_widget, "restore_session": _restore_session}
