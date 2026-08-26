"""Tab 8: Neighborhood Enrichment."""

from __future__ import annotations
from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_tqdm_progress, qt_tqdm_context, make_progress_bar, combo_value_kwargs
from xenium_viewer.utils.plot_output import recorded_save_code, safe_stem
from xenium_viewer.utils.prov_graph import ARTIFACT, TERMINAL
from xenium_viewer.utils.steps import Step, StepError, coerce
from xenium_viewer.utils.step_templates import (
    Preview, builtin_spec, builtin_text, step_template as _resolved,
)

TEMPLATE_ID = "spatial.nhood"

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# The source below is what the viewer executes *and* what the notebook records.
# It runs on ``adata_norm`` with the spatial graph the ``spatial_neighbors``
# step built on that same object; the old recorded cell called
# ``sq.gr.nhood_enrichment(adata, ...)`` on raw, graph-less counts.
_NHOOD_TEMPLATE = builtin_text("spatial.nhood")


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    ne_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
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
    ne_export_button = PushButton(label="Export Z-scores CSV...", enabled=False)
    ne_status = StatusProxy(ctx.viewer)
    ne_progress = make_progress_bar()

    from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
    from xenium_viewer.utils.spatial_analysis import make_nhood_enrichment_plot

    def _nhood_preview() -> Preview:
        """What "Run Nhood Enrichment" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane. ``n_neighs`` is not here: it belongs to
        the ``spatial_neighbors`` step this one depends on, not to this template.
        """
        clustering_key = ne_clustering_widget.value
        return Preview(
            list(builtin_spec(TEMPLATE_ID).blocks),
            {
                "cluster_key": clustering_key,
                "uns_key": f"{clustering_key}_nhood_enrichment",
                "n_perms": coerce(ne_perms_slider.value),
                "seed": 42,
            },
        )

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _nhood_preview

    def on_run_nhood():
        ne_status.value = "Running neighborhood enrichment..."
        ne_run_button.enabled = False

        blocks, params, _ = _nhood_preview()
        clustering_key = params["cluster_key"]
        n_perms = params["n_perms"]
        n_neighs = ne_neighs_slider.value
        state["_ne_params"] = {"n_perms": n_perms, "n_neighs": n_neighs}
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        # The clustering must exist as a node, and in adata.obs, before a step
        # can declare it as a dependency and read it. Both on the GUI thread.
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata, ctx.clusterings[clustering_key], clustering_key)

        step = Step(
            id=f"nhood:{clustering_key}",
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=[f"clustering:{clustering_key}", "spatial_neighbors"],
            kind=ARTIFACT,
            label=f"Neighborhood enrichment: {clustering_key}",
            outputs=["adata_norm"],
        )

        _progress = [None]  # filled after worker is created

        @thread_worker
        def _run():
            ctx.ensure_spatial_neighbors(n_neighs)   # implies ensure_normalized()
            try:
                with qt_tqdm_context(_progress[0], "Enrichment permutations: "):
                    adata_norm = ctx.run_step(step)["adata_norm"]
            except StepError as e:
                return {'zscore': np.array([]), 'count': np.array([]),
                        'clusters': [], 'warning': str(e)}
            uns = adata_norm.uns[f'{clustering_key}_nhood_enrichment']
            return {
                'zscore': np.array(uns['zscore']),
                'count': np.array(uns['count']),
                'clusters': list(adata_norm.obs[clustering_key].cat.categories.astype(str)),
                'warning': None,
                '_adata_norm': adata_norm,
                '_cluster_key': clustering_key,
            }

        worker = _run()
        _progress[0], state['_progress_timer'] = attach_tqdm_progress(
            worker,
            lambda m: setattr(ne_status, 'value', m),
            "Enrichment permutations: ",
            progress_bar=ne_progress,
        )
        worker.returned.connect(_on_nhood_ready)
        worker.start()

    def _on_nhood_ready(result):
        state["nhood_result"] = result
        ne_run_button.enabled = True

        # Persist to adata.uns immediately
        from xenium_viewer.utils.adata_persistence import save_nhood_to_adata
        save_nhood_to_adata(ctx, result)

        warning = result.get('warning')
        if warning:
            ne_status.value = f"Nhood: {warning}"
            ne_results_text.setPlainText(warning)
            ne_plot_button.enabled = False
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

        # Recording happened inside ctx.run_step(), which recorded the source it ran.

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
            _ne_ck = result.get('_cluster_key', '')
            paths = ctx.show_plot(
                fig, f"nhood_enrichment_{safe_stem(_ne_ck)}",
                title=f"Neighborhood enrichment: {_ne_ck}")
            if groups:
                ne_status.value = (f"Heatmap displayed (clusters: {', '.join(groups)}) "
                                   f"— saved to {', '.join(paths)}")
            else:
                ne_status.value = f"Heatmap displayed — saved to {', '.join(paths)}"

            _ne_mode = ne_mode_widget.value
            _saves = recorded_save_code(ctx.recorded_plot_paths(paths))
            # TERMINAL, still on the legacy recorder — plot nodes are the E4
            # view/analysis split. It reads adata_norm because that is where the
            # nhood step now puts the result. The savefig lines name the files
            # the viewer actually wrote, not a bare relative guess.
            ctx.record_node(
                f"plot:nhood:{_ne_ck}",
                f"\n# Nhood enrichment heatmap (mode={_ne_mode})\n"
                f"sq.pl.nhood_enrichment(adata_norm, cluster_key=\"{_ne_ck}\", "
                f"mode=\"{_ne_mode}\")\n"
                f"fig = plt.gcf()\n"
                f"{_saves}",
                deps=[f"nhood:{_ne_ck}"],
                kind=TERMINAL,
                label="Nhood enrichment heatmap",
            )
        except Exception as e:
            ne_status.value = f"Plot error: {e}"

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
        _ne_ck = result.get('_cluster_key', '')
        _fname = os.path.basename(path)
        ctx.record_node(
            f"export:nhood_zscore:{_ne_ck}",
            f"\n# Export nhood z-scores\n"
            f"_cats = adata_norm.obs[\"{_ne_ck}\"].cat.categories\n"
            f"pd.DataFrame(adata_norm.uns[\"{_ne_ck}_nhood_enrichment\"][\"zscore\"], "
            f"index=_cats, columns=_cats).to_csv(\"{_fname}\")",
            deps=[f"nhood:{_ne_ck}"],
            kind=TERMINAL,
            label="Export nhood z-scores",
        )

    ne_run_button.clicked.connect(on_run_nhood)
    ne_plot_button.clicked.connect(on_show_nhood_plot)
    ne_export_button.clicked.connect(on_export_nhood)

    widget = make_tab(
        ne_clustering_widget,
        ne_perms_slider,
        ne_neighs_slider,
        ne_run_button,
        ne_progress,
        ne_mode_widget,
        ne_results_text,
        ne_plot_button,
        ne_export_button,
    )

    def _restore_session(session):
        # Check state first (populated from adata.uns at startup), fall back to session
        nh = state.get("nhood_result") or session.get("nhood_result")
        if nh is not None:
            state["nhood_result"] = nh
            ne_plot_button.enabled = True
            ne_export_button.enabled = True
            n = len(nh.get('clusters', []))
            print(f"  Restored nhood enrichment ({n} clusters)")

    return widget, {"ne_clustering_widget": ne_clustering_widget, "restore_session": _restore_session}
