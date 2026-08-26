"""Tab: Distance to Annotation.

For each Xenium cell, computes the minimum Euclidean distance (µm) to the
boundary of a selected annotation type, then visualises the distribution per
cell-type cluster as violin / box plots.  Optionally colours cells on the
napari canvas by their distance value.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, PushButton, SpinBox
from xenium_viewer.utils.coloring import AVAILABLE_COLORMAPS
from qtpy.QtWidgets import (
    QTextEdit, QFileDialog,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import StatusProxy, combo_value_kwargs, make_tab
from xenium_viewer.utils.plot_output import safe_stem
from xenium_viewer.utils.prov_graph import NOTE

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Selectors ─────────────────────────────────────────────────────────────
    clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.annot_dist_clustering_widget = clustering_widget

    def _annot_choices():
        from xenium_viewer.utils.annotation_utils import get_annotation_types
        return get_annotation_types(ctx.annotation_layer) or ["(none)"]

    annot_widget = ComboBox(label="Annotation type", choices=_annot_choices())
    refresh_btn = PushButton(label="Refresh annotation types")

    # ── Parameters ────────────────────────────────────────────────────────────
    plot_type_widget = ComboBox(
        label="Plot type", choices=["violin", "box", "strip"], value="violin",
    )
    max_dist_spin = SpinBox(
        label="Max distance to show (µm, 0=all)", min=0, max=100000, value=0,
    )

    dist_colormap_widget = ComboBox(
        label="Distance colormap", choices=AVAILABLE_COLORMAPS, value="plasma",
    )

    # ── Controls ──────────────────────────────────────────────────────────────
    run_btn = PushButton(label="Run Distance Analysis", enabled=True)
    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(200)
    plot_btn = PushButton(label="Show Plot", enabled=False)
    export_btn = PushButton(label="Export CSV...", enabled=False)
    color_btn = PushButton(label="Colour cells by distance", enabled=False)
    clear_color_btn = PushButton(label="Clear distance colouring", enabled=False)

    status = StatusProxy(ctx.viewer)

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def _refresh_annot():
        annot_widget.choices = _annot_choices()

    def _on_run():
        from xenium_viewer.utils.annotation_utils import compute_distance_to_annotation

        clustering_key = clustering_widget.value
        annot_type = annot_widget.value
        if not clustering_key or annot_type in (None, "(none)"):
            results_text.setPlainText("Select a clustering and annotation type.")
            return
        if ctx.annotation_layer is None or not ctx.annotation_layer.data:
            results_text.setPlainText("No annotations defined.")
            return

        run_btn.enabled = False
        status.value = f"Computing distances to '{annot_type}'..."

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        centroids_um = np.asarray(_adata.obsm['spatial'], dtype=np.float64)  # (N,2) xy µm

        @thread_worker
        def _run():
            distances = compute_distance_to_annotation(
                centroids_um, ctx.annotation_layer, annot_type, ctx.pixel_size
            )
            return distances

        def _on_done(distances):
            run_btn.enabled = True
            state["annot_dist_distances"] = distances
            state["annot_dist_annot_type"] = annot_type
            state["annot_dist_clustering_key"] = clustering_key

            valid = distances[~np.isnan(distances)]
            if len(valid) == 0:
                results_text.setPlainText(
                    f"No cells could be measured against annotation '{annot_type}'.\n"
                    "Ensure the annotation layer contains polygons of this type."
                )
                return

            lines = [
                f"Distance to '{annot_type}' boundary (µm)",
                f"Clustering: {clustering_key}",
                f"Cells measured: {len(valid):,}",
                f"Overall: min={valid.min():.1f}  median={np.median(valid):.1f}  "
                f"mean={valid.mean():.1f}  max={valid.max():.1f}",
                "",
            ]

            # Per-cluster summary
            cluster_series = ctx.clusterings.get(clustering_key)
            if cluster_series is not None:
                if 'cell_id' in _adata.obs.columns:
                    cell_ids = _adata.obs['cell_id'].values
                    aligned = cluster_series.reindex(cell_ids).values
                else:
                    aligned = cluster_series.reindex(_adata.obs_names).values

                unique_clusters = sorted(set(str(c) for c in aligned if str(c) != 'nan'))
                lines.append("Per-cluster median distance (µm):")
                for c in unique_clusters:
                    mask = np.array([str(v) == c for v in aligned])
                    d_c = distances[mask]
                    valid_c = d_c[~np.isnan(d_c)]
                    if len(valid_c) > 0:
                        lines.append(f"  {c}: {np.median(valid_c):.1f} µm  (n={len(valid_c):,})")

            results_text.setPlainText("\n".join(lines))
            status.value = f"Distance analysis done: {len(valid):,} cells measured"
            plot_btn.enabled = True
            export_btn.enabled = True
            color_btn.enabled = True

        worker = _run()
        worker.returned.connect(_on_done)
        worker.start()

    def _on_show_plot():
        distances = state.get("annot_dist_distances")
        clustering_key = state.get("annot_dist_clustering_key")
        annot_type = state.get("annot_dist_annot_type")
        if distances is None or clustering_key is None:
            return

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        cluster_series = ctx.clusterings.get(clustering_key)
        if cluster_series is None:
            return

        if 'cell_id' in _adata.obs.columns:
            cell_ids = _adata.obs['cell_id'].values
            aligned = cluster_series.reindex(cell_ids).values
        else:
            aligned = cluster_series.reindex(_adata.obs_names).values

        labels = ctx.get_labels_for(clustering_key)
        unique_clusters = sorted(set(str(c) for c in aligned if str(c) != 'nan'))

        max_d = max_dist_spin.value
        plot_type = plot_type_widget.value

        ctx.apply_plot_font_size()

        import matplotlib.pyplot as plt
        import pandas as pd

        # Build DataFrame for plotting
        rows = []
        for i, (c, d) in enumerate(zip(aligned, distances)):
            if np.isnan(d):
                continue
            if max_d > 0 and d > max_d:
                continue
            c_str = str(c)
            display = labels.get(c_str, labels.get(c, c_str)) if labels else c_str
            rows.append({"cluster": str(display), "distance_um": d})

        if not rows:
            status.value = "No data to plot."
            return

        df = pd.DataFrame(rows)
        cluster_order = [
            str(labels.get(c, labels.get(str(c), c))) if labels else str(c)
            for c in unique_clusters
        ]
        cluster_order = [c for c in cluster_order if c in df['cluster'].values]

        fig, ax = plt.subplots(figsize=(max(8, len(cluster_order) * 0.6), 5))

        if plot_type == "violin":
            import seaborn as sns
            sns.violinplot(
                data=df, x="cluster", y="distance_um", order=cluster_order,
                ax=ax, cut=0, inner="quartile", palette="tab20",
            )
        elif plot_type == "box":
            import seaborn as sns
            sns.boxplot(
                data=df, x="cluster", y="distance_um", order=cluster_order,
                ax=ax, palette="tab20", flierprops={"markersize": 2},
            )
        else:  # strip
            import seaborn as sns
            sns.stripplot(
                data=df, x="cluster", y="distance_um", order=cluster_order,
                ax=ax, size=2, alpha=0.5, palette="tab20",
            )

        ax.set_xlabel(clustering_key)
        ax.set_ylabel("Distance to annotation boundary (µm)")
        ax.set_title(f"Distance to '{annot_type}' boundary by {clustering_key}")
        ax.tick_params(axis="x", labelrotation=45)
        for tick in ax.get_xticklabels():
            tick.set_horizontalalignment("right")
        fig.tight_layout()
        # ``plt.show()`` here was *blocking*, and this figure was never recorded.
        stem = f"annot_distance_{safe_stem(annot_type)}_{safe_stem(clustering_key)}"
        paths = ctx.show_plot(
            fig, stem,
            title=f"Distance to '{annot_type}' by {clustering_key}")
        status.value = f"Distance plot displayed — saved to {', '.join(paths)}"
        # NOTE for the same reason as the annotation-nhood plot: the distances
        # are measured against shapes drawn in the viewer, which the notebook
        # cannot reach, so there is no code to record — only the fact that this
        # figure exists and where it went.
        ctx.record_node(
            "viewer:annot_distance_plot",
            f"\n# Distance to '{annot_type}' boundary by {clustering_key} "
            f"({plot_type} plot)\n"
            f"# Measured against annotation shapes drawn in the viewer, which\n"
            f"# this notebook cannot reach. Figure written to:\n"
            + "".join(f"#   {p}\n" for p in ctx.recorded_plot_paths(paths)),
            deps=["preamble"],
            kind=NOTE,
            label="Annotation distance plot",
        )

    def _on_export():
        distances = state.get("annot_dist_distances")
        clustering_key = state.get("annot_dist_clustering_key")
        annot_type = state.get("annot_dist_annot_type")
        if distances is None:
            return

        path, _ = QFileDialog.getSaveFileName(
            None, "Export Distance CSV", "annotation_distances.csv", "CSV Files (*.csv)"
        )
        if not path:
            return

        import pandas as pd
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        cluster_series = ctx.clusterings.get(clustering_key)
        if 'cell_id' in _adata.obs.columns:
            cell_ids = _adata.obs['cell_id'].values
            aligned = cluster_series.reindex(cell_ids).values if cluster_series is not None else [None] * len(distances)
        else:
            aligned = cluster_series.reindex(_adata.obs_names).values if cluster_series is not None else [None] * len(distances)

        centroids_um = np.asarray(_adata.obsm['spatial'], dtype=np.float64)
        df = pd.DataFrame({
            "x_um": centroids_um[:, 0],
            "y_um": centroids_um[:, 1],
            "cluster": [str(c) for c in aligned],
            f"dist_to_{annot_type}_um": distances,
        })
        df.to_csv(path, index=False)
        status.value = f"Exported {len(df):,} rows to {path}"

    def _on_colour_cells():
        distances = state.get("annot_dist_distances")
        annot_type = state.get("annot_dist_annot_type", "annotation")
        if distances is None or ctx.cell_labels_layer is None:
            return

        import matplotlib.pyplot as plt

        # Normalise distances to [0, 1], treating NaN as transparent
        valid_mask = ~np.isnan(distances)
        dist_norm = np.zeros(len(distances), dtype=np.float32)
        if valid_mask.any():
            vmin = distances[valid_mask].min()
            vmax = distances[valid_mask].max()
            if vmax > vmin:
                dist_norm[valid_mask] = (distances[valid_mask] - vmin) / (vmax - vmin)
            else:
                dist_norm[valid_mask] = 1.0

        # Map to RGBA via chosen colormap; NaN cells get alpha=0
        cmap = plt.get_cmap(dist_colormap_widget.value)
        rgba_obs = cmap(dist_norm).astype(np.float32)
        rgba_obs[~valid_mask, 3] = 0.0  # transparent for unmeasured cells

        # Build label-indexed colour array using the same label_to_obs mapping
        cm = ctx.color_manager
        color_arr = cm._empty_color_array()
        lto = cm.label_to_obs
        valid_labels = np.where(lto >= 0)[0]
        color_arr[valid_labels] = rgba_obs[lto[valid_labels]]

        # Save existing colormap so we can restore it on clear
        state["annot_dist_prev_colormap"] = ctx.cell_labels_layer.colormap

        colormap = cm.build_direct_label_colormap(color_arr)
        ctx.cell_labels_layer.colormap = colormap
        clear_color_btn.enabled = True
        status.value = f"Cells coloured by distance to '{annot_type}' ({dist_colormap_widget.value})"

    def _clear_distance_layer():
        prev = state.pop("annot_dist_prev_colormap", None)
        if prev is not None and ctx.cell_labels_layer is not None:
            ctx.cell_labels_layer.colormap = prev
        clear_color_btn.enabled = False

    # ── Connect ───────────────────────────────────────────────────────────────
    refresh_btn.changed.connect(_refresh_annot)
    run_btn.changed.connect(_on_run)
    plot_btn.changed.connect(_on_show_plot)
    export_btn.changed.connect(_on_export)
    color_btn.changed.connect(_on_colour_cells)
    clear_color_btn.changed.connect(_clear_distance_layer)

    # ── Session restore ────────────────────────────────────────────────────────
    def _restore_session(session):
        pass

    # ── Build tab layout ──────────────────────────────────────────────────────
    # Was a hand-rolled QVBoxLayout + QScrollArea doing exactly what make_tab
    # does. That duplication is what kept this tab's labels invisible after the
    # fix landed in make_tab, so it goes through the helper like every other tab.
    scroll = make_tab(
        clustering_widget,
        annot_widget,
        refresh_btn,
        plot_type_widget,
        max_dist_spin,
        run_btn,
        results_text,
        plot_btn,
        export_btn,
        dist_colormap_widget,
        color_btn,
        clear_color_btn,
    )

    return scroll, {"restore_session": _restore_session}
