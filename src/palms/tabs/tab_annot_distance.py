"""Tab: Distance to Annotation.

For each Xenium cell, computes the minimum Euclidean distance (µm) to the
boundary of a selected annotation type, then visualises the distribution per
cell-type cluster as violin / box plots.  Optionally colours cells on the
napari canvas by their distance value.

The distances land in ``adata.obs`` under ``dist_to_<type>_um``, which is what
makes them a result the notebook can go on to use rather than a number the
viewer kept to itself. The napari colouring below reads the same column.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, PushButton, SpinBox
from palms.utils.coloring import AVAILABLE_COLORMAPS
from qtpy.QtWidgets import (
    QTextEdit, QFileDialog,
)
from napari.qt.threading import thread_worker
from palms.tabs._helpers import (
    StatusProxy, annotation_polygons_preview, combo_value_kwargs, make_tab,
)
from palms.utils.plot_output import safe_stem
from palms.utils.prov_graph import ARTIFACT, TERMINAL
from palms.utils.steps import Step, StepError, coerce
from palms.utils.step_templates import (
    Preview, builtin_spec, step_template as _resolved,
)

# The filename the save dialog proposes, and what the Templates pane shows as
# the sample path before a dialog has returned one.
_EXPORT_FILENAME = "annotation_distances.csv"

TEMPLATE_ID = "annot.distance"
PLOT_TEMPLATE_ID = "annot.distance_plot"
EXPORT_TEMPLATE_ID = "annot.export_distance"

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Selectors ─────────────────────────────────────────────────────────────
    clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.annot_dist_clustering_widget = clustering_widget

    def _annot_choices():
        from palms.utils.annotation_utils import get_annotation_types
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

    def _obs_key(annot_type) -> str:
        """Where the distances land in ``adata.obs``. One column per annotation
        type, so measuring against a second type does not erase the first."""
        return f"dist_to_{safe_stem(annot_type)}_um"

    def _plot_stem(annot_type, clustering_key) -> str:
        return f"annot_distance_{safe_stem(annot_type)}_{safe_stem(clustering_key)}"

    # ── Preview providers ─────────────────────────────────────────────────────

    def _distance_preview() -> Preview:
        annot_type = annot_widget.value
        if annot_type in (None, "(none)"):
            annot_type = ""
        return Preview(
            list(builtin_spec(TEMPLATE_ID).blocks),
            {"annotation_type": annot_type, "obs_key": _obs_key(annot_type)},
        )

    def _distance_plot_blocks(relabel: bool, clip: bool) -> list[str]:
        """The plot type *is* a block, and so are the two optional steps before
        it — which is why block selection lives here and not in the registry."""
        blocks = ["head"]
        if relabel:
            blocks.append("relabel")
        if clip:
            blocks.append("clip")
        return blocks + [f"plot.{plot_type_widget.value}", "save"]

    def _distance_plot_preview() -> Preview:
        annot_type = state.get("annot_dist_annot_type") or annot_widget.value
        if annot_type in (None, "(none)"):
            annot_type = ""
        clustering_key = (state.get("annot_dist_clustering_key")
                          or clustering_widget.value or "")
        categories = ctx.get_labels_for(clustering_key) or None
        max_dist = max_dist_spin.value

        params = {
            "obs_key": _obs_key(annot_type),
            "cluster_key": clustering_key,
            "annotation_type": annot_type,
            "paths": ctx.plot_paths(_plot_stem(annot_type, clustering_key)),
        }
        if categories is not None:
            params["categories"] = {str(k): str(v) for k, v in categories.items()}
        if max_dist > 0:
            params["max_dist"] = float(coerce(max_dist))
        return Preview(
            _distance_plot_blocks(categories is not None, max_dist > 0), params)

    def _distance_export_preview(path: str = None) -> Preview:
        """What "Export CSV" would run: the last measured type, and where to.

        ``path`` only exists once the save dialog has returned, so the Templates
        pane is shown the filename that dialog would propose and told, in the
        header, that it is the one value not yet settled.
        """
        annot_type = state.get("annot_dist_annot_type") or annot_widget.value
        if annot_type in (None, "(none)"):
            annot_type = ""
        cluster_key = (state.get("annot_dist_clustering_key")
                       or clustering_widget.value or "")
        return Preview(
            list(builtin_spec(EXPORT_TEMPLATE_ID).blocks),
            {"annotation_type": annot_type, "obs_key": _obs_key(annot_type),
             "cluster_key": cluster_key,
             "path": os.fspath(path) if path else _EXPORT_FILENAME},
            note="" if path else "path chosen on save",
        )

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _distance_preview
    ctx.state.setdefault("template_preview", {})[PLOT_TEMPLATE_ID] = _distance_plot_preview
    ctx.state.setdefault("template_preview", {})[EXPORT_TEMPLATE_ID] = _distance_export_preview

    def _on_run():
        clustering_key = clustering_widget.value
        annot_type = annot_widget.value
        if not clustering_key or annot_type in (None, "(none)"):
            results_text.setPlainText("Select a clustering and annotation type.")
            return
        if ctx.ensure_annotations(annotation_polygons_preview(ctx)) is None:
            results_text.setPlainText(
                "No typed annotations have been drawn — use the Annotations tab."
            )
            return

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        # The clustering is a dependency of the *plot*, which groups by it, so
        # it has to be a node and an obs column before either step runs.
        from palms.utils.gene_analysis import add_clustering_to_obs
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata, ctx.clusterings[clustering_key],
                              clustering_key)

        blocks, params, _ = _distance_preview()
        step = Step(
            id=f"annot_distance:{annot_type}",
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=[ctx.cell_root(), "annotations"],
            kind=ARTIFACT,
            label=f"Distance to '{annot_type}'",
            outputs=["annot_distances"],
        )

        run_btn.enabled = False
        status.value = f"Computing distances to '{annot_type}'..."

        @thread_worker
        def _run():
            return ctx.run_step(step)["annot_distances"]

        def _failed(exc):
            run_btn.enabled = True
            results_text.setPlainText(f"Distance analysis failed: {exc}")
            status.value = f"Distance analysis failed: {exc}"

        def _on_done(distances):
            run_btn.enabled = True
            distances = np.asarray(distances, dtype=np.float64)
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
        worker.errored.connect(_failed)
        worker.start()

    def _on_show_plot():
        distances = state.get("annot_dist_distances")
        clustering_key = state.get("annot_dist_clustering_key")
        annot_type = state.get("annot_dist_annot_type")
        if distances is None or clustering_key is None:
            return

        ctx.apply_plot_font_size()
        blocks, params, _ = _distance_plot_preview()
        step = Step(
            id=f"plot:annot_distance:{annot_type}:{clustering_key}",
            **_resolved(PLOT_TEMPLATE_ID, blocks),
            params=params,
            deps=[f"annot_distance:{annot_type}",
                  f"clustering:{clustering_key}"],
            kind=TERMINAL,
            label=f"Distance to '{annot_type}': {clustering_key}",
            outputs=["fig"],
        )
        # A TERMINAL now, not a NOTE — the distances are a recorded step, so the
        # figure is drawable from them. Seaborn rather than sc.pl.violin because
        # the tab offers box and strip too, and the recorded cell has to be the
        # cell that ran; the DataFrame it plots comes out of sc.get.obs_df.
        try:
            fig = ctx.run_step(step)["fig"]
        except StepError as e:
            status.value = f"Plot error: {e}"
            return

        paths = ctx.show_plot(
            fig, _plot_stem(annot_type, clustering_key),
            title=f"Distance to '{annot_type}' by {clustering_key}",
            save=False, paths=params["paths"])
        status.value = f"Distance plot displayed — saved to {', '.join(paths)}"

    def _on_export():
        distances = state.get("annot_dist_distances")
        clustering_key = state.get("annot_dist_clustering_key")
        annot_type = state.get("annot_dist_annot_type")
        if distances is None or clustering_key is None or annot_type is None:
            return

        path, _ = QFileDialog.getSaveFileName(
            None, "Export Distance CSV", _EXPORT_FILENAME, "CSV Files (*.csv)"
        )
        if not path:
            return

        # Written *by* the recorded code rather than beside it: the notebook's
        # cell is the statement that produced the file the user has. The full
        # path is recorded, not the basename — a cell that writes somewhere
        # other than where the export went would be a lie about what ran.
        blocks, params, _ = _distance_export_preview(path)
        try:
            ctx.run_step(Step(
                id=f"export:annot_distance:{annot_type}:{clustering_key}",
                **_resolved(EXPORT_TEMPLATE_ID, blocks),
                params=params,
                deps=[f"annot_distance:{annot_type}",
                      f"clustering:{clustering_key}"],
                kind=TERMINAL,
                label=f"Export distance to '{annot_type}'",
            ))
        except StepError as e:
            status.value = f"Export failed: {e}"
            return
        status.value = f"Exported {len(distances):,} rows to {path}"

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
