"""Tab: Annotation Neighbourhood Enrichment.

Runs squidpy neighbourhood enrichment on real Xenium cells combined with
virtual "cells" sampled from user-drawn annotation polygons.  Annotation
types appear as additional rows/columns in the Z-score heatmap alongside
real cell-type clusters.

Every button here runs a ``Step``, so the recorded code is the executed code.
This tab recorded nothing but a ``NOTE`` until 2026-09-01, on the stated
grounds that the notebook cannot reach a napari shapes layer. It does not have
to: ``annot.polygons`` inlines the drawn shapes as literals, exactly as
``roi.polygons`` has done for ROIs since the ROI migration.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, PushButton, Slider, SpinBox
from qtpy.QtWidgets import (
    QVBoxLayout, QLabel, QTextEdit, QFileDialog,
    QCheckBox, QGroupBox,
)
from napari.qt.threading import thread_worker
from palms.tabs._helpers import (
    StatusProxy, annotation_polygons_preview, attach_tqdm_progress,
    qt_tqdm_context, make_progress_bar, combo_value_kwargs, make_tab,
)
from palms.utils.plot_output import safe_stem
from palms.utils.prov_graph import ARTIFACT, TERMINAL
from palms.utils.steps import Step, StepError, coerce
from palms.utils.step_templates import (
    Preview, builtin_spec, step_template as _resolved,
)

#: This tab runs four steps: the drawn shapes, the virtual cells sampled inside
#: them, the enrichment over both, and the heatmap.
POLYGONS_TEMPLATE_ID = "annot.polygons"
VIRTUAL_CELLS_TEMPLATE_ID = "annot.virtual_cells"
TEMPLATE_ID = "annot.nhood"
PLOT_TEMPLATE_ID = "annot.nhood_plot"
EXPORT_TEMPLATE_ID = "annot.export_nhood"

# The filename the save dialog proposes, and what the Templates pane shows as
# the sample path before a dialog has returned one.
_EXPORT_FILENAME = "annot_nhood_zscore.csv"

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Clustering selector ───────────────────────────────────────────────────
    clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.annot_nhood_clustering_widget = clustering_widget

    # ── Annotation type checkboxes ────────────────────────────────────────────
    annot_group_label = QLabel("Include annotation types as virtual cells:")
    annot_group = QGroupBox()
    annot_group_layout = QVBoxLayout()
    annot_group_layout.setContentsMargins(4, 4, 4, 4)
    annot_group.setLayout(annot_group_layout)

    _annot_checkboxes: list[QCheckBox] = []

    def _rebuild_annot_checkboxes():
        for cb in _annot_checkboxes:
            annot_group_layout.removeWidget(cb)
            cb.deleteLater()
        _annot_checkboxes.clear()

        from palms.utils.annotation_utils import get_annotation_types
        types = get_annotation_types(ctx.annotation_layer)
        if not types:
            placeholder = QLabel("  (no annotation types defined — use the Annotations tab)")
            annot_group_layout.addWidget(placeholder)
            _annot_checkboxes.append(placeholder)  # type: ignore[arg-type]
            return
        for t in types:
            cb = QCheckBox(t)
            cb.setChecked(True)
            annot_group_layout.addWidget(cb)
            _annot_checkboxes.append(cb)

    _rebuild_annot_checkboxes()

    refresh_annot_btn = PushButton(label="Refresh annotation types")

    # ── Parameters ────────────────────────────────────────────────────────────
    density_spin = SpinBox(label="Grid density (µm²/virtual cell)", min=10, max=10000, value=100)
    perms_slider = Slider(label="Permutations", min=100, max=1000, value=1000)
    neighs_slider = Slider(label="Neighbours", min=3, max=20, value=6)

    # ── Controls ──────────────────────────────────────────────────────────────
    run_btn = PushButton(label="Run Annotation Nhood Enrichment", enabled=True)
    mode_widget = ComboBox(label="Display mode", choices=["zscore", "count"], value="zscore")
    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(250)
    plot_btn = PushButton(label="Show Heatmap", enabled=False)
    export_btn = PushButton(label="Export Z-scores CSV...", enabled=False)
    status = StatusProxy(ctx.viewer)
    annot_nhood_progress = make_progress_bar()

    from palms.utils.gene_analysis import add_clustering_to_obs

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def _refresh_annot():
        _rebuild_annot_checkboxes()

    def _get_selected_annot_types() -> list[str]:
        return [
            cb.text() for cb in _annot_checkboxes
            if isinstance(cb, QCheckBox) and cb.isChecked()
        ]

    # ── Preview providers ─────────────────────────────────────────────────────
    # One expression per step of what this tab would run as the widgets stand,
    # called by the run below *and* by the Templates tab's preview pane.

    def _annot_polygons_preview() -> Preview:
        """The drawn shapes. Delegates, because the distance tab needs the same
        expression and two copies of "what has been drawn" would drift."""
        return annotation_polygons_preview(ctx)

    def _virtual_cells_preview() -> Preview:
        return Preview(
            list(builtin_spec(VIRTUAL_CELLS_TEMPLATE_ID).blocks),
            {
                "types": _get_selected_annot_types(),
                "density_um2": float(coerce(density_spin.value)),
            },
        )

    def _nhood_preview() -> Preview:
        clustering_key = clustering_widget.value
        return Preview(
            list(builtin_spec(TEMPLATE_ID).blocks),
            {
                "cluster_key": clustering_key,
                "uns_key": f"{clustering_key}_nhood_enrichment",
                "n_neighs": coerce(neighs_slider.value),
                "n_perms": coerce(perms_slider.value),
                "seed": 42,
            },
        )

    def _nhood_plot_preview() -> Preview:
        result = state.get("annot_nhood_result") or {}
        clustering_key = result.get("_cluster_key") or clustering_widget.value
        annot_types = result.get("_annot_types") or _get_selected_annot_types()
        return Preview(
            list(builtin_spec(PLOT_TEMPLATE_ID).blocks),
            {
                "cluster_key": clustering_key,
                "mode": mode_widget.value,
                "title": ("Annotation Nhood Enrichment\n"
                          f"(annotation types: {', '.join(annot_types)})"),
                "paths": ctx.plot_paths(_plot_stem(clustering_key)),
            },
        )

    def _nhood_export_preview(path: str = None) -> Preview:
        """What "Export Z-scores CSV" would run, and where to.

        ``path`` only exists once the save dialog has returned, so the Templates
        pane is shown the filename that dialog would propose and told, in the
        header, that it is the one value not yet settled.
        """
        result = state.get("annot_nhood_result") or {}
        clustering_key = result.get("_cluster_key") or clustering_widget.value
        return Preview(
            list(builtin_spec(EXPORT_TEMPLATE_ID).blocks),
            {
                "cluster_key": clustering_key,
                "uns_key": f"{clustering_key}_nhood_enrichment",
                "path": os.fspath(path) if path else _EXPORT_FILENAME,
            },
            note="" if path else "path chosen on save",
        )

    def _plot_stem(clustering_key) -> str:
        return f"annot_nhood_enrichment_{safe_stem(clustering_key)}"

    ctx.state.setdefault("template_preview", {})[POLYGONS_TEMPLATE_ID] = _annot_polygons_preview
    ctx.state.setdefault("template_preview", {})[VIRTUAL_CELLS_TEMPLATE_ID] = _virtual_cells_preview
    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _nhood_preview
    ctx.state.setdefault("template_preview", {})[PLOT_TEMPLATE_ID] = _nhood_plot_preview
    ctx.state.setdefault("template_preview", {})[EXPORT_TEMPLATE_ID] = _nhood_export_preview

    def _on_run():
        clustering_key = clustering_widget.value
        if clustering_key is None or clustering_key not in ctx.clusterings:
            results_text.setPlainText("No clustering selected.")
            return

        annot_types = _get_selected_annot_types()
        if not annot_types:
            results_text.setPlainText("No annotation types selected.")
            return

        # The shapes, first: the two analysis steps declare this node as a
        # dependency, and a run with nothing drawn has no region to sample.
        if ctx.ensure_annotations(_annot_polygons_preview()) is None:
            results_text.setPlainText(
                "No typed annotations have been drawn — use the Annotations tab."
            )
            return

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        # The clustering must exist as a node, and in adata.obs, before a step
        # can declare it as a dependency and read it. Both on the GUI thread.
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata, ctx.clusterings[clustering_key],
                              clustering_key)

        virtual_blocks, virtual_params, _ = _virtual_cells_preview()
        nhood_blocks, nhood_params, _ = _nhood_preview()

        virtual_step = Step(
            id="annot_virtual_cells",
            **_resolved(VIRTUAL_CELLS_TEMPLATE_ID, virtual_blocks),
            params=virtual_params,
            deps=["annotations"],
            kind=ARTIFACT,
            label=f"Annotation virtual cells ({', '.join(annot_types)})",
            outputs=["annot_virtual_cells"],
        )
        nhood_step = Step(
            id=f"annot_nhood:{clustering_key}",
            **_resolved(TEMPLATE_ID, nhood_blocks),
            params=nhood_params,
            deps=["normalize", f"clustering:{clustering_key}",
                  "annot_virtual_cells"],
            kind=ARTIFACT,
            label=f"Annotation nhood enrichment: {clustering_key}",
            outputs=["adata_annot"],
        )

        run_btn.enabled = False
        status.value = "Building augmented dataset with annotation virtual cells..."
        _progress = [None]

        @thread_worker
        def _run():
            ctx.ensure_normalized()
            try:
                virtual = ctx.run_step(virtual_step)["annot_virtual_cells"]
                if len(virtual) == 0:
                    return {"warning": "No virtual cells fell inside the selected "
                                       "annotations — try a finer grid density.",
                            "zscore": np.array([]), "count": np.array([]),
                            "clusters": []}
                with qt_tqdm_context(_progress[0], "Enrichment permutations: "):
                    adata_annot = ctx.run_step(nhood_step)["adata_annot"]
            except StepError as e:
                return {"warning": str(e), "zscore": np.array([]),
                        "count": np.array([]), "clusters": []}

            uns = adata_annot.uns[nhood_params["uns_key"]]
            return {
                "zscore": np.array(uns["zscore"]),
                "count": np.array(uns["count"]),
                "clusters": list(
                    adata_annot.obs[clustering_key].cat.categories.astype(str)),
                "warning": None,
                "_adata_annot": adata_annot,
                "_cluster_key": clustering_key,
                "_annot_types": annot_types,
                "_n_virtual": len(virtual),
            }

        worker = _run()
        _progress[0], state['_annot_nhood_progress_timer'] = attach_tqdm_progress(
            worker,
            lambda m: setattr(status, 'value', m),
            "Enrichment permutations: ",
            progress_bar=annot_nhood_progress,
        )
        worker.returned.connect(_on_done)
        worker.start()

    def _on_done(result):
        state["annot_nhood_result"] = result
        run_btn.enabled = True

        warning = result.get('warning')
        if warning:
            status.value = f"Annot nhood: {warning}"
            results_text.setPlainText(warning)
            plot_btn.enabled = False
            export_btn.enabled = False
            return

        zscore = result['zscore']
        clusters = result['clusters']
        n = len(clusters)
        annot_types = result.get('_annot_types', [])

        lines = [
            f"Annotation neighbourhood enrichment: {n}x{n} matrix",
            f"Clusters: {', '.join(clusters)}",
            f"Annotation types included: {', '.join(annot_types)}",
            "",
            "Note: virtual cells have zero gene expression; Z-scores reflect spatial",
            "proximity only (not gene-expression-derived similarity).",
            "",
        ]

        if zscore.size > 0:
            # Show top enrichments involving annotation types
            ann_set = set(annot_types)
            pairs = []
            for i, c1 in enumerate(clusters):
                for j, c2 in enumerate(clusters):
                    if i != j and (c1 in ann_set or c2 in ann_set):
                        pairs.append((c1, c2, zscore[i, j]))
            pairs.sort(key=lambda x: x[2], reverse=True)
            if pairs:
                lines.append("Top enrichments involving annotation types:")
                for c1, c2, z in pairs[:10]:
                    lines.append(f"  {c1} <-> {c2}: {z:.2f}")

        results_text.setPlainText("\n".join(lines))
        status.value = f"Annotation nhood enrichment done: {n} groups"
        plot_btn.enabled = True
        export_btn.enabled = True

    def _on_show_plot():
        result = state.get("annot_nhood_result")
        if result is None:
            return
        clustering_key = result.get("_cluster_key", "")
        ctx.apply_plot_font_size()

        blocks, params, _ = _nhood_plot_preview()
        step = Step(
            id=f"plot:annot_nhood:{clustering_key}",
            **_resolved(PLOT_TEMPLATE_ID, blocks),
            params=params,
            deps=[f"annot_nhood:{clustering_key}"],
            kind=TERMINAL,
            label=f"Annotation nhood heatmap: {clustering_key}",
            outputs=["fig"],
        )
        # A TERMINAL now, not a NOTE. The NOTE was right while the figure came
        # from a viewer-only computation the notebook could not express — a
        # terminal that saved a figure nothing had drawn would have replayed as
        # an empty file. The enrichment is a recorded step now, so the heatmap
        # is drawable from it, and the cell that draws it is the cell that ran.
        try:
            fig = ctx.run_step(step)["fig"]
        except StepError as e:
            status.value = f"Plot error: {e}"
            return
        paths = ctx.show_plot(
            fig, _plot_stem(clustering_key),
            title=f"Annotation nhood enrichment: {clustering_key}",
            save=False, paths=params["paths"])
        status.value = f"Heatmap displayed — saved to {', '.join(paths)}"

    def _on_export():
        result = state.get("annot_nhood_result")
        if result is None:
            return
        zscore = result.get('zscore')
        clustering_key = result.get('_cluster_key')
        if zscore is None or zscore.size == 0 or clustering_key is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Z-scores CSV", _EXPORT_FILENAME, "CSV Files (*.csv)"
        )
        if not path:
            return
        # Written *by* the recorded code rather than beside it: the notebook's
        # cell is the statement that produced the file the user has. The full
        # path is recorded, not the basename — a cell that writes somewhere
        # other than where the export went would be a lie about what ran.
        blocks, params, _ = _nhood_export_preview(path)
        try:
            ctx.run_step(Step(
                id=f"export:annot_nhood:{clustering_key}",
                **_resolved(EXPORT_TEMPLATE_ID, blocks),
                params=params,
                deps=[f"annot_nhood:{clustering_key}"],
                kind=TERMINAL,
                label=f"Export annotation nhood z-scores: {clustering_key}",
            ))
        except StepError as e:
            status.value = f"Export failed: {e}"
            return
        status.value = f"Exported Z-scores to {path}"

    # ── Connect ───────────────────────────────────────────────────────────────
    refresh_annot_btn.changed.connect(_refresh_annot)
    run_btn.changed.connect(_on_run)
    plot_btn.changed.connect(_on_show_plot)
    export_btn.changed.connect(_on_export)

    # ── Session restore ────────────────────────────────────────────────────────
    def _restore_session(session):
        pass  # No persistent state for this tab

    # ── Build tab layout ──────────────────────────────────────────────────────
    # Was a hand-rolled QVBoxLayout + QScrollArea doing exactly what make_tab
    # does. That duplication is what kept this tab's labels invisible after the
    # fix landed in make_tab, so it goes through the helper like every other tab.
    scroll = make_tab(
        clustering_widget,
        annot_group_label,
        annot_group,
        refresh_annot_btn,
        density_spin,
        perms_slider,
        neighs_slider,
        run_btn,
        annot_nhood_progress,
        mode_widget,
        results_text,
        plot_btn,
        export_btn,
    )

    return scroll, {"restore_session": _restore_session}
