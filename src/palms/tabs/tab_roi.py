"""Tab 4: ROI Analysis + ROI DEG."""

from __future__ import annotations
from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton
from qtpy.QtWidgets import QTextEdit, QFileDialog, QLabel as QtLabel
from napari.qt.threading import thread_worker
from palms.tabs._helpers import make_tab, StatusProxy, make_progress_bar
from palms.utils.plot_output import batch_dir, plot_formats, save_figure
from palms.utils.prov_graph import ARTIFACT, SETUP, TERMINAL
from palms.utils.steps import Step, StepError, coerce
from palms.utils.step_templates import (
    Preview, builtin_assemble, builtin_spec, builtin_text,
    step_template as _resolved,
)

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


# The drawn polygons are inlined as literals so the notebook reproduces them
# without the viewer's zarr cache. ``rois`` is a SETUP node: it binds a constant.
_ROIS_TEMPLATE = builtin_text("roi.polygons")


# ROI DEG, executed and recorded from one string. The old recorded cell called
# ``palms.utils.gene_analysis.compute_roi_deg`` — so the notebook needed
# this package installed and the reader could not see what it did. It is plain
# spatialdata + scanpy, so it is now written out in full: membership comes from
# ``sd.polygon_query`` against a points element whose µm→pixel scale is a
# declared transformation, not a hand-rolled ``contains_xy`` loop.






# Per-region expression of one gene. This used to record two comment lines
# saying the numbers were "shown in the viewer" — a cell that replays as a
# silent no-op, which ``allow_errors=False`` can never catch. The same
# ``polygon_query`` membership test as the DEG step above, then the statistics
# the tab prints.










ROI_EXPR_TEMPLATE_ID = "roi.expression"
ROI_DEG_TEMPLATE_ID = "roi.deg"
ROI_EXPORT_TEMPLATE_ID = "roi.export_expression"
ROI_POLYGONS_TEMPLATE_ID = "roi.polygons"


def _roi_blocks(filtered: bool) -> list[str]:
    """Both ROI templates share the same two filter injection points.

    The cluster mask is built once before the loop and applied inside it, so a
    filtered run selects two blocks rather than one.
    """
    return (["head"] + (["filter"] if filtered else [])
            + ["loop_head"] + (["loop_filter"] if filtered else []) + ["tail"])


def _roi_expr_template(filtered: bool) -> str:
    return builtin_assemble(ROI_EXPR_TEMPLATE_ID, _roi_blocks(filtered))


def _roi_deg_template(filtered: bool) -> str:
    return builtin_assemble(ROI_DEG_TEMPLATE_ID, _roi_blocks(filtered))


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    roi_calc_button = PushButton(label="Calculate Expression", enabled=True)
    roi_export_button = PushButton(label="Export CSV", enabled=False)
    roi_text = QTextEdit()
    roi_text.setReadOnly(True)
    roi_text.setFontFamily("monospace")
    roi_text.setMaximumHeight(300)

    status_label = StatusProxy(ctx.viewer)
    roi_deg_status = StatusProxy(ctx.viewer)
    roi_deg_progress = make_progress_bar()

    def _rois_preview() -> Preview:
        """The drawn polygons, as the literal the ``rois`` step would bind.

        Shown untruncated in the Templates pane. A shortened polygon list would
        read more easily and would no longer be the string that gets executed,
        which is the one property that pane has.
        """
        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        return Preview(
            list(builtin_spec(ROI_POLYGONS_TEMPLATE_ID).blocks),
            {"polygons": [np.round(np.asarray(p), 2).tolist() for p in polygons]},
        )

    ctx.state.setdefault(
        "template_preview", {})[ROI_POLYGONS_TEMPLATE_ID] = _rois_preview

    def _record_rois():
        """Bind and record ``roi_polygons`` from the drawn shapes."""
        blocks, params, _ = _rois_preview()
        if not params["polygons"]:
            return
        ctx.record_preamble()
        ctx.run_step(Step(
            id="rois",
            **_resolved(ROI_POLYGONS_TEMPLATE_ID, blocks),
            params=params,
            deps=["preamble"],
            kind=SETUP,
            label=f"ROI polygons ({len(params['polygons'])})",
            outputs=["roi_polygons"],
        ))

    def _cluster_filter_params(use_filter: bool) -> dict:
        """The two filter params, or nothing. Shared by the expression and DEG
        previews, which inject the cluster mask at the same two points.

        Read-only: recording the clustering node and mirroring it onto ``obs``
        are the run's job, not a side effect of drawing a preview pane.
        """
        if not use_filter:
            return {}
        clustering_key = ctx.clustering_widget.value
        if not clustering_key:
            return {}
        return {
            "clustering": clustering_key,
            "selected": sorted({str(i) for i in ctx.get_selected_cluster_ids()}),
        }

    def _roi_expr_preview() -> Preview:
        """What "Calculate Expression" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane. "Filter by cluster" selects two blocks
        as well as filling two params, so both halves have to travel together.
        """
        params = {"gene": ctx.gene_widget.value,
                  "pixel_size": coerce(ctx.pixel_size)}
        filter_params = _cluster_filter_params(ctx.filter_check.value)
        params.update(filter_params)
        return Preview(_roi_blocks(bool(filter_params)), params)

    ctx.state.setdefault(
        "template_preview", {})[ROI_EXPR_TEMPLATE_ID] = _roi_expr_preview

    def on_calculate_roi():
        gene = ctx.gene_widget.value
        if gene is None:
            roi_text.setPlainText("No gene selected.")
            return

        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        if len(polygons) == 0:
            roi_text.setPlainText("No ROI polygons drawn.\nUse the Shapes layer to draw polygons.")
            return

        # The polygons must be bound before the step, which reads them.
        _record_rois()

        blocks, params, _ = _roi_expr_preview()
        deps = ["rois"]
        filter_desc = ""
        clustering_key = params.get("clustering")
        if clustering_key is not None:
            ctx.record_clustering(clustering_key)
            from palms.utils.gene_analysis import add_clustering_to_obs
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps.append(f"clustering:{clustering_key}")
            filter_desc = f" ({clustering_key} clusters: {params['selected']})"

        try:
            out = ctx.run_step(Step(
                id=f"roi_expression:{gene}",
                **_resolved(ROI_EXPR_TEMPLATE_ID, blocks),
                params=params,
                deps=deps,
                kind=ARTIFACT,
                label=f"ROI expression: {gene}",
                outputs=["roi_expr_cells", "roi_expr_stats", "roi_expr_tests"],
            ))
        except StepError as e:
            roi_text.setPlainText(str(e))
            status_label.value = f"ROI expression failed: {e}"
            roi_export_button.enabled = False
            return

        cells, stats_df, tests = (out["roi_expr_cells"], out["roi_expr_stats"],
                                  out["roi_expr_tests"])
        state["roi_expr_cells"] = cells
        state["roi_gene"] = gene

        lines = [f"Gene: {gene}{filter_desc}", ""]
        for region_id, row in stats_df.iterrows():
            if row["count"] == 0:
                lines.append(f"Region {region_id}: 0 cells")
                continue
            lines.append(
                f"Region {region_id}: {int(row['count'])} cells, "
                f"mean={row['mean']:.2f}, median={row['median']:.2f}, "
                f"std={row['std']:.2f}, min={row['min']:.0f}, max={row['max']:.0f}"
            )
        if len(tests):
            lines += ["", "── Pairwise Welch's t-tests ──"]
            corrected = len(tests) > 1
            for row in tests.itertuples():
                sig = " *" if (row.p_adj if corrected else row.p) < 0.05 else ""
                adj = f", p_adj(BH)={row.p_adj:.2e}" if corrected else ""
                lines.append(
                    f"  Region {row.region_1} vs {row.region_2}: "
                    f"t={row.t:.3f}, p={row.p:.2e}{adj}{sig}"
                )
            if corrected:
                lines.append(f"  ({len(tests)} comparisons, "
                             f"Benjamini-Hochberg correction)")

        roi_text.setPlainText("\n".join(lines))
        roi_export_button.enabled = len(cells) > 0

    def _roi_export_preview(path: str = None) -> Preview:
        """What "Export CSV" would run: the last calculated gene, and where to.

        ``path`` only exists once the save dialog has returned, so the Templates
        pane is shown the filename that dialog would propose and told, in the
        header, that it is the one value not yet settled.
        """
        gene = state.get("roi_gene", "gene")
        return Preview(
            list(builtin_spec(ROI_EXPORT_TEMPLATE_ID).blocks),
            {"gene": gene, "path": os.fspath(path) if path else f"roi_{gene}.csv"},
            note="" if path else "path chosen on save",
        )

    ctx.state.setdefault(
        "template_preview", {})[ROI_EXPORT_TEMPLATE_ID] = _roi_export_preview

    def on_export_csv():
        cells = state.get("roi_expr_cells")
        gene = state.get("roi_gene", "gene")
        if cells is None or not len(cells):
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ROI Data", f"roi_{gene}.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        # Written *by* the recorded code rather than beside it: the notebook's
        # cell is the statement that produced the file the user has. The full
        # path is recorded, not the basename — a cell that writes somewhere
        # other than where the export went would be a lie about what ran.
        blocks, params, _ = _roi_export_preview(path)
        try:
            ctx.run_step(Step(
                id="export:roi_expression",
                **_resolved(ROI_EXPORT_TEMPLATE_ID, blocks),
                params=params,
                deps=[f"roi_expression:{gene}"],
                kind=TERMINAL,
                label="Export ROI expression",
            ))
        except StepError as e:
            status_label.value = f"Export failed: {e}"
            return
        status_label.value = f"Exported {len(cells)} cells to {path}"

    roi_calc_button.clicked.connect(on_calculate_roi)
    roi_export_button.clicked.connect(on_export_csv)

    # ── ROI DEG widgets ──────────────────────────────────────────────────
    roi_deg_method_widget = ComboBox(
        label="Method", choices=["wilcoxon", "t-test"], value="wilcoxon",
    )
    roi_deg_filter_check = CheckBox(label="Filter by cluster", value=False)
    roi_deg_button = PushButton(label="Run ROI DEG", enabled=True)
    roi_deg_text = QTextEdit()
    roi_deg_text.setReadOnly(True)
    roi_deg_text.setFontFamily("monospace")
    roi_deg_text.setMaximumHeight(250)
    roi_deg_export_button = PushButton(label="Export DEG CSV...", enabled=False)
    roi_volcano_button = PushButton(label="Save Volcano Plot(s)...", enabled=False)

    def _roi_deg_preview() -> Preview:
        """What "Run ROI DEG" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane. Its own filter checkbox, not the ROI
        expression one — the two analyses are filtered independently.
        """
        params = {
            "method": roi_deg_method_widget.value,
            "pixel_size": coerce(ctx.pixel_size),
        }
        filter_params = _cluster_filter_params(roi_deg_filter_check.value)
        params.update(filter_params)
        return Preview(_roi_blocks(bool(filter_params)), params)

    ctx.state.setdefault(
        "template_preview", {})[ROI_DEG_TEMPLATE_ID] = _roi_deg_preview

    def on_roi_deg():
        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        if len(polygons) < 2:
            roi_deg_status.value = "Need at least 2 ROI polygons drawn"
            return

        roi_deg_status.value = "Running differential expression..."
        roi_deg_button.enabled = False
        gen = ctx.dataset_generation

        # The polygons must be bound before the DEG step, which reads them.
        _record_rois()

        blocks, params, _ = _roi_deg_preview()
        deps = ["rois"]
        clustering_key = params.get("clustering")
        if clustering_key is not None:
            ctx.record_clustering(clustering_key)
            from palms.utils.gene_analysis import add_clustering_to_obs
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps.append(f"clustering:{clustering_key}")

        step = Step(
            id="roi_deg",
            **_resolved(ROI_DEG_TEMPLATE_ID, blocks),
            params=params,
            deps=deps,
            kind=ARTIFACT,
            label="ROI differential expression",
            outputs=["roi_deg_df", "roi_adata"],
        )

        @thread_worker
        def _run():
            try:
                out = ctx.run_step(step)
            except StepError as e:
                import pandas as _pd
                return _pd.DataFrame(), None, str(e)
            return out["roi_deg_df"], out["roi_adata"], None

        worker = _run()
        worker.returned.connect(lambda r: _on_roi_deg_ready(r, gen))
        roi_deg_progress.setVisible(True)
        worker.finished.connect(lambda: roi_deg_progress.setVisible(False))
        worker.start()

    def _on_roi_deg_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        df, adata_norm, error = result
        state["roi_deg_df"] = df
        state["roi_deg_adata_norm"] = adata_norm
        roi_deg_button.enabled = True
        # Recording happened inside ctx.run_step(), which recorded the source it ran.
        if error:
            roi_deg_text.setPlainText(error)
            roi_deg_status.value = f"DEG failed: {error}"
            roi_deg_export_button.enabled = False
            roi_volcano_button.enabled = False
            return
        if df.empty:
            roi_deg_text.setPlainText("No significant results or insufficient cells in ROIs.")
            roi_deg_status.value = "DEG: no results"
            roi_deg_export_button.enabled = False
            roi_volcano_button.enabled = False
            return
        preview = df.head(50).to_string(index=False)
        roi_deg_text.setPlainText(preview)
        roi_deg_status.value = f"DEG complete: {len(df)} gene-group results"
        roi_deg_export_button.enabled = True
        roi_volcano_button.enabled = adata_norm is not None
        from palms.utils.adata_persistence import save_roi_deg_to_sdata
        save_roi_deg_to_sdata(ctx, df)

    def on_export_roi_deg():
        df = state.get("roi_deg_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ROI DEG Results", "roi_deg_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        roi_deg_status.value = f"Exported {len(df)} rows to {path}"
        ctx.record_node(
            "export:roi_deg",
            f"\n# Export ROI DEG results\n"
            f"roi_deg_df.to_csv(\"{os.path.basename(path)}\", index=False)",
            deps=["roi_deg"],
            kind=TERMINAL,
            label="Export ROI DEG results",
        )

    def on_roi_generate_volcanos():
        adata_norm = state.get("roi_deg_adata_norm")
        if adata_norm is None:
            roi_deg_status.value = "Run ROI DEG first"
            return
        default_dir = batch_dir(ctx.data_path, "roi_volcano")
        default_dir.mkdir(parents=True, exist_ok=True)
        output_dir = QFileDialog.getExistingDirectory(
            None, "Select output directory for ROI volcano plots",
            str(default_dir))
        if not output_dir:
            return
        roi_volcano_button.enabled = False
        roi_deg_status.value = "Generating ROI volcano plots..."
        method = roi_deg_method_widget.value
        gen = ctx.dataset_generation
        formats = plot_formats(state)
        _recorded_dir = ctx.recorded_plot_paths([output_dir])[0]
        ctx.record_node(
            "plot:roi_volcano",
            f"\n# ROI pairwise volcano plots (method={method})\n"
            f"import itertools\n"
            f"from pathlib import Path\n"
            f"from palms.utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"roi_volcano_dir = Path(r\"{_recorded_dir}\"); "
            f"roi_volcano_dir.mkdir(parents=True, exist_ok=True)\n"
            f"_groups = sorted(roi_adata_norm.obs['roi_region'].cat.categories.tolist())\n"
            f"for _a, _b in itertools.combinations(_groups, 2):\n"
            f"    _df = run_pairwise_deg(roi_adata_norm, 'roi_region', str(_a), str(_b), method=\"{method}\")\n"
            f"    _vfig = make_volcano_plot(_df, str(_a), str(_b), lfc_thresh=1.0, pval_thresh=0.01)\n"
            f"    _stem = 'roi_volcano_' + str(_a).replace(' ', '_') + '_vs_' + str(_b).replace(' ', '_')\n"
            f"    for _ext in {formats!r}:\n"
            f"        _vfig.savefig(roi_volcano_dir / f'{{_stem}}.{{_ext}}', dpi=300)\n"
            f"    plt.close(_vfig)",
            deps=["roi_deg"],
            kind=TERMINAL,
            label="ROI pairwise volcano plots",
        )

        @thread_worker
        def _run():
            from pathlib import Path
            import itertools as _it
            from palms.utils.gene_analysis import run_pairwise_deg, make_volcano_plot
            import matplotlib.pyplot as _plt

            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            groups = sorted(adata_norm.obs['roi_region'].cat.categories.tolist())
            pairs = list(_it.combinations(groups, 2))
            total = len(pairs)
            for i, (a, b) in enumerate(pairs):
                yield f"Volcano plot {i + 1}/{total}: {a} vs {b}"
                df = run_pairwise_deg(adata_norm, 'roi_region', str(a), str(b), method=method)
                fig = make_volcano_plot(df, str(a), str(b), lfc_thresh=1.0, pval_thresh=0.01)
                safe_a = str(a).replace(' ', '_')
                safe_b = str(b).replace(' ', '_')
                save_figure(fig, [out / f'roi_volcano_{safe_a}_vs_{safe_b}.{ext}'
                                  for ext in formats])
                _plt.close(fig)
            return total, output_dir

        worker = _run()
        worker.yielded.connect(lambda msg: setattr(roi_deg_status, 'value', msg) if ctx.dataset_generation == gen else None)
        worker.returned.connect(lambda result: _on_roi_volcanos_done(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def _on_roi_volcanos_done(result):
        count, out_dir = result
        roi_volcano_button.enabled = True
        roi_deg_status.value = f"{count} ROI volcano plot(s) saved to {out_dir}"

    roi_deg_button.clicked.connect(on_roi_deg)
    roi_deg_export_button.clicked.connect(on_export_roi_deg)
    roi_volcano_button.clicked.connect(on_roi_generate_volcanos)

    # ── Layout ───────────────────────────────────────────────────────────
    roi_deg_header = QtLabel("── Differential Expression Between ROIs ──")
    roi_deg_header.setStyleSheet("font-weight: bold; margin-top: 10px;")

    widget = make_tab(
        roi_calc_button,
        roi_text,
        roi_export_button,
        roi_deg_header,
        roi_deg_method_widget,
        roi_deg_filter_check,
        roi_deg_button,
        roi_deg_progress,
        roi_deg_text,
        roi_deg_export_button,
        roi_volcano_button,
    )

    def _restore_session(session):
        rois = session.get("rois", [])
        if rois and ctx.roi_layer is not None:
            ctx.roi_layer.data = rois
            print(f"  Restored {len(rois)} ROI polygons")

        rd = session.get("roi_deg_df")
        if rd is not None:
            state["roi_deg_df"] = rd
            roi_deg_export_button.enabled = True
            print(f"  Restored ROI DEG ({len(rd)} rows)")

    return widget, {"restore_session": _restore_session}
