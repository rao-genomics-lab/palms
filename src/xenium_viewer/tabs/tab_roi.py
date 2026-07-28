"""Tab 4: ROI Analysis + ROI DEG."""

from __future__ import annotations
from typing import TYPE_CHECKING

import os

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton
from qtpy.QtWidgets import QTextEdit, QFileDialog, QLabel as QtLabel
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, make_progress_bar
from xenium_viewer.utils.prov_graph import ARTIFACT, SETUP, TERMINAL
from xenium_viewer.utils.steps import Step, StepError, coerce

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# The drawn polygons are inlined as literals so the notebook reproduces them
# without the viewer's zarr cache. ``rois`` is a SETUP node: it binds a constant.
_ROIS_TEMPLATE = """
# ROI polygons drawn in the viewer (Nx2 arrays, pixel coords, (y, x) order)
roi_polygons = [np.array(_p) for _p in $polygons]"""


# ROI DEG, executed and recorded from one string. The old recorded cell called
# ``xenium_viewer.utils.gene_analysis.compute_roi_deg`` — so the notebook needed
# this package installed and the reader could not see what it did. It is plain
# shapely + scanpy, so it is now written out in full. (E3 replaces the
# point-in-polygon block with ``spatialdata.polygon_query``.)
_ROI_DEG_HEAD = """
# ROI differential expression (method=$method)
from shapely import contains_xy
from shapely.geometry import Polygon

centroids_yx = adata.obsm['spatial'][:, ::-1] / $pixel_size   # µm→px, xy→yx
roi_region = np.full(adata.n_obs, '', dtype=object)"""

_ROI_DEG_FILTER = """
# Cluster filter: cells must be inside an ROI *and* in the selected clusters
cluster_mask = adata.obs[$clustering].astype(str).isin($selected).to_numpy()"""

_ROI_DEG_LOOP_HEAD = """
for _i, _poly_yx in enumerate(roi_polygons):
    _poly = Polygon(_poly_yx[:, ::-1])
    if not _poly.is_valid:
        _poly = _poly.buffer(0)
    _inside = contains_xy(_poly, centroids_yx[:, 1], centroids_yx[:, 0])"""

_ROI_DEG_LOOP_FILTER = """
    _inside = _inside & cluster_mask"""

_ROI_DEG_TAIL = """
    roi_region[_inside] = f'Region {_i + 1}'

roi_adata = adata[roi_region != ''].copy()
roi_adata.obs['roi_region'] = pd.Categorical(roi_region[roi_region != ''])
sc.pp.normalize_total(roi_adata, target_sum=1e4)
sc.pp.log1p(roi_adata)
sc.tl.rank_genes_groups(
    roi_adata, 'roi_region', method=$method, reference='rest', key_added=$method,
)
roi_deg_df = sc.get.rank_genes_groups_df(roi_adata, group=None, key=$method)"""


def _roi_deg_template(filtered: bool) -> str:
    parts = [_ROI_DEG_HEAD]
    if filtered:
        parts.append(_ROI_DEG_FILTER)
    parts.append(_ROI_DEG_LOOP_HEAD)
    if filtered:
        parts.append(_ROI_DEG_LOOP_FILTER)
    parts.append(_ROI_DEG_TAIL)
    return "".join(parts)


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

    def _record_rois():
        """Bind and record ``roi_polygons`` from the drawn shapes."""
        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        if len(polygons) == 0:
            return
        ctx.record_preamble()
        ctx.run_step(Step(
            id="rois",
            template=_ROIS_TEMPLATE,
            params={"polygons": [
                np.round(np.asarray(p), 2).tolist() for p in polygons
            ]},
            deps=["preamble"],
            kind=SETUP,
            label=f"ROI polygons ({len(polygons)})",
            outputs=["roi_polygons"],
        ))

    def on_calculate_roi():
        from shapely.geometry import Polygon as ShapelyPolygon
        from shapely import contains_xy

        gene = ctx.gene_widget.value
        if gene is None:
            roi_text.setPlainText("No gene selected.")
            return

        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        if len(polygons) == 0:
            roi_text.setPlainText("No ROI polygons drawn.\nUse the Shapes layer to draw polygons.")
            return

        adata = ctx.color_manager.adata
        gene_idx = adata.var_names.get_loc(gene)
        X = adata.X
        if hasattr(X, "toarray"):
            expr = np.asarray(X[:, gene_idx].toarray()).ravel().astype(np.float32)
        else:
            expr = np.asarray(X[:, gene_idx]).ravel().astype(np.float32)

        use_filter = ctx.filter_check.value
        cluster_mask = None
        filter_desc = ""
        if use_filter:
            clustering_key = ctx.clustering_widget.value
            selected_ids = ctx.get_selected_cluster_ids()
            cluster_series = ctx.clusterings[clustering_key]
            if 'cell_id' in adata.obs.columns:
                cell_ids_arr = adata.obs['cell_id'].values
                clusters_aligned = cluster_series.reindex(cell_ids_arr)
            else:
                clusters_aligned = cluster_series.reindex(adata.obs_names)
            cluster_mask = ctx.make_cluster_mask(clusters_aligned.values, selected_ids)
            filter_desc = f" ({clustering_key} clusters: {sorted(selected_ids)})"

        from scipy import stats
        from itertools import combinations

        lines = [f"Gene: {gene}{filter_desc}", ""]
        roi_results = []
        region_exprs = []

        for i, poly_yx in enumerate(polygons):
            poly_xy = poly_yx[:, ::-1]
            shapely_poly = ShapelyPolygon(poly_xy)
            if not shapely_poly.is_valid:
                shapely_poly = shapely_poly.buffer(0)

            inside = contains_xy(shapely_poly, ctx.centroids_yx[:, 1], ctx.centroids_yx[:, 0])
            if cluster_mask is not None:
                inside = inside & cluster_mask
            inside_idx = np.where(inside)[0]
            n_cells = len(inside_idx)

            if n_cells == 0:
                lines.append(f"Region {i+1}: 0 cells")
                region_exprs.append((i + 1, np.array([], dtype=np.float32)))
            else:
                region_expr = expr[inside_idx]
                lines.append(
                    f"Region {i+1}: {n_cells} cells, "
                    f"mean={region_expr.mean():.2f}, "
                    f"median={np.median(region_expr):.2f}, "
                    f"std={region_expr.std():.2f}, "
                    f"min={region_expr.min():.0f}, "
                    f"max={region_expr.max():.0f}"
                )
                region_exprs.append((i + 1, region_expr))

            for idx in inside_idx:
                x_um = ctx.centroids_yx[idx, 1] * ctx.pixel_size
                y_um = ctx.centroids_yx[idx, 0] * ctx.pixel_size
                cell_id = adata.obs['cell_id'].values[idx] if 'cell_id' in adata.obs.columns else str(idx)
                roi_results.append({
                    "region_id": i + 1,
                    "cell_id": cell_id,
                    "x_centroid_um": x_um,
                    "y_centroid_um": y_um,
                    "expression": expr[idx],
                })

        # Significance testing
        testable = [(r, e) for r, e in region_exprs if len(e) >= 2]
        pairs = list(combinations(testable, 2))
        if pairs:
            lines.append("")
            lines.append("── Pairwise Welch's t-tests ──")
            raw_pvals = []
            pair_labels = []
            for (r1, e1), (r2, e2) in pairs:
                t_stat, p_val = stats.ttest_ind(e1, e2, equal_var=False)
                raw_pvals.append(p_val)
                pair_labels.append((r1, r2, t_stat, p_val))

            n_tests = len(raw_pvals)
            if n_tests > 1:
                sorted_idx = np.argsort(raw_pvals)
                adjusted = np.empty(n_tests, dtype=np.float64)
                for rank_pos, orig_idx in enumerate(sorted_idx):
                    adjusted[orig_idx] = raw_pvals[orig_idx] * n_tests / (rank_pos + 1)
                adjusted_sorted = adjusted[sorted_idx]
                for j in range(n_tests - 2, -1, -1):
                    adjusted_sorted[j] = min(adjusted_sorted[j], adjusted_sorted[j + 1])
                adjusted[sorted_idx] = adjusted_sorted
                adjusted = np.minimum(adjusted, 1.0)

                for k, (r1, r2, t_stat, p_raw) in enumerate(pair_labels):
                    p_adj = adjusted[k]
                    sig = " *" if p_adj < 0.05 else ""
                    lines.append(
                        f"  Region {r1} vs {r2}: t={t_stat:.3f}, "
                        f"p={p_raw:.2e}, p_adj(BH)={p_adj:.2e}{sig}"
                    )
                lines.append(f"  ({n_tests} comparisons, Benjamini-Hochberg correction)")
            else:
                r1, r2, t_stat, p_val = pair_labels[0]
                sig = " *" if p_val < 0.05 else ""
                lines.append(
                    f"  Region {r1} vs {r2}: t={t_stat:.3f}, p={p_val:.2e}{sig}"
                )

        roi_text.setPlainText("\n".join(lines))
        state["roi_results"] = roi_results
        state["roi_gene"] = gene
        n_regions = len(polygons)
        _record_rois()
        ctx.record_node(
            f"roi_expression:{gene}",
            f"\n# ROI expression analysis: gene '{gene}' across {n_regions} region(s){filter_desc}\n"
            f"# Per-region mean expression is shown in the viewer; ROI polygons are in roi_polygons.",
            deps=["rois"],
            kind=TERMINAL,
            label=f"ROI expression: {gene}",
        )
        roi_export_button.enabled = len(roi_results) > 0

    def on_export_csv():
        results = state.get("roi_results", [])
        if not results:
            return
        import csv
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ROI Data", f"roi_{state.get('roi_gene', 'gene')}.csv",
            "CSV Files (*.csv)",
        )
        if not path:
            return
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["region_id", "cell_id", "x_centroid_um", "y_centroid_um", "expression"])
            writer.writeheader()
            writer.writerows(results)
        status_label.value = f"Exported {len(results)} cells to {path}"
        _record_rois()
        ctx.record_node(
            "export:roi_expression",
            f"\n# Export ROI per-cell expression\n"
            f"# (region_id, cell_id, centroid, expression for each cell inside an ROI)\n"
            f"# saved from the viewer to {os.path.basename(path)}",
            deps=["rois"],
            kind=TERMINAL,
            label="Export ROI expression",
        )

    roi_calc_button.clicked.connect(on_calculate_roi)
    roi_export_button.clicked.connect(on_export_csv)

    # ── ROI DEG widgets ──────────────────────────────────────────────────
    roi_deg_method_widget = ComboBox(
        label="DEG Method", choices=["wilcoxon", "t-test"], value="wilcoxon",
    )
    roi_deg_filter_check = CheckBox(label="Filter by cluster", value=False)
    roi_deg_button = PushButton(label="Run ROI DEG", enabled=True)
    roi_deg_text = QTextEdit()
    roi_deg_text.setReadOnly(True)
    roi_deg_text.setFontFamily("monospace")
    roi_deg_text.setMaximumHeight(250)
    roi_deg_export_button = PushButton(label="Export DEG CSV...", enabled=False)
    roi_volcano_button = PushButton(label="Save Volcano Plot(s)...", enabled=False)

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

        use_filter = roi_deg_filter_check.value
        deps = ["rois"]
        params = {
            "method": roi_deg_method_widget.value,
            "pixel_size": coerce(ctx.pixel_size),
        }
        clustering_key = None
        if use_filter:
            clustering_key = ctx.clustering_widget.value
            selected_ids = ctx.get_selected_cluster_ids()
            ctx.record_clustering(clustering_key)
            from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            params["clustering"] = clustering_key
            params["selected"] = sorted({str(i) for i in selected_ids})
            deps.append(f"clustering:{clustering_key}")

        step = Step(
            id="roi_deg",
            template=_roi_deg_template(use_filter),
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
        from xenium_viewer.utils.adata_persistence import save_roi_deg_to_sdata
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
        output_dir = QFileDialog.getExistingDirectory(None, "Select output directory for ROI volcano plots")
        if not output_dir:
            return
        roi_volcano_button.enabled = False
        roi_deg_status.value = "Generating ROI volcano plots..."
        method = roi_deg_method_widget.value
        gen = ctx.dataset_generation
        ctx.record_node(
            "plot:roi_volcano",
            f"\n# ROI pairwise volcano plots (method={method})\n"
            f"import itertools\n"
            f"from pathlib import Path\n"
            f"from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"roi_volcano_dir = Path(\"{os.path.basename(output_dir)}\"); "
            f"roi_volcano_dir.mkdir(parents=True, exist_ok=True)\n"
            f"_groups = sorted(roi_adata_norm.obs['roi_region'].cat.categories.tolist())\n"
            f"for _a, _b in itertools.combinations(_groups, 2):\n"
            f"    _df = run_pairwise_deg(roi_adata_norm, 'roi_region', str(_a), str(_b), method=\"{method}\")\n"
            f"    _vfig = make_volcano_plot(_df, str(_a), str(_b), lfc_thresh=1.0, pval_thresh=0.01)\n"
            f"    _name = 'roi_volcano_' + str(_a).replace(' ', '_') + '_vs_' + str(_b).replace(' ', '_') + '.png'\n"
            f"    _vfig.savefig(roi_volcano_dir / _name, dpi=300)\n"
            f"    plt.close(_vfig)",
            deps=["roi_deg"],
            kind=TERMINAL,
            label="ROI pairwise volcano plots",
        )

        @thread_worker
        def _run():
            from pathlib import Path
            import itertools as _it
            from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot
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
                fig.savefig(out / f'roi_volcano_{safe_a}_vs_{safe_b}.png', dpi=300)
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
