"""Tab 4: ROI Analysis + ROI DEG."""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton
from qtpy.QtWidgets import QTextEdit, QFileDialog, QLabel as QtLabel
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


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
        ctx.record_code(
            f"\n# ROI expression analysis\n"
            f"# gene={gene}, {n_regions} ROI regions{filter_desc}"
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
        ctx.record_code(
            f"\n# Export ROI results\n"
            f"# roi_{state.get('roi_gene', 'gene')}.csv -> \"{path}\""
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

    from utils.gene_analysis import compute_roi_deg

    def on_roi_deg():
        polygons = ctx.roi_layer.data if ctx.roi_layer is not None else []
        if len(polygons) < 2:
            roi_deg_status.value = "Need at least 2 ROI polygons drawn"
            return

        roi_deg_status.value = "Running differential expression..."
        roi_deg_button.enabled = False
        gen = ctx.dataset_generation

        use_filter = roi_deg_filter_check.value
        cluster_mask = None
        if use_filter:
            clustering_key = ctx.clustering_widget.value
            selected_ids = ctx.get_selected_cluster_ids()
            cluster_series = ctx.clusterings[clustering_key]
            _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
            if 'cell_id' in _adata.obs.columns:
                clusters_aligned = cluster_series.reindex(_adata.obs['cell_id'].values)
            else:
                clusters_aligned = cluster_series.reindex(_adata.obs_names)
            cluster_mask = ctx.make_cluster_mask(clusters_aligned.values, selected_ids)

        method = roi_deg_method_widget.value
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        @thread_worker
        def _run():
            return compute_roi_deg(
                _adata, ctx.centroids_yx, polygons, ctx.pixel_size,
                cluster_mask=cluster_mask, method=method,
            )

        _use_filter = use_filter
        _selected_ids = selected_ids if use_filter else None

        worker = _run()
        worker.returned.connect(lambda df: _on_roi_deg_ready(df, gen, _use_filter, _selected_ids))
        worker.start()

    def _on_roi_deg_ready(result, _gen, _use_filter, _selected_ids):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        df, adata_norm = result
        state["roi_deg_df"] = df
        state["roi_deg_adata_norm"] = adata_norm
        roi_deg_button.enabled = True
        if not df.empty:
            filter_line = ""
            if _use_filter:
                clustering_key = ctx.clustering_widget.value
                filter_line = (
                    f"\n# cluster filter (clustering='{clustering_key}', "
                    f"clusters={sorted(_selected_ids)})\n"
                    f"# cells must be inside an ROI AND in the selected clusters\n"
                    f"cluster_series = adata.obs[{clustering_key!r}]\n"
                    f"clusters_aligned = cluster_series.reindex(adata.obs['cell_id'].values)\n"
                    f"cluster_mask = np.isin(clusters_aligned.values, {sorted(_selected_ids)})\n"
                )
            ctx.record_code(
                f"\n# ROI differential expression\n"
                f"import json, numpy as np\n"
                f"from utils.adata_persistence import load_rois_from_sdata\n"
                f"from utils.gene_analysis import compute_roi_deg\n"
                f"pixel_size = float(json.load(open(data_path / 'experiment.xenium'))['pixel_size'])\n"
                f"centroids_yx = adata.obsm['spatial'][:, ::-1] / pixel_size  # µm→px, xy→yx\n"
                f"roi_polygons = load_rois_from_sdata(sdata)  # Nx2 yx arrays\n"
                f"{filter_line}"
                f"roi_deg_df, roi_adata_norm = compute_roi_deg(\n"
                f"    adata, centroids_yx, roi_polygons, pixel_size,\n"
                f"    method={roi_deg_method_widget.value!r},\n"
                f"    cluster_mask={'cluster_mask' if _use_filter else 'None'},\n"
                f")"
            )
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
        from utils.adata_persistence import save_roi_deg_to_sdata
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
        ctx.record_code(
            f"\n# Export ROI DEG results\n"
            f"# roi_deg_results.csv -> \"{path}\""
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
        ctx.record_code(
            f"\n# Generate ROI pairwise volcano plots\n"
            f"from utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"import itertools\n"
            f"roi_volcano_dir = \"{output_dir}\"\n"
            f"# Uses roi_adata_norm from DEG step, method=\"{method}\""
        )

        @thread_worker
        def _run():
            from pathlib import Path
            import itertools as _it
            from utils.gene_analysis import run_pairwise_deg, make_volcano_plot
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
