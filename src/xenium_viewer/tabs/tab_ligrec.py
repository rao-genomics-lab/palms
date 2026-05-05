"""Tab 7: Ligand-Receptor analysis."""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, Slider
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QFileDialog, QCheckBox, QGridLayout, QGroupBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_tqdm_progress, qt_tqdm_context

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

from xenium_viewer.utils.gene_analysis import get_normalized_adata, add_clustering_to_obs
from xenium_viewer.utils.spatial_analysis import compute_spatial_neighbors, run_ligrec, make_ligrec_plot


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    lr_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        value=ctx.clustering_names[0] if ctx.clustering_names else None,
    )
    ctx.lr_clustering_widget = lr_clustering_widget

    lr_perms_slider = Slider(label="Permutations", min=100, max=1000, value=1000)
    lr_neighs_slider = Slider(label="N neighbors", min=3, max=20, value=6)
    lr_run_button = PushButton(label="Run L-R Analysis", enabled=True)

    # Interaction database filter checkboxes
    lr_ds_group = QGroupBox("Interaction datasets")
    lr_ds_layout = QGridLayout()
    lr_ds_omnipath = QCheckBox("OmniPath"); lr_ds_omnipath.setChecked(True)
    lr_ds_ligrecextra = QCheckBox("LigRecExtra"); lr_ds_ligrecextra.setChecked(True)
    lr_ds_pathwayextra = QCheckBox("PathwayExtra"); lr_ds_pathwayextra.setChecked(True)
    lr_ds_kinaseextra = QCheckBox("KinaseExtra"); lr_ds_kinaseextra.setChecked(True)
    lr_ds_layout.addWidget(lr_ds_omnipath, 0, 0)
    lr_ds_layout.addWidget(lr_ds_ligrecextra, 0, 1)
    lr_ds_layout.addWidget(lr_ds_pathwayextra, 1, 0)
    lr_ds_layout.addWidget(lr_ds_kinaseextra, 1, 1)
    lr_ds_group.setLayout(lr_ds_layout)
    lr_cpdb_only = QCheckBox("CellPhoneDB only")
    lr_cpdb_only.setChecked(False)

    lr_status = StatusProxy(ctx.viewer)

    lr_results_text = QTextEdit()
    lr_results_text.setReadOnly(True)
    lr_results_text.setFontFamily("monospace")
    lr_results_text.setMaximumHeight(250)

    lr_pval_widget = ComboBox(
        label="P-value threshold", choices=["0.001", "0.005", "0.01", "0.05"],
        value="0.05",
    )
    lr_plot_button = PushButton(label="Show L-R Plot", enabled=False)
    lr_export_means_button = PushButton(label="Export Means CSV...", enabled=False)
    lr_export_pvals_button = PushButton(label="Export P-values CSV...", enabled=False)

    def on_run_ligrec():
        lr_status.value = "Running L-R analysis..."
        lr_run_button.enabled = False

        clustering_key = lr_clustering_widget.value
        n_perms = lr_perms_slider.value
        n_neighs = lr_neighs_slider.value
        # Build interactions description for code recording
        ds_names = []
        if lr_ds_omnipath.isChecked(): ds_names.append("OmniPath")
        if lr_ds_ligrecextra.isChecked(): ds_names.append("LigRecExtra")
        if lr_ds_pathwayextra.isChecked(): ds_names.append("PathwayExtra")
        if lr_ds_kinaseextra.isChecked(): ds_names.append("KinaseExtra")
        interactions_desc = ", ".join(ds_names) if ds_names else "none"
        if lr_cpdb_only.isChecked():
            interactions_desc += ", CellPhoneDB_only=True"

        state["_lr_params"] = {
            "clustering_key": clustering_key,
            "n_perms": n_perms,
            "n_neighs": n_neighs,
            "interactions_desc": interactions_desc,
        }
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        from omnipath.constants import InteractionDataset
        include = []
        if lr_ds_omnipath.isChecked():
            include.append(InteractionDataset.OMNIPATH)
        if lr_ds_ligrecextra.isChecked():
            include.append(InteractionDataset.LIGREC_EXTRA)
        if lr_ds_pathwayextra.isChecked():
            include.append(InteractionDataset.PATHWAY_EXTRA)
        if lr_ds_kinaseextra.isChecked():
            include.append(InteractionDataset.KINASE_EXTRA)
        interactions_params = {"include": tuple(include)} if include else {}
        if lr_cpdb_only.isChecked():
            interactions_params["resources"] = "CellPhoneDB"

        _progress = [None]  # filled after worker is created

        @thread_worker
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, ctx.clusterings[clustering_key], clustering_key)
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()
            compute_spatial_neighbors(adata_norm, n_neighs=n_neighs)
            with qt_tqdm_context(_progress[0], "L-R permutations: "):
                result = run_ligrec(adata_norm, clustering_key, n_perms=n_perms,
                                    interactions_params=interactions_params)
            return result

        worker = _run()
        _progress[0], state['_progress_timer'] = attach_tqdm_progress(
            worker,
            lambda m: setattr(lr_status, 'value', m),
            "L-R permutations: ",
        )
        worker.returned.connect(_on_ligrec_ready)
        worker.start()

    def _on_ligrec_ready(result):
        state["ligrec_result"] = result
        lr_run_button.enabled = True

        # Persist to adata.uns immediately
        from xenium_viewer.utils.adata_persistence import save_ligrec_to_adata
        save_ligrec_to_adata(ctx, result)

        warning = result.get('warning')
        means = result['means']
        pvalues = result['pvalues']

        if warning:
            lr_status.value = f"L-R: {warning}"
            lr_plot_button.enabled = False
            lr_export_means_button.enabled = False
            lr_export_pvals_button.enabled = False
            return

        n_interactions = means.shape[0]
        pval_thresh = float(lr_pval_widget.value)
        n_sig = (pvalues < pval_thresh).sum().sum() if not pvalues.empty else 0

        lines = [
            f"L-R interactions found: {n_interactions}",
            f"Significant (p < {pval_thresh}): {n_sig}",
            "",
        ]
        if not means.empty:
            lines.append("Top interactions by mean expression:")
            top_means = means.max(axis=1).sort_values(ascending=False).head(20)
            for idx, val in top_means.items():
                lines.append(f"  {idx}: {val:.4f}")

        lr_results_text.setPlainText("\n".join(lines))
        lr_status.value = f"L-R done: {n_interactions} interactions, {n_sig} significant"
        lr_plot_button.enabled = n_interactions > 0
        lr_export_means_button.enabled = not means.empty
        lr_export_pvals_button.enabled = not pvalues.empty

        _lr_p = state.get("_lr_params", {})
        _lr_ck = _lr_p.get("clustering_key", "")
        _lr_np = _lr_p.get("n_perms", 1000)
        _lr_nn = _lr_p.get("n_neighs", 6)
        ctx.record_clustering(_lr_ck)
        ctx.record_spatial_neighbors(_lr_nn)
        _lr_idesc = _lr_p.get("interactions_desc", "")
        ctx.record_code(
            f"\n# Ligand-receptor analysis (n_perms={_lr_np})\n"
            f"# interactions: {_lr_idesc}\n"
            f"sq.gr.ligrec(\n"
            f"    adata, cluster_key=\"{_lr_ck}\", n_perms={_lr_np},\n"
            f"    threshold=0.01, seed=42,\n"
            f"    transmitter_params={{\"categories\": \"ligand\"}},\n"
            f"    receiver_params={{\"categories\": \"receptor\"}},\n"
            f")"
        )

    def on_show_lr_plot():
        result = state.get("ligrec_result")
        if result is None:
            return
        pval_thresh = float(lr_pval_widget.value)
        groups = ctx.get_cluster_filter()
        lr_ck = state.get("_lr_params", {}).get("clustering_key", "")
        labels = ctx.get_labels_for(lr_ck)
        import matplotlib.pyplot as _plt
        ctx.apply_plot_font_size()
        try:
            fig = make_ligrec_plot(
                result, pvalue_threshold=pval_thresh,
                source_groups=groups, target_groups=groups,
                cluster_labels=labels,
            )
            state["ligrec_fig"] = fig
            _plt.show(block=False)
            path = ctx.auto_save_plot(fig, "ligrec")
            if groups:
                lr_status.value = f"L-R plot displayed (clusters: {', '.join(groups)}) — saved to {path}"
            else:
                lr_status.value = f"L-R plot displayed — saved to {path}"

            _lr_ck = state.get("_lr_params", {}).get("clustering_key", "")
            _lr_fmt = ctx.state.get("plot_format", "svg")
            ctx.record_code(
                f"\n# L-R dotplot (pvalue_threshold={pval_thresh})\n"
                f"sq.pl.ligrec(adata, cluster_key=\"{_lr_ck}\", "
                f"pvalue_threshold={pval_thresh}"
                + (f", source_groups={groups}, target_groups={groups}" if groups else "")
                + f")\nplt.show()\n"
                + f"fig.savefig(\"ligrec.{_lr_fmt}\", dpi=300, bbox_inches='tight')"
            )
        except Exception as e:
            lr_status.value = f"Plot error: {e}"

    def on_export_lr_means():
        result = state.get("ligrec_result")
        if result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export L-R Means", "ligrec_means.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        result['means'].to_csv(path)
        lr_status.value = f"Means exported to {path}"
        ctx.record_code(f"\n# Export L-R means\n# ligrec_means.csv -> \"{path}\"")

    def on_export_lr_pvals():
        result = state.get("ligrec_result")
        if result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export L-R P-values", "ligrec_pvalues.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        result['pvalues'].to_csv(path)
        lr_status.value = f"P-values exported to {path}"
        ctx.record_code(f"\n# Export L-R p-values\n# ligrec_pvalues.csv -> \"{path}\"")

    lr_run_button.clicked.connect(on_run_ligrec)
    lr_plot_button.clicked.connect(on_show_lr_plot)
    lr_export_means_button.clicked.connect(on_export_lr_means)
    lr_export_pvals_button.clicked.connect(on_export_lr_pvals)

    lr_export_btn_row = QWidget()
    lr_export_btn_layout = QHBoxLayout()
    lr_export_btn_layout.setContentsMargins(0, 0, 0, 0)
    lr_export_btn_layout.addWidget(lr_export_means_button.native)
    lr_export_btn_layout.addWidget(lr_export_pvals_button.native)
    lr_export_btn_row.setLayout(lr_export_btn_layout)

    widget = make_tab(
        lr_clustering_widget,
        lr_perms_slider,
        lr_neighs_slider,
        lr_ds_group,
        lr_cpdb_only,
        lr_run_button,
        lr_pval_widget,
        lr_plot_button,
        lr_export_btn_row,
    )

    def _restore_session(session):
        import pandas as pd
        # Check state first (populated from adata.uns at startup), fall back to session
        lr = state.get("ligrec_result")
        if lr is None:
            lm = session.get("ligrec_means")
            lp = session.get("ligrec_pvalues")
            if lm is not None or lp is not None:
                lr = {
                    "means": lm if lm is not None else pd.DataFrame(),
                    "pvalues": lp if lp is not None else pd.DataFrame(),
                    "warning": None,
                }
        if lr is not None:
            state["ligrec_result"] = lr
            lm = lr.get("means")
            lp = lr.get("pvalues")
            lr_export_means_button.enabled = lm is not None and not lm.empty
            lr_export_pvals_button.enabled = lp is not None and not lp.empty
            lr_plot_button.enabled = lm is not None and not lm.empty
            print(f"  Restored L-R results")

    return widget, {"lr_clustering_widget": lr_clustering_widget, "restore_session": _restore_session}
