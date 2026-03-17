"""Tab 6: Gene Analysis — rank genes, dotplot, volcanos."""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext

from utils.gene_analysis import (
    get_normalized_adata, add_clustering_to_obs, run_rank_genes,
    make_rank_genes_dotplot, make_rank_genes_plot, generate_all_volcano_plots,
)


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    ga_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        value=ctx.clustering_names[0] if ctx.clustering_names else None,
    )
    ctx.ga_clustering_widget = ga_clustering_widget

    ga_method_widget = ComboBox(
        label="Method", choices=["wilcoxon", "t-test", "logreg"], value="wilcoxon",
    )
    ga_n_genes_slider = Slider(label="Top N genes", min=5, max=50, value=25)
    ga_run_button = PushButton(label="Run Rank Genes", enabled=True)

    ga_dotplot_n_slider = Slider(label="Genes per cluster", min=3, max=20, value=5)
    ga_dendro_check = CheckBox(label="Dendrogram", value=True)
    ga_dotplot_button = PushButton(label="Show Dotplot", enabled=False)
    ga_edit_labels_button = PushButton(label="Edit Cluster Labels...", enabled=False)
    ga_save_dotplot_button = PushButton(label="Save Dotplot as PNG...", enabled=False)

    ga_rank_plot_button = PushButton(label="Show Rank Genes Plot", enabled=False)

    ga_results_text = QTextEdit()
    ga_results_text.setReadOnly(True)
    ga_results_text.setFontFamily("monospace")
    ga_results_text.setMaximumHeight(300)
    ga_export_button = PushButton(label="Export Full Results CSV...", enabled=False)
    ga_volcano_button = PushButton(label="Generate All Volcano Plots...", enabled=False)
    ga_status = StatusProxy(ctx.viewer)

    def on_run_rank_genes():
        ga_status.value = "Running rank genes (normalizing + computing)..."
        ga_run_button.enabled = False

        clustering_key = ga_clustering_widget.value
        method = ga_method_widget.value
        n_genes = ga_n_genes_slider.value
        state["_rg_method"] = method
        state["_rg_n_genes"] = n_genes

        if not ctx.clusterings or not clustering_key or clustering_key not in ctx.clusterings:
            ga_status.value = "Error: clustering data not available. Please wait for data to finish loading."
            ga_run_button.enabled = True
            return

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        _clustering_series = ctx.clusterings[clustering_key]

        @thread_worker(connect={"returned": _on_rank_genes_ready})
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, _clustering_series, clustering_key)
            df = run_rank_genes(adata_norm, clustering_key, method=method, n_genes=n_genes)
            return df, adata_norm, clustering_key
        _run()

    def _on_rank_genes_ready(result):
        df, adata_norm, clustering_key = result
        state["rank_genes_df"] = df
        state["rank_genes_adata_norm"] = adata_norm
        state["rank_genes_groupby"] = clustering_key
        ga_run_button.enabled = True
        ga_dotplot_button.enabled = True
        ga_rank_plot_button.enabled = True
        ga_edit_labels_button.enabled = True
        ga_export_button.enabled = True
        ga_volcano_button.enabled = True
        preview = df.head(50).to_string(index=False)
        ga_results_text.setPlainText(preview)
        ga_status.value = f"Rank genes done: {len(df)} results ({clustering_key}, {ga_method_widget.value})"

        _rg_method = state.get("_rg_method", "wilcoxon")
        _rg_n = state.get("_rg_n_genes", 25)
        ctx.record_clustering(clustering_key)
        ctx.record_code(
            f"\n# Rank genes: method={_rg_method}, groupby={clustering_key}, n_genes={_rg_n}\n"
            f"sc.tl.rank_genes_groups(adata, groupby=\"{clustering_key}\", "
            f"method=\"{_rg_method}\", n_genes={_rg_n})\n"
            f"rank_df = sc.get.rank_genes_groups_df(adata, group=None)"
        )

    def on_show_dotplot():
        adata_norm = state.get("rank_genes_adata_norm")
        groupby = state.get("rank_genes_groupby")
        if adata_norm is None or groupby is None:
            ga_status.value = "Run rank genes first"
            return
        ga_status.value = "Generating dotplot..."
        ga_dotplot_button.enabled = False

        n_genes = ga_dotplot_n_slider.value
        dendro = ga_dendro_check.value
        labels = ctx.get_labels_for(groupby)

        @thread_worker(connect={"returned": _on_dotplot_ready})
        def _run():
            ctx.apply_plot_font_size()
            fig = make_rank_genes_dotplot(
                adata_norm, groupby, n_genes=n_genes,
                cluster_labels=labels, dendrogram=dendro,
            )
            return fig
        _run()

    def _on_dotplot_ready(fig):
        state["dotplot_fig"] = fig
        ga_dotplot_button.enabled = True
        ga_save_dotplot_button.enabled = True
        import matplotlib.pyplot as _plt
        _plt.show(block=False)
        ga_status.value = "Dotplot displayed"

        _dp_n = ga_dotplot_n_slider.value
        _dp_dendro = ga_dendro_check.value
        _dp_groupby = state.get("rank_genes_groupby", "")
        ctx.record_code(
            f"\n# Dotplot (n_genes={_dp_n}, dendrogram={_dp_dendro})\n"
            + (f"sc.tl.dendrogram(adata, groupby=\"{_dp_groupby}\")\n" if _dp_dendro else "")
            + f"sc.pl.rank_genes_groups_dotplot(adata, n_genes={_dp_n}, "
            f"dendrogram={_dp_dendro})\nplt.show()"
        )

    def _open_label_editor():
        clustering_key = ga_clustering_widget.value
        if ctx.build_label_editor_dialog(clustering_key):
            ga_status.value = f"Labels updated for {clustering_key}"

    def on_save_dotplot():
        fig = state.get("dotplot_fig")
        if fig is None:
            return
        path = ctx.get_plot_save_path("Save Dotplot", "dotplot")
        if not path:
            return
        fig.savefig(path, dpi=300, bbox_inches='tight')
        ga_status.value = f"Dotplot saved to {path}"
        ctx.record_code(f"\n# Save dotplot\n# dotplot -> \"{path}\"")

    def on_show_rank_plot():
        adata_norm = state.get("rank_genes_adata_norm")
        if adata_norm is None:
            ga_status.value = "Run rank genes first"
            return
        import matplotlib.pyplot as _plt
        ctx.apply_plot_font_size()
        _rp_n = ga_n_genes_slider.value
        groupby = state.get("rank_genes_groupby", "")
        labels = ctx.get_labels_for(groupby)
        fig = make_rank_genes_plot(adata_norm, n_genes=_rp_n, cluster_labels=labels)
        _plt.show(block=False)
        ga_status.value = "Rank genes plot displayed"
        ctx.record_code(
            f"\n# Rank genes panel plot (n_genes={_rp_n})\n"
            f"sc.pl.rank_genes_groups(adata, n_genes={_rp_n})\nplt.show()"
        )

    def on_export_rank_genes():
        df = state.get("rank_genes_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Rank Genes Results", "rank_genes_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        ga_status.value = f"Exported {len(df)} rows to {path}"
        ctx.record_code(f"\n# Export rank genes results\n# rank_genes_results.csv -> \"{path}\"")

    def on_generate_volcanos():
        adata_norm = state.get("rank_genes_adata_norm")
        groupby = state.get("rank_genes_groupby")
        if adata_norm is None or groupby is None:
            ga_status.value = "Run rank genes first"
            return
        output_dir = QFileDialog.getExistingDirectory(None, "Select output directory for volcano plots")
        if not output_dir:
            return
        ga_volcano_button.enabled = False
        ga_status.value = "Generating volcano plots..."
        method = state.get("_rg_method", "wilcoxon")
        ctx.record_code(
            f"\n# Generate pairwise volcano plots\n"
            f"from utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"import itertools\n"
            f"volcano_dir = \"{output_dir}\"\n"
            f"# groupby={groupby}, method=\"{method}\""
        )

        @thread_worker(connect={"yielded": lambda msg: setattr(ga_status, 'value', msg),
                                "returned": _on_volcanos_done})
        def _run():
            from pathlib import Path
            import itertools as _it
            from utils.gene_analysis import run_pairwise_deg, make_volcano_plot
            import matplotlib.pyplot as _plt

            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            groups = sorted(
                [g for g in adata_norm.obs[groupby].cat.categories if str(g) != '-1'],
                key=lambda x: (int(x) if str(x).lstrip('-').isdigit() else 0, str(x)),
            )
            pairs = list(_it.combinations(groups, 2))
            total = len(pairs)
            for i, (a, b) in enumerate(pairs):
                yield f"Volcano plot {i + 1}/{total}: {a} vs {b}"
                df = run_pairwise_deg(adata_norm, groupby, str(a), str(b), method=method)
                fig = make_volcano_plot(df, str(a), str(b))
                fig.savefig(out / f'volcano_{a}_vs_{b}.png', dpi=300)
                _plt.close(fig)
            return total, output_dir
        _run()

    def _on_volcanos_done(result):
        count, out_dir = result
        ga_volcano_button.enabled = True
        ga_status.value = f"{count} volcano plots saved to {out_dir}"

    ga_run_button.clicked.connect(on_run_rank_genes)
    ga_dotplot_button.clicked.connect(on_show_dotplot)
    ga_edit_labels_button.clicked.connect(_open_label_editor)
    ga_save_dotplot_button.clicked.connect(on_save_dotplot)
    ga_rank_plot_button.clicked.connect(on_show_rank_plot)
    ga_export_button.clicked.connect(on_export_rank_genes)
    ga_volcano_button.clicked.connect(on_generate_volcanos)

    # Layout
    ga_dotplot_btn_row = QWidget()
    ga_dotplot_btn_layout = QHBoxLayout()
    ga_dotplot_btn_layout.setContentsMargins(0, 0, 0, 0)
    ga_dotplot_btn_layout.addWidget(ga_dotplot_button.native)
    ga_dotplot_btn_layout.addWidget(ga_edit_labels_button.native)
    ga_dotplot_btn_layout.addWidget(ga_save_dotplot_button.native)
    ga_dotplot_btn_row.setLayout(ga_dotplot_btn_layout)

    widget = make_tab(
        ga_clustering_widget,
        ga_method_widget,
        ga_n_genes_slider,
        ga_run_button,
        ga_dotplot_n_slider,
        ga_dendro_check,
        ga_dotplot_btn_row,
        ga_rank_plot_button,
        ga_export_button,
        ga_volcano_button,
    )

    def _restore_session(session):
        rg = session.get("rank_genes_df")
        if rg is not None:
            state["rank_genes_df"] = rg
            ga_export_button.enabled = True
            ga_volcano_button.enabled = True
            print(f"  Restored rank genes ({len(rg)} rows)")

    return widget, {"ga_clustering_widget": ga_clustering_widget, "restore_session": _restore_session}
