from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from magicgui.widgets import ComboBox, PushButton, CheckBox
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    gene_a_widget = ComboBox(label="Gene A", choices=ctx.gene_names,
                             **combo_value_kwargs(ctx.gene_names, index=0))
    gene_b_widget = ComboBox(label="Gene B", choices=ctx.gene_names,
                             **combo_value_kwargs(ctx.gene_names, index=1))
    norm_widget    = ComboBox(
        label="Normalisation",
        choices=["Raw counts", "Fraction of total", "Log1p(CPM)"],
        value="Log1p(CPM)",
    )
    filter_check  = CheckBox(label="Filter by current cluster selection", value=False)
    plot_button   = PushButton(label="Plot Correlation")
    status        = StatusProxy(ctx.viewer)

    def _get_expr(gene_name, norm):
        """Return 1-D float32 expression for gene_name under the requested normalisation."""
        if norm == "Log1p(CPM)":
            from xenium_viewer.utils.gene_analysis import get_normalized_adata
            a = get_normalized_adata(ctx.adata)
        else:
            a = ctx.adata
        gene_idx = a.var_names.get_loc(gene_name)
        X = a.X
        col = X[:, gene_idx]
        if hasattr(col, "toarray"):
            col = col.toarray()
        return np.asarray(col).ravel().astype(np.float32)

    def on_plot():
        from napari.qt.threading import thread_worker
        gene_a = gene_a_widget.value
        gene_b = gene_b_widget.value
        norm   = norm_widget.value
        status.value = f"Computing correlation: {gene_a} vs {gene_b}…"
        gen = ctx.dataset_generation

        # Build optional cell mask
        mask = None
        if filter_check.value:
            cluster_filter = ctx.get_cluster_filter()
            if cluster_filter:
                key = state.get("active_clustering_name") or (
                    ctx.clustering_widget.value if ctx.clustering_widget else None)
                if key:
                    cluster_values, _ = ctx.get_cluster_ids_per_obs(key)
                    selected_ids = ctx.get_selected_cluster_ids()
                    mask = ctx.make_cluster_mask(cluster_values, selected_ids)

        @thread_worker
        def compute():
            from scipy.stats import pearsonr, spearmanr
            x = _get_expr(gene_a, norm)
            y = _get_expr(gene_b, norm)

            if norm == "Fraction of total":
                raw_X = ctx.adata.X
                totals = np.asarray(raw_X.sum(axis=1)).ravel().astype(np.float64)
                totals[totals == 0] = 1
                x = (x / totals).astype(np.float32)
                y = (y / totals).astype(np.float32)

            if mask is not None:
                x, y = x[mask], y[mask]
            pr, pp = pearsonr(x, y)
            sr, sp = spearmanr(x, y)
            return gene_a, gene_b, x, y, pr, pp, sr, sp, norm

        def _on_result(result):
            if ctx.dataset_generation != gen:
                return
            gene_a_r, gene_b_r, x, y, pr, pp, sr, sp, norm_r = result
            import matplotlib.pyplot as plt
            ctx.apply_plot_font_size()

            norm_label = {
                "Raw counts":      "raw counts",
                "Fraction of total": "fraction of total",
                "Log1p(CPM)":      "log1p(CPM)",
            }[norm_r]

            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(x, y, s=1, alpha=0.3, rasterized=True, color="#1f77b4")
            ax.set_xlabel(f"{gene_a_r} [{norm_label}]")
            ax.set_ylabel(f"{gene_b_r} [{norm_label}]")
            n_label = f"n={len(x):,}" + (" (filtered)" if mask is not None else "")
            ax.set_title(f"{gene_a_r} vs {gene_b_r}  [{n_label}]")

            def _fmt_p(p):
                return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"

            ann = (
                f"Pearson  r = {pr:.3f}, p = {_fmt_p(pp)}\n"
                f"Spearman ρ = {sr:.3f}, p = {_fmt_p(sp)}"
            )
            ax.text(0.03, 0.97, ann, transform=ax.transAxes,
                    va="top", ha="left", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

            fig.tight_layout()
            state["corr_fig"] = fig
            state["corr_genes"] = (gene_a_r, gene_b_r)
            path = ctx.auto_save_plot(fig, "gene_correlation")
            status.value = (
                f"{gene_a_r}/{gene_b_r} [{norm_label}] — "
                f"Pearson r={pr:.3f} (p={_fmt_p(pp)}), "
                f"Spearman ρ={sr:.3f} (p={_fmt_p(sp)}) — saved to {path}"
            )

            # Build reproducible code snippet for each normalisation mode
            if norm_r == "Raw counts":
                norm_code = (
                    f"x = adata[:, '{gene_a_r}'].X\n"
                    f"y = adata[:, '{gene_b_r}'].X\n"
                    f"if hasattr(x, 'toarray'): x = x.toarray()\n"
                    f"if hasattr(y, 'toarray'): y = y.toarray()\n"
                    f"x, y = x.ravel().astype('float32'), y.ravel().astype('float32')\n"
                )
            elif norm_r == "Fraction of total":
                norm_code = (
                    f"totals = np.asarray(adata.X.sum(axis=1)).ravel().astype('float64')\n"
                    f"totals[totals == 0] = 1\n"
                    f"x = adata[:, '{gene_a_r}'].X\n"
                    f"y = adata[:, '{gene_b_r}'].X\n"
                    f"if hasattr(x, 'toarray'): x = x.toarray()\n"
                    f"if hasattr(y, 'toarray'): y = y.toarray()\n"
                    f"x = (x.ravel() / totals).astype('float32')\n"
                    f"y = (y.ravel() / totals).astype('float32')\n"
                )
            else:  # Log1p(CPM)
                norm_code = (
                    f"import scanpy as sc\n"
                    f"_norm = adata.copy()\n"
                    f"sc.pp.normalize_total(_norm, target_sum=1e4)\n"
                    f"sc.pp.log1p(_norm)\n"
                    f"x = np.asarray(_norm[:, '{gene_a_r}'].X).ravel().astype('float32')\n"
                    f"y = np.asarray(_norm[:, '{gene_b_r}'].X).ravel().astype('float32')\n"
                )

            _gc_fmt = ctx.state.get("plot_format", "svg")
            ctx.record_node(
                "plot:gene_correlation",
                f"\n# Gene correlation ({norm_r}): {gene_a_r} vs {gene_b_r}\n"
                f"from scipy.stats import pearsonr, spearmanr\n"
                f"import matplotlib.pyplot as plt\n"
                f"import numpy as np\n"
                + norm_code +
                f"pr, pp = pearsonr(x, y)\n"
                f"sr, sp = spearmanr(x, y)\n"
                f"fig, ax = plt.subplots(figsize=(5, 5))\n"
                f"ax.scatter(x, y, s=1, alpha=0.3)\n"
                f"ax.set_xlabel('{gene_a_r} [{norm_label}]')\n"
                f"ax.set_ylabel('{gene_b_r} [{norm_label}]')\n"
                f"print(f'Pearson r={{pr:.3f}}, p={{pp:.3e}}')\n"
                f"print(f'Spearman rho={{sr:.3f}}, p={{sp:.3e}}')\n"
                f"plt.tight_layout(); plt.show()\n"
                f"fig.savefig(\"gene_correlation.{_gc_fmt}\", dpi=300, bbox_inches='tight')",
                deps=["preamble"],
                kind=TERMINAL,
                label="Gene correlation plot",
            )
            fig.show()

        worker = compute()
        worker.returned.connect(_on_result)
        worker.start()

    plot_button.clicked.connect(on_plot)

    widget = make_tab(gene_a_widget, gene_b_widget, norm_widget, filter_check, plot_button)
    return widget, {}
