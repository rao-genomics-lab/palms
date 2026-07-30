from __future__ import annotations
from typing import TYPE_CHECKING
import os

from magicgui.widgets import ComboBox, PushButton, CheckBox
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.steps import Step
from xenium_viewer.utils.step_templates import builtin_assemble

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# The whole figure is now one templated step, so the scatter the viewer shows is
# the scatter the notebook draws. Previously the GUI built the figure by hand and
# a *parallel* code string was recorded that omitted the annotation box, the n in
# the title and the cluster filter entirely.
#
# Expression is pulled with ``sc.get.obs_df`` rather than indexing ``.X`` and
# calling ``.toarray()`` by hand — same values, and it is the idiomatic accessor.




# The "(filtered)" note used to be spliced in by ``str.replace``-ing a fake
# ``$n_suffix`` token out of the tail before Step saw it. It cannot be a real
# param — it sits *inside* an f-string, where ``repr('')`` would render as ``''``
# and break the literal — so the whole title line is a variant instead. That
# also removes a trap: a stray ``$n_suffix`` surviving into a template is a hard
# ``StepError`` from ``Template.substitute``, with a message that names a param
# no call site has ever declared.


_NORM_LABELS = {
    "Raw counts": "raw counts",
    "Fraction of total": "fraction of total",
    "Log1p(CPM)": "log1p(CPM)",
}


TEMPLATE_ID = "genes.correlation"

#: Widget label -> the ``expr.*`` block in the .tmpl that reads expression that
#: way. Three genuinely different expressions, not three parameterisations of
#: one — raw counts and fraction-of-total read ``adata``, log1p(CPM) reads
#: ``adata_norm``, which is why the step's deps differ too.
_EXPR_BLOCK = {
    "Raw counts": "expr.raw",
    "Fraction of total": "expr.fraction",
    "Log1p(CPM)": "expr.log1p_cpm",
}


def _gene_corr_blocks(norm: str, filtered: bool) -> list[str]:
    return (["head", _EXPR_BLOCK[norm]] + (["filter"] if filtered else [])
            + ["stats", "title.filtered" if filtered else "title.plain", "tail"])


def _gene_corr_template(norm: str, filtered: bool) -> str:
    return builtin_assemble(TEMPLATE_ID, _gene_corr_blocks(norm, filtered))


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

    def on_plot():
        from napari.qt.threading import thread_worker
        gene_a = gene_a_widget.value
        gene_b = gene_b_widget.value
        norm   = norm_widget.value
        norm_label = _NORM_LABELS[norm]
        status.value = f"Computing correlation: {gene_a} vs {gene_b}…"
        gen = ctx.dataset_generation

        # The cluster filter reaches the notebook as an explicit
        # ``obs[key].isin([...])`` rather than as an opaque boolean array.
        clustering_key = None
        selected = None
        if filter_check.value and ctx.get_cluster_filter():
            clustering_key = state.get("active_clustering_name") or (
                ctx.clustering_widget.value if ctx.clustering_widget else None)
            if clustering_key:
                ids = ctx.get_selected_cluster_ids()
                selected = sorted({str(i) for i in ids}) if ids else None
        filtered = clustering_key is not None and selected is not None

        fmt = state.get("plot_format", "svg")
        plots_dir = os.path.join(ctx.data_path, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        path = os.path.join(plots_dir, f"gene_correlation.{fmt}")

        params = {
            "gene_a": gene_a, "gene_b": gene_b,
            "norm_label": norm_label,
            "xlabel": f"{gene_a} [{norm_label}]",
            "ylabel": f"{gene_b} [{norm_label}]",
            "title_prefix": f"{gene_a} vs {gene_b}",
            "path": path,
        }
        deps = ["preamble"]
        if norm == "Log1p(CPM)":
            deps = ["normalize"]
        if filtered:
            params["clustering"] = clustering_key
            params["selected"] = selected
            ctx.record_clustering(clustering_key)
            from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps = deps + [f"clustering:{clustering_key}"]

        step = Step(
            id="plot:gene_correlation",
            template=_gene_corr_template(norm, filtered),
            params=params,
            deps=deps,
            kind=TERMINAL,
            label="Gene correlation plot",
            outputs=["fig", "x", "pr", "pp", "sr", "sp"],
        )

        @thread_worker
        def compute():
            if norm == "Log1p(CPM)":
                ctx.ensure_normalized()
            else:
                ctx.record_preamble()
            return ctx.run_step(step)

        def _on_result(out):
            if ctx.dataset_generation != gen:
                return
            import matplotlib.pyplot as plt
            state["corr_fig"] = out["fig"]
            state["corr_genes"] = (gene_a, gene_b)
            plt.show(block=False)
            out["fig"].show()

            def _fmt_p(p):
                return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"
            status.value = (
                f"{gene_a}/{gene_b} [{norm_label}] — "
                f"Pearson r={out['pr']:.3f} (p={_fmt_p(out['pp'])}), "
                f"Spearman ρ={out['sr']:.3f} (p={_fmt_p(out['sp'])}) — saved to {path}"
            )

        def _on_error(e):
            status.value = f"Correlation failed: {e}"

        worker = compute()
        worker.returned.connect(_on_result)
        worker.errored.connect(_on_error)
        worker.start()

    plot_button.clicked.connect(on_plot)

    widget = make_tab(gene_a_widget, gene_b_widget, norm_widget, filter_check, plot_button)
    return widget, {}
