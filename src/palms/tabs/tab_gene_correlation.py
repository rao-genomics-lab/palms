from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, CheckBox
from palms.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from palms.utils.plot_output import safe_stem
from palms.utils.prov_graph import TERMINAL
from palms.utils.steps import Step
from palms.utils.step_templates import (
    Preview, builtin_assemble, step_template as _resolved,
)

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


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
    # Published on ctx so refresh_gene_choices can re-populate it: the
    # choices are bound to var_names once, here, and a gene filter or a
    # segmentation swap can shrink that panel underneath them.
    ctx.corr_gene_a_widget = gene_a_widget
    ctx.corr_gene_b_widget = gene_b_widget
    norm_widget    = ComboBox(
        label="Normalisation",
        choices=["Raw counts", "Fraction of total", "Log1p(CPM)"],
        value="Log1p(CPM)",
    )
    filter_check  = CheckBox(label="Filter by current cluster selection", value=False)
    plot_button   = PushButton(label="Plot Correlation")
    status        = StatusProxy(ctx.viewer)

    def _stem(gene_a: str, gene_b: str) -> str:
        """Keyed by the gene pair — a second correlation used to overwrite the
        first, because every run wrote ``plots/gene_correlation.<fmt>``."""
        return f"gene_correlation_{safe_stem(gene_a)}_{safe_stem(gene_b)}"

    def _gene_corr_preview() -> Preview:
        """What "Plot Correlation" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane. Both halves matter here: the
        normalisation combo selects which ``expr.*`` block reads the expression
        — three genuinely different expressions, not one parameterised — so a
        params-only preview would show the wrong statement entirely.

        Read-only, deliberately. The run creates ``plots/`` before writing into
        it; opening the Templates tab must not create a directory as a side
        effect of drawing a pane.
        """
        gene_a = gene_a_widget.value
        gene_b = gene_b_widget.value
        norm_label = _NORM_LABELS[norm_widget.value]

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

        params = {
            "gene_a": gene_a, "gene_b": gene_b,
            "norm_label": norm_label,
            "xlabel": f"{gene_a} [{norm_label}]",
            "ylabel": f"{gene_b} [{norm_label}]",
            "title_prefix": f"{gene_a} vs {gene_b}",
            "paths": ctx.plot_paths(_stem(gene_a, gene_b)),
        }
        if filtered:
            params["clustering"] = clustering_key
            params["selected"] = selected
        return Preview(_gene_corr_blocks(norm_widget.value, filtered), params)

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _gene_corr_preview

    def on_plot():
        from napari.qt.threading import thread_worker
        norm = norm_widget.value
        blocks, params, _ = _gene_corr_preview()
        gene_a, gene_b = params["gene_a"], params["gene_b"]
        norm_label = params["norm_label"]
        clustering_key = params.get("clustering")
        status.value = f"Computing correlation: {gene_a} vs {gene_b}…"
        gen = ctx.dataset_generation
        ctx.apply_plot_font_size()

        deps = ["preamble"]
        if norm == "Log1p(CPM)":
            deps = ["normalize"]
        if clustering_key is not None:
            ctx.record_clustering(clustering_key)
            from palms.utils.gene_analysis import add_clustering_to_obs
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps = deps + [f"clustering:{clustering_key}"]

        step = Step(
            id="plot:gene_correlation",
            **_resolved(TEMPLATE_ID, blocks),
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
            state["corr_fig"] = out["fig"]
            state["corr_genes"] = (gene_a, gene_b)
            # The template wrote the files; show_plot only publishes them.
            paths = ctx.show_plot(out["fig"], _stem(gene_a, gene_b),
                                  title=f"{gene_a} vs {gene_b} [{norm_label}]",
                                  save=False, paths=params["paths"])

            def _fmt_p(p):
                return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"
            status.value = (
                f"{gene_a}/{gene_b} [{norm_label}] — "
                f"Pearson r={out['pr']:.3f} (p={_fmt_p(out['pp'])}), "
                f"Spearman ρ={out['sr']:.3f} (p={_fmt_p(out['sp'])}) — "
                f"saved to {', '.join(paths)}"
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
