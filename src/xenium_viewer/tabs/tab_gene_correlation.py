from __future__ import annotations
from typing import TYPE_CHECKING
import os

from magicgui.widgets import ComboBox, PushButton, CheckBox
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.steps import Step

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# The whole figure is now one templated step, so the scatter the viewer shows is
# the scatter the notebook draws. Previously the GUI built the figure by hand and
# a *parallel* code string was recorded that omitted the annotation box, the n in
# the title and the cluster filter entirely.
#
# Expression is pulled with ``sc.get.obs_df`` rather than indexing ``.X`` and
# calling ``.toarray()`` by hand — same values, and it is the idiomatic accessor.
_GENE_CORR_HEAD = """
# Gene correlation ($norm_label): $gene_a vs $gene_b
from scipy.stats import pearsonr, spearmanr"""

_GENE_CORR_EXPR = {
    "Raw counts": """
_expr = sc.get.obs_df(adata, keys=[$gene_a, $gene_b])
x = _expr[$gene_a].to_numpy(dtype='float32')
y = _expr[$gene_b].to_numpy(dtype='float32')""",
    "Fraction of total": """
_expr = sc.get.obs_df(adata, keys=[$gene_a, $gene_b])
_totals = np.asarray(adata.X.sum(axis=1)).ravel().astype('float64')
_totals[_totals == 0] = 1
x = (_expr[$gene_a].to_numpy() / _totals).astype('float32')
y = (_expr[$gene_b].to_numpy() / _totals).astype('float32')""",
    "Log1p(CPM)": """
_expr = sc.get.obs_df(adata_norm, keys=[$gene_a, $gene_b])
x = _expr[$gene_a].to_numpy(dtype='float32')
y = _expr[$gene_b].to_numpy(dtype='float32')""",
}

_GENE_CORR_FILTER = """
_sel = adata.obs[$clustering].astype(str).isin($selected).to_numpy()
x, y = x[_sel], y[_sel]"""

_GENE_CORR_TAIL = """
pr, pp = pearsonr(x, y)
sr, sp = spearmanr(x, y)
_p = lambda p: f'{p:.2e}' if p < 0.001 else f'{p:.4f}'

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(x, y, s=1, alpha=0.3, rasterized=True, color='#1f77b4')
ax.set_xlabel($xlabel)
ax.set_ylabel($ylabel)
ax.set_title($title_prefix + f'  [n={len(x):,}$n_suffix]')
ax.text(
    0.03, 0.97,
    f'Pearson  r = {pr:.3f}, p = {_p(pp)}\\nSpearman \\u03c1 = {sr:.3f}, p = {_p(sp)}',
    transform=ax.transAxes, va='top', ha='left', fontsize=9,
    bbox={'boxstyle': 'round,pad=0.3', 'fc': 'white', 'alpha': 0.7},
)
fig.tight_layout()
fig.savefig($path, dpi=300, bbox_inches='tight')"""

_NORM_LABELS = {
    "Raw counts": "raw counts",
    "Fraction of total": "fraction of total",
    "Log1p(CPM)": "log1p(CPM)",
}


def _gene_corr_template(norm: str, filtered: bool) -> str:
    parts = [_GENE_CORR_HEAD, _GENE_CORR_EXPR[norm]]
    if filtered:
        parts.append(_GENE_CORR_FILTER)
    parts.append(_GENE_CORR_TAIL.replace(
        "$n_suffix", " (filtered)" if filtered else "",
    ))
    return "".join(parts)


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
