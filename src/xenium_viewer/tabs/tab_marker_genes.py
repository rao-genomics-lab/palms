"""Tab: Marker Genes — dotplot/heatmap/matrixplot/tracksplot/correlation
from a user-supplied marker gene dict."""

from __future__ import annotations
from typing import TYPE_CHECKING

import json

from magicgui.widgets import ComboBox, PushButton
from qtpy.QtWidgets import QTextEdit, QFileDialog, QLabel as QtLabel
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.steps import Step, StepError
from xenium_viewer.utils.step_templates import builtin_assemble

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

_PLACEHOLDER = '''\
{
  "Cell type A": ["Gene1", "Gene2", "Gene3"],
  "Cell type B": ["Gene4", "Gene5"]
}'''


# This tab recorded nothing at all before, despite being five plain scanpy
# plotting calls. It is now templated like the rest: the marker dict and the
# display labels reach the notebook as literals, so a replay reproduces the
# figure rather than a differently-grouped approximation of it.


# sc.pl.correlation_matrix is the one that needs its statistic computed first,
# and the one that ignores var_names.

# Raster formats get an explicit dpi; SVG does not, where it means nothing.
# This used to splice ``, dpi=150`` in by ``str.replace``-ing a fake
# ``$dpi_kwarg`` token out of the tail before Step saw it — a token that is not
# a param, and whose escape into a template would raise ``StepError`` naming a
# param no call site declares. Two whole-line variants instead.


TEMPLATE_ID = "genes.marker_plot"

#: Plot types the tab offers; each is one ``call.*`` block in the .tmpl file.
MARKER_PLOTS = ("dotplot", "heatmap", "matrixplot", "tracksplot",
                "correlation_matrix")


def _marker_plot_blocks(plot_name: str, relabel: bool, dpi: bool) -> list[str]:
    return (["head"] + (["relabel"] if relabel else [])
            + [f"call.{plot_name}"] + ["save.dpi" if dpi else "save.plain"])


def _marker_plot_template(plot_name: str, relabel: bool, dpi: bool) -> str:
    return builtin_assemble(
        TEMPLATE_ID, _marker_plot_blocks(plot_name, relabel, dpi))


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Clustering selector (registered on ctx for refresh_clustering_choices) ──
    mg_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.mg_clustering_widget = mg_clustering_widget

    # ── Marker dict input ────────────────────────────────────────────────
    input_label = QtLabel("Marker genes dict (JSON):")
    marker_text = QTextEdit()
    marker_text.setPlaceholderText(_PLACEHOLDER)
    marker_text.setFontFamily("monospace")
    marker_text.setMinimumHeight(180)

    # ── Format selector ──────────────────────────────────────────────────
    fmt_widget = ComboBox(label="Save format", choices=["PNG", "SVG"], value="PNG")

    # ── Action buttons ───────────────────────────────────────────────────
    dot_button    = PushButton(label="Dotplot")
    heat_button   = PushButton(label="Heatmap")
    matrix_button = PushButton(label="Matrix plot")
    tracks_button = PushButton(label="Tracks plot")
    corr_button   = PushButton(label="Correlation matrix")

    status_label = StatusProxy(ctx.viewer)

    # ─────────────────────────────────────────────────────────────────────
    def _parse_marker_dict():
        """Parse and validate the JSON marker dict from the text widget.
        Returns the dict, or None on error (sets status message).
        """
        text = marker_text.toPlainText().strip()
        if not text:
            status_label.value = "Enter a marker genes dict first"
            return None
        try:
            marker_dict = json.loads(text)
        except json.JSONDecodeError as e:
            status_label.value = f"JSON parse error: {e}"
            return None
        if not isinstance(marker_dict, dict):
            status_label.value = "Input must be a JSON object (dict)"
            return None

        # Flatten to check gene names
        all_genes = [g for genes in marker_dict.values() for g in genes]
        var_set = set(ctx.adata.var_names)
        unknown = [g for g in all_genes if g not in var_set]
        if unknown:
            status_label.value = (
                f"Unknown genes (ignored): {', '.join(unknown[:5])}"
                + (" ..." if len(unknown) > 5 else "")
            )
            # Filter out unknown genes but continue
            marker_dict = {
                k: [g for g in v if g in var_set]
                for k, v in marker_dict.items()
                if any(g in var_set for g in v)
            }
            if not marker_dict:
                status_label.value = "No valid genes found in marker dict"
                return None

        # Persist the (possibly cleaned) dict as JSON
        state["marker_genes_json"] = json.dumps(marker_dict, indent=2)
        return marker_dict

    def _pick_save_path(plot_name: str) -> str | None:
        fmt = fmt_widget.value.lower()
        path, _ = QFileDialog.getSaveFileName(
            None, f"Save {plot_name}",
            f"{plot_name}.{fmt}",
            f"{fmt_widget.value} Files (*.{fmt})",
        )
        return path or None

    # ── Generic runner ───────────────────────────────────────────────────
    def _run_plot(plot_name: str):
        """Validate inputs, ask for save path, then run the plot as a Step."""
        marker_dict = _parse_marker_dict()
        if marker_dict is None:
            return
        path = _pick_save_path(plot_name)
        if not path:
            return

        clustering_key = mg_clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            status_label.value = "Select a valid clustering first"
            return

        for btn in (dot_button, heat_button, matrix_button, tracks_button, corr_button):
            btn.enabled = False
        status_label.value = f"Generating {plot_name}..."
        gen = ctx.dataset_generation

        _adata  = ctx.adata
        _labels = ctx.get_labels_for(clustering_key)
        _fmt    = fmt_widget.value.lower()

        # The clustering must exist as a node, and in adata.obs, before a step
        # can declare it as a dependency and read it. Both on the GUI thread.
        from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata,
                              ctx.clusterings[clustering_key], clustering_key)

        # Display names, resolved here so they render as a literal list.
        categories = None
        if _labels:
            categories = []
            for c in _adata.obs[clustering_key].cat.categories:
                try:
                    label = _labels.get(int(c), _labels.get(c, c))
                except (ValueError, TypeError):
                    label = _labels.get(c, c)
                categories.append(str(label))

        params = {
            "plot_name": plot_name,
            "groupby": clustering_key,
            "markers": marker_dict,
            "path": path,
        }
        if categories is not None:
            params["categories"] = categories

        step = Step(
            id=f"plot:markers:{plot_name}:{clustering_key}",
            template=_marker_plot_template(
                plot_name, relabel=categories is not None, dpi=_fmt == "png",
            ),
            params=params,
            deps=["normalize", f"clustering:{clustering_key}"],
            kind=TERMINAL,
            label=f"Marker {plot_name}: {clustering_key}",
        )

        @thread_worker
        def _run():
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            ctx.ensure_normalized()
            plt.close('all')
            try:
                ctx.run_step(step)
            except StepError as e:
                return str(e)
            plt.close('all')
            return path

        worker = _run()
        worker.returned.connect(
            lambda p: _on_done(plot_name, p) if ctx.dataset_generation == gen else None
        )
        worker.start()

    def _on_done(plot_name: str, path: str):
        for btn in (dot_button, heat_button, matrix_button, tracks_button, corr_button):
            btn.enabled = True
        status_label.value = f"{plot_name} saved: {path}"

    dot_button.clicked.connect(lambda: _run_plot("dotplot"))
    heat_button.clicked.connect(lambda: _run_plot("heatmap"))
    matrix_button.clicked.connect(lambda: _run_plot("matrixplot"))
    tracks_button.clicked.connect(lambda: _run_plot("tracksplot"))
    corr_button.clicked.connect(lambda: _run_plot("correlation_matrix"))

    # ── Layout ───────────────────────────────────────────────────────────
    widget = make_tab(
        mg_clustering_widget,
        input_label,
        marker_text,
        fmt_widget,
        dot_button,
        heat_button,
        matrix_button,
        tracks_button,
        corr_button,
    )

    # ── Session persistence ───────────────────────────────────────────────
    def _restore_session(session):
        json_str = session.get("marker_genes_json")
        if json_str:
            state["marker_genes_json"] = json_str
            marker_text.setPlainText(json_str)
            print(f"  Restored marker genes dict")

    return widget, {"restore_session": _restore_session}
