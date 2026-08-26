"""Tab: Marker Genes — dotplot/heatmap/matrixplot/tracksplot/correlation
from a user-supplied marker gene dict."""

from __future__ import annotations
from typing import TYPE_CHECKING

import json

from magicgui.widgets import ComboBox, PushButton
from qtpy.QtWidgets import QTextEdit, QLabel as QtLabel
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, combo_value_kwargs
from xenium_viewer.utils.plot_output import safe_stem
from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.steps import Step
from xenium_viewer.utils.step_templates import (
    Preview, builtin_assemble, step_template as _resolved,
)

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

# The save block used to come in two variants — one with ``dpi=150``, one
# without — chosen by the tab's PNG/SVG combo. There is no such choice any more:
# a plot is written in every format Preferences asks for, so the block loops over
# a ``paths`` list and passes dpi unconditionally (PDF and SVG ignore it). That
# is also what lets the notebook name the files the GUI actually wrote.
#
# The dpi variants themselves replaced an earlier trick — splicing ``, dpi=150``
# in by ``str.replace``-ing a fake ``$dpi_kwarg`` token out of the tail before
# Step saw it, hiding the template from ``Template.substitute``'s check. Neither
# idiom should come back; ``tests/test_template_placeholders.py`` guards it.


TEMPLATE_ID = "genes.marker_plot"

#: Plot types the tab offers; each is one ``call.*`` block in the .tmpl file.
MARKER_PLOTS = ("dotplot", "heatmap", "matrixplot", "tracksplot",
                "correlation_matrix")


def _marker_plot_blocks(plot_name: str, relabel: bool) -> list[str]:
    return (["head"] + (["relabel"] if relabel else [])
            + [f"call.{plot_name}", "save"])


def _marker_plot_template(plot_name: str, relabel: bool) -> str:
    """The *shipped* text for these options. Tests pin this.

    ``builtin_assemble`` is right here and wrong at the call site: a test that
    asserts what the shipped template says must not read the developer's own
    overrides. The run below goes through ``step_template`` instead — reading
    builtin text there meant a user could edit this template, validate it, save
    it, and watch nothing change.
    """
    return builtin_assemble(
        TEMPLATE_ID, _marker_plot_blocks(plot_name, relabel))


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

    # No format selector: Preferences → Plot format governs every figure the
    # viewer writes, so a per-tab combo could only disagree with it.

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

    def _stem(plot_name: str, clustering_key) -> str:
        """Keyed by clustering, so a second run does not overwrite the first."""
        return f"{plot_name}_{safe_stem(clustering_key or 'none')}"

    def _display_categories(clustering_key: str) -> dict | None:
        """Original cluster id -> display name, or None if unnamed.

        Reads ``adata.obs`` when the clustering has been mirrored there and the
        stored series otherwise, so the preview can answer before a run has put
        the column in place.
        """
        labels = ctx.get_labels_for(clustering_key)
        if not labels:
            return None
        obs = ctx.adata.obs if ctx.adata is not None else None
        if obs is not None and clustering_key in obs:
            column = obs[clustering_key]
            source = (column.cat.categories if hasattr(column, "cat")
                      else sorted(column.unique()))
        else:
            series = ctx.clusterings.get(clustering_key)
            if series is None:
                return None
            source = sorted(series.unique())
        # A mapping, not a list: the template merges clusters that share a
        # display name, and it needs the original key to map from. A list also
        # silently mis-aligns if the category order ever shifts.
        display = {}
        for c in source:
            try:
                label = labels.get(int(c), labels.get(c, c))
            except (ValueError, TypeError):
                label = labels.get(c, c)
            display[str(c)] = str(label)
        return display

    def _marker_preview(plot_name: str = None, marker_dict: dict = None) -> Preview:
        """What one of the five plot buttons would run, as the widgets stand.

        One expression of the current settings, called by ``_run_plot`` with the
        values it has validated, and by the Templates tab's preview pane with
        none — five buttons share this callback, so the pane shows the last one
        run (a dotplot until something else has been). The blocks belong with the
        params: the plot type *is* a block, so a params-only preview would show
        the wrong call.

        The destination is no longer a dialog answer: every figure goes to
        ``<dataset>/plots/`` in the configured formats, so the preview can show
        the real paths rather than a filename a dialog might propose.
        """
        plot_name = plot_name or state.get("_mg_last_plot") or MARKER_PLOTS[0]
        if marker_dict is None:
            # Best effort and side-effect free: the pane must render whatever is
            # in the box, without the status messages and state writes that
            # _parse_marker_dict owes the user when they press a button.
            try:
                marker_dict = json.loads(marker_text.toPlainText().strip())
            except (json.JSONDecodeError, TypeError):
                marker_dict = None
            if not isinstance(marker_dict, dict):
                marker_dict = json.loads(_PLACEHOLDER)

        clustering_key = mg_clustering_widget.value
        categories = _display_categories(clustering_key) if clustering_key else None

        params = {
            "plot_name": plot_name,
            "groupby": clustering_key,
            "markers": marker_dict,
            "paths": ctx.plot_paths(_stem(plot_name, clustering_key)),
        }
        if categories is not None:
            params["categories"] = categories
        return Preview(
            _marker_plot_blocks(plot_name, relabel=categories is not None),
            params,
        )

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _marker_preview

    # ── Generic runner ───────────────────────────────────────────────────
    def _run_plot(plot_name: str):
        """Validate inputs, then run the plot as a Step."""
        marker_dict = _parse_marker_dict()
        if marker_dict is None:
            return

        clustering_key = mg_clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            status_label.value = "Select a valid clustering first"
            return

        for btn in (dot_button, heat_button, matrix_button, tracks_button, corr_button):
            btn.enabled = False
        status_label.value = f"Generating {plot_name}..."
        gen = ctx.dataset_generation
        state["_mg_last_plot"] = plot_name

        _adata = ctx.adata

        # The clustering must exist as a node, and in adata.obs, before a step
        # can declare it as a dependency and read it. Both on the GUI thread.
        # Also before the params are built: the display categories come off the
        # obs column this puts there.
        from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
        ctx.record_clustering(clustering_key)
        add_clustering_to_obs(_adata, _adata,
                              ctx.clusterings[clustering_key], clustering_key)

        ctx.apply_plot_font_size()
        blocks, params, _ = _marker_preview(plot_name, marker_dict)

        step = Step(
            id=f"plot:markers:{plot_name}:{clustering_key}",
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=["normalize", f"clustering:{clustering_key}"],
            kind=TERMINAL,
            label=f"Marker {plot_name}: {clustering_key}",
            outputs=["fig"],
        )

        # No ``matplotlib.use('Agg')`` here. It used to be set inside this
        # worker and never restored, which silently disabled every later
        # ``plt.show`` in the session — figures from other tabs simply stopped
        # appearing. Display now goes through the Plots dock's own canvas, so
        # the process-wide backend is nobody's business.
        @thread_worker
        def _run():
            ctx.ensure_normalized()
            return ctx.run_step(step)

        def _done(out):
            if ctx.dataset_generation != gen:
                return
            _enable()
            paths = ctx.show_plot(
                out["fig"], _stem(plot_name, clustering_key),
                title=f"{plot_name}: {clustering_key}",
                save=False, paths=params["paths"])
            status_label.value = f"{plot_name} saved: {', '.join(paths)}"

        def _failed(exc):
            _enable()
            status_label.value = f"{plot_name} failed: {exc}"

        worker = _run()
        worker.returned.connect(_done)
        worker.errored.connect(_failed)
        worker.start()

    def _enable():
        for btn in (dot_button, heat_button, matrix_button, tracks_button, corr_button):
            btn.enabled = True

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
