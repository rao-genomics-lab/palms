"""Tab: Marker Genes — dotplot/heatmap/matrixplot/tracksplot/correlation
from a user-supplied marker gene dict."""

from __future__ import annotations
from typing import TYPE_CHECKING

import json

from magicgui.widgets import ComboBox, PushButton
from qtpy.QtWidgets import QTextEdit, QFileDialog, QLabel as QtLabel
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext

_PLACEHOLDER = '''\
{
  "Cell type A": ["Gene1", "Gene2", "Gene3"],
  "Cell type B": ["Gene4", "Gene5"]
}'''


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Clustering selector (registered on ctx for refresh_clustering_choices) ──
    mg_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names or [""],
        value=ctx.clustering_names[0] if ctx.clustering_names else "",
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

    def _get_adata_norm():
        """Return (adata_norm, groupby_key) with clustering added to obs."""
        from utils.gene_analysis import get_normalized_adata, add_clustering_to_obs

        clustering_key = mg_clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            return None, None

        _adata = ctx.adata
        adata_norm = get_normalized_adata(_adata)
        add_clustering_to_obs(adata_norm, ctx.clusterings[clustering_key],
                              _adata, key_name=clustering_key)

        # Apply cluster labels as category names
        labels = ctx.get_labels_for(clustering_key)
        if labels:
            cats = adata_norm.obs[clustering_key].cat.categories
            new_cats = []
            for c in cats:
                try:
                    label = labels.get(int(c), labels.get(c, c))
                except (ValueError, TypeError):
                    label = labels.get(c, c)
                new_cats.append(str(label))
            adata_norm.obs[clustering_key] = (
                adata_norm.obs[clustering_key].cat.rename_categories(new_cats)
            )
        return adata_norm, clustering_key

    def _pick_save_path(plot_name: str) -> str | None:
        fmt = fmt_widget.value.lower()
        path, _ = QFileDialog.getSaveFileName(
            None, f"Save {plot_name}",
            f"{plot_name}.{fmt}",
            f"{fmt_widget.value} Files (*.{fmt})",
        )
        return path or None

    # ── Generic runner ───────────────────────────────────────────────────
    def _run_plot(plot_name: str, plot_fn):
        """Validate inputs, ask for save path, run plot_fn in thread."""
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

        _adata    = ctx.adata
        _ckey     = clustering_key
        _cser     = ctx.clusterings[clustering_key].copy()
        _labels   = ctx.get_labels_for(clustering_key)
        _mdict    = marker_dict
        _fmt      = fmt_widget.value.lower()

        @thread_worker
        def _run():
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from utils.gene_analysis import get_normalized_adata, add_clustering_to_obs

            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(adata_norm, _adata, _cser, key_name=_ckey)

            if _labels:
                cats = adata_norm.obs[_ckey].cat.categories
                new_cats = []
                for c in cats:
                    try:
                        label = _labels.get(int(c), _labels.get(c, c))
                    except (ValueError, TypeError):
                        label = _labels.get(c, c)
                    new_cats.append(str(label))
                adata_norm.obs[_ckey] = (
                    adata_norm.obs[_ckey].cat.rename_categories(new_cats)
                )

            # plot_fn renders to the current figure; capture with gcf()
            plt.close('all')
            plot_fn(adata_norm, _mdict, _ckey)
            fig = plt.gcf()
            save_kwargs = dict(bbox_inches='tight')
            if _fmt == 'png':
                save_kwargs['dpi'] = 150
            fig.savefig(path, **save_kwargs)
            plt.close(fig)
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

    # ── Plot functions (each runs inside a thread via _run_plot) ─────────
    # These render to plt's current figure; _run_plot captures with plt.gcf().
    def _dotplot_fn(adata_norm, marker_dict, groupby):
        import scanpy as sc
        sc.pl.dotplot(adata_norm, var_names=marker_dict, groupby=groupby, show=False)

    def _heatmap_fn(adata_norm, marker_dict, groupby):
        import scanpy as sc
        sc.pl.heatmap(adata_norm, var_names=marker_dict, groupby=groupby, show=False)

    def _matrixplot_fn(adata_norm, marker_dict, groupby):
        import scanpy as sc
        sc.pl.matrixplot(adata_norm, var_names=marker_dict, groupby=groupby, show=False)

    def _tracksplot_fn(adata_norm, marker_dict, groupby):
        import scanpy as sc
        sc.pl.tracksplot(adata_norm, var_names=marker_dict, groupby=groupby, show=False)

    def _corrplot_fn(adata_norm, marker_dict, groupby):
        import scanpy as sc
        sc.tl.correlation_matrix(adata_norm, groupby)
        sc.pl.correlation_matrix(adata_norm, groupby, show=False)

    dot_button.clicked.connect(
        lambda: _run_plot("dotplot", _dotplot_fn))
    heat_button.clicked.connect(
        lambda: _run_plot("heatmap", _heatmap_fn))
    matrix_button.clicked.connect(
        lambda: _run_plot("matrixplot", _matrixplot_fn))
    tracks_button.clicked.connect(
        lambda: _run_plot("tracksplot", _tracksplot_fn))
    corr_button.clicked.connect(
        lambda: _run_plot("correlation_matrix", _corrplot_fn))

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
