"""Tab 3: UMAP — the linked scatter window, plus publication UMAP figures.

Two figures, one template. **By gene** (issue #34) draws one panel per selected
gene, each with its own colour scale — the thing the napari UMAP window cannot
do, since a Points layer has no colour bar and no way to show several genes at
once. **By cluster** replaces the old "Save UMAP Plot…", which built a throwaway
AnnData by hand, wrote it through a file dialog under ``matplotlib.use('Agg')``
(leaking that backend for the rest of the session), and recorded a cell that
*recomputed* the embedding the viewer had never used.

Both go through ``ctx.show_plot``, so they appear in the Plots dock and land in
``<dataset>/plots/`` like every other figure.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, Slider, SpinBox
from qtpy.QtWidgets import QLabel, QListWidget
from napari.qt.threading import thread_worker

from xenium_viewer.tabs._helpers import make_tab, StatusProxy
from xenium_viewer.utils.coloring import AVAILABLE_COLORMAPS
from xenium_viewer.utils.plot_output import safe_stem
from xenium_viewer.utils.prov_graph import NOTE, TERMINAL
from xenium_viewer.utils.step_templates import Preview, step_template as _resolved
from xenium_viewer.utils.steps import Step, coerce

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


TEMPLATE_ID = "umap.plot"

#: The issue asks for "up to 10? 15?" panels. Fifteen: beyond that the panels are
#: too small to read at any figure size scanpy will produce.
MAX_GENES = 15

#: Where ``loader.load_umap`` reads Xenium's own embedding from. The template
#: reads the same file, so the notebook draws the coordinates the viewer drew.
UMAP_PROJECTION = ("analysis", "umap", "gene_expression_2_components",
                   "projection.csv")


def has_xenium_umap(data_path) -> bool:
    """Whether this dataset ships Xenium's UMAP, or the step must recompute one.

    A Crop Dataset export has no ``analysis/`` folder — the case
    ``loader.load_umap`` already returns ``None`` for — and the viewer then
    falls back to whatever is in ``obsm['X_umap']``. The notebook has no such
    fallback, so it computes its own embedding and says so.
    """
    if data_path is None:
        return False
    from pathlib import Path
    return Path(data_path).joinpath(*UMAP_PROJECTION).exists()


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    state.setdefault("umap_genes", [])

    # ── The linked napari scatter window (unchanged) ──────────────────────
    show_umap_button = PushButton(label="Show UMAP Window", enabled=True)
    umap_size_slider = Slider(label="Point size", min=1, max=50, value=15)

    # ── Gene picker ───────────────────────────────────────────────────────
    gene_widget = ComboBox(
        label="Gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    add_gene_button = PushButton(label="Add Gene")
    remove_gene_button = PushButton(label="Remove Selected")
    clear_genes_button = PushButton(label="Clear All")

    gene_list_qt = QListWidget()
    gene_list_qt.setMaximumHeight(140)

    hint_label_qt = QLabel(
        f"Select up to {MAX_GENES} genes — one panel each, with a colour scale.")
    hint_label_qt.setWordWrap(True)
    hint_label_qt.setStyleSheet("color: gray;")

    cmap_widget = ComboBox(label="Colormap", choices=AVAILABLE_COLORMAPS,
                           value="viridis")
    ncols_widget = SpinBox(label="Columns", min=1, max=6, value=3)
    plot_genes_button = PushButton(label="Plot UMAP by gene")
    plot_clusters_button = PushButton(label="Plot UMAP by cluster")

    status_label = StatusProxy(ctx.viewer)

    # ── Napari window callbacks ───────────────────────────────────────────
    def on_show_umap():
        ctx.umap_viewer.show()
        ctx.record_node(
            "viewer:umap_window",
            "\n# Show UMAP scatter plot (viewer window; coords from Xenium analysis output)",
            deps=["preamble"],
            kind=NOTE,
            label="Show UMAP window",
        )

    def on_umap_size_change(value):
        ctx.umap_viewer.set_point_size(value / 100.0)

    # ── Gene list callbacks ───────────────────────────────────────────────
    def on_add_gene():
        gene = gene_widget.value
        genes = state["umap_genes"]
        if gene is None:
            return
        if gene in genes:
            status_label.value = f"'{gene}' already selected"
            return
        if len(genes) >= MAX_GENES:
            status_label.value = f"Maximum {MAX_GENES} genes reached"
            return
        genes.append(gene)
        gene_list_qt.addItem(gene)

    def on_remove_gene():
        selected = gene_list_qt.currentRow()
        if selected >= 0:
            gene_list_qt.takeItem(selected)
            state["umap_genes"].pop(selected)

    def on_clear_genes():
        state["umap_genes"].clear()
        gene_list_qt.clear()

    # ── Preview provider ──────────────────────────────────────────────────
    def _display_categories(clustering_key: str) -> dict | None:
        """Original cluster id -> display name, or None if unnamed."""
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

    def _umap_preview(mode: str = None) -> Preview:
        """What the two buttons would run, as the widgets stand.

        The blocks travel with the params because block selection *is* what the
        widgets mean here: which button was pressed picks the colouring, and
        whether the dataset ships an ``analysis/`` folder picks the embedding.
        A params-only provider would pin the pane to one assembly while the
        numbers tracked the widgets.
        """
        mode = mode or state.get("_umap_last_mode") or "genes"
        embed = "embed.xenium" if has_xenium_umap(ctx.data_path) else "embed.recompute"

        if mode == "clusters":
            key = ctx.clustering_widget.value if ctx.clustering_widget else None
            key = key or "clustering"
            categories = _display_categories(key)
            blocks = [embed]
            params = {"color": [key], "paths": ctx.plot_paths(_stem_for_clusters(key))}
            if categories is not None:
                blocks.append("relabel")
                params["groupby"] = key
                params["categories"] = categories
            blocks += ["color.clusters", "save"]
            return Preview(blocks, params)

        genes = list(state["umap_genes"]) or [
            g for g in ctx.gene_names[:1]] or ["GENE"]
        params = {
            "color": genes,
            "cmap": cmap_widget.value,
            "ncols": coerce(ncols_widget.value),
            "paths": ctx.plot_paths(_stem_for_genes(genes)),
        }
        return Preview([embed, "color.genes", "save"], params)

    state.setdefault("template_preview", {})[TEMPLATE_ID] = _umap_preview

    def _stem_for_genes(genes) -> str:
        """``umap_EPCAM_KRT5``, truncated so a 15-gene name stays a filename."""
        joined = "_".join(safe_stem(g) for g in genes[:4])
        suffix = f"_and_{len(genes) - 4}_more" if len(genes) > 4 else ""
        return f"umap_{joined}{suffix}"

    def _stem_for_clusters(key: str) -> str:
        return f"umap_{safe_stem(key)}"

    # ── Runners ───────────────────────────────────────────────────────────
    def _run(mode: str, node_id: str, deps: list, label: str, stem: str):
        blocks, params, _ = _umap_preview(mode)
        state["_umap_last_mode"] = mode
        step = Step(
            id=node_id,
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=deps,
            kind=TERMINAL,
            label=label,
            outputs=["fig"],
        )
        plot_genes_button.enabled = False
        plot_clusters_button.enabled = False
        gen = ctx.dataset_generation

        @thread_worker
        def _work():
            ctx.ensure_normalized()
            return ctx.run_step(step)

        def _done(out):
            if ctx.dataset_generation != gen:
                return
            plot_genes_button.enabled = True
            plot_clusters_button.enabled = True
            # The template wrote the files, so show_plot only publishes them.
            paths = ctx.show_plot(out["fig"], stem, title=label,
                                  save=False, paths=params["paths"])
            status_label.value = f"{label} — saved to {', '.join(paths)}"

        def _failed(exc):
            plot_genes_button.enabled = True
            plot_clusters_button.enabled = True
            status_label.value = f"UMAP plot failed: {exc}"

        worker = _work()
        worker.returned.connect(_done)
        worker.errored.connect(_failed)
        worker.start()

    def on_plot_genes():
        genes = list(state["umap_genes"])
        if not genes:
            status_label.value = "Add at least one gene first"
            return
        ctx.apply_plot_font_size()
        _run("genes", f"plot:umap_genes:{'_'.join(genes)}",
             deps=["normalize"],
             label=f"UMAP: {', '.join(genes)}",
             stem=_stem_for_genes(genes))

    def on_plot_clusters():
        key = ctx.clustering_widget.value if ctx.clustering_widget else None
        if not key or key not in ctx.clusterings:
            status_label.value = "No clustering selected"
            return
        # The column must exist in obs and as a node before a step can read it
        # and declare it as a dependency. Both on the GUI thread, and before the
        # params are built — the display categories come off that column.
        from xenium_viewer.utils.gene_analysis import add_clustering_to_obs
        ctx.record_clustering(key)
        add_clustering_to_obs(ctx.adata, ctx.adata, ctx.clusterings[key], key)
        ctx.apply_plot_font_size()
        _run("clusters", f"plot:umap:{key}",
             deps=["normalize", f"clustering:{key}"],
             label=f"UMAP: {key}",
             stem=_stem_for_clusters(key))

    show_umap_button.clicked.connect(on_show_umap)
    umap_size_slider.changed.connect(on_umap_size_change)
    add_gene_button.clicked.connect(on_add_gene)
    remove_gene_button.clicked.connect(on_remove_gene)
    clear_genes_button.clicked.connect(on_clear_genes)
    plot_genes_button.clicked.connect(on_plot_genes)
    plot_clusters_button.clicked.connect(on_plot_clusters)

    widget = make_tab(
        show_umap_button,
        umap_size_slider,
        hint_label_qt,
        gene_widget,
        add_gene_button,
        gene_list_qt,
        remove_gene_button,
        clear_genes_button,
        cmap_widget,
        ncols_widget,
        plot_genes_button,
        plot_clusters_button,
    )

    def _restore_session(session):
        genes = session.get("umap_genes")
        if genes:
            state["umap_genes"] = [g for g in genes if g in set(ctx.gene_names)]
            gene_list_qt.clear()
            for gene in state["umap_genes"]:
                gene_list_qt.addItem(gene)

    return widget, {"restore_session": _restore_session}
