"""Tab 3: UMAP — show UMAP window + point size slider + save scanpy UMAP plot."""

from __future__ import annotations
from typing import TYPE_CHECKING

import os

from magicgui.widgets import PushButton, Slider, ComboBox
from qtpy.QtWidgets import QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy
from xenium_viewer.utils.prov_graph import NOTE, TERMINAL

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    show_umap_button = PushButton(label="Show UMAP Window", enabled=True)
    umap_size_slider = Slider(label="UMAP pt size", min=1, max=50, value=15)
    umap_save_format = ComboBox(label="Save format", choices=["PNG", "SVG"], value="PNG")
    save_umap_button = PushButton(label="Save UMAP Plot...", enabled=True)
    status_label = StatusProxy(ctx.viewer)

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

    def on_save_umap():
        clustering_key = ctx.clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            status_label.value = "No clustering selected"
            return

        fmt = umap_save_format.value.lower()
        path, _ = QFileDialog.getSaveFileName(
            None, "Save UMAP Plot",
            f"umap_{clustering_key}.{fmt}",
            f"{umap_save_format.value} Files (*.{fmt})",
        )
        if not path:
            return

        save_umap_button.enabled = False
        status_label.value = "Saving UMAP plot..."
        gen = ctx.dataset_generation

        _adata          = ctx.adata
        _clustering_key = clustering_key
        _clustering_ser = ctx.clusterings[clustering_key].copy()
        _labels         = ctx.get_labels_for(clustering_key)
        _fmt            = fmt

        @thread_worker
        def _run():
            import pandas as pd
            import anndata as ad
            import scanpy as sc
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from xenium_viewer.utils.gene_analysis import add_clustering_to_obs

            # Minimal AnnData — only what sc.pl.umap needs
            plot_adata = ad.AnnData(obs=pd.DataFrame(index=_adata.obs_names))
            plot_adata.obsm['X_umap'] = _adata.obsm['X_umap'].copy()

            add_clustering_to_obs(plot_adata, _adata, _clustering_ser,
                                  key_name=_clustering_key)

            # Rename categories to cluster labels if available
            if _labels:
                cats = plot_adata.obs[_clustering_key].cat.categories
                new_cats = []
                for c in cats:
                    try:
                        label = _labels.get(int(c), _labels.get(c, c))
                    except (ValueError, TypeError):
                        label = _labels.get(c, c)
                    new_cats.append(str(label))
                plot_adata.obs[_clustering_key] = (
                    plot_adata.obs[_clustering_key].cat.rename_categories(new_cats)
                )

            sc.pl.umap(
                plot_adata,
                color=_clustering_key,
                legend_loc='on data',
                show=False,
                title=f'UMAP — {_clustering_key}',
            )
            fig = plt.gcf()
            save_kwargs = dict(bbox_inches='tight')
            if _fmt == 'png':
                save_kwargs['dpi'] = 150
            fig.savefig(path, **save_kwargs)
            plt.close(fig)
            return path

        worker = _run()
        worker.returned.connect(
            lambda p: _on_save_done(p) if ctx.dataset_generation == gen else None
        )
        worker.start()

    def _on_save_done(path):
        save_umap_button.enabled = True
        status_label.value = f"UMAP plot saved: {path}"
        _ck = ctx.clustering_widget.value
        _fmt = os.path.splitext(path)[1].lstrip('.') or 'png'
        ctx.record_clustering(_ck)
        ctx.record_node(
            "plot:umap",
            f"\n# UMAP embedding + plot colored by '{_ck}'\n"
            f"# (the viewer uses Xenium-provided UMAP coords; here we recompute so it replays)\n"
            f"sc.pp.neighbors(adata)\n"
            f"sc.tl.umap(adata, random_state=0)\n"
            f"sc.pl.umap(adata, color='{_ck}', legend_loc='on data', show=False)\n"
            f"plt.gcf().savefig('umap.{_fmt}', bbox_inches='tight')",
            deps=[f"clustering:{_ck}"],
            kind=TERMINAL,
            label="UMAP plot",
        )

    show_umap_button.clicked.connect(on_show_umap)
    umap_size_slider.changed.connect(on_umap_size_change)
    save_umap_button.clicked.connect(on_save_umap)

    widget = make_tab(show_umap_button, umap_size_slider, umap_save_format, save_umap_button)
    return widget, {}
