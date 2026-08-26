"""Tab: Annotation Neighbourhood Enrichment.

Runs squidpy neighbourhood enrichment on real Xenium cells combined with
virtual "cells" sampled from user-drawn annotation polygons.  Annotation
types appear as additional rows/columns in the Z-score heatmap alongside
real cell-type clusters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, PushButton, Slider, SpinBox
from qtpy.QtWidgets import (
    QVBoxLayout, QLabel, QTextEdit, QFileDialog,
    QCheckBox, QGroupBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import StatusProxy, attach_tqdm_progress, qt_tqdm_context, make_progress_bar, combo_value_kwargs, make_tab
from xenium_viewer.utils.plot_output import safe_stem
from xenium_viewer.utils.prov_graph import NOTE

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    # ── Clustering selector ───────────────────────────────────────────────────
    clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.annot_nhood_clustering_widget = clustering_widget

    # ── Annotation type checkboxes ────────────────────────────────────────────
    annot_group_label = QLabel("Include annotation types as virtual cells:")
    annot_group = QGroupBox()
    annot_group_layout = QVBoxLayout()
    annot_group_layout.setContentsMargins(4, 4, 4, 4)
    annot_group.setLayout(annot_group_layout)

    _annot_checkboxes: list[QCheckBox] = []

    def _rebuild_annot_checkboxes():
        for cb in _annot_checkboxes:
            annot_group_layout.removeWidget(cb)
            cb.deleteLater()
        _annot_checkboxes.clear()

        from xenium_viewer.utils.annotation_utils import get_annotation_types
        types = get_annotation_types(ctx.annotation_layer)
        if not types:
            placeholder = QLabel("  (no annotation types defined — use the Annotations tab)")
            annot_group_layout.addWidget(placeholder)
            _annot_checkboxes.append(placeholder)  # type: ignore[arg-type]
            return
        for t in types:
            cb = QCheckBox(t)
            cb.setChecked(True)
            annot_group_layout.addWidget(cb)
            _annot_checkboxes.append(cb)

    _rebuild_annot_checkboxes()

    refresh_annot_btn = PushButton(label="Refresh annotation types")

    # ── Parameters ────────────────────────────────────────────────────────────
    density_spin = SpinBox(label="Grid density (µm²/virtual cell)", min=10, max=10000, value=100)
    perms_slider = Slider(label="Permutations", min=100, max=1000, value=1000)
    neighs_slider = Slider(label="Neighbours", min=3, max=20, value=6)

    # ── Controls ──────────────────────────────────────────────────────────────
    run_btn = PushButton(label="Run Annotation Nhood Enrichment", enabled=True)
    mode_widget = ComboBox(label="Display mode", choices=["zscore", "count"], value="zscore")
    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(250)
    plot_btn = PushButton(label="Show Heatmap", enabled=False)
    export_btn = PushButton(label="Export Z-scores CSV...", enabled=False)
    status = StatusProxy(ctx.viewer)
    annot_nhood_progress = make_progress_bar()

    from xenium_viewer.utils.gene_analysis import get_normalized_adata, add_clustering_to_obs
    from xenium_viewer.utils.spatial_analysis import (
        compute_spatial_neighbors, run_nhood_enrichment, make_nhood_enrichment_plot,
    )

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def _refresh_annot():
        _rebuild_annot_checkboxes()

    def _get_selected_annot_types() -> list[str]:
        return [
            cb.text() for cb in _annot_checkboxes
            if isinstance(cb, QCheckBox) and cb.isChecked()
        ]

    def _on_run():
        from xenium_viewer.utils.annotation_utils import sample_annotation_centroids
        import anndata
        import scipy.sparse as sp
        import pandas as pd

        clustering_key = clustering_widget.value
        if clustering_key is None or clustering_key not in ctx.clusterings:
            results_text.setPlainText("No clustering selected.")
            return

        annot_types = _get_selected_annot_types()
        if not annot_types:
            results_text.setPlainText("No annotation types selected.")
            return

        n_perms = perms_slider.value
        n_neighs = neighs_slider.value
        density = density_spin.value

        run_btn.enabled = False
        status.value = "Building augmented dataset with annotation virtual cells..."

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        _progress = [None]

        @thread_worker
        def _run():
            adata_norm = get_normalized_adata(_adata)
            add_clustering_to_obs(
                adata_norm, _adata, ctx.clusterings[clustering_key], clustering_key
            )
            adata_norm.obsm['spatial'] = _adata.obsm['spatial'].copy()

            # Sample virtual cells for each annotation type
            virtual_blocks = []
            for atype in annot_types:
                centroids = sample_annotation_centroids(
                    ctx.annotation_layer, atype, ctx.pixel_size, density_um2=density
                )
                if len(centroids) == 0:
                    continue
                virtual_blocks.append((atype, centroids))

            if not virtual_blocks:
                return {"warning": "No virtual cells could be sampled (no annotation polygons?)",
                        "zscore": np.array([]), "count": np.array([]), "clusters": []}

            # Build augmented AnnData
            n_real = adata_norm.n_obs
            n_genes = adata_norm.n_vars
            all_spatial = [adata_norm.obsm['spatial']]
            all_labels = list(adata_norm.obs[clustering_key].astype(str))

            for atype, centroids in virtual_blocks:
                all_spatial.append(centroids)
                all_labels += [atype] * len(centroids)

            spatial_aug = np.vstack(all_spatial)
            n_virtual = len(spatial_aug) - n_real

            X_real = adata_norm.X
            X_virtual = sp.csr_matrix((n_virtual, n_genes))

            adata_aug = anndata.AnnData(
                X=sp.vstack([X_real, X_virtual]),
                var=adata_norm.var.copy(),
            )
            adata_aug.obsm['spatial'] = spatial_aug

            # All categories including annotation types
            all_cats = sorted(set(all_labels))
            adata_aug.obs[clustering_key] = pd.Categorical(all_labels, categories=all_cats)

            compute_spatial_neighbors(adata_aug, n_neighs=n_neighs)
            with qt_tqdm_context(_progress[0], "Enrichment permutations: "):
                result = run_nhood_enrichment(adata_aug, clustering_key, n_perms=n_perms)

            result['_adata_norm'] = adata_aug
            result['_cluster_key'] = clustering_key
            result['_annot_types'] = annot_types
            return result

        worker = _run()
        _progress[0], state['_annot_nhood_progress_timer'] = attach_tqdm_progress(
            worker,
            lambda m: setattr(status, 'value', m),
            "Enrichment permutations: ",
            progress_bar=annot_nhood_progress,
        )
        worker.returned.connect(_on_done)
        worker.start()

    def _on_done(result):
        state["annot_nhood_result"] = result
        run_btn.enabled = True

        warning = result.get('warning')
        if warning:
            status.value = f"Annot nhood: {warning}"
            results_text.setPlainText(warning)
            plot_btn.enabled = False
            export_btn.enabled = False
            return

        zscore = result['zscore']
        clusters = result['clusters']
        n = len(clusters)
        annot_types = result.get('_annot_types', [])

        lines = [
            f"Annotation neighbourhood enrichment: {n}x{n} matrix",
            f"Clusters: {', '.join(clusters)}",
            f"Annotation types included: {', '.join(annot_types)}",
            "",
            "Note: virtual cells have zero gene expression; Z-scores reflect spatial",
            "proximity only (not gene-expression-derived similarity).",
            "",
        ]

        if zscore.size > 0:
            # Show top enrichments involving annotation types
            ann_set = set(annot_types)
            pairs = []
            for i, c1 in enumerate(clusters):
                for j, c2 in enumerate(clusters):
                    if i != j and (c1 in ann_set or c2 in ann_set):
                        pairs.append((c1, c2, zscore[i, j]))
            pairs.sort(key=lambda x: x[2], reverse=True)
            if pairs:
                lines.append("Top enrichments involving annotation types:")
                for c1, c2, z in pairs[:10]:
                    lines.append(f"  {c1} <-> {c2}: {z:.2f}")

        results_text.setPlainText("\n".join(lines))
        status.value = f"Annotation nhood enrichment done: {n} groups"
        plot_btn.enabled = True
        export_btn.enabled = True

    def _on_show_plot():
        result = state.get("annot_nhood_result")
        if result is None:
            return
        ctx.apply_plot_font_size()
        fig = make_nhood_enrichment_plot(result, mode=mode_widget.value)
        annot_types = result.get('_annot_types', [])
        fig.suptitle(f"Annotation Nhood Enrichment\n(annotation types: {', '.join(annot_types)})",
                     fontsize=10)
        fig.tight_layout()
        # ``plt.show()`` here was *blocking* — the only plot in the app that
        # froze the viewer until its window was closed.
        cluster_key = result.get('_cluster_key', '')
        paths = ctx.show_plot(
            fig, f"annot_nhood_enrichment_{safe_stem(cluster_key)}",
            title=f"Annotation nhood enrichment: {cluster_key}")
        status.value = f"Heatmap displayed — saved to {', '.join(paths)}"
        # This tab recorded nothing at all before, so the figure existed with no
        # trace in the notebook. It is a NOTE rather than a TERMINAL because the
        # enrichment is computed over virtual cells sampled from a napari shapes
        # layer the notebook has no access to: there is no code to record yet,
        # and a cell that ran `plt.gcf().savefig(...)` with no preceding plot
        # would replay as a silent no-op writing an empty figure. NOTE says
        # "not code, and not meant to be", renders as markdown, and is counted
        # apart from the comment-only punch list.
        ctx.record_node(
            "viewer:annot_nhood_plot",
            f"\n# Annotation neighbourhood enrichment heatmap "
            f"(mode={mode_widget.value}, annotation types: {', '.join(annot_types)})\n"
            f"# Computed over virtual cells sampled from annotation shapes drawn\n"
            f"# in the viewer, which this notebook cannot reach. Figure written to:\n"
            + "".join(f"#   {p}\n" for p in ctx.recorded_plot_paths(paths)),
            deps=["preamble"],
            kind=NOTE,
            label="Annotation nhood heatmap",
        )

    def _on_export():
        result = state.get("annot_nhood_result")
        if result is None:
            return
        zscore = result.get('zscore')
        clusters = result.get('clusters', [])
        if zscore is None or zscore.size == 0:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Z-scores CSV", "annot_nhood_zscore.csv", "CSV Files (*.csv)"
        )
        if not path:
            return
        import pandas as pd
        df = pd.DataFrame(zscore, index=clusters, columns=clusters)
        df.to_csv(path)
        status.value = f"Exported Z-scores to {path}"

    # ── Connect ───────────────────────────────────────────────────────────────
    refresh_annot_btn.changed.connect(_refresh_annot)
    run_btn.changed.connect(_on_run)
    plot_btn.changed.connect(_on_show_plot)
    export_btn.changed.connect(_on_export)

    # ── Session restore ────────────────────────────────────────────────────────
    def _restore_session(session):
        pass  # No persistent state for this tab

    # ── Build tab layout ──────────────────────────────────────────────────────
    # Was a hand-rolled QVBoxLayout + QScrollArea doing exactly what make_tab
    # does. That duplication is what kept this tab's labels invisible after the
    # fix landed in make_tab, so it goes through the helper like every other tab.
    scroll = make_tab(
        clustering_widget,
        annot_group_label,
        annot_group,
        refresh_annot_btn,
        density_spin,
        perms_slider,
        neighs_slider,
        run_btn,
        annot_nhood_progress,
        mode_widget,
        results_text,
        plot_btn,
        export_btn,
    )

    return scroll, {"restore_session": _restore_session}
