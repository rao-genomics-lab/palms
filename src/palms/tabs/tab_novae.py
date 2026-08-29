"""Tab: Spatial Domains — Novae zero-shot domain inference."""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit
from napari.qt.threading import thread_worker
from palms.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    species_widget = ComboBox(
        label="Species",
        choices=["human", "mouse"],
        value="human",
    )

    n_domains_slider = Slider(label="N domains (0=auto)", min=0, max=30, value=0)

    level_slider = Slider(label="Level", min=1, max=15, value=7)

    run_button = PushButton(label="Run Novae Domains", enabled=True)

    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(200)

    novae_status = StatusProxy(ctx.viewer)
    novae_progress = make_progress_bar()

    def on_run_novae():
        if ctx.adata is None:
            novae_status.value = "No AnnData loaded"
            return

        run_button.enabled = False
        results_text.setPlainText("Running Novae domain inference...")

        species = species_widget.value
        n_domains = n_domains_slider.value
        level = level_slider.value
        _adata = ctx.adata
        gen = ctx.dataset_generation

        @thread_worker
        def _run():
            try:
                import novae
            except ImportError:
                raise ImportError(
                    "novae is not installed. Run: pip install novae"
                )

            _adata_copy = _adata.copy()
            novae.spatial_neighbors(_adata_copy)

            model = novae.Novae.from_pretrained(f"MICS-Lab/novae-{species}-0")
            model.compute_representations(_adata_copy, zero_shot=True)

            n_dom = n_domains if n_domains > 0 else None
            model.assign_domains(_adata_copy, n_domains=n_dom, level=level)

            # novae names the column "novae_domain" or "novae_domains_N" depending on version
            domain_cols = [c for c in _adata_copy.obs.columns if c.startswith("novae_domain")]
            if not domain_cols:
                raise KeyError("No 'novae_domain*' column found in adata.obs after assign_domains()")
            # prefer exact match, otherwise take the last one added
            col = "novae_domain" if "novae_domain" in domain_cols else domain_cols[-1]
            series = _adata_copy.obs[col].copy()
            # Re-index by cell barcode so coloring lookups (reindex(cell_ids)) work
            if 'cell_id' in _adata_copy.obs.columns:
                series.index = _adata_copy.obs['cell_id'].values
            series.name = "novae_domains"
            return series, series.nunique()

        worker = _run()
        _timer, _ = attach_spinner(
            worker,
            lambda m: setattr(novae_status, "value", m),
            "Running Novae...",
            progress_bar=novae_progress,
        )
        state["_novae_spinner_timer"] = _timer  # keep reference to prevent GC

        worker.returned.connect(lambda result: _on_novae_ready(result, gen))
        worker.errored.connect(_on_novae_error)
        worker.start()

    def _on_novae_ready(result, gen):
        run_button.enabled = True
        if ctx.dataset_generation != gen:
            return  # dataset reloaded while worker ran

        series, n_unique = result
        key = "novae_domains"

        # Store in clusterings (invalidate stale color cache for this key first)
        ctx.color_manager.invalidate_cluster_cache(key)
        ctx.clusterings[key] = series
        if "custom_clusterings" not in state:
            state["custom_clusterings"] = {}
        state["custom_clusterings"][key] = series

        # Refresh dropdowns
        ctx.refresh_clustering_choices()

        # Auto-apply coloring in background thread
        @thread_worker
        def _apply_colors():
            color_arr, cluster_to_color = ctx.color_manager.get_cluster_colors(series)
            cluster_ids_per_obs, label_to_cluster = ctx.get_cluster_ids_per_obs(key)
            colormap = ctx.color_manager.build_direct_label_colormap(color_arr)
            return colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs

        def _on_colors_ready(color_result):
            colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs = color_result
            state["cluster_to_color"] = cluster_to_color
            state["label_to_cluster"] = label_to_cluster
            state["active_clustering_name"] = key
            if ctx.cell_labels_layer is not None:
                ctx.cell_labels_layer.colormap = colormap
                ctx.cell_labels_layer.refresh()
            ctx.umap_viewer.color_by_cluster(
                key, color_arr, ctx.label_to_obs,
                cluster_ids_per_obs=cluster_ids_per_obs,
            )
            novae_status.value = f"Novae domains applied: {n_unique} domains"

        color_worker = _apply_colors()
        color_worker.returned.connect(_on_colors_ready)
        color_worker.start()

        # Save to adata.obs
        from palms.utils.adata_persistence import save_clustering_to_adata
        save_clustering_to_adata(ctx, key, series)

        # Record code
        species = species_widget.value
        n_domains = n_domains_slider.value
        level = level_slider.value
        n_dom_arg = n_domains if n_domains > 0 else None
        ctx.record_preamble()
        # Recorded as ``clustering:novae_domains`` — the id of the artifact it
        # produces. Under the old ``novae`` id nothing could declare a dependency
        # on these domains, so any tab that analysed them fell through to the
        # generic clustering fallback and recorded a loader for a CSV that does
        # not exist. The column is renamed here too: novae writes
        # ``novae_domain`` (or ``novae_domains_N``), the viewer stores it as
        # ``novae_domains``, and the recorded cell claimed the former.
        ctx.record_node(
            "clustering:novae_domains",
            f"\n# Novae spatial domain inference\n"
            f"import novae\n"
            f"adata_novae = adata.copy()\n"
            f"novae.spatial_neighbors(adata_novae)\n"
            f"model = novae.Novae.from_pretrained('MICS-Lab/novae-{species}-0')\n"
            f"model.compute_representations(adata_novae, zero_shot=True)\n"
            f"model.assign_domains(adata_novae, n_domains={n_dom_arg!r}, level={level})\n"
            f"# novae's column name varies by version: 'novae_domain', else the last added\n"
            f"_dom = [c for c in adata_novae.obs.columns if c.startswith('novae_domain')]\n"
            f"_col = 'novae_domain' if 'novae_domain' in _dom else _dom[-1]\n"
            f"adata.obs['novae_domains'] = pd.Categorical(adata_novae.obs[_col].values)",
            deps=["preamble"],
            label="Novae spatial domains",
        )

        results_text.setPlainText(
            f"Novae domain inference complete\n"
            f"  Domains: {n_unique}\n"
            f"  Species model: {species}\n"
            f"  Level: {level}\n"
            f"  N domains requested: {'auto' if n_domains == 0 else n_domains}\n"
            f"\nResult stored as: novae_domains\n"
            f"Use Cell Coloring tab to re-apply or switch colorings."
        )

    def _on_novae_error(exc):
        run_button.enabled = True
        novae_status.value = f"Novae error: {exc}"
        results_text.setPlainText(f"Error running Novae:\n{exc}")

    run_button.clicked.connect(on_run_novae)

    widget = make_tab(
        species_widget,
        n_domains_slider,
        level_slider,
        run_button,
        novae_progress,
        results_text,
    )

    return widget, {}
