"""Tab 0: Clustering — Leiden, import/export, label editor."""

from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path

from magicgui.widgets import CheckBox, PushButton, Slider, FloatSpinBox
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

from xenium_viewer.utils.gene_analysis import get_normalized_adata


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    leiden_n_neighbors = Slider(label="n_neighbors", min=5, max=50, value=15)
    leiden_n_pcs = Slider(label="n_pcs", min=10, max=50, value=40)
    leiden_resolution = FloatSpinBox(label="resolution", min=0.1, max=5.0, step=0.1, value=1.0)
    leiden_hvg_check = CheckBox(label="Use HVGs only", value=False)
    leiden_n_hvgs = Slider(label="n_top_genes", min=500, max=4000, value=2000, enabled=False)
    leiden_scale_check = CheckBox(label="Scale (max_value=10)", value=False)

    def _on_hvg_toggle(val):
        leiden_n_hvgs.enabled = val
    leiden_hvg_check.changed.connect(_on_hvg_toggle)

    leiden_run_button = PushButton(label="Run Leiden Clustering", enabled=True)
    leiden_import_button = PushButton(label="Import Clustering...", enabled=True)
    leiden_export_button = PushButton(label="Export Clustering...", enabled=True)

    leiden_status_text = QTextEdit()
    leiden_status_text.setReadOnly(True)
    leiden_status_text.setFontFamily("monospace")
    leiden_status_text.setMaximumHeight(150)

    leiden_status = StatusProxy(ctx.viewer)

    def _on_leiden_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        series, n_clusters, resolution, n_neighbors, n_pcs, use_hvg, do_scale, n_hvgs = result
        key = f"leiden_r{resolution}"
        ctx.clusterings[key] = series
        state["custom_clusterings"][key] = series
        ctx.refresh_clustering_choices()

        leiden_status_text.setPlainText(
            f"Leiden clustering complete\n"
            f"  Key: {key}\n"
            f"  Clusters: {n_clusters}\n"
            f"  n_neighbors: {n_neighbors}\n"
            f"  n_pcs: {n_pcs}\n"
            f"  resolution: {resolution}\n"
            f"  HVGs: {n_hvgs if use_hvg else 'all genes'}\n"
            f"  Scaled: {'yes (max=10)' if do_scale else 'no'}"
        )
        leiden_status.value = f"Leiden done: {n_clusters} clusters ({key})"
        leiden_run_button.enabled = True

        from xenium_viewer.utils.adata_persistence import save_clustering_to_adata
        save_clustering_to_adata(ctx, key, series)

        if use_hvg or do_scale:
            code_lines = [
                "\n# Custom preprocessing for Leiden",
                "adata_leiden = adata.copy()",
                "sc.pp.normalize_total(adata_leiden, target_sum=1e4)",
                "sc.pp.log1p(adata_leiden)",
            ]
            if use_hvg:
                code_lines.append(f"sc.pp.highly_variable_genes(adata_leiden, n_top_genes={n_hvgs}, flavor='seurat')")
                code_lines.append("adata_leiden = adata_leiden[:, adata_leiden.var.highly_variable].copy()")
            if do_scale:
                code_lines.append("sc.pp.scale(adata_leiden, max_value=10)")
            code_lines.append("sc.pp.pca(adata_leiden)")
            code_lines.append(f"sc.pp.neighbors(adata_leiden, n_neighbors={n_neighbors}, n_pcs={n_pcs})")
            code_lines.append(f'sc.tl.leiden(adata_leiden, resolution={resolution}, key_added="{key}")')
            ctx.record_code("\n".join(code_lines), tag=f"leiden_{key}")
        else:
            ctx.record_normalize()
            ctx.record_code(
                f"\n# Leiden clustering (n_neighbors={n_neighbors}, n_pcs={n_pcs}, resolution={resolution})\n"
                f"sc.pp.neighbors(adata, n_neighbors={n_neighbors}, n_pcs={n_pcs})\n"
                f"sc.tl.leiden(adata, resolution={resolution}, key_added=\"{key}\")",
                tag=f"leiden_{key}"
            )

    def on_run_leiden():
        n_neighbors = leiden_n_neighbors.value
        n_pcs = leiden_n_pcs.value
        resolution = leiden_resolution.value
        use_hvg = leiden_hvg_check.value
        do_scale = leiden_scale_check.value
        n_hvgs = leiden_n_hvgs.value
        leiden_run_button.enabled = False
        leiden_status.value = "Running Leiden clustering..."
        gen = ctx.dataset_generation

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata

        @thread_worker
        def _run():
            import scanpy as sc
            yield "Preparing data..."
            if use_hvg or do_scale:
                adata_work = _adata.copy()
                sc.pp.normalize_total(adata_work, target_sum=1e4)
                sc.pp.log1p(adata_work)
                if use_hvg:
                    sc.pp.highly_variable_genes(adata_work, n_top_genes=n_hvgs, flavor="seurat")
                    adata_work = adata_work[:, adata_work.var.highly_variable].copy()
                if do_scale:
                    sc.pp.scale(adata_work, max_value=10)
                sc.pp.pca(adata_work)
            else:
                adata_work = get_normalized_adata(_adata)
            yield "Computing neighbors..."
            sc.pp.neighbors(adata_work, n_neighbors=n_neighbors, n_pcs=n_pcs)
            yield "Running Leiden algorithm..."
            sc.tl.leiden(adata_work, resolution=resolution, key_added="leiden")
            import pandas as pd
            cell_ids = _adata.obs['cell_id'].values if 'cell_id' in _adata.obs.columns else _adata.obs_names
            series = pd.Series(
                adata_work.obs['leiden'].astype(int).values,
                index=cell_ids,
                name="leiden",
            )
            n_clusters = series.nunique()
            return series, n_clusters, resolution, n_neighbors, n_pcs, use_hvg, do_scale, n_hvgs

        worker = _run()
        worker.returned.connect(lambda result: _on_leiden_ready(result, gen))
        timer, update_msg = attach_spinner(
            worker,
            lambda m: setattr(leiden_status, 'value', m),
            "Preparing data...",
        )
        state['_spinner_timer'] = timer  # prevent GC
        worker.yielded.connect(update_msg)
        worker.start()

    leiden_run_button.clicked.connect(on_run_leiden)

    # ── Import / Export callbacks ─────────────────────────────────────────
    def _on_import_clustering():
        import pandas as pd
        path, _ = QFileDialog.getOpenFileName(
            None, "Import Clustering", "",
            "CSV/TSV Files (*.csv *.tsv *.txt);;All Files (*)",
        )
        if not path:
            return
        df = pd.read_csv(path, sep=None, engine='python')
        if 'cell_id' in df.columns and 'group' in df.columns:
            series = pd.Series(df['group'].values, index=df['cell_id'].values)
        else:
            series = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0].values)
        name = Path(path).stem
        ctx.clusterings[name] = series
        state["custom_clusterings"][name] = series
        ctx.refresh_clustering_choices()
        ctx.clustering_widget.value = name
        leiden_status.value = f"Imported '{name}' ({series.nunique()} groups, {len(series)} cells)"

        from xenium_viewer.utils.adata_persistence import save_clustering_to_adata
        save_clustering_to_adata(ctx, name, series)
        ctx.record_code(
            f"\n# Import clustering from file\n"
            f"# name={name}, {series.nunique()} groups, {len(series)} cells\n"
            f"# source: \"{path}\""
        )

    def _on_export_clustering():
        import pandas as pd
        clustering_key = ctx.clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            leiden_status.value = "No clustering selected"
            return
        series = ctx.clusterings[clustering_key]
        labels = ctx.get_active_labels()
        if labels:
            mapped = series.map(lambda x: labels.get(x, labels.get(str(x),
                                labels.get(int(x) if str(x).lstrip('-').isdigit() else x, x))))
        else:
            mapped = series
        df = pd.DataFrame({'cell_id': mapped.index, 'group': mapped.values})
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Clustering", f"{clustering_key}.csv",
            "CSV Files (*.csv);;TSV Files (*.tsv)",
        )
        if not path:
            return
        sep = '\t' if path.endswith('.tsv') else ','
        df.to_csv(path, index=False, sep=sep)
        leiden_status.value = f"Exported {len(df)} cells to {path}"
        ctx.record_code(
            f"\n# Export clustering\n"
            f"# {clustering_key} -> \"{path}\""
        )

    leiden_import_button.clicked.connect(_on_import_clustering)
    leiden_export_button.clicked.connect(_on_export_clustering)

    # ── Layout ───────────────────────────────────────────────────────────
    leiden_io_row = QWidget()
    leiden_io_layout = QHBoxLayout()
    leiden_io_layout.setContentsMargins(0, 0, 0, 0)
    leiden_io_layout.addWidget(leiden_import_button.native)
    leiden_io_layout.addWidget(leiden_export_button.native)
    leiden_io_row.setLayout(leiden_io_layout)

    widget = make_tab(
        leiden_n_neighbors,
        leiden_n_pcs,
        leiden_resolution,
        leiden_hvg_check,
        leiden_n_hvgs,
        leiden_scale_check,
        leiden_run_button,
        leiden_status_text,
        leiden_io_row,
    )

    def _restore_session(session):
        # Custom clusterings are now loaded from adata.obs at startup;
        # sync them into state["custom_clusterings"] for compatibility
        from xenium_viewer.utils.adata_persistence import load_custom_clusterings_from_adata
        cc = load_custom_clusterings_from_adata(ctx.adata)
        if cc:
            for name, series in cc.items():
                state["custom_clusterings"][name] = series
            ctx.refresh_clustering_choices()
            print(f"  Restored {len(cc)} custom clustering(s) from adata.obs")

        cl = session.get("cluster_labels")
        if cl and isinstance(cl, dict):
            state["cluster_labels"] = cl
            n_clusterings = len(cl)
            n_labels = sum(len(v) for v in cl.values() if isinstance(v, dict))
            print(f"  Restored cluster labels: {n_labels} labels across {n_clusterings} clustering(s)")

    return widget, {"restore_session": _restore_session}
