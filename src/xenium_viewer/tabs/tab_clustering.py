"""Tab 0: Clustering — Leiden, import/export, label editor."""

from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path

import os

from magicgui.widgets import CheckBox, PushButton, Slider, FloatSpinBox
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar
from xenium_viewer.utils.prov_graph import ARTIFACT, TERMINAL
from xenium_viewer.utils.steps import Step, coerce

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


# Leiden runs as a single Step: the source below is what the viewer executes
# *and* what the notebook records — there is no second expression of this
# pipeline to drift from (see utils/steps.py).
#
# It starts from ``adata_norm`` (the log-normalised copy bound by the
# "normalize" step) rather than re-normalising ``adata`` itself, so the graph
# carries a real ``normalize -> clustering`` edge and the notebook normalises
# once. This is safe only because ``normalize`` binds a *copy*: when it mutated
# ``adata`` in place, any step that copied ``adata`` risked normalising twice.
#
# It still works on its own ``adata_leiden`` copy, so that neighbours/leiden
# (and any HVG subsetting) don't mutate the shared ``adata_norm`` that
# rank-genes and the other expression analyses read.
_LEIDEN_TEMPLATE_HEAD = """
# Leiden clustering ($key)
adata_leiden = adata_norm.copy()"""

_LEIDEN_TEMPLATE_HVG = """
sc.pp.highly_variable_genes(adata_leiden, n_top_genes=$n_top_genes, flavor='seurat')
adata_leiden = adata_leiden[:, adata_leiden.var.highly_variable].copy()"""

_LEIDEN_TEMPLATE_SCALE = """
sc.pp.scale(adata_leiden, max_value=10)"""

# PCA is recomputed only when the gene set or scaling changed; otherwise the
# X_pca carried over from `normalize` is exactly what we would recompute.
_LEIDEN_TEMPLATE_PCA = """
sc.pp.pca(adata_leiden)"""

_LEIDEN_TEMPLATE_TAIL = """
sc.pp.neighbors(adata_leiden, n_neighbors=$n_neighbors, n_pcs=$n_pcs)
sc.tl.leiden(
    adata_leiden, resolution=$resolution, key_added=$key,
    flavor='igraph', n_iterations=2, random_state=$random_state,
)
adata.obs[$key] = adata_leiden.obs[$key].values"""


def _leiden_template(use_hvg: bool, do_scale: bool) -> str:
    """Assemble the Leiden template for the selected preprocessing options."""
    parts = [_LEIDEN_TEMPLATE_HEAD]
    if use_hvg:
        parts.append(_LEIDEN_TEMPLATE_HVG)
    if do_scale:
        parts.append(_LEIDEN_TEMPLATE_SCALE)
    if use_hvg or do_scale:
        parts.append(_LEIDEN_TEMPLATE_PCA)
    parts.append(_LEIDEN_TEMPLATE_TAIL)
    return "".join(parts)


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
    leiden_progress = make_progress_bar()

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

        # Recording happened inside ctx.run_step(), which recorded the very
        # source it executed — there is nothing to re-describe here.

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

        key = f"leiden_r{resolution}"
        step = Step(
            id=f"clustering:{key}",
            template=_leiden_template(use_hvg, do_scale),
            params={
                "key": key,
                "resolution": coerce(resolution),
                "n_neighbors": coerce(n_neighbors),
                "n_pcs": coerce(n_pcs),
                "n_top_genes": coerce(n_hvgs),
                "random_state": 0,
            },
            deps=["normalize"],
            kind=ARTIFACT,
            label=f"Clustering: {key}",
        )

        @thread_worker
        def _run():
            import pandas as pd
            yield "Normalizing..."
            # Binds adata_norm and records the "normalize" node this step
            # declares as its dependency. Idempotent per adata.
            ctx.ensure_normalized()

            yield "Running Leiden clustering..."
            # One call: this executes exactly the source it records.
            ctx.run_step(step)

            labels = _adata.obs[key]
            cell_ids = (_adata.obs['cell_id'].values
                        if 'cell_id' in _adata.obs.columns else _adata.obs_names)
            series = pd.Series(
                labels.astype(int).values, index=cell_ids, name="leiden",
            )
            n_clusters = series.nunique()
            return series, n_clusters, resolution, n_neighbors, n_pcs, use_hvg, do_scale, n_hvgs

        worker = _run()
        worker.returned.connect(lambda result: _on_leiden_ready(result, gen))
        timer, update_msg = attach_spinner(
            worker,
            lambda m: setattr(leiden_status, 'value', m),
            "Preparing data...",
            progress_bar=leiden_progress,
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
        ctx.record_node(
            f"clustering:{name}",
            f"\n# Import clustering '{name}' from file "
            f"({series.nunique()} groups, {len(series)} cells)\n"
            f"_imp = pd.read_csv(r\"{path}\", sep=None, engine=\"python\")\n"
            f"_idx = \"cell_id\" if \"cell_id\" in _imp.columns else _imp.columns[0]\n"
            f"_col = \"group\" if \"group\" in _imp.columns else _imp.columns[1]\n"
            f"adata.obs[\"{name}\"] = pd.Categorical("
            f"_imp.set_index(_idx)[_col].astype(str).reindex(adata.obs_names).values)",
            deps=["preamble"],
            label=f"Import clustering: {name}",
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
        ctx.record_clustering(clustering_key)
        ctx.record_node(
            f"export:clustering:{clustering_key}",
            f"\n# Export clustering '{clustering_key}'\n"
            f"pd.DataFrame({{\"cell_id\": adata.obs_names, "
            f"\"group\": adata.obs[\"{clustering_key}\"].values}})"
            f".to_csv(\"{os.path.basename(path)}\", index=False, sep={sep!r})",
            deps=[f"clustering:{clustering_key}"],
            kind=TERMINAL,
            label=f"Export clustering: {clustering_key}",
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
        leiden_progress,
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
