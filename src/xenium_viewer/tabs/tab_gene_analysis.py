"""Tab 6: Gene Analysis — rank genes, dotplot, volcanos."""

from __future__ import annotations
import json
import os
from pathlib import Path
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar, combo_value_kwargs
from xenium_viewer.utils.prov_graph import TERMINAL

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

from xenium_viewer.utils.gene_analysis import (
    add_clustering_to_obs,
    make_rank_genes_dotplot, make_rank_genes_plot, generate_all_volcano_plots,
    run_celltypist_annotation,
    load_reference_h5ad, get_annotation_columns, run_label_transfer,
    run_llm_annotation,
)
from xenium_viewer.utils.steps import Step, coerce
from xenium_viewer.utils.step_templates import (
    Preview, builtin_spec, builtin_text, step_template as _resolved,
)

RANK_GENES_TEMPLATE_ID = "genes.rank_genes"

# Executed *and* recorded — see utils/steps.py. Reads the clustering from
# adata.obs (where the clustering cell puts it) and ranks on the normalized
# copy bound by the "normalize" step, which is what the viewer has always
# actually done; the previous recorded cell ranked on raw `adata`.
_RANK_GENES_TEMPLATE = builtin_text("genes.rank_genes")


def dotplot_code(groupby: str, n_genes: int, dendrogram: bool, fmt: str) -> str:
    """The recorded dotplot cell.

    ``adata_norm``, not ``adata``: the ranked genes live on the normalised copy
    (see ``_RANK_GENES_TEMPLATE``). Recorded against ``adata`` this cell died on
    replay with ``KeyError: 'rank_genes_groups'`` — found by
    ``scripts/verify_notebook.py``, not by reading it. Module-level so a test can
    execute the exact string the viewer records.
    """
    return (
        f"\n# Dotplot (n_genes={n_genes}, dendrogram={dendrogram})\n"
        + (f"sc.tl.dendrogram(adata_norm, groupby=\"{groupby}\")\n"
           if dendrogram else "")
        + f"dotplot = sc.pl.rank_genes_groups_dotplot(\n"
        f"    adata_norm, groupby=\"{groupby}\", n_genes={n_genes},\n"
        f"    dendrogram={dendrogram}, show=False, return_fig=True,\n"
        f")\n"
        f"dotplot.savefig(\"dotplot.{fmt}\", dpi=300, bbox_inches='tight')"
    )


_APP_DIR = Path(__file__).resolve().parent.parent
_REF_DIR = _APP_DIR / "reference_datasets"


def _load_dataset_registry() -> dict:
    """Discover reference datasets from .metadata.json sidecar files."""
    registry = {}
    for meta_path in sorted(_REF_DIR.glob("*.metadata.json")):
        try:
            with open(meta_path) as f:
                meta = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue
        h5ad_path = _REF_DIR / meta["filename"]
        if not h5ad_path.exists():
            continue  # skip if h5ad not yet downloaded
        display = meta.get("display_name", meta["filename"])
        registry[display] = {
            "path": str(h5ad_path),
            "default_col": meta.get("default_col"),
            "metadata": meta,
        }
    registry["Browse..."] = {"path": None, "default_col": None, "metadata": None}
    return registry


_REFERENCE_DATASETS = _load_dataset_registry()

# Check if celltypist is available
try:
    import celltypist
    from celltypist import models as ct_models
    _HAS_CELLTYPIST = True
except ImportError:
    _HAS_CELLTYPIST = False


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    ga_clustering_widget = ComboBox(
        label="Clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.ga_clustering_widget = ga_clustering_widget

    ga_method_widget = ComboBox(
        label="Method", choices=["wilcoxon", "t-test", "logreg"], value="wilcoxon",
    )
    ga_n_genes_slider = Slider(label="Top N genes", min=5, max=50, value=25)
    ga_run_button = PushButton(label="Run Rank Genes", enabled=True)

    ga_dotplot_n_slider = Slider(label="Genes per cluster", min=3, max=20, value=5)
    ga_dendro_check = CheckBox(label="Dendrogram", value=True)
    ga_dotplot_button = PushButton(label="Show Dotplot", enabled=False)
    ga_edit_labels_button = PushButton(label="Edit Cluster Labels...", enabled=False)
    ga_reset_labels_button = PushButton(label="Reset Labels", enabled=False)

    ga_rank_plot_button = PushButton(label="Show Rank Genes Plot", enabled=False)

    ga_results_text = QTextEdit()
    ga_results_text.setReadOnly(True)
    ga_results_text.setFontFamily("monospace")
    ga_results_text.setMaximumHeight(300)
    ga_export_button = PushButton(label="Export Full Results CSV...", enabled=False)
    ga_volcano_button = PushButton(label="Generate All Volcano Plots...", enabled=False)

    # -- CellTypist annotation widgets --
    from qtpy.QtWidgets import QLabel
    ga_ct_separator = QLabel("── CellTypist Annotation ──")
    ga_ct_separator.setStyleSheet("font-weight: bold; margin-top: 8px;")
    ga_ct_model_widget = ComboBox(label="CellTypist Model", choices=[])
    from magicgui.widgets import FloatSlider
    ga_ct_conf_slider = FloatSlider(label="Min confidence", min=0.0, max=1.0, value=0.5, step=0.05)
    ga_ct_download_button = PushButton(label="Download Models")
    ga_ct_annotate_button = PushButton(label="Annotate with CellTypist", enabled=False)

    if not _HAS_CELLTYPIST:
        ga_ct_model_widget.enabled = False
        ga_ct_conf_slider.enabled = False
        ga_ct_download_button.enabled = False
        ga_ct_annotate_button.enabled = False
        ga_ct_separator.setText("── CellTypist (not installed) ──")
    else:
        # Populate with locally available models
        try:
            _local_models = sorted(Path(ct_models.models_path).glob("*.pkl"))
            _model_names = [m.name for m in _local_models]
            if _model_names:
                ga_ct_model_widget.choices = _model_names
                ga_ct_annotate_button.enabled = True
        except Exception:
            pass

    ga_status = StatusProxy(ctx.viewer)
    ga_progress = make_progress_bar()

    def _rank_genes_preview() -> Preview:
        """What "Run Rank Genes" would run with the widgets as they stand.

        One expression of the current settings, called by the run below and by
        the Templates tab's preview pane. Every block, always: this template
        declares a single assembly, so there is nothing for the widgets to
        select between.
        """
        return Preview(
            list(builtin_spec(RANK_GENES_TEMPLATE_ID).blocks),
            {
                "groupby": ga_clustering_widget.value,
                "method": ga_method_widget.value,
                "n_genes": coerce(ga_n_genes_slider.value),
            },
        )

    ctx.state.setdefault(
        "template_preview", {})[RANK_GENES_TEMPLATE_ID] = _rank_genes_preview

    def on_run_rank_genes():
        ga_status.value = "Running rank genes (normalizing + computing)..."
        ga_run_button.enabled = False

        blocks, params, _ = _rank_genes_preview()
        clustering_key = params["groupby"]
        method = params["method"]
        n_genes = params["n_genes"]
        state["_rg_method"] = method
        state["_rg_n_genes"] = n_genes

        if not ctx.clusterings or not clustering_key or clustering_key not in ctx.clusterings:
            ga_status.value = "Error: clustering data not available. Please wait for data to finish loading."
            ga_run_button.enabled = True
            return

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        _clustering_series = ctx.clusterings[clustering_key]

        # The step declares clustering:<key> as a dependency, so the node must
        # exist first. Recorded on the GUI thread — it writes analysis.py and
        # refreshes the Notebook tab.
        ctx.record_clustering(clustering_key)
        # Mirror the clustering onto adata.obs, which is where the recorded
        # clustering cell puts it and where the rank-genes step reads it from.
        add_clustering_to_obs(_adata, _adata, _clustering_series, clustering_key)

        step = Step(
            id=f"rank_genes:{clustering_key}",
            **_resolved(RANK_GENES_TEMPLATE_ID, blocks),
            params=params,
            deps=["normalize", f"clustering:{clustering_key}"],
            label=f"Rank genes: {clustering_key}",
            outputs=["rank_df", "adata_norm"],
        )

        @thread_worker
        def _run():
            ctx.ensure_normalized()          # binds adata_norm, records "normalize"
            outputs = ctx.run_step(step)     # executes exactly what it records
            return outputs["rank_df"], outputs["adata_norm"], clustering_key

        worker = _run()
        timer, _ = attach_spinner(
            worker,
            lambda m: setattr(ga_status, 'value', m),
            "Running rank genes...",
            progress_bar=ga_progress,
        )
        state['_spinner_timer'] = timer  # prevent GC
        worker.returned.connect(_on_rank_genes_ready)
        worker.start()

    def _on_rank_genes_ready(result):
        df, adata_norm, clustering_key = result
        state["rank_genes_df"] = df
        state["rank_genes_adata_norm"] = adata_norm
        state["rank_genes_groupby"] = clustering_key
        from xenium_viewer.utils.adata_persistence import save_rank_genes_to_adata
        save_rank_genes_to_adata(ctx, df, adata_norm, clustering_key)
        ga_run_button.enabled = True
        ga_dotplot_button.enabled = True
        ga_rank_plot_button.enabled = True
        ga_edit_labels_button.enabled = True
        ga_reset_labels_button.enabled = True
        ga_export_button.enabled = True
        ga_volcano_button.enabled = True
        llm_annotate_button.enabled = True
        preview = df.head(50).to_string(index=False)
        ga_results_text.setPlainText(preview)
        ga_status.value = f"Rank genes done: {len(df)} results ({clustering_key}, {ga_method_widget.value})"

        # Recording happened inside ctx.run_step(), which recorded the source it ran.

    def on_show_dotplot():
        adata_norm = state.get("rank_genes_adata_norm")
        groupby = state.get("rank_genes_groupby")
        if adata_norm is None or groupby is None:
            ga_status.value = "Run rank genes first"
            return
        ga_status.value = "Generating dotplot..."
        ga_dotplot_button.enabled = False

        n_genes = ga_dotplot_n_slider.value
        dendro = ga_dendro_check.value
        labels = ctx.get_labels_for(groupby)

        @thread_worker(connect={"returned": _on_dotplot_ready})
        def _run():
            ctx.apply_plot_font_size()
            fig = make_rank_genes_dotplot(
                adata_norm, groupby, n_genes=n_genes,
                cluster_labels=labels, dendrogram=dendro,
            )
            return fig
        _run()

    def _on_dotplot_ready(fig):
        state["dotplot_fig"] = fig
        ga_dotplot_button.enabled = True
        import matplotlib.pyplot as _plt
        _plt.show(block=False)
        path = ctx.auto_save_plot(fig, "dotplot")
        ga_status.value = f"Dotplot displayed — saved to {path}"

        _dp_n = ga_dotplot_n_slider.value
        _dp_dendro = ga_dendro_check.value
        _dp_groupby = state.get("rank_genes_groupby", "")
        _dp_fmt = ctx.state.get("plot_format", "svg")
        ctx.record_node(
            f"plot:dotplot:{_dp_groupby}",
            dotplot_code(_dp_groupby, _dp_n, _dp_dendro, _dp_fmt),
            deps=[f"rank_genes:{_dp_groupby}"],
            kind=TERMINAL,
            label="Rank-genes dotplot",
        )

    def _reset_labels():
        clustering_key = ga_clustering_widget.value
        all_labels = state.get("cluster_labels", {})
        if clustering_key in all_labels:
            del all_labels[clustering_key]
        # Drop the annotation node so it no longer appears in the notebook.
        graph = state.get("prov_graph")
        if graph is not None and f"annotation:{clustering_key}" in graph:
            try:
                graph.remove(f"annotation:{clustering_key}")
            except ValueError:
                pass
        ga_status.value = f"Labels reset for {clustering_key}"

    def _open_label_editor():
        clustering_key = ga_clustering_widget.value
        if ctx.build_label_editor_dialog(clustering_key):
            ga_status.value = f"Labels updated for {clustering_key}"

    def on_show_rank_plot():
        adata_norm = state.get("rank_genes_adata_norm")
        if adata_norm is None:
            ga_status.value = "Run rank genes first"
            return
        import matplotlib.pyplot as _plt
        ctx.apply_plot_font_size()
        _rp_n = ga_n_genes_slider.value
        groupby = state.get("rank_genes_groupby", "")
        labels = ctx.get_labels_for(groupby)
        fig = make_rank_genes_plot(adata_norm, n_genes=_rp_n, cluster_labels=labels)
        _plt.show(block=False)
        ga_status.value = "Rank genes plot displayed"
        ctx.record_node(
            f"plot:rank_panel:{groupby}",
            f"\n# Rank genes panel plot (n_genes={_rp_n})\n"
            f"sc.pl.rank_genes_groups(adata_norm, n_genes={_rp_n})",
            deps=[f"rank_genes:{groupby}"],
            kind=TERMINAL,
            label="Rank-genes panel plot",
        )

    def on_export_rank_genes():
        df = state.get("rank_genes_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Rank Genes Results", "rank_genes_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        ga_status.value = f"Exported {len(df)} rows to {path}"
        _rg_groupby = state.get("rank_genes_groupby", "")
        ctx.record_node(
            f"export:rank_genes:{_rg_groupby}",
            f"\n# Export rank genes results\n"
            f"rank_df.to_csv(\"{os.path.basename(path)}\", index=False)",
            deps=[f"rank_genes:{_rg_groupby}"],
            kind=TERMINAL,
            label="Export rank genes results",
        )

    def on_generate_volcanos():
        adata_norm = state.get("rank_genes_adata_norm")
        groupby = state.get("rank_genes_groupby")
        if adata_norm is None or groupby is None:
            ga_status.value = "Run rank genes first"
            return
        output_dir = QFileDialog.getExistingDirectory(None, "Select output directory for volcano plots")
        if not output_dir:
            return
        ga_volcano_button.enabled = False
        ga_status.value = "Generating volcano plots..."
        method = state.get("_rg_method", "wilcoxon")
        ctx.record_node(
            f"plot:volcano:{groupby}",
            f"\n# Pairwise volcano plots (groupby={groupby}, method={method})\n"
            f"import itertools\n"
            f"from pathlib import Path\n"
            f"from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"volcano_dir = Path(\"{os.path.basename(output_dir)}\"); "
            f"volcano_dir.mkdir(parents=True, exist_ok=True)\n"
            f"_groups = [g for g in adata_norm.obs[\"{groupby}\"].cat.categories if str(g) != '-1']\n"
            f"for _a, _b in itertools.combinations(_groups, 2):\n"
            f"    _df = run_pairwise_deg(adata_norm, \"{groupby}\", str(_a), str(_b), method=\"{method}\")\n"
            f"    _vfig = make_volcano_plot(_df, str(_a), str(_b))\n"
            f"    _vfig.savefig(volcano_dir / f'volcano_{{_a}}_vs_{{_b}}.png', dpi=300)\n"
            f"    plt.close(_vfig)",
            deps=[f"rank_genes:{groupby}"],
            kind=TERMINAL,
            label="Pairwise volcano plots",
        )

        @thread_worker(connect={"yielded": lambda msg: setattr(ga_status, 'value', msg),
                                "returned": _on_volcanos_done})
        def _run():
            from pathlib import Path
            import itertools as _it
            from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot
            import matplotlib.pyplot as _plt

            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            groups = sorted(
                [g for g in adata_norm.obs[groupby].cat.categories if str(g) != '-1'],
                key=lambda x: (int(x) if str(x).lstrip('-').isdigit() else 0, str(x)),
            )
            pairs = list(_it.combinations(groups, 2))
            total = len(pairs)
            for i, (a, b) in enumerate(pairs):
                yield f"Volcano plot {i + 1}/{total}: {a} vs {b}"
                df = run_pairwise_deg(adata_norm, groupby, str(a), str(b), method=method)
                fig = make_volcano_plot(df, str(a), str(b))
                fig.savefig(out / f'volcano_{a}_vs_{b}.png', dpi=300)
                _plt.close(fig)
            return total, output_dir
        _run()

    def _on_volcanos_done(result):
        count, out_dir = result
        ga_volcano_button.enabled = True
        ga_status.value = f"{count} volcano plots saved to {out_dir}"

    # -- CellTypist handlers --
    def on_download_models():
        if not _HAS_CELLTYPIST:
            ga_status.value = "CellTypist is not installed (pip install celltypist)"
            return
        ga_ct_download_button.enabled = False
        ga_status.value = "Downloading CellTypist models..."

        @thread_worker
        def _download():
            ct_models.download_models()
            return sorted(m.name for m in Path(ct_models.models_path).glob("*.pkl"))

        def _on_download_done(model_names):
            ga_ct_download_button.enabled = True
            if model_names:
                ga_ct_model_widget.choices = model_names
                ga_ct_annotate_button.enabled = True
                ga_status.value = f"Downloaded {len(model_names)} CellTypist models"
            else:
                ga_status.value = "No CellTypist models found after download"

        worker = _download()
        worker.returned.connect(_on_download_done)
        worker.start()

    def on_celltypist_annotate():
        if not _HAS_CELLTYPIST:
            ga_status.value = "CellTypist is not installed"
            return
        model_name = ga_ct_model_widget.value
        if not model_name:
            ga_status.value = "Select a CellTypist model first"
            return
        clustering_key = ga_clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            ga_status.value = "Select a valid clustering first"
            return

        ga_ct_annotate_button.enabled = False
        conf_threshold = ga_ct_conf_slider.value
        ga_status.value = f"Running CellTypist ({model_name}, conf≥{conf_threshold:.2f})..."

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        _clustering_series = ctx.clusterings[clustering_key]

        @thread_worker
        def _run():
            # The same adata_norm every other analysis uses, rather than a second
            # normalised copy from get_normalized_adata's id()-keyed cache. The
            # CellTypist call itself is not a Step: what gets recorded is the
            # resolved cluster→label mapping, not the model run (see below).
            adata_norm = ctx.ensure_normalized()
            cell_predictions, cell_confidence = run_celltypist_annotation(adata_norm, model_name)
            return cell_predictions, cell_confidence, clustering_key

        def _on_celltypist_ready(result):
            cell_predictions, cell_confidence, clust_key = result
            _clustering_series_local = ctx.clusterings[clust_key]
            _orig_adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
            ga_ct_annotate_button.enabled = True

            # Build per-obs_name cluster assignment (same alignment as add_clustering_to_obs)
            if 'cell_id' in _orig_adata.obs.columns:
                cell_ids = _orig_adata.obs['cell_id'].values
                aligned_clusters = _clustering_series_local.reindex(cell_ids)
            else:
                aligned_clusters = _clustering_series_local.reindex(_orig_adata.obs_names)
            # Index by obs_names so it matches cell_predictions index
            aligned_clusters.index = _orig_adata.obs_names

            # Filter out low-confidence predictions
            high_conf_mask = cell_confidence >= conf_threshold
            filtered_predictions = cell_predictions[high_conf_mask]
            n_total = len(cell_predictions)
            n_passed = len(filtered_predictions)

            # Majority vote: map per-cell predictions to per-cluster labels
            labels = {}
            for cluster_id in aligned_clusters.dropna().unique():
                mask = aligned_clusters == cluster_id
                obs_in_cluster = aligned_clusters.index[mask]
                preds = filtered_predictions.reindex(obs_in_cluster).dropna()
                if len(preds) == 0:
                    labels[cluster_id] = "Unknown"
                else:
                    labels[cluster_id] = preds.value_counts().idxmax()

            # Store in the same system as manual labels
            if "cluster_labels" not in state:
                state["cluster_labels"] = {}
            state["cluster_labels"][clust_key] = labels
            from xenium_viewer.utils.adata_persistence import save_cluster_labels_to_sdata
            save_cluster_labels_to_sdata(ctx, clust_key, labels)

            # Summary
            label_counts = {}
            for lbl in labels.values():
                label_counts[lbl] = label_counts.get(lbl, 0) + 1
            summary_parts = [f"{lbl} ({n})" for lbl, n in sorted(label_counts.items(), key=lambda x: -x[1])]
            summary = ", ".join(summary_parts)
            pct_passed = 100 * n_passed / n_total if n_total > 0 else 0
            ga_status.value = (
                f"{len(labels)} clusters annotated ({n_passed}/{n_total} cells "
                f"passed conf≥{conf_threshold:.2f}, {pct_passed:.0f}%): {summary}"
            )

            # Code recording — emit the resolved cluster→label mapping as real,
            # runnable code (deterministic), with the CellTypist derivation as a
            # provenance comment above it.
            ctx.record_clustering(clust_key)
            _ct_map = {str(k): v for k, v in labels.items()}
            ctx.record_node(
                f"annotation:{clust_key}",
                f"\n# Cell-type annotation of '{clust_key}' (CellTypist, "
                f"confidence >= {conf_threshold})\n"
                f"# Derived with: celltypist.annotate(adata, "
                f"model=celltypist.models.Model.load('{model_name}'), "
                f"majority_voting=False),\n"
                f"# then majority-voted per cluster among cells passing the threshold.\n"
                f"annotation_map = {_ct_map!r}\n"
                f"adata.obs[\"{clust_key}_annotated\"] = ("
                f"adata.obs[\"{clust_key}\"].astype(str).map(annotation_map)"
                f".astype(\"category\"))",
                deps=[f"clustering:{clust_key}"],
                label=f"Cell-type annotation: {clust_key}",
            )

        worker = _run()
        timer, _ = attach_spinner(worker, lambda m: setattr(ga_status, 'value', m), "Running CellTypist...")
        state['_ct_spinner_timer'] = timer
        worker.returned.connect(_on_celltypist_ready)
        worker.start()

    ga_ct_download_button.clicked.connect(on_download_models)
    ga_ct_annotate_button.clicked.connect(on_celltypist_annotate)

    # -- LLM Annotation widgets --
    llm_separator = QLabel("── LLM Annotation ──")
    llm_separator.setStyleSheet("font-weight: bold; margin-top: 8px;")
    llm_provider_widget = ComboBox(
        label="LLM Provider",
        choices=["Claude (claude)", "Gemini (gemini)", "Codex (codex)"],
        value="Claude (claude)",
    )
    llm_annotate_button = PushButton(label="Annotate with LLM", enabled=False)

    def on_llm_annotate():
        rank_df = state.get("rank_genes_df")
        if rank_df is None or rank_df.empty:
            ga_status.value = "Run Rank Genes first"
            return
        clustering_key = state.get("rank_genes_groupby")
        if not clustering_key:
            ga_status.value = "Run Rank Genes first"
            return

        # Parse CLI name from combo choice
        provider_str = llm_provider_widget.value
        cli = provider_str.split("(")[-1].rstrip(")")

        llm_annotate_button.enabled = False
        ga_status.value = f"Running LLM annotation ({cli})..."

        @thread_worker
        def _run():
            return run_llm_annotation(rank_df, cli, n_genes=10)

        def _on_llm_ready(labels):
            llm_annotate_button.enabled = True
            if "cluster_labels" not in state:
                state["cluster_labels"] = {}
            state["cluster_labels"][clustering_key] = labels
            from xenium_viewer.utils.adata_persistence import save_cluster_labels_to_sdata
            save_cluster_labels_to_sdata(ctx, clustering_key, labels)

            label_counts = {}
            for lbl in labels.values():
                label_counts[lbl] = label_counts.get(lbl, 0) + 1
            summary_parts = [f"{lbl} ({n})" for lbl, n in sorted(label_counts.items(), key=lambda x: -x[1])]
            summary = ", ".join(summary_parts)
            ga_status.value = f"{len(labels)} clusters annotated via {cli}: {summary}"

            ctx.record_clustering(clustering_key)
            _llm_map = {str(k): v for k, v in labels.items()}
            ctx.record_node(
                f"annotation:{clustering_key}",
                f"\n# Cell-type annotation of '{clustering_key}' (LLM via {cli})\n"
                f"# Derived with: run_llm_annotation(rank_df, cli='{cli}', n_genes=10)\n"
                f"annotation_map = {_llm_map!r}\n"
                f"adata.obs[\"{clustering_key}_annotated\"] = ("
                f"adata.obs[\"{clustering_key}\"].astype(str).map(annotation_map)"
                f".astype(\"category\"))",
                deps=[f"clustering:{clustering_key}"],
                label=f"Cell-type annotation: {clustering_key}",
            )

        def _on_llm_error(exc):
            llm_annotate_button.enabled = True
            ga_status.value = f"LLM annotation failed: {exc}"

        worker = _run()
        timer, _ = attach_spinner(worker, lambda m: setattr(ga_status, 'value', m), f"Running {cli}...")
        state['_llm_spinner_timer'] = timer
        worker.returned.connect(_on_llm_ready)
        worker.errored.connect(_on_llm_error)
        worker.start()

    llm_annotate_button.clicked.connect(on_llm_annotate)

    # -- Label Transfer (sc.tl.ingest) widgets --
    lt_separator = QLabel("── Label Transfer (sc.tl.ingest) ──")
    lt_separator.setStyleSheet("font-weight: bold; margin-top: 8px;")
    lt_ref_widget = ComboBox(label="Reference Dataset", choices=list(_REFERENCE_DATASETS.keys()))
    lt_col_widget = ComboBox(label="Annotation Column", choices=[])
    lt_load_ref_button = PushButton(label="Load Reference")
    lt_transfer_button = PushButton(label="Run Label Transfer", enabled=False)

    state["_lt_custom_paths"] = {}
    state["_lt_ref_cache"] = {}

    def on_lt_ref_changed(value):
        if value == "Browse...":
            path, _ = QFileDialog.getOpenFileName(
                None, "Select Reference h5ad", "", "AnnData Files (*.h5ad)",
            )
            if path:
                name = Path(path).stem
                display = f"{name} (custom)"
                state["_lt_custom_paths"][display] = path
                # Add to choices before Browse...
                choices = list(lt_ref_widget.choices)
                if display not in choices:
                    choices.insert(len(choices) - 1, display)
                    lt_ref_widget.choices = choices
                lt_ref_widget.value = display
            else:
                # User cancelled — revert to first choice
                lt_ref_widget.value = list(_REFERENCE_DATASETS.keys())[0]
        elif value in _REFERENCE_DATASETS:
            meta = _REFERENCE_DATASETS[value].get("metadata")
            if meta:
                authors = meta.get("authors", "")
                journal = meta.get("journal", "")
                year = meta.get("year", "")
                platform = meta.get("platform", "")
                tissue = meta.get("tissue", "")
                ga_status.value = f"Paper: {authors} ({journal}, {year}) | Platform: {platform} | Tissue: {tissue}"

    lt_ref_widget.changed.connect(on_lt_ref_changed)

    def _resolve_ref_path():
        """Resolve the path for the currently selected reference dataset."""
        value = lt_ref_widget.value
        if value in _REFERENCE_DATASETS:
            return _REFERENCE_DATASETS[value]["path"]
        return state["_lt_custom_paths"].get(value)

    def _get_default_col():
        value = lt_ref_widget.value
        if value in _REFERENCE_DATASETS:
            return _REFERENCE_DATASETS[value]["default_col"]
        return None

    def on_load_reference():
        path = _resolve_ref_path()
        if not path:
            ga_status.value = "No reference dataset selected"
            return
        lt_load_ref_button.enabled = False
        ga_status.value = "Loading reference dataset..."

        cache_key = path

        @thread_worker
        def _run():
            if cache_key in state["_lt_ref_cache"]:
                return state["_lt_ref_cache"][cache_key]
            ref = load_reference_h5ad(path)
            state["_lt_ref_cache"][cache_key] = ref
            return ref

        def _on_ref_loaded(ref):
            lt_load_ref_button.enabled = True
            cols = get_annotation_columns(ref)
            if not cols:
                ga_status.value = "No suitable annotation columns found in reference"
                return
            lt_col_widget.choices = cols
            default = _get_default_col()
            if default and default in cols:
                lt_col_widget.value = default
            lt_transfer_button.enabled = True
            ga_status.value = (
                f"Loaded: {ref.n_obs} cells × {ref.n_vars} genes, "
                f"{len(cols)} annotation columns"
            )

        worker = _run()
        timer, _ = attach_spinner(worker, lambda m: setattr(ga_status, 'value', m), "Loading reference...")
        state['_lt_spinner_timer'] = timer
        worker.returned.connect(_on_ref_loaded)
        worker.start()

    def on_label_transfer():
        clustering_key = ga_clustering_widget.value
        if not clustering_key or clustering_key not in ctx.clusterings:
            ga_status.value = "Select a valid clustering first"
            return
        ref_path = _resolve_ref_path()
        if not ref_path or ref_path not in state["_lt_ref_cache"]:
            ga_status.value = "Load a reference dataset first"
            return
        annotation_col = lt_col_widget.value
        if not annotation_col:
            ga_status.value = "Select an annotation column"
            return

        lt_transfer_button.enabled = False
        ga_status.value = f"Running label transfer ({annotation_col})..."

        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        ref_adata = state["_lt_ref_cache"][ref_path]
        _clustering_series = ctx.clusterings[clustering_key]

        @thread_worker
        def _run():
            predictions, n_common = run_label_transfer(_adata, ref_adata, annotation_col)
            return predictions, n_common, clustering_key

        def _on_lt_ready(result):
            predictions, n_common, clust_key = result
            _clustering_series_local = ctx.clusterings[clust_key]
            _orig_adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
            lt_transfer_button.enabled = True

            # Align clustering to obs_names (same pattern as CellTypist)
            if 'cell_id' in _orig_adata.obs.columns:
                cell_ids = _orig_adata.obs['cell_id'].values
                aligned_clusters = _clustering_series_local.reindex(cell_ids)
            else:
                aligned_clusters = _clustering_series_local.reindex(_orig_adata.obs_names)
            aligned_clusters.index = _orig_adata.obs_names

            # Majority vote per cluster (no confidence filtering for ingest)
            labels = {}
            for cluster_id in aligned_clusters.dropna().unique():
                mask = aligned_clusters == cluster_id
                obs_in_cluster = aligned_clusters.index[mask]
                preds = predictions.reindex(obs_in_cluster).dropna()
                if len(preds) == 0:
                    labels[cluster_id] = "Unknown"
                else:
                    labels[cluster_id] = preds.value_counts().idxmax()

            # Store labels
            if "cluster_labels" not in state:
                state["cluster_labels"] = {}
            state["cluster_labels"][clust_key] = labels
            from xenium_viewer.utils.adata_persistence import save_cluster_labels_to_sdata
            save_cluster_labels_to_sdata(ctx, clust_key, labels)

            # Summary
            label_counts = {}
            for lbl in labels.values():
                label_counts[lbl] = label_counts.get(lbl, 0) + 1
            summary_parts = [f"{lbl} ({n})" for lbl, n in sorted(label_counts.items(), key=lambda x: -x[1])]
            summary = ", ".join(summary_parts)
            ga_status.value = (
                f"{len(labels)} clusters annotated ({n_common} common genes): {summary}"
            )

            # Code recording — deterministic mapping application, with the ingest
            # derivation kept as a provenance comment (ingest/UMAP is stochastic).
            ctx.record_clustering(clust_key)
            _lt_map = {str(k): v for k, v in labels.items()}
            ctx.record_node(
                f"annotation:{clust_key}",
                f"\n# Cell-type annotation of '{clust_key}' (label transfer via sc.tl.ingest)\n"
                f"# Reference: {ref_path} | column: {annotation_col} | common genes: {n_common}\n"
                f"# Derived by ingesting the Xenium cells into the reference embedding:\n"
                f"#   ref = sc.read_h5ad('{ref_path}')\n"
                f"#   genes = sorted(set(adata.var_names) & set(ref.var_names))\n"
                f"#   ref, xen = ref[:, genes].copy(), adata[:, genes].copy()\n"
                f"#   sc.pp.normalize_total(ref, target_sum=1e4); sc.pp.log1p(ref)\n"
                f"#   sc.pp.highly_variable_genes(ref); sc.pp.pca(ref); "
                f"sc.pp.neighbors(ref); sc.tl.umap(ref)\n"
                f"#   sc.pp.normalize_total(xen, target_sum=1e4); sc.pp.log1p(xen)\n"
                f"#   sc.tl.ingest(xen, ref, obs='{annotation_col}')  # majority-vote per cluster\n"
                f"annotation_map = {_lt_map!r}\n"
                f"adata.obs[\"{clust_key}_annotated\"] = ("
                f"adata.obs[\"{clust_key}\"].astype(str).map(annotation_map)"
                f".astype(\"category\"))",
                deps=[f"clustering:{clust_key}"],
                label=f"Cell-type annotation: {clust_key}",
            )

        worker = _run()
        timer, _ = attach_spinner(worker, lambda m: setattr(ga_status, 'value', m), "Running label transfer...")
        state['_lt_transfer_timer'] = timer
        worker.returned.connect(_on_lt_ready)
        worker.start()

    lt_load_ref_button.clicked.connect(on_load_reference)
    lt_transfer_button.clicked.connect(on_label_transfer)

    ga_run_button.clicked.connect(on_run_rank_genes)
    ga_dotplot_button.clicked.connect(on_show_dotplot)
    ga_edit_labels_button.clicked.connect(_open_label_editor)
    ga_reset_labels_button.clicked.connect(_reset_labels)
    ga_rank_plot_button.clicked.connect(on_show_rank_plot)
    ga_export_button.clicked.connect(on_export_rank_genes)
    ga_volcano_button.clicked.connect(on_generate_volcanos)

    # Layout
    ga_dotplot_btn_row = QWidget()
    ga_dotplot_btn_layout = QHBoxLayout()
    ga_dotplot_btn_layout.setContentsMargins(0, 0, 0, 0)
    ga_dotplot_btn_layout.addWidget(ga_dotplot_button.native)
    ga_dotplot_btn_layout.addWidget(ga_edit_labels_button.native)
    ga_dotplot_btn_layout.addWidget(ga_reset_labels_button.native)
    ga_dotplot_btn_row.setLayout(ga_dotplot_btn_layout)

    # CellTypist button row
    ga_ct_btn_row = QWidget()
    ga_ct_btn_layout = QHBoxLayout()
    ga_ct_btn_layout.setContentsMargins(0, 0, 0, 0)
    ga_ct_btn_layout.addWidget(ga_ct_download_button.native)
    ga_ct_btn_layout.addWidget(ga_ct_annotate_button.native)
    ga_ct_btn_row.setLayout(ga_ct_btn_layout)

    # Label Transfer button row
    lt_btn_row = QWidget()
    lt_btn_layout = QHBoxLayout()
    lt_btn_layout.setContentsMargins(0, 0, 0, 0)
    lt_btn_layout.addWidget(lt_load_ref_button.native)
    lt_btn_layout.addWidget(lt_transfer_button.native)
    lt_btn_row.setLayout(lt_btn_layout)

    widget = make_tab(
        ga_clustering_widget,
        ga_method_widget,
        ga_n_genes_slider,
        ga_run_button,
        ga_progress,
        ga_dotplot_n_slider,
        ga_dendro_check,
        ga_dotplot_btn_row,
        ga_rank_plot_button,
        ga_export_button,
        ga_volcano_button,
        ga_ct_separator,
        ga_ct_model_widget,
        ga_ct_conf_slider,
        ga_ct_btn_row,
        llm_separator,
        llm_provider_widget,
        llm_annotate_button,
        lt_separator,
        lt_ref_widget,
        lt_col_widget,
        lt_btn_row,
    )

    def _restore_session(session):
        # Check state first (populated from adata.uns at startup), fall back to session
        rg = state.get("rank_genes_df") if state.get("rank_genes_df") is not None else session.get("rank_genes_df")
        if rg is not None:
            state["rank_genes_df"] = rg
            rg_norm = state.get("rank_genes_adata_norm") or session.get("rank_genes_adata_norm")
            state["rank_genes_adata_norm"] = rg_norm
            state["rank_genes_groupby"] = (
                state.get("rank_genes_groupby") or session.get("rank_genes_groupby")
            )
            ga_dotplot_button.enabled = rg_norm is not None
            ga_rank_plot_button.enabled = True
            ga_edit_labels_button.enabled = True
            ga_reset_labels_button.enabled = True
            ga_export_button.enabled = True
            ga_volcano_button.enabled = rg_norm is not None
            llm_annotate_button.enabled = True
            preview = rg.head(50).to_string(index=False)
            ga_results_text.setPlainText(preview)
            print(f"  Restored rank genes ({len(rg)} rows)")

    return widget, {"ga_clustering_widget": ga_clustering_widget, "restore_session": _restore_session}
