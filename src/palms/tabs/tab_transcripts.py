"""Tab 2: Transcript overlay — multi-gene points display, and density binning.

The **points overlay** is display: it reads the per-gene feather index
(`utils/transcript_index.py`, built by `palms-preprocess`) because drawing a
gene has to feel instant, and nothing about which points are on screen belongs
in a notebook.

The **density heatmap** is analysis, and is recorded. It reads
`sdata.points['transcripts']` instead of that index — the index is a viewer
artifact a replayed notebook would not have — which costs a few seconds the
first time a gene is binned and nothing thereafter, since the fetch is its own
step and the bin-size knob only re-runs the histogram.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import QListWidget, QHBoxLayout, QWidget, QLabel
from napari.qt.threading import thread_worker
from palms.tabs._helpers import make_tab, StatusProxy
from palms.utils.prov_graph import ARTIFACT, NOTE
from palms.utils.steps import Step, coerce
from palms.utils.step_templates import (
    Preview, builtin_spec, step_template as _resolved,
)

GENE_TEMPLATE_ID = "transcripts.gene"
DENSITY_TEMPLATE_ID = "transcripts.density"

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    transcript_gene_widget = ComboBox(
        label="Transcript gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    add_gene_button = PushButton(label="Add Gene", enabled=True)
    remove_gene_button = PushButton(label="Remove Selected", enabled=True)
    clear_genes_button = PushButton(label="Clear All", enabled=True)

    gene_list_qt = QListWidget()
    gene_list_qt.setMaximumHeight(150)

    legend_label_qt = QLabel("")
    legend_label_qt.setWordWrap(True)

    transcript_check = CheckBox(label="Show transcripts", value=False)
    qv_slider = Slider(label="Min QV", min=0, max=40, value=20)
    apply_transcripts_button = PushButton(label="Apply Transcripts", enabled=True)

    status_label = StatusProxy(ctx.viewer)

    def _update_legend():
        genes = state["transcript_genes"]
        if not genes:
            legend_label_qt.setText("")
            return
        color_names = [
            "Yellow", "Cyan", "Magenta", "Orange", "Green",
            "Sky Blue", "Red", "Violet", "Pink", "Brown",
        ]
        parts = []
        for i, g in enumerate(genes):
            parts.append(f"{color_names[i % len(color_names)]}: {g}")
        legend_label_qt.setText(" | ".join(parts))

    def on_add_gene():
        gene = transcript_gene_widget.value
        genes = state["transcript_genes"]
        if gene in genes:
            status_label.value = f"'{gene}' already in list"
            return
        if len(genes) >= 10:
            status_label.value = "Maximum 10 genes reached"
            return
        genes.append(gene)
        gene_list_qt.addItem(gene)
        _update_legend()

    def on_remove_gene():
        selected = gene_list_qt.currentRow()
        if selected >= 0:
            gene_list_qt.takeItem(selected)
            state["transcript_genes"].pop(selected)
            _update_legend()

    def on_clear_genes():
        state["transcript_genes"].clear()
        gene_list_qt.clear()
        _update_legend()

    def _on_transcripts_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        points, colors = result
        ctx.transcript_layer.data = points
        ctx.transcript_layer.face_color = colors
        ctx.transcript_layer.visible = True
        genes = state["transcript_genes"]
        status_label.value = (
            f"Transcripts: {', '.join(genes)} ({len(points):,} spots total)"
        )
        ctx.record_node(
            "viewer:transcript_overlay",
            f"\n# Display transcript overlay (viewer)\n"
            f"# genes={genes}, qv_threshold={qv_slider.value}, {len(points):,} spots total",
            deps=["preamble"], kind=NOTE, label="Transcript overlay",
        )
        apply_transcripts_button.enabled = True

    def on_apply_transcripts():
        genes = state["transcript_genes"]
        if transcript_check.value and genes:
            status_label.value = f"Loading transcripts for {len(genes)} gene(s)..."
            apply_transcripts_button.enabled = False
            gen = ctx.dataset_generation

            @thread_worker
            def fetch():
                return ctx.transcript_loader.get_multi_gene_points(genes)

            worker = fetch()
            worker.returned.connect(lambda result: _on_transcripts_ready(result, gen))
            worker.start()
        elif transcript_check.value and not genes:
            status_label.value = "No genes in list — add genes first"
        else:
            ctx.transcript_layer.visible = False
            status_label.value = "Transcripts hidden"

    add_gene_button.clicked.connect(on_add_gene)
    remove_gene_button.clicked.connect(on_remove_gene)
    clear_genes_button.clicked.connect(on_clear_genes)
    apply_transcripts_button.clicked.connect(on_apply_transcripts)

    # Button row
    btn_row = QWidget()
    btn_layout = QHBoxLayout()
    btn_layout.setContentsMargins(0, 0, 0, 0)
    btn_layout.addWidget(add_gene_button.native)
    btn_layout.addWidget(remove_gene_button.native)
    btn_layout.addWidget(clear_genes_button.native)
    btn_row.setLayout(btn_layout)

    # ── Transcript Density section ────────────────────────────────────────────
    density_separator = QLabel("── Transcript Density ──")

    density_gene_widget = ComboBox(
        label="Gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    bin_size_slider = Slider(label="Bin size (µm)", min=10, max=500, value=50)
    cluster_filter_density_check = CheckBox(label="Filter by selected clusters", value=False)
    normalise_cells_check = CheckBox(label="Normalise by cells per bin", value=False)
    compute_density_button = PushButton(label="Compute Density")
    density_status = QLabel("")

    def _cluster_filter_params() -> dict:
        """The two filter params, or nothing at all.

        Read-only, like tab_roi's namesake: recording the clustering node and
        mirroring it onto obs are the run's job, not a side effect of drawing a
        preview pane.
        """
        if not cluster_filter_density_check.value:
            return {}
        key = (ctx.state.get("active_clustering_name")
               or getattr(ctx.clustering_widget, "value", None))
        if not key or key not in (ctx.clusterings or {}):
            return {}
        return {
            "clustering": key,
            "selected": sorted({str(i) for i in ctx.get_selected_cluster_ids()}),
        }

    # ── Preview providers ─────────────────────────────────────────────────────

    def _transcripts_gene_preview() -> Preview:
        return Preview(
            list(builtin_spec(GENE_TEMPLATE_ID).blocks),
            {"gene": density_gene_widget.value,
             # The QV slider is provenance-only for the point overlay, whose
             # feather files were filtered when palms-preprocess built them.
             # Here it is a real filter, applied as the transcripts are read.
             "min_qv": int(coerce(qv_slider.value))},
        )

    def _density_blocks(filtered: bool, normalised: bool) -> list[str]:
        """Both switches are blocks as well as params: "filter by cluster" adds
        a mask, "normalise" adds a second histogram."""
        blocks = ["head"]
        if filtered:
            blocks.append("filter")
        blocks.append("main")
        if normalised:
            blocks.append("normalise")
        return blocks

    def _density_preview() -> Preview:
        params = {
            "gene": density_gene_widget.value,
            "bin_size_um": float(coerce(bin_size_slider.value)),
            "pixel_size": float(coerce(ctx.pixel_size)),
        }
        filter_params = _cluster_filter_params()
        params.update(filter_params)
        return Preview(
            _density_blocks(bool(filter_params), normalise_cells_check.value),
            params,
        )

    ctx.state.setdefault("template_preview", {})[GENE_TEMPLATE_ID] = _transcripts_gene_preview
    ctx.state.setdefault("template_preview", {})[DENSITY_TEMPLATE_ID] = _density_preview

    def _on_density_ready(result):
        density, bin_px, normalised, needs_clustering = result
        compute_density_button.enabled = True
        if density is None:
            return
        ctx.transcript_bins_layer.data = density.astype(np.float32)
        ctx.transcript_bins_layer.scale = [bin_px, bin_px]
        nonzero = density[density > 0]
        if nonzero.size > 0:
            ctx.transcript_bins_layer.contrast_limits = [0, float(np.percentile(nonzero, 99))]
        else:
            ctx.transcript_bins_layer.contrast_limits = [0, 1]
        ctx.transcript_bins_layer.visible = True
        if normalised:
            density_status.setText(
                f"Done — {int(nonzero.size):,} non-zero bins, normalised by cells")
        else:
            density_status.setText(f"Done — {int(density.sum()):,} transcripts binned")

    def on_compute_density():
        gene = density_gene_widget.value
        if gene is None:
            return
        if "transcripts" not in getattr(ctx.sdata, "points", {}):
            density_status.setText("No transcripts element — cannot compute density")
            return
        if "morphology_focus" not in getattr(ctx.sdata, "images", {}):
            density_status.setText("No morphology image — cannot compute density")
            return

        gene_blocks, gene_params, _ = _transcripts_gene_preview()
        blocks, params, _ = _density_preview()
        if cluster_filter_density_check.value and "clustering" not in params:
            density_status.setText("No clustering applied — filter skipped")
            return

        # The clustering must exist as a node, and in adata.obs, before the step
        # can declare it as a dependency and read it. Both on the GUI thread.
        deps = [f"transcripts:{gene}"]
        clustering_key = params.get("clustering")
        if clustering_key:
            from palms.utils.gene_analysis import add_clustering_to_obs
            ctx.record_clustering(clustering_key)
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps.append(f"clustering:{clustering_key}")

        gene_step = Step(
            id=f"transcripts:{gene}",
            **_resolved(GENE_TEMPLATE_ID, gene_blocks),
            params=gene_params,
            deps=["preamble"],
            kind=ARTIFACT,
            label=f"Transcripts: {gene}",
            outputs=["transcript_points"],
        )
        density_step = Step(
            id=f"transcript_density:{gene}",
            **_resolved(DENSITY_TEMPLATE_ID, blocks),
            params=params,
            deps=deps,
            kind=ARTIFACT,
            label=f"Transcript density: {gene}",
            outputs=["transcript_density"],
        )

        # Fetching a gene's transcripts is the slow half — seconds, against
        # milliseconds for the histogram — so it is its own step and runs only
        # when the gene changes. Sliding the bin size re-runs the cheap half.
        fetch_key = (gene, gene_params["min_qv"], ctx.dataset_generation)
        need_fetch = state.get("_transcript_points_key") != fetch_key

        compute_density_button.enabled = False
        density_status.setText(
            f"Fetching transcripts for {gene}..." if need_fetch
            else f"Computing density for {gene}...")
        normalised = normalise_cells_check.value
        bin_px = params["bin_size_um"] / params["pixel_size"]

        @thread_worker
        def _run():
            ctx.record_preamble()
            if need_fetch:
                ctx.run_step(gene_step)
                state["_transcript_points_key"] = fetch_key
            density = ctx.run_step(density_step)["transcript_density"]
            return density, bin_px, normalised, clustering_key

        def _failed(exc):
            compute_density_button.enabled = True
            state.pop("_transcript_points_key", None)
            density_status.setText(f"Density failed: {exc}")

        worker = _run()
        worker.returned.connect(_on_density_ready)
        worker.errored.connect(_failed)
        worker.start()

    compute_density_button.clicked.connect(on_compute_density)

    widget = make_tab(
        transcript_gene_widget,
        transcript_check,
        qv_slider,
        apply_transcripts_button,
        gene_list_qt,
        btn_row,
        legend_label_qt,
        density_separator,
        density_gene_widget,
        bin_size_slider,
        cluster_filter_density_check,
        normalise_cells_check,
        compute_density_button,
        density_status,
    )
    return widget, {}
