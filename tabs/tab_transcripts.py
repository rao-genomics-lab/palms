"""Tab 2: Transcript overlay — multi-gene points display."""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import QListWidget, QHBoxLayout, QWidget, QLabel
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


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
        ctx.record_code(
            f"\n# Display transcript overlay\n"
            f"# genes={genes}, {len(points):,} spots total"
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
        label="Density gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    bin_size_slider = Slider(label="Bin size (µm)", min=10, max=500, value=50)
    cluster_filter_density_check = CheckBox(label="Filter by selected clusters", value=False)
    normalise_cells_check = CheckBox(label="Normalise by cells per bin", value=False)
    compute_density_button = PushButton(label="Compute Density")
    density_status = QLabel("")

    def _compute_density_worker(gene, bin_size_um, use_cluster_filter, normalise):
        pts = ctx.transcript_loader.get_points_array(gene)  # (N, 2) [y, x] pixels
        H, W = ctx.morph_full_shape_yx
        bin_px = bin_size_um / ctx.pixel_size
        n_bins_y = max(1, int(H / bin_px))
        n_bins_x = max(1, int(W / bin_px))

        active_centroids = None
        if use_cluster_filter and ctx.state.get("label_to_cluster") is not None:
            selected_ids = ctx.translate_selected_ids_to_int(ctx.get_selected_cluster_ids())
            if selected_ids:
                selected_set = set(selected_ids)
                label_to_cluster = ctx.state["label_to_cluster"]
                all_labels = np.arange(len(label_to_cluster))
                selected_labels = all_labels[np.isin(label_to_cluster, list(selected_set))]
                obs_indices = ctx.label_to_obs[selected_labels]
                valid = obs_indices >= 0
                active_centroids = ctx.centroids_yx[obs_indices[valid]]  # (K, 2) [y, x]
                cy_bin = np.clip((active_centroids[:, 0] / bin_px).astype(int), 0, n_bins_y - 1)
                cx_bin = np.clip((active_centroids[:, 1] / bin_px).astype(int), 0, n_bins_x - 1)
                coverage = np.zeros((n_bins_y, n_bins_x), dtype=bool)
                coverage[cy_bin, cx_bin] = True
                ty_bin = np.clip((pts[:, 0] / bin_px).astype(int), 0, n_bins_y - 1)
                tx_bin = np.clip((pts[:, 1] / bin_px).astype(int), 0, n_bins_x - 1)
                mask = coverage[ty_bin, tx_bin]
                pts = pts[mask]

        hist, _, _ = np.histogram2d(
            pts[:, 0], pts[:, 1],
            bins=[n_bins_y, n_bins_x],
            range=[[0, H], [0, W]],
        )

        if normalise:
            centroids_for_count = active_centroids if active_centroids is not None else ctx.centroids_yx
            cell_count, _, _ = np.histogram2d(
                centroids_for_count[:, 0], centroids_for_count[:, 1],
                bins=[n_bins_y, n_bins_x],
                range=[[0, H], [0, W]],
            )
            hist = np.where(cell_count > 0, hist / cell_count, 0).astype(np.float32)

        return hist.astype(np.float32), bin_px, normalise

    def _on_density_ready(result):
        hist, bin_px, normalised = result
        ctx.transcript_bins_layer.data = hist
        ctx.transcript_bins_layer.scale = [bin_px, bin_px]
        nonzero = hist[hist > 0]
        if nonzero.size > 0:
            ctx.transcript_bins_layer.contrast_limits = [0, float(np.percentile(nonzero, 99))]
        else:
            ctx.transcript_bins_layer.contrast_limits = [0, 1]
        ctx.transcript_bins_layer.visible = True
        if normalised:
            density_status.setText(f"Done — {int(hist[hist > 0].size):,} non-zero bins, normalised by cells")
        else:
            density_status.setText(f"Done — {int(hist.sum()):,} transcripts binned")
        compute_density_button.enabled = True

    def on_compute_density():
        gene = density_gene_widget.value
        if gene is None:
            return
        if ctx.morph_full_shape_yx is None:
            density_status.setText("No morphology data — cannot compute density")
            return
        compute_density_button.enabled = False
        density_status.setText(f"Computing density for {gene}...")
        bin_um = bin_size_slider.value
        use_filter = cluster_filter_density_check.value
        normalise = normalise_cells_check.value

        @thread_worker
        def _run():
            return _compute_density_worker(gene, bin_um, use_filter, normalise)

        worker = _run()
        worker.returned.connect(_on_density_ready)
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
