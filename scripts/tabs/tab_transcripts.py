"""Tab 2: Transcript overlay — multi-gene points display."""

from __future__ import annotations
from typing import TYPE_CHECKING

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

    def _on_transcripts_ready(result):
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

            @thread_worker(connect={"returned": _on_transcripts_ready})
            def fetch():
                return ctx.transcript_loader.get_multi_gene_points(genes)
            fetch()
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

    widget = make_tab(
        transcript_gene_widget,
        transcript_check,
        qv_slider,
        apply_transcripts_button,
        gene_list_qt,
        btn_row,
        legend_label_qt,
    )
    return widget, {}
