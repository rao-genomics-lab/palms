"""Tab: Notebook — Jupyter-style code cells derived from the provenance graph.

Cells are *derived* from ``state["prov_graph"]``: one cell per node, in
topological order, with a staleness badge when an upstream input changed. The
graph is the source of truth, so re-running a step revises its cell in place and
the whole notebook stays dependency-ordered — including across sessions. Users
can still add/edit/run free-form cells; "Export .ipynb" writes a standalone,
replayable notebook (graph cells + any user cells).
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPlainTextEdit, QLabel,
    QPushButton, QFrame, QScrollArea, QFileDialog,
)
from qtpy.QtGui import QFont, QSyntaxHighlighter, QTextCharFormat, QColor
from qtpy.QtCore import Qt

from palms.tabs._helpers import toolbar_row
from palms.utils.prov_graph import NOTE, TEMPLATE_HAND_EDITED

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


# ── Folding user edits back into the graph ───────────────────────────────────

def reconcile_edits(graph, cells) -> list[str]:
    """Fold user edits of graph cells back into *graph*. Returns the ids touched.

    The edited text is written straight onto the node **without being
    executed**, so from that moment the node's ``code`` is no longer the source
    that produced the artifact the viewer is showing. Every other node in the
    graph carries exactly that guarantee, and nothing in the exported notebook
    distinguishes this one — so it is marked ``hand-edited``. ``template_hash``
    is cleared because the code no longer derives from any template.

    Module-level rather than a closure so it can be tested against a real graph
    without building the tab (which needs a napari viewer).
    """
    touched: list[str] = []
    if graph is None:
        return touched
    for cell in cells:
        if cell.node_id is None or not cell.edited_by_user:
            continue
        node = graph.get(cell.node_id)
        if node is None:
            continue
        node.code = cell.get_code()
        node.template_origin = TEMPLATE_HAND_EDITED
        node.template_hash = None
        cell.edited_by_user = False
        touched.append(node.id)
    return touched


# ── Syntax highlighting ──────────────────────────────────────────────────────

class PythonSyntaxHighlighter(QSyntaxHighlighter):
    """Lightweight Python syntax highlighter using pygments."""

    def __init__(self, parent=None):
        super().__init__(parent)
        from pygments.lexers import PythonLexer
        self._lexer = PythonLexer()
        self._formats = {}
        self._build_formats()

    def _build_formats(self):
        from pygments.token import Token
        color_map = {
            Token.Keyword: ("#0000FF", False),
            Token.Name.Builtin: ("#008080", False),
            Token.Literal.String: ("#008000", False),
            Token.Comment: ("#808080", True),
            Token.Literal.Number: ("#FF8000", False),
            Token.Operator: ("#AA22FF", False),
            Token.Name.Decorator: ("#AA22FF", True),
        }
        for tok, (color, italic) in color_map.items():
            fmt = QTextCharFormat()
            fmt.setForeground(QColor(color))
            if italic:
                fmt.setFontItalic(True)
            self._formats[tok] = fmt

    def highlightBlock(self, text):
        from pygments import lex
        index = 0
        for token_type, value in lex(text, self._lexer):
            length = len(value)
            tt = token_type
            while tt not in self._formats and tt.parent:
                tt = tt.parent
            if tt in self._formats:
                self.setFormat(index, length, self._formats[tt])
            index += length


# ── Cell widget ──────────────────────────────────────────────────────────────

class NotebookCell(QFrame):
    """A single notebook cell with code editor and output area.

    ``node_id`` links the cell to a provenance-graph node (None for free-form
    user cells). Graph cells are refreshed from the graph unless the user has
    edited them.
    """

    def __init__(self, cell_number, code="", on_run=None, on_delete=None,
                 node_id=None, node_label=None, stale=False, parent=None):
        super().__init__(parent)
        self.cell_number = cell_number
        self.node_id = node_id
        self.edited_by_user = False
        self._loading = True
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setLineWidth(1)

        layout = QVBoxLayout()
        layout.setContentsMargins(6, 4, 6, 4)

        # Header row
        header = QHBoxLayout()
        self.label = QLabel()
        font = QFont("monospace")
        font.setBold(True)
        self.label.setFont(font)
        header.addWidget(self.label)
        self.stale_label = QLabel("⚠ stale — input changed; re-run in the viewer")
        self.stale_label.setStyleSheet("color: #CC6600;")
        self.stale_label.setVisible(False)
        header.addWidget(self.stale_label)
        header.addStretch()

        self.run_btn = QPushButton("Run")
        self.run_btn.setFixedWidth(50)
        self.delete_btn = QPushButton("Delete")
        self.delete_btn.setFixedWidth(60)
        header.addWidget(self.run_btn)
        header.addWidget(self.delete_btn)
        layout.addLayout(header)

        # Code editor
        self.code_edit = QPlainTextEdit()
        self.code_edit.setFont(QFont("monospace", 10))
        self.code_edit.setTabStopDistance(28)
        self.code_edit.setPlainText(code)
        self.code_edit.setMinimumHeight(60)
        self.code_edit.setMaximumHeight(300)
        self._auto_resize()
        self.code_edit.textChanged.connect(self._on_text_changed)
        layout.addWidget(self.code_edit)

        # Syntax highlighter
        self._highlighter = PythonSyntaxHighlighter(self.code_edit.document())

        # Output area (hidden initially)
        self.output_area = QWidget()
        self.output_layout = QVBoxLayout()
        self.output_layout.setContentsMargins(8, 4, 8, 4)
        self.output_area.setLayout(self.output_layout)
        self.output_area.hide()
        layout.addWidget(self.output_area)

        self.setLayout(layout)

        # Connect signals
        if on_run:
            self.run_btn.clicked.connect(lambda: on_run(self))
        if on_delete:
            self.delete_btn.clicked.connect(lambda: on_delete(self))

        self._loading = False
        self.set_meta(node_label, stale)

    def _on_text_changed(self):
        self._auto_resize()
        if not self._loading:
            self.edited_by_user = True

    def _auto_resize(self):
        doc_height = self.code_edit.document().size().toSize().height() + 10
        self.code_edit.setFixedHeight(max(60, min(300, doc_height)))

    def get_code(self):
        return self.code_edit.toPlainText()

    def set_code(self, code):
        self._loading = True
        self.code_edit.setPlainText(code)
        self._loading = False

    def set_meta(self, node_label, stale):
        """Update the header label + staleness badge (graph cells)."""
        self.stale = bool(stale)
        if self.node_id is not None:
            self.label.setText(node_label or self.node_id)
        else:
            self.label.setText(f"In [{self.cell_number}]:")
        self.stale_label.setVisible(self.stale)

    def clear_output(self):
        while self.output_layout.count():
            item = self.output_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self.output_area.hide()

    def display_result(self, result):
        self.clear_output()
        has_output = False
        mono = QFont("monospace", 9)

        # Return value
        if result.result_repr and result.result_repr != "None":
            lbl = QLabel(f"Out: {result.result_repr}")
            lbl.setFont(mono)
            lbl.setWordWrap(True)
            lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.output_layout.addWidget(lbl)
            has_output = True

        # Stdout
        if result.stdout.strip():
            lbl = QLabel(result.stdout.rstrip())
            lbl.setFont(mono)
            lbl.setWordWrap(True)
            lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.output_layout.addWidget(lbl)
            has_output = True

        # Figures
        for pixmap in result.figures:
            fig_lbl = QLabel()
            fig_lbl.setPixmap(pixmap)
            fig_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft)
            self.output_layout.addWidget(fig_lbl)
            has_output = True

        # Error
        if result.error:
            err_lbl = QLabel(result.error.rstrip())
            err_lbl.setFont(mono)
            err_lbl.setStyleSheet("color: red;")
            err_lbl.setWordWrap(True)
            err_lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.output_layout.addWidget(err_lbl)
            has_output = True
        elif result.stderr.strip():
            err_lbl = QLabel(result.stderr.rstrip())
            err_lbl.setFont(mono)
            err_lbl.setStyleSheet("color: #CC6600;")
            err_lbl.setWordWrap(True)
            err_lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            self.output_layout.addWidget(err_lbl)
            has_output = True

        if has_output:
            self.output_area.show()


# ── Tab builder ──────────────────────────────────────────────────────────────

def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    _engine = [None]  # mutable container for lazy init
    cells = []
    cell_counter = [0]

    # ── Layout ────────────────────────────────────────────────────────────
    outer = QWidget()
    outer_layout = QVBoxLayout()
    outer_layout.setContentsMargins(4, 4, 4, 4)

    # Toolbar. Pinned above the cells but horizontally scrollable — six buttons
    # in a plain row cannot shrink below their labels, which made this tab the
    # floor for the whole control dock's width.
    sync_btn = QPushButton("Sync Graph")
    add_btn = QPushButton("+ Cell")
    run_all_btn = QPushButton("Run All")
    clear_btn = QPushButton("Clear Outputs")
    dag_btn = QPushButton("Show DAG")
    export_btn = QPushButton("Export .ipynb")
    outer_layout.addWidget(toolbar_row(
        sync_btn, add_btn, run_all_btn, clear_btn, dag_btn, export_btn
    ))

    # Scrollable cell area
    cell_container = QWidget()
    cell_layout = QVBoxLayout()
    cell_layout.setContentsMargins(0, 0, 0, 0)
    cell_layout.addStretch()
    cell_container.setLayout(cell_layout)

    scroll = QScrollArea()
    scroll.setWidget(cell_container)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.Shape.NoFrame)
    outer_layout.addWidget(scroll)
    outer.setLayout(outer_layout)

    # ── Helpers ───────────────────────────────────────────────────────────
    def _get_engine():
        if _engine[0] is None:
            from palms.utils.notebook_engine import NotebookEngine
            _engine[0] = NotebookEngine(ctx.viewer)
        return _engine[0]

    def _on_run_cell(cell):
        engine = _get_engine()
        result = engine.run_cell(cell.get_code())
        cell.display_result(result)

    def _on_delete_cell(cell):
        if cell in cells:
            cells.remove(cell)
        cell_layout.removeWidget(cell)
        cell.deleteLater()

    def _add_cell(code="", node_id=None, node_label=None, stale=False):
        cell_counter[0] += 1
        cell = NotebookCell(
            cell_counter[0], code,
            on_run=_on_run_cell,
            on_delete=_on_delete_cell,
            node_id=node_id, node_label=node_label, stale=stale,
        )
        cells.append(cell)
        # Insert before the trailing stretch
        cell_layout.insertWidget(cell_layout.count() - 1, cell)
        return cell

    def _reorder(order):
        """Graph cells in topological order first, then free-form user cells."""
        pos = {nid: i for i, nid in enumerate(order)}
        graph_cells = sorted(
            (c for c in cells if c.node_id is not None),
            key=lambda c: pos.get(c.node_id, 1 << 30),
        )
        user_cells = [c for c in cells if c.node_id is None]
        ordered = graph_cells + user_cells
        for c in ordered:
            cell_layout.removeWidget(c)
        for c in ordered:
            cell_layout.insertWidget(cell_layout.count() - 1, c)
        cells[:] = ordered

    def _sync_from_graph():
        """Render/refresh cells from the provenance graph (source of truth)."""
        graph = state.get("prov_graph")
        if graph is None:
            return
        try:
            order = graph.topo_sort()
        except Exception:
            return
        existing = {c.node_id: c for c in cells if c.node_id is not None}
        seen = set()
        for nid in order:
            node = graph.get(nid)
            if node is None:
                continue
            # A note is viewer state, not code — say so in its header, since
            # here (unlike the exported notebook) it still shows as a cell.
            label = node.label
            if node.kind == NOTE:
                label = f"{label or nid} — viewer state, not code"
            cell = existing.get(nid)
            if cell is not None:
                if not cell.edited_by_user and cell.get_code() != node.code:
                    cell.set_code(node.code)
                cell.set_meta(label, node.stale)
            else:
                _add_cell(code=node.code, node_id=nid,
                          node_label=label, stale=node.stale)
            seen.add(nid)
        # Drop cells whose node was removed from the graph
        for c in list(cells):
            if c.node_id is not None and c.node_id not in seen:
                _on_delete_cell(c)
        _reorder(order)

    def _run_all():
        engine = _get_engine()
        for cell in cells:
            result = engine.run_cell(cell.get_code())
            cell.display_result(result)

    def _clear_all_outputs():
        for cell in cells:
            cell.clear_output()

    def _reconcile_edits():
        reconcile_edits(state.get("prov_graph"), cells)

    def _export_cells():
        """(cell_type, source) list for export: graph cells + user cells."""
        from palms.utils import notebook_export
        _reconcile_edits()
        graph = state.get("prov_graph")
        out = notebook_export.graph_to_cells(graph) if graph is not None else []
        for c in cells:
            if c.node_id is None and c is not welcome_cell:
                src = c.get_code().strip()
                if src:
                    out.append(("code", src))
        return out

    def _on_show_dag():
        from palms.utils.dag_view import render_dag
        graph = state.get("prov_graph")
        if graph is None or len(graph) == 0:
            ctx.set_status("Provenance graph is empty — record some steps first")
            return
        try:
            fig = render_dag(graph)
            paths = ctx.show_plot(fig, "provenance_dag",
                                  title=f"Provenance DAG ({len(graph)} nodes)")
            ctx.set_status(f"Provenance DAG ({len(graph)} nodes) — "
                           f"saved to {', '.join(paths)}")
        except Exception as e:
            ctx.set_status(f"DAG render failed: {e}")

    def _on_export_ipynb():
        from palms.utils import notebook_export
        default = str(ctx.data_path / "analysis_notebook.ipynb")
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Notebook", default, "Jupyter Notebook (*.ipynb)",
        )
        if not path:
            return
        try:
            notebook_export.write_notebook(_export_cells(), path)
            ctx.set_status(f"Notebook exported to {path}")
        except Exception as e:
            ctx.set_status(f"Notebook export failed: {e}")

    # ── Connect buttons ──────────────────────────────────────────────────
    sync_btn.clicked.connect(_sync_from_graph)
    add_btn.clicked.connect(lambda: _add_cell())
    run_all_btn.clicked.connect(_run_all)
    clear_btn.clicked.connect(_clear_all_outputs)
    dag_btn.clicked.connect(_on_show_dag)
    export_btn.clicked.connect(_on_export_ipynb)

    # ── Register auto-sync callback (called by record_node) ──────────────
    state["_notebook_sync_fn"] = _sync_from_graph

    # ── Welcome cell ─────────────────────────────────────────────────────
    welcome_cell = _add_cell(
        code=(
            "# Available objects in this session\n"
            "print('Available variables:')\n"
            "print('  adata        - AnnData object')\n"
            "print('  sdata        - SpatialData container')\n"
            "print('  viewer       - napari Viewer')\n"
            "print('  ctx          - ViewerContext (all state)')\n"
            "print('  clusterings  - dict of clustering assignments')\n"
            "print('  color_manager - CellColorManager')\n"
            "print('  gene_names   - list of gene names')\n"
            "print('  data_path    - Path to dataset')"
        ),
    )

    # ── Initial render from the graph ────────────────────────────────────
    _sync_from_graph()

    return outer, {
        "restore_session": lambda session: _sync_from_graph(),
        "get_export_cells": _export_cells,
    }
