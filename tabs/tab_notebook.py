"""Tab: Notebook — Jupyter-style code cells with inline output."""
from __future__ import annotations
from typing import TYPE_CHECKING

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPlainTextEdit, QLabel,
    QPushButton, QFrame, QScrollArea,
)
from qtpy.QtGui import QFont, QSyntaxHighlighter, QTextCharFormat, QColor
from qtpy.QtCore import Qt

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


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
    """A single notebook cell with code editor and output area."""

    def __init__(self, cell_number, code="", on_run=None, on_delete=None, parent=None):
        super().__init__(parent)
        self.cell_number = cell_number
        self.from_journal = False
        self.setFrameShape(QFrame.StyledPanel)
        self.setLineWidth(1)

        layout = QVBoxLayout()
        layout.setContentsMargins(6, 4, 6, 4)

        # Header row
        header = QHBoxLayout()
        self.label = QLabel(f"In [{cell_number}]:")
        font = QFont("monospace")
        font.setBold(True)
        self.label.setFont(font)
        header.addWidget(self.label)
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
        self.code_edit.textChanged.connect(self._auto_resize)
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

    def _auto_resize(self):
        doc_height = self.code_edit.document().size().toSize().height() + 10
        self.code_edit.setFixedHeight(max(60, min(300, doc_height)))

    def get_code(self):
        return self.code_edit.toPlainText()

    def set_code(self, code):
        self.code_edit.setPlainText(code)

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
            lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.output_layout.addWidget(lbl)
            has_output = True

        # Stdout
        if result.stdout.strip():
            lbl = QLabel(result.stdout.rstrip())
            lbl.setFont(mono)
            lbl.setWordWrap(True)
            lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.output_layout.addWidget(lbl)
            has_output = True

        # Figures
        for pixmap in result.figures:
            fig_lbl = QLabel()
            fig_lbl.setPixmap(pixmap)
            fig_lbl.setAlignment(Qt.AlignLeft)
            self.output_layout.addWidget(fig_lbl)
            has_output = True

        # Error
        if result.error:
            err_lbl = QLabel(result.error.rstrip())
            err_lbl.setFont(mono)
            err_lbl.setStyleSheet("color: red;")
            err_lbl.setWordWrap(True)
            err_lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.output_layout.addWidget(err_lbl)
            has_output = True
        elif result.stderr.strip():
            err_lbl = QLabel(result.stderr.rstrip())
            err_lbl.setFont(mono)
            err_lbl.setStyleSheet("color: #CC6600;")
            err_lbl.setWordWrap(True)
            err_lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
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
    journal_cell_count = [0]

    # ── Layout ────────────────────────────────────────────────────────────
    outer = QWidget()
    outer_layout = QVBoxLayout()
    outer_layout.setContentsMargins(4, 4, 4, 4)

    # Toolbar
    toolbar = QHBoxLayout()
    sync_btn = QPushButton("Sync Journal")
    add_btn = QPushButton("+ Cell")
    run_all_btn = QPushButton("Run All")
    clear_btn = QPushButton("Clear Outputs")
    for btn in [sync_btn, add_btn, run_all_btn, clear_btn]:
        toolbar.addWidget(btn)
    outer_layout.addLayout(toolbar)

    # Scrollable cell area
    cell_container = QWidget()
    cell_layout = QVBoxLayout()
    cell_layout.setContentsMargins(0, 0, 0, 0)
    cell_layout.addStretch()
    cell_container.setLayout(cell_layout)

    scroll = QScrollArea()
    scroll.setWidget(cell_container)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.NoFrame)
    outer_layout.addWidget(scroll)
    outer.setLayout(outer_layout)

    # ── Helpers ───────────────────────────────────────────────────────────
    def _get_engine():
        if _engine[0] is None:
            from utils.notebook_engine import NotebookEngine
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

    def _add_cell(code="", from_journal=False):
        cell_counter[0] += 1
        cell = NotebookCell(
            cell_counter[0], code,
            on_run=_on_run_cell,
            on_delete=_on_delete_cell,
        )
        cell.from_journal = from_journal
        cells.append(cell)
        # Insert before the stretch
        cell_layout.insertWidget(cell_layout.count() - 1, cell)
        return cell

    def _sync_journal():
        journal = state.get("code_journal", [])
        new_count = len(journal)
        old_count = journal_cell_count[0]
        for i in range(old_count, new_count):
            _add_cell(code=journal[i], from_journal=True)
        journal_cell_count[0] = new_count

    def _run_all():
        engine = _get_engine()
        for cell in cells:
            result = engine.run_cell(cell.get_code())
            cell.display_result(result)

    def _clear_all_outputs():
        for cell in cells:
            cell.clear_output()

    # ── Connect buttons ──────────────────────────────────────────────────
    sync_btn.clicked.connect(_sync_journal)
    add_btn.clicked.connect(lambda: _add_cell())
    run_all_btn.clicked.connect(_run_all)
    clear_btn.clicked.connect(_clear_all_outputs)

    # ── Register auto-sync callback ──────────────────────────────────────
    state["_notebook_sync_fn"] = _sync_journal

    # ── Welcome cell ─────────────────────────────────────────────────────
    _add_cell(
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

    # ── Initial sync ─────────────────────────────────────────────────────
    _sync_journal()

    return outer, {}
