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

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPlainTextEdit, QLabel,
    QPushButton, QFrame, QScrollArea, QFileDialog, QMessageBox,
)
from qtpy.QtGui import QFont, QSyntaxHighlighter, QTextCharFormat, QColor
from qtpy.QtCore import Qt

from palms.tabs._helpers import PROV_GRAPH_SIDECAR, toolbar_row
from palms.utils import stale_results
from palms.utils.prov_graph import NOTE, TEMPLATE_HAND_EDITED
from palms.utils.reporting import get_logger

log = get_logger(__name__)

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


# ── Pruning stale nodes ──────────────────────────────────────────────────────

@dataclass(frozen=True)
class PrunePlan:
    """Which stale nodes can be dropped, and which are held by a fresh one."""

    #: Ids to remove, in reverse topological order — a parent is only removed
    #: after everything depending on it, which is the order ``ProvGraph.remove``
    #: requires.
    remove: tuple[str, ...] = ()
    #: (stale id, the non-stale nodes that still depend on it). Kept, because
    #: removing one would break a step that is current.
    blocked: tuple[tuple[str, tuple[str, ...]], ...] = ()

    @property
    def is_empty(self) -> bool:
        return not self.remove


def plan_prune(graph) -> PrunePlan:
    """Every stale node that can be removed without orphaning a fresh one.

    The stale set is **not** dependency-closed. ``upsert`` clears ``stale`` on
    the node it re-records while flagging that node's descendants, so a step run
    *after* a stale one is fresh and depends on it — and ``ProvGraph.remove``
    refuses any node another still names. Both halves are handled here: the
    removable set excludes anything reachable from a fresh node, and what
    survives is emitted in reverse topological order so each removal is legal
    when it happens.
    """
    if graph is None or not len(graph):
        return PrunePlan()

    stale = {n.id for n in graph.nodes() if n.stale}
    if not stale:
        return PrunePlan()

    # Walk up from every fresh node: anything it can reach must stay.
    held: dict[str, set[str]] = {}
    for node in graph.nodes():
        if node.id in stale:
            continue
        seen: set[str] = set()
        stack = list(node.deps)
        while stack:
            dep = stack.pop()
            if dep in seen:
                continue
            seen.add(dep)
            if dep in stale:
                held.setdefault(dep, set()).add(node.id)
            parent = graph.get(dep)
            if parent is not None:
                stack.extend(parent.deps)

    try:
        order = graph.topo_sort()
    except Exception:
        # A cycle means the graph is already broken; pruning is not the repair.
        return PrunePlan()

    remove = tuple(reversed([nid for nid in order
                             if nid in stale and nid not in held]))
    blocked = tuple(sorted(
        (nid, tuple(sorted(holders))) for nid, holders in held.items()))
    return PrunePlan(remove=remove, blocked=blocked)


def describe_prune(plan: PrunePlan, orphans: dict[str, tuple[str, ...]],
                   inventory_known: bool = True) -> str:
    """The confirm-dialog body.

    *orphans* maps an id to the results it left in the store. ``inventory_known``
    is False when no scan has happened this session, in which case the dialog has
    to say that it cannot tell rather than implying there is nothing there — an
    empty ``orphans`` means two very different things.
    """
    lines = [f"Remove {len(plan.remove)} stale step(s) from the provenance graph:"]
    lines += [f"  • {nid}" for nid in plan.remove]
    if not inventory_known:
        lines += [
            "",
            "Whether these steps left results in the dataset was not checked —",
            "run Tools → Dataset → Scan Dataset first if you want to know.",
        ]
    if orphans:
        lines += [
            "",
            "⚠ These steps still have results stored in the dataset. Dropping the",
            "  step leaves the result with nothing to explain it — clear them first",
            "  with Tools → Dataset → Select Stale Results if that is what you want:",
        ]
        for nid in plan.remove:
            for key in orphans.get(nid, ()):
                lines.append(f"      {nid} → {key}")
    if plan.blocked:
        lines += ["", "Kept — a step that is still current depends on them:"]
        for nid, holders in plan.blocked:
            lines.append(f"  • {nid} (required by {', '.join(holders)})")
    lines += [
        "",
        "The notebook loses these steps permanently. A copy of the graph is saved",
        f"to {PROV_GRAPH_SIDECAR.replace('.json', '.backup_<time>.json')} first.",
    ]
    return "\n".join(lines)


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
    prune_btn = QPushButton("Drop Stale Nodes...")
    export_btn = QPushButton("Export .ipynb")
    outer_layout.addWidget(toolbar_row(
        sync_btn, add_btn, run_all_btn, clear_btn, dag_btn, prune_btn, export_btn
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

    def _persist_pruned_graph(graph) -> str:
        """Write the shrunken graph to *both* records that store it.

        The sidecar alone is not enough, and the failure is silent. Two guards —
        ``app._load_prov_graph_items`` and ``session._build_session_attrs`` —
        refuse a graph smaller than the one already stored, on the stated
        assumption that nothing in the GUI removes nodes. Writing only the
        sidecar therefore has the pruned steps restored from the session attr at
        the next exit *and* the next launch, with a printed line about a "partial
        graph" as the only clue. Writing both leaves the two the same size, so
        neither guard fires and neither needs weakening.

        Returns "" on success, or a message naming what could not be written —
        in which case the guards will restore the pre-prune graph, which is a
        safe way to fail but has to be said out loud.
        """
        from palms.utils import zarr_safe

        ctx.save_prov_graph()
        cache = None
        if ctx.sdata is not None and getattr(ctx.sdata, "path", None):
            cache = zarr_safe.cache_path_of(ctx.sdata)
        elif ctx.data_path is not None:
            candidate = Path(ctx.data_path) / "sdata_cached.zarr"
            cache = candidate if candidate.is_dir() else None
        if cache is None:
            return ("no zarr cache for this dataset, so only the sidecar was "
                    "updated")
        try:
            with zarr_safe.safe_group_update(cache, "viewer_session") as (group, _):
                group.attrs["prov_graph"] = graph.to_list()
        except Exception as e:
            log.warning("could not write the pruned provenance graph to the "
                        "session attrs: %s", e)
            return (f"the session copy could not be updated ({e}); the pruned "
                    "steps will come back on the next launch")
        return ""

    def _backup_graph(graph) -> str:
        """Snapshot the graph before pruning. Returns the path, or "" if it failed.

        ``prov_graph.backup_*.json`` is already protected from deletion by
        ``store_inventory``'s prov-graph prefix rule, so the copy cannot be
        swept away by the Dataset tab afterwards.
        """
        try:
            from palms.utils.adata_persistence import sidecar_write_path
            from palms.utils.zarr_safe import atomic_json
            stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            name = PROV_GRAPH_SIDECAR.replace(".json", f".backup_{stamp}.json")
            path = sidecar_write_path(ctx, name)
            atomic_json(path, graph.to_list())
            return str(path)
        except Exception as e:
            log.warning("could not back up the provenance graph: %s", e)
            return ""

    def _stale_artifacts(graph) -> dict:
        """Which pruned steps still have results in the store, id → keys.

        Reuses the inventory the Dataset tab already scanned rather than walking
        the dataset itself: ``build_inventory`` sizes every file under a store
        that is routinely tens of gigabytes, and doing that from a button
        callback would freeze the GUI for as long as it took — to decorate a
        dialog. If no scan has happened this session the dialog says so instead,
        which is a smaller cost than the freeze and is honest about what it does
        not know.
        """
        sections = state.get("_dataset_sections")
        if not sections:
            return {}
        try:
            return {nid: keys
                    for nid, keys in stale_results.select_stale(graph, sections).matched}
        except Exception as e:
            log.warning("could not resolve stored results for the prune dialog: %s", e)
            return {}

    def _on_prune_stale():
        graph = state.get("prov_graph")
        plan = plan_prune(graph)
        if plan.is_empty:
            if plan.blocked:
                ctx.set_status(
                    f"{len(plan.blocked)} stale step(s), all still required by "
                    "steps that are current — nothing to drop.")
            else:
                ctx.set_status("Nothing in the provenance graph is stale.")
            return

        body = describe_prune(plan, _stale_artifacts(graph),
                              inventory_known=bool(state.get("_dataset_sections")))
        if QMessageBox.question(
            None, "Drop stale nodes from the notebook", body,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        ) != QMessageBox.StandardButton.Yes:
            ctx.set_status("Prune cancelled.")
            return

        backup = _backup_graph(graph)
        dropped, failed = [], []
        for node_id in plan.remove:
            try:
                graph.remove(node_id)
                dropped.append(node_id)
            except ValueError as e:
                # Only reachable if the plan and the graph disagree, which is a
                # bug — name it rather than pruning around it.
                failed.append(f"{node_id}: {e}")
        problem = _persist_pruned_graph(graph) if dropped else ""
        _sync_from_graph()

        msg = f"Dropped {len(dropped)} stale step(s)."
        if backup:
            msg += f" Backup: {backup}."
        if failed:
            msg += f" Could not drop: {'; '.join(failed)}."
        if problem:
            msg += f" ⚠ {problem}."
        ctx.set_status(msg)
        log.info("prune: %s", msg)

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
    prune_btn.clicked.connect(_on_prune_stale)
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
