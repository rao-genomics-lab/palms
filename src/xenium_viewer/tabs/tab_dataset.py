"""Tools → Dataset: what this dataset holds on disk, and what may be removed.

The Cache tab answers "is the store healthy". This one answers the two questions
it cannot: *where did the space go*, and *how do I get rid of something I no
longer want*. Both need a per-item view, so the model lives in
:mod:`utils.store_inventory` — pure, filesystem-only and tested on its own — and
this file is the tree plus the executor.

The original 10x output is listed here, greyed out. Showing it is the point: the
guarantee that the viewer never deletes it is visible rather than asserted, and
every path this tab hands to the filesystem goes through
``store_inventory.assert_node_deletable`` first.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from magicgui.widgets import PushButton
from napari.qt.threading import thread_worker
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QLabel, QTextEdit, QTreeWidget, QTreeWidgetItem

from xenium_viewer.tabs._helpers import (
    StatusProxy,
    make_progress_bar,
    make_tab,
)
from xenium_viewer.utils import store_inventory
from xenium_viewer.utils.cache_repair import human_bytes

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def _cache_path(ctx) -> "Path | None":
    """Same resolution as the Cache tab: the live store, else the usual name."""
    if ctx.sdata is not None and getattr(ctx.sdata, "path", None):
        return Path(ctx.sdata.path)
    if ctx.data_path is not None:
        candidate = Path(ctx.data_path) / "sdata_cached.zarr"
        if candidate.exists():
            return candidate
    return None


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    status = StatusProxy(ctx.viewer)
    progress = make_progress_bar()

    header = QLabel()
    header.setWordWrap(True)

    tree = QTreeWidget()
    tree.setColumnCount(3)
    tree.setHeaderLabels(["Item", "Size", "Detail"])
    # Check state is the selection mechanism; a row highlight on top of it only
    # invites the reading that the highlighted row is the one about to go.
    tree.setSelectionMode(QTreeWidget.NoSelection)
    tree.setColumnWidth(0, 240)
    tree.setColumnWidth(1, 80)
    tree.setMinimumHeight(320)
    tree.setUniformRowHeights(True)

    report_text = QTextEdit()
    report_text.setReadOnly(True)
    report_text.setFontFamily("monospace")
    report_text.setMinimumHeight(120)
    report_text.setVisible(False)

    scan_btn = PushButton(label="Scan Dataset")
    expand_btn = PushButton(label="Expand All")
    collapse_btn = PushButton(label="Collapse All")

    def _set_busy(busy: bool):
        progress.setVisible(busy)
        scan_btn.enabled = not busy

    def _no_cache_message() -> str:
        if ctx.no_cache:
            return ("Running with --no-cache: nothing is persisted, so there is "
                    "nothing here to inspect or remove.")
        return "No zarr cache found for this dataset."

    # ── Building the tree ────────────────────────────────────────────────
    def _checked_keys() -> set[str]:
        """Keys currently ticked, so a Refresh does not lose a selection."""
        found: set[str] = set()
        stack = [tree.topLevelItem(i) for i in range(tree.topLevelItemCount())]
        while stack:
            item = stack.pop()
            if item is None:
                continue
            key = item.data(0, Qt.UserRole)
            if key and item.checkState(0) == Qt.Checked:
                found.add(key)
            stack += [item.child(i) for i in range(item.childCount())]
        return found

    def _add_node(node, parent_item, restore: set[str]) -> QTreeWidgetItem:
        item = QTreeWidgetItem(parent_item)
        item.setText(0, _display_name(node))
        item.setText(1, human_bytes(node.size_bytes) if node.size_bytes else "")
        item.setText(2, node.detail or node.blocked_reason)
        item.setData(0, Qt.UserRole, node.key)
        if node.deletable:
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(
                0, Qt.Checked if node.key in restore else Qt.Unchecked)
            if node.recoverable == store_inventory.RECOVER_NONE:
                item.setToolTip(0, "Not recoverable — no copy is kept.")
        else:
            # Disabled rather than merely uncheckable: the reason is in the
            # Detail column, and a greyed row reads as "not yours to remove".
            item.setFlags(item.flags() & ~Qt.ItemIsUserCheckable)
            item.setDisabled(True)
            if node.blocked_reason:
                item.setToolTip(0, node.blocked_reason)
        return item

    def _populate(sections, restore: set[str]):
        tree.clear()
        for section in sections:
            section_item = QTreeWidgetItem(tree)
            section_item.setText(0, section.title)
            section_item.setText(1, human_bytes(section.total_bytes)
                                 if section.total_bytes else "")
            section_item.setText(2, section.note)
            section_item.setFirstColumnSpanned(False)
            font = section_item.font(0)
            font.setBold(True)
            section_item.setFont(0, font)
            items: dict[str, QTreeWidgetItem] = {}
            for node in section.nodes:
                parent_item = items.get(node.parent, section_item)
                items[node.key] = _add_node(node, parent_item, restore)
            section_item.setExpanded(True)

    def _refresh_header(sections=None):
        cache_path = _cache_path(ctx)
        lines = [f"<b>{ctx.data_path}</b>"]
        if cache_path is None:
            lines.append(_no_cache_message())
        elif sections is not None:
            total = sum(s.total_bytes for s in sections)
            viewer_total = sum(s.total_bytes for s in sections
                               if s.title != "Original Xenium output")
            lines.append(
                f"{human_bytes(total)} on disk · "
                f"{human_bytes(viewer_total)} of it written by the viewer"
            )
        else:
            lines.append("Press <b>Scan Dataset</b> to measure what is on disk.")
        lines.append("Tick what you want removed. The original Xenium output is "
                     "listed read-only and can never be selected.")
        header.setText("<br>".join(lines))

    def _on_scan():
        if ctx.data_path is None:
            header.setText("No dataset directory for this session.")
            return
        data_path = Path(ctx.data_path)
        cache_path = _cache_path(ctx)
        restore = _checked_keys()
        _set_busy(True)
        status.value = "Scanning dataset..."

        @thread_worker
        def _run():
            return store_inventory.build_inventory(data_path, cache_path)

        def _done(sections):
            state["_dataset_sections"] = sections
            _populate(sections, restore)
            _refresh_header(sections)
            _after_refresh(sections)
            status.value = "Dataset scanned."
            _set_busy(False)

        _start(_run, _done, "Scan failed")

    def _after_refresh(sections):
        """Hook the deletion half fills in; a no-op for the read-only tree."""
        return None

    # Nothing scans at build time, deliberately: the walk covers the whole
    # dataset directory, and charging every launch for a tab most sessions never
    # open is the same mistake as loading a 30 GB store to show a size.

    def _start(make_worker, on_done, error_prefix):
        worker = make_worker()

        def _failed(exc):
            # napari's `errored` emits the exception itself, not an exc_info
            # triple, so the traceback has to be rebuilt from it.
            import traceback
            detail = "".join(traceback.format_exception(
                type(exc), exc, exc.__traceback__))
            report_text.setVisible(True)
            report_text.setPlainText(f"{error_prefix}: {exc}\n\n{detail}")
            status.value = f"{error_prefix}: {exc}"
            _set_busy(False)

        worker.returned.connect(on_done)
        worker.errored.connect(_failed)
        # Keep a reference so the worker is not garbage-collected mid-run.
        state["_dataset_worker"] = worker
        worker.start()

    scan_btn.clicked.connect(_on_scan)
    expand_btn.clicked.connect(tree.expandAll)
    collapse_btn.clicked.connect(tree.collapseAll)

    # Not disabled under --no-cache: this session persists nothing, but earlier
    # ones left sidecars, a transcript cache and backups on disk, and seeing —
    # and reclaiming — those is exactly what this tab is for.
    _refresh_header()

    widget = make_tab(
        header, scan_btn, progress, tree,
        expand_btn, collapse_btn, report_text,
    )

    def _restore_session(session):
        return None

    return widget, {"restore_session": _restore_session}


def _display_name(node) -> str:
    """Row label. Table contents are tagged, since 'obs' vs 'uns' matters."""
    if node.kind in (store_inventory.OBS, store_inventory.UNS,
                     store_inventory.OBSM):
        return f"{node.kind.lower()}  {node.name}"
    return node.name
