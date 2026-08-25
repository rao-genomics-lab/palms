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

import gc
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from magicgui.widgets import PushButton
from napari.qt.threading import thread_worker
from qtpy.QtCore import Qt
from qtpy.QtGui import QBrush, QPalette
from qtpy.QtWidgets import (
    QApplication,
    QLabel,
    QMessageBox,
    QTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
)

from xenium_viewer.tabs._helpers import (
    StatusProxy,
    make_progress_bar,
    make_tab,
)
from xenium_viewer.utils import store_inventory, zarr_safe
from xenium_viewer.utils.adata_persistence import (
    CLUSTERING_PREFIX,
    _persist_table,
)
from xenium_viewer.utils.cache_repair import human_bytes
from xenium_viewer.utils.reporting import get_logger, report_write_failure

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

log = get_logger(__name__)


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
    tree.setSelectionMode(QTreeWidget.SelectionMode.NoSelection)
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
    delete_btn = PushButton(label="Delete Selected...", enabled=False)
    trash_btn = PushButton(label="Empty Trash", enabled=False)

    def _set_busy(busy: bool):
        progress.setVisible(busy)
        scan_btn.enabled = not busy
        for button in (delete_btn, trash_btn):
            button.enabled = not busy and bool(state.get("_dataset_sections"))

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
            key = item.data(0, Qt.ItemDataRole.UserRole)
            if key and item.checkState(0) == Qt.CheckState.Checked:
                found.add(key)
            stack += [item.child(i) for i in range(item.childCount())]
        return found

    def _populate(sections, restore: set[str]):
        _populate_tree(tree, sections, restore)

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
        delete_btn.enabled = True
        trash_btn.enabled = any(
            n.key == "trash:all" for s in sections for n in s.nodes)

    # ── Deleting ─────────────────────────────────────────────────────────
    def _confirm_and_delete(keys, title):
        sections = state.get("_dataset_sections")
        if not sections:
            status.value = "Scan the dataset first."
            return
        try:
            plan = store_inventory.plan_deletion(sections, keys)
        except store_inventory.NotDeletable as exc:
            # Only reachable if a blocked row somehow became tickable, which is
            # a bug rather than a user error — say so instead of deleting.
            report_text.setVisible(True)
            report_text.setPlainText(f"Refused: {exc}")
            status.value = f"Refused: {exc}"
            return
        if plan.is_empty:
            status.value = "Nothing selected."
            return
        answer = QMessageBox.question(
            None, title, store_inventory.describe_plan(plan),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, QMessageBox.StandardButton.Cancel,
        )
        if answer != QMessageBox.StandardButton.Yes:
            status.value = "Deletion cancelled."
            return

        _set_busy(True)
        status.value = f"Deleting {len(plan.nodes)} item(s)..."

        @thread_worker
        def _run():
            return _apply_deletion(ctx, plan)

        def _done(result: "DeletionResult"):
            report_text.setVisible(True)
            report_text.setPlainText(result.summary())
            status.value = (f"Removed {len(result.removed)} item(s), "
                            f"{human_bytes(result.bytes_freed)} reclaimed."
                            if result.removed else "Nothing was removed.")
            _set_busy(False)
            _on_scan()
            # Last, because reload_dataset() tears down this very widget.
            if result.needs_reload:
                _offer_reload(result)

        _start(_run, _done, "Deletion failed")

    def _on_delete():
        _confirm_and_delete(sorted(_checked_keys()), "Delete dataset components")

    def _on_empty_trash():
        _confirm_and_delete(["trash:all"], "Empty the cache trash")

    def _offer_reload(result):
        """Element deletions need the dataset rebuilt to leave the viewer sane.

        Offered from the worker's `returned` slot, i.e. on the main thread:
        reload_dataset() destroys every tab widget, including the one whose
        callback is running.
        """
        reload_dataset = getattr(ctx, "reload_dataset", None)
        if reload_dataset is None:
            report_text.append(
                "\n\nReopen this dataset (File → Open Dataset) so the viewer "
                "stops showing what was just removed.")
            return
        answer = QMessageBox.question(
            None, "Reload dataset?",
            f"{len(result.removed)} item(s) were removed from the store.\n\n"
            "Reload the dataset now? The viewer still has the old elements "
            "loaded until you do.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.Yes,
        )
        if answer != QMessageBox.StandardButton.Yes:
            report_text.append(
                "\n\nNot reloaded. Reopen the dataset when convenient.")
            return
        status.value = "Reloading dataset..."
        reload_dataset()

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
    # Wrapped in lambdas because magicgui inspects the signature of what it
    # connects, and a Qt builtin slot has none to find.
    expand_btn.clicked.connect(lambda: tree.expandAll())
    collapse_btn.clicked.connect(lambda: tree.collapseAll())
    delete_btn.clicked.connect(_on_delete)
    trash_btn.clicked.connect(_on_empty_trash)

    # Not disabled under --no-cache: this session persists nothing, but earlier
    # ones left sidecars, a transcript cache and backups on disk, and seeing —
    # and reclaiming — those is exactly what this tab is for.
    _refresh_header()

    widget = make_tab(
        header, scan_btn, progress, tree,
        expand_btn, collapse_btn, delete_btn, trash_btn, report_text,
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


def _dim_brush():
    """The palette's disabled-text colour, so dimming follows the napari theme."""
    app = QApplication.instance()
    if app is None:
        return None
    return QBrush(app.palette().color(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text))


def _add_node(tree_parent, node, restore: set[str]) -> QTreeWidgetItem:
    item = QTreeWidgetItem(tree_parent)
    item.setText(0, _display_name(node))
    item.setText(1, human_bytes(node.size_bytes) if node.size_bytes else "")
    item.setText(2, node.detail or node.blocked_reason)
    item.setData(0, Qt.ItemDataRole.UserRole, node.key)
    if node.deletable:
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(0, Qt.CheckState.Checked if node.key in restore else Qt.CheckState.Unchecked)
        if node.recoverable == store_inventory.RECOVER_NONE:
            item.setToolTip(0, "Not recoverable — no copy is kept.")
    else:
        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsUserCheckable)
        # Dimmed, never setDisabled(True): Qt propagates a disabled item down
        # its whole subtree, so blocking `group:tables` and `tables/table` also
        # greyed out every clustering underneath them — the one thing this tab
        # exists to let you delete. Group headings keep normal text; only real
        # blocked rows are dimmed, since for them it means "not yours to remove".
        if node.kind != store_inventory.GROUP:
            brush = _dim_brush()
            if brush is not None:
                for column in range(3):
                    item.setForeground(column, brush)
        if node.blocked_reason:
            item.setToolTip(0, node.blocked_reason)
    return item


def _populate_tree(tree, sections, restore: set[str]) -> None:
    """Render *sections* as a three-level tree, restoring previous check state."""
    tree.clear()
    for section in sections:
        section_item = QTreeWidgetItem(tree)
        section_item.setText(0, section.title)
        section_item.setText(1, human_bytes(section.total_bytes)
                             if section.total_bytes else "")
        section_item.setText(2, section.note)
        font = section_item.font(0)
        font.setBold(True)
        section_item.setFont(0, font)

        items: dict[str, QTreeWidgetItem] = {}
        for node in section.nodes:
            parent_item = items.get(node.parent, section_item)
            items[node.key] = _add_node(parent_item, node, restore)

        # Open the path to anything actionable. Clusterings sit three levels
        # down, under a blocked table, and a user should not have to go hunting
        # for the rows they came here to tick.
        for node in section.nodes:
            if not node.deletable:
                continue
            item = items[node.key].parent()
            while item is not None:
                item.setExpanded(True)
                item = item.parent()
        section_item.setExpanded(True)


# ── The executor ─────────────────────────────────────────────────────────────

@dataclass
class DeletionResult:
    removed: list = field(default_factory=list)          # node names
    failed: list = field(default_factory=list)           # (name, reason)
    bytes_freed: int = 0
    needs_reload: bool = False

    def summary(self) -> str:
        lines: list[str] = []
        if self.removed:
            lines.append(f"Removed {len(self.removed)} item(s), "
                         f"{human_bytes(self.bytes_freed)} reclaimed:")
            lines += [f"    {name}" for name in self.removed]
        if self.failed:
            if lines:
                lines.append("")
            lines.append(f"✗ {len(self.failed)} item(s) could not be removed:")
            lines += [f"    {name}: {reason}" for name, reason in self.failed]
        return "\n".join(lines) or "Nothing was removed."


def _run_node(result: DeletionResult, node, action) -> bool:
    """Apply one node's deletion, keeping the batch alive if it fails.

    Per-node rather than per-batch because several of these are irreversible: a
    partially applied batch has to be describable afterwards, not rolled back.
    """
    try:
        action()
    except Exception as exc:
        report_write_failure(exc, f"delete {node.name}")
        result.failed.append((node.name, _explain(exc)))
        return False
    result.removed.append(node.name)
    result.bytes_freed += node.size_bytes or 0
    return True


def _explain(exc: Exception) -> str:
    if isinstance(exc, zarr_safe.ZarrSafeError):
        return ("in use — reload the dataset (Tools → Dataset offers this after "
                "a delete) and try again")
    return str(exc)


def _remove_tree(path: Path) -> None:
    """The one place this module removes anything from the filesystem.

    Every caller must have vetted the node with ``assert_node_deletable`` first;
    a test parses this file and fails if a removal appears anywhere else, or in a
    function that skipped the check. A symlink is unlinked rather than followed —
    deleting its target could reach outside a root even when the link itself is
    inside one.
    """
    path = Path(path)
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _remove_path(node, roots) -> None:
    """Remove the file or directory a node stands for."""
    store_inventory.assert_node_deletable(node, roots)
    _remove_tree(Path(node.path))


def _delete_table_entries(ctx, nodes, roots, result: DeletionResult) -> None:
    """Drop obs/uns/obsm keys, then persist the table exactly once.

    One write for the whole batch: ``_persist_table`` rewrites the entire table,
    so per-column writes would be N full rewrites of the hottest element in the
    store for no benefit.
    """
    adata = ctx.adata
    if adata is None:
        for node in nodes:
            result.failed.append((node.name, "no table is loaded"))
        return

    staged = []
    for node in nodes:
        try:
            store_inventory.assert_node_deletable(node, roots)
            # A column can be on disk without being in memory (a failed
            # restore, a --no-cache session, an external write). It still counts
            # as deleted: _persist_table rewrites the whole table from memory.
            if node.kind == store_inventory.OBS:
                if node.name in adata.obs.columns:
                    del adata.obs[node.name]
            elif node.kind == store_inventory.UNS:
                adata.uns.pop(node.name, None)
            elif node.kind == store_inventory.OBSM:
                if node.name in adata.obsm:
                    del adata.obsm[node.name]
        except Exception as exc:
            report_write_failure(exc, f"delete {node.name}")
            result.failed.append((node.name, _explain(exc)))
        else:
            staged.append(node)

    if not staged:
        return
    try:
        _persist_table(ctx)
    except Exception as exc:
        report_write_failure(exc, "clustering/analysis data")
        for node in staged:
            result.failed.append((node.name, _explain(exc)))
        return

    for node in staged:
        result.removed.append(node.name)
        result.bytes_freed += node.size_bytes or 0
    _forget_clusterings(ctx, staged)


def _forget_clusterings(ctx, nodes) -> None:
    """Drop deleted clusterings from the in-memory dicts the combos read.

    ``refresh_clustering_choices`` reads ``ctx.clusterings``, *not*
    ``adata.obs`` — it is a plain dict populated once at load and mutated by
    hand at each producer. Without this the deleted clustering stays in every
    combo and every ``ctx.clusterings[key]`` lookup still resolves against the
    cached Series, so the column is gone from disk and still colouring cells.
    """
    removed = [n.name[len(CLUSTERING_PREFIX):] for n in nodes
               if n.kind == store_inventory.OBS
               and n.name.startswith(CLUSTERING_PREFIX)]
    if not removed:
        return
    for name in removed:
        if isinstance(getattr(ctx, "clusterings", None), dict):
            ctx.clusterings.pop(name, None)
        custom = ctx.state.get("custom_clusterings")
        if isinstance(custom, dict):
            custom.pop(name, None)
        labels = ctx.state.get("cluster_labels")
        if isinstance(labels, dict):
            labels.pop(name, None)
    refresh = getattr(ctx, "refresh_clustering_choices", None)
    if callable(refresh):
        refresh()


def _delete_session_node(ctx, node, cache_path, roots) -> None:
    """Remove a session group/array/attr *and* clear its in-memory mirror.

    ``save_session`` rebuilds ``viewer_session`` from ctx.state / ctx.he_state /
    ctx.arms_state at exit and carries unrecognised attrs forward from what is
    already stored. Deleting only the disk copy therefore puts it straight back
    on the next clean exit — the mirror image of what
    ``tab_cache._recover_session_attrs`` has to do when recovering.
    """
    if cache_path is None:
        raise RuntimeError("no cache to delete session state from")
    name = node.key.rpartition("/")[2]
    with zarr_safe.safe_group_update(cache_path, "viewer_session") as (group, stage):
        if node.key.startswith("session:attr/"):
            _drop_attrs(group, [name])
        else:
            store_inventory.assert_node_deletable(node, roots)
            # The group yielded here is backed by staging, seeded from the live
            # copy, and swapped in on a clean exit — so removing it from staging
            # is the delete, and an exception leaves the live group untouched.
            _remove_tree(Path(stage) / name)
            # Retire the attrs the group is the storage for as well, or
            # _build_session_attrs carries them forward from prev_attrs.
            _drop_attrs(group, [k for k in dict(group.attrs)
                                if _attr_belongs_to(k, name)])

    for holder, key in store_inventory.session_memory_keys(node.key):
        target = getattr(ctx, holder, None)
        if isinstance(target, dict):
            target.pop(key, None)


def _drop_attrs(group, names) -> None:
    for name in names:
        try:
            del group.attrs[name]
        except KeyError:
            pass


_SESSION_ATTR_PREFIXES = {
    "he": ("he_", "flip_v", "flip_h"),
    "arms": ("arms_",),
}


def _attr_belongs_to(attr: str, group_name: str) -> bool:
    prefixes = _SESSION_ATTR_PREFIXES.get(group_name)
    return bool(prefixes) and attr.startswith(prefixes)


def _drop_layer(viewer, layer) -> None:
    if viewer is None or layer is None:
        return
    try:
        viewer.layers.remove(layer)
    except Exception:
        pass


def _release_layers(ctx, element_names) -> None:
    """Let go of anything holding a lazily-loaded element, before deleting it.

    ``safe_delete_element`` refuses to rename an element whose files still back
    a live dask graph (``_assert_not_dask_backed``), and a napari image layer is
    exactly such a reader. This is the teardown
    ``tab_external_images.on_remove`` already does before its own delete: pull
    the layer, close the tif, then collect.
    """
    viewer = getattr(ctx, "viewer", None)
    wanted = set(element_names)

    for store in (getattr(ctx, "external_images_state", None) or [],
                  getattr(ctx, "patch_overlays_state", None) or []):
        for entry in list(store):
            if not isinstance(entry, dict):
                continue
            if entry.get("element_name") not in wanted:
                continue
            disconnect = entry.get("affine_disconnect")
            if callable(disconnect):
                try:
                    disconnect()
                except Exception:
                    pass
            for key in ("layer_ref", "xenium_lm_layer", "image_lm_layer"):
                _drop_layer(viewer, entry.get(key))
                entry[key] = None
            tif = entry.get("tif")
            if tif is not None:
                try:
                    tif.close()
                except Exception:
                    pass
                entry["tif"] = None
            try:
                store.remove(entry)
            except ValueError:
                pass

    # H&E and ARMS keep their layer on a state dict rather than in a list.
    for holder_name, element in (("he_state", "he_image"),
                                 ("arms_state", "arms_he_image")):
        if element not in wanted:
            continue
        holder = getattr(ctx, holder_name, None)
        if not isinstance(holder, dict):
            continue
        for key in ("he_layer", "he_lm_layer", "xenium_lm_layer", "shapes_layer"):
            _drop_layer(viewer, holder.get(key))
            if key in holder:
                holder[key] = None
        tif = holder.get("he_tif")
        if tif is not None:
            try:
                tif.close()
            except Exception:
                pass
            holder["he_tif"] = None

    gc.collect()


def _apply_deletion(ctx, plan) -> DeletionResult:
    """Apply *plan*, in its own kind order, reporting per node.

    Runs on a worker: it can rmtree a multi-gigabyte backup. Nothing here
    touches Qt, and the reload offer is left to the caller's `returned` slot.
    """
    result = DeletionResult()
    cache_path = _cache_path(ctx)
    data_path = Path(ctx.data_path) if ctx.data_path is not None else None
    roots = store_inventory.deletable_roots(data_path, cache_path)

    table_nodes = plan.of_kinds(*store_inventory.TABLE_KINDS)
    if table_nodes:
        _delete_table_entries(ctx, table_nodes, roots, result)

    for node in plan.of_kinds(store_inventory.SESSION):
        _run_node(result, node,
                  lambda n=node: _delete_session_node(ctx, n, cache_path, roots))

    elements = plan.of_kinds(store_inventory.ELEMENT)
    if elements:
        _release_layers(ctx, [n.name for n in elements])
        for node in elements:
            def _delete(n=node):
                store_inventory.assert_node_deletable(n, roots)
                if ctx.sdata is None or n.name not in ctx.sdata:
                    raise RuntimeError("not in the loaded dataset any more")
                zarr_safe.safe_delete_element(ctx.sdata, n.name)
            if _run_node(result, node, _delete):
                result.needs_reload = True

    for node in plan.of_kinds(store_inventory.SIDECAR, store_inventory.DERIVED,
                              store_inventory.TRASH, store_inventory.BACKUP):
        _run_node(result, node, lambda n=node: _remove_path(n, roots))

    return result
