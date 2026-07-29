"""Tab: Cache — inspect and repair the zarr store.

The viewer's cache used to be a black box: when it broke, the loader silently
moved it aside (or deleted it) and rebuilt from the raw Xenium output, which on
a real dataset means 30 GB and a long wait. This tab exposes what the loader now
does automatically, so it can also be run deliberately:

* **Verify** — read-only health check, safe at any time.
* **Re-consolidate** — rebuilds the root metadata index. Fixes the most common
  corruption (an element present on disk but missing from the index) without
  touching data.
* **Recover from backup** — pull elements out of a ``.xv_trash`` copy or a
  previous cache the loader kept aside.
* **Rebuild** — force the from-raw rebuild, preserving user-generated data.

Every button runs in a ``thread_worker`` and every mutation goes through
``zarr_safe.store_lock``, so nothing here can race ``_persist_table``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton
from napari.qt.threading import thread_worker
from qtpy.QtWidgets import QLabel, QMessageBox, QTextEdit

from xenium_viewer.tabs._helpers import (
    StatusProxy, make_progress_bar, make_tab,
)
from xenium_viewer.utils import cache_repair
from xenium_viewer.utils.cache_repair import human_bytes

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def _cache_path(ctx) -> "Path | None":
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

    report_text = QTextEdit()
    report_text.setReadOnly(True)
    report_text.setFontFamily("monospace")
    report_text.setMinimumHeight(220)

    verify_btn = PushButton(label="Verify (read-only)")
    consolidate_btn = PushButton(label="Re-consolidate Metadata")
    backup_widget = ComboBox(label="Backup", choices=[])
    recover_btn = PushButton(label="Recover from Backup...", enabled=False)
    rebuild_btn = PushButton(label="Force Rebuild + Restore...")
    log_btn = PushButton(label="Open Log File")
    copy_btn = PushButton(label="Copy Report")

    mutating = (consolidate_btn, recover_btn, rebuild_btn)

    def _set_busy(busy: bool):
        progress.setVisible(busy)
        for button in (verify_btn, *mutating):
            button.enabled = not busy
        if not busy:
            recover_btn.enabled = bool(backup_widget.choices)

    def _no_cache_message() -> str:
        if ctx.no_cache:
            return ("Running with --no-cache: nothing is persisted, so there is "
                    "no cache to inspect.")
        return "No zarr cache found for this dataset."

    # ── Header: size, free space, manifest, failures this session ────────
    def _refresh_header(report=None, described=None):
        from xenium_viewer.utils.reporting import configured_log_path, failure_summary

        cache_path = _cache_path(ctx)
        if cache_path is None:
            header.setText(_no_cache_message())
            return
        lines = [f"<b>{cache_path}</b>"]
        if described:
            lines.append(
                f"Size {human_bytes(described['size_bytes'])} · "
                f"{human_bytes(described['free_bytes'])} free on this disk"
            )
        manifest = (report.manifest if report is not None
                    else cache_repair.read_manifest(cache_path))
        if manifest:
            built = manifest.get("built_at", "?")[:19].replace("T", " ")
            lines.append(f"Built {built} · spatialdata "
                         f"{manifest.get('spatialdata_version', '?')}")
        else:
            lines.append("No build manifest — freshness falls back to file "
                         "timestamps, which a copy also changes.")
        if report is not None:
            lines.append("Status: <b>healthy</b>" if report.ok
                         else "Status: <b>needs attention</b> (see below)")
        lines.append(failure_summary())
        log_file = configured_log_path()
        if log_file is not None:
            lines.append(f"Log: {log_file}")
        header.setText("<br>".join(lines))

    def _refresh_backups(report=None):
        cache_path = _cache_path(ctx)
        if cache_path is None:
            backup_widget.choices = []
            recover_btn.enabled = False
            return
        if report is None:
            report = cache_repair.verify(cache_path)
        options: list[str] = []
        state["_cache_backup_map"] = {}
        for element, copies in sorted(report.trash_available.items()):
            label = f"{element}  (previous version)"
            options.append(label)
            state["_cache_backup_map"][label] = ("trash", element, copies[0])
        for backup in report.backups:
            label = f"{backup.name}  (whole cache)"
            options.append(label)
            state["_cache_backup_map"][label] = ("cache", None, backup)
        backup_widget.choices = options
        recover_btn.enabled = bool(options)

    # ── Verify ───────────────────────────────────────────────────────────
    def _on_verify():
        cache_path = _cache_path(ctx)
        if cache_path is None:
            report_text.setPlainText(_no_cache_message())
            _refresh_header()
            return
        _set_busy(True)
        status.value = "Checking cache..."

        @thread_worker
        def _run():
            report = cache_repair.verify(cache_path)
            return report, cache_repair.describe_store(cache_path)

        def _done(result):
            report, described = result
            report_text.setPlainText(report.summary())
            _refresh_header(report, described)
            _refresh_backups(report)
            status.value = ("Cache is healthy." if report.ok
                            else "Cache needs attention — see the Cache tab.")
            _set_busy(False)

        _start(_run, _done, "Verify failed")

    # ── Re-consolidate ───────────────────────────────────────────────────
    def _on_consolidate():
        cache_path = _cache_path(ctx)
        if cache_path is None:
            return
        _set_busy(True)
        status.value = "Re-consolidating cache metadata..."

        @thread_worker
        def _run():
            report = cache_repair.verify(cache_path)
            result = cache_repair.repair(cache_path, report, level=cache_repair.AUTO)
            return result, cache_repair.verify(cache_path)

        def _done(payload):
            result, report = payload
            lines = ["Repair actions:"] + [f"  • {a}" for a in result.actions]
            if result.failures:
                lines += ["", "Could not fix:"] + [f"  • {f}" for f in result.failures]
            lines += ["", report.summary()]
            report_text.setPlainText("\n".join(lines))
            _refresh_header(report)
            _refresh_backups(report)
            status.value = ("Cache metadata re-consolidated."
                            if not result.failures else "Repair finished with warnings.")
            _set_busy(False)

        _start(_run, _done, "Repair failed")

    # ── Recover from a backup ────────────────────────────────────────────
    def _on_recover():
        cache_path = _cache_path(ctx)
        choice = backup_widget.value
        mapping = state.get("_cache_backup_map", {})
        if cache_path is None or not choice or choice not in mapping:
            return
        kind, element, source = mapping[choice]

        if kind == "trash":
            question = (f"Restore the previous version of {element} from a backup?\n\n"
                        "The current copy is kept aside, so this is reversible.")
        else:
            question = (
                f"Copy user-generated data out of {source.name} into the live cache?\n\n"
                "Elements already present in the live cache are left alone."
            )
        if QMessageBox.question(
            None, "Recover from Backup", question,
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel,
        ) != QMessageBox.Yes:
            return

        _set_busy(True)
        status.value = "Recovering from backup..."

        @thread_worker
        def _run():
            if kind == "trash":
                cache_repair.restore_from_trash(cache_path, element, source)
                return [f"restored {element}"]
            return _recover_whole_cache(ctx, cache_path, source)

        def _done(actions):
            report = cache_repair.verify(cache_path)
            report_text.setPlainText(
                "\n".join(["Recovered:"] + [f"  • {a}" for a in actions]
                          + ["", report.summary()])
            )
            _refresh_header(report)
            _refresh_backups(report)
            status.value = (f"Recovered {len(actions)} item(s). Reopen the dataset "
                            "to see restored data.")
            _set_busy(False)

        _start(_run, _done, "Recovery failed")

    # ── Force rebuild ────────────────────────────────────────────────────
    def _on_rebuild():
        cache_path = _cache_path(ctx)
        if cache_path is None:
            return
        described = cache_repair.describe_store(cache_path)
        size = described["size_bytes"]
        free = described["free_bytes"]
        warning = ""
        if free < size * 1.1:
            warning = (f"\n\nWARNING: only {human_bytes(free)} free, and the rebuild "
                       f"needs roughly {human_bytes(size)} alongside the copy it "
                       "keeps. It will refuse to start rather than fill the disk.")
        if QMessageBox.question(
            None, "Force Rebuild",
            f"Rebuild the cache from the raw Xenium files?\n\n"
            f"The current cache ({human_bytes(size)}) is moved aside, not deleted, "
            f"and user-generated data is restored into the new one.\n\n"
            f"This re-reads the whole dataset and can take a long time."
            f"{warning}\n\nThe viewer must be restarted afterwards.",
            QMessageBox.Yes | QMessageBox.Cancel, QMessageBox.Cancel,
        ) != QMessageBox.Yes:
            return

        # The rebuild replaces the object every layer and manager is holding, so
        # doing it live would leave the session pointing at a freed store.
        # Staging it for the next launch is both safer and simpler than trying
        # to rewire the running viewer.
        from datetime import datetime
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        moved = cache_path.with_name(f"sdata_cached_prev_{stamp}.zarr")
        try:
            import shutil
            shutil.move(str(cache_path), str(moved))
        except Exception as e:
            QMessageBox.critical(None, "Force Rebuild",
                                 f"Could not move the existing cache aside:\n{e}")
            return
        report_text.setPlainText(
            f"The cache has been moved to:\n  {moved}\n\n"
            "Restart the viewer on this dataset to rebuild it. You will be asked "
            "whether to restore your data — the moved cache is the source it "
            "restores from, and nothing has been deleted."
        )
        status.value = "Cache moved aside — restart the viewer to rebuild."
        _refresh_header()
        _refresh_backups()

    # ── Small actions ────────────────────────────────────────────────────
    def _on_open_log():
        from xenium_viewer.utils.reporting import open_log_file
        if not open_log_file():
            status.value = "No log file for this session."

    def _on_copy():
        try:
            from qtpy.QtWidgets import QApplication
            QApplication.clipboard().setText(report_text.toPlainText())
            status.value = "Report copied to clipboard."
        except Exception as e:
            status.value = f"Could not copy: {e}"

    # ── Worker plumbing ──────────────────────────────────────────────────
    def _start(make_worker, on_done, error_prefix):
        worker = make_worker()

        def _failed(exc_info):
            report_text.setPlainText(f"{error_prefix}: {exc_info[1]}")
            status.value = f"{error_prefix}: {exc_info[1]}"
            _set_busy(False)

        worker.returned.connect(on_done)
        worker.errored.connect(_failed)
        # Keep a reference so the worker is not garbage-collected mid-run.
        state["_cache_worker"] = worker
        worker.start()

    verify_btn.clicked.connect(_on_verify)
    consolidate_btn.clicked.connect(_on_consolidate)
    recover_btn.clicked.connect(_on_recover)
    rebuild_btn.clicked.connect(_on_rebuild)
    log_btn.clicked.connect(_on_open_log)
    copy_btn.clicked.connect(_on_copy)

    _refresh_header()
    _refresh_backups()
    if ctx.no_cache or _cache_path(ctx) is None:
        for button in (verify_btn, *mutating):
            button.enabled = False

    widget = make_tab(
        header,
        verify_btn,
        progress,
        report_text,
        consolidate_btn,
        backup_widget,
        recover_btn,
        rebuild_btn,
        log_btn,
        copy_btn,
    )

    def _restore_session(session):
        return None

    return widget, {"restore_session": _restore_session}


def _recover_whole_cache(ctx, cache_path: Path, backup: Path) -> list[str]:
    """Copy user-generated elements and sidecars out of *backup* into the live cache.

    Only fills gaps: an element already in the live cache is left alone, because
    it is the newer of the two.
    """
    import shutil
    import warnings

    import spatialdata

    from xenium_viewer.loader import _detect_user_data
    from xenium_viewer.utils.adata_persistence import sidecar_dir
    from xenium_viewer.utils.zarr_safe import safe_write_element

    actions: list[str] = []
    user_data = _detect_user_data(backup)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        old = spatialdata.read_zarr(str(backup))
        live = ctx.sdata

        for key in user_data["shapes"] + user_data["images"]:
            if key in live:
                continue
            try:
                safe_write_element(live, key, old[key])
                actions.append(key)
            except Exception as e:
                actions.append(f"{key} — FAILED: {e}")

    # Sidecars live beside the store now; copy any the live dataset lacks.
    destination = sidecar_dir(Path(ctx.data_path), create=True)
    for name in user_data["sidecars"]:
        for source in (backup / name, backup.parent / "viewer_cache" / name):
            if source.exists() and not (destination / name).exists():
                try:
                    shutil.copy2(source, destination / name)
                    actions.append(name)
                except Exception as e:
                    actions.append(f"{name} — FAILED: {e}")
                break

    if not actions:
        actions.append("nothing to recover — the live cache already has it all")
    return actions
