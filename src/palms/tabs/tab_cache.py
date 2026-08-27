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

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton
from napari.qt.threading import thread_worker
from qtpy.QtWidgets import QLabel, QMessageBox, QTextEdit

from palms.tabs._helpers import (
    StatusProxy, make_progress_bar, make_tab,
)
from palms.utils import cache_repair
from palms.utils.cache_repair import human_bytes

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def _cache_path(ctx) -> "Path | None":
    if ctx.sdata is not None and getattr(ctx.sdata, "path", None):
        return Path(ctx.sdata.path)
    if ctx.data_path is not None:
        candidate = Path(ctx.data_path) / "sdata_cached.zarr"
        if candidate.exists():
            return candidate
    return None


def _rebuild_blocked_reason(ctx) -> "str | None":
    """Why Force Rebuild must not run here, or ``None`` if it may.

    Pure and Qt-free so it can be tested: ``_on_rebuild`` itself opens a modal,
    which is why the precondition that used to be missing lives out here rather
    than inline. A Crop Dataset export has no raw 10x files, so the rebuild it
    stages for the next launch cannot happen — and the staging step renames the
    only copy of the data aside first.
    """
    from palms.loader import has_raw_xenium_source

    data_path = getattr(ctx, "data_path", None)
    if data_path is None:
        return None
    if has_raw_xenium_source(Path(data_path)):
        return None
    return (
        "This dataset has no raw Xenium output (no cells.zarr.zip, no "
        "cell_feature_matrix, no morphology_focus) — most likely a Crop Dataset "
        "export. Its zarr cache is the only copy of the data, so there is "
        "nothing to rebuild from. Use Verify, Re-consolidate or Recover from "
        "Backup instead."
    )


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

    rebuild_blocked = _rebuild_blocked_reason(ctx)
    if rebuild_blocked:
        rebuild_btn.enabled = False
        rebuild_btn.tooltip = rebuild_blocked

    def _set_busy(busy: bool):
        progress.setVisible(busy)
        for button in (verify_btn, *mutating):
            button.enabled = not busy
        if not busy:
            recover_btn.enabled = bool(backup_widget.choices)
            # Re-enabling everything after a worker finishes must not undo the
            # precondition — this is the button that renames the cache aside.
            if rebuild_blocked:
                rebuild_btn.enabled = False

    def _no_cache_message() -> str:
        if ctx.no_cache:
            return ("Running with --no-cache: nothing is persisted, so there is "
                    "no cache to inspect.")
        return "No zarr cache found for this dataset."

    # ── Header: size, free space, manifest, failures this session ────────
    def _refresh_header(report=None, described=None):
        from palms.utils.reporting import configured_log_path, failure_summary

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
        elif not rebuild_blocked:
            lines.append("No build manifest — freshness falls back to file "
                         "timestamps, which a copy also changes.")
        if rebuild_blocked:
            # Worth saying up front rather than only on a disabled button: for
            # this dataset the cache is the data, not a derivative of it.
            lines.append("<b>Cache-only dataset</b> — no raw Xenium output "
                         "beside it, so this cache is the only copy of the data "
                         "and is never rebuilt.")
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
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, QMessageBox.StandardButton.Cancel,
        ) != QMessageBox.StandardButton.Yes:
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
            _set_busy(False)

            recovered = [a for a in actions if "FAILED" not in a
                         and not a.startswith("nothing to recover")]
            if not recovered:
                status.value = "Nothing was recovered."
                return
            status.value = f"Recovered {len(recovered)} item(s)."
            _offer_reload(len(recovered))

        _start(_run, _done, "Recovery failed")

    # ── Force rebuild ────────────────────────────────────────────────────
    def _on_rebuild():
        cache_path = _cache_path(ctx)
        if cache_path is None:
            return
        blocked = _rebuild_blocked_reason(ctx)
        if blocked:
            # The button is already disabled; this is the guard that matters,
            # since what follows renames the store before anything checks that
            # the rebuild it stages is even possible.
            QMessageBox.information(None, "Force Rebuild", blocked)
            report_text.setPlainText(blocked)
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
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, QMessageBox.StandardButton.Cancel,
        ) != QMessageBox.StandardButton.Yes:
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

    def _offer_reload(n: int):
        """Recovered data is on disk but not in memory — offer to reload.

        The elements were written straight into the zarr store, so the live
        SpatialData, the napari layers and every tab's widgets know nothing
        about them. Reloading rebuilds all of it from disk.
        """
        reload_dataset = getattr(ctx, "reload_dataset", None)
        if reload_dataset is None:
            report_text.append(
                "\n\nRecovered data is on disk. Reopen this dataset "
                "(File → Open Dataset) to see it."
            )
            return
        if QMessageBox.question(
            None, "Reload Dataset",
            f"Recovered {n} item(s) into the cache.\n\n"
            "They are on disk but not yet loaded — reload the dataset now to "
            "see them?\n\nThis rebuilds the layers and tabs from disk and "
            "takes a few seconds.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.Yes,
        ) != QMessageBox.StandardButton.Yes:
            report_text.append(
                "\n\nRecovered data is on disk. Reopen this dataset "
                "(File → Open Dataset) when you want to see it."
            )
            return
        status.value = "Reloading dataset..."
        # Deliberately synchronous: this tears down and rebuilds every widget,
        # including the one this callback belongs to, so it must not run while
        # a worker holds a reference to the old tab.
        reload_dataset()

    # ── Small actions ────────────────────────────────────────────────────
    def _on_open_log():
        from palms.utils.reporting import open_log_file
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

        def _failed(exc):
            # napari's `errored` emits the exception itself, not an exc_info
            # triple. Indexing it raised TypeError, which then replaced the real
            # error in the traceback.
            import traceback
            detail = "".join(traceback.format_exception(
                type(exc), exc, exc.__traceback__))
            report_text.setPlainText(f"{error_prefix}: {exc}\n\n{detail}")
            status.value = f"{error_prefix}: {exc}"
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
    """Salvage user-generated data out of *backup* into the live cache.

    Works entirely at the filesystem level and **never opens the backup as a
    SpatialData**. A cache worth recovering from is one that failed to open, so
    anything that starts with ``read_zarr`` on it cannot work: the first attempt
    at this did exactly that and died in ``_read_table`` on a corrupt table,
    taking the recoverable shapes and images down with it.

    Element directories are self-contained, and obs columns are individual zarr
    arrays, so a broken root index or an unreadable table does not condemn what
    sits beside it. Only gaps are filled — anything already in the live cache is
    the newer copy and is left alone.
    """
    from palms.loader import (
        _USER_IMAGE_KEYS, _USER_OBS_PREFIXES, _USER_SHAPE_KEYS, _is_user_element,
    )
    from palms.utils.adata_persistence import sidecar_dir
    from palms.utils.zarr_safe import safe_import_element

    actions: list[str] = []

    # ── Whole elements: ROIs, annotations, landmarks, tiles, images ──────
    for element in cache_repair.salvageable_elements(backup):
        etype, _, name = element.partition("/")
        keys = _USER_IMAGE_KEYS if etype == "images" else _USER_SHAPE_KEYS
        if etype not in ("shapes", "images") or not _is_user_element(name, keys):
            continue
        if (cache_path / element).exists():
            continue                      # live copy is newer
        try:
            safe_import_element(cache_path, element, backup / element)
            actions.append(element)
        except Exception as e:
            actions.append(f"{element} — FAILED: {e}")

    # ── Clusterings and CNV scores, read column by column ────────────────
    # The most valuable thing in a broken cache, and the part a whole-store
    # read can never reach.
    actions += _recover_obs_columns(ctx, backup, _USER_OBS_PREFIXES)

    # ── Registration state ───────────────────────────────────────────────
    # Recovering he_image and its landmarks without this is half a job: the
    # image element would exist while the session still says no H&E is loaded,
    # so nothing would display it.
    actions += _recover_session_attrs(ctx, cache_path, backup)

    # ── Sidecars (CNV caches, DEG tables) ────────────────────────────────
    destination = sidecar_dir(Path(ctx.data_path), create=True)
    for pattern in ("*.h5ad", "*.parquet", "cnv_*_result.json"):
        for source in sorted(list(backup.glob(pattern))
                             + list((backup.parent / "viewer_cache").glob(pattern))):
            target = destination / source.name
            if target.exists():
                continue
            try:
                shutil.copy2(source, target)
                actions.append(source.name)
            except Exception as e:
                actions.append(f"{source.name} — FAILED: {e}")

    if not actions:
        actions.append("nothing to recover — the live cache already has it all")
    return actions


_SESSION_KEYS = (
    "he_filename", "he_path", "he_shape_yx", "flip_v", "flip_h",
    "arms_he_filename", "arms_he_path", "arms_he_shape_yx",
    "arms_affine_3x3", "arms_flip_v", "arms_flip_h",
    "arms_geojson_path", "arms_csv_path",
    "cluster_labels", "marker_genes_json", "prov_graph",
    "external_images_ui", "patch_overlays_ui",
)


def _recover_session_attrs(ctx, cache_path: Path, backup: Path) -> list[str]:
    """Merge registration and UI state from a backup's viewer_session.

    Writes to disk **and** hydrates ``ctx.he_state`` / ``ctx.arms_state``. The
    in-memory half is not optional: reloading the dataset saves the current
    session first, and ``save_session`` deletes the ``he``/``arms`` groups and
    rewrites them from ``he_state``. With an empty ``he_state`` that erased the
    affine we had just recovered — the images came back but unaligned, which is
    exactly the symptom that exposed this.

    Only fills keys the live session lacks. Uses plain zarr on both sides, so a
    backup whose *store* is unreadable still gives up its session.
    """
    import numpy as np
    import zarr

    from palms.utils.zarr_safe import safe_group_update

    source = backup / "viewer_session"
    if not source.is_dir():
        return []
    try:
        old = zarr.open_group(str(source), mode="r", use_consolidated=False)
        old_attrs = dict(old.attrs)
    except Exception as e:
        return [f"viewer_session — FAILED: {e}"]

    live_group = cache_path / "viewer_session"
    live_attrs: dict = {}
    if live_group.is_dir():
        try:
            live_attrs = dict(zarr.open_group(
                str(live_group), mode="r", use_consolidated=False).attrs)
        except Exception:
            live_attrs = {}

    def _empty(value) -> bool:
        return value is None or (isinstance(value, (list, dict, str)) and len(value) == 0)

    missing = {k: old_attrs[k] for k in _SESSION_KEYS
               if not _empty(old_attrs.get(k)) and _empty(live_attrs.get(k))}
    affines = [name for name in ("he", "arms")
               if (source / name).is_dir() and not (live_group / name / "affine_3x3").exists()]
    if not missing and not affines:
        return []

    actions: list[str] = []
    try:
        with safe_group_update(cache_path, "viewer_session") as (session, stage):
            for name in affines:
                target = stage / name
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source / name, target)
                actions.append(f"{name} registration affine")
            if missing:
                session.attrs.update(missing)
                actions.append(f"session state ({', '.join(sorted(missing))})")
    except Exception as e:
        return [f"viewer_session — FAILED: {e}"]

    # ── Hydrate the live state, or the reload's save will undo all of it ──
    def _array(group_name: str, array_name: str):
        path = source / group_name / array_name
        if not path.exists():
            return None
        try:
            return np.array(zarr.open_array(str(path), mode="r")[:], dtype=np.float64)
        except Exception:
            return None

    he_state = getattr(ctx, "he_state", None)
    if isinstance(he_state, dict):
        for key, attr in (("he_filename", "he_filename"), ("he_path", "he_path"),
                          ("he_shape_yx", "he_shape_yx"),
                          ("flip_v", "flip_v"), ("flip_h", "flip_h")):
            if _empty(he_state.get(key)) and not _empty(old_attrs.get(attr)):
                he_state[key] = old_attrs[attr]
        for key, array_name in (("affine_3x3", "affine_3x3"),
                                ("coarse_affine", "coarse_affine")):
            if he_state.get(key) is None:
                recovered = _array("he", array_name)
                if recovered is not None:
                    he_state[key] = recovered

    arms_state = getattr(ctx, "arms_state", None)
    if isinstance(arms_state, dict):
        for key, attr in (("he_filename", "arms_he_filename"), ("he_path", "arms_he_path"),
                          ("he_shape_yx", "arms_he_shape_yx"),
                          ("flip_v", "arms_flip_v"), ("flip_h", "arms_flip_h"),
                          ("geojson_path", "arms_geojson_path"),
                          ("csv_path", "arms_csv_path")):
            if _empty(arms_state.get(key)) and not _empty(old_attrs.get(attr)):
                arms_state[key] = old_attrs[attr]
        if arms_state.get("affine_3x3") is None:
            recovered = _array("arms", "affine_3x3")
            if recovered is None and not _empty(old_attrs.get("arms_affine_3x3")):
                recovered = np.asarray(old_attrs["arms_affine_3x3"], dtype=np.float64)
            if recovered is not None:
                arms_state["affine_3x3"] = recovered

    return actions


def _recover_obs_columns(ctx, backup: Path, prefixes: tuple) -> list[str]:
    """Pull clustering / CNV obs columns out of a backup table into the live one."""
    import pandas as pd

    from palms.utils.adata_persistence import _persist_table

    live_adata = getattr(ctx, "adata", None)
    if live_adata is None:
        return []

    columns = cache_repair.read_obs_columns(backup, prefixes)
    if not columns:
        return []

    live_index = (live_adata.obs["cell_id"].astype(str).to_numpy()
                  if "cell_id" in live_adata.obs.columns
                  else live_adata.obs_names.to_numpy().astype(str))

    actions: list[str] = []
    added = False
    for name, (index, values) in columns.items():
        if name in live_adata.obs.columns:
            continue
        try:
            series = pd.Series(list(values), index=[str(i) for i in index])
            aligned = series.reindex(live_index)
            if aligned.isna().all():
                actions.append(f"{name} — skipped (no matching cells)")
                continue
            live_adata.obs[name] = pd.Categorical(aligned.astype(str).values)
            actions.append(f"obs/{name}")
            added = True
        except Exception as e:
            actions.append(f"obs/{name} — FAILED: {e}")

    if added:
        _persist_table(ctx)
    return actions
