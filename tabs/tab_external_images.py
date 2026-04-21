"""
External Images tab — loads arbitrary multichannel OME-TIFF/TIFF/SVS files as
additional napari layers. Channels appear as individual sub-layers with
napari's native per-channel controls; this tab adds minimal group-level
widgets (group opacity, show/hide all, affine mirror, remove).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from qtpy.QtWidgets import (
    QCheckBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QListWidget,
    QListWidgetItem, QMessageBox, QPushButton, QSlider, QVBoxLayout, QWidget,
    QComboBox,
)
from qtpy.QtCore import Qt
from napari.qt.threading import thread_worker

from tabs._helpers import make_tab
from utils.registration import load_multichannel_pyramid
from utils.affine_linking import (
    link_affine, list_transformable_layers, find_layer_by_name,
)
from utils.adata_persistence import (
    save_external_image_to_sdata, load_external_images_from_sdata,
    save_overlay_affine_to_sdata, _slugify,
)

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


def build_tab(ctx: "ViewerContext"):
    viewer = ctx.viewer

    # ── UI ───────────────────────────────────────────────────────────────
    root = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)

    add_btn = QPushButton("Add image…")
    layout.addWidget(add_btn)

    list_widget = QListWidget()
    list_widget.setMinimumHeight(120)
    layout.addWidget(list_widget)

    # Per-selection panel
    panel = QGroupBox("Selected image")
    panel_layout = QVBoxLayout()

    status_label = QLabel("—")
    status_label.setWordWrap(True)
    panel_layout.addWidget(status_label)

    opacity_row = QHBoxLayout()
    opacity_row.addWidget(QLabel("Opacity:"))
    opacity_slider = QSlider(Qt.Horizontal)
    opacity_slider.setRange(0, 100)
    opacity_slider.setValue(100)
    opacity_row.addWidget(opacity_slider)
    panel_layout.addLayout(opacity_row)

    vis_row = QHBoxLayout()
    show_all_btn = QPushButton("Show all channels")
    hide_all_btn = QPushButton("Hide all channels")
    vis_row.addWidget(show_all_btn)
    vis_row.addWidget(hide_all_btn)
    panel_layout.addLayout(vis_row)

    affine_row = QHBoxLayout()
    affine_row.addWidget(QLabel("Apply transform from:"))
    affine_combo = QComboBox()
    affine_row.addWidget(affine_combo, 1)
    panel_layout.addLayout(affine_row)

    remove_btn = QPushButton("Remove image")
    panel_layout.addWidget(remove_btn)
    panel.setLayout(panel_layout)
    layout.addWidget(panel)

    root.setLayout(layout)
    tab_widget = make_tab(root)

    # ── Helpers ──────────────────────────────────────────────────────────

    def _current_entry():
        row = list_widget.currentRow()
        if row < 0 or row >= len(ctx.external_images_state):
            return None
        return ctx.external_images_state[row]

    def _refresh_affine_choices():
        entry = _current_entry()
        excluded = entry["sub_layer_refs"] if entry else ()
        layers = list_transformable_layers(viewer, exclude=excluded)
        affine_combo.blockSignals(True)
        affine_combo.clear()
        affine_combo.addItem("(none)", None)
        for lyr in layers:
            affine_combo.addItem(lyr.name, lyr.name)
        # Restore current selection if still present
        if entry and entry.get("affine_source_name"):
            idx = affine_combo.findData(entry["affine_source_name"])
            if idx >= 0:
                affine_combo.setCurrentIndex(idx)
        affine_combo.blockSignals(False)

    def _update_panel():
        entry = _current_entry()
        if entry is None:
            status_label.setText("—")
            opacity_slider.setEnabled(False)
            show_all_btn.setEnabled(False)
            hide_all_btn.setEnabled(False)
            affine_combo.setEnabled(False)
            remove_btn.setEnabled(False)
            return
        opacity_slider.setEnabled(True)
        show_all_btn.setEnabled(True)
        hide_all_btn.setEnabled(True)
        affine_combo.setEnabled(True)
        remove_btn.setEnabled(True)

        n = len(entry.get("sub_layer_refs", []))
        src_path = entry.get("path") or entry.get("element_name", "?")
        status_label.setText(f"{Path(src_path).name} — {n} channel(s)")

        opacity_slider.blockSignals(True)
        opacity_slider.setValue(int(entry.get("opacity", 1.0) * 100))
        opacity_slider.blockSignals(False)

        _refresh_affine_choices()

    def _apply_opacity(value: int):
        entry = _current_entry()
        if entry is None:
            return
        op = value / 100.0
        entry["opacity"] = op
        for lyr in entry["sub_layer_refs"]:
            try:
                lyr.opacity = op
            except Exception:
                pass

    def _set_visible(vis: bool):
        entry = _current_entry()
        if entry is None:
            return
        for lyr in entry["sub_layer_refs"]:
            try:
                lyr.visible = vis
            except Exception:
                pass

    def _apply_affine_source(_ix: int):
        entry = _current_entry()
        if entry is None:
            return
        name = affine_combo.currentData()
        # Disconnect previous
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
            entry["affine_disconnect"] = None
        entry["affine_source_name"] = name
        if not name:
            return
        source = find_layer_by_name(viewer, name)
        if source is None:
            return
        # Mirror affine onto every sub-layer — group a single disconnect
        disconnects = []
        for lyr in entry["sub_layer_refs"]:
            try:
                disconnects.append(link_affine(lyr, source, viewer=viewer))
            except Exception as e:
                print(f"  Warning: could not link external image affine: {e}")

        def _disconnect_all():
            for d in disconnects:
                try:
                    d()
                except Exception:
                    pass

        entry["affine_disconnect"] = _disconnect_all

        # Persist affine to sdata
        if entry["sub_layer_refs"]:
            try:
                save_overlay_affine_to_sdata(
                    ctx, entry["element_name"],
                    entry["sub_layer_refs"][0].affine.affine_matrix,
                )
            except Exception:
                pass

    def _build_layers(pyramid, channel_axis, channel_names, display_name,
                      opacity=1.0):
        """Create napari image layers for one external image. Returns refs list."""
        if channel_axis is None:
            lyr = viewer.add_image(
                pyramid,
                name=f"{display_name} :: {channel_names[0] if channel_names else 'C0'}",
                blending="additive",
                opacity=opacity,
            )
            return [lyr]
        try:
            layers = viewer.add_image(
                pyramid,
                channel_axis=channel_axis,
                name=[f"{display_name} :: {n}" for n in channel_names],
                blending="additive",
                opacity=opacity,
            )
        except TypeError:
            # Older napari: name must be single string
            layers = viewer.add_image(
                pyramid, channel_axis=channel_axis, name=display_name,
                blending="additive", opacity=opacity,
            )
        if not isinstance(layers, list):
            layers = [layers]
        return layers

    def _register_entry(entry):
        ctx.external_images_state.append(entry)
        item = QListWidgetItem(
            f"{Path(entry.get('path') or entry.get('element_name', '?')).name}"
            f"  ({len(entry['sub_layer_refs'])}ch)"
        )
        list_widget.addItem(item)
        list_widget.setCurrentRow(list_widget.count() - 1)

    def _on_loaded(result, display_name, path):
        pyramid, tif, channel_axis, channel_names = result
        layers = _build_layers(pyramid, channel_axis, channel_names, display_name)
        element_name = f"ext_{_slugify(Path(path).stem)}"
        entry = {
            "element_name": element_name,
            "path": str(path),
            "tif": tif,
            "pyramid": pyramid,
            "channel_axis": channel_axis,
            "channel_names": list(channel_names),
            "sub_layer_refs": layers,
            "affine_source_name": None,
            "affine_disconnect": None,
            "opacity": 1.0,
        }
        _register_entry(entry)
        save_external_image_to_sdata(
            ctx, element_name, pyramid, channel_axis, channel_names,
        )
        ctx.record_code(
            f"\n# Load multichannel image\n"
            f"from utils.registration import load_multichannel_pyramid\n"
            f"ext_pyramid, ext_tif, channel_axis, channel_names = "
            f"load_multichannel_pyramid(\"{path}\")"
        )
        ctx.set_status(f"Loaded {Path(path).name} ({len(channel_names)} channels)")

    def on_add():
        default_dir = str(ctx.data_path) if ctx.data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Add multichannel image", default_dir,
            "Image Files (*.ome.tif *.tif *.tiff *.svs);;All Files (*)",
        )
        if not path:
            return
        ctx.set_status(f"Loading {Path(path).name}…")
        gen = ctx.dataset_generation
        display_name = Path(path).stem

        @thread_worker
        def _task():
            return load_multichannel_pyramid(path)

        worker = _task()
        worker.returned.connect(
            lambda result: _on_loaded(result, display_name, path)
            if ctx.dataset_generation == gen else None
        )
        worker.errored.connect(
            lambda e: QMessageBox.critical(None, "Load failed", str(e))
        )
        worker.start()

    def on_remove():
        row = list_widget.currentRow()
        if row < 0 or row >= len(ctx.external_images_state):
            return
        entry = ctx.external_images_state.pop(row)
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
        for lyr in entry.get("sub_layer_refs", []):
            try:
                viewer.layers.remove(lyr)
            except Exception:
                pass
        try:
            tif = entry.get("tif")
            if tif is not None:
                tif.close()
        except Exception:
            pass
        # Remove from sdata too
        try:
            element = entry.get("element_name")
            if element and ctx.sdata is not None and element in ctx.sdata:
                ctx.sdata.delete_element_from_disk(element)
        except Exception as e:
            print(f"  Warning: could not delete {element} from sdata: {e}")
        list_widget.takeItem(row)
        _update_panel()

    # ── Wire signals ─────────────────────────────────────────────────────
    add_btn.clicked.connect(on_add)
    list_widget.currentRowChanged.connect(lambda _i: _update_panel())
    opacity_slider.valueChanged.connect(_apply_opacity)
    show_all_btn.clicked.connect(lambda: _set_visible(True))
    hide_all_btn.clicked.connect(lambda: _set_visible(False))
    affine_combo.currentIndexChanged.connect(_apply_affine_source)
    remove_btn.clicked.connect(on_remove)

    # Refresh affine choices whenever layers change
    def _on_layers_changed(_event=None):
        _refresh_affine_choices()

    try:
        viewer.layers.events.inserted.connect(_on_layers_changed)
        viewer.layers.events.removed.connect(_on_layers_changed)
    except Exception:
        pass

    _update_panel()

    # ── Session restore ──────────────────────────────────────────────────
    def restore_session(session):
        """Re-hydrate external images from sdata.images[ext_*]."""
        entries = load_external_images_from_sdata(ctx.sdata)
        ui_list = (session or {}).get("external_images_ui") or []
        ui_by_name = {u.get("element_name"): u for u in ui_list if isinstance(u, dict)}

        for meta in entries:
            element_name = meta["element_name"]
            pyramid = meta["pyramid"]
            channel_names = meta["channel_names"]
            # sdata-stored pyramids are in (c, y, x) layout.
            channel_axis = 0 if pyramid[0].ndim >= 3 else None
            ui = ui_by_name.get(element_name, {})
            display_name = Path(ui.get("path", element_name)).stem or element_name
            layers = _build_layers(
                pyramid, channel_axis, channel_names, display_name,
                opacity=float(ui.get("opacity", 1.0)),
            )
            entry = {
                "element_name": element_name,
                "path": ui.get("path") or element_name,
                "tif": None,
                "pyramid": pyramid,
                "channel_axis": channel_axis,
                "channel_names": channel_names,
                "sub_layer_refs": layers,
                "affine_source_name": ui.get("affine_source_name"),
                "affine_disconnect": None,
                "opacity": float(ui.get("opacity", 1.0)),
            }
            # Apply affine from sdata (authoritative), fall back to session attrs
            saved_affine = meta.get("affine_matrix")  # from sdata transformations
            if saved_affine is None:
                saved_affine = ui.get("affine_matrix")  # fallback: session attrs
            if saved_affine is not None:
                for lyr in layers:
                    try:
                        lyr.affine = np.array(saved_affine, dtype=np.float64)
                    except Exception:
                        pass

            _register_entry(entry)
            # Re-link affine if named source is present
            if entry["affine_source_name"]:
                src = find_layer_by_name(viewer, entry["affine_source_name"])
                if src is not None:
                    disconnects = []
                    for lyr in entry["sub_layer_refs"]:
                        try:
                            disconnects.append(link_affine(lyr, src, viewer=viewer))
                        except Exception:
                            pass
                    if disconnects:
                        def _disc_all(dd=disconnects):
                            for d in dd:
                                try: d()
                                except Exception: pass
                        entry["affine_disconnect"] = _disc_all

        # Deferred affine linking: source layers (e.g. H&E) may not exist
        # yet because they load asynchronously. Watch for layer insertions
        # and link once the named source appears.
        _pending = [
            e for e in ctx.external_images_state
            if e.get("affine_source_name") and e.get("affine_disconnect") is None
        ]

        def _on_layer_inserted_for_pending(event=None):
            still_pending = []
            for e in _pending:
                if e.get("affine_disconnect") is not None:
                    continue
                src = find_layer_by_name(viewer, e["affine_source_name"])
                if src is None:
                    still_pending.append(e)
                    continue
                disconnects = []
                for lyr in e.get("sub_layer_refs", []):
                    try:
                        disconnects.append(link_affine(lyr, src, viewer=viewer))
                    except Exception:
                        pass
                if disconnects:
                    def _disc_all(dd=disconnects):
                        for d in dd:
                            try: d()
                            except Exception: pass
                    e["affine_disconnect"] = _disc_all
                    # Persist to sdata
                    refs = e.get("sub_layer_refs", [])
                    if refs:
                        try:
                            save_overlay_affine_to_sdata(
                                ctx, e["element_name"],
                                refs[0].affine.affine_matrix,
                            )
                        except Exception:
                            pass
            _pending[:] = still_pending
            if not _pending:
                try:
                    viewer.layers.events.inserted.disconnect(
                        _on_layer_inserted_for_pending
                    )
                except Exception:
                    pass

        if _pending:
            try:
                viewer.layers.events.inserted.connect(
                    _on_layer_inserted_for_pending
                )
            except Exception:
                pass

    return tab_widget, {"restore_session": restore_session}
