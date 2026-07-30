"""
External Images tab — loads arbitrary multichannel OME-TIFF/TIFF/SVS files as
a single RGB composite napari layer per image. Per-channel visibility, color,
and contrast are controlled from the tab widget. Includes landmark-based
registration for independent alignment with the Xenium image.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from qtpy.QtWidgets import (
    QCheckBox, QColorDialog, QFileDialog, QGridLayout, QGroupBox, QHBoxLayout,
    QLabel, QListWidget, QListWidgetItem, QMessageBox, QPushButton, QScrollArea,
    QSlider, QVBoxLayout, QWidget, QComboBox,
)
from qtpy.QtCore import Qt, QTimer
from qtpy.QtGui import QColor
from superqt import QDoubleRangeSlider
from napari.qt.threading import thread_worker

from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.zarr_safe import safe_delete_element
from xenium_viewer.tabs._helpers import make_tab
from xenium_viewer.utils.registration import load_multichannel_pyramid, compute_landmark_affine
from xenium_viewer.utils.composite import (
    build_composite_pyramid, default_channel_colors, auto_contrast,
)
from xenium_viewer.utils.affine_linking import (
    link_affine, list_transformable_layers, find_layer_by_name,
)
from xenium_viewer.utils.adata_persistence import (
    save_external_image_to_sdata, load_external_images_from_sdata,
    save_overlay_affine_to_sdata, save_landmarks_to_sdata,
    load_landmarks_from_sdata, _slugify,
)

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: "ViewerContext"):
    viewer = ctx.viewer

    # ── UI ───────────────────────────────────────────────────────────────
    root = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)

    add_btn = QPushButton("Add image…")
    layout.addWidget(add_btn)

    list_widget = QListWidget()
    list_widget.setMinimumHeight(80)
    layout.addWidget(list_widget)

    # Per-selection panel
    panel = QGroupBox("Selected image")
    panel_layout = QVBoxLayout()

    status_label = QLabel("—")
    status_label.setWordWrap(True)
    panel_layout.addWidget(status_label)

    # Opacity
    opacity_row = QHBoxLayout()
    opacity_row.addWidget(QLabel("Opacity:"))
    opacity_slider = QSlider(Qt.Horizontal)
    opacity_slider.setRange(0, 100)
    opacity_slider.setValue(100)
    opacity_row.addWidget(opacity_slider)
    panel_layout.addLayout(opacity_row)

    # Channel controls (scrollable)
    ch_group = QGroupBox("Channels")
    ch_group_layout = QVBoxLayout()
    ch_vis_row = QHBoxLayout()
    show_all_btn = QPushButton("All on")
    hide_all_btn = QPushButton("All off")
    ch_vis_row.addWidget(show_all_btn)
    ch_vis_row.addWidget(hide_all_btn)
    ch_group_layout.addLayout(ch_vis_row)
    ch_scroll = QScrollArea()
    ch_scroll.setWidgetResizable(True)
    ch_scroll.setMinimumHeight(150)
    ch_scroll.setMaximumHeight(250)
    ch_scroll_inner = QWidget()
    ch_grid = QGridLayout()
    ch_grid.setContentsMargins(2, 2, 2, 2)
    ch_scroll_inner.setLayout(ch_grid)
    ch_scroll.setWidget(ch_scroll_inner)
    ch_group_layout.addWidget(ch_scroll)
    ch_group.setLayout(ch_group_layout)
    panel_layout.addWidget(ch_group)

    # Registration section
    reg_group = QGroupBox("Registration")
    reg_layout = QVBoxLayout()
    flip_row = QHBoxLayout()
    flip_v_chk = QCheckBox("Flip V")
    flip_h_chk = QCheckBox("Flip H")
    flip_row.addWidget(flip_v_chk)
    flip_row.addWidget(flip_h_chk)
    reg_layout.addLayout(flip_row)
    lm_row = QHBoxLayout()
    add_xen_lm_btn = QPushButton("Add Xenium LM")
    add_img_lm_btn = QPushButton("Add Image LM")
    clear_lm_btn = QPushButton("Clear")
    lm_row.addWidget(add_xen_lm_btn)
    lm_row.addWidget(add_img_lm_btn)
    lm_row.addWidget(clear_lm_btn)
    reg_layout.addLayout(lm_row)
    register_btn = QPushButton("Compute Registration")
    register_btn.setEnabled(False)
    reg_layout.addWidget(register_btn)
    reg_status_label = QLabel("")
    reg_status_label.setWordWrap(True)
    reg_layout.addWidget(reg_status_label)
    # Alternative: mirror existing affine
    affine_row = QHBoxLayout()
    affine_row.addWidget(QLabel("Or apply transform from:"))
    affine_combo = QComboBox()
    affine_row.addWidget(affine_combo, 1)
    reg_layout.addLayout(affine_row)
    reg_group.setLayout(reg_layout)
    panel_layout.addWidget(reg_group)

    remove_btn = QPushButton("Remove image")
    panel_layout.addWidget(remove_btn)
    panel.setLayout(panel_layout)
    layout.addWidget(panel)

    root.setLayout(layout)
    tab_widget = make_tab(root)

    # Debounce timer for contrast slider changes
    _debounce_timer = QTimer()
    _debounce_timer.setSingleShot(True)
    _debounce_timer.setInterval(100)

    # ── Helpers ──────────────────────────────────────────────────────────

    def _current_entry():
        row = list_widget.currentRow()
        if row < 0 or row >= len(ctx.external_images_state):
            return None
        return ctx.external_images_state[row]

    def _refresh_affine_choices():
        entry = _current_entry()
        lyr = entry["layer_ref"] if entry else None
        excluded = [lyr] if lyr else []
        # Also exclude landmark layers
        if entry:
            for k in ("xenium_lm_layer", "image_lm_layer"):
                l = entry.get(k)
                if l is not None:
                    excluded.append(l)
        layers = list_transformable_layers(viewer, exclude=excluded)
        affine_combo.blockSignals(True)
        affine_combo.clear()
        affine_combo.addItem("(none)", None)
        for lyr in layers:
            affine_combo.addItem(lyr.name, lyr.name)
        if entry and entry.get("affine_source_name"):
            idx = affine_combo.findData(entry["affine_source_name"])
            if idx >= 0:
                affine_combo.setCurrentIndex(idx)
        affine_combo.blockSignals(False)

    # ── Composite rendering ─────────────────────────────────────────────

    def _update_composite(entry):
        """Rebuild the dask composite and push to napari layer."""
        if entry is None or entry.get("layer_ref") is None:
            return
        pyramid = entry["pyramid"]
        ch_axis = entry["channel_axis"]
        if ch_axis is None:
            return  # single-channel, nothing to composite
        comp = build_composite_pyramid(pyramid, entry["channel_states"],
                                       channel_axis=ch_axis)
        entry["layer_ref"].data = comp

    def _build_composite_layer(pyramid, channel_axis, channel_states,
                               display_name, opacity=1.0):
        """Create a single RGB composite napari layer. Returns layer ref."""
        if channel_axis is None:
            # Single-channel — just add directly
            lyr = viewer.add_image(
                pyramid, name=display_name,
                blending="translucent", opacity=opacity,
            )
            return lyr
        comp = build_composite_pyramid(pyramid, channel_states,
                                       channel_axis=channel_axis)
        lyr = viewer.add_image(
            comp, rgb=True, name=display_name,
            blending="translucent", opacity=opacity,
        )
        return lyr

    # ── Channel control grid ────────────────────────────────────────────

    def _clear_channel_grid():
        """Remove all widgets from the channel grid layout."""
        while ch_grid.count():
            item = ch_grid.takeAt(0)
            w = item.widget()
            if w:
                w.setParent(None)
                w.deleteLater()

    def _rebuild_channel_grid(entry):
        """Populate channel grid with checkboxes, color buttons, sliders."""
        _clear_channel_grid()
        if entry is None:
            return
        states = entry.get("channel_states")
        if not states:
            return

        ch_widgets = []
        for i, cs in enumerate(states):
            row = i
            # Checkbox
            cb = QCheckBox(entry["channel_names"][i] if i < len(entry["channel_names"]) else f"C{i}")
            cb.setChecked(cs.get("visible", True))
            ch_grid.addWidget(cb, row, 0)

            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(20, 20)
            r, g, b = [int(c * 255) for c in cs["color"][:3]]
            color_btn.setStyleSheet(f"background-color: rgb({r},{g},{b}); border: 1px solid gray;")
            ch_grid.addWidget(color_btn, row, 1)

            # Contrast range slider
            lo, hi = cs["clim"]
            data_min = cs.get("data_min", 0.0)
            data_max = cs.get("data_max", 65535.0)
            range_slider = QDoubleRangeSlider(Qt.Horizontal)
            range_slider.setRange(data_min, data_max)
            range_slider.setValue((lo, hi))
            ch_grid.addWidget(range_slider, row, 2)

            ch_widgets.append((cb, color_btn, range_slider))

            # Wire signals
            def _on_vis_toggled(checked, idx=i):
                e = _current_entry()
                if e is not None and idx < len(e["channel_states"]):
                    e["channel_states"][idx]["visible"] = checked
                    _update_composite(e)

            def _on_color_clicked(_, idx=i, btn=color_btn):
                e = _current_entry()
                if e is None or idx >= len(e["channel_states"]):
                    return
                cur = e["channel_states"][idx]["color"]
                initial = QColor(int(cur[0]*255), int(cur[1]*255), int(cur[2]*255))
                color = QColorDialog.getColor(initial, None, f"Channel {idx} color")
                if color.isValid():
                    e["channel_states"][idx]["color"] = [
                        color.redF(), color.greenF(), color.blueF(),
                    ]
                    btn.setStyleSheet(
                        f"background-color: rgb({color.red()},{color.green()},{color.blue()}); "
                        f"border: 1px solid gray;"
                    )
                    _update_composite(e)

            def _on_contrast_changed(value, idx=i):
                e = _current_entry()
                if e is not None and idx < len(e["channel_states"]):
                    e["channel_states"][idx]["clim"] = [value[0], value[1]]
                    # Debounce — disconnect previous timeout handler safely
                    _debounce_timer.stop()
                    try:
                        _debounce_timer.timeout.disconnect()
                    except (TypeError, RuntimeError):
                        pass
                    _debounce_timer.timeout.connect(lambda: _update_composite(e))
                    _debounce_timer.start()

            cb.toggled.connect(_on_vis_toggled)
            color_btn.clicked.connect(_on_color_clicked)
            range_slider.valueChanged.connect(_on_contrast_changed)

        entry["_ch_widgets"] = ch_widgets

    # ── Affine / registration ───────────────────────────────────────────

    def _build_flip_affine(entry):
        shape = entry.get("image_shape_yx")
        if shape is None:
            return np.eye(3)
        h, w = shape
        M = np.eye(3, dtype=np.float64)
        if entry.get("flip_v", False):
            M = np.array([[-1, 0, h - 1], [0, 1, 0], [0, 0, 1]], dtype=np.float64) @ M
        if entry.get("flip_h", False):
            M = np.array([[1, 0, 0], [0, -1, w - 1], [0, 0, 1]], dtype=np.float64) @ M
        return M

    def _apply_entry_affine(entry):
        """Combine fine + flip and apply to layer + landmark layer."""
        flip = _build_flip_affine(entry)
        fine = entry.get("affine_3x3")
        combined = fine @ flip if fine is not None else flip
        lyr = entry.get("layer_ref")
        if lyr is not None:
            lyr.affine = combined
        img_lm = entry.get("image_lm_layer")
        if img_lm is not None:
            img_lm.affine = combined
        # Persist to sdata
        save_overlay_affine_to_sdata(ctx, entry["element_name"], combined)

    def _create_landmark_layers(entry):
        """Create Xenium + Image landmark point layers for registration."""
        if entry.get("xenium_lm_layer") is not None:
            return
        display = Path(entry.get("path") or entry.get("element_name", "Ext")).stem
        xen_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name=f"{display} Xenium LM", size=30, face_color="cyan",
            symbol="cross", border_color="cyan",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )
        img_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name=f"{display} Image LM", size=30, face_color="red",
            symbol="cross", border_color="red",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )

        def _check_lm_count(*_args):
            n = min(len(xen_lm.data), len(img_lm.data))
            register_btn.setEnabled(n >= 3)

        xen_lm.events.data.connect(_check_lm_count)
        img_lm.events.data.connect(_check_lm_count)
        entry["xenium_lm_layer"] = xen_lm
        entry["image_lm_layer"] = img_lm
        # Apply current affine to image landmark layer
        lyr = entry.get("layer_ref")
        if lyr is not None:
            try:
                img_lm.affine = lyr.affine
            except Exception:
                pass

    def _remove_landmark_layers(entry):
        for key in ("xenium_lm_layer", "image_lm_layer"):
            lyr = entry.get(key)
            if lyr is not None:
                try:
                    viewer.layers.remove(lyr)
                except Exception:
                    pass
            entry[key] = None

    # ── Panel update ────────────────────────────────────────────────────

    def _update_panel():
        entry = _current_entry()
        has_entry = entry is not None
        for w in (opacity_slider, show_all_btn, hide_all_btn, affine_combo,
                  remove_btn, flip_v_chk, flip_h_chk, add_xen_lm_btn,
                  add_img_lm_btn, clear_lm_btn):
            w.setEnabled(has_entry)
        register_btn.setEnabled(False)
        if not has_entry:
            status_label.setText("—")
            reg_status_label.setText("")
            _clear_channel_grid()
            return

        n = len(entry.get("channel_names", []))
        src_path = entry.get("path") or entry.get("element_name", "?")
        status_label.setText(f"{Path(src_path).name} — {n} channel(s)")

        opacity_slider.blockSignals(True)
        opacity_slider.setValue(int(entry.get("opacity", 1.0) * 100))
        opacity_slider.blockSignals(False)

        flip_v_chk.blockSignals(True)
        flip_v_chk.setChecked(entry.get("flip_v", False))
        flip_v_chk.blockSignals(False)
        flip_h_chk.blockSignals(True)
        flip_h_chk.setChecked(entry.get("flip_h", False))
        flip_h_chk.blockSignals(False)

        # Check landmark count
        xen_lm = entry.get("xenium_lm_layer")
        img_lm = entry.get("image_lm_layer")
        if xen_lm is not None and img_lm is not None:
            n_lm = min(len(xen_lm.data), len(img_lm.data))
            register_btn.setEnabled(n_lm >= 3)

        _rebuild_channel_grid(entry)
        _refresh_affine_choices()

    # ── Signal handlers ─────────────────────────────────────────────────

    def _apply_opacity(value: int):
        entry = _current_entry()
        if entry is None:
            return
        op = value / 100.0
        entry["opacity"] = op
        lyr = entry.get("layer_ref")
        if lyr is not None:
            lyr.opacity = op

    def _set_all_channels(vis: bool):
        entry = _current_entry()
        if entry is None:
            return
        for cs in entry.get("channel_states", []):
            cs["visible"] = vis
        _update_composite(entry)
        # Update checkboxes
        widgets = entry.get("_ch_widgets", [])
        for cb, _, _ in widgets:
            cb.blockSignals(True)
            cb.setChecked(vis)
            cb.blockSignals(False)

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
        lyr = entry.get("layer_ref")
        if lyr is None:
            return
        try:
            disconnect = link_affine(lyr, source, viewer=viewer)
            entry["affine_disconnect"] = disconnect
            # Also link image landmark layer if present
            img_lm = entry.get("image_lm_layer")
            if img_lm is not None:
                img_lm.affine = lyr.affine
            save_overlay_affine_to_sdata(
                ctx, entry["element_name"], lyr.affine.affine_matrix,
            )
        except Exception as e:
            print(f"  Warning: could not link external image affine: {e}")

    def on_flip_changed(_value=None):
        entry = _current_entry()
        if entry is None:
            return
        entry["flip_v"] = flip_v_chk.isChecked()
        entry["flip_h"] = flip_h_chk.isChecked()
        _apply_entry_affine(entry)

    def on_add_xenium_lm():
        entry = _current_entry()
        if entry is None:
            return
        _create_landmark_layers(entry)
        lm = entry["xenium_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.setText("Click on a feature in the Xenium image")

    def on_add_image_lm():
        entry = _current_entry()
        if entry is None:
            return
        _create_landmark_layers(entry)
        lm = entry["image_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.setText("Click on the same feature in the external image")

    def on_clear_lm():
        entry = _current_entry()
        if entry is None:
            return
        for key in ("xenium_lm_layer", "image_lm_layer"):
            lm = entry.get(key)
            if lm is not None:
                lm.selected_data = set()
                lm.data = np.empty((0, 2), dtype=np.float64)
        entry["affine_3x3"] = None
        _apply_entry_affine(entry)
        register_btn.setEnabled(False)
        reg_status_label.setText("Landmarks cleared")
        # Clear from sdata
        ename = entry["element_name"]
        save_landmarks_to_sdata(ctx, f"{ename}_xenium_lm", None)
        save_landmarks_to_sdata(ctx, f"{ename}_image_lm", None)

    def on_register():
        entry = _current_entry()
        if entry is None:
            return
        xen_lm = entry.get("xenium_lm_layer")
        img_lm = entry.get("image_lm_layer")
        if xen_lm is None or img_lm is None:
            return
        xen_pts = np.asarray(xen_lm.data, dtype=np.float64)
        img_pts = np.asarray(img_lm.data, dtype=np.float64)
        n = min(len(xen_pts), len(img_pts))
        if n < 3:
            reg_status_label.setText("Need at least 3 paired landmarks")
            return
        xen_pts = xen_pts[:n]
        img_pts = img_pts[:n]
        affine, residuals = compute_landmark_affine(xen_pts, img_pts)
        entry["affine_3x3"] = affine
        # Disconnect any linked affine — manual registration takes over
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
            entry["affine_disconnect"] = None
        entry["affine_source_name"] = None
        affine_combo.blockSignals(True)
        affine_combo.setCurrentIndex(0)
        affine_combo.blockSignals(False)
        _apply_entry_affine(entry)
        # Display residuals
        pixel_size = ctx.pixel_size
        lines = [f"Registered: {n} landmarks, similarity transform"]
        lines.append(f"Mean residual: {residuals.mean():.1f} px"
                     f" ({residuals.mean() * pixel_size:.1f} um)")
        lines.append(f"Max residual: {residuals.max():.1f} px"
                     f" ({residuals.max() * pixel_size:.1f} um)")
        scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
        lines.append(f"Scale: {scale:.4f}")
        reg_status_label.setText("\n".join(lines))
        # Persist landmarks + affine
        ename = entry["element_name"]
        save_landmarks_to_sdata(ctx, f"{ename}_xenium_lm", xen_pts)
        save_landmarks_to_sdata(ctx, f"{ename}_image_lm", img_pts)
        ctx.record_node(
            f"extimg:register:{ename}",
            f"\n# External image landmark registration ({ename})\n"
            f"from xenium_viewer.utils.registration import compute_landmark_affine\n"
            f"ext_xen_pts = np.array({xen_pts.tolist()})\n"
            f"ext_img_pts = np.array({img_pts.tolist()})\n"
            f"ext_affine, ext_residuals = compute_landmark_affine(ext_xen_pts, ext_img_pts)",
            deps=["preamble"], kind=TERMINAL, label=f"External image registration: {ename}",
        )

    # ── Entry creation / registration ───────────────────────────────────

    def _register_entry(entry):
        ctx.external_images_state.append(entry)
        n_ch = len(entry.get("channel_names", []))
        item = QListWidgetItem(
            f"{Path(entry.get('path') or entry.get('element_name', '?')).name}"
            f"  ({n_ch}ch)"
        )
        list_widget.addItem(item)
        list_widget.setCurrentRow(list_widget.count() - 1)

    def _on_loaded(result, display_name, path):
        pyramid, tif, channel_axis, channel_names = result
        # Build initial channel states
        colors = default_channel_colors(channel_names)
        clims = auto_contrast(pyramid, channel_axis=channel_axis or 0)
        # Compute data range from clims for slider bounds
        channel_states = []
        for i in range(len(channel_names)):
            lo, hi = clims[i]
            # Use a wider range for slider: allow going beyond auto-contrast
            data_min = 0.0
            data_max = max(hi * 1.5, 65535.0) if pyramid[0].dtype == np.uint16 else max(hi * 1.5, 255.0)
            channel_states.append({
                "visible": True,
                "color": colors[i],
                "clim": [lo, hi],
                "data_min": data_min,
                "data_max": data_max,
            })

        layer = _build_composite_layer(
            pyramid, channel_axis, channel_states, display_name,
        )
        # Determine image shape for flip affine
        base = pyramid[0]
        if channel_axis is not None:
            spatial = [s for i, s in enumerate(base.shape) if i != channel_axis]
        else:
            spatial = list(base.shape[:2])
        image_shape_yx = tuple(spatial[:2])

        element_name = f"ext_{_slugify(Path(path).stem)}"
        entry = {
            "element_name": element_name,
            "path": str(path),
            "tif": tif,
            "pyramid": pyramid,
            "channel_axis": channel_axis,
            "channel_names": list(channel_names),
            "channel_states": channel_states,
            "layer_ref": layer,
            "affine_source_name": None,
            "affine_disconnect": None,
            "opacity": 1.0,
            # Registration state
            "flip_v": False,
            "flip_h": False,
            "affine_3x3": None,
            "image_shape_yx": image_shape_yx,
            "xenium_lm_layer": None,
            "image_lm_layer": None,
        }
        _register_entry(entry)
        save_external_image_to_sdata(
            ctx, element_name, pyramid, channel_axis, channel_names,
        )
        ctx.record_node(
            f"extimg:load:{element_name}",
            f"\n# Load multichannel image ({element_name})\n"
            f"from xenium_viewer.utils.registration import load_multichannel_pyramid\n"
            f"ext_pyramid, ext_tif, channel_axis, channel_names = "
            f"load_multichannel_pyramid(\"{path}\")",
            deps=["preamble"], kind=TERMINAL, label=f"Load external image: {element_name}",
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
        # Disconnect affine
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
        # Remove composite layer
        lyr = entry.get("layer_ref")
        if lyr is not None:
            try:
                viewer.layers.remove(lyr)
            except Exception:
                pass
        # Remove landmark layers
        _remove_landmark_layers(entry)
        # Close file handle
        try:
            tif = entry.get("tif")
            if tif is not None:
                tif.close()
        except Exception:
            pass
        # Remove from sdata
        try:
            element = entry.get("element_name")
            if element and ctx.sdata is not None and element in ctx.sdata:
                safe_delete_element(ctx.sdata, element)
            # Also remove landmarks from sdata
            for suffix in ("_xenium_lm", "_image_lm"):
                lm_name = f"{element}{suffix}"
                if lm_name in ctx.sdata:
                    safe_delete_element(ctx.sdata, lm_name)
        except Exception as e:
            print(f"  Warning: could not delete from sdata: {e}")
        list_widget.takeItem(row)
        _update_panel()

    # ── Wire signals ─────────────────────────────────────────────────────
    add_btn.clicked.connect(on_add)
    list_widget.currentRowChanged.connect(lambda _i: _update_panel())
    opacity_slider.valueChanged.connect(_apply_opacity)
    show_all_btn.clicked.connect(lambda: _set_all_channels(True))
    hide_all_btn.clicked.connect(lambda: _set_all_channels(False))
    affine_combo.currentIndexChanged.connect(_apply_affine_source)
    flip_v_chk.toggled.connect(on_flip_changed)
    flip_h_chk.toggled.connect(on_flip_changed)
    add_xen_lm_btn.clicked.connect(on_add_xenium_lm)
    add_img_lm_btn.clicked.connect(on_add_image_lm)
    clear_lm_btn.clicked.connect(on_clear_lm)
    register_btn.clicked.connect(on_register)
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
            channel_axis = 0 if pyramid[0].ndim >= 3 else None
            ui = ui_by_name.get(element_name, {})
            display_name = Path(ui.get("path", element_name)).stem or element_name

            # Restore or rebuild channel states
            saved_ch = ui.get("channel_states")
            if saved_ch and len(saved_ch) == len(channel_names):
                channel_states = saved_ch
            else:
                colors = default_channel_colors(channel_names)
                clims = auto_contrast(pyramid, channel_axis=channel_axis or 0)
                channel_states = []
                for i in range(len(channel_names)):
                    lo, hi = clims[i]
                    data_max = 65535.0 if pyramid[0].dtype == np.uint16 else 255.0
                    channel_states.append({
                        "visible": True, "color": colors[i],
                        "clim": [lo, hi], "data_min": 0.0,
                        "data_max": max(hi * 1.5, data_max),
                    })

            layer = _build_composite_layer(
                pyramid, channel_axis, channel_states, display_name,
                opacity=float(ui.get("opacity", 1.0)),
            )

            # Image shape for flip affine
            base = pyramid[0]
            if channel_axis is not None:
                spatial = [s for i, s in enumerate(base.shape) if i != channel_axis]
            else:
                spatial = list(base.shape[:2])
            image_shape_yx = tuple(spatial[:2])

            entry = {
                "element_name": element_name,
                "path": ui.get("path") or element_name,
                "tif": None,
                "pyramid": pyramid,
                "channel_axis": channel_axis,
                "channel_names": channel_names,
                "channel_states": channel_states,
                "layer_ref": layer,
                "affine_source_name": ui.get("affine_source_name"),
                "affine_disconnect": None,
                "opacity": float(ui.get("opacity", 1.0)),
                "flip_v": bool(ui.get("flip_v", False)),
                "flip_h": bool(ui.get("flip_h", False)),
                "affine_3x3": None,
                "image_shape_yx": image_shape_yx,
                "xenium_lm_layer": None,
                "image_lm_layer": None,
            }

            # Apply affine from sdata (authoritative), fall back to session attrs
            saved_affine = meta.get("affine_matrix")
            if saved_affine is None:
                saved_affine = ui.get("affine_matrix")
            if saved_affine is not None:
                try:
                    layer.affine = np.array(saved_affine, dtype=np.float64)
                except Exception:
                    pass

            _register_entry(entry)

            # Restore landmarks if available
            xen_lm_data = load_landmarks_from_sdata(
                ctx.sdata, f"{element_name}_xenium_lm")
            img_lm_data = load_landmarks_from_sdata(
                ctx.sdata, f"{element_name}_image_lm")
            if xen_lm_data is not None or img_lm_data is not None:
                _create_landmark_layers(entry)
                if xen_lm_data is not None and entry["xenium_lm_layer"] is not None:
                    entry["xenium_lm_layer"].data = xen_lm_data
                if img_lm_data is not None and entry["image_lm_layer"] is not None:
                    entry["image_lm_layer"].data = img_lm_data

            # Re-link affine if named source is present
            if entry["affine_source_name"]:
                src = find_layer_by_name(viewer, entry["affine_source_name"])
                if src is not None:
                    try:
                        disconnect = link_affine(layer, src, viewer=viewer)
                        entry["affine_disconnect"] = disconnect
                    except Exception:
                        pass

        # Deferred affine linking for async-loaded source layers
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
                lyr = e.get("layer_ref")
                if src is not None and lyr is not None:
                    try:
                        e["affine_disconnect"] = link_affine(lyr, src, viewer=viewer)
                        save_overlay_affine_to_sdata(
                            ctx, e["element_name"], lyr.affine.affine_matrix,
                        )
                    except Exception:
                        pass
                else:
                    still_pending.append(e)
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
