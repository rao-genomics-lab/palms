"""
Patch Overlays tab — loads phikon patch-cluster outputs and subclone
prediction CSVs. Patches render as a napari Shapes layer of coloured
rectangles (one polygon per patch) so fills and outlines follow zoom.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import numpy as np

from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.zarr_safe import safe_delete_element

from qtpy.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QGridLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem,
    QMessageBox, QPushButton, QScrollArea, QSlider, QVBoxLayout, QWidget,
    QSpinBox,
)
from qtpy.QtCore import Qt

from xenium_viewer.tabs._helpers import make_tab
from xenium_viewer.utils.coloring import PATCH_PALETTES
from xenium_viewer.utils.affine_linking import (
    link_affine, list_transformable_layers, find_layer_by_name,
)
from xenium_viewer.utils.adata_persistence import (
    save_patch_overlay_to_sdata, load_patch_overlays_from_sdata,
    save_overlay_affine_to_sdata, _slugify,
)
from xenium_viewer.utils.patch_overlay_io import (
    PatchOverlayData, estimate_stride, infer_patch_size_from_path,
    load_phikon_folder, load_subclone_csv,
)

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def _confirm_patch_size(parent, inferred: Optional[int], stride: Optional[int]) -> Optional[int]:
    """Modal dialog asking the user to confirm / override patch size. Returns px or None."""
    dlg = QDialog(parent)
    dlg.setWindowTitle("Confirm patch size")
    layout = QVBoxLayout()
    info = [
        f"Inferred patch size: {inferred if inferred else '(unknown)'} px",
        f"Estimated grid stride: {stride if stride else '(unknown)'} px",
        "",
        "Override patch size in pixels if needed:",
    ]
    for line in info:
        layout.addWidget(QLabel(line))
    spin = QSpinBox()
    spin.setRange(1, 8192)
    spin.setValue(inferred or stride or 128)
    layout.addWidget(spin)
    btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    btns.accepted.connect(dlg.accept)
    btns.rejected.connect(dlg.reject)
    layout.addWidget(btns)
    dlg.setLayout(layout)
    if dlg.exec_() == QDialog.Accepted:
        return int(spin.value())
    return None


def _build_rectangles_yx(coords_xy: np.ndarray, patch_size: int) -> list:
    """Build a list of Nx(4,2) rectangle polygon arrays in napari (y, x) order.

    Each rectangle has corners TL → TR → BR → BL.
    """
    xs = coords_xy[:, 0].astype(np.float64)
    ys = coords_xy[:, 1].astype(np.float64)
    s = float(patch_size)
    # (N, 4, 2) stacked
    corners = np.empty((len(coords_xy), 4, 2), dtype=np.float64)
    corners[:, 0, 0] = ys
    corners[:, 0, 1] = xs
    corners[:, 1, 0] = ys
    corners[:, 1, 1] = xs + s
    corners[:, 2, 0] = ys + s
    corners[:, 2, 1] = xs + s
    corners[:, 3, 0] = ys + s
    corners[:, 3, 1] = xs
    return [corners[i] for i in range(len(coords_xy))]


def _base_face_colors(cluster_ids: np.ndarray, palette: np.ndarray) -> np.ndarray:
    """Map cluster ids (int64) → RGBA rows via palette, 0-normalised so the
    smallest cluster ID always maps to palette[0]."""
    ids = np.asarray(cluster_ids, dtype=np.int64)
    out = np.zeros((ids.size, 4), dtype=np.float32)
    valid = ids >= 0
    if valid.any():
        normed = ids[valid] - ids[valid].min()
        out[valid] = palette[normed % len(palette)]
    return out


def _compute_display_colors(entry: dict) -> tuple[np.ndarray, np.ndarray]:
    """Return (face, edge) RGBA arrays given entry state (hidden ids, confidence,
    outline-only, palette, active cluster column).
    """
    palette = PATCH_PALETTES.get(entry["palette_name"], PATCH_PALETTES["tab20"])
    vals = entry["cluster_columns"][entry["active_cluster_column"]]
    base = _base_face_colors(vals, palette)
    face = base.copy()
    edge = base.copy()
    # visibility mask (True = visible)
    visible = np.ones(len(vals), dtype=bool)
    hidden = entry.get("hidden_cluster_ids") or set()
    if hidden:
        for cid in hidden:
            visible &= vals != cid
    thr = entry.get("confidence_threshold", 0.0)
    conf = entry.get("confidence")
    if conf is not None and thr > 0:
        visible &= conf >= thr

    if entry.get("outline_only", False):
        face[:, 3] = 0.0
        edge[:, 3] = 1.0
    else:
        edge[:, 3] = 0.0

    face[~visible, 3] = 0.0
    edge[~visible, 3] = 0.0
    return face, edge


def build_tab(ctx: "ViewerContext"):
    viewer = ctx.viewer

    # ── UI ───────────────────────────────────────────────────────────────
    root = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)

    add_phikon_btn = QPushButton("Add phikon clustering…")
    add_subclone_btn = QPushButton("Add subclone predictions…")
    add_row = QHBoxLayout()
    add_row.addWidget(add_phikon_btn)
    add_row.addWidget(add_subclone_btn)
    layout.addLayout(add_row)

    list_widget = QListWidget()
    list_widget.setMinimumHeight(120)
    layout.addWidget(list_widget)

    panel = QGroupBox("Selected overlay")
    panel_layout = QVBoxLayout()

    status_label = QLabel("—")
    status_label.setWordWrap(True)
    panel_layout.addWidget(status_label)

    col_row = QHBoxLayout()
    col_row.addWidget(QLabel("Cluster column:"))
    cluster_col_combo = QComboBox()
    col_row.addWidget(cluster_col_combo, 1)
    panel_layout.addLayout(col_row)

    pal_row = QHBoxLayout()
    pal_row.addWidget(QLabel("Palette:"))
    palette_combo = QComboBox()
    palette_combo.addItems(["tab10", "tab20", "glasbey_dark", "Set1", "Set3", "ARMS (Set1+Set2+Dark2)"])
    pal_row.addWidget(palette_combo, 1)
    panel_layout.addLayout(pal_row)

    aff_row = QHBoxLayout()
    aff_row.addWidget(QLabel("Apply transform from:"))
    affine_combo = QComboBox()
    aff_row.addWidget(affine_combo, 1)
    panel_layout.addLayout(aff_row)

    outline_chk = QCheckBox("Outline only")
    panel_layout.addWidget(outline_chk)

    edge_row = QHBoxLayout()
    edge_row.addWidget(QLabel("Edge width:"))
    edge_slider = QSlider(Qt.Horizontal)
    edge_slider.setRange(0, 20)
    edge_slider.setValue(2)
    edge_row.addWidget(edge_slider)
    edge_value_label = QLabel("2")
    edge_row.addWidget(edge_value_label)
    panel_layout.addLayout(edge_row)

    op_row = QHBoxLayout()
    op_row.addWidget(QLabel("Opacity:"))
    opacity_slider = QSlider(Qt.Horizontal)
    opacity_slider.setRange(0, 100)
    opacity_slider.setValue(80)
    op_row.addWidget(opacity_slider)
    panel_layout.addLayout(op_row)

    conf_row = QHBoxLayout()
    conf_label = QLabel("Confidence ≥ 0.00:")
    conf_row.addWidget(conf_label)
    conf_slider = QSlider(Qt.Horizontal)
    conf_slider.setRange(0, 100)
    conf_slider.setValue(0)
    conf_row.addWidget(conf_slider)
    panel_layout.addLayout(conf_row)

    panel_layout.addWidget(QLabel("Visible clusters:"))
    filter_scroll = QScrollArea()
    filter_scroll.setWidgetResizable(True)
    filter_scroll.setMinimumHeight(100)
    filter_scroll.setMaximumHeight(200)
    filter_container = QWidget()
    filter_grid = QGridLayout()
    filter_container.setLayout(filter_grid)
    filter_scroll.setWidget(filter_container)
    panel_layout.addWidget(filter_scroll)

    sel_row = QHBoxLayout()
    select_all_btn = QPushButton("Select all")
    deselect_all_btn = QPushButton("Deselect all")
    sel_row.addWidget(select_all_btn)
    sel_row.addWidget(deselect_all_btn)
    panel_layout.addLayout(sel_row)

    remove_btn = QPushButton("Remove overlay")
    panel_layout.addWidget(remove_btn)

    panel.setLayout(panel_layout)
    layout.addWidget(panel)

    root.setLayout(layout)
    tab_widget = make_tab(root)

    # ── Helpers ──────────────────────────────────────────────────────────

    def _current_entry() -> Optional[dict]:
        row = list_widget.currentRow()
        if row < 0 or row >= len(ctx.patch_overlays_state):
            return None
        return ctx.patch_overlays_state[row]

    def _apply_colors(entry):
        """Recompute face + edge colour arrays and push them to the shapes layer."""
        face, edge = _compute_display_colors(entry)
        lyr = entry["shapes_layer"]
        try:
            lyr.face_color = face
            lyr.edge_color = edge
            lyr.edge_width = int(entry.get("edge_width", 2))
        except Exception as e:
            print(f"  Warning: could not update patch colors: {e}")

    def _on_cluster_cb_toggled(entry, cid, checked):
        hidden = entry.setdefault("hidden_cluster_ids", set())
        if checked:
            hidden.discard(int(cid))
        else:
            hidden.add(int(cid))
        _apply_colors(entry)

    def _rebuild_cluster_filter_grid(entry):
        while filter_grid.count():
            item = filter_grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        entry["cluster_checkboxes"] = {}

        vals = entry["cluster_columns"][entry["active_cluster_column"]]
        uniq = sorted(set(int(v) for v in vals if v >= 0))
        cols = 3
        hidden = set(entry.get("hidden_cluster_ids", []))
        for i, cid in enumerate(uniq):
            cb = QCheckBox(str(cid))
            cb.setChecked(cid not in hidden)
            cb.toggled.connect(
                lambda chk, _e=entry, _c=cid: _on_cluster_cb_toggled(_e, _c, chk)
            )
            filter_grid.addWidget(cb, i // cols, i % cols)
            entry["cluster_checkboxes"][cid] = cb

    def _refresh_affine_choices():
        entry = _current_entry()
        excluded = [entry["shapes_layer"]] if entry and entry.get("shapes_layer") else ()
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

    def _update_panel():
        entry = _current_entry()
        enabled = entry is not None
        for w in (cluster_col_combo, palette_combo, affine_combo, outline_chk,
                  edge_slider, opacity_slider, conf_slider, remove_btn,
                  select_all_btn, deselect_all_btn):
            w.setEnabled(enabled)
        if not enabled:
            status_label.setText("—")
            return

        n = len(entry["coords_xy"])
        src = entry.get("source_path", "?")
        kind = entry.get("source_kind", "?")
        status_label.setText(
            f"{Path(src).name} ({kind}) — {n} patches, size {entry['patch_size_px']}px"
        )

        cluster_col_combo.blockSignals(True)
        cluster_col_combo.clear()
        for col in entry["cluster_columns"].keys():
            cluster_col_combo.addItem(col, col)
        idx = cluster_col_combo.findData(entry["active_cluster_column"])
        if idx >= 0:
            cluster_col_combo.setCurrentIndex(idx)
        cluster_col_combo.blockSignals(False)
        cluster_col_combo.setEnabled(len(entry["cluster_columns"]) > 1)

        palette_combo.blockSignals(True)
        idx = palette_combo.findText(entry.get("palette_name", "tab20"))
        if idx >= 0:
            palette_combo.setCurrentIndex(idx)
        palette_combo.blockSignals(False)

        outline_chk.blockSignals(True)
        outline_chk.setChecked(entry.get("outline_only", False))
        outline_chk.blockSignals(False)
        edge_slider.blockSignals(True)
        edge_slider.setValue(int(entry.get("edge_width", 2)))
        edge_slider.blockSignals(False)
        edge_value_label.setText(str(int(entry.get("edge_width", 2))))

        opacity_slider.blockSignals(True)
        opacity_slider.setValue(int(entry.get("opacity", 0.8) * 100))
        opacity_slider.blockSignals(False)

        has_conf = entry.get("confidence") is not None
        conf_slider.setEnabled(has_conf)
        conf_label.setEnabled(has_conf)
        thr = entry.get("confidence_threshold", 0.0)
        conf_slider.blockSignals(True)
        conf_slider.setValue(int(thr * 100))
        conf_slider.blockSignals(False)
        conf_label.setText(f"Confidence ≥ {thr:.2f}:")

        _rebuild_cluster_filter_grid(entry)
        _refresh_affine_choices()

    def _create_shapes_layer(coords_xy, patch_size, face, edge, edge_width,
                             layer_name, opacity):
        rects = _build_rectangles_yx(coords_xy, patch_size)
        return viewer.add_shapes(
            rects,
            shape_type="polygon",
            face_color=face,
            edge_color=edge,
            edge_width=int(edge_width),
            name=layer_name,
            opacity=float(opacity),
        )

    def _add_overlay_from_data(data: PatchOverlayData):
        patch_size = data.patch_size or 0
        stride = estimate_stride(data.coords_xy)
        if not patch_size:
            patch_size = stride or 128
        confirmed = _confirm_patch_size(root, patch_size, stride)
        if confirmed is None:
            return
        patch_size = confirmed

        entry = {
            "element_name": f"patch_{_slugify(Path(data.source_path).stem)}",
            "source_path": data.source_path,
            "source_kind": data.source_kind,
            "coords_xy": data.coords_xy,
            "cluster_columns": data.cluster_columns,
            "active_cluster_column": data.active_cluster_column,
            "confidence": data.confidence,
            "confidence_threshold": 0.0,
            "patch_size_px": int(patch_size),
            "palette_name": "ARMS (Set1+Set2+Dark2)" if data.source_kind == "subclone" else "tab20",
            "affine_source_name": None,
            "affine_disconnect": None,
            "outline_only": False,
            "edge_width": 2,
            "opacity": 0.8,
            "hidden_cluster_ids": set(),
            "cluster_checkboxes": {},
            "shapes_layer": None,
        }
        face, edge = _compute_display_colors(entry)
        layer_name = f"patches :: {Path(data.source_path).stem}"
        ctx.set_status(f"Creating {len(data.coords_xy)} patch shapes…")
        try:
            shapes_layer = _create_shapes_layer(
                data.coords_xy, patch_size, face, edge,
                entry["edge_width"], layer_name, entry["opacity"],
            )
        except Exception as e:
            QMessageBox.critical(None, "Layer creation failed", str(e))
            return
        entry["shapes_layer"] = shapes_layer

        ctx.patch_overlays_state.append(entry)
        list_widget.addItem(QListWidgetItem(
            f"{data.source_kind}: {Path(data.source_path).name} "
            f"({len(data.coords_xy)}×{patch_size}px)"
        ))
        list_widget.setCurrentRow(list_widget.count() - 1)

        save_patch_overlay_to_sdata(
            ctx, entry["element_name"], data.coords_xy, patch_size,
            data.cluster_columns, data.confidence,
        )
        ctx.record_node(
            "viewer:patch_overlay",
            f"\n# Load patch overlay ({data.source_kind}) — viewer overlay\n"
            f"# source: {data.source_path}\n"
            f"# patch_size: {patch_size} px, N={len(data.coords_xy)}",
            deps=["preamble"],
            kind=TERMINAL,
            label="Patch overlay",
        )
        ctx.set_status(
            f"Loaded {len(data.coords_xy)} patches from {Path(data.source_path).name}"
        )

    def on_add_phikon():
        default_dir = str(ctx.data_path) if ctx.data_path else ""
        folder = QFileDialog.getExistingDirectory(
            None, "Choose phikon results folder", default_dir,
        )
        if not folder:
            return
        try:
            data = load_phikon_folder(Path(folder))
        except Exception as e:
            QMessageBox.critical(None, "Load failed", str(e))
            return
        _add_overlay_from_data(data)

    def on_add_subclone():
        default_dir = str(ctx.data_path) if ctx.data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Choose subclone predictions CSV", default_dir,
            "CSV Files (*.csv);;All Files (*)",
        )
        if not path:
            return
        try:
            data = load_subclone_csv(Path(path))
        except Exception as e:
            QMessageBox.critical(None, "Load failed", str(e))
            return
        _add_overlay_from_data(data)

    def on_cluster_col_changed(_ix):
        entry = _current_entry()
        if entry is None:
            return
        col = cluster_col_combo.currentData()
        if not col or col == entry.get("active_cluster_column"):
            return
        entry["active_cluster_column"] = col
        entry["hidden_cluster_ids"] = set()
        _rebuild_cluster_filter_grid(entry)
        _apply_colors(entry)

    def on_palette_changed(_ix):
        entry = _current_entry()
        if entry is None:
            return
        entry["palette_name"] = palette_combo.currentText()
        _apply_colors(entry)

    def on_outline_changed(chk):
        entry = _current_entry()
        if entry is None:
            return
        entry["outline_only"] = bool(chk)
        _apply_colors(entry)

    def on_edge_width(value):
        entry = _current_entry()
        edge_value_label.setText(str(int(value)))
        if entry is None:
            return
        entry["edge_width"] = int(value)
        try:
            entry["shapes_layer"].edge_width = int(value)
        except Exception:
            pass

    def on_opacity(value):
        entry = _current_entry()
        if entry is None:
            return
        entry["opacity"] = value / 100.0
        try:
            entry["shapes_layer"].opacity = value / 100.0
        except Exception:
            pass

    def on_conf_changed(value):
        entry = _current_entry()
        if entry is None:
            return
        thr = value / 100.0
        entry["confidence_threshold"] = thr
        conf_label.setText(f"Confidence ≥ {thr:.2f}:")
        _apply_colors(entry)

    def on_affine_changed(_ix):
        entry = _current_entry()
        if entry is None:
            return
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
            entry["affine_disconnect"] = None
        name = affine_combo.currentData()
        entry["affine_source_name"] = name
        if not name:
            return
        source = find_layer_by_name(viewer, name)
        if source is None:
            return
        try:
            entry["affine_disconnect"] = link_affine(
                entry["shapes_layer"], source, viewer=viewer,
            )
            # Persist to sdata
            save_overlay_affine_to_sdata(
                ctx, entry["element_name"],
                entry["shapes_layer"].affine.affine_matrix,
            )
        except Exception as e:
            print(f"  Warning: could not link patch affine: {e}")

    def on_select_all():
        entry = _current_entry()
        if entry is None:
            return
        entry["hidden_cluster_ids"] = set()
        for cb in entry["cluster_checkboxes"].values():
            cb.blockSignals(True)
            cb.setChecked(True)
            cb.blockSignals(False)
        _apply_colors(entry)

    def on_deselect_all():
        entry = _current_entry()
        if entry is None:
            return
        entry["hidden_cluster_ids"] = set(entry["cluster_checkboxes"].keys())
        for cb in entry["cluster_checkboxes"].values():
            cb.blockSignals(True)
            cb.setChecked(False)
            cb.blockSignals(False)
        _apply_colors(entry)

    def on_remove():
        row = list_widget.currentRow()
        if row < 0 or row >= len(ctx.patch_overlays_state):
            return
        entry = ctx.patch_overlays_state.pop(row)
        cb = entry.get("affine_disconnect")
        if cb is not None:
            try:
                cb()
            except Exception:
                pass
        try:
            viewer.layers.remove(entry["shapes_layer"])
        except Exception:
            pass
        try:
            element = entry.get("element_name")
            if element and ctx.sdata is not None and element in ctx.sdata:
                safe_delete_element(ctx.sdata, element)
        except Exception as e:
            from xenium_viewer.utils.adata_persistence import _maybe_show_permission_dialog
            _maybe_show_permission_dialog(e, f"delete '{element}' from zarr cache")
            print(f"  Warning: could not delete {element} from sdata: {e}")
        list_widget.takeItem(row)
        _update_panel()

    # ── Wire signals ─────────────────────────────────────────────────────
    add_phikon_btn.clicked.connect(on_add_phikon)
    add_subclone_btn.clicked.connect(on_add_subclone)
    list_widget.currentRowChanged.connect(lambda _i: _update_panel())
    cluster_col_combo.currentIndexChanged.connect(on_cluster_col_changed)
    palette_combo.currentIndexChanged.connect(on_palette_changed)
    outline_chk.toggled.connect(on_outline_changed)
    edge_slider.valueChanged.connect(on_edge_width)
    opacity_slider.valueChanged.connect(on_opacity)
    conf_slider.valueChanged.connect(on_conf_changed)
    affine_combo.currentIndexChanged.connect(on_affine_changed)
    select_all_btn.clicked.connect(on_select_all)
    deselect_all_btn.clicked.connect(on_deselect_all)
    remove_btn.clicked.connect(on_remove)

    try:
        viewer.layers.events.inserted.connect(lambda _e: _refresh_affine_choices())
        viewer.layers.events.removed.connect(lambda _e: _refresh_affine_choices())
    except Exception:
        pass

    _update_panel()

    # ── Session restore ──────────────────────────────────────────────────
    def restore_session(session):
        entries = load_patch_overlays_from_sdata(ctx.sdata)
        ui_list = (session or {}).get("patch_overlays_ui") or []
        ui_by_name = {u.get("element_name"): u for u in ui_list if isinstance(u, dict)}
        for meta in entries:
            element_name = meta["element_name"]
            ui = ui_by_name.get(element_name, {})
            active = ui.get("active_cluster_column") or next(
                iter(meta["cluster_columns"]), None,
            )
            if active is None:
                continue
            source_kind = ui.get("source_kind") or (
                "phikon" if "phikon_cluster" in meta["cluster_columns"]
                else "subclone"
            )
            patch_size = int(ui.get("patch_size_px") or meta.get("patch_size") or 0)
            if patch_size <= 0:
                patch_size = 128

            entry = {
                "element_name": element_name,
                "source_path": ui.get("source_path") or element_name,
                "source_kind": source_kind,
                "coords_xy": meta["coords_xy"],
                "cluster_columns": meta["cluster_columns"],
                "active_cluster_column": active,
                "confidence": meta.get("confidence"),
                "confidence_threshold": float(ui.get("confidence_threshold", 0.0)),
                "patch_size_px": patch_size,
                "palette_name": ui.get("palette_name", "tab20"),
                "affine_source_name": ui.get("affine_source_name"),
                "affine_disconnect": None,
                "outline_only": bool(ui.get("outline_only", False)),
                "edge_width": int(ui.get("edge_width", 2)),
                "opacity": float(ui.get("opacity", 0.8)),
                "hidden_cluster_ids": set(ui.get("hidden_cluster_ids", []) or []),
                "cluster_checkboxes": {},
                "shapes_layer": None,
            }
            face, edge = _compute_display_colors(entry)
            layer_name = f"patches :: {Path(entry['source_path']).stem}"
            try:
                shapes_layer = _create_shapes_layer(
                    entry["coords_xy"], patch_size, face, edge,
                    entry["edge_width"], layer_name, entry["opacity"],
                )
            except Exception as e:
                print(f"  Warning: could not restore patch overlay {element_name}: {e}")
                continue
            entry["shapes_layer"] = shapes_layer

            # Apply affine from sdata (authoritative), fall back to session attrs
            saved_affine = meta.get("affine_matrix")  # from sdata transformations
            if saved_affine is None:
                saved_affine = ui.get("affine_matrix")  # fallback: session attrs
            if saved_affine is not None:
                try:
                    shapes_layer.affine = np.array(saved_affine, dtype=np.float64)
                except Exception:
                    pass

            ctx.patch_overlays_state.append(entry)
            list_widget.addItem(QListWidgetItem(
                f"{source_kind}: {Path(entry['source_path']).name} "
                f"({len(entry['coords_xy'])}×{patch_size}px)"
            ))
            if entry["affine_source_name"]:
                src = find_layer_by_name(viewer, entry["affine_source_name"])
                if src is not None:
                    try:
                        entry["affine_disconnect"] = link_affine(
                            shapes_layer, src, viewer=viewer,
                        )
                    except Exception:
                        pass

        # Deferred affine linking: source layers (e.g. H&E) may not exist
        # yet because they load asynchronously. Watch for layer insertions
        # and link once the named source appears.
        _pending = [
            e for e in ctx.patch_overlays_state
            if e.get("affine_source_name") and e.get("affine_disconnect") is None
        ]

        def _on_layer_inserted_for_pending(event=None):
            still_pending = []
            for e in _pending:
                if e.get("affine_disconnect") is not None:
                    continue  # already linked
                src = find_layer_by_name(viewer, e["affine_source_name"])
                lyr = e.get("shapes_layer")
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

        if ctx.patch_overlays_state:
            list_widget.setCurrentRow(0)

    return tab_widget, {"restore_session": restore_session}
