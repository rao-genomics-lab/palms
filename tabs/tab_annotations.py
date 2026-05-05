"""Tab: Annotation Management.

Lets users draw named tissue annotations (bone, adipocyte, vessel, etc.) on the
napari Annotations shapes layer, assign type labels, pick per-type colours,
and import/export GeoJSON.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import numpy as np
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QTableWidget, QTableWidgetItem, QHeaderView,
    QFileDialog, QMessageBox,
)

from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext

# Default palette for auto-assigning colours to new annotation types
_DEFAULT_PALETTE = [
    (1.0, 1.0, 0.0),   # yellow
    (1.0, 0.5, 0.0),   # orange
    (0.0, 1.0, 1.0),   # cyan
    (1.0, 0.0, 1.0),   # magenta
    (0.0, 1.0, 0.0),   # green
    (0.5, 0.0, 1.0),   # purple
    (1.0, 0.0, 0.0),   # red
    (0.0, 0.5, 1.0),   # sky blue
]


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    layer = ctx.annotation_layer
    status = StatusProxy(ctx.viewer)

    # ── Type entry ────────────────────────────────────────────────────────────
    type_label = QLabel("Annotation type:")
    type_edit = QLineEdit()
    type_edit.setPlaceholderText("e.g. bone, adipocyte, vessel")
    assign_btn = QPushButton("Assign to selected shapes")

    type_row = QWidget()
    type_row_layout = QHBoxLayout()
    type_row_layout.setContentsMargins(0, 0, 0, 0)
    type_row_layout.addWidget(type_label)
    type_row_layout.addWidget(type_edit)
    type_row_layout.addWidget(assign_btn)
    type_row.setLayout(type_row_layout)

    # ── Type list table ───────────────────────────────────────────────────────
    table_label = QLabel("Annotation types (click colour to change):")
    type_table = QTableWidget(0, 3)
    type_table.setHorizontalHeaderLabels(["Type", "Count", "Colour"])
    type_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
    type_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
    type_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
    type_table.setMaximumHeight(200)
    type_table.setEditTriggers(QTableWidget.NoEditTriggers)
    type_table.setSelectionBehavior(QTableWidget.SelectRows)

    # ── Action buttons ────────────────────────────────────────────────────────
    btn_row1 = QWidget()
    btn_row1_layout = QHBoxLayout()
    btn_row1_layout.setContentsMargins(0, 0, 0, 0)
    delete_btn = QPushButton("Delete selected shapes")
    clear_btn = QPushButton("Clear all annotations")
    btn_row1_layout.addWidget(delete_btn)
    btn_row1_layout.addWidget(clear_btn)
    btn_row1.setLayout(btn_row1_layout)

    btn_row2 = QWidget()
    btn_row2_layout = QHBoxLayout()
    btn_row2_layout.setContentsMargins(0, 0, 0, 0)
    import_btn = QPushButton("Import GeoJSON...")
    export_btn = QPushButton("Export GeoJSON...")
    btn_row2_layout.addWidget(import_btn)
    btn_row2_layout.addWidget(export_btn)
    btn_row2.setLayout(btn_row2_layout)

    status_label = QLabel("")
    status_label.setWordWrap(True)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _get_types() -> list[str]:
        return list(layer.properties.get("annotation_type", []))

    def _get_colors() -> list:
        """Return per-shape edge colors as list of RGBA tuples."""
        ec = layer.edge_color  # (N, 4) ndarray or list
        if hasattr(ec, '__len__') and len(ec) > 0:
            return [tuple(float(v) for v in row) for row in ec]
        return []

    def _next_color(type_name: str) -> tuple:
        types_seen = list(state.get("annot_types", {}).keys())
        if type_name in types_seen:
            return state["annot_types"][type_name]
        idx = len(types_seen) % len(_DEFAULT_PALETTE)
        return _DEFAULT_PALETTE[idx]

    def _is_valid_type(t) -> bool:
        """Return True only for non-empty, non-NaN type strings."""
        if t is None:
            return False
        try:
            import math
            if isinstance(t, float) and math.isnan(t):
                return False
        except (TypeError, ValueError):
            pass
        return bool(str(t).strip())

    def _refresh_table():
        types = [str(t) if _is_valid_type(t) else "" for t in _get_types()]
        colors = _get_colors()
        # Count shapes per type (skip unassigned)
        from collections import Counter
        counts = Counter(t for t in types if t)

        # Build unique type→colour mapping from current shapes
        type_color: dict[str, tuple] = dict(state.get("annot_types", {}))
        for i, t in enumerate(types):
            if t and t not in type_color:
                if i < len(colors):
                    type_color[t] = colors[i]

        type_table.setRowCount(0)
        for typename, count in sorted(counts.items()):
            row = type_table.rowCount()
            type_table.insertRow(row)
            type_table.setItem(row, 0, QTableWidgetItem(typename))
            type_table.setItem(row, 1, QTableWidgetItem(str(count)))

            rgba = type_color.get(typename, (1.0, 1.0, 0.0, 1.0))
            r, g, b = int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
            color_item = QTableWidgetItem()
            color_item.setBackground(
                __import__('qtpy.QtGui', fromlist=['QColor']).QColor(r, g, b)
            )
            type_table.setItem(row, 2, color_item)

    def _assign_type():
        typename = type_edit.text().strip()
        if not typename:
            status_label.setText("Enter an annotation type name first.")
            return
        if layer is None:
            return

        selected = sorted(layer.selected_data)
        if not selected:
            status_label.setText("Select shapes on the Annotations layer first.")
            return

        # Ensure annot_types has a colour for this type
        if "annot_types" not in state:
            state["annot_types"] = {}
        if typename not in state["annot_types"]:
            state["annot_types"][typename] = _next_color(typename) + (1.0,)

        # Update properties
        current = list(layer.properties.get("annotation_type", [""] * len(layer.data)))
        while len(current) < len(layer.data):
            current.append("")
        for idx in selected:
            current[idx] = typename
        layer.properties = {"annotation_type": current}
        layer.refresh_colors(update_color_mapping=False)

        # Update edge colour for affected shapes
        rgba = state["annot_types"][typename]
        colors = list(layer.edge_color) if len(layer.edge_color) else []
        while len(colors) < len(layer.data):
            colors.append(np.array([1., 1., 0., 1.]))
        for idx in selected:
            colors[idx] = np.array([rgba[0], rgba[1], rgba[2], 1.0])
        layer.edge_color = colors

        _persist()
        _refresh_table()
        status_label.setText(f"Assigned type '{typename}' to {len(selected)} shape(s).")

    def _on_table_cell_clicked(row, col):
        if col != 2:
            return
        typename_item = type_table.item(row, 0)
        if typename_item is None:
            return
        typename = typename_item.text()

        from qtpy.QtWidgets import QColorDialog
        from qtpy.QtGui import QColor
        current = state.get("annot_types", {}).get(typename, (1.0, 1.0, 0.0, 1.0))
        qc = QColor(int(current[0]*255), int(current[1]*255), int(current[2]*255))
        chosen = QColorDialog.getColor(qc, type_table, f"Colour for '{typename}'")
        if not chosen.isValid():
            return

        rgba = (chosen.redF(), chosen.greenF(), chosen.blueF(), 1.0)
        if "annot_types" not in state:
            state["annot_types"] = {}
        state["annot_types"][typename] = rgba

        # Re-apply colour to all shapes of this type
        types = _get_types()
        colors = list(layer.edge_color) if len(layer.edge_color) else []
        while len(colors) < len(layer.data):
            colors.append(np.array([1., 1., 0., 1.]))
        for i, t in enumerate(types):
            if t == typename:
                colors[i] = np.array([rgba[0], rgba[1], rgba[2], 1.0])
        layer.edge_color = colors

        _persist()
        _refresh_table()

    def _delete_selected():
        selected = sorted(layer.selected_data, reverse=True)
        if not selected:
            status_label.setText("No shapes selected.")
            return
        for idx in selected:
            layer.remove(idx)
        _persist()
        _refresh_table()
        status_label.setText(f"Deleted {len(selected)} shape(s).")

    def _clear_all():
        if not layer.data:
            return
        reply = QMessageBox.question(
            None, "Clear annotations",
            "Remove all annotations? This cannot be undone.",
            QMessageBox.Yes | QMessageBox.Cancel,
        )
        if reply != QMessageBox.Yes:
            return
        layer.data = []
        layer.properties = {"annotation_type": []}
        state.pop("annot_types", None)
        _persist()
        _refresh_table()
        status_label.setText("All annotations cleared.")

    def _import_geojson():
        path, _ = QFileDialog.getOpenFileName(
            None, "Import Annotations GeoJSON", "", "GeoJSON Files (*.geojson *.json)"
        )
        if not path:
            return
        try:
            _do_import_geojson(path)
            _persist()
            _refresh_table()
            status_label.setText(f"Imported annotations from {path}")
        except Exception as e:
            QMessageBox.critical(None, "Import failed", str(e))

    def _do_import_geojson(path: str):
        from shapely.geometry import shape as shapely_shape
        with open(path) as f:
            fc = json.load(f)
        features = fc if isinstance(fc, list) else fc.get("features", [])
        if not features and "geometry" in fc:
            features = [fc]

        existing_data = list(layer.data)
        existing_types = list(layer.properties.get("annotation_type", []))

        if "annot_types" not in state:
            state["annot_types"] = {}

        for feat in features:
            geom = feat.get("geometry") if isinstance(feat, dict) else feat
            props = feat.get("properties") or {} if isinstance(feat, dict) else {}
            typename = str(props.get("annotation_type", props.get("name", "imported")))

            shp = shapely_shape(geom)
            # Handle Polygon and MultiPolygon
            polys = []
            if shp.geom_type == "Polygon":
                polys = [shp]
            elif shp.geom_type == "MultiPolygon":
                polys = list(shp.geoms)
            else:
                continue

            if typename not in state["annot_types"]:
                state["annot_types"][typename] = _next_color(typename) + (1.0,)

            for poly in polys:
                coords_xy = np.array(poly.exterior.coords[:-1], dtype=np.float64)
                coords_yx = coords_xy[:, ::-1]  # xy → yx for napari
                existing_data.append(coords_yx)
                existing_types.append(typename)

        layer.data = existing_data
        layer.properties = {"annotation_type": existing_types}

        # Apply colours
        colors = list(layer.edge_color) if len(layer.edge_color) else []
        while len(colors) < len(existing_data):
            colors.append(np.array([1., 1., 0., 1.]))
        for i, t in enumerate(existing_types):
            if t in state["annot_types"]:
                rgba = state["annot_types"][t]
                colors[i] = np.array([rgba[0], rgba[1], rgba[2], 1.0])
        layer.edge_color = colors

    def _export_geojson():
        if not layer.data:
            status_label.setText("No annotations to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export Annotations GeoJSON", "annotations.geojson",
            "GeoJSON Files (*.geojson *.json)"
        )
        if not path:
            return
        try:
            types = _get_types()
            while len(types) < len(layer.data):
                types.append("")
            features = []
            for arr, t in zip(layer.data, types):
                xy = arr[:, ::-1].tolist()  # yx → xy
                xy.append(xy[0])  # close ring
                feat = {
                    "type": "Feature",
                    "geometry": {"type": "Polygon", "coordinates": [xy]},
                    "properties": {"annotation_type": t},
                }
                features.append(feat)
            fc = {"type": "FeatureCollection", "features": features}
            with open(path, "w") as f:
                json.dump(fc, f)
            status_label.setText(f"Exported {len(features)} annotation(s) to {path}")
        except Exception as e:
            QMessageBox.critical(None, "Export failed", str(e))

    def _persist():
        from utils.adata_persistence import save_annotations_to_sdata
        save_annotations_to_sdata(ctx)

    # ── Connect signals ────────────────────────────────────────────────────────
    assign_btn.clicked.connect(_assign_type)
    type_table.cellClicked.connect(_on_table_cell_clicked)
    delete_btn.clicked.connect(_delete_selected)
    clear_btn.clicked.connect(_clear_all)
    import_btn.clicked.connect(_import_geojson)
    export_btn.clicked.connect(_export_geojson)

    # Auto-refresh table when shapes are added/removed
    if layer is not None:
        layer.events.data.connect(lambda _: _refresh_table())

    # ── Session restore ────────────────────────────────────────────────────────
    def _restore_session(session):
        from utils.adata_persistence import load_annotations_from_sdata
        shapes, types = load_annotations_from_sdata(ctx.sdata)
        if not shapes:
            return
        layer.data = shapes
        layer.properties = {"annotation_type": types}

        # Rebuild annot_types colour map
        if "annot_types" not in state:
            state["annot_types"] = {}
        for t in set(types):
            if t and t not in state["annot_types"]:
                state["annot_types"][t] = _next_color(t) + (1.0,)

        # Apply edge colours
        colors = []
        for t in types:
            rgba = state["annot_types"].get(t, (1.0, 1.0, 0.0, 1.0))
            colors.append(np.array([rgba[0], rgba[1], rgba[2], 1.0]))
        if colors:
            layer.edge_color = colors

        _refresh_table()

    # ── Build tab ─────────────────────────────────────────────────────────────
    container = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)
    layout.addWidget(type_row)
    layout.addWidget(table_label)
    layout.addWidget(type_table)
    layout.addWidget(btn_row1)
    layout.addWidget(btn_row2)
    layout.addWidget(status_label)
    layout.addStretch()
    container.setLayout(layout)

    from qtpy.QtWidgets import QScrollArea
    scroll = QScrollArea()
    scroll.setWidget(container)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.NoFrame)

    return scroll, {"restore_session": _restore_session}
