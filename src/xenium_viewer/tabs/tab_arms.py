"""Tab 10: ARMS Overlay — H&E, landmarks, GeoJSON tiles, tile DEG."""

from __future__ import annotations
import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QFileDialog,
    QCheckBox, QGridLayout, QScrollArea,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

from xenium_viewer.utils.registration import (
    load_he_pyramid, compute_landmark_affine, save_landmarks, load_landmarks,
)
from xenium_viewer.utils.gene_analysis import compute_arms_tile_deg

import matplotlib.pyplot as _plt
_set1 = _plt.get_cmap("Set1")
_set2 = _plt.get_cmap("Set2")
ARMS_CLUSTER_PALETTE = np.array(
    [list(_set1(i)[:3]) + [0.6] for i in range(9)] +
    [list(_set2(i)[:3]) + [0.6] for i in range(8)],
    dtype=np.float32,
)


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    arms_state = ctx.arms_state

    viewer = ctx.viewer
    sdata = ctx.sdata
    data_path = ctx.data_path
    no_cache = ctx.no_cache
    pixel_size = ctx.pixel_size
    centroids_yx = ctx.centroids_yx

    arms_status_label = StatusProxy(viewer)

    arms_load_he_button = PushButton(label="Load ARMS H&E Image...", enabled=True)
    arms_flip_v = CheckBox(label="Flip vertically", value=False)
    arms_flip_h = CheckBox(label="Flip horizontally", value=False)
    arms_opacity_slider = Slider(label="H&E opacity", min=0, max=100, value=70)
    arms_opacity_slider.enabled = False

    arms_add_xenium_lm_button = PushButton(label="Add Xenium Landmark", enabled=False)
    arms_add_he_lm_button = PushButton(label="Add ARMS H&E Landmark", enabled=False)
    arms_clear_lm_button = PushButton(label="Clear All", enabled=False)
    arms_register_button = PushButton(label="Compute Registration", enabled=False)
    arms_save_lm_button = PushButton(label="Save Landmarks...", enabled=False)
    arms_load_lm_button = PushButton(label="Load Landmarks...", enabled=True)

    arms_residuals_qt = QTextEdit()
    arms_residuals_qt.setReadOnly(True)
    arms_residuals_qt.setFontFamily("monospace")
    arms_residuals_qt.setMaximumHeight(180)

    arms_load_geojson_button = PushButton(label="Load GeoJSON + CSV...", enabled=False)
    arms_tile_opacity_slider = Slider(label="Tile opacity", min=0, max=100, value=50)
    arms_tile_opacity_slider.enabled = False
    arms_outline_only_check = CheckBox(label="Outline only", value=False, enabled=False)
    arms_edge_width_slider = Slider(label="Tile edge width", min=1, max=100, value=20)
    arms_edge_width_slider.enabled = False

    # ARMS Tile DEG widgets
    arms_deg_method = ComboBox(label="DEG Method", choices=["wilcoxon", "t-test"], value="wilcoxon")
    arms_deg_filter_check = CheckBox(label="Filter by cluster", value=False)
    arms_deg_button = PushButton(label="Run ARMS Tile DEG", enabled=False)
    arms_deg_text = QTextEdit()
    arms_deg_text.setReadOnly(True)
    arms_deg_text.setFontFamily("monospace")
    arms_deg_text.setMaximumHeight(250)
    arms_deg_export_button = PushButton(label="Export ARMS DEG CSV...", enabled=False)
    arms_volcano_button = PushButton(label="Generate ARMS Volcano Plots...", enabled=False)

    # Cluster filter
    arms_cluster_filter_container = QWidget()
    arms_cluster_filter_grid = QGridLayout()
    arms_cluster_filter_grid.setContentsMargins(0, 0, 0, 0)
    arms_cluster_filter_container.setLayout(arms_cluster_filter_grid)

    arms_cluster_scroll = QScrollArea()
    arms_cluster_scroll.setWidget(arms_cluster_filter_container)
    arms_cluster_scroll.setWidgetResizable(True)
    arms_cluster_scroll.setMaximumHeight(120)

    arms_select_all_btn = PushButton(label="Select All", enabled=False)
    arms_deselect_all_btn = PushButton(label="Deselect All", enabled=False)

    # ── Flip / Affine helpers ────────────────────────────────────────────
    def _build_arms_flip_affine():
        shape = arms_state.get("he_shape_yx")
        if shape is None:
            return np.eye(3)
        h, w = shape
        M = np.eye(3)
        if arms_flip_v.value:
            M = np.array([[-1, 0, h - 1], [0, 1, 0], [0, 0, 1]], dtype=np.float64) @ M
        if arms_flip_h.value:
            M = np.array([[1, 0, 0], [0, -1, w - 1], [0, 0, 1]], dtype=np.float64) @ M
        return M

    def _apply_arms_affine():
        flip = _build_arms_flip_affine()
        fine = arms_state["affine_3x3"]
        if fine is not None:
            combined = fine @ flip
        else:
            combined = flip
        if arms_state["he_layer"] is not None:
            arms_state["he_layer"].affine = combined
        if arms_state["he_lm_layer"] is not None:
            arms_state["he_lm_layer"].affine = combined
        if arms_state["shapes_layer"] is not None:
            arms_state["shapes_layer"].affine = combined

    def on_arms_flip_changed(_value=None):
        _apply_arms_affine()
        arms_state["flip_v"] = arms_flip_v.value
        arms_state["flip_h"] = arms_flip_h.value
        _save_arms_affine_to_sdata()
        flips = []
        if arms_flip_v.value: flips.append("V")
        if arms_flip_h.value: flips.append("H")
        if flips:
            arms_status_label.value = f"ARMS flip applied: {'+'.join(flips)}"
        ctx.record_code(
            f"\n# ARMS H&E image flip\n"
            f"# flip_vertical={arms_flip_v.value}, flip_horizontal={arms_flip_h.value}"
        )

    arms_flip_v.changed.connect(on_arms_flip_changed)
    arms_flip_h.changed.connect(on_arms_flip_changed)

    # ── Landmarks ────────────────────────────────────────────────────────
    def _check_arms_landmark_count(*_args):
        xen = arms_state["xenium_lm_layer"]
        he = arms_state["he_lm_layer"]
        if xen is not None and he is not None:
            n = min(len(xen.data), len(he.data))
            arms_register_button.enabled = n >= 3
            arms_save_lm_button.enabled = n >= 1

    def _create_arms_landmark_layers():
        if arms_state["xenium_lm_layer"] is not None:
            return
        xen_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="ARMS Xenium Landmarks", size=30, face_color="cyan",
            symbol="cross", border_color="cyan",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )
        he_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="ARMS H&E Landmarks", size=30, face_color="orange",
            symbol="cross", border_color="orange",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )
        xen_lm.events.data.connect(_check_arms_landmark_count)
        he_lm.events.data.connect(_check_arms_landmark_count)
        arms_state["xenium_lm_layer"] = xen_lm
        arms_state["he_lm_layer"] = he_lm
        arms_add_xenium_lm_button.enabled = True
        arms_add_he_lm_button.enabled = True
        arms_clear_lm_button.enabled = True

    # ── Load H&E ─────────────────────────────────────────────────────────
    def _on_arms_he_loaded(result):
        (pyramid, tif), path = result
        if arms_state["he_layer"] is not None:
            try:
                viewer.layers.remove(arms_state["he_layer"])
            except ValueError:
                pass
        arms_state["he_tif"] = tif
        arms_state["he_filename"] = Path(path).name
        arms_state["he_path"] = str(path)
        base = pyramid[0]
        arms_state["he_shape_yx"] = (base.shape[0], base.shape[1])
        he_layer = viewer.add_image(
            pyramid, name=f"ARMS H&E ({Path(path).name})",
            rgb=True, blending="translucent",
            opacity=arms_opacity_slider.value / 100.0,
        )
        arms_state["he_layer"] = he_layer
        arms_state["affine_3x3"] = None
        _apply_arms_affine()
        _create_arms_landmark_layers()
        _save_arms_he_to_sdata(pyramid, Path(path).name)
        arms_opacity_slider.enabled = True
        arms_load_he_button.enabled = True
        arms_load_geojson_button.enabled = True
        shape_str = "x".join(str(s) for s in pyramid[0].shape)
        arms_status_label.value = f"ARMS H&E loaded: {Path(path).name} ({shape_str}, {len(pyramid)} levels)"
        ctx.record_code(
            f"\n# Load ARMS H&E image\n"
            f"from xenium_viewer.utils.registration import load_he_pyramid\n"
            f"arms_pyramid, arms_tif = load_he_pyramid(\"{path}\")"
        )

    def on_arms_load_he():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load ARMS H&E Image", default_dir,
            "Image Files (*.ome.tif *.tif *.tiff *.svs);;All Files (*)",
        )
        if not path:
            return
        arms_status_label.value = "Loading ARMS H&E..."
        arms_load_he_button.enabled = False
        gen = ctx.dataset_generation

        @thread_worker
        def load_task():
            return load_he_pyramid(path), path

        worker = load_task()
        worker.returned.connect(lambda result: _on_arms_he_loaded(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def on_arms_he_opacity(value):
        if arms_state["he_layer"] is not None:
            arms_state["he_layer"].opacity = value / 100.0
    arms_opacity_slider.changed.connect(on_arms_he_opacity)

    def on_arms_add_xenium_lm():
        lm = arms_state["xenium_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            arms_status_label.value = "Click on a feature in the Xenium image"

    def on_arms_add_he_lm():
        lm = arms_state["he_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            arms_status_label.value = "Click on the same feature in the ARMS H&E image"

    def on_arms_clear_lm():
        for key in ("xenium_lm_layer", "he_lm_layer"):
            lm = arms_state[key]
            if lm is not None:
                lm.selected_data = set()
                lm.data = np.empty((0, 2), dtype=np.float64)
        arms_state["affine_3x3"] = None
        _apply_arms_affine()
        _save_arms_affine_to_sdata()
        arms_residuals_qt.clear()
        arms_status_label.value = "ARMS landmarks cleared"
        arms_register_button.enabled = False
        arms_save_lm_button.enabled = False
        from xenium_viewer.utils.adata_persistence import save_landmarks_to_sdata
        save_landmarks_to_sdata(ctx, 'arms_xenium_landmarks', None)
        save_landmarks_to_sdata(ctx, 'arms_he_landmarks', None)

    def on_arms_register():
        xen_pts = arms_state["xenium_lm_layer"].data
        he_pts = arms_state["he_lm_layer"].data
        n = min(len(xen_pts), len(he_pts))
        if n < 3:
            arms_status_label.value = "Need at least 3 paired landmarks"
            return
        xen_pts = np.asarray(xen_pts[:n], dtype=np.float64)
        he_pts = np.asarray(he_pts[:n], dtype=np.float64)
        affine, residuals = compute_landmark_affine(xen_pts, he_pts)
        arms_state["affine_3x3"] = affine
        _apply_arms_affine()
        lines = [f"ARMS Registration: {n} landmarks, similarity transform"]
        lines.append(f"Mean residual: {residuals.mean():.1f} px ({residuals.mean() * pixel_size:.1f} um)")
        lines.append(f"Max  residual: {residuals.max():.1f} px ({residuals.max() * pixel_size:.1f} um)")
        lines.append("")
        for i, r in enumerate(residuals):
            lines.append(f"  Landmark {i+1}: {r:.1f} px ({r * pixel_size:.1f} um)")
        scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
        lines.append(f"\nScale factor: {scale:.4f}")
        arms_residuals_qt.setPlainText("\n".join(lines))
        arms_status_label.value = f"ARMS registered ({n} landmarks, mean residual {residuals.mean():.1f} px)"
        _save_arms_affine_to_sdata()
        from xenium_viewer.utils.adata_persistence import save_landmarks_to_sdata
        save_landmarks_to_sdata(ctx, 'arms_xenium_landmarks', xen_pts)
        save_landmarks_to_sdata(ctx, 'arms_he_landmarks', he_pts)
        ctx.record_code(
            f"\n# ARMS landmark registration\n"
            f"from xenium_viewer.utils.registration import compute_landmark_affine\n"
            f"arms_xen_pts = np.array({xen_pts.tolist()})\n"
            f"arms_he_pts = np.array({he_pts.tolist()})\n"
            f"arms_affine, arms_residuals = compute_landmark_affine(arms_xen_pts, arms_he_pts)"
        )

    def on_arms_save_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save ARMS Landmarks", default_dir + "/arms_landmarks.json", "JSON Files (*.json)",
        )
        if not path:
            return
        xen_pts = np.asarray(arms_state["xenium_lm_layer"].data, dtype=np.float64)
        he_pts = np.asarray(arms_state["he_lm_layer"].data, dtype=np.float64)
        save_landmarks(
            path, xen_pts, he_pts,
            affine=arms_state["affine_3x3"], he_filename=arms_state.get("he_filename"),
        )
        arms_status_label.value = f"Landmarks saved to {Path(path).name}"
        ctx.record_code(f"\n# Save ARMS landmarks\n# landmarks -> \"{path}\"")

    def on_arms_load_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load ARMS Landmarks", default_dir, "JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        lm_data = load_landmarks(path)
        _create_arms_landmark_layers()
        arms_state["xenium_lm_layer"].data = lm_data["xenium_landmarks_yx"]
        arms_state["he_lm_layer"].data = lm_data["he_landmarks_yx"]
        if "affine_3x3_yx" in lm_data:
            affine = lm_data["affine_3x3_yx"]
            arms_state["affine_3x3"] = affine
            _apply_arms_affine()
            scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
            arms_residuals_qt.setPlainText(f"Loaded affine (scale={scale:.4f})")
        if "he_filename" in lm_data:
            arms_state["he_filename"] = lm_data["he_filename"]
        n = min(len(lm_data["xenium_landmarks_yx"]), len(lm_data["he_landmarks_yx"]))
        arms_status_label.value = f"Loaded {n} landmarks from {Path(path).name}"
        ctx.record_code(
            f"\n# Load ARMS landmarks from file\n"
            f"from xenium_viewer.utils.registration import load_landmarks\n"
            f"landmarks = load_landmarks(\"{path}\")"
        )

    # ── Save H&E / affine to sdata ──────────────────────────────────────
    def _save_arms_he_to_sdata(pyramid, he_filename):
        if sdata is None or no_cache:
            return
        try:
            from spatialdata.models import Image2DModel
            base = np.asarray(pyramid[0])
            if base.ndim == 3 and base.shape[-1] in (3, 4):
                base_cyx = np.transpose(base, (2, 0, 1))
            else:
                base_cyx = base
            parsed = Image2DModel.parse(
                base_cyx.astype(np.uint8), dims=("c", "y", "x"),
                scale_factors=[2, 2, 2, 2], chunks=(3, 1024, 1024),
            )
            if "arms_he_image" in sdata.images:
                del sdata.images["arms_he_image"]
            sdata.images["arms_he_image"] = parsed
            sdata.write_element("arms_he_image", overwrite=True)
            zarr_p = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_p), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            store["viewer_session"].attrs["arms_he_filename"] = he_filename
            store["viewer_session"].attrs["arms_he_shape_yx"] = list(base.shape[:2])
            print(f"  ARMS H&E image saved to sdata zarr cache ({base_cyx.shape})")
        except Exception as e:
            print(f"  Warning: could not save ARMS H&E to sdata: {e}")

    def _save_arms_affine_to_sdata():
        if sdata is None or no_cache or "arms_he_image" not in sdata.images:
            return
        try:
            from spatialdata.transformations import Affine as SdAffine, set_transformation
            flip = _build_arms_flip_affine()
            fine = arms_state["affine_3x3"]
            if fine is not None:
                combined = fine @ flip
            else:
                combined = flip
            sd_affine = SdAffine(combined, input_axes=("y", "x"), output_axes=("y", "x"))
            set_transformation(sdata.images["arms_he_image"], sd_affine, "global")
            sdata.write_transformations("arms_he_image")
            zarr_p = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_p), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            sess = store["viewer_session"]
            sess.attrs["arms_flip_v"] = bool(arms_state.get("flip_v", False))
            sess.attrs["arms_flip_h"] = bool(arms_state.get("flip_h", False))
            if fine is not None:
                sess.attrs["arms_affine_3x3"] = fine.tolist()
            sess.attrs["arms_he_filename"] = arms_state.get("he_filename")
            sess.attrs["arms_he_path"] = arms_state.get("he_path")
            sess.attrs["arms_he_shape_yx"] = (
                list(arms_state["he_shape_yx"]) if arms_state.get("he_shape_yx") else None
            )
            sess.attrs["arms_geojson_path"] = arms_state.get("geojson_path")
            sess.attrs["arms_csv_path"] = arms_state.get("csv_path")
        except Exception as e:
            print(f"  Warning: could not save ARMS affine: {e}")

    # ── GeoJSON/CSV loading ──────────────────────────────────────────────
    def _load_geojson_csv(geojson_path, csv_path):
        import csv as csv_mod
        import re as _re

        try:
            with open(geojson_path) as f:
                geojson = json.load(f)
            features = geojson.get("features", [])
            tile_polygons = {}
            for feat in features:
                name = feat.get("properties", {}).get("name")
                if name is None:
                    continue
                coords_xy = feat["geometry"]["coordinates"]
                if feat["geometry"]["type"] == "MultiPolygon":
                    ring = coords_xy[0][0]
                elif feat["geometry"]["type"] == "Polygon":
                    ring = coords_xy[0]
                else:
                    continue
                arr = np.array(ring, dtype=np.float64)
                tile_polygons[name] = arr[:, ::-1]

            tile_clusters = {}
            with open(csv_path) as f:
                reader = csv_mod.DictReader(f)
                for row in reader:
                    tile_name = row.get("tile", row.get("sample", "")).strip()
                    cluster_val = row.get("cluster", "").strip()
                    if tile_name and cluster_val:
                        try:
                            tile_clusters[tile_name] = int(cluster_val)
                        except ValueError:
                            pass

            def _normalize_tile(name):
                n = name.strip().lower()
                n = _re.sub(r'^plate', 'p', n)
                return n

            csv_norm = {}
            for raw_name, cid in tile_clusters.items():
                csv_norm[_normalize_tile(raw_name)] = cid

            polygon_data = []
            face_colors = []
            tile_names = []
            cluster_ids = []
            matched = 0
            for name, poly_yx in tile_polygons.items():
                cid = tile_clusters.get(name)
                if cid is None:
                    cid = csv_norm.get(_normalize_tile(name))
                if cid is None:
                    continue
                matched += 1
                polygon_data.append(poly_yx)
                tile_names.append(name)
                cluster_ids.append(cid)
                idx = max(0, min(cid - 1, len(ARMS_CLUSTER_PALETTE) - 1))
                face_colors.append(ARMS_CLUSTER_PALETTE[idx])

            if not polygon_data:
                geojson_sample = list(tile_polygons.keys())[:5]
                csv_sample = list(tile_clusters.keys())[:5]
                arms_status_label.value = "No matching tiles found between GeoJSON and CSV"
                arms_residuals_qt.setPlainText(
                    f"No matching tiles found.\n\n"
                    f"GeoJSON tile names (sample): {geojson_sample}\n"
                    f"CSV tile names (sample): {csv_sample}\n\n"
                    f"GeoJSON has {len(tile_polygons)} tiles, CSV has {len(tile_clusters)} entries."
                )
                return False

            if arms_state["shapes_layer"] is not None:
                try:
                    viewer.layers.remove(arms_state["shapes_layer"])
                except ValueError:
                    pass

            cluster_ids_arr = np.array(cluster_ids, dtype=int)
            if arms_outline_only_check.value:
                init_face = np.zeros((len(cluster_ids_arr), 4), dtype=np.float32)
                init_edge = np.array(face_colors)
            else:
                init_face = np.array(face_colors)
                init_edge = np.ones((len(cluster_ids_arr), 4), dtype=np.float32)
            shapes_layer = viewer.add_shapes(
                polygon_data, shape_type="polygon",
                face_color=init_face, edge_color=init_edge,
                edge_width=20, name="ARMS Tiles",
                opacity=arms_tile_opacity_slider.value / 100.0,
            )
            arms_state["shapes_layer"] = shapes_layer
            arms_state["tile_names"] = tile_names
            arms_state["cluster_ids"] = cluster_ids_arr
            arms_state["geojson_path"] = geojson_path
            arms_state["csv_path"] = csv_path

            _apply_arms_affine()
            arms_tile_opacity_slider.enabled = True
            arms_outline_only_check.enabled = True
            arms_edge_width_slider.enabled = True
            arms_deg_button.enabled = True

            # Populate cluster filter
            unique_clusters = sorted(set(cluster_ids))
            while arms_cluster_filter_grid.count():
                item = arms_cluster_filter_grid.takeAt(0)
                w = item.widget()
                if w:
                    w.deleteLater()
            arms_state["cluster_checkboxes"].clear()

            cols = 3
            for i, cid in enumerate(unique_clusters):
                cb = QCheckBox(f"C{cid}")
                cb.setChecked(True)
                arms_cluster_filter_grid.addWidget(cb, i // cols, i % cols)
                arms_state["cluster_checkboxes"][cid] = cb

            arms_select_all_btn.enabled = True
            arms_deselect_all_btn.enabled = True

            legend_lines = [f"Loaded {matched} tiles from {len(tile_polygons)} GeoJSON features"]
            legend_lines.append(f"Clusters found: {unique_clusters}")
            legend_lines.append("")
            legend_lines.append("Cluster legend:")
            for cid in unique_clusters:
                count = cluster_ids.count(cid)
                legend_lines.append(f"  Cluster {cid} ({count} tiles)")
            arms_residuals_qt.setPlainText("\n".join(legend_lines))
            arms_status_label.value = f"ARMS tiles loaded: {matched} tiles, {len(unique_clusters)} clusters"
            return True

        except Exception as e:
            arms_status_label.value = f"Error loading GeoJSON/CSV: {e}"
            import traceback
            traceback.print_exc()
            return False

    def on_arms_load_geojson():
        default_dir = str(data_path) if data_path else ""
        geojson_path, _ = QFileDialog.getOpenFileName(
            None, "Load GeoJSON (tile boundaries)", default_dir,
            "GeoJSON Files (*.geojson *.json);;All Files (*)",
        )
        if not geojson_path:
            return
        csv_path, _ = QFileDialog.getOpenFileName(
            None, "Load CSV (tile cluster IDs)", str(Path(geojson_path).parent),
            "CSV Files (*.csv);;All Files (*)",
        )
        if not csv_path:
            return
        arms_status_label.value = "Loading GeoJSON + CSV..."
        ok = _load_geojson_csv(geojson_path, csv_path)
        if ok and not no_cache:
            try:
                zarr_p = data_path / "sdata_cached.zarr"
                if zarr_p.exists():
                    import zarr as zarr_mod
                    store = zarr_mod.open_group(str(zarr_p), mode="r+", use_consolidated=False)
                    if "viewer_session" not in store:
                        store.create_group("viewer_session")
                    store["viewer_session"].attrs["arms_geojson_path"] = geojson_path
                    store["viewer_session"].attrs["arms_csv_path"] = csv_path
            except Exception:
                pass
            from xenium_viewer.utils.adata_persistence import save_arms_tiles_to_sdata
            polys = list(arms_state["shapes_layer"].data) if arms_state["shapes_layer"] is not None else []
            save_arms_tiles_to_sdata(
                ctx, polys, arms_state.get("tile_names", []), arms_state.get("cluster_ids", []),
            )
        ctx.record_code(
            f"\n# Load ARMS tile boundaries + cluster assignments\n"
            f"import json\n"
            f"arms_geojson_path = \"{geojson_path}\"\n"
            f"arms_csv_path = \"{csv_path}\""
        )

    def on_arms_tile_opacity(value):
        if arms_state["shapes_layer"] is not None:
            arms_state["shapes_layer"].opacity = value / 100.0
    arms_tile_opacity_slider.changed.connect(on_arms_tile_opacity)

    def on_arms_edge_width(value):
        layer = arms_state["shapes_layer"]
        if layer is not None:
            layer.edge_width = value
            layer.refresh()
    arms_edge_width_slider.changed.connect(on_arms_edge_width)

    def _apply_tile_display_mode():
        layer = arms_state["shapes_layer"]
        cluster_ids_arr = arms_state.get("cluster_ids")
        if layer is None or cluster_ids_arr is None:
            return
        cluster_colors = np.array([
            ARMS_CLUSTER_PALETTE[max(0, min(int(cid) - 1, len(ARMS_CLUSTER_PALETTE) - 1))]
            for cid in cluster_ids_arr
        ])
        if arms_outline_only_check.value:
            layer.face_color = np.zeros((len(cluster_ids_arr), 4), dtype=np.float32)
            layer.edge_color = cluster_colors
        else:
            layer.face_color = cluster_colors
            layer.edge_color = np.ones((len(cluster_ids_arr), 4), dtype=np.float32)  # white

    arms_outline_only_check.changed.connect(lambda _: _apply_tile_display_mode())

    # ── ARMS Tile DEG ────────────────────────────────────────────────────
    def on_arms_deg():
        if arms_state["shapes_layer"] is None or arms_state["cluster_ids"] is None:
            arms_status_label.value = "Load ARMS tiles first"
            return

        tile_data = list(arms_state["shapes_layer"].data)
        cluster_ids_arr = arms_state["cluster_ids"].copy()

        if arms_state["cluster_checkboxes"]:
            selected = {cid for cid, cb in arms_state["cluster_checkboxes"].items() if cb.isChecked()}
            if len(selected) < 2:
                arms_status_label.value = "Select at least 2 ARMS clusters for DEG"
                return
            keep_mask = np.array([cid in selected for cid in cluster_ids_arr])
            tile_data = [tile_data[i] for i in range(len(tile_data)) if keep_mask[i]]
            cluster_ids_arr = cluster_ids_arr[keep_mask]

        affine_mat = arms_state["shapes_layer"].affine.affine_matrix

        transformed_polys = []
        for poly_yx in tile_data:
            ones = np.ones((len(poly_yx), 1))
            coords_h = np.hstack([poly_yx, ones])
            transformed = (affine_mat @ coords_h.T).T[:, :2]
            transformed_polys.append(transformed)

        arms_status_label.value = "Running ARMS Tile DEG..."
        arms_deg_button.enabled = False
        method = arms_deg_method.value
        _adata = ctx.adata if ctx.adata is not None else ctx.color_manager.adata
        gen = ctx.dataset_generation

        cluster_mask = None
        if arms_deg_filter_check.value:
            clustering_key = ctx.clustering_widget.value
            selected_ids = ctx.get_selected_cluster_ids()
            cluster_series = ctx.clusterings[clustering_key]
            if 'cell_id' in _adata.obs.columns:
                clusters_aligned = cluster_series.reindex(_adata.obs['cell_id'].values)
            else:
                clusters_aligned = cluster_series.reindex(_adata.obs_names)
            cluster_mask = ctx.make_cluster_mask(clusters_aligned.values, selected_ids)

        selected_clusters = sorted(cid for cid, cb in arms_state["cluster_checkboxes"].items() if cb.isChecked()) if arms_state["cluster_checkboxes"] else "all"

        # Build code snippet for recorder
        _fine = arms_state.get("affine_3x3")
        _flip_v = arms_flip_v.value
        _flip_h = arms_flip_h.value
        _shape = arms_state.get("he_shape_yx")  # (h, w) or None
        _fine_repr = f"np.array({np.asarray(_fine, dtype=np.float64).tolist()!r})" if _fine is not None else "np.eye(3)"
        _shape_repr = f"({int(_shape[0])}, {int(_shape[1])})" if _shape else "(0, 0)"
        _sel_clusters_repr = repr(list(selected_clusters)) if selected_clusters != "all" else repr([])

        _xenium_filter_line = ""
        if arms_deg_filter_check.value:
            _ck = clustering_key
            _si = sorted(selected_ids)
            _xenium_filter_line = (
                f"# Xenium cluster filter (clustering={_ck!r}, clusters={_si})\n"
                f"cluster_series = adata.obs[{_ck!r}]\n"
                f"clusters_aligned = cluster_series.reindex(adata.obs['cell_id'].values)\n"
                f"cluster_mask = np.isin(clusters_aligned.values, {_si})\n"
            )

        ctx.record_code(
            f"\n# ARMS Tile DEG analysis\n"
            f"import json, numpy as np\n"
            f"from xenium_viewer.utils.adata_persistence import load_arms_tiles_from_sdata\n"
            f"from xenium_viewer.utils.gene_analysis import compute_arms_tile_deg\n"
            f"pixel_size = float(json.load(open(data_path / 'experiment.xenium'))['pixel_size'])\n"
            f"centroids_yx = adata.obsm['spatial'][:, ::-1] / pixel_size  # µm→px, xy→yx\n"
            f"arms_polys_yx, arms_tile_names, arms_cluster_ids = load_arms_tiles_from_sdata(sdata)\n"
            f"# Apply ARMS registration affine (fine @ flip)\n"
            f"h, w = {_shape_repr}\n"
            f"M = np.eye(3)\n"
            f"if {_flip_v}:  # flip_v\n"
            f"    M = np.array([[-1, 0, h-1], [0, 1, 0], [0, 0, 1]], dtype=float) @ M\n"
            f"if {_flip_h}:  # flip_h\n"
            f"    M = np.array([[1, 0, 0], [0, -1, w-1], [0, 0, 1]], dtype=float) @ M\n"
            f"affine_3x3 = {_fine_repr}\n"
            f"combined = affine_3x3 @ M\n"
            f"transformed_polys = []\n"
            f"for poly_yx in arms_polys_yx:\n"
            f"    ones = np.ones((len(poly_yx), 1))\n"
            f"    transformed_polys.append((combined @ np.hstack([poly_yx, ones]).T).T[:, :2])\n"
            f"# Select tile clusters: {list(selected_clusters) if selected_clusters != 'all' else 'all'}\n"
            f"_sel = {_sel_clusters_repr}\n"
            f"if _sel:\n"
            f"    _keep = np.isin(arms_cluster_ids, _sel)\n"
            f"    transformed_polys = [transformed_polys[i] for i in range(len(transformed_polys)) if _keep[i]]\n"
            f"    arms_cluster_ids = arms_cluster_ids[_keep]\n"
            f"{_xenium_filter_line}"
            f"arms_deg_df, arms_summary, arms_adata_norm = compute_arms_tile_deg(\n"
            f"    adata, centroids_yx, transformed_polys, arms_cluster_ids,\n"
            f"    method={method!r},\n"
            f"    cluster_mask={'cluster_mask' if arms_deg_filter_check.value else 'None'},\n"
            f")"
        )

        @thread_worker
        def _run():
            return compute_arms_tile_deg(
                _adata, centroids_yx, transformed_polys, cluster_ids_arr,
                method=method, cluster_mask=cluster_mask,
            )

        worker = _run()
        worker.returned.connect(lambda result: _on_arms_deg_ready(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def _on_arms_deg_ready(result):
        deg_df, summary, adata_norm = result
        state["arms_tile_deg_df"] = deg_df
        arms_deg_button.enabled = True
        if deg_df.empty:
            summary_str = ", ".join(f"C{k}: {v} cells" for k, v in sorted(summary.items()))
            arms_deg_text.setPlainText(
                f"No results (need ≥2 clusters with ≥10 cells each).\n"
                f"Cells per cluster: {summary_str or 'none'}"
            )
            arms_status_label.value = "ARMS DEG: no results"
            arms_deg_export_button.enabled = False
            arms_volcano_button.enabled = False
            return
        state["arms_deg_adata_norm"] = adata_norm
        summary_str = ", ".join(f"C{k}: {v}" for k, v in sorted(summary.items()))
        preview = deg_df.head(50).to_string(index=False)
        arms_deg_text.setPlainText(f"Cells per cluster: {summary_str}\n\n{preview}")
        arms_status_label.value = f"ARMS DEG complete: {len(deg_df)} gene-group results"
        arms_deg_export_button.enabled = True
        arms_volcano_button.enabled = True
        from xenium_viewer.utils.adata_persistence import save_arms_tile_deg_to_sdata
        save_arms_tile_deg_to_sdata(ctx, deg_df)

    def on_export_arms_deg():
        df = state.get("arms_tile_deg_df")
        if df is None or df.empty:
            return
        path, _ = QFileDialog.getSaveFileName(
            None, "Export ARMS DEG Results", "arms_tile_deg_results.csv", "CSV Files (*.csv)",
        )
        if not path:
            return
        df.to_csv(path, index=False)
        arms_status_label.value = f"Exported {len(df)} rows to {path}"
        ctx.record_code(f"\n# Export ARMS DEG results\narms_deg_df.to_csv(\"{path}\", index=False)")

    def on_arms_generate_volcanos():
        adata_norm = state.get("arms_deg_adata_norm")
        if adata_norm is None:
            arms_status_label.value = "Run ARMS Tile DEG first"
            return
        output_dir = QFileDialog.getExistingDirectory(None, "Select output directory for ARMS volcano plots")
        if not output_dir:
            return
        arms_volcano_button.enabled = False
        arms_status_label.value = "Generating ARMS volcano plots..."
        method = arms_deg_method.value
        gen = ctx.dataset_generation
        ctx.record_code(
            f"\n# Generate ARMS pairwise volcano plots\n"
            f"from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot\n"
            f"import itertools\n"
            f"arms_volcano_dir = \"{output_dir}\"\n"
            f"# Uses arms_adata_norm from DEG step, method=\"{method}\""
        )

        @thread_worker
        def _run():
            from pathlib import Path
            import itertools as _it
            from xenium_viewer.utils.gene_analysis import run_pairwise_deg, make_volcano_plot
            import matplotlib.pyplot as _plt

            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            groups = sorted(
                [g for g in adata_norm.obs['arms_cluster'].cat.categories if str(g) != '-1'],
                key=lambda x: (int(x) if str(x).lstrip('-').isdigit() else 0, str(x)),
            )
            pairs = list(_it.combinations(groups, 2))
            total = len(pairs)
            for i, (a, b) in enumerate(pairs):
                yield f"Volcano plot {i + 1}/{total}: ARMS C{a} vs C{b}"
                df = run_pairwise_deg(adata_norm, 'arms_cluster', str(a), str(b), method=method)
                fig = make_volcano_plot(df, f"ARMS C{a}", f"ARMS C{b}")
                fig.savefig(out / f'volcano_ARMS_C{a}_vs_ARMS_C{b}.png', dpi=300)
                _plt.close(fig)
            return total, output_dir

        worker = _run()
        worker.yielded.connect(lambda msg: setattr(arms_status_label, 'value', msg) if ctx.dataset_generation == gen else None)
        worker.returned.connect(lambda result: _on_arms_volcanos_done(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def _on_arms_volcanos_done(result):
        count, out_dir = result
        arms_volcano_button.enabled = True
        arms_status_label.value = f"{count} ARMS volcano plots saved to {out_dir}"

    def _on_arms_select_all():
        for cb in arms_state["cluster_checkboxes"].values():
            cb.setChecked(True)

    def _on_arms_deselect_all():
        for cb in arms_state["cluster_checkboxes"].values():
            cb.setChecked(False)

    arms_select_all_btn.clicked.connect(_on_arms_select_all)
    arms_deselect_all_btn.clicked.connect(_on_arms_deselect_all)
    arms_deg_button.clicked.connect(on_arms_deg)
    arms_deg_export_button.clicked.connect(on_export_arms_deg)
    arms_volcano_button.clicked.connect(on_arms_generate_volcanos)
    arms_load_he_button.clicked.connect(on_arms_load_he)
    arms_add_xenium_lm_button.clicked.connect(on_arms_add_xenium_lm)
    arms_add_he_lm_button.clicked.connect(on_arms_add_he_lm)
    arms_clear_lm_button.clicked.connect(on_arms_clear_lm)
    arms_register_button.clicked.connect(on_arms_register)
    arms_save_lm_button.clicked.connect(on_arms_save_landmarks)
    arms_load_lm_button.clicked.connect(on_arms_load_landmarks)
    arms_load_geojson_button.clicked.connect(on_arms_load_geojson)

    # Layout rows
    arms_flip_row = QWidget()
    arms_flip_layout = QHBoxLayout()
    arms_flip_layout.setContentsMargins(0, 0, 0, 0)
    arms_flip_layout.addWidget(arms_flip_v.native)
    arms_flip_layout.addWidget(arms_flip_h.native)
    arms_flip_row.setLayout(arms_flip_layout)

    arms_lm_btn_row = QWidget()
    arms_lm_layout = QHBoxLayout()
    arms_lm_layout.setContentsMargins(0, 0, 0, 0)
    arms_lm_layout.addWidget(arms_add_xenium_lm_button.native)
    arms_lm_layout.addWidget(arms_add_he_lm_button.native)
    arms_lm_layout.addWidget(arms_clear_lm_button.native)
    arms_lm_btn_row.setLayout(arms_lm_layout)

    arms_io_btn_row = QWidget()
    arms_io_btn_layout = QHBoxLayout()
    arms_io_btn_layout.setContentsMargins(0, 0, 0, 0)
    arms_io_btn_layout.addWidget(arms_save_lm_button.native)
    arms_io_btn_layout.addWidget(arms_load_lm_button.native)
    arms_io_btn_row.setLayout(arms_io_btn_layout)

    arms_cluster_btn_row = QWidget()
    arms_cluster_btn_layout = QHBoxLayout()
    arms_cluster_btn_layout.setContentsMargins(0, 0, 0, 0)
    arms_cluster_btn_layout.addWidget(arms_select_all_btn.native)
    arms_cluster_btn_layout.addWidget(arms_deselect_all_btn.native)
    arms_cluster_btn_row.setLayout(arms_cluster_btn_layout)

    widget = make_tab(
        arms_load_he_button,
        arms_flip_row,
        arms_opacity_slider,
        arms_lm_btn_row,
        arms_register_button,
        arms_io_btn_row,
        arms_residuals_qt,
        arms_load_geojson_button,
        arms_tile_opacity_slider,
        arms_outline_only_check,
        arms_edge_width_slider,
        arms_cluster_btn_row,
        arms_cluster_scroll,
        arms_deg_filter_check,
        arms_deg_method,
        arms_deg_button,
        arms_deg_text,
        arms_deg_export_button,
        arms_volcano_button,
    )

    # ── Session restore ──────────────────────────────────────────────────
    def _on_arms_restored(pyramid_rgb, session_arms_data):
        if arms_state["he_layer"] is not None:
            try:
                viewer.layers.remove(arms_state["he_layer"])
            except ValueError:
                pass
        he_filename = session_arms_data.get("he_filename", "ARMS H&E")
        arms_state["he_tif"] = None
        arms_state["he_filename"] = he_filename
        arms_state["he_path"] = session_arms_data.get("he_path")
        base = pyramid_rgb[0]
        arms_state["he_shape_yx"] = session_arms_data.get("he_shape_yx") or (base.shape[0], base.shape[1])
        arms_state["geojson_path"] = session_arms_data.get("geojson_path")
        arms_state["csv_path"] = session_arms_data.get("csv_path")
        he_layer = viewer.add_image(
            pyramid_rgb, name=f"ARMS H&E ({he_filename})",
            rgb=True, blending="translucent",
            opacity=arms_opacity_slider.value / 100.0,
        )
        arms_state["he_layer"] = he_layer
        arms_state["affine_3x3"] = session_arms_data.get("affine_3x3")
        _apply_arms_affine()
        _save_arms_affine_to_sdata()
        _create_arms_landmark_layers()
        xen_lm = session_arms_data.get("xenium_landmarks")
        he_lm = session_arms_data.get("he_landmarks")
        if xen_lm is not None and arms_state["xenium_lm_layer"] is not None:
            arms_state["xenium_lm_layer"].data = xen_lm
        if he_lm is not None and arms_state["he_lm_layer"] is not None:
            arms_state["he_lm_layer"].data = he_lm
        arms_opacity_slider.enabled = True
        arms_load_he_button.enabled = True
        arms_load_geojson_button.enabled = True

        sdata_tiles = session_arms_data.get("arms_tiles_sdata")
        if sdata_tiles and sdata_tiles[0]:
            # Tiles loaded from sdata.shapes — skip file re-load
            polys_yx, t_names, c_ids = sdata_tiles
            face_colors = []
            for cid in c_ids:
                idx = max(0, min(int(cid) - 1, len(ARMS_CLUSTER_PALETTE) - 1))
                face_colors.append(ARMS_CLUSTER_PALETTE[idx])
            if arms_state["shapes_layer"] is not None:
                try:
                    viewer.layers.remove(arms_state["shapes_layer"])
                except ValueError:
                    pass
            c_ids_arr = np.array(c_ids, dtype=int)
            if arms_outline_only_check.value:
                init_face = np.zeros((len(c_ids_arr), 4), dtype=np.float32)
                init_edge = np.array(face_colors)
            else:
                init_face = np.array(face_colors)
                init_edge = np.ones((len(c_ids_arr), 4), dtype=np.float32)
            shapes_layer = viewer.add_shapes(
                polys_yx, shape_type="polygon",
                face_color=init_face, edge_color=init_edge,
                edge_width=20, name="ARMS Tiles",
                opacity=arms_tile_opacity_slider.value / 100.0,
            )
            arms_state["shapes_layer"] = shapes_layer
            arms_state["tile_names"] = t_names
            arms_state["cluster_ids"] = c_ids_arr
            _apply_arms_affine()
            arms_tile_opacity_slider.enabled = True
            arms_outline_only_check.enabled = True
            arms_edge_width_slider.enabled = True
            arms_deg_button.enabled = True
            unique_clusters = sorted(set(int(c) for c in c_ids))
            while arms_cluster_filter_grid.count():
                item = arms_cluster_filter_grid.takeAt(0)
                w = item.widget()
                if w:
                    w.deleteLater()
            arms_state["cluster_checkboxes"].clear()
            cols = 3
            from qtpy.QtWidgets import QCheckBox as _QCheckBox
            for i, cid in enumerate(unique_clusters):
                cb = _QCheckBox(f"C{cid}")
                cb.setChecked(True)
                arms_cluster_filter_grid.addWidget(cb, i // cols, i % cols)
                arms_state["cluster_checkboxes"][cid] = cb
            arms_select_all_btn.enabled = True
            arms_deselect_all_btn.enabled = True
            print(f"  Restored ARMS tiles from sdata ({len(polys_yx)} tiles)")
        else:
            geojson_path = session_arms_data.get("geojson_path")
            csv_path = session_arms_data.get("csv_path")
            if geojson_path and csv_path:
                if Path(geojson_path).exists() and Path(csv_path).exists():
                    _load_geojson_csv(geojson_path, csv_path)
                else:
                    missing = []
                    if not Path(geojson_path).exists():
                        missing.append(f"GeoJSON: {geojson_path}")
                    if not Path(csv_path).exists():
                        missing.append(f"CSV: {csv_path}")
                    print(f"  Warning: ARMS tile files moved/missing: {', '.join(missing)}")
                    arms_status_label.value = "ARMS H&E restored (tile files not found)"

        has_affine = arms_state["affine_3x3"] is not None
        if arms_status_label.value.startswith("Restoring"):
            arms_status_label.value = f"ARMS restored: {he_filename}" + (" (with registration)" if has_affine else "")
        print(f"  Restored ARMS from cache: {he_filename}" + (" with registration" if has_affine else ""))

    def _restore_session(session):
        # ARMS tile DEG
        atd = session.get("arms_tile_deg_df")
        if atd is not None:
            state["arms_tile_deg_df"] = atd
            arms_deg_export_button.enabled = True
            arms_deg_text.setPlainText(atd.head(50).to_string(index=False))
            print(f"  Restored ARMS Tile DEG ({len(atd)} rows)")

        # ARMS H&E from sdata cache
        if sdata is not None and "arms_he_image" in sdata.images:
            arms_flip_v.value = session.get("arms_flip_v", False)
            arms_flip_h.value = session.get("arms_flip_h", False)
            _session_arms_data = {
                "affine_3x3": session.get("arms_affine_3x3"),
                "xenium_landmarks": session.get("arms_xenium_landmarks"),
                "he_landmarks": session.get("arms_he_landmarks"),
                "he_filename": session.get("arms_he_filename", "ARMS H&E"),
                "he_path": session.get("arms_he_path"),
                "he_shape_yx": session.get("arms_he_shape_yx"),
                "geojson_path": session.get("arms_geojson_path"),
                "csv_path": session.get("arms_csv_path"),
                "arms_tiles_sdata": session.get("arms_tiles_sdata"),
            }
            arms_status_label.value = "Restoring ARMS H&E from cache..."
            gen = ctx.dataset_generation

            from xenium_viewer.tabs.tab_he_registration import _extract_dt_scales

            @thread_worker
            def _load_arms_from_sdata():
                import dask.array as da
                arms_dt = sdata.images["arms_he_image"]
                pyramid = _extract_dt_scales(arms_dt)
                pyramid_rgb = []
                for arr in pyramid:
                    if not isinstance(arr, da.Array):
                        arr = da.from_array(arr)
                    if arr.ndim == 3 and arr.shape[0] in (3, 4):
                        arr = da.transpose(arr, (1, 2, 0))
                    pyramid_rgb.append(arr)
                return pyramid_rgb

            worker = _load_arms_from_sdata()
            worker.returned.connect(
                lambda result: _on_arms_restored(result, _session_arms_data)
                if ctx.dataset_generation == gen else None
            )
            worker.start()
        elif session.get("arms_he_filename"):
            print(f"  Warning: ARMS H&E image not found in sdata cache, skipping ARMS restore")

    return widget, {
        "restore_session": _restore_session,
        "arms_flip_v": arms_flip_v,
        "arms_flip_h": arms_flip_h,
    }
