"""Tab 5: H&E Registration — load, flip, coarse align, landmarks, register."""

from __future__ import annotations
import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import CheckBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from tabs._helpers import make_tab, StatusProxy

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext

from utils.registration import (
    load_he_pyramid, compute_landmark_affine, save_landmarks, load_landmarks,
    extract_tissue_mask_fluorescence, extract_tissue_mask_he, compute_coarse_affine,
)


def build_tab(ctx: ViewerContext) -> tuple:
    he_state = ctx.he_state

    he_load_button = PushButton(label="Load H&E Image...", enabled=True)
    he_flip_v = CheckBox(label="Flip vertically", value=False)
    he_flip_h = CheckBox(label="Flip horizontally", value=False)
    he_opacity_slider = Slider(label="H&E opacity", min=0, max=100, value=70)
    he_opacity_slider.enabled = False

    coarse_align_button = PushButton(label="Coarse Align", enabled=False)

    add_xenium_lm_button = PushButton(label="Add Xenium Landmark", enabled=False)
    add_he_lm_button = PushButton(label="Add H&E Landmark", enabled=False)
    clear_lm_button = PushButton(label="Clear All", enabled=False)

    register_button = PushButton(label="Compute Registration", enabled=False)

    reg_residuals_qt = QTextEdit()
    reg_residuals_qt.setReadOnly(True)
    reg_residuals_qt.setFontFamily("monospace")
    reg_residuals_qt.setMaximumHeight(150)

    save_lm_button = PushButton(label="Save Landmarks...", enabled=False)
    load_lm_button = PushButton(label="Load Landmarks...", enabled=True)
    save_affine_button = PushButton(label="Save Affine...", enabled=False)

    he_status_label = StatusProxy(ctx.viewer)
    reg_status_label = StatusProxy(ctx.viewer)

    viewer = ctx.viewer
    sdata = ctx.sdata
    data_path = ctx.data_path
    no_cache = ctx.no_cache
    pixel_size = ctx.pixel_size
    morph_thumb = ctx.morph_thumb
    morph_full_shape_yx = ctx.morph_full_shape_yx

    def _build_flip_affine():
        shape = he_state.get("he_shape_yx")
        if shape is None:
            return np.eye(3)
        h, w = shape
        M = np.eye(3)
        if he_flip_v.value:
            M = np.array([[-1, 0, h - 1], [0, 1, 0], [0, 0, 1]], dtype=np.float64) @ M
        if he_flip_h.value:
            M = np.array([[1, 0, 0], [0, -1, w - 1], [0, 0, 1]], dtype=np.float64) @ M
        return M

    def _apply_he_affine():
        flip = _build_flip_affine()
        fine = he_state["affine_3x3"]
        coarse = he_state["coarse_affine"]
        if fine is not None:
            combined = fine @ flip
        elif coarse is not None:
            combined = coarse @ flip
        else:
            combined = flip
        if he_state["he_layer"] is not None:
            he_state["he_layer"].affine = combined
        if he_state["he_lm_layer"] is not None:
            he_state["he_lm_layer"].affine = combined

    def on_flip_changed(_value=None):
        _apply_he_affine()
        he_state["flip_v"] = he_flip_v.value
        he_state["flip_h"] = he_flip_h.value
        _save_he_affine_to_sdata()
        flips = []
        if he_flip_v.value: flips.append("V")
        if he_flip_h.value: flips.append("H")
        if flips:
            he_status_label.value = f"Flip applied: {'+'.join(flips)}"

    he_flip_v.changed.connect(on_flip_changed)
    he_flip_h.changed.connect(on_flip_changed)

    def _check_landmark_count(*_args):
        xen = he_state["xenium_lm_layer"]
        he = he_state["he_lm_layer"]
        if xen is not None and he is not None:
            n = min(len(xen.data), len(he.data))
            register_button.enabled = n >= 3
            save_lm_button.enabled = n >= 1

    def _create_landmark_layers():
        if he_state["xenium_lm_layer"] is not None:
            return
        xen_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="Xenium Landmarks", size=30, face_color="cyan",
            symbol="cross", border_color="cyan",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )
        he_lm = viewer.add_points(
            np.empty((0, 2), dtype=np.float64),
            name="H&E Landmarks", size=30, face_color="red",
            symbol="cross", border_color="red",
            border_width=0.1, border_width_is_relative=True, opacity=1.0,
        )
        xen_lm.events.data.connect(_check_landmark_count)
        he_lm.events.data.connect(_check_landmark_count)
        he_state["xenium_lm_layer"] = xen_lm
        he_state["he_lm_layer"] = he_lm
        add_xenium_lm_button.enabled = True
        add_he_lm_button.enabled = True
        clear_lm_button.enabled = True

    def _save_he_to_sdata(pyramid, he_filename):
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
            if "he_image" in sdata.images:
                del sdata.images["he_image"]
            sdata.images["he_image"] = parsed
            sdata.write_element("he_image", overwrite=True)
            zarr_path = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_path), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            store["viewer_session"].attrs["he_filename"] = he_filename
            store["viewer_session"].attrs["he_shape_yx"] = list(base.shape[:2])
            print(f"  H&E image saved to sdata zarr cache ({base_cyx.shape})")
        except Exception as e:
            print(f"  Warning: could not save H&E to sdata: {e}")

    def _save_he_affine_to_sdata():
        if sdata is None or no_cache or "he_image" not in sdata.images:
            return
        try:
            from spatialdata.transformations import Affine as SdAffine, set_transformation
            flip = _build_flip_affine()
            fine = he_state["affine_3x3"]
            coarse = he_state["coarse_affine"]
            if fine is not None:
                combined = fine @ flip
            elif coarse is not None:
                combined = coarse @ flip
            else:
                combined = flip
            sd_affine = SdAffine(combined, input_axes=("y", "x"), output_axes=("y", "x"))
            set_transformation(sdata.images["he_image"], sd_affine, "global")
            sdata.write_transformations("he_image")
            zarr_path = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_path), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            sess = store["viewer_session"]
            sess.attrs["flip_v"] = bool(he_state.get("flip_v", False))
            sess.attrs["flip_h"] = bool(he_state.get("flip_h", False))
            if fine is not None:
                sess.attrs["affine_3x3"] = fine.tolist()
            if coarse is not None:
                sess.attrs["coarse_affine"] = coarse.tolist()
        except Exception as e:
            print(f"  Warning: could not save H&E affine: {e}")

    def _on_he_loaded(result):
        (pyramid, tif), path = result
        if he_state["he_layer"] is not None:
            try:
                viewer.layers.remove(he_state["he_layer"])
            except ValueError:
                pass
        he_state["he_tif"] = tif
        he_state["he_filename"] = Path(path).name
        he_state["he_path"] = str(path)
        base = pyramid[0]
        he_state["he_shape_yx"] = (base.shape[0], base.shape[1])
        he_layer = viewer.add_image(
            pyramid, name=f"H&E ({Path(path).name})",
            rgb=True, blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        he_state["he_layer"] = he_layer
        he_state["affine_3x3"] = None
        he_state["coarse_affine"] = None
        _apply_he_affine()
        _create_landmark_layers()
        _save_he_to_sdata(pyramid, Path(path).name)
        he_opacity_slider.enabled = True
        he_load_button.enabled = True
        coarse_align_button.enabled = morph_thumb is not None
        shape_str = "x".join(str(s) for s in pyramid[0].shape)
        he_status_label.value = f"H&E loaded: {Path(path).name} ({shape_str}, {len(pyramid)} levels)"
        ctx.record_code(
            f"\n# Load H&E image\n"
            f"from utils.registration import load_he_pyramid\n"
            f"he_pyramid, he_tif = load_he_pyramid(\"{path}\")"
        )

    def on_load_he():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load H&E Image", default_dir,
            "Image Files (*.ome.tif *.tif *.tiff *.svs);;All Files (*)",
        )
        if not path:
            return
        he_status_label.value = "Loading H&E..."
        he_load_button.enabled = False
        gen = ctx.dataset_generation

        @thread_worker
        def load_task():
            return load_he_pyramid(path), path

        worker = load_task()
        worker.returned.connect(lambda result: _on_he_loaded(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def on_he_opacity(value):
        if he_state["he_layer"] is not None:
            he_state["he_layer"].opacity = value / 100.0

    def _on_coarse_done(coarse_affine):
        he_state["coarse_affine"] = coarse_affine
        he_state["affine_3x3"] = None
        _apply_he_affine()
        _save_he_affine_to_sdata()
        coarse_align_button.enabled = True
        scale = np.sqrt(coarse_affine[0, 0]**2 + coarse_affine[0, 1]**2)
        reg_status_label.value = f"Coarse aligned (scale={scale:.4f}). Place landmarks to refine."
        ctx.record_code(f"\n# Coarse H&E alignment (tissue outlines)\n# scale={scale:.4f}")
        reg_residuals_qt.setPlainText(
            f"Coarse tissue-outline alignment applied.\n"
            f"Scale: {scale:.4f}\n"
            f"Place >= 3 matching landmarks, then click 'Compute Registration'\n"
            f"to refine alignment."
        )

    def on_coarse_align():
        if he_state["he_layer"] is None:
            reg_status_label.value = "Load H&E image first"
            return
        if morph_thumb is None:
            reg_status_label.value = "No morphology data available"
            return
        reg_status_label.value = "Computing coarse alignment..."
        coarse_align_button.enabled = False
        gen = ctx.dataset_generation

        @thread_worker
        def _compute_coarse():
            he_pyramid = he_state["he_layer"].data
            he_low = np.asarray(he_pyramid[-1])
            morph_mask = extract_tissue_mask_fluorescence(morph_thumb)
            he_mask = extract_tissue_mask_he(he_low)
            target_ds = morph_full_shape_yx[0] / morph_thumb.shape[1]
            he_full_shape = he_state["he_shape_yx"]
            source_ds = he_full_shape[0] / he_low.shape[0]
            coarse_affine_yx = compute_coarse_affine(
                target_mask=morph_mask, source_mask=he_mask,
                target_downsample=target_ds, source_downsample=source_ds,
            )
            return coarse_affine_yx

        worker = _compute_coarse()
        worker.returned.connect(lambda result: _on_coarse_done(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def on_add_xenium_lm():
        lm = he_state["xenium_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.value = "Click on a feature in the Xenium image"

    def on_add_he_lm():
        lm = he_state["he_lm_layer"]
        if lm is not None:
            viewer.layers.selection.active = lm
            lm.mode = "add"
            reg_status_label.value = "Click on the same feature in the H&E image"

    def on_clear_lm():
        for key in ("xenium_lm_layer", "he_lm_layer"):
            lm = he_state[key]
            if lm is not None:
                lm.selected_data = set()
                lm.data = np.empty((0, 2), dtype=np.float64)
        he_state["affine_3x3"] = None
        he_state["coarse_affine"] = None
        _apply_he_affine()
        reg_residuals_qt.clear()
        reg_status_label.value = "Landmarks cleared"
        register_button.enabled = False
        save_lm_button.enabled = False
        save_affine_button.enabled = False

    def on_register():
        xen_pts = he_state["xenium_lm_layer"].data
        he_pts = he_state["he_lm_layer"].data
        n = min(len(xen_pts), len(he_pts))
        if n < 3:
            reg_status_label.value = "Need at least 3 paired landmarks"
            return
        xen_pts = np.asarray(xen_pts[:n], dtype=np.float64)
        he_pts = np.asarray(he_pts[:n], dtype=np.float64)
        affine, residuals = compute_landmark_affine(xen_pts, he_pts)
        he_state["affine_3x3"] = affine
        _apply_he_affine()
        lines = [f"Registration: {n} landmarks, similarity transform"]
        lines.append(f"Mean residual: {residuals.mean():.1f} px ({residuals.mean() * pixel_size:.1f} um)")
        lines.append(f"Max  residual: {residuals.max():.1f} px ({residuals.max() * pixel_size:.1f} um)")
        lines.append("")
        for i, r in enumerate(residuals):
            lines.append(f"  Landmark {i+1}: {r:.1f} px ({r * pixel_size:.1f} um)")
        scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
        lines.append(f"\nScale factor: {scale:.4f}")
        reg_residuals_qt.setPlainText("\n".join(lines))
        reg_status_label.value = f"Registered ({n} landmarks, mean residual {residuals.mean():.1f} px)"
        ctx.record_code(
            f"\n# H&E landmark registration\n"
            f"from utils.registration import compute_landmark_affine\n"
            f"he_xen_pts = np.array({xen_pts.tolist()})\n"
            f"he_he_pts = np.array({he_pts.tolist()})\n"
            f"he_affine, he_residuals = compute_landmark_affine(he_xen_pts, he_he_pts)"
        )
        _save_he_affine_to_sdata()
        save_affine_button.enabled = True

    def on_save_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Landmarks", default_dir + "/landmarks.json", "JSON Files (*.json)",
        )
        if not path:
            return
        xen_pts = np.asarray(he_state["xenium_lm_layer"].data, dtype=np.float64)
        he_pts = np.asarray(he_state["he_lm_layer"].data, dtype=np.float64)
        save_landmarks(
            path, xen_pts, he_pts,
            affine=he_state["affine_3x3"], he_filename=he_state["he_filename"],
        )
        reg_status_label.value = f"Landmarks saved to {Path(path).name}"
        ctx.record_code(f"\n# Save landmarks\n# landmarks -> \"{path}\"")

    def on_load_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getOpenFileName(
            None, "Load Landmarks", default_dir, "JSON Files (*.json);;All Files (*)",
        )
        if not path:
            return
        data = load_landmarks(path)
        _create_landmark_layers()
        he_state["xenium_lm_layer"].data = data["xenium_landmarks_yx"]
        he_state["he_lm_layer"].data = data["he_landmarks_yx"]
        if "affine_3x3_yx" in data:
            affine = data["affine_3x3_yx"]
            he_state["affine_3x3"] = affine
            _apply_he_affine()
            save_affine_button.enabled = True
            scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
            reg_residuals_qt.setPlainText(f"Loaded affine (scale={scale:.4f})")
        if "he_filename" in data:
            he_state["he_filename"] = data["he_filename"]
        n = min(len(data["xenium_landmarks_yx"]), len(data["he_landmarks_yx"]))
        reg_status_label.value = f"Loaded {n} landmarks from {Path(path).name}"
        ctx.record_code(
            f"\n# Load landmarks from file\n"
            f"from utils.registration import load_landmarks\n"
            f"landmarks = load_landmarks(\"{path}\")"
        )

    def on_save_affine():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Affine", default_dir + "/he_affine.json", "JSON Files (*.json)",
        )
        if not path:
            return
        affine = he_state["affine_3x3"]
        if affine is None:
            return
        with open(path, "w") as f:
            json.dump({"affine_3x3_yx": affine.tolist()}, f, indent=2)
        reg_status_label.value = f"Affine saved to {Path(path).name}"
        ctx.record_code(f"\n# Save affine transform\n# affine -> \"{path}\"")

    # Wire events
    he_load_button.clicked.connect(on_load_he)
    he_opacity_slider.changed.connect(on_he_opacity)
    coarse_align_button.clicked.connect(on_coarse_align)
    add_xenium_lm_button.clicked.connect(on_add_xenium_lm)
    add_he_lm_button.clicked.connect(on_add_he_lm)
    clear_lm_button.clicked.connect(on_clear_lm)
    register_button.clicked.connect(on_register)
    save_lm_button.clicked.connect(on_save_landmarks)
    load_lm_button.clicked.connect(on_load_landmarks)
    save_affine_button.clicked.connect(on_save_affine)

    # Rows
    lm_btn_row = QWidget()
    lm_btn_layout = QHBoxLayout()
    lm_btn_layout.setContentsMargins(0, 0, 0, 0)
    lm_btn_layout.addWidget(add_xenium_lm_button.native)
    lm_btn_layout.addWidget(add_he_lm_button.native)
    lm_btn_layout.addWidget(clear_lm_button.native)
    lm_btn_row.setLayout(lm_btn_layout)

    io_btn_row = QWidget()
    io_btn_layout = QHBoxLayout()
    io_btn_layout.setContentsMargins(0, 0, 0, 0)
    io_btn_layout.addWidget(save_lm_button.native)
    io_btn_layout.addWidget(load_lm_button.native)
    io_btn_layout.addWidget(save_affine_button.native)
    io_btn_row.setLayout(io_btn_layout)

    flip_row = QWidget()
    flip_layout = QHBoxLayout()
    flip_layout.setContentsMargins(0, 0, 0, 0)
    flip_layout.addWidget(he_flip_v.native)
    flip_layout.addWidget(he_flip_h.native)
    flip_row.setLayout(flip_layout)

    widget = make_tab(
        he_load_button,
        flip_row,
        he_opacity_slider,
        coarse_align_button,
        lm_btn_row,
        register_button,
        io_btn_row,
    )

    # ── Session restore helpers (referenced by orchestrator) ─────────────
    # We need to expose these for H&E session restore
    def _on_he_restored_from_sdata(pyramid_rgb, session_he_data):
        if he_state["he_layer"] is not None:
            try:
                viewer.layers.remove(he_state["he_layer"])
            except ValueError:
                pass
        he_filename = session_he_data.get("he_filename", "H&E")
        he_state["he_tif"] = None
        he_state["he_filename"] = he_filename
        he_state["he_path"] = None
        base = pyramid_rgb[0]
        he_state["he_shape_yx"] = session_he_data.get("he_shape_yx") or (base.shape[0], base.shape[1])
        he_layer = viewer.add_image(
            pyramid_rgb, name=f"H&E ({he_filename})",
            rgb=True, blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        he_state["he_layer"] = he_layer
        he_state["affine_3x3"] = session_he_data.get("affine_3x3")
        he_state["coarse_affine"] = session_he_data.get("coarse_affine")
        _apply_he_affine()
        _create_landmark_layers()
        xen_lm = session_he_data.get("xenium_landmarks")
        he_lm_data = session_he_data.get("he_landmarks")
        if xen_lm is not None and he_state["xenium_lm_layer"] is not None:
            he_state["xenium_lm_layer"].data = xen_lm
        if he_lm_data is not None and he_state["he_lm_layer"] is not None:
            he_state["he_lm_layer"].data = he_lm_data
        he_opacity_slider.enabled = True
        he_load_button.enabled = True
        coarse_align_button.enabled = morph_thumb is not None
        save_affine_button.enabled = he_state["affine_3x3"] is not None
        has_affine = he_state["affine_3x3"] is not None or he_state["coarse_affine"] is not None
        he_status_label.value = f"H&E restored: {he_filename}" + (" (with registration)" if has_affine else "")
        print(f"  Restored H&E from cache: {he_filename}" + (" with registration" if has_affine else ""))

    def _restore_session(session):
        from tabs._helpers import StatusProxy as _SP  # noqa: avoid circular

        if sdata is not None and "he_image" in sdata.images:
            he_flip_v.value = session.get("flip_v", False)
            he_flip_h.value = session.get("flip_h", False)
            _session_he_data = {
                "affine_3x3": session.get("affine_3x3"),
                "coarse_affine": session.get("coarse_affine"),
                "xenium_landmarks": session.get("xenium_landmarks"),
                "he_landmarks": session.get("he_landmarks"),
                "he_filename": session.get("he_filename", "H&E"),
                "he_shape_yx": session.get("he_shape_yx"),
            }
            he_status_label.value = "Restoring H&E from cache..."
            gen = ctx.dataset_generation

            # Import _extract_dt_scales at runtime from the main module
            from tabs.tab_he_registration import _extract_dt_scales

            @thread_worker
            def _load_he_from_sdata():
                he_dt = sdata.images["he_image"]
                pyramid = _extract_dt_scales(he_dt)
                pyramid_rgb = []
                for arr in pyramid:
                    computed = arr.compute() if hasattr(arr, 'compute') else np.asarray(arr)
                    if computed.ndim == 3 and computed.shape[0] in (3, 4):
                        computed = np.transpose(computed, (1, 2, 0))
                    pyramid_rgb.append(computed)
                return pyramid_rgb

            worker = _load_he_from_sdata()
            worker.returned.connect(
                lambda result: _on_he_restored_from_sdata(result, _session_he_data)
                if ctx.dataset_generation == gen else None
            )
            worker.start()
        elif session.get("he_filename"):
            print(f"  Warning: H&E image not found in sdata cache, skipping H&E restore")

    return widget, {
        "restore_session": _restore_session,
        "he_flip_v": he_flip_v,
        "he_flip_h": he_flip_h,
        "he_opacity_slider": he_opacity_slider,
        "save_affine_button": save_affine_button,
        "coarse_align_button": coarse_align_button,
        "create_landmark_layers": _create_landmark_layers,
        "apply_he_affine": _apply_he_affine,
    }


def _extract_dt_scales(dt):
    """Extract an ordered list of dask arrays from a spatialdata DataTree."""
    import re

    def _sort_key(name):
        nums = re.findall(r'\d+', name)
        return int(nums[0]) if nums else 0

    scales = []
    for name in sorted(dt.children.keys(), key=_sort_key):
        child = dt.children[name]
        ds = getattr(child, 'ds', None)
        if ds is None:
            continue
        if 'image' in ds:
            scales.append(ds['image'].data)
        elif ds.data_vars:
            first = next(iter(ds.data_vars))
            scales.append(ds[first].data)
    return scales
