"""Tab 5: H&E Registration — load, flip, coarse align, landmarks, register."""

from __future__ import annotations
import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import CheckBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from xenium_viewer.utils.units import px_affine_to_world
from xenium_viewer.tabs._helpers import make_tab, StatusProxy
from xenium_viewer.utils.prov_graph import TERMINAL
from xenium_viewer.utils.zarr_safe import safe_write_element

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

from xenium_viewer.utils.registration import (
    load_he_pyramid, compute_landmark_affine, save_landmarks, load_landmarks,
    extract_tissue_mask_fluorescence, extract_tissue_mask_he, compute_coarse_affine,
    describe_pyramid, parse_rgb_image_for_store,
)
from xenium_viewer.utils.reporting import report_write_failure


def build_tab(ctx: ViewerContext) -> tuple:
    he_state = ctx.he_state

    he_load_button = PushButton(label="Load H&E Image...", enabled=True)
    he_flip_v = CheckBox(label="Flip vertically", value=False)
    he_flip_h = CheckBox(label="Flip horizontally", value=False)
    he_opacity_slider = Slider(label="Opacity", min=0, max=100, value=70)
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
        # `combined` is in Xenium pixels — that is what the store and the crop
        # export expect. napari applies a layer's affine *after* its scale, so
        # the layer's copy has to be in world (µm) units. See utils/units.py.
        world = px_affine_to_world(combined, ctx.pixel_size)
        if he_state["he_layer"] is not None:
            he_state["he_layer"].affine = world
        if he_state["he_lm_layer"] is not None:
            he_state["he_lm_layer"].affine = world

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
        ctx.record_node(
            "he:flip",
            f"\n# H&E image flip, applied before registration\n"
            f"he_flip_vertical = {he_flip_v.value}\n"
            f"he_flip_horizontal = {he_flip_h.value}",
            deps=["preamble"], kind=TERMINAL, label="H&E flip",
        )

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
            parsed, shape_yx = parse_rgb_image_for_store(pyramid[0])
            # replace_backed: _on_he_loaded removed the old layer before calling
            # us, so the only thing still holding the stored element is the
            # sdata binding this write replaces. Without it, loading a second
            # H&E over one restored from the cache was refused and silently
            # never persisted.
            safe_write_element(sdata, "he_image", parsed, replace_backed=True)
            zarr_path = data_path / "sdata_cached.zarr"
            import zarr as zarr_mod
            store = zarr_mod.open_group(str(zarr_path), mode="r+", use_consolidated=False)
            if "viewer_session" not in store:
                store.create_group("viewer_session")
            store["viewer_session"].attrs["he_filename"] = he_filename
            store["viewer_session"].attrs["he_shape_yx"] = list(shape_yx)
            print(f"  H&E image saved to sdata zarr cache ({shape_yx[0]}x{shape_yx[1]})")
        except Exception as e:
            report_write_failure(e, "H&E image")

    def _save_he_affine_to_sdata():
        if sdata is None or no_cache or "he_image" not in sdata.images:
            return
        try:
            from spatialdata.transformations import Affine as SdAffine, set_transformation
            from xenium_viewer.utils.adata_persistence import _load_affine_from_sdata_element
            flip = _build_flip_affine()
            fine = he_state["affine_3x3"]
            coarse = he_state["coarse_affine"]
            if fine is not None:
                combined = fine @ flip
            elif coarse is not None:
                combined = coarse @ flip
            else:
                combined = flip
            # See the twin guard in tab_arms._save_arms_affine_to_sdata: an
            # identity built from an empty he_state must never overwrite a
            # registration the element already carries.
            if fine is None and coarse is None and np.allclose(combined, np.eye(3), atol=1e-6):
                if _load_affine_from_sdata_element(sdata, "he_image") is not None:
                    return
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
        print(f"  {describe_pyramid(pyramid, f'H&E {Path(path).name}')}")
        he_status_label.value = f"H&E loaded: {Path(path).name} ({shape_str}, {len(pyramid)} levels)"
        ctx.record_node(
            "he:load",
            f"\n# Load H&E image\n"
            f"from xenium_viewer.utils.registration import load_he_pyramid\n"
            f"he_pyramid, he_tif = load_he_pyramid(\"{path}\")",
            deps=["preamble"], kind=TERMINAL, label="Load H&E image",
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
        ctx.record_node(
            "he:coarse_align",
            f"\n# Coarse H&E alignment to tissue outlines (scale={scale:.4f}).\n"
            f"# Computed from the morphology thumbnail; the matrix it produced:\n"
            f"he_coarse_affine = np.array({np.asarray(coarse_affine).tolist()})",
            deps=["preamble"], kind=TERMINAL, label="H&E coarse align",
        )
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
        from xenium_viewer.utils.adata_persistence import save_landmarks_to_sdata
        save_landmarks_to_sdata(ctx, 'he_xenium_landmarks', None)
        save_landmarks_to_sdata(ctx, 'he_he_landmarks', None)

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
        ctx.record_node(
            "he:landmark_register",
            f"\n# H&E landmark registration\n"
            f"from xenium_viewer.utils.registration import compute_landmark_affine\n"
            f"he_xen_pts = np.array({xen_pts.tolist()})\n"
            f"he_he_pts = np.array({he_pts.tolist()})\n"
            f"he_affine, he_residuals = compute_landmark_affine(he_xen_pts, he_he_pts)",
            deps=["preamble"], kind=TERMINAL, label="H&E landmark registration",
        )
        _save_he_affine_to_sdata()
        from xenium_viewer.utils.adata_persistence import save_landmarks_to_sdata
        save_landmarks_to_sdata(ctx, 'he_xenium_landmarks', np.asarray(xen_pts))
        save_landmarks_to_sdata(ctx, 'he_he_landmarks', np.asarray(he_pts))

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
        # The points are inlined rather than referenced: landmarks can be saved
        # before a registration is computed, so ``he_xen_pts`` may not exist.
        affine = he_state["affine_3x3"]
        ctx.record_node(
            "he:save_landmarks",
            f"\n# Save H&E landmarks to {Path(path).name}\n"
            f"from xenium_viewer.utils.registration import save_landmarks\n"
            f"save_landmarks(\n"
            f"    r\"{path}\",\n"
            f"    np.array({xen_pts.tolist()}),\n"
            f"    np.array({he_pts.tolist()}),\n"
            f"    affine={None if affine is None else f'np.array({np.asarray(affine).tolist()})'},\n"
            f"    he_filename={he_state['he_filename']!r},\n"
            f")",
            deps=["preamble"], kind=TERMINAL, label="Save H&E landmarks",
        )

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
            scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
            reg_residuals_qt.setPlainText(f"Loaded affine (scale={scale:.4f})")
        if "he_filename" in data:
            he_state["he_filename"] = data["he_filename"]
        n = min(len(data["xenium_landmarks_yx"]), len(data["he_landmarks_yx"]))
        reg_status_label.value = f"Loaded {n} landmarks from {Path(path).name}"
        ctx.record_node(
            "he:load_landmarks",
            f"\n# Load H&E landmarks from file\n"
            f"from xenium_viewer.utils.registration import load_landmarks\n"
            f"landmarks = load_landmarks(\"{path}\")",
            deps=["preamble"], kind=TERMINAL, label="Load H&E landmarks",
        )

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
        # Filled at :361 but never laid out. Tab-HE-Registration.md has always
        # documented a "Residuals (read-only)" control; this is what makes that
        # true. Per-landmark residuals are the only way to tell which landmark
        # is dragging the fit, and the status bar shows the mean alone.
        reg_residuals_qt,
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
        # Resolve the placement before the flip widgets are read: on a store with
        # no viewer_session the element's own transform is the only record of the
        # registration, and it already contains the flip. See utils/registration_seed.
        from xenium_viewer.utils.registration_seed import seed_registration
        reg = seed_registration(
            session_he_data, sdata, "he_image",
            affine_key="affine_3x3", coarse_key="coarse_affine",
            flip_v_key="flip_v", flip_h_key="flip_h", shape_key="he_shape_yx",
            element_shape_yx=(base.shape[0], base.shape[1]),
        )
        he_state["he_shape_yx"] = reg.shape_yx or (base.shape[0], base.shape[1])
        # Blocked: on_flip_changed writes the element transform and records a
        # `he:flip` provenance node. Restoring a value the user chose earlier is
        # not the user choosing it again.
        with he_flip_v.changed.blocked(), he_flip_h.changed.blocked():
            he_flip_v.value = reg.flip_v
            he_flip_h.value = reg.flip_h
        he_layer = viewer.add_image(
            pyramid_rgb, name=f"H&E ({he_filename})",
            rgb=True, blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        he_state["he_layer"] = he_layer
        he_state["affine_3x3"] = reg.fine
        he_state["coarse_affine"] = reg.coarse
        he_state["flip_v"] = reg.flip_v
        he_state["flip_h"] = reg.flip_h
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
        has_affine = he_state["affine_3x3"] is not None or he_state["coarse_affine"] is not None
        how = ""
        if has_affine:
            how = (" (registration from the image itself)" if reg.source == "element"
                   else " (with registration)")
        he_status_label.value = f"H&E restored: {he_filename}{how}"
        print(f"  Restored H&E from cache: {he_filename}"
              + (f" with registration ({reg.source})" if has_affine else ""))

    def _restore_session(session):
        from xenium_viewer.tabs._helpers import StatusProxy as _SP  # noqa: avoid circular

        if sdata is not None and "he_image" in sdata.images:
            # The flips are set in _on_he_restored_from_sdata, from the same
            # record the affine comes from — setting them here too would pair a
            # session flip with an element affine that already contains one.
            _session_he_data = {
                "affine_3x3": session.get("affine_3x3"),
                "coarse_affine": session.get("coarse_affine"),
                "xenium_landmarks": session.get("xenium_landmarks"),
                "he_landmarks": session.get("he_landmarks"),
                "he_filename": session.get("he_filename", "H&E"),
                "he_shape_yx": session.get("he_shape_yx"),
                "flip_v": session.get("flip_v", False),
                "flip_h": session.get("flip_h", False),
            }
            he_status_label.value = "Restoring H&E from cache..."
            gen = ctx.dataset_generation

            # Import _extract_dt_scales at runtime from the main module
            from xenium_viewer.tabs.tab_he_registration import _extract_dt_scales

            @thread_worker
            def _load_he_from_sdata():
                # Lazy on purpose: the eager version computed every level,
                # scale0 included, into dense numpy on every launch. The ARMS
                # copy of this was fixed in 9cad210 and this one was missed.
                # napari fetches only the tiles it draws from a dask multiscale.
                import dask.array as da
                he_dt = sdata.images["he_image"]
                pyramid = _extract_dt_scales(he_dt)
                pyramid_rgb = []
                for arr in pyramid:
                    if not isinstance(arr, da.Array):
                        arr = da.from_array(arr)
                    if arr.ndim == 3 and arr.shape[0] in (3, 4):
                        arr = da.transpose(arr, (1, 2, 0))
                    pyramid_rgb.append(arr)
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
