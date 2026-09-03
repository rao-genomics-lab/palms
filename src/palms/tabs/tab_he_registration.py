"""Tab 5: H&E Registration — load, flip, coarse align, landmarks, register.

**The frame rule.** ``coarse`` and ``fine`` are both maps from the *flipped* H&E
frame to Xenium pixels. That is what ``_apply_he_affine`` composes (``X @ flip``)
and what docs/Tab-HE-Registration.md describes — but neither fit site honoured it
until 2026-09-02: the coarse mask came off the unflipped pyramid and the landmark
points came off ``he_lm_layer.data``, which is unflipped layer-data coordinates.
Any registration made with a Flip ticked was therefore flipped twice. It matters
beyond that latent bug, because automatic mirror detection routes *through* the
flip: ``compute_landmark_affine`` fits a similarity, which has no reflection, so
a mirror baked into ``coarse`` would be silently discarded the moment the user
pressed Compute Registration.
"""

from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import CheckBox, PushButton, Slider
from qtpy.QtWidgets import QTextEdit, QHBoxLayout, QWidget, QFileDialog
from napari.qt.threading import thread_worker
from palms.utils.units import px_affine_to_world
from palms.tabs._helpers import make_tab, StatusProxy
from palms.utils.prov_graph import ARTIFACT, TERMINAL
from palms.utils.step_templates import Preview, step_template
from palms.utils.steps import Step, coerce
from palms.utils.zarr_safe import safe_write_element

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext

from palms.utils.registration import (
    save_landmarks, load_landmarks,
    describe_pyramid, parse_rgb_image_for_store, flip_matrix, pyramid_levels,
)
from palms.utils.reporting import report_write_failure

# The five steps this tab runs. Named once so the Templates tab, the preview
# providers and the run sites cannot disagree about which template is which.
LOAD_TEMPLATE = "he.load"
FLIP_TEMPLATE = "he.flip"
COARSE_TEMPLATE = "he.coarse_align"
NUCLEI_TEMPLATE = "he.nuclei_align"
LANDMARK_TEMPLATE = "he.landmark_align"


def build_tab(ctx: ViewerContext) -> tuple:
    he_state = ctx.he_state

    he_load_button = PushButton(label="Load H&E Image...", enabled=True)
    he_flip_v = CheckBox(label="Flip vertically", value=False)
    he_flip_h = CheckBox(label="Flip horizontally", value=False)
    he_opacity_slider = Slider(label="Opacity", min=0, max=100, value=70)
    he_opacity_slider.enabled = False

    coarse_align_button = PushButton(label="Coarse Align", enabled=False)
    nuclei_align_button = PushButton(label="Fine Align (nuclei)", enabled=False)

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
    # Only to decide whether Coarse Align can be offered at all. The search
    # derives its own thumbnail inside the template, from `sdata`, so that the
    # notebook's cell reads the same element rather than a value the viewer
    # happened to compute at launch.
    morph_thumb = getattr(ctx, "morph_thumb", None)

    def _step_progress(prefix: str, status):
        """Adapt ``StepExecutor``'s per-statement callback to the status bar.

        Both fits are long -- the nuclei one runs for ten minutes on a full
        slide over a network share -- so silence is not an option. The reporting
        rides on ``run(progress=)``, which fires before each top-level statement
        of the very source being recorded, so the readout cannot claim a stage
        the cell does not contain; a hand-maintained list of stage names beside
        the template could, and the previous generator-based version did exactly
        that in a second copy of the orchestration.

        Bounced to the GUI thread: steps run in a napari worker and a status
        write is a Qt write.
        """
        from superqt.utils import ensure_main_thread

        @ensure_main_thread
        def _show(text):
            status.value = text

        def _report(index, total, label):
            _show(f"{prefix} ({index}/{total}): {label}")

        return _report

    def _flip_step() -> Step:
        blocks, params, _ = _flip_preview()
        return Step(
            id="he:flip", **step_template(FLIP_TEMPLATE, blocks),
            params=params,
            deps=["preamble"], kind=ARTIFACT, label="H&E flip",
            outputs=["he_flip_vertical", "he_flip_horizontal"],
        )

    def _record_flip():
        """Run and record the flip declaration both fits consume.

        An ``ARTIFACT`` rather than a terminal, and that is the point: the two
        fits ``deps`` on it, so changing a flip checkbox marks them stale --
        which is what the GUI has always *done* (``on_flip_changed`` drops the
        coarse affine) without the graph ever saying so.
        """
        ctx.run_step(_flip_step())

    def _build_flip_affine():
        shape = he_state.get("he_shape_yx")
        if shape is None:
            return np.eye(3)
        return flip_matrix(shape, he_flip_v.value, he_flip_h.value)

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

    def _flip_points(pts_yx):
        """H&E landmark data (unflipped layer coordinates) in the flipped frame.

        See the frame rule in the module docstring: the fit has to happen in the
        frame the affine is composed with, or the flip is applied twice.
        """
        pts = np.asarray(pts_yx, dtype=np.float64)
        if len(pts) == 0:
            return pts
        flip = _build_flip_affine()
        homo = np.hstack([pts, np.ones((len(pts), 1))])
        return (flip @ homo.T).T[:, :2]

    def on_flip_changed(_value=None):
        he_state["flip_v"] = he_flip_v.value
        he_state["flip_h"] = he_flip_h.value
        # A coarse affine was fitted for the orientation that has just changed,
        # so it no longer describes this image. Cleared *before* the layer is
        # re-placed, or the overlay jumps to a transform that is already void.
        dropped_coarse = he_state.get("coarse_affine") is not None
        if dropped_coarse:
            he_state["coarse_affine"] = None
        _apply_he_affine()
        _save_he_affine_to_sdata()
        flips = []
        if he_flip_v.value: flips.append("V")
        if he_flip_h.value: flips.append("H")
        parts = ([f"Flip applied: {'+'.join(flips)}"] if flips else [])
        if dropped_coarse:
            parts.append("coarse alignment cleared — run Coarse Align again")
        if parts:
            he_status_label.value = " — ".join(parts)
        _record_flip()

    he_flip_v.changed.connect(on_flip_changed)
    he_flip_h.changed.connect(on_flip_changed)

    def _check_landmark_count(*_args):
        xen = he_state.get("xenium_lm_layer")
        he = he_state.get("he_lm_layer")
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
            from palms.utils.adata_persistence import _load_affine_from_sdata_element
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
            if he_state.get("he_pixel_size_um"):
                sess.attrs["he_pixel_size_um"] = float(he_state["he_pixel_size_um"])
            if fine is not None:
                sess.attrs["affine_3x3"] = fine.tolist()
                if he_state.get("affine_source"):
                    sess.attrs["affine_source"] = he_state["affine_source"]
            if coarse is not None:
                sess.attrs["coarse_affine"] = coarse.tolist()
        except Exception as e:
            print(f"  Warning: could not save H&E affine: {e}")

    def _he_load_step(path: str) -> Step:
        """The step that binds ``he_pyramid`` and ``he_px_um``.

        ``from_file`` whenever the original image is still on this machine,
        which is the only variant a notebook replayed against the raw Xenium
        output can run; ``from_store`` is the fallback for a session whose H&E
        survives only in the viewer's cache. The block is chosen here rather
        than in the template because which of the two is available is runtime
        state, exactly like every other block selection in this codebase.
        """
        blocks, params, _ = _he_load_preview(path)
        return Step(
            id="he:load", **step_template(LOAD_TEMPLATE, blocks),
            params=params,
            deps=["preamble"], kind=ARTIFACT, label="Load H&E image",
            outputs=["he_pyramid", "he_px_um"],
        )

    def _he_load_preview(path=None) -> Preview:
        """Which image the load step would read, and how it would reach it.

        *path* is passed by the run site at the moment the user picks a file --
        before ``he_state`` knows about it -- and defaults to whatever is on
        record, which is what the Templates pane asks for.
        """
        path = path or he_state.get("he_path")
        return Preview(
            blocks=["from_file"] if path else ["from_store"],
            params={"path": str(path) if path else None,
                    "px_um": coerce(he_state.get("he_pixel_size_um"))},
            note=("" if path else
                  "no file path is on record for this H&E, so the cell reads "
                  "the copy in the viewer's cache"),
        )

    def _flip_preview() -> Preview:
        return Preview(blocks=["main"],
                       params={"flip_v": bool(he_flip_v.value),
                               "flip_h": bool(he_flip_h.value)})

    def _ensure_he_loaded():
        """Bind the H&E in the executor namespace, running ``he:load`` if needed.

        Beside ``ctx.ensure_normalized`` and ``ctx.ensure_spatial_neighbors``,
        and for the same reason: a fit that consumes ``he_pyramid`` must be able
        to say what produced it. A session restored from the cache has an H&E on
        the canvas that no step in *this* session ever bound, so without this the
        first Coarse Align after a restore would fail on a NameError from a
        template that is perfectly correct.
        """
        executor = ctx.executor
        if executor is not None and "he_pyramid" in executor.names():
            return
        ctx.run_step(_he_load_step(he_state.get("he_path")))

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
        # Read by the recorded step, not here: it is the step's `he_px_um`
        # output. Stored in the session because a cache-restored H&E has no
        # TiffFile to read it back from.
        he_state["he_pixel_size_um"] = ctx.executor.get("he_px_um")
        he_layer = viewer.add_image(
            pyramid, name=f"H&E ({Path(path).name})",
            rgb=True, blending="translucent",
            opacity=he_opacity_slider.value / 100.0,
        )
        he_state["he_layer"] = he_layer
        he_state["affine_3x3"] = None
        he_state["coarse_affine"] = None
        _apply_he_affine()
        _refresh_nuclei_button()
        _create_landmark_layers()
        _save_he_to_sdata(pyramid, Path(path).name)
        he_opacity_slider.enabled = True
        he_load_button.enabled = True
        coarse_align_button.enabled = morph_thumb is not None
        shape_str = "x".join(str(s) for s in pyramid[0].shape)
        print(f"  {describe_pyramid(pyramid, f'H&E {Path(path).name}')}")
        he_status_label.value = f"H&E loaded: {Path(path).name} ({shape_str}, {len(pyramid)} levels)"
        # No recording here: `load_task` ran the step, which recorded itself.

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
            # Run the recorded step rather than calling the loader beside it:
            # the pyramid the viewer displays is then the object the notebook's
            # cell produces, not a second read of the same file that happens to
            # agree. `he_tif` comes back out of the namespace because the tab
            # needs the open handle to write the image into the store.
            out = ctx.run_step(_he_load_step(path))
            return (out["he_pyramid"], ctx.executor.get("he_tif")), path

        worker = load_task()
        worker.returned.connect(lambda result: _on_he_loaded(result) if ctx.dataset_generation == gen else None)
        worker.start()

    def on_he_opacity(value):
        if he_state["he_layer"] is not None:
            he_state["he_layer"].opacity = value / 100.0

    def _on_coarse_done(result):
        coarse_align_button.enabled = True
        if result is None:
            return
        if result.mirrored:
            # The reflection is carried as an explicit flip, not inside the
            # matrix, so that the landmark refinement inherits it — a similarity
            # fit cannot reproduce a reflection.
            #
            # *Toggling* rather than setting: the search ran on the image as
            # already flipped, so what it found is one further reflection on top
            # of that. Composed with the flips in force it is exactly this — H on
            # nothing is H, H on H is nothing, H on V is both, H on both is V.
            #
            # Blocked because this is the tool choosing the orientation, not the
            # user: on_flip_changed would clear the very coarse affine that is
            # about to be applied.
            with he_flip_v.changed.blocked(), he_flip_h.changed.blocked():
                he_flip_h.value = not he_flip_h.value
            he_state["flip_h"] = he_flip_h.value
            _record_flip()
        he_state["coarse_affine"] = result.affine_3x3_yx
        he_state["affine_3x3"] = None
        _apply_he_affine()
        _save_he_affine_to_sdata()
        _refresh_nuclei_button()
        confidence = "" if result.confident else " — LOW CONFIDENCE, check the overlay"
        reg_status_label.value = (
            f"Coarse aligned (scale={result.scale:.4f}, "
            f"{result.rotation_deg:.1f} deg{', mirrored' if result.mirrored else ''}, "
            f"match {result.score:.2f}){confidence}. Place landmarks to refine."
        )
        # No recording here: `_compute_coarse` ran the step, which recorded the
        # search it performed rather than the matrix that came out of it.
        reg_residuals_qt.setPlainText(
            "Coarse alignment applied.\n"
            + result.summary()
            + "\n\nPlace >= 3 matching landmarks, then click 'Compute Registration'\n"
              "to refine alignment."
        )

    def _coarse_preview() -> Preview:
        """What "Coarse Align" would run as things stand.

        The run site builds its Step from this, so the pane in Tools ->
        Templates and the button are one expression of the settings rather than
        two. There is only one param -- everything else the search reads is
        derived inside the template from `sdata` and the steps it depends on --
        but going through the provider is the property, not the saving.
        """
        return Preview(blocks=["main"],
                       params={"pixel_size": coerce(ctx.pixel_size)})

    def _coarse_step() -> Step:
        blocks, params, _ = _coarse_preview()
        return Step(
            id="he:coarse_align", **step_template(COARSE_TEMPLATE, blocks),
            params=params,
            # Both, and not "preamble": the search reads the image `he:load`
            # binds and works in the frame `he:flip` declares, so re-loading a
            # different H&E or re-ticking a flip has to mark this stale.
            deps=["he:load", "he:flip"], kind=ARTIFACT,
            label="H&E coarse align",
            outputs=["he_coarse_fit", "he_coarse_affine"],
        )

    def _on_coarse_failed(message):
        coarse_align_button.enabled = True
        reg_status_label.value = f"Coarse align failed: {message}"
        reg_residuals_qt.setPlainText(f"Coarse alignment could not run.\n{message}")

    def on_coarse_align():
        if he_state["he_layer"] is None:
            reg_status_label.value = "Load H&E image first"
            return
        if morph_thumb is None:
            reg_status_label.value = "No morphology data available"
            return
        reg_status_label.value = "Computing coarse alignment (global rotation search)..."
        coarse_align_button.enabled = False
        gen = ctx.dataset_generation

        @thread_worker
        def _compute_coarse():
            # Everything the fit reads -- the morphology thumbnail, the H&E
            # level, the flip, the scale prior -- is derived inside the template
            # from `sdata` and the two steps this one depends on. The tab's only
            # contribution is the dataset's pixel size, as a param.
            _ensure_he_loaded()
            _record_flip()
            out = ctx.run_step(_coarse_step(), progress=_step_progress(
                "Coarse align", reg_status_label))
            return out["he_coarse_fit"]

        worker = _compute_coarse()
        worker.returned.connect(
            lambda result: _on_coarse_done(result) if ctx.dataset_generation == gen else None)
        worker.errored.connect(
            lambda exc: _on_coarse_failed(str(exc)) if ctx.dataset_generation == gen else None)
        worker.start()

    def _nuclei_labels_element():
        """The finest level of ``nucleus_labels``, or None if the store has none.

        The nuclei fit is refused rather than approximated when it is missing: a
        cell-boundary raster is not a nuclear one, and matching haematoxylin
        peaks against cell centroids would return a confident transform that is
        wrong by however far a nucleus sits from its cell's centre.
        """
        if sdata is None or "nucleus_labels" not in getattr(sdata, "labels", {}):
            return None
        levels = pyramid_levels(sdata.labels["nucleus_labels"])
        return levels[0] if levels else None

    def _refresh_nuclei_button():
        seed = he_state["affine_3x3"] if he_state["affine_3x3"] is not None \
            else he_state["coarse_affine"]
        nuclei_align_button.enabled = bool(
            he_state["he_layer"] is not None and seed is not None
            and _nuclei_labels_element() is not None)

    def _on_nuclei_done(result):
        nuclei_align_button.enabled = True
        if result is None:
            return
        he_state["affine_3x3"] = result.affine_3x3_yx
        he_state["affine_source"] = "nuclei"
        _apply_he_affine()
        _save_he_affine_to_sdata()
        confidence = "" if result.confident else " — LOW CONFIDENCE, check the overlay"
        reg_status_label.value = (
            f"Fine aligned on {result.n_matched:,} nuclei "
            f"(median residual {result.median_residual_um:.2f} um, "
            f"moved {result.seed_shift_um:.1f} um){confidence}"
        )
        reg_residuals_qt.setPlainText(
            "Automatic fine registration from nuclei.\n" + result.summary()
            + "\n\nPlace landmarks and press 'Compute Registration' only if this\n"
              "needs overriding — it replaces this fit."
        )
        # No recording here: `_compute_nuclei` ran the step, which recorded the
        # detection and the fit rather than the matrix they produced.

    def _on_nuclei_failed(message):
        nuclei_align_button.enabled = True
        reg_status_label.value = f"Fine align failed: {message}"
        reg_residuals_qt.setPlainText(f"Automatic fine registration could not run.\n{message}")

    def _has_step(node_id: str, binding: str) -> bool:
        """Is *node_id* a step this session can actually name as a dependency?

        Both halves are needed and they fail apart. A node restored from the
        session says a step *ran*, in some earlier session, but binds nothing in
        this one -- a template naming its output would raise ``NameError``. A
        binding with no node is the mirror case, and ``upsert`` rejects a
        dependency on a node that does not exist.
        """
        graph = ctx.state.get("prov_graph")
        executor = getattr(ctx, "executor", None)
        return bool(graph is not None and graph.get(node_id) is not None
                    and executor is not None and binding in executor.names())

    def _nuclei_seed():
        """The seed to refine: ``(blocks, extra_params, dependency, note)``.

        Which transform is on hand is runtime state, so the choice is made here
        rather than in the template -- and it *is* the dependency edge, which is
        why the three come back together: a fit refined from the landmark
        transform depends on the landmark step, one refined from the coarse
        search depends on that, and one refined from a transform this session
        merely restored depends on neither and has to inline it.

        Note what is deliberately absent: seeding a nuclei fit from a *previous
        nuclei fit*. That is what the GUI used to do, and it cannot be recorded
        -- the step would depend on itself, which is a cycle, and no notebook can
        express "run this again on its own output". Re-running therefore starts
        from the same place the first run did, which is also what makes the
        recorded step the one that replays.
        """
        fine = he_state.get("affine_3x3")
        if (fine is not None and he_state.get("affine_source") == "landmarks"
                and _has_step("he:landmark_register", "he_affine")):
            return ["seed_fine", "main"], {}, "he:landmark_register", ""
        if _has_step("he:coarse_align", "he_coarse_affine"):
            return ["seed_coarse", "main"], {}, "he:coarse_align", ""
        seed = fine if fine is not None else he_state.get("coarse_affine")
        if seed is None:
            # Only the preview reaches this: the button refuses to run without a
            # starting transform. Shown as the identity so the pane can still
            # say what the cell would look like, with the header saying which
            # value is not settled -- the alternative is a blank pane whenever
            # the tab is opened before Coarse Align has been pressed.
            return (["seed_restored", "main"], {"seed": np.eye(3).tolist()},
                    None, "no alignment yet; the seed shown is the identity")
        return (["seed_restored", "main"],
                {"seed": coerce(np.asarray(seed, dtype=float))}, None, "")

    def _nuclei_preview() -> Preview:
        blocks, extra, _, note = _nuclei_seed()
        return Preview(blocks=blocks,
                       params={"pixel_size": coerce(ctx.pixel_size), **extra},
                       note=note)

    def _nuclei_step() -> Step:
        blocks, params, _ = _nuclei_preview()
        _, _, seed_dep, _ = _nuclei_seed()
        deps = ["he:load", "he:flip"] + ([seed_dep] if seed_dep else [])
        return Step(
            id="he:nuclei_register", **step_template(NUCLEI_TEMPLATE, blocks),
            params=params,
            deps=deps, kind=ARTIFACT,
            label="H&E fine align (nuclei)",
            outputs=["he_nuclei_fit", "he_affine"],
        )

    def on_nuclei_align():
        if _nuclei_labels_element() is None:
            reg_status_label.value = "No nucleus_labels in this dataset"
            return
        seed = he_state["affine_3x3"] if he_state["affine_3x3"] is not None \
            else he_state["coarse_affine"]
        if seed is None:
            reg_status_label.value = "Run Coarse Align first — the nuclei fit needs a starting transform"
            return
        nuclei_align_button.enabled = False
        reg_status_label.value = "Fine align: matching nuclei to nuclear masks..."
        gen = ctx.dataset_generation

        @thread_worker
        def _compute_nuclei():
            _ensure_he_loaded()
            _record_flip()
            out = ctx.run_step(_nuclei_step(), progress=_step_progress(
                "Fine align", reg_status_label))
            return out["he_nuclei_fit"]

        worker = _compute_nuclei()
        worker.returned.connect(
            lambda result: _on_nuclei_done(result) if ctx.dataset_generation == gen else None)
        worker.errored.connect(
            lambda exc: _on_nuclei_failed(str(exc)) if ctx.dataset_generation == gen else None)
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
        _refresh_nuclei_button()
        reg_residuals_qt.clear()
        reg_status_label.value = "Landmarks cleared"
        register_button.enabled = False
        save_lm_button.enabled = False
        from palms.utils.adata_persistence import save_landmarks_to_sdata
        save_landmarks_to_sdata(ctx, 'he_xenium_landmarks', None)
        save_landmarks_to_sdata(ctx, 'he_he_landmarks', None)

    def _landmark_step() -> Step:
        blocks, params, _ = _landmark_preview()
        return Step(
            id="he:landmark_register",
            **step_template(LANDMARK_TEMPLATE, blocks),
            params=params,
            # The landmarks are inlined, so the fit needs nothing the preamble
            # does not already provide -- but it is still an ARTIFACT, because
            # the nuclei fit can be seeded from the transform it binds.
            deps=["preamble"], kind=ARTIFACT,
            label="H&E landmark registration",
            outputs=["he_affine", "he_residuals"],
        )

    def _landmark_preview() -> Preview:
        """The points as they stand, flipped into the frame the fit works in.

        Read-only, and it has to stay that way: drawing a preview must not move
        a landmark or touch the layer. Empty layers render as empty arrays,
        which is the honest picture of a button that cannot yet be pressed.
        """
        xen = he_state.get("xenium_lm_layer")
        he = he_state.get("he_lm_layer")
        xen_pts = np.asarray(getattr(xen, "data", np.empty((0, 2))), dtype=float)
        he_pts = np.asarray(getattr(he, "data", np.empty((0, 2))), dtype=float)
        n = min(len(xen_pts), len(he_pts))
        return Preview(
            blocks=["main"],
            params={"xenium_points": coerce(xen_pts[:n]),
                    "he_points": coerce(_flip_points(he_pts[:n]))},
        )

    def on_register():
        xen_pts = he_state["xenium_lm_layer"].data
        he_pts = he_state["he_lm_layer"].data
        n = min(len(xen_pts), len(he_pts))
        if n < 3:
            reg_status_label.value = "Need at least 3 paired landmarks"
            return
        xen_pts = np.asarray(xen_pts[:n], dtype=np.float64)
        # No arguments: the step reads the layers through the same provider the
        # Templates pane does, so the points that are fitted are the points the
        # pane showed. The H&E points are flipped in there -- they are read in
        # unflipped layer-data coordinates while the affine is composed as
        # `fine @ flip`, so the fit has to happen in the flipped frame. See the
        # frame rule in the module docstring.
        out = ctx.run_step(_landmark_step())
        affine, residuals = out["he_affine"], out["he_residuals"]
        he_state["affine_3x3"] = affine
        he_state["affine_source"] = "landmarks"
        _apply_he_affine()
        _refresh_nuclei_button()
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
        # No recording here: the fit above ran through `ctx.run_step`.
        _save_he_affine_to_sdata()
        from palms.utils.adata_persistence import save_landmarks_to_sdata
        # Persist the points as *clicked* — the layer's own data coordinates —
        # so restoring them puts each marker back where the user put it. The
        # flip is applied at fit time, not baked into the stored geometry.
        save_landmarks_to_sdata(ctx, 'he_xenium_landmarks', np.asarray(xen_pts))
        save_landmarks_to_sdata(ctx, 'he_he_landmarks',
                                np.asarray(he_state["he_lm_layer"].data[:n], dtype=np.float64))

    def on_save_landmarks():
        default_dir = str(data_path) if data_path else ""
        path, _ = QFileDialog.getSaveFileName(
            None, "Save Landmarks", default_dir + "/landmarks.json", "JSON Files (*.json)",
        )
        if not path:
            return
        xen_pts = np.asarray(he_state["xenium_lm_layer"].data, dtype=np.float64)
        # Saved in the frame the affine was fitted in, so the file is
        # self-consistent: `affine @ he_pts == xen_pts`, which is what
        # scripts/compare_he_registration.py and any other reader assume.
        he_pts = _flip_points(np.asarray(he_state["he_lm_layer"].data, dtype=np.float64))
        save_landmarks(
            path, xen_pts, he_pts,
            affine=he_state["affine_3x3"], he_filename=he_state["he_filename"],
            flip_v=he_flip_v.value, flip_h=he_flip_h.value,
        )
        reg_status_label.value = f"Landmarks saved to {Path(path).name}"
        # The points are inlined rather than referenced: landmarks can be saved
        # before a registration is computed, so ``he_xen_pts`` may not exist.
        affine = he_state["affine_3x3"]
        ctx.record_node(
            "he:save_landmarks",
            f"\n# Save H&E landmarks to {Path(path).name}\n"
            f"from palms.utils.registration import save_landmarks\n"
            f"save_landmarks(\n"
            f"    r\"{path}\",\n"
            f"    np.array({xen_pts.tolist()}),\n"
            f"    np.array({he_pts.tolist()}),\n"
            f"    affine={None if affine is None else f'np.array({np.asarray(affine).tolist()})'},\n"
            f"    he_filename={he_state['he_filename']!r},\n"
            f"    flip_v={he_flip_v.value}, flip_h={he_flip_h.value},\n"
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
        # The file stores the H&E points in the orientation they were fitted in;
        # the layer wants its own (unflipped) data coordinates. Adopt the file's
        # flips first, then undo them — both flips are involutions, so the same
        # transform converts each way.
        with he_flip_v.changed.blocked(), he_flip_h.changed.blocked():
            he_flip_v.value = bool(data.get("flip_v", False))
            he_flip_h.value = bool(data.get("flip_h", False))
        he_state["flip_v"] = he_flip_v.value
        he_state["flip_h"] = he_flip_h.value
        he_state["xenium_lm_layer"].data = data["xenium_landmarks_yx"]
        he_state["he_lm_layer"].data = _flip_points(data["he_landmarks_yx"])
        if "affine_3x3_yx" in data:
            affine = data["affine_3x3_yx"]
            he_state["affine_3x3"] = affine
            he_state["affine_source"] = "landmarks"
            _apply_he_affine()
            _refresh_nuclei_button()
            scale = np.sqrt(affine[0, 0]**2 + affine[0, 1]**2)
            reg_residuals_qt.setPlainText(f"Loaded affine (scale={scale:.4f})")
        if "he_filename" in data:
            he_state["he_filename"] = data["he_filename"]
        n = min(len(data["xenium_landmarks_yx"]), len(data["he_landmarks_yx"]))
        reg_status_label.value = f"Loaded {n} landmarks from {Path(path).name}"
        ctx.record_node(
            "he:load_landmarks",
            f"\n# Load H&E landmarks from file\n"
            f"from palms.utils.registration import load_landmarks\n"
            f"landmarks = load_landmarks(\"{path}\")",
            deps=["preamble"], kind=TERMINAL, label="Load H&E landmarks",
        )

    # Wire events
    he_load_button.clicked.connect(on_load_he)
    he_opacity_slider.changed.connect(on_he_opacity)
    coarse_align_button.clicked.connect(on_coarse_align)
    nuclei_align_button.clicked.connect(on_nuclei_align)
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

    # Registered so Tools -> Templates can show what each button would run, with
    # this tab's current settings rather than the template's sample values.
    ctx.state.setdefault("template_preview", {})[LOAD_TEMPLATE] = _he_load_preview
    ctx.state.setdefault("template_preview", {})[FLIP_TEMPLATE] = _flip_preview
    ctx.state.setdefault("template_preview", {})[COARSE_TEMPLATE] = _coarse_preview
    ctx.state.setdefault("template_preview", {})[NUCLEI_TEMPLATE] = _nuclei_preview
    ctx.state.setdefault(
        "template_preview", {})[LANDMARK_TEMPLATE] = _landmark_preview

    widget = make_tab(
        he_load_button,
        flip_row,
        he_opacity_slider,
        coarse_align_button,
        nuclei_align_button,
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
        # The pyramid on the canvas came from the cache, but the *path* is still
        # the honest record of where this H&E is, and it is what lets `he:load`
        # record a cell a notebook can replay against the raw output. Dropping it
        # here also silently erased it from the session on the next save, since
        # `_build_session_attrs` writes back whatever `he_state` holds.
        he_state["he_path"] = session_he_data.get("he_path")
        base = pyramid_rgb[0]
        # Resolve the placement before the flip widgets are read: on a store with
        # no viewer_session the element's own transform is the only record of the
        # registration, and it already contains the flip. See utils/registration_seed.
        from palms.utils.registration_seed import seed_registration
        reg = seed_registration(
            session_he_data, sdata, "he_image",
            affine_key="affine_3x3", coarse_key="coarse_affine",
            flip_v_key="flip_v", flip_h_key="flip_h", shape_key="he_shape_yx",
            element_shape_yx=(base.shape[0], base.shape[1]),
        )
        he_state["he_shape_yx"] = reg.shape_yx or (base.shape[0], base.shape[1])
        he_state["he_pixel_size_um"] = session_he_data.get("he_pixel_size_um")
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
        he_state["affine_source"] = session_he_data.get("affine_source")
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
        _refresh_nuclei_button()
        has_affine = he_state["affine_3x3"] is not None or he_state["coarse_affine"] is not None
        how = ""
        if has_affine:
            how = (" (registration from the image itself)" if reg.source == "element"
                   else " (with registration)")
        he_status_label.value = f"H&E restored: {he_filename}{how}"
        print(f"  Restored H&E from cache: {he_filename}"
              + (f" with registration ({reg.source})" if has_affine else ""))

    def _restore_session(session):
        if sdata is not None and "he_image" in sdata.images:
            # The flips are set in _on_he_restored_from_sdata, from the same
            # record the affine comes from — setting them here too would pair a
            # session flip with an element affine that already contains one.
            _session_he_data = {
                "affine_3x3": session.get("affine_3x3"),
                "affine_source": session.get("affine_source"),
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

            @thread_worker
            def _load_he_from_sdata():
                # Lazy on purpose: the eager version computed every level,
                # scale0 included, into dense numpy on every launch. The ARMS
                # copy of this was fixed in 9cad210 and this one was missed.
                # napari fetches only the tiles it draws from a dask multiscale.
                import dask.array as da
                he_dt = sdata.images["he_image"]
                pyramid = pyramid_levels(he_dt)
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
        "nuclei_align_button": nuclei_align_button,
        "create_landmark_layers": _create_landmark_layers,
        "apply_he_affine": _apply_he_affine,
    }

