"""Resolving an overlay's coordinate frame from live viewer state.

The element's stored transformation is not a trustworthy authority: a registered
H&E whose affine lives in the session reads as identity, and so do `arms_tiles`
and every `patch_*` overlay, which are placed by a *linked* layer's affine at
display time and never get one written to disk. Reading the identity and slicing
on it is what put a strip of the wrong picture into a real export.

No Qt and no napari here — the layers are stubs exposing the one attribute the
capture reads, which is also the point: `crop_state` must stay describable as
plain data so `crop_overlays` can stay pure.
"""

from __future__ import annotations

import numpy as np
import pytest

from palms.utils.crop_state import capture_overlay_frames


class _Affine:
    def __init__(self, m):
        self.affine_matrix = np.asarray(m, dtype=np.float64)


class _Layer:
    def __init__(self, name, affine=None):
        self.name = name
        self.affine = _Affine(affine if affine is not None else np.eye(3))


class _Ctx:
    """The handful of ViewerContext fields the capture reads."""

    def __init__(self, he_state=None, arms_state=None,
                 external_images_state=None, patch_overlays_state=None,
                 pixel_size=1.0):
        self.he_state = he_state or {}
        self.arms_state = arms_state or {}
        self.external_images_state = external_images_state or []
        self.patch_overlays_state = patch_overlays_state or []
        self.pixel_size = pixel_size
        self.sdata = None          # _morphology_shape returns None, which is fine


_REG = np.array([[-2.3535, 0.0, 34396.65],
                 [0.0, -2.3535, 53689.95],
                 [0.0, 0.0, 1.0]])


def test_the_live_layer_affine_is_what_gets_captured():
    ctx = _Ctx(he_state={"he_layer": _Layer("H&E (slide.tif)", _REG)})

    frames = capture_overlay_frames(ctx).frames

    np.testing.assert_allclose(frames["he_image"], _REG)


def test_the_stored_state_affine_is_the_fallback():
    """Before the layer exists — or when it was never given one."""
    ctx = _Ctx(he_state={"he_layer": None, "affine_3x3": _REG})

    np.testing.assert_allclose(capture_overlay_frames(ctx).frames["he_image"], _REG)


def test_an_identity_layer_yields_no_frame_at_all():
    """"Unregistered" and "registered to the identity" are the same thing.

    Reporting `eye(3)` would look like a resolved frame and defeat the
    credibility check downstream; absence lets the element have its say.
    """
    ctx = _Ctx(he_state={"he_layer": _Layer("H&E", np.eye(3))})

    assert "he_image" not in capture_overlay_frames(ctx).frames


def test_arms_tiles_are_declared_to_be_in_the_arms_image_frame():
    """The defect: tiles hold ARMS/SVS pixels while declaring identity on disk."""
    ctx = _Ctx(arms_state={"he_layer": _Layer("ARMS H&E (x.svs)", _REG)})

    got = capture_overlay_frames(ctx)

    assert got.companions["arms_tiles"] == "arms_he_image"
    assert got.companions["arms_he_landmarks"] == "arms_he_image"
    np.testing.assert_allclose(got.frames["arms_he_image"], _REG)


def test_he_landmarks_name_their_companion_even_with_no_registration():
    """The companion carries the slice origin, which matters at identity too."""
    assert capture_overlay_frames(_Ctx()).companions["he_he_landmarks"] == "he_image"


def test_a_patch_overlay_is_resolved_to_the_element_it_is_linked_to():
    """The link is by *layer name* — "H&E (None)" is a title, not an element id.

    Resolving it through the layer object rather than parsing the string is what
    keeps this working when the filename is None, which is exactly the case a
    real dataset had.
    """
    he_layer = _Layer("H&E (None)", _REG)
    ctx = _Ctx(
        he_state={"he_layer": he_layer},
        patch_overlays_state=[{"element_name": "patch_HE_R2",
                               "affine_source_name": "H&E (None)"}],
    )

    got = capture_overlay_frames(ctx)

    assert got.companions["patch_HE_R2"] == "he_image", (
        "a patch overlay linked to the H&E is in H&E pixels; read as Xenium data it "
        "is clipped to nothing, which is how three patch overlays were lost"
    )


def test_a_patch_overlay_linked_to_something_unnameable_still_gets_its_matrix():
    """No companion means no slice origin, but placement is better than nothing."""
    ctx = _Ctx(patch_overlays_state=[{
        "element_name": "patch_orphan",
        "affine_source_name": "some layer that is not an element",
        "shapes_layer": _Layer("patches", _REG),
    }])

    got = capture_overlay_frames(ctx)

    assert "patch_orphan" not in got.companions
    np.testing.assert_allclose(got.frames["patch_orphan"], _REG)


def test_external_images_and_their_landmarks():
    ctx = _Ctx(external_images_state=[{
        "element_name": "ext_phenocycler",
        "layer_ref": _Layer("PhenoCycler", _REG),
    }])

    got = capture_overlay_frames(ctx)

    np.testing.assert_allclose(got.frames["ext_phenocycler"], _REG)
    assert got.companions["ext_phenocycler_image_lm"] == "ext_phenocycler"


def test_a_serialised_affine_is_accepted_when_the_layer_is_gone():
    """Session UI rows keep `affine_matrix`; a torn-down layer must not lose it."""
    ctx = _Ctx(external_images_state=[{
        "element_name": "ext_x", "layer_ref": None,
        "affine_matrix": _REG.tolist(),
    }])

    np.testing.assert_allclose(capture_overlay_frames(ctx).frames["ext_x"], _REG)


def test_an_empty_context_is_not_an_error():
    got = capture_overlay_frames(_Ctx())
    assert got.frames == {}
    assert got.xenium_shape_yx is None


def test_the_napari_and_spatialdata_affine_conventions_agree():
    """Both are 3x3 (row, col) applied left-multiplied.

    Asserted rather than assumed: if napari ever changed this, every overlay would
    be silently transposed and nothing would raise.
    """
    napari = pytest.importorskip("napari")
    m = np.array([[2.0, 0.5, 3.0], [0.1, 3.0, 7.0], [0.0, 0.0, 1.0]])
    a = napari.utils.transforms.Affine(affine_matrix=m)

    np.testing.assert_allclose(a.affine_matrix, m)
    np.testing.assert_allclose(a([2.0, 3.0]), (m @ np.array([2.0, 3.0, 1.0]))[:2])


# ── the world/pixel boundary ─────────────────────────────────────────────────

def test_a_layer_affine_is_converted_out_of_world_units():
    """`layer.affine` is in micrometres; everything downstream is in pixels.

    Since the scale bar work every layer carries `scale = (pixel_size,)*ndim`, and
    napari applies a layer's affine *after* its scale — so `layer.affine` is a
    world-unit matrix. `crop_translation` and `overlay_pixel_bbox` are both in
    pixels. Passing the raw matrix through is a silent misplacement by a factor of
    1/pixel_size: on the real H&E registration below that is 27000 x 42000 px, a
    crop of the wrong part of the slide.

    Every other affine-persisting call site already goes through
    `units.layer_affine_px`; this is the one that has to as well.
    """
    from palms.utils.units import px_affine_to_world

    px = 0.2125
    reg_px = np.array([[-2.3535, 0.0, 34396.65],
                       [0.0, -2.3535, 53689.95],
                       [0.0, 0.0, 1.0]])
    # A layer placed exactly as the viewer places it: pixels -> world at the boundary.
    layer = _Layer("H&E (slide.tif)", px_affine_to_world(reg_px, px))

    got = capture_overlay_frames(
        _Ctx(he_state={"he_layer": layer}, pixel_size=px)).frames["he_image"]

    np.testing.assert_allclose(got, reg_px, atol=1e-6)


def test_a_stored_state_affine_is_already_in_pixels_and_is_not_converted():
    """The asymmetry: only the *layer* is in world units.

    `he_state["affine_3x3"]` is the pixel-space matrix the store holds, so
    converting it too would break it in the opposite direction.
    """
    reg_px = np.array([[2.0, 0.0, 100.0], [0.0, 2.0, 200.0], [0.0, 0.0, 1.0]])

    got = capture_overlay_frames(
        _Ctx(he_state={"he_layer": None, "affine_3x3": reg_px},
             pixel_size=0.2125)).frames["he_image"]

    np.testing.assert_allclose(got, reg_px)
