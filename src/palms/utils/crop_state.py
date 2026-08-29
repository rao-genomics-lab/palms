"""The live viewer state a crop needs, captured on the GUI thread.

``crop_overlays`` decides where an overlay goes from *the element's own stored
transformation*. That is wrong often enough to have shipped a broken export: a
registered H&E whose affine lives in the session and was never mirrored onto the
element reads as identity, and so do ``arms_tiles`` and every ``patch_*`` overlay,
which are placed by a linked layer's affine at display time and never get one of
their own written to disk.

The authority that is always right is the one the user is looking at. A layer's
``affine`` is what puts the overlay on screen; if the export reproduces that, the
export matches what the user cropped.

Two constraints shape this module:

* ``crop_overlays`` stays pure — no ``ViewerContext``, no Qt, no napari — because
  that is what lets its geometry be tested against synthetic elements. So the
  resolution happens here and is passed *in* as plain data.
* ``crop_and_export`` runs on a ``QThread``. Reading ``layer.affine`` from a worker
  is a live-object read across threads, so :func:`capture_overlay_frames` is called
  from the GUI thread before the worker starts, and the worker receives a snapshot.
  This is the same discipline ``app._snapshot_layers`` already follows for
  ``save_session``.
"""

from __future__ import annotations

import logging
from typing import NamedTuple

import numpy as np

log = logging.getLogger(__name__)


class OverlayFrames(NamedTuple):
    """Everything the crop needs to know about where overlays live.

    ``frames`` maps element name -> 3x3 (y, x) affine into Xenium pixel space.
    ``companions`` maps a shapes element name -> the raster element whose pixel
    frame its geometry is written in. ``xenium_shape_yx`` is the morphology grid,
    the only shape against which an identity frame can be checked.
    """

    frames: dict
    companions: dict
    xenium_shape_yx: tuple | None


def _layer_affine(layer, pixel_size: float):
    """The 3x3 (y, x) affine of a napari layer **in Xenium pixels**, or None.

    ``layer.affine`` is in *world* units — micrometres, since every layer carries
    ``scale = (pixel_size, pixel_size)`` and napari applies a layer's affine after
    its scale. Everything downstream of here is in pixels: ``crop_translation`` is
    a pixel offset and ``overlay_pixel_bbox`` maps a pixel bounding box. Handing
    the raw world matrix on is a silent misplacement by a factor of
    ``1 / pixel_size`` — measured at 27000 x 42000 px on a real H&E registration,
    which is a crop of the wrong part of the slide, not a slightly wrong one.

    So the conversion happens here, through ``units.layer_affine_px``, the one
    function that owns that boundary. Identity comes back as ``None`` rather than
    ``eye(3)``: "this layer is unregistered" and "this layer is registered to the
    identity" are the same thing, and the caller wants to fall through to the
    element.
    """
    if layer is None:
        return None
    try:
        from palms.utils.units import layer_affine_px

        m = layer_affine_px(layer, pixel_size)
    except Exception:                                  # pragma: no cover - defensive
        return None
    if m is None or np.allclose(m, np.eye(3), atol=1e-6):
        return None
    return np.asarray(m, dtype=np.float64)


def _state_affine(state, pixel_size: float):
    """The placement recorded in ``he_state`` / ``arms_state``.

    Prefers the live layer, because that is what is on screen; falls back to the
    stored ``affine_3x3``. The stored copy is **already in pixels** — that is the
    frame the store and this module use — so only the layer needs converting,
    which is exactly the asymmetry ``_layer_affine`` encodes.
    """
    if not state:
        return None
    live = _layer_affine(state.get("he_layer"), pixel_size)
    if live is not None:
        return live
    for key in ("affine_3x3", "coarse_affine"):
        m = state.get(key)
        if m is not None:
            return np.asarray(m, dtype=np.float64)
    return None


def capture_overlay_frames(ctx) -> OverlayFrames:
    """Resolve every overlay's frame from live viewer state. GUI thread only.

    Precedence per element is **live layer -> stored element transform ->
    identity**, and only the first is captured here; ``crop_export`` applies the
    other two, because they need the element and this function deliberately does
    not touch the store.
    """
    # Every layer affine read below is converted to pixels against this.
    pixel_size = float(getattr(ctx, "pixel_size", 1.0) or 1.0)

    frames: dict = {}
    companions: dict = {}

    # ── H&E and its landmarks ────────────────────────────────────────────────
    he = _state_affine(getattr(ctx, "he_state", None), pixel_size)
    if he is not None:
        frames["he_image"] = he
    # The landmark set is in H&E pixels whether or not we resolved an affine —
    # naming the companion is what carries the slice origin across (see
    # crop_overlays.crop_vector_overlay), and that is useful even at identity.
    companions["he_he_landmarks"] = "he_image"

    # ── ARMS: image, its landmarks, and the tiles ────────────────────────────
    arms = _state_affine(getattr(ctx, "arms_state", None), pixel_size)
    if arms is not None:
        frames["arms_he_image"] = arms
    companions["arms_he_landmarks"] = "arms_he_image"
    # The tiles are drawn with the ARMS image's affine (tab_arms assigns the same
    # matrix to shapes_layer), so they are in ARMS pixels — not Xenium, which is
    # what the element's identity transform would otherwise say.
    companions["arms_tiles"] = "arms_he_image"

    # ── External images ──────────────────────────────────────────────────────
    for entry in (getattr(ctx, "external_images_state", None) or []):
        name = entry.get("element_name")
        if not name:
            continue
        m = _layer_affine(entry.get("layer_ref"), pixel_size)
        if m is None and entry.get("affine_matrix") is not None:
            m = np.asarray(entry["affine_matrix"], dtype=np.float64)
        if m is not None:
            frames[name] = m
        lm = f"{name}_image_lm"
        companions[lm] = name

    # ── Patch overlays ───────────────────────────────────────────────────────
    # A patch overlay's geometry is in its *source image's* pixels. The link is
    # by layer name, so resolve that name back to an element rather than parsing
    # it — "H&E (None)" is a layer title, not an element id.
    layer_to_element = {}
    for state, element in ((getattr(ctx, "he_state", None), "he_image"),
                           (getattr(ctx, "arms_state", None), "arms_he_image")):
        layer = (state or {}).get("he_layer")
        if layer is not None:
            layer_to_element[getattr(layer, "name", None)] = element
    for entry in (getattr(ctx, "external_images_state", None) or []):
        layer = entry.get("layer_ref")
        if layer is not None and entry.get("element_name"):
            layer_to_element[getattr(layer, "name", None)] = entry["element_name"]

    for entry in (getattr(ctx, "patch_overlays_state", None) or []):
        name = entry.get("element_name")
        if not name:
            continue
        source_element = layer_to_element.get(entry.get("affine_source_name"))
        if source_element is not None:
            companions[name] = source_element
        else:
            # Linked to something we cannot name as an element — keep the matrix
            # so the geometry is at least placed, without a slice origin.
            m = _layer_affine(entry.get("shapes_layer"), pixel_size)
            if m is None and entry.get("affine_matrix") is not None:
                m = np.asarray(entry["affine_matrix"], dtype=np.float64)
            if m is not None:
                frames[name] = m

    return OverlayFrames(frames=frames, companions=companions,
                         xenium_shape_yx=_morphology_shape(ctx))


def _morphology_shape(ctx) -> tuple | None:
    """The morphology grid, ``(y, x)`` — the yardstick for a credible identity."""
    try:
        from palms.utils.crop_export import _extract_dt_scales

        scales = _extract_dt_scales(ctx.sdata.images["morphology_focus"])
        if not scales:
            return None
        shape = scales[0].shape
        return tuple(int(v) for v in shape[-2:])
    except Exception:                                  # pragma: no cover - defensive
        log.debug("could not read the morphology shape", exc_info=True)
        return None
