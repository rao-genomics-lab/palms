"""Carry registered overlays and drawn regions through a dataset crop.

``crop_export`` builds the exported store from five core elements. Everything a
user *added* — a registered H&E, an ARMS section, external images, patch
overlays, ROIs, annotations, ARMS tiles, registration landmarks — lived only in
the source store, so a crop silently threw away the work that made the dataset
worth cropping. This module is the part that brings it along.

Pure functions over spatialdata elements: no ``ViewerContext``, no Qt, no napari.
That is what lets them be tested against synthetic elements, which matters here
because the failure mode is geometric — an overlay that lands in the wrong place
still *looks* like a successful export.

The one rule everything else follows
------------------------------------
**An element's transformation says which space it is in.** An element with a
non-identity affine holds coordinates in its own pixel space and is placed into
Xenium morphology pixel space by that affine (``adata_persistence`` writes it
with ``input_axes=("y","x") → output_axes=("y","x")``). An element with no
transformation is already in Xenium pixel space.

So the crop, which is a pure translation ``T`` in Xenium pixel space, composes
one of two ways:

* identity element — move the geometry itself by ``-(col_min, row_min)``, and
  clip it to the crop region, because the geometry *is* the data;
* affine element — leave its own coordinates alone and rewrite the transform to
  ``T ∘ A``, because its coordinates mean nothing in Xenium space by themselves.

Deriving this from the transformation rather than from a list of element names
is deliberate: the loader's own comments record that a fixed name list has
already gone stale once, and new overlay kinds are named per file (``ext_``,
``patch_``) precisely so they cannot be enumerated up front.

The exception the disk forces: image-space landmarks
---------------------------------------------------
The rule above would be the whole story if every element declared its frame. It
does not. ``save_overlay_affine_to_sdata`` is only ever called with an image or
patch element's own name, and ``_save_he_affine_to_sdata`` writes only to
``images["he_image"]`` — so **no landmark element carries a transformation on
disk**, including the ones whose coordinates are H&E pixels
(``he_he_landmarks``, ``arms_he_landmarks``, ``*_image_lm``). Read by the rule
alone they look like Xenium-space elements and would be translated by the crop
origin, quietly corrupting the registration they exist to reconstruct.

They are recognised by name and passed through verbatim, which is also exactly
what the source store holds — the viewer applies the *image's* affine to the
landmark layer at display time and reads the raw coordinates back off disk.

Why landmarks are not filtered
------------------------------
It is tempting to drop registration landmarks that fall outside the crop. Doing
so silently *changes the registration*: the affine is a least-squares fit over
the landmark set, so a subset fits a different transform, and the exported
dataset would re-derive an affine that no longer matches the one it shipped
with. Every landmark is kept, translated into the crop frame if it was in Xenium
space. A landmark outside the cropped image is harmless — it is a reference
point, not data — and keeping it means a re-fit in the exported dataset
reproduces the original registration exactly.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)

#: Landmark sets whose coordinates are the *overlay's* pixels even though nothing
#: on disk says so. See "The exception the disk forces" above. Recognised by name
#: because there is no other signal; the suffix covers external images, whose
#: element names are per-file.
_IMAGE_SPACE_SHAPE_NAMES = frozenset({"he_he_landmarks", "arms_he_landmarks"})
_IMAGE_SPACE_SHAPE_SUFFIXES = ("_image_lm",)

#: Landmark sets in Xenium pixel space. Translated with the crop but never
#: clipped — dropping one re-fits the registration; see the module docstring.
_XENIUM_LANDMARK_NAMES = frozenset({"he_xenium_landmarks", "arms_xenium_landmarks"})
_XENIUM_LANDMARK_SUFFIXES = ("_xenium_lm",)


def is_image_space_shape(name: str) -> bool:
    """True when *name* holds coordinates in an overlay's pixels, not Xenium's."""
    return (name in _IMAGE_SPACE_SHAPE_NAMES
            or name.endswith(_IMAGE_SPACE_SHAPE_SUFFIXES))


def is_landmark_shape(name: str) -> bool:
    """True for a registration landmark set in either frame."""
    return (is_image_space_shape(name)
            or name in _XENIUM_LANDMARK_NAMES
            or name.endswith(_XENIUM_LANDMARK_SUFFIXES))


class OverlayCropError(Exception):
    """An overlay could not be carried through the crop."""


# ── transformations ──────────────────────────────────────────────────────────

def element_affine(element) -> np.ndarray:
    """The 3x3 (y, x) affine mapping *element*'s own pixels into Xenium space.

    Identity when the element carries no transformation, which is how an element
    declares "my coordinates are already Xenium pixels".
    """
    from spatialdata.transformations import get_transformation

    try:
        t = get_transformation(element, "global")
    except (KeyError, ValueError):
        return np.eye(3)
    if t is None:
        return np.eye(3)
    try:
        return np.asarray(
            t.to_affine_matrix(input_axes=("y", "x"), output_axes=("y", "x")),
            dtype=np.float64,
        )
    except (AttributeError, ValueError) as e:      # pragma: no cover - defensive
        log.warning("could not read a transformation, treating it as identity: %s", e)
        return np.eye(3)


def crop_translation(row_min: int, col_min: int) -> np.ndarray:
    """``T``: the crop's own translation in Xenium pixel space, as a (y, x) affine."""
    t = np.eye(3)
    t[0, 2] = -float(row_min)
    t[1, 2] = -float(col_min)
    return t


def set_element_affine(element, affine_yx: np.ndarray):
    """Attach *affine_yx* to *element* as its ``global`` transformation.

    An identity affine is written as ``Identity()`` rather than as a degenerate
    ``Affine``, so an element that is genuinely in Xenium space keeps saying so
    and the rule in the module docstring stays readable off disk.
    """
    from spatialdata.transformations import Affine, Identity, set_transformation

    m = np.asarray(affine_yx, dtype=np.float64)
    if np.allclose(m, np.eye(3)):
        set_transformation(element, Identity(), "global")
    else:
        set_transformation(
            element, Affine(m, input_axes=("y", "x"), output_axes=("y", "x")), "global",
        )
    return element


# ── rasters ──────────────────────────────────────────────────────────────────

def overlay_pixel_bbox(affine_yx: np.ndarray, crop_bbox: tuple, shape_yx: tuple) -> tuple | None:
    """Where *crop_bbox* lands in an overlay's own pixel space.

    ``crop_bbox`` is ``(row_min, row_max, col_min, col_max)`` in Xenium pixels.
    Returns the axis-aligned bbox of the crop's four corners mapped through
    ``affine_yx⁻¹``, clipped to ``shape_yx``, or ``None`` when the overlay does
    not reach the crop at all.

    Registration fits a *similarity* transform (rotation, uniform scale,
    translation), so a rotated overlay gives an axis-aligned box larger than the
    crop itself. That is correct rather than sloppy: the extra margin is carried
    along and the transform still places every pixel where it belongs.
    """
    row_min, row_max, col_min, col_max = crop_bbox
    corners = np.array([
        [row_min, col_min, 1.0], [row_min, col_max, 1.0],
        [row_max, col_min, 1.0], [row_max, col_max, 1.0],
    ]).T
    try:
        inv = np.linalg.inv(np.asarray(affine_yx, dtype=np.float64))
    except np.linalg.LinAlgError:
        raise OverlayCropError("overlay transformation is singular and cannot be inverted")

    mapped = (inv @ corners)[:2]                       # 2xN, rows then cols
    r0 = int(np.floor(mapped[0].min()))
    r1 = int(np.ceil(mapped[0].max()))
    c0 = int(np.floor(mapped[1].min()))
    c1 = int(np.ceil(mapped[1].max()))

    r0, c0 = max(0, r0), max(0, c0)
    r1, c1 = min(int(shape_yx[0]), r1), min(int(shape_yx[1]), c1)
    if r1 <= r0 or c1 <= c0:
        return None
    return r0, r1, c0, c1


def crop_raster_overlay(element, crop_bbox: tuple, scale_factors_fn):
    """Crop a registered image/label overlay to the region the crop covers.

    Returns ``(new_element, is_labels)`` or ``None`` if the overlay does not
    overlap the crop.

    Stays lazy throughout: the level-0 dask array is *sliced*, never computed.
    ``crop_export`` carries long comments about why — materialising a
    full-resolution raster here has caused a real OOM — and an overlay is no
    smaller than the morphology image it was registered against.
    """
    from spatialdata.models import Image2DModel, Labels2DModel
    from xenium_viewer.utils.crop_export import _extract_dt_scales

    scales = _extract_dt_scales(element)
    if not scales:
        raise OverlayCropError("overlay has no readable image data")
    full = scales[0]

    is_labels = full.ndim == 2
    shape_yx = full.shape if is_labels else full.shape[-2:]

    affine = element_affine(element)
    own = overlay_pixel_bbox(affine, crop_bbox, shape_yx)
    if own is None:
        return None
    r0, r1, c0, c1 = own

    cropped = full[r0:r1, c0:c1] if is_labels else full[:, r0:r1, c0:c1]
    factors = scale_factors_fn(cropped.shape[-2], cropped.shape[-1])

    # The overlay was sliced, so its own origin moved to (r0, c0); undo that
    # before applying the registration, then apply the crop's own translation.
    # Reading right to left: own-pixel -> unsliced own-pixel -> Xenium -> crop.
    origin = np.eye(3)
    origin[0, 2] = float(r0)
    origin[1, 2] = float(c0)
    row_min, _, col_min, _ = crop_bbox
    combined = crop_translation(row_min, col_min) @ affine @ origin

    model = Labels2DModel if is_labels else Image2DModel
    dims = ("y", "x") if is_labels else ("c", "y", "x")
    parsed = model.parse(cropped, dims=dims, scale_factors=factors)
    return set_element_affine(parsed, combined), is_labels


# ── vectors ──────────────────────────────────────────────────────────────────

def crop_vector_overlay(element, name: str, crop_bbox: tuple, crop_polygon_xy):
    """Carry a shapes element through the crop, or ``None`` if nothing survives.

    Three cases, and each is a different question about what the coordinates mean:

    * **image-space landmarks** — passed through untouched (see the module
      docstring: nothing on disk marks their frame, so the rule cannot see them);
    * **Xenium-space landmarks** — translated into the crop frame but never
      clipped, because dropping one re-fits the registration;
    * **everything else is data** — clipped to the drawn region and translated.
      When the element carries its own affine its geometry is not in Xenium
      pixels at all, so the region is mapped into *its* frame first and the
      coordinates are left alone.

    That last case is not hypothetical tidiness: ``cell_circles`` is stored in
    microns with a 1/pixel_size scale, and patch overlays in their source image's
    pixels. Carrying either one whole would ship an export whose circles or
    patches cover cells the table no longer contains.

    Every non-geometry column survives — ``arms_tiles`` carries its
    ``cluster_id``, patch overlays their per-patch cluster columns — because
    geometry without them means nothing.
    """
    import geopandas as gpd
    from shapely import make_valid
    from shapely.affinity import translate as shapely_translate
    from spatialdata.models import ShapesModel

    gdf = element
    if len(gdf) == 0:
        return None

    affine = element_affine(gdf)
    row_min, _, col_min, _ = crop_bbox
    identity = np.allclose(affine, np.eye(3))

    if is_image_space_shape(name):
        return set_element_affine(ShapesModel.parse(gdf.copy()), affine)

    out = gdf.copy()

    if identity:
        out["geometry"] = out["geometry"].apply(
            lambda g: shapely_translate(g, xoff=-float(col_min), yoff=-float(row_min))
        )
        new_affine = np.eye(3)
        region = crop_polygon_xy
        if region is not None:
            region = shapely_translate(region, xoff=-float(col_min), yoff=-float(row_min))
    else:
        # Own-frame geometry: leave the coordinates, compose the transform, and
        # bring the crop region *down* into this element's frame to clip there.
        new_affine = crop_translation(row_min, col_min) @ affine
        region = _region_in_own_frame(crop_polygon_xy, affine)

    if not is_landmark_shape(name) and region is not None:
        if not region.is_valid:
            # make_valid, not buffer(0): buffer(0) silently *deletes* a lobe of a
            # self-intersecting polygon rather than repairing it.
            region = make_valid(region)
        out["geometry"] = out["geometry"].apply(
            lambda g: (g if g.is_valid else make_valid(g)).intersection(region)
        )
        # An explicit mask rather than `~is_empty & notna()`: geopandas warns that
        # notna()'s treatment of empty geometry changed, and "which rows survived a
        # clip" is not a question to answer through an API in the middle of
        # changing what it means.
        keep = np.array([g is not None and not g.is_empty for g in out["geometry"]], dtype=bool)
        out = out[keep]
        if len(out) == 0:
            return None

    out = gpd.GeoDataFrame(out, geometry="geometry")
    return set_element_affine(ShapesModel.parse(out), new_affine)


def _region_in_own_frame(crop_polygon_xy, affine_yx: np.ndarray):
    """The crop region expressed in an element's own coordinates.

    ``affine_yx`` is (y, x); shapely wants (x, y), so the rotation block is
    transposed about the anti-diagonal rather than simply reused — getting this
    backwards is a silent 90-degree error, not a crash.
    """
    from shapely.affinity import affine_transform

    if crop_polygon_xy is None:
        return None
    try:
        inv = np.linalg.inv(np.asarray(affine_yx, dtype=np.float64))
    except np.linalg.LinAlgError:
        raise OverlayCropError("overlay transformation is singular and cannot be inverted")

    # (y, x) matrix -> (x, y) matrix: swap both rows and columns of the 2x2 part,
    # and swap the translation components.
    a, b, ty = inv[0, 0], inv[0, 1], inv[0, 2]
    c, d, tx = inv[1, 0], inv[1, 1], inv[1, 2]
    # shapely's order is (a, b, d, e, xoff, yoff) for x' = a*x + b*y + xoff
    return affine_transform(crop_polygon_xy, (d, c, b, a, tx, ty))


# ── the whole set ────────────────────────────────────────────────────────────

def user_overlay_names(sdata) -> dict:
    """``{"images": [...], "labels": [...], "shapes": [...]}`` of carryable elements.

    Uses ``loader._is_user_element`` so this and the "does the cache hold user
    data" check can never disagree about what counts as the user's work.
    ``cell_circles`` is added because it is real Xenium output that the crop was
    dropping — it is not *user* data, so the loader's predicate rightly ignores
    it, but it should still survive a crop.
    """
    from xenium_viewer import loader

    out = {"images": [], "labels": [], "shapes": []}
    if sdata is None:
        return out
    for group, keys in (
        ("images", loader._USER_IMAGE_KEYS),
        ("labels", []),
        ("shapes", loader._USER_SHAPE_KEYS),
    ):
        try:
            present = list(getattr(sdata, group).keys())
        except (AttributeError, KeyError):            # pragma: no cover - defensive
            continue
        for elem_name in present:
            if elem_name in ("morphology_focus", "cell_labels", "nucleus_labels"):
                continue
            if loader._is_user_element(elem_name, list(keys)) or elem_name == "cell_circles":
                out[group].append(elem_name)
    return out
