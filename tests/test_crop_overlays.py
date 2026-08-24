"""Carrying registered overlays and drawn regions through a dataset crop.

The failure mode here is geometric, not exceptional: an overlay placed at the
wrong offset still exports cleanly, still opens, and still *looks* like a
registered H&E — it is simply on the wrong tissue. So these tests assert
placement, by round-tripping a landmark pixel through the exported transform and
checking where it lands, rather than asserting that a function returned
something.

No dataset and no GUI: every element is built here, which is also the only way
to test this at all — a rotated, scaled, registered overlay is exactly what no
local dataset happens to have.
"""

from __future__ import annotations

import numpy as np
import pytest

import geopandas as gpd
from shapely.geometry import Point, Polygon

from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Affine, Identity, set_transformation

from xenium_viewer.utils import crop_overlays
from xenium_viewer.utils.crop_export import _safe_scale_factors


# ── fixtures ─────────────────────────────────────────────────────────────────

def _similarity(theta_deg: float, scale: float, ty: float, tx: float) -> np.ndarray:
    """A (y, x) similarity affine — what utils/registration.py actually fits."""
    t = np.deg2rad(theta_deg)
    c, s = np.cos(t) * scale, np.sin(t) * scale
    return np.array([[c, -s, ty],
                     [s,  c, tx],
                     [0,  0,  1.0]])


def _image(shape_yx=(400, 500), channels=3, affine=None):
    import dask.array as da
    arr = da.zeros((channels, *shape_yx), dtype=np.uint8, chunks=(channels, 128, 128))
    el = Image2DModel.parse(arr, dims=("c", "y", "x"), scale_factors=[2])
    set_transformation(
        el,
        Identity() if affine is None
        else Affine(affine, input_axes=("y", "x"), output_axes=("y", "x")),
        "global",
    )
    return el


def _labels(shape_yx=(400, 500), affine=None):
    import dask.array as da
    arr = da.zeros(shape_yx, dtype=np.uint32, chunks=(128, 128))
    el = Labels2DModel.parse(arr, dims=("y", "x"), scale_factors=[2])
    set_transformation(
        el,
        Identity() if affine is None
        else Affine(affine, input_axes=("y", "x"), output_axes=("y", "x")),
        "global",
    )
    return el


def _shapes(gdf, affine=None):
    el = ShapesModel.parse(gdf)
    set_transformation(
        el,
        Identity() if affine is None
        else Affine(affine, input_axes=("y", "x"), output_axes=("y", "x")),
        "global",
    )
    return el


CROP = (100, 300, 150, 450)          # row_min, row_max, col_min, col_max


def _place(element, own_yx):
    """Where own-pixel *own_yx* lands, per the element's stored transformation."""
    a = crop_overlays.element_affine(element)
    y, x = own_yx
    out = a @ np.array([y, x, 1.0])
    return out[0], out[1]


# ── rasters ──────────────────────────────────────────────────────────────────

def test_registered_overlay_lands_in_the_same_place_after_the_crop():
    """The whole point: a pixel on a tissue feature stays on that feature.

    A rotated, scaled, translated H&E is cropped alongside the core data. A
    point in the overlay's own pixels must end up exactly ``(row_min, col_min)``
    lower in the exported frame than it was in the source frame — that, and
    nothing weaker, is what "the overlay still lines up" means.
    """
    affine = _similarity(theta_deg=17.0, scale=0.8, ty=60.0, tx=-30.0)
    src = _image(affine=affine)

    own = (120.0, 200.0)
    before = _place(src, own)

    out, is_labels = crop_overlays.crop_raster_overlay(src, CROP, _safe_scale_factors)
    assert not is_labels

    # The export sliced the overlay, so that same tissue point now sits at
    # own-pixel (own - slice origin). Ask the module where the slice started —
    # then map that *specific* pixel forward and check where it lands. Solving
    # for the pixel with the inverse and mapping it back would only assert
    # A·A⁻¹ = I, which is true of any affine at all.
    r0, _, c0, _ = crop_overlays.overlay_pixel_bbox(affine, CROP, (400, 500))
    own_after = np.array([own[0] - r0, own[1] - c0, 1.0])
    landed = crop_overlays.element_affine(out) @ own_after

    assert landed[0] == pytest.approx(before[0] - CROP[0], abs=1e-9)
    assert landed[1] == pytest.approx(before[1] - CROP[2], abs=1e-9)


def test_identity_overlay_is_translated_by_the_crop_origin():
    """An unregistered raster is already in Xenium pixels, so T alone applies."""
    src = _image(affine=None)
    out, _ = crop_overlays.crop_raster_overlay(src, CROP, _safe_scale_factors)
    a = crop_overlays.element_affine(out)
    # own-pixel (r0, c0) of the slice is Xenium (100, 150) -> crop-frame (0, 0)
    assert (a @ np.array([0.0, 0.0, 1.0]))[:2] == pytest.approx([0.0, 0.0])


def test_label_overlay_round_trips_as_labels():
    out = crop_overlays.crop_raster_overlay(_labels(), CROP, _safe_scale_factors)
    assert out is not None
    _, is_labels = out
    assert is_labels, "a 2-D raster must come back as a Labels element, not an Image"


def test_overlay_that_misses_the_crop_is_skipped_not_raised():
    """A stray registered image far from the crop must not fail the export."""
    far = _similarity(0.0, 1.0, ty=100_000.0, tx=100_000.0)
    assert crop_overlays.crop_raster_overlay(_image(affine=far), CROP, _safe_scale_factors) is None


def test_raster_crop_never_materialises_the_array():
    """crop_export's OOM history is the reason this stays lazy end to end."""
    import dask.array as da
    out, _ = crop_overlays.crop_raster_overlay(
        _image(shape_yx=(2000, 2000)), CROP, _safe_scale_factors)
    scales = [out[k].ds["image"].data for k in out.children] if hasattr(out, "children") else []
    assert scales, "expected a multiscale DataTree"
    for arr in scales:
        assert isinstance(arr, da.Array), "an eager .compute() crept into the overlay crop"


def test_singular_transform_is_reported_not_silently_wrong():
    bad = np.array([[1.0, 1.0, 0.0], [2.0, 2.0, 0.0], [0.0, 0.0, 1.0]])   # rank-deficient
    with pytest.raises(crop_overlays.OverlayCropError):
        crop_overlays.crop_raster_overlay(_image(affine=bad), CROP, _safe_scale_factors)


# ── vectors ──────────────────────────────────────────────────────────────────

def _crop_polygon():
    r0, r1, c0, c1 = CROP
    return Polygon([(c0, r0), (c1, r0), (c1, r1), (c0, r1)])


def test_roi_inside_the_crop_is_translated_into_the_crop_frame():
    roi = Polygon([(200, 150), (300, 150), (300, 250), (200, 250)])   # xy, inside
    el = _shapes(gpd.GeoDataFrame({"geometry": [roi]}))
    out = crop_overlays.crop_vector_overlay(el, "rois", CROP, _crop_polygon())
    minx, miny, maxx, maxy = out.geometry.iloc[0].bounds
    assert (minx, miny) == pytest.approx((200 - CROP[2], 150 - CROP[0]))
    assert (maxx, maxy) == pytest.approx((300 - CROP[2], 250 - CROP[0]))


def test_roi_straddling_the_boundary_is_clipped_to_the_region():
    straddle = Polygon([(100, 150), (250, 150), (250, 250), (100, 250)])  # x starts left of c0=150
    el = _shapes(gpd.GeoDataFrame({"geometry": [straddle]}))
    out = crop_overlays.crop_vector_overlay(el, "rois", CROP, _crop_polygon())
    minx = out.geometry.iloc[0].bounds[0]
    assert minx == pytest.approx(0.0), "the part outside the crop should have been clipped away"


def test_roi_entirely_outside_the_crop_yields_nothing():
    outside = Polygon([(1000, 1000), (1100, 1000), (1100, 1100), (1000, 1100)])
    el = _shapes(gpd.GeoDataFrame({"geometry": [outside]}))
    assert crop_overlays.crop_vector_overlay(el, "rois", CROP, _crop_polygon()) is None


def test_tile_attributes_survive_the_crop():
    """Geometry without its cluster_id is geometry with no meaning."""
    tiles = gpd.GeoDataFrame({
        "geometry": [Polygon([(200, 150), (260, 150), (260, 210), (200, 210)]),
                     Polygon([(300, 200), (360, 200), (360, 260), (300, 260)])],
        "tile_name": ["A1", "B2"],
        "cluster_id": np.array([3, 7], dtype=np.int32),
    })
    out = crop_overlays.crop_vector_overlay(_shapes(tiles), "arms_tiles", CROP, _crop_polygon())
    assert list(out["tile_name"]) == ["A1", "B2"]
    assert list(out["cluster_id"]) == [3, 7]


def test_self_intersecting_region_is_repaired_not_lobotomised():
    """buffer(0) drops a lobe of a bowtie; make_valid keeps both."""
    bowtie = Polygon([(200, 160), (300, 260), (300, 160), (200, 260)])
    el = _shapes(gpd.GeoDataFrame({"geometry": [bowtie]}))
    out = crop_overlays.crop_vector_overlay(el, "rois", CROP, _crop_polygon())
    assert out is not None and out.geometry.iloc[0].area > 0


# ── landmarks: the pairing rule ──────────────────────────────────────────────

def _landmarks(points_xy, affine=None):
    return _shapes(
        gpd.GeoDataFrame({"geometry": [Point(x, y) for x, y in points_xy], "radius": 1.0}),
        affine=affine,
    )


def test_landmark_pairs_stay_the_same_length_even_when_some_fall_outside():
    """Dropping an out-of-crop landmark would silently re-fit the registration.

    ``he_xenium_landmarks[i]`` pairs with ``he_he_landmarks[i]``. The affine is a
    least-squares fit over the whole set, so filtering it — however carefully —
    yields a *different* transform than the one the export ships. Both sets must
    come through whole.
    """
    xen = _landmarks([(200, 160), (5, 5), (400, 280), (9000, 9000)])   # two outside the crop
    he = _landmarks([(10, 10), (20, 20), (30, 30), (40, 40)],
                    affine=_similarity(5.0, 1.2, 3.0, 4.0))

    out_xen = crop_overlays.crop_vector_overlay(xen, "he_xenium_landmarks", CROP, _crop_polygon())
    out_he = crop_overlays.crop_vector_overlay(he, "he_he_landmarks", CROP, _crop_polygon())

    assert len(out_xen) == 4, "a Xenium landmark outside the crop must not be dropped"
    assert len(out_he) == len(out_xen), "landmark sets must stay paired index-for-index"


def test_xenium_landmarks_move_into_the_crop_frame():
    out = crop_overlays.crop_vector_overlay(
        _landmarks([(200, 160)]), "he_xenium_landmarks", CROP, _crop_polygon())
    p = out.geometry.iloc[0]
    assert (p.x, p.y) == pytest.approx((200 - CROP[2], 160 - CROP[0]))


@pytest.mark.parametrize("name", ["he_he_landmarks", "arms_he_landmarks", "phenocycler_image_lm"])
def test_image_space_landmarks_are_left_exactly_alone(name):
    """The case the transformation rule cannot see, and would get wrong.

    These hold H&E/overlay pixels, but **nothing on disk says so**:
    ``save_overlay_affine_to_sdata`` is only ever called with an image or patch
    element's name, and ``_save_he_affine_to_sdata`` writes only to
    ``images["he_image"]``. So they arrive with an identity transformation and
    read exactly like Xenium-space geometry. Translating them by the crop origin
    would corrupt the registration they exist to reconstruct — silently, since
    the export would still succeed and still open.

    Note the fixture passes *no* affine, matching what the store actually holds.
    """
    src = _landmarks([(10, 10), (20, 20)])
    out = crop_overlays.crop_vector_overlay(src, name, CROP, _crop_polygon())

    assert [(p.x, p.y) for p in out.geometry] == [(10, 10), (20, 20)]
    assert np.allclose(crop_overlays.element_affine(out), np.eye(3))


def test_xenium_and_image_landmark_sets_are_told_apart():
    """The naming rule itself, since it is the only signal available."""
    for xen in ("he_xenium_landmarks", "arms_xenium_landmarks", "phenocycler_xenium_lm"):
        assert not crop_overlays.is_image_space_shape(xen)
    for img in ("he_he_landmarks", "arms_he_landmarks", "phenocycler_image_lm"):
        assert crop_overlays.is_image_space_shape(img)


# ── which elements are carried ───────────────────────────────────────────────

def test_the_carry_list_comes_from_the_loader_predicate():
    """Not a second hand-written list — the loader's is the single source."""
    from xenium_viewer import loader

    class _FakeSData:
        images = {"morphology_focus": 1, "he_image": 1, "arms_he_image": 1, "ext_phenocycler": 1}
        labels = {"cell_labels": 1, "nucleus_labels": 1}
        shapes = {"rois": 1, "annotations": 1, "arms_tiles": 1, "cell_circles": 1,
                  "he_xenium_landmarks": 1, "phenocycler_image_lm": 1, "patch_tumour": 1}

    got = crop_overlays.user_overlay_names(_FakeSData())

    assert set(got["images"]) == {"he_image", "arms_he_image", "ext_phenocycler"}
    assert got["labels"] == [], "the two core rasters are handled by crop_export itself"
    assert set(got["shapes"]) == {
        "rois", "annotations", "arms_tiles", "cell_circles",
        "he_xenium_landmarks", "phenocycler_image_lm", "patch_tumour",
    }
    # The core elements must never appear, whatever the predicate says.
    for group in got.values():
        assert "morphology_focus" not in group and "cell_labels" not in group

    # And the predicate really is the loader's, not a copy.
    assert loader._is_user_element("ext_anything", loader._USER_IMAGE_KEYS)
    assert not loader._is_user_element("morphology_focus", loader._USER_IMAGE_KEYS)


# ── own-frame clipping ───────────────────────────────────────────────────────

PIXEL_SIZE = 0.2125


def test_micron_space_circles_are_clipped_in_their_own_frame():
    """``cell_circles`` is stored in microns with a 1/pixel_size scale.

    Carrying it whole would ship an export whose circles cover cells the table no
    longer contains, and translating it by the crop's *pixel* origin would move
    every circle by a factor of ~4.7. Neither raises; both are wrong. So the crop
    region is brought down into micron space and the clip happens there, leaving
    the coordinates untouched.
    """
    scale = np.diag([1.0 / PIXEL_SIZE, 1.0 / PIXEL_SIZE, 1.0])   # microns -> Xenium px

    inside_px, outside_px = (200.0, 250.0), (50.0, 60.0)         # (row, col)
    circles = gpd.GeoDataFrame({
        "geometry": [Point(inside_px[1] * PIXEL_SIZE, inside_px[0] * PIXEL_SIZE),
                     Point(outside_px[1] * PIXEL_SIZE, outside_px[0] * PIXEL_SIZE)],
        "radius": 3.0,
    })
    out = crop_overlays.crop_vector_overlay(
        _shapes(circles, affine=scale), "cell_circles", CROP, _crop_polygon())

    assert len(out) == 1, "the circle outside the crop should have been dropped"
    kept = out.geometry.iloc[0]
    assert (kept.x, kept.y) == pytest.approx(
        (inside_px[1] * PIXEL_SIZE, inside_px[0] * PIXEL_SIZE)
    ), "micron coordinates must not be translated by a pixel offset"

    landed = crop_overlays.element_affine(out) @ np.array([kept.y, kept.x, 1.0])
    assert landed[0] == pytest.approx(inside_px[0] - CROP[0])
    assert landed[1] == pytest.approx(inside_px[1] - CROP[2])


def test_region_mapped_into_a_rotated_frame_is_not_transposed():
    """(y, x) affine vs shapely's (x, y): getting it backwards is a silent rotation.

    Map the crop region into a rotated element's frame, then map a point of that
    region back out by hand. It has to land on the original region.
    """
    affine = _similarity(theta_deg=30.0, scale=1.5, ty=20.0, tx=-40.0)
    region_own = crop_overlays._region_in_own_frame(_crop_polygon(), affine)

    ox, oy = region_own.exterior.coords[0]
    back = affine @ np.array([oy, ox, 1.0])                      # own (y, x) -> Xenium
    assert _crop_polygon().buffer(1e-6).contains(Point(back[1], back[0])), (
        "a corner of the mapped region did not map back onto the crop region — "
        "the (y, x) -> (x, y) conversion is transposed"
    )


def test_patch_overlay_outside_the_crop_is_dropped():
    """Patch overlays are data in their source image's pixels, not landmarks."""
    affine = _similarity(0.0, 1.0, ty=0.0, tx=0.0)
    affine[0, 2] = 5.0                                            # a nudge, still own-frame
    patches = gpd.GeoDataFrame({
        "geometry": [Point(250.0, 200.0), Point(9000.0, 9000.0)],
        "radius": 8.0,
        "cluster": np.array([1, 2], dtype=np.int64),
    })
    out = crop_overlays.crop_vector_overlay(
        _shapes(patches, affine=affine), "patch_tumour", CROP, _crop_polygon())
    assert len(out) == 1
    assert list(out["cluster"]) == [1], "the surviving patch kept its cluster column"


# ── disk round-trip ──────────────────────────────────────────────────────────

def test_cropped_overlays_survive_a_write_and_reopen(tmp_path):
    """In-memory placement is worth nothing if the transform does not persist.

    The export's whole promise is that the *written* dataset opens with its
    alignment intact, so assert that across an actual zarr round-trip rather
    than against the objects still in hand.
    """
    import spatialdata as sd

    affine = _similarity(theta_deg=12.0, scale=0.9, ty=25.0, tx=-15.0)
    img, _ = crop_overlays.crop_raster_overlay(
        _image(affine=affine), CROP, _safe_scale_factors)
    rois = crop_overlays.crop_vector_overlay(
        _shapes(gpd.GeoDataFrame({
            "geometry": [Polygon([(200, 150), (300, 150), (300, 250), (200, 250)])]})),
        "rois", CROP, _crop_polygon())
    lms = crop_overlays.crop_vector_overlay(
        _landmarks([(10, 10), (20, 20)]), "he_he_landmarks", CROP, _crop_polygon())

    before = crop_overlays.element_affine(img)

    out = tmp_path / "export.zarr"
    sd.SpatialData(images={"he_image": img},
                   shapes={"rois": rois, "he_he_landmarks": lms}).write(out)
    reopened = sd.read_zarr(out)

    assert np.allclose(crop_overlays.element_affine(reopened.images["he_image"]), before), \
        "the composed registration did not survive the write"
    assert (reopened.shapes["rois"].geometry.iloc[0].bounds[0]
            == pytest.approx(200 - CROP[2]))
    assert [(p.x, p.y) for p in reopened.shapes["he_he_landmarks"].geometry] == [(10, 10), (20, 20)]
