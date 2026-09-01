"""Annotation geometry: the polygon repair must keep the whole shape.

Both annotation helpers used ``poly.buffer(0)`` to repair a self-intersecting
polygon. That does not repair it — on a bowtie it *deletes* a lobe, so virtual
cells were sampled from, and distances measured to, half the region the user
drew, with no warning. ``roi.polygons.tmpl`` already carried the rule and the
reason; these tests are what stops the idiom coming back on the python side.
"""

import numpy as np
import pytest

from palms.utils.annotation_utils import (
    compute_distance_to_annotation,
    get_annotation_types,
    sample_annotation_centroids,
)


class FakeShapesLayer:
    """The two attributes the annotation helpers read off a napari Shapes layer."""

    def __init__(self, shapes, types):
        self.data = [np.asarray(s, dtype=np.float64) for s in shapes]
        self.properties = {"annotation_type": list(types)}


# A bowtie in napari yx-pixel coordinates: two 50 µm² lobes meeting at (5, 5).
# buffer(0) keeps one of them; make_valid keeps both.
BOWTIE = [(0, 0), (10, 10), (10, 0), (0, 10)]
LOWER_LOBE_CENTRE = (5.0, 2.0)   # xy-µm, inside the lobe buffer(0) drops
UPPER_LOBE_CENTRE = (5.0, 8.0)   # xy-µm, inside the lobe buffer(0) keeps


@pytest.fixture
def bowtie_layer():
    return FakeShapesLayer([BOWTIE], ["tumour"])


def _contains(points, xy, tol=1.5):
    """Is any sampled point within `tol` µm of `xy`?"""
    if len(points) == 0:
        return False
    return bool((np.abs(points - np.asarray(xy)).max(axis=1) <= tol).any())


def test_sampling_covers_both_lobes_of_a_self_intersecting_polygon(bowtie_layer):
    pts = sample_annotation_centroids(bowtie_layer, "tumour", pixel_size=1.0, density_um2=1.0)

    assert _contains(pts, UPPER_LOBE_CENTRE)
    # The regression: buffer(0) returns a single 25 µm² triangle, so nothing was
    # ever sampled from this half of the annotation.
    assert _contains(pts, LOWER_LOBE_CENTRE)


def test_sampled_area_is_the_whole_bowtie_not_half_of_it(bowtie_layer):
    """One virtual cell per µm², so the count estimates the area: 50, not 25."""
    pts = sample_annotation_centroids(bowtie_layer, "tumour", pixel_size=1.0, density_um2=1.0)

    assert 40 <= len(pts) <= 60, f"sampled {len(pts)} points; a half-bowtie gives ~25"


def test_distance_is_measured_to_both_lobes(bowtie_layer):
    """A cell in the dropped lobe is *inside* the annotation, so it sits within a
    few µm of a boundary — not the ~5 µm+ it would read if that lobe were gone."""
    centroids = np.array([LOWER_LOBE_CENTRE, UPPER_LOBE_CENTRE], dtype=np.float64)

    d = compute_distance_to_annotation(centroids, bowtie_layer, "tumour", pixel_size=1.0)

    assert d[1] == pytest.approx(2.0, abs=0.5)   # the kept lobe, unchanged
    assert d[0] == pytest.approx(2.0, abs=0.5)   # the dropped lobe


def test_a_line_drawn_with_the_polygon_tool_is_ignored():
    """make_valid answers a collinear ring with a LineString, which has bounds and
    a boundary and would otherwise pass for a region. buffer(0) returned an empty
    polygon here, and being ignored is the behaviour worth keeping."""
    layer = FakeShapesLayer([[(0, 0), (5, 5), (10, 10)]], ["scratch"])

    assert len(sample_annotation_centroids(layer, "scratch", pixel_size=1.0)) == 0
    d = compute_distance_to_annotation(
        np.array([[1.0, 1.0]]), layer, "scratch", pixel_size=1.0)
    assert np.isnan(d).all()


def test_pixel_size_scales_the_geometry():
    """The helpers take yx-pixels and work in xy-microns."""
    layer = FakeShapesLayer([[(0, 0), (0, 10), (10, 10), (10, 0)]], ["a"])

    at_1 = sample_annotation_centroids(layer, "a", pixel_size=1.0, density_um2=1.0)
    at_2 = sample_annotation_centroids(layer, "a", pixel_size=2.0, density_um2=1.0)

    assert len(at_2) == pytest.approx(4 * len(at_1), rel=0.15)  # 2x linear → 4x area


def test_missing_type_returns_empty_and_nan():
    layer = FakeShapesLayer([[(0, 0), (0, 10), (10, 10), (10, 0)]], ["a"])

    assert len(sample_annotation_centroids(layer, "absent", pixel_size=1.0)) == 0
    d = compute_distance_to_annotation(
        np.array([[1.0, 1.0], [2.0, 2.0]]), layer, "absent", pixel_size=1.0)
    assert d.shape == (2,) and np.isnan(d).all()


def test_types_shorter_than_shapes_is_padded_not_zipped_short():
    """A shape drawn but not yet assigned a type must not silently take the next
    shape's label — layer.properties can lag layer.data by one entry."""
    layer = FakeShapesLayer(
        [[(0, 0), (0, 10), (10, 10), (10, 0)], [(20, 20), (20, 30), (30, 30), (30, 20)]],
        ["a"],  # the second shape has no type yet
    )

    pts = sample_annotation_centroids(layer, "a", pixel_size=1.0, density_um2=4.0)

    assert len(pts) > 0
    assert pts[:, 0].max() < 15, "the untyped second shape was sampled as type 'a'"


def test_get_annotation_types_ignores_blanks():
    layer = FakeShapesLayer([[(0, 0)]] * 4, ["b", "", "a", "  "])

    assert get_annotation_types(layer) == ["a", "b"]
