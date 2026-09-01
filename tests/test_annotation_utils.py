"""Reading the annotation layer: what becomes an annotation, and what does not.

The geometry these tests used to cover moved into the ``annot.*`` templates
(see ``tests/test_annotation_steps.py``, which asserts the bowtie repair on the
code that now runs). What is left on the python side is the boundary between
the napari layer and the recorded step: which shapes are inlined, with which
type, at what precision. Every one of those choices ends up as a literal in the
notebook, so it is worth pinning here rather than only through a rendered
template.
"""

import numpy as np
import pytest

from palms.tabs._helpers import annotation_polygons_preview
from palms.utils.annotation_utils import get_annotation_types


class FakeShapesLayer:
    """The two attributes the annotation readers touch on a napari Shapes layer."""

    def __init__(self, shapes, types):
        self.data = [np.asarray(s, dtype=np.float64) for s in shapes]
        self.properties = {"annotation_type": list(types)}


class FakeCtx:
    def __init__(self, layer, pixel_size=0.2125):
        self.annotation_layer = layer
        self.pixel_size = pixel_size


SQUARE = [(0, 0), (0, 10), (10, 10), (10, 0)]
OTHER = [(20, 20), (20, 30), (30, 30), (30, 20)]


def test_an_untyped_shape_is_not_an_annotation():
    """A shape drawn but not yet labelled is not an annotation of type "" — the
    Annotations tab has simply not been told what it is. Inlining it would put
    an unnamed group in the notebook's results."""
    ctx = FakeCtx(FakeShapesLayer([SQUARE, OTHER], ["tumour", "  "]))

    _blocks, params, _ = annotation_polygons_preview(ctx)

    assert params["types"] == ["tumour"]
    assert len(params["polygons"]) == 1


def test_types_shorter_than_shapes_does_not_shift_the_labels():
    """layer.properties can lag layer.data by an entry. Zipping them short would
    silently give the second shape the third one's label."""
    ctx = FakeCtx(FakeShapesLayer([SQUARE, OTHER], ["tumour"]))

    _blocks, params, _ = annotation_polygons_preview(ctx)

    assert params["types"] == ["tumour"]
    assert params["polygons"][0][0] == [0.0, 0.0]     # the square, not the other


def test_the_params_are_plain_literals():
    """They are rendered into the recorded cell, and Step validates them with
    ast.literal_eval(repr(v)) == v — a numpy array would not survive it."""
    ctx = FakeCtx(FakeShapesLayer([SQUARE], ["tumour"]))

    _blocks, params, _ = annotation_polygons_preview(ctx)

    import ast
    assert ast.literal_eval(repr(params)) == params
    assert isinstance(params["polygons"], list)
    assert isinstance(params["polygons"][0][0][0], float)
    assert isinstance(params["pixel_size"], float)


def test_coordinates_are_rounded_so_the_recorded_cell_stays_readable():
    ctx = FakeCtx(FakeShapesLayer([[(0.123456, 1.987654)] * 3], ["tumour"]))

    _blocks, params, _ = annotation_polygons_preview(ctx)

    assert params["polygons"][0][0] == [0.12, 1.99]


def test_the_pixel_size_travels_with_the_shapes():
    """The template scales the drawn pixels to microns, so the conversion is
    recorded with the coordinates it applies to rather than assumed."""
    ctx = FakeCtx(FakeShapesLayer([SQUARE], ["tumour"]), pixel_size=0.5)

    _blocks, params, _ = annotation_polygons_preview(ctx)

    assert params["pixel_size"] == 0.5


def test_nothing_drawn_yields_no_polygons():
    """What ``ctx.ensure_annotations`` checks before running an analysis over an
    empty region — and what a preview must answer for a tab built before the
    user has drawn anything."""
    for layer in (None, FakeShapesLayer([], [])):
        _blocks, params, _ = annotation_polygons_preview(FakeCtx(layer))
        assert params["polygons"] == []


def test_get_annotation_types_ignores_blanks():
    layer = FakeShapesLayer([[(0, 0)]] * 4, ["b", "", "a", "  "])

    assert get_annotation_types(layer) == ["a", "b"]


def test_get_annotation_types_of_no_layer():
    assert get_annotation_types(None) == []
