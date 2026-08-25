"""Where an overlay's placement comes from, and the mixture that must never happen.

The defect these pin: a cropped export (and any store built by `xenium-build-cache`,
or whose session node was deleted) has no `viewer_session`, so the H&E and ARMS
restorers found no affine and added the image at the origin. The element carries the
registration; nothing read it back.

The subtler half is *how* it is read back. The element stores ``fine @ flip`` while the
session stores ``fine`` and re-derives the flip. Taking the affine from one and the flip
from the other applies the flip twice — an image that sits somewhere wrong, with no
error anywhere.
"""

import numpy as np
import pytest

from xenium_viewer.utils.registration_seed import seed_registration


class _FakeSdata:
    """Just enough of a SpatialData for `_load_affine_from_sdata_element`."""

    def __init__(self, transforms):
        self._t = transforms

    def __contains__(self, name):
        return name in self._t

    def __getitem__(self, name):
        return self._t[name]


def _element_with(affine_yx):
    """A bare object carrying a spatialdata 'global' transformation."""
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import Affine, Identity, set_transformation

    el = Image2DModel.parse(np.zeros((3, 4, 4), dtype=np.uint8), dims=("c", "y", "x"))
    if affine_yx is None:
        set_transformation(el, Identity(), "global")
    else:
        set_transformation(
            el, Affine(affine_yx, input_axes=("y", "x"), output_axes=("y", "x")), "global")
    return el


_KEYS = dict(affine_key="affine_3x3", coarse_key="coarse_affine",
             flip_v_key="flip_v", flip_h_key="flip_h", shape_key="he_shape_yx")

_REG = np.array([[-2.3535, 0.0, 34396.65],
                 [0.0, -2.3535, 53689.95],
                 [0.0, 0.0, 1.0]])


def test_the_session_wins_when_it_has_an_affine():
    """An ordinary store is unaffected — this path must not change behaviour."""
    sdata = _FakeSdata({"he_image": _element_with(_REG)})
    session = {"affine_3x3": np.eye(3) * 2, "flip_v": True, "he_shape_yx": (100, 200)}

    reg = seed_registration(session, sdata, "he_image", **_KEYS,
                            element_shape_yx=(4, 4))

    assert reg.source == "session"
    np.testing.assert_allclose(reg.fine, np.eye(3) * 2)
    assert reg.flip_v is True, "a session flip must survive when the session has the affine"
    assert reg.shape_yx == (100, 200)


def test_the_element_places_itself_when_the_session_has_nothing():
    """The crop-export case: elements on disk, no viewer_session to read."""
    sdata = _FakeSdata({"he_image": _element_with(_REG)})

    reg = seed_registration({}, sdata, "he_image", **_KEYS, element_shape_yx=(2254, 16371))

    assert reg.source == "element"
    np.testing.assert_allclose(reg.fine, _REG)
    assert reg.shape_yx == (2254, 16371), "the element's own shape, not the session's"


def test_a_session_flip_is_never_paired_with_an_element_affine():
    """The double-flip. The element transform already contains the flip.

    A session that remembers `flip_v=True` alongside an element affine that was
    written as ``fine @ flip`` would otherwise yield ``flip @ fine @ flip``.
    """
    sdata = _FakeSdata({"he_image": _element_with(_REG)})
    session = {"flip_v": True, "flip_h": True, "he_shape_yx": (13690, 23092)}

    reg = seed_registration(session, sdata, "he_image", **_KEYS,
                            element_shape_yx=(2254, 16371))

    assert reg.source == "element"
    assert reg.flip_v is False and reg.flip_h is False, (
        "the element's transform already contains the flip; re-deriving one from the "
        "session applies it twice"
    )
    assert reg.shape_yx == (2254, 16371), (
        "a flip built against the session's shape would be off by the difference "
        "between the source raster and the exported one"
    )


def test_a_coarse_affine_alone_still_counts_as_the_session_having_one():
    """Coarse-only is a real placement — the element must not override it."""
    sdata = _FakeSdata({"he_image": _element_with(_REG)})
    coarse = np.eye(3) * 3

    reg = seed_registration({"coarse_affine": coarse, "flip_h": True},
                            sdata, "he_image", **_KEYS, element_shape_yx=(4, 4))

    assert reg.source == "session"
    assert reg.fine is None
    np.testing.assert_allclose(reg.coarse, coarse)
    assert reg.flip_h is True


def test_an_identity_element_is_not_a_registration():
    """`_load_affine_from_sdata_element` returns None for identity, and so must this.

    Otherwise every unregistered overlay would come back claiming a placement.
    """
    sdata = _FakeSdata({"he_image": _element_with(None)})

    reg = seed_registration({"flip_v": True}, sdata, "he_image", **_KEYS,
                            element_shape_yx=(4, 4))

    assert reg.fine is None and reg.coarse is None
    assert reg.source == "session"
    assert reg.flip_v is True, "with no affine anywhere the session's flip is all there is"


def test_a_missing_element_is_not_an_error():
    reg = seed_registration({}, _FakeSdata({}), "he_image", **_KEYS,
                            element_shape_yx=(4, 4))
    assert reg.fine is None and reg.source == "session"


def test_no_sdata_at_all_is_not_an_error():
    reg = seed_registration({}, None, "he_image", **_KEYS, element_shape_yx=(4, 4))
    assert reg.fine is None


@pytest.mark.parametrize("session", [None, {}])
def test_an_absent_session_is_the_same_as_an_empty_one(session):
    sdata = _FakeSdata({"he_image": _element_with(_REG)})
    reg = seed_registration(session, sdata, "he_image", **_KEYS, element_shape_yx=(4, 4))
    assert reg.source == "element"


def test_arms_has_no_coarse_key_and_that_is_not_a_crash():
    """ARMS calls this with `coarse_key=None` — it has no coarse alignment."""
    sdata = _FakeSdata({"arms_he_image": _element_with(_REG)})

    reg = seed_registration({}, sdata, "arms_he_image",
                            affine_key="affine_3x3", coarse_key=None,
                            flip_v_key="flip_v", flip_h_key="flip_h",
                            shape_key="he_shape_yx", element_shape_yx=(10, 10))

    assert reg.coarse is None
    np.testing.assert_allclose(reg.fine, _REG)
