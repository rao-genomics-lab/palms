"""Unit tests for registration landmark math and landmark JSON I/O.

`compute_landmark_affine` fits a similarity transform; `save_landmarks`/
`load_landmarks` are a JSON round-trip. Both are pure (numpy/skimage/json).

Run with:  pytest tests/test_registration.py
"""
from __future__ import annotations

import numpy as np
import pytest

from palms.utils.registration import (
    as_rgb_yxc, compute_landmark_affine, extract_tissue_mask_he,
    load_landmarks, nuclear_density_he, rgb_pyramid_yxc, save_landmarks,
)


def _similarity_yx(pts_yx, scale, theta_deg, translation_yx):
    """Apply a known similarity (scaled rotation + translation) in (y, x)."""
    t = np.deg2rad(theta_deg)
    R = scale * np.array([[np.cos(t), -np.sin(t)],
                          [np.sin(t),  np.cos(t)]])
    return (R @ pts_yx.T).T + np.asarray(translation_yx, dtype=float)


def test_compute_landmark_affine_recovers_known_similarity():
    he = np.array([[0.0, 0.0], [0.0, 10.0], [10.0, 0.0], [5.0, 7.0]])
    xenium = _similarity_yx(he, scale=1.5, theta_deg=30.0, translation_yx=(3.0, -2.0))

    affine, residuals = compute_landmark_affine(xenium, he)

    assert affine.shape == (3, 3)
    assert residuals.shape == (len(he),)
    assert residuals.max() < 1e-6, "an exact similarity must fit with ~zero residual"

    # The returned affine should map H&E (y, x) points onto the Xenium points.
    he_homo = np.hstack([he, np.ones((len(he), 1))])
    recovered = (affine @ he_homo.T).T[:, :2]
    assert np.allclose(recovered, xenium, atol=1e-6)


def test_landmarks_json_roundtrip_full(tmp_path):
    xenium = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    he = np.array([[0.5, 0.5], [1.5, 2.5], [7.0, 8.0]])
    affine = np.eye(3)
    path = tmp_path / "landmarks.json"

    save_landmarks(path, xenium, he, affine=affine, he_filename="slide_he.ome.tif")
    data = load_landmarks(path)

    assert np.allclose(data["xenium_landmarks_yx"], xenium)
    assert np.allclose(data["he_landmarks_yx"], he)
    assert np.allclose(data["affine_3x3_yx"], affine)
    assert data["he_filename"] == "slide_he.ome.tif"


def test_landmarks_json_roundtrip_optional_fields_absent(tmp_path):
    xenium = np.array([[1.0, 2.0], [3.0, 4.0]])
    he = np.array([[0.0, 0.0], [1.0, 1.0]])
    path = tmp_path / "landmarks_min.json"

    save_landmarks(path, xenium, he)  # no affine, no filename
    data = load_landmarks(path)

    assert np.allclose(data["xenium_landmarks_yx"], xenium)
    assert np.allclose(data["he_landmarks_yx"], he)
    assert "affine_3x3_yx" not in data
    assert "he_filename" not in data


# ── The channel-layout contract ──────────────────────────────────────────────
#
# `he.load` binds one name, `he_pyramid`, from two blocks that read different
# things: a TIFF, which is (Y, X, C), and `sdata.images['he_image']`, which is an
# `Image2DModel` and therefore (C, Y, X). Every consumer assumes channel-last, so
# for a session restored from the cache Coarse Align either died inside
# `cv2.cvtColor` ("Invalid number of channels ... 'scn' is 1", because OpenCV
# builds a 3-D single-channel Mat when the trailing dim is not a channel count)
# or -- with a pixel size on record, which takes the other branch -- silently
# fitted a transform 16x out of scale. Measured on a real prostate section:
# 10.0826 against a truth of 0.6298.

def test_as_rgb_yxc_normalises_every_layout():
    rng = np.random.default_rng(0)
    yxc = (rng.random((16, 24, 3)) * 255).astype(np.uint8)
    cyx = np.ascontiguousarray(np.transpose(yxc, (2, 0, 1)))

    assert np.array_equal(as_rgb_yxc(cyx), yxc)
    assert as_rgb_yxc(yxc) is yxc               # already channel-last: untouched
    gray = yxc[..., 0]
    assert as_rgb_yxc(gray) is gray             # 2-D: nothing to do


def test_as_rgb_yxc_keeps_a_dask_array_lazy():
    """The restore path leaves the pyramid lazy so napari fetches only the tiles
    it draws. A helper that materialised it would pull a whole slide into RAM."""
    import dask.array as da

    cyx = da.zeros((3, 64, 96), chunks=(3, 32, 32), dtype="uint8")
    out = as_rgb_yxc(cyx)

    assert isinstance(out, da.Array)
    assert out.shape == (64, 96, 3)


def test_rgb_pyramid_yxc_normalises_every_level():
    levels = [np.zeros((3, 64, 96), np.uint8), np.zeros((3, 32, 48), np.uint8)]
    assert [level.shape for level in rgb_pyramid_yxc(levels)] == [(64, 96, 3), (32, 48, 3)]


def test_the_he_mask_and_density_do_not_depend_on_the_layout():
    """The property that was false, and the reason the bug had a silent half."""
    rng = np.random.default_rng(1)
    yxc = (rng.random((48, 64, 3)) * 255).astype(np.uint8)
    cyx = np.ascontiguousarray(np.transpose(yxc, (2, 0, 1)))

    assert np.array_equal(extract_tissue_mask_he(yxc), extract_tissue_mask_he(cyx))
    assert np.array_equal(nuclear_density_he(yxc), nuclear_density_he(cyx))


def test_nuclear_density_he_still_accepts_a_grayscale_image():
    gray = np.full((8, 8), 51, np.uint8)
    assert np.allclose(nuclear_density_he(gray), 1.0 - 51 / 255.0)


@pytest.mark.parametrize("bad", [
    np.zeros((5, 16, 24), np.uint8),   # neither a channel count nor an image
    np.zeros((16,), np.uint8),         # 1-D
])
def test_a_non_rgb_image_is_refused_by_name(bad):
    """Named here rather than left to OpenCV, whose message for this reports a
    channel count ('scn' is 1) that nothing in the call ever passed in."""
    with pytest.raises(ValueError, match="RGB"):
        extract_tissue_mask_he(bad)


def test_nothing_outside_registration_hand_rolls_the_channel_transpose():
    """One definition of "channel-last", in the idiom of ``flip_matrix``.

    Four modules had grown their own copy of ``shape[0] in (3, 4)`` followed by a
    transpose, and they did not agree: the two tabs applied it for the napari
    layer while ``he.load``'s ``from_store`` block did not, so the layer and the
    recorded step were looking at differently shaped arrays.

    ``detect_he_nuclei`` keeps its own test, and the exemption is named here with
    its reason: it needs the layout *before* slicing a tile out of a lazy array,
    which a helper returning a normalised array cannot give it.
    """
    import ast as _ast
    from pathlib import Path as _Path

    src = _Path(__file__).resolve().parents[1] / "src" / "palms"
    # The whole of registration.py: it is the module that owns the layout, and
    # `parse_rgb_image_for_store` (the reverse, for the (c, y, x) store) and
    # `load_multichannel_pyramid` (which locates a channel axis) belong there too.
    allowed_modules = {"utils/registration.py"}
    allowed = {("utils/nuclei_registration.py", "detect_he_nuclei")}

    def _is_channel_test(node):
        if not (isinstance(node, _ast.Compare) and len(node.ops) == 1
                and isinstance(node.ops[0], _ast.In)):
            return False
        right = node.comparators[0]
        if not isinstance(right, (_ast.Tuple, _ast.List)):
            return False
        values = [e.value for e in right.elts if isinstance(e, _ast.Constant)]
        return values == [3, 4] and "shape" in _ast.dump(node.left)

    offenders = []
    for path in sorted(src.rglob("*.py")):
        tree = _ast.parse(path.read_text())
        funcs = [n for n in _ast.walk(tree)
                 if isinstance(n, (_ast.FunctionDef, _ast.AsyncFunctionDef))]
        for node in _ast.walk(tree):
            if not _is_channel_test(node):
                continue
            owner = max((f for f in funcs
                         if any(sub is node for sub in _ast.walk(f))),
                        key=lambda f: f.lineno, default=None)
            key = (path.relative_to(src).as_posix(), owner.name if owner else "<module>")
            if key[0] not in allowed_modules and key not in allowed:
                offenders.append(f"{key[0]}:{node.lineno} in {key[1]}")

    assert not offenders, (
        "use registration.as_rgb_yxc / rgb_pyramid_yxc instead of a hand-rolled "
        "channel-first test: " + ", ".join(offenders)
    )
