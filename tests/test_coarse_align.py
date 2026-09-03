"""Properties of `Coarse Align` — the tissue masks and the transform search.

Written against transforms and images built here, so the right answer is known
by construction rather than remembered from a dataset that cannot ship in the
repo. The real measurement is `scripts/score_coarse_align.py` against the two
datasets that carry ground truth; this is the gate that stops the failures it
found from coming back.

Every assertion below corresponds to something that was actually broken:

  * the fluorescence mask thresholded *within* tissue and kept a third of the
    section, which made the true transform score 0.216 where a wrong one scored
    0.493 — no search can survive an objective maximised in the wrong place;
  * rotation was estimated from a moment axis and searched only at 0/90/180/270
    degrees off it, so an error of 35 degrees was unreachable;
  * reflection was in neither the search nor the similarity fit, so a mirrored
    H&E could not be aligned at all;
  * a fixed refinement window narrower than the scale grid's own spacing left
    the answer unreachable even when the grid straddled it.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

cv2 = pytest.importorskip("cv2")

from palms.utils.affine_compare import angle_difference_deg, decompose, to_yx
from palms.utils.registration import (
    SCALE_BAND_FOR, build_alignment_fields, compute_coarse_affine,
    compute_landmark_affine, extract_tissue_mask, he_pixel_size_um,
    mask_area_scale, pick_level,
)

# Small fields and a coarse grid: these tests are about behaviour, not precision,
# and the search cost is quadratic in the things being trimmed here.
FAST = dict(rotation_step=5.0, n_scales=3, scale_band=0.25)


def _blobs(shape=(240, 200), n=45, seed=0, radius=(6, 18), support=None):
    """A synthetic tissue-like field: structure at the architecture scale."""
    rng = np.random.default_rng(seed)
    field = np.zeros(shape, dtype=np.float32)
    ys, xs = np.mgrid[0:shape[0], 0:shape[1]]
    for _ in range(n):
        cy, cx = rng.uniform(0, shape[0]), rng.uniform(0, shape[1])
        r = rng.uniform(*radius)
        field += np.exp(-((ys - cy) ** 2 + (xs - cx) ** 2) / (2 * r * r)).astype(np.float32)
    field /= field.max()
    if support is not None:
        field *= support
    return field


def _similarity_xy(scale, deg, tx, ty):
    t = np.radians(deg)
    return np.array([[scale * np.cos(t), -scale * np.sin(t), tx],
                     [scale * np.sin(t), scale * np.cos(t), ty],
                     [0.0, 0.0, 1.0]])


def _target_from(source, scale, deg, out_shape, mirror=False):
    """Render *source* into a target frame through a known transform.

    Returns ``(target, M_xy)`` where ``M_xy`` maps source pixels to target
    pixels — exactly what ``compute_coarse_affine`` has to recover.
    """
    sh, sw = source.shape
    M = _similarity_xy(scale, deg, 0.0, 0.0)
    if mirror:
        M = M @ np.array([[-1.0, 0.0, sw - 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    # Centre the source's centre on the target's centre.
    c = M @ np.array([sw / 2.0, sh / 2.0, 1.0])
    M[0, 2] += out_shape[1] / 2.0 - c[0]
    M[1, 2] += out_shape[0] / 2.0 - c[1]
    target = cv2.warpAffine(source, M[:2, :], (out_shape[1], out_shape[0]))
    return target, M


def _run(target, source, scale_prior, **kw):
    opts = dict(FAST, mirror=False)
    opts.update(kw)
    return compute_coarse_affine(target, source, scale_prior=scale_prior, **opts)


# ── the search ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("scale,deg", [(1.0, 37.0), (1.3, -112.0), (0.75, 4.0)])
def test_a_known_similarity_comes_back(scale, deg):
    source = _blobs(seed=1)
    target, truth = _target_from(source, scale, deg, (260, 300))

    result = _run(target, source, scale_prior=scale)

    want = decompose(to_yx(truth))
    assert result.scale == pytest.approx(want["scale"], rel=0.02)
    assert abs(angle_difference_deg(result.rotation_deg, want["rotation_deg"])) < 1.0
    assert not result.mirrored
    assert result.confident


def test_a_wrong_scale_prior_is_recovered_from():
    """The prior is a starting point, not the answer.

    A fixed refinement window narrower than the scale grid's spacing could not
    cross from one grid point to the next, so a prior that was merely *near* the
    truth pinned the result to itself.
    """
    source = _blobs(seed=2)
    target, truth = _target_from(source, 1.20, 22.0, (280, 280))

    result = _run(target, source, scale_prior=1.20 * 1.12)   # 12% out

    assert result.scale == pytest.approx(decompose(to_yx(truth))["scale"], rel=0.02)
    assert abs(angle_difference_deg(result.rotation_deg, 22.0)) < 1.0


def test_a_rotation_between_the_right_angles_is_found():
    """The failure that started this: a moment axis plus 0/90/180/270.

    On the pancreas reference the moment axis was 35 degrees out, and none of
    the four hypotheses could reach the truth. Any angle must be reachable.
    """
    source = _blobs(seed=3)
    target, _ = _target_from(source, 1.0, 53.0, (300, 300))

    result = _run(target, source, scale_prior=1.0)

    assert abs(angle_difference_deg(result.rotation_deg, 53.0)) < 1.0


def test_a_near_isotropic_outline_is_still_recovered():
    """A round section has no principal axis; its interior still fixes the angle."""
    shape = (240, 240)
    ys, xs = np.mgrid[0:shape[0], 0:shape[1]]
    disc = (((ys - 120) ** 2 + (xs - 120) ** 2) < 110 ** 2).astype(np.float32)
    source = _blobs(shape=shape, seed=4, support=disc)
    target, _ = _target_from(source, 1.0, -66.0, (280, 280))

    result = _run(target, source, scale_prior=1.0)

    assert abs(angle_difference_deg(result.rotation_deg, -66.0)) < 1.5


def test_a_mirrored_source_is_detected_and_factored_out():
    """A reflection is reported, and left for the caller's flip to supply.

    ``compute_landmark_affine`` fits a similarity, which has no reflection, so a
    mirror baked into the returned matrix would be discarded by the very next
    step of the workflow.
    """
    source = _blobs(seed=5)
    sh, sw = source.shape
    target, _ = _target_from(source, 1.1, 15.0, (300, 300), mirror=True)

    result = compute_coarse_affine(target, source, scale_prior=1.1,
                                   source_shape_yx=(sh, sw), mirror=True, **FAST)

    assert result.mirrored
    assert np.linalg.det(result.affine_3x3_yx[:2, :2]) > 0, "the reflection must be factored out"
    # Composed with the flip the caller ticks, it maps the source onto the target.
    flip = np.array([[1, 0, 0], [0, -1, sw - 1], [0, 0, 1]], dtype=float)
    combined = result.affine_3x3_yx @ flip
    assert np.linalg.det(combined[:2, :2]) < 0


def test_unrelated_images_are_reported_as_not_confident():
    """No match is an answer. A confident wrong transform is not."""
    result = _run(_blobs(seed=6), _blobs(seed=7, n=90, radius=(2, 5)), scale_prior=1.0)
    assert not result.confident
    assert "LOW CONFIDENCE" in result.summary()


# ── the masks ─────────────────────────────────────────────────────────────────

def test_triangle_keeps_both_tissue_populations_where_otsu_splits_them():
    """The pancreas failure in miniature.

    A fluorescence max-projection has one background mode and *several* tissue
    modes — bright acinar, dim fibrotic. Otsu cuts between the two brightest,
    which threw away two thirds of the section and left the true transform
    scoring below a wrong one.
    """
    image = np.zeros((200, 300), dtype=np.uint8)
    image[40:160, 20:150] = 60      # dim tissue
    image[40:160, 150:280] = 220    # bright tissue
    rng = np.random.default_rng(0)
    image = np.clip(image + rng.normal(0, 4, image.shape), 0, 255).astype(np.uint8)

    otsu = extract_tissue_mask(image, method="otsu", outline=True)
    triangle = extract_tissue_mask(image, method="triangle", outline=True)

    dim = (slice(60, 140), slice(40, 130))
    bright = (slice(60, 140), slice(170, 260))
    assert (otsu[bright] > 0).mean() > 0.9
    assert (otsu[dim] > 0).mean() < 0.1, "Otsu is expected to cut within the tissue"
    assert (triangle[dim] > 0).mean() > 0.9
    assert (triangle[bright] > 0).mean() > 0.9


def test_every_scale_source_has_a_band():
    """The tab indexes SCALE_BAND_FOR by the source it just decided on.

    A source with no entry is not a bad default, it is a KeyError inside the
    worker thread — Coarse Align would simply stop with no transform.
    """
    assert set(SCALE_BAND_FOR) >= {"metadata", "tissue-area", "unknown"}
    assert all(band >= 0 for band in SCALE_BAND_FOR.values())
    tab = Path(__file__).resolve().parent.parent / "src/palms/tabs/tab_he_registration.py"
    for source in re.findall(r'scale_source = "([a-z-]+)"', tab.read_text()):
        assert source in SCALE_BAND_FOR, f"tab_he_registration uses {source!r}"


def test_mask_area_scale_reads_the_downsample_factors():
    small = np.zeros((100, 100), np.uint8); small[25:75, 25:75] = 255      # area 2500
    big = np.zeros((100, 100), np.uint8); big[10:90, 10:90] = 255          # area 6400
    # sqrt(6400/2500) = 1.6, times the ratio of the two downsamples.
    assert mask_area_scale(big, small, 4.0, 2.0) == pytest.approx(1.6 * 2.0)
    with pytest.raises(ValueError):
        mask_area_scale(np.zeros((10, 10), np.uint8), small)


def test_pick_level_chooses_by_size_not_by_position():
    """`pyramid[-1]` is not a size — it is however many levels a file carries.

    The same pancreas H&E bottoms out at 860x466 read from its OME-TIFF and at
    1718x931 read back from the zarr cache, so the identical dataset was aligned
    at two different resolutions depending on whether the session was restored.
    """
    pyramid = [np.zeros((1718, 931, 3)), np.zeros((859, 465, 3)),
               np.zeros((429, 232, 3)), np.zeros((214, 116, 3))]
    assert pick_level(pyramid, min_long_side=384).shape[0] == 429
    assert pick_level(pyramid, min_long_side=1000).shape[0] == 1718


def test_he_pixel_size_um_is_read_when_declared_and_none_when_not(tmp_path):
    tifffile = pytest.importorskip("tifffile")
    image = np.zeros((32, 32, 3), dtype=np.uint8)

    declared = tmp_path / "declared.ome.tif"
    tifffile.imwrite(declared, image, photometric="rgb",
                     metadata={"PhysicalSizeX": 0.27377, "PhysicalSizeXUnit": "µm"})
    with tifffile.TiffFile(declared) as tif:
        assert he_pixel_size_um(tif) == pytest.approx(0.27377, rel=1e-4)

    plain = tmp_path / "plain.tif"
    tifffile.imwrite(plain, image, photometric="rgb")
    with tifffile.TiffFile(plain) as tif:
        assert he_pixel_size_um(tif) is None


def test_build_alignment_fields_does_not_blur():
    """The smoothing belongs to the search, at its own working resolution.

    Blurring here means blurring in thumbnail pixels, which makes the amount of
    smoothing an accident of which pyramid level the file happened to carry.
    """
    morph = np.zeros((2, 40, 40), dtype=np.uint16)
    morph[0, 20, 20] = 4000
    he = np.full((40, 40, 3), 255, dtype=np.uint8)
    he[20, 20, 1] = 0

    target, source = build_alignment_fields(morph, he)

    assert target[20, 20] == target.max() and target[20, 19] == 0
    assert source[20, 20] == pytest.approx(1.0) and source[20, 19] == pytest.approx(0.0)


# ── the frame rule ────────────────────────────────────────────────────────────

def test_landmarks_must_be_fitted_in_the_flipped_frame():
    """`fine @ flip` is only correct if `fine` was fitted on flipped points.

    The tab reads H&E landmarks in unflipped layer-data coordinates and composes
    the result as ``fine @ flip``, so without this the flip is applied twice.
    Both reference datasets have their flips off, which is why it stayed latent —
    and it has to be right for automatic mirror detection to be usable at all.
    """
    width = 500
    flip = np.array([[1, 0, 0], [0, -1, width - 1], [0, 0, 1]], dtype=float)
    he_unflipped = np.array([[10.0, 20.0], [300.0, 60.0], [120.0, 400.0]])
    he_flipped = (flip @ np.hstack([he_unflipped, np.ones((3, 1))]).T).T[:, :2]

    truth = to_yx(_similarity_xy(1.4, 25.0, 33.0, -12.0))
    xenium = (truth @ np.hstack([he_flipped, np.ones((3, 1))]).T).T[:, :2]

    right, _ = compute_landmark_affine(xenium, he_flipped)
    np.testing.assert_allclose(right @ flip, truth @ flip, atol=1e-6)

    wrong, _ = compute_landmark_affine(xenium, he_unflipped)
    assert not np.allclose(wrong @ flip, truth @ flip, atol=1.0), (
        "fitting on unflipped points and composing with the flip must not agree "
        "with the truth — that is the double-flip this rule exists to prevent")


def test_an_extra_mirror_is_the_same_as_toggling_the_horizontal_flip():
    """Why `_on_coarse_done` toggles `Flip horizontally` instead of setting it.

    The search runs on the image *as already flipped*, so a reflection it finds
    is one further mirror on top of the flips in force. Composed with them it is
    always the same thing: H on nothing is H, H on H is nothing, H on V is both,
    H on both is V. That identity needs the two flips to commute and each to be
    its own inverse, which they are — they negate different axes.
    """
    h, w = 400, 500
    V = np.array([[-1, 0, h - 1], [0, 1, 0], [0, 0, 1]], dtype=float)
    H = np.array([[1, 0, 0], [0, -1, w - 1], [0, 0, 1]], dtype=float)
    I = np.eye(3)

    np.testing.assert_allclose(H @ H, I, atol=1e-9)
    np.testing.assert_allclose(V @ V, I, atol=1e-9)
    np.testing.assert_allclose(H @ V, V @ H, atol=1e-9)

    # (flip_v, flip_h) before -> after, and the composed transform it must equal.
    for (fv, fh), expected in (((False, False), H), ((False, True), I),
                               ((True, False), H @ V), ((True, True), V)):
        current = (H if fh else I) @ (V if fv else I)
        toggled = (H if not fh else I) @ (V if fv else I)
        np.testing.assert_allclose(H @ current, expected, atol=1e-9)
        np.testing.assert_allclose(toggled, expected, atol=1e-9)
