"""Properties of the automatic fine registration — nuclei against nuclear masks.

Synthetic, in the idiom of `test_coarse_align.py`: the transform is built here so
the right answer is known by construction rather than remembered from a dataset
that cannot ship in the repo. The real measurement is
`scripts/score_nuclei_align.py` against the two datasets that carry an
independent transform; this is the gate that stops the failures found while
building it from coming back.

Each assertion corresponds to something that was measured, or to an invariant a
plausible edit would break:

  * ICP cannot solve this problem — from a seed one nucleus spacing out, nearest
    neighbour matching locks onto the wrong nuclei and converges, confidently, to
    the wrong answer (15.5 um -> 13.1 um on the pancreas). The annealed fit must
    keep recovering from a seed several spacings out;
  * the fit must never introduce a reflection, because a mirrored H&E is carried
    as an explicit flip and a reflection here would silently undo it;
  * `flip_matrix` is the single definition of the flip, and the frame rule says
    coarse and fine are both maps *from* the flipped frame;
  * a label that first appears in a later row block must not drop out of the
    centroid accumulators;
  * a nucleus on a detector tile boundary must be found once, not twice and not
    never;
  * unrelated point sets must come back not confident, with enrichment ~1.
"""

from __future__ import annotations

import numpy as np
import pytest

cv2 = pytest.importorskip("cv2")
pytest.importorskip("skimage")
pytest.importorskip("scipy")

from palms.utils.nuclei_registration import (
    MIN_ENRICHMENT, _weighted_similarity, detect_he_nuclei, finest_level,
    fit_nuclei_similarity, haematoxylin_od, nucleus_centroids, register_he_nuclei,
)
from palms.utils.registration import flip_matrix

PIXEL_SIZE = 0.2125          # um per Xenium pixel, as on the reference datasets
SPACING_UM = 6.0             # median nucleus separation, as on the pancreas


def _similarity_yx(scale, deg, ty, tx):
    """A similarity in napari (y, x), the convention every transform here uses."""
    t = np.radians(deg)
    m = np.eye(3)
    m[:2, :2] = scale * np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])
    m[:2, 2] = (ty, tx)
    return m


def _nuclei(n=4000, seed=0, spacing_um=SPACING_UM):
    """An irregular point cloud with a realistic nearest-neighbour spacing.

    Uniform-random, **not a jittered lattice**. A lattice is periodic, so it
    looks the same shifted by one spacing and no matcher can resolve which
    offset is right: on one the fit settles 1.6 spacings out and reports a
    plausible residual. Real nuclei are irregular, which is what makes the
    correspondence unique — and it is what the enrichment statistic's Poisson
    baseline assumes, so the fixture and the metric agree by construction.
    """
    rng = np.random.default_rng(seed)
    spacing_px = spacing_um / PIXEL_SIZE
    # Density from the median nearest-neighbour distance of a Poisson field.
    lam = np.log(2.0) / (np.pi * spacing_px ** 2)
    side = np.sqrt(n / lam)
    return rng.uniform(0, side, (n, 2))


def _apply(m, pts):
    return (np.asarray(m, float) @ np.hstack([pts, np.ones((len(pts), 1))]).T).T[:, :2]


def _nudge(matrix, shift_um, deg=0.0):
    """A transform displaced from *matrix* by roughly *shift_um* of translation."""
    return _similarity_yx(1.0, deg, shift_um / PIXEL_SIZE, shift_um / PIXEL_SIZE) @ matrix


# ── The fit ───────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("scale,deg", [(1.0, 0.0), (1.2891, -89.88), (2.3532, 180.0)])
def test_a_known_similarity_comes_back(scale, deg):
    truth = _similarity_yx(scale, deg, 300.0, -120.0)
    target = _nuclei(3000, seed=1)
    source = _apply(np.linalg.inv(truth), target)

    result = fit_nuclei_similarity(target, source, _nudge(truth, 8.0),
                                   pixel_size_um=PIXEL_SIZE)

    assert result.confident
    residual = np.linalg.norm(_apply(result.affine_3x3_yx, source) - target, axis=1)
    assert np.median(residual) * PIXEL_SIZE < 0.05


def test_it_recovers_from_a_seed_several_nucleus_spacings_out():
    """The property ICP does not have.

    A seed 15 um out with nuclei 6 um apart puts the nearest target to most
    source points at the *wrong* nucleus, and those wrong correspondences agree
    with each other. Measured on the pancreas: trimmed ICP drove its own
    matched-pair RMSE from 5.10 um to 1.06 um while leaving the transform 13.1 um
    wrong. The annealed fit has no such basin limit.
    """
    truth = _similarity_yx(1.2891, -89.88, 300.0, -120.0)
    target = _nuclei(4000, seed=2)
    source = _apply(np.linalg.inv(truth), target)

    result = fit_nuclei_similarity(target, source, _nudge(truth, 15.0),
                                   pixel_size_um=PIXEL_SIZE)

    residual = np.linalg.norm(_apply(result.affine_3x3_yx, source) - target, axis=1)
    assert np.median(residual) * PIXEL_SIZE < 0.1
    assert result.seed_shift_um > 5.0        # it really did move that far


def test_nearest_neighbour_matching_is_what_fails_here():
    """The seed the previous test recovers from does defeat ICP, as claimed.

    Without this, "it is not ICP" is an assertion about a rejected alternative
    with nothing holding it true.
    """
    truth = _similarity_yx(1.2891, -89.88, 300.0, -120.0)
    target = _nuclei(4000, seed=2)
    source = _apply(np.linalg.inv(truth), target)
    seed = _nudge(truth, 15.0)

    from scipy.spatial import cKDTree
    tree = cKDTree(target)
    matrix = seed.copy()
    for _ in range(30):                       # trimmed ICP, generously iterated
        moved = _apply(matrix, source)
        dist, idx = tree.query(moved, workers=-1)
        keep = dist < np.percentile(dist, 60)
        matrix = _weighted_similarity(source[keep], target[idx[keep]],
                                      np.ones(int(keep.sum())))

    icp = np.median(np.linalg.norm(_apply(matrix, source) - target, axis=1)) * PIXEL_SIZE
    soft = fit_nuclei_similarity(target, source, seed, pixel_size_um=PIXEL_SIZE)
    soft_residual = np.median(
        np.linalg.norm(_apply(soft.affine_3x3_yx, source) - target, axis=1)) * PIXEL_SIZE
    assert icp > 1.0                          # ICP is still a whole nucleus out
    assert soft_residual < icp / 10


def test_unrelated_point_sets_are_reported_as_not_confident():
    """Overlapping on purpose — the hard case, where the fit gets to try.

    Two clouds far apart is the easy null: nothing matches at any transform. Two
    clouds of the same density *in the same place* let the optimiser look for
    structure that is not there, which is what an H&E of a different section
    looks like once Coarse Align has put it roughly on top.
    """
    target = _nuclei(3000, seed=3)
    source = _nuclei(3000, seed=99)

    result = fit_nuclei_similarity(target, source, np.eye(3), pixel_size_um=PIXEL_SIZE)

    assert not result.confident
    assert result.enrichment < MIN_ENRICHMENT


def test_a_source_with_no_targets_at_all_stays_at_the_seed():
    """It must not collapse the cloud onto the origin.

    Every weight underflowing to zero used to leave each point's virtual target
    at (0, 0) — a division by an epsilon — and the least-squares answer to "put
    every point at the origin" is scale 0. Measured before the fix: scale 6.8e-28
    and an *enrichment of 26x*, reported confident, because a cloud collapsed to a
    point lands on top of some nucleus. The failure mode is the dangerous kind:
    an H&E of the wrong section returning a confident transform.
    """
    target = _nuclei(3000, seed=13)
    source = _nuclei(3000, seed=14) + 1e5     # nothing within any sigma

    result = fit_nuclei_similarity(target, source, np.eye(3), pixel_size_um=PIXEL_SIZE)

    assert np.allclose(result.affine_3x3_yx, np.eye(3), atol=1e-6)
    assert result.scale == pytest.approx(1.0)
    assert not result.confident


def test_enrichment_is_about_one_for_a_point_set_with_nothing_to_match():
    """It must be a *ratio to chance*, not a raw matched fraction.

    Nuclei 6 um apart put a target within 1 um of an arbitrary point about 9% of
    the time, and a denser section more often still, so a raw fraction is not
    comparable between datasets — 31% on the pancreas and 45% on crop_6 are the
    same quality of fit.
    """
    target = _nuclei(4000, seed=4)
    source = _nuclei(4000, seed=100)

    result = fit_nuclei_similarity(target, source, np.eye(3), pixel_size_um=PIXEL_SIZE)

    # 1.11-1.24 over three unrelated pairings here; 0.6-1.8 on the two real
    # datasets under a deliberately corrupted transform.
    assert 0.5 < result.enrichment < 2.0


def test_too_few_nuclei_is_refused_rather_than_fitted():
    target = _nuclei(3000, seed=5)
    with pytest.raises(ValueError, match="not enough nuclei"):
        fit_nuclei_similarity(target, target[:10], np.eye(3), pixel_size_um=PIXEL_SIZE)


# ── No reflection, ever ───────────────────────────────────────────────────────

def test_the_fit_never_introduces_a_reflection():
    """A mirrored H&E is carried as an explicit flip; see the frame rule.

    If this fit were allowed a reflection it would absorb the mirror the tab is
    separately composing, and the image would be flipped twice.
    """
    target = _nuclei(2000, seed=6)
    mirrored = target * np.array([1.0, -1.0])        # a reflected source

    matrix = _weighted_similarity(mirrored, target, np.ones(len(target)))

    assert np.linalg.det(matrix[:2, :2]) > 0

    result = fit_nuclei_similarity(target, mirrored, np.eye(3), pixel_size_um=PIXEL_SIZE)
    assert np.linalg.det(result.affine_3x3_yx[:2, :2]) > 0
    assert not result.confident                      # and it does not pretend to fit


def test_weighted_similarity_is_a_similarity_and_honours_its_weights():
    truth = _similarity_yx(1.7, 33.0, 12.0, -4.0)
    src = _nuclei(500, seed=7)
    dst = _apply(truth, src)
    dst[:50] += 5000.0                                # gross outliers
    weights = np.ones(len(src))
    weights[:50] = 1e-9

    matrix = _weighted_similarity(src, dst, weights)

    assert np.allclose(matrix, truth, atol=1e-6)
    block = matrix[:2, :2]
    assert np.allclose(block @ block.T, np.eye(2) * (block[0, 0] ** 2 + block[0, 1] ** 2))


# ── The flip is defined once ──────────────────────────────────────────────────

def test_flip_matrix_is_an_involution_and_composes_both_axes():
    shape = (400, 250)
    for flip_v, flip_h in [(True, False), (False, True), (True, True)]:
        matrix = flip_matrix(shape, flip_v, flip_h)
        assert np.allclose(matrix @ matrix, np.eye(3))
    assert np.allclose(flip_matrix(shape, False, False), np.eye(3))
    corner = _apply(flip_matrix(shape, True, True), np.array([[0.0, 0.0]]))
    assert np.allclose(corner, [[shape[0] - 1, shape[1] - 1]])


def test_the_flip_has_exactly_one_definition():
    """A source guard: four hand-written copies of this matrix is how it drifts.

    The flip is load-bearing — ``compute_coarse_affine`` reports a detected
    reflection by asking the caller to tick one — so a second copy that composed
    the axes in the other order would misplace a mirrored H&E with nothing to
    catch it.
    """
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent / "src" / "palms"
    offenders = []
    for path in list(root.rglob("*.py")) + list(
            (root.parent.parent / "scripts").glob("*.py")):
        if path.name == "registration.py":
            continue
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            # The signature of a hand-rolled flip: a 3x3 literal whose leading
            # entry is -1 and whose last column carries a `- 1`.
            if not isinstance(node, ast.List) or len(node.elts) != 3:
                continue
            rows = [e for e in node.elts if isinstance(e, ast.List) and len(e.elts) == 3]
            if len(rows) != 3:
                continue
            text = ast.unparse(node)
            if "-1" in text and "- 1" in text and text.count("0") >= 4:
                offenders.append(f"{path.name}: {text}")
    assert not offenders, (
        "hand-written flip matrices found; use registration.flip_matrix:\n  "
        + "\n  ".join(offenders))


def test_the_fit_must_be_given_the_flipped_points():
    """The frame rule, and why the reflection cannot live in the matrix.

    ``coarse`` and ``fine`` are both maps *from* the flipped H&E frame. So when
    an H&E needs a flip, the map from the flipped frame is a proper rotation and
    the map from the *unflipped* one is a reflection — which
    ``_weighted_similarity`` refuses by construction. Handing the fit unflipped
    points therefore does not merely give a worse answer, it gives one that
    cannot be right, and it says so.
    """
    shape = (4000.0, 3600.0)
    truth = _similarity_yx(1.0, 0.0, 500.0, 400.0)      # from the *flipped* frame
    target = _nuclei(3000, seed=8)
    flipped_source = _apply(np.linalg.inv(truth), target)
    flip = flip_matrix(shape, True, False)
    unflipped = _apply(flip, flipped_source)            # what the layer holds

    good = fit_nuclei_similarity(target, flipped_source, _nudge(truth, 4.0),
                                 pixel_size_um=PIXEL_SIZE, image_shape_yx=shape)
    residual = np.linalg.norm(_apply(good.affine_3x3_yx, flipped_source) - target, axis=1)
    assert good.confident
    assert np.median(residual) * PIXEL_SIZE < 0.05

    bad = fit_nuclei_similarity(target, unflipped, _nudge(truth, 4.0),
                                pixel_size_um=PIXEL_SIZE, image_shape_yx=shape)
    assert not bad.confident


# ── The point sets ────────────────────────────────────────────────────────────

def _label_raster(shape=(200, 260), n=400, seed=0, radius=3, min_sep=30.0):
    """A label image of round nuclei, plus their exact centroids.

    ``min_sep`` keeps the nuclei further apart than the detector's own
    ``_MIN_SEPARATION_UM``; packed tighter, ``peak_local_max`` merges neighbours
    and the fixture measures the packing rather than the detector.
    """
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    labels = np.zeros(shape, dtype=np.uint32)
    margin = radius + 2
    # Rejection-sample the centres first, then rasterise each into its own small
    # window: an mgrid over the whole image per nucleus is minutes at the field
    # size these tests need.
    candidates = rng.uniform([margin, margin],
                             [shape[0] - margin, shape[1] - margin], (40 * n, 2))
    keep = []
    tree = None
    for point in candidates:
        if len(keep) >= n:
            break
        if keep:
            if tree is None or len(keep) % 64 == 0:
                tree = cKDTree(np.array(keep))
            if tree.query(point)[0] < min_sep:
                continue
            if np.min(np.hypot(*(np.array(keep[-64:]) - point).T)) < min_sep:
                continue
        keep.append(tuple(point))
    centres = []
    for label, (cy, cx) in enumerate(keep, start=1):
        y0, x0 = int(cy) - radius - 1, int(cx) - radius - 1
        ys, xs = np.mgrid[y0:y0 + 2 * radius + 3, x0:x0 + 2 * radius + 3]
        disc = (ys - cy) ** 2 + (xs - cx) ** 2 <= radius ** 2
        labels[y0:y0 + 2 * radius + 3, x0:x0 + 2 * radius + 3][disc] = label
        centres.append((ys[disc].mean(), xs[disc].mean()))
    return labels, np.array(centres)


def _he_with_nuclei(shape=(300, 400), centres=(), radius=4.0):
    """A pink field with dark-purple discs where the nuclei are."""
    density = np.zeros(shape, dtype=np.float32)
    reach = int(np.ceil(4 * radius))
    for cy, cx in centres:
        y0, y1 = max(0, int(cy) - reach), min(shape[0], int(cy) + reach + 1)
        x0, x1 = max(0, int(cx) - reach), min(shape[1], int(cx) + reach + 1)
        if y0 >= y1 or x0 >= x1:
            continue
        ys, xs = np.mgrid[y0:y1, x0:x1]
        density[y0:y1, x0:x1] += np.exp(
            -((ys - cy) ** 2 + (xs - cx) ** 2) / (2 * radius ** 2)).astype(np.float32)
    density = np.clip(density, 0, 1)
    rgb = np.empty(shape + (3,), dtype=np.uint8)
    rgb[..., 0] = 235 - 120 * density        # haematoxylin darkens all channels,
    rgb[..., 1] = 170 - 140 * density        # green most
    rgb[..., 2] = 225 - 80 * density
    return rgb


def test_detect_he_nuclei_finds_the_nuclei_that_are_there():
    centres = np.array([(40.0, 60.0), (120.0, 200.0), (220.0, 330.0), (70.0, 300.0)])
    rgb = _he_with_nuclei(centres=centres)

    found = detect_he_nuclei(rgb, pixel_size_um=0.2735)

    assert len(found) == len(centres)
    order = np.lexsort((found[:, 1], found[:, 0]))
    ref = centres[np.lexsort((centres[:, 1], centres[:, 0]))]
    assert np.allclose(found[order][:, :2], ref, atol=0.6)


def test_a_nucleus_on_a_tile_boundary_is_found_once():
    """Tiles overlap and each keeps only its own interior, so neither duplicates
    a nucleus nor loses the one that straddles the seam."""
    centres = np.array([(64.0, 64.0), (64.0, 128.0), (128.0, 64.0), (30.0, 30.0)])
    rgb = _he_with_nuclei(shape=(200, 200), centres=centres)

    whole = detect_he_nuclei(rgb, pixel_size_um=0.2735, tile=4096)
    tiled = detect_he_nuclei(rgb, pixel_size_um=0.2735, tile=64)

    assert len(tiled) == len(whole) == len(centres)
    order_w = np.lexsort((whole[:, 1], whole[:, 0]))
    order_t = np.lexsort((tiled[:, 1], tiled[:, 0]))
    assert np.allclose(whole[order_w][:, :2], tiled[order_t][:, :2], atol=1e-6)


def test_the_detector_is_physical_not_pixel_based():
    """The same tissue at two samplings must give the same positions in microns.

    Coarse Align's original defect in miniature: a smoothing stated in pixels
    makes the amount of smoothing an accident of which level the file carries.
    """
    centres = np.array([(60.0, 80.0), (160.0, 240.0), (260.0, 360.0)])
    fine = _he_with_nuclei(shape=(320, 420), centres=centres, radius=6.0)
    coarse = np.ascontiguousarray(fine[::2, ::2])

    found_fine = detect_he_nuclei(fine, pixel_size_um=0.25)
    found_coarse = detect_he_nuclei(coarse, pixel_size_um=0.50, downsample=2.0)

    assert len(found_fine) == len(found_coarse) == len(centres)
    a = found_fine[np.lexsort((found_fine[:, 1], found_fine[:, 0]))][:, :2]
    b = found_coarse[np.lexsort((found_coarse[:, 1], found_coarse[:, 0]))][:, :2]
    assert np.abs(a - b).max() < 2.0


def test_haematoxylin_od_accepts_both_channel_orders():
    rgb = _he_with_nuclei(shape=(64, 64), centres=[(32.0, 32.0)])
    assert np.allclose(haematoxylin_od(rgb),
                       haematoxylin_od(np.transpose(rgb, (2, 0, 1))))


def test_background_only_tiles_are_skipped():
    white = np.full((256, 256, 3), 250, dtype=np.uint8)
    assert len(detect_he_nuclei(white, pixel_size_um=0.2735)) == 0


# ── End to end ────────────────────────────────────────────────────────────────

def test_register_he_nuclei_recovers_a_known_transform():
    """The orchestrator, from a label raster and an RGB image to a matrix."""
    labels, centres = _label_raster(shape=(2400, 2400), n=1600, seed=11, radius=8)
    truth = _similarity_yx(1.0, 0.0, 20.0, 15.0)
    he_centres = _apply(np.linalg.inv(truth), centres)
    shape = (int(he_centres[:, 0].max()) + 20, int(he_centres[:, 1].max()) + 20)
    rgb = _he_with_nuclei(shape=shape, centres=he_centres, radius=6.0)

    result = register_he_nuclei(labels, rgb, _nudge(truth, 1.5),
                                pixel_size_um=PIXEL_SIZE,
                                he_pixel_size_um=PIXEL_SIZE)

    assert result.confident
    assert result.seed_shift_um > 0.5
    residual = np.linalg.norm(
        _apply(result.affine_3x3_yx, he_centres) - centres, axis=1) * PIXEL_SIZE
    assert np.median(residual) < 0.3


def test_the_he_pixel_size_falls_back_to_the_seed_scale():
    """An H&E that declares no pixel size still gets a physical detector.

    The seed's scale is H&E px -> Xenium px, which is a measurement of exactly
    the quantity the detector needs, and a seed is always required.
    """
    labels, centres = _label_raster(shape=(2400, 2400), n=1600, seed=12, radius=8)
    truth = _similarity_yx(1.0, 0.0, 20.0, 15.0)
    he_centres = _apply(np.linalg.inv(truth), centres)
    shape = (int(he_centres[:, 0].max()) + 20, int(he_centres[:, 1].max()) + 20)
    rgb = _he_with_nuclei(shape=shape, centres=he_centres, radius=6.0)

    declared = register_he_nuclei(labels, rgb, _nudge(truth, 1.5),
                                  pixel_size_um=PIXEL_SIZE,
                                  he_pixel_size_um=PIXEL_SIZE)
    derived = register_he_nuclei(labels, rgb, _nudge(truth, 1.5),
                                 pixel_size_um=PIXEL_SIZE)      # scale is 1.0

    assert np.allclose(declared.affine_3x3_yx, derived.affine_3x3_yx)


def test_register_he_nuclei_applies_the_flip_to_its_detections():
    """The orchestrator's half of the frame rule, through the real code path.

    The H&E on disk is unflipped; the transform is composed as ``fine @ flip``.
    So ``register_he_nuclei`` must map its detections into the flipped frame
    before fitting — detecting on the flipped *array* instead would be a 1.2 GB
    copy for the identical answer, which is why it does it to the points.
    """
    labels, centres = _label_raster(shape=(2400, 2400), n=1600, seed=15, radius=8)
    truth = _similarity_yx(1.0, 0.0, 20.0, 15.0)        # from the flipped frame
    flipped_he = _apply(np.linalg.inv(truth), centres)
    shape = (float(int(flipped_he[:, 0].max()) + 40),
             float(int(flipped_he[:, 1].max()) + 40))
    unflipped_he = _apply(flip_matrix(shape, True, False), flipped_he)
    rgb = _he_with_nuclei(shape=(int(shape[0]), int(shape[1])),
                          centres=unflipped_he, radius=6.0)

    result = register_he_nuclei(labels, rgb, _nudge(truth, 1.5),
                                pixel_size_um=PIXEL_SIZE, he_pixel_size_um=PIXEL_SIZE,
                                he_shape_yx=shape, flip_v=True, flip_h=False)

    assert result.confident
    residual = np.linalg.norm(
        _apply(result.affine_3x3_yx, flipped_he) - centres, axis=1) * PIXEL_SIZE
    assert np.median(residual) < 0.3


def test_a_pyramid_is_reduced_to_its_finest_level():
    """The layer hands over a whole pyramid; the accuracy is in level 0.

    Reported from the GUI as "list indices must be integers or slices, not
    tuple": the tab tested ``isinstance(data, (list, tuple))``, and napari wraps
    a multiscale layer in its own ``MultiScaleData``, which is neither. That
    object answers ``.ndim`` and ``.shape`` for level 0, so the whole pyramid
    passed every guard and only failed on the first tile slice, deep in the
    detector. A ``Sequence`` test separates every pyramid container from every
    array — numpy and dask arrays are not Sequences — so the rule lives here
    rather than at one call site.
    """
    import collections.abc

    base = np.zeros((8, 8, 3), np.uint8)
    half = np.zeros((4, 4, 3), np.uint8)

    assert finest_level(base) is base                       # a plain array
    assert finest_level([base, half]) is base               # a list pyramid
    assert finest_level((base, half)) is base               # a tuple pyramid

    class _MultiScaleLike(collections.abc.Sequence):
        """Stands in for napari's MultiScaleData without importing napari."""

        def __init__(self, levels):
            self._levels = list(levels)

        def __getitem__(self, i):
            return self._levels[i]

        def __len__(self):
            return len(self._levels)

        @property
        def ndim(self):
            return self._levels[0].ndim

        @property
        def shape(self):
            return self._levels[0].shape

    assert finest_level(_MultiScaleLike([base, half])) is base

    # And the detector works when handed one, which is what actually broke.
    centres = np.array([(40.0, 60.0), (120.0, 200.0)])
    rgb = _he_with_nuclei(shape=(200, 260), centres=centres)
    pyramid = _MultiScaleLike([rgb, np.ascontiguousarray(rgb[::2, ::2])])
    assert len(detect_he_nuclei(pyramid, pixel_size_um=0.2735)) == len(centres)
