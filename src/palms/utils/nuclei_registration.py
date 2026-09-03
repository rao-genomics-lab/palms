"""Automatic fine H&E registration, by matching nuclei to nuclei.

What the manual landmark step is actually doing is matching nuclei in the H&E
against Xenium's nuclear masks — that is what makes it accurate, and a landmark
is only a hand-placed sample of that correspondence. This module does it
directly, over every nucleus in the section rather than the six a user clicks:

1. ``nucleus_centroids`` — area centroids of ``labels/nucleus_labels``, which is
   Xenium's own segmentation of the DAPI channel. Nothing is inferred here; the
   nuclei are already delineated.
2. ``detect_he_nuclei`` — haematoxylin optical density from the H&E (a real
   colour deconvolution, ``skimage.color.rgb2hed``), smoothed to a nuclear scale
   and reduced to sub-pixel local maxima.
3. ``bracket_seed`` then ``fit_nuclei_similarity`` — the similarity that best
   puts one point set on the other, starting from ``compute_coarse_affine``'s
   transform.

The fit is two passes, and both are load-bearing
------------------------------------------------
``bracket_seed`` searches a grid in *scale and rotation* with the translation for
each hypothesis read off a coincidence peak; the anneal below then refines what
is left. Neither does the other's job. A scale error displaces each point in
proportion to its distance from the centre, so on a 10.6 x 6.3 mm section a 1.5%
error is under a micron in the middle and 127 um at the corners — and no sigma
sees both: small enough to resolve nuclei in the middle and the corners match at
random, large enough to reach the corners and the density is flat with no
gradient to follow. Conversely the grid is far too coarse to land on an answer.
Measured on the prostate section that exposed this, whose coarse seed was 1.5% of
scale out: the anneal alone reached enrichment 1.4x and said LOW CONFIDENCE; the
bracket alone landed 1.3 um away; together, 9.6x and confident. Raising the
anneal's starting sigma to 144 um instead — with the target set thinned so the
neighbour budget could not silently truncate it — did not help, and cost twice
the time. The blind spot is structural, not a matter of range.

Why the anneal is not ICP
-------------------------
It is annealed soft assignment (the rigid/similarity case of coherent point
drift, with the neighbour sum truncated), and plain ICP is not an option rather
than a stylistic alternative. Nuclei in the pancreas reference sit a median
6.1 um apart and Coarse Align lands 15 um out, so the nearest target to a source
point is usually *the wrong nucleus* — and those wrong correspondences are
mutually consistent, being a whole-field shift of about one nucleus spacing.
Measured: trimmed ICP from the coarse seed drives its own matched-pair RMSE from
5.10 um to 1.06 um while moving the transform from 15.5 um of disagreement to
13.1 um. It converges, confidently, to the wrong answer.

Soft assignment has no such basin limit. At large sigma each source point is
drawn toward a weighted mean of many targets, which is a smooth field with the
whole nuclear architecture in it — the same quantity Coarse Align correlates,
but expressed on points. Annealing sigma from 20 um to 1 um walks that
continuously down to per-nucleus correspondence: on the very same point sets that
left ICP at 13.1 um, it reaches 1.02 um. (Both figures are from the first
detections tried, at s1; the resolution paragraph below is what takes the fit
from there to 0.70 um.)

Resolution is the accuracy, not the point count
-----------------------------------------------
Both point sets must be built at **full resolution**. Measured on the pancreas
reference, holding the detector's parameters fixed *in microns* so only the
sampling changes:

    H&E detections   nucleus_labels   disagreement with 10x
    s1  (0.55 um/px)  s1 (0.43 um/px)        1.023 um
    s1, 2.2x as many  s1                     1.021 um
    s0  (0.27 um/px)  s1                     0.878 um
    s0                s0 (0.21 um/px)        0.682 um

Detecting 2.2x as many nuclei at s1 bought 0.002 um; the same detector run one
level finer bought 0.145 um, and the label centroids finer again bought 0.196 um.
A coarser pixel does not merely add noise that 10^5 points average away — it
biases each peak, and a bias does not average away. So the level is not tunable
here: use what the file has.

Reading the confidence
----------------------
``enrichment`` is the number to look at, and roughly 1x means the fit found
nothing. Measured on three datasets: 17.2x (pancreas), 22.0x (crop_6) and 9.6x
(a 10.6 x 6.3 mm prostate section), against 1.0-1.8x for every deliberately wrong
transform tried. When it comes back low the usual cause is the *seed*, not the
detections and not the section — check ``bracket_shift_um``, which says how much
of the move was correcting Coarse Align's scale and rotation.

What the numbers mean
---------------------
0.93 um is the disagreement between PALMS's *landmark* fit and 10x's shipped
``he_imagealignment.csv`` on the pancreas reference; this module reaches 0.70 um
with nothing placed by hand, from a coarse seed 15.5 um out, in 99 s. But all
three estimates agree with each other to about a micron in every pairing, while
the nuclei fit's own reproducibility is 0.02 um — eight independent random halves
of the detections agree with each other to 0.010-0.035 um. So the residual micron
is not this estimator's noise, and **10x's matrix is not a ruler fine enough to
measure it**: "beat 0.93 um against 10x" is a bar this passes, not evidence about
which of the two is closer to the section.

``scripts/score_nuclei_align.py`` reports the check that is evidence, and it is
not agreement with another estimate: fit on one half of the *section*, then score
every candidate transform by how well it puts the *other* half's nuclei onto the
nuclear masks — data that fit has never seen, scored the same way for all of
them. On both datasets the held-out automatic fit wins, including against the
reference it is being compared to:

                          held-out median residual, um
                        automatic   10x/landmark   coarse
    pancreas  lower       1.472        1.479        2.078
              upper       1.139        1.240        2.095
    crop_6    lower       0.982        1.299        2.109
              upper       1.012        1.217        2.066
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from palms.utils.registration import flip_matrix

#: Nuclear scale, in microns. The haematoxylin field is smoothed by
#: ``_SMOOTH_UM`` and peaks closer than ``_MIN_SEPARATION_UM`` are one nucleus.
#: Stated in microns rather than pixels so a 0.27 um/px H&E and a 0.50 um/px one
#: are smoothed by the same physical amount — the mistake Coarse Align made in
#: thumbnail pixels, where the smoothing became an accident of the file's levels.
_SMOOTH_UM = 0.88
_MIN_SEPARATION_UM = 1.64

#: Haematoxylin optical density below which a local maximum is not a nucleus.
#: Flat across the range that matters: 0.03, 0.04 and 0.05 give 162909, 126609
#: and 89456 detections on the pancreas reference and move the fitted transform
#: by 0.005 um. The detector is limited by where it puts a peak, not how many it
#: finds, so this is set low enough to keep faint nuclei and no lower.
_HAEM_OD_FLOOR = 0.04

#: Tiles whose darkest pixel is above this are slide background, and are skipped
#: before the deconvolution — the expensive step — rather than after it.
_BACKGROUND_U8 = 200

#: How far the bracket searches around the seed, in fractional scale and degrees,
#: and how finely. The residual after Coarse Align is a *similarity* error, and
#: its scale and rotation parts are what the anneal cannot see — see
#: ``bracket_seed``. Measured seed errors: 0.05% and 0.24% of scale on the two
#: reference datasets, 1.49% and 0.52 deg on the prostate section that made this
#: necessary. The band is set to twice the worst of those.
BRACKET_SCALE_BAND = 0.03
BRACKET_ROTATION_BAND_DEG = 1.0
_BRACKET_SCALE_STEP = 0.0025
_BRACKET_ROTATION_STEP = 0.1
#: Source points the bracket scores each hypothesis on. It is looking for a
#: displacement field, not a per-nucleus match, so a few thousand is plenty and
#: the cost is linear in this — 4,000 points took 147 s over the grid, 2,000 take
#: a quarter of that at the same answer.
_BRACKET_SAMPLE = 2000
#: Radius the bracket's offset histogram covers, in microns. It must exceed the
#: seed's translation error — 24.8 um on the prostate section — and the cost goes
#: as its square.
_BRACKET_RADIUS_UM = 60.0
_BRACKET_BIN_UM = 1.5
#: Fraction of the half-diagonal at which a scale or rotation error produces its
#: *typical* displacement. Measured from the width of the peak itself: on both
#: reference sections it falls to half its height 0.15% of scale from the optimum,
#: and 6.4 um (one nuclear spacing) divided by 0.0015 puts the typical radius at
#: 0.7 of the half-diagonal.
_TYPICAL_RADIUS_FRACTION = 0.7

#: Annealing schedule for the soft assignment, in microns. The first value is how
#: far the *bracketed* seed may be wrong, not the raw one: the bracket is what
#: crosses a large similarity error, and 20 um then covers what is left (1.3 um
#: on the prostate section). The last is the scale at which a source point sees
#: one nucleus.
SIGMA_START_UM = 20.0
SIGMA_END_UM = 1.0
_SIGMA_STEPS = 10
_ITERS_PER_SIGMA = 6

#: Weight of the uniform outlier component, as a fraction of the truncated
#: neighbour sum evaluated at 3 sigma. A detection with no nucleus near it — H&E
#: that extends past the Xenium field, a fold, a pigment granule — contributes
#: almost nothing rather than dragging the fit toward whatever is nearest.
_OUTLIER = 0.3

#: A match is a source point whose nearest target is within this, in microns.
#: Reported, never used to select: the fit itself never thresholds a distance.
MATCH_RADIUS_UM = 1.0

#: Below these the fit is not worth presenting as an answer. ``enrichment`` is
#: how many times more source points land within ``MATCH_RADIUS_UM`` of a target
#: than a Poisson field of the same target density would put there by chance, so
#: it is ~1.0 for two unrelated point sets whatever their density — which is what
#: makes it comparable across datasets, where the raw matched fraction is not
#: (31% on the pancreas and 45% on crop_6 are the same quality of fit at
#: different nuclear densities).
#:
#: The floor is set from the gap, which is wide. Measured on the two reference
#: datasets: 17.2x and 22.0x after fitting, against 1.6x and 1.8x at the coarse
#: seed the fit started from, 1.4x/1.4x for the fitted transform rotated 5 deg,
#: and 1.0x/0.8x for a fit seeded 640 um out, which does not converge. Nothing
#: observed sits between 1.8 and 17. Like Coarse Align's floors this labels
#: rather than refuses — the transform is applied either way.
MIN_ENRICHMENT = 2.5


@dataclass
class NucleiAlignResult:
    """What the nuclei fit found, and how much to believe it."""

    affine_3x3_yx: np.ndarray
    n_target: int
    n_source: int
    n_matched: int
    matched_fraction: float
    median_residual_um: float
    enrichment: float
    #: How far the fit moved the H&E from the seed, in microns, averaged over the
    #: image. A fit that barely moved found nothing to improve; one that moved
    #: much further than the seed's own error found something else.
    seed_shift_um: float
    #: How far the (scale, rotation) bracket moved the seed before the anneal ran,
    #: and how strong the coincidence it settled on was. A large bracket shift
    #: means Coarse Align's scale or rotation was materially out — 60 um on the
    #: prostate section, against 0 to a few microns where it was already right.
    bracket_shift_um: float
    bracket_peak: float
    scale: float
    rotation_deg: float
    confident: bool = True
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"Nuclei matched: {self.n_matched:,} of {self.n_source:,} H&E "
            f"detections within {MATCH_RADIUS_UM:g} um "
            f"({self.matched_fraction * 100:.0f}%), against "
            f"{self.n_target:,} nucleus masks",
            f"Enrichment over chance: {self.enrichment:.1f}x",
            f"Median residual: {self.median_residual_um:.2f} um",
            f"Moved from the seed by: {self.seed_shift_um:.1f} um "
            f"({self.bracket_shift_um:.1f} um of it by the scale/rotation bracket)",
            f"Scale: {self.scale:.4f}",
            f"Rotation: {self.rotation_deg:.2f} deg",
        ]
        if not self.confident:
            lines.append(
                "LOW CONFIDENCE — the H&E nuclei did not settle onto the nuclear "
                "masks. Check the overlay; the coarse alignment may be too far "
                "out to refine, or the two images may not be the same section.")
        lines.extend(self.notes)
        return "\n".join(lines)


# ─── Point sets ───────────────────────────────────────────────────────────────

def finest_level(image):
    """The full-resolution level of anything pyramid-like, or *image* itself.

    A pyramid is a ``Sequence``; numpy and dask arrays are not, so this separates
    them without importing napari or guessing from ``.ndim``. It exists because
    an ``isinstance(x, (list, tuple))`` test in the tab silently let napari's own
    ``MultiScaleData`` through as if it were a single array — it answers ``.ndim``
    and ``.shape`` for level 0, so the whole pyramid travelled all the way into
    the detector and failed there with "list indices must be integers or slices,
    not tuple". The rule belongs here, where every caller gets it.
    """
    import collections.abc

    if isinstance(image, collections.abc.Sequence) and len(image):
        return image[0]
    return image


def nucleus_centroids(labels, downsample: float = 1.0,
                      block_rows: int = 1024) -> np.ndarray:
    """Area centroids of every label in a 2-D label raster, as (N, 2) in (y, x).

    ``downsample`` is full_res / this_level, so the result is always in
    full-resolution pixels whichever level was passed.

    Accumulated with ``np.bincount`` over row blocks rather than
    ``scipy.ndimage.center_of_mass``: the full-resolution raster of the pancreas
    reference is 13770x34155 uint32 (1.9 GB) and holds 141,810 labels, which this
    reduces in 17 s without ever holding more than a block.
    """
    labels = finest_level(labels)
    h, w = int(labels.shape[-2]), int(labels.shape[-1])
    n = None
    count = sum_y = sum_x = None
    cols = np.arange(w, dtype=np.float64)
    for y0 in range(0, h, block_rows):
        block = np.asarray(labels[y0:y0 + block_rows])
        rows = block.shape[0]
        flat = block.ravel()
        if n is None:
            # One pass to size the accumulators, so a label id larger than the
            # first block's maximum cannot silently drop out of the bincounts.
            n = int(np.asarray(labels[::max(1, h // 64)]).max()) + 1
            n = max(n, int(flat.max()) + 1)
            count = np.zeros(n, dtype=np.int64)
            sum_y = np.zeros(n, dtype=np.float64)
            sum_x = np.zeros(n, dtype=np.float64)
        if int(flat.max()) >= n:
            grow = int(flat.max()) + 1
            count = np.pad(count, (0, grow - n))
            sum_y = np.pad(sum_y, (0, grow - n))
            sum_x = np.pad(sum_x, (0, grow - n))
            n = grow
        count += np.bincount(flat, minlength=n)
        sum_y += np.bincount(
            flat, weights=np.repeat(np.arange(y0, y0 + rows, dtype=np.float64), w),
            minlength=n)
        sum_x += np.bincount(flat, weights=np.tile(cols, rows), minlength=n)
    keep = np.zeros(n, dtype=bool)
    keep[1:] = count[1:] > 0          # label 0 is background
    return np.column_stack([sum_y[keep] / count[keep],
                            sum_x[keep] / count[keep]]) * float(downsample)


def haematoxylin_od(rgb) -> np.ndarray:
    """Haematoxylin optical density of an RGB tile, as float32.

    A real colour deconvolution (``skimage.color.rgb2hed``), not the ``1 - green``
    proxy ``nuclear_density_he`` uses. That proxy is right for Coarse Align, which
    reduces the image to a few hundred pixels across and only needs tissue
    architecture; here the peak of this field *is* the estimate, and eosin
    bleeding into it moves that peak.
    """
    from skimage.color import rgb2hed

    arr = np.asarray(rgb)
    if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
        arr = np.transpose(arr, (1, 2, 0))
    arr = arr[..., :3]
    if arr.dtype != np.float32 and arr.dtype != np.float64:
        arr = arr.astype(np.float32) / 255.0
    return rgb2hed(arr)[..., 0].astype(np.float32)


def detect_he_nuclei(image, pixel_size_um: float, downsample: float = 1.0,
                     od_floor: float = _HAEM_OD_FLOOR, tile: int = 2048) -> np.ndarray:
    """Sub-pixel nucleus positions in an H&E, as (N, 3): y, x, haematoxylin OD.

    ``image`` is (Y, X, C) or (C, Y, X) and may be a dask array — it is read one
    tile at a time. ``pixel_size_um`` is *this level's* pixel size, which is what
    makes the smoothing and the minimum separation physical; ``downsample`` is
    full_res / this_level, so the returned points are in full-resolution pixels.

    Tiles overlap by the smoothing support and each keeps only the peaks whose
    refined position lands in its own interior, so a nucleus on a tile boundary
    is found once, by exactly one tile, and with its neighbourhood intact.
    """
    import cv2
    from skimage.feature import peak_local_max

    arr = finest_level(image)
    channel_first = (getattr(arr, "ndim", 0) == 3 and arr.shape[0] in (3, 4)
                     and arr.shape[-1] not in (3, 4))
    h, w = (arr.shape[1], arr.shape[2]) if channel_first else (arr.shape[0], arr.shape[1])

    sigma_px = max(_SMOOTH_UM / float(pixel_size_um), 0.6)
    min_sep = max(round(_MIN_SEPARATION_UM / float(pixel_size_um)), 1)
    pad = int(np.ceil(4 * sigma_px)) + min_sep

    found = []
    for y0 in range(0, h, tile):
        for x0 in range(0, w, tile):
            ya, yb = max(0, y0 - pad), min(h, y0 + tile + pad)
            xa, xb = max(0, x0 - pad), min(w, x0 + tile + pad)
            block = (arr[:, ya:yb, xa:xb] if channel_first
                     else arr[ya:yb, xa:xb])
            block = np.asarray(block)
            if channel_first:
                block = np.transpose(block, (1, 2, 0))
            if block.size == 0 or block[..., :3].min() > _BACKGROUND_U8:
                continue
            smooth = cv2.GaussianBlur(haematoxylin_od(block), (0, 0), sigma_px)
            peaks = peak_local_max(smooth, min_distance=min_sep,
                                   threshold_abs=od_floor)
            if len(peaks) == 0:
                continue
            yy, xx = peaks[:, 0], peaks[:, 1]
            inner = ((yy > 0) & (yy < smooth.shape[0] - 1)
                     & (xx > 0) & (xx < smooth.shape[1] - 1))
            yy, xx = yy[inner], xx[inner]
            if len(yy) == 0:
                continue
            centre = smooth[yy, xx]
            # Sub-pixel peak, from the quadratic through the three samples across
            # it in each axis. Without it every position is quantised to the pixel
            # grid, which at 0.27 um/px is a third of the accuracy on offer.
            dy = _parabolic(smooth[yy - 1, xx], centre, smooth[yy + 1, xx])
            dx = _parabolic(smooth[yy, xx - 1], centre, smooth[yy, xx + 1])
            py, px = yy + dy + ya, xx + dx + xa
            keep = ((py >= y0) & (py < min(h, y0 + tile))
                    & (px >= x0) & (px < min(w, x0 + tile)))
            if keep.any():
                found.append(np.column_stack([py[keep], px[keep], centre[keep]]))
    if not found:
        return np.empty((0, 3), dtype=np.float64)
    out = np.vstack(found)
    out[:, :2] *= float(downsample)
    return out


def _parabolic(left, centre, right):
    """Offset of the vertex of the parabola through three equally spaced samples."""
    denom = left - 2.0 * centre + right
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = 0.5 * (left - right) / denom
    return np.clip(np.nan_to_num(delta), -1.0, 1.0)


# ─── The fit ──────────────────────────────────────────────────────────────────

def _apply(m_yx, pts_yx):
    pts = np.asarray(pts_yx, dtype=np.float64)
    return (np.asarray(m_yx, dtype=np.float64)
            @ np.hstack([pts, np.ones((len(pts), 1))]).T).T[:, :2]


def _weighted_similarity(src, dst, weights):
    """The weighted-least-squares similarity mapping *src* onto *dst*.

    Umeyama, and deliberately **rotation only** — the reflection branch flips a
    singular value rather than accepting ``det(R) < 0``. A mirrored H&E is
    carried as an explicit flip by Coarse Align precisely because a similarity
    cannot hold a reflection, and this fit runs in that already-flipped frame, so
    letting one in here would silently undo the flip the tab is composing.
    """
    w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    mu_s = (w[:, None] * src).sum(0)
    mu_d = (w[:, None] * dst).sum(0)
    s_c, d_c = src - mu_s, dst - mu_d
    cov = (w[:, None] * d_c).T @ s_c
    u, sv, vt = np.linalg.svd(cov)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u = u.copy()
        u[:, -1] *= -1
        rot = u @ vt
        sv = sv.copy()
        sv[-1] *= -1
    var = (w[:, None] * s_c ** 2).sum()
    scale = sv.sum() / var if var > 0 else 1.0
    out = np.eye(3)
    out[:2, :2] = scale * rot
    out[:2, 2] = mu_d - scale * rot @ mu_s
    return out


def _neighbours_for(sigma_px, spacing_px):
    """How many nearest targets a source point needs at this sigma.

    Enough to cover 3 sigma: for a Poisson field whose median nearest-neighbour
    distance is *spacing*, the expected count inside radius r is
    ``ln(2) (r/spacing)^2``. Truncating tighter biases the weighted mean toward
    the source point; carrying more is wasted work, and at the end of the anneal
    it is 10x the work.
    """
    k = int(np.ceil(9.0 * np.log(2.0) * (sigma_px / max(spacing_px, 1e-6)) ** 2))
    return int(np.clip(k, 6, 64))


def bracket_seed(tree, source_yx, seed_3x3, pixel_size_um: float, image_shape_yx,
                 spacing_px: float = 1.0,
                 scale_band: float = BRACKET_SCALE_BAND,
                 rotation_band_deg: float = BRACKET_ROTATION_BAND_DEG,
                 sample: int = _BRACKET_SAMPLE):
    """Search (scale, rotation) around a seed, translation free per hypothesis.

    Returns ``(matrix, peak)`` — the best seed found, and how many times more
    nuclei coincide at its best translation than chance would put there.

    **Why the anneal needs this.** Soft assignment pulls each point toward the
    local nuclear density, which is a translation-like force: it fixes an error
    that is roughly the same vector everywhere. A *scale* error is not — it
    displaces each point in proportion to its distance from the centre, so on a
    10.6 x 6.3 mm section a 1.5% error is under a micron at the middle and 127 um
    at the corners. There is no sigma that sees both: small enough to resolve
    nuclei in the middle and the corners are matched at random; large enough to
    reach the corners and the density is flat, so there is no gradient to follow.
    Measured on that section — the prostate dataset that exposed this — the anneal
    moved 15.9 um out of the 60 um it needed and reported LOW CONFIDENCE, and
    raising the starting sigma to 144 um (with the target set thinned so the
    neighbour budget could not silently truncate it) did not help: it made it
    worse. The blind spot is structural, not a matter of range.

    A grid has no blind spot, and scale and rotation are only two dimensions.
    Translation is not searched — it is read off the peak of the offset
    histogram, exactly as ``compute_coarse_affine`` takes each hypothesis's
    translation from phase correlation.

    Scoring on the *peak* rather than on a match count is what makes the grid
    work at all: a hypothesis with the right scale but the wrong translation
    matches nothing, so a match count is flat across the grid until translation
    happens to be right. Measured: holding the centre fixed and scoring matches,
    the grid found +1.00% and stopped 28.4 um from the answer; solving the
    translation per hypothesis it found +1.50% and stopped 1.3 um away.
    """
    import cv2

    source = np.asarray(source_yx, dtype=np.float64)
    if len(source) > sample:
        step = max(1, len(source) // sample)
        source = source[::step][:sample]
    centre = np.array([image_shape_yx[0] / 2.0, image_shape_yx[1] / 2.0])
    seed = np.asarray(seed_3x3, dtype=np.float64)
    radius_px = _BRACKET_RADIUS_UM / pixel_size_um
    bins = int(2 * _BRACKET_RADIUS_UM / _BRACKET_BIN_UM)
    span = [[-_BRACKET_RADIUS_UM, _BRACKET_RADIUS_UM]] * 2

    # How far the grid can miss by, and therefore how much the peak has to be
    # smeared for a neighbouring grid point to still see it.
    #
    # This is the one thing that has to scale with the section, and it is why the
    # grid step can stay fixed. The peak is *narrow*: measured on both reference
    # datasets it falls to half height 0.15% of scale from the optimum, which is
    # finer than the 0.25% step — the search found it by luck, and on a section
    # twice as long the peak would be 0.07% wide and every grid point would sit in
    # background with the argmax deciding on noise. Smoothing the offset histogram
    # by exactly the displacement half a step produces makes the peak as wide as
    # the grid is coarse, at any section size. It costs nothing in accuracy: the
    # bracket only has to land inside the anneal's capture range, not on the
    # answer.
    half_diagonal = 0.5 * float(np.hypot(*image_shape_yx)) * pixel_size_um
    miss_um = (0.5 * _BRACKET_SCALE_STEP * _TYPICAL_RADIUS_FRACTION * half_diagonal)
    smooth_bins = float(np.clip(miss_um / _BRACKET_BIN_UM, 1.0, 12.0))

    def candidate(d_scale, d_deg):
        t = np.radians(d_deg)
        a = np.eye(3)
        a[:2, :2] = (1 + d_scale) * np.array([[np.cos(t), np.sin(t)],
                                              [-np.sin(t), np.cos(t)]])
        a[:2, 2] = centre - a[:2, :2] @ centre
        return seed @ a

    # A fixed-k query rather than query_ball_point: the latter returns a list of
    # lists, and the Python loop to turn that into offsets is pure overhead. k is
    # sized from the target density so it covers the radius, and anything beyond
    # it is masked out, so the answer is unchanged — measured identical transforms
    # on both reference datasets, at 30 s for the 525 hypotheses against 47 s.
    neighbours = int(np.clip(2 * np.log(2.0) * (radius_px / max(spacing_px, 1e-6)) ** 2,
                             16, 256))

    def score(matrix):
        moved = _apply(matrix, source)
        dist, idx = tree.query(moved, k=neighbours, workers=-1)
        inside = dist < radius_px
        if not inside.any():
            return -1.0, (0.0, 0.0)
        offsets = (tree.data[idx] - moved[:, None, :])[inside] * pixel_size_um
        hist, _, _ = np.histogram2d(offsets[:, 0], offsets[:, 1],
                                    bins=bins, range=span)
        hist = cv2.GaussianBlur(hist.astype(np.float32), (0, 0), smooth_bins)
        chance = len(offsets) * _BRACKET_BIN_UM ** 2 / (np.pi * _BRACKET_RADIUS_UM ** 2)
        peak = np.unravel_index(int(np.argmax(hist)), hist.shape)
        shift = (-_BRACKET_RADIUS_UM + (peak[0] + 0.5) * _BRACKET_BIN_UM,
                 -_BRACKET_RADIUS_UM + (peak[1] + 0.5) * _BRACKET_BIN_UM)
        return float(hist[peak] / max(chance, 1e-12)), shift

    scales = np.arange(-scale_band, scale_band + 1e-9, _BRACKET_SCALE_STEP)
    degrees = np.arange(-rotation_band_deg, rotation_band_deg + 1e-9,
                        _BRACKET_ROTATION_STEP)
    best = (-1.0, 0.0, 0.0, (0.0, 0.0))
    for d_scale in scales:
        for d_deg in degrees:
            peak, shift = score(candidate(d_scale, d_deg))
            if peak > best[0]:
                best = (peak, float(d_scale), float(d_deg), shift)

    peak, d_scale, d_deg, shift = best
    matrix = candidate(d_scale, d_deg)
    matrix = matrix.copy()
    matrix[:2, 2] += np.asarray(shift) / pixel_size_um
    return matrix, peak


def fit_nuclei_similarity(target_yx, source_yx, seed_3x3, pixel_size_um: float,
                          sigma_start_um: float = SIGMA_START_UM,
                          sigma_end_um: float = SIGMA_END_UM,
                          steps: int = _SIGMA_STEPS,
                          iters: int = _ITERS_PER_SIGMA,
                          image_shape_yx=None,
                          bracket: bool = True) -> NucleiAlignResult:
    """Fit the similarity that puts *source_yx* onto *target_yx*, from *seed_3x3*.

    Both point sets are in the pixel frame the seed maps between: *source_yx* in
    H&E pixels (already flipped, if the H&E is), *target_yx* in Xenium morphology
    pixels. ``pixel_size_um`` is the Xenium pixel size, which is the frame every
    reported distance is in.

    The fit is a ``bracket_seed`` pass followed by annealed soft assignment; see
    that function for why the bracket is needed and the module docstring for why
    the anneal is not ICP. ``bracket=False`` runs the anneal alone.
    """
    from scipy.spatial import cKDTree

    target = np.asarray(target_yx, dtype=np.float64)
    source = np.asarray(source_yx, dtype=np.float64)
    if len(target) < 50 or len(source) < 50:
        raise ValueError(
            f"not enough nuclei to register — {len(source)} in the H&E, "
            f"{len(target)} nuclear masks")

    tree = cKDTree(target)
    # The target's own nearest-neighbour spacing, which sets how many neighbours
    # each sigma needs. Sampled rather than measured over every point: it is a
    # scale, and 20k points give it to well under a percent.
    sample = target if len(target) <= 20000 else target[
        np.linspace(0, len(target) - 1, 20000).astype(int)]
    spacing_px = float(np.median(cKDTree(target).query(sample, k=2, workers=-1)[0][:, 1]))

    seed = np.asarray(seed_3x3, dtype=np.float64)
    if image_shape_yx is None:
        image_shape_yx = (float(source[:, 0].max()), float(source[:, 1].max()))
    # Cross the part of the seed's error the anneal cannot see, before annealing.
    bracket_peak = 0.0
    bracket_shift = 0.0
    if bracket:
        bracketed, bracket_peak = bracket_seed(
            tree, source, seed, pixel_size_um, image_shape_yx, spacing_px)
        bracket_shift = _seed_shift_um(bracketed, seed, image_shape_yx, pixel_size_um)
        seed = bracketed

    affine = seed.copy()
    sigmas = np.geomspace(sigma_start_um / pixel_size_um,
                          sigma_end_um / pixel_size_um, steps)
    for sigma in sigmas:
        k = _neighbours_for(sigma, spacing_px)
        floor = _OUTLIER * k * np.exp(-4.5)     # as if 0.3k targets sat at 3 sigma
        for _ in range(iters):
            moved = _apply(affine, source)
            dist, idx = tree.query(moved, k=k, workers=-1)
            weight = np.exp(-0.5 * (dist / sigma) ** 2)
            total = weight.sum(axis=1)
            # The target each source point is pulled toward: the mean of its
            # neighbours weighted by how well each explains it. At large sigma
            # that is the local nuclear density; at small sigma it is one nucleus.
            #
            # A point with no neighbour inside any meaningful weight is pulled
            # toward *itself*, so it exerts no force. Dividing by a small epsilon
            # instead puts its virtual target at the origin, and when that happens
            # to every point — an H&E of a different section, a seed hundreds of
            # microns out — the least-squares answer is to collapse the whole
            # source cloud onto the origin. Measured before this line: scale
            # 6.8e-28 and an *enrichment of 26x*, because every point then lands
            # on top of some nucleus. It reported itself confident.
            live = total > 1e-12
            virtual = np.where(
                live[:, None],
                (weight[:, :, None] * target[idx]).sum(axis=1)
                / np.where(live, total, 1.0)[:, None],
                moved)
            confidence = total / (total + floor)
            affine = _weighted_similarity(source, virtual, confidence + 1e-9)

    moved = _apply(affine, source)
    dist, _ = tree.query(moved, workers=-1)
    dist_um = dist * pixel_size_um
    matched = dist_um < MATCH_RADIUS_UM
    n_matched = int(matched.sum())

    # Chance is not negligible and depends on the density, so it is divided out
    # rather than assumed small: for a Poisson target field, the probability that
    # a point unrelated to it has a target within r is 1 - exp(-lambda pi r^2),
    # with lambda from the same median-spacing relation _neighbours_for uses.
    lam = np.log(2.0) / (np.pi * spacing_px ** 2)
    chance = 1.0 - np.exp(-lam * np.pi * (MATCH_RADIUS_UM / pixel_size_um) ** 2)
    enrichment = float((n_matched / len(source)) / chance) if chance > 0 else float("inf")

    shift = _seed_shift_um(affine, np.asarray(seed_3x3, dtype=np.float64),
                           image_shape_yx, pixel_size_um)

    scale = float(np.hypot(affine[0, 0], affine[0, 1]))
    rotation = float(np.degrees(np.arctan2(affine[0, 1], affine[0, 0])))
    notes = []
    if n_matched < 200:
        notes.append(f"Only {n_matched} nuclei matched; the fit rests on very few.")
    return NucleiAlignResult(
        affine_3x3_yx=affine,
        n_target=len(target), n_source=len(source), n_matched=n_matched,
        matched_fraction=n_matched / len(source),
        median_residual_um=float(np.median(dist_um[matched])) if n_matched else float("nan"),
        enrichment=enrichment,
        seed_shift_um=shift,
        bracket_shift_um=bracket_shift, bracket_peak=bracket_peak,
        scale=scale, rotation_deg=rotation,
        confident=bool(enrichment >= MIN_ENRICHMENT and n_matched >= 200),
        notes=notes,
    )


def _seed_shift_um(fitted, seed, image_shape_yx, pixel_size_um, grid=20):
    """Mean distance, in microns, between where two transforms put the H&E."""
    from palms.utils.affine_compare import disagreement_um

    return float(disagreement_um(fitted, seed, image_shape_yx,
                                 pixel_size=pixel_size_um, grid=grid)["mean"])


# ─── Orchestration ────────────────────────────────────────────────────────────

def register_he_nuclei_steps(nucleus_labels, he_image, seed_3x3, pixel_size_um: float,
                             he_pixel_size_um: float | None = None,
                             label_downsample: float = 1.0,
                             he_downsample: float = 1.0,
                             he_shape_yx=None,
                             flip_v: bool = False, flip_h: bool = False):
    """Register an H&E to Xenium by matching its nuclei to the nuclear masks.

    A **generator**, yielding ``(stage, fraction_done)`` before each of the three
    stages and returning the ``NucleiAlignResult`` (as ``StopIteration.value``).
    ``register_he_nuclei`` is the blocking wrapper; the GUI drives this one
    directly, because a napari ``thread_worker`` reports progress by yielding and
    a callback invoked from the worker thread would be touching Qt from off the
    main thread. One mechanism, so the two cannot drift.

    ``nucleus_labels`` and ``he_image`` are single 2-D / 3-D arrays — pass the
    **finest level** of each; see the module docstring for the measurement that
    says so. ``seed_3x3`` is the starting transform in napari (y, x), mapping
    *flipped* H&E pixels to Xenium pixels, which is exactly what
    ``compute_coarse_affine`` returns and what the tab composes.

    ``he_pixel_size_um`` is only used to make the detector physical; when the H&E
    declares none it is derived from the seed's scale, which is a measurement of
    the same quantity and is always available because a seed is required.

    The detection runs on the H&E *unflipped* and the points are put into the
    flipped frame afterwards. A flip is a relabelling of coordinates and the
    haematoxylin field is unchanged by it, so flipping the array instead would be
    a 1.2 GB copy for nothing.
    """
    seed = np.asarray(seed_3x3, dtype=np.float64)
    if he_pixel_size_um is None:
        # scale is H&E px -> Xenium px, so px * scale * um-per-Xenium-px is the
        # H&E pixel in microns.
        he_pixel_size_um = float(np.hypot(seed[0, 0], seed[0, 1])) * pixel_size_um

    yield ("reading the nuclear masks", 0.0)
    target = nucleus_centroids(nucleus_labels, downsample=label_downsample)

    yield (f"finding nuclei in the H&E ({len(target):,} masks to match)", 0.25)
    detections = detect_he_nuclei(he_image,
                                  pixel_size_um=he_pixel_size_um * he_downsample,
                                  downsample=he_downsample)
    source = detections[:, :2]

    if he_shape_yx is None:
        base = finest_level(he_image)
        channel_first = (getattr(base, "ndim", 0) == 3
                         and base.shape[0] in (3, 4)
                         and base.shape[-1] not in (3, 4))
        spatial = ((base.shape[1], base.shape[2]) if channel_first
                   else (base.shape[0], base.shape[1]))
        he_shape_yx = (spatial[0] * he_downsample, spatial[1] * he_downsample)

    if flip_v or flip_h:
        source = _apply(flip_matrix(he_shape_yx, flip_v, flip_h), source)

    yield (f"matching {len(source):,} H&E nuclei to {len(target):,} masks", 0.75)
    result = fit_nuclei_similarity(target, source, seed, pixel_size_um=pixel_size_um,
                                   image_shape_yx=he_shape_yx)
    yield ("done", 1.0)
    return result


def register_he_nuclei(*args, **kwargs) -> NucleiAlignResult:
    """Blocking form of :func:`register_he_nuclei_steps`."""
    steps = register_he_nuclei_steps(*args, **kwargs)
    while True:
        try:
            next(steps)
        except StopIteration as done:
            return done.value
