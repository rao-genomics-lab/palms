"""Comparing two estimates of the same H&E registration.

Two conventions meet here and getting either backwards is the whole risk: PALMS
writes a 3x3 in napari ``(y, x)``, 10x ships ``<sample>_he_imagealignment.csv``
in ``(x, y)``, and both map H&E pixels to morphology pixels. A swap done wrong
reports a large disagreement between transforms that agree, or — worse — a small
one between transforms that do not, and either reads as a statement about the
registration rather than about the arithmetic.

The number that means something is not a difference of matrix entries, which is
uninterpretable. It is the disagreement *in image space*: push a grid of H&E
pixels through both transforms and report the distance distribution in microns.

Used by ``scripts/compare_he_registration.py`` (landmark fit against 10x) and
``scripts/score_coarse_align.py`` (coarse align against either). The two report
the same quantity on purpose, so a coarse figure and a landmark figure can be
put side by side.
"""

from __future__ import annotations

import numpy as np

#: The exchange matrix. ``P M P`` is the same map read in the other axis order.
P = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)


def to_yx(m_xy: np.ndarray) -> np.ndarray:
    """A 3x3 affine in (x, y) expressed in napari's (y, x)."""
    return P @ np.asarray(m_xy, dtype=float) @ P


def decompose(m_yx: np.ndarray) -> dict:
    """Scale, rotation and translation of a similarity, read in (y, x).

    The swap turns an (x, y) rotation by theta into
    ``[[s cos, s sin], [-s sin, s cos]]`` in (y, x), so theta comes off the
    *top* row. Reading it off the bottom row instead returns ``-theta``: the
    magnitudes and every pairwise difference survive that, which is why it
    looked right against real data, and only a transform built from a known
    angle exposes it.
    """
    m_yx = np.asarray(m_yx, dtype=float)
    return {
        "scale": float(np.hypot(m_yx[0, 0], m_yx[0, 1])),
        "rotation_deg": float(np.degrees(np.arctan2(m_yx[0, 1], m_yx[0, 0]))),
        "translation_yx": m_yx[:2, 2].tolist(),
        "mirrored": bool(np.linalg.det(m_yx[:2, :2]) < 0),
    }


def apply(m_yx: np.ndarray, pts_yx: np.ndarray) -> np.ndarray:
    pts_yx = np.asarray(pts_yx, dtype=float)
    homo = np.hstack([pts_yx, np.ones((len(pts_yx), 1))])
    return (np.asarray(m_yx, dtype=float) @ homo.T).T[:, :2]


def angle_difference_deg(a_deg: float, b_deg: float) -> float:
    """Signed difference between two angles, wrapped to (-180, 180].

    +179.95 and -179.97 are 0.08 degrees apart, not 359.9. Reporting the
    unwrapped difference turns a correct 180-degree match into a headline
    failure.
    """
    return float((a_deg - b_deg + 180.0) % 360.0 - 180.0)


def grid_over(shape_yx, n: int = 50) -> np.ndarray:
    """An ``n`` x ``n`` grid of (y, x) points spanning an image."""
    h, w = shape_yx
    gy, gx = np.meshgrid(np.linspace(0, h, n), np.linspace(0, w, n), indexing="ij")
    return np.column_stack([gy.ravel(), gx.ravel()])


def disagreement_um(a_yx, b_yx, shape_yx, pixel_size: float = 0.2125,
                    grid: int = 50) -> dict:
    """How far apart two transforms place the same H&E pixels, in microns."""
    pts = grid_over(shape_yx, grid)
    d = np.linalg.norm(apply(a_yx, pts) - apply(b_yx, pts), axis=1) * pixel_size
    return {
        "mean": float(d.mean()),
        "median": float(np.median(d)),
        "p95": float(np.percentile(d, 95)),
        "max": float(d.max()),
        "grid_points": int(d.size),
    }
