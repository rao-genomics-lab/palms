#!/usr/bin/env python
"""Compare PALMS's landmark H&E registration against 10x's own alignment matrix.

The residual PALMS reports is *self-referential*: it measures how well a
similarity fits the very points that were clicked, so it cannot tell you the
registration is right, only that it is internally consistent. 10x ships
`<sample>_he_imagealignment.csv` for this dataset — an independent estimate of
the same transform — so the disagreement between the two is an external number.

Conventions, both verified rather than assumed:
  * PALMS  : 3x3 in napari (y, x), H&E pixels -> morphology pixels.
             registration.compute_landmark_affine, saved as `affine_3x3_yx`.
  * 10x    : 3x3 in (x, y). Direction confirmed by pushing the H&E corners
             through it: they land in x -850..34640, y -3518..15750, matching the
             morphology image's 34155 x 13770. So it is H&E -> morphology too.
  * Swap   : M_yx = P @ M_xy @ P, with P the exchange matrix.

The number that matters is not a difference of matrix entries, which is not
interpretable. It is the disagreement in image space: push a grid of H&E pixels
through both transforms and report the distance distribution in microns.

Measured once, on Xenium_V1_human_Pancreas_FFPE with 3 landmarks (2026-09-02):
scale differed by 0.0173% and rotation by 0.0083 deg (-89.8796 vs -89.8880); the self-referential
residual was 0.97 um mean and the independent disagreement 0.93 um mean / 1.81 um
max over the whole H&E. **The two agreeing is the point.** A landmark residual
alone cannot detect systematically mis-clicked pairs, because it is fitted to
those pairs; that it matches the independent number to within 0.05 um is what
says nothing systematic is hiding in it.

Usage:
    python scripts/compare_he_registration.py \
        <dataset>/landmarks.json <sample>_he_imagealignment.csv \
        --he-shape 27502 14896 --pixel-size 0.2125 --out report.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

P = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)


def to_yx(m_xy: np.ndarray) -> np.ndarray:
    """A 3x3 affine in (x, y) expressed in napari's (y, x)."""
    return P @ m_xy @ P


def decompose(m_yx: np.ndarray) -> dict:
    """Scale, rotation and translation of a similarity, read in (y, x).

    The swap turns an (x, y) rotation by theta into
    ``[[s cos, s sin], [-s sin, s cos]]`` in (y, x), so theta comes off the
    *top* row. Reading it off the bottom row instead returns ``-theta``: the
    magnitudes and every pairwise difference survive that, which is why it
    looked right against real data, and only a transform built from a known
    angle exposes it.
    """
    return {
        "scale": float(np.hypot(m_yx[0, 0], m_yx[0, 1])),
        "rotation_deg": float(np.degrees(np.arctan2(m_yx[0, 1], m_yx[0, 0]))),
        "translation_yx": m_yx[:2, 2].tolist(),
    }


def apply(m_yx: np.ndarray, pts_yx: np.ndarray) -> np.ndarray:
    homo = np.hstack([pts_yx, np.ones((len(pts_yx), 1))])
    return (m_yx @ homo.T).T[:, :2]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("landmarks", type=Path, help="landmarks.json from Save Landmarks...")
    ap.add_argument("alignment", type=Path, help="10x <sample>_he_imagealignment.csv")
    ap.add_argument("--he-shape", type=int, nargs=2, metavar=("H", "W"), required=True,
                    help="H&E image shape in pixels (rows cols)")
    ap.add_argument("--pixel-size", type=float, default=0.2125, help="um per morphology px")
    ap.add_argument("--grid", type=int, default=50, help="grid points per axis")
    ap.add_argument("--out", type=Path, default=None, help="write the report as JSON")
    args = ap.parse_args()

    lm = json.loads(args.landmarks.read_text())
    if "affine_3x3_yx" not in lm:
        raise SystemExit("no affine in the landmarks file — press Compute Registration "
                         "before Save Landmarks...")
    ours = np.array(lm["affine_3x3_yx"], dtype=float)
    theirs = to_yx(np.loadtxt(args.alignment, delimiter=","))

    # A grid over the H&E, in its own pixels, (y, x).
    h, w = args.he_shape
    gy, gx = np.meshgrid(np.linspace(0, h, args.grid),
                         np.linspace(0, w, args.grid), indexing="ij")
    grid = np.column_stack([gy.ravel(), gx.ravel()])

    d_px = np.linalg.norm(apply(ours, grid) - apply(theirs, grid), axis=1)
    d_um = d_px * args.pixel_size

    # The self-referential number, for contrast.
    xen = np.array(lm["xenium_landmarks_yx"], dtype=float)
    he = np.array(lm["he_landmarks_yx"], dtype=float)
    own = np.linalg.norm(apply(ours, he) - xen, axis=1) * args.pixel_size

    report = {
        "n_landmarks": len(xen),
        "palms": decompose(ours),
        "tenx": decompose(theirs),
        "landmark_residual_um": {
            "mean": float(own.mean()), "max": float(own.max()),
            "per_landmark": own.tolist(),
        },
        "disagreement_with_10x_um": {
            "mean": float(d_um.mean()), "median": float(np.median(d_um)),
            "p95": float(np.percentile(d_um, 95)), "max": float(d_um.max()),
            "grid_points": int(d_um.size),
        },
    }

    print(f"Landmarks: {report['n_landmarks']}")
    print("\n                        scale    rotation      translation (y, x)")
    for who, key in (("PALMS (landmarks)", "palms"), ("10x  (shipped)   ", "tenx")):
        d = report[key]
        print(f"  {who}  {d['scale']:.5f}  {d['rotation_deg']:9.4f}   "
              f"({d['translation_yx'][0]:11.2f}, {d['translation_yx'][1]:10.2f})")
    ds, ts = report["palms"]["scale"], report["tenx"]["scale"]
    print(f"\n  scale differs by {abs(ds - ts) / ts * 100:.4f}%, "
          f"rotation by {abs(report['palms']['rotation_deg'] - report['tenx']['rotation_deg']):.4f} deg")

    r = report["landmark_residual_um"]
    print(f"\nSelf-referential (fit against the clicked points):")
    print(f"  mean {r['mean']:.2f} um, max {r['max']:.2f} um")
    d = report["disagreement_with_10x_um"]
    print(f"\nIndependent (vs 10x's matrix, over {d['grid_points']:,} points across the H&E):")
    print(f"  mean {d['mean']:.2f} um, median {d['median']:.2f} um, "
          f"p95 {d['p95']:.2f} um, max {d['max']:.2f} um")

    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
