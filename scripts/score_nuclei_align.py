#!/usr/bin/env python
"""Score PALMS's automatic nuclei registration against a known transform.

Runs the whole automatic chain exactly as the GUI does — Coarse Align, then the
nuclei fit seeded by it — and reports how far each places the H&E from a
reference: 10x's shipped `<sample>_he_imagealignment.csv` where the dataset has
one, otherwise the landmark fit stored in the session. The metric is the one
`scripts/compare_he_registration.py` and `scripts/score_coarse_align.py` report,
so a coarse figure, a landmark figure and an automatic figure read side by side.

Headless: reads the zarr store directly, with no SpatialData load and no napari.

**The reference is not a ruler fine enough for this.** On the pancreas the three
estimates agree with each other to about a micron in every pairing, while the
nuclei fit's own reproducibility is 0.02 um. So `--holdout` is the check that
means something and it is on by default: fit on one half of the *section*, then
score the matched-pair residual on the nuclei of the other half, which that fit
has never seen. It compares like with like — every transform is scored on the
same held-out nuclei, and the one that puts them on the nuclear masks best wins.

Measured 2026-09-02:

    Xenium_V1_human_Pancreas_FFPE   coarse 15.5 um -> nuclei 0.70 um vs 10x
                                    (landmark fit, by hand: 0.93 um)
                                    held-out residual  nuclei 1.06 um,
                                    10x 1.36 um, landmark 1.35 um
    demo_data/crop_6                coarse 17.5 um -> nuclei 1.31 um vs landmarks

Usage:
    python scripts/score_nuclei_align.py <dataset_dir> \
        [--alignment <sample>_he_imagealignment.csv] [--out report.json]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from palms.utils.affine_compare import (  # noqa: E402
    to_yx, decompose, disagreement_um, angle_difference_deg,
)
from palms.utils.registration import (  # noqa: E402
    SCALE_BAND_FOR, build_alignment_fields, compute_coarse_affine,
    extract_tissue_mask_fluorescence, extract_tissue_mask_he, flip_matrix,
    mask_area_scale,
)
from palms.utils.nuclei_registration import (  # noqa: E402
    detect_he_nuclei, fit_nuclei_similarity, nucleus_centroids,
)


def _levels(group):
    """Multiscale level names, coarsest last."""
    return sorted(group.array_keys() if hasattr(group, "array_keys") else group,
                  key=lambda k: int("".join(c for c in k if c.isdigit()) or 0))


def _pick(group, min_long_side=384):
    """The coarsest level whose longest spatial side is still >= min_long_side."""
    names = _levels(group)
    chosen = names[0]
    for name in names:
        if max(sorted(group[name].shape)[-2:]) >= min_long_side:
            chosen = name
    return chosen


def _apply(m, pts):
    return (np.asarray(m, float) @ np.hstack([pts, np.ones((len(pts), 1))]).T).T[:, :2]


def _held_out_residual(matrix, points, tree, pixel_size, radius_um=3.0):
    """Median distance from a held-out detection to its nearest nuclear mask.

    Capped at *radius_um* so a detection with no counterpart at all — H&E outside
    the Xenium field, a fold — is excluded rather than dominating the median. The
    same cap is applied to every transform, so the comparison stays like for like.
    """
    d, _ = tree.query(_apply(matrix, points), workers=-1)
    d = d * pixel_size
    keep = d < radius_um
    return float(np.median(d[keep])), float(keep.mean()), int(keep.sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", type=Path, help="dataset directory (holding sdata_cached.zarr)")
    ap.add_argument("--alignment", type=Path, default=None,
                    help="10x <sample>_he_imagealignment.csv to score against; "
                         "without it the session's landmark affine is the reference")
    ap.add_argument("--landmarks", type=Path, default=None,
                    help="a landmarks.json from Save Landmarks... to score against. "
                         "Worth reaching for: the session's fine affine is "
                         "overwritten by whichever method ran last, this one "
                         "included, but a saved landmark file is a record of a "
                         "registration made another way and it survives")
    ap.add_argument("--pixel-size", type=float, default=None,
                    help="um per morphology pixel (default: read from the store)")
    ap.add_argument("--grid", type=int, default=50, help="grid points per axis")
    ap.add_argument("--seed", choices=("coarse", "reference"), default="coarse",
                    help="what to seed the nuclei fit with (default: coarse, as the GUI does)")
    ap.add_argument("--no-holdout", action="store_true",
                    help="skip the spatially held-out comparison, which doubles the runtime")
    ap.add_argument("--out", type=Path, default=None, help="write the report as JSON")
    args = ap.parse_args()

    import zarr
    from scipy.spatial import cKDTree

    store = args.dataset / "sdata_cached.zarr"
    if not store.exists():
        store = args.dataset
    g = zarr.open_group(str(store), mode="r")
    if "images" not in g or "he_image" not in g["images"]:
        raise SystemExit(f"no he_image in {store} — register an H&E first")
    if "labels" not in g or "nucleus_labels" not in g["labels"]:
        raise SystemExit(f"no nucleus_labels in {store} — nothing to match the H&E against")

    morph_group, he_group = g["images/morphology_focus"], g["images/he_image"]
    mk, hk = _pick(morph_group), _pick(he_group)
    morph = np.asarray(morph_group[mk])
    he_low = np.transpose(np.asarray(he_group[hk]), (1, 2, 0))
    morph_full = morph_group[_levels(morph_group)[0]].shape
    he_full = he_group[_levels(he_group)[0]].shape
    target_ds = morph_full[-2] / morph.shape[-2]
    source_ds = he_full[-2] / he_low.shape[0]
    he_shape_yx = (he_full[-2], he_full[-1])

    session = dict(g["viewer_session"].attrs) if "viewer_session" in g else {}
    pixel_size = args.pixel_size or session.get("pixel_size") or 0.2125
    flip_v, flip_h = bool(session.get("flip_v")), bool(session.get("flip_h"))

    if args.alignment is not None:
        reference = to_yx(np.loadtxt(args.alignment, delimiter=","))
        reference_source = f"10x {args.alignment.name}"
    elif args.landmarks is not None:
        from palms.utils.registration import load_landmarks
        stored = load_landmarks(args.landmarks)
        if "affine_3x3_yx" not in stored:
            raise SystemExit(f"{args.landmarks} holds no affine — it was saved "
                             "before Compute Registration was pressed")
        reference = stored["affine_3x3_yx"]
        reference_source = f"landmark fit from {args.landmarks.name}"
    elif session.get("affine_3x3") is not None:
        # The session's fine affine is written by whichever method ran last, and
        # this script's own method is one of them. Scoring against a stored run
        # of itself is not a measurement: on demo_data/crop_6 it duly reported
        # 0.000 um of disagreement, which reads as a perfect result and says
        # nothing at all. `affine_source` exists to make that detectable.
        if session.get("affine_source") == "nuclei":
            raise SystemExit(
                "the stored registration was made by this same nuclei fit, so "
                "scoring against it would be circular.\n"
                "Pass --alignment <sample>_he_imagealignment.csv or --landmarks "
                "landmarks.json for an independent reference.\n"
                "(--no-holdout off, the held-out check below needs no reference "
                "at all — run with --seed coarse and read that section.)")
        if session.get("affine_source") is None:
            print("  NOTE: the stored affine does not say which method produced "
                  "it (written before affine_source existed). If it came from "
                  "this fit, the comparison below is circular — the held-out "
                  "check is not.")
        reference = np.array(session["affine_3x3"], dtype=float)
        reference_source = "session landmark fit"
    else:
        raise SystemExit("no reference: pass --alignment, or register landmarks first")

    # ── Coarse Align, exactly as the tab runs it ──────────────────────────────
    he_px_um = session.get("he_pixel_size_um")
    if he_px_um:
        scale_prior, scale_source = he_px_um / pixel_size, "metadata"
    else:
        scale_prior = mask_area_scale(
            extract_tissue_mask_fluorescence(morph), extract_tissue_mask_he(he_low),
            target_ds, source_ds)
        scale_source = "tissue-area"
    target_field, source_field = build_alignment_fields(morph, he_low)
    t0 = time.perf_counter()
    coarse_result = compute_coarse_affine(
        target_field, source_field,
        target_downsample=target_ds, source_downsample=source_ds,
        scale_prior=scale_prior, scale_band=SCALE_BAND_FOR[scale_source],
        scale_source=scale_source, source_shape_yx=he_shape_yx)
    coarse_seconds = time.perf_counter() - t0

    # The reference carries whatever flip the session has; the coarse result is
    # expressed in the flipped frame, so compose the same flip before comparing.
    coarse_flip = flip_matrix(he_shape_yx, flip_v,
                              flip_h or coarse_result.mirrored)
    coarse = coarse_result.affine_3x3_yx @ coarse_flip

    # ── The nuclei fit ────────────────────────────────────────────────────────
    print(f"{args.dataset}\n  levels: morphology {mk}, H&E {hk} (coarse search)")
    print(f"  building point sets at full resolution "
          f"({_levels(he_group)[0]}: {he_full[1]}x{he_full[2]}, "
          f"{_levels(morph_group)[0]}: {morph_full[-2]}x{morph_full[-1]})", flush=True)
    t0 = time.perf_counter()
    target_pts = nucleus_centroids(g[f"labels/nucleus_labels/{_levels(g['labels/nucleus_labels'])[0]}"])
    t_masks = time.perf_counter() - t0
    seed_scale = float(np.hypot(coarse_result.affine_3x3_yx[0, 0],
                                coarse_result.affine_3x3_yx[0, 1]))
    he_pixel = (he_px_um or seed_scale * pixel_size)
    t0 = time.perf_counter()
    detections = detect_he_nuclei(g[f"images/he_image/{_levels(he_group)[0]}"],
                                  pixel_size_um=he_pixel)
    t_detect = time.perf_counter() - t0
    source_pts = detections[:, :2]
    if flip_v or flip_h:
        source_pts = _apply(flip_matrix(he_shape_yx, flip_v, flip_h), source_pts)
    print(f"  {len(target_pts):,} nuclear masks ({t_masks:.0f}s), "
          f"{len(source_pts):,} H&E nuclei ({t_detect:.0f}s)", flush=True)

    seed = coarse if args.seed == "coarse" else reference
    t0 = time.perf_counter()
    fit = fit_nuclei_similarity(target_pts, source_pts, seed, pixel_size_um=pixel_size,
                                image_shape_yx=he_shape_yx)
    t_fit = time.perf_counter() - t0
    fine = fit.affine_3x3_yx

    ours, theirs = decompose(fine), decompose(reference)
    coarse_d = disagreement_um(coarse, reference, he_shape_yx, pixel_size, args.grid)
    fine_d = disagreement_um(fine, reference, he_shape_yx, pixel_size, args.grid)

    report = {
        "dataset": str(args.dataset),
        "reference": reference_source,
        "seed": args.seed,
        "point_sets": {"nuclear_masks": int(len(target_pts)),
                       "he_detections": int(len(source_pts)),
                       "he_pixel_size_um": float(he_pixel)},
        "coarse": dict(decompose(coarse), score=coarse_result.score,
                       confident=coarse_result.confident,
                       seconds=round(coarse_seconds, 1),
                       disagreement_um=coarse_d),
        "nuclei": dict(ours, matched=fit.n_matched,
                       matched_fraction=fit.matched_fraction,
                       median_residual_um=fit.median_residual_um,
                       enrichment=fit.enrichment,
                       seed_shift_um=fit.seed_shift_um,
                       confident=fit.confident,
                       seconds=round(t_masks + t_detect + t_fit, 1),
                       disagreement_um=fine_d),
        "reference_transform": theirs,
        "error": {
            "scale_pct": float((ours["scale"] / theirs["scale"] - 1) * 100),
            "rotation_deg": angle_difference_deg(ours["rotation_deg"], theirs["rotation_deg"]),
        },
    }

    print(f"  reference: {reference_source}\n")
    print(f"                    scale     rotation    disagreement with the reference")
    for label, mat, dd in (("Coarse Align", coarse, coarse_d), ("Nuclei fit  ", fine, fine_d)):
        dec = decompose(mat)
        print(f"  {label} {dec['scale']:9.4f}  {dec['rotation_deg']:9.3f}    "
              f"mean {dd['mean']:7.3f} um, p95 {dd['p95']:7.3f}, max {dd['max']:7.3f}")
    print(f"  reference    {theirs['scale']:9.4f}  {theirs['rotation_deg']:9.3f}")
    print(f"\n  {fit.summary()}")
    print(f"  ({t_masks:.0f}s masks + {t_detect:.0f}s detection + {t_fit:.0f}s fit)")

    # ── The independent check ─────────────────────────────────────────────────
    if not args.no_holdout:
        tree = cKDTree(target_pts)
        split = np.median(source_pts[:, 0])
        halves = {"lower half": source_pts[:, 0] >= split,
                  "upper half": source_pts[:, 0] < split}
        report["holdout"] = {}
        print("\n  Held-out check — fit on one half of the section, score the other:")
        for name, mask in halves.items():
            held = fit_nuclei_similarity(target_pts, source_pts[~mask], seed,
                                         pixel_size_um=pixel_size,
                                         image_shape_yx=he_shape_yx)
            rows = {}
            for label, mat in (("nuclei (held out)", held.affine_3x3_yx),
                               ("reference", reference),
                               ("coarse", coarse)):
                med, frac, n = _held_out_residual(mat, source_pts[mask], tree, pixel_size)
                rows[label] = {"median_residual_um": med, "matched_fraction": frac, "n": n}
            report["holdout"][name] = rows
            best = min(rows, key=lambda k: rows[k]["median_residual_um"])
            print(f"    scoring the {name} ({int(mask.sum()):,} nuclei it never saw):")
            for label, r in rows.items():
                mark = "  <-- best" if label == best else ""
                print(f"      {label:20s} median residual {r['median_residual_um']:.3f} um"
                      f"   ({r['matched_fraction'] * 100:.0f}% within 3 um){mark}")

    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
