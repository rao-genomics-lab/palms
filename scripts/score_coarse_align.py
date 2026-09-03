#!/usr/bin/env python
"""Score PALMS's Coarse Align against an independently known transform.

`Coarse Align` produces the starting point for landmark registration, and "it
rarely lands close enough" is not a number. This makes it one: run the coarse
search exactly as the GUI does, then report how far its transform places the
H&E from a reference — 10x's shipped `<sample>_he_imagealignment.csv` where the
dataset has one, otherwise the landmark fit stored in the session.

The metric is the same one `scripts/compare_he_registration.py` reports for the
landmark fit — the distance distribution, in microns, between where the two
transforms put a grid of H&E pixels — so a coarse figure and a landmark figure
can be read side by side.

Headless: reads the zarr store directly, with no SpatialData load and no napari,
the idiom `cache_repair.verify` and `scripts/verify_notebook.py` already use. It
runs over ssh on a dataset that is too big to open.

Measured 2026-09-02, against the implementation this replaced:

    Xenium_V1_human_Pancreas_FFPE   scale 0.5515 -> 1.2749  (truth 1.2891)
                                    rotation -55.0 -> -89.85 deg (truth -89.88)
                                    disagreement  n/a -> 45 um mean, 88 um max
    demo_data/crop_6                scale 1.1702 -> 2.3521  (truth 2.3532)
                                    rotation -174.83 -> -179.95 deg (truth -179.97)
                                    disagreement  n/a -> 5 um mean, 6 um max

Usage:
    python scripts/score_coarse_align.py <dataset_dir> \
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


def _levels(group):
    """Multiscale level names, coarsest last."""
    return sorted(group.array_keys() if hasattr(group, "array_keys") else group,
                  key=lambda k: int("".join(c for c in k if c.isdigit()) or 0))


def _pick(group, min_long_side=384):
    """The coarsest level whose longest spatial side is still >= min_long_side."""
    names = _levels(group)
    chosen = names[0]
    for name in names:
        shape = group[name].shape
        if max(sorted(shape)[-2:]) >= min_long_side:
            chosen = name
    return chosen


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", type=Path, help="dataset directory (holding sdata_cached.zarr)")
    ap.add_argument("--alignment", type=Path, default=None,
                    help="10x <sample>_he_imagealignment.csv to score against; "
                         "without it the session's landmark affine is the reference")
    ap.add_argument("--pixel-size", type=float, default=None,
                    help="um per morphology pixel (default: read from the store)")
    ap.add_argument("--grid", type=int, default=50, help="grid points per axis")
    ap.add_argument("--no-mirror", action="store_true", help="skip reflection hypotheses")
    ap.add_argument("--out", type=Path, default=None, help="write the report as JSON")
    args = ap.parse_args()

    import zarr

    store = args.dataset / "sdata_cached.zarr"
    if not store.exists():
        store = args.dataset
    g = zarr.open_group(str(store), mode="r")
    if "images" not in g or "he_image" not in g["images"]:
        raise SystemExit(f"no he_image in {store} — register an H&E first")

    morph_group, he_group = g["images/morphology_focus"], g["images/he_image"]
    mk, hk = _pick(morph_group), _pick(he_group)
    morph = np.asarray(morph_group[mk])                              # (C, Y, X)
    he = np.transpose(np.asarray(he_group[hk]), (1, 2, 0))           # (Y, X, C)
    morph_full = morph_group[_levels(morph_group)[0]].shape
    he_full = he_group[_levels(he_group)[0]].shape
    target_ds = morph_full[-2] / morph.shape[-2]
    source_ds = he_full[-2] / he.shape[0]
    he_shape_yx = (he_full[-2], he_full[-1])

    session = dict(g["viewer_session"].attrs) if "viewer_session" in g else {}
    pixel_size = args.pixel_size or session.get("pixel_size") or 0.2125

    # --- the reference, and what makes it independent ---
    if args.alignment is not None:
        reference = to_yx(np.loadtxt(args.alignment, delimiter=","))
        reference_source = f"10x {args.alignment.name}"
    elif session.get("affine_3x3") is not None:
        reference = np.array(session["affine_3x3"], dtype=float)
        reference_source = ("session nuclei fit"
                            if session.get("affine_source") == "nuclei"
                            else "session landmark fit")
    else:
        raise SystemExit("no reference: pass --alignment, or register landmarks first")

    # --- the scale prior, exactly as the tab derives it ---
    he_px_um = session.get("he_pixel_size_um")
    if he_px_um:
        scale_prior, scale_source = he_px_um / pixel_size, "metadata"
    else:
        scale_prior = mask_area_scale(
            extract_tissue_mask_fluorescence(morph), extract_tissue_mask_he(he),
            target_ds, source_ds)
        scale_source = "tissue-area"

    target_field, source_field = build_alignment_fields(morph, he)
    t0 = time.perf_counter()
    result = compute_coarse_affine(
        target_field, source_field,
        target_downsample=target_ds, source_downsample=source_ds,
        scale_prior=scale_prior, scale_band=SCALE_BAND_FOR[scale_source],
        scale_source=scale_source, mirror=not args.no_mirror,
        source_shape_yx=he_shape_yx)
    elapsed = time.perf_counter() - t0

    # The reference includes whatever flip the session carries; the coarse result
    # is expressed in the flipped frame, so compose the same flip onto it before
    # comparing. Otherwise a mirrored dataset reads as a total failure.
    coarse = result.affine_3x3_yx @ flip_matrix(
        he_shape_yx, bool(session.get("flip_v")),
        bool(session.get("flip_h")) or result.mirrored)

    ours, theirs = decompose(coarse), decompose(reference)
    report = {
        "dataset": str(args.dataset),
        "reference": reference_source,
        "levels": {"morphology": mk, "he": hk,
                   "target_downsample": target_ds, "source_downsample": source_ds},
        "scale_prior": {"value": float(scale_prior), "source": scale_source,
                        "error_pct": float((scale_prior / theirs["scale"] - 1) * 100)},
        "coarse": dict(ours, score=result.score, runner_up_score=result.runner_up_score,
                       mirrored=result.mirrored, confident=result.confident,
                       seconds=round(elapsed, 1)),
        "reference_transform": theirs,
        "error": {
            "scale_pct": float((ours["scale"] / theirs["scale"] - 1) * 100),
            "rotation_deg": angle_difference_deg(ours["rotation_deg"], theirs["rotation_deg"]),
        },
        "disagreement_um": disagreement_um(coarse, reference, he_shape_yx,
                                           pixel_size=pixel_size, grid=args.grid),
    }

    print(f"{args.dataset}")
    print(f"  levels: morphology {mk} {morph.shape}, H&E {hk} {he.shape}")
    print(f"  reference: {reference_source}")
    print(f"  scale prior: {scale_prior:.4f} from {scale_source} "
          f"({report['scale_prior']['error_pct']:+.2f}%)")
    print(f"\n                    scale     rotation   mirrored")
    print(f"  Coarse Align   {ours['scale']:9.4f}  {ours['rotation_deg']:9.2f}   "
          f"{str(result.mirrored):>5}   (match {result.score:.3f}, "
          f"next distinct {result.runner_up_score:.3f}, {elapsed:.1f}s)")
    print(f"  reference      {theirs['scale']:9.4f}  {theirs['rotation_deg']:9.2f}   "
          f"{str(theirs['mirrored']):>5}")
    print(f"\n  error: scale {report['error']['scale_pct']:+.2f}%, "
          f"rotation {report['error']['rotation_deg']:+.3f} deg")
    d = report["disagreement_um"]
    print(f"  disagreement over {d['grid_points']:,} points across the H&E: "
          f"mean {d['mean']:.0f} um, median {d['median']:.0f} um, "
          f"p95 {d['p95']:.0f} um, max {d['max']:.0f} um")
    if not result.confident:
        print("\n  LOW CONFIDENCE — the search did not find a clear match.")

    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
