#!/usr/bin/env python3
"""
build_custom_segmentation.py — Stage 2: Build custom segmentation assets for Xenium Viewer.

Reads the intermediate files produced by extract_seurat_segmentation.R and builds:
  custom_labels.zarr/      Multi-scale integer label raster (same pixel grid as native labels)
  custom_segmentation.h5ad AnnData with counts + obsm['spatial'] (xy µm) + obs['cell_id']
  custom_segmentation.json Metadata JSON consumed by the viewer

Usage:
  conda activate xenium_viewer
  python scripts/build_custom_segmentation.py <xenium_dir> <stage1_dir> [--out <output_dir>]

Arguments:
  xenium_dir   Path to Xenium output directory (contains cells.parquet, sdata_cached.zarr, etc.)
  stage1_dir   Directory containing Stage 1 outputs (segmentation_polygons.csv, counts.mtx, etc.)
  --out        Output directory (default: <xenium_dir>/custom_segmentation/)
  --scales     Number of pyramid levels to write in the zarr (default: 4)
  --chunk      Zarr chunk size in pixels (default: 1024)
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import zarr
from scipy.io import mmread
from scipy.ndimage import zoom
from shapely.geometry import Polygon
from shapely import contains_xy


def _read_pixel_size(xenium_dir: Path) -> float:
    """Read pixel size (µm/pixel) from experiment.xenium or fall back to Xenium default."""
    import json as _json
    exp = xenium_dir / "experiment.xenium"
    if exp.exists():
        with open(exp) as f:
            meta = _json.load(f)
        px = meta.get("pixel_size", None)
        if px:
            return float(px)
    print("  Warning: could not read pixel_size from experiment.xenium; using 0.2125 µm/px")
    return 0.2125


def _get_label_image_shape(xenium_dir: Path) -> tuple[int, int]:
    """Return (H, W) of the full-resolution label raster from the sdata zarr cache."""
    zarr_path = xenium_dir / "sdata_cached.zarr"
    if zarr_path.exists():
        z = zarr.open(str(zarr_path), mode="r")
        for lbl_key in ["cell_labels", "nucleus_labels"]:
            if "labels" in z and lbl_key in z["labels"]:
                lbl = z["labels"][lbl_key]
                # Keys are '0','1',... — '0' is finest resolution
                if "0" in lbl:
                    arr = lbl["0"]
                    return arr.shape[-2], arr.shape[-1]  # (H, W)
    # Fall back: read from cells.parquet centroids
    cells_path = xenium_dir / "cells.parquet"
    if cells_path.exists():
        df = pd.read_parquet(cells_path, columns=["x_centroid", "y_centroid"])
        # Infer from centroid range + margin
        pixel_size = _read_pixel_size(xenium_dir)
        max_x_px = int(df["x_centroid"].max() / pixel_size) + 100
        max_y_px = int(df["y_centroid"].max() / pixel_size) + 100
        print(f"  Warning: inferred shape ({max_y_px}, {max_x_px}) from centroid range")
        return max_y_px, max_x_px
    raise RuntimeError(
        "Cannot determine label image shape. Ensure sdata_cached.zarr exists or run the viewer once."
    )


def _validate_coordinate_system(
    polys_df: pd.DataFrame, xenium_dir: Path, pixel_size: float
) -> str:
    """
    Determine whether Seurat polygon coordinates are in pixel or micron units
    by comparing the coordinate range against the known image extent.

    For Xenium (pixel_size ≈ 0.2125 µm/px) the two ranges differ ~4.7×, making
    them unambiguous regardless of whether barcodes match between segmentations.

    Returns 'pixel' or 'micron'.
    """
    H, W = _get_label_image_shape(xenium_dir)
    max_extent_px = max(H, W)
    max_extent_um = max_extent_px * pixel_size

    coord_max = max(polys_df["x"].max(), polys_df["y"].max())
    print(f"  Coordinate range max: {coord_max:.1f}  "
          f"(image extent: {max_extent_px:.0f} px = {max_extent_um:.0f} µm)")

    # Midpoint between the two extents (in log space) as the decision threshold
    threshold = (max_extent_um * max_extent_px) ** 0.5
    if coord_max > threshold:
        print("  → Seurat coordinates are in pixel units")
        return "pixel"
    else:
        print("  → Seurat coordinates are in micron units")
        return "micron"


def _rasterize_chunk(args):
    """Worker function for multiprocessing rasterization of a spatial chunk."""
    shapes, out_shape, dtype = args
    from rasterio.features import rasterize as _rasterize
    return _rasterize(shapes, out_shape=out_shape, dtype=dtype, fill=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xenium_dir", type=Path)
    parser.add_argument("stage1_dir", type=Path)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--scales", type=int, default=4,
                        help="Number of pyramid levels (default: 4)")
    parser.add_argument("--chunk", type=int, default=1024,
                        help="Zarr chunk size in pixels (default: 1024)")
    parser.add_argument("--no-resume", action="store_true",
                        help="Ignore cached intermediate files and recompute from scratch")
    parser.add_argument("--coord-unit", choices=["auto", "pixel", "micron"], default="auto",
                        help="Coordinate unit of Seurat polygons (default: auto-detect)")
    args = parser.parse_args()

    xenium_dir: Path = args.xenium_dir
    stage1_dir: Path = args.stage1_dir
    out_dir: Path = args.out or (xenium_dir / "custom_segmentation")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Xenium dir : {xenium_dir}")
    print(f"Stage1 dir : {stage1_dir}")
    print(f"Output dir : {out_dir}")

    # ── Read pixel size ───────────────────────────────────────────────────────
    pixel_size = _read_pixel_size(xenium_dir)
    print(f"Pixel size : {pixel_size} µm/px")

    # ── Intermediate cache paths ───────────────────────────────────────────────
    cache_raster    = out_dir / "_cache_label_raster.npy"
    cache_centroids = out_dir / "_cache_centroids_px.npy"
    cache_barcodes  = out_dir / "_cache_barcodes.txt"
    cache_coord_unit = out_dir / "_cache_coord_unit.txt"

    def _cache_exists():
        return all(p.exists() for p in
                   [cache_raster, cache_centroids, cache_barcodes, cache_coord_unit])

    if not args.no_resume and _cache_exists():
        print("Resuming from cached intermediate files (use --no-resume to recompute)...")
        t0 = time.perf_counter()
        label_raster = np.load(str(cache_raster))
        centroids_px = np.load(str(cache_centroids))
        barcodes     = cache_barcodes.read_text().splitlines()
        coord_unit   = cache_coord_unit.read_text().strip()
        n_cells      = len(barcodes)
        label_ints   = np.arange(1, n_cells + 1, dtype=np.int32)
        H, W = label_raster.shape
        print(f"  Loaded: {n_cells:,} cells, raster ({H}×{W}), coord_unit={coord_unit} "
              f"({time.perf_counter()-t0:.1f}s)")
    else:
        # ── Load polygon CSV ──────────────────────────────────────────────────
        t0 = time.perf_counter()
        polys_path = stage1_dir / "segmentation_polygons.csv"
        print(f"Loading polygon vertices from {polys_path} ...")
        polys_df = pd.read_csv(polys_path)
        print(f"  {len(polys_df):,} vertex rows, {polys_df['cell_id'].nunique():,} unique cells "
              f"({time.perf_counter()-t0:.1f}s)")

        # ── Validate coordinate system ────────────────────────────────────────
        if args.coord_unit != "auto":
            coord_unit = args.coord_unit
            print(f"Coordinate unit: {coord_unit} (overridden via --coord-unit)")
        else:
            print("Validating coordinate system...")
            coord_unit = _validate_coordinate_system(polys_df, xenium_dir, pixel_size)
        if coord_unit == "micron":
            polys_df["x"] = polys_df["x"] / pixel_size
            polys_df["y"] = polys_df["y"] / pixel_size

        # ── Build shapely polygons + centroids ────────────────────────────────
        t0 = time.perf_counter()
        print("Building Shapely polygons and computing centroids...")
        grouped = polys_df.groupby("cell_id", sort=False)
        barcodes = list(grouped.groups.keys())
        n_cells = len(barcodes)
        label_ints = np.arange(1, n_cells + 1, dtype=np.int32)

        shapely_polys = []
        centroids_px = np.empty((n_cells, 2), dtype=np.float64)
        for i, bc in enumerate(barcodes):
            sub = grouped.get_group(bc)
            coords = sub[["x", "y"]].values
            if len(coords) < 3:
                poly = None
            else:
                poly = Polygon(coords)
                if not poly.is_valid:
                    poly = poly.buffer(0)
            shapely_polys.append(poly)
            if poly is not None and not poly.is_empty:
                c = poly.centroid
                centroids_px[i] = [c.x, c.y]
            else:
                centroids_px[i] = [sub["x"].mean(), sub["y"].mean()]
        print(f"  {n_cells:,} polygons built ({time.perf_counter()-t0:.1f}s)")

        # ── Determine raster dimensions ───────────────────────────────────────
        H, W = _get_label_image_shape(xenium_dir)
        print(f"Label raster shape: ({H}, {W})")

        # ── Rasterize polygons ────────────────────────────────────────────────
        try:
            from rasterio.features import rasterize
        except ImportError:
            print("ERROR: rasterio is required. Install with: pip install rasterio")
            sys.exit(1)

        t0 = time.perf_counter()
        print(f"Rasterizing {n_cells:,} polygons into ({H}×{W}) raster...")
        shapes = [
            (poly.__geo_interface__, int(lbl))
            for poly, lbl in zip(shapely_polys, label_ints)
            if poly is not None and not poly.is_empty
        ]
        label_raster = rasterize(
            shapes,
            out_shape=(H, W),
            dtype=np.int32,
            fill=0,
        )
        print(f"  Rasterization done ({time.perf_counter()-t0:.1f}s). "
              f"Non-zero pixels: {(label_raster > 0).sum():,}")

        # ── Save intermediates ────────────────────────────────────────────────
        print("Saving intermediate files for resume...")
        np.save(str(cache_raster), label_raster)
        np.save(str(cache_centroids), centroids_px)
        cache_barcodes.write_text("\n".join(str(b) for b in barcodes))
        cache_coord_unit.write_text(coord_unit)
        print("  Intermediates saved.")

    # ── Write multi-scale zarr ────────────────────────────────────────────────
    zarr_path = out_dir / "custom_labels.zarr"
    print(f"Writing multi-scale zarr to {zarr_path} ...")
    if zarr_path.exists():
        import shutil
        shutil.rmtree(zarr_path)

    store = zarr.open_group(str(zarr_path), mode="w")
    chunk = args.chunk
    current = label_raster
    _zarr_major = int(zarr.__version__.split(".")[0])
    for level in range(args.scales):
        key = str(level)
        if _zarr_major >= 3:
            arr = store.create_array(
                key, shape=current.shape,
                chunks=(chunk, chunk), dtype=np.int32,
            )
            arr[:] = current
        else:
            from numcodecs import Blosc
            store.create_dataset(key, data=current,
                                 chunks=(chunk, chunk), dtype=np.int32,
                                 compressor=Blosc(cname="lz4", clevel=5))
        print(f"  scale {level}: {current.shape}")
        if level < args.scales - 1:
            # Downsample 2× using nearest-neighbour (preserve label integers)
            current = current[::2, ::2]
    print(f"  Zarr written.")

    # ── Load counts matrix ────────────────────────────────────────────────────
    t0 = time.perf_counter()
    print("Loading counts matrix (Matrix Market)...")
    counts_path = stage1_dir / "counts.mtx"
    genes_path  = stage1_dir / "genes.txt"
    barcodes_path = stage1_dir / "barcodes.txt"

    counts_gxc = mmread(str(counts_path)).tocsr()  # genes × cells
    genes_list    = Path(genes_path).read_text().splitlines()
    barcodes_list = Path(barcodes_path).read_text().splitlines()
    print(f"  {len(genes_list)} genes × {len(barcodes_list)} cells ({time.perf_counter()-t0:.1f}s)")

    # Reorder to match our barcode order
    bc_to_col = {bc: i for i, bc in enumerate(barcodes_list)}
    col_order = np.array([bc_to_col[bc] for bc in barcodes if bc in bc_to_col])
    missing = [bc for bc in barcodes if bc not in bc_to_col]
    if missing:
        print(f"  Warning: {len(missing)} barcodes in polygons not found in counts matrix")
    counts_cxg = counts_gxc[:, col_order].T.tocsr()  # cells × genes

    # ── Load metadata ─────────────────────────────────────────────────────────
    meta_path = stage1_dir / "cell_metadata.csv"
    if meta_path.exists():
        meta = pd.read_csv(meta_path, index_col="cell_id")
        # Reorder to match barcode order
        meta = meta.reindex(barcodes)
    else:
        meta = pd.DataFrame(index=barcodes)

    # ── Build AnnData ─────────────────────────────────────────────────────────
    print("Building AnnData...")
    obs = meta.copy()
    obs.index.name = "barcode"   # avoid collision: index holds barcode strings,
    obs["cell_id"] = label_ints  # column holds integer label matching the raster

    var = pd.DataFrame(index=genes_list)
    var.index.name = "gene"

    # centroids_px is (N,2) xy-pixels; convert to xy-microns for obsm['spatial']
    spatial_um = centroids_px * pixel_size  # (N,2) xy µm

    # Convert ArrowStringArray columns/categories to plain object dtype (pandas 3.0 compat).
    # Temporarily disable infer_string so pandas doesn't re-wrap strings as ArrowStringArray.
    obs = obs.copy()
    _old_infer = pd.options.future.infer_string
    pd.options.future.infer_string = False
    try:
        def _fix_df(df):
            if pd.api.types.is_string_dtype(df.index):
                df.index = pd.Index(df.index.to_numpy(dtype=object))
            for col in df.columns:
                s = df[col]
                if isinstance(s.dtype, pd.CategoricalDtype):
                    cat = s.cat
                    if pd.api.types.is_string_dtype(cat.categories):
                        new_cats = cat.categories.astype(object)
                        df[col] = s.cat.rename_categories(dict(zip(cat.categories, new_cats)))
                elif pd.api.types.is_string_dtype(s):
                    df[col] = s.to_numpy(dtype=object)
            return df
        obs = _fix_df(obs)
        var = _fix_df(var)
    finally:
        pd.options.future.infer_string = _old_infer

    adata = anndata.AnnData(X=counts_cxg, obs=obs, var=var)
    adata.obsm["spatial"] = spatial_um.astype(np.float32)
    adata.uns["custom_segmentation"] = {
        "source": "seurat_rds",
        "pixel_size": pixel_size,
        "n_cells": n_cells,
        "coord_unit_original": coord_unit,
    }

    h5ad_path = out_dir / "custom_segmentation.h5ad"
    adata.write_h5ad(str(h5ad_path))
    print(f"  AnnData saved: {h5ad_path} ({adata.n_obs:,} cells × {adata.n_vars} genes)")

    # ── Write metadata JSON ───────────────────────────────────────────────────
    meta_json = {
        "pixel_size": pixel_size,
        "n_cells": n_cells,
        "label_raster": "custom_labels.zarr",
        "adata": "custom_segmentation.h5ad",
        "coord_unit_original": coord_unit,
        "raster_shape": [H, W],
        "n_scales": args.scales,
    }
    json_path = out_dir / "custom_segmentation.json"
    json_path.write_text(json.dumps(meta_json, indent=2))
    print(f"  Metadata JSON saved: {json_path}")

    print("\nStage 2 complete.")
    print(f"Output directory: {out_dir}")
    print("Load in the viewer via the 'Segmentation' tab → 'Load Custom Segmentation...'")
    print(f"  Select: {h5ad_path}")


if __name__ == "__main__":
    main()
