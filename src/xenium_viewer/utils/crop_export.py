"""
Crop the loaded dataset to a user-drawn polygon and write it out as a
standalone, independently-openable xenium-viewer data directory.

Images/labels are cropped to the polygon's pixel bounding box. Cells and
transcripts are filtered to the exact polygon via true point-in-polygon
tests (mirrors the technique used by tabs/tab_roi.py), so the exported
table/points are precise even though the raster extends slightly beyond
the drawn shape.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Callable, TYPE_CHECKING

import numpy as np
import pandas as pd
from shapely import contains_xy
from shapely.affinity import scale as shapely_scale
from shapely.geometry import Polygon

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

_PYRAMID_FACTORS = [2, 2, 2, 2, 2]
_MIN_LEVEL_SIZE = 8  # stop building pyramid levels once a dim would drop below this


class CropExportError(Exception):
    """User-facing crop/export failure (empty region, invalid geometry, etc.)."""


def _safe_scale_factors(height: int, width: int) -> list:
    """Pyramid scale factors that always yield at least one downsampled level.

    Image2DModel/Labels2DModel.parse return a bare DataArray (not a DataTree)
    when scale_factors is falsy, but the app's multiscale-loading code
    (_extract_dt_scales) assumes a DataTree on reopen — so always build at
    least one level, even for small crops.
    """
    factors = []
    h, w = height, width
    for f in _PYRAMID_FACTORS:
        h, w = h // f, w // f
        if h < _MIN_LEVEL_SIZE or w < _MIN_LEVEL_SIZE:
            break
        factors.append(f)
    return factors or [2]


def _extract_dt_scales(dt) -> list:
    """Extract ordered dask arrays (highest to lowest res) from a spatialdata
    multiscale DataTree.

    Duplicated from app.py/tab_segmentation.py: utils modules shouldn't
    import from app.py (which imports tab modules at call time).
    """
    def _sort_key(name):
        nums = re.findall(r"\d+", name)
        return int(nums[0]) if nums else 0

    scales = []
    for name in sorted(dt.children.keys(), key=_sort_key):
        child = dt.children[name]
        ds = getattr(child, "ds", None)
        if ds is None:
            continue
        if "image" in ds:
            scales.append(ds["image"].data)
        elif ds.data_vars:
            first = next(iter(ds.data_vars))
            scales.append(ds[first].data)
    return scales


def _polygon_pixel_bbox(polygon_yx: np.ndarray, shape_yx: tuple) -> tuple:
    """Return (row_min, row_max, col_min, col_max), clipped to shape_yx (exclusive upper bounds)."""
    row_min = max(0, int(np.floor(polygon_yx[:, 0].min())))
    row_max = min(shape_yx[0], int(np.ceil(polygon_yx[:, 0].max())))
    col_min = max(0, int(np.floor(polygon_yx[:, 1].min())))
    col_max = min(shape_yx[1], int(np.ceil(polygon_yx[:, 1].max())))
    return row_min, row_max, col_min, col_max


def crop_and_export(
    ctx: "ViewerContext",
    polygon_yx: np.ndarray,
    output_dir: Path,
    name: str,
    progress_cb: Callable[[int, str], None] | None = None,
) -> Path:
    """Crop the loaded dataset to *polygon_yx* and write a standalone
    xenium-viewer data directory at ``output_dir / name``.

    Read-only with respect to *ctx* — safe to call from a background thread.
    Writes to a hidden staging directory first and only moves it into place
    on full success, so a failure never leaves a partial dataset on disk.

    Returns the final output path.
    """
    def _progress(pct, msg):
        if progress_cb is not None:
            progress_cb(pct, msg)

    if ctx.segmentation_source != "xenium":
        raise CropExportError(
            "Crop Dataset requires native Xenium segmentation. "
            "Revert custom segmentation before cropping."
        )
    if ctx.sdata is None or ctx.adata is None:
        raise CropExportError("No dataset is loaded.")
    if ctx.data_path is None:
        raise CropExportError("Current dataset has no source directory (experiment.xenium unavailable).")

    polygon_yx = np.asarray(polygon_yx, dtype=np.float64)
    if len(polygon_yx) < 3:
        raise CropExportError("A crop region needs at least 3 points.")

    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, Labels2DModel, TableModel, PointsModel
    from spatialdata.transformations import Identity, Scale
    from xenium_viewer.loader import _convert_arrow_strings
    from xenium_viewer.preprocess import preprocess

    _progress(0, f"Validating region '{name}'...")

    poly_xy = Polygon(polygon_yx[:, ::-1])  # napari yx -> shapely xy, pixel space
    if not poly_xy.is_valid:
        poly_xy = poly_xy.buffer(0)
    if poly_xy.is_empty:
        raise CropExportError("Drawn region is degenerate (self-intersecting with no valid area).")

    full_shape_yx = ctx.morph_full_shape_yx
    if full_shape_yx is None:
        raise CropExportError("Morphology image shape is unavailable.")

    row_min, row_max, col_min, col_max = _polygon_pixel_bbox(polygon_yx, full_shape_yx)
    if row_max <= row_min or col_max <= col_min:
        raise CropExportError("Drawn region has zero area within the image bounds.")

    pixel_size = ctx.pixel_size
    origin_x_um = col_min * pixel_size
    origin_y_um = row_min * pixel_size

    # ── Cell membership: bbox pre-filter, then true point-in-polygon ─────────
    _progress(5, "Selecting cells inside region...")
    centroids_yx = ctx.centroids_yx
    bbox_pref = (
        (centroids_yx[:, 0] >= row_min) & (centroids_yx[:, 0] < row_max)
        & (centroids_yx[:, 1] >= col_min) & (centroids_yx[:, 1] < col_max)
    )
    bbox_idx = np.where(bbox_pref)[0]
    if len(bbox_idx) == 0:
        raise CropExportError("No cells inside this region.")
    inside = contains_xy(poly_xy, centroids_yx[bbox_idx, 1], centroids_yx[bbox_idx, 0])
    kept_idx = np.sort(bbox_idx[inside])
    if len(kept_idx) == 0:
        raise CropExportError("No cells inside this region.")

    # Map obs rows -> cell_labels raster pixel values via ctx.label_to_obs
    # (label_value -> obs_row_index), inverted. Xenium datasets vary in which
    # obs column actually holds this integer link (e.g. 'cell_labels' vs
    # 'cell_id' — 'cell_id' can be a non-numeric barcode string), so we reuse
    # the app's own label_to_obs mapping (loader.get_label_to_obs_mapping)
    # rather than re-guessing a column name.
    if ctx.label_to_obs is None:
        raise CropExportError("Cell-label mapping is unavailable.")
    label_to_obs = ctx.label_to_obs
    obs_to_label = np.zeros(len(ctx.adata), dtype=np.int64)
    valid = label_to_obs >= 0
    obs_to_label[label_to_obs[valid]] = np.arange(len(label_to_obs))[valid]

    adata_cropped = ctx.adata[kept_idx].copy()
    kept_cell_ids = obs_to_label[kept_idx]
    adata_cropped.obsm["spatial"] = adata_cropped.obsm["spatial"] - np.array([origin_x_um, origin_y_um])

    # ── Crop morphology image (bbox, full-res, rebuild pyramid) ──────────────
    _progress(20, "Cropping morphology image...")
    morph_scales = _extract_dt_scales(ctx.sdata.images["morphology_focus"])
    morph_full = morph_scales[0]
    cropped_img = np.asarray(morph_full[:, row_min:row_max, col_min:col_max].compute())
    scale_factors = _safe_scale_factors(cropped_img.shape[-2], cropped_img.shape[-1])
    img_element = Image2DModel.parse(
        cropped_img, dims=("c", "y", "x"),
        transformations={"global": Identity()},
        scale_factors=scale_factors,
    )

    # ── Crop cell_labels (bbox, then zero out non-kept cells) ────────────────
    _progress(35, "Cropping cell labels...")
    cl_scales = _extract_dt_scales(ctx.sdata.labels["cell_labels"])
    cropped_cl = np.asarray(cl_scales[0][row_min:row_max, col_min:col_max].compute())
    cropped_cl = cropped_cl.copy()
    cropped_cl[~np.isin(cropped_cl, kept_cell_ids)] = 0
    cl_element = Labels2DModel.parse(
        cropped_cl, dims=("y", "x"),
        transformations={"global": Identity()},
        scale_factors=scale_factors,
    )

    # ── Crop nucleus_labels (mask only if its ID space matches cell_id) ──────
    _progress(45, "Cropping nucleus labels...")
    nl_scales = _extract_dt_scales(ctx.sdata.labels["nucleus_labels"])
    cropped_nl = np.asarray(nl_scales[0][row_min:row_max, col_min:col_max].compute())
    cropped_nl = cropped_nl.copy()
    unique_nl = np.unique(cropped_nl)
    unique_nl = unique_nl[unique_nl > 0]
    overlap = np.isin(unique_nl, kept_cell_ids).mean() if len(unique_nl) else 1.0
    if overlap > 0.5:
        cropped_nl[~np.isin(cropped_nl, kept_cell_ids)] = 0
    else:
        print(
            f"Warning: nucleus_labels values do not appear to share cell_id's ID space "
            f"({overlap:.0%} overlap) — leaving nucleus_labels crop unmasked."
        )
    nl_element = Labels2DModel.parse(
        cropped_nl, dims=("y", "x"),
        transformations={"global": Identity()},
        scale_factors=scale_factors,
    )

    # ── Filter transcripts to the true polygon (in microns) ──────────────────
    # sdata.points['transcripts'] uses spatialdata's canonical column names
    # ("x", "y", "z", already in microns), not the raw Xenium parquet names
    # ("x_location", "y_location", "z_location") — those are only used in the
    # exported transcripts.parquet (below), to match what preprocess.py and
    # TranscriptLoader expect.
    _progress(60, "Filtering transcripts...")
    poly_xy_um = shapely_scale(poly_xy, xfact=pixel_size, yfact=pixel_size, origin=(0, 0))
    transcripts_ddf = ctx.sdata.points["transcripts"]

    def _filter_partition(pdf):
        if len(pdf) == 0:
            return pdf
        mask = contains_xy(poly_xy_um, pdf["x"].to_numpy(), pdf["y"].to_numpy())
        return pdf[mask]

    filtered_df = transcripts_ddf.map_partitions(_filter_partition).compute()
    filtered_df = filtered_df.copy()
    filtered_df["x"] = filtered_df["x"] - origin_x_um
    filtered_df["y"] = filtered_df["y"] - origin_y_um

    points_coords = {"x": "x", "y": "y"}
    if "z" in filtered_df.columns:
        points_coords["z"] = "z"
    points_element = PointsModel.parse(
        filtered_df,
        coordinates=points_coords,
        feature_key="feature_name",
        instance_key="cell_id",
        transformations={"global": Scale([1.0 / pixel_size, 1.0 / pixel_size], axes=("x", "y"))},
    )

    # ── Table ──────────────────────────────────────────────────────────────
    # Set/overwrite a 'cell_labels' obs column holding the raster pixel values
    # (kept_cell_ids) and link the table to it directly, regardless of how the
    # source dataset's own region/instance_key were named — this guarantees
    # get_label_to_obs_mapping() finds a consistent link when this exported
    # dataset is reopened.
    adata_cropped.obs["cell_labels"] = kept_cell_ids
    adata_cropped.obs["region"] = pd.Categorical(["cell_labels"] * adata_cropped.n_obs)
    table_element = TableModel.parse(
        adata_cropped, region="cell_labels", region_key="region", instance_key="cell_labels",
        overwrite_metadata=True,
    )

    # ── Assemble and write to a staging directory ────────────────────────────
    _progress(75, "Writing dataset...")
    new_sdata = SpatialData(
        images={"morphology_focus": img_element},
        labels={"cell_labels": cl_element, "nucleus_labels": nl_element},
        points={"transcripts": points_element},
        tables={"table": table_element},
    )
    _convert_arrow_strings(new_sdata)

    output_dir = Path(output_dir)
    final_dir = output_dir / name
    staging_dir = output_dir / f".{name}__crop_tmp"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True)

    try:
        # Copy experiment.xenium first so its mtime predates the zarr write —
        # loader.py's cache-freshness check requires the zarr be >= as fresh.
        shutil.copy(ctx.data_path / "experiment.xenium", staging_dir / "experiment.xenium")

        new_sdata.write(str(staging_dir / "sdata_cached.zarr"), overwrite=True)

        _progress(90, "Writing transcripts.parquet...")
        parquet_path = staging_dir / "transcripts.parquet"
        raw_schema_df = filtered_df.rename(columns={"x": "x_location", "y": "y_location", "z": "z_location"})
        raw_schema_df.attrs = {}   # spatialdata's transform metadata isn't JSON-serializable
        raw_schema_df.to_parquet(parquet_path, index=False)

        _progress(95, "Building transcript cache...")
        preprocess(parquet_path=parquet_path, cache_dir=staging_dir / "transcript_cache")
    except Exception:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise

    if final_dir.exists():
        shutil.rmtree(final_dir)
    shutil.move(str(staging_dir), str(final_dir))

    _progress(100, f"Done: {final_dir}")
    return final_dir
