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

import dask
import dask.array as da
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from shapely import contains_xy
from shapely.affinity import scale as shapely_scale
from shapely.geometry import Polygon

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

_PYRAMID_FACTORS = [2, 2, 2, 2, 2]
_MIN_LEVEL_SIZE = 8  # stop building pyramid levels once a dim would drop below this

# Transcript partitions can be far bigger than an image/label chunk (~150MB+
# compressed, all columns, vs ~30-130MB for an image/label chunk) — dask's
# default scheduler runs up to one task per CPU core concurrently, so a crop
# that overlaps most of the source table's partitions (e.g. a crop spanning
# most of the slide's width, which barely benefits from the `filters=`
# pruning above) can still pull several GB into memory at once even though no
# single partition is ever fully collected. Measured: `PointsModel.parse`'s
# own internal index-monotonicity check alone hit an ArrowMemoryError under a
# 24GB cap with the default (~40-way) scheduler. Capping workers here is a
# blunt, deliberately conservative safety net for the transcript path only —
# trades write speed for a hard bound on concurrent partition memory.
_TRANSCRIPT_WRITE_WORKERS = 2


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


def _carry_over_clusterings(ctx: "ViewerContext", adata_cropped) -> None:
    """Write ctx.clusterings (subset to adata_cropped's cells) into adata_cropped.obs.

    Uses the same 'clustering_<name>' / 'cluster_labels_<name>' column
    convention as save_clustering_to_adata / save_cluster_labels_to_sdata
    (utils/adata_persistence.py), so load_custom_clusterings_from_adata picks
    them up automatically the next time this exported dataset is opened — no
    analysis/clustering/ folder needed. Covers both built-in Xenium
    clusterings (graphclust, kmeans_*) and custom ones (e.g. Leiden), since
    neither is guaranteed to already be written into ctx.adata.obs.
    """
    from xenium_viewer.utils.adata_persistence import CLUSTERING_PREFIX, CLUSTER_LABELS_PREFIX

    if not ctx.clusterings:
        return

    has_cell_id = "cell_id" in adata_cropped.obs.columns
    if has_cell_id:
        cell_id_to_idx = pd.Series(adata_cropped.obs.index, index=adata_cropped.obs["cell_id"].values)

    cluster_labels = ctx.state.get("cluster_labels", {}) if ctx.state else {}

    for name, series in ctx.clusterings.items():
        col = f"{CLUSTERING_PREFIX}{name}"
        if has_cell_id:
            aligned = series.rename(cell_id_to_idx).reindex(adata_cropped.obs.index)
        else:
            aligned = series.reindex(adata_cropped.obs.index)
        adata_cropped.obs[col] = pd.Categorical(aligned)

        label_dict = cluster_labels.get(name)
        if label_dict:
            str_map = {str(k): str(v) for k, v in label_dict.items()}
            adata_cropped.obs[f"{CLUSTER_LABELS_PREFIX}{name}"] = (
                adata_cropped.obs[col].astype(str).map(str_map).fillna("").astype(object)
            )


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
    # Kept lazy (a dask array, not `.compute()`'d here): for a crop spanning a
    # large fraction of the slide, materializing the full-res crop as a dense
    # numpy array before it ever reaches disk can be tens of GB (measured: a
    # 27043x51144x4ch uint16 crop is ~11GB for the image alone, ~5.5GB each for
    # the two label rasters below, plus np.isin/np.unique working buffers on
    # top) — on top of whatever the running viewer already holds, that's a real
    # OOM risk, and it's the only place in the codebase that eagerly
    # materializes an image/label array this size (loader.py's initial cache
    # build writes the full, larger, uncropped image the same dask-backed way).
    # Image2DModel.parse/Labels2DModel.parse accept a dask array directly, so
    # `new_sdata.write(...)` below computes and writes these chunk-by-chunk
    # instead of all at once.
    _progress(20, "Cropping morphology image...")
    morph_scales = _extract_dt_scales(ctx.sdata.images["morphology_focus"])
    morph_full = morph_scales[0]
    cropped_img = morph_full[:, row_min:row_max, col_min:col_max]
    scale_factors = _safe_scale_factors(cropped_img.shape[-2], cropped_img.shape[-1])
    img_element = Image2DModel.parse(
        cropped_img, dims=("c", "y", "x"),
        transformations={"global": Identity()},
        scale_factors=scale_factors,
    )

    # ── Crop cell_labels (bbox, then zero out non-kept cells) ────────────────
    _progress(35, "Cropping cell labels...")
    cl_scales = _extract_dt_scales(ctx.sdata.labels["cell_labels"])
    cropped_cl = cl_scales[0][row_min:row_max, col_min:col_max]
    masked_cl = da.where(da.isin(cropped_cl, kept_cell_ids), cropped_cl, 0)
    cl_element = Labels2DModel.parse(
        masked_cl, dims=("y", "x"),
        transformations={"global": Identity()},
        scale_factors=scale_factors,
    )

    # ── Crop nucleus_labels (bbox, then zero out nuclei outside kept cells) ───
    # nucleus_labels pixel values are their own independent numbering, not the
    # same per-cell IDs as cell_labels/cell_id (verified: sampling both
    # rasters at the same spatial location gives different, unrelated
    # numbers, even though both value ranges span roughly 1..n_obs). So
    # instead of matching ID numbers, find which nucleus IDs spatially
    # overlap a kept cell's footprint in the already-masked cell-label crop
    # (any overlapping pixel counts — a nucleus straddling a kept/non-kept
    # cell boundary is conservatively kept rather than silently dropped) and
    # zero out the rest, mirroring the cell_labels masking above. The only
    # eager `.compute()` in this whole section: its result is just the
    # distinct nucleus IDs (a few thousand values), not the underlying array.
    _progress(45, "Cropping nucleus labels...")
    nl_scales = _extract_dt_scales(ctx.sdata.labels["nucleus_labels"])
    cropped_nl = nl_scales[0][row_min:row_max, col_min:col_max]
    nl_where_kept_cell = da.where(masked_cl > 0, cropped_nl, 0)
    kept_nucleus_ids = da.unique(nl_where_kept_cell).compute()
    kept_nucleus_ids = kept_nucleus_ids[kept_nucleus_ids > 0]
    masked_nl = da.where(da.isin(cropped_nl, kept_nucleus_ids), cropped_nl, 0)
    nl_element = Labels2DModel.parse(
        masked_nl, dims=("y", "x"),
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

    # ctx.sdata.points["transcripts"] is spatialdata's dask dataframe over the
    # on-disk points.parquet, opened with a plain dask.dataframe.read_parquet(path)
    # and no predicate (spatialdata._io.io_points._read_points). A boolean-mask
    # bbox filter applied on top of that dataframe does *not* get pushed down
    # into the parquet reader, and even where an explicit `filters=` kwarg does
    # get real row-group pruning (confirmed for a bbox narrow in x), a crop
    # that spans nearly the full width — e.g. splitting a two-section slide by
    # height, which is what this tissue's own naming ("PCA043_PCA044") and a
    # real crash both point at — barely prunes anything: the on-disk
    # partitions are banded by x, not y, so a full-width bbox still overlaps
    # every partition. Measured directly: for that crop shape, collecting the
    # filtered result into one pandas DataFrame (the previous version of this
    # code) pulled ~100M+ rows into memory and hit an ArrowMemoryError under a
    # 24GB test cap. So this stays fully lazy end-to-end, like the image/label
    # crops above — never call `.compute()` on the filtered dataframe itself;
    # let `new_sdata.write()` and `to_parquet()` below stream it partition by
    # partition. `filters=` is still applied as a first-pass narrowing where it
    # helps (small/x-narrow crops), it just isn't relied on for correctness or
    # memory safety.
    x_min_um, x_max_um = col_min * pixel_size, col_max * pixel_size
    y_min_um, y_max_um = row_min * pixel_size, row_max * pixel_size
    bbox_filters = [
        ("x", ">=", x_min_um), ("x", "<=", x_max_um),
        ("y", ">=", y_min_um), ("y", "<=", y_max_um),
    ]
    try:
        import dask.dataframe as dd
        points_parquet_path = Path(ctx.sdata.path) / "points" / "transcripts" / "points.parquet"
        transcripts_ddf = dd.read_parquet(points_parquet_path, filters=bbox_filters)
    except Exception:
        transcripts_ddf = ctx.sdata.points["transcripts"]
        transcripts_ddf = transcripts_ddf[
            (transcripts_ddf["x"] >= x_min_um) & (transcripts_ddf["x"] <= x_max_um)
            & (transcripts_ddf["y"] >= y_min_um) & (transcripts_ddf["y"] <= y_max_um)
        ]

    def _filter_partition(pdf):
        if len(pdf) == 0:
            return pdf
        mask = contains_xy(poly_xy_um, pdf["x"].to_numpy(), pdf["y"].to_numpy())
        return pdf[mask]

    filtered_ddf = transcripts_ddf.map_partitions(_filter_partition)
    # Each source partition keeps its own index (e.g. restarting at 0), so
    # without this the partitions have duplicate index labels — spatialdata's
    # parquet writer (dask-expr `assign`) fails on that with "cannot reindex
    # on an axis with duplicate labels". reset_index on a dask dataframe is a
    # lightweight per-partition-length computation, not a full materialization.
    filtered_ddf = filtered_ddf.reset_index(drop=True)
    filtered_ddf["x"] = filtered_ddf["x"] - origin_x_um
    filtered_ddf["y"] = filtered_ddf["y"] - origin_y_um
    points_coords = {"x": "x", "y": "y"}
    if "z" in filtered_ddf.columns:
        points_coords["z"] = "z"
    # Everything from here through the write() calls below stays under a
    # concurrency-capped scheduler (see _TRANSCRIPT_WRITE_WORKERS) — including
    # `.cat.as_known()` just below, which triggers its own eager `.compute()`
    # and hit the same unbounded-concurrency ArrowMemoryError as the rest of
    # this section when it ran under the default (~40-way) scheduler.
    with dask.config.set(scheduler="threads", num_workers=_TRANSCRIPT_WRITE_WORKERS):
        # feature_name is dictionary-encoded with ~5000+ possible gene names
        # (this panel: 5001 predesigned + 100 custom, well past 127) — read
        # lazily/filtered, each partition only knows the categories *it*
        # happens to contain, so pyarrow infers a different index width per
        # partition (int8 vs int16) and `to_parquet` fails with a schema
        # mismatch across partitions once two disagree. spatialdata's own
        # `write_points()` hits the identical issue and fixes it the same
        # way: force one *known*, shared category set before any write. This
        # computes only the distinct gene names actually present (a few
        # thousand strings), not the data.
        if filtered_ddf["feature_name"].dtype == "category" and not filtered_ddf["feature_name"].cat.known:
            filtered_ddf["feature_name"] = filtered_ddf["feature_name"].cat.as_known()

        points_element = PointsModel.parse(
            filtered_ddf,
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

    _progress(70, "Carrying over clusterings...")
    _carry_over_clusterings(ctx, adata_cropped)

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

        # Image/label pyramid chunks are small enough that the default (~40-way)
        # scheduler already proved safe for them (measured); the points element
        # embedded in the same `new_sdata` is not, so the whole write stays
        # under the capped scheduler rather than trying to split images/labels
        # out into a separate write just to give them back full concurrency.
        with dask.config.set(scheduler="threads", num_workers=_TRANSCRIPT_WRITE_WORKERS):
            new_sdata.write(str(staging_dir / "sdata_cached.zarr"), overwrite=True)

            _progress(90, "Writing transcripts.parquet...")
            parquet_path = staging_dir / "transcripts.parquet"
            # Re-derives the same filtered dask dataframe rather than reusing a
            # materialized result from the write above (there isn't one to reuse
            # — filtered_ddf was never `.compute()`'d) — this re-reads and
            # re-filters the source transcripts a second time, which is wasted
            # work for a large crop, but keeps this write bounded the same way:
            # each partition is materialized, written, and dropped in turn.
            #
            # This has to be `transcripts.parquet` as a single *file*, not a
            # directory — preprocess.py and TranscriptLoader read it as plain
            # Xenium output, matching the real 10x format. dask's own
            # `.to_parquet()` always writes a directory of per-partition files
            # (like spatialdata's own points.parquet), so it can't be used
            # here; write one shared `pq.ParquetWriter` across partitions
            # instead — same bounded-memory streaming, single physical file.
            raw_schema_ddf = filtered_ddf.rename(columns={"x": "x_location", "y": "y_location", "z": "z_location"})
            raw_schema_ddf.attrs = {}   # spatialdata's transform metadata isn't JSON-serializable
            schema = pa.Table.from_pandas(raw_schema_ddf._meta, preserve_index=False).schema
            writer = pq.ParquetWriter(parquet_path, schema)
            try:
                for partition in raw_schema_ddf.partitions:
                    pdf = partition.compute()
                    if len(pdf) == 0:
                        continue
                    writer.write_table(pa.Table.from_pandas(pdf, preserve_index=False).cast(schema))
            finally:
                writer.close()

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
