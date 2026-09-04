"""
Crop the loaded dataset to a user-drawn polygon and write it out as a
standalone, independently-openable palms data directory.

Images/labels are cropped to the polygon's pixel bounding box. Cells and
transcripts are filtered to the exact polygon via true point-in-polygon
tests (mirrors the technique used by tabs/tab_roi.py), so the exported
table/points are precise even though the raster extends slightly beyond
the drawn shape.
"""

from __future__ import annotations

import logging
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

from palms.utils import sdata_write, xenium_specs

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext

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
# 24GB cap with the default (~40-way) scheduler.
#
# Shared with the cache build in `loader.py`, which has exactly the same problem
# for exactly the same reason — the two drifting apart is how one write path ends
# up bounded and the other does not.
_TRANSCRIPT_WRITE_WORKERS = sdata_write.WRITE_WORKERS


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
    from palms.utils.adata_persistence import CLUSTERING_PREFIX, CLUSTER_LABELS_PREFIX

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
        # `any(...)`, not just truthiness: a clustering whose names were all
        # cleared still has a dict, and writing it produced a
        # `cluster_labels_<name>` column of empty strings — which reads back as a
        # named clustering with 24 blank names.
        if label_dict and any(str(v).strip() for v in label_dict.values()):
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


def crop_export_note(name: str, output_dir, include_overlays: bool, notes) -> tuple:
    """The `crop_export:<name>` provenance node, as ``(id, code, label)``.

    One function because the node is built twice — once into the graph the export
    carries, before the session is written, and once into the source dataset's own
    graph after the worker returns. Two constructions of the same string is the
    drift `run_step` exists to prevent, and here the two would be in different
    files.
    """
    joined = "; ".join(notes) if notes else ""
    carried = ("with every registered overlay and drawn region, cropped to "
               "the same region" if include_overlays else "core elements only")
    code = (
        f"\n# Crop & Export dataset '{name}' (writes a standalone dataset)\n"
        f"# from palms.utils.crop_export import crop_and_export\n"
        f"# crop_and_export(ctx, polygon_yx=<region drawn in the \"Crop Regions\" layer>,\n"
        f"#                 output_dir=Path({str(output_dir)!r}), name={name!r},\n"
        f"#                 include_overlays={include_overlays!r})\n"
        f"# -> wrote standalone dataset to \"{Path(output_dir) / name}\" ({carried})"
        + (f"\n# -> not carried: {joined}" if joined else "")
    )
    return f"crop_export:{name}", code, f"Crop & export: {name}"


def _unresolved_note(elem_name, shape_yx, xenium_shape_yx) -> str:
    """Why an overlay was not sliced, in terms the user can act on."""
    return (
        f"{elem_name}: frame unknown — the element declares an identity transform and "
        f"its extent ({shape_yx[0]}x{shape_yx[1]}) is not the morphology grid "
        f"({xenium_shape_yx[0]}x{xenium_shape_yx[1]}), and no registration for it was "
        "found in the viewer. Not sliced, because slicing it would take a rectangle of "
        "the wrong picture. Register it in the source dataset and re-export."
    ) if xenium_shape_yx else (
        f"{elem_name}: frame unknown and the morphology grid could not be read. Not sliced."
    )


def _collect_overlays(ctx, crop_bbox, poly_xy, progress,
                      viewer_state=None) -> tuple[dict, dict, dict, list]:
    """Crop every registered overlay and drawn region the store holds.

    Returns ``(images, labels, shapes, notes)`` ready to hand to ``SpatialData``.
    A single overlay that cannot be carried is skipped and named in *notes*
    rather than failing the export: the core dataset is the thing the user asked
    for, and losing it over one unregistered stray image would be the worse
    trade. The report is what keeps that from being silent.

    ``viewer_state`` is an :class:`~palms.utils.crop_state.OverlayFrames`
    captured on the GUI thread. Where it names a frame, that frame wins over the
    element's stored transformation — see ``crop_state`` for why the element is
    not trustworthy on its own.

    Rasters are done before shapes because a shapes element written in a raster's
    pixels needs that raster's slice origin, which only exists once it is cropped.
    """
    from palms.utils import crop_overlays
    from palms.utils.crop_state import OverlayFrames, capture_overlay_frames

    if viewer_state is None:
        viewer_state = capture_overlay_frames(ctx)
    if not isinstance(viewer_state, OverlayFrames):    # pragma: no cover - defensive
        viewer_state = OverlayFrames(frames={}, companions={}, xenium_shape_yx=None)

    frames = viewer_state.frames
    companions = viewer_state.companions
    xenium_shape = viewer_state.xenium_shape_yx

    images, labels, shapes, notes = {}, {}, {}, []
    origins: dict = {}
    names = crop_overlays.user_overlay_names(ctx.sdata)
    total = sum(len(v) for v in names.values())
    if not total:
        return images, labels, shapes, notes

    def _frame_for(elem_name, element):
        """live layer -> stored element transform -> identity."""
        if elem_name in frames:
            return np.asarray(frames[elem_name], dtype=np.float64)
        return crop_overlays.element_affine(element)

    done = 0
    for group in ("images", "labels", "shapes"):
        for elem_name in names[group]:
            done += 1
            progress(int(76 + 6 * done / total), f"Cropping overlay '{elem_name}'...")
            try:
                element = ctx.sdata[elem_name]
                if group == "shapes":
                    companion = companions.get(elem_name)
                    if companion is not None and companion in origins:
                        frame = crop_overlays.Frame(
                            _frame_for(companion, ctx.sdata[companion]),
                            origins[companion])
                    elif companion is not None and companion in frames:
                        frame = crop_overlays.Frame(
                            np.asarray(frames[companion], dtype=np.float64), None)
                    elif elem_name in frames:
                        frame = crop_overlays.Frame(
                            np.asarray(frames[elem_name], dtype=np.float64), None)
                    else:
                        frame = None
                    out = crop_overlays.crop_vector_overlay(
                        element, elem_name, crop_bbox, poly_xy, frame=frame)
                    if out is None:
                        notes.append(f"{elem_name}: nothing inside the crop")
                        continue
                    shapes[elem_name] = out
                else:
                    frame_affine = _frame_for(elem_name, element)
                    shape_yx = _raster_shape_yx(element)
                    if shape_yx is not None and not crop_overlays.frame_is_credible(
                            frame_affine, shape_yx, xenium_shape):
                        # Not a failure — an honest "I do not know where this goes".
                        # Slicing on a guess is the one irreversible move here.
                        notes.append(_unresolved_note(elem_name, shape_yx, xenium_shape))
                        continue
                    out = crop_overlays.crop_raster_overlay(
                        element, crop_bbox, _safe_scale_factors,
                        frame_affine=frame_affine)
                    if out is None:
                        notes.append(f"{elem_name}: does not overlap the crop")
                        continue
                    (labels if out.is_labels else images)[elem_name] = out.element
                    origins[elem_name] = out.origin_yx
            except Exception as e:
                log.warning("could not carry overlay %r through the crop: %s", elem_name, e)
                notes.append(f"{elem_name}: skipped ({e})")
    return images, labels, shapes, notes


def _write_session(ctx, staging_dir, name, output_dir, include_overlays,
                   overlay_notes, new_sdata, adata_cropped) -> None:
    """Give the export a `viewer_session` and its own copy of the graph.

    Failing here must not fail the export — the dataset itself is written and
    valid, and the session is recoverable by opening it once. So this reports
    and returns rather than raising.
    """
    from palms.utils import crop_session

    try:
        carried = set(new_sdata.images) | set(new_sdata.labels) | set(new_sdata.shapes)

        labels = {}
        try:
            from palms.utils.adata_persistence import load_cluster_labels_from_sdata
            labels = load_cluster_labels_from_sdata(new_sdata) or {}
        except Exception:
            log.debug("could not read cluster labels for the export session", exc_info=True)

        graph_items = []
        graph = (getattr(ctx, "state", None) or {}).get("prov_graph")
        if graph is not None and len(graph):
            export_path = str(Path(output_dir) / name)
            graph_items = crop_session.rewrite_graph_paths(
                graph.to_list(), str(ctx.data_path), export_path)
            # The rewrite repoints every recorded path at the export, including
            # ones under directories the export does not carry (10x's
            # analysis/clustering CSVs are the case that bites). Existence is
            # checked against the staging directory — what is actually here —
            # rather than against a list of what an export is believed to hold.
            graph_items, dangling = crop_session.repair_missing_reads(
                graph_items, export_path, staging_dir,
                # Written after this call, by design: see the comment on the
                # _write_session call site.
                also_expected=("transcripts.parquet", "transcript_cache"))
            for node_id, missing in dangling:
                rel = Path(missing).name
                overlay_notes.append(
                    f"provenance: '{node_id}' reads {rel}, which this export does "
                    "not contain — that cell will fail on replay"
                )
            node_id, code, label = crop_export_note(
                name, output_dir, include_overlays, overlay_notes)
            graph_items = [g for g in graph_items if g.get("id") != node_id]
            graph_items.append({
                "id": node_id, "code": code, "deps": ["preamble"], "kind": "note",
                "label": label, "params": {}, "stale": False,
                "seq": len(graph_items) + 1,
                "template_id": None, "template_origin": "builtin", "template_hash": None,
            })

        attrs = crop_session.build_session_attrs(
            ctx,
            carried_elements=carried,
            cluster_labels=labels,
            graph_items=graph_items,
            he_shape_yx=_raster_shape_yx(new_sdata.images.get("he_image")),
            arms_shape_yx=_raster_shape_yx(new_sdata.images.get("arms_he_image")),
            roi_count=len(new_sdata.shapes["rois"]) if "rois" in new_sdata.shapes else 0,
        )
        crop_session.write_export_session(staging_dir, attrs, graph_items)
    except Exception as e:
        log.warning("could not write the export's session for %r: %s", name, e)
        overlay_notes.append(
            f"session/provenance not written ({e}) — the exported data is complete; "
            "open it once and close it to write one"
        )


def _raster_shape_yx(element):
    """``(y, x)`` of a raster element's full-resolution level, or None."""
    if element is None:
        return None
    try:
        scales = _extract_dt_scales(element)
        if not scales:
            return None
        shape = scales[0].shape
        return tuple(int(v) for v in (shape if len(shape) == 2 else shape[-2:]))
    except Exception:                                  # pragma: no cover - defensive
        return None


def crop_and_export(
    ctx: "ViewerContext",
    polygon_yx: np.ndarray,
    output_dir: Path,
    name: str,
    progress_cb: Callable[[int, str], None] | None = None,
    include_overlays: bool = True,
    viewer_state=None,
) -> Path:
    """Crop the loaded dataset to *polygon_yx* and write a standalone
    palms data directory at ``output_dir / name``.

    With *include_overlays* (the default), every registered overlay and drawn
    region travels with the crop — see ``utils/crop_overlays.py``.

    *viewer_state* is an ``OverlayFrames`` captured on the GUI thread by
    ``crop_state.capture_overlay_frames``. Passing it is how a caller running on
    a worker thread avoids reading napari layers from that thread; omitting it
    makes this function capture its own, which is right for a direct or test
    call but not from inside ``_CropExportWorker``.

    Read-only with respect to *ctx* — safe to call from a background thread.
    Writes to a hidden staging directory first and only moves it into place
    on full success, so a failure never leaves a partial dataset on disk.

    Returns ``(final_path, overlay_notes)`` — the second is one line per overlay
    that could *not* be carried, empty when everything travelled. It is returned
    rather than only logged because a partly-populated export has to be
    describable to the person who asked for it.
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
    from palms.loader import _convert_arrow_strings, write_manifest
    from palms.utils.zarr_safe import atomic_json
    from palms.preprocess import preprocess

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

    # ── Registered overlays and drawn regions ────────────────────────────────
    # A crop that keeps only the five core elements throws away the registered
    # H&E, the ARMS tiles, the ROIs — the work that made the dataset worth
    # cropping in the first place. See utils/crop_overlays.py for how each is
    # carried across, and utils/crop_state.py for how its coordinate frame is
    # decided — the element's own transformation is a fallback, not the
    # authority, because a registration can live only in the viewer.
    overlay_images, overlay_labels, overlay_shapes, overlay_notes = {}, {}, {}, []
    if include_overlays:
        _progress(76, "Cropping registered overlays...")
        overlay_images, overlay_labels, overlay_shapes, overlay_notes = _collect_overlays(
            ctx, (row_min, row_max, col_min, col_max), poly_xy, _progress,
            viewer_state=viewer_state,
        )

    # ── Assemble and write to a staging directory ────────────────────────────
    _progress(82, "Writing dataset...")
    new_sdata = SpatialData(
        images={"morphology_focus": img_element, **overlay_images},
        labels={"cell_labels": cl_element, "nucleus_labels": nl_element, **overlay_labels},
        points={"transcripts": points_element},
        shapes=overlay_shapes,
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
        # Write experiment.xenium first so its mtime predates the zarr write —
        # loader.py's cache-freshness check requires the zarr be >= as fresh.
        #
        # Restated, never copied: every quantity in the parent's file describes
        # the parent (measured on demo_data/crop_6: num_cells 299769 for a
        # 76,577-cell crop), and this is the one Xenium-format file the export
        # writes, so it is the one a reader would trust. It also gives the crop
        # a source hash of its own — loader's freshness fingerprint is a sha256
        # of this file, so a copied one hashed identically to the parent's.
        # utils/xenium_specs.py owns which keys are carried, restated, or
        # demoted into the embedded copy of the run's own file.
        atomic_json(
            staging_dir / "experiment.xenium",
            xenium_specs.crop_specs(
                xenium_specs.read_specs(ctx.data_path / "experiment.xenium"),
                stats=xenium_specs.crop_table_stats(adata_cropped.obs),
                # The drawn polygon in microns — already built above for the
                # transcript filter, and the crop's analogue of region_area.
                region_area_um2=poly_xy_um.area,
                source_path=ctx.data_path,
            ),
        )

        # Image/label pyramid chunks are small enough that the default (~40-way)
        # scheduler already proved safe for them (measured); the points element
        # embedded in the same `new_sdata` is not, so the whole write stays
        # under the capped scheduler rather than trying to split images/labels
        # out into a separate write just to give them back full concurrency.
        with dask.config.set(scheduler="threads", num_workers=_TRANSCRIPT_WRITE_WORKERS):
            new_sdata.write(str(staging_dir / "sdata_cached.zarr"), overwrite=True)

            # The session and the provenance graph, so the export opens as a
            # dataset rather than as a bag of elements. Written after the store
            # exists and before the manifest, so a failure here leaves staging
            # to be discarded rather than a half-labelled export in place.
            _progress(86, "Writing session and provenance...")
            _write_session(ctx, staging_dir, name, output_dir, include_overlays,
                           overlay_notes, new_sdata, adata_cropped)

            # Stamp the manifest here rather than leaving it to the first load.
            # Without one, freshness falls back to comparing experiment.xenium's
            # mtime against the cache directory's — and an export that is copied,
            # unzipped or synced has those timestamps reordered, which used to
            # send the loader into a rebuild with no raw files to rebuild from.
            # `cache_only` says the same thing declaratively, for readers that
            # would otherwise have to infer it from absent files.
            write_manifest(
                staging_dir / "sdata_cached.zarr",
                staging_dir / "experiment.xenium",
                extra={"cache_only": True, "derived_from": str(ctx.data_path)},
            )

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
    if overlay_notes:
        log.info("crop %r: %d overlay(s) not carried: %s",
                 name, len(overlay_notes), "; ".join(overlay_notes))
    return final_dir, overlay_notes
