"""
AnnData persistence — save/load analysis results in native adata locations.

Phase 1 of the SpatialData storage refactoring. Moves clusterings,
nhood enrichment, co-occurrence, L-R results, and UMAP coordinates
from custom viewer_session/ storage into adata.obs/obsm/uns, persisted
via sdata.write_element("table", overwrite=True).
"""

from __future__ import annotations

import threading
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData
    from spatialdata import SpatialData
    from xenium_viewer.utils.viewer_context import ViewerContext

_persist_lock = threading.Lock()

# Prefix for custom clustering columns in adata.obs
CLUSTERING_PREFIX = "clustering_"

# SpatialData element keys for custom segmentation
CUSTOM_LABELS_KEY = "custom_cell_labels"
CUSTOM_TABLE_KEY = "custom_table"

_permission_dialog_shown = False


def _maybe_show_permission_dialog(e: Exception, operation: str = "data") -> None:
    """Show a QMessageBox if e is a read-only / permission-denied error."""
    global _permission_dialog_shown
    is_perm = (
        isinstance(e, PermissionError)
        or (isinstance(e, OSError) and getattr(e, "errno", None) in (13, 30))
        or "Permission denied" in str(e)
        or "read-only" in str(e).lower()
        or "Read-only" in str(e)
    )
    if not is_perm or _permission_dialog_shown:
        return
    _permission_dialog_shown = True
    try:
        from qtpy.QtWidgets import QApplication, QMessageBox
        if QApplication.instance() is None:
            return
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("Write Permission Error")
        msg.setText(f"Could not save {operation} — permission denied.")
        msg.setInformativeText(
            "The dataset folder may be read-only (e.g. a shared or mounted drive).\n\n"
            "To enable saving, copy the dataset to a writable location and reopen it, "
            "or launch with --no-cache to skip zarr persistence entirely."
        )
        msg.exec_()
    except Exception:
        pass


def _convert_adata_arrow_strings(adata) -> None:
    """Convert ArrowStringArray columns in an AnnData object to object dtype.

    pandas 3.0 with future.infer_string=True uses PyArrow-backed string arrays
    by default, but anndata's zarr writer has no serializer for ArrowStringArray.
    """
    old_infer = pd.options.future.infer_string
    pd.options.future.infer_string = False
    try:
        for attr in ["obs", "var"]:
            df = getattr(adata, attr).copy()
            if pd.api.types.is_string_dtype(df.index):
                df.index = pd.Index(df.index.to_numpy(dtype=object))
            for col in df.columns:
                if isinstance(df[col].dtype, pd.CategoricalDtype):
                    cat = df[col].cat
                    if pd.api.types.is_string_dtype(cat.categories):
                        new_cats = cat.categories.astype(object)
                        df[col] = df[col].cat.rename_categories(
                            dict(zip(cat.categories, new_cats))
                        )
                elif pd.api.types.is_string_dtype(df[col]):
                    df[col] = df[col].to_numpy(dtype=object)
            setattr(adata, attr, df)
    finally:
        pd.options.future.infer_string = old_infer


def _convert_arrow_strings(sdata) -> None:
    """Convert ArrowStringArray columns in sdata's main AnnData table to object dtype."""
    _convert_adata_arrow_strings(sdata["table"])


def _persist_custom_table(ctx: ViewerContext) -> None:
    """Write custom adata to sdata.tables['custom_table']. Caller must hold _persist_lock."""
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        from spatialdata.models import TableModel
        adata_copy = ctx.adata.copy()
        _convert_adata_arrow_strings(adata_copy)
        adata_copy = TableModel.parse(
            adata_copy,
            region=CUSTOM_LABELS_KEY,
            region_key="region",
            instance_key="cell_id",
        )
        sdata = ctx.sdata
        if CUSTOM_TABLE_KEY in sdata.tables:
            sdata.delete_element_from_disk(CUSTOM_TABLE_KEY)
        sdata.tables[CUSTOM_TABLE_KEY] = adata_copy
        sdata.write_element(CUSTOM_TABLE_KEY)
    except Exception as e:
        _maybe_show_permission_dialog(e, "clustering data")
        print(f"Warning: could not persist custom table to sdata: {e}")


def _persist_table(ctx: ViewerContext) -> None:
    """Write adata back to the zarr-backed sdata store.

    When custom segmentation is active, delegates to _persist_custom_table
    so all per-action saves (clustering, rank genes, etc.) go to the right table.
    No-op if running in no-cache mode or sdata has no backing store.
    Thread-safe (serialized via lock).
    """
    if ctx.no_cache:
        return
    if getattr(ctx, "segmentation_source", "xenium") == "custom":
        with _persist_lock:
            _persist_custom_table(ctx)
        return
    sdata = ctx.sdata
    if sdata is None or sdata.path is None:
        return
    with _persist_lock:
        try:
            _convert_arrow_strings(sdata)
            # spatialdata refuses write_element(overwrite=True) when the
            # element is backed by the same zarr store.  The documented
            # workaround is delete-then-write.
            sdata.delete_element_from_disk("table")
            sdata.write_element("table")
        except Exception as e:
            _maybe_show_permission_dialog(e, "clustering/analysis data")
            print(f"Warning: could not persist adata table: {e}")


# Prefix for cluster label columns in adata.obs
CLUSTER_LABELS_PREFIX = "cluster_labels_"

# ── Save functions ────────────────────────────────────────────────────────

def save_cluster_labels_to_sdata(ctx: "ViewerContext", clustering_key: str, label_dict: dict) -> None:
    """Map cluster IDs → labels and persist as adata.obs['cluster_labels_<key>'].

    Each cell gets the label for its cluster ID. This column can be read back
    in a standalone Python session without the viewer.

    The source clustering column may be stored as `clustering_key` directly
    (built-in Xenium clusterings) or as `clustering_<clustering_key>` (custom
    clusterings added via save_clustering_to_adata).
    """
    if ctx.no_cache or ctx.sdata is None:
        return
    adata = ctx.adata
    if adata is None:
        return
    # Resolve which obs column holds the cluster IDs for this key
    if clustering_key in adata.obs.columns:
        src_col = clustering_key
    elif f"{CLUSTERING_PREFIX}{clustering_key}" in adata.obs.columns:
        src_col = f"{CLUSTERING_PREFIX}{clustering_key}"
    else:
        return
    try:
        str_map = {str(k): str(v) for k, v in label_dict.items()}
        obs_key = f"{CLUSTER_LABELS_PREFIX}{clustering_key}"
        adata.obs[obs_key] = (
            adata.obs[src_col].astype(str).map(str_map).fillna("").astype(object)
        )
        _persist_table(ctx)
    except Exception as e:
        _maybe_show_permission_dialog(e, "cluster labels")
        print(f"Warning: could not save cluster labels to sdata: {e}")


def load_cluster_labels_from_sdata(sdata) -> dict:
    """Read cluster_labels_* columns from adata.obs and reconstruct labels dict.

    Returns {clustering_key: {cluster_id: label}} or {} if none found.
    Works with both xenium 'table' and custom 'custom_table'.
    """
    if sdata is None:
        return {}
    table_key = "custom_table" if "custom_table" in (getattr(sdata, "tables", {}) or {}) else "table"
    if table_key not in (getattr(sdata, "tables", {}) or {}):
        return {}
    try:
        adata = sdata.tables[table_key]
        result = {}
        for col in adata.obs.columns:
            if not col.startswith(CLUSTER_LABELS_PREFIX):
                continue
            clustering_key = col[len(CLUSTER_LABELS_PREFIX):]
            # Find which obs column holds the cluster IDs
            if clustering_key in adata.obs.columns:
                src_col = clustering_key
            elif f"{CLUSTERING_PREFIX}{clustering_key}" in adata.obs.columns:
                src_col = f"{CLUSTERING_PREFIX}{clustering_key}"
            else:
                continue
            pairs = adata.obs[[src_col, col]].drop_duplicates()
            label_dict = {}
            for _, row in pairs.iterrows():
                cid = row[src_col]
                try:
                    cid = int(cid)
                except (ValueError, TypeError):
                    pass
                label_dict[cid] = str(row[col])
            result[clustering_key] = label_dict
        return result
    except Exception as e:
        print(f"Warning: could not load cluster labels from sdata: {e}")
        return {}


def save_clustering_to_adata(ctx: ViewerContext, name: str, series: pd.Series) -> None:
    """Write a clustering Series to adata.obs and persist.

    The series is indexed by cell_id (barcode), which may differ from
    adata.obs.index. We align via the 'cell_id' column if present.
    """
    col = f"{CLUSTERING_PREFIX}{name}"
    adata = ctx.adata

    if "cell_id" in adata.obs.columns:
        # Build cell_id -> obs_index mapping, then align the series
        cell_id_to_idx = pd.Series(adata.obs.index, index=adata.obs["cell_id"].values)
        aligned = series.rename(cell_id_to_idx).reindex(adata.obs.index)
    else:
        aligned = series.reindex(adata.obs.index)

    adata.obs[col] = pd.Categorical(aligned)
    _persist_table(ctx)


def save_nhood_to_adata(ctx: ViewerContext, result: dict) -> None:
    """Write nhood enrichment result to adata.uns and persist."""
    if result is None or result.get("warning") is not None:
        return
    ctx.adata.uns["nhood_enrichment"] = {
        "zscore": np.asarray(result["zscore"], dtype=np.float64),
        "count": np.asarray(result["count"], dtype=np.float64),
        "clusters": list(result["clusters"]),
    }
    _persist_table(ctx)


def save_co_occurrence_to_adata(ctx: ViewerContext, result: dict) -> None:
    """Write co-occurrence result to adata.uns and persist."""
    if result is None or result.get("warning") is not None:
        return
    ctx.adata.uns["co_occurrence"] = {
        "occ": np.asarray(result["occ"], dtype=np.float64),
        "interval": np.asarray(result["interval"], dtype=np.float64),
        "clusters": list(result["clusters"]),
    }
    _persist_table(ctx)


def save_ligrec_to_adata(ctx: ViewerContext, result: dict) -> None:
    """Write L-R interaction result to adata.uns and persist."""
    if result is None or result.get("warning") is not None:
        return
    means = result.get("means")
    pvalues = result.get("pvalues")
    if means is None or means.empty:
        return
    ctx.adata.uns["ligrec"] = {
        "means": means,
        "pvalues": pvalues if pvalues is not None else pd.DataFrame(),
    }
    _persist_table(ctx)


# ── Phase 2 save/load: rank genes, adata_norm, ROI polygons ──────────────

def _persist_adata_norm(ctx: ViewerContext, adata_norm) -> None:
    """Write adata_norm as h5ad alongside the zarr cache.

    Saves to <zarr_path>/adata_norm_cache.h5ad. Using a plain h5ad rather
    than sdata.tables avoids spatialdata's table validation requirements
    (region/region_key/instance_key), which adata_norm doesn't satisfy.
    No-op if running in no-cache mode or sdata has no backing store.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        norm_path = Path(ctx.sdata.path) / "adata_norm_cache.h5ad"
        _convert_adata_arrow_strings(adata_norm)
        adata_norm.write_h5ad(norm_path)
    except Exception as e:
        _maybe_show_permission_dialog(e, "normalized expression cache")
        print(f"Warning: could not persist adata_norm: {e}")


def save_custom_seg_to_sdata(ctx: ViewerContext, new_adata, scales: list) -> None:
    """Persist custom segmentation labels and AnnData to the sdata zarr cache.

    Writes:
      - sdata.labels["custom_cell_labels"] — multiscale label raster
      - sdata.tables["custom_table"]       — AnnData annotated via TableModel

    No-op if running in no-cache mode or sdata has no backing store.
    Thread-safe (serialized via _persist_lock).
    scales: list of arrays (highest→lowest resolution) as loaded from source zarr.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    with _persist_lock:
        # ── Save label raster ────────────────────────────────────────────
        try:
            from spatialdata.models import LabelsModel
            arr = scales[0]
            if hasattr(arr, "compute"):
                arr = arr.compute()
            arr = np.asarray(arr)
            n = len(scales)
            scale_factors = [2] * (n - 1) if n > 1 else None
            parsed = LabelsModel.parse(arr, dims=("y", "x"), scale_factors=scale_factors)
            sdata = ctx.sdata
            if CUSTOM_LABELS_KEY in sdata.labels:
                sdata.delete_element_from_disk(CUSTOM_LABELS_KEY)
            sdata.labels[CUSTOM_LABELS_KEY] = parsed
            sdata.write_element(CUSTOM_LABELS_KEY)
            print(f"Custom labels saved to sdata.labels['{CUSTOM_LABELS_KEY}']")
        except Exception as e:
            _maybe_show_permission_dialog(e, "custom segmentation labels")
            print(f"Warning: could not save custom labels to sdata: {e}")
        # ── Save AnnData as sdata.tables["custom_table"] ─────────────────
        _persist_custom_table(ctx)


def load_custom_seg_from_sdata(sdata) -> tuple:
    """Load custom segmentation from sdata if previously persisted.

    Returns (adata, scales) where scales is a list of dask arrays suitable
    for viewer.add_labels(), or (None, None) if not fully cached.
    """
    if sdata is None or CUSTOM_LABELS_KEY not in sdata.labels:
        return None, None
    if CUSTOM_TABLE_KEY not in sdata.tables:
        return None, None
    try:
        import re
        adata = sdata.tables[CUSTOM_TABLE_KEY]
        dt = sdata.labels[CUSTOM_LABELS_KEY]

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
                scales.append(ds[next(iter(ds.data_vars))].data)

        if not scales:
            return None, None
        return adata, scales
    except Exception as e:
        print(f"Warning: could not load custom seg from sdata: {e}")
        return None, None


def save_rank_genes_to_adata(ctx: ViewerContext, df, adata_norm, groupby: str) -> None:
    """Write rank genes results to adata.uns and persist.

    Copies rank_genes_groups from adata_norm.uns into main adata.uns (so
    sc.get.rank_genes_groups_df can reconstruct the DataFrame on restore),
    and stores adata_norm as sdata.tables['adata_norm'] for dotplot/volcano.
    """
    if 'rank_genes_groups' in adata_norm.uns:
        ctx.adata.uns['rank_genes_groups'] = adata_norm.uns['rank_genes_groups']
    ctx.adata.uns['rank_genes_groupby'] = groupby
    _persist_table(ctx)
    _persist_adata_norm(ctx, adata_norm)


def load_rank_genes_from_adata(adata, sdata) -> tuple:
    """Read rank genes results from adata.uns and adata_norm from h5ad cache.

    Returns (df, adata_norm, groupby) or (None, None, None) if not stored.
    adata_norm is loaded from <zarr_path>/adata_norm_cache.h5ad if present.
    """
    rgg = adata.uns.get('rank_genes_groups')
    if rgg is None:
        return None, None, None

    groupby = adata.uns.get('rank_genes_groupby')

    import scanpy as sc
    try:
        df = sc.get.rank_genes_groups_df(adata, group=None)
    except Exception as e:
        print(f"Warning: could not reconstruct rank genes DataFrame: {e}")
        return None, None, None

    adata_norm = None
    if sdata is not None and sdata.path is not None:
        norm_path = Path(sdata.path) / "adata_norm_cache.h5ad"
        if norm_path.exists():
            try:
                adata_norm = sc.read_h5ad(norm_path)
            except Exception as e:
                print(f"Warning: could not load adata_norm cache: {e}")

    return df, adata_norm, groupby


def save_cnv_results_to_adata(ctx: ViewerContext, result: dict) -> None:
    """Write CNV inference results to adata.obs/uns and cache adata_cnv as h5ad.

    Mirrors save_rank_genes_to_adata's shape: the derived per-cell score and
    run metadata go into the main adata (so they round-trip on reopen), and
    the full CNV-profile AnnData (with obsm['X_cnv'], gene positions, etc.)
    is cached alongside the zarr store so the chromosome heatmap can be
    regenerated without recomputation.
    """
    adata = ctx.adata
    score = result["cnv_score"]
    if "cell_id" in adata.obs.columns:
        cell_id_to_idx = pd.Series(adata.obs.index, index=adata.obs["cell_id"].values)
        aligned = score.rename(cell_id_to_idx).reindex(adata.obs.index)
    else:
        aligned = score.reindex(adata.obs.index)
    adata.obs["cnv_score"] = aligned.astype(np.float64)

    adata.uns["cnv_run_info"] = {
        "reference_obs_key": result["reference_obs_key"],
        "reference_clustering_name": result.get("reference_clustering_name", ""),
        "reference_categories": list(result["reference_categories"]),
        "cluster_key": result["cluster_key"],
        "cluster_keys": list(result.get("cluster_keys", [result["cluster_key"]])),
        "n_genes_total": int(result["n_genes_total"]),
        "n_genes_mapped": int(result["n_genes_mapped"]),
        "n_windows": int(result["n_windows"]),
        "params": dict(result["params"]),
    }
    _persist_table(ctx)
    _persist_cnv_adata(ctx, result["adata_cnv"])


def _persist_cnv_adata(ctx: ViewerContext, adata_cnv) -> None:
    """Write adata_cnv as h5ad alongside the zarr cache (see _persist_adata_norm)."""
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        cnv_path = Path(ctx.sdata.path) / "adata_cnv_cache.h5ad"
        _convert_adata_arrow_strings(adata_cnv)
        adata_cnv.write_h5ad(cnv_path)
    except Exception as e:
        _maybe_show_permission_dialog(e, "CNV profile cache")
        print(f"Warning: could not persist adata_cnv: {e}")


def load_cnv_results_from_adata(adata: AnnData, sdata) -> "dict | None":
    """Read CNV run info from adata.uns/obs and adata_cnv from the h5ad cache.

    Returns None if no CNV run has been saved.
    """
    info = adata.uns.get("cnv_run_info")
    if info is None:
        return None

    score = None
    if "cnv_score" in adata.obs.columns:
        if "cell_id" in adata.obs.columns:
            score = pd.Series(
                adata.obs["cnv_score"].values,
                index=adata.obs["cell_id"].values,
                name="cnv_score",
            )
        else:
            score = adata.obs["cnv_score"].rename("cnv_score")

    adata_cnv = None
    if sdata is not None and sdata.path is not None:
        cnv_path = Path(sdata.path) / "adata_cnv_cache.h5ad"
        if cnv_path.exists():
            try:
                import scanpy as sc
                adata_cnv = sc.read_h5ad(cnv_path)
            except Exception as e:
                print(f"Warning: could not load adata_cnv cache: {e}")

    cluster_key = info.get("cluster_key")
    cluster_keys = list(info.get("cluster_keys") or ([cluster_key] if cluster_key else []))
    params = dict(info.get("params", {}))
    # Recompute the core-params signature so a subsequent run under the same
    # parameters can keep accumulating resolutions (see _cnv_signature in tab_cnv).
    signature = (
        info.get("reference_clustering_name", ""),
        tuple(info.get("reference_categories", [])),
        params.get("n_neighbors"),
        params.get("smoothing_neighbors"),
        params.get("window_size"),
        params.get("step"),
        params.get("lfc_clip"),
    )
    return {
        "reference_obs_key": info.get("reference_obs_key"),
        "reference_clustering_name": info.get("reference_clustering_name", ""),
        "reference_categories": list(info.get("reference_categories", [])),
        "cluster_key": cluster_key,
        "cluster_keys": cluster_keys,
        "signature": signature,
        "n_genes_total": info.get("n_genes_total"),
        "n_genes_mapped": info.get("n_genes_mapped"),
        "n_windows": info.get("n_windows"),
        "params": params,
        "cnv_score": score,
        "adata_cnv": adata_cnv,
    }


# ── DEG result caches (sidecar parquets at zarr root) ─────────────────────────

ROI_DEG_CACHE = "roi_deg_cache.parquet"
ARMS_DEG_CACHE = "arms_tile_deg_cache.parquet"


def save_roi_deg_to_sdata(ctx: "ViewerContext", df) -> None:
    """Persist ROI DEG DataFrame to <zarr_path>/roi_deg_cache.parquet."""
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    if df is None or df.empty:
        return
    try:
        cache_path = Path(ctx.sdata.path) / ROI_DEG_CACHE
        df.to_parquet(cache_path, index=False)
    except Exception as e:
        _maybe_show_permission_dialog(e, "ROI DEG results")
        print(f"Warning: could not save ROI DEG to sdata: {e}")


def load_roi_deg_from_sdata(sdata) -> "pd.DataFrame | None":
    """Load ROI DEG DataFrame from <zarr_path>/roi_deg_cache.parquet."""
    if sdata is None or sdata.path is None:
        return None
    cache_path = Path(sdata.path) / ROI_DEG_CACHE
    if not cache_path.exists():
        return None
    try:
        return pd.read_parquet(cache_path)
    except Exception as e:
        print(f"Warning: could not load ROI DEG from sdata: {e}")
        return None


def save_arms_tile_deg_to_sdata(ctx: "ViewerContext", df) -> None:
    """Persist ARMS tile DEG DataFrame to <zarr_path>/arms_tile_deg_cache.parquet."""
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    if df is None or df.empty:
        return
    try:
        cache_path = Path(ctx.sdata.path) / ARMS_DEG_CACHE
        df.to_parquet(cache_path, index=False)
    except Exception as e:
        _maybe_show_permission_dialog(e, "ARMS tile DEG results")
        print(f"Warning: could not save ARMS tile DEG to sdata: {e}")


def load_arms_tile_deg_from_sdata(sdata) -> "pd.DataFrame | None":
    """Load ARMS tile DEG DataFrame from <zarr_path>/arms_tile_deg_cache.parquet."""
    if sdata is None or sdata.path is None:
        return None
    cache_path = Path(sdata.path) / ARMS_DEG_CACHE
    if not cache_path.exists():
        return None
    try:
        return pd.read_parquet(cache_path)
    except Exception as e:
        print(f"Warning: could not load ARMS tile DEG from sdata: {e}")
        return None


def save_rois_to_sdata(ctx: ViewerContext, roi_data: list) -> None:
    """Save ROI polygons to sdata.shapes['rois'] as a GeoDataFrame.

    roi_data: list of Nx2 float64 arrays in yx (napari pixel) coordinates.
    Converts to shapely Polygons in xy coords for GeoDataFrame storage.
    No-op if no-cache mode or sdata has no backing store.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return

    from shapely.geometry import Polygon
    import geopandas as gpd
    from spatialdata.models import ShapesModel

    try:
        if 'rois' in ctx.sdata:
            ctx.sdata.delete_element_from_disk('rois')

        if not roi_data:
            return

        polys = [Polygon(arr[:, ::-1]) for arr in roi_data]  # yx → xy
        gdf = ShapesModel.parse(gpd.GeoDataFrame(geometry=polys))
        ctx.sdata['rois'] = gdf
        ctx.sdata.write_element('rois')
    except Exception as e:
        _maybe_show_permission_dialog(e, "ROI shapes")
        print(f"Warning: could not save ROIs to sdata: {e}")


def load_rois_from_sdata(sdata) -> list:
    """Load ROI polygons from sdata.shapes['rois'].

    Returns list of Nx2 float64 numpy arrays in yx (napari pixel) coords.
    Returns [] if no ROIs are stored in sdata.
    """
    if sdata is None or 'rois' not in sdata:
        return []

    try:
        gdf = sdata['rois']
        rois = []
        for geom in gdf.geometry:
            if geom is None or geom.is_empty:
                continue
            coords = np.array(geom.exterior.coords[:-1], dtype=np.float64)  # xy, no closing pt
            rois.append(coords[:, ::-1])  # xy → yx
        return rois
    except Exception as e:
        print(f"Warning: could not load ROIs from sdata: {e}")
        return []


def save_annotations_to_sdata(ctx: "ViewerContext") -> None:
    """Save annotation layer shapes to sdata.shapes['annotations'] as a GeoDataFrame.

    Stores shapely Polygon geometries with an 'annotation_type' column.
    Coordinates are converted from yx-pixels (napari) to xy-pixels for shapely.
    No-op if no-cache mode, sdata has no backing store, or no annotations exist.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    if ctx.annotation_layer is None:
        return

    from shapely.geometry import Polygon
    import geopandas as gpd
    from spatialdata.models import ShapesModel

    try:
        if 'annotations' in ctx.sdata:
            ctx.sdata.delete_element_from_disk('annotations')

        layer = ctx.annotation_layer
        shapes = layer.data
        raw_types = list(layer.properties.get("annotation_type", []))
        # Normalise: convert NaN / non-string entries to empty string
        import math
        types = [
            str(t) if (t is not None and not (isinstance(t, float) and math.isnan(t))) else ""
            for t in raw_types
        ]
        # Pad types list if shorter than shapes (can happen on initial assignment)
        while len(types) < len(shapes):
            types.append("")

        if not shapes:
            return

        polys = [Polygon(arr[:, ::-1]) for arr in shapes]  # yx → xy
        gdf = gpd.GeoDataFrame(
            {"geometry": polys, "annotation_type": types},
        )
        gdf = ShapesModel.parse(gdf)
        ctx.sdata['annotations'] = gdf
        ctx.sdata.write_element('annotations')
    except Exception as e:
        _maybe_show_permission_dialog(e, "annotations")
        print(f"Warning: could not save annotations to sdata: {e}")


def load_annotations_from_sdata(sdata) -> tuple[list, list]:
    """Load annotation shapes from sdata.shapes['annotations'].

    Returns (list_of_yx_arrays, list_of_type_strings).
    Arrays are Nx2 float64 in yx (napari pixel) coordinates.
    Returns ([], []) if no annotations are stored.
    """
    if sdata is None or 'annotations' not in sdata:
        return [], []

    try:
        gdf = sdata['annotations']
        shapes = []
        types = []
        ann_col = gdf['annotation_type'] if 'annotation_type' in gdf.columns else None
        for i, geom in enumerate(gdf.geometry):
            if geom is None or geom.is_empty:
                continue
            coords = np.array(geom.exterior.coords[:-1], dtype=np.float64)  # xy, no closing pt
            shapes.append(coords[:, ::-1])  # xy → yx
            types.append(ann_col.iloc[i] if ann_col is not None else "")
        return shapes, types
    except Exception as e:
        print(f"Warning: could not load annotations from sdata: {e}")
        return [], []


# ── Phase 3: Landmark and ARMS tile sdata persistence ─────────────────────

def _make_landmark_gdf(points_yx):
    """Build a spatialdata-compatible GeoDataFrame of Point shapes from Nx2 yx array."""
    from shapely.geometry import Point
    import geopandas as gpd
    from spatialdata.models import ShapesModel
    pts = np.asarray(points_yx, dtype=np.float64)
    gdf = gpd.GeoDataFrame(
        {"geometry": [Point(float(x), float(y)) for y, x in pts], "radius": 1.0},
    )
    return ShapesModel.parse(gdf)


def save_landmarks_to_sdata(ctx, element_name: str, points_yx) -> None:
    """Save Nx2 yx landmark array to sdata.shapes[element_name] as Point GeoDataFrame.

    points_yx: Nx2 numpy array in yx (napari) coords, or None/empty to clear.
    Converts to shapely Points in xy coords before storing.
    No-op if no-cache mode or sdata has no backing store.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        if element_name in ctx.sdata:
            ctx.sdata.delete_element_from_disk(element_name)
        if points_yx is None or len(points_yx) == 0:
            return
        ctx.sdata[element_name] = _make_landmark_gdf(points_yx)
        ctx.sdata.write_element(element_name)
    except Exception as e:
        _maybe_show_permission_dialog(e, "registration landmarks")
        print(f"Warning: could not save {element_name} to sdata: {e}")


def load_landmarks_from_sdata(sdata, element_name: str):
    """Load landmarks from sdata.shapes[element_name].

    Returns Nx2 float64 array in yx (napari) coords, or None if not present.
    """
    if sdata is None or element_name not in sdata:
        return None
    try:
        gdf = sdata[element_name]
        return np.array([[p.y, p.x] for p in gdf.geometry], dtype=np.float64)
    except Exception as e:
        print(f"Warning: could not load {element_name} from sdata: {e}")
        return None


def save_arms_tiles_to_sdata(ctx, tile_polygons_yx, tile_names, cluster_ids) -> None:
    """Save ARMS tile polygons with metadata to sdata.shapes['arms_tiles'].

    tile_polygons_yx: list of Nx2 yx numpy arrays (napari pixel coords).
    tile_names: list of str tile name identifiers.
    cluster_ids: array-like of int cluster IDs.
    No-op if no-cache mode or sdata has no backing store.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    from shapely.geometry import Polygon
    import geopandas as gpd
    from spatialdata.models import ShapesModel
    try:
        if 'arms_tiles' in ctx.sdata:
            ctx.sdata.delete_element_from_disk('arms_tiles')
        if not tile_polygons_yx:
            return
        polys = [Polygon(arr[:, ::-1]) for arr in tile_polygons_yx]  # yx → xy
        gdf = ShapesModel.parse(gpd.GeoDataFrame({
            'geometry': polys,
            'tile_name': list(tile_names),
            'cluster_id': np.asarray(cluster_ids, dtype=np.int32),
        }))
        ctx.sdata['arms_tiles'] = gdf
        ctx.sdata.write_element('arms_tiles')
    except Exception as e:
        _maybe_show_permission_dialog(e, "ARMS tiles")
        print(f"Warning: could not save ARMS tiles to sdata: {e}")


def load_arms_tiles_from_sdata(sdata) -> tuple:
    """Load ARMS tile polygons from sdata.shapes['arms_tiles'].

    Returns (polygons_yx, tile_names, cluster_ids) where polygons_yx is a list
    of Nx2 yx float64 arrays in napari pixel coords. Returns ([], [], []) if
    no tiles are stored.
    """
    if sdata is None or 'arms_tiles' not in sdata:
        return [], [], []
    try:
        gdf = sdata['arms_tiles']
        polys = [np.array(p.exterior.coords[:-1])[:, ::-1] for p in gdf.geometry]  # xy → yx
        return polys, list(gdf['tile_name']), np.array(gdf['cluster_id'])
    except Exception as e:
        print(f"Warning: could not load ARMS tiles from sdata: {e}")
        return [], [], []


# ── Load functions ────────────────────────────────────────────────────────

def load_custom_clusterings_from_adata(adata: AnnData) -> dict:
    """Read custom clustering columns from adata.obs.

    Returns dict mapping clustering name -> Series indexed by cell_id
    (to match ctx.clusterings convention).
    """
    clusterings = {}
    has_cell_id = "cell_id" in adata.obs.columns
    for col in adata.obs.columns:
        if col.startswith(CLUSTERING_PREFIX):
            name = col[len(CLUSTERING_PREFIX):]
            raw = adata.obs[col].dropna()
            if len(raw) == 0:
                continue
            if has_cell_id:
                # Re-index from adata.obs.index → cell_id values
                cell_ids = adata.obs.loc[raw.index, "cell_id"].values
                series = pd.Series(raw.values, index=cell_ids, name=name)
            else:
                series = raw
            clusterings[name] = series
    return clusterings


def load_analysis_results_from_adata(adata: AnnData) -> dict:
    """Read nhood/co-occ/ligrec results from adata.uns.

    Returns dict with keys nhood_result, co_result, ligrec_result (each may be None).
    """
    results = {
        "nhood_result": None,
        "co_result": None,
        "ligrec_result": None,
    }

    nh = adata.uns.get("nhood_enrichment")
    if nh is not None and isinstance(nh, dict):
        results["nhood_result"] = {
            "zscore": np.asarray(nh["zscore"]),
            "count": np.asarray(nh["count"]),
            "clusters": list(nh["clusters"]),
            "warning": None,
        }

    co = adata.uns.get("co_occurrence")
    if co is not None and isinstance(co, dict):
        results["co_result"] = {
            "occ": np.asarray(co["occ"]),
            "interval": np.asarray(co["interval"]),
            "clusters": list(co["clusters"]),
            "warning": None,
        }

    lr = adata.uns.get("ligrec")
    if lr is not None and isinstance(lr, dict):
        results["ligrec_result"] = {
            "means": lr.get("means", pd.DataFrame()),
            "pvalues": lr.get("pvalues", pd.DataFrame()),
            "warning": None,
        }

    return results


# ── Migration ─────────────────────────────────────────────────────────────

def migrate_landmarks_to_sdata(zarr_path: "Path", sdata: "SpatialData", session: dict) -> None:
    """One-time migration of zarr-stored landmarks and ARMS tile files → sdata.shapes.

    Reads landmarks from the zarr session arrays and ARMS tiles from the original
    GeoJSON/CSV files (if they still exist), writes them to sdata.shapes, and marks
    the migration done so it never re-runs.

    No-op if already migrated, if sdata has no backing store, or if no data to migrate.
    """
    if sdata is None or sdata.path is None:
        return

    import zarr as _zarr
    try:
        store = _zarr.open_group(str(zarr_path), mode="r+", use_consolidated=False)
    except Exception:
        return
    if "viewer_session" not in store:
        return
    vs = store["viewer_session"]
    if dict(vs.attrs).get("migrated_landmarks_to_sdata"):
        return

    from shapely.geometry import Point, Polygon
    import geopandas as gpd

    migrated_any = False

    # ── Landmarks ────────────────────────────────────────────────────────
    # Source 1: zarr session arrays (snapshot-captured at close time)
    # Source 2: landmarks.json in the dataset folder (saved via Save Landmarks button)
    # Source 2 takes priority for H&E since it's more reliably populated.

    # Try to load from landmarks.json first (zarr_path.parent == data_path)
    lm_json_path = zarr_path.parent / "landmarks.json"
    lm_json = {}
    if lm_json_path.exists():
        try:
            import json as _json
            with open(lm_json_path) as f:
                lm_json = _json.load(f)
            print(f"  Found landmarks.json with {len(lm_json.get('xenium_landmarks_yx', []))} H&E landmark pairs")
        except Exception as e:
            print(f"  Warning: could not read landmarks.json: {e}")

    # Build a unified lookup: prefer json, fall back to zarr session
    def _get_pts(json_key, session_key):
        pts = lm_json.get(json_key) or session.get(session_key)
        return np.asarray(pts, dtype=np.float64) if pts is not None and len(pts) > 0 else None

    _landmark_map = [
        ('he_xenium_landmarks',   'xenium_landmarks_yx',   'xenium_landmarks'),
        ('he_he_landmarks',       'he_landmarks_yx',       'he_landmarks'),
        ('arms_xenium_landmarks', 'arms_xenium_landmarks_yx', 'arms_xenium_landmarks'),
        ('arms_he_landmarks',     'arms_he_landmarks_yx',  'arms_he_landmarks'),
    ]
    for sdata_name, json_key, session_key in _landmark_map:
        if sdata_name in sdata:
            continue  # already present — skip
        pts = _get_pts(json_key, session_key)
        if pts is None:
            continue
        try:
            sdata[sdata_name] = _make_landmark_gdf(pts)
            sdata.write_element(sdata_name)
            migrated_any = True
            print(f"  Migrated {len(pts)} landmarks → sdata.shapes['{sdata_name}']")
        except Exception as e:
            print(f"  Warning: could not migrate {sdata_name}: {e}")

    # ── ARMS tiles — re-parse from GeoJSON + CSV if files still exist ────
    if 'arms_tiles' not in sdata:
        attrs = dict(vs.attrs)
        geojson_path = attrs.get("arms_geojson_path")
        csv_path = attrs.get("arms_csv_path")
        if geojson_path and csv_path:
            gj = Path(geojson_path)
            cp = Path(csv_path)
            if gj.exists() and cp.exists():
                try:
                    import json as _json, csv as _csv, re as _re

                    with open(gj) as f:
                        geojson = _json.load(f)
                    tile_polygons_xy = {}
                    for feat in geojson.get("features", []):
                        name = feat.get("properties", {}).get("name")
                        if name is None:
                            continue
                        coords = feat["geometry"]["coordinates"]
                        gtype = feat["geometry"]["type"]
                        ring = coords[0][0] if gtype == "MultiPolygon" else (coords[0] if gtype == "Polygon" else None)
                        if ring is None:
                            continue
                        tile_polygons_xy[name] = np.array(ring, dtype=np.float64)

                    tile_clusters = {}
                    with open(cp) as f:
                        for row in _csv.DictReader(f):
                            tname = row.get("tile", row.get("sample", "")).strip()
                            cval = row.get("cluster", "").strip()
                            if tname and cval:
                                try:
                                    tile_clusters[tname] = int(cval)
                                except ValueError:
                                    pass

                    def _norm(n):
                        n = n.strip().lower()
                        return _re.sub(r'^plate', 'p', n)

                    csv_norm = {_norm(k): v for k, v in tile_clusters.items()}

                    polys, t_names, c_ids = [], [], []
                    for name, arr_xy in tile_polygons_xy.items():
                        cid = tile_clusters.get(name) or csv_norm.get(_norm(name))
                        if cid is None:
                            continue
                        polys.append(Polygon(arr_xy))  # coords already in xy
                        t_names.append(name)
                        c_ids.append(cid)

                    if polys:
                        from spatialdata.models import ShapesModel as _SM
                        gdf = _SM.parse(gpd.GeoDataFrame({
                            'geometry': polys,
                            'tile_name': t_names,
                            'cluster_id': np.asarray(c_ids, dtype=np.int32),
                        }))
                        sdata['arms_tiles'] = gdf
                        sdata.write_element('arms_tiles')
                        migrated_any = True
                        print(f"  Migrated {len(polys)} ARMS tiles → sdata.shapes['arms_tiles']")
                except Exception as e:
                    print(f"  Warning: could not migrate ARMS tiles: {e}")

    # Mark migration done (even if there was nothing to migrate)
    try:
        vs.attrs["migrated_landmarks_to_sdata"] = True
    except Exception:
        pass



def migrate_old_session_to_adata(
    zarr_path: Path,
    sdata: SpatialData,
    adata: AnnData,
) -> None:
    """One-time migration of old viewer_session data into adata.

    Reads clusterings, nhood, co-occurrence, and ligrec from the old
    viewer_session/ zarr/parquet storage, writes them into adata.obs/uns,
    and persists via sdata.write_element. Sets a migration marker to
    prevent re-running.
    """
    import zarr

    try:
        store = zarr.open_group(str(zarr_path), mode="r+", use_consolidated=False)
    except Exception:
        return

    if "viewer_session" not in store:
        return

    session = store["viewer_session"]
    attrs = dict(session.attrs)

    # Already migrated?
    if attrs.get("migrated_to_adata"):
        return

    migrated_any = False

    # ── Clusterings ───────────────────────────────────────────────────
    clust_dir = Path(zarr_path) / "viewer_session" / "clusterings"
    custom_names = attrs.get("custom_clustering_names", [])
    for name in custom_names:
        col = f"{CLUSTERING_PREFIX}{name}"
        if col in adata.obs.columns:
            continue  # already present
        pq = clust_dir / f"{name}.parquet"
        if pq.exists():
            df = pd.read_parquet(pq)
            series = df.iloc[:, 0]  # indexed by cell_id (barcode)
            # Align cell_id-indexed series to adata.obs.index
            if "cell_id" in adata.obs.columns:
                cell_id_to_idx = pd.Series(adata.obs.index, index=adata.obs["cell_id"].values)
                aligned = series.rename(cell_id_to_idx).reindex(adata.obs.index)
            else:
                aligned = series.reindex(adata.obs.index)
            adata.obs[col] = pd.Categorical(aligned)
            migrated_any = True
            print(f"  Migrated clustering '{name}' to adata.obs")

    # ── Nhood enrichment ──────────────────────────────────────────────
    if attrs.get("has_nhood") and "nhood" in session and "nhood_enrichment" not in adata.uns:
        nhood_group = session["nhood"]
        zscore = np.array(nhood_group["zscore"]) if "zscore" in nhood_group else None
        count = np.array(nhood_group["count"]) if "count" in nhood_group else None
        clusters = list(dict(nhood_group.attrs).get("clusters", []))
        if zscore is not None and zscore.size > 0:
            adata.uns["nhood_enrichment"] = {
                "zscore": zscore,
                "count": count if count is not None else np.zeros_like(zscore),
                "clusters": clusters,
            }
            migrated_any = True
            print("  Migrated nhood enrichment to adata.uns")

    # ── Co-occurrence ─────────────────────────────────────────────────
    if attrs.get("has_co") and "co_occ" in session and "co_occurrence" not in adata.uns:
        co_group = session["co_occ"]
        occ = np.array(co_group["occ"]) if "occ" in co_group else None
        interval = np.array(co_group["interval"]) if "interval" in co_group else None
        clusters = list(dict(co_group.attrs).get("clusters", []))
        if occ is not None and occ.size > 0:
            adata.uns["co_occurrence"] = {
                "occ": occ,
                "interval": interval if interval is not None else np.array([]),
                "clusters": clusters,
            }
            migrated_any = True
            print("  Migrated co-occurrence to adata.uns")

    # ── Ligrec ────────────────────────────────────────────────────────
    if attrs.get("has_ligrec") and "ligrec" not in adata.uns:
        session_dir = Path(zarr_path) / "viewer_session"
        means_path = session_dir / "ligrec_means.parquet"
        pvals_path = session_dir / "ligrec_pvalues.parquet"
        means = pd.read_parquet(means_path) if means_path.exists() else None
        pvals = pd.read_parquet(pvals_path) if pvals_path.exists() else None
        if means is not None:
            adata.uns["ligrec"] = {
                "means": means,
                "pvalues": pvals if pvals is not None else pd.DataFrame(),
            }
            migrated_any = True
            print("  Migrated L-R results to adata.uns")

    # ── Persist Phase 1 and mark ──────────────────────────────────────
    if migrated_any:
        try:
            _convert_arrow_strings(sdata)
            sdata.delete_element_from_disk("table")
            sdata.write_element("table")
            print("  Migration persisted to zarr store")
        except Exception as e:
            print(f"  Warning: migration persist failed: {e}")

    # Mark migration done (even if no data to migrate)
    try:
        session.attrs["migrated_to_adata"] = True
    except Exception:
        pass

    # ── Phase 2 migration: rank genes parquet/h5ad → adata.uns + sdata ────
    if not attrs.get("migrated_rank_genes_to_adata"):
        session_dir = Path(zarr_path) / "viewer_session"
        rg_parquet = session_dir / "rank_genes.parquet"
        rg_h5ad = session_dir / "rank_genes_adata_norm.h5ad"
        groupby = attrs.get("rank_genes_groupby")

        if attrs.get("has_rank_genes") and rg_parquet.exists() and "rank_genes_groups" not in adata.uns:
            # Migrate rank_genes_groups from h5ad into main adata.uns
            if rg_h5ad.exists():
                try:
                    import scanpy as sc
                    adata_norm = sc.read_h5ad(rg_h5ad)
                    if 'rank_genes_groups' in adata_norm.uns:
                        adata.uns['rank_genes_groups'] = adata_norm.uns['rank_genes_groups']
                    if groupby:
                        adata.uns['rank_genes_groupby'] = groupby
                    # Persist updated main adata
                    _convert_arrow_strings(sdata)
                    sdata.delete_element_from_disk("table")
                    sdata.write_element("table")
                    # Migrate adata_norm to h5ad cache alongside the zarr store
                    norm_cache = Path(zarr_path) / "adata_norm_cache.h5ad"
                    if not norm_cache.exists():
                        try:
                            _convert_adata_arrow_strings(adata_norm)
                            adata_norm.write_h5ad(norm_cache)
                            print("  Migrated rank genes and adata_norm to zarr cache")
                        except Exception as e:
                            print(f"  Warning: could not migrate adata_norm: {e}")
                    # Clean up old files
                    try:
                        rg_parquet.unlink(missing_ok=True)
                        rg_h5ad.unlink(missing_ok=True)
                    except Exception:
                        pass
                except Exception as e:
                    print(f"  Warning: rank genes migration failed: {e}")

        try:
            session.attrs["migrated_rank_genes_to_adata"] = True
        except Exception:
            pass

    # ── Phase N migration: ROI DEG + ARMS DEG parquets → zarr root ───────
    if not attrs.get("migrated_deg_to_sdata"):
        session_dir = Path(zarr_path) / "viewer_session"
        for old_name, new_name in [
            ("roi_deg.parquet", ROI_DEG_CACHE),
            ("arms_tile_deg.parquet", ARMS_DEG_CACHE),
        ]:
            old_path = session_dir / old_name
            new_path = Path(zarr_path) / new_name
            if old_path.exists() and not new_path.exists():
                try:
                    import shutil
                    shutil.copy2(old_path, new_path)
                    print(f"  Migrated {old_name} → zarr root")
                except Exception as e:
                    print(f"  Warning: could not migrate {old_name}: {e}")
        try:
            session.attrs["migrated_deg_to_sdata"] = True
        except Exception:
            pass


# ── External multichannel image persistence ────────────────────────────────

def _slugify(text: str) -> str:
    """Make a filesystem-safe slug from a free-form name."""
    import re
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", str(text)).strip("_")
    return s or "image"


def save_external_image_to_sdata(
    ctx, element_name: str, pyramid, channel_axis, channel_names,
) -> None:
    """Persist an external multichannel image into ``sdata.images[element_name]``.

    Uses ``Image2DModel.parse`` like the H&E / ARMS loaders so the image
    travels with the zarr store. The base level of ``pyramid`` is written;
    downstream levels are reconstructed via ``scale_factors`` (powers of two).
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        from spatialdata.models import Image2DModel

        # Remove existing element if possible (may fail if backed by active files)
        if element_name in ctx.sdata:
            try:
                del ctx.sdata[element_name]
            except Exception:
                pass

        import dask.array as da

        base = pyramid[0]
        # Wrap via map_blocks to detach from ZarrTiffStore backing.
        # ZarrTiffStore (tifffile) lacks the .root attribute that zarr v3 /
        # spatialdata expects; map_blocks creates a new graph without exposing
        # the store to spatialdata's introspection.
        if hasattr(base, "dask"):
            base = da.map_blocks(lambda x: x, base, dtype=base.dtype)
        # Build axes string — Image2DModel wants ('c', 'y', 'x') or ('y', 'x', 'c').
        if channel_axis is None:
            base = base[None, ...]
            dims = ("c", "y", "x")
        else:
            if channel_axis != 0:
                base = da.moveaxis(base, channel_axis, 0)
            dims = ("c", "y", "x")

        n_levels = len(pyramid)
        scale_factors = [2] * max(0, n_levels - 1) if n_levels > 1 else None

        parsed = Image2DModel.parse(
            base, dims=dims, scale_factors=scale_factors,
            c_coords=list(channel_names) if channel_names else None,
        )
        # Remove in-memory ref to avoid conflicts, then write with overwrite
        if element_name in ctx.sdata:
            try:
                del ctx.sdata[element_name]
            except Exception:
                pass
        ctx.sdata[element_name] = parsed
        ctx.sdata.write_element(element_name, overwrite=True)
    except Exception as e:
        _maybe_show_permission_dialog(e, f"external image '{element_name}'")
        print(f"Warning: could not save external image {element_name} to sdata: {e}")


def load_external_images_from_sdata(sdata, prefix: str = "ext_") -> list:
    """Return entries for every ``sdata.images`` element whose name starts with ``prefix``.

    Each entry is ``{"element_name": str, "pyramid": list[dask], "channel_names": list[str]}``.
    """
    out = []
    if sdata is None:
        return out
    try:
        for name in list(sdata.images.keys()):
            if not name.startswith(prefix):
                continue
            try:
                img = sdata.images[name]
                pyramid, channel_names = _extract_image_pyramid(img)
                out.append({
                    "element_name": name,
                    "pyramid": pyramid,
                    "channel_names": channel_names,
                    "affine_matrix": _load_affine_from_sdata_element(sdata, name),
                })
            except Exception as e:
                print(f"Warning: could not read external image {name}: {e}")
    except Exception:
        pass
    return out


def _extract_image_pyramid(img):
    """Return (list_of_dask_levels, channel_names) from a spatialdata image element."""
    import dask.array as da
    try:
        from datatree import DataTree  # optional
    except Exception:
        DataTree = None
    import xarray as xr

    pyramid = []
    channel_names: list[str] = []

    if DataTree is not None and isinstance(img, DataTree):
        scales = sorted(img.children.keys())
        for s in scales:
            node = img[s]
            arr = node[list(node.data_vars)[0]] if hasattr(node, "data_vars") else node.to_array()
            pyramid.append(da.asarray(arr.data))
            if not channel_names and "c" in arr.coords:
                channel_names = [str(c) for c in arr.coords["c"].values]
    elif isinstance(img, xr.DataArray):
        pyramid.append(da.asarray(img.data))
        if "c" in img.coords:
            channel_names = [str(c) for c in img.coords["c"].values]
    else:
        pyramid.append(da.asarray(img))

    if not channel_names and pyramid:
        n = pyramid[0].shape[0] if pyramid[0].ndim >= 3 else 1
        channel_names = [f"C{i}" for i in range(n)]

    return pyramid, channel_names


# ── Patch overlay persistence ───────────────────────────────────────────────

def save_patch_overlay_to_sdata(
    ctx,
    element_name: str,
    coords_xy,
    patch_size: int,
    cluster_columns: dict,
    confidence=None,
) -> None:
    """Store a patch-grid overlay as a Point GeoDataFrame in ``sdata.shapes``.

    Columns: ``radius`` (= patch_size/2) plus one integer column per cluster
    labelling in ``cluster_columns``. Optional ``confidence`` is included if
    provided. Top-left patch (x, y) coordinates are encoded as shapely Points.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    try:
        from shapely.geometry import Point
        import geopandas as gpd
        from spatialdata.models import ShapesModel

        if element_name in ctx.sdata:
            ctx.sdata.delete_element_from_disk(element_name)

        coords = np.asarray(coords_xy, dtype=np.float64)
        if coords.size == 0:
            return

        data = {
            "geometry": [Point(float(x), float(y)) for x, y in coords],
            "radius": float(patch_size) / 2.0,
            "patch_size": np.full(len(coords), int(patch_size), dtype=np.int32),
        }
        for col, vals in cluster_columns.items():
            data[col] = np.asarray(vals, dtype=np.int64)
        if confidence is not None:
            data["confidence"] = np.asarray(confidence, dtype=np.float32)

        gdf = ShapesModel.parse(gpd.GeoDataFrame(data))
        ctx.sdata[element_name] = gdf
        ctx.sdata.write_element(element_name)
    except Exception as e:
        _maybe_show_permission_dialog(e, f"patch overlay '{element_name}'")
        print(f"Warning: could not save patch overlay {element_name} to sdata: {e}")


def save_overlay_affine_to_sdata(ctx, element_name: str, affine_matrix) -> None:
    """Write an affine transformation to an sdata element (image or shape).

    Uses ``spatialdata.transformations.set_transformation`` +
    ``sdata.write_transformations`` — the same pattern as
    ``_save_he_affine_to_sdata`` in ``tab_he_registration.py``.
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    if element_name not in ctx.sdata:
        return
    try:
        from spatialdata.transformations import (
            Affine as SdAffine, set_transformation,
        )
        m = np.asarray(affine_matrix, dtype=np.float64)
        # SpatialData expects a 3×3 matrix with input/output axes (y, x).
        if m.shape == (3, 3):
            sd_affine = SdAffine(m, input_axes=("y", "x"), output_axes=("y", "x"))
        else:
            return  # unexpected shape — skip
        elem = ctx.sdata[element_name]
        set_transformation(elem, sd_affine, "global")
        ctx.sdata.write_transformations(element_name)
    except Exception as e:
        _maybe_show_permission_dialog(e, f"affine transform for '{element_name}'")
        print(f"Warning: could not save affine for {element_name}: {e}")


def _load_affine_from_sdata_element(sdata, element_name: str):
    """Return the 3×3 affine matrix stored on an sdata element, or None."""
    if sdata is None or element_name not in sdata:
        return None
    try:
        from spatialdata.transformations import get_transformation, Affine as SdAffine
        elem = sdata[element_name]
        t = get_transformation(elem, "global")
        if isinstance(t, SdAffine):
            m = np.asarray(t.to_affine_matrix(input_axes=("y", "x"),
                                               output_axes=("y", "x")),
                           dtype=np.float64)
            if not np.allclose(m, np.eye(3), atol=1e-6):
                return m
    except Exception:
        pass
    return None


def load_patch_overlays_from_sdata(sdata, prefix: str = "patch_") -> list:
    """Return entries for every ``sdata.shapes`` element whose name starts with ``prefix``.

    Each entry is a dict::

        {
          "element_name": str,
          "coords_xy": np.ndarray (N, 2),
          "patch_size": int,
          "cluster_columns": {col: np.ndarray},
          "confidence": np.ndarray | None,
        }
    """
    out = []
    if sdata is None:
        return out
    try:
        shape_names = list(sdata.shapes.keys())
    except Exception:
        return out
    for name in shape_names:
        if not name.startswith(prefix):
            continue
        try:
            gdf = sdata.shapes[name]
            coords = np.array(
                [[float(p.x), float(p.y)] for p in gdf.geometry],
                dtype=np.float64,
            )
            # Patch size is stored per-point; all rows share the same value.
            if "patch_size" in gdf.columns and len(gdf) > 0:
                patch_size = int(gdf["patch_size"].iloc[0])
            elif "radius" in gdf.columns and len(gdf) > 0:
                patch_size = int(round(float(gdf["radius"].iloc[0]) * 2))
            else:
                patch_size = 0

            reserved = {"geometry", "radius", "patch_size", "confidence"}
            cluster_columns = {
                col: np.asarray(gdf[col].values, dtype=np.int64)
                for col in gdf.columns if col not in reserved
            }
            confidence = None
            if "confidence" in gdf.columns:
                confidence = np.asarray(gdf["confidence"].values, dtype=np.float32)

            out.append({
                "element_name": name,
                "coords_xy": coords,
                "patch_size": patch_size,
                "cluster_columns": cluster_columns,
                "confidence": confidence,
                "affine_matrix": _load_affine_from_sdata_element(sdata, name),
            })
        except Exception as e:
            print(f"Warning: could not load patch overlay {name}: {e}")
    return out
