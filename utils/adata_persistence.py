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
    from utils.viewer_context import ViewerContext

_persist_lock = threading.Lock()

# Prefix for custom clustering columns in adata.obs
CLUSTERING_PREFIX = "clustering_"


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


def _persist_table(ctx: ViewerContext) -> None:
    """Write adata back to the zarr-backed sdata store.

    No-op if running in no-cache mode or sdata has no backing store.
    Thread-safe (serialized via lock).
    """
    if ctx.no_cache:
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
            print(f"Warning: could not persist adata table: {e}")


# ── Save functions ────────────────────────────────────────────────────────

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
    """Write adata_norm to sdata.tables['adata_norm'].

    No-op if running in no-cache mode or sdata has no backing store.
    Thread-safe (serialized via lock).
    """
    if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
        return
    with _persist_lock:
        try:
            _convert_adata_arrow_strings(adata_norm)
            if 'adata_norm' in ctx.sdata:
                ctx.sdata.delete_element_from_disk('adata_norm')
            ctx.sdata['adata_norm'] = adata_norm
            ctx.sdata.write_element('adata_norm')
        except Exception as e:
            print(f"Warning: could not persist adata_norm: {e}")


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
    """Read rank genes results from adata.uns and adata_norm from sdata.tables.

    Returns (df, adata_norm, groupby) or (None, None, None) if not stored.
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
    if sdata is not None:
        try:
            adata_norm = sdata['adata_norm'] if 'adata_norm' in sdata else None
        except Exception:
            pass

    return df, adata_norm, groupby


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

    try:
        if 'rois' in ctx.sdata:
            ctx.sdata.delete_element_from_disk('rois')

        if not roi_data:
            return

        polys = [Polygon(arr[:, ::-1]) for arr in roi_data]  # yx → xy
        gdf = gpd.GeoDataFrame(geometry=polys)
        ctx.sdata['rois'] = gdf
        ctx.sdata.write_element('rois')
    except Exception as e:
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
                    # Migrate adata_norm to sdata.tables['adata_norm']
                    if 'adata_norm' not in sdata:
                        try:
                            _convert_adata_arrow_strings(adata_norm)
                            sdata['adata_norm'] = adata_norm
                            sdata.write_element('adata_norm')
                            print("  Migrated rank genes and adata_norm to sdata")
                        except Exception as e:
                            print(f"  Warning: could not migrate adata_norm to sdata: {e}")
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
