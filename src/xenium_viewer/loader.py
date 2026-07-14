"""
SpatialData loading for Xenium 3.x output.

Usage (standalone test):
    python scripts/01_load_sdata.py

Returns a SpatialData object with:
  - images:  morphology_focus (multiscale, 4-channel, CYX)
  - labels:  cell_labels, nucleus_labels
  - points:  transcripts (dask-backed, lazy)
  - tables:  table (AnnData 318K cells x 480+ genes)
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ─── Channel names for the 4 morphology_focus planes ───────────────────────
CHANNEL_NAMES = [
    "DAPI",
    "ATP1A1-CD45-E-Cadherin",
    "18S",
    "AlphaSMA-Vimentin",
]

# ─── User-generated element keys ────────────────────────────────────────────
_USER_SHAPE_KEYS = [
    "rois", "he_xenium_landmarks", "he_he_landmarks",
    "arms_xenium_landmarks", "arms_he_landmarks", "arms_tiles", "annotations",
]
_USER_IMAGE_KEYS = ["he_image", "arms_he_image"]
_USER_UNS_KEYS = ["nhood_enrichment", "co_occurrence", "ligrec", "rank_genes_groups"]
_SIDECAR_FILES = ["roi_deg_cache.parquet", "arms_tile_deg_cache.parquet", "adata_norm_cache.h5ad"]

_SHAPE_LABELS = {
    "rois": "ROIs",
    "he_xenium_landmarks": "H&E landmarks (Xenium side)",
    "he_he_landmarks": "H&E landmarks (H&E image side)",
    "arms_xenium_landmarks": "ARMS landmarks (Xenium side)",
    "arms_he_landmarks": "ARMS landmarks (ARMS image side)",
    "arms_tiles": "ARMS tiles",
    "annotations": "Annotations",
}
_IMAGE_LABELS = {"he_image": "H&E image", "arms_he_image": "ARMS image"}


def _detect_user_data(cache_path: Path) -> dict:
    """Fast filesystem scan for user-added elements without loading sdata.

    Returns a dict with keys: shapes, images, clusterings, uns_keys,
    has_obsm_umap, sidecars, has_viewer_session.
    """
    found: dict = {
        "shapes": [],
        "images": [],
        "clusterings": [],
        "uns_keys": [],
        "has_obsm_umap": False,
        "sidecars": [],
        "has_viewer_session": False,
    }
    for key in _USER_SHAPE_KEYS:
        if (cache_path / "shapes" / key).exists():
            found["shapes"].append(key)
    for key in _USER_IMAGE_KEYS:
        if (cache_path / "images" / key).exists():
            found["images"].append(key)
    obs_dir = cache_path / "tables" / "table" / "obs"
    if obs_dir.exists():
        for item in obs_dir.iterdir():
            if item.name.startswith("clustering_"):
                found["clusterings"].append(item.name)
    uns_dir = cache_path / "tables" / "table" / "uns"
    if uns_dir.exists():
        for key in _USER_UNS_KEYS:
            if (uns_dir / key).exists():
                found["uns_keys"].append(key)
    obsm_dir = cache_path / "tables" / "table" / "obsm"
    if obsm_dir.exists() and (obsm_dir / "X_umap").exists():
        found["has_obsm_umap"] = True
    for fname in _SIDECAR_FILES:
        if (cache_path / fname).exists():
            found["sidecars"].append(fname)
    if (cache_path / "viewer_session").exists():
        found["has_viewer_session"] = True
    return found


def _has_any_user_data(user_data: dict) -> bool:
    return bool(
        user_data["shapes"] or user_data["images"] or
        user_data["clusterings"] or user_data["uns_keys"] or
        user_data["sidecars"] or user_data["has_viewer_session"]
    )


def _format_user_data_message(user_data: dict) -> str:
    lines = []
    for key in user_data["shapes"]:
        lines.append(f"  • {_SHAPE_LABELS.get(key, key)}")
    for key in user_data["images"]:
        lines.append(f"  • {_IMAGE_LABELS.get(key, key)}")
    if user_data["clusterings"]:
        n = len(user_data["clusterings"])
        lines.append(f"  • {n} custom clustering{'s' if n > 1 else ''}")
    if user_data["uns_keys"]:
        labels = {
            "nhood_enrichment": "Neighborhood enrichment results",
            "co_occurrence": "Co-occurrence results",
            "ligrec": "Ligand-receptor results",
            "rank_genes_groups": "Rank genes results",
        }
        for key in user_data["uns_keys"]:
            lines.append(f"  • {labels.get(key, key)}")
    _sidecar_labels = {
        "roi_deg_cache.parquet": "ROI DEG results",
        "arms_tile_deg_cache.parquet": "ARMS tile DEG results",
        "adata_norm_cache.h5ad": "Normalized expression cache",
    }
    for fname in user_data["sidecars"]:
        lines.append(f"  • {_sidecar_labels.get(fname, fname)}")
    return "\n".join(lines)


def _ask_rebuild_preference(user_data: dict) -> str:
    """Show a Qt dialog asking what to do when a stale cache has user data.

    Returns 'restore', 'rebuild', or 'keep'.
    Falls back to 'rebuild' if called from a non-main thread or if Qt is absent.
    """
    import threading
    if threading.current_thread() is not threading.main_thread():
        print(
            "Warning: zarr cache is stale and contains user data, but running in a "
            "background thread — rebuilding without restoring."
        )
        return "rebuild"
    try:
        from qtpy.QtWidgets import QApplication, QMessageBox
    except Exception:
        print("Warning: zarr cache is stale and contains user data (no Qt for dialog). Rebuilding.")
        return "rebuild"
    app = QApplication.instance()
    if app is None:
        print("Warning: zarr cache is stale and contains user data (no Qt app). Rebuilding.")
        return "rebuild"

    summary = _format_user_data_message(user_data)
    msg = QMessageBox()
    msg.setWindowTitle("Zarr Cache Needs Rebuilding")
    msg.setText(
        "experiment.xenium has been modified since the cache was last built.\n\n"
        "Your cache contains user-generated data:\n"
        f"{summary}\n\n"
        "What would you like to do?"
    )
    msg.setIcon(QMessageBox.Warning)
    restore_btn = msg.addButton("Rebuild and restore my data", QMessageBox.AcceptRole)
    msg.addButton("Rebuild without restoring", QMessageBox.DestructiveRole)
    keep_btn = msg.addButton("Keep existing cache", QMessageBox.RejectRole)
    msg.setDefaultButton(restore_btn)
    msg.exec_()

    clicked = msg.clickedButton()
    if clicked is restore_btn:
        return "restore"
    elif clicked is keep_btn:
        return "keep"
    else:
        return "rebuild"


def _restore_user_elements(old_sdata, sdata, user_data: dict) -> list[str]:
    """Merge user-added elements from old_sdata into sdata in memory.

    Returns a list of successfully restored element names.
    """
    restored = []
    for key in user_data["shapes"]:
        try:
            if key in old_sdata.shapes:
                sdata[key] = old_sdata.shapes[key]
                restored.append(_SHAPE_LABELS.get(key, key))
        except Exception as e:
            print(f"  Warning: could not restore shape '{key}': {e}")
    for key in user_data["images"]:
        try:
            if key in old_sdata.images:
                sdata[key] = old_sdata.images[key]
                restored.append(_IMAGE_LABELS.get(key, key))
        except Exception as e:
            print(f"  Warning: could not restore image '{key}': {e}")
    if "table" in old_sdata.tables and "table" in sdata.tables:
        old_adata = old_sdata["table"]
        new_adata = sdata["table"]
        user_obs_cols = [
            c for c in old_adata.obs.columns
            if c.startswith("clustering_") or c.startswith("cluster_labels_")
        ]
        for col in user_obs_cols:
            try:
                new_adata.obs[col] = old_adata.obs[col].reindex(new_adata.obs.index)
                if col.startswith("clustering_"):
                    restored.append(col)
            except Exception as e:
                print(f"  Warning: could not restore obs column '{col}': {e}")
        for key in user_data["uns_keys"]:
            if key in old_adata.uns:
                try:
                    new_adata.uns[key] = old_adata.uns[key]
                    restored.append(f"uns/{key}")
                except Exception as e:
                    print(f"  Warning: could not restore uns['{key}']: {e}")
        if user_data["has_obsm_umap"] and "X_umap" in old_adata.obsm:
            try:
                if len(old_adata.obsm["X_umap"]) == len(new_adata.obs):
                    new_adata.obsm["X_umap"] = old_adata.obsm["X_umap"]
                    restored.append("UMAP coordinates")
            except Exception as e:
                print(f"  Warning: could not restore UMAP coordinates: {e}")
    return restored


def _copy_sidecars_and_session(backup_path: Path, cache_path: Path, user_data: dict) -> None:
    """Copy sidecar files and viewer_session zarr group from backup to new cache."""
    import shutil
    for fname in user_data.get("sidecars", []):
        src = backup_path / fname
        dst = cache_path / fname
        if src.exists() and not dst.exists():
            try:
                shutil.copy2(str(src), str(dst))
            except Exception as e:
                print(f"  Warning: could not restore {fname}: {e}")
    if user_data.get("has_viewer_session"):
        src = backup_path / "viewer_session"
        dst = cache_path / "viewer_session"
        if src.exists() and not dst.exists():
            try:
                shutil.copytree(str(src), str(dst))
            except Exception as e:
                print(f"  Warning: could not restore viewer_session: {e}")


def _convert_arrow_strings(sdata):
    """
    Convert ArrowStringArray columns in sdata's AnnData table to object dtype.

    pandas 3.0 with future.infer_string=True uses PyArrow-backed string arrays
    by default, but anndata's zarr writer has no serializer for ArrowStringArray.
    We must temporarily disable infer_string and rebuild obs/var DataFrames with
    plain object dtype so sdata.write() can serialize the table.
    """
    adata = sdata["table"]
    old_infer = pd.options.future.infer_string
    pd.options.future.infer_string = False
    try:
        for attr in ["obs", "var"]:
            df = getattr(adata, attr).copy()
            if pd.api.types.is_string_dtype(df.index):
                df.index = pd.Index(df.index.to_numpy(dtype=object))
            for col in df.columns:
                if isinstance(df[col].dtype, pd.CategoricalDtype):
                    # Rebuild categorical with object-dtype categories
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


def load_sdata(
    path: Path,
    build_pyramid: bool = True,
    n_jobs: int = 8,
    use_cache: bool = True,
):
    """
    Load the Xenium 3.x output as a SpatialData object.

    Parameters
    ----------
    path : Path
        Root directory of the Xenium output.
    build_pyramid : bool
        If True, build a 5-level software image pyramid for the morphology_focus
        TIFFs (which have no internal pyramid). Required for smooth napari
        pan/zoom performance.
    n_jobs : int
        Number of threads for spatialdata_io.
    use_cache : bool
        If True, use/create a zarr cache for faster subsequent loads.

    Returns
    -------
    spatialdata.SpatialData
    """
    import shutil
    from datetime import datetime

    cache_path = path / "sdata_cached.zarr"
    experiment_path = path / "experiment.xenium"
    backup_path = None   # set when we move the old cache for a user-restore
    user_data = None     # populated when we detect user data in the old cache

    # A directory with no raw Xenium files (e.g. one exported by the Crop
    # Dataset tool — just experiment.xenium + sdata_cached.zarr + transcripts)
    # can only ever be loaded from its zarr cache; there's nothing to rebuild
    # from. Use the cache even if the caller asked for use_cache=False (e.g.
    # launched with --no-cache), rather than failing on a missing cells.zarr.zip.
    if not use_cache and cache_path.exists() and not (path / "cells.zarr.zip").exists():
        print(
            "No raw Xenium files found in this directory (likely a Crop Dataset "
            "export) — loading from the zarr cache regardless of --no-cache."
        )
        use_cache = True

    # Try loading from zarr cache if it exists and is fresh
    if use_cache and cache_path.exists():
        cache_fresh = True
        if experiment_path.exists():
            cache_mtime = cache_path.stat().st_mtime
            exp_mtime = experiment_path.stat().st_mtime
            if exp_mtime > cache_mtime:
                cache_fresh = False

        if cache_fresh:
            import spatialdata
            print(f"Loading SpatialData from zarr cache: {cache_path}")
            try:
                sdata = spatialdata.read_zarr(str(cache_path))
                print("SpatialData loaded from cache.")
                print(sdata)
                return sdata
            except Exception as e:
                # Cache is unreadable — preserve it so the user can recover data,
                # then rebuild from raw Xenium files.
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                corrupt_dest = cache_path.with_name(f"sdata_cached_corrupt_{timestamp}.zarr")
                try:
                    shutil.move(str(cache_path), str(corrupt_dest))
                    print(
                        f"Warning: zarr cache is corrupt ({e}).\n"
                        f"Corrupt cache preserved at:\n  {corrupt_dest}\n"
                        "You may be able to recover data from it manually.\n"
                        "Rebuilding cache from raw Xenium files..."
                    )
                except Exception:
                    shutil.rmtree(cache_path, ignore_errors=True)
                    print(f"Warning: zarr cache is corrupt ({e}). Deleted, rebuilding...")
        else:
            # Stale cache: check whether it contains user-generated data.
            user_data = _detect_user_data(cache_path)
            if _has_any_user_data(user_data):
                preference = _ask_rebuild_preference(user_data)
                if preference == "keep":
                    import spatialdata
                    print("Using existing cache (stale, kept at user request).")
                    sdata = spatialdata.read_zarr(str(cache_path))
                    print(sdata)
                    return sdata
                elif preference == "restore":
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    backup_path = cache_path.with_name(f"sdata_cached_backup_{timestamp}.zarr")
                    shutil.move(str(cache_path), str(backup_path))
                    print(
                        f"Old cache backed up to:\n  {backup_path}\n"
                        "Rebuilding cache and restoring user data..."
                    )
                else:
                    print("Rebuilding cache without restoring user data...")
            else:
                print("Zarr cache is stale (experiment.xenium is newer). Rebuilding...")

    from spatialdata_io import xenium

    image_models_kwargs = {}
    if build_pyramid:
        image_models_kwargs = {"scale_factors": [2, 2, 2, 2, 2]}

    print(f"Loading SpatialData from {path} ...")
    sdata = xenium(
        path=path,
        cells_boundaries=False,
        nucleus_boundaries=False,
        cells_labels=True,
        nucleus_labels=True,
        transcripts=True,
        morphology_focus=True,
        morphology_mip=False,   # not in Xenium 3.x output
        cells_table=True,
        n_jobs=n_jobs,
        image_models_kwargs=image_models_kwargs,
    )
    print("SpatialData loaded successfully.")
    print(sdata)

    # Restore user elements from backup into the fresh sdata (before writing).
    if backup_path is not None:
        import spatialdata as _sd
        try:
            print("Reading user data from backup...")
            old_sdata = _sd.read_zarr(str(backup_path))
            restored = _restore_user_elements(old_sdata, sdata, user_data)
            if restored:
                print(f"  Restored: {', '.join(restored)}")
            else:
                print("  Nothing to restore.")
        except Exception as e:
            print(
                f"Warning: could not read backup for restoration ({e}).\n"
                f"Backup preserved at:\n  {backup_path}"
            )
            backup_path = None  # keep backup; don't delete later

    # Write zarr cache for next time
    if use_cache:
        try:
            _convert_arrow_strings(sdata)
            print(f"Writing zarr cache to {cache_path} ...")
            sdata.write(str(cache_path), overwrite=True)
            print("Zarr cache written.")
            # Copy sidecar files and viewer_session from backup into the new cache.
            if backup_path is not None:
                _copy_sidecars_and_session(backup_path, cache_path, user_data)
                shutil.rmtree(str(backup_path), ignore_errors=True)
                print("Cache rebuild and data restoration complete.")
        except Exception as e:
            shutil.rmtree(cache_path, ignore_errors=True)
            print(f"Warning: could not write zarr cache: {e}")
            if backup_path is not None:
                print(f"Backup preserved at:\n  {backup_path}")

    return sdata


def load_umap(path: Path):
    """
    Load precomputed UMAP coordinates.

    Returns a DataFrame with columns ['UMAP_1', 'UMAP_2'] indexed by cell barcode.
    Note: the UMAP has 91 fewer cells than the AnnData — handled with reindex.
    """
    analysis_path = path / "analysis"
    umap_path = analysis_path / "umap" / "gene_expression_2_components" / "projection.csv"
    if not umap_path.exists():
        raise FileNotFoundError(f"UMAP projection not found at {umap_path}")
    umap_df = pd.read_csv(umap_path, index_col=0)
    umap_df.columns = ["UMAP_1", "UMAP_2"]
    print(f"Loaded UMAP: {umap_df.shape[0]} cells")
    return umap_df


def load_clusterings(path: Path):
    """
    Load all cluster assignments from analysis/clustering/.

    Returns a dict: {clustering_name -> pd.Series(cluster_id, index=cell_barcode)}
    """
    analysis_path = path / "analysis"
    clustering_root = analysis_path / "clustering"
    clusterings = {}
    for subdir in sorted(clustering_root.iterdir()):
        if not subdir.is_dir():
            continue
        csv_path = subdir / "clusters.csv"
        if not csv_path.exists():
            continue
        # strip leading "gene_expression_" prefix for display
        name = subdir.name.replace("gene_expression_", "")
        df = pd.read_csv(csv_path, index_col=0)
        clusterings[name] = df.iloc[:, 0]  # first column = cluster id
    print(f"Loaded {len(clusterings)} clusterings: {list(clusterings.keys())}")
    return clusterings


def get_label_to_obs_mapping(sdata):
    """
    Build a mapping from integer label value -> AnnData obs index position.

    Xenium cell labels are stored as integer raster masks where pixel value k
    corresponds to cell k. The AnnData obs DataFrame has a 'cell_id' or
    similar column linking to label values.

    Returns
    -------
    np.ndarray of shape (max_label + 1,) where arr[k] = row index in adata.obs
    or -1 if no cell has label k.
    """
    adata = sdata["table"]
    # spatialdata_io stores the label value in obs under the region column
    # Look for a column that links obs rows to label values
    region_col = None
    for col in ["cell_labels", "cell_id", "label"]:
        if col in adata.obs.columns:
            region_col = col
            break

    if region_col is None:
        # Fall back: assume obs are ordered 1..N matching label 1..N
        print("Warning: no label-obs link column found; assuming sequential ordering")
        n = len(adata.obs)
        label_to_obs = np.arange(-1, n, dtype=np.int32)  # index 0 = -1 (background)
        return label_to_obs

    label_values = adata.obs[region_col].values.astype(np.int32)
    max_label = int(label_values.max())
    label_to_obs = np.full(max_label + 1, -1, dtype=np.int32)
    for obs_idx, lv in enumerate(label_values):
        if lv > 0:
            label_to_obs[lv] = obs_idx
    return label_to_obs


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir",
                        help="Path to Xenium output directory")
    parser.add_argument("--no-cache", action="store_true",
                        help="Skip zarr cache, load from raw output")
    args = parser.parse_args()

    data_path = Path(args.data_dir)

    sdata = load_sdata(data_path, use_cache=not args.no_cache)
    umap_df = load_umap(data_path)
    clusterings = load_clusterings(data_path)

    # Print summary
    print("\n=== SpatialData summary ===")
    print(f"Images:  {list(sdata.images.keys())}")
    print(f"Labels:  {list(sdata.labels.keys())}")
    if hasattr(sdata, 'shapes') and sdata.shapes:
        print(f"Shapes:  {list(sdata.shapes.keys())}")
    print(f"Points:  {list(sdata.points.keys())}")
    print(f"Tables:  {list(sdata.tables.keys())}")

    adata = sdata["table"]
    print(f"\nAnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"UMAP:    {umap_df.shape[0]} cells")
    print(f"Clusterings: {list(clusterings.keys())}")


if __name__ == "__main__":
    main()
