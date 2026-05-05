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
    cache_path = path / "sdata_cached.zarr"
    experiment_path = path / "experiment.xenium"

    # Try loading from zarr cache if it exists and is fresh
    if use_cache and cache_path.exists():
        cache_fresh = True
        if experiment_path.exists():
            cache_mtime = cache_path.stat().st_mtime
            exp_mtime = experiment_path.stat().st_mtime
            if exp_mtime > cache_mtime:
                print("Zarr cache is stale (experiment.xenium is newer). Rebuilding...")
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
                import shutil
                print(f"Warning: zarr cache is corrupt ({e}). Deleting and rebuilding...")
                shutil.rmtree(cache_path, ignore_errors=True)

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

    # Write zarr cache for next time
    if use_cache:
        try:
            _convert_arrow_strings(sdata)
            print(f"Writing zarr cache to {cache_path} ...")
            sdata.write(str(cache_path), overwrite=True)
            print("Zarr cache written.")
        except Exception as e:
            import shutil
            shutil.rmtree(cache_path, ignore_errors=True)
            print(f"Warning: could not write zarr cache: {e}")

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
