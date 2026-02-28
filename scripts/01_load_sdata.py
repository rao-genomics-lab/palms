"""
SpatialData loading for Xenium 3.x output.

Usage (standalone test):
    python scripts/01_load_sdata.py

Returns a SpatialData object with:
  - images:  morphology_focus (multiscale, 4-channel, CYX)
  - labels:  cell_labels, nucleus_labels
  - shapes:  cell_boundaries, nucleus_boundaries
  - points:  transcripts (dask-backed, lazy)
  - tables:  table (AnnData 318K cells × 480+ genes)
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA_PATH = Path(__file__).parent.parent  # /media/.../output-XETG...
ANALYSIS_PATH = DATA_PATH / "analysis"

# ─── Channel names for the 4 morphology_focus planes ───────────────────────
CHANNEL_NAMES = [
    "DAPI",
    "ATP1A1-CD45-E-Cadherin",
    "18S",
    "AlphaSMA-Vimentin",
]


def load_sdata(
    path: Path = DATA_PATH,
    build_pyramid: bool = True,
    n_jobs: int = 8,
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

    Returns
    -------
    spatialdata.SpatialData
    """
    from spatialdata_io import xenium

    image_models_kwargs = {}
    if build_pyramid:
        image_models_kwargs = {"scale_factors": [2, 2, 2, 2, 2]}

    print(f"Loading SpatialData from {path} ...")
    sdata = xenium(
        path=path,
        cells_boundaries=True,
        nucleus_boundaries=True,
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
    return sdata


def load_umap(sdata=None):
    """
    Load precomputed UMAP coordinates.

    Returns a DataFrame with columns ['UMAP_1', 'UMAP_2'] indexed by cell barcode.
    Note: the UMAP has 91 fewer cells than the AnnData — handled with reindex.
    """
    umap_path = ANALYSIS_PATH / "umap" / "gene_expression_2_components" / "projection.csv"
    if not umap_path.exists():
        raise FileNotFoundError(f"UMAP projection not found at {umap_path}")
    umap_df = pd.read_csv(umap_path, index_col=0)
    umap_df.columns = ["UMAP_1", "UMAP_2"]
    print(f"Loaded UMAP: {umap_df.shape[0]} cells")
    return umap_df


def load_clusterings():
    """
    Load all cluster assignments from analysis/clustering/.

    Returns a dict: {clustering_name -> pd.Series(cluster_id, index=cell_barcode)}
    """
    clustering_root = ANALYSIS_PATH / "clustering"
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


if __name__ == "__main__":
    sdata = load_sdata()
    umap_df = load_umap()
    clusterings = load_clusterings()

    # Print summary
    print("\n=== SpatialData summary ===")
    print(f"Images:  {list(sdata.images.keys())}")
    print(f"Labels:  {list(sdata.labels.keys())}")
    print(f"Shapes:  {list(sdata.shapes.keys())}")
    print(f"Points:  {list(sdata.points.keys())}")
    print(f"Tables:  {list(sdata.tables.keys())}")

    adata = sdata["table"]
    print(f"\nAnnData: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"UMAP:    {umap_df.shape[0]} cells")
    print(f"Clusterings: {list(clusterings.keys())}")
