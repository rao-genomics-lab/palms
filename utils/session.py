"""
Session persistence for the Xenium viewer.

Saves and loads viewer state (ROIs, H&E registration, ARMS overlay,
rank genes, ROI DEG, cluster labels) to/from a zarr store so the
session can be restored on the next launch.

NOTE: Clusterings, nhood enrichment, co-occurrence, L-R results, and
UMAP coordinates are now persisted in adata.obs/obsm/uns via
utils/adata_persistence.py (Phase 1 SpatialData refactoring).

Storage layout inside sdata_cached.zarr/viewer_session/:
  rois/0, rois/1, ...         — zarr arrays (Nx2 float64) for each ROI polygon
  he/affine_3x3               — zarr array 3x3
  he/coarse_affine             — zarr array 3x3
  he/xenium_landmarks          — zarr array Nx2
  he/he_landmarks              — zarr array Nx2
  arms/affine_3x3              — zarr array 3x3
  arms/xenium_landmarks        — zarr array Nx2
  arms/he_landmarks            — zarr array Nx2
  rank_genes.parquet           — DataFrame
  roi_deg.parquet              — DataFrame
  (group attrs contain JSON metadata)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import zarr


def _write_array(group, name, data):
    """Write a numpy array to a zarr group (compatible with zarr v2 and v3)."""
    arr = np.asarray(data, dtype=np.float64)
    ds = group.create_array(name, shape=arr.shape, dtype=arr.dtype)
    ds[:] = arr


def save_session(
    zarr_path: Path,
    state: dict,
    he_state: dict,
    snapshot: dict,
):
    """
    Save viewer session state to zarr store.

    Parameters
    ----------
    zarr_path : Path
        Path to sdata_cached.zarr directory.
    state : dict
        Viewer state dict (_state from _build_control_panel).
    he_state : dict
        H&E registration state (_he_state from _build_control_panel).
        Must include 'flip_v' and 'flip_h' bool keys.
    snapshot : dict
        Data captured from napari layers before Qt teardown.
        Keys: 'roi_data' (list of arrays), 'xenium_landmarks' (array|None),
        'he_landmarks' (array|None).
    """
    try:
        store = zarr.open_group(str(zarr_path), mode="r+", use_consolidated=False)

        # Clean up parquet files from previous session (not tracked by zarr)
        session_dir = Path(zarr_path) / "viewer_session"
        if session_dir.exists():
            for pq in session_dir.glob("*.parquet"):
                pq.unlink()
            # Clean up custom clustering subdirectory
            clust_dir = session_dir / "clusterings"
            if clust_dir.exists():
                for pq in clust_dir.glob("*.parquet"):
                    pq.unlink()

        # Preserve real-time-saved ARMS attrs before overwrite wipes them
        _prev_arms_attrs = {}
        if "viewer_session" in store:
            prev = store["viewer_session"].attrs
            for key in ("arms_he_filename", "arms_he_path", "arms_he_shape_yx",
                        "arms_affine_3x3", "arms_flip_v", "arms_flip_h",
                        "arms_geojson_path", "arms_csv_path"):
                if key in prev and prev[key] is not None:
                    _prev_arms_attrs[key] = prev[key]

        session = store.create_group("viewer_session", overwrite=True)
        attrs = {}

        # ── ROIs ──────────────────────────────────────────────────────────
        # ROIs are now persisted to sdata.shapes['rois'] by save_rois_to_sdata()
        # called from 02_xenium_viewer.py at exit. Keep roi_count attr for
        # legacy reads and to avoid breaking old zarr session structures.
        roi_data = snapshot.get("roi_data", [])
        attrs["roi_count"] = len(roi_data)

        # ── H&E registration ─────────────────────────────────────────────
        he_group = session.create_group("he")

        if he_state.get("affine_3x3") is not None:
            _write_array(he_group, "affine_3x3", he_state["affine_3x3"])

        if he_state.get("coarse_affine") is not None:
            _write_array(he_group, "coarse_affine", he_state["coarse_affine"])

        # Landmarks from snapshot (captured before Qt teardown)
        xen_lm = snapshot.get("xenium_landmarks")
        he_lm = snapshot.get("he_landmarks")
        if xen_lm is not None:
            _write_array(he_group, "xenium_landmarks", xen_lm)
        if he_lm is not None:
            _write_array(he_group, "he_landmarks", he_lm)

        attrs["he_filename"] = he_state.get("he_filename")
        attrs["he_path"] = he_state.get("he_path")
        attrs["he_shape_yx"] = (
            list(he_state["he_shape_yx"]) if he_state.get("he_shape_yx") else None
        )
        attrs["flip_v"] = bool(he_state.get("flip_v", False))
        attrs["flip_h"] = bool(he_state.get("flip_h", False))

        # ── ARMS overlay ─────────────────────────────────────────────────
        arms_state = snapshot.get("arms_state", {})
        arms_group = session.create_group("arms")

        if arms_state.get("affine_3x3") is not None:
            _write_array(arms_group, "affine_3x3", arms_state["affine_3x3"])

        arms_xen_lm = snapshot.get("arms_xenium_landmarks")
        arms_he_lm = snapshot.get("arms_he_landmarks")
        if arms_xen_lm is not None:
            _write_array(arms_group, "xenium_landmarks", arms_xen_lm)
        if arms_he_lm is not None:
            _write_array(arms_group, "he_landmarks", arms_he_lm)

        attrs["arms_he_filename"] = arms_state.get("he_filename") or _prev_arms_attrs.get("arms_he_filename")
        attrs["arms_he_path"] = arms_state.get("he_path") or _prev_arms_attrs.get("arms_he_path")
        attrs["arms_he_shape_yx"] = (
            list(arms_state["he_shape_yx"]) if arms_state.get("he_shape_yx")
            else _prev_arms_attrs.get("arms_he_shape_yx")
        )
        attrs["arms_flip_v"] = bool(arms_state.get("flip_v", _prev_arms_attrs.get("arms_flip_v", False)))
        attrs["arms_flip_h"] = bool(arms_state.get("flip_h", _prev_arms_attrs.get("arms_flip_h", False)))
        if arms_state.get("affine_3x3") is not None:
            attrs["arms_affine_3x3"] = np.asarray(arms_state["affine_3x3"], dtype=np.float64).tolist()
        elif "arms_affine_3x3" in _prev_arms_attrs:
            attrs["arms_affine_3x3"] = _prev_arms_attrs["arms_affine_3x3"]
        attrs["arms_geojson_path"] = arms_state.get("geojson_path") or _prev_arms_attrs.get("arms_geojson_path")
        attrs["arms_csv_path"] = arms_state.get("csv_path") or _prev_arms_attrs.get("arms_csv_path")

        # ── Cluster labels (per-clustering nested dict) ────────────────────
        cluster_labels = state.get("cluster_labels")
        if cluster_labels and isinstance(cluster_labels, dict):
            # Nested dict: {clustering_name: {cluster_id: label}}
            serialized = {}
            for clust_name, label_dict in cluster_labels.items():
                if isinstance(label_dict, dict):
                    serialized[clust_name] = {str(k): v for k, v in label_dict.items()}
            attrs["cluster_labels"] = serialized if serialized else None
        else:
            attrs["cluster_labels"] = None

        # (Custom clusterings are now saved to adata.obs via adata_persistence)

        # ── Analysis DataFrames ───────────────────────────────────────────
        session_dir = Path(zarr_path) / "viewer_session"

        # rank_genes — now persisted to adata.uns + sdata.tables['adata_norm']
        # by save_rank_genes_to_adata(). Keep attrs for legacy reads.
        attrs["has_rank_genes"] = state.get("rank_genes_df") is not None
        attrs["rank_genes_groupby"] = state.get("rank_genes_groupby")

        # roi_deg_df
        roi_deg_df = state.get("roi_deg_df")
        if roi_deg_df is not None and not roi_deg_df.empty:
            roi_deg_df.to_parquet(session_dir / "roi_deg.parquet")
            attrs["has_roi_deg"] = True
        else:
            attrs["has_roi_deg"] = False

        # arms_tile_deg_df
        arms_tile_deg_df = state.get("arms_tile_deg_df")
        if arms_tile_deg_df is not None and not arms_tile_deg_df.empty:
            arms_tile_deg_df.to_parquet(session_dir / "arms_tile_deg.parquet")
            attrs["has_arms_tile_deg"] = True
        else:
            attrs["has_arms_tile_deg"] = False

        # (ligrec, nhood, co-occurrence are now saved to adata.uns via adata_persistence)

        # Write all attrs at once
        session.attrs.update(attrs)

        # Summary
        parts = []
        if attrs["roi_count"] > 0:
            parts.append(f"{attrs['roi_count']} ROIs")
        if attrs.get("he_filename"):
            parts.append(f"H&E ({attrs['he_filename']})")
        if attrs.get("cluster_labels"):
            parts.append(f"{len(attrs['cluster_labels'])} cluster labels")
        if attrs["has_rank_genes"]:
            parts.append("rank genes")
        if attrs["has_roi_deg"]:
            parts.append("ROI DEG")
        if attrs.get("has_arms_tile_deg"):
            parts.append("ARMS Tile DEG")
        if attrs.get("arms_he_filename"):
            parts.append(f"ARMS ({attrs['arms_he_filename']})")
        summary = ", ".join(parts) if parts else "empty session"
        print(f"Session saved: {summary}")

    except Exception as e:
        print(f"Warning: could not save session: {e}")



def save_rank_genes_incremental(zarr_path: Path, df, adata_norm, groupby: str) -> None:
    """Deprecated stub — rank genes are now saved via save_rank_genes_to_adata().

    Callers in tab_gene_analysis have been updated to call that function directly.
    This stub is retained for safety in case of any remaining references.
    """
    pass


def load_session(zarr_path: Path) -> Optional[dict]:
    """
    Load viewer session state from zarr store.

    Parameters
    ----------
    zarr_path : Path
        Path to sdata_cached.zarr directory.

    Returns
    -------
    dict or None
        Session data dict, or None if no session found.
    """
    import pandas as pd

    try:
        store = zarr.open_group(str(zarr_path), mode="r", use_consolidated=False)
    except Exception:
        return None

    if "viewer_session" not in store:
        return None

    session = store["viewer_session"]
    attrs = dict(session.attrs)

    result = {
        "rois": [],
        "he_filename": attrs.get("he_filename"),
        "he_path": attrs.get("he_path"),
        "he_shape_yx": tuple(attrs["he_shape_yx"]) if attrs.get("he_shape_yx") else None,
        "affine_3x3": None,
        "coarse_affine": None,
        "xenium_landmarks": None,
        "he_landmarks": None,
        "flip_v": attrs.get("flip_v", False),
        "flip_h": attrs.get("flip_h", False),
        "cluster_labels": None,
        "rank_genes_df": None,
        "rank_genes_adata_norm": None,
        "rank_genes_groupby": None,
        "roi_deg_df": None,
        "arms_tile_deg_df": None,
        "ligrec_means": None,
        "ligrec_pvalues": None,
        "nhood_result": None,
        "co_result": None,
        # ARMS overlay
        "arms_he_filename": attrs.get("arms_he_filename"),
        "arms_he_path": attrs.get("arms_he_path"),
        "arms_he_shape_yx": tuple(attrs["arms_he_shape_yx"]) if attrs.get("arms_he_shape_yx") else None,
        "arms_affine_3x3": None,
        "arms_xenium_landmarks": None,
        "arms_he_landmarks": None,
        "arms_flip_v": attrs.get("arms_flip_v", False),
        "arms_flip_h": attrs.get("arms_flip_h", False),
        "arms_geojson_path": attrs.get("arms_geojson_path"),
        "arms_csv_path": attrs.get("arms_csv_path"),
    }

    # ── ROIs ──────────────────────────────────────────────────────────
    roi_count = attrs.get("roi_count", 0)
    if "rois" in session:
        rois_group = session["rois"]
        for i in range(roi_count):
            key = str(i)
            if key in rois_group:
                result["rois"].append(np.array(rois_group[key]))

    # ── H&E registration ─────────────────────────────────────────────
    # Affine matrices: prefer attrs (saved in real-time) over zarr arrays (snapshot)
    if "affine_3x3" in attrs and attrs["affine_3x3"] is not None:
        result["affine_3x3"] = np.array(attrs["affine_3x3"], dtype=np.float64)
    if "coarse_affine" in attrs and attrs["coarse_affine"] is not None:
        result["coarse_affine"] = np.array(attrs["coarse_affine"], dtype=np.float64)

    # Fall back to zarr arrays if attrs don't have them
    if "he" in session:
        he_group = session["he"]
        if result["affine_3x3"] is None and "affine_3x3" in he_group:
            result["affine_3x3"] = np.array(he_group["affine_3x3"])
        if result["coarse_affine"] is None and "coarse_affine" in he_group:
            result["coarse_affine"] = np.array(he_group["coarse_affine"])
        if "xenium_landmarks" in he_group:
            result["xenium_landmarks"] = np.array(he_group["xenium_landmarks"])
        if "he_landmarks" in he_group:
            result["he_landmarks"] = np.array(he_group["he_landmarks"])

    # ── Cluster labels (per-clustering nested dict) ────────────────────
    cl = attrs.get("cluster_labels")
    if cl:
        # Detect old flat format (str(int) -> label) vs new nested format
        first_val = next(iter(cl.values()), None) if cl else None
        if isinstance(first_val, dict):
            # New nested format: {clustering_name: {cluster_id: label}}
            restored = {}
            for clust_name, label_dict in cl.items():
                inner = {}
                for k, v in label_dict.items():
                    # Try converting keys back to int
                    try:
                        inner[int(k)] = v
                    except (ValueError, TypeError):
                        inner[k] = v
                restored[clust_name] = inner
            result["cluster_labels"] = restored
        else:
            # Old flat format: {str(int): label} — wrap under a generic key
            result["cluster_labels"] = {"_legacy": {int(k): v for k, v in cl.items()}}

    # (Custom clusterings are now loaded from adata.obs via adata_persistence)

    # ── Analysis DataFrames (rank genes, ROI DEG, ARMS tile DEG remain here) ──
    session_dir = Path(zarr_path) / "viewer_session"

    result["rank_genes_groupby"] = attrs.get("rank_genes_groupby")
    if attrs.get("has_rank_genes"):
        p = session_dir / "rank_genes.parquet"
        if p.exists():
            result["rank_genes_df"] = pd.read_parquet(p)
        p = session_dir / "rank_genes_adata_norm.h5ad"
        if p.exists():
            import scanpy as sc
            result["rank_genes_adata_norm"] = sc.read_h5ad(p)

    if attrs.get("has_roi_deg"):
        p = session_dir / "roi_deg.parquet"
        if p.exists():
            result["roi_deg_df"] = pd.read_parquet(p)

    if attrs.get("has_arms_tile_deg"):
        p = session_dir / "arms_tile_deg.parquet"
        if p.exists():
            result["arms_tile_deg_df"] = pd.read_parquet(p)

    # (ligrec, nhood, co-occurrence are now loaded from adata.uns via adata_persistence)

    # ── ARMS overlay ─────────────────────────────────────────────────
    # Prefer attrs (real-time saved) over zarr arrays (snapshot)
    if "arms_affine_3x3" in attrs and attrs["arms_affine_3x3"] is not None:
        result["arms_affine_3x3"] = np.array(attrs["arms_affine_3x3"], dtype=np.float64)

    if "arms" in session:
        arms_group = session["arms"]
        if result["arms_affine_3x3"] is None and "affine_3x3" in arms_group:
            result["arms_affine_3x3"] = np.array(arms_group["affine_3x3"])
        if "xenium_landmarks" in arms_group:
            result["arms_xenium_landmarks"] = np.array(arms_group["xenium_landmarks"])
        if "he_landmarks" in arms_group:
            result["arms_he_landmarks"] = np.array(arms_group["he_landmarks"])

    return result
