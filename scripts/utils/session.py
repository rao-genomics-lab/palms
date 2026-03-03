"""
Session persistence for the Xenium viewer.

Saves and loads viewer state (ROIs, H&E registration, analysis results,
cluster labels) to/from a zarr store so the session can be restored on
the next launch.

Storage layout inside sdata_cached.zarr/viewer_session/:
  rois/0, rois/1, ...         — zarr arrays (Nx2 float64) for each ROI polygon
  he/affine_3x3               — zarr array 3x3
  he/coarse_affine             — zarr array 3x3
  he/xenium_landmarks          — zarr array Nx2
  he/he_landmarks              — zarr array Nx2
  rank_genes.parquet           — DataFrame
  roi_deg.parquet              — DataFrame
  ligrec_means.parquet         — DataFrame
  ligrec_pvalues.parquet       — DataFrame
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

        session = store.create_group("viewer_session", overwrite=True)
        attrs = {}

        # ── ROIs ──────────────────────────────────────────────────────────
        rois_group = session.create_group("rois")
        roi_data = snapshot.get("roi_data", [])
        for i, poly in enumerate(roi_data):
            _write_array(rois_group, str(i), poly)
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

        # ── Cluster labels ────────────────────────────────────────────────
        cluster_labels = state.get("cluster_labels")
        if cluster_labels:
            # Convert int keys to str for JSON serialization
            attrs["cluster_labels"] = {str(k): v for k, v in cluster_labels.items()}
        else:
            attrs["cluster_labels"] = None

        # ── Analysis DataFrames ───────────────────────────────────────────
        session_dir = Path(zarr_path) / "viewer_session"

        # rank_genes_df
        rank_genes_df = state.get("rank_genes_df")
        if rank_genes_df is not None and not rank_genes_df.empty:
            rank_genes_df.to_parquet(session_dir / "rank_genes.parquet")
            attrs["has_rank_genes"] = True
        else:
            attrs["has_rank_genes"] = False

        # roi_deg_df
        roi_deg_df = state.get("roi_deg_df")
        if roi_deg_df is not None and not roi_deg_df.empty:
            roi_deg_df.to_parquet(session_dir / "roi_deg.parquet")
            attrs["has_roi_deg"] = True
        else:
            attrs["has_roi_deg"] = False

        # ligrec results
        ligrec_result = state.get("ligrec_result")
        has_ligrec = False
        if ligrec_result is not None:
            means = ligrec_result.get("means")
            pvalues = ligrec_result.get("pvalues")
            if means is not None and not means.empty:
                means.to_parquet(session_dir / "ligrec_means.parquet")
                has_ligrec = True
            if pvalues is not None and not pvalues.empty:
                pvalues.to_parquet(session_dir / "ligrec_pvalues.parquet")
        attrs["has_ligrec"] = has_ligrec

        # nhood enrichment results
        nhood_result = state.get("nhood_result")
        has_nhood = False
        if nhood_result is not None and nhood_result.get("warning") is None:
            zscore = nhood_result.get("zscore")
            count = nhood_result.get("count")
            clusters = nhood_result.get("clusters", [])
            if zscore is not None and zscore.size > 0:
                nhood_group = session.create_group("nhood")
                _write_array(nhood_group, "zscore", zscore)
                _write_array(nhood_group, "count", count)
                nhood_group.attrs["clusters"] = clusters
                has_nhood = True
        attrs["has_nhood"] = has_nhood

        # co-occurrence results
        co_result = state.get("co_result")
        has_co = False
        if co_result is not None and co_result.get("warning") is None:
            occ = co_result.get("occ")
            interval_arr = co_result.get("interval")
            clusters = co_result.get("clusters", [])
            if occ is not None and occ.size > 0:
                co_group = session.create_group("co_occ")
                _write_array(co_group, "occ", occ)
                _write_array(co_group, "interval", interval_arr)
                co_group.attrs["clusters"] = clusters
                has_co = True
        attrs["has_co"] = has_co

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
        if attrs["has_ligrec"]:
            parts.append("L-R results")
        if attrs.get("has_nhood"):
            parts.append("nhood enrichment")
        if attrs.get("has_co"):
            parts.append("co-occurrence")
        summary = ", ".join(parts) if parts else "empty session"
        print(f"Session saved: {summary}")

    except Exception as e:
        print(f"Warning: could not save session: {e}")


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
        "roi_deg_df": None,
        "ligrec_means": None,
        "ligrec_pvalues": None,
        "nhood_result": None,
        "co_result": None,
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

    # ── Cluster labels ────────────────────────────────────────────────
    cl = attrs.get("cluster_labels")
    if cl:
        # Convert string keys back to int
        result["cluster_labels"] = {int(k): v for k, v in cl.items()}

    # ── Analysis DataFrames ───────────────────────────────────────────
    session_dir = Path(zarr_path) / "viewer_session"

    if attrs.get("has_rank_genes"):
        p = session_dir / "rank_genes.parquet"
        if p.exists():
            result["rank_genes_df"] = pd.read_parquet(p)

    if attrs.get("has_roi_deg"):
        p = session_dir / "roi_deg.parquet"
        if p.exists():
            result["roi_deg_df"] = pd.read_parquet(p)

    if attrs.get("has_ligrec"):
        p = session_dir / "ligrec_means.parquet"
        if p.exists():
            result["ligrec_means"] = pd.read_parquet(p)
        p = session_dir / "ligrec_pvalues.parquet"
        if p.exists():
            result["ligrec_pvalues"] = pd.read_parquet(p)

    # ── Nhood enrichment ─────────────────────────────────────────────
    if attrs.get("has_nhood") and "nhood" in session:
        nhood_group = session["nhood"]
        zscore = np.array(nhood_group["zscore"]) if "zscore" in nhood_group else np.array([])
        count = np.array(nhood_group["count"]) if "count" in nhood_group else np.array([])
        clusters = list(dict(nhood_group.attrs).get("clusters", []))
        result["nhood_result"] = {
            "zscore": zscore,
            "count": count,
            "clusters": clusters,
            "warning": None,
        }

    # ── Co-occurrence ────────────────────────────────────────────────
    if attrs.get("has_co") and "co_occ" in session:
        co_group = session["co_occ"]
        occ = np.array(co_group["occ"]) if "occ" in co_group else np.array([])
        interval_arr = np.array(co_group["interval"]) if "interval" in co_group else np.array([])
        clusters = list(dict(co_group.attrs).get("clusters", []))
        result["co_result"] = {
            "occ": occ,
            "interval": interval_arr,
            "clusters": clusters,
            "warning": None,
        }

    return result
