"""
Session persistence for the PALMS.

Saves and loads viewer state (ROIs, H&E registration, ARMS overlay,
rank genes, ROI DEG, cluster labels) to/from a zarr store so the
session can be restored on the next launch.

NOTE: Clusterings, nhood enrichment, co-occurrence, L-R results, and
UMAP coordinates are now persisted in adata.obs/obsm/uns via
utils/adata_persistence.py (Phase 1 SpatialData refactoring).

Storage layout inside sdata_cached.zarr/viewer_session/:
  he/affine_3x3               — zarr array 3x3
  he/coarse_affine             — zarr array 3x3
  arms/affine_3x3              — zarr array 3x3
  roi_deg.parquet              — DataFrame
  (group attrs contain JSON metadata)

Items now stored directly in sdata (via adata_persistence.py):
  sdata.shapes['rois']                — ROI polygons (Phase 2)
  sdata.shapes['he_xenium_landmarks'] — H&E Xenium landmark points (Phase 3)
  sdata.shapes['he_he_landmarks']     — H&E image landmark points (Phase 3)
  sdata.shapes['arms_xenium_landmarks'] — ARMS Xenium landmark points (Phase 3)
  sdata.shapes['arms_he_landmarks']   — ARMS image landmark points (Phase 3)
  sdata.shapes['arms_tiles']          — ARMS tile polygons + metadata (Phase 3)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import zarr

from palms.utils.reporting import get_logger, report_write_failure
from palms.utils.zarr_safe import safe_group_update

log = get_logger(__name__)

# Attrs written by other code paths and never recomputed here. Carrying every
# unknown key forward by default (rather than an allow-list) is what keeps them
# alive: the previous version rebuilt attrs from scratch, so each of these was
# silently wiped on every clean exit and its migration re-ran at the next
# launch — including two that themselves rewrote the whole cell table.
_PRESERVED_ON_SAVE = (
    "migrated_to_adata",
    "migrated_landmarks_to_sdata",
    "migrated_rank_genes_to_adata",
    "migrated_deg_to_sdata",
)

# Keys deliberately dropped rather than carried forward. Empty today; kept as
# the explicit place to retire a key, so the default stays "preserve".
_TRANSIENT_ATTR_KEYS: frozenset[str] = frozenset()


def _write_array(group, name, data):
    """Write a numpy array to a zarr group (compatible with zarr v2 and v3)."""
    arr = np.asarray(data, dtype=np.float64)
    ds = group.create_array(name, shape=arr.shape, dtype=arr.dtype)
    ds[:] = arr


def _read_prev_attrs(zarr_path: Path) -> dict:
    """Existing viewer_session attrs, or {} if absent/unreadable."""
    try:
        store = zarr.open_group(str(zarr_path), mode="r", use_consolidated=False)
        if "viewer_session" not in store:
            return {}
        return dict(store["viewer_session"].attrs)
    except Exception:
        return {}


def _json_safe(attrs: dict) -> tuple[dict, list[str]]:
    """Split *attrs* into what zarr can store and what it cannot.

    A single non-serialisable value used to abort the entire save — after the
    group had already been destroyed. Dropping the offender and reporting it
    keeps the rest of the session.
    """
    safe: dict = {}
    dropped: list[str] = []
    for key, value in attrs.items():
        try:
            json.dumps(value)
        except (TypeError, ValueError):
            dropped.append(key)
        else:
            safe[key] = value
    return safe, dropped


def _build_session_attrs(state: dict, he_state: dict, snapshot: dict,
                         prev_attrs: dict) -> dict:
    """Assemble the session attrs. Pure — no I/O, so it can be tested directly.

    Starts from *prev_attrs* so anything this function does not compute (the
    migration markers, and any key a future version adds) survives. Computed
    values always win, including explicit ``None`` — that is how clearing the
    H&E image actually clears it.
    """
    attrs = {k: v for k, v in prev_attrs.items() if k not in _TRANSIENT_ATTR_KEYS}

    # ── ROIs (stored in sdata.shapes['rois']; count kept for legacy reads) ──
    attrs["roi_count"] = len(snapshot.get("roi_data", []))

    # ── H&E registration ─────────────────────────────────────────────────
    attrs["he_filename"] = he_state.get("he_filename")
    attrs["he_path"] = he_state.get("he_path")
    attrs["he_shape_yx"] = (
        list(he_state["he_shape_yx"]) if he_state.get("he_shape_yx") else None
    )
    attrs["flip_v"] = bool(he_state.get("flip_v", False))
    attrs["flip_h"] = bool(he_state.get("flip_h", False))
    # The H&E's own declared pixel size, kept because a cache-restored H&E has
    # no TiffFile to read it from and it is the best scale prior coarse
    # alignment can get.
    attrs["he_pixel_size_um"] = (
        float(he_state["he_pixel_size_um"]) if he_state.get("he_pixel_size_um")
        else prev_attrs.get("he_pixel_size_um")
    )

    # ── ARMS overlay ─────────────────────────────────────────────────────
    # These are also written in real time by tab_arms, so an absent value in
    # arms_state means "unchanged", not "cleared" — hence the explicit fallback.
    arms_state = snapshot.get("arms_state", {})
    attrs["arms_he_filename"] = arms_state.get("he_filename") or prev_attrs.get("arms_he_filename")
    attrs["arms_he_path"] = arms_state.get("he_path") or prev_attrs.get("arms_he_path")
    attrs["arms_he_shape_yx"] = (
        list(arms_state["he_shape_yx"]) if arms_state.get("he_shape_yx")
        else prev_attrs.get("arms_he_shape_yx")
    )
    attrs["arms_flip_v"] = bool(arms_state.get("flip_v", prev_attrs.get("arms_flip_v", False)))
    attrs["arms_flip_h"] = bool(arms_state.get("flip_h", prev_attrs.get("arms_flip_h", False)))
    if arms_state.get("affine_3x3") is not None:
        attrs["arms_affine_3x3"] = np.asarray(
            arms_state["affine_3x3"], dtype=np.float64).tolist()
    attrs["arms_geojson_path"] = arms_state.get("geojson_path") or prev_attrs.get("arms_geojson_path")
    attrs["arms_csv_path"] = arms_state.get("csv_path") or prev_attrs.get("arms_csv_path")

    # ── Cluster labels (per-clustering nested dict) ──────────────────────
    cluster_labels = state.get("cluster_labels")
    serialized = {}
    if cluster_labels and isinstance(cluster_labels, dict):
        for clust_name, label_dict in cluster_labels.items():
            if isinstance(label_dict, dict):
                serialized[clust_name] = {str(k): v for k, v in label_dict.items()}
    attrs["cluster_labels"] = serialized or None

    # ── Analysis DataFrames (all now persisted into sdata/adata) ─────────
    attrs["has_rank_genes"] = state.get("rank_genes_df") is not None
    attrs["rank_genes_groupby"] = state.get("rank_genes_groupby")
    attrs["has_roi_deg"] = False
    attrs["has_arms_tile_deg"] = False

    attrs["marker_genes_json"] = state.get("marker_genes_json")
    attrs["umap_genes"] = list(state.get("umap_genes") or [])
    attrs["segmentation_source"] = state.get("segmentation_source", "xenium")

    # ── Reproducible-code provenance graph ───────────────────────────────
    # Never shrinks *here*. A smaller in-memory graph means it is not the
    # session's — a restore that did not happen, a launch that seeded only the
    # preamble. Saving it anyway is how a 13-node analysis became a 1-node stub:
    # the viewer came up empty, and its exit wrote that emptiness over the only
    # remaining copy.
    #
    # Deliberate removal exists — the Notebook tab's "Drop Stale Nodes" — but it
    # never relies on this path. It writes the pruned graph straight to the attr
    # itself (tab_notebook._persist_pruned_graph), so by the time an exit reaches
    # here `previous` is already the pruned copy and the sizes agree. Which is
    # what lets this stay a flat "never shrinks" rather than needing to tell a
    # prune apart from a failed restore, a distinction it has no way to make.
    prov = state.get("prov_graph")
    try:
        items = prov.to_list() if prov is not None and len(prov) else None
    except Exception:
        items = None
    previous = prev_attrs.get("prov_graph") or []
    if items is not None and len(items) < len(previous):
        log.warning(
            "keeping the stored provenance graph (%d nodes); the session holds "
            "only %d, which would lose recorded steps", len(previous), len(items),
        )
        attrs["prov_graph"] = previous
    elif items is None and previous:
        log.warning("keeping the stored provenance graph (%d nodes); the session "
                    "holds none", len(previous))
        attrs["prov_graph"] = previous
    else:
        attrs["prov_graph"] = items

    # ── External images / patch overlays UI residuals ────────────────────
    # An *empty* list means "none are loaded right now" — which is equally true
    # before restore has run and immediately after a cache recovery, so it must
    # not overwrite what is stored. Falling back on empty rather than only on
    # None is safe because restore is driven by the sdata elements, with these
    # attrs used only for contrast/opacity/affine: an entry left behind for a
    # removed image is simply never looked up.
    ext_ui = snapshot.get("external_images_ui") or prev_attrs.get("external_images_ui")
    attrs["external_images_ui"] = ext_ui or []

    patch_ui = snapshot.get("patch_overlays_ui") or prev_attrs.get("patch_overlays_ui")
    attrs["patch_overlays_ui"] = patch_ui or []

    return attrs


def _session_summary(attrs: dict) -> str:
    parts = []
    if attrs.get("roi_count"):
        parts.append(f"{attrs['roi_count']} ROIs")
    if attrs.get("he_filename"):
        parts.append(f"H&E ({attrs['he_filename']})")
    if attrs.get("cluster_labels"):
        parts.append(f"{len(attrs['cluster_labels'])} cluster labels")
    if attrs.get("has_rank_genes"):
        parts.append("rank genes")
    if attrs.get("arms_he_filename"):
        parts.append(f"ARMS ({attrs['arms_he_filename']})")
    return ", ".join(parts) if parts else "empty session"


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
        # Everything that can fail is done *before* the live group is touched:
        # the previous version destroyed viewer_session with
        # create_group(overwrite=True) and only wrote the replacement ~110 lines
        # later, so one non-serialisable value left an empty group and a printed
        # warning the user — who had already closed the window — never saw.
        prev_attrs = _read_prev_attrs(zarr_path)
        attrs = _build_session_attrs(state, he_state, snapshot, prev_attrs)
        attrs, dropped = _json_safe(attrs)
        if dropped:
            log.warning("session keys could not be serialized and were dropped: %s",
                        ", ".join(sorted(dropped)))

        arms_state = snapshot.get("arms_state", {})
        with safe_group_update(Path(zarr_path), "viewer_session") as (session, stage):
            # Parquet sidecars from previous sessions. Removed from *staging*,
            # so a failure below leaves the originals in place.
            for pattern in ("*.parquet", "clusterings/*.parquet"):
                for pq in stage.glob(pattern):
                    pq.unlink()

            # The affine subgroups are fully rewritten each save, so clear them
            # rather than merging into the seeded copy — otherwise a cleared
            # registration would leave its old affine behind.
            for sub in ("he", "arms"):
                if sub in session:
                    del session[sub]

            he_group = session.create_group("he")
            if he_state.get("affine_3x3") is not None:
                _write_array(he_group, "affine_3x3", he_state["affine_3x3"])
            if he_state.get("coarse_affine") is not None:
                _write_array(he_group, "coarse_affine", he_state["coarse_affine"])

            arms_group = session.create_group("arms")
            if arms_state.get("affine_3x3") is not None:
                _write_array(arms_group, "affine_3x3", arms_state["affine_3x3"])

            # Landmarks live in sdata.shapes now (save_landmarks_to_sdata).
            session.attrs.update(attrs)

        print(f"Session saved: {_session_summary(attrs)}")

    except Exception as e:
        report_write_failure(e, "session state")



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
        "he_pixel_size_um": attrs.get("he_pixel_size_um"),
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
        "marker_genes_json": attrs.get("marker_genes_json"),
        "umap_genes": list(attrs.get("umap_genes") or []),
        "segmentation_source": attrs.get("segmentation_source", "xenium"),
        "external_images_ui": attrs.get("external_images_ui") or [],
        "patch_overlays_ui": attrs.get("patch_overlays_ui") or [],
        "prov_graph": attrs.get("prov_graph"),
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

    # ── Analysis DataFrames ───────────────────────────────────────────────
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

    # roi_deg and arms_tile_deg are now loaded from sdata via
    # load_roi_deg_from_sdata / load_arms_tile_deg_from_sdata (called in 02_palms.py).

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
