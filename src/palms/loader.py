"""
SpatialData loading for Xenium 3.x output.

Also the ``palms-build-cache`` console script (see ``main`` at the bottom of
this file), which runs the same load headlessly so the slow first read of a
dataset can happen without holding a napari window open.

Returns a SpatialData object with:
  - images:  morphology_focus (multiscale, 4-channel, CYX)
  - labels:  cell_labels, nucleus_labels
  - points:  transcripts (dask-backed, lazy)
  - tables:  table (AnnData 318K cells x 480+ genes)
"""

import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from palms.utils import mem_probe, sdata_write

# ─── Channel names for the 4 morphology_focus planes ───────────────────────
CHANNEL_NAMES = [
    "DAPI",
    "ATP1A1-CD45-E-Cadherin",
    "18S",
    "AlphaSMA-Vimentin",
]

# Pyramid levels built for morphology_focus: scale0 plus one per factor.
PYRAMID_SCALE_FACTORS = (2, 2, 2, 2, 2)

# ─── User-generated element keys ────────────────────────────────────────────
_USER_SHAPE_KEYS = [
    "rois", "he_xenium_landmarks", "he_he_landmarks",
    "arms_xenium_landmarks", "arms_he_landmarks", "arms_tiles", "annotations",
]
_USER_IMAGE_KEYS = ["he_image", "arms_he_image"]

# External images and patch overlays are named per file, so a fixed list cannot
# see them: a dataset with a registered PhenoCycler image and its landmarks
# reported "no user data" and could be rebuilt over without a prompt.
_USER_SHAPE_SUFFIXES = ("_xenium_lm", "_image_lm")
_USER_ELEMENT_PREFIXES = ("ext_", "patch_")


def _is_user_element(name: str, keys: list) -> bool:
    return (name in keys
            or name.startswith(_USER_ELEMENT_PREFIXES)
            or name.endswith(_USER_SHAPE_SUFFIXES))


_USER_UNS_KEYS = ["nhood_enrichment", "co_occurrence", "ligrec", "rank_genes_groups"]

# A prefix, not a name: the rank-genes step writes `rank_genes_<clustering>` so
# ranking a second clustering does not overwrite the first, and `rank_genes_groupby`
# names the most recent. A fixed list would have stopped counting a ranking as
# user data the moment the keying landed — and "no user data" is what lets a
# cache be rebuilt with no dialog.
_USER_UNS_PREFIXES = ("rank_genes",)


def _is_user_uns(name: str) -> bool:
    return name in _USER_UNS_KEYS or name.startswith(_USER_UNS_PREFIXES)

# Globs, not literal names: the previous fixed list omitted the CNV caches, so a
# cache whose only user data was a multi-hour CopyKAT run reported "no user
# data" and was rebuilt with no dialog at all.
_SIDECAR_PATTERNS = [
    "roi_deg_cache.parquet",
    "arms_tile_deg_cache.parquet",
    "adata_norm_cache.h5ad",
    "adata_cnv_cache_*.h5ad",
    "cnv_*_result.json",
]

# obs columns that exist only because the user ran something.
_USER_OBS_PREFIXES = ("clustering_", "cluster_labels_", "cnv_score", "copykat_leiden_res")


class CacheLoadAborted(RuntimeError):
    """The user chose to quit rather than rebuild a cache."""


class NoRawSourceError(RuntimeError):
    """A rebuild was requested for a dataset that has nothing to rebuild from.

    Raised instead of letting ``spatialdata_io.xenium`` fail on a missing
    ``cell_feature_matrix.h5`` — the h5py error names a file the user never had
    and says nothing about what to do next.
    """


# The files ``spatialdata_io.xenium`` actually reads. A Crop Dataset export has
# none of them: it ships ``experiment.xenium`` (which the app requires) plus the
# zarr cache and derived transcripts, and that is the whole dataset.
_RAW_PALMS_MARKERS = (
    "cells.zarr.zip",
    "cell_feature_matrix.h5",
    "cell_feature_matrix",
    "morphology_focus",
)


def has_raw_xenium_source(path: Path) -> bool:
    """Can this directory be rebuilt from raw 10x output?

    Deliberately conservative — ``True`` unless *none* of the markers is present.
    A dataset missing only some of them is broken raw output, and should keep
    raising whatever error describes that, rather than being reclassified as a
    Crop Dataset export.

    This is the single definition. It used to be an inline
    ``not (path / "cells.zarr.zip").exists()`` guarding only the ``--no-cache``
    override, which is why every rebuild path walked straight into a read that
    could not work.
    """
    path = Path(path)
    return any((path / marker).exists() for marker in _RAW_PALMS_MARKERS)


def cache_only_declared(path: Path) -> Optional[bool]:
    """What the cache's own manifest says about being cache-only, if anything.

    ``None`` means the manifest is absent, unreadable, or silent — every cache
    built before the Crop Dataset tool started stamping it, and every raw
    dataset. Never raises: a manifest that cannot be parsed is a manifest that
    says nothing, and the caller falls back to ``has_raw_xenium_source``.
    """
    import json

    from palms.utils.zarr_safe import MANIFEST_FILE

    try:
        manifest = json.loads((Path(path) / "sdata_cached.zarr" / MANIFEST_FILE).read_text())
    except (OSError, ValueError):
        return None
    value = manifest.get("cache_only")
    return bool(value) if isinstance(value, bool) else None


def is_cache_only(path: Path) -> bool:
    """Is this dataset's zarr store the source of record rather than a cache?

    The **declaration wins over the inference.** ``has_raw_xenium_source`` asks
    a different question — "are raw files present?" — and only answers this one
    by accident, because a Crop Dataset export happens to write none of them.
    An export that writes raw-shaped files (the raw-format export) would flip
    that inference, reopen the rebuild paths on a dataset whose cache is the
    only copy of the data, and silently revert the recorded preamble to
    ``xenium(data_path)``, which reads the raw half and drops every derived
    layer the crop carried. The export already stamps ``cache_only`` for
    exactly this reason; until now nothing read it.

    Only a ``True`` stamp overrides. A ``False`` one is read but not trusted to
    *unset* cache-only, because being wrong in that direction sends the loader
    down a rebuild path on a dataset whose only copy is the cache — the stamp
    may add certainty, never remove protection.
    """
    if cache_only_declared(path):
        return True
    return not has_raw_xenium_source(path)


def cache_only_reason(path: Path) -> str:
    """Why ``is_cache_only`` said yes — for messages that must not misdescribe it.

    Enumerating absent files at a dataset that *declares* itself cache-only
    while holding raw-shaped files would be a message contradicted by the
    directory it describes.
    """
    if cache_only_declared(path):
        return "its cache manifest declares it one (a Crop Dataset export)"
    return ("it holds no raw Xenium output — no cells.zarr.zip, no "
            "cell_feature_matrix, no morphology_focus")


def _no_raw_source_message(path: Path, reason: str) -> str:
    return (
        f"{reason}\n\n"
        f"{path} holds no raw Xenium output — no cells.zarr.zip, no "
        "cell_feature_matrix, no morphology_focus. It is most likely a Crop "
        "Dataset export, whose zarr cache is the only copy of the data, so "
        "there is nothing to rebuild from.\n\n"
        "Nothing has been moved or deleted. If the cache itself is damaged, "
        "run `palms-build-cache --check` on this directory to see what is "
        "wrong, and look for a sdata_cached_backup_*.zarr or "
        "sdata_cached_prev_*.zarr sibling to recover from — an earlier forced "
        "rebuild may have left one."
    )


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
    shapes_dir = cache_path / "shapes"
    if shapes_dir.is_dir():
        for entry in sorted(shapes_dir.iterdir()):
            if entry.is_dir() and _is_user_element(entry.name, _USER_SHAPE_KEYS):
                found["shapes"].append(entry.name)
    images_dir = cache_path / "images"
    if images_dir.is_dir():
        for entry in sorted(images_dir.iterdir()):
            if entry.is_dir() and _is_user_element(entry.name, _USER_IMAGE_KEYS):
                found["images"].append(entry.name)
    obs_dir = cache_path / "tables" / "table" / "obs"
    if obs_dir.exists():
        for item in obs_dir.iterdir():
            if item.name.startswith("clustering_"):
                found["clusterings"].append(item.name)
    uns_dir = cache_path / "tables" / "table" / "uns"
    if uns_dir.exists():
        for item in sorted(uns_dir.iterdir()):
            if not item.name.startswith(".") and _is_user_uns(item.name):
                found["uns_keys"].append(item.name)
    obsm_dir = cache_path / "tables" / "table" / "obsm"
    if obsm_dir.exists() and (obsm_dir / "X_umap").exists():
        found["has_obsm_umap"] = True
    # Sidecars now live beside the store, where a cache rebuild cannot touch
    # them; the in-store location is still scanned so pre-existing datasets
    # still count them as data at stake.
    sidecar_home = cache_path.parent / "viewer_cache"
    for pattern in _SIDECAR_PATTERNS:
        found["sidecars"].extend(sorted(p.name for p in cache_path.glob(pattern)))
        found["sidecars"].extend(sorted(p.name for p in sidecar_home.glob(pattern)))
    found["sidecars"] = sorted(set(found["sidecars"]))
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
            "rank_genes_groupby": "Rank genes results",
        }
        seen_uns: set[str] = set()
        for key in user_data["uns_keys"]:
            # Every rank_genes_<clustering> slot is one line — a dataset with
            # four rankings should not print four indistinguishable bullets.
            label = labels.get(key) or (
                "Rank genes results" if key.startswith("rank_genes") else key)
            if label in seen_uns:
                continue
            seen_uns.add(label)
            lines.append(f"  • {label}")
    _sidecar_labels = {
        "roi_deg_cache.parquet": "ROI DEG results",
        "arms_tile_deg_cache.parquet": "ARMS tile DEG results",
        "adata_norm_cache.h5ad": "Normalized expression cache",
    }
    seen: set[str] = set()
    for fname in user_data["sidecars"]:
        if fname.startswith("adata_cnv_cache_") or fname.startswith("cnv_"):
            backend = "CopyKAT" if "copykat" in fname else "inferCNV"
            label = f"{backend} CNV results (hours of compute)"
        else:
            label = _sidecar_labels.get(fname, fname)
        # Several files map to one label (the h5ad and its result JSON).
        if label in seen:
            continue
        seen.add(label)
        lines.append(f"  • {label}")
    return "\n".join(lines)


def _stale_preference(user_data: dict, certain: bool,
                      on_stale: Optional[str] = None) -> Optional[str]:
    """Decide what to do about a stale cache. Pure apart from its printing.

    Returns 'restore', 'rebuild', 'keep', or ``None`` for "say nothing more and
    rebuild" — the pre-manifest case, where there is nothing to lose.

    ``on_stale`` is the caller answering in advance (``palms-build-cache
    --on-stale``), which is the only way to authorise a rebuild from a terminal:
    with no dialog every branch below either keeps the cache or asks. It is
    checked first so it covers the branches that never reach a prompt, and it is
    deliberately not a default — an explicit opt-in is not the same thing as a
    silent one.
    """
    if on_stale is not None:
        print(f"Zarr cache is stale; proceeding with '{on_stale}' as instructed.")
        return on_stale

    if not certain and not _has_any_user_data(user_data):
        # Pre-manifest cache, mtime-only signal, nothing to lose:
        # rebuilding is cheap insurance and stamps a manifest.
        print("Zarr cache may be stale (no manifest; experiment.xenium is "
              "newer). Rebuilding...")
        return None

    if not certain:
        print("Zarr cache may be stale — the check is a file timestamp, "
              "which a copy or re-download also changes.")
        return _ask_rebuild_preference(user_data, certain=False)

    preference = (_ask_rebuild_preference(user_data)
                  if _has_any_user_data(user_data) else "rebuild")
    if preference == "rebuild":
        print("Zarr cache is stale (experiment.xenium has changed). Rebuilding...")
    return preference


def _ask_rebuild_preference(user_data: dict, certain: bool = True) -> str:
    """Ask what to do about a stale cache. Returns 'restore', 'rebuild' or 'keep'.

    Without a way to prompt this now returns 'keep' rather than 'rebuild'. The
    old default silently discarded a cache — potentially 30 GB and hours of CNV
    compute — whenever the dialog was unavailable, including on any background
    thread. Keeping is recoverable; rebuilding is not.
    """
    message_box = _qt_message_box()
    if message_box is None:
        print("Warning: zarr cache looks stale and contains user data, but there "
              "is no way to ask. Keeping the existing cache — delete it by hand "
              "to force a rebuild.")
        return "keep"

    summary = _format_user_data_message(user_data)
    reason = (
        "experiment.xenium has changed since the cache was built."
        if certain else
        "experiment.xenium is newer than the cache. This cache predates content\n"
        "checking, so the comparison is a file timestamp — which copying or\n"
        "re-downloading the dataset also changes, even when nothing differs."
    )
    msg = message_box()
    msg.setWindowTitle("Zarr Cache May Need Rebuilding")
    msg.setText(
        f"{reason}\n\n"
        "Your cache contains user-generated data:\n"
        f"{summary}\n\n"
        "Rebuilding re-reads the raw Xenium files and can take a long time.\n"
        "What would you like to do?"
    )
    msg.setIcon(message_box.Warning)
    restore_btn = msg.addButton("Rebuild and restore my data", message_box.AcceptRole)
    msg.addButton("Rebuild without restoring", message_box.DestructiveRole)
    keep_btn = msg.addButton("Keep existing cache", message_box.RejectRole)
    msg.setDefaultButton(keep_btn if not certain else restore_btn)
    msg.exec()

    clicked = msg.clickedButton()
    if clicked is restore_btn:
        return "restore"
    elif clicked is keep_btn:
        return "keep"
    else:
        return "rebuild"


def _qt_message_box():
    """Return a usable QMessageBox class, or None when we cannot prompt."""
    import threading
    if threading.current_thread() is not threading.main_thread():
        return None
    try:
        from qtpy.QtWidgets import QApplication, QMessageBox
    except Exception:
        return None
    return QMessageBox if QApplication.instance() is not None else None


def _ask_corrupt_cache(error: Exception, report, user_data: dict,
                       preset: Optional[str] = None) -> str:
    """Ask what to do about a cache that will not open even after repair.

    Returns 'restore', 'rebuild' or 'quit'. Never returns a destructive default:
    with no way to prompt we raise instead, because the previous behaviour —
    silently rebuilding 30 GB — is the thing being fixed.

    ``preset`` answers in advance for a headless caller, but only for the two
    answers that make sense here. 'keep' does not: this cache cannot be opened,
    so keeping it is not a way to carry on — it falls through to the raise
    below, which is what tells the user their cache is broken.
    """
    if preset in ("restore", "rebuild"):
        return preset

    summary = _format_user_data_message(user_data)
    detail = (
        f"The zarr cache could not be opened:\n  {error}\n\n"
        f"{report.summary()}\n"
    )
    message_box = _qt_message_box()
    if message_box is None:
        raise CacheLoadAborted(
            f"{detail}\nThe cache contains user-generated data:\n{summary}\n\n"
            "Refusing to rebuild it without confirmation. Re-run with a GUI to "
            "choose, or move the cache aside yourself to force a rebuild."
            if _has_any_user_data(user_data) else
            f"{detail}\nRe-run with a GUI to choose how to proceed."
        )

    msg = message_box()
    msg.setWindowTitle("Zarr Cache Could Not Be Opened")
    msg.setText(
        "The zarr cache could not be opened, and automatic repair did not fix it.\n\n"
        + (f"Your cache contains user-generated data:\n{summary}\n\n"
           if _has_any_user_data(user_data) else "")
        + "Rebuilding re-reads the raw Xenium files and can take a long time.\n"
          "The existing cache will be kept aside either way — nothing is deleted."
    )
    msg.setDetailedText(detail)
    msg.setIcon(message_box.Warning)
    restore_btn = msg.addButton("Rebuild and restore my data", message_box.AcceptRole)
    msg.addButton("Rebuild without restoring", message_box.DestructiveRole)
    quit_btn = msg.addButton("Quit", message_box.RejectRole)
    msg.setDefaultButton(restore_btn)
    msg.exec()

    clicked = msg.clickedButton()
    if clicked is restore_btn:
        return "restore"
    if clicked is quit_btn:
        return "quit"
    return "rebuild"


def _reopen_written_cache(sdata, cache_path: Path):
    """Swap the freshly built SpatialData for the one we just wrote to disk.

    The object returned by `xenium()` reaches its rasters through lazy dask
    graphs rooted in the OME-TIFF, and its pyramid levels above `scale0` exist
    only as a chained `coarsen().mean()`. The level *napari shows first* is the
    smallest one (`_data_level = len(data) - 1`), so adding the layer walks that
    chain all the way down: on a full slide `scale5` is 20 MiB of pixels
    standing on ~13,400 tasks over all 780 `scale0` chunks, promoted to float64
    on the way. That is tens of GB of pressure for a thumbnail, three layers
    over, and it is why building and displaying in one session was killed while
    a restart against the same cache is fine.

    Re-reading costs a few seconds and is lazy; every level is then a direct
    read of an array that is already on disk. It also settles `sdata.path`
    without a separate assignment — the write pointed it at the staging
    directory, which no longer exists, and every later element write resolves
    against it (`zarr_safe._store_root`), so a stale value sent the first ROI or
    clustering into a resurrected `.sdata_cached__building.zarr`.

    Never raises: a cache that was just written should open, but if it does not,
    the caller is better off with the in-memory object than with an exception
    that loses a build. In that case the path is set the old way, and the layers
    are built from the chain — slowly, and only warned about by
    `app._warn_if_pyramid_is_not_stored`; nothing caps it.
    """
    try:
        reopened = _open_cache(cache_path)
    except Exception as e:
        sdata.path = cache_path
        _print_and_log(
            f"Warning: could not re-open the cache just written ({e}). "
            "Continuing with the in-memory dataset — building the layers will "
            "use much more memory than usual."
        )
        return sdata
    return reopened


def _open_cache(cache_path: Path):
    """Open the cache, repairing first and escalating if that is not enough.

    Returns the SpatialData, or raises the last read error. Repair is attempted
    before *and* after a plain read because the common corruption — an element
    present on disk but absent from the consolidated metadata — makes the read
    fail while costing nothing to fix.
    """
    import spatialdata

    from palms.utils import cache_repair

    report = cache_repair.verify(cache_path)
    if not report.ok:
        result = cache_repair.repair(cache_path, report, level=cache_repair.AUTO)
        for action in result.actions:
            print(f"  Cache: {action}")

    try:
        return spatialdata.read_zarr(str(cache_path))
    except Exception as first_error:
        # Still broken. If a previous version of a missing element is sitting
        # in the trash, putting it back is strictly better than discarding
        # 30 GB.
        report = cache_repair.verify(cache_path)
        if report.missing_on_disk and report.repairable:
            print("  Cache: attempting to restore missing elements from backups...")
            result = cache_repair.repair(cache_path, report, level=cache_repair.FULL)
            for action in result.actions:
                print(f"  Cache: {action}")
            return spatialdata.read_zarr(str(cache_path))
        raise first_error


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

        # Deny-list, not allow-list. The previous version restored only
        # clustering_*/cluster_labels_* obs and four named uns keys, so CNV
        # scores, CNV run metadata and the rank-genes groupby were silently
        # dropped even when the user explicitly chose "restore my data".
        # Anything a freshly-built table does not already have is, by
        # definition, something the user's session added.
        for col in old_adata.obs.columns:
            if col in new_adata.obs.columns and not col.startswith(_USER_OBS_PREFIXES):
                continue
            try:
                new_adata.obs[col] = old_adata.obs[col].reindex(new_adata.obs.index)
                if col.startswith("clustering_"):
                    restored.append(col)
            except Exception as e:
                print(f"  Warning: could not restore obs column '{col}': {e}")

        for key in old_adata.uns:
            if key in new_adata.uns:
                continue
            try:
                new_adata.uns[key] = old_adata.uns[key]
                restored.append(f"uns/{key}")
            except Exception as e:
                print(f"  Warning: could not restore uns['{key}']: {e}")

        for key in old_adata.obsm:
            if key in new_adata.obsm:
                continue
            try:
                if len(old_adata.obsm[key]) == len(new_adata.obs):
                    new_adata.obsm[key] = old_adata.obsm[key]
                    restored.append("UMAP coordinates" if key == "X_umap" else f"obsm/{key}")
            except Exception as e:
                print(f"  Warning: could not restore obsm['{key}']: {e}")
    return restored


# ─── Cache freshness ────────────────────────────────────────────────────────

def _source_fingerprint(experiment_path: Path) -> dict:
    """Content hash of experiment.xenium — a few KB of JSON, so this is free."""
    import hashlib
    data = experiment_path.read_bytes()
    return {
        "source": experiment_path.name,
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "source_size": len(data),
    }


def write_manifest(cache_path: Path, experiment_path: Path,
                   extra: Optional[dict] = None) -> None:
    """Record what this cache was built from, so staleness can be exact.

    ``extra`` is merged in for facts the fingerprint cannot carry — Crop Dataset
    exports stamp ``cache_only``/``derived_from`` so a reader knows the cache is
    the source of record rather than a derivative of files sitting beside it.
    """
    from datetime import datetime, timezone

    from palms.utils.zarr_safe import MANIFEST_FILE, atomic_json

    manifest = {"built_at": datetime.now(timezone.utc).isoformat()}
    if extra:
        manifest.update(extra)
    try:
        import spatialdata
        import zarr as _zarr
        manifest["spatialdata_version"] = spatialdata.__version__
        manifest["zarr_version"] = _zarr.__version__
    except Exception:
        pass
    try:
        manifest.update(_source_fingerprint(experiment_path))
    except OSError:
        return          # no source to fingerprint; leave the cache unstamped
    try:
        atomic_json(cache_path / MANIFEST_FILE, manifest)
    except OSError as e:
        print(f"  Warning: could not write cache manifest: {e}")


def _is_cache_stale(cache_path: Path, experiment_path: Path) -> tuple[bool, bool]:
    """Return ``(stale, certain)``.

    The old rule compared ``experiment.xenium``'s mtime against the *directory*
    mtime of the cache. Directory mtime only moves when a direct child is
    added or removed, while ``rsync``/``cp -p``/re-downloading the dataset bumps
    the source — so perfectly good caches were condemned. When a manifest is
    present the answer comes from a content hash and is exact; without one
    (every cache built before this change) the mtime comparison survives only as
    an uncertain hint, and the caller must ask rather than rebuild silently.
    """
    from palms.utils.cache_repair import read_manifest

    if not experiment_path.exists():
        return False, True

    manifest = read_manifest(cache_path)
    if manifest and manifest.get("source_sha256"):
        try:
            current = _source_fingerprint(experiment_path)
        except OSError:
            return False, True
        return current["source_sha256"] != manifest["source_sha256"], True

    return experiment_path.stat().st_mtime > cache_path.stat().st_mtime, False


def _copy_sidecars_and_session(backup_path: Path, cache_path: Path, user_data: dict) -> None:
    """Copy sidecar files and viewer_session zarr group from backup to new cache."""
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


def _print_and_log(line: str) -> None:
    """Progress detail that belongs both on screen and in the dataset log."""
    from palms.utils import reporting
    print(f"  {line}")
    reporting.get_logger().info(line)


def _retile_morphology_focus(sdata, path: Path) -> bool:
    """Re-read morphology_focus from the OME-TIFF's own 1024px tiles.

    ``spatialdata_io`` reads that image with ``dask_image.imread``, which gives
    one dask chunk per *full channel page* — 5.93 GB each on a typical slide, and
    every one of them has to be decoded and held before a single output tile can
    be written. The file is tiled, so reading it that way instead removes ~24 GB
    from the cache build. The pyramid is still built by spatialdata from this
    array, so **every written byte is unchanged**; only the chunking of the read
    differs. See ``utils/raster_io.py`` for the measurements, including why the
    TIFF's own pyramid levels are *not* reused.

    Returns True if the swap happened. Declining leaves spatialdata_io's element
    untouched — a performance regression at worst, never a correctness one — so
    this never raises.
    """
    from spatialdata.models import Image2DModel

    from palms.utils import raster_io, reporting

    if "morphology_focus" not in sdata.images:
        return False
    reference = sdata.images["morphology_focus"]
    try:
        level0 = raster_io.tiled_morphology_image(Path(path), reference)
        scale_factors = raster_io.pyramid_scale_factors(reference)
        if level0 is None or scale_factors is None:
            reporting.get_logger().info(
                "morphology_focus: no usable tiled source; reading it as before.")
            return False
        sdata.images["morphology_focus"] = Image2DModel.parse(
            level0,
            dims=("c", "y", "x"),
            transformations=raster_io.reference_transformations(reference),
            chunks=raster_io.DEFAULT_CHUNKS,
            # Taken from the element being replaced, not from a constant: this
            # has to reproduce whatever levels spatialdata_io actually built,
            # and it fills in its own default scale_factors when we pass none.
            scale_factors=scale_factors or None,
            c_coords=raster_io.reference_channels(reference),
            rgb=None,
        )
    except Exception as exc:
        reporting.get_logger().warning(
            "Could not re-read morphology_focus from its tiles (%s); "
            "reading it as before.", exc)
        return False

    print(f"  morphology_focus: reading the OME-TIFF's own "
          f"{level0.chunksize[-2:]} tiles instead of whole pages.")
    return True


ON_STALE_CHOICES = ("keep", "rebuild", "restore")


def load_sdata(
    path: Path,
    build_pyramid: bool = True,
    n_jobs: int = 8,
    use_cache: bool = True,
    on_stale: Optional[str] = None,
):
    """
    Load the Xenium 3.x output as a SpatialData object.

    Parameters
    ----------
    path : Path
        Root directory of the Xenium output.
    build_pyramid : bool
        If True, build a 5-level image pyramid for morphology_focus. Required for
        smooth napari pan/zoom performance. The source OME-TIFFs do carry their
        own pyramid, but it is computed with a different filter than spatialdata's
        (measured ~15% per-pixel difference — see ``utils/raster_io.py``), so only
        their full-resolution level is read and the rest is recomputed.
    n_jobs : int
        Number of threads for spatialdata_io.
    use_cache : bool
        If True, use/create a zarr cache for faster subsequent loads.
    on_stale : {'keep', 'rebuild', 'restore'}, optional
        Answer the stale/unopenable-cache question in advance instead of
        prompting. ``None`` (the default, and what the GUI passes) prompts if a
        dialog is available and otherwise keeps the cache. Used by
        ``palms-build-cache --on-stale``, which is the only way to authorise a
        rebuild from a terminal.

    Returns
    -------
    spatialdata.SpatialData
    """
    from datetime import datetime

    if on_stale is not None and on_stale not in ON_STALE_CHOICES:
        raise ValueError(
            f"on_stale must be one of {ON_STALE_CHOICES} or None, got {on_stale!r}"
        )

    cache_path = path / "sdata_cached.zarr"
    experiment_path = path / "experiment.xenium"
    backup_path = None   # set when we move the old cache for a user-restore
    user_data = None     # populated when we detect user data in the old cache

    # A directory with no raw Xenium files (e.g. one exported by the Crop
    # Dataset tool — just experiment.xenium + sdata_cached.zarr + transcripts)
    # can only ever be loaded from its zarr cache; there's nothing to rebuild
    # from. Use the cache even if the caller asked for use_cache=False (e.g.
    # launched with --no-cache), rather than failing on a missing cells.zarr.zip.
    cache_only = is_cache_only(path)
    if not use_cache and cache_path.exists() and cache_only:
        print(
            "No raw Xenium files found in this directory (likely a Crop Dataset "
            "export) — loading from the zarr cache regardless of --no-cache."
        )
        use_cache = True

    if cache_only and not cache_path.exists():
        raise NoRawSourceError(_no_raw_source_message(
            path, "There is no zarr cache to load and no raw Xenium output to build one from."
        ))

    # Try loading from zarr cache if it exists and is fresh
    if use_cache and cache_path.exists():
        # A cache-only dataset is never "stale": the cache *is* the source of
        # record, and experiment.xenium is a copy carried along beside it. Asking
        # the staleness question at all is what armed the trap — the mtime
        # fallback (no manifest, a copy or a touch reorders the timestamps) sent
        # crop exports into a rebuild branch with nothing to rebuild from.
        if cache_only:
            stale, certain = False, True
            print("No raw Xenium output beside this cache (likely a Crop Dataset "
                  "export) — the cache is the source of record, so it is never "
                  "rebuilt.")
        else:
            stale, certain = _is_cache_stale(cache_path, experiment_path)
        user_data = _detect_user_data(cache_path)
        preference = None

        if stale:
            preference = _stale_preference(user_data, certain, on_stale)

        if not stale or preference == "keep":
            print(f"Loading SpatialData from zarr cache: {cache_path}")
            try:
                sdata = _open_cache(cache_path)
                print("SpatialData loaded from cache.")
                if preference == "keep":
                    # They chose to keep it, so stop asking: stamp the manifest
                    # against the source as it stands now.
                    write_manifest(cache_path, experiment_path)
                elif not certain:
                    write_manifest(cache_path, experiment_path)
                elif cache_only:
                    # Stamp exports made before this change, so the next reader
                    # gets the fact from the manifest rather than re-deriving it.
                    write_manifest(cache_path, experiment_path,
                                   extra={"cache_only": True})
            except Exception as e:
                from palms.utils import cache_repair
                report = cache_repair.verify(cache_path)
                if cache_only:
                    # Every answer _ask_corrupt_cache can give other than
                    # 'quit' moves the cache aside and then rebuilds. Here that
                    # renames away the only copy of the data and fails anyway.
                    raise NoRawSourceError(_no_raw_source_message(
                        path,
                        f"The zarr cache could not be opened:\n  {e}\n\n{report.summary()}",
                    )) from e
                choice = _ask_corrupt_cache(e, report, user_data, preset=on_stale)
                if choice == "quit":
                    raise CacheLoadAborted(
                        "Cancelled at the user's request; the cache was left untouched."
                    ) from e
                preference = choice
            else:
                # The summary is diagnostics, and it stays *outside* the block
                # above on purpose: a repr that cannot render is not a cache
                # that cannot be opened. It used to be inside, and when
                # spatialdata's ``__repr__`` raised — it walks each points
                # element's dask graph for its backing files, which a
                # non-default parquet reader can make unparseable — a healthy
                # store was routed to _ask_corrupt_cache and the user was
                # offered a rebuild, with verify printing "✓ Cache is healthy"
                # two lines later.
                try:
                    print(sdata)
                except Exception as exc:
                    from palms.utils import reporting
                    reporting.get_logger().warning(
                        "The SpatialData summary could not be rendered (%s); "
                        "the cache itself opened fine.", exc)
                return sdata

        if preference == "restore":
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_path = cache_path.with_name(f"sdata_cached_backup_{timestamp}.zarr")
            shutil.move(str(cache_path), str(backup_path))
            print(f"Old cache moved aside to:\n  {backup_path}\n"
                  "Rebuilding cache and restoring user data...")
        elif preference == "rebuild":
            # Never overwrite in place and never delete: the old rebuild path
            # wrote over the live cache and rmtree'd it if the write failed,
            # which is unrecoverable. Moving aside costs a rename.
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            kept_aside = cache_path.with_name(f"sdata_cached_prev_{timestamp}.zarr")
            shutil.move(str(cache_path), str(kept_aside))
            print(f"Previous cache kept at:\n  {kept_aside}\n"
                  "Rebuilding without restoring user data...")

    # The one line that makes every route safe, including any added later: below
    # this point the only way to produce an sdata is to read the raw output.
    # Reaching here without it used to surface as h5py failing to open a
    # cell_feature_matrix.h5 the user never had.
    if cache_only:
        raise NoRawSourceError(_no_raw_source_message(
            path, "A cache rebuild was requested, but it cannot be done."
        ))

    from spatialdata_io import xenium

    image_models_kwargs = {}
    if build_pyramid:
        image_models_kwargs = {"scale_factors": list(PYRAMID_SCALE_FACTORS)}

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

    # Unconditional, not gated on build_pyramid: the whole-page read happens
    # either way (spatialdata_io fills in its own default scale_factors when we
    # pass none, so build_pyramid=False does not actually mean "single scale").
    _retile_morphology_focus(sdata, path)

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

    # Write zarr cache for next time — into a staging directory, then rename.
    # Writing straight to cache_path with overwrite=True meant a failure part
    # way through destroyed the only copy, and the old error handler then
    # rmtree'd whatever was left.
    if use_cache:
        staging = cache_path.with_name(".sdata_cached__building.zarr")
        shutil.rmtree(staging, ignore_errors=True)
        try:
            _check_free_space(path, backup_path)
            _convert_arrow_strings(sdata)
            print(f"Writing zarr cache to {cache_path} ...")
            # Element by element rather than `sdata.write()`: it caps dask's
            # concurrency only for the elements that need it (labels, points),
            # and it is what lets memory be released — and reported — between
            # elements.
            sdata_write.write_sdata(
                sdata, staging,
                progress_cb=lambda pct, msg: print(f"  [{pct}%] {msg}"),
                # To the terminal as well as the log: a cache build runs for
                # tens of minutes, and "how much memory is this using" is a
                # question people ask *while* it is happening.
                log=_print_and_log,
            )
            write_manifest(staging, experiment_path)
            if backup_path is not None:
                _copy_sidecars_and_session(backup_path, staging, user_data)
            if cache_path.exists():
                # Should be unreachable — every path that gets here moved the
                # old cache aside already. Move rather than delete anyway: the
                # whole point is that nothing removes a cache irreversibly.
                displaced = cache_path.with_name(
                    f"sdata_cached_prev_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zarr")
                shutil.move(str(cache_path), str(displaced))
                print(f"  Existing cache moved to {displaced.name}")
            os.rename(staging, cache_path)
            sdata = _reopen_written_cache(sdata, cache_path)
            # After the rebinding, not inside the helper: the caller's reference
            # is what keeps the build-time object alive, so trimming any earlier
            # trims nothing.
            mem_probe.release()
            print("Zarr cache written.")
            if backup_path is not None:
                shutil.rmtree(str(backup_path), ignore_errors=True)
                print("Cache rebuild and data restoration complete.")
        except Exception as e:
            shutil.rmtree(staging, ignore_errors=True)
            print(f"Warning: could not write zarr cache: {e}")
            if backup_path is not None:
                print(f"Your previous cache is preserved at:\n  {backup_path}")

    return sdata


def _check_free_space(path: Path, backup_path: Optional[Path]) -> None:
    """Refuse to start a rebuild that cannot finish.

    A rebuild writes a full second copy before the old one is released, so on a
    nearly-full disk it can fail part way — which used to take the original with
    it. Estimating from the preserved copy is the best signal available.
    """
    reference = backup_path
    if reference is None:
        candidates = sorted(path.glob("sdata_cached_prev_*.zarr"), reverse=True)
        reference = candidates[0] if candidates else None
    if reference is None or not reference.exists():
        return
    from palms.utils.cache_repair import describe_store, human_bytes

    needed = describe_store(reference)["size_bytes"]
    free = shutil.disk_usage(path).free
    if free < needed * 1.1:
        raise OSError(
            f"not enough free space to rebuild the cache: about "
            f"{human_bytes(needed)} needed, {human_bytes(free)} free"
        )


def load_umap(path: Path):
    """
    Load precomputed UMAP coordinates.

    Returns a DataFrame with columns ['UMAP_1', 'UMAP_2'] indexed by cell barcode,
    or None if no 'analysis/' folder exists (e.g. a Crop Dataset export, which
    has no raw Xenium analysis outputs — the UMAP tab falls back to whatever
    coordinates are already embedded in adata.obsm['X_umap'], if any).
    Note: the UMAP has 91 fewer cells than the AnnData — handled with reindex.
    """
    analysis_path = path / "analysis"
    umap_path = analysis_path / "umap" / "gene_expression_2_components" / "projection.csv"
    if not umap_path.exists():
        print(f"No UMAP projection found at {umap_path} (no 'analysis/' folder in this dataset).")
        return None
    umap_df = pd.read_csv(umap_path, index_col=0)
    umap_df.columns = ["UMAP_1", "UMAP_2"]
    print(f"Loaded UMAP: {umap_df.shape[0]} cells")
    return umap_df


#: Directory-name prefix 10x uses under ``analysis/clustering/``.
_TENX_CLUSTERING_DIR_PREFIX = "gene_expression_"

#: What 10x's own clustering names look like. Only consulted when the CSVs are
#: not on disk to be asked — see :func:`is_native_clustering`.
_TENX_CLUSTERING_NAME = re.compile(r"^(graphclust|kmeans_\d+_clusters)$")


def native_clustering_names(path) -> set:
    """The clusterings 10x shipped with this dataset, read from disk.

    Empty for a dataset with no ``analysis/`` folder, which is not the same
    statement as "it has none" — a Crop Dataset export carries the obs columns
    into its cache without copying the CSVs they came from.
    """
    root = Path(path) / "analysis" / "clustering"
    if not root.is_dir():
        return set()
    return {d.name.replace(_TENX_CLUSTERING_DIR_PREFIX, "")
            for d in root.iterdir() if (d / "clusters.csv").exists()}


def is_native_clustering(key: str, path=None) -> bool:
    """Whether *key* is 10x's own clustering rather than one the viewer derived.

    The distinction matters wherever the two must not be conflated: a native
    clustering is an **input** that arrived with the dataset, so no recorded
    step produces it and a replayed notebook is not expected to. A
    viewer-derived clustering with no step behind it is a defect
    (``tests/test_clustering_recording.py``), and the two look identical in
    ``obs`` — both are ``clustering_<key>`` columns.

    Disk wins where it can answer: if this dataset has an ``analysis/clustering``
    folder, membership in it is definitive. Only when there is none does the
    name decide, because a Crop Dataset export drops that folder while keeping
    the columns, and refusing to answer there would make every crop export look
    like it had lost ten analyses.
    """
    if path is not None:
        on_disk = native_clustering_names(path)
        if on_disk:
            return key in on_disk
    return bool(_TENX_CLUSTERING_NAME.match(key))


def load_clusterings(path: Path):
    """
    Load all cluster assignments from analysis/clustering/.

    Returns a dict: {clustering_name -> pd.Series(cluster_id, index=cell_barcode)}.
    Returns {} if no 'analysis/' folder exists (e.g. a Crop Dataset export).
    """
    analysis_path = path / "analysis"
    clustering_root = analysis_path / "clustering"
    clusterings = {}
    if not clustering_root.exists():
        print(f"No clustering results found at {clustering_root} (no 'analysis/' folder in this dataset).")
        return clusterings
    for subdir in sorted(clustering_root.iterdir()):
        if not subdir.is_dir():
            continue
        csv_path = subdir / "clusters.csv"
        if not csv_path.exists():
            continue
        # strip leading "gene_expression_" prefix for display
        name = subdir.name.replace(_TENX_CLUSTERING_DIR_PREFIX, "")
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


def _build_parser():
    """The ``palms-build-cache`` argument parser, separated so it can be tested.

    The description is spelled out rather than taken from ``__doc__``: the
    module docstring describes the library, and putting it here is how the
    ``--help`` text came to advertise a script deleted in the 2026 refactor.
    """
    import argparse
    parser = argparse.ArgumentParser(
        prog="palms-build-cache",
        description=(
            "Build the SpatialData zarr cache (sdata_cached.zarr/) for a Xenium "
            "output directory, without starting the GUI. The viewer builds this "
            "on first launch anyway; running it here lets the slow first read "
            "happen over ssh or overnight. Complements palms-preprocess, which "
            "builds the separate per-gene transcript cache."
        ),
    )
    parser.add_argument("data_dir",
                        help="Path to Xenium output directory")
    parser.add_argument("--on-stale", choices=("ask",) + ON_STALE_CHOICES,
                        default="ask",
                        help="What to do when the existing cache is stale or will "
                             "not open. 'ask' (default) prompts if a GUI is "
                             "available and otherwise keeps the cache untouched, "
                             "so a rebuild from a terminal has to be asked for: "
                             "'rebuild' discards user data in the cache, 'restore' "
                             "rebuilds and carries it over, 'keep' loads it as-is.")
    parser.add_argument("--no-pyramid", action="store_true",
                        help="Do not build the morphology_focus image pyramid "
                             "(faster, but pan/zoom in the viewer will stutter)")
    parser.add_argument("--n-jobs", type=int, default=8,
                        help="Threads for spatialdata_io (default: 8)")
    parser.add_argument("--check", action="store_true",
                        help="Report the cache's status and exit without building "
                             "anything. Exits non-zero if it is missing, stale or "
                             "does not verify.")
    parser.add_argument("--no-cache", action="store_true",
                        help="Read the raw output without writing a cache — a dry "
                             "run of the read path, not a way to build one")
    return parser


def _check_cache(data_path: Path) -> int:
    """Print what is known about a cache without opening it. Returns an exit code.

    Every call here is read-only and filesystem-level (the same property
    ``cache_repair.verify`` relies on), so this works on a store too broken for
    zarr to open — which is exactly when someone runs it.
    """
    from palms.utils import cache_repair

    cache_path = data_path / "sdata_cached.zarr"
    experiment_path = data_path / "experiment.xenium"

    cache_only = is_cache_only(data_path)

    print(f"Dataset: {data_path}")
    if not cache_path.exists():
        print(f"No zarr cache at {cache_path}.")
        if cache_only:
            print("There is no raw Xenium output here either, so one cannot be "
                  "built. If this was a Crop Dataset export, its cache is the "
                  "only copy of the data — look for a sdata_cached_*.zarr "
                  "sibling before assuming it is lost.")
        else:
            print("Run this command without --check to build it.")
        return 1

    print(f"Cache:   {cache_path}")
    ok = True

    stale, certain = _is_cache_stale(cache_path, experiment_path)
    if cache_only:
        # Not "unknown": the answer is known and is "not applicable". Reporting
        # this cache as stale would be advising a rebuild that cannot happen.
        print("Freshness: n/a — cache-only dataset (no raw Xenium output beside "
              "it, likely a Crop Dataset export). The cache is the source of "
              "record and is never rebuilt.")
    elif not experiment_path.exists():
        print("Freshness: unknown — no experiment.xenium to compare against.")
    elif stale and certain:
        print("Freshness: STALE — experiment.xenium has changed since the build.")
        ok = False
    elif stale:
        print("Freshness: possibly stale — no manifest, so this is an mtime "
              "comparison, which copying the dataset also trips.")
        ok = False
    else:
        print("Freshness: up to date.")

    report = cache_repair.verify(cache_path)
    print("\nIntegrity:")
    print(report.summary())
    if not report.ok:
        ok = False

    user_data = _detect_user_data(cache_path)
    if _has_any_user_data(user_data):
        # Not a duplicate of the sidecar list above: this is the "what would
        # --on-stale rebuild throw away" view, the same one the GUI dialog shows.
        print("\nUser-generated data in this cache:")
        print(_format_user_data_message(user_data))
    else:
        print("\nNo user-generated data in this cache.")

    return 0 if ok else 1


def main():
    parser = _build_parser()
    args = parser.parse_args()

    data_path = Path(args.data_dir)
    if not data_path.is_dir():
        sys.exit(f"Error: {data_path} is not a directory.")

    if args.check:
        sys.exit(_check_cache(data_path))

    sdata = load_sdata(
        data_path,
        build_pyramid=not args.no_pyramid,
        n_jobs=args.n_jobs,
        use_cache=not args.no_cache,
        on_stale=None if args.on_stale == "ask" else args.on_stale,
    )
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
    # load_umap returns None when there is no 'analysis/' folder — a 10x bundle
    # whose analysis.tar.gz was never extracted, or a Crop Dataset export. The
    # cache is already written by this point, so the summary must not turn that
    # into a traceback. load_clusterings returns {} in the same situation.
    if umap_df is not None:
        print(f"UMAP:    {umap_df.shape[0]} cells")
    else:
        print("UMAP:    not present (no 'analysis/' folder)")
    print(f"Clusterings: {list(clusterings.keys())}")


if __name__ == "__main__":
    main()
