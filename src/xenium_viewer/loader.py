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
import shutil
import sys
from pathlib import Path
from typing import Optional

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
        }
        for key in user_data["uns_keys"]:
            lines.append(f"  • {labels.get(key, key)}")
    _sidecar_labels = {
        "roi_deg_cache.parquet": "ROI DEG results",
        "arms_tile_deg_cache.parquet": "ARMS tile DEG results",
        "adata_norm_cache.h5ad": "Normalized expression cache",
    }
    for fname in user_data["sidecars"]:
        if fname.startswith("adata_cnv_cache_") or fname.startswith("cnv_"):
            backend = "CopyKAT" if "copykat" in fname else "inferCNV"
            label = f"{backend} CNV results (hours of compute)"
        else:
            label = _sidecar_labels.get(fname, fname)
        lines.append(f"  • {label}")
    return "\n".join(lines)


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
    msg.exec_()

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


def _ask_corrupt_cache(error: Exception, report, user_data: dict) -> str:
    """Ask what to do about a cache that will not open even after repair.

    Returns 'restore', 'rebuild' or 'quit'. Never returns a destructive default:
    with no way to prompt we raise instead, because the previous behaviour —
    silently rebuilding 30 GB — is the thing being fixed.
    """
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
    msg.exec_()

    clicked = msg.clickedButton()
    if clicked is restore_btn:
        return "restore"
    if clicked is quit_btn:
        return "quit"
    return "rebuild"


def _open_cache(cache_path: Path):
    """Open the cache, repairing first and escalating if that is not enough.

    Returns the SpatialData, or raises the last read error. Repair is attempted
    before *and* after a plain read because the common corruption — an element
    present on disk but absent from the consolidated metadata — makes the read
    fail while costing nothing to fix.
    """
    import spatialdata

    from xenium_viewer.utils import cache_repair

    report = cache_repair.verify(cache_path)
    if not report.ok:
        result = cache_repair.repair(cache_path, report, level=cache_repair.AUTO)
        for action in result.actions:
            print(f"  Cache: {action}")

    try:
        return spatialdata.read_zarr(str(cache_path))
    except Exception as first_error:
        # Still broken. If a previous version of a missing element is sitting in
        # the trash, putting it back is strictly better than discarding 30 GB.
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


def write_manifest(cache_path: Path, experiment_path: Path) -> None:
    """Record what this cache was built from, so staleness can be exact."""
    from datetime import datetime, timezone

    from xenium_viewer.utils.zarr_safe import MANIFEST_FILE, atomic_json

    manifest = {"built_at": datetime.now(timezone.utc).isoformat()}
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
    from xenium_viewer.utils.cache_repair import read_manifest

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
        stale, certain = _is_cache_stale(cache_path, experiment_path)
        user_data = _detect_user_data(cache_path)
        preference = None

        if stale:
            if not certain and not _has_any_user_data(user_data):
                # Pre-manifest cache, mtime-only signal, nothing to lose:
                # rebuilding is cheap insurance and stamps a manifest.
                print("Zarr cache may be stale (no manifest; experiment.xenium is "
                      "newer). Rebuilding...")
            elif not certain:
                print("Zarr cache may be stale — the check is a file timestamp, "
                      "which a copy or re-download also changes.")
                preference = _ask_rebuild_preference(user_data, certain=False)
            else:
                preference = (_ask_rebuild_preference(user_data)
                              if _has_any_user_data(user_data) else "rebuild")
                if preference == "rebuild":
                    print("Zarr cache is stale (experiment.xenium has changed). "
                          "Rebuilding...")

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
                print(sdata)
                return sdata
            except Exception as e:
                from xenium_viewer.utils import cache_repair
                choice = _ask_corrupt_cache(e, cache_repair.verify(cache_path), user_data)
                if choice == "quit":
                    raise CacheLoadAborted(
                        "Cancelled at the user's request; the cache was left untouched."
                    ) from e
                preference = choice

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
            sdata.write(str(staging))
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
    from xenium_viewer.utils.cache_repair import describe_store, human_bytes

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
