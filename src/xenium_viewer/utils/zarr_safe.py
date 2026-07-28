"""Crash-safe writes into the SpatialData zarr cache.

The viewer used to persist elements with ``delete_element_from_disk`` followed by
``write_element``. That is not a metadata operation: spatialdata's
``delete_element_from_disk`` does ``del root[element_type][element_name]``, which
recursively unlinks — **the bytes are gone before the replacement starts being
written**. spatialdata's own docstring says so ("data loss may occur if the
execution is interrupted during writing"). Since ``_persist_table`` ran on every
analysis action, every clustering, DEG run and label edit opened a window in
which killing the app, filling the disk, or any exception left the store
structurally invalid — and the loader then discarded the whole cache.

This module replaces that with stage-then-swap:

1. Write the new element into a *throwaway sibling store* under
   ``<cache>/.xv_staging/``. The live store is untouched, so a failure here
   costs nothing.
2. Under the store lock, drop a journal file, then swap with two ``os.rename``
   calls — microseconds, and same-filesystem because staging lives inside the
   cache.
3. Consolidate, then retire the journal and prune the trash.

The old element is *moved to* ``<cache>/.xv_trash/`` rather than deleted, so the
previous good copy is available to the repair tool.

Three facts about spatialdata 0.8 / zarr 3.1 make this work, all verified
against the installed packages rather than assumed:

* ``SpatialData.write(..., update_sdata_path=False)`` leaves ``tmp.path is None``
  at write time, so the "target is inside a store in use" guard never fires and
  staging can live inside the cache.
* ``read_zarr`` skips dot-prefixed members, so the staging, trash and journal
  directories are invisible to readers.
* Element *content* is read from the element's own ``zarr.json``, not from the
  root consolidated metadata, so a swapped-in element reads correctly even
  before consolidation. That makes the post-rename window benign.

Recovery infers what happened purely from the filesystem, so the journal is
written once and never updated mid-flight. See :func:`recover_pending`.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import threading
import time
import uuid
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Optional

log = logging.getLogger(__name__)

# Element type directories in a SpatialData zarr store.
ELEMENT_TYPES = ("images", "labels", "points", "shapes", "tables")

# All dot-prefixed so ``read_zarr`` ignores them (spatialdata _io/io_zarr.py).
STAGING_DIR = ".xv_staging"
TRASH_DIR = ".xv_trash"
JOURNAL_DIR = ".xv_journal"
LOCK_FILE = ".xv.lock"
MANIFEST_FILE = ".xv_manifest.json"

INTERNAL_NAMES = frozenset(
    {STAGING_DIR, TRASH_DIR, JOURNAL_DIR, LOCK_FILE, MANIFEST_FILE}
)

# Don't keep a backup copy of an element larger than this — the point of the
# trash is cheap insurance on the 320 MB table, not a second copy of a 20 GB
# image pyramid.
DEFAULT_MAX_TRASH_BYTES = 2 * 1024 ** 3


class ZarrSafeError(RuntimeError):
    """A safe-write operation could not be completed."""


# ── locking ──────────────────────────────────────────────────────────────────

# Reentrant: save_custom_seg_to_sdata holds the lock and then calls
# _persist_custom_table, which takes it again. With a plain Lock that only works
# by careful nesting; an RLock makes it robust.
_STORE_LOCK = threading.RLock()

_flock_state: dict[str, Any] = {"depth": 0, "fd": None, "path": None}

try:  # pragma: no cover - platform dependent
    import fcntl
except ImportError:  # Windows
    fcntl = None  # type: ignore[assignment]


@contextmanager
def store_lock(cache_path: Optional[Path] = None, timeout: float = 120.0) -> Iterator[None]:
    """Serialise writes to *cache_path*, in-process and across processes.

    The in-process ``RLock`` is what actually prevents the interleaved
    ``rename`` + ``consolidate_metadata`` that used to drop a live element from
    the root metadata. The ``flock`` on top guards against a second viewer
    instance on the same dataset; it is advisory and the kernel releases it if
    the process dies, so there is no stale-lock problem.
    """
    _STORE_LOCK.acquire()
    took_flock = False
    try:
        if fcntl is not None and cache_path is not None:
            took_flock = _acquire_flock(Path(cache_path), timeout)
        yield
    finally:
        if took_flock:
            _release_flock()
        _STORE_LOCK.release()


def _acquire_flock(cache_path: Path, timeout: float) -> bool:
    """Take the cross-process lock, or re-enter it if already held."""
    if _flock_state["depth"] > 0:
        _flock_state["depth"] += 1
        return True
    lock_path = cache_path / LOCK_FILE
    try:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(lock_path), os.O_RDWR | os.O_CREAT, 0o644)
    except OSError as exc:
        # A read-only dataset directory is a legitimate configuration; the
        # in-process lock still applies and the write itself will fail with a
        # clearer error than we could produce here.
        log.debug("could not open cache lock file %s: %s", lock_path, exc)
        return False
    deadline = time.monotonic() + timeout
    while True:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            break
        except OSError:
            if time.monotonic() >= deadline:
                os.close(fd)
                raise ZarrSafeError(
                    f"another process is writing to {cache_path} "
                    f"(waited {timeout:.0f}s for the cache lock)"
                ) from None
            time.sleep(0.05)
    _flock_state.update(depth=1, fd=fd, path=lock_path)
    return True


def _release_flock() -> None:
    _flock_state["depth"] -= 1
    if _flock_state["depth"] > 0:
        return
    fd = _flock_state["fd"]
    if fd is not None:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)
    _flock_state.update(depth=0, fd=None, path=None)


# ── primitives ───────────────────────────────────────────────────────────────

def consolidate(cache_path: Path) -> None:
    """Rebuild the root consolidated metadata, suppressing sidecar warnings.

    The viewer drops non-zarr files into the store root (h5ad caches, parquet
    DEG results, CopyKAT JSON). zarr emits one ``ZarrUserWarning`` per
    unrecognised member while walking the hierarchy. spatialdata suppresses
    these in its own consolidation path; the viewer's bare call in ``app.py``
    did not, which is the most likely source of the reported "several warnings".
    """
    import zarr
    import zarr.errors

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=zarr.errors.ZarrUserWarning)
        zarr.consolidate_metadata(str(cache_path))


def atomic_json(path: Path, obj: Any) -> None:
    """Write JSON via a temp file + ``os.replace``, so readers never see a partial."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.{uuid.uuid4().hex}.tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    os.replace(tmp, path)


def cache_path_of(sdata) -> Path:
    if getattr(sdata, "path", None) is None:
        raise ZarrSafeError("SpatialData object is not backed by a zarr store")
    return Path(sdata.path)


def element_type_of(sdata, name: str, element: Any = None) -> str:
    """Which of images/labels/points/shapes/tables *name* belongs to."""
    for etype in ELEMENT_TYPES:
        if name in getattr(sdata, etype, {}):
            return etype
    if element is not None:
        return _element_type_from_model(element)
    # Fall back to what is on disk, so an element present only in the store
    # (e.g. deleting something the in-memory object never held) still resolves.
    for etype in ELEMENT_TYPES:
        if (cache_path_of(sdata) / etype / name).exists():
            return etype
    raise ZarrSafeError(f"could not determine element type for {name!r}")


def _element_type_from_model(element: Any) -> str:
    import anndata

    if isinstance(element, anndata.AnnData):
        return "tables"
    from spatialdata.models import (
        Image2DModel, Image3DModel, Labels2DModel, Labels3DModel,
        PointsModel, ShapesModel, get_model,
    )
    model = get_model(element)
    if model in (Image2DModel, Image3DModel):
        return "images"
    if model in (Labels2DModel, Labels3DModel):
        return "labels"
    if model is PointsModel:
        return "points"
    if model is ShapesModel:
        return "shapes"
    raise ZarrSafeError(f"unsupported element model {model!r}")


def _assert_not_dask_backed(sdata, live: Path, name: str) -> None:
    """Refuse to move an element whose files back a live dask graph.

    Tables are read eagerly so this never trips for them, but images, labels and
    points are lazily backed — renaming their directory out from under a live
    napari layer would break it. ``delete_element_from_disk`` guards the same
    way; doing it explicitly here keeps it a property of the design rather than
    an accident of the call sites.
    """
    try:
        from spatialdata._io._utils import _backed_elements_contained_in_path
    except ImportError:  # pragma: no cover - upstream layout change
        return
    if any(_backed_elements_contained_in_path(path=live, object=sdata)):
        raise ZarrSafeError(
            f"cannot safely rewrite {name!r}: its files back a lazily-loaded "
            "element still in use. Load it into memory first."
        )


def _journal_path(cache_path: Path, uid: str) -> Path:
    return cache_path / JOURNAL_DIR / f"{uid}.json"


def _dir_size(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for fname in files:
            try:
                total += os.path.getsize(os.path.join(root, fname))
            except OSError:
                pass
    return total


def prune_trash(
    cache_path: Path,
    etype: str,
    name: str,
    keep: int = 1,
    max_bytes: int = DEFAULT_MAX_TRASH_BYTES,
) -> None:
    """Keep at most *keep* recent backups of an element, and drop huge ones."""
    trash_dir = cache_path / TRASH_DIR / etype
    if not trash_dir.exists():
        return
    entries = sorted(
        (p for p in trash_dir.iterdir() if p.name.startswith(f"{name}.")),
        key=lambda p: p.name,
        reverse=True,
    )
    for index, entry in enumerate(entries):
        drop = index >= keep
        if not drop and _dir_size(entry) > max_bytes:
            drop = True  # too big to be worth a safety copy
        if drop:
            shutil.rmtree(entry, ignore_errors=True)


# ── the safe write ───────────────────────────────────────────────────────────

def safe_write_element(
    sdata,
    name: str,
    element: Any = None,
    *,
    keep_backup: bool = True,
    max_backup_bytes: int = DEFAULT_MAX_TRASH_BYTES,
) -> None:
    """Write *element* into the store under *name*, without a loss window.

    Replaces ``sdata.delete_element_from_disk(name); sdata.write_element(name)``.
    Also updates the in-memory ``sdata`` so it stays consistent with disk.
    """
    from spatialdata import SpatialData

    cache_path = cache_path_of(sdata)
    if element is None:
        element = sdata[name]
    etype = element_type_of(sdata, name, element)

    live = cache_path / etype / name
    _assert_not_dask_backed(sdata, live, name)

    uid = uuid.uuid4().hex[:12]
    stage_store = cache_path / STAGING_DIR / f"{name}.{uid}.zarr"
    stage_element = stage_store / etype / name
    stamp = time.strftime("%Y%m%d_%H%M%S")
    trash = cache_path / TRASH_DIR / etype / f"{name}.{stamp}.{uid}"

    # 1. Stage. Deliberately outside the lock: this is the slow part (a full
    #    serialize of the element) and it cannot affect the live store, so
    #    holding the lock here would stall the GUI for no safety benefit.
    stage_store.parent.mkdir(parents=True, exist_ok=True)
    try:
        SpatialData(**{etype: {name: element}}).write(
            stage_store,
            overwrite=True,
            consolidate_metadata=False,   # root zarr.json is rewritten once, after the swap
            update_sdata_path=False,      # keeps tmp.path None, so the in-store guard is skipped
        )
    except Exception:
        shutil.rmtree(stage_store, ignore_errors=True)
        raise

    journal = _journal_path(cache_path, uid)
    record = {
        "version": 1,
        "op": "write_element",
        "element_type": etype,
        "name": name,
        "live": str(live.relative_to(cache_path)),
        "stage": str(stage_element.relative_to(cache_path)),
        "stage_store": str(stage_store.relative_to(cache_path)),
        "trash": str(trash.relative_to(cache_path)),
        "created": stamp,
    }
    journaled = False
    try:
        with store_lock(cache_path):
            # 2. Journal before the first destructive step, so recovery can tell
            #    a half-swap from an abandoned staging directory.
            atomic_json(journal, record)
            journaled = True

            # 3. Swap: two renames, microseconds, no I/O.
            if live.exists():
                trash.parent.mkdir(parents=True, exist_ok=True)
                os.rename(live, trash)
            live.parent.mkdir(parents=True, exist_ok=True)
            os.rename(stage_element, live)

            # 4. Commit.
            _set_in_memory(sdata, etype, name, element)
            consolidate(cache_path)
            journal.unlink(missing_ok=True)
            shutil.rmtree(stage_store, ignore_errors=True)
            if keep_backup:
                prune_trash(cache_path, etype, name, keep=1, max_bytes=max_backup_bytes)
            else:
                shutil.rmtree(trash, ignore_errors=True)
    except Exception:
        # An exception we can see is recoverable *now*: a running GUI must not
        # be left with an inconsistent store. Never delete staging while a
        # journal stands — it may hold the only copy of the new data, and
        # dropping it would downgrade a roll-forward into a roll-back.
        if journaled:
            _heal(cache_path, record, journal, stage_store)
        else:
            shutil.rmtree(stage_store, ignore_errors=True)
        raise


def _heal(cache_path: Path, record: dict, journal: Path, stage_store: Path) -> None:
    """Best-effort in-process replay after a failed swap, then clear the journal.

    If the replay itself fails the journal is deliberately left behind, so
    :func:`recover_pending` retries at the next launch.
    """
    try:
        _replay(cache_path, record)
        consolidate(cache_path)
    except Exception as exc:  # pragma: no cover - defensive
        log.warning(
            "could not heal an interrupted write of %s; leaving the journal for "
            "startup recovery: %s", record.get("name"), exc, exc_info=True,
        )
        return
    journal.unlink(missing_ok=True)
    shutil.rmtree(stage_store, ignore_errors=True)


def safe_delete_element(sdata, name: str, *, keep_backup: bool = True) -> None:
    """Remove *name* from the store, keeping a backup copy in the trash."""
    cache_path = cache_path_of(sdata)
    etype = element_type_of(sdata, name)
    live = cache_path / etype / name
    if not live.exists():
        _drop_in_memory(sdata, etype, name)
        return
    _assert_not_dask_backed(sdata, live, name)

    uid = uuid.uuid4().hex[:12]
    stamp = time.strftime("%Y%m%d_%H%M%S")
    trash = cache_path / TRASH_DIR / etype / f"{name}.{stamp}.{uid}"
    journal = _journal_path(cache_path, uid)

    with store_lock(cache_path):
        atomic_json(journal, {
            "version": 1,
            "op": "delete_element",
            "element_type": etype,
            "name": name,
            "live": str(live.relative_to(cache_path)),
            "trash": str(trash.relative_to(cache_path)),
            "created": stamp,
        })
        trash.parent.mkdir(parents=True, exist_ok=True)
        os.rename(live, trash)
        _drop_in_memory(sdata, etype, name)
        consolidate(cache_path)
        journal.unlink(missing_ok=True)
        if keep_backup:
            prune_trash(cache_path, etype, name, keep=1)
        else:
            shutil.rmtree(trash, ignore_errors=True)


def _set_in_memory(sdata, etype: str, name: str, element: Any) -> None:
    try:
        getattr(sdata, etype)[name] = element
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("could not update in-memory %s/%s: %s", etype, name, exc)


def _drop_in_memory(sdata, etype: str, name: str) -> None:
    try:
        container = getattr(sdata, etype)
        if name in container:
            del container[name]
    except Exception as exc:  # pragma: no cover - defensive
        log.debug("could not drop in-memory %s/%s: %s", etype, name, exc)


# ── plain zarr groups (viewer_session) ───────────────────────────────────────

@contextmanager
def safe_group_update(cache_path: Path, group_name: str) -> Iterator[Any]:
    """Edit a plain zarr group under the store root without a loss window.

    ``viewer_session`` is not a SpatialData element, so :func:`safe_write_element`
    does not apply — but it had the same shape of bug, and worse: ``save_session``
    destroyed the group with ``create_group(overwrite=True)`` and only wrote the
    replacement ~110 lines later.

    Yields a zarr group backed by staging, seeded from the current contents so
    callers that update only part of it keep working. On a clean exit the
    staging group is swapped in; on an exception the live group is untouched.
    """
    import zarr

    cache_path = Path(cache_path)
    uid = uuid.uuid4().hex[:12]
    stamp = time.strftime("%Y%m%d_%H%M%S")
    live = cache_path / group_name
    stage = cache_path / STAGING_DIR / f"{group_name}.{uid}"
    trash = cache_path / TRASH_DIR / "_groups" / f"{group_name}.{stamp}.{uid}"
    journal = _journal_path(cache_path, uid)

    stage.parent.mkdir(parents=True, exist_ok=True)
    if live.exists():
        shutil.copytree(live, stage)
    else:
        stage.mkdir(parents=True)

    try:
        group = zarr.open_group(str(stage), mode="a", use_consolidated=False)
        yield group
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise

    try:
        with store_lock(cache_path):
            atomic_json(journal, {
                "version": 1,
                "op": "group_update",
                "name": group_name,
                "live": str(live.relative_to(cache_path)),
                "stage": str(stage.relative_to(cache_path)),
                "trash": str(trash.relative_to(cache_path)),
                "created": stamp,
            })
            if live.exists():
                trash.parent.mkdir(parents=True, exist_ok=True)
                os.rename(live, trash)
            os.rename(stage, live)
            consolidate(cache_path)
            journal.unlink(missing_ok=True)
            prune_trash(cache_path, "_groups", group_name, keep=1)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise


# ── recovery ─────────────────────────────────────────────────────────────────

def recover_pending(cache_path: Path) -> list[str]:
    """Finish or undo any interrupted safe write, and clear debris.

    Called at startup before the store is opened. The phase is inferred from
    the filesystem rather than from journal bookkeeping, so a journal written
    once is enough:

    ==================================== ==========================
    State                                Action
    ==================================== ==========================
    journal, live present, stage present crashed while staging → drop staging
    journal, live *missing*, stage       **roll forward**
    journal, live missing, stage gone,   **roll back** from trash
    trash present
    journal, live present, stage gone    crashed before cleanup → just clean
    no journal, orphan staging/.partial  garbage-collect
    ==================================== ==========================

    Returns a human-readable list of what it did (empty if nothing was pending).
    """
    cache_path = Path(cache_path)
    actions: list[str] = []
    journal_dir = cache_path / JOURNAL_DIR
    changed = False

    if journal_dir.exists():
        for journal in sorted(journal_dir.glob("*.json")):
            try:
                record = json.loads(journal.read_text())
            except (OSError, ValueError) as exc:
                log.warning("discarding unreadable journal %s: %s", journal, exc)
                journal.unlink(missing_ok=True)
                continue
            try:
                acted = _replay(cache_path, record)
            except Exception as exc:  # pragma: no cover - defensive
                log.warning("could not replay journal %s: %s", journal, exc, exc_info=True)
                continue
            if acted:
                actions.append(acted)
                changed = True
            journal.unlink(missing_ok=True)
            stage_store = record.get("stage_store") or record.get("stage")
            if stage_store:
                shutil.rmtree(cache_path / stage_store, ignore_errors=True)

    # Debris left by hard kills: abandoned staging trees and zarr's own
    # partial-write temp files. The trash is deliberately NOT cleaned here —
    # it is the recovery source for the repair tool.
    staging_dir = cache_path / STAGING_DIR
    if staging_dir.exists():
        for entry in staging_dir.iterdir():
            shutil.rmtree(entry, ignore_errors=True) if entry.is_dir() else entry.unlink(missing_ok=True)
            actions.append(f"removed abandoned staging {entry.name}")
    for partial in cache_path.rglob("*.partial"):
        try:
            partial.unlink()
            actions.append(f"removed partial write {partial.name}")
        except OSError:
            pass

    if changed:
        try:
            consolidate(cache_path)
        except Exception as exc:  # pragma: no cover - defensive
            log.warning("could not consolidate after recovery: %s", exc, exc_info=True)
    return actions


def _replay(cache_path: Path, record: dict) -> Optional[str]:
    """Apply one journal record. Returns a description, or None if nothing to do."""
    name = record.get("name", "?")
    live = cache_path / record["live"]
    stage = cache_path / record["stage"] if record.get("stage") else None
    trash = cache_path / record["trash"] if record.get("trash") else None
    op = record.get("op", "write_element")

    if op == "delete_element":
        # The rename is the whole operation; if live is gone it succeeded.
        if live.exists() and trash is not None and trash.exists():
            shutil.rmtree(trash, ignore_errors=True)
            return None
        return None

    if live.exists():
        # Either we crashed while staging, or after the swap but before cleanup.
        # Both are already-correct states.
        return None

    if stage is not None and stage.exists():
        live.parent.mkdir(parents=True, exist_ok=True)
        os.rename(stage, live)
        return f"recovered {name}: completed an interrupted write"

    if trash is not None and trash.exists():
        live.parent.mkdir(parents=True, exist_ok=True)
        os.rename(trash, live)
        return f"recovered {name}: rolled back an interrupted write"

    return f"WARNING: {name} was lost by an interrupted write and has no backup"


def list_trash(cache_path: Path) -> dict[str, list[Path]]:
    """Map ``<type>/<name>`` to the backup copies available for recovery."""
    cache_path = Path(cache_path)
    trash_root = cache_path / TRASH_DIR
    found: dict[str, list[Path]] = {}
    if not trash_root.exists():
        return found
    for type_dir in sorted(trash_root.iterdir()):
        if not type_dir.is_dir():
            continue
        for entry in sorted(type_dir.iterdir(), reverse=True):
            # "<name>.<stamp>.<uid>" — the name may itself contain dots.
            base = entry.name.rsplit(".", 2)[0]
            found.setdefault(f"{type_dir.name}/{base}", []).append(entry)
    return found
