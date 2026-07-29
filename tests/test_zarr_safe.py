"""Crash-safe zarr writes: the store must survive an interrupted write.

The reported bug: the viewer persisted the cell table with
``delete_element_from_disk`` then ``write_element``, on every analysis action.
That unlinks the on-disk table *before* the replacement exists, so any
interruption left the store structurally invalid — and the loader then threw the
whole 30 GB cache away.

``test_live_store_untouched_during_staging`` and the ``_fail_nth_rename`` tests
below are the regression tests for that. They interrupt a write by making the
staging write, or a specific ``os.rename``, raise — with ``OSError`` for a
failure the running app can see and heal, and ``KeyboardInterrupt`` (a
``BaseException``, so no handler runs) for what a ``kill -9`` actually leaves
on disk.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_zarr_safe.py
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")

from xenium_viewer.utils import zarr_safe  # noqa: E402
from xenium_viewer.utils.zarr_safe import (  # noqa: E402
    JOURNAL_DIR, STAGING_DIR, TRASH_DIR, ZarrSafeError,
    consolidate, list_trash, recover_pending, safe_delete_element,
    safe_write_element, store_lock,
)


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _reread(cache):
    from spatialdata import read_zarr
    return read_zarr(str(cache))


def _marker(sdata):
    return list(sdata["table"].obs["marker"])[0]


# ── the happy path ───────────────────────────────────────────────────────────

def test_safe_write_element_roundtrip(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))

    assert _marker(_reread(cache)) == "NEW"
    # the other element is untouched and still enumerable
    assert sorted(_reread(cache).labels) == ["lab"]
    # in-memory object followed the write
    assert _marker(tiny_sdata) == "NEW"


def test_no_journal_or_staging_survives_a_clean_write(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))

    assert list((cache / JOURNAL_DIR).glob("*.json")) == []
    assert list((cache / STAGING_DIR).iterdir()) == []


def test_previous_version_is_kept_in_the_trash(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))

    backups = list_trash(cache)
    assert "tables/table" in backups
    assert len(backups["tables/table"]) == 1


def test_trash_is_pruned_to_one_copy(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    for value in ("V1", "V2", "V3"):
        safe_write_element(tiny_sdata, "table", make_table(value))
    assert len(list_trash(cache)["tables/table"]) == 1


def test_oversized_backups_are_not_kept(tiny_sdata, make_table):
    """A 20 GB image pyramid must not get a shadow copy on a 96%-full disk."""
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"), max_backup_bytes=0)
    assert list_trash(cache).get("tables/table", []) == []


def test_writing_the_first_element_of_a_type_creates_its_group(tiny_sdata):
    """A store with no shapes has no ``shapes/`` zarr group.

    Renaming into a plain directory leaves the element physically present but
    invisible to ``read_zarr``, because the consolidation walk never descends
    into a non-group. Caught by the ROI persistence test, pinned here.
    """
    pytest.importorskip("geopandas")
    import geopandas as gpd
    from shapely.geometry import Polygon
    from spatialdata.models import ShapesModel

    cache = Path(tiny_sdata.path)
    assert not (cache / "shapes").exists()

    gdf = ShapesModel.parse(gpd.GeoDataFrame(geometry=[Polygon([(0, 0), (4, 0), (4, 4)])]))
    safe_write_element(tiny_sdata, "rois", gdf)

    assert (cache / "shapes" / "zarr.json").exists()
    assert sorted(_reread(cache).shapes) == ["rois"]


def test_safe_delete_element(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "extra", make_table("X"))
    assert "extra" in _reread(cache).tables

    safe_delete_element(tiny_sdata, "extra")
    assert "extra" not in _reread(cache).tables
    assert "extra" not in tiny_sdata.tables


# ── the regression tests ─────────────────────────────────────────────────────

def test_live_store_untouched_during_staging(tiny_sdata, make_table, monkeypatch):
    """The reported bug, in miniature: a failed write must cost nothing.

    Under delete-then-write, a failure here left no table on disk at all.
    """
    from spatialdata import SpatialData

    cache = Path(tiny_sdata.path)
    before = sorted(p.name for p in (cache / "tables" / "table").rglob("*"))

    real_write = SpatialData.write

    def boom(self, *args, **kwargs):
        real_write(self, *args, **kwargs)
        raise RuntimeError("disk full")

    monkeypatch.setattr(SpatialData, "write", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    assert _marker(_reread(cache)) == "OLD"
    assert sorted(p.name for p in (cache / "tables" / "table").rglob("*")) == before
    assert list((cache / STAGING_DIR).iterdir()) == []


def _fail_nth_rename(monkeypatch, n, exc=None):
    """Make the nth os.rename raise.

    ``OSError`` models a failure the running app can see and heal from (disk
    full, EPERM). ``KeyboardInterrupt`` — a ``BaseException``, so it sails past
    ``except Exception`` — models an actual kill: no cleanup runs, and the crash
    state is left on disk for startup recovery. Both matter, and they take
    different paths.
    """
    calls = {"n": 0}
    real = os.rename
    error = exc or OSError("simulated interruption")

    def flaky(src, dst):
        calls["n"] += 1
        if calls["n"] == n:
            raise error
        return real(src, dst)

    monkeypatch.setattr(zarr_safe.os, "rename", flaky)
    return calls


def test_interrupted_before_swap_leaves_the_store_correct(tiny_sdata, make_table, monkeypatch):
    """Failed before the first rename: nothing changed."""
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 1)
    with pytest.raises(OSError):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    assert _marker(_reread(cache)) == "OLD"
    recover_pending(cache)
    assert _marker(_reread(cache)) == "OLD"


def test_failure_between_renames_self_heals_in_process(tiny_sdata, make_table, monkeypatch):
    """A failure the app can see must never leave the GUI on a broken store.

    The swap already removed the live element, so the error handler completes
    the write rather than unwinding it — the new data is the only copy that
    exists at that instant.
    """
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2)
    with pytest.raises(OSError):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    assert _marker(_reread(cache)) == "NEW"
    assert list((cache / JOURNAL_DIR).glob("*.json")) == []


def test_kill_between_renames_rolls_forward_at_startup(tiny_sdata, make_table, monkeypatch):
    """The genuinely dangerous instant, with no chance to clean up.

    KeyboardInterrupt bypasses the handler, so this leaves exactly the on-disk
    state a ``kill -9`` would.
    """
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2, exc=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    # The store is genuinely broken at this point — that is the premise.
    assert not (cache / "tables" / "table").exists()

    actions = recover_pending(cache)
    assert any("completed an interrupted write" in a for a in actions), actions
    assert _marker(_reread(cache)) == "NEW"


def test_kill_between_renames_rolls_back_when_staging_is_gone(
    tiny_sdata, make_table, monkeypatch,
):
    """Same instant, but the staging tree did not survive: fall back to the trash."""
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2, exc=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    shutil.rmtree(cache / STAGING_DIR)
    actions = recover_pending(cache)
    assert any("rolled back" in a for a in actions), actions
    assert _marker(_reread(cache)) == "OLD"


def test_store_is_readable_after_the_swap_before_consolidation(tiny_sdata, make_table, monkeypatch):
    """Element content comes from the element's own metadata, not the root.

    This is why the post-rename window is benign, and it is load-bearing for the
    whole design — pin it so an upstream change is caught here.
    """
    cache = Path(tiny_sdata.path)
    monkeypatch.setattr(zarr_safe, "consolidate", lambda p: (_ for _ in ()).throw(
        RuntimeError("consolidate interrupted")))
    with pytest.raises(RuntimeError, match="consolidate interrupted"):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    assert _marker(_reread(cache)) == "NEW"


def test_recovery_is_idempotent(tiny_sdata, make_table, monkeypatch):
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2, exc=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    first = recover_pending(cache)
    second = recover_pending(cache)
    assert first and second == []
    assert _marker(_reread(cache)) == "NEW"


def test_recovery_reports_unrecoverable_loss_instead_of_hiding_it(tiny_sdata, make_table, monkeypatch):
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2, exc=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    shutil.rmtree(cache / STAGING_DIR)
    shutil.rmtree(cache / TRASH_DIR)
    actions = recover_pending(cache)
    assert any("lost by an interrupted write" in a for a in actions), actions


def test_recovery_clears_debris(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / STAGING_DIR).mkdir(parents=True, exist_ok=True)
    (cache / STAGING_DIR / "abandoned.zarr").mkdir()
    (cache / "tables" / "zarr.json.deadbeef.partial").write_text("{}")

    actions = recover_pending(cache)
    assert any("abandoned" in a for a in actions)
    assert any("partial" in a for a in actions)
    assert not (cache / "tables" / "zarr.json.deadbeef.partial").exists()


def test_unreadable_journal_is_discarded_not_fatal(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / JOURNAL_DIR).mkdir(parents=True, exist_ok=True)
    (cache / JOURNAL_DIR / "garbage.json").write_text("{not json")

    assert recover_pending(cache) == []
    assert list((cache / JOURNAL_DIR).glob("*.json")) == []


# ── consolidation and store hygiene ──────────────────────────────────────────

def test_consolidate_suppresses_sidecar_warnings(tiny_sdata):
    """The viewer drops non-zarr files in the store root; zarr warns per file.

    app.py called consolidate_metadata bare, so those surfaced to the user —
    the most likely source of the reported "several warnings".
    """
    cache = Path(tiny_sdata.path)
    (cache / "adata_norm_cache.h5ad").write_bytes(b"not zarr")
    (cache / "roi_deg_cache.parquet").write_bytes(b"not zarr")
    (cache / "cnv_copykat_result.json").write_text("{}")

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        consolidate(cache)


def test_sidecars_survive_a_swap(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    sidecar = cache / "adata_norm_cache.h5ad"
    sidecar.write_bytes(b"payload")

    safe_write_element(tiny_sdata, "table", make_table("NEW"))
    assert sidecar.read_bytes() == b"payload"


def test_internal_dirs_are_invisible_to_read_zarr(tiny_sdata, make_table):
    """Staging/trash/journal are dot-prefixed so the reader skips them."""
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))
    assert (cache / TRASH_DIR).exists()

    sdata = _reread(cache)
    assert sorted(sdata.tables) == ["table"]
    assert sorted(sdata.labels) == ["lab"]


def test_refuses_to_move_a_dask_backed_element(tiny_sdata):
    """Renaming files out from under a lazily-loaded layer would break it."""
    cache = Path(tiny_sdata.path)
    lazy = _reread(cache)          # labels come back dask-backed
    with pytest.raises(ZarrSafeError, match="lazily-loaded"):
        safe_write_element(lazy, "lab", lazy["lab"])


# ── grouped (non-element) updates ────────────────────────────────────────────

def test_safe_group_update_commits(tiny_sdata):
    cache = Path(tiny_sdata.path)
    with zarr_safe.safe_group_update(cache, "viewer_session") as group:
        group.attrs["hello"] = "world"

    import zarr
    reopened = zarr.open_group(str(cache / "viewer_session"), mode="r", use_consolidated=False)
    assert reopened.attrs["hello"] == "world"


def test_safe_group_update_leaves_the_group_intact_on_failure(tiny_sdata):
    """save_session destroyed viewer_session ~110 lines before rewriting it."""
    import zarr

    cache = Path(tiny_sdata.path)
    with zarr_safe.safe_group_update(cache, "viewer_session") as group:
        group.attrs["keep"] = "original"

    with pytest.raises(RuntimeError, match="boom"):
        with zarr_safe.safe_group_update(cache, "viewer_session") as group:
            group.attrs["keep"] = "clobbered"
            raise RuntimeError("boom")

    reopened = zarr.open_group(str(cache / "viewer_session"), mode="r", use_consolidated=False)
    assert reopened.attrs["keep"] == "original"


# ── locking ──────────────────────────────────────────────────────────────────

def test_store_lock_is_reentrant(tiny_sdata):
    """save_custom_seg_to_sdata takes the lock then calls _persist_custom_table."""
    cache = Path(tiny_sdata.path)
    with store_lock(cache):
        with store_lock(cache):
            pass


def test_write_takes_the_lock_while_swapping(tiny_sdata, make_table):
    """A concurrent consolidate must not observe a half-swapped tree."""
    seen = []
    real_consolidate = zarr_safe.consolidate

    def watched(path):
        seen.append(zarr_safe._STORE_LOCK._is_owned())
        return real_consolidate(path)

    zarr_safe.consolidate = watched
    try:
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    finally:
        zarr_safe.consolidate = real_consolidate
    assert seen == [True]


def test_journal_records_relative_paths(tiny_sdata, make_table, monkeypatch):
    """So a cache that gets moved can still be recovered."""
    cache = Path(tiny_sdata.path)
    _fail_nth_rename(monkeypatch, 2, exc=KeyboardInterrupt())
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    journal = next((cache / JOURNAL_DIR).glob("*.json"))
    record = json.loads(journal.read_text())
    for key in ("live", "stage", "trash"):
        assert not Path(record[key]).is_absolute()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
