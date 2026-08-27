"""Cache verification and repair.

The loader used to equate "read_zarr raised" with "this cache is worthless" and
move it aside (or delete it). The most common failure is a consolidated-metadata
entry that no longer matches disk, which re-consolidating fixes in a second —
versus a 30 GB rebuild from the raw Xenium output.

``verify`` must work on a store too broken for zarr to open, and must never
mutate anything.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_cache_repair.py
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

from palms.utils import cache_repair  # noqa: E402
from palms.utils.cache_repair import (  # noqa: E402
    AUTO, FULL, human_bytes, repair, verify,
)
from palms.utils.zarr_safe import (  # noqa: E402
    JOURNAL_DIR, STAGING_DIR, consolidate, safe_write_element,
)


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _reread(cache):
    from spatialdata import read_zarr
    return read_zarr(str(cache))


# ── a healthy store ──────────────────────────────────────────────────────────

def test_a_fresh_cache_is_healthy(tiny_sdata):
    report = verify(Path(tiny_sdata.path))
    assert report.ok
    assert report.repairable
    assert sorted(report.on_disk) == ["labels/lab", "tables/table"]
    assert "healthy" in report.summary()


def test_verify_reports_a_missing_cache(tmp_path):
    report = verify(tmp_path / "nope.zarr")
    assert not report.exists
    assert not report.ok
    assert "No cache" in report.summary()


def test_verify_never_mutates(tiny_sdata):
    cache = Path(tiny_sdata.path)
    before = sorted(str(p.relative_to(cache)) for p in cache.rglob("*"))
    verify(cache)
    assert sorted(str(p.relative_to(cache)) for p in cache.rglob("*")) == before


def test_sidecars_are_reported_but_are_not_faults(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"hours of compute")
    (cache / "roi_deg_cache.parquet").write_bytes(b"x")

    report = verify(cache)
    assert report.ok
    assert "adata_cnv_cache_copykat.h5ad" in report.sidecars


# ── the fixable class: on disk, absent from metadata ─────────────────────────

def test_verify_detects_metadata_that_lost_an_element(tiny_sdata):
    """Exactly the state an interleaved consolidate used to produce."""
    cache = Path(tiny_sdata.path)
    stash = cache.parent / "stash"
    os.rename(cache / "tables" / "table", stash)
    consolidate(cache)                       # metadata now omits the table
    os.rename(stash, cache / "tables" / "table")

    report = verify(cache)
    assert report.missing_in_meta == ["tables/table"]
    assert report.missing_on_disk == []
    assert not report.ok
    assert report.repairable


def test_repair_reconsolidates_and_the_store_opens_again(tiny_sdata):
    cache = Path(tiny_sdata.path)
    stash = cache.parent / "stash"
    os.rename(cache / "tables" / "table", stash)
    consolidate(cache)
    os.rename(stash, cache / "tables" / "table")
    assert "table" not in _reread(cache).tables      # broken as far as callers see

    result = repair(cache)
    assert result.changed and not result.failures
    assert "table" in _reread(cache).tables
    assert verify(cache).ok


def test_repair_is_idempotent(tiny_sdata):
    cache = Path(tiny_sdata.path)
    repair(cache)
    second = repair(cache)
    assert not second.failures
    assert verify(cache).ok


# ── the serious class: in metadata, absent from disk ─────────────────────────

def test_verify_detects_an_element_missing_from_disk(tiny_sdata):
    cache = Path(tiny_sdata.path)
    shutil.rmtree(cache / "tables" / "table")

    report = verify(cache)
    assert report.missing_on_disk == ["tables/table"]
    assert not report.ok
    assert not report.repairable          # nothing to restore from
    assert "missing from disk" in report.summary()


def test_full_repair_restores_a_lost_element_from_the_trash(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))   # seeds .xv_trash
    shutil.rmtree(cache / "tables" / "table")

    report = verify(cache)
    assert report.missing_on_disk == ["tables/table"]
    assert report.repairable                                     # a backup exists

    result = repair(cache, report, level=FULL)
    assert "restored tables/table from backup" in result.actions
    assert list(_reread(cache)["table"].obs["marker"])[0] == "OLD"   # the prior version


def test_auto_repair_does_not_touch_data(tiny_sdata, make_table):
    """AUTO fixes bookkeeping only; restoring data needs an explicit choice."""
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))
    shutil.rmtree(cache / "tables" / "table")

    repair(cache, level=AUTO)
    assert not (cache / "tables" / "table").exists()


# ── stores too broken for zarr to open ───────────────────────────────────────

def test_verify_survives_a_truncated_root_metadata(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / "zarr.json").write_text("{")

    report = verify(cache)                 # must not raise
    assert not report.readable_metadata
    assert report.metadata_error
    assert not report.repairable
    assert sorted(report.on_disk) == ["labels/lab", "tables/table"]


def test_repair_rebuilds_a_truncated_root_metadata(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / "zarr.json").write_text("{")

    repair(cache)
    assert "table" in _reread(cache).tables


def test_verify_survives_a_missing_root_metadata(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / "zarr.json").unlink()
    report = verify(cache)
    assert not report.readable_metadata
    assert "no zarr.json" in report.metadata_error


# ── strays and debris ────────────────────────────────────────────────────────

def test_stray_groups_are_detected_and_dropped(tiny_sdata):
    """Generalises the hard-coded tables/adata_norm cleanup from app.py."""
    cache = Path(tiny_sdata.path)
    stray = cache / "tables" / "adata_norm"
    stray.mkdir()
    (stray / "zarr.json").write_text(json.dumps({"node_type": "group", "attributes": {}}))

    report = verify(cache)
    assert report.stray_elements == ["tables/adata_norm"]
    assert "tables/adata_norm" not in report.on_disk    # not counted as an element

    repair(cache, report)
    assert not stray.exists()


def test_debris_is_reported_and_cleared(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / STAGING_DIR).mkdir(parents=True, exist_ok=True)
    (cache / STAGING_DIR / "abandoned.zarr").mkdir()
    (cache / "tables" / "zarr.json.abc.partial").write_text("{}")

    report = verify(cache)
    assert len(report.debris) == 2
    assert not report.ok

    repair(cache, report)
    assert verify(cache).debris == []


def test_pending_writes_are_reported(tiny_sdata, make_table, monkeypatch):
    from palms.utils import zarr_safe

    cache = Path(tiny_sdata.path)
    real = os.rename
    calls = {"n": 0}

    def flaky(src, dst):
        calls["n"] += 1
        if calls["n"] == 2:
            raise KeyboardInterrupt()
        return real(src, dst)

    monkeypatch.setattr(zarr_safe.os, "rename", flaky)
    with pytest.raises(KeyboardInterrupt):
        safe_write_element(tiny_sdata, "table", make_table("NEW"))
    monkeypatch.undo()

    report = verify(cache)
    assert report.pending_ops and "table" in report.pending_ops[0]

    repair(cache, report)
    assert verify(cache).pending_ops == []
    assert list(_reread(cache)["table"].obs["marker"])[0] == "NEW"


def test_unreadable_journal_is_reported_not_fatal(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / JOURNAL_DIR).mkdir(parents=True, exist_ok=True)
    (cache / JOURNAL_DIR / "bad.json").write_text("{oops")

    report = verify(cache)
    assert any("unreadable journal" in op for op in report.pending_ops)


# ── incidental helpers ───────────────────────────────────────────────────────

def test_backups_are_found(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache.parent / "sdata_cached_corrupt_20260101_000000.zarr").mkdir()
    assert len(verify(cache).backups) == 1


def test_describe_store_excludes_internal_dirs(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))   # creates .xv_trash
    described = cache_repair.describe_store(cache)
    assert described["size_bytes"] > 0
    assert described["free_bytes"] > 0


def test_human_bytes():
    assert human_bytes(512) == "512 B"
    assert human_bytes(30 * 1024 ** 3) == "30.0 GB"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
