"""The Cache tab: the repair actions a user can reach without a terminal.

The cache used to be a black box — when it broke, the loader moved it aside and
rebuilt from raw, and the user's only signal was a long wait. These tests drive
the tab's callbacks against a real on-disk store, headless.

They deliberately exercise the *callbacks*, not the widgets: the wiring between
a QPushButton and a function is not where bugs live, but "does re-consolidate
actually fix a broken store" is.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_tab_cache.py
"""
from __future__ import annotations

import os
import shutil
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
pytest.importorskip("qtpy")
np = pytest.importorskip("numpy")

from xenium_viewer.tabs import tab_cache  # noqa: E402
from xenium_viewer.utils import cache_repair  # noqa: E402
from xenium_viewer.utils.zarr_safe import consolidate, safe_write_element  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _ctx(sdata=None, data_path=None, no_cache=False):
    return SimpleNamespace(
        sdata=sdata, viewer=None, state={}, no_cache=no_cache,
        data_path=data_path or (Path(sdata.path).parent if sdata else None),
    )


def _reread(cache):
    from spatialdata import read_zarr
    return read_zarr(str(cache))


def _break_metadata(cache: Path):
    """Drop the table from the consolidated index without touching the data."""
    stash = cache.parent / "stash"
    os.rename(cache / "tables" / "table", stash)
    consolidate(cache)
    os.rename(stash, cache / "tables" / "table")


# ── locating the cache ───────────────────────────────────────────────────────

def test_cache_path_from_sdata(tiny_sdata):
    assert tab_cache._cache_path(_ctx(tiny_sdata)) == Path(tiny_sdata.path)


def test_cache_path_from_data_path(tiny_sdata):
    """A dataset whose sdata failed to load still has a cache to repair."""
    data_path = Path(tiny_sdata.path).parent
    assert tab_cache._cache_path(_ctx(None, data_path)) == Path(tiny_sdata.path)


def test_cache_path_is_none_when_there_is_no_cache(tmp_path):
    assert tab_cache._cache_path(_ctx(None, tmp_path)) is None


# ── building the tab ─────────────────────────────────────────────────────────

def test_build_tab_returns_a_widget_and_exports(tiny_sdata, qapp):
    widget, exports = tab_cache.build_tab(_ctx(tiny_sdata))
    assert widget is not None
    assert callable(exports["restore_session"])


def test_build_tab_survives_no_cache_mode(tmp_path, qapp):
    """--no-cache means nothing is persisted; the tab must still build."""
    widget, _ = tab_cache.build_tab(_ctx(None, tmp_path, no_cache=True))
    assert widget is not None


# ── the actions, driven directly ─────────────────────────────────────────────

def test_verify_reports_a_healthy_cache(tiny_sdata):
    report = cache_repair.verify(Path(tiny_sdata.path))
    assert report.ok
    assert "healthy" in report.summary()


def test_reconsolidate_fixes_a_store_that_will_not_open(tiny_sdata):
    """What the "Re-consolidate Metadata" button does, end to end."""
    cache = Path(tiny_sdata.path)
    _break_metadata(cache)
    assert "table" not in _reread(cache).tables          # broken

    report = cache_repair.verify(cache)
    result = cache_repair.repair(cache, report, level=cache_repair.AUTO)

    assert result.changed and not result.failures
    assert "table" in _reread(cache).tables
    assert cache_repair.verify(cache).ok


def test_recover_lists_a_previous_version_after_a_write(tiny_sdata, make_table):
    """The backup combo is populated from .xv_trash."""
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))

    report = cache_repair.verify(cache)
    assert "tables/table" in report.trash_available


def test_recover_from_trash_puts_the_previous_version_back(tiny_sdata, make_table):
    cache = Path(tiny_sdata.path)
    safe_write_element(tiny_sdata, "table", make_table("NEW"))
    assert list(_reread(cache)["table"].obs["marker"])[0] == "NEW"

    backup = cache_repair.verify(cache).trash_available["tables/table"][0]
    cache_repair.restore_from_trash(cache, "tables/table", backup)

    assert list(_reread(cache)["table"].obs["marker"])[0] == "OLD"


def test_recovering_a_whole_cache_fills_only_the_gaps(tiny_sdata, tmp_path):
    """Elements already in the live cache are newer — leave them alone."""
    import geopandas as gpd
    from shapely.geometry import Polygon
    from spatialdata.models import ShapesModel

    cache = Path(tiny_sdata.path)
    # A "previous cache" holding an ROI the live one lacks.
    backup = cache.parent / "sdata_cached_prev_20260101_000000.zarr"
    shutil.copytree(cache, backup)
    gdf = ShapesModel.parse(gpd.GeoDataFrame(geometry=[Polygon([(0, 0), (4, 0), (4, 4)])]))
    backup_sdata = _reread(backup)
    safe_write_element(backup_sdata, "rois", gdf)

    ctx = _ctx(tiny_sdata)
    actions = tab_cache._recover_whole_cache(ctx, cache, backup)

    assert "rois" in actions
    assert "rois" in _reread(cache).shapes


def test_recovering_when_nothing_is_missing_says_so(tiny_sdata):
    cache = Path(tiny_sdata.path)
    backup = cache.parent / "sdata_cached_prev_20260101_000000.zarr"
    shutil.copytree(cache, backup)

    actions = tab_cache._recover_whole_cache(_ctx(tiny_sdata), cache, backup)
    assert actions == ["nothing to recover — the live cache already has it all"]


def test_recovering_a_cache_pulls_back_sidecars(tiny_sdata, tmp_path):
    """CopyKAT results are the most expensive thing a rebuild can lose."""
    from xenium_viewer.utils.adata_persistence import sidecar_dir

    cache = Path(tiny_sdata.path)
    backup = cache.parent / "sdata_cached_corrupt_20260101_000000.zarr"
    shutil.copytree(cache, backup)
    (backup / "adata_cnv_cache_copykat.h5ad").write_bytes(b"hours of compute")

    actions = tab_cache._recover_whole_cache(_ctx(tiny_sdata), cache, backup)

    assert "adata_cnv_cache_copykat.h5ad" in actions
    restored = sidecar_dir(Path(tiny_sdata.path).parent) / "adata_cnv_cache_copykat.h5ad"
    assert restored.read_bytes() == b"hours of compute"


def test_recovery_uses_the_safe_write_path(tiny_sdata):
    """Recovery must not be able to corrupt the cache it is repairing."""
    source = Path(tab_cache.__file__).read_text()
    assert "safe_write_element" in source
    assert "write_element(" not in source.replace("safe_write_element(", "")


# ── the guard the whole feature exists for ───────────────────────────────────

def test_the_tab_never_deletes_the_cache():
    """Force Rebuild moves the cache aside; nothing here may remove one."""
    import re

    source = Path(tab_cache.__file__).read_text()
    offenders = [
        line.strip() for line in source.splitlines()
        if re.search(r"rmtree\s*\(", line)
    ]
    assert offenders == [], (
        "the Cache tab must never delete a cache — move it aside: "
        + "; ".join(offenders)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
