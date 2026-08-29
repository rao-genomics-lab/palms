"""Sidecar analysis outputs belong beside the zarr store, not inside it.

h5ad caches, parquet DEG tables and the CopyKAT worker's JSON are not zarr
nodes. Written into the store root, zarr's hierarchy walk hits each one, fails
to open it, and emits a ZarrUserWarning — the most likely source of the "several
warnings" in the bug report. Living inside the store also meant a cache rebuild
destroyed them, including multi-hour CopyKAT runs.

Readers must still find the old location, because existing datasets have files
there and nothing is migrated eagerly.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_sidecar_location.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")

from palms.utils.adata_persistence import (  # noqa: E402
    ARMS_DEG_CACHE, ROI_DEG_CACHE, SIDECAR_DIRNAME, find_sidecar, glob_sidecars,
    load_arms_tile_deg_from_sdata, load_roi_deg_from_sdata,
    save_arms_tile_deg_to_sdata, save_roi_deg_to_sdata, sidecar_dir,
    sidecar_write_path,
)


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _ctx(sdata):
    return SimpleNamespace(
        sdata=sdata, adata=sdata["table"], no_cache=False,
        data_path=Path(sdata.path).parent, segmentation_source="xenium",
    )


def _df():
    return pd.DataFrame({"group": ["a", "b"], "score": [1.0, 2.0]})


# ── where things get written ─────────────────────────────────────────────────

def test_sidecars_are_written_beside_the_store_not_inside_it(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)

    save_roi_deg_to_sdata(ctx, _df())

    assert (cache.parent / SIDECAR_DIRNAME / ROI_DEG_CACHE).exists()
    assert not (cache / ROI_DEG_CACHE).exists()


def test_arms_deg_too(tiny_sdata):
    cache = Path(tiny_sdata.path)
    save_arms_tile_deg_to_sdata(_ctx(tiny_sdata), _df())
    assert (cache.parent / SIDECAR_DIRNAME / ARMS_DEG_CACHE).exists()
    assert not (cache / ARMS_DEG_CACHE).exists()


def test_sidecar_write_path_accepts_a_bare_sdata(tiny_sdata):
    """Not every caller has a ViewerContext."""
    path = sidecar_write_path(tiny_sdata, "x.h5ad")
    assert path.parent.name == SIDECAR_DIRNAME
    assert path.parent.parent == Path(tiny_sdata.path).parent


# ── reading back ─────────────────────────────────────────────────────────────

def test_roundtrip_through_the_new_location(tiny_sdata):
    ctx = _ctx(tiny_sdata)
    save_roi_deg_to_sdata(ctx, _df())
    pd.testing.assert_frame_equal(load_roi_deg_from_sdata(tiny_sdata), _df())


def test_legacy_files_in_the_store_root_are_still_found(tiny_sdata):
    """Existing datasets have their sidecars inside the store."""
    cache = Path(tiny_sdata.path)
    _df().to_parquet(cache / ROI_DEG_CACHE, index=False)

    assert find_sidecar(tiny_sdata, ROI_DEG_CACHE) == cache / ROI_DEG_CACHE
    pd.testing.assert_frame_equal(load_roi_deg_from_sdata(tiny_sdata), _df())


def test_the_new_location_wins_over_the_legacy_one(tiny_sdata):
    cache = Path(tiny_sdata.path)
    _df().to_parquet(cache / ROI_DEG_CACHE, index=False)
    save_roi_deg_to_sdata(_ctx(tiny_sdata), _df().assign(score=[9.0, 9.0]))

    assert list(load_roi_deg_from_sdata(tiny_sdata)["score"]) == [9.0, 9.0]


def test_missing_sidecar_reads_as_none(tiny_sdata):
    assert find_sidecar(tiny_sdata, "nope.parquet") is None
    assert load_roi_deg_from_sdata(tiny_sdata) is None
    assert load_arms_tile_deg_from_sdata(tiny_sdata) is None


def test_glob_finds_both_locations_with_the_new_one_winning(tiny_sdata):
    cache = Path(tiny_sdata.path)
    (cache / "cnv_infercnv_result.json").write_text("{}")
    new_home = sidecar_dir(cache.parent, create=True)
    (new_home / "cnv_copykat_result.json").write_text("{}")
    (new_home / "cnv_infercnv_result.json").write_text("{}")

    found = glob_sidecars(tiny_sdata, "cnv_*_result.json")
    assert [p.name for p in found] == ["cnv_copykat_result.json", "cnv_infercnv_result.json"]
    assert all(p.parent.name == SIDECAR_DIRNAME for p in found)


# ── why this matters ─────────────────────────────────────────────────────────

def test_new_sidecars_do_not_make_consolidation_warn(tiny_sdata):
    """The store root is where the ZarrUserWarnings came from."""
    from palms.utils.zarr_safe import consolidate

    ctx = _ctx(tiny_sdata)
    save_roi_deg_to_sdata(ctx, _df())
    sidecar_write_path(ctx, "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        consolidate(Path(tiny_sdata.path))


def test_a_cache_rebuild_cannot_destroy_them(tiny_sdata):
    """They now live outside the directory a rebuild replaces."""
    import shutil

    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    save_roi_deg_to_sdata(ctx, _df())

    shutil.rmtree(cache)                       # what a rebuild does to the store
    assert (cache.parent / SIDECAR_DIRNAME / ROI_DEG_CACHE).exists()


def test_legacy_sidecars_still_count_as_user_data_at_stake(tmp_path):
    """A rebuild must still warn about CNV results left in the old location."""
    from palms.loader import _detect_user_data, _has_any_user_data

    cache = tmp_path / "sdata_cached.zarr"
    cache.mkdir()
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    assert _has_any_user_data(_detect_user_data(cache))


def test_new_location_sidecars_are_reported_too(tmp_path):
    from palms.loader import _detect_user_data

    cache = tmp_path / "sdata_cached.zarr"
    cache.mkdir()
    home = tmp_path / SIDECAR_DIRNAME
    home.mkdir()
    (home / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")

    assert "adata_cnv_cache_copykat.h5ad" in _detect_user_data(cache)["sidecars"]


def test_no_writer_targets_the_store_root_any_more():
    """Source guard: sidecar writes must go through the sidecar helpers."""
    import re

    src = Path(__file__).resolve().parent.parent / "src" / "palms"
    pattern = re.compile(
        r"Path\((?:ctx\.)?sdata\.path\)\s*/\s*[\"'f]?[^\"')]*"
        r"(?:\.h5ad|\.parquet|_result\.json)")
    offenders = []
    for path in sorted(src.rglob("*.py")):
        text = path.read_text()
        for match in pattern.finditer(text):
            line = text[: match.start()].count("\n") + 1
            offenders.append(f"{path.relative_to(src)}:{line}")
    assert offenders == [], (
        "sidecar files must not be written into the zarr store root — use "
        "sidecar_write_path(): " + ", ".join(offenders)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
