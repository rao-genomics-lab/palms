"""When the loader is allowed to throw a cache away.

The reported complaint was that it did so too easily. Three separate paths did:

* an unreadable cache was renamed aside — or ``rmtree``'d if the rename failed —
  with no repair attempt and no check for user data;
* staleness compared ``experiment.xenium``'s mtime against the cache *directory*
  mtime, so ``rsync``/``cp -p``/a re-download condemned a perfectly good cache;
* the sidecar list omitted the CNV caches, so a cache whose only user data was a
  multi-hour CopyKAT run reported "no user data" and was rebuilt with no dialog.

These are pure-function tests; none of them builds a 30 GB anything.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_loader_policy.py
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")

from xenium_viewer import loader  # noqa: E402
from xenium_viewer.loader import (  # noqa: E402
    _detect_user_data, _format_user_data_message, _has_any_user_data,
    _is_cache_stale, _restore_user_elements, write_manifest,
)


@pytest.fixture
def fake_cache(tmp_path):
    """A directory shaped like a cache, without the cost of building one."""
    cache = tmp_path / "sdata_cached.zarr"
    (cache / "tables" / "table" / "obs").mkdir(parents=True)
    (cache / "tables" / "table" / "uns").mkdir(parents=True)
    (cache / "tables" / "table" / "obsm").mkdir(parents=True)
    experiment = tmp_path / "experiment.xenium"
    experiment.write_text('{"pixel_size": 0.2125}')
    return cache, experiment


# ── the silent-rebuild case ──────────────────────────────────────────────────

def test_copykat_results_count_as_user_data(fake_cache):
    """The worst case in the report: hours of compute, rebuilt with no dialog."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    (cache / "cnv_copykat_result.json").write_text("{}")

    user_data = _detect_user_data(cache)
    assert _has_any_user_data(user_data)
    assert "adata_cnv_cache_copykat.h5ad" in user_data["sidecars"]

    message = _format_user_data_message(user_data)
    assert "CopyKAT" in message and "hours of compute" in message


def test_infercnv_results_count_as_user_data(fake_cache):
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_infercnv.h5ad").write_bytes(b"x")
    assert _has_any_user_data(_detect_user_data(cache))
    assert "inferCNV" in _format_user_data_message(_detect_user_data(cache))


def test_an_untouched_cache_has_no_user_data(fake_cache):
    cache, _ = fake_cache
    assert not _has_any_user_data(_detect_user_data(cache))


def test_clusterings_and_rois_still_count(fake_cache):
    cache, _ = fake_cache
    (cache / "tables" / "table" / "obs" / "clustering_leiden_r1.0").mkdir()
    (cache / "shapes" / "rois").mkdir(parents=True)

    user_data = _detect_user_data(cache)
    assert _has_any_user_data(user_data)
    assert user_data["clusterings"] == ["clustering_leiden_r1.0"]
    assert "ROIs" in _format_user_data_message(user_data)


# ── staleness ────────────────────────────────────────────────────────────────

def test_without_a_manifest_staleness_is_uncertain(fake_cache):
    cache, experiment = fake_cache
    os.utime(experiment, (time.time() + 100, time.time() + 100))

    stale, certain = _is_cache_stale(cache, experiment)
    assert stale and not certain      # uncertain ⇒ the caller must ask


def test_a_manifest_makes_a_touched_source_a_non_event(fake_cache):
    """Copying a dataset bumps the mtime without changing a byte."""
    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    os.utime(experiment, (time.time() + 100, time.time() + 100))

    stale, certain = _is_cache_stale(cache, experiment)
    assert not stale and certain


def test_a_manifest_still_detects_a_real_change(fake_cache):
    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    experiment.write_text('{"pixel_size": 0.4250}')

    stale, certain = _is_cache_stale(cache, experiment)
    assert stale and certain


def test_a_missing_source_is_never_stale(fake_cache):
    """A Crop Dataset export has no experiment.xenium to compare against."""
    cache, experiment = fake_cache
    experiment.unlink()
    assert _is_cache_stale(cache, experiment) == (False, True)


def test_manifest_records_the_versions_it_was_built_with(fake_cache):
    from xenium_viewer.utils.cache_repair import read_manifest

    cache, experiment = fake_cache
    write_manifest(cache, experiment)
    manifest = read_manifest(cache)
    assert manifest["source"] == "experiment.xenium"
    assert len(manifest["source_sha256"]) == 64
    assert "spatialdata_version" in manifest


def test_the_manifest_is_invisible_to_readers(tiny_sdata):
    """It lives in the store root, so it must not look like an element."""
    import warnings

    from spatialdata import read_zarr

    cache = Path(tiny_sdata.path)
    experiment = cache.parent / "experiment.xenium"
    experiment.write_text("{}")
    write_manifest(cache, experiment)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sorted(read_zarr(str(cache)).tables) == ["table"]


# ── restoring user data on rebuild ───────────────────────────────────────────

def _adata_pair():
    import anndata
    import pandas as pd

    old = anndata.AnnData(np.ones((4, 2), dtype="float32"))
    old.obs_names = [f"c{i}" for i in range(4)]
    old.obs["clustering_leiden_r1.0"] = pd.Categorical(["0", "0", "1", "1"])
    old.obs["cnv_score_copykat"] = [0.1, 0.2, 0.3, 0.4]
    old.obs["copykat_leiden_res0.2"] = pd.Categorical(["a", "a", "b", "b"])
    old.obs["total_counts"] = [1, 2, 3, 4]          # rebuilt anyway
    old.uns["rank_genes_groups"] = {"names": ["g1"]}
    old.uns["cnv_runs"] = {"copykat": {"resolution": 0.2}}
    old.uns["rank_genes_groupby"] = "clustering_leiden_r1.0"
    old.obsm["X_umap"] = np.zeros((4, 2))
    old.obsm["X_cnv_umap"] = np.ones((4, 2))

    new = anndata.AnnData(np.ones((4, 2), dtype="float32"))
    new.obs_names = [f"c{i}" for i in range(4)]
    new.obs["total_counts"] = [9, 9, 9, 9]
    return old, new


class _FakeSdata:
    def __init__(self, table):
        self.tables = {"table": table}
        self.shapes: dict = {}
        self.images: dict = {}

    def __getitem__(self, key):
        return self.tables[key]

    def __setitem__(self, key, value):
        self.tables[key] = value


def test_restore_covers_cnv_and_rank_genes_keys():
    """Regression: these were dropped even on 'Rebuild and restore my data'."""
    old, new = _adata_pair()
    _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                           {"shapes": [], "images": [], "uns_keys": [],
                            "has_obsm_umap": True})

    assert "clustering_leiden_r1.0" in new.obs
    assert "cnv_score_copykat" in new.obs
    assert "copykat_leiden_res0.2" in new.obs
    assert new.uns["cnv_runs"] == {"copykat": {"resolution": 0.2}}
    assert new.uns["rank_genes_groupby"] == "clustering_leiden_r1.0"
    assert "X_umap" in new.obsm and "X_cnv_umap" in new.obsm


def test_restore_does_not_clobber_freshly_built_columns():
    """A rebuilt total_counts must win over the stale one."""
    old, new = _adata_pair()
    _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                           {"shapes": [], "images": [], "uns_keys": [],
                            "has_obsm_umap": True})
    assert list(new.obs["total_counts"]) == [9, 9, 9, 9]


def test_restore_reports_what_it_moved():
    old, new = _adata_pair()
    restored = _restore_user_elements(_FakeSdata(old), _FakeSdata(new),
                                      {"shapes": [], "images": [], "uns_keys": [],
                                       "has_obsm_umap": True})
    assert "clustering_leiden_r1.0" in restored
    assert "UMAP coordinates" in restored


# ── never discard without asking ─────────────────────────────────────────────

def test_stale_with_user_data_and_no_dialog_keeps_the_cache(monkeypatch, fake_cache):
    """The old default here was 'rebuild' — a silent 30 GB discard."""
    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    assert loader._ask_rebuild_preference(_detect_user_data(cache)) == "keep"


def test_an_unopenable_cache_with_no_dialog_raises_rather_than_rebuilding(
    monkeypatch, fake_cache,
):
    from xenium_viewer.utils import cache_repair

    cache, _ = fake_cache
    (cache / "adata_cnv_cache_copykat.h5ad").write_bytes(b"x")
    monkeypatch.setattr(loader, "_qt_message_box", lambda: None)

    with pytest.raises(loader.CacheLoadAborted, match="Refusing to rebuild"):
        loader._ask_corrupt_cache(
            OSError("bad store"), cache_repair.verify(cache), _detect_user_data(cache),
        )


def test_the_loader_never_deletes_a_cache_directory():
    """Every destructive branch must be a rename. Source-level guard."""
    import re

    source = (Path(__file__).resolve().parent.parent
              / "src" / "xenium_viewer" / "loader.py").read_text()
    offenders = [
        line.strip() for line in source.splitlines()
        if re.search(r"rmtree\(\s*(str\()?cache_path", line)
    ]
    assert offenders == [], (
        "loader.py must never rmtree the live cache — move it aside instead: "
        + "; ".join(offenders)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
