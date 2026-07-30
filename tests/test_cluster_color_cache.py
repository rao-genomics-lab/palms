"""Regression: re-running a clustering under an existing key must not keep the
previous run's colors.

``CellColorManager.get_cluster_colors`` caches on the series' ``name``, which is
the clustering key. Producers that *replace* the series behind an existing key —
re-running Leiden at the same resolution, re-importing a file, a new CNV run —
got the old cached color array back, while the legend and cluster filter were
rebuilt from the new assignment. The result on screen was a clustering that
looked partially overwritten: cells whose new cluster id happened to match their
old one were right, the rest kept stale colors.

Run standalone:   python tests/test_cluster_color_cache.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")

from xenium_viewer.utils.coloring import CellColorManager  # noqa: E402


N_CELLS = 12
KEY = "leiden_r1.0"


def _manager():
    """A manager over 12 cells whose napari labels are 1..12 (0 = background)."""
    adata = anndata.AnnData(np.zeros((N_CELLS, 3), dtype="float32"))
    adata.obs["cell_id"] = [f"cell{i}" for i in range(N_CELLS)]
    label_to_obs = np.full(N_CELLS + 1, -1, dtype=np.int64)
    label_to_obs[1:] = np.arange(N_CELLS)
    return CellColorManager(adata, label_to_obs)


def _series(values, name=KEY):
    return pd.Series(values, index=[f"cell{i}" for i in range(N_CELLS)], name=name)


def _first_run():
    return _series([0] * 6 + [1] * 6)


def _second_run():
    """A genuinely different partition of the same cells, same key."""
    return _series([0] * 3 + [1] * 3 + [2] * 6)


def test_colors_are_cached_within_a_run():
    """The cache exists for a reason — don't regress it into a no-op."""
    mgr = _manager()
    first, _ = mgr.get_cluster_colors(_first_run())
    again, _ = mgr.get_cluster_colors(_first_run())
    assert again is first          # same object: served from the cache


def test_recomputing_the_same_key_returns_stale_colors_without_invalidation():
    """Pins the mechanism the bug rode on, so the fix cannot be silently undone.

    This is the *buggy* behaviour, reproduced deliberately: without an explicit
    invalidation the manager cannot tell that the series behind the key changed.
    """
    mgr = _manager()
    before, _ = mgr.get_cluster_colors(_first_run())
    after, _ = mgr.get_cluster_colors(_second_run())
    assert np.array_equal(before, after)          # stale — the point of the test


def test_invalidating_the_key_recomputes_the_colors():
    mgr = _manager()
    before, before_map = mgr.get_cluster_colors(_first_run())
    before = before.copy()

    mgr.invalidate_cluster_cache(KEY)
    after, after_map = mgr.get_cluster_colors(_second_run())

    assert not np.array_equal(before, after)
    assert len(before_map) == 2 and len(after_map) == 3
    # cells 3-5 moved from cluster 0 to cluster 1, so their color must change
    assert not np.array_equal(before[4], after[4])


def test_invalidating_without_a_name_clears_everything():
    mgr = _manager()
    mgr.get_cluster_colors(_first_run())
    mgr.get_cluster_colors(_series([0] * 12, name="other"))
    assert set(mgr._cluster_cache) == {KEY, "other"}

    mgr.invalidate_cluster_cache()
    assert mgr._cluster_cache == {}


def test_invalidating_an_absent_key_is_a_no_op():
    mgr = _manager()
    mgr.invalidate_cluster_cache("never-computed")   # must not raise


def test_different_keys_do_not_share_a_cache_entry():
    """Leiden named every series "leiden" regardless of resolution, so one entry
    served every run. The key is now the series name."""
    mgr = _manager()
    a, _ = mgr.get_cluster_colors(_series([0] * 6 + [1] * 6, name="leiden_r0.5"))
    b, _ = mgr.get_cluster_colors(_series([0] * 3 + [1] * 9, name="leiden_r1.0"))
    assert not np.array_equal(a, b)


def test_every_producer_that_replaces_a_clustering_invalidates_it():
    """A source-level guard: a tab that rebinds ``ctx.clusterings[key]`` and
    never invalidates is the bug returning under a different name."""
    import re

    tabs = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer" / "tabs"
    assignment = re.compile(r"ctx\.clusterings\[([A-Za-z_][\w]*)\]\s*=")
    offenders = []
    for path in sorted(tabs.glob("tab_*.py")):
        source = path.read_text()
        for var in set(assignment.findall(source)):
            if f"invalidate_cluster_cache({var})" not in source:
                offenders.append(f"{path.name}: ctx.clusterings[{var}] = ...")
    assert offenders == [], (
        "these rebind a clustering without invalidating its cached colors: "
        + "; ".join(offenders)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
