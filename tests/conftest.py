"""Shared fixtures.

``tiny_sdata`` builds a real on-disk SpatialData zarr store small enough to
write in milliseconds, so the persistence paths — which had no coverage at all —
can be exercised in CI rather than only by hand on a 30 GB cache.
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# Qt picks its platform plugin when QApplication is constructed, and on a
# headless box the default (xcb) does not fail — it *aborts the process*. CI
# reported nothing but "Aborted (core dumped)" with no test name, and a bare
# `pytest` over ssh did the same. Choosing offscreen up front turns that into a
# working run; an explicit QT_QPA_PLATFORM or a real display still wins, so this
# does not change what happens on a desktop.
if not os.environ.get("QT_QPA_PLATFORM") and not os.environ.get("DISPLAY"):
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
# Same shape of problem: matplotlib's default backend needs a display.
os.environ.setdefault("MPLBACKEND", "Agg")

# Analysis templates can be overridden from ~/.config/xenium-viewer/templates/.
# An empty search path disables that for the whole suite: a developer who has
# customised a template must not get different test results from CI, and the
# tests that pin template text are asserting what the package *ships*, not what
# this machine happens to run. Set deliberately rather than with setdefault —
# an inherited value would reintroduce exactly the divergence this prevents.
#
# NOTE: this also means the code path taken when the variable is *unset* — which
# is every real user — runs nowhere in this suite by default. That gap shipped a
# viewer that could not start: `search_path()` and `user_template_dir()`
# delegated to each other and recursed forever, and no test noticed because none
# of them ever had the variable unset. `tests/test_template_overrides.py` now has
# a `no_env` fixture that deletes it (redirecting the platform config dir), and
# anything reachable at launch should be covered there too.
os.environ["XENIUM_VIEWER_TEMPLATE_PATH"] = ""


@pytest.fixture(scope="session")
def qapp():
    """A QApplication, so widget-building code can be exercised headless.

    Run the suite with QT_QPA_PLATFORM=offscreen; the instance is shared because
    Qt allows only one per process.
    """
    pytest.importorskip("qtpy")
    from qtpy.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def make_table():
    """Factory for a valid SpatialData table carrying a marker value."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    from spatialdata.models import TableModel

    def _make(marker: str = "OLD", n_obs: int = 6):
        adata = anndata.AnnData(np.ones((n_obs, 3), dtype="float32"))
        adata.obs["region"] = pd.Categorical(["lab"] * n_obs)
        adata.obs["instance_id"] = list(range(n_obs))
        adata.obs["marker"] = [marker] * n_obs
        return TableModel.parse(
            adata, region="lab", region_key="region", instance_key="instance_id",
        )

    return _make


@pytest.fixture
def tiny_sdata(tmp_path, make_table):
    """A written-to-disk SpatialData store: one labels element and one table."""
    pytest.importorskip("spatialdata")
    np = pytest.importorskip("numpy")
    from spatialdata import SpatialData, read_zarr
    from spatialdata.models import Labels2DModel

    cache = tmp_path / "sdata_cached.zarr"
    labels = Labels2DModel.parse(np.arange(16, dtype=np.int32).reshape(4, 4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SpatialData(labels={"lab": labels}, tables={"table": make_table("OLD")}).write(cache)
        return read_zarr(cache)


@pytest.fixture(scope="session")
def replay_adata():
    """A small AnnData with structure to find: two populations plus coordinates.

    Used by the notebook-replay test, which runs the real analysis steps over it
    and then re-runs the exported notebook. Counts are drawn from a seeded RNG so
    the *input* is identical in-process and on replay — any divergence in the
    output is then attributable to the recorded code.
    """
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")

    def _make(n_obs: int = 200, n_vars: int = 60):
        rng = np.random.default_rng(0)
        counts = rng.poisson(3, size=(n_obs, n_vars)).astype("float32")
        half = n_obs // 2
        counts[:half, : n_vars // 3] += 12          # population A markers
        counts[half:, n_vars // 3: 2 * n_vars // 3] += 12   # population B markers
        adata = anndata.AnnData(counts)
        adata.obs_names = [f"cell{i}" for i in range(n_obs)]
        adata.var_names = [f"gene{i}" for i in range(n_vars)]
        # Two spatially separated blobs, so the neighbour graph is meaningful.
        centers = np.where(np.arange(n_obs) < half, 0.0, 400.0)
        adata.obsm["spatial"] = np.column_stack([
            centers + rng.normal(0, 40, n_obs),
            rng.normal(200, 40, n_obs),
        ])
        return adata

    return _make


@pytest.fixture
def marker_of():
    """Read the marker column back out of a store or SpatialData."""
    def _marker(sdata_or_path):
        import warnings as _w
        from spatialdata import read_zarr
        if isinstance(sdata_or_path, (str, Path)):
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                sdata_or_path = read_zarr(str(sdata_or_path))
        return list(sdata_or_path["table"].obs["marker"])[0]
    return _marker
