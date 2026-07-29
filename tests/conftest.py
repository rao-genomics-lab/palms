"""Shared fixtures.

``tiny_sdata`` builds a real on-disk SpatialData zarr store small enough to
write in milliseconds, so the persistence paths — which had no coverage at all —
can be exercised in CI rather than only by hand on a 30 GB cache.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


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
