"""``write_sdata`` must stay interchangeable with ``SpatialData.write()``.

The point of writing element by element is progress reporting and releasing
memory between elements — not a different store. So the property under test is
equivalence: same elements, same order, same bytes. If a spatialdata upgrade
changes what ``write()`` produces, these fail rather than the cache quietly
becoming something the viewer reads differently.
"""
from __future__ import annotations

import warnings

import pytest


@pytest.fixture
def sample_sdata(make_table):
    """Two rasters with pyramids, points and a table — one of each write path."""
    pytest.importorskip("spatialdata")
    np = pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    pd = pytest.importorskip("pandas")
    from spatialdata import SpatialData
    from spatialdata.models import Image2DModel, Labels2DModel, PointsModel
    from spatialdata.transformations import Identity

    def _build():
        rng = np.random.default_rng(0)
        img = da.from_array(
            rng.integers(0, 500, (2, 64, 48), dtype=np.uint16), chunks=(1, 16, 16))
        lab = da.from_array(
            rng.integers(0, 20, (64, 48), dtype=np.uint32), chunks=(16, 16))
        pts = pd.DataFrame({
            "x": rng.random(100) * 48,
            "y": rng.random(100) * 64,
            "feature_name": pd.Categorical(rng.choice(list("abcd"), 100)),
        })
        return SpatialData(
            images={"im": Image2DModel.parse(
                img, dims=("c", "y", "x"), transformations={"global": Identity()},
                scale_factors=[2], c_coords=["a", "b"])},
            labels={"lab": Labels2DModel.parse(
                lab, dims=("y", "x"), transformations={"global": Identity()},
                scale_factors=[2])},
            points={"pts": PointsModel.parse(
                pts, coordinates={"x": "x", "y": "y"}, feature_key="feature_name",
                transformations={"global": Identity()})},
            tables={"table": make_table("OLD")},
        )

    return _build


def _write_both(sample_sdata, tmp_path, **kwargs):
    from xenium_viewer.utils import sdata_write

    reference, produced = tmp_path / "ref.zarr", tmp_path / "new.zarr"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sample_sdata().write(reference)
        sdata_write.write_sdata(sample_sdata(), produced,
                                log=lambda line: None, **kwargs)
    return reference, produced


def test_matches_spatialdata_write(sample_sdata, tmp_path):
    np = pytest.importorskip("numpy")
    from spatialdata import read_zarr
    from xenium_viewer.utils import sdata_write

    reference, produced = _write_both(sample_sdata, tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref, new = read_zarr(str(reference)), read_zarr(str(produced))

    assert sdata_write.element_names(ref) == sdata_write.element_names(new)

    for element_type, name, element in ref.gen_elements():
        other = new[name]
        if element_type in ("images", "labels"):
            for level in element:
                a = element[level]["image"].data
                b = other[level]["image"].data
                assert (a.shape, a.dtype, a.chunks) == (b.shape, b.dtype, b.chunks), \
                    f"{name}/{level}"
                assert np.array_equal(a.compute(), b.compute()), f"{name}/{level}"
        elif element_type == "points":
            assert element.compute().equals(other.compute())
        else:
            assert np.array_equal(element.X, other.X)
            assert element.obs.equals(other.obs)


def test_reports_progress_for_every_element(sample_sdata, tmp_path):
    from xenium_viewer.utils import sdata_write

    seen = []
    _write_both(sample_sdata, tmp_path,
                progress_cb=lambda pct, msg: seen.append((pct, msg)))

    names = sdata_write.element_names(sample_sdata())
    # One call per element plus the final one; percentages never go backwards.
    assert len(seen) == len(names) + 1
    assert [pct for pct, _ in seen] == sorted(pct for pct, _ in seen)
    for name, (_, message) in zip(names, seen):
        assert name in message


def test_sets_path_and_consolidates_like_write(sample_sdata, tmp_path):
    from spatialdata import read_zarr
    from xenium_viewer.utils import sdata_write

    sdata = sample_sdata()
    assert sdata.path is None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sdata_write.write_sdata(sdata, tmp_path / "out.zarr", log=lambda line: None)
        assert sdata.path == tmp_path / "out.zarr"
        assert sdata.has_consolidated_metadata()
        read_zarr(str(tmp_path / "out.zarr"))   # opens without repair


def test_failure_restores_the_original_path(sample_sdata, tmp_path, monkeypatch):
    """A half-written store must not leave the object pointing at it.

    The caller writes to a staging directory and renames on success; if a
    failure left ``sdata.path`` on the abandoned staging path, the next
    ``write_element`` anywhere in the viewer would write into a directory that
    is about to be deleted.
    """
    from spatialdata import SpatialData
    from xenium_viewer.utils import sdata_write

    sdata = sample_sdata()
    sdata.path = tmp_path / "original.zarr"

    def boom(self, name, *args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(SpatialData, "write_element", boom)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(RuntimeError, match="disk full"):
            sdata_write.write_sdata(sdata, tmp_path / "staging.zarr",
                                    log=lambda line: None)
    assert sdata.path == tmp_path / "original.zarr"


def test_writing_points_sdata_at_the_destination(sample_sdata, tmp_path):
    """``sdata.path`` follows the write — which is what makes the loader's
    post-rename fixup necessary; see the source guard below."""
    from xenium_viewer.utils import sdata_write, zarr_safe

    staging = tmp_path / ".sdata_cached__building.zarr"
    sdata = sample_sdata()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sdata_write.write_sdata(sdata, staging, log=lambda line: None)

    assert sdata.path == staging
    assert zarr_safe.cache_path_of(sdata) == staging


def test_loader_repoints_sdata_after_renaming_staging_into_place():
    """A rebuild must not leave the object pointing at the staging directory.

    ``load_sdata`` writes to ``.sdata_cached__building.zarr`` and renames it into
    place, but the write sets ``sdata.path`` to the *staging* path, which by then
    no longer exists. Every later element write resolves against that attribute
    (``zarr_safe.cache_path_of``), so a stale value sends the first ROI or
    clustering saved after a rebuild into a resurrected staging directory instead
    of the cache — silently, and only noticed much later.

    A source guard rather than a behavioural test: reproducing it needs a real
    Xenium dataset, which CI does not have.
    """
    import ast
    import inspect
    import textwrap

    from xenium_viewer import loader

    tree = ast.parse(textwrap.dedent(inspect.getsource(loader.load_sdata)))
    renamed_at = repointed_at = None
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "rename"
                and any(getattr(a, "id", None) == "staging" for a in node.args)):
            renamed_at = node.lineno
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Attribute) and t.attr == "path"
                        and getattr(t.value, "id", None) == "sdata"
                        for t in node.targets)):
            repointed_at = node.lineno

    assert renamed_at is not None, "load_sdata no longer renames the staging dir"
    assert repointed_at is not None, (
        "load_sdata renames the staging directory into place but never resets "
        "sdata.path, so later element writes go to the staging path"
    )
    assert repointed_at > renamed_at, "sdata.path is reset before the rename"


def test_element_names_matches_gen_elements(sample_sdata):
    """The order elements are written in is spatialdata's, not ours."""
    from xenium_viewer.utils import sdata_write

    sdata = sample_sdata()
    assert sdata_write.element_names(sdata) == [
        name for _, name, _ in sdata.gen_elements()
    ]


def test_write_workers_is_shared_with_crop_export():
    """One bounded write path, not two that can drift apart."""
    from xenium_viewer.utils import crop_export, sdata_write

    assert crop_export._TRANSCRIPT_WRITE_WORKERS == sdata_write.WRITE_WORKERS


def test_images_keep_the_default_scheduler_width():
    """Images are the one raster type that does not need throttling.

    Measured on a real slide: capping them moved the peak 8.71 GB -> 8.31 GB and
    made the raster write several times slower. Labels are capped because their
    pyramid goes through float64 and two rechunks per level; points because a
    transcript partition is ~150 MB. If someone adds "images" here, they should
    have a measurement saying why.
    """
    from xenium_viewer.utils import sdata_write

    assert set(sdata_write.ELEMENT_WORKERS) == {"labels", "points"}


@pytest.mark.parametrize("element_type,expected", [
    ("points", 2), ("labels", 4), ("images", None), ("tables", None),
])
def test_scheduler_is_bounded_only_where_declared(element_type, expected):
    import dask
    from xenium_viewer.utils import sdata_write

    with sdata_write._scheduler_for(element_type):
        assert dask.config.get("num_workers", None) == expected
