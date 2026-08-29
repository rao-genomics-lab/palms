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
    from palms.utils import sdata_write

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
    from palms.utils import sdata_write

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
    from palms.utils import sdata_write

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
    from palms.utils import sdata_write

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
    from palms.utils import sdata_write

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
    from palms.utils import sdata_write, zarr_safe

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

    from palms import loader

    tree = ast.parse(textwrap.dedent(inspect.getsource(loader.load_sdata)))
    renamed_at = repointed_at = None
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "rename"
                and any(getattr(a, "id", None) == "staging" for a in node.args)):
            renamed_at = node.lineno
        # Either form settles the path: assigning it outright, or rebinding
        # `sdata` to the re-opened cache, which carries its own path.
        if isinstance(node, ast.Assign):
            if any(isinstance(t, ast.Attribute) and t.attr == "path"
                   and getattr(t.value, "id", None) == "sdata"
                   for t in node.targets):
                repointed_at = node.lineno
            elif (any(getattr(t, "id", None) == "sdata" for t in node.targets)
                  and isinstance(node.value, ast.Call)
                  and getattr(node.value.func, "id", None)
                  == "_reopen_written_cache"):
                repointed_at = node.lineno

    assert renamed_at is not None, "load_sdata no longer renames the staging dir"
    assert repointed_at is not None, (
        "load_sdata renames the staging directory into place but never resets "
        "sdata.path, so later element writes go to the staging path"
    )
    assert repointed_at > renamed_at, "sdata.path is reset before the rename"

    # The re-open helper is the one carrying that guarantee now, so it has to
    # keep the fallback that sets the path when the re-open fails.
    helper = ast.parse(
        textwrap.dedent(inspect.getsource(loader._reopen_written_cache)))
    assert any(
        isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Attribute) and t.attr == "path"
                and getattr(t.value, "id", None) == "sdata"
                for t in node.targets)
        for node in ast.walk(helper)
    ), ("_reopen_written_cache no longer sets sdata.path on the fallback path, "
        "so a failed re-open leaves the object pointing at the staging dir")


def _levels(element):
    """The dask array behind every scale of a multiscale element."""
    return [element[name]["image"].data
            for name in sorted(element.children, key=lambda n: int(n[5:]))]


def test_reopen_returns_arrays_that_are_read_rather_than_recomputed(
        sample_sdata, tmp_path):
    """The whole point of the swap: no level may still be a coarsen chain.

    A freshly parsed multiscale element holds only ``scale0`` as data; every
    smaller level is a lazy ``coarsen().mean()`` over the one below it. napari
    displays the *smallest* level first, so on a full slide adding the layer
    walked the entire pyramid — tens of GB — for a thumbnail. After the re-open
    each level is an array on disk, which shows up as one task per chunk.
    """
    pytest.importorskip("spatialdata")
    from palms import loader
    from palms.utils import sdata_write

    cache = tmp_path / "sdata_cached.zarr"
    built = sample_sdata()

    # The premise. If this stops holding, the fix is no longer needed and this
    # test is the one that should say so.
    assert any(len(arr.dask) > arr.npartitions
               for arr in _levels(built.images["im"])), (
        "expected the parsed pyramid to be a lazy chain, not stored levels")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sdata_write.write_sdata(built, cache, log=lambda line: None)
        reopened = loader._reopen_written_cache(built, cache)

    # A stored array is one task per chunk, plus the single bookkeeping key
    # `da.from_zarr` adds. A computed one carries the whole chain below it.
    for name, element in (("im", reopened.images["im"]),
                          ("lab", reopened.labels["lab"])):
        for i, arr in enumerate(_levels(element)):
            assert len(arr.dask) <= arr.npartitions + 1, (
                f"{name} scale{i} is still computed ({len(arr.dask)} tasks for "
                f"{arr.npartitions} chunks) — the build graph was returned")

    from pathlib import Path
    assert Path(reopened.path) == cache


def test_reopen_falls_back_to_the_in_memory_dataset(sample_sdata, tmp_path,
                                                    monkeypatch):
    """A cache that will not re-open must not lose the build."""
    pytest.importorskip("spatialdata")
    from pathlib import Path

    from palms import loader

    def _boom(_path):
        raise RuntimeError("nope")

    monkeypatch.setattr(loader, "_open_cache", _boom)

    built = sample_sdata()
    cache = tmp_path / "sdata_cached.zarr"
    returned = loader._reopen_written_cache(built, cache)

    assert returned is built
    assert Path(returned.path) == cache


def test_element_names_matches_gen_elements(sample_sdata):
    """The order elements are written in is spatialdata's, not ours."""
    from palms.utils import sdata_write

    sdata = sample_sdata()
    assert sdata_write.element_names(sdata) == [
        name for _, name, _ in sdata.gen_elements()
    ]


def test_write_workers_is_shared_with_crop_export():
    """One bounded write path, not two that can drift apart."""
    from palms.utils import crop_export, sdata_write

    assert crop_export._TRANSCRIPT_WRITE_WORKERS == sdata_write.WRITE_WORKERS


def test_images_keep_the_default_scheduler_width():
    """Images are the one raster type that does not need throttling.

    Measured on a real slide: capping them moved the peak 8.71 GB -> 8.31 GB and
    made the raster write several times slower. Labels are capped because their
    pyramid goes through float64 and two rechunks per level; points because a
    transcript partition is ~150 MB. If someone adds "images" here, they should
    have a measurement saying why.
    """
    from palms.utils import sdata_write

    assert set(sdata_write.ELEMENT_WORKERS) == {"labels", "points"}


@pytest.mark.parametrize("element_type,expected", [
    ("points", 2), ("labels", 4), ("images", None), ("tables", None),
])
def test_scheduler_is_bounded_only_where_declared(element_type, expected):
    import dask
    from palms.utils import sdata_write

    with sdata_write._scheduler_for(element_type):
        assert dask.config.get("num_workers", None) == expected
