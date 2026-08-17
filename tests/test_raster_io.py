"""Reading morphology_focus from the OME-TIFF's tiles instead of whole pages.

The swap is a pure read-path change: same array values, same parse call, same
written bytes — only the dask chunking of the source differs. So the properties
worth testing are that the data really is identical, and that every disagreement
makes the reader **decline**, leaving the old path in place.

The declines matter more than they look. They are the branch that never runs on
a developer's machine, and the one deciding whether an older Xenium layout still
loads at all.
"""
from __future__ import annotations

import pytest

from xenium_viewer.utils import raster_io


@pytest.fixture
def write_pyramid_tiff(tmp_path):
    """Write a tiled, multi-resolution OME-TIFF; return its path and level 0."""
    tifffile = pytest.importorskip("tifffile")
    np = pytest.importorskip("numpy")

    def _write(shape=(2, 256, 192), levels=6, dtype="uint16", tile=32,
               name="morphology_focus.ome.tif"):
        rng = np.random.default_rng(0)
        base = rng.integers(0, 400, shape, dtype=dtype)
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        opts = {"tile": (tile, tile), "photometric": "minisblack"}
        with tifffile.TiffWriter(path, ome=True, bigtiff=True) as tif:
            tif.write(base, subifds=levels - 1, **opts)
            current = base
            for _ in range(levels - 1):
                current = current[:, ::2, ::2]
                tif.write(current, subfiletype=1, **opts)
        return path, base

    return _write


def _reference(array, c_coords=("a", "b"), chunks=None):
    """A stand-in for the element spatialdata_io hands us.

    ``chunks=array.shape`` on the source is the point: that is what
    ``dask_image.imread`` produces — one chunk per full channel page.
    """
    da = pytest.importorskip("dask.array")
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import Identity

    return Image2DModel.parse(
        da.from_array(array, chunks=array.shape), dims=("c", "y", "x"),
        transformations={"global": Identity()}, c_coords=list(c_coords),
        chunks=chunks or raster_io.DEFAULT_CHUNKS, scale_factors=[2, 2],
    )


def test_reads_every_level_tiled(write_pyramid_tiff):
    np = pytest.importorskip("numpy")
    path, base = write_pyramid_tiff()

    levels = raster_io.open_ome_tiff_pyramid(path)

    assert len(levels) == 6
    assert np.array_equal(levels[0].compute(), base)
    # The whole point: chunks come from the TIFF's tiles, not one per page.
    assert levels[0].chunksize[-2:] == (32, 32)


def test_level_keys_sort_numerically(write_pyramid_tiff):
    """A lexical sort would place level 10 between 1 and 2."""
    path, _ = write_pyramid_tiff(shape=(1, 2048, 2048), levels=11, tile=16)

    levels = raster_io.open_ome_tiff_pyramid(path)

    shapes = [level.shape[-1] for level in levels]
    assert shapes == sorted(shapes, reverse=True)


def test_tiled_image_is_identical_to_the_reference(write_pyramid_tiff):
    np = pytest.importorskip("numpy")
    path, base = write_pyramid_tiff()
    reference = _reference(base)

    level0 = raster_io.tiled_morphology_image(path.parent, reference)

    assert level0 is not None
    assert np.array_equal(level0.compute(), base)
    assert np.array_equal(
        level0.compute(), reference["scale0"]["image"].data.compute())
    # ...but chunked to the tiles rather than to one chunk per page.
    assert level0.chunksize < reference["scale0"]["image"].data.chunksize


def test_parsed_pyramid_is_unchanged_by_the_swap(write_pyramid_tiff):
    """The reason only level 0 is taken: every written byte must stay the same."""
    np = pytest.importorskip("numpy")
    from spatialdata.models import Image2DModel

    path, base = write_pyramid_tiff()
    reference = _reference(base, chunks=(1, 64, 64))
    level0 = raster_io.tiled_morphology_image(path.parent, reference)

    reparsed = Image2DModel.parse(
        level0, dims=("c", "y", "x"),
        transformations=raster_io.reference_transformations(reference),
        chunks=(1, 64, 64),
        scale_factors=raster_io.pyramid_scale_factors(reference),
        c_coords=raster_io.reference_channels(reference), rgb=None,
    )

    assert sorted(reparsed.children) == sorted(reference.children)
    for level in reference:
        old = reference[level]["image"].data
        new = reparsed[level]["image"].data
        assert (old.shape, old.dtype, old.chunks) == (new.shape, new.dtype, new.chunks)
        assert np.array_equal(old.compute(), new.compute()), level


def test_carries_over_channels_and_transformations(write_pyramid_tiff):
    from spatialdata.transformations import Scale, get_transformation

    path, base = write_pyramid_tiff()
    reference = _reference(base, c_coords=("DAPI", "18S"))
    from spatialdata.transformations import set_transformation
    set_transformation(reference, Scale([2.0, 2.0], axes=("x", "y")), "global")

    assert raster_io.reference_channels(reference) == ["DAPI", "18S"]
    carried = raster_io.reference_transformations(reference)
    assert carried.keys() == get_transformation(reference, get_all=True).keys()
    assert isinstance(carried["global"], Scale)


@pytest.mark.parametrize("scale_factors,expected", [
    ([2, 2], [2, 2]), ([2], [2]), (None, []),
])
def test_scale_factors_are_derived_from_the_reference(scale_factors, expected):
    """The replacement must rebuild whatever levels the original actually had.

    ``spatialdata_io`` fills in its own default ``scale_factors`` when the caller
    passes none, so a caller-side constant is not the authority on how many
    levels exist.
    """
    da = pytest.importorskip("dask.array")
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import Identity

    reference = Image2DModel.parse(
        da.zeros((2, 256, 192), dtype="uint16", chunks=(1, 64, 64)),
        dims=("c", "y", "x"), transformations={"global": Identity()},
        c_coords=["a", "b"], scale_factors=scale_factors,
    )
    assert raster_io.pyramid_scale_factors(reference) == expected


def test_scale_factors_declines_on_a_non_halving_pyramid():
    """A chain ``scale_factors`` cannot express must not be guessed at."""
    da = pytest.importorskip("dask.array")
    from spatialdata.models.pyramids_utils import dask_arrays_to_datatree

    odd = dask_arrays_to_datatree(
        [da.zeros((2, 300, 300), dtype="uint16"),
         da.zeros((2, 100, 100), dtype="uint16")],
        dims=("c", "y", "x"), channels=["a", "b"],
    )
    assert raster_io.pyramid_scale_factors(odd) is None


def test_scale_factors_sorts_levels_numerically():
    """scale10 must not sort between scale1 and scale2."""
    da = pytest.importorskip("dask.array")
    from spatialdata.models.pyramids_utils import dask_arrays_to_datatree

    deep = dask_arrays_to_datatree(
        [da.zeros((1, 4096 >> i, 4096 >> i), dtype="uint16") for i in range(11)],
        dims=("c", "y", "x"), channels=["a"],
    )
    assert raster_io.pyramid_scale_factors(deep) == [2] * 10


# ── Declining: each of these must leave spatialdata_io's element alone ───────

def test_declines_when_there_is_no_tiff(tmp_path):
    np = pytest.importorskip("numpy")
    reference = _reference(np.zeros((2, 64, 48), dtype="uint16"))
    assert raster_io.tiled_morphology_image(tmp_path, reference) is None


def test_declines_when_the_file_is_not_readable(tmp_path):
    np = pytest.importorskip("numpy")
    (tmp_path / "morphology_focus.ome.tif").write_bytes(b"not a tiff")
    reference = _reference(np.zeros((2, 64, 48), dtype="uint16"))
    assert raster_io.tiled_morphology_image(tmp_path, reference) is None


def test_reads_a_tiff_with_no_sub_resolutions(write_pyramid_tiff):
    """Older single-level layouts open as a bare array, not a group.

    The tiling is the whole benefit, so these must still be read this way rather
    than falling back to a whole-page decode.
    """
    np = pytest.importorskip("numpy")
    path, base = write_pyramid_tiff(levels=1)

    levels = raster_io.open_ome_tiff_pyramid(path)
    assert len(levels) == 1
    assert levels[0].chunksize[-2:] == (32, 32)
    assert np.array_equal(levels[0].compute(), base)

    assert raster_io.tiled_morphology_image(
        path.parent, _reference(base)) is not None


def test_channels_are_carried_through_even_when_unnamed(write_pyramid_tiff):
    """A reference with no channel *names* still has positional coords.

    Those get carried through unchanged, which is the point: the swap must
    reproduce the reference's channel identity, not invent one.
    """
    path, base = write_pyramid_tiff()
    nameless = _reference(base)["scale0"]["image"].drop_vars("c")

    assert raster_io.reference_channels(nameless) == [0, 1]
    assert raster_io.tiled_morphology_image(path.parent, nameless) is not None


def test_declines_when_the_reference_is_not_an_array():
    """The guard that stops a nonsense reference producing relabelled channels."""
    assert raster_io.reference_channels(object()) is None


def test_declines_when_the_reference_disagrees(write_pyramid_tiff):
    """Shape or dtype mismatch means we are not looking at the same image."""
    np = pytest.importorskip("numpy")
    path, base = write_pyramid_tiff()

    assert raster_io.tiled_morphology_image(
        path.parent, _reference(np.zeros((2, 128, 96), dtype="uint16"))) is None
    assert raster_io.tiled_morphology_image(
        path.parent, _reference(base.astype("uint8"))) is None


def test_finds_the_tiff_in_the_morphology_focus_directory(write_pyramid_tiff):
    path, base = write_pyramid_tiff(name="morphology_focus/ch0000_dapi.ome.tif")

    assert raster_io.find_morphology_tiff(path.parent.parent) == path
    assert raster_io.tiled_morphology_image(
        path.parent.parent, _reference(base)) is not None


def test_ignores_apple_double_files(write_pyramid_tiff):
    """``._``-prefixed siblings appear on any NAS a Mac has touched."""
    path, _ = write_pyramid_tiff(name="morphology_focus/ch0000_dapi.ome.tif")
    (path.parent / "._ch0000_dapi.ome.tif").write_bytes(b"resource fork")

    assert raster_io.find_morphology_tiff(path.parent.parent) == path


# ── The loader-side swap ─────────────────────────────────────────────────────

def test_loader_swap_never_raises(monkeypatch, tmp_path):
    """A broken read must degrade to the old path, not take the load down."""
    import xenium_viewer.loader as loader

    class FakeSdata:
        images = {"morphology_focus": object()}

    def boom(*args, **kwargs):
        raise OSError("nope")

    monkeypatch.setattr(raster_io, "tiled_morphology_image", boom)
    assert loader._retile_morphology_focus(FakeSdata(), tmp_path) is False


def test_loader_swap_skips_datasets_without_the_image(tmp_path):
    import xenium_viewer.loader as loader

    class FakeSdata:
        images: dict = {}

    assert loader._retile_morphology_focus(FakeSdata(), tmp_path) is False


def test_loader_swap_replaces_the_element(write_pyramid_tiff):
    np = pytest.importorskip("numpy")
    import xenium_viewer.loader as loader

    path, base = write_pyramid_tiff(shape=(2, 256, 192), levels=6)

    class FakeSdata:
        def __init__(self):
            self.images = {"morphology_focus": _reference(base)}

    sdata = FakeSdata()
    before = sdata.images["morphology_focus"]

    assert loader._retile_morphology_focus(sdata, path.parent) is True
    after = sdata.images["morphology_focus"]

    assert after is not before
    assert np.array_equal(after["scale0"]["image"].data.compute(), base)
    # Same level structure as the element it replaced, derived from it.
    assert sorted(after.children) == sorted(before.children)
    for level in before:
        old, new = before[level]["image"].data, after[level]["image"].data
        assert (old.shape, old.dtype, old.chunks) == (new.shape, new.dtype, new.chunks)
        assert np.array_equal(old.compute(), new.compute()), level
