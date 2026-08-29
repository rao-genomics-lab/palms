"""Read the morphology OME-TIFF as the tiled image it already is.

``spatialdata_io.xenium()`` reads ``morphology_focus`` through
``dask_image.imread.imread``, which yields **one dask chunk per full channel
page**. Measured on a 57887x51217 4-channel Xenium slide::

    dask_image.imread(...)  ->  chunksize (1, 57887, 51217) = 5.93 GB per chunk

``Image2DModel.parse(chunks=(1, 4096, 4096))`` then rechunks that into ~195
tiles per channel, but dask cannot produce *any* tile without decoding the whole
5.93 GB page, and must hold it until every consumer is done — including the
chained pyramid levels above it. Four channels under the default scheduler means
all four pages resident at once: **a ~24 GB floor before a byte is written**, and
it is the single largest term in the cache build's memory use.

The file itself has none of that problem — it is tiled 1024x1024::

    tifffile.imread(ch0000_dapi.ome.tif, aszarr=True)
      level 0  (4, 57887, 51217) uint16  chunks (1, 1024, 1024)
      level 1  (4, 28943, 25608) uint16  chunks (1, 1024, 1024)
      ...                                       8 levels

so :func:`tiled_morphology_image` reads level 0 through that route instead. The
OME metadata stitches all four ``ch000N_*.ome.tif`` files into one 4-channel
series, which is why one filename gives the whole image.

**Only level 0 is taken, and the pyramid is still computed by spatialdata.** The
file's own levels 1+ were the obvious thing to reuse — they exist, they are
tiled, and their shapes match ``scale_factors=[2, 2, 2, 2, 2]`` exactly — but
measured against the real slide they are *not* what spatialdata produces:

    10x level1 vs coarsen-mean   mean|diff| = 8.63   (level mean value 57.6)
    10x level1 vs decimation     mean|diff| = 5.98
    10x level1 vs max-pool       mean|diff| = 14.74

10x uses some other filter, and none of those is a rounding difference — it is a
~15% change to every pixel the viewer draws when zoomed out. Reading level 0 and
recomputing the rest removes the same 24 GB while leaving every written byte
identical to what the old path produced, so that is what this does.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

# Chunks for the parsed element and the written store. 1024 -> 4096 is a clean
# 4x4 gather (one output chunk needs 16 source tiles, ~33 MB), so the on-disk
# layout is unchanged from what spatialdata_io produced.
DEFAULT_CHUNKS = (1, 4096, 4096)


def level_is_computed(level) -> bool:
    """True if *level* is built from a coarsen chain rather than read from disk.

    A stored level needs one task per chunk (plus a bookkeeping key); a level
    standing on a chain of ``coarsen().mean()`` carries the tasks of every level
    below it, so touching it materialises the largest one. That is the whole of
    the first-load crash ``loader._reopen_written_cache`` exists to prevent, and
    the same test applies to any pyramid from any source — an H&E TIFF with no
    internal levels included — so it lives here rather than inline in ``app.py``.
    """
    dask_graph = getattr(level, "dask", None)
    if dask_graph is None:
        return False
    return len(dask_graph) > level.npartitions + 1


def find_morphology_tiff(data_path: Path) -> Optional[Path]:
    """The OME-TIFF holding morphology_focus, across Xenium output layouts."""
    focus_dir = data_path / "morphology_focus"
    if focus_dir.is_dir():
        # Any one file resolves the whole multi-file OME series through its
        # metadata; sorted() only makes the choice deterministic.
        tiffs = sorted(
            p for p in focus_dir.iterdir()
            if p.name.endswith(".ome.tif") and not p.name.startswith("._")
        )
        if tiffs:
            return tiffs[0]
    flat = data_path / "morphology_focus.ome.tif"
    return flat if flat.is_file() else None


def open_ome_tiff_pyramid(tif_path: Path) -> list:
    """All resolution levels of *tif_path* as dask arrays, chunked as on disk.

    Raises whatever tifffile/zarr raise; callers wanting a fallback should use
    :func:`tiled_morphology_image`.
    """
    import dask.array as da
    import tifffile
    import zarr

    opened = zarr.open(tifffile.imread(str(tif_path), aszarr=True), mode="r")
    if not hasattr(opened, "array_keys"):
        # A TIFF with no sub-resolutions opens as a bare array rather than a
        # group. Still worth reading this way — the tiling is what matters, and
        # the older single-file Xenium layouts land here.
        return [da.from_zarr(opened)]
    # Level keys are "0", "1", ... — sort numerically, or level 10 would sort
    # between 1 and 2.
    keys = sorted(opened.array_keys(), key=int)
    return [da.from_zarr(opened[k]) for k in keys]


def tiled_morphology_image(data_path: Path, reference):
    """Full-resolution morphology_focus as a tile-chunked dask array.

    Returns ``None`` — meaning "keep what spatialdata_io built" — if there is no
    readable OME-TIFF, or if its full-resolution level does not match *reference*
    in both shape and dtype. Declining is always safe: the caller keeps the
    existing element and the build behaves exactly as it did before.

    Parameters
    ----------
    reference
        The element ``spatialdata_io`` produced. Used to confirm we are looking
        at the same image before swapping anything.
    """
    tif_path = find_morphology_tiff(Path(data_path))
    if tif_path is None:
        return None

    try:
        levels = open_ome_tiff_pyramid(tif_path)
    except Exception:
        return None
    if not levels:
        return None

    scale0 = reference_scale0(reference)
    if scale0 is None:
        return None
    if levels[0].shape != tuple(scale0.shape) or levels[0].dtype != scale0.dtype:
        return None
    # Channel names come from the reference, so no names means no swap: parsing
    # without them would silently relabel the channels 0, 1, 2, 3.
    if reference_channels(reference) is None:
        return None
    return levels[0]


def pyramid_scale_factors(reference) -> Optional[list]:
    """The ``scale_factors`` that would rebuild *reference*'s level structure.

    Derived from the reference rather than from a constant, so the re-parsed
    element has the same levels as the one it replaces no matter who configured
    them — ``spatialdata_io`` fills in its own default ``scale_factors`` whenever
    the caller does not pass one, so a caller-side constant is not the authority.

    Returns ``[]`` for a single-scale reference, and ``None`` if the levels are
    not a clean chain of halvings (which ``scale_factors`` cannot express).
    """
    if not hasattr(reference, "children"):
        return []
    names = sorted(reference.children, key=lambda k: int(k.removeprefix("scale")))
    shapes = [reference[name]["image"].shape for name in names]
    factors = []
    for previous, current in zip(shapes, shapes[1:]):
        if current[:-2] != previous[:-2]:
            return None
        if current[-2:] != (previous[-2] // 2, previous[-1] // 2):
            return None
        factors.append(2)
    return factors


def reference_scale0(reference):
    """Full-resolution DataArray of *reference*, whether DataTree or DataArray."""
    try:
        if hasattr(reference, "children") and "scale0" in reference:
            return reference["scale0"]["image"]
        return reference
    except Exception:
        return None


def reference_channels(reference) -> Optional[list]:
    """Channel names carried by *reference*, or None if it has none."""
    scale0 = reference_scale0(reference)
    try:
        return list(scale0.coords["c"].values)
    except Exception:
        return None


def reference_transformations(reference) -> dict:
    """Transformations of *reference*, defaulting to a global identity."""
    from spatialdata.transformations import Identity
    from spatialdata.transformations._utils import _get_transformations_xarray

    scale0 = reference_scale0(reference)
    try:
        transformations = _get_transformations_xarray(scale0)
    except Exception:
        transformations = None
    return dict(transformations) if transformations else {"global": Identity()}
