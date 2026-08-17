"""Write a SpatialData store one element at a time.

``SpatialData.write()`` is a single opaque call: it loops over
``gen_elements()`` internally, so there is no point at which a caller can
report progress, drop references, or measure anything. For a Xenium cache
build that loop runs for tens of minutes and moves tens of GB, and "RAM climbs
and never comes back" was not answerable without a per-element breakdown.

This module does the same work with the loop turned inside out, using only
public API:

* an empty :class:`SpatialData` writes the store shell (root group, attrs,
  consolidated metadata) — this is also what performs the
  "can I safely write here" checks;
* the real object is pointed at that store and each element is written with
  the public :meth:`SpatialData.write_element`;
* between elements we collect garbage and ``malloc_trim``, then log RSS.

Deliberately *not* replicating ``write()``'s up-front ``_validate_all_elements()``:
it only changes when a bad element is discovered, not whether, and every caller
here writes to a staging directory that is renamed into place on success — so a
partial store is discarded rather than seen. Avoiding it keeps this module off
spatialdata's private surface entirely.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Optional

import dask
from spatialdata import SpatialData

from xenium_viewer.utils import mem_probe

ProgressCb = Callable[[int, str], None]

# Concurrency cap for writing a points element.
#
# dask's default threaded scheduler runs one task per core — 40 on this box — and
# a transcript partition is ~150 MB uncompressed with all columns, so the default
# width can pull several GB in at once even though no single partition is ever
# fully collected. The crop path measured an ArrowMemoryError under a 24 GB cap
# at the default width, from ``PointsModel.parse``'s own eager monotonicity
# check alone. This trades write throughput for a hard bound.
WRITE_WORKERS = 2

# Labels need a cap too, for a different reason than points, and images do not.
#
# Both raster pyramids are built lazily, but by different code. Images (which
# have a ``c`` dim) go through ``multiscale_spatial_image``'s XARRAY_COARSEN;
# labels go through ``spatialdata.models.pyramids_utils.to_multiscale``, which
# does an explicit ``prev.astype(float)`` and **two rechunks per level** (before
# and after ``ome_zarr.dask_utils.resize``'s ``map_blocks``). A rechunk is a
# many-to-many gather, and the float64 intermediates are 134 MB a chunk, so 40 of
# them in flight is how the label write — not the far larger image write — was
# the one that died: ``Unable to allocate 177. MiB for an array with shape
# (3620, 6404) and data type float64``.
#
# Images are deliberately left at the default width. Measured on the real slide,
# capping them to 2 moved the peak only 8.71 GB -> 8.31 GB while making the
# raster write several times slower. What made the image write expensive was
# never its concurrency but its *chunking* — a source read as one chunk per full
# channel page had to materialise 5.93 GB before it could write anything, which
# is what ``utils/raster_io.py`` fixes.
LABEL_WORKERS = 4

ELEMENT_WORKERS = {"labels": LABEL_WORKERS, "points": WRITE_WORKERS}


def element_names(sdata: SpatialData) -> list[str]:
    """Element names in the order ``SpatialData.write()`` would write them."""
    return [name for _, name, _ in sdata.gen_elements()]


@contextmanager
def _scheduler_for(element_type: str):
    """Bound dask's concurrency for element types that need it (see above)."""
    workers = ELEMENT_WORKERS.get(element_type)
    if workers is None:
        yield
    else:
        with dask.config.set(scheduler="threads", num_workers=workers):
            yield


def write_sdata(
    sdata: SpatialData,
    path: str | Path,
    progress_cb: Optional[ProgressCb] = None,
    log: Optional[Callable[[str], None]] = None,
    start_pct: int = 0,
    end_pct: int = 100,
) -> None:
    """Write *sdata* to *path*, one element at a time.

    Equivalent to ``sdata.write(path)`` — same store, same element order, same
    ``sdata.path`` afterwards — but reports progress per element and releases
    memory between them.

    Parameters
    ----------
    progress_cb
        Called as ``(percent, message)`` before each element, with *percent*
        interpolated between *start_pct* and *end_pct*.
    log
        Where the per-element memory line goes. Defaults to the package logger.
    """
    path = Path(path)
    # Resolved before the first write: gen_elements() walks the live object, and
    # writing an element is not the moment to discover the sequence has changed.
    plan = [(etype, name) for etype, name, _ in sdata.gen_elements()]
    span = max(0, end_pct - start_pct)

    # An empty SpatialData writes the root group and attrs, and performs the
    # same "is this path safe to write to" validation that a full write() does.
    SpatialData(attrs=dict(sdata.attrs)).write(path)

    previous_path = sdata.path
    sdata.path = path
    try:
        for i, (element_type, name) in enumerate(plan):
            if progress_cb is not None:
                pct = start_pct + (span * i // max(1, len(plan)))
                progress_cb(pct, f"Writing {name} ({i + 1}/{len(plan)})...")
            with _scheduler_for(element_type):
                sdata.write_element(name)
            # Drop what the element's compute left behind *and* hand the pages
            # back to the OS: without malloc_trim, glibc keeps freed blocks in
            # its arenas and RSS only ever ratchets upwards, which is what made
            # this loop look like it leaked.
            mem_probe.release()
            mem_probe.log_memory(f"wrote {name}", log=log)
    except BaseException:
        sdata.path = previous_path
        raise

    sdata.write_consolidated_metadata()
    if progress_cb is not None:
        progress_cb(end_pct, "Wrote all elements.")
