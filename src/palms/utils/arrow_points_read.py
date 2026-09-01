"""Read the cache's points elements through pyarrow's filesystem, not fsspec.

A compatibility patch, in the idiom of the others listed under "Known
Compatibility Patches" in CLAUDE.md — and the only one that reaches into a
*private* module path, so it is written to fail quietly and to announce its own
obsolescence.

Why it exists
-------------
``spatialdata._io.io_points._read_points`` calls ``dask.dataframe.read_parquet``
with no ``filesystem=``, so the points element comes back as a
``ReadParquetFSSpec``. In dask, predicate pushdown into the parquet reader is
``_filter_passthrough``, and that is ``True`` on exactly one class:

===========================  =====================
class                        ``_filter_passthrough``
===========================  =====================
``ReadParquet`` (base)       ``False``
``ReadParquetPyarrowFS``     ``True``
``ReadParquetFSSpec``        ``False``
===========================  =====================

So a boolean mask over ``sdata.points['transcripts']`` can never become a
parquet filter, and the transcript-density step reads all 128.7 M rows of a
real Xenium cache to answer a question about one gene. Passing
``filesystem="arrow"`` yields the class that does push the filter down:
measured 5.44 s → 1.32 s. See ``docs/transcript-read-pushdown.md`` for the
measurement, and why sorting the cache (the obvious alternative) does nothing.

What it does not do
-------------------
It is a context manager rather than a process-wide patch, applied only around
the cache read in ``loader._open_cache``. The pyarrow reader returns identical
data — same rows, dtypes, partition count and transformations — but emits
partitions in a different **order** and keeps ``__null_dask_index__`` as the
index name. Nothing in PALMS depends on points row order (``crop_export`` does
its own filtered read and then ``reset_index(drop=True)``; the density step
histograms; the overlay uses the feather index), but a process-wide patch would
impose that ordering on every points element anywhere in the process, including
one a caller loaded through stock spatialdata for its own reasons.

Removing it
-----------
``tests/test_arrow_points_read.py`` fails as soon as stock spatialdata passes a
``filesystem`` of its own — that is the signal to delete this module, its call
site, and the ``== True`` comment in ``transcripts.gene.tmpl``. The wrapper also
declines to act when a ``filesystem`` is already given, so an upstream fix takes
effect even before anyone gets round to the deletion.
"""

from __future__ import annotations

import contextlib
import logging

logger = logging.getLogger(__name__)


def _is_local(path) -> bool:
    """Only local stores get the pyarrow filesystem.

    ``_read_points`` prefixes remote stores with ``simplecache::``, an fsspec
    chained URL that the pyarrow filesystem cannot parse — reading one that way
    would turn a working remote dataset into a crash.
    """
    text = str(path)
    return not text.startswith("simplecache::") and "://" not in text


@contextlib.contextmanager
def arrow_points_reader():
    """Make points reads inside this block use pyarrow's filesystem.

    Yields ``True`` if the patch went on, ``False`` if the private module path
    it needs has moved — in which case the reads are stock and simply slower,
    which is the right way for a private-API patch to fail.
    """
    try:
        import spatialdata._io.io_points as io_points

        stock = io_points.read_parquet
    except Exception as exc:  # pragma: no cover - depends on spatialdata layout
        logger.info(
            "points reads stay on the fsspec reader: %s. Transcript queries will "
            "be slower; see palms/utils/arrow_points_read.py", exc,
        )
        yield False
        return

    if not callable(stock):  # pragma: no cover - defensive
        yield False
        return

    def _read(path, *args, **kwargs):
        # A caller that names a filesystem outranks us — including a future
        # spatialdata that passes one itself, which is what makes this wrapper
        # inert rather than harmful once upstream ships the fix.
        if kwargs.get("filesystem") is None and _is_local(path):
            try:
                return stock(path, *args, filesystem="arrow", **kwargs)
            except Exception as exc:
                logger.warning(
                    "pyarrow filesystem read failed (%s); falling back to fsspec",
                    exc,
                )
        return stock(path, *args, **kwargs)

    io_points.read_parquet = _read
    try:
        yield True
    finally:
        io_points.read_parquet = stock
