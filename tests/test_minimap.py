"""The overview thumbnail is read for the minimap, and for nothing else.

``_populate_viewer`` used to ``.compute()`` the bottom morphology level on the
main thread at every load, under a comment saying it was "for coarse tissue
alignment". That stopped being true when ``he.coarse_align.tmpl`` shipped: the
fit derives its own thumbnail from ``sdata``, so the recorded cell reads the
element rather than a launch-time value, and the H&E tab was left reading the
array only to decide whether to enable a button.

What actually reads those pixels is the minimap, one channel of them, into a
200x160 pixmap. Measured on the cached pancreas section (``morphology_focus/s5``,
4x430x1067 uint16) the load-time read cost 0.66 s of the main thread cold and
0.03 s warm -- small, but paid on every launch, larger over a network share, and
unbounded when the pyramid is a coarsen chain rather than stored bytes (the
first-load crash: 23.1 GB, killed session, see ``level_is_computed``).

Two properties are pinned here:

* ``raster_io.overview_thumbnail`` returns the DAPI *plane* of a stored level,
  and ``None`` for a level that is a computation.
* no tab reads ``morph_thumb`` -- the H&E tab must ask ``sdata`` whether the
  dataset has morphology, which is the same question the template asks, and is
  answerable before the background read has finished.

Run with:  pytest tests/test_minimap.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

da = pytest.importorskip("dask.array")

from palms.utils.raster_io import overview_thumbnail  # noqa: E402

SRC = Path(__file__).resolve().parent.parent / "src" / "palms"


def _stored(shape, chunks):
    """A level read from disk: one task per chunk."""
    data = np.arange(int(np.prod(shape)), dtype=np.uint16).reshape(shape)
    return da.from_array(data, chunks=chunks)


# ── what the minimap gets ────────────────────────────────────────────────────

def test_overview_thumbnail_takes_the_dapi_plane():
    """Channel 0 of (C, Y, X), as 2-D — not the whole stack."""
    level = _stored((4, 16, 24), chunks=(1, 8, 8))
    thumb = overview_thumbnail(level)

    assert thumb is not None
    assert thumb.shape == (16, 24), "the minimap draws one plane, not four"
    assert np.array_equal(thumb, np.asarray(level[0]))


def test_overview_thumbnail_passes_a_two_dimensional_level_through():
    """A single-channel element is already the plane."""
    level = _stored((16, 24), chunks=(8, 8))
    thumb = overview_thumbnail(level)

    assert thumb is not None
    assert thumb.shape == (16, 24)


def test_overview_thumbnail_refuses_a_computed_level():
    """A coarsen chain is not read to fill a 200x160 pixmap.

    This is the ``--no-cache`` case, and the one that used to cost 23.1 GB and
    the session. Refusing costs a minimap; reading costs the run.
    """
    stored = _stored((4, 32, 32), chunks=(1, 16, 16))
    chained = da.coarsen(np.mean, stored, {1: 2, 2: 2})

    assert overview_thumbnail(chained) is None
    assert overview_thumbnail(None) is None


def test_overview_thumbnail_refuses_a_shape_it_does_not_understand():
    """Anything but (C, Y, X) or (Y, X) is refused, not indexed into on a guess."""
    assert overview_thumbnail(_stored((2, 4, 8, 8), chunks=(1, 1, 4, 4))) is None


# ── the dead coupling must not come back ─────────────────────────────────────

def test_no_tab_reads_the_morphology_thumbnail():
    """A tab asks ``sdata`` what the dataset has, not what the viewer computed.

    ``tab_he_registration`` gated Coarse Align on ``ctx.morph_thumb is not
    None`` — a display array, for a step that reads the element itself. With the
    read moved into the background that gate would also have been *wrong* for as
    long as the read took.
    """
    offenders = [
        path.relative_to(SRC).as_posix()
        for path in sorted((SRC / "tabs").rglob("*.py"))
        if "morph_thumb" in path.read_text()
    ]
    assert not offenders, (
        "these tabs read the minimap's display array; ask sdata instead: "
        + ", ".join(offenders)
    )
