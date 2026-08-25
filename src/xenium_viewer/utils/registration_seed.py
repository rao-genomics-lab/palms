"""Where an overlay's placement comes from when a session is missing or partial.

The viewer places a registered overlay from *session state* — `he_state["affine_3x3"]`
composed with a flip derived from `flip_v`/`flip_h` and the image's shape. The element's
own spatialdata transformation is written alongside it but never read back at restore
time. That works right up until a store has no `viewer_session` to read: a crop export, a
cache built by `xenium-build-cache`, a session node deleted in Tools -> Dataset, or a
recovered cache. Then the overlay is added at the origin, unregistered, with no error.

So the element has to be able to place itself. `seed_registration` is the one place that
decides, and it exists as a single function because the decision cannot be half-applied.

Why it is all-or-nothing
------------------------
The two records do not store the same thing. `_save_he_affine_to_sdata` writes
``fine @ flip`` onto the element — the flip is *baked in* — while the session stores
``fine`` alone and re-derives the flip from `flip_v`/`flip_h` plus `he_shape_yx`. Mixing
them applies the flip twice: a session that says `flip_v=True` combined with an element
affine that already contains that flip yields ``flip @ fine @ flip``, which is not a
misplacement anyone would notice as a bug — the image simply sits somewhere wrong.

Hence: the element is used only when the session offers *no* affine at all, and when it
is used the flips are forced False and the shape is taken from the element. The caller
gets back a complete, self-consistent set or the session's own values, never a blend.

What the caller should know
---------------------------
When the element is the source, the flip checkboxes read unticked. That is correct rather
than a compromise — the stored raster already has the flip inside its transform, so
"flip" from that point on means "flip this raster *again*".
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class Registration(NamedTuple):
    """A complete placement for one overlay. Never a mix of two sources."""

    fine: np.ndarray | None
    coarse: np.ndarray | None
    flip_v: bool
    flip_h: bool
    shape_yx: tuple | None
    #: Which record this came from — ``"session"`` or ``"element"``. Callers use it
    #: for the status line; nothing branches on it.
    source: str


def seed_registration(session, sdata, element_name: str,
                      affine_key: str, coarse_key: str | None,
                      flip_v_key: str, flip_h_key: str, shape_key: str,
                      element_shape_yx=None) -> Registration:
    """Resolve one overlay's placement: session first, element second, never both.

    ``element_shape_yx`` is the *stored* raster's ``(y, x)`` — needed because an
    exported crop's overlay has a different shape from the one the session recorded,
    and a flip built against the wrong shape is off by the difference.
    """
    from xenium_viewer.utils.adata_persistence import _load_affine_from_sdata_element

    session = session or {}
    fine = session.get(affine_key)
    coarse = session.get(coarse_key) if coarse_key else None

    if fine is not None or coarse is not None:
        shape = session.get(shape_key) or element_shape_yx
        return Registration(
            fine=fine, coarse=coarse,
            flip_v=bool(session.get(flip_v_key, False)),
            flip_h=bool(session.get(flip_h_key, False)),
            shape_yx=tuple(shape) if shape else None,
            source="session",
        )

    stored = _load_affine_from_sdata_element(sdata, element_name)
    if stored is None:
        # Nothing anywhere. Keep whatever the session said about flips and shape;
        # with no affine the layer is placed by the flip alone, which is what the
        # viewer has always done for an unregistered overlay.
        shape = session.get(shape_key) or element_shape_yx
        return Registration(
            fine=None, coarse=None,
            flip_v=bool(session.get(flip_v_key, False)),
            flip_h=bool(session.get(flip_h_key, False)),
            shape_yx=tuple(shape) if shape else None,
            source="session",
        )

    return Registration(
        fine=stored, coarse=None, flip_v=False, flip_h=False,
        shape_yx=tuple(element_shape_yx) if element_shape_yx else None,
        source="element",
    )
