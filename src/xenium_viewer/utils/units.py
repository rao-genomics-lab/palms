"""Physical units for the napari canvas — one boundary, stated once.

The viewer's scale bar read in pixels because napari measures the *world*, and
every layer was added with the default scale of 1, so one world unit was one
Xenium morphology pixel. Giving each layer ``scale = (pixel_size, pixel_size)``
and ``units = ("um", "um")`` makes one world unit a micrometre and the scale bar
correct, in µm or mm as the zoom warrants.

That change is not free, and the cost is concentrated in one place.

Why the affines have to be converted
------------------------------------
napari composes a layer's transform as ``world = affine(scale(data))`` — the
affine is applied **after** the scale, so its translation is in *world* units.
Verified directly: a layer whose affine translates by ``(100, 50)`` puts
``data(0, 0)`` at world ``(100, 50)`` both before and after ``scale`` is set.

Every affine this codebase stores is in **Xenium pixels**: registration fits one
from landmark pixel coordinates, and ``adata_persistence`` writes it to the zarr
element in the same frame — which the crop export then composes with, and which
is the dataset's coordinate contract on disk. Switch the world to micrometres
without touching those and every registered overlay silently shifts by a factor
of ``1 / pixel_size`` — roughly 4.7×, with nothing raised and nothing logged.

So the rule is: **stored affines stay in pixels; napari layer affines are in
world units; convert at that boundary and nowhere else.** :func:`px_affine_to_world`
and :func:`world_affine_to_px` are that boundary.

The conversion is a similarity conjugation, ``A' = S A S⁻¹`` with
``S = diag(pixel_size, pixel_size)``, which for an affine reduces to scaling the
translation column and leaving the linear part alone — a rotation is the same
rotation whatever the unit, but an offset of 100 pixels is an offset of 21.25 µm.
It is written as the conjugation anyway, because that is the reason it is right.

What does *not* need converting
-------------------------------
Copying an affine from one layer to another (``utils/affine_linking.py``,
``img_lm.affine = lyr.affine``) is unit-agnostic: both layers carry the same
``scale``, so a world affine stays a world affine. Those sites are correct
untouched, and adding a conversion to them would be a double-scaling bug.

Layer **data** stays in Xenium pixels throughout. Nothing that reads
``layer.data`` — the crop export, the ROI tab, the ARMS tile ingest — is affected
by any of this.
"""

from __future__ import annotations

import numpy as np

#: The unit string napari's pint-backed ``layer.units`` understands.
#:
#: Note that a *scaled* unit does not work: ``layer.units = "0.2125 um"`` is
#: accepted and then silently discards the magnitude, leaving ``Unit('micrometer')``
#: — one pixel labelled as one micrometre. The magnitude has to live in ``scale``.
MICRON = "um"


def scale_for(pixel_size: float) -> tuple:
    """The ``layer.scale`` that puts a layer's pixel data into micrometres."""
    p = float(pixel_size)
    if not np.isfinite(p) or p <= 0:
        raise ValueError(f"pixel_size must be a positive finite number, got {pixel_size!r}")
    return (p, p)


def _conjugate(affine, pixel_size: float, forward: bool) -> np.ndarray:
    m = np.asarray(affine, dtype=np.float64)
    if m.shape != (3, 3):
        raise ValueError(f"expected a 3x3 (y, x) affine, got shape {m.shape}")
    p = float(pixel_size)
    if not np.isfinite(p) or p <= 0:
        raise ValueError(f"pixel_size must be a positive finite number, got {pixel_size!r}")

    s = np.diag([p, p, 1.0])
    s_inv = np.diag([1.0 / p, 1.0 / p, 1.0])
    return s @ m @ s_inv if forward else s_inv @ m @ s


def px_affine_to_world(affine_px, pixel_size: float) -> np.ndarray:
    """A stored pixel-space affine, expressed for ``layer.affine``."""
    return _conjugate(affine_px, pixel_size, forward=True)


def world_affine_to_px(affine_world, pixel_size: float) -> np.ndarray:
    """A ``layer.affine`` read back, expressed for storage and for geometry."""
    return _conjugate(affine_world, pixel_size, forward=False)


def layer_affine_px(layer, pixel_size: float) -> np.ndarray:
    """The pixel-space affine of *layer*, or identity if it has none.

    The one function to use when persisting a layer's registration or applying
    it to pixel coordinates, so that no call site has to remember which frame
    ``layer.affine`` is in.
    """
    try:
        m = np.asarray(layer.affine.affine_matrix, dtype=np.float64)
    except AttributeError:
        return np.eye(3)
    if m.shape[0] > 3:          # napari pads to the viewer's dimensionality
        m = m[-3:, -3:]
    return world_affine_to_px(m, pixel_size)


def apply_to_layer(layer, pixel_size: float) -> None:
    """Put one napari layer into micrometres.

    Applied to every layer as it is inserted (see ``app.py``) rather than at each
    ``add_image`` / ``add_shapes`` call site: there are more than twenty of those
    spread across eight modules, and a layer added by a tab written later would
    otherwise be silently left in pixels — misplaced relative to everything else,
    which is a worse failure than a wrong scale bar.
    """
    ndim = getattr(layer, "ndim", 2)
    p = float(pixel_size)
    try:
        layer.scale = tuple([p] * ndim)
        layer.units = tuple([MICRON] * ndim)
    except (ValueError, TypeError, AttributeError):
        # A layer that will not take a scale is left alone rather than half-set;
        # a partially converted layer is misplaced, and misplacement is the thing
        # this module exists to prevent.
        return
