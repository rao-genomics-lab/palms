"""Physical units for the napari canvas — one boundary, stated once.

The viewer's scale bar read in pixels because napari measures the *world*, and
every layer was added with the default scale of 1, so one world unit was one
Xenium morphology pixel. Giving each layer ``scale = (pixel_size, pixel_size)``
and ``units = ("um", "um")`` makes one world unit a micrometre and the scale bar
correct, in µm or mm as the zoom warrants.

That change is not free, and the cost is concentrated in one place. It is also not
a choice — the display-level route it replaces is the obvious one, and napari
removed it.

Why there is no display-level fix
---------------------------------
The natural approach is to tell the scale bar the pixel size and change nothing
else. napari <= 0.5 supported exactly that: ``viewer.scale_bar.unit`` held a
``pint.Quantity``, and the bar's reading is computed as ``self._unit *
desired_length`` (``_vispy/overlays/scale_bar.py``), so ``"0.2125 um"`` carried
its own magnitude and did the whole job.

That attribute is now a deprecated no-op — its setter warns and does nothing, and
napari's own test asserts it reads back ``None``; it disappears in napari 0.9.0.
The unit is derived from the layers instead: ``unit = viewer.layers.units[-1]``
followed by ``self._unit = unit * 1`` — a ``pint.Unit``, which carries a dimension
but no magnitude, promoted to a Quantity of magnitude *one*.

So the conversion factor has nowhere left to live except the world coordinates
themselves, which is what makes ``layer.scale`` the only lever. (The same
narrowing one level down is why ``layer.units = "0.2125 um"`` silently keeps only
``micrometer``, see :data:`MICRON`. And ``viewer.layers.units`` is no help either:
it relabels axes and never converts scales, and on layers still in pixels it
raises, because ``pixel`` is dimensionless while ``um`` is a length.)
``tests/test_units.py`` pins all three, so if napari restores a display-level
route this module can be deleted rather than quietly kept.

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


def apply_to_layer(layer, pixel_size: float) -> bool:
    """Put one napari layer into micrometres. True if it took.

    Applied to every layer as it is inserted (see ``app.py``) rather than at each
    ``add_image`` / ``add_shapes`` call site: there are more than twenty of those
    spread across eight modules, and a layer added by a tab written later would
    otherwise be silently left in pixels — misplaced relative to everything else,
    which is a worse failure than a wrong scale bar.

    The return value exists so the caller can tell a layer that was converted from
    one that refused. That distinction used to be invisible, and it is the whole
    difference between napari's units warning being noise and being a real report
    — see :func:`quiet_insertion`.
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
        return False
    return True


# ── napari's transient units warning ─────────────────────────────────────────
#
# Adding a layer emits, once per layer after the first:
#
#     Inconsistent units across layers; units will not be used for rendering.
#
# It is spurious here, and the reason is a timing window we do not control. A new
# layer arrives carrying napari's default (pixel) units while every existing layer
# is already in micrometres. napari's own canvas handler runs on the same
# `inserted` event and triggers a draw before ours has stamped the newcomer:
#
#     on_draw -> _update_scenegraph -> add_layer_visual_mapping -> _update_world_units
#
# which reads `viewer.layers.extent.units`, finds a disagreement, and warns. By the
# time the next draw happens the layer is stamped and everything agrees — measured:
# `extent.units` settles to micrometres and the scale bar renders correctly. So the
# warning describes a state that has already stopped being true.
#
# Connection order does not help (tried: `position="first"` on `inserted` changes
# nothing, because the draw is not ordered by it), and setting units at the ~20
# `add_*` call sites is the design this module exists to avoid.
#
# So the message is dropped, but *only* for the span of one insertion — from
# `inserting` to the end of our own `inserted` handler. A genuine mismatch outside
# that window still reaches the user untouched. And the one real failure this could
# have masked — a layer that refuses the scale, which really would leave the world
# inconsistent — is reported explicitly instead, naming the layer, which is a better
# message than the one being suppressed.

_TRANSIENT_UNITS_WARNING = "Inconsistent units across layers"


class _InsertionQuiet:
    """Drops napari's transient units warning while a layer is being inserted.

    The patch exists **only for the span of the window** and the previous value is
    restored exactly. A version of this that installed a permanent wrapper and
    gated it on a depth counter worked, but it captured whatever ``show_warning``
    happened to be bound at first use and kept wrapping it forever — which made it
    order-dependent with anything else that rebinds the name, and it showed up as a
    test that passed alone and failed in a suite. A wrapper that outlives the reason
    for it is the wrong shape.
    """

    def __init__(self):
        self._depth = 0
        self._module = None
        self._saved = None

    def __enter__(self):
        if self._depth == 0:
            self._patch()
        self._depth += 1
        return self

    def _patch(self) -> None:
        try:
            import napari._vispy.canvas as canvas_mod
        except Exception:                              # pragma: no cover - defensive
            return
        original = getattr(canvas_mod, "show_warning", None)
        if not callable(original):
            # napari moved it. Say nothing and warn about nothing: a noisier console
            # is not worth a crash at startup.
            return

        def filtered(message, *args, **kwargs):
            if _TRANSIENT_UNITS_WARNING in str(message):
                return None
            return original(message, *args, **kwargs)

        self._module, self._saved = canvas_mod, original
        canvas_mod.show_warning = filtered

    def __exit__(self, *exc):
        # Never leave the window open: a stuck counter would suppress the warning for
        # the rest of the session, which is the failure this approach exists to avoid.
        self._depth = max(0, self._depth - 1)
        if self._depth == 0 and self._module is not None:
            self._module.show_warning = self._saved
            self._module = self._saved = None
        return False


def quiet_insertion() -> _InsertionQuiet:
    """The suppression window. One shared instance, so nesting counts correctly."""
    return _QUIET


_QUIET = _InsertionQuiet()
