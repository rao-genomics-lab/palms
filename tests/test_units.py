"""Micrometre world units, and the registration they could silently break.

Putting the canvas into micrometres means giving every layer a ``scale``. napari
composes ``world = affine(scale(data))`` — the affine lands *after* the scale, so
its translation is in world units, while every affine this codebase stores is in
Xenium pixels. Get that boundary wrong and every registered overlay shifts by
``1 / pixel_size`` (~4.7x) with nothing raised, nothing logged, and an export that
still opens.

So the assertions here are about *placement*: where a known pixel ends up. A test
that only checked "the scale bar says µm" would pass with every overlay in the
wrong place.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from xenium_viewer.utils import units

PX = 0.2125          # the Xenium default, and the fallback in transcript_index


def _similarity(theta_deg, scale, ty, tx):
    t = np.deg2rad(theta_deg)
    c, s = np.cos(t) * scale, np.sin(t) * scale
    return np.array([[c, -s, ty], [s, c, tx], [0, 0, 1.0]])


# ── the conversion itself ────────────────────────────────────────────────────

def test_only_the_translation_changes_units():
    """A rotation is the same rotation in any unit; an offset is not."""
    a = _similarity(30.0, 1.4, ty=100.0, tx=-50.0)
    w = units.px_affine_to_world(a, PX)

    assert np.allclose(w[:2, :2], a[:2, :2]), "the linear part must be unit-free"
    assert w[0, 2] == pytest.approx(100.0 * PX)
    assert w[1, 2] == pytest.approx(-50.0 * PX)


def test_round_trip_is_exact_enough_to_persist():
    a = _similarity(11.0, 0.77, ty=1234.5, tx=-987.6)
    back = units.world_affine_to_px(units.px_affine_to_world(a, PX), PX)
    assert np.allclose(back, a, atol=1e-9)


def test_identity_stays_identity():
    """A layer with no registration must not acquire one from a unit change."""
    assert np.allclose(units.px_affine_to_world(np.eye(3), PX), np.eye(3))
    assert np.allclose(units.world_affine_to_px(np.eye(3), PX), np.eye(3))


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_a_nonsense_pixel_size_is_refused(bad):
    """Silently producing a degenerate transform is the failure to avoid."""
    with pytest.raises(ValueError):
        units.scale_for(bad)
    with pytest.raises(ValueError):
        units.px_affine_to_world(np.eye(3), bad)


def test_a_non_3x3_affine_is_refused():
    with pytest.raises(ValueError):
        units.px_affine_to_world(np.eye(4), PX)


# ── the napari behaviour the conversion exists for ───────────────────────────

def _layer(**kw):
    """A bare napari Image layer — deliberately *not* added to a Viewer.

    Everything asserted here (``data_to_world``, ``scale``, ``units``, ``affine``,
    ``extent``) is layer-level transform maths that needs no canvas. Building a
    Viewer would drag in vispy and an OpenGL context, which CI does not have —
    the suite would fail there with ``GLXBadFBConfig`` for reasons that have
    nothing to do with units.
    """
    from napari.layers import Image
    return Image(np.zeros((8, 8)), **kw)


def test_napari_applies_affine_after_scale():
    """The premise. If this ever stops holding, everything below is wrong.

    Pinned as a test rather than left as a comment because the whole conversion
    exists only because of it, and it is an upstream behaviour we do not control.
    """
    a = np.array([[1, 0, 100.0], [0, 1, 50.0], [0, 0, 1]])
    layer = _layer(affine=a)
    assert layer.data_to_world((0, 0)) == pytest.approx((100.0, 50.0))

    layer.scale = (PX, PX)
    assert layer.data_to_world((0, 0)) == pytest.approx((100.0, 50.0)), (
        "napari applied the affine before the scale — the unit conversion in "
        "utils/units.py is built on the opposite assumption"
    )


def test_a_registered_overlay_does_not_move_when_the_world_becomes_microns():
    """The regression this whole change risks, stated as an assertion.

    An H&E pixel sat on some tissue feature at Xenium pixel P. After the switch to
    micrometres it must sit at ``P * pixel_size`` — the *same tissue*, relabelled —
    and not at P, which is where it lands if the affine is left in pixels.
    """
    affine_px = _similarity(23.0, 0.6, ty=800.0, tx=-250.0)
    own = np.array([120.0, 200.0, 1.0])
    expected_px = affine_px @ own

    layer = _layer()
    units.apply_to_layer(layer, PX)
    layer.affine = units.px_affine_to_world(affine_px, PX)

    got = layer.data_to_world((own[0], own[1]))
    assert got == pytest.approx((expected_px[0] * PX, expected_px[1] * PX)), (
        "the overlay moved; a stored pixel affine was applied as if it were µm"
    )


def test_unscaled_and_scaled_layers_agree_on_where_a_point_is():
    """Registered and unregistered layers must land in the same world.

    A cell centroid at pixel (y, x) on the labels layer and the H&E pixel that
    was registered onto it have to coincide, or the overlay is misaligned even
    though each layer alone looks right.
    """
    p_px = np.array([300.0, 450.0])

    plain = _layer()
    units.apply_to_layer(plain, PX)

    affine_px = np.array([[1, 0, 100.0], [0, 1, 200.0], [0, 0, 1]])
    own = p_px - np.array([100.0, 200.0])           # the own-pixel mapping to p_px
    overlay = _layer()
    units.apply_to_layer(overlay, PX)
    overlay.affine = units.px_affine_to_world(affine_px, PX)

    assert plain.data_to_world(tuple(p_px)) == pytest.approx(
        overlay.data_to_world(tuple(own))
    )


def test_layer_affine_px_reads_back_what_was_written():
    """Persisting a registration must recover pixels, whatever napari padded it to."""
    affine_px = _similarity(7.0, 1.3, ty=42.0, tx=-17.0)
    layer = _layer()
    units.apply_to_layer(layer, PX)
    layer.affine = units.px_affine_to_world(affine_px, PX)

    assert np.allclose(units.layer_affine_px(layer, PX), affine_px, atol=1e-9)


def test_a_layer_with_no_affine_reads_back_as_identity():
    layer = _layer()
    units.apply_to_layer(layer, PX)
    assert np.allclose(units.layer_affine_px(layer, PX), np.eye(3), atol=1e-9)


def test_units_reach_the_layer_so_the_scale_bar_can_read_them():
    """``layer.units`` is what napari 0.8 uses; ScaleBarOverlay has no `unit`."""
    from napari.components.overlays import ScaleBarOverlay

    assert "unit" not in ScaleBarOverlay.model_fields, (
        "napari grew back scale_bar.unit — a scaled unit string would be simpler "
        "than scaling every layer, but check it keeps its magnitude first"
    )

    layer = _layer()
    units.apply_to_layer(layer, PX)
    assert all(str(u) == "micrometer" for u in layer.units)
    # napari's world extent spans pixel *centres*, so an 8-px axis ends at 7.
    assert layer.extent.world[1][-1] == pytest.approx(7 * PX)


def test_a_scaled_unit_string_silently_loses_its_magnitude():
    """Why the magnitude lives in `scale` and not in the unit.

    ``"0.2125 um"`` is *accepted* and then discarded, leaving one pixel labelled
    as one micrometre — a wrong scale bar with no error anywhere. Pinned so
    nobody re-introduces the napari-0.5-era shortcut.
    """
    layer = _layer()
    layer.units = "0.2125 um"
    assert all(str(u) == "micrometer" for u in layer.units)
    assert layer.extent.world[1][-1] == pytest.approx(7.0), (
        "napari started honouring a scaled unit string; units.py could be simplified"
    )


# ── the display-level route, and why it is not available ─────────────────────
#
# The obvious fix is to tell the scale bar the pixel size and change nothing
# else. That is exactly what napari <= 0.5 supported, and it is the first thing
# anyone will reach for. These two tests record that napari removed it, so the
# question is answered by the suite rather than by memory — and so that the day
# napari brings it back, they fail and prompt simplifying utils/units.py.


def test_the_scale_bar_can_no_longer_be_told_a_unit():
    """``scale_bar.unit`` is a deprecated no-op, removed in napari 0.9.0.

    It used to hold a ``pint.Quantity``, so ``"0.2125 um"`` carried its magnitude
    and did the whole job. The unit is now derived from the layers instead
    (``_vispy/overlays/scale_bar.py``: ``unit = viewer.layers.units[-1]`` then
    ``self._unit = unit * 1``) — a ``pint.Unit`` promoted to a Quantity of
    magnitude *one*. There is nowhere left to put a scale factor except the world
    coordinates, which is the entire reason utils/units.py exists.

    Setting it does not raise; it warns and does nothing.
    """
    from napari.components.overlays import ScaleBarOverlay

    bar = ScaleBarOverlay()
    with pytest.warns(FutureWarning, match="no longer has any effect"):
        bar.unit = "0.2125 um"
    with pytest.warns(FutureWarning, match="always returns None"):
        assert bar.unit is None


def test_relabelling_the_layer_list_does_not_convert_anything():
    """The other display-level candidate: ``viewer.layers.units``.

    It sets labels, never scales — and on layers still in pixels it does not even
    do that, because ``pixel`` is dimensionless and ``um`` is a length.
    """
    from napari.components.layerlist import LayerList

    layers = LayerList()
    layers.append(_layer())
    assert all(str(u) == "pixel" for u in layers.units)

    with pytest.raises(ValueError, match="dimensionality"):
        layers.units = ("um", "um")


# ── the minimap, which reads the camera ──────────────────────────────────────

def test_minimap_maps_a_click_to_world_not_to_pixels():
    """The minimap converts between the camera and the morphology shape.

    The camera is in world units and the shape is in pixels, so without the
    pixel size the click would land ~4.7x off — inside the tissue, plausibly,
    which is what makes it easy to miss.
    """
    from xenium_viewer.utils.minimap_widget import MinimapWidget

    full_shape = (4000, 5000)
    scaled = MinimapWidget.__new__(MinimapWidget)
    scaled._morph_full_shape_yx = full_shape
    scaled._pixel_size = PX
    scaled._world_shape_yx = (full_shape[0] * PX, full_shape[1] * PX)

    assert scaled._world_shape_yx == pytest.approx((4000 * PX, 5000 * PX))
    # An unscaled viewer (pixel_size=1.0) must behave exactly as before.
    assert (full_shape[0] * 1.0, full_shape[1] * 1.0) == full_shape


# ── the console noise the conversion used to cost ────────────────────────────
#
# Every layer added after the first emitted:
#
#     Inconsistent units across layers; units will not be used for rendering.
#
# ~20 times on a real dataset load. It is spurious: the new layer briefly carries
# napari's default pixel units while the rest are in µm, napari's canvas handler
# triggers a draw inside that window, and by the next draw everything agrees.

# ── the suppression mechanism, canvas-free (this is what guards it in CI) ────

def _canvas_module():
    """napari's canvas module — importable without a GL context."""
    return pytest.importorskip("napari._vispy.canvas")


def test_the_target_message_is_dropped_inside_the_window(monkeypatch):
    canvas_mod = _canvas_module()
    seen = []
    monkeypatch.setattr(canvas_mod, "show_warning",
                        lambda m, *a, **k: seen.append(str(m)), raising=False)

    with units.quiet_insertion():
        canvas_mod.show_warning("Inconsistent units across layers; units will not be used.")

    assert seen == [], "the transient message must not reach the user mid-insertion"


def test_other_messages_are_never_dropped(monkeypatch):
    """Suppression is one message, not a quiet mode."""
    canvas_mod = _canvas_module()
    seen = []
    monkeypatch.setattr(canvas_mod, "show_warning",
                        lambda m, *a, **k: seen.append(str(m)), raising=False)

    with units.quiet_insertion():
        canvas_mod.show_warning("something else entirely went wrong")

    assert seen == ["something else entirely went wrong"]


def test_the_target_message_outside_the_window_gets_through(monkeypatch):
    """The window is the whole safety property: a real mismatch still reports."""
    canvas_mod = _canvas_module()
    seen = []
    monkeypatch.setattr(canvas_mod, "show_warning",
                        lambda m, *a, **k: seen.append(str(m)), raising=False)

    with units.quiet_insertion():
        pass
    canvas_mod.show_warning("Inconsistent units across layers; units will not be used.")

    assert len(seen) == 1, "outside an insertion this message is real information"


def test_inserting_opens_the_window_and_inserted_closes_it():
    """The wiring, canvas-free — this is what guards the fix in CI.

    Connecting only `inserted` is the bug: napari's canvas draws before it runs,
    which is exactly when the spurious warning is emitted. So it is not enough that
    a callback exists on each event; firing `inserting` must actually open the
    suppression window, and `inserted` must close it again.
    """
    from napari.layers import Image

    from xenium_viewer.app import _install_unit_scaling

    class _Emitter:
        def __init__(self): self.callbacks = []
        def connect(self, cb, **kw): self.callbacks.append(cb)
        def fire(self, event=None):
            for cb in self.callbacks:
                cb(event)

    class _Events:
        def __init__(self): self.inserting, self.inserted = _Emitter(), _Emitter()

    class _Layers(list):
        def __init__(self): super().__init__(); self.events = _Events()

    class _Viewer:
        class scale_bar: visible = False; colored = True
        def __init__(self): self.layers = _Layers()

    class _Event:
        def __init__(self, value): self.value = value

    v = _Viewer()
    _install_unit_scaling(v, PX)
    quiet = units.quiet_insertion()

    assert quiet._depth == 0
    v.layers.events.inserting.fire(None)
    assert quiet._depth == 1, (
        "firing `inserting` did not open the suppression window — napari's canvas "
        "draws inside this gap, and that is when it warns"
    )

    layer = Image(np.zeros((4, 4), dtype=np.uint8))
    v.layers.events.inserted.fire(_Event(layer))
    assert quiet._depth == 0, "`inserted` must close the window it opened"
    assert str(layer.units[0]) == "micrometer", "and stamp the layer on the way"


def _has_real_display() -> bool:
    """Whether a GPU-backed napari canvas can be built here.

    Deliberately a check on the *environment* rather than a probe, because probing
    is not safe: CI constructs the `Viewer` fine and only fails on the first draw,
    deep in vispy (`glGetParameter` returns None), and the run then ends in a
    segfault at teardown — exit 139. There is nothing to catch that reliably, so
    the Viewer is not built at all where there is no display.

    This matches the project's standing position that the napari GUI proper has no
    automated coverage; what these two tests add is desktop-only confirmation of
    behaviour that is *also* covered canvas-free above.
    """
    if os.environ.get("QT_QPA_PLATFORM", "").startswith("offscreen"):
        return False
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


requires_display = pytest.mark.skipif(
    not _has_real_display(),
    reason="needs a GPU-backed napari canvas; covered canvas-free by the tests above",
)


@pytest.fixture
def viewer():
    """A real offscreen-window napari Viewer, torn down per test."""
    napari = pytest.importorskip("napari")
    v = napari.Viewer(show=False)
    try:
        yield v
    finally:
        v.close()


@requires_display
def test_adding_layers_does_not_warn_about_inconsistent_units(viewer, monkeypatch):
    """The whole point: a normal dataset load must not produce that message."""
    import napari._vispy.canvas as canvas_mod
    from xenium_viewer.app import _install_unit_scaling

    _install_unit_scaling(viewer, PX)

    seen = []
    real = canvas_mod.show_warning
    monkeypatch.setattr(canvas_mod, "show_warning", lambda m, *a, **k: seen.append(str(m)))

    for name in ("morphology_focus", "cell_labels", "nucleus_labels", "he_image"):
        viewer.add_image(np.zeros((16, 16), dtype=np.uint8), name=name)

    noisy = [m for m in seen if "Inconsistent units" in m]
    assert noisy == [], (
        f"adding layers emitted {len(noisy)} spurious units warning(s); the world is "
        "consistent by the next draw, so the message describes a state that is over"
    )
    assert viewer.layers.extent.units is not None, "units must still actually agree"


@requires_display
def test_a_real_mismatch_outside_an_insertion_still_warns(viewer, monkeypatch):
    """Suppression is a window, not a mute.

    If it silenced the message outright, a layer genuinely stuck in pixels would
    become invisible — which is a worse outcome than the noise.
    """
    import napari._vispy.canvas as canvas_mod
    from xenium_viewer.app import _install_unit_scaling

    _install_unit_scaling(viewer, PX)
    viewer.add_image(np.zeros((16, 16), dtype=np.uint8), name="a")
    viewer.add_image(np.zeros((16, 16), dtype=np.uint8), name="b")

    seen = []
    monkeypatch.setattr(canvas_mod, "show_warning", lambda m, *a, **k: seen.append(str(m)))

    viewer.layers["b"].units = ("pixel", "pixel")          # a genuine disagreement
    viewer.window._qt_viewer.canvas._update_world_units()

    assert any("Inconsistent units" in m for m in seen), (
        "a mismatch outside the insertion window must reach the user"
    )


def test_the_suppression_window_always_closes(monkeypatch):
    """A stuck window would mute the warning for the rest of the session."""
    canvas_mod = _canvas_module()
    sentinel = lambda m, *a, **k: None
    monkeypatch.setattr(canvas_mod, "show_warning", sentinel, raising=False)

    quiet = units.quiet_insertion()
    before = quiet._depth
    try:
        with quiet:
            raise RuntimeError("boom")
    except RuntimeError:
        pass

    assert quiet._depth == before, "the depth must unwind even when the body raises"
    assert canvas_mod.show_warning is sentinel, (
        "the patch must be lifted too — a wrapper left behind outlives its reason "
        "and becomes order-dependent with anything else that rebinds the name"
    )


# ── the failure the suppression must not hide ────────────────────────────────

def test_apply_to_layer_reports_whether_it_took():
    from napari.layers import Image

    layer = Image(np.zeros((4, 4), dtype=np.uint8))
    assert units.apply_to_layer(layer, PX) is True
    assert str(layer.units[0]) == "micrometer"


def test_a_layer_that_refuses_the_scale_is_reported_not_swallowed(monkeypatch):
    """The one case napari's warning was right about, said better.

    A layer left in pixels really is misplaced relative to every other layer, so
    this must be surfaced — naming the layer, which napari's message does not.
    """
    from xenium_viewer.utils import reporting

    class _Stubborn:
        name = "wont_scale"
        ndim = 2

        @property
        def scale(self):
            return (1.0, 1.0)

        @scale.setter
        def scale(self, value):
            raise ValueError("this layer does not do scales")

    assert units.apply_to_layer(_Stubborn(), PX) is False

    notified = []
    monkeypatch.setattr(reporting, "_notify", lambda m: notified.append(m))
    reporting._layer_scaling_failures.discard("wont_scale")
    reporting.report_layer_scaling_failure("wont_scale")

    assert notified and "wont_scale" in notified[0]
    assert "pixel coordinates" in notified[0]
    assert "wont_scale" in reporting.layer_scaling_failures()

    notified.clear()
    reporting.report_layer_scaling_failure("wont_scale")
    assert notified == [], "once per layer per session, not once per redraw"
