"""Cells -> Cell Coloring -> "White background" must survive a theme change.

Reported from a cluster desktop: ticking the box did nothing. The handler set
``canvas.bgcolor``, which napari rewrites from the theme on every theme or
theme-canvas event, so anything that fired one -- a theme setting sync, the
``system`` theme -- undid the checkbox silently. ``background_color_override``
is the attribute napari keeps across those events. On napari 0.9 the handler
did not even get that far: 0.9 dropped ``VispyCanvas.bgcolor`` (the override
moved to the public ``viewer.canvas`` model), so it raised ``AttributeError``.
Assertions therefore read the *drawn* colour from the vispy scene canvas, which
both versions keep, rather than either version's API.

Builds the real tab against a stub context, as the rest of the tab is not under
test here.
"""
from __future__ import annotations

import gc
import types

import pytest


def drawn(viewer) -> str:
    """The colour vispy clears the canvas to, in either napari version."""
    return viewer.window._qt_viewer.canvas._scene_canvas.bgcolor.hex


@pytest.fixture
def tab(qapp):
    napari = pytest.importorskip("napari")
    from magicgui.widgets import CheckBox
    from palms.tabs import tab_cell_coloring

    viewer = napari.Viewer(show=False)
    recorded = []

    class Ctx(types.SimpleNamespace):
        def __getattr__(self, name):           # every other ctx callable: no-op
            return lambda *a, **kw: None

    ctx = Ctx(viewer=viewer, state={"cluster_checkboxes": {}},
              gene_names=["A"], clustering_names=[])
    ctx.record_node = lambda node_id, *a, **kw: recorded.append(node_id)
    widget, _ = tab_cell_coloring.build_tab(ctx)
    box = next(o for o in gc.get_objects()
               if isinstance(o, CheckBox) and o.label == "White background"
               and o.native.window() is widget.window())
    try:
        yield viewer, box, recorded
    finally:
        viewer.close()


def test_ticking_turns_the_canvas_white(tab):
    viewer, box, recorded = tab
    box.value = True
    assert drawn(viewer) == "#ffffff"
    assert recorded == ["viewer:background"]


def test_white_survives_a_theme_change(tab):
    viewer, box, _ = tab
    box.value = True
    viewer.theme = "light"
    viewer.theme = "dark"
    assert drawn(viewer) == "#ffffff"


def test_unticking_hands_the_canvas_back_to_the_theme(tab):
    from napari.utils.theme import get_theme
    from vispy.color import Color

    viewer, box, _ = tab
    box.value = True
    box.value = False
    viewer.theme = "light"
    viewer.theme = "dark"     # dark, because the light theme's canvas is white
    assert drawn(viewer) == Color(get_theme("dark").canvas.as_hex()).hex
    assert drawn(viewer) != "#ffffff"
