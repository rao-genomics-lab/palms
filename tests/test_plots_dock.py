"""The Plots dock has to survive being closed.

Reported from a real session: *"the plots dock disappears when I undock it and
try to move it around ... the view/hide toggle for this is not working"*. One
cause, two faces.

napari's dock title bar has an "x", and it does not hide the dock. It calls
``destroyOnClose`` -> ``Window.remove_dock_widget``, which reparents the inner
widget to ``None`` and ``deleteLater()``s the dock. After one click the dock
object was a dangling C++ pointer and the gallery was an orphan; the View-menu
toggle then called ``setVisible`` on that pointer, raised ``RuntimeError`` into a
bare ``except: pass``, and looked broken.

The panel always survives — Python holds it — so the dock is treated as
disposable and rebuilt around the panel whenever it is needed.

Most of this runs against a stub main window; one test uses a real
``napari.Viewer``, because the close-button shim is specific to napari's dock.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xenium_viewer.tabs._helpers import (                        # noqa: E402
    dock_is_alive, ensure_plots_dock, reveal_plots_dock,
)


def _drain():
    """Run Qt's deferred deletions, which is what actually frees a dock."""
    from qtpy.QtCore import QCoreApplication, QEvent
    from qtpy.QtWidgets import QApplication
    QApplication.processEvents()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    QApplication.processEvents()


@pytest.fixture
def stub_viewer(qapp):
    """Just enough of ``viewer.window`` to add and remove a dock widget."""
    from types import SimpleNamespace
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import QDockWidget, QMainWindow

    window = QMainWindow()

    def add_dock_widget(widget, name="", area="bottom"):
        dock = QDockWidget(name, window)
        dock.setWidget(widget)
        window.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
        return dock

    def remove_dock_widget(dock):
        # napari's own teardown, which is the part that bites: the inner widget
        # is detached *and* the dock is scheduled for deletion.
        if dock.widget() is not None:
            dock.widget().setParent(None)
        window.removeDockWidget(dock)
        dock.deleteLater()

    return SimpleNamespace(
        window=SimpleNamespace(add_dock_widget=add_dock_widget,
                               remove_dock_widget=remove_dock_widget,
                               _qt_window=window),
        _window=window,
    )


@pytest.fixture
def app_state(qapp):
    from xenium_viewer.utils.plots_panel import PlotsPanel
    panel = PlotsPanel()
    return {"plots_panel": panel, "plots_dock": None, "plots_action": None}


# ── liveness ─────────────────────────────────────────────────────────────────

def test_dock_is_alive_sees_through_a_deleted_dock(stub_viewer, app_state):
    dock = ensure_plots_dock(stub_viewer, app_state)
    assert dock_is_alive(dock)

    stub_viewer.window.remove_dock_widget(dock)
    _drain()
    assert not dock_is_alive(dock), (
        "a deleted dock must not read as alive — calling into it is what made "
        "the View-menu toggle a silent no-op"
    )


def test_dock_is_alive_tolerates_none():
    assert not dock_is_alive(None)


# ── recreation ───────────────────────────────────────────────────────────────

def test_the_dock_is_rebuilt_after_it_is_destroyed(stub_viewer, app_state):
    first = ensure_plots_dock(stub_viewer, app_state)
    stub_viewer.window.remove_dock_widget(first)
    _drain()

    second = ensure_plots_dock(stub_viewer, app_state)
    assert dock_is_alive(second)
    assert second is not first
    assert app_state["plots_dock"] is second


def test_the_gallery_survives_the_dock_being_destroyed(stub_viewer, app_state):
    """The figures are the point; losing them to a stray click is not on."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panel = app_state["plots_panel"]
    dock = ensure_plots_dock(stub_viewer, app_state)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    panel.add_figure(fig, "kept", [])

    stub_viewer.window.remove_dock_widget(dock)
    _drain()
    assert panel.count() == 1, "the panel must outlive its dock"

    rebuilt = ensure_plots_dock(stub_viewer, app_state)
    assert rebuilt.widget() is panel, "the surviving panel is re-docked"
    assert panel.count() == 1
    plt.close("all")


def test_without_a_panel_there_is_nothing_to_dock(stub_viewer):
    """Before the first dataset finishes loading."""
    assert ensure_plots_dock(stub_viewer, {"plots_panel": None}) is None


def test_an_existing_live_dock_is_reused(stub_viewer, app_state):
    first = ensure_plots_dock(stub_viewer, app_state)
    assert ensure_plots_dock(stub_viewer, app_state) is first


# ── reveal ───────────────────────────────────────────────────────────────────

def test_reveal_shows_a_hidden_dock(stub_viewer, app_state):
    dock = ensure_plots_dock(stub_viewer, app_state)
    dock.setVisible(False)
    stub_viewer._window.show()

    assert reveal_plots_dock(stub_viewer, app_state) is dock
    assert dock.isVisible()


def test_off_screen_geometry_is_recognised(qapp):
    """The predicate on its own, on a widget whose geometry we fully control.

    Deliberately not by dragging a real window off the desktop: a window manager
    clamps that, so such a test passes offscreen (CI) and fails on a desktop.
    """
    from qtpy.QtGui import QGuiApplication
    from qtpy.QtWidgets import QWidget

    from xenium_viewer.tabs._helpers import _is_on_a_screen

    stray = QWidget()
    stray.setGeometry(-100000, -100000, 100, 100)
    assert not _is_on_a_screen(stray)

    centre = QGuiApplication.primaryScreen().availableGeometry().center()
    reachable = QWidget()
    reachable.setGeometry(centre.x() - 50, centre.y() - 50, 100, 100)
    assert _is_on_a_screen(reachable)


def test_reveal_rescues_a_dock_dragged_off_the_desktop(
        stub_viewer, app_state, monkeypatch):
    """A floating dock moved somewhere unreachable looks like a dead menu item.

    Re-showing it where it is changes nothing the user can see, so it is docked
    back into the main window instead. The "unreachable" judgement is stubbed
    here — see the test above for the predicate itself — because a window
    manager will not let a test put a real window out of reach.
    """
    import xenium_viewer.tabs._helpers as helpers

    stub_viewer._window.show()
    dock = ensure_plots_dock(stub_viewer, app_state)
    dock.setFloating(True)
    _drain()
    monkeypatch.setattr(helpers, "_is_on_a_screen", lambda widget: False)

    reveal_plots_dock(stub_viewer, app_state)
    assert not dock.isFloating(), "an unreachable dock should come home"
    assert dock.isVisible()


def test_reveal_leaves_a_reachable_floating_dock_floating(
        stub_viewer, app_state, monkeypatch):
    """Floating on purpose is a legitimate choice; do not undo it."""
    import xenium_viewer.tabs._helpers as helpers

    stub_viewer._window.show()
    dock = ensure_plots_dock(stub_viewer, app_state)
    dock.setFloating(True)
    _drain()
    monkeypatch.setattr(helpers, "_is_on_a_screen", lambda widget: True)

    reveal_plots_dock(stub_viewer, app_state)
    assert dock.isFloating()


# ── the napari-specific half ─────────────────────────────────────────────────

def test_the_close_button_hides_rather_than_destroys(qapp):
    """napari's title-bar "x" calls ``destroyOnClose``; ours must only hide.

    Closing a gallery should not throw away the figures in it. This is the one
    test that needs a real napari dock, because the method being shadowed is
    napari's.
    """
    napari = pytest.importorskip("napari")
    from xenium_viewer.utils.plots_panel import PlotsPanel

    viewer = napari.Viewer(show=False)
    try:
        state = {"plots_panel": PlotsPanel(), "plots_dock": None,
                 "plots_action": None}
        dock = ensure_plots_dock(viewer, state)
        assert dock.destroyOnClose.__self__ is dock, "shim not installed"

        dock.destroyOnClose()            # exactly what the "x" invokes
        _drain()

        assert dock_is_alive(dock), "the close button must not destroy the dock"
        assert not dock.isVisible()
        assert dock.widget() is state["plots_panel"], "panel still attached"
    finally:
        viewer.close()
