"""Drive the built application through its real widgets.

``Rig`` is the handle a caller gets on a running viewer: navigate to a tab, find
a widget by the label the user reads, click it, wait for the worker it started.
It was written for ``scripts/capture_screenshots.py`` and lives here so the
end-to-end tests can use the same one — two copies of "click the button labelled
X and wait for it to finish" would disagree the first time either learned
something about a widget the other did not.

Nothing here imports napari. The Qt imports are at module scope because every
caller has already built a ``QApplication`` by the time it constructs a ``Rig``.
"""
from __future__ import annotations

import time

import numpy as np
from qtpy.QtWidgets import QApplication, QPushButton, QTabWidget, QWidget

# Everything a working session restores that is not the dataset itself is hidden
# by default: registration landmarks, an ARMS scan, patch overlays and prediction
# rasters sit outside the Xenium extent or on top of it, and none of them is what
# a reference image of the viewer should be showing. An allow-list rather than a
# block-list, because a session can restore arbitrarily named overlays and an
# unrecognised one must not end up in a published image.
BASE_LAYERS = ("cell_labels", "morphology_focus")
FRAME_LAYER = "cell_labels"


def process_events(pause=0.08):
    QApplication.processEvents()
    time.sleep(pause)
    QApplication.processEvents()


# The name the capture script has always used, kept so its call sites read the
# same. New code should prefer ``process_events``.
_process_events = process_events


class Rig:
    """The handle a shot's setup gets: the app, plus the few verbs it needs.

    ``strict`` decides what a timeout means. A screenshot of a step that is
    still running is still a picture of that step, so the capture script sets
    ``strict=False`` and a timeout only prints; a test wants a failure, so the
    default raises ``TimeoutError``.
    """

    def __init__(self, viewer, ctx, dock, panel, strict: bool = True):
        self.viewer = viewer
        self.ctx = ctx
        self.dock = dock
        self.panel = panel
        self.strict = strict

    @classmethod
    def from_app(cls, viewer, ctx, app_state: dict, strict: bool = True) -> "Rig":
        """Build a rig from what ``palms.app.create_viewer`` returns.

        The dock's own widget is the outer ``QTabWidget`` of the five groups —
        derived here rather than passed, so a caller cannot get the pair wrong.
        """
        dock = app_state["dock_widget"]
        if dock is None:
            raise RuntimeError("no Controls dock — was a dataset loaded?")
        return cls(viewer, ctx, dock, dock.widget(), strict=strict)

    # ── navigation ──────────────────────────────────────────────────────
    def navigate(self, outer, inner):
        self.panel.setCurrentIndex(outer)
        process_events(0.05)
        sub = self.panel.currentWidget()
        if sub is not None and hasattr(sub, "setCurrentIndex"):
            sub.setCurrentIndex(inner)
        process_events(0.08)

    def tab(self, group: str, name: str):
        """Navigate by the labels the user reads: ``rig.tab("Cells", "Clustering")``.

        Indices are what the capture script uses, because its shot table is
        ordered by them anyway. A test should not encode them: reordering a tab
        would silently point every assertion at a different page instead of
        failing.
        """
        outer = self._index_of(self.panel, group)
        self.panel.setCurrentIndex(outer)
        process_events(0.05)
        sub = self.panel.currentWidget()
        if not isinstance(sub, QTabWidget):
            raise LookupError(f"tab group {group!r} holds no sub-tabs")
        sub.setCurrentIndex(self._index_of(sub, name))
        process_events(0.08)
        return self.page

    @staticmethod
    def _index_of(tabs: QTabWidget, label: str) -> int:
        want = Rig._norm(label)
        for i in range(tabs.count()):
            if Rig._norm(tabs.tabText(i)) == want:
                return i
        have = ", ".join(tabs.tabText(i) for i in range(tabs.count()))
        raise LookupError(f"no tab labelled {label!r} (have: {have})")

    @property
    def page(self):
        sub = self.panel.currentWidget()
        return sub.currentWidget() if hasattr(sub, "currentWidget") else sub

    # ── widgets ─────────────────────────────────────────────────────────
    @staticmethod
    def _norm(text):
        return "".join(ch for ch in str(text).lower() if ch.isalnum())

    def mg(self, label):
        """A magicgui widget on the current page, by its label.

        magicgui stores a back-reference on the Qt widget it wraps
        (``native._magic_widget``), so the whole tab can be driven through the
        real widget objects — ``.value = x``, ``.native.click()`` — without the
        app having to export them and without matching on Qt layout structure.
        """
        want = self._norm(label)
        for w in self.page.findChildren(QWidget):
            mw = getattr(w, "_magic_widget", None)
            if mw is not None and self._norm(getattr(mw, "label", "")) == want:
                return mw
        raise LookupError(f"no magicgui widget labelled {label!r} on this page")

    def qbtn(self, text):
        """A plain QPushButton on the current page, by its text."""
        want = self._norm(text)
        for b in self.page.findChildren(QPushButton):
            if self._norm(b.text()) == want:
                return b
        raise LookupError(f"no QPushButton {text!r} on this page")

    def click(self, label):
        try:
            self.mg(label).native.click()
        except LookupError:
            self.qbtn(label).click()
        process_events(0.2)

    def _timed_out(self, what, timeout, strict):
        strict = self.strict if strict is None else strict
        message = f"timed out after {timeout}s waiting for {what}"
        if strict:
            raise TimeoutError(message)
        print(f"    ! {message}")
        return False

    def wait_idle(self, label, timeout=600, strict=None):
        """Wait for a button that disables itself while its worker runs."""
        try:
            widget = self.mg(label).native
        except LookupError:
            widget = self.qbtn(label)
        deadline = time.time() + timeout
        # give the worker a moment to actually disable it
        for _ in range(10):
            process_events(0.1)
            if not widget.isEnabled():
                break
        while time.time() < deadline:
            process_events(0.25)
            if widget.isEnabled():
                return True
        return self._timed_out(repr(label), timeout, strict)

    def wait_for(self, cond, timeout=600, what="", strict=None):
        deadline = time.time() + timeout
        while time.time() < deadline:
            process_events(0.25)
            if cond():
                return True
        return self._timed_out(what, timeout, strict)

    # ── layers and camera ───────────────────────────────────────────────
    def show_only(self, *prefixes):
        prefixes = prefixes or BASE_LAYERS
        for layer in self.viewer.layers:
            layer.visible = layer.name.startswith(tuple(prefixes))

    def show_also(self, *prefixes):
        for layer in self.viewer.layers:
            if layer.name.startswith(tuple(prefixes)):
                layer.visible = True

    def select(self, name):
        if name in self.viewer.layers:
            self.viewer.layers.selection = {self.viewer.layers[name]}

    def frame(self, name=FRAME_LAYER, margin=0.95, zoom_factor=1.0):
        """Put the camera on one layer's own extent.

        Not viewer.reset_view(): napari's fit_to_view measures
        layers._extent_world_augmented, which ignores `visible`, so a hidden
        ARMS scan far outside the tissue would still set the frame.
        """
        if name not in self.viewer.layers:
            return
        lo, hi = self.viewer.layers[name].extent.world
        size = np.maximum(hi - lo, 1.0)
        self.viewer.camera.center = (0.0, *((lo + hi) / 2.0))
        base = float(np.min(np.array(self.viewer._canvas_size) / size))
        self.viewer.camera.zoom = margin * base * zoom_factor

    def zoom_in(self, factor=4.0, name=FRAME_LAYER, offset=(0.0, 0.0)):
        self.frame(name)
        cy, cx = self.viewer.camera.center[1:]
        lo, hi = self.viewer.layers[name].extent.world
        span = hi - lo
        self.viewer.camera.center = (0.0, cy + offset[0] * span[0], cx + offset[1] * span[1])
        self.viewer.camera.zoom *= factor
