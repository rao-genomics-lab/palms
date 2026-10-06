"""Diagnose why the napari canvas background does not change on this machine.

Standalone: imports napari only, needs no dataset and no PALMS install. Run it in
the same conda env PALMS runs in, on the machine where "White background" does
nothing (e.g. a cluster desktop):

    python scripts/diagnose_canvas_bg.py

It prints the stack (napari / vispy / Qt / GL renderer, napari theme), then
sets the canvas white twice -- once through ``VispyCanvas.bgcolor`` (what PALMS
used to do; gone in napari 0.9) and once through ``background_color_override``
(what it does now, on ``viewer.canvas`` from 0.9) -- reads the colour back, grabs the rendered pixels, and ends with one
verdict per route:

- ``works``        the canvas reports white and renders white;
- ``reset``        something set the colour back after we set it;
- ``not rendered`` the canvas reports white but the pixels are not -- the GL
                   stack (remote X, VirtualGL, software Mesa) is not showing
                   the clear colour, which no code change in PALMS can fix.
"""
from __future__ import annotations

import platform
import sys

import numpy as np


def _versions():
    import napari
    import qtpy
    import vispy

    lines = [
        f"python   {sys.version.split()[0]} ({platform.platform()})",
        f"napari   {napari.__version__}",
        f"vispy    {vispy.__version__}",
        f"qt       {qtpy.API_NAME} {qtpy.QT_VERSION}",
    ]
    try:
        import OpenGL
        lines.append(f"pyopengl {OpenGL.__version__}")
    except ImportError:
        lines.append("pyopengl (not importable)")
    try:
        from napari.settings import get_settings
        lines.append(f"theme    {get_settings().appearance.theme}")
    except Exception as e:  # noqa: BLE001 -- diagnostic, report and carry on
        lines.append(f"theme    (unreadable: {e})")
    return lines


def _gl_info(canvas):
    try:
        from vispy.gloo import gl
        canvas._scene_canvas.set_current()
        return (f"GL renderer {gl.glGetParameter(gl.GL_RENDERER)} | "
                f"vendor {gl.glGetParameter(gl.GL_VENDOR)} | "
                f"version {gl.glGetParameter(gl.GL_VERSION)}")
    except Exception as e:  # noqa: BLE001
        return f"GL info unavailable: {e}"


def _is_white(px) -> bool:
    return bool(np.all(np.asarray(px)[:3] >= 250))


def _pixels(viewer, canvas):
    """Centre pixel from napari's screenshot and from Qt's framebuffer grab."""
    shot = viewer.screenshot(canvas_only=True, flash=False)
    h, w = shot.shape[:2]
    napari_px = tuple(int(v) for v in shot[h // 2, w // 2])
    try:
        img = canvas.native.grabFramebuffer()
        c = img.pixelColor(img.width() // 2, img.height() // 2)
        qt_px = (c.red(), c.green(), c.blue(), c.alpha())
    except Exception as e:  # noqa: BLE001
        qt_px = f"unavailable: {e}"
    return napari_px, qt_px


def main():
    import napari
    from qtpy.QtCore import QTimer
    from qtpy.QtWidgets import QApplication

    print("\n".join(_versions()))
    viewer = napari.Viewer(show=True, title="canvas background diagnostic")
    # An empty viewer shows napari's welcome screen instead of the canvas, so
    # give it one empty layer; it draws nothing over the background.
    viewer.add_points(name="empty")
    canvas = viewer.window._qt_viewer.canvas
    print(_gl_info(canvas))

    # napari 0.9 moved the override to the public ``viewer.canvas`` model and
    # dropped ``VispyCanvas.bgcolor``; 0.8 has no ``viewer.canvas``. What is
    # *drawn* is the vispy scene canvas's colour in both, so read that.
    model = getattr(viewer, "canvas", None)
    if model is None or not hasattr(model, "background_color_override"):
        model = canvas

    def drawn():
        return canvas._scene_canvas.bgcolor.hex

    print(f"initial  drawn={drawn()} override={model.background_color_override}")

    verdicts = {}

    def check(route):
        read_back = drawn()
        napari_px, qt_px = _pixels(viewer, canvas)
        print(f"[{route}] read back {read_back}; screenshot px {napari_px}; "
              f"framebuffer px {qt_px}")
        if read_back.lower() != "#ffffff":
            verdicts[route] = "reset"
        elif not _is_white(napari_px):
            verdicts[route] = "not rendered"
        else:
            verdicts[route] = "works"

    def set_bgcolor():
        if not hasattr(canvas, "bgcolor"):
            verdicts["bgcolor"] = ("unsupported: this napari has no "
                                   "VispyCanvas.bgcolor (removed in 0.9)")
            return set_override()
        canvas.bgcolor = (1, 1, 1, 1)
        QTimer.singleShot(1000, lambda: (check("bgcolor"), reset()))

    def reset():
        canvas.bgcolor = (0, 0, 0, 1)
        QTimer.singleShot(500, set_override)

    def set_override():
        model.background_color_override = "white"
        QTimer.singleShot(1000, lambda: (check("override"), finish()))

    def finish():
        for route, verdict in verdicts.items():
            print(f"VERDICT {route}: {verdict}")
        viewer.close()
        QApplication.instance().quit()

    QTimer.singleShot(1000, set_bgcolor)
    napari.run()


if __name__ == "__main__":
    main()
