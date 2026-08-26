"""Put a matplotlib figure on screen without going through ``plt.show``.

Two things live here: a figure → ``QPixmap`` thumbnail (which the Notebook tab
has done since it was written, and which the Plots dock now also needs), and a
plain window holding a live canvas plus matplotlib's navigation toolbar.

The point of the second one is that it does **not** use pyplot's interactive
backend. Two workers called ``matplotlib.use('Agg')`` globally and never
restored it, so after a user saved a UMAP or ran a marker plot every later
``plt.show(block=False)`` in the session was a silent no-op — figures simply
stopped appearing, with no error. Embedding our own ``FigureCanvasQTAgg``
sidesteps the process-wide backend entirely: the figure is drawn because we drew
it, not because pyplot happened to be pointed at Qt.
"""

from __future__ import annotations

#: Thumbnail width in the Notebook tab and the Plots dock.
THUMBNAIL_WIDTH = 580


def fig_to_pixmap(fig, max_width: int = THUMBNAIL_WIDTH):
    """Render *fig* to a ``QPixmap``, scaled down to *max_width*.

    Uses the Agg canvas explicitly rather than whatever backend pyplot is on,
    so the result does not depend on process-wide state.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from qtpy.QtCore import Qt
    from qtpy.QtGui import QImage, QPixmap

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = canvas.buffer_rgba()
    w, h = canvas.get_width_height()
    qimg = QImage(bytes(buf), w, h, 4 * w, QImage.Format.Format_RGBA8888)
    pixmap = QPixmap.fromImage(qimg)
    if pixmap.width() > max_width:
        pixmap = pixmap.scaledToWidth(
            max_width, Qt.TransformationMode.SmoothTransformation)
    return pixmap


def open_figure_window(fig, title: str, parent=None):
    """Show *fig* in its own resizable window with pan/zoom/save.

    Returns the window, which the caller must keep a reference to — a top-level
    Qt widget with no parent is garbage-collected the moment it goes out of
    scope, and closes itself in front of the user.
    """
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg, NavigationToolbar2QT,
    )
    from qtpy.QtWidgets import QVBoxLayout, QWidget

    window = QWidget(parent)
    window.setWindowTitle(title)
    canvas = FigureCanvasQTAgg(fig)
    layout = QVBoxLayout()
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(NavigationToolbar2QT(canvas, window))
    layout.addWidget(canvas)
    window.setLayout(layout)

    width, height = fig.get_size_inches() * fig.dpi
    window.resize(int(min(width, 1400)) + 20, int(min(height, 1000)) + 60)
    window.show()
    canvas.draw_idle()
    return window


def detach_from_pyplot(fig) -> None:
    """Drop pyplot's reference to *fig* while keeping the figure usable.

    Figures built with ``plt.subplots`` stay in pyplot's registry until closed,
    and the viewer produces one per user action — so without this a long session
    accumulates every figure it ever drew. ``plt.close`` destroys the *manager*,
    not the ``Figure``: the object stays alive as long as we hold it, and can be
    given a new canvas, which is exactly what ``open_figure_window`` does.
    """
    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception:
        # Never let bookkeeping lose a figure the user just asked for.
        pass
