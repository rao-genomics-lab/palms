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


def to_figure(obj):
    """Return the matplotlib ``Figure`` behind *obj*.

    Not every "figure" in this codebase is one. ``sc.pl.rank_genes_groups_dotplot
    (return_fig=True)`` hands back a scanpy ``DotPlot``, which quacks *almost*
    like a Figure — it has ``savefig``, which is why the old one-format
    ``auto_save_plot`` never noticed — but has no canvas, so anything that draws
    it fails with ``AttributeError: 'DotPlot' object has no attribute
    'set_canvas'``.

    Two details make this less obvious than it looks:

    * a ``BasePlot``'s ``.fig`` is ``None`` until the plot is actually built, so
      reading the attribute is not enough. ``get_axes()`` builds it if needed and
      is idempotent (it guards on ``ax_dict is None``);
    * ``BasePlot.savefig`` writes ``plt.gcf()``, not ``self.fig``. Resolving to
      the concrete Figure and saving *that* also removes the chance of writing
      whichever figure happens to be current instead.

    Raises ``TypeError`` naming the type, rather than failing later inside Qt.
    """
    from matplotlib.figure import Figure

    if isinstance(obj, Figure):
        return obj
    if hasattr(obj, "get_axes") and not isinstance(obj, Figure):
        # scanpy BasePlot: build the figure if it has not been built yet.
        axes = obj.get_axes()
        figure = getattr(obj, "fig", None)
        if isinstance(figure, Figure):
            return figure
        if isinstance(axes, dict) and axes:
            return next(iter(axes.values())).figure
    figure = getattr(obj, "figure", None)      # an Axes
    if isinstance(figure, Figure):
        return figure
    raise TypeError(
        f"expected a matplotlib Figure or a scanpy plot object, got "
        f"{type(obj).__module__}.{type(obj).__name__}"
    )


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
