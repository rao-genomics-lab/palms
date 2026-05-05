"""Notebook execution engine — runs code in napari's IPython namespace."""
from __future__ import annotations

import io
import sys
from dataclasses import dataclass, field


@dataclass
class CellResult:
    stdout: str = ""
    stderr: str = ""
    error: str | None = None
    figures: list = field(default_factory=list)  # list of QPixmap
    result_repr: str = ""
    success: bool = True


class NotebookEngine:
    """Execute code cells in napari's in-process IPython kernel."""

    def __init__(self, viewer):
        console = viewer.window._qt_viewer.console
        kernel = console.kernel_manager.kernel
        self.shell = kernel.shell

    def run_cell(self, code: str) -> CellResult:
        import matplotlib.pyplot as plt
        from qtpy.QtCore import Qt

        # 1. Snapshot existing figures
        pre_figs = set(plt.get_fignums())

        # 2. Intercept plt.show to prevent popup windows
        _orig_show = plt.show
        plt.show = lambda *a, **kw: None

        # 3. Capture stdout/stderr
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = stdout_buf
        sys.stderr = stderr_buf

        try:
            result = self.shell.run_cell(code, store_history=True)
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            plt.show = _orig_show

        # 4. Collect new figures
        new_fig_nums = sorted(set(plt.get_fignums()) - pre_figs)
        pixmaps = []
        for num in new_fig_nums:
            fig = plt.figure(num)
            pixmaps.append(self._fig_to_pixmap(fig))
            plt.close(fig)

        # 5. Build result
        error_str = None
        if result.error_in_exec is not None:
            import traceback
            error_str = "".join(traceback.format_exception(
                type(result.error_in_exec), result.error_in_exec,
                result.error_in_exec.__traceback__,
            ))
        elif result.error_before_exec is not None:
            import traceback
            error_str = "".join(traceback.format_exception(
                type(result.error_before_exec), result.error_before_exec,
                result.error_before_exec.__traceback__,
            ))

        result_repr = ""
        if result.result is not None:
            try:
                result_repr = repr(result.result)
            except Exception:
                result_repr = "<repr failed>"

        return CellResult(
            stdout=stdout_buf.getvalue(),
            stderr=stderr_buf.getvalue(),
            error=error_str,
            figures=pixmaps,
            result_repr=result_repr,
            success=result.success,
        )

    def _fig_to_pixmap(self, fig, max_width=580):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from qtpy.QtGui import QImage, QPixmap
        from qtpy.QtCore import Qt

        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        buf = canvas.buffer_rgba()
        w, h = canvas.get_width_height()
        qimg = QImage(bytes(buf), w, h, 4 * w, QImage.Format_RGBA8888)
        pixmap = QPixmap.fromImage(qimg)
        if pixmap.width() > max_width:
            pixmap = pixmap.scaledToWidth(max_width, Qt.SmoothTransformation)
        return pixmap
