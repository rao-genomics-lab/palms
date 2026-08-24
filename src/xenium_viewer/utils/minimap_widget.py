"""MinimapWidget — floating canvas overlay showing tissue thumbnail + viewport rect."""

from __future__ import annotations

import numpy as np
from qtpy.QtWidgets import QWidget
from qtpy.QtGui import QImage, QPixmap, QPainter, QColor, QPen
from qtpy.QtCore import Qt, QRect, QEvent


class MinimapWidget(QWidget):
    """Floating minimap overlay in the top-right corner of the napari canvas.

    Shows the DAPI channel of the morphology thumbnail with a white rectangle
    indicating the current camera viewport. Clicking navigates the camera.
    """

    _WIDGET_W = 200
    _WIDGET_H = 160
    _MARGIN = 10

    def __init__(self, viewer, morph_thumb: np.ndarray,
                 morph_full_shape_yx: tuple, canvas_native: QWidget,
                 pixel_size: float = 1.0):
        super().__init__(canvas_native)
        self._viewer = viewer
        self._morph_full_shape_yx = morph_full_shape_yx  # (H, W) in data pixels
        # The camera reports *world* coordinates, and the viewer's world is
        # micrometres (utils/units.py), while the morphology shape is in pixels.
        # Everything below works in world units so the two cannot be mixed up;
        # the default of 1.0 keeps an unscaled viewer behaving as it always did.
        self._pixel_size = float(pixel_size)
        self._world_shape_yx = (morph_full_shape_yx[0] * self._pixel_size,
                                morph_full_shape_yx[1] * self._pixel_size)
        self._canvas_native = canvas_native

        # Build grayscale QPixmap from DAPI channel (morph_thumb[0])
        dapi = morph_thumb[0].astype(np.float32)
        p99 = np.percentile(dapi, 99)
        if p99 > 0:
            dapi = np.clip(dapi / p99 * 255, 0, 255).astype(np.uint8)
        else:
            dapi = np.zeros_like(dapi, dtype=np.uint8)

        # Scale thumb to fixed widget size
        thumb_h, thumb_w = dapi.shape
        scale = min(self._WIDGET_W / thumb_w, self._WIDGET_H / thumb_h)
        scaled_w = int(thumb_w * scale)
        scaled_h = int(thumb_h * scale)

        # Convert to QPixmap
        img = QImage(dapi.data, thumb_w, thumb_h, thumb_w, QImage.Format_Grayscale8)
        pixmap = QPixmap.fromImage(img).scaled(
            scaled_w, scaled_h, Qt.KeepAspectRatio, Qt.SmoothTransformation
        )
        self._pixmap = pixmap
        self._thumb_display_size = (scaled_w, scaled_h)

        # Widget setup
        self.setFixedSize(self._WIDGET_W, self._WIDGET_H)
        self.setStyleSheet("background: rgba(0, 0, 0, 160);")
        self.setAttribute(Qt.WA_TranslucentBackground, False)
        self.setWindowFlags(Qt.SubWindow)

        # Connect camera events
        viewer.camera.events.center.connect(self._on_camera_changed)
        viewer.camera.events.zoom.connect(self._on_camera_changed)

        # Install event filter on canvas to track resize
        canvas_native.installEventFilter(self)

        self._reposition()

    # ── Camera events ─────────────────────────────────────────────────────

    def _on_camera_changed(self, event=None):
        self.update()

    # ── Qt event handling ─────────────────────────────────────────────────

    def eventFilter(self, obj, event):
        if obj is self._canvas_native and event.type() == QEvent.Resize:
            self._reposition()
        return False

    def _reposition(self):
        cw = self._canvas_native.width()
        self.move(cw - self._WIDGET_W - self._MARGIN, self._MARGIN)
        self.raise_()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor(0, 0, 0, 160))

        # Draw thumbnail centered in widget
        sw, sh = self._thumb_display_size
        ox = (self._WIDGET_W - sw) // 2
        oy = (self._WIDGET_H - sh) // 2
        painter.drawPixmap(ox, oy, self._pixmap)

        # Compute viewport rectangle in data coords
        try:
            center = self._viewer.camera.center  # (z, y, x) in world units
            zoom = self._viewer.camera.zoom      # canvas pixels per world unit
            canvas_h = self._canvas_native.height()
            canvas_w = self._canvas_native.width()

            half_h_world = canvas_h / (2.0 * zoom)
            half_w_world = canvas_w / (2.0 * zoom)

            full_h, full_w = self._world_shape_yx
            scale_x = sw / full_w
            scale_y = sh / full_h

            # Viewport in minimap pixel coords
            cy_world = center[-2] if len(center) >= 2 else center[0]
            cx_world = center[-1] if len(center) >= 1 else center[0]

            rect_y = int((cy_world - half_h_world) * scale_y) + oy
            rect_x = int((cx_world - half_w_world) * scale_x) + ox
            rect_h = int(2 * half_h_world * scale_y)
            rect_w = int(2 * half_w_world * scale_x)

            # Clip to minimap bounds
            rect_x = max(ox, min(rect_x, ox + sw))
            rect_y = max(oy, min(rect_y, oy + sh))
            rect_x2 = max(ox, min(rect_x + rect_w, ox + sw))
            rect_y2 = max(oy, min(rect_y + rect_h, oy + sh))

            pen = QPen(QColor(255, 255, 255), 1)
            painter.setPen(pen)
            painter.drawRect(QRect(rect_x, rect_y, rect_x2 - rect_x, rect_y2 - rect_y))
        except Exception:
            pass

        painter.end()

    def mousePressEvent(self, event):
        """Click to navigate camera to that point in the tissue."""
        sw, sh = self._thumb_display_size
        ox = (self._WIDGET_W - sw) // 2
        oy = (self._WIDGET_H - sh) // 2

        click_x = event.x() - ox
        click_y = event.y() - oy

        if click_x < 0 or click_y < 0 or click_x > sw or click_y > sh:
            return

        full_h, full_w = self._world_shape_yx
        scale_x = sw / full_w
        scale_y = sh / full_h

        world_x = click_x / scale_x
        world_y = click_y / scale_y

        # napari camera.center is (z, y, x), in world units
        center = list(self._viewer.camera.center)
        center[-2] = world_y
        center[-1] = world_x
        self._viewer.camera.center = tuple(center)
