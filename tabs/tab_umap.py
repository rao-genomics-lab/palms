"""Tab 3: UMAP — show UMAP window + point size slider."""

from __future__ import annotations
from typing import TYPE_CHECKING

from magicgui.widgets import PushButton, Slider
from tabs._helpers import make_tab

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    show_umap_button = PushButton(label="Show UMAP Window", enabled=True)
    umap_size_slider = Slider(label="UMAP pt size", min=1, max=50, value=15)

    def on_show_umap():
        ctx.umap_viewer.show()

    def on_umap_size_change(value):
        ctx.umap_viewer.set_point_size(value / 100.0)

    show_umap_button.clicked.connect(on_show_umap)
    umap_size_slider.changed.connect(on_umap_size_change)

    widget = make_tab(show_umap_button, umap_size_slider)
    return widget, {}
