"""The Plots dock — every figure the viewer makes, in one place.

Issue #35: a dotplot appeared in a floating matplotlib window, a CNV heatmap
appeared nowhere at all (its button said "Show Heatmap" and only wrote a file),
and an annotation-distance plot opened a *blocking* window. Four display modes,
none of them discoverable. This panel is the one answer: a figure is added here,
it is on screen, and it says where it was written.

Qt-only. It never imports pyplot except through
:mod:`palms.utils.fig_render`, so nothing here depends on which backend
the process happens to be using.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget,
)

from palms.utils.fig_render import (
    detach_from_pyplot, fig_to_pixmap, open_figure_window,
)

#: Thumbnails are sized for a bottom dock, which is wide and short.
THUMBNAIL_WIDTH = 280

#: How many figures the gallery keeps. Beyond this the oldest is dropped and its
#: figure closed — a long session otherwise holds every figure it ever drew, and
#: a rank-genes heatmap is not small.
MAX_ENTRIES = 20


@dataclass
class PlotEntry:
    """One figure in the gallery."""
    fig: Any
    title: str
    paths: list = field(default_factory=list)
    card: Any = None


class PlotsPanel(QWidget):
    """Scrollable gallery of the session's figures, newest first."""

    def __init__(self, parent=None):
        super().__init__(parent)
        # Imported here rather than at module scope: ``tabs._helpers`` imports
        # from ``utils``, so a top-level import the other way round would close
        # a cycle at interpreter start.
        from palms.tabs._helpers import scrollable, toolbar_row

        self._entries: list[PlotEntry] = []
        self._windows: list[QWidget] = []   # open figure windows, kept alive

        self._cards = QWidget()
        self._cards_layout = QVBoxLayout()
        self._cards_layout.setContentsMargins(4, 4, 4, 4)
        self._cards_layout.setSpacing(6)
        self._cards_layout.addStretch()
        self._cards.setLayout(self._cards_layout)

        self._empty = QLabel("No plots yet — figures from the analysis tabs "
                             "appear here.")
        self._empty.setWordWrap(True)
        self._empty.setStyleSheet("color: gray;")
        self._cards_layout.insertWidget(0, self._empty)

        clear_button = QPushButton("Clear all")
        clear_button.clicked.connect(self.clear)

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        # A button row outside a scroll area cannot shrink below its labels and
        # would set the dock's minimum width all by itself — the defect that
        # pinned the Controls dock at 536px. ``toolbar_row`` keeps it
        # pinned but horizontally scrollable.
        layout.addWidget(toolbar_row(clear_button))
        layout.addWidget(scrollable(self._cards))
        self.setLayout(layout)

    # ── public API ───────────────────────────────────────────────────────────

    def add_figure(self, fig, title: str, paths=None) -> PlotEntry:
        """Add *fig* to the top of the gallery and return its entry."""
        entry = PlotEntry(fig=fig, title=title, paths=list(paths or []))
        entry.card = self._build_card(entry)
        self._entries.insert(0, entry)
        self._cards_layout.insertWidget(0, entry.card)
        self._empty.setVisible(False)
        self._evict()
        # Hand the figure to us alone: pyplot's registry would otherwise keep
        # every figure of the session alive behind our back.
        detach_from_pyplot(fig)
        return entry

    def clear(self) -> None:
        for entry in list(self._entries):
            self._remove(entry)

    def count(self) -> int:
        return len(self._entries)

    @property
    def entries(self) -> list[PlotEntry]:
        return list(self._entries)

    # ── card construction ────────────────────────────────────────────────────

    def _build_card(self, entry: PlotEntry) -> QWidget:
        card = QFrame()
        card.setFrameShape(QFrame.Shape.StyledPanel)
        row = QHBoxLayout()
        row.setContentsMargins(6, 6, 6, 6)

        thumb = QLabel()
        thumb.setPixmap(fig_to_pixmap(entry.fig, THUMBNAIL_WIDTH))
        thumb.setAlignment(Qt.AlignmentFlag.AlignTop)
        thumb.setToolTip("Click Open for a full-size, pannable view")
        row.addWidget(thumb)

        side = QVBoxLayout()
        name = QLabel(f"<b>{entry.title}</b>")
        name.setWordWrap(True)
        side.addWidget(name)

        if entry.paths:
            where = QLabel("\n".join(entry.paths))
            where.setWordWrap(True)
            where.setStyleSheet("color: gray; font-size: 10px;")
            where.setTextInteractionFlags(
                Qt.TextInteractionFlag.TextSelectableByMouse)
            side.addWidget(where)

        buttons = QHBoxLayout()
        open_button = QPushButton("Open")
        open_button.clicked.connect(lambda: self._open(entry))
        save_button = QPushButton("Save as…")
        save_button.clicked.connect(lambda: self._save_as(entry))
        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda: self._remove(entry))
        for button in (open_button, save_button, remove_button):
            buttons.addWidget(button)
        buttons.addStretch()
        side.addLayout(buttons)
        side.addStretch()
        row.addLayout(side)
        card.setLayout(row)
        return card

    # ── actions ──────────────────────────────────────────────────────────────

    def _open(self, entry: PlotEntry) -> None:
        window = open_figure_window(entry.fig, entry.title, parent=None)
        self._windows.append(window)
        window.destroyed.connect(
            lambda *_: self._windows.remove(window)
            if window in self._windows else None)

    def _save_as(self, entry: PlotEntry) -> None:
        from qtpy.QtWidgets import QFileDialog
        from palms.utils.plot_output import safe_stem
        suggestion = entry.paths[0] if entry.paths else f"{safe_stem(entry.title)}.png"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot As", suggestion,
            "Images (*.png *.pdf *.svg);;All Files (*)",
        )
        if path:
            entry.fig.savefig(path, dpi=300, bbox_inches="tight")

    def _remove(self, entry: PlotEntry) -> None:
        if entry not in self._entries:
            return
        self._entries.remove(entry)
        if entry.card is not None:
            self._cards_layout.removeWidget(entry.card)
            entry.card.setParent(None)
            entry.card.deleteLater()
        detach_from_pyplot(entry.fig)
        if not self._entries:
            self._empty.setVisible(True)

    def _evict(self) -> None:
        while len(self._entries) > MAX_ENTRIES:
            self._remove(self._entries[-1])
