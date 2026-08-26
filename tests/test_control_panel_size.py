"""No control-panel page may pin the size of the Xenium Controls dock.

The panel is a ``QTabWidget`` of ``QTabWidget``s, and a stacked widget's minimum
size is the maximum over *all* its pages — the hidden ones included. So one page
that reports a large minimum becomes the floor for the whole dock, and the
separator between the dock and the canvas stops responding: there is nothing
left to give. ``_helpers.make_tab`` wraps page content in a ``QScrollArea``
precisely to stop that, since a scroll area reports a fixed ~68px minimum
whatever it holds.

Two pages had bypassed the wrapper and reintroduced the failure on the width
axis. Measured before the fix, with the same stubs used here: Notebook 528x107
(its six-button toolbar sat outside its own scroll area), Templates 389x468 (no
scroll area at all, and up to ~609 wide once its fourth button is shown),
combining to a 536x534 panel. Afterwards: 76x109, 68x68, and a 131x175 panel.

The bound below is deliberately loose — ``minimumSizeHint`` depends on the font,
and the point is to catch a page reporting hundreds of pixels, not to pin an
exact number.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

# Well above the ~120px a scroll-wrapped page reports, far below the hundreds
# that a bare button row or splitter costs.
MAX_MINIMUM = 250


@pytest.fixture
def notebook_tab(qapp):
    from xenium_viewer.tabs.tab_notebook import build_tab
    ctx = SimpleNamespace(state={}, viewer=None, data_path=None,
                          set_status=lambda *a, **k: None)
    widget, _ = build_tab(ctx)
    return widget


@pytest.fixture
def templates_tab(qapp):
    from xenium_viewer.tabs.tab_templates import build_tab
    widget, _ = build_tab(SimpleNamespace(state={}))
    return widget


def _hint(widget):
    hint = widget.minimumSizeHint()
    return hint.width(), hint.height()


def test_the_notebook_page_does_not_pin_the_dock(notebook_tab):
    width, height = _hint(notebook_tab)
    assert width <= MAX_MINIMUM, f"Notebook page minimum width {width}px"
    assert height <= MAX_MINIMUM, f"Notebook page minimum height {height}px"


def test_the_templates_page_does_not_pin_the_dock(templates_tab):
    width, height = _hint(templates_tab)
    assert width <= MAX_MINIMUM, f"Templates page minimum width {width}px"
    assert height <= MAX_MINIMUM, f"Templates page minimum height {height}px"


def test_the_templates_page_stays_small_with_every_button_shown(templates_tab):
    """The fourth button appears only while a review is pending.

    It carries the longest label in the tab, so a size measured with it hidden
    would miss the worst case — which is the state the user is actually in when
    they go to the tab to resolve an upgrade.
    """
    from qtpy.QtWidgets import QPushButton

    for button in templates_tab.findChildren(QPushButton):
        button.setVisible(True)

    width, height = _hint(templates_tab)
    assert width <= MAX_MINIMUM, f"Templates page minimum width {width}px"
    assert height <= MAX_MINIMUM, f"Templates page minimum height {height}px"


def test_the_assembled_panel_can_shrink(qapp, notebook_tab, templates_tab):
    """The aggregate is what the dock inherits, so measure it the way app.py nests it."""
    from qtpy.QtWidgets import QTabWidget

    group = QTabWidget()
    group.setTabPosition(QTabWidget.TabPosition.South)
    group.addTab(notebook_tab, "Notebook")
    group.addTab(templates_tab, "Templates")
    panel = QTabWidget()
    panel.addTab(group, "Tools")

    width, height = _hint(panel)
    assert width <= MAX_MINIMUM, f"Control panel minimum width {width}px"
    assert height <= MAX_MINIMUM, f"Control panel minimum height {height}px"


def test_a_page_with_wide_content_reports_a_small_minimum(qapp):
    """The mechanism itself: this is why every page must be wrapped."""
    from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

    from xenium_viewer.tabs._helpers import scrollable

    wide = QWidget()
    layout = QHBoxLayout()
    for _ in range(6):
        layout.addWidget(QPushButton("a button with a fairly long label"))
    wide.setLayout(layout)
    assert wide.minimumSizeHint().width() > 500

    assert scrollable(wide).minimumSizeHint().width() <= MAX_MINIMUM
