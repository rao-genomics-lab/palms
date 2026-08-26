"""The Plots dock: it must not pin the window, and it must not hoard figures.

Two properties, each aimed at a defect this codebase has already had once.

The size one is the March/August dock-resize bug, which recurred because a page
was added without a scroll area: a widget in a dock that reports a large minimum
takes that space away from the canvas permanently. A second dock has exactly the
same failure mode, and a gallery of 280px thumbnails plus a three-button row is
precisely the shape that causes it.

The eviction one is new but predictable: the viewer produces a figure per user
action, matplotlib figures are not small, and a session lasts hours.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

# Same loose bound as tests/test_control_panel_size.py: well above the ~120px a
# scroll-wrapped widget reports, far below the hundreds a bare button row costs.
MAX_MINIMUM = 250


def _figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1], [1, 0])
    return fig


@pytest.fixture
def panel(qapp):
    from xenium_viewer.utils.plots_panel import PlotsPanel
    return PlotsPanel()


def test_an_empty_panel_does_not_pin_the_dock(panel):
    hint = panel.minimumSizeHint()
    assert hint.width() <= MAX_MINIMUM, f"minimum width {hint.width()}px"
    assert hint.height() <= MAX_MINIMUM, f"minimum height {hint.height()}px"


def test_a_populated_panel_does_not_pin_the_dock_either(panel):
    """The thumbnails are what would do it — 280px each, stacked."""
    for i in range(3):
        panel.add_figure(_figure(), f"plot {i}", ["plots/x.png"])
    hint = panel.minimumSizeHint()
    assert hint.width() <= MAX_MINIMUM, f"minimum width {hint.width()}px"
    assert hint.height() <= MAX_MINIMUM, f"minimum height {hint.height()}px"


def test_figures_arrive_newest_first(panel):
    panel.add_figure(_figure(), "first", [])
    panel.add_figure(_figure(), "second", [])
    assert [e.title for e in panel.entries] == ["second", "first"]


def test_the_gallery_is_capped_and_drops_the_oldest(panel):
    from xenium_viewer.utils.plots_panel import MAX_ENTRIES

    for i in range(MAX_ENTRIES + 3):
        panel.add_figure(_figure(), f"plot {i}", [])

    assert panel.count() == MAX_ENTRIES
    titles = [e.title for e in panel.entries]
    assert titles[0] == f"plot {MAX_ENTRIES + 2}"       # newest kept
    assert "plot 0" not in titles                       # oldest dropped


def test_adding_a_figure_hands_it_over_from_pyplot(panel):
    """pyplot keeps every figure it made until closed. The panel holds the one
    reference that matters, so a long session does not accumulate the rest."""
    import matplotlib.pyplot as plt

    plt.close("all")
    fig = _figure()
    assert plt.get_fignums()                    # pyplot is holding it
    panel.add_figure(fig, "handed over", [])
    assert not plt.get_fignums(), "pyplot still holds the figure"
    # …and the Figure itself is still usable — that is the whole trick.
    assert fig.axes


def test_clearing_empties_the_gallery(panel):
    for i in range(3):
        panel.add_figure(_figure(), f"plot {i}", [])
    panel.clear()
    assert panel.count() == 0
