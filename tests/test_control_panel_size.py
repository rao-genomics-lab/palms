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
    from palms.tabs.tab_notebook import build_tab
    ctx = SimpleNamespace(state={}, viewer=None, data_path=None,
                          set_status=lambda *a, **k: None)
    widget, _ = build_tab(ctx)
    return widget


@pytest.fixture
def templates_tab(qapp):
    from palms.tabs.tab_templates import build_tab
    widget, _ = build_tab(SimpleNamespace(state={}))
    return widget


@pytest.fixture
def umap_tab(qapp):
    """The UMAP page, rebuilt for issue #34.

    Worth measuring: it gained a ``QListWidget`` and five buttons, which is the
    shape that pinned the dock twice before.
    """
    from palms.tabs.tab_umap import build_tab
    ctx = SimpleNamespace(
        state={}, viewer=None, adata=None, data_path=None,
        gene_names=[f"Gene{i}" for i in range(20)],
        clusterings={}, clustering_widget=None,
        get_labels_for=lambda key: {},
        record_node=lambda *a, **k: None,
        record_clustering=lambda *a, **k: None,
        umap_viewer=None, dataset_generation=0,
        plot_paths=lambda stem: [f"plots/{stem}.png"],
        apply_plot_font_size=lambda: None,
    )
    widget, _ = build_tab(ctx)
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


def test_the_umap_page_does_not_pin_the_dock(umap_tab):
    width, height = _hint(umap_tab)
    assert width <= MAX_MINIMUM, f"UMAP page minimum width {width}px"
    assert height <= MAX_MINIMUM, f"UMAP page minimum height {height}px"


def test_the_assembled_panel_can_shrink(qapp, notebook_tab, templates_tab, umap_tab):
    """The aggregate is what the dock inherits, so measure it the way app.py nests it."""
    from qtpy.QtWidgets import QTabWidget

    group = QTabWidget()
    group.setTabPosition(QTabWidget.TabPosition.South)
    group.addTab(notebook_tab, "Notebook")
    group.addTab(templates_tab, "Templates")
    panel = QTabWidget()
    panel.addTab(group, "Tools")
    cells = QTabWidget()
    cells.setTabPosition(QTabWidget.TabPosition.South)
    cells.addTab(umap_tab, "UMAP")
    panel.addTab(cells, "Cells")

    width, height = _hint(panel)
    assert width <= MAX_MINIMUM, f"Control panel minimum width {width}px"
    assert height <= MAX_MINIMUM, f"Control panel minimum height {height}px"


def test_make_tab_renders_a_widgets_label(qapp):
    """The caption must survive the trip through ``make_tab``.

    ``.native`` is the bare control; magicgui keeps the caption in a wrapper only
    a ``Container`` builds, so ``layout.addWidget(w.native)`` dropped every one.
    The whole control panel rendered as anonymous sliders — Clustering showed
    ``15``, ``40``, ``1.00``, ``2``, ``2000`` with nothing to say which was
    which — and no test noticed, because none of them rendered a tab. Found by
    screenshotting the running viewer over the ``--mcp`` bridge.
    """
    from qtpy.QtWidgets import QLabel
    from magicgui.widgets import Slider

    from palms.tabs._helpers import make_tab

    page = make_tab(Slider(label="probe_label", min=0, max=10, value=5))
    texts = [w.text() for w in page.findChildren(QLabel)]
    assert "probe_label" in texts, f"label dropped by make_tab; found {texts}"


def test_make_tab_does_not_label_a_widget_that_paints_its_own_text(qapp):
    """A ``CheckBox``/``PushButton`` carries its text *on* the control.

    magicgui's own ``Container`` skips ``ButtonWidget`` for this reason, and
    ``labelled()`` mirrors it. Pinned because the failure is silent and ugly
    rather than loud: every checkbox in the panel would read its text twice.
    """
    from qtpy.QtWidgets import QCheckBox, QLabel

    from magicgui.widgets import CheckBox

    from palms.tabs._helpers import make_tab

    page = make_tab(CheckBox(label="Use HVGs only", value=False))
    on_control = [w.text() for w in page.findChildren(QCheckBox)]
    captions = [w.text() for w in page.findChildren(QLabel) if w.text()]
    assert on_control == ["Use HVGs only"]
    assert captions == [], f"checkbox text duplicated as a caption: {captions}"


def test_make_tab_invents_no_label_for_an_unlabelled_widget(qapp):
    """An unlabelled widget must stay unlabelled.

    magicgui reports ``label == ''`` when ``label=`` was never passed — it does
    not derive one from a variable name — so wrapping unconditionally is safe.
    Asserted rather than assumed: if that ever changed, every bare control in
    the panel would sprout a caption nobody wrote.
    """
    from qtpy.QtWidgets import QLabel
    from magicgui.widgets import Slider

    from palms.tabs._helpers import make_tab

    page = make_tab(Slider(min=0, max=10, value=5))
    captions = [w.text() for w in page.findChildren(QLabel) if w.text()]
    assert captions == [], f"invented a caption: {captions}"


def test_a_page_with_wide_content_reports_a_small_minimum(qapp):
    """The mechanism itself: this is why every page must be wrapped."""
    from qtpy.QtWidgets import QHBoxLayout, QPushButton, QWidget

    from palms.tabs._helpers import scrollable

    wide = QWidget()
    layout = QHBoxLayout()
    for _ in range(6):
        layout.addWidget(QPushButton("a button with a fairly long label"))
    wide.setLayout(layout)
    assert wide.minimumSizeHint().width() > 500

    assert scrollable(wide).minimumSizeHint().width() <= MAX_MINIMUM
