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

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib                                                # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

# Same loose bound as tests/test_control_panel_size.py: well above the ~120px a
# scroll-wrapped widget reports, far below the hundreds a bare button row costs.
MAX_MINIMUM = 250


def _figure():
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


# ── not everything that looks like a Figure is one ───────────────────────────

def test_a_scanpy_plot_object_is_resolved_to_its_figure():
    """The bug this section exists for.

    ``sc.pl.rank_genes_groups_dotplot(return_fig=True)`` returns a ``DotPlot``,
    which has ``savefig`` — so the old save-only path never noticed — but no
    canvas, so drawing it raised ``AttributeError: 'DotPlot' object has no
    attribute 'set_canvas'`` the moment the Plots dock tried to make a
    thumbnail. Reported from a real session on the dotplot button.

    ``DotPlot.fig`` is ``None`` until the plot is built, so reading the
    attribute is not enough either — hence ``get_axes()``.
    """
    sc = pytest.importorskip("scanpy")
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    import warnings
    from matplotlib.figure import Figure

    from xenium_viewer.utils.fig_render import to_figure

    n = 60
    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(3, (n, 10)).astype("float32"))
    adata.obs_names = [f"c{i}" for i in range(n)]
    adata.var_names = [f"Gene{i}" for i in range(10)]
    adata.obs["leiden"] = pd.Categorical([str(i % 3) for i in range(n)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        dotplot = sc.pl.rank_genes_groups_dotplot(
            adata, groupby="leiden", n_genes=2, show=False, return_fig=True)

    assert not isinstance(dotplot, Figure), (
        "fixture no longer reproduces the case — scanpy now returns a Figure"
    )
    assert getattr(dotplot, "fig", None) is None, (
        "fixture no longer reproduces the case — the figure is built eagerly"
    )
    assert isinstance(to_figure(dotplot), Figure)
    plt.close("all")


def test_a_scanpy_plot_object_can_be_shown_and_saved(panel, tmp_path):
    """End to end: the whole show_plot path, on the object that broke it."""
    sc = pytest.importorskip("scanpy")
    ad = pytest.importorskip("anndata")
    pd = pytest.importorskip("pandas")
    import types
    import warnings

    from xenium_viewer.tabs._helpers import create_shared_helpers
    from xenium_viewer.utils.viewer_context import ViewerContext

    n = 60
    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(3, (n, 10)).astype("float32"))
    adata.obs_names = [f"c{i}" for i in range(n)]
    adata.var_names = [f"Gene{i}" for i in range(10)]
    adata.obs["leiden"] = pd.Categorical([str(i % 3) for i in range(n)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        dotplot = sc.pl.rank_genes_groups_dotplot(
            adata, groupby="leiden", n_genes=2, show=False, return_fig=True)

    ctx = ViewerContext(data_path=tmp_path, state={},
                        viewer=types.SimpleNamespace(status=""))
    ctx.plots_panel = panel
    create_shared_helpers(ctx)

    paths = ctx.show_plot(dotplot, "dotplot_leiden", title="Dotplot")
    assert len(paths) == 2 and all(Path(p).stat().st_size > 0 for p in paths)
    assert panel.count() == 1
    plt.close("all")


def test_a_thing_that_is_not_a_figure_at_all_says_so():
    """Better a named TypeError here than an AttributeError inside Qt."""
    from xenium_viewer.utils.fig_render import to_figure

    with pytest.raises(TypeError, match="expected a matplotlib Figure"):
        to_figure(object())
