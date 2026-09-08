"""A launch builds the whole UI, and says so if anything went wrong on the way.

This is the coverage gap the suite exists for: every tab is unit-tested against a
``SimpleNamespace`` stub, and none had ever been built against a real
``ViewerContext`` — let alone alongside the other 28, which is the only
configuration a user ever sees.
"""
from __future__ import annotations

import pytest

# app.py's addTab calls are the authoritative count: 4 Cells + 4 Genes +
# 7 Spatial + 4 Images + 10 Tools.
TAB_GROUPS = {
    "Cells": ["Clustering", "Coloring", "Transcripts", "UMAP"],
    "Genes": ["Rank Genes", "Markers", "Correlation", "CNV"],
    "Spatial": ["ROI DEG", "Lig-Rec", "Nhood Enrich", "Co-occur", "Domains",
                "Annot Nhood", "Annot Dist"],
    "Images": ["H&E", "ARMS", "Ext Images", "Patches"],
    "Tools": ["Annotations", "QC", "Preprocess", "Segmentation", "Crop Dataset",
              "Publish", "Notebook", "Dataset", "Cache", "Templates"],
}

EXPECTED_LAYERS = {
    "cell_labels", "nucleus_labels", "morphology_focus",
    "ROIs", "Annotations", "Crop Regions", "transcripts",
}


def test_the_layers_the_app_adds_are_all_there(launched):
    names = {layer.name for layer in launched.viewer.layers}
    assert EXPECTED_LAYERS <= names, f"missing: {sorted(EXPECTED_LAYERS - names)}"
    assert launched.viewer.layers["cell_labels"].visible
    assert launched.ctx.cell_labels_layer is not None


def test_every_tab_is_built(launched):
    """All 29 pages, in five groups, from one real context.

    Not a size check (``test_control_panel_size.py`` does that) — a *build*
    check. A tab whose ``build_tab`` raises on real data would take the whole
    dock down with it, and nothing before this test could see that.
    """
    panel = launched.app["dock_widget"].widget()
    groups = [panel.tabText(i) for i in range(panel.count())]
    assert groups == list(TAB_GROUPS)

    for i, group in enumerate(groups):
        sub = panel.widget(i)
        assert [sub.tabText(j) for j in range(sub.count())] == TAB_GROUPS[group]
        for j in range(sub.count()):
            assert sub.widget(j) is not None

    assert sum(len(v) for v in TAB_GROUPS.values()) == 29


def test_the_plots_dock_exists_and_starts_hidden(launched):
    dock = launched.app["plots_dock"]
    assert dock is not None
    assert not dock.isVisible()
    assert launched.ctx.plots_panel is not None


def test_navigating_to_a_tab_by_name_finds_its_widgets(launched):
    """The rig can reach a tab's real widgets — the premise of every other test."""
    launched.rig.tab("Cells", "Clustering")
    assert launched.rig.mg("Resolution") is not None
    launched.rig.tab("Tools", "QC")
    assert launched.rig.page is not None


def test_nothing_reported_a_failure_during_startup(launched):
    """The app swallows these on purpose; a test must not.

    ``reporting`` collects what would otherwise only be a napari toast: failed
    writes, steps that could not be recorded, customised templates that were
    rejected. ``_load_dataset`` resets the registry at the start of each launch,
    and every test in this module is read-only, so what is left here is what the
    launch itself produced.
    """
    from palms.utils import reporting
    assert reporting.failures() == []
    assert reporting.failure_summary() == "No write failures this session."


def test_the_restored_session_is_the_one_on_disk(launched, e2e_source):
    """The fixture's own session comes back — including both QC lineages."""
    if e2e_source.name != "crop_7":
        pytest.skip("assertions below are about the committed crop_7 fixture")

    assert launched.ctx.adata.n_obs == 151
    assert launched.ctx.pixel_size == pytest.approx(0.2125)

    graph = launched.graph
    assert len(graph) == 31
    for node_id in ("preamble", "normalize", "normalize:qc", "qc_filter",
                    "spatial_neighbors", "spatial_neighbors:qc",
                    "clustering:leiden_igraph_r1.0", "cnv:infercnv"):
        assert node_id in graph, node_id

    # The barrier is what keeps the exported notebook running the unfiltered
    # clusterings on unfiltered cells; it has to survive a round trip.
    assert graph.get("qc_filter").barrier is True

    # Reverted before the fixture was saved, so the node is present with no
    # filter in force — the state `clear_qc_filter` leaves when something still
    # depends on the node.
    assert launched.ctx.state.get("qc_filter") is None
    assert launched.ctx.cell_root() == "preamble"
