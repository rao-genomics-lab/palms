"""ROI expression must record code that computes the numbers it shows.

This is the node the Tier-2 verification flagged: ``roi_expression:A1cf``
recorded two comment lines saying the per-region means were "shown in the
viewer". That cell executes cleanly and does nothing, so ``allow_errors=False``
could never catch it — the notebook "passed" while carrying no ROI analysis at
all, and the CSV export cell said only that a file had been saved.

Both are now real steps, so these tests execute the rendered source and check
the numbers against the same computation done by hand.

Importing the tab pulls in Qt/napari:
    QT_QPA_PLATFORM=offscreen pytest tests/test_roi_expression_step.py
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
sd = pytest.importorskip("spatialdata")
anndata = pytest.importorskip("anndata")
pytest.importorskip("shapely")
pytest.importorskip("qtpy")

from palms.utils.prov_graph import ARTIFACT, TERMINAL  # noqa: E402
from palms.utils.steps import Step, StepExecutor, check_step  # noqa: E402
from palms.tabs.tab_roi import (  # noqa: E402
    _ROIS_TEMPLATE, _roi_expr_template,
)

GENE = "GENE1"

# Two 10×10 squares, in pixel (y, x) order as napari's shapes layer stores them.
POLYGONS = [
    [[0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0]],
    [[20.0, 20.0], [20.0, 30.0], [30.0, 30.0], [30.0, 20.0]],
]


def _adata():
    """Six cells: four in region 1, two in region 2, one outside both."""
    coords_um = np.array([
        [2.0, 2.0], [4.0, 4.0], [6.0, 3.0], [3.0, 7.0],   # region 1
        [22.0, 22.0], [27.0, 26.0],                        # region 2
        [50.0, 50.0],                                      # outside
    ])
    expression = np.array([1.0, 3.0, 5.0, 7.0, 10.0, 20.0, 99.0], dtype="float32")

    adata = anndata.AnnData(np.column_stack([expression, np.zeros(7, "float32")]))
    adata.var_names = [GENE, "GENE2"]
    adata.obs_names = [f"cell{i}" for i in range(7)]
    adata.obs["cell_id"] = [f"cell{i}" for i in range(7)]
    adata.obs["leiden"] = pd.Categorical(["0", "0", "1", "1", "0", "1", "0"])
    adata.obsm["spatial"] = coords_um          # µm, (x, y)
    return adata


def _run(filtered: bool = False, polygons=POLYGONS, adata=None, **extra):
    """Execute the rois + roi_expression steps and return the executor."""
    ex = StepExecutor(namespace={"sc": sc, "sd": sd, "np": np, "pd": pd,
                                 "adata": adata if adata is not None else _adata()})
    ex.run(Step(id="rois", template=_ROIS_TEMPLATE,
                params={"polygons": polygons}, outputs=["roi_polygons"]))
    params = {"gene": GENE, "pixel_size": 1.0}
    if filtered:
        params.update(clustering="leiden", selected=["0"])
    params.update(extra)
    ex.run(Step(
        id=f"roi_expression:{GENE}", template=_roi_expr_template(filtered),
        params=params, deps=["rois"], kind=ARTIFACT,
        outputs=["roi_expr_cells", "roi_expr_stats", "roi_expr_tests"],
    ))
    return ex


# ── the cell is code, and it is self-contained ───────────────────────────────

@pytest.mark.parametrize("filtered", [False, True])
def test_the_template_is_executable_and_self_contained(filtered):
    step = Step(id="roi_expression:X", template=_roi_expr_template(filtered),
                params={"gene": "X", "pixel_size": 0.2125,
                        "clustering": "leiden", "selected": ["0"]})
    ast.parse(step.render())
    missing = check_step(
        step, available={"sc", "sd", "np", "pd", "adata", "roi_polygons"})
    assert missing == set()


def test_the_cell_is_not_a_comment(filtered=False):
    """The exact defect: a node whose body is prose replays as a no-op."""
    step = Step(id="roi_expression:X", template=_roi_expr_template(filtered),
                params={"gene": "X", "pixel_size": 1.0})
    assert ast.parse(step.render()).body


# ── the numbers ──────────────────────────────────────────────────────────────

def test_each_region_gets_the_cells_inside_it():
    cells = _run().get("roi_expr_cells")
    assert list(cells.groupby("region_id").size()) == [4, 2]
    assert set(cells["cell_id"]) == {f"cell{i}" for i in range(6)}   # not cell6


def test_the_statistics_match_the_expression_of_the_cells_inside():
    stats = _run().get("roi_expr_stats")
    region1 = np.array([1.0, 3.0, 5.0, 7.0])

    assert stats.loc[1, "count"] == 4
    assert stats.loc[1, "mean"] == pytest.approx(region1.mean())
    assert stats.loc[1, "median"] == pytest.approx(np.median(region1))
    assert stats.loc[1, "min"] == pytest.approx(1.0)
    assert stats.loc[1, "max"] == pytest.approx(7.0)
    assert stats.loc[2, "mean"] == pytest.approx(15.0)


def test_the_frame_stays_region_grouped_in_cell_order():
    """The per-region frames became one `sc.get.obs_df` call plus a filter.

    That builds the table in *cell* order, so without the stable sort the
    exported CSV silently changed row order for every existing user. Columns,
    order and dtypes are the export's contract.
    """
    cells = _run().get("roi_expr_cells")

    assert list(cells.columns) == ["region_id", "cell_id", "x_centroid_um",
                                   "y_centroid_um", "expression"]
    assert list(cells["region_id"]) == [1, 1, 1, 1, 2, 2]
    assert list(cells["cell_id"]) == [f"cell{i}" for i in range(6)]
    assert list(cells.index) == list(range(6))          # reset, not inherited


def test_the_coordinates_are_the_micron_ones_the_export_promises():
    cells = _run().get("roi_expr_cells")
    first = cells.loc[cells["cell_id"] == "cell0"].iloc[0]
    assert (first["x_centroid_um"], first["y_centroid_um"]) == (2.0, 2.0)


def test_membership_is_what_the_shapely_loop_computed():
    """The migration's whole claim, over a *non-convex* ROI.

    The trap this pins: ``polygon_query`` on the SpatialData as the viewer loads
    it delegates to ``bounding_box_query`` (the table is annotated by a labels
    raster, not by points), which would return the ROI's bounding box and pass
    every convex-ROI test silently. An L-shape distinguishes the two.
    """
    from shapely import contains_xy
    from shapely.geometry import Polygon

    rng = np.random.default_rng(0)
    coords = rng.uniform(0, 40, size=(400, 2))
    adata = anndata.AnnData(rng.random((400, 2), dtype="float32"))
    adata.var_names = [GENE, "GENE2"]
    adata.obs_names = [f"cell{i}" for i in range(400)]
    adata.obs["cell_id"] = list(adata.obs_names)
    adata.obsm["spatial"] = coords

    # An L, in napari (y, x) order. Its bounding box is the full 30x30 square.
    l_shape_yx = [[0.0, 0.0], [0.0, 30.0], [10.0, 30.0],
                  [10.0, 10.0], [30.0, 10.0], [30.0, 0.0]]

    cells = _run(polygons=[l_shape_yx], adata=adata).get("roi_expr_cells")

    poly_xy = Polygon(np.asarray(l_shape_yx)[:, ::-1])
    legacy = contains_xy(poly_xy, coords[:, 0], coords[:, 1])

    assert set(cells["cell_id"]) == {f"cell{i}" for i in np.flatnonzero(legacy)}
    # ...and the bounding box would have been a strictly larger answer, so the
    # equality above is load-bearing rather than trivially satisfied.
    assert legacy.sum() < ((coords < 30).all(axis=1)).sum()


def test_an_roi_holding_no_cells_is_handled():
    """``polygon_query`` returns ``None``, not an empty frame.

    The boolean mask it replaced had no equivalent case, so nothing in the old
    template exercised this path.
    """
    import spatialdata as _sd
    from shapely.geometry import Polygon

    adata = _adata()
    points = _sd.models.PointsModel.parse(
        pd.DataFrame({"x": adata.obsm["spatial"][:, 0],
                      "y": adata.obsm["spatial"][:, 1],
                      "cell_index": np.arange(adata.n_obs)}),
        coordinates={"x": "x", "y": "y"},
    )
    far = Polygon([[900, 900], [900, 910], [910, 910], [910, 900]])
    assert _sd.polygon_query(points, far, target_coordinate_system="global") is None


def test_a_self_intersecting_roi_keeps_both_lobes():
    """``buffer(0)`` silently deleted one lobe; ``make_valid`` repairs the shape.

    A figure-eight is the shape a user draws by crossing the polygon tool over
    its own edge, and it used to lose half its cells with no warning.
    """
    from shapely import make_valid
    from shapely.geometry import Polygon

    adata = _adata()
    # One cell in each lobe of a bow-tie spanning (0,0)-(10,10).
    adata.obsm["spatial"] = np.array([
        [5.0, 8.0], [5.0, 2.0],
        [50.0, 50.0], [51.0, 51.0], [52.0, 52.0], [53.0, 53.0], [54.0, 54.0],
    ])
    figure_eight = [[0.0, 0.0], [10.0, 10.0], [10.0, 0.0], [0.0, 10.0]]

    cells = _run(polygons=[figure_eight], adata=adata).get("roi_expr_cells")
    assert set(cells["cell_id"]) == {"cell0", "cell1"}

    # The regression itself, stated directly.
    bow_tie = Polygon(np.asarray(figure_eight)[:, ::-1])
    assert make_valid(bow_tie).area == pytest.approx(2 * bow_tie.buffer(0).area)


def test_an_empty_region_is_reported_rather_than_dropped():
    """The tab prints 'Region N: 0 cells'; the frame has to carry that row."""
    empty_far_away = [[900.0, 900.0], [900.0, 910.0], [910.0, 910.0], [910.0, 900.0]]
    stats = _run(polygons=POLYGONS + [empty_far_away]).get("roi_expr_stats")

    assert list(stats.index) == [1, 2, 3]
    assert stats.loc[3, "count"] == 0
    assert np.isnan(stats.loc[3, "mean"])


def test_the_cluster_filter_keeps_only_the_selected_clusters():
    cells = _run(filtered=True).get("roi_expr_cells")
    assert set(cells["cell_id"]) == {"cell0", "cell1", "cell4"}


# ── the significance testing the tab prints ──────────────────────────────────

def test_the_pairwise_test_is_welchs_t_test():
    from scipy import stats as scipy_stats

    tests = _run().get("roi_expr_tests")
    expected = scipy_stats.ttest_ind(np.array([1.0, 3.0, 5.0, 7.0]),
                                     np.array([10.0, 20.0]), equal_var=False)

    assert len(tests) == 1
    assert tests.loc[0, "t"] == pytest.approx(expected.statistic)
    assert tests.loc[0, "p"] == pytest.approx(expected.pvalue)


def test_a_single_comparison_is_not_multiplicity_corrected():
    tests = _run().get("roi_expr_tests")
    assert tests.loc[0, "p_adj"] == pytest.approx(tests.loc[0, "p"])


def test_several_comparisons_are_benjamini_hochberg_corrected():
    from scipy import stats as scipy_stats

    third = [[20.0, 0.0], [20.0, 10.0], [30.0, 10.0], [30.0, 0.0]]

    # A region holding a single cell has no variance, so it is not testable and
    # must not appear as a comparison.
    lonely = _adata()
    lonely.obsm["spatial"] = np.vstack([lonely.obsm["spatial"][:6], [[5.0, 25.0]]])
    tests = _run(polygons=POLYGONS + [third], adata=lonely).get("roi_expr_tests")
    assert len(tests) == 1

    # With three populated regions there are three comparisons.
    adata = _adata()
    adata.obsm["spatial"] = np.array([
        [2.0, 2.0], [4.0, 4.0],
        [22.0, 22.0], [27.0, 26.0],
        [5.0, 22.0], [8.0, 27.0],
        [50.0, 50.0],
    ])
    tests = _run(polygons=POLYGONS + [third], adata=adata).get("roi_expr_tests")

    assert len(tests) == 3
    expected = scipy_stats.false_discovery_control(tests["p"], method="bh")
    assert np.allclose(tests["p_adj"], expected)
    assert (tests["p_adj"] >= tests["p"] - 1e-12).all()


# ── the export writes the file it records ────────────────────────────────────

def test_the_export_cell_writes_the_csv(tmp_path):
    ex = _run()
    out = tmp_path / "roi_GENE1.csv"
    ex.run(Step(
        id="export:roi_expression",
        template="\n# Export ROI per-cell expression of $gene\n"
                 "roi_expr_cells.to_csv($path, index=False)",
        params={"gene": GENE, "path": str(out)},
        deps=[f"roi_expression:{GENE}"], kind=TERMINAL,
    ))

    written = pd.read_csv(out)
    assert list(written.columns) == ["region_id", "cell_id", "x_centroid_um",
                                     "y_centroid_um", "expression"]
    assert len(written) == 6


def test_the_tab_no_longer_records_prose_for_either_node():
    """Source guard: both nodes were comment-only, and both must stay code."""
    text = (Path(__file__).resolve().parent.parent / "src" / "palms"
            / "tabs" / "tab_roi.py").read_text()
    for prose in ("mean expression is shown in the viewer",
                  "saved from the viewer to"):
        assert prose not in text, prose


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
