"""The annotation tabs' migration to ``run_step``.

Both tabs produced results and recorded a ``NOTE`` — a comment cell that
replays as a silent no-op — on the stated grounds that the notebook cannot
reach a napari shapes layer. It does not have to: ``annot.polygons`` inlines
the drawn shapes as literals, the way ``roi.polygons`` has for ROIs since the
ROI migration. These tests execute the real templates against a synthetic
AnnData and then replay the recorded graph in a clean namespace, which is the
reproducibility claim reduced to something CI can assert.

Importing a tab module pulls in Qt/napari, so run headless:
    QT_QPA_PLATFORM=offscreen pytest tests/test_annotation_steps.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
sq = pytest.importorskip("squidpy")
sd = pytest.importorskip("spatialdata")
gpd = pytest.importorskip("geopandas")
pytest.importorskip("shapely")
pytest.importorskip("seaborn")
plt = pytest.importorskip("matplotlib.pyplot")

from palms.tabs._helpers import _NORMALIZE_TEMPLATE  # noqa: E402
from palms.utils.step_templates import (  # noqa: E402
    builtin_assemble, builtin_spec, builtin_text,
)
from palms.utils.steps import Step, StepExecutor, check_step  # noqa: E402

CLUSTER_KEY = "leiden_r1.0"
PIXEL_SIZE = 0.5

#: A square and a bowtie, in napari's (y, x) pixels. The bowtie is the shape
#: that told us buffer(0) was wrong; it is here so the template keeps both lobes.
SQUARE_PX = [[0, 0], [0, 100], [100, 100], [100, 0]]        # 50x50 um
BOWTIE_PX = [[120, 120], [160, 160], [160, 120], [120, 160]]


def _adata(n_obs: int = 200, n_vars: int = 20):
    rng = np.random.default_rng(0)
    a = anndata.AnnData(rng.poisson(3, size=(n_obs, n_vars)).astype("float32"))
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_vars)]
    a.obsm["spatial"] = rng.uniform(0, 100, size=(n_obs, 2))
    half = n_obs // 2
    a.obs[CLUSTER_KEY] = pd.Categorical(["0"] * half + ["1"] * (n_obs - half))
    return a


def _ns(adata=None):
    """Mirrors EXECUTOR_BASE_NAMES, which is what the notebook preamble binds."""
    return {"sc": sc, "sq": sq, "sd": sd, "pd": pd, "np": np, "plt": plt,
            "Path": Path, "adata": adata if adata is not None else _adata()}


def _blocks(template_id):
    return list(builtin_spec(template_id).blocks)


def _polygons_step(types=("tumour", "stroma")):
    return Step(
        id="annotations", template=builtin_text("annot.polygons"),
        params={"polygons": [SQUARE_PX, BOWTIE_PX], "types": list(types),
                "pixel_size": PIXEL_SIZE},
        kind="setup", outputs=["annotations"],
    )


def _virtual_cells_step(types=("tumour", "stroma"), density_um2=25.0):
    return Step(
        id="annot_virtual_cells", template=builtin_text("annot.virtual_cells"),
        params={"types": list(types), "density_um2": density_um2},
        deps=["annotations"], outputs=["annot_virtual_cells"],
    )


def _nhood_step(n_perms=100):
    return Step(
        id=f"annot_nhood:{CLUSTER_KEY}", template=builtin_text("annot.nhood"),
        params={"cluster_key": CLUSTER_KEY,
                "uns_key": f"{CLUSTER_KEY}_nhood_enrichment",
                "n_neighs": 6, "n_perms": n_perms, "seed": 42},
        deps=["normalize", "annot_virtual_cells"], outputs=["adata_annot"],
    )


def _distance_step(annotation_type="tumour"):
    return Step(
        id=f"annot_distance:{annotation_type}",
        template=builtin_text("annot.distance"),
        params={"annotation_type": annotation_type,
                "obs_key": f"dist_to_{annotation_type}_um"},
        deps=["annotations"], outputs=["annot_distances"],
    )


def _normalize_step():
    return Step(id="normalize", template=_NORMALIZE_TEMPLATE, kind="setup",
                outputs=["adata_norm"])


def _run(steps, adata=None):
    ex = StepExecutor(namespace=_ns(adata))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for step in steps:
            ex.run(step)
    return ex


def _replay(ex, adata=None):
    """Execute every recorded cell in dependency order, as the notebook does."""
    ns = _ns(adata)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for node_id in ex.graph.topo_sort():
            exec(compile(ex.graph.get(node_id).code, "<replay>", "exec"), ns)  # noqa: S102
    return ns


# ── the polygons step ────────────────────────────────────────────────────────

def test_every_annotation_template_is_self_contained():
    """Rule (a): a template may only reach for EXECUTOR_BASE_NAMES plus the
    names its declared dependencies bind. A name in neither replays as a
    NameError in a clean kernel, which is where nobody is watching."""
    base = set(_ns())
    for step, bound in (
        (_polygons_step(), set()),
        (_virtual_cells_step(), {"annotations"}),
        (_nhood_step(), {"adata_norm", "annot_virtual_cells"}),
        (_distance_step(), {"annotations"}),
    ):
        assert check_step(step, base | bound) == set(), step.id


def test_the_drawn_shapes_become_one_geometry_per_type():
    ex = _run([_polygons_step()])
    annotations = ex.ns["annotations"]

    assert list(annotations.index) == ["stroma", "tumour"]
    # 100x100 px at 0.5 um/px is 50x50 um.
    assert annotations.loc["tumour", "geometry"].area == pytest.approx(2500.0)


def test_a_self_intersecting_annotation_keeps_both_lobes():
    """The buffer(0) defect, at the template level: the bowtie's two 100 um^2
    lobes must both survive the repair, or half the region a user drew is
    silently dropped before anything is sampled from it."""
    ex = _run([_polygons_step()])

    assert ex.ns["annotations"].loc["stroma", "geometry"].area == pytest.approx(200.0)


def test_pixels_are_converted_once():
    """The whole pipeline works in microns because this step converts, and
    nothing downstream converts again — a second scale would be invisible in
    every result except as wrong numbers."""
    ex = _run([_polygons_step()])
    minx, miny, maxx, maxy = ex.ns["annotations"].total_bounds

    assert (minx, miny) == pytest.approx((0.0, 0.0))
    assert (maxx, maxy) == pytest.approx((80.0, 80.0))   # 160 px x 0.5 um/px


# ── virtual cells ────────────────────────────────────────────────────────────

def test_virtual_cells_land_inside_their_annotation_and_are_labelled():
    ex = _run([_polygons_step(), _virtual_cells_step()])
    cells = ex.ns["annot_virtual_cells"]

    assert set(cells["annotation_type"]) == {"tumour", "stroma"}
    tumour = cells[cells["annotation_type"] == "tumour"]
    assert tumour["x_um"].between(0, 50).all() and tumour["y_um"].between(0, 50).all()


def test_grid_density_sets_the_number_of_virtual_cells():
    """One cell per density_um2, so the 2500 um^2 square gets ~100 at 25 um^2."""
    ex = _run([_polygons_step(), _virtual_cells_step(density_um2=25.0)])
    cells = ex.ns["annot_virtual_cells"]
    n_tumour = (cells["annotation_type"] == "tumour").sum()

    assert 90 <= n_tumour <= 110

    finer = _run([_polygons_step(), _virtual_cells_step(density_um2=100.0)])
    assert (finer.ns["annot_virtual_cells"]["annotation_type"] == "tumour").sum() < n_tumour


def test_sampling_is_deterministic_so_it_needs_no_seed():
    """The tracker carried "plus a seed so the sampling replays" as remaining
    work. There is nothing to seed: the grid is a lattice over the bounds."""
    first = _run([_polygons_step(), _virtual_cells_step()])
    second = _run([_polygons_step(), _virtual_cells_step()])

    np.testing.assert_array_equal(
        first.ns["annot_virtual_cells"][["x_um", "y_um"]].to_numpy(),
        second.ns["annot_virtual_cells"][["x_um", "y_um"]].to_numpy(),
    )


def test_only_the_selected_types_are_sampled():
    ex = _run([_polygons_step(), _virtual_cells_step(types=("tumour",))])

    assert set(ex.ns["annot_virtual_cells"]["annotation_type"]) == {"tumour"}


# ── neighbourhood enrichment ─────────────────────────────────────────────────

def test_annotation_types_join_the_clustering_as_extra_categories():
    adata = _adata()
    ex = _run([_normalize_step(), _polygons_step(), _virtual_cells_step(),
               _nhood_step()], adata)
    adata_annot = ex.ns["adata_annot"]

    assert list(adata_annot.obs[CLUSTER_KEY].cat.categories) == [
        "0", "1", "stroma", "tumour"]
    assert adata_annot.n_obs == adata.n_obs + len(ex.ns["annot_virtual_cells"])
    zscore = np.asarray(adata_annot.uns[f"{CLUSTER_KEY}_nhood_enrichment"]["zscore"])
    assert zscore.shape == (4, 4)


def test_virtual_cells_carry_no_expression():
    """They are places, not cells. Anything else would put invented counts into
    an object the notebook goes on to analyse."""
    ex = _run([_normalize_step(), _polygons_step(), _virtual_cells_step(),
               _nhood_step()], _adata())
    adata_annot = ex.ns["adata_annot"]
    virtual = adata_annot[adata_annot.obs_names.str.startswith("annot-")]

    assert virtual.X.sum() == 0


def test_the_graph_is_built_on_the_augmented_object():
    """Not reused from the spatial_neighbors step: adding virtual cells changes
    who every real cell's neighbours are, which is the analysis."""
    ex = _run([_normalize_step(), _polygons_step(), _virtual_cells_step(),
               _nhood_step()], _adata())

    assert ex.ns["adata_annot"].obsp["spatial_connectivities"].shape[0] == (
        ex.ns["adata_annot"].n_obs)


# ── distance ─────────────────────────────────────────────────────────────────

def test_distances_are_microns_to_the_nearest_boundary():
    adata = _adata()
    ex = _run([_polygons_step(), _distance_step()], adata)
    distances = ex.ns["annot_distances"]

    assert distances.shape == (adata.n_obs,)
    assert np.isfinite(distances).all()
    # Cells are spread over 0-100 um and the square spans 0-50 um, so the
    # farthest any of them can be from its boundary is the far corner.
    assert distances.max() < 100.0


def test_the_distances_land_in_obs_under_a_per_type_key():
    adata = _adata()
    _run([_polygons_step(), _distance_step("tumour"), _distance_step("stroma")],
         adata)

    assert "dist_to_tumour_um" in adata.obs
    assert "dist_to_stroma_um" in adata.obs
    assert not adata.obs["dist_to_tumour_um"].equals(adata.obs["dist_to_stroma_um"])


# ── figures ──────────────────────────────────────────────────────────────────

def test_the_nhood_heatmap_writes_the_files_it_names(tmp_path):
    paths = [str(tmp_path / "plots" / "annot_nhood.png")]
    ex = _run([_normalize_step(), _polygons_step(), _virtual_cells_step(),
               _nhood_step()], _adata())
    ex.run(Step(
        id=f"plot:annot_nhood:{CLUSTER_KEY}",
        template=builtin_assemble("annot.nhood_plot", _blocks("annot.nhood_plot")),
        params={"cluster_key": CLUSTER_KEY, "mode": "zscore",
                "title": "Annotation Nhood Enrichment", "paths": paths},
        kind="terminal", deps=[f"annot_nhood:{CLUSTER_KEY}"], outputs=["fig"],
    ))

    assert Path(paths[0]).exists()


@pytest.mark.parametrize("plot_type", ["violin", "box", "strip"])
def test_every_distance_plot_type_renders(tmp_path, plot_type):
    adata = _adata()
    paths = [str(tmp_path / "plots" / f"annot_distance_{plot_type}.png")]
    ex = _run([_polygons_step(), _distance_step()], adata)
    ex.run(Step(
        id=f"plot:annot_distance:tumour:{CLUSTER_KEY}",
        template=builtin_assemble(
            "annot.distance_plot", ["head", "clip", f"plot.{plot_type}", "save"]),
        params={"obs_key": "dist_to_tumour_um", "cluster_key": CLUSTER_KEY,
                "annotation_type": "tumour", "paths": paths, "max_dist": 60.0},
        kind="terminal", deps=["annot_distance:tumour"], outputs=["fig"],
    ))

    assert Path(paths[0]).exists()


def test_the_clip_block_drops_the_far_tail_rather_than_flattening_it():
    """"Max distance to show" must not clip values onto the limit: the drawn
    quartiles have to be the quartiles of what is shown."""
    adata = _adata()
    ex = _run([_polygons_step(), _distance_step()], adata)
    ns = dict(ex.ns)
    exec(compile(builtin_assemble("annot.distance_plot", ["head", "clip"]).replace(  # noqa: S102
        "$obs_key", "'dist_to_tumour_um'").replace(
        "$cluster_key", f"'{CLUSTER_KEY}'").replace(
        "$annotation_type", "tumour").replace("$max_dist", "20.0"),
        "<clip>", "exec"), ns)

    assert (ns["_dist"]["dist_to_tumour_um"] <= 20.0).all()
    assert len(ns["_dist"]) < adata.n_obs


# ── replay ───────────────────────────────────────────────────────────────────

def test_the_recorded_graph_replays_to_the_same_numbers():
    """The claim, end to end: run the steps, then execute the recorded cells in
    a clean namespace and compare. Nothing in the graph may depend on the GUI —
    the drawn shapes reach the replay as literals in the ``annotations`` cell."""
    adata = _adata()
    ex = _run([_normalize_step(), _polygons_step(), _virtual_cells_step(),
               _nhood_step(), _distance_step()], adata)

    replayed = _replay(ex, _adata())

    np.testing.assert_allclose(
        np.asarray(ex.ns["adata_annot"].uns[
            f"{CLUSTER_KEY}_nhood_enrichment"]["zscore"], dtype=float),
        np.asarray(replayed["adata_annot"].uns[
            f"{CLUSTER_KEY}_nhood_enrichment"]["zscore"], dtype=float),
    )
    np.testing.assert_allclose(ex.ns["annot_distances"],
                               replayed["annot_distances"])
    np.testing.assert_array_equal(
        ex.ns["annot_virtual_cells"][["x_um", "y_um"]].to_numpy(),
        replayed["annot_virtual_cells"][["x_um", "y_um"]].to_numpy(),
    )
