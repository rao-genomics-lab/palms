"""The transcript density heatmap, as recorded code.

This was the last node in the app recording prose — a `TERMINAL` whose cell was
a comment, which replays as a silent no-op: the notebook "passes" with the
analysis simply missing. It was held back because the viewer computed the
histogram from its own per-gene feather index (`utils/transcript_index.py`,
built by `palms-preprocess`), which a notebook replaying from raw Xenium output
does not have.

It reads `sdata.points['transcripts']` now, in two steps — the fetch, which is
the slow half, and the histogram, which is the half the bin-size slider moves.
The tests below build a small synthetic SpatialData with the same shape as a
real one: transcripts in microns carrying a `Scale` to a pixel-frame image.
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
plt = pytest.importorskip("matplotlib.pyplot")

from palms.utils.step_templates import builtin_assemble, builtin_text  # noqa: E402
from palms.utils.steps import Step, StepExecutor, check_step  # noqa: E402

PIXEL_SIZE = 0.5
#: 100 x 200 pixels, i.e. 50 x 100 microns.
IMAGE_HW = (100, 200)
CLUSTER_KEY = "leiden_r1.0"
GENE = "GeneA"


def _sdata():
    """A miniature Xenium: transcripts in microns, image in pixels, one table.

    The `Scale` on the points is the thing under test as much as the histogram
    is — it is what the template asks spatialdata to apply, instead of dividing
    by the pixel size itself.
    """
    from spatialdata.models import Image2DModel, PointsModel, TableModel

    rng = np.random.default_rng(0)

    # Two genes. GeneA sits in the left half of the image, GeneB in the right,
    # so a wrong axis order or a missed filter shows up as a shifted picture.
    rows = []
    for gene, x_range in ((GENE, (0.0, 20.0)), ("GeneB", (30.0, 50.0))):
        for i in range(200):
            rows.append({
                "x": rng.uniform(*x_range),
                "y": rng.uniform(0.0, 50.0),
                "feature_name": gene,
                "cell_id": f"cell{i % 20}",
                # Every twentieth call is below the quality floor, and every
                # twenty-fifth is a control probe rather than a gene.
                "qv": 10.0 if i % 20 == 0 else 40.0,
                "is_gene": i % 25 != 0,
            })
    points = PointsModel.parse(
        pd.DataFrame(rows),
        coordinates={"x": "x", "y": "y"},
        transformations={"global": sd.transformations.Scale(
            [1 / PIXEL_SIZE, 1 / PIXEL_SIZE], axes=("x", "y"))},
    )

    image = Image2DModel.parse(
        np.zeros((1, *IMAGE_HW), dtype="uint16"), dims=("c", "y", "x"))

    n = 20
    adata = anndata.AnnData(rng.poisson(3, size=(n, 4)).astype("float32"))
    adata.obs_names = [str(i) for i in range(n)]
    adata.var_names = [f"Gene{i}" for i in range(4)]
    adata.obs["cell_id"] = [f"cell{i}" for i in range(n)]
    adata.obs[CLUSTER_KEY] = pd.Categorical(["0"] * (n // 2) + ["1"] * (n - n // 2))
    adata.obs["region"] = pd.Categorical(["cells"] * n)
    adata.obs["instance_id"] = np.arange(n)
    adata.obsm["spatial"] = np.column_stack([
        rng.uniform(0, 50, n), rng.uniform(0, 50, n)])   # microns
    adata = TableModel.parse(adata, region="cells", region_key="region",
                             instance_key="instance_id")

    return sd.SpatialData(
        points={"transcripts": points},
        images={"morphology_focus": image},
        tables={"table": adata},
    )


def _ns(sdata):
    """Mirrors EXECUTOR_BASE_NAMES, which is what the notebook preamble binds."""
    return {"sc": sc, "sq": sq, "sd": sd, "pd": pd, "np": np, "plt": plt,
            "Path": Path, "data_path": Path("/tmp/xv-not-written"),
            "sdata": sdata, "adata": sdata.tables["table"]}


def _gene_step(gene=GENE, min_qv=20):
    return Step(id=f"transcripts:{gene}", template=builtin_text("transcripts.gene"),
                params={"gene": gene, "min_qv": min_qv},
                outputs=["transcript_points"])


def _density_step(blocks=("head", "main"), **params):
    full = {"gene": GENE, "bin_size_um": 10.0, "pixel_size": PIXEL_SIZE}
    full.update(params)
    return Step(
        id=f"transcript_density:{full['gene']}",
        template=builtin_assemble("transcripts.density", list(blocks)),
        params=full, deps=[f"transcripts:{full['gene']}"],
        outputs=["transcript_density"],
    )


def _run(steps, sdata=None):
    sdata = sdata if sdata is not None else _sdata()
    ex = StepExecutor(namespace=_ns(sdata))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for step in steps:
            ex.run(step)
    return ex


def test_both_templates_are_self_contained():
    """Rule (a): only EXECUTOR_BASE_NAMES plus what a declared dependency binds.
    A name in neither replays as a NameError in a clean kernel."""
    base = set(_ns(_sdata()))

    assert check_step(_gene_step(), base) == set()
    for blocks in (("head", "main"), ("head", "filter", "main"),
                   ("head", "main", "normalise"),
                   ("head", "filter", "main", "normalise")):
        step = _density_step(blocks, clustering=CLUSTER_KEY, selected=["0"])
        assert check_step(step, base | {"transcript_points"}) == set(), blocks


#: Of the 200 rows per gene the fixture writes, i%20==0 is below the quality
#: floor (10) and i%25==0 is a control probe (8); i=0 and i=100 are both, so 16
#: rows are dropped and 184 survive.
KEPT = 184


def test_the_fetch_returns_one_gene_in_the_image_frame():
    ex = _run([_gene_step()])
    points = ex.ns["transcript_points"]

    assert len(points) == KEPT
    assert set(points.columns) >= {"x", "y", "cell_id"}
    # Microns x2 (pixel size 0.5) — GeneA spans 0-20 um, so 0-40 px.
    assert points["x"].max() <= 40.0 + 1e-6
    assert points["x"].max() > 20.0, "the declared Scale was not applied"


def test_the_histogram_matches_the_image_grid():
    ex = _run([_gene_step(), _density_step()])
    density = ex.ns["transcript_density"]

    # 10 um bins over a 100x200 px image at 0.5 um/px = 20 px bins -> 5 x 10.
    assert density.shape == (5, 10)
    assert density.sum() == KEPT


def test_the_quality_floor_and_the_control_probes_are_filtered_out():
    """The two filters palms-preprocess bakes into its feather files. The
    density read those files before it read the points element, so dropping
    either here would silently raise every count."""
    assert len(_run([_gene_step()]).ns["transcript_points"]) == KEPT
    # Below every qv in the fixture: only the control probes are left out.
    assert len(_run([_gene_step(min_qv=0)]).ns["transcript_points"]) == 192
    # Above every qv: nothing survives.
    assert len(_run([_gene_step(min_qv=50)]).ns["transcript_points"]) == 0


def test_transcripts_land_where_they_were_drawn():
    """A transposed histogram or a dropped transform would still sum to 200.

    GeneA spans x 0-20 um, i.e. 0-40 px, i.e. the first two of the ten 20 px
    columns. The image is wider than the transcripts, so the rest stays empty.
    """
    ex = _run([_gene_step(), _density_step()])
    density = ex.ns["transcript_density"]

    assert density[:, :2].sum() == KEPT
    assert density[:, 2:].sum() == 0


def test_a_second_gene_is_binned_separately():
    """GeneB spans x 30-50 um — 60-100 px, columns 3 and 4 — so a filter that
    silently matched everything would light up columns 0 and 1 as well."""
    ex = _run([_gene_step("GeneB"), _density_step(gene="GeneB")])
    density = ex.ns["transcript_density"]

    assert density[:, 3:5].sum() == KEPT
    assert density[:, :3].sum() == 0
    assert density[:, 5:].sum() == 0


def test_bin_size_sets_the_grid():
    coarse = _run([_gene_step(), _density_step(bin_size_um=25.0)])

    assert coarse.ns["transcript_density"].shape == (2, 4)
    assert coarse.ns["transcript_density"].sum() == KEPT


def test_the_cluster_filter_keeps_only_transcripts_from_selected_cells():
    """Per cell, via the assignment on the points element — the viewer used to
    fall back to a bin-level approximation when its cache lacked cell ids."""
    ex = _run([_gene_step(),
               _density_step(("head", "filter", "main"),
                             clustering=CLUSTER_KEY, selected=["0"])])

    total = _run([_gene_step(), _density_step()]).ns["transcript_density"].sum()
    filtered = ex.ns["transcript_density"].sum()
    assert 0 < filtered < total

    # Exactly the transcripts assigned to a cluster "0" cell — cells 0-9. Counted
    # off the fetched points rather than assumed to be half: the quality floor
    # does not drop rows evenly across cells.
    points = _run([_gene_step()]).ns["transcript_points"]
    expected = points["cell_id"].isin({f"cell{i}" for i in range(10)}).sum()
    assert filtered == expected


def test_normalising_divides_by_the_cells_in_each_bin():
    plain = _run([_gene_step(), _density_step()]).ns["transcript_density"]
    normalised = _run([_gene_step(),
                       _density_step(("head", "main", "normalise"))]).ns["transcript_density"]

    assert normalised.shape == plain.shape
    assert normalised.sum() != plain.sum()
    # Empty bins stay 0 rather than becoming inf or nan — the whole point of the
    # np.where, and what the napari contrast limits depend on.
    assert np.isfinite(normalised).all()
    assert (normalised[plain == 0] == 0).all()


def test_the_recorded_cells_replay_to_the_same_histogram():
    """The claim, end to end: run the steps, then execute the recorded cells in
    a clean namespace and compare."""
    sdata = _sdata()
    ex = _run([_gene_step(), _density_step(("head", "main", "normalise"))], sdata)

    ns = _ns(_sdata())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for node_id in ex.graph.topo_sort():
            exec(compile(ex.graph.get(node_id).code, "<replay>", "exec"), ns)  # noqa: S102

    np.testing.assert_array_equal(ex.ns["transcript_density"],
                                  ns["transcript_density"])
