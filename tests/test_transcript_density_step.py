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


# ── the preview tier ─────────────────────────────────────────────────────────
# The density heatmap can be previewed from the per-gene feather index while
# the user hunts for a gene, because the two routes were measured to agree
# exactly. These tests are what turn that one-off measurement into a gate: if
# the preview and the recorded step ever stop agreeing, the preview is a lie
# about what Compute Density will draw.

def _write_feather_cache(tmp_path, sdata):
    """The fixture's transcripts, as `palms-preprocess` would leave them.

    Same two filters (`is_gene`, `qv >= MIN_QV`), same column names, one file
    per gene — including a `cell_id`, which caches built before it was kept do
    not have. `test_a_cache_without_cell_ids_...` writes the other shape.
    """
    from palms.preprocess import MIN_QV, _write_sentinel

    cache = tmp_path / "transcript_cache"
    cache.mkdir()
    df = sdata.points["transcripts"].compute()
    df = df[df["is_gene"].astype(bool) & (df["qv"] >= MIN_QV)]
    for gene, sub in df.groupby("feature_name", observed=True):
        sub.rename(columns={"x": "x_location", "y": "y_location"})[
            ["x_location", "y_location", "qv", "feature_name", "cell_id"]
        ].reset_index(drop=True).to_feather(cache / f"{gene}.feather")
    # The parquet the cache claims to have come from, and the stamp that says
    # so — written with the real function, so the test cannot drift from what
    # `palms-preprocess` actually leaves behind.
    parquet = tmp_path / "transcripts.parquet"
    plain = df.reset_index(drop=True)
    plain.attrs = {}          # the element's Scale rides along and is not JSON
    plain.to_parquet(parquet)
    _write_sentinel(cache, parquet)
    return cache


def _loader(cache_dir):
    from palms.utils.transcript_index import TranscriptLoader

    return TranscriptLoader(cache_dir=cache_dir,
                            parquet_path=cache_dir.parent / "transcripts.parquet",
                            min_qv=20, pixel_size=PIXEL_SIZE)


def _preview_density(loader, sdata, blocks, **params):
    """Run the preview the way the tab does: feather rows, same template text."""
    from palms.utils.transcript_index import points_for_preview

    step = _density_step(blocks, **params)
    points, reason = points_for_preview(
        loader, sdata, params.get("gene", GENE), params.pop("min_qv", 20),
        need_cell_id="cell_id" in step.render())
    assert points is not None, reason
    ex = StepExecutor(namespace=_ns(sdata))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ex, ex.preview(step, bindings={"transcript_points": points})


@pytest.mark.parametrize("blocks", [
    ("head", "main"),
    ("head", "filter", "main"),
    ("head", "main", "normalise"),
    ("head", "filter", "main", "normalise"),
])
def test_the_feather_preview_bins_identically_to_the_recorded_step(tmp_path, blocks):
    """The claim the whole preview tier rests on, for every assembly.

    Measured once by hand on a real dataset — 211,154 rows matching
    `TranscriptLoader.load_gene` exactly and the histogram equal bin for bin,
    max absolute difference 0. Asserted here so it stays true.
    """
    sdata = _sdata()
    cache = _write_feather_cache(tmp_path, sdata)
    extra = {"clustering": CLUSTER_KEY, "selected": ["0"]} if "filter" in blocks else {}

    recorded = _run([_gene_step(), _density_step(blocks, **extra)], sdata)
    _, preview = _preview_density(_loader(cache), sdata, blocks, **extra)

    np.testing.assert_array_equal(recorded.ns["transcript_density"],
                                  preview["transcript_density"])


def test_a_preview_of_the_density_records_nothing(tmp_path):
    """Browsing genes must not leave a trace in the notebook."""
    sdata = _sdata()
    cache = _write_feather_cache(tmp_path, sdata)
    ex, _ = _preview_density(_loader(cache), sdata, ("head", "main"))

    assert list(ex.graph.nodes()) == []
    assert "transcript_points" not in ex.names()
    assert "transcript_density" not in ex.names()


def test_a_preview_still_works_after_a_recorded_run(tmp_path):
    """The regression a real dataset caught and the synthetic tests did not.

    `transcript_points` is bound in the shared namespace by the recorded fetch,
    and an earlier version of `preview` refused to shadow any bound name — so
    the preview raised, and stopped working, from the first Compute Density
    onwards. Substituting that result is the entire mechanism; it is the *base*
    names a preview may not swap.
    """
    sdata = _sdata()
    cache = _write_feather_cache(tmp_path, sdata)
    from palms.utils.transcript_index import points_for_preview

    ex = _run([_gene_step(), _density_step(("head", "main"))], sdata)
    assert "transcript_points" in ex.names()          # the recorded run bound it

    points, reason = points_for_preview(_loader(cache), sdata, GENE, 20,
                                        need_cell_id=False)
    assert points is not None, reason
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        preview = ex.preview(_density_step(("head", "main")),
                             bindings={"transcript_points": points})

    np.testing.assert_array_equal(ex.ns["transcript_density"],
                                  preview["transcript_density"])


def test_the_preview_takes_its_frame_from_the_element_not_the_pixel_size(tmp_path):
    """The scale comes off the transcripts element, the way the template asks
    spatialdata to apply it — not from dividing by a pixel size passed in."""
    from palms.utils.transcript_index import points_for_preview

    sdata = _sdata()
    cache = _write_feather_cache(tmp_path, sdata)
    loader = _loader(cache)
    loader.pixel_size = PIXEL_SIZE * 7      # nonsense; must not be consulted

    points, reason = points_for_preview(loader, sdata, GENE, 20, need_cell_id=False)
    assert points is not None, reason
    recorded = _run([_gene_step()], sdata).ns["transcript_points"]
    np.testing.assert_allclose(np.sort(points["x"].to_numpy()),
                               np.sort(recorded["x"].to_numpy()))


def test_a_quality_floor_below_the_caches_is_refused(tmp_path):
    """The feather files are filtered at MIN_QV, so the slider cannot go under
    it — returning the qv>=20 subset for a qv>=10 question is the wrong picture."""
    from palms.utils.transcript_index import points_for_preview

    sdata = _sdata()
    loader = _loader(_write_feather_cache(tmp_path, sdata))

    points, reason = points_for_preview(loader, sdata, GENE, 10, need_cell_id=False)
    assert points is None and "Min QV 20" in reason

    # Above the floor it is a real filter, applied to the index's own column.
    points, reason = points_for_preview(loader, sdata, GENE, 30, need_cell_id=False)
    assert points is not None, reason
    assert len(points) == len(_run([_gene_step(min_qv=30)], sdata).ns["transcript_points"])


def test_a_cache_without_cell_ids_refuses_only_the_filtered_preview(tmp_path):
    """Caches built before cell ids were kept are still on disk and still load.
    Only the cluster filter reads the column, so only it may be refused."""
    from palms.utils.transcript_index import points_for_preview

    sdata = _sdata()
    cache = _write_feather_cache(tmp_path, sdata)
    for path in cache.glob("*.feather"):
        pd.read_feather(path).drop(columns=["cell_id"]).to_feather(path)
    loader = _loader(cache)

    points, reason = points_for_preview(loader, sdata, GENE, 20, need_cell_id=True)
    assert points is None and "palms-preprocess" in reason

    points, reason = points_for_preview(loader, sdata, GENE, 20, need_cell_id=False)
    assert points is not None, reason


def test_an_index_built_from_another_parquet_is_refused(tmp_path):
    """The loader globs `*.feather` and never checks the sentinel — fine for the
    overlay, not here. A cache built from a different transcripts.parquet would
    preview rows the recorded step would not reproduce, which is the one way a
    preview could be a lie about the result rather than merely unavailable."""
    from palms.utils.transcript_index import points_for_preview

    sdata = _sdata()
    loader = _loader(_write_feather_cache(tmp_path, sdata))

    # The parquet changed after the index was built.
    (tmp_path / "transcripts.parquet").write_bytes(b"different bytes entirely")
    points, reason = points_for_preview(loader, sdata, GENE, 20, need_cell_id=False)
    assert points is None and "transcripts.parquet" in reason

    # And an index with no stamp at all is refused the same way.
    (loader.cache_dir / ".complete").unlink()
    points, reason = points_for_preview(loader, sdata, GENE, 20, need_cell_id=False)
    assert points is None


def test_a_gene_outside_the_cache_is_not_a_parquet_scan(tmp_path):
    """The loader's fallback is a ~22 s scan that also drops cell_id. A preview
    that takes 22 s is not a preview, so it must refuse instead."""
    from palms.utils import transcript_index

    sdata = _sdata()
    loader = _loader(_write_feather_cache(tmp_path, sdata))
    (loader.cache_dir / f"{GENE}.feather").unlink()
    loader._cached_genes = None

    def _boom(*a, **k):
        raise AssertionError("the preview fell through to the parquet scan")

    loader._load_from_parquet = _boom
    points, reason = transcript_index.points_for_preview(
        loader, sdata, GENE, 20, need_cell_id=False)
    assert points is None and "not in the transcript index" in reason
