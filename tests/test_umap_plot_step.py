"""The UMAP figure step (issue #34), and the join that nearly broke it.

``embed.xenium`` reads Xenium's own UMAP so the notebook reproduces the figure
that was on screen rather than a freshly computed, equally valid, *different*
layout. Getting the coordinates is the whole job, and there is exactly one way
to get it wrong:

    spatialdata_io indexes the table **positionally** — obs_names are '0', '1',
    '2', … — and keeps the barcode in ``obs['cell_id']``. ``projection.csv`` is
    indexed by barcode. Reindexing on ``obs_names`` therefore matches *nothing*:
    every coordinate comes back NaN.

That failure is close to invisible. It does not raise; ``sc.pl.umap`` happily
draws an empty panel, and only the PDF backend eventually complains ("need at
least one array to concatenate") because a collection of all-NaN points has no
extent. It was found by running the step against a real dataset, not by reading
it, so it is pinned here against a fixture shaped the way the real table is.

No Qt, no dataset: a synthetic AnnData plus a projection.csv on disk.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
ad = pytest.importorskip("anndata")

import matplotlib                                                # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

from xenium_viewer.utils.step_templates import (                 # noqa: E402
    builtin_assemble, builtin_text,
)
from xenium_viewer.utils.steps import Step, StepExecutor         # noqa: E402

N = 60
GENES = ["Gene0", "Gene1", "Gene2"]


@pytest.fixture
def dataset(tmp_path):
    """A table shaped like the real one, plus Xenium's projection.csv.

    Positional ``obs_names`` and a ``cell_id`` column: that pairing is the whole
    point of the fixture, so a version that joined on the index would pass.
    """
    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(3, (N, 5)).astype("float32"))
    adata.obs_names = [str(i) for i in range(N)]
    adata.var_names = [f"Gene{i}" for i in range(5)]
    barcodes = [f"cell{i}-1" for i in range(N)]
    adata.obs["cell_id"] = barcodes
    adata.obs["leiden"] = pd.Categorical([str(i % 3) for i in range(N)])

    projection = tmp_path / "analysis" / "umap" / "gene_expression_2_components"
    projection.mkdir(parents=True)
    coords = rng.normal(size=(N, 2))
    # Xenium omits a handful of cells from its UMAP — 91 on the dataset this was
    # measured against — so the join is deliberately not total.
    frame = pd.DataFrame(coords[:-3], index=barcodes[:-3],
                         columns=["UMAP-1", "UMAP-2"])
    frame.to_csv(projection / "projection.csv")
    return adata, tmp_path


def _executor(adata, data_path):
    namespace = {"sc": sc, "pd": pd, "np": np, "plt": plt, "Path": Path,
                 "data_path": data_path, "adata": adata}
    executor = StepExecutor(namespace=namespace)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        executor.run(Step(id="normalize", template=builtin_text("normalize"),
                          kind="setup", outputs=["adata_norm"]))
    return executor


def _run(executor, blocks, params):
    step = Step(id="plot:umap", template=builtin_assemble("umap.plot", blocks),
                params=params, deps=["normalize"], kind="terminal",
                outputs=["fig"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return step, executor.run(step)


def test_the_xenium_embedding_is_joined_on_cell_id_not_on_obs_names(dataset, tmp_path):
    """The regression itself: obs_names are positional, the join key is not."""
    adata, data_path = dataset
    executor = _executor(adata, data_path)
    _run(executor, ["embed.xenium", "color.genes", "save"],
         {"color": GENES[:1], "cmap": "viridis", "ncols": 1,
          "paths": [str(tmp_path / "plots" / "u.png")]})

    coords = executor.ns["adata_norm"].obsm["X_umap"]
    finite = int(np.isfinite(coords[:, 0]).sum())
    assert finite == N - 3, (
        f"only {finite}/{N} cells got coordinates — the join matched the wrong "
        f"key (obs_names are positional; the barcode is in obs['cell_id'])"
    )
    plt.close("all")


def test_the_template_names_the_cell_id_column():
    """Asserted on the text too, so the reason survives a refactor of the test."""
    source = builtin_assemble("umap.plot", ["embed.xenium"])
    assert "obs['cell_id']" in source
    assert "reindex(adata_norm.obs_names)" not in source


def test_a_pdf_is_written_and_not_only_a_png(dataset, tmp_path):
    """The all-NaN join surfaced *only* in the PDF backend.

    ``sc.pl.umap`` draws an empty panel without complaint and the PNG saves
    fine; ``backend_pdf`` is what raised. Writing both formats is therefore not
    just issue #35's ask — it is what makes an empty figure loud.
    """
    adata, data_path = dataset
    executor = _executor(adata, data_path)
    paths = [str(tmp_path / "plots" / f"u.{ext}") for ext in ("png", "pdf")]
    step, _ = _run(executor, ["embed.xenium", "color.genes", "save"],
                   {"color": GENES, "cmap": "viridis", "ncols": 3,
                    "paths": paths})
    for path in paths:
        assert Path(path).stat().st_size > 0
    plt.close("all")


def test_one_gene_gives_one_panel_and_several_give_a_grid(dataset, tmp_path):
    """Issue #34's actual ask, and each panel carries its own colour bar."""
    adata, data_path = dataset
    executor = _executor(adata, data_path)

    _, out = _run(executor, ["embed.xenium", "color.genes", "save"],
                  {"color": GENES[:1], "cmap": "viridis", "ncols": 3,
                   "paths": [str(tmp_path / "plots" / "one.png")]})
    single = out["fig"]
    assert len(single.axes) == 2, "one panel plus its colour bar"
    plt.close("all")

    _, out = _run(executor, ["embed.xenium", "color.genes", "save"],
                  {"color": GENES, "cmap": "viridis", "ncols": 3,
                   "paths": [str(tmp_path / "plots" / "three.png")]})
    assert len(out["fig"].axes) == 6, "three panels, three colour bars"
    plt.close("all")


def test_the_cluster_assembly_relabels_from_the_display_names(dataset, tmp_path):
    adata, data_path = dataset
    executor = _executor(adata, data_path)
    step, out = _run(
        executor, ["embed.xenium", "relabel", "color.clusters", "save"],
        {"color": ["leiden"], "groupby": "leiden",
         "categories": {"0": "Tumour", "1": "Stroma", "2": "Immune"},
         "paths": [str(tmp_path / "plots" / "c.png")]})

    assert "{'0': 'Tumour', '1': 'Stroma', '2': 'Immune'}" in step.render()
    assert list(executor.ns["adata_norm"].obs["leiden"].cat.categories) == [
        "Tumour", "Stroma", "Immune"], "cluster order should survive relabelling"
    plt.close("all")


def test_two_clusters_may_share_a_display_name(dataset, tmp_path):
    """Reported from a real session: *ValueError: Categorical categories must
    be unique*, raised by ``.cat.rename_categories`` at statement 4 of 6.

    Giving two clusters one label is an ordinary thing to want — it means "these
    are the same cell type". ``.map()`` merges them; ``rename_categories``
    refuses, because it can only rename one-for-one.
    """
    adata, data_path = dataset
    executor = _executor(adata, data_path)
    step, _ = _run(
        executor, ["embed.xenium", "relabel", "color.clusters", "save"],
        {"color": ["leiden"], "groupby": "leiden",
         "categories": {"0": "Tumour", "1": "Stroma", "2": "Tumour"},
         "paths": [str(tmp_path / "plots" / "merged.png")]})

    column = executor.ns["adata_norm"].obs["leiden"]
    assert list(column.cat.categories) == ["Tumour", "Stroma"]
    # Cells from both source clusters really are in the merged category.
    assert (column == "Tumour").sum() == (adata.obs["leiden"].isin(["0", "2"])).sum()
    assert Path(step.params["paths"][0]).stat().st_size > 0
    plt.close("all")


def test_the_relabel_block_does_not_use_rename_categories():
    """Asserted on the text, so the reason survives a refactor of the test.

    Comment lines are stripped first — the block explains *why* it avoids
    ``rename_categories``, and that explanation must not trip the check.
    """
    source = builtin_assemble("umap.plot", ["relabel"])
    code = "\n".join(line for line in source.splitlines()
                     if not line.lstrip().startswith("#"))
    assert "rename_categories" not in code
    assert ".map(" in code


def test_the_fallback_computes_its_own_embedding_and_says_so(dataset, tmp_path):
    """A Crop Dataset export has no ``analysis/`` folder to read."""
    adata, data_path = dataset
    executor = _executor(adata, data_path)
    step, out = _run(executor, ["embed.recompute", "color.genes", "save"],
                     {"color": GENES[:1], "cmap": "viridis", "ncols": 1,
                      "paths": [str(tmp_path / "plots" / "r.png")]})

    source = step.render()
    assert "sc.tl.umap" in source and "projection.csv" not in source
    assert np.isfinite(executor.ns["adata_norm"].obsm["X_umap"]).all()
    plt.close("all")


def test_the_tab_picks_the_embedding_block_from_the_dataset(dataset, tmp_path):
    """``has_xenium_umap`` is what chooses between the two embed blocks."""
    pytest.importorskip("qtpy")
    from xenium_viewer.tabs.tab_umap import has_xenium_umap

    _adata, data_path = dataset
    assert has_xenium_umap(data_path)
    assert not has_xenium_umap(tmp_path / "nope")
    assert not has_xenium_umap(None)
