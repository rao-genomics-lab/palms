"""The rest of the E2 migration: spatial, marker, correlation and ROI steps.

Each of these tabs previously ran one expression and recorded a different one.
The tests below execute the *real* templates from the tab modules against a
small synthetic AnnData, then replay the recorded graph in a clean namespace and
compare — which is the reproducibility claim reduced to something CI can assert.

Importing a tab module pulls in Qt/napari, so run headless:
    QT_QPA_PLATFORM=offscreen pytest tests/test_spatial_roi_steps.py
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
pytest.importorskip("shapely")
plt = pytest.importorskip("matplotlib.pyplot")

from xenium_viewer.utils.steps import Step, StepExecutor, check_step  # noqa: E402
from xenium_viewer.tabs._helpers import (  # noqa: E402
    _NORMALIZE_TEMPLATE, _SPATIAL_NEIGHBORS_TEMPLATE,
)
from xenium_viewer.tabs.tab_nhood import _NHOOD_TEMPLATE  # noqa: E402
from xenium_viewer.tabs.tab_co_occurrence import _COOCCUR_TEMPLATE  # noqa: E402
from xenium_viewer.tabs.tab_ligrec import _ligrec_template  # noqa: E402
from xenium_viewer.tabs.tab_marker_genes import _marker_plot_template  # noqa: E402
from xenium_viewer.tabs.tab_gene_correlation import _gene_corr_template  # noqa: E402
from xenium_viewer.tabs.tab_roi import _ROIS_TEMPLATE, _roi_deg_template  # noqa: E402


CLUSTER_KEY = "leiden_r1.0"


def _adata(n_obs: int = 120, n_vars: int = 30):
    """Two spatially separated populations, so the spatial stats have signal."""
    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(n_obs, n_vars)).astype("float32")
    half = n_obs // 2
    counts[:half, : n_vars // 3] += 12
    a = anndata.AnnData(counts)
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_vars)]
    xy = rng.uniform(0, 100, size=(n_obs, 2))
    xy[:half, 0] += 200          # push population A away from B
    a.obsm["spatial"] = xy
    a.obs[CLUSTER_KEY] = pd.Categorical(["0"] * half + ["1"] * (n_obs - half))
    return a


def _ns(adata=None):
    return {"sc": sc, "sq": sq, "pd": pd, "np": np, "plt": plt,
            "adata": adata if adata is not None else _adata()}


def _executor(adata=None):
    return StepExecutor(namespace=_ns(adata))


def _normalize_step():
    return Step(id="normalize", template=_NORMALIZE_TEMPLATE, kind="setup",
                outputs=["adata_norm"])


def _clustering_step():
    """Stand-in for the real ``clustering:<key>`` node, which the steps under
    test declare as a dependency (the graph rejects unknown deps at record
    time). Mirrors what ``record_clustering`` emits: a categorical in obs."""
    return Step(
        id=f"clustering:{CLUSTER_KEY}",
        template="\nadata.obs[$key] = pd.Categorical($labels)",
        params={"key": CLUSTER_KEY, "labels": ["0"] * 60 + ["1"] * 60},
    )


def _neighbors_step(n_neighs=6):
    return Step(id="spatial_neighbors", template=_SPATIAL_NEIGHBORS_TEMPLATE,
                params={"n_neighs": n_neighs}, deps=["normalize"])


def _run(steps, adata=None):
    ex = _executor(adata)
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


# ── spatial neighbours ───────────────────────────────────────────────────────

def test_spatial_neighbors_builds_the_graph_on_the_object_consumers_use():
    """Regression: the old node built the graph on `adata`, but every consumer
    was handed `adata_norm` — so a replay ran nhood against no graph at all."""
    assert "sq.gr.spatial_neighbors(adata_norm" in _SPATIAL_NEIGHBORS_TEMPLATE
    ex = _run([_normalize_step(), _neighbors_step()])
    assert "spatial_connectivities" in ex.ns["adata_norm"].obsp


def test_spatial_neighbors_records_the_k_it_used():
    ex = _run([_normalize_step(), _neighbors_step(n_neighs=9)])
    assert "n_neighs=9" in ex.graph.get("spatial_neighbors").code


# ── neighbourhood enrichment ─────────────────────────────────────────────────

def _nhood_step(n_perms=20):
    return Step(
        id=f"nhood:{CLUSTER_KEY}", template=_NHOOD_TEMPLATE,
        params={"cluster_key": CLUSTER_KEY,
                "uns_key": f"{CLUSTER_KEY}_nhood_enrichment",
                "n_perms": n_perms, "seed": 42},
        deps=[f"clustering:{CLUSTER_KEY}", "spatial_neighbors"],
        outputs=["adata_norm"],
    )


def test_nhood_runs_on_the_normalized_copy_not_raw_counts():
    assert "sq.gr.nhood_enrichment(\n    adata_norm" in _NHOOD_TEMPLATE
    assert "nhood_enrichment(adata," not in _NHOOD_TEMPLATE.replace(" ", "")


def test_nhood_replays_to_the_same_zscores():
    ex = _run([_normalize_step(), _clustering_step(), _neighbors_step(), _nhood_step()])
    gui = np.asarray(ex.ns["nhood_zscore"])
    replayed = np.asarray(_replay(ex)["nhood_zscore"])
    np.testing.assert_allclose(gui, replayed)


def test_nhood_is_self_contained():
    available = {"sc", "sq", "pd", "np", "adata", "adata_norm"}
    assert check_step(_nhood_step(), available=available) == set()


# ── co-occurrence ────────────────────────────────────────────────────────────

def _cooccur_step(interval=10):
    return Step(
        id=f"cooccur:{CLUSTER_KEY}", template=_COOCCUR_TEMPLATE,
        params={"cluster_key": CLUSTER_KEY, "interval": interval},
        deps=[f"clustering:{CLUSTER_KEY}"], outputs=["adata_norm"],
    )


def test_cooccur_runs_on_the_normalized_copy_and_needs_no_hand_built_spatial():
    """The old cell rebuilt obsm['spatial'] from x_centroid/y_centroid columns
    that the Xenium table does not have under those names."""
    assert "sq.gr.co_occurrence(adata_norm" in _COOCCUR_TEMPLATE
    assert "x_centroid" not in _COOCCUR_TEMPLATE


def test_cooccur_replays_to_the_same_matrix():
    ex = _run([_normalize_step(), _clustering_step(), _cooccur_step()])
    key = f"{CLUSTER_KEY}_co_occurrence"
    gui = np.asarray(ex.ns["adata_norm"].uns[key]["occ"])
    replayed = np.asarray(_replay(ex)["adata_norm"].uns[key]["occ"])
    np.testing.assert_allclose(gui, replayed)


# ── ligand-receptor ──────────────────────────────────────────────────────────

def _ligrec_step(include=None, resources=None):
    params = {"cluster_key": CLUSTER_KEY, "n_perms": 10,
              "threshold": 0.01, "seed": 42}
    if include:
        params["include"] = include
    if resources:
        params["resources"] = resources
    return Step(
        id=f"ligrec:{CLUSTER_KEY}",
        template=_ligrec_template(bool(include), resources is not None),
        params=params, deps=[f"clustering:{CLUSTER_KEY}", "spatial_neighbors"],
        outputs=["ligrec_res"],
    )


def test_ligrec_records_the_interaction_databases_as_code_not_prose():
    """The checkbox selection used to survive only as a `# interactions: ...`
    comment, so a replay silently fell back to omnipath's defaults."""
    source = _ligrec_step(include=["OMNIPATH", "KINASE_EXTRA"]).render()
    assert "InteractionDataset[_n] for _n in ['OMNIPATH', 'KINASE_EXTRA']" in source
    assert "from omnipath.constants import InteractionDataset" in source


def test_ligrec_records_the_cellphonedb_restriction():
    source = _ligrec_step(resources="CellPhoneDB").render()
    assert "interactions_params['resources'] = 'CellPhoneDB'" in source


def test_ligrec_runs_on_the_normalized_copy_with_use_raw_false():
    source = _ligrec_step().render()
    assert "adata_norm, cluster_key=" in source
    assert "use_raw=False" in source
    assert "copy=True" in source


def test_ligrec_omits_the_databases_block_when_nothing_is_selected():
    source = _ligrec_step().render()
    assert "interactions_params = {}" in source
    assert "InteractionDataset[" not in source


# ── marker gene plots ────────────────────────────────────────────────────────

MARKERS = {"A": ["gene0", "gene1"], "B": ["gene2"]}


def _marker_step(plot_name="dotplot", relabel=False, tmp_path=None):
    params = {"plot_name": plot_name, "groupby": CLUSTER_KEY,
              "markers": MARKERS, "path": str(tmp_path / f"{plot_name}.png")}
    if relabel:
        params["categories"] = ["Tumour", "Stroma"]
    return Step(
        id=f"plot:markers:{plot_name}:{CLUSTER_KEY}",
        template=_marker_plot_template(plot_name, relabel, dpi=True),
        params=params, deps=["normalize", f"clustering:{CLUSTER_KEY}"],
        kind="terminal",
    )


@pytest.mark.parametrize("plot_name", [
    "dotplot", "heatmap", "matrixplot", "tracksplot", "correlation_matrix",
])
def test_marker_plot_templates_execute_and_write_a_file(plot_name, tmp_path):
    """This tab recorded nothing at all before, despite being plain scanpy."""
    step = _marker_step(plot_name, tmp_path=tmp_path)
    ex = _run([_normalize_step(), _clustering_step(), step])
    assert Path(step.params["path"]).exists()
    assert ex.graph.get(step.id).code == step.render()
    plt.close("all")


def test_marker_plot_carries_display_labels_into_the_recorded_source(tmp_path):
    source = _marker_step(relabel=True, tmp_path=tmp_path).render()
    assert "rename_categories(['Tumour', 'Stroma'])" in source


def test_marker_dict_survives_as_a_dict_literal(tmp_path):
    """string.Template is used precisely so `{...}` in params survives."""
    source = _marker_step(tmp_path=tmp_path).render()
    assert "var_names={'A': ['gene0', 'gene1'], 'B': ['gene2']}" in source


# ── gene correlation ─────────────────────────────────────────────────────────

def _gene_corr_step(norm="Log1p(CPM)", filtered=False, tmp_path=None):
    params = {"gene_a": "gene0", "gene_b": "gene1", "norm_label": "log1p(CPM)",
              "xlabel": "gene0", "ylabel": "gene1", "title_prefix": "gene0 vs gene1",
              "path": str(tmp_path / "corr.png")}
    if filtered:
        params["clustering"] = CLUSTER_KEY
        params["selected"] = ["0"]
    return Step(
        id="plot:gene_correlation",
        template=_gene_corr_template(norm, filtered),
        params=params, deps=["normalize"], kind="terminal",
        outputs=["fig", "x", "pr", "pp", "sr", "sp"],
    )


@pytest.mark.parametrize("norm", ["Raw counts", "Fraction of total", "Log1p(CPM)"])
def test_gene_correlation_executes_for_every_normalisation(norm, tmp_path):
    step = _gene_corr_step(norm, tmp_path=tmp_path)
    ex = _run([_normalize_step(), step])
    assert Path(step.params["path"]).exists()
    assert -1.0 <= ex.ns["pr"] <= 1.0
    plt.close("all")


def test_gene_correlation_filter_is_recorded_as_code(tmp_path):
    """The cluster filter used to be applied in the GUI and omitted entirely
    from the recorded cell, so the notebook correlated all cells."""
    step = _gene_corr_step(filtered=True, tmp_path=tmp_path)
    ex = _run([_normalize_step(), step])
    assert "isin(['0'])" in ex.graph.get(step.id).code
    assert len(ex.ns["x"]) == 60          # only the first population
    plt.close("all")


def test_gene_correlation_figure_is_the_one_the_viewer_shows(tmp_path):
    """The figure is an output of the step, not rebuilt alongside it."""
    ex = _run([_normalize_step(), _gene_corr_step(tmp_path=tmp_path)])
    assert ex.ns["fig"] is plt.gcf()
    plt.close("all")


# ── ROI DEG ──────────────────────────────────────────────────────────────────

def _roi_polygons():
    # Two boxes in napari (y, x) pixel order, one over each population.
    return [
        [[0.0, 200.0], [0.0, 320.0], [120.0, 320.0], [120.0, 200.0]],
        [[0.0, -10.0], [0.0, 110.0], [120.0, 110.0], [120.0, -10.0]],
    ]


def _rois_step():
    return Step(id="rois", template=_ROIS_TEMPLATE,
                params={"polygons": _roi_polygons()}, kind="setup",
                outputs=["roi_polygons"])


def _roi_deg_step(filtered=False):
    params = {"method": "wilcoxon", "pixel_size": 1.0}
    if filtered:
        params["clustering"] = CLUSTER_KEY
        params["selected"] = ["0", "1"]
    return Step(
        id="roi_deg", template=_roi_deg_template(filtered),
        params=params, deps=["rois"], outputs=["roi_deg_df", "roi_adata"],
    )


def test_roi_deg_no_longer_imports_the_viewer_package():
    """The old recorded cell called `xenium_viewer.utils.gene_analysis`, so the
    notebook was not standalone scverse code. It is plain shapely + scanpy now."""
    source = _roi_deg_template(False)
    assert "xenium_viewer" not in source
    assert "from shapely import contains_xy" in source
    assert "sc.tl.rank_genes_groups(" in source


def test_roi_deg_matches_compute_roi_deg():
    """The step must reproduce what the viewer has always actually computed."""
    from xenium_viewer.utils.gene_analysis import compute_roi_deg

    adata = _adata()
    polygons = [np.array(p) for p in _roi_polygons()]
    centroids_yx = adata.obsm["spatial"][:, ::-1] / 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy, _ = compute_roi_deg(adata, centroids_yx, polygons, 1.0,
                                    method="wilcoxon")
    ex = _run([_rois_step(), _roi_deg_step()], adata=_adata())
    pd.testing.assert_frame_equal(ex.ns["roi_deg_df"], legacy)


def test_roi_deg_replays_to_the_same_frame():
    ex = _run([_rois_step(), _roi_deg_step()])
    pd.testing.assert_frame_equal(ex.ns["roi_deg_df"], _replay(ex)["roi_deg_df"])


def test_roi_deg_cluster_filter_is_recorded_as_code():
    ex = _run([_rois_step(), _roi_deg_step(filtered=True)])
    code = ex.graph.get("roi_deg").code
    assert "cluster_mask = adata.obs['leiden_r1.0'].astype(str).isin(['0', '1'])" in code
    assert "_inside = _inside & cluster_mask" in code


def test_roi_polygons_round_trip_as_literals():
    ex = _run([_rois_step()])
    polygons = ex.ns["roi_polygons"]
    assert len(polygons) == 2
    np.testing.assert_allclose(polygons[0], np.array(_roi_polygons()[0]))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
