"""Every clustering the viewer shows must have a node whose code produces it.

Tabs that analyse a clustering declare ``deps=["clustering:<key>"]``, so a node
with that id has to exist. ``record_clustering`` is the backstop that creates
one when no producer has. It used to create it by pointing ``pd.read_csv`` at
``analysis/clustering/<key>/clusters.csv`` **whether or not that file existed** —
true only for the clusterings 10x ships, and false for every one the viewer
derives (Leiden, CNV, Novae, an import). The first replay of a real session
against its dataset died there with ``FileNotFoundError``, three cells in.

The fix is in two halves, and both are tested here:

- the producers record the code that actually made the column (source guard
  below), so the backstop is rarely reached at all;
- when it *is* reached — a column from a session recorded before its producer
  did — it emits a *reload* from the viewer's cache, says so in the cell, and,
  as asserted below, actually runs.

Qt is imported by the tab helpers, so these need
``QT_QPA_PLATFORM=offscreen``.
"""
from __future__ import annotations

import ast
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")
pytest.importorskip("qtpy")

SRC = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer"


@pytest.fixture
def ctx(tmp_path, qapp):
    """A ViewerContext with the shared helpers bound, and nothing else."""
    from xenium_viewer.tabs._helpers import create_shared_helpers
    from xenium_viewer.utils.viewer_context import ViewerContext

    context = ViewerContext(
        data_path=tmp_path, state={"record_code": True, "code_journal": []},
    )
    create_shared_helpers(context)
    return context


def _write_10x_clustering(data_path: Path, key: str, prefixed: bool = True):
    """Write an ``analysis/clustering`` CSV the way the Xenium output does."""
    name = f"gene_expression_{key}" if prefixed else key
    directory = data_path / "analysis" / "clustering" / name
    directory.mkdir(parents=True)
    csv = directory / "clusters.csv"
    pd.DataFrame({"Barcode": ["0", "1"], "Cluster": [1, 2]}).to_csv(csv, index=False)
    return csv


def _write_cache_with_clustering(data_path: Path, key: str, labels):
    """A minimal viewer cache carrying ``clustering_<key>`` in the table obs."""
    pytest.importorskip("spatialdata")
    from spatialdata import SpatialData
    from spatialdata.models import Labels2DModel, TableModel

    n_obs = len(labels)
    table = anndata.AnnData(np.ones((n_obs, 3), dtype="float32"))
    table.obs["region"] = pd.Categorical(["lab"] * n_obs)
    table.obs["instance_id"] = list(range(n_obs))
    table.obs[f"clustering_{key}"] = pd.Categorical([str(v) for v in labels])
    table = TableModel.parse(
        table, region="lab", region_key="region", instance_key="instance_id",
    )
    labels_elem = Labels2DModel.parse(np.arange(16, dtype=np.int32).reshape(4, 4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        SpatialData(labels={"lab": labels_elem}, tables={"table": table}).write(
            data_path / "sdata_cached.zarr",
        )


def _code(ctx, key: str) -> str:
    return ctx.state["prov_graph"].get(f"clustering:{key}").code


# ── the backstop ─────────────────────────────────────────────────────────────

def test_a_10x_clustering_is_recorded_as_the_csv_read_that_produced_it(ctx):
    """The dataset really does ship this file, so reading it *is* the real code."""
    csv = _write_10x_clustering(ctx.data_path, "graphclust")
    ctx.record_clustering("graphclust")

    code = _code(ctx, "graphclust")
    assert "pd.read_csv" in code
    assert str(csv) in code


def test_the_unprefixed_csv_layout_is_found_too(ctx):
    csv = _write_10x_clustering(ctx.data_path, "kmeans_2", prefixed=False)
    ctx.record_clustering("kmeans_2")
    assert str(csv) in _code(ctx, "kmeans_2")


def test_a_viewer_derived_clustering_never_records_a_read_of_a_missing_file(ctx):
    """Regression: the defect that broke the first real replay.

    No ``analysis/clustering/cnv_leiden_res0.5/`` exists — the CNV tab made this
    column — so recording a ``read_csv`` of that path produced a notebook that
    could not run.
    """
    ctx.record_clustering("cnv_leiden_res0.5")

    code = _code(ctx, "cnv_leiden_res0.5")
    assert "read_csv" not in code
    assert "analysis/clustering" not in code
    assert ast.parse(code).body, "the fallback must be executable, not a comment"


def test_the_reload_fallback_says_it_reloads_rather_than_recomputes(ctx):
    """An honest cell, on the CopyKAT precedent.

    Reloading stored labels cannot reproduce them, and a reader must be able to
    tell the two apart at a glance rather than by tracing where the data came from.
    """
    ctx.record_clustering("cnv_leiden_res0.5")
    code = _code(ctx, "cnv_leiden_res0.5").lower()
    assert "reload" in code
    assert "does not recompute" in code


def test_the_reload_fallback_actually_runs_against_the_cache(ctx):
    """Executability is the whole point, so it is executed rather than inspected."""
    labels = ["a", "b", "a", "b", "c", "c"]
    _write_cache_with_clustering(ctx.data_path, "cnv_leiden_res0.5", labels)
    ctx.record_clustering("cnv_leiden_res0.5")

    adata = anndata.AnnData(np.ones((len(labels), 2), dtype="float32"))
    namespace = {"pd": pd, "adata": adata}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(compile(_code(ctx, "cnv_leiden_res0.5"), "<recorded>", "exec"), namespace)  # noqa: S102

    assert list(namespace["adata"].obs["cnv_leiden_res0.5"]) == labels


def test_a_producers_code_is_never_overwritten_by_the_backstop(ctx):
    """The producer knows how the column was made; the backstop only guesses."""
    ctx.record_preamble()
    ctx.record_node("clustering:leiden_igraph_r1.0", "sc.tl.leiden(adata_leiden)",
                    deps=["preamble"])
    ctx.record_clustering("leiden_igraph_r1.0")

    assert _code(ctx, "leiden_igraph_r1.0") == "sc.tl.leiden(adata_leiden)"


def test_the_backstop_does_not_flag_dependents_stale_on_a_second_call(ctx):
    ctx.record_clustering("graphclust")
    ctx.record_node("rank_genes:graphclust", "sc.tl.rank_genes_groups(adata_norm)",
                    deps=["clustering:graphclust"])
    ctx.record_clustering("graphclust")

    assert not ctx.state["prov_graph"].get("rank_genes:graphclust").stale


# ── persisting the graph when it changes ─────────────────────────────────────

def _sidecar(ctx) -> Path:
    from xenium_viewer.tabs._helpers import PROV_GRAPH_SIDECAR
    from xenium_viewer.utils.adata_persistence import sidecar_dir
    return sidecar_dir(ctx.data_path) / PROV_GRAPH_SIDECAR


def test_the_graph_reaches_disk_as_soon_as_a_step_is_recorded(ctx):
    """Measured failure: the artifacts were persisted eagerly, the code lazily.

    ``save_clustering_to_adata`` writes the column immediately; the graph used
    to be written only by ``save_session`` (dataset switch / viewer exit). On a
    real session that gap left a 16-minute-old three-node graph on disk beside a
    table holding two Leiden clusterings and a rank-genes result — so the
    verification replayed a session the user had not run.
    """
    import json

    ctx.record_preamble()
    ctx.record_node("clustering:leiden_igraph_r1.0", "sc.tl.leiden(adata_leiden)",
                    deps=["preamble"])

    items = json.loads(_sidecar(ctx).read_text())
    assert [item["id"] for item in items] == ["preamble", "clustering:leiden_igraph_r1.0"]


def test_the_persisted_graph_reloads_into_an_equivalent_graph(ctx):
    import json
    from xenium_viewer.utils.prov_graph import ProvGraph

    ctx.record_preamble()
    ctx.record_node("clustering:k", "sc.tl.leiden(adata_leiden)", deps=["preamble"],
                    label="Clustering: k")

    restored = ProvGraph.from_list(json.loads(_sidecar(ctx).read_text()))
    assert restored.topo_sort() == ctx.state["prov_graph"].topo_sort()
    assert restored.get("clustering:k").deps == ["preamble"]
    assert restored.get("clustering:k").label == "Clustering: k"


def test_a_revised_node_is_persisted_too(ctx):
    """Re-running a step revises its node; disk must follow, not lag a version."""
    import json

    ctx.record_preamble()
    ctx.record_node("clustering:k", "sc.tl.leiden(adata_leiden, resolution=0.5)",
                    deps=["preamble"])
    ctx.record_node("clustering:k", "sc.tl.leiden(adata_leiden, resolution=1.0)",
                    deps=["preamble"])

    items = {item["id"]: item for item in json.loads(_sidecar(ctx).read_text())}
    assert "resolution=1.0" in items["clustering:k"]["code"]


def test_the_sidecar_wins_over_a_stale_session_attr():
    """Load precedence: the eagerly-written file beats the exit-written attr."""
    import json
    from xenium_viewer.app import _load_prov_graph_items
    from xenium_viewer.utils.adata_persistence import sidecar_dir
    from xenium_viewer.tabs._helpers import PROV_GRAPH_SIDECAR
    import tempfile

    data_path = Path(tempfile.mkdtemp())
    stale = [{"id": "preamble", "code": "import scanpy as sc", "kind": "setup"}]
    fresh = stale + [{"id": "clustering:k", "code": "sc.tl.leiden(adata_leiden)",
                      "deps": ["preamble"]}]
    directory = sidecar_dir(data_path, create=True)
    (directory / PROV_GRAPH_SIDECAR).write_text(json.dumps(fresh))

    items = _load_prov_graph_items(data_path, {"prov_graph": stale})
    assert [item["id"] for item in items] == ["preamble", "clustering:k"]


def test_a_dataset_without_the_sidecar_still_restores_from_the_attr():
    """Every cache written before this existed has only the attr."""
    from xenium_viewer.app import _load_prov_graph_items
    import tempfile

    stale = [{"id": "preamble", "code": "import scanpy as sc", "kind": "setup"}]
    items = _load_prov_graph_items(Path(tempfile.mkdtemp()), {"prov_graph": stale})
    assert [item["id"] for item in items] == ["preamble"]


def test_an_unreadable_sidecar_falls_back_rather_than_losing_the_graph():
    from xenium_viewer.app import _load_prov_graph_items
    from xenium_viewer.utils.adata_persistence import sidecar_dir
    from xenium_viewer.tabs._helpers import PROV_GRAPH_SIDECAR
    import tempfile

    data_path = Path(tempfile.mkdtemp())
    (sidecar_dir(data_path, create=True) / PROV_GRAPH_SIDECAR).write_text("{tru")
    attr = [{"id": "preamble", "code": "import scanpy as sc", "kind": "setup"}]

    assert _load_prov_graph_items(data_path, {"prov_graph": attr}) == attr


def test_the_sidecar_lives_outside_the_zarr_store(ctx):
    """Files inside the store make zarr's hierarchy walk warn, and a cache
    rebuild would delete them — the same rule the h5ad/parquet sidecars follow."""
    ctx.record_preamble()
    assert "sdata_cached.zarr" not in str(_sidecar(ctx))
    assert _sidecar(ctx).exists()


# ── the producers ────────────────────────────────────────────────────────────

def test_every_producer_of_a_clustering_column_records_a_clustering_node():
    """Source guard: persisting a clustering and recording it must travel together.

    ``save_clustering_to_adata`` is what puts a column into ``adata.obs`` and on
    disk. A module that calls it without recording a ``clustering:<key>`` node
    leaves the artifact visible in the GUI with no code behind it, which is
    exactly how the CNV and Novae columns ended up on the backstop path.
    """
    producers = [
        path for path in sorted((SRC / "tabs").glob("*.py"))
        if "save_clustering_to_adata(" in path.read_text()
    ]
    assert producers, "expected to find the tabs that persist clusterings"

    missing = [
        path.name for path in producers
        if "clustering:" not in path.read_text()
    ]
    assert missing == [], (
        f"{missing} persist a clustering column but record no 'clustering:<key>' "
        f"node, so nothing in the notebook produces it"
    )


def test_novae_records_the_column_name_the_viewer_actually_stores():
    """The recorded cell claimed ``novae_domain``; the viewer stores ``novae_domains``.

    A notebook that binds a different column than the one every downstream cell
    names fails at the dependent, not here — so it is pinned at the source.
    """
    source = (SRC / "tabs" / "tab_novae.py").read_text()
    assert 'clustering:novae_domains' in source
    assert "adata.obs['novae_domains'] = " in source


def test_cnv_publishes_its_clustering_onto_the_main_table():
    """Both backends leave labels on their own object; the column has to be moved.

    inferCNV binds ``cnv_clusters``; CopyKAT leaves them on the subsampled
    ``adata_copykat``. Publishing them is real code, and it is what the viewer
    does via ``save_clustering_to_adata``.
    """
    source = (SRC / "tabs" / "tab_cnv.py").read_text()
    assert "_record_cnv_clustering_node" in source
    assert "cnv_clusters.reindex(_ids)" in source      # inferCNV
    assert '_ck.reindex(_ids)' in source               # CopyKAT
    assert 'var = "adata_copykat"' in source


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
