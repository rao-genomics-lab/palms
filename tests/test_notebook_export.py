"""Unit tests for notebook_export: graph -> .ipynb -> graph round-trip.

Requires `nbformat` (present in the runtime env). No napari/zarr.

Run with:  pytest tests/test_notebook_export.py
"""
from __future__ import annotations

import pytest

from palms.utils.prov_graph import ProvGraph, SETUP
from palms.utils import notebook_export


def _graph():
    g = ProvGraph()
    g.upsert("preamble", "import scanpy as sc", kind=SETUP, label="Setup")
    g.upsert("normalize", "sc.pp.normalize_total(adata)", deps=["preamble"],
             label="Normalize")
    return g


def test_write_read_roundtrip_preserves_cells(tmp_path):
    g = _graph()
    expected = notebook_export.graph_to_cells(g)  # [(cell_type, source), ...]
    path = tmp_path / "analysis.ipynb"

    notebook_export.write_graph_notebook(g, path)
    assert path.exists()

    got = notebook_export.read_notebook(path)
    # _build_notebook strips leading/trailing newlines; our sources have none.
    assert got == [(t, (s or "").strip("\n")) for t, s in expected]

    # A labeled node yields a markdown header directly before its code cell.
    assert got[0][0] == "markdown" and "Setup" in got[0][1]
    assert ("code", "import scanpy as sc") in got


def test_written_notebook_is_valid_nbformat(tmp_path):
    import nbformat

    g = _graph()
    path = tmp_path / "nb.ipynb"
    notebook_export.write_graph_notebook(g, path)

    nb = nbformat.read(str(path), as_version=4)
    nbformat.validate(nb)  # raises on an invalid notebook
    assert nb.metadata["kernelspec"]["name"] == "python3"


# ── the verifier's notebook ──────────────────────────────────────────────────

def _verify_notebook_module():
    """Load ``scripts/verify_notebook.py``, which is a script, not a package."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parent.parent / "scripts" / "verify_notebook.py"
    spec = importlib.util.spec_from_file_location("_verify_notebook", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _customised_graph():
    from palms.utils.prov_graph import TEMPLATE_BLENDED

    g = _graph()
    g.upsert("clustering:leiden_x", "sc.tl.leiden(adata, resolution=1.0 * 2)",
             deps=["normalize"], label="Clustering: leiden_x",
             template_id="clustering.leiden", template_origin=TEMPLATE_BLENDED,
             template_hash="deadbeef")
    return g


def test_the_verifiers_notebook_carries_the_customisation_banner(tmp_path):
    """The verifier writes a notebook too, and it must not be the quiet one.

    ``build_notebook`` takes its cells from ``prov_graph.graph_to_cells``,
    because that is the one carrying ``node_id`` — which is what turns
    nbclient's "cell 4 failed" into a named step. But the banner is added by
    ``notebook_export.graph_to_cells``, so taking the other path silently
    dropped it: a customised session's ``--work-dir`` notebook was the single
    artifact that did not say it was customised, while the viewer's own export
    and the ``--out`` report both did. It is the artifact most likely to be kept
    and forwarded.
    """
    vn = _verify_notebook_module()
    graph = _customised_graph()

    nb_path, node_ids = vn.build_notebook(graph, tmp_path)
    cells = notebook_export.read_notebook(nb_path)

    kind, source = cells[0]
    assert kind == "markdown"
    assert "did not use the shipped templates" in source
    assert "clustering:leiden_x" in source

    # node_ids indexes cells, so a cell inserted at the front must shift it.
    assert len(node_ids) == len(cells), "node_ids no longer lines up with cells"
    assert node_ids[0] is None, "the banner belongs to no node"
    assert node_ids[-1] is None, "the injected dump cell belongs to no node"
    assert "clustering:leiden_x" in node_ids


def test_a_stock_run_gets_no_banner_from_the_verifier(tmp_path):
    """The ordinary case stays uncluttered, and the alignment still holds."""
    vn = _verify_notebook_module()

    nb_path, node_ids = vn.build_notebook(_graph(), tmp_path)
    cells = notebook_export.read_notebook(nb_path)

    assert "did not use the shipped templates" not in cells[0][1]
    assert len(node_ids) == len(cells)
    assert node_ids[-1] is None


# ── 10x-native clusterings are inputs, not results ───────────────────────────
# A clustering 10x shipped with the dataset is not produced by any recorded
# step, so a replayed notebook is not expected to reproduce it and its absence
# is not a divergence. A clustering the *viewer* persisted with no step behind
# it still is: that is the defect tests/test_clustering_recording.py guards
# against, and the two are indistinguishable in obs — both are
# `clustering_<key>` columns.

def _obs_frame(pd, **columns):
    return pd.DataFrame(columns, index=[f"cell{i}" for i in range(4)])


def test_native_clusterings_are_read_from_disk_when_the_dataset_has_them(tmp_path):
    from palms import loader

    root = tmp_path / "analysis" / "clustering"
    (root / "gene_expression_graphclust").mkdir(parents=True)
    (root / "gene_expression_graphclust" / "clusters.csv").write_text("Barcode,Cluster\n")
    (root / "gene_expression_my_own").mkdir(parents=True)
    (root / "gene_expression_my_own" / "clusters.csv").write_text("Barcode,Cluster\n")

    assert loader.native_clustering_names(tmp_path) == {"graphclust", "my_own"}
    # Disk wins over the naming rule in both directions.
    assert loader.is_native_clustering("my_own", tmp_path) is True
    assert loader.is_native_clustering("kmeans_5_clusters", tmp_path) is False


def test_the_name_decides_only_when_the_dataset_has_no_analysis_folder(tmp_path):
    """A crop export drops analysis/ but keeps the obs columns.

    Refusing to answer there would make every crop export look as though it had
    lost ten analyses.
    """
    from palms import loader

    assert loader.native_clustering_names(tmp_path) == set()
    assert loader.is_native_clustering("graphclust", tmp_path) is True
    assert loader.is_native_clustering("kmeans_10_clusters", tmp_path) is True
    assert loader.is_native_clustering("leiden_igraph_r1.0", tmp_path) is False
    assert loader.is_native_clustering("cnv_leiden_res0.2", tmp_path) is False
    # Near-misses must not be swept in.
    assert loader.is_native_clustering("my_kmeans", tmp_path) is False
    assert loader.is_native_clustering("kmeans_clusters", tmp_path) is False


def test_a_native_clustering_absent_from_the_replay_is_an_input_not_a_failure(tmp_path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("sklearn")
    vn = _verify_notebook_module()

    viewer = _obs_frame(pd, **{
        "clustering_graphclust": ["1", "2", "1", "2"],
        "clustering_leiden_r1.0": ["a", "b", "a", "b"],
    })
    replay = _obs_frame(pd, **{"leiden_r1.0": ["a", "b", "a", "b"]})

    by_key = {e["clustering"]: e for e in
              vn.compare_clusterings(viewer, replay, tmp_path)}
    assert by_key["graphclust"]["status"] == "input"
    assert by_key["leiden_r1.0"]["status"] == "ok"
    assert by_key["leiden_r1.0"]["ari"] == 1.0


def test_a_viewer_derived_clustering_with_no_step_is_still_a_failure(tmp_path):
    """The guard that must survive treating native clusterings as inputs."""
    pd = pytest.importorskip("pandas")
    pytest.importorskip("sklearn")
    vn = _verify_notebook_module()

    viewer = _obs_frame(pd, **{"clustering_leiden_r1.0": ["a", "b", "a", "b"]})
    replay = _obs_frame(pd, **{"unrelated": ["x", "y", "x", "y"]})

    entry = vn.compare_clusterings(viewer, replay, tmp_path)[0]
    assert entry["status"] == "not_in_replay"
    assert entry["status"] not in ("ok", "input"), "must still count as a failure"
