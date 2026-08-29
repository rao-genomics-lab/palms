"""Unit tests for notebook_export: graph -> .ipynb -> graph round-trip.

Requires `nbformat` (present in the runtime env). No napari/zarr.

Run with:  pytest tests/test_notebook_export.py
"""
from __future__ import annotations

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
