"""Unit tests for notebook_export: graph -> .ipynb -> graph round-trip.

Requires `nbformat` (present in the runtime env). No napari/zarr.

Run with:  pytest tests/test_notebook_export.py
"""
from __future__ import annotations

from xenium_viewer.utils.prov_graph import ProvGraph, SETUP
from xenium_viewer.utils import notebook_export


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
