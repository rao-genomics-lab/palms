"""
Export the reproducible-code provenance graph as a Jupyter notebook (.ipynb).

The notebook body is *derived* from the provenance graph (topologically sorted,
one cell per node — see ``prov_graph.graph_to_cells``), so the exported notebook
always respects dependencies regardless of the order actions were recorded. The
notebook is code-only (no stored outputs); outputs regenerate when it is run.

``nbformat`` is imported lazily so importing this module never hard-requires it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional


_KERNELSPEC = {"name": "python3", "display_name": "Python 3", "language": "python"}


def _build_notebook(cells: Iterable[tuple]):
    """Build an nbformat NotebookNode from (cell_type, source) pairs."""
    import nbformat
    from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell

    nb = new_notebook()
    nb.metadata["kernelspec"] = dict(_KERNELSPEC)
    for cell_type, source in cells:
        src = (source or "").strip("\n")
        if cell_type == "markdown":
            nb.cells.append(new_markdown_cell(src))
        else:
            nb.cells.append(new_code_cell(src))
    return nb


def graph_to_cells(graph, include_terminals: bool = True) -> list[tuple]:
    """Topologically-ordered [(cell_type, source), ...] derived from the graph."""
    from xenium_viewer.utils.prov_graph import graph_to_cells as _g2c
    return [(c.cell_type, c.source) for c in _g2c(graph, include_terminals=include_terminals)]


def write_notebook(cells: Iterable[tuple], path: str | Path) -> None:
    """Write an explicit list of (cell_type, source) cells to *path* as .ipynb.

    Used by the Notebook tab, which exports graph-derived cells plus any
    user-authored cells.
    """
    import nbformat
    nb = _build_notebook(cells)
    nbformat.write(nb, str(path))


def write_graph_notebook(graph, path: str | Path,
                         include_terminals: bool = True) -> None:
    """Convenience: derive cells from *graph* and write them to *path*."""
    write_notebook(graph_to_cells(graph, include_terminals=include_terminals), path)


def read_notebook(path: str | Path) -> list[tuple]:
    """Read a .ipynb back into a list of (cell_type, source) pairs."""
    import nbformat
    nb = nbformat.read(str(path), as_version=4)
    return [(c.cell_type, c.source) for c in nb.cells]
