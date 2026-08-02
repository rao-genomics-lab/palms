"""
Export the reproducible-code provenance graph as a Jupyter notebook (.ipynb).

The notebook body is *derived* from the provenance graph (topologically sorted,
one cell per node — see ``prov_graph.graph_to_cells``), so the exported notebook
always respects dependencies regardless of the order actions were recorded. The
notebook is code-only (no stored outputs); outputs regenerate when it is run.

``nbformat`` is imported lazily so importing this module never hard-requires it.

:func:`execute_notebook` runs an exported notebook in a clean kernel. It is the
measurement half of the reproducibility claim — the viewer can only assert that
the recorded code *is* the executed code; whether replaying it reproduces the
result has to be run to be known (``tests/test_notebook_replay.py`` and
``scripts/verify_notebook.py``).
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterable, Optional


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


def customisation_banner(graph) -> str | None:
    """A markdown note naming steps that did not use the shipped template.

    ``code`` always records what actually ran, so the notebook replays correctly
    whether or not a template was customised. What the *reader* cannot tell from
    the source — a customised template renders to code that looks entirely
    ordinary — is that it is not the stock pipeline. Someone opening this file a
    year later, or reviewing it alongside a paper, needs that stated once at the
    top rather than inferred.

    Returns ``None`` for a fully stock run, so the ordinary notebook is not
    cluttered with a banner saying nothing happened.
    """
    from xenium_viewer.utils.prov_graph import TEMPLATE_BUILTIN

    rows = []
    for node_id in graph.topo_sort():
        node = graph.get(node_id)
        if node is None or node.template_origin == TEMPLATE_BUILTIN:
            continue
        detail = node.template_id or "—"
        short = (node.template_hash or "")[:12]
        rows.append(f"| `{node_id}` | `{detail}` | {node.template_origin} | `{short}` |")
    if not rows:
        return None

    return "\n".join([
        "## ⚠ This analysis did not use the shipped templates",
        "",
        "The cells below are the code that actually ran, so this notebook "
        "replays exactly what the viewer did. But some steps were generated "
        "from **customised** analysis templates rather than the ones shipped "
        "with Xenium Viewer, which the source alone cannot show.",
        "",
        "| step | template | origin | template hash |",
        "|---|---|---|---|",
        *rows,
        "",
        "`user` / `user+builtin` mean a customised template; `hand-edited` "
        "means the cell was typed into the Notebook tab and **not re-executed**, "
        "so it is the only kind here that may not describe what produced the "
        "result.",
    ])


def graph_to_cells(graph, include_terminals: bool = True) -> list[tuple]:
    """Topologically-ordered [(cell_type, source), ...] derived from the graph.

    A customisation banner is prepended when any step used a non-shipped
    template. It goes here rather than in ``prov_graph.graph_to_cells`` so the
    flat ``analysis.py`` is unaffected, and it is a single markdown cell so the
    verbatim-source property the replay test asserts over *code* cells still
    holds.
    """
    from xenium_viewer.utils.prov_graph import graph_to_cells as _g2c
    cells = [(c.cell_type, c.source)
             for c in _g2c(graph, include_terminals=include_terminals)]
    banner = customisation_banner(graph)
    return ([("markdown", banner)] + cells) if banner else cells


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


# ── execution ────────────────────────────────────────────────────────────────

@contextmanager
def _own_interpreter_kernel(name: str = "xv-replay"):
    """Yield a kernel name that launches *this* interpreter.

    The installed ``python3`` kernelspec belongs to whichever environment
    registered it last — on a conda/micromamba box that is routinely the base
    env, which has no scanpy. Replaying into it would report a failure that is
    an artifact of kernel discovery rather than of the recorded code. So a
    throwaway kernelspec pointing at ``sys.executable`` is written to a temp
    dir and prepended to ``JUPYTER_PATH`` for the duration.
    """
    with tempfile.TemporaryDirectory(prefix="xv-kernel-") as tmp:
        spec_dir = Path(tmp) / "kernels" / name
        spec_dir.mkdir(parents=True)
        (spec_dir / "kernel.json").write_text(json.dumps({
            "argv": [sys.executable, "-m", "ipykernel_launcher",
                     "-f", "{connection_file}"],
            "display_name": f"Python ({Path(sys.prefix).name})",
            "language": "python",
        }))
        previous = os.environ.get("JUPYTER_PATH")
        os.environ["JUPYTER_PATH"] = (
            f"{tmp}{os.pathsep}{previous}" if previous else tmp
        )
        try:
            yield name
        finally:
            if previous is None:
                os.environ.pop("JUPYTER_PATH", None)
            else:
                os.environ["JUPYTER_PATH"] = previous


def execute_notebook(
    path: str | Path,
    cwd: Optional[str | Path] = None,
    timeout: int = 1800,
    on_cell_start: Optional[Callable] = None,
    on_cell_executed: Optional[Callable] = None,
    on_cell_error: Optional[Callable] = None,
):
    """Execute the notebook at *path* in a fresh kernel and return it.

    Errors are *not* allowed: the first failing cell raises
    ``nbclient.exceptions.CellExecutionError``. A notebook that cannot run is a
    reproducibility failure, so it must not be swallowed into a stored traceback.

    Note for callers wiring up the hooks: nbclient calls ``on_cell_executed``
    for a *failing* cell as well, just before it raises — so treating that hook
    as "this cell succeeded" mis-reports the one cell you most want named. Use
    ``on_cell_error`` for the failure.

    The executed notebook (with outputs) is returned rather than written back;
    callers that want it on disk write it themselves.
    """
    import nbformat
    from nbclient import NotebookClient

    path = Path(path)
    nb = nbformat.read(str(path), as_version=4)
    run_dir = Path(cwd) if cwd is not None else path.parent
    with _own_interpreter_kernel() as kernel_name:
        client = NotebookClient(
            nb,
            kernel_name=kernel_name,
            timeout=timeout,
            allow_errors=False,
            resources={"metadata": {"path": str(run_dir)}},
            on_cell_start=on_cell_start,
            on_cell_executed=on_cell_executed,
            on_cell_error=on_cell_error,
        )
        client.execute()
    return nb
