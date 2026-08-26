"""Where a figure goes: the plots directory, the file formats, the names.

Every figure the viewer produces is written to ``<data_path>/plots/`` in each
configured format. Before this module that was five hand-rolled
``os.makedirs(os.path.join(ctx.data_path, "plots"))`` calls, four different save
policies and one plot (the rank-genes panel) that was never written at all — so
"where did that figure go?" had no single answer, which is half of issue #35.

Deliberately pure Python: no Qt, no napari, no pyplot import at module level.
The figure argument is duck-typed on ``savefig`` so this can be unit-tested
headless, in the style of :mod:`xenium_viewer.utils.store_inventory`.
"""

from __future__ import annotations

import re
from pathlib import Path

#: Written for every plot unless the user narrows it in Preferences → Plot format.
#: Two formats rather than one because the two audiences differ: a PNG to look at
#: and paste into a slide, a PDF to hand to a journal.
DEFAULT_FORMATS = ("png", "pdf")

#: What Preferences offers, and the only values ``plot_formats`` will return.
KNOWN_FORMATS = ("png", "pdf", "svg")

#: Directory name under the dataset. Listed in
#: ``store_inventory.DERIVED_IN_DATA_DIR`` as viewer-created and therefore
#: deletable, so it must stay in step with that entry.
PLOTS_DIRNAME = "plots"

#: Raster dpi. PDF/SVG ignore it; passing it unconditionally keeps one call.
SAVE_DPI = 300


def plot_formats(state) -> list[str]:
    """The formats to write, from the viewer's state dict.

    Unknown entries are dropped rather than passed to ``savefig``, which would
    raise from inside a worker thread and lose the figure. An empty or missing
    setting falls back to the default rather than writing nothing: a plot the
    user asked for and cannot find afterwards is the failure mode this whole
    module exists to remove.
    """
    requested = state.get("plot_formats") if state else None
    if isinstance(requested, str):          # tolerate the old scalar setting
        requested = [requested]
    formats = [f for f in (requested or ()) if f in KNOWN_FORMATS]
    return formats or list(DEFAULT_FORMATS)


def primary_format(state) -> str:
    """The format a recorded code snippet names when it can only name one."""
    return plot_formats(state)[0]


def plots_dir(data_path, create: bool = False) -> Path:
    """``<data_path>/plots``. The one definition of it.

    Mirrors :func:`adata_persistence.sidecar_dir`, including the ``create``
    flag, so a reader who knows one knows the other.
    """
    directory = Path(data_path) / PLOTS_DIRNAME
    if create:
        directory.mkdir(parents=True, exist_ok=True)
    return directory


def safe_stem(text: str) -> str:
    """Make *text* safe to embed in a filename.

    A clustering key reaches a filename verbatim, and those carry spaces,
    slashes and parentheses (``leiden (res 1.0)``). This is the rule
    ``tab_cnv`` already applied to its CNV heatmap names, lifted here so every
    plot name is built the same way.
    """
    cleaned = re.sub(r"[^0-9A-Za-z._-]+", "_", str(text)).strip("_")
    return cleaned or "plot"


def save_paths(data_path, stem: str, formats=None, state=None) -> list[Path]:
    """The files ``stem`` would be written to, in order.

    Pass *formats* explicitly, or *state* to read the user's preference.
    """
    if formats is None:
        formats = plot_formats(state)
    directory = plots_dir(data_path)
    return [directory / f"{stem}.{ext}" for ext in formats]


def batch_dir(data_path, name: str) -> Path:
    """``<data_path>/plots/<name>`` — where a *set* of figures goes.

    The three pairwise-volcano generators produce N×N figures at once. Those
    stay out of the Plots gallery (they would bury every other plot) but they
    belong in the same tree as everything else, so the default the directory
    chooser opens on is here rather than the user's home.
    """
    return plots_dir(data_path) / safe_stem(name)


def save_figure(fig, paths, dpi: int = SAVE_DPI) -> list[str]:
    """Write *fig* to every path in *paths*, creating their directory.

    Returns the paths as strings, for the status line and for the recorded
    ``savefig`` call — which quotes what was actually written rather than a
    bare relative guess, as every hand-written ``plot:*`` node used to.
    """
    written = []
    for path in paths:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        written.append(str(path))
    return written


def recorded_save_code(paths, fig_expr: str = "fig", dpi: int = SAVE_DPI) -> str:
    """The lines a recorded cell needs to write *paths*, ``plots/`` included.

    One implementation for the half-dozen tabs that still build their cell as a
    string. The ``mkdir`` is not decoration: the paths are relative to the
    dataset directory, and a replayed notebook has no ``plots/`` until something
    makes one — ``savefig`` does not.

    ``Path`` is bound by the preamble and declared in ``EXECUTOR_BASE_NAMES``,
    so the line is valid in both the executor and the exported notebook.
    """
    lines = []
    for path in paths:
        text = str(path)
        lines.append(f"Path({text!r}).parent.mkdir(parents=True, exist_ok=True)")
        lines.append(
            f"{fig_expr}.savefig({text!r}, dpi={dpi}, bbox_inches='tight')")
    return "\n".join(lines)


def recorded_paths(data_path, paths) -> list[str]:
    """The same files as a notebook should name them: relative to the dataset.

    An exported notebook replays from the raw Xenium output with ``data_path``
    bound, so ``plots/dotplot.png`` resolves to the same file the GUI wrote
    while an absolute path would pin the notebook to one machine.
    """
    root = Path(data_path)
    out = []
    for path in paths:
        path = Path(path)
        try:
            out.append(str(path.relative_to(root)))
        except ValueError:
            out.append(str(path))
    return out
