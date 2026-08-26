"""One route on screen, one route to disk — enforced in the source.

Issue #35 was not one bug. Eighteen sites produced a figure and no two agreed:
four display modes (floating ``plt.show(block=False)``, *blocking* ``plt.show()``,
invisible ``Agg``-only, one embedded ``QPixmap``) and four save policies. It got
that way one tab at a time, and it would again, so the invariant is checked
statically rather than remembered:

    A tab may not display a figure itself, may not set the matplotlib backend,
    and may not write one to disk except through ``ctx.show_plot`` or its Step
    template.

The ``matplotlib.use`` half is the one that bit hardest. Two workers called it
globally and never restored it, so after a user saved a UMAP or ran any marker
plot, every later ``plt.show`` in the session became a silent no-op — figures
from unrelated tabs simply stopped appearing, with no error anywhere.

Pure ast: no Qt, no imports of the tabs.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

TABS = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer" / "tabs"
UTILS = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer" / "utils"

#: The modules allowed to call ``savefig``, each with the reason.
_SAVE_ALLOWED = {
    # What ctx.show_plot uses; the definition of the policy.
    "plot_output.py",
    # The Plots dock's own "Save as…" button — a copy to a place the *user*
    # picked, which is display-side, not a second output policy.
    "plots_panel.py",
}

def _tab_sources():
    return sorted(p for p in TABS.glob("tab_*.py"))


def _calls(tree: ast.AST):
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            yield node


def _dotted(func: ast.AST) -> str:
    """``plt.show`` for an Attribute chain, ``show`` for a bare Name."""
    parts = []
    while isinstance(func, ast.Attribute):
        parts.append(func.attr)
        func = func.value
    if isinstance(func, ast.Name):
        parts.append(func.id)
    return ".".join(reversed(parts))


@pytest.mark.parametrize("path", _tab_sources(), ids=lambda p: p.name)
def test_no_tab_displays_a_figure_itself(path: Path):
    """``plt.show`` is what made display depend on process-wide backend state.

    Display goes through the Plots dock, which draws into its own
    ``FigureCanvasQTAgg`` — so it works regardless of what pyplot is pointed at.
    """
    tree = ast.parse(path.read_text())
    offenders = [
        node.lineno for node in _calls(tree)
        if _dotted(node.func).split(".")[-1] == "show"
        and _dotted(node.func).split(".")[0] in {"plt", "_plt", "pyplot"}
    ]
    assert not offenders, (
        f"{path.name} calls plt.show at line(s) {offenders} — use "
        f"ctx.show_plot(fig, stem) so the figure lands in the Plots dock"
    )


@pytest.mark.parametrize("path", _tab_sources(), ids=lambda p: p.name)
def test_no_tab_sets_the_matplotlib_backend(path: Path):
    """Setting it is process-wide and was never restored.

    ``cnv_copykat_worker.py`` still calls it, correctly: that is a separate
    detached process with no GUI, so it is not a tab and not covered here.
    """
    tree = ast.parse(path.read_text())
    offenders = [node.lineno for node in _calls(tree)
                 if _dotted(node.func) in {"matplotlib.use", "mpl.use"}]
    assert not offenders, (
        f"{path.name} calls matplotlib.use at line(s) {offenders} — it is "
        f"process-wide and silently disables every later plt.show in the session"
    )


@pytest.mark.parametrize("path", _tab_sources(), ids=lambda p: p.name)
def test_a_tab_saves_a_figure_only_through_the_shared_helper(path: Path):
    """Every ``savefig`` in a tab must come from ``plot_output.save_figure``.

    Hand-rolled saves are how four different policies arose, and how the path a
    ``plot:*`` node recorded came to name a file that was never written.
    Recorded *source strings* are exempt — those are code for the notebook to
    run, not a call this module makes — which the ast walk gets for free, since
    a string constant holds no Call node.
    """
    tree = ast.parse(path.read_text())
    offenders = [node.lineno for node in _calls(tree)
                 if _dotted(node.func).split(".")[-1] == "savefig"]
    assert not offenders, (
        f"{path.name} calls savefig directly at line(s) {offenders} — use "
        f"ctx.show_plot(fig, stem), or plot_output.save_figure for a batch"
    )


def test_only_plot_output_writes_a_figure_among_the_utils():
    """The other half of the same rule, for the modules the tabs call into.

    A figure *factory* returns a Figure and leaves saving to its caller; that is
    the contract ``make_cnv_heatmap`` already documented. This asserts the rest
    of ``utils`` keeps it.
    """
    offenders = {}
    for path in sorted(UTILS.glob("*.py")):
        if path.name in _SAVE_ALLOWED:
            continue
        lines = [node.lineno for node in _calls(ast.parse(path.read_text()))
                 if _dotted(node.func).split(".")[-1] == "savefig"]
        if lines:
            offenders[path.name] = lines
    assert not offenders, (
        f"savefig outside utils/plot_output.py: {offenders} — a figure factory "
        f"returns the Figure and lets ctx.show_plot write it"
    )


def test_the_plots_directory_has_exactly_one_definition():
    """Five modules used to build ``<data_path>/plots`` by hand."""
    package = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer"
    offenders = {}
    for path in sorted(package.rglob("*.py")):
        if path.name in {
            "plot_output.py",       # the definition itself
            "cnv_copykat_worker.py",  # a detached process; cannot import ctx
            "rename_dataset.py",    # names a CopyKAT sentinel, not a plot path
        }:
            continue
        text = path.read_text()
        hits = [i + 1 for i, line in enumerate(text.splitlines())
                if ('/ "plots"' in line or "'plots'" in line
                    or '"plots")' in line)
                and not line.lstrip().startswith("#")]
        if hits:
            offenders[str(path.relative_to(package))] = hits
    assert not offenders, (
        f"plots directory built by hand in {offenders} — use "
        f"plot_output.plots_dir(data_path)"
    )
