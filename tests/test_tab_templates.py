"""The Templates tab shows what a step would run, before it runs.

The value of this view rests entirely on the preview being the *real* rendered
source rather than a reconstruction of it — a preview assembled by a second code
path would be exactly the drift the Step system exists to eliminate, just moved
somewhere less visible. So the tests here check that the preview equals
``Step.render()`` of the same template and params, and that a tab supplying live
widget values is actually consulted.

The tab is built against a stand-in context: the real ``ViewerContext`` needs a
napari viewer, and this view reads only ``ctx.state``.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from xenium_viewer.utils.step_templates import builtin_ids, builtin_spec
from xenium_viewer.utils.steps import Step


@pytest.fixture
def ctx():
    return SimpleNamespace(state={})


@pytest.fixture
def tab(qapp, ctx):
    from xenium_viewer.tabs.tab_templates import build_tab
    widget, exports = build_tab(ctx)
    return widget, exports


def _tree(widget):
    from qtpy.QtWidgets import QTreeWidget
    return widget.findChild(QTreeWidget)


def test_every_registered_template_is_listed(tab):
    from qtpy.QtCore import Qt

    widget, _ = tab
    tree = _tree(widget)
    listed = set()
    for i in range(tree.topLevelItemCount()):
        group = tree.topLevelItem(i)
        for j in range(group.childCount()):
            listed.add(group.child(j).data(0, Qt.UserRole))
    assert listed == set(builtin_ids())


def test_templates_are_grouped_by_owner(tab):
    widget, _ = tab
    tree = _tree(widget)
    groups = {tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())}
    assert {"Clustering", "Genes", "Spatial", "ROI", "Setup"} == groups


def test_the_preview_is_the_real_rendered_source(ctx):
    """Not a reconstruction: identical to what StepExecutor would compile."""
    from xenium_viewer.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    assembly = spec.assemblies[0]
    shown = _preview(ctx, spec, assembly)

    expected = Step(id="x", template=spec.assemble(assembly),
                    params=spec.synth_params()).render()
    assert shown.endswith(expected)
    assert shown.startswith("# preview — sample values")


def test_a_live_tab_supplies_the_real_widget_values(ctx):
    """When the owning tab registers a provider, the preview uses it."""
    from xenium_viewer.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview_params"] = {
        "clustering.leiden": lambda: dict(spec.synth_params(), resolution=0.42),
    }
    shown = _preview(ctx, spec, spec.assemblies[0])

    assert "resolution=0.42" in shown
    assert shown.startswith("# preview — current widget values")


def test_a_broken_provider_does_not_break_the_view(ctx):
    """A half-built tab must degrade to sample values, not raise into the GUI."""
    from xenium_viewer.tabs.tab_templates import _preview

    spec = builtin_spec("clustering.leiden")
    ctx.state["template_preview_params"] = {
        "clustering.leiden": lambda: (_ for _ in ()).throw(RuntimeError("not ready")),
    }
    shown = _preview(ctx, spec, spec.assemblies[0])
    assert "sc.tl.leiden(" in shown


def test_the_clustering_tab_registers_a_preview_provider():
    """Source guard: the wiring that makes the preview real, not illustrative.

    ``_leiden_params`` is the single expression of "the current settings" — the
    run and the preview both call it. If the tab stopped registering it the
    preview would silently fall back to sample values and look fine.
    """
    import ast
    import inspect

    from xenium_viewer.tabs import tab_clustering

    source = inspect.getsource(tab_clustering)
    assert "template_preview_params" in source, (
        "tab_clustering no longer registers its preview provider"
    )
    tree = ast.parse(source)
    names = {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert "_leiden_params" in names


def test_the_contract_block_reports_the_declared_interface(tab):
    from xenium_viewer.tabs.tab_templates import _contract_text

    spec = builtin_spec("roi.expression")
    text = _contract_text(spec)
    assert "roi_expr_cells" in text          # a declared output
    assert "pixel_size" in text              # a required param
    assert f"blocks    : {', '.join(spec.blocks)}" in text


def test_a_frozen_block_says_so(tab):
    """The Arrow shim is not editable; the view has to make that visible."""
    from xenium_viewer.tabs.tab_templates import _contract_text

    text = _contract_text(builtin_spec("genes.cnv_infercnv"))
    assert "frozen    : arrow_shim" in text


def test_selecting_a_template_fills_all_three_panes(tab):
    from qtpy.QtWidgets import QPlainTextEdit

    widget, _ = tab
    tree = _tree(widget)
    group = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
                 if tree.topLevelItem(i).text(0) == "Clustering")
    tree.setCurrentItem(group.child(0))

    panes = widget.findChildren(QPlainTextEdit)
    filled = [p for p in panes if p.toPlainText().strip()]
    assert len(filled) == 3, "contract, source and preview should all populate"
    assert all(p.isReadOnly() for p in panes), "this view is read-only"
