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


def _select(widget, group_name="Clustering", index=0):
    tree = _tree(widget)
    group = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
                 if tree.topLevelItem(i).text(0) == group_name)
    tree.setCurrentItem(group.child(index))
    return group.child(index)


def test_selecting_a_template_fills_every_pane(tab):
    from qtpy.QtWidgets import QPlainTextEdit

    widget, _ = tab
    _select(widget)

    panes = widget.findChildren(QPlainTextEdit)
    assert len(panes) == 4, "contract, default, yours, preview"
    assert all(p.toPlainText().strip() for p in panes)


def test_exactly_one_pane_is_editable(tab):
    """The default must stay read-only, or there is no reference to diff against."""
    from qtpy.QtWidgets import QPlainTextEdit

    widget, _ = tab
    _select(widget)
    panes = widget.findChildren(QPlainTextEdit)
    writable = [p for p in panes if not p.isReadOnly()]
    assert len(writable) == 1, "only 'Yours' may be edited"


# ── editing ──────────────────────────────────────────────────────────────────

@pytest.fixture
def editable(qapp, tmp_path, monkeypatch):
    """A tab whose overrides land in tmp_path rather than the real config dir.

    One environment variable is enough because saving resolves its destination
    from the same search path reading uses — a write cannot land somewhere the
    reader does not look. Before that, these tests wrote into the developer's
    actual ``~/.config``.
    """
    from xenium_viewer.utils.step_templates import (
        TEMPLATE_PATH_ENV, clear_cache, user_template_dir,
    )

    monkeypatch.setenv(TEMPLATE_PATH_ENV, str(tmp_path))
    clear_cache()
    assert user_template_dir() == tmp_path, "saves must be redirected too"

    from xenium_viewer.tabs.tab_templates import build_tab
    widget, _ = build_tab(SimpleNamespace(state={}))
    yield widget, tmp_path
    clear_cache()


def _panes(widget):
    from qtpy.QtWidgets import QPlainTextEdit
    return widget.findChildren(QPlainTextEdit)


def _editor(widget):
    return next(p for p in _panes(widget) if not p.isReadOnly())


def _button(widget, text):
    from qtpy.QtWidgets import QPushButton
    return next(b for b in widget.findChildren(QPushButton)
                if b.text().replace("&", "") == text)


def _problems(widget):
    from qtpy.QtWidgets import QListWidget
    lst = widget.findChild(QListWidget)
    return [lst.item(i).text() for i in range(lst.count())]


@pytest.mark.parametrize("template_id", builtin_ids())
def test_the_editor_round_trips_every_template_body(template_id):
    """Parsing an edit back must agree with how the loader splits a file.

    Every template, because the bug this caught was specific to one: an
    annotation on the marker line of a frozen block came back as a *different*
    block name, so editing the CNV Arrow shim silently created a new block
    instead of modifying the protected one.
    """
    from xenium_viewer.tabs.tab_templates import blocks_to_text, text_to_blocks

    spec = builtin_spec(template_id)
    assert text_to_blocks(blocks_to_text(spec)) == {
        name: block.text for name, block in spec.blocks.items()
    }


def test_validate_reports_a_broken_edit_without_saving(editable):
    widget, tmp_path = editable
    _select(widget)
    _editor(widget).setPlainText("#--- block head\nadata_leiden = (")

    _button(widget, "Validate").click()
    assert any("not valid Python" in p for p in _problems(widget))
    assert list(tmp_path.glob("*.tmpl")) == [], "Validate must not write anything"


def test_saving_a_valid_edit_activates_it(editable):
    from xenium_viewer.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    spec = builtin_spec("clustering.leiden")
    from xenium_viewer.tabs.tab_templates import blocks_to_text
    edited = blocks_to_text(spec).replace("max_value=10", "max_value=99")
    _editor(widget).setPlainText(edited)

    _button(widget, "Save  Activate").click()

    assert (tmp_path / "clustering.leiden.tmpl").is_file()
    resolved = resolve("clustering.leiden")
    assert resolved.is_customised
    assert resolved.changed_blocks() == ["scale"], (
        "only the block that actually differs should be written"
    )


def test_saving_a_broken_edit_keeps_the_text_but_does_not_activate(editable):
    """Never refuse the write; gate activation instead."""
    from xenium_viewer.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    from xenium_viewer.tabs.tab_templates import blocks_to_text
    spec = builtin_spec("clustering.leiden")
    _editor(widget).setPlainText(
        blocks_to_text(spec).replace("sc.pp.scale(adata_leiden, max_value=10)",
                                     "sc.pp.scale(adata_leiden,"))

    _button(widget, "Save  Activate").click()

    assert (tmp_path / "clustering.leiden.tmpl").is_file(), "the edit is kept"
    resolved = resolve("clustering.leiden")
    assert resolved.rejected
    assert not resolved.is_customised, "the shipped template runs instead"


def test_reverting_removes_the_override(editable):
    from xenium_viewer.utils.step_templates import resolve

    widget, tmp_path = editable
    _select(widget)
    from xenium_viewer.tabs.tab_templates import blocks_to_text
    spec = builtin_spec("clustering.leiden")
    _editor(widget).setPlainText(
        blocks_to_text(spec).replace("max_value=10", "max_value=99"))
    _button(widget, "Save  Activate").click()
    assert resolve("clustering.leiden").is_customised

    _button(widget, "Revert to default").click()
    assert list(tmp_path.glob("*.tmpl")) == []
    assert not resolve("clustering.leiden").is_customised


def test_saving_an_unchanged_edit_removes_rather_than_writes_an_inert_file(editable):
    widget, tmp_path = editable
    _select(widget)
    _button(widget, "Save  Activate").click()
    assert list(tmp_path.glob("*.tmpl")) == [], (
        "a file that says nothing would show as 'customised' forever"
    )


def test_a_frozen_block_edit_is_reported(editable):
    widget, _ = editable
    _select(widget, "Genes")
    from qtpy.QtCore import Qt
    tree = _tree(widget)
    cnv = next(
        tree.topLevelItem(g).child(c)
        for g in range(tree.topLevelItemCount())
        for c in range(tree.topLevelItem(g).childCount())
        if tree.topLevelItem(g).child(c).data(0, Qt.UserRole) == "genes.cnv_infercnv"
    )
    tree.setCurrentItem(cnv)

    editor = _editor(widget)
    editor.setPlainText(editor.toPlainText().replace(
        "_old_infer = pd.options.future.infer_string", "_old_infer = False"))
    _button(widget, "Validate").click()
    assert any("workaround" in p for p in _problems(widget))
