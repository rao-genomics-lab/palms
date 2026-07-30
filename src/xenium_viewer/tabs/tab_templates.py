"""Tab: Templates — what each analysis button is actually going to run.

The source a step executes has always been recoverable *after* the fact, from
the Notebook tab. What there was no way to see is what a step *would* run before
running it — the text lived in private module constants, and the parameters only
became visible once a ten-minute analysis had already produced a cell.

This view is read-only on purpose. It shows, per registered template: the
contract (params, required namespace names, declared outputs), the shipped
source of every block, and a **live preview** of exactly the string that would
be handed to ``exec`` and recorded, rendered with the parameters currently set
in the owning tab.

The preview is not a reconstruction. It calls ``Step.render()`` — the same
method ``StepExecutor.run`` calls — so what is displayed here is what would run,
by the same construction that makes the recorded code equal the executed code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from qtpy.QtWidgets import (
    QHBoxLayout, QLabel, QPlainTextEdit, QSplitter, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

from xenium_viewer.tabs.tab_notebook import PythonSyntaxHighlighter
from xenium_viewer.utils.step_templates import builtin_ids, builtin_spec
from xenium_viewer.utils.steps import Step, StepError

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext

#: Tree grouping. A template's id prefix says which part of the app owns it;
#: the two setup steps have no prefix because every analysis depends on them.
_GROUPS = {
    "clustering": "Clustering",
    "genes": "Genes",
    "spatial": "Spatial",
    "roi": "ROI",
}
_SETUP_GROUP = "Setup"


def _group_of(template_id: str) -> str:
    return _GROUPS.get(template_id.split(".", 1)[0], _SETUP_GROUP)


def _mono(widget) -> None:
    font = QFont("monospace")
    font.setStyleHint(QFont.TypeWriter)
    widget.setFont(font)
    widget.setTabStopDistance(28)


def _contract_text(spec) -> str:
    """The template's contract, as a short readable block."""
    lines = []
    if spec.doc:
        lines += [spec.doc, ""]
    required = sorted(spec.required_params)
    optional = sorted(spec.param_names - spec.required_params)
    lines.append(f"params    : {', '.join(required) or '(none)'}")
    if optional:
        lines.append(f"optional  : {', '.join(optional)}")
    lines.append(f"requires  : {', '.join(sorted(spec.requires)) or '(none)'}")
    lines.append(f"outputs   : {', '.join(spec.outputs) or '(none)'}")
    lines.append(f"blocks    : {', '.join(spec.blocks)}")
    lines.append(f"assemblies: {len(spec.assemblies)}")
    frozen = [b.name for b in spec.blocks.values() if not b.editable]
    if frozen:
        lines.append(f"frozen    : {', '.join(frozen)}  (version workaround)")
    return "\n".join(lines)


def _preview(ctx: ViewerContext, spec, assembly) -> str:
    """The exact source this template would execute, for *assembly*.

    Uses the owning tab's current widget values when it has registered a
    ``preview_params`` callable, and falls back to synthesised literals of the
    right type otherwise — so the preview is real for a live tab and still
    illustrative for one the user has not opened.
    """
    providers = ctx.state.get("template_preview_params", {})
    provider = providers.get(spec.id)
    params, source = None, "sample values"
    if callable(provider):
        try:
            params = provider()
            source = "current widget values"
        except Exception:                      # a half-built tab must not break the view
            params = None
    if params is None:
        params = spec.synth_params()
    else:
        # A provider may legitimately omit optional params; fill the rest so the
        # template still renders rather than raising at the user.
        merged = spec.synth_params()
        merged.update(params)
        params = {k: v for k, v in merged.items() if f"${k}" in spec.assemble(assembly)}

    step = Step(id=f"preview:{spec.id}", template=spec.assemble(assembly),
                params=params, template_id=spec.id)
    try:
        return f"# preview — {source}\n{step.render()}"
    except StepError as exc:
        return f"# could not render this template:\n# {exc}"


def build_tab(ctx: ViewerContext) -> tuple:
    tree = QTreeWidget()
    tree.setHeaderLabels(["Template", "Blocks"])
    tree.setColumnWidth(0, 260)

    contract = QPlainTextEdit()
    contract.setReadOnly(True)
    contract.setMaximumHeight(150)
    _mono(contract)

    source = QPlainTextEdit()
    source.setReadOnly(True)
    _mono(source)
    source_hl = PythonSyntaxHighlighter(source.document())

    preview = QPlainTextEdit()
    preview.setReadOnly(True)
    _mono(preview)
    preview_hl = PythonSyntaxHighlighter(preview.document())
    # Held so Qt does not garbage-collect the highlighters with the local names.
    ctx.state["_template_highlighters"] = (source_hl, preview_hl)

    def _populate():
        tree.clear()
        groups: dict[str, QTreeWidgetItem] = {}
        for template_id in builtin_ids():
            spec = builtin_spec(template_id)
            name = _group_of(template_id)
            parent = groups.get(name)
            if parent is None:
                parent = QTreeWidgetItem(tree, [name, ""])
                parent.setExpanded(True)
                groups[name] = parent
            item = QTreeWidgetItem(parent, [template_id, str(len(spec.blocks))])
            item.setData(0, Qt.UserRole, template_id)

    def _on_selection():
        items = tree.selectedItems()
        if not items:
            return
        template_id = items[0].data(0, Qt.UserRole)
        if not template_id:                    # a group header
            contract.setPlainText("")
            source.setPlainText("")
            preview.setPlainText("")
            return
        spec = builtin_spec(template_id)
        contract.setPlainText(_contract_text(spec))

        blocks = []
        for block in spec.blocks.values():
            mark = "" if block.editable else "   [frozen]"
            blocks.append(f"# ── block: {block.name}{mark}"
                          f"{block.text.rstrip()}\n")
        source.setPlainText("\n".join(blocks).strip("\n"))

        # The first declared assembly is the one with no optional blocks, which
        # is the most representative thing to show without a widget to choose.
        preview.setPlainText(_preview(ctx, spec, spec.assemblies[0]))

    tree.itemSelectionChanged.connect(_on_selection)
    _populate()

    right = QWidget()
    right_layout = QVBoxLayout()
    right_layout.setContentsMargins(0, 0, 0, 0)
    right_layout.addWidget(QLabel("Contract"))
    right_layout.addWidget(contract)

    panes = QSplitter(Qt.Vertical)
    for label, editor in (("Shipped source (read-only)", source),
                          ("Preview — what would run", preview)):
        holder = QWidget()
        holder_layout = QVBoxLayout()
        holder_layout.setContentsMargins(0, 0, 0, 0)
        holder_layout.addWidget(QLabel(label))
        holder_layout.addWidget(editor)
        holder.setLayout(holder_layout)
        panes.addWidget(holder)
    right_layout.addWidget(panes, 1)
    right.setLayout(right_layout)

    split = QSplitter(Qt.Horizontal)
    split.addWidget(tree)
    split.addWidget(right)
    split.setStretchFactor(1, 3)

    widget = QWidget()
    layout = QHBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)
    layout.addWidget(split)
    widget.setLayout(layout)
    return widget, {}
