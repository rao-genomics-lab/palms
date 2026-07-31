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

import difflib
from dataclasses import replace
from typing import TYPE_CHECKING

from qtpy.QtCore import Qt
from qtpy.QtGui import QFont
from qtpy.QtWidgets import (
    QHBoxLayout, QLabel, QListWidget, QPlainTextEdit, QPushButton, QSplitter,
    QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from xenium_viewer.tabs.tab_notebook import PythonSyntaxHighlighter
from xenium_viewer.utils.step_templates import (
    ERROR, EXECUTOR_BASE_NAMES, builtin_ids, builtin_spec, parse_template,
    resolve, validate,
)
from xenium_viewer.utils.step_templates.overrides import (
    remove_override, save_override,
)
from xenium_viewer.utils.step_templates.spec import BLOCK_MARKER, BlockSpec
from xenium_viewer.utils.step_templates.validate import Problem
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


def blocks_to_text(spec) -> str:
    """A template's blocks as one editable ``.tmpl`` body.

    The user edits the whole thing rather than one block at a time: the block
    markers are already how the file is written, and a per-block editor would
    hide the boundaries that matter (which block a statement belongs to decides
    whether it runs at all).

    Nothing is annotated onto the marker line — not even "frozen". Anything
    after the block name is *part of the name*, so a helpful note there comes
    back as a different block on the round trip. Which blocks are frozen is
    stated in the Contract pane instead.
    """
    parts = [
        f"{BLOCK_MARKER}{block.name}\n{block.text.strip(chr(10))}"
        for block in spec.blocks.values()
    ]
    return "\n\n".join(parts) + "\n"


def text_to_blocks(text: str) -> dict:
    """Parse an edited body back into ``{block name: text}``.

    Mirrors the loader's own splitter so the editor and the file format cannot
    disagree about where a block starts.
    """
    spec = parse_template(f"# id: _editor\n# schema-version: 1\n\n{text}")
    return {name: block.text for name, block in spec.blocks.items()}


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
    """The exact source this template would execute right now.

    Uses the owning tab's :class:`Preview` when it has registered one — *both*
    the blocks it would select and the params it would pass, so the pane tracks
    the widgets in shape as well as in value. Falls back to *assembly* and
    synthesised literals of the right type otherwise, so the preview is real for
    a live tab and still illustrative for one the user has not opened.

    Blocks have to come from the provider rather than from a fixed assembly.
    The call site owns block selection by design (the branch structure *is* what
    the widgets mean), so a params-only provider left the preview pinned to the
    first declared assembly: untick "use HVGs" and the numbers moved while the
    code shape did not.
    """
    providers = ctx.state.get("template_preview", {})
    provider = providers.get(spec.id)
    params, note, source = None, "", "sample values"
    if callable(provider):
        try:
            blocks, params, note = provider()
            assembly = list(blocks)
            source = "current widget values"
        except Exception:                      # a half-built tab must not break the view
            params = None
    try:
        template = spec.assemble(assembly)
    except KeyError as exc:                    # a provider naming a block an override dropped
        return f"# could not assemble this template:\n# {exc}"

    if params is None:
        params = spec.synth_params()
    else:
        # A provider may legitimately omit optional params; fill the rest so the
        # template still renders rather than raising at the user.
        merged = spec.synth_params()
        merged.update(params)
        params = {k: v for k, v in merged.items() if f"${k}" in template}

    step = Step(id=f"preview:{spec.id}", template=template,
                params=params, template_id=spec.id)
    header = f"# preview — {source}" + (f" ({note})" if note else "")
    try:
        return f"{header}\n{step.render()}"
    except StepError as exc:
        return f"# could not render this template:\n# {exc}"


def _diff_text(resolved) -> str:
    """Yours vs the new shipped text, for the blocks that moved.

    Two-way, not three-way. A real merge would put conflict markers into Python
    source that the user then has to un-mangle by hand, in a language where a
    stray ``<<<<<<<`` is a syntax error rather than a visible annotation. Seeing
    both versions and choosing is the honest amount of help.
    """
    if not resolved.stale_blocks:
        return ("# The shipped template's contract changed (schema-version).\n"
                "# Re-check the parameters in the Contract pane above.")
    chunks = []
    for name in resolved.stale_blocks:
        # splitlines() without keepends, and lineterm="": block text has no
        # trailing newline, so keepends left the final '-' and '+' lines
        # concatenated on one line — in the pane whose whole job is showing the
        # user which line changed.
        diff = difflib.unified_diff(
            resolved.spec.blocks[name].text.splitlines(),
            resolved.builtin.blocks[name].text.splitlines(),
            fromfile=f"yours/{name}", tofile=f"new default/{name}",
            n=3, lineterm="",
        )
        chunks.append("\n".join(diff).rstrip())
    return "\n\n".join(chunks)


def check_edit(template_id: str, text: str) -> list:
    """Validate an edited body without saving it. Returns Problems.

    The whole point of a Validate button is that it answers in milliseconds
    against no data — a template edit that can only be checked by running a
    ten-minute analysis would not get checked.
    """
    builtin = builtin_spec(template_id)
    try:
        blocks = text_to_blocks(text)
    except Exception as exc:
        return [Problem(f"could not parse the edited template: {exc}")]

    # Merged onto the shipped blocks, not substituted for them: an override is
    # per block, so a body that omits a block inherits it rather than deleting
    # it. Replacing the whole dict reported a cascade of "no such block" errors
    # for every block the user simply had not changed.
    merged_blocks = dict(builtin.blocks)
    for name, text in blocks.items():
        editable = builtin.blocks[name].editable if name in builtin.blocks else True
        merged_blocks[name] = BlockSpec(name=name, text=text, editable=editable)

    return validate(replace(builtin, blocks=merged_blocks), builtin=builtin,
                    available=EXECUTOR_BASE_NAMES | builtin.requires)


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

    editor = QPlainTextEdit()
    _mono(editor)

    preview = QPlainTextEdit()
    preview.setReadOnly(True)
    _mono(preview)
    # Held so Qt does not garbage-collect the highlighters with the local names.
    ctx.state["_template_highlighters"] = tuple(
        PythonSyntaxHighlighter(w.document()) for w in (source, editor, preview)
    )

    problems_list = QListWidget()
    problems_list.setMaximumHeight(110)

    validate_button = QPushButton("Validate")
    save_button = QPushButton("Save && Activate")
    revert_button = QPushButton("Revert to default")
    # Only meaningful while a review is pending; shown by _on_selection.
    take_button = QPushButton("Take new default for changed blocks")
    take_button.setVisible(False)
    status = QLabel("")
    status.setWordWrap(True)

    current = {"id": None}

    def _badge(template_id: str) -> str:
        resolved = resolve(template_id)
        if resolved.rejected:
            return "✕ not used"
        if resolved.needs_review:
            return "⚠ review"
        if resolved.is_customised:
            return f"● customised ({len(resolved.changed_blocks())})"
        return ""

    def _populate(keep: str = None):
        tree.blockSignals(True)
        tree.clear()
        groups: dict[str, QTreeWidgetItem] = {}
        chosen = None
        for template_id in builtin_ids():
            spec = builtin_spec(template_id)
            name = _group_of(template_id)
            parent = groups.get(name)
            if parent is None:
                parent = QTreeWidgetItem(tree, [name, ""])
                parent.setExpanded(True)
                groups[name] = parent
            item = QTreeWidgetItem(
                parent, [template_id, _badge(template_id) or str(len(spec.blocks))])
            item.setData(0, Qt.UserRole, template_id)
            if template_id == keep:
                chosen = item
        tree.blockSignals(False)
        if chosen is not None:
            tree.setCurrentItem(chosen)

    def _show_problems(problems) -> int:
        problems_list.clear()
        errors = 0
        for problem in problems:
            problems_list.addItem(str(problem))
            errors += problem.severity == ERROR
        return errors

    def _on_selection():
        items = tree.selectedItems()
        if not items:
            return
        template_id = items[0].data(0, Qt.UserRole)
        current["id"] = template_id
        for widget in (validate_button, save_button, revert_button):
            widget.setEnabled(bool(template_id))
        if not template_id:                    # a group header
            for pane in (contract, source, editor, preview):
                pane.setPlainText("")
            problems_list.clear()
            status.setText("")
            return

        builtin = builtin_spec(template_id)
        resolved = resolve(template_id)
        contract.setPlainText(_contract_text(builtin))
        source.setPlainText(blocks_to_text(builtin))
        editor.setPlainText(blocks_to_text(resolved.spec))
        _show_problems(resolved.problems)
        take_button.setVisible(resolved.needs_review)
        if resolved.rejected:
            status.setText(
                "⚠ This customised template was NOT used — the shipped version "
                "ran instead. Fix the problems below and save again."
            )
        elif resolved.needs_review:
            what = []
            if resolved.stale_blocks:
                what.append(f"block(s) {', '.join(resolved.stale_blocks)} changed")
            if resolved.schema_moved:
                what.append("its contract changed")
            status.setText(
                f"⚠ Your customisation is still active, but the shipped template "
                f"has moved on since you forked it ({'; '.join(what)}). The "
                f"update may be a fix your edit is now shadowing — see the diff "
                f"below, then either take the new default or save again to "
                f"confirm you have reviewed it."
            )
        elif resolved.is_customised:
            status.setText(
                f"● Customised: {', '.join(resolved.changed_blocks())}. "
                f"Every other block still tracks the shipped template."
            )
        else:
            status.setText("Shipped template, unmodified.")

        if resolved.needs_review:
            preview_label.setText("Diff — yours vs the new default")
            preview.setPlainText(_diff_text(resolved))
        else:
            preview_label.setText("Preview — what would run")
            preview.setPlainText(
                _preview(ctx, resolved.spec, resolved.spec.assemblies[0]))

    def _on_validate():
        template_id = current["id"]
        if not template_id:
            return
        problems = check_edit(template_id, editor.toPlainText())
        errors = _show_problems(problems)
        status.setText(
            "No problems — safe to save and activate." if not problems
            else f"{errors} error(s), {len(problems) - errors} warning(s)."
        )

    def _on_save():
        template_id = current["id"]
        if not template_id:
            return
        text = editor.toPlainText()
        problems = check_edit(template_id, text)
        errors = _show_problems(problems)
        try:
            blocks = text_to_blocks(text)
        except Exception as exc:
            status.setText(f"Could not parse the edit, so nothing was saved: {exc}")
            return

        # Saved either way. Refusing to write would send the user to an external
        # editor and out of the loop that gives them feedback; what is gated is
        # *activation*, and an invalid file on disk is simply rejected by the
        # resolver, which falls back to the shipped template and says so.
        path = save_override(template_id, blocks)
        _populate(keep=template_id)
        if errors:
            status.setText(
                f"Saved to {path}, but NOT active: {errors} error(s) below. "
                f"The shipped template will run until they are fixed."
            )
        elif path is None:
            status.setText("Matches the shipped template — customisation removed.")
        else:
            status.setText(f"Saved and active: {path}")

    def _on_revert():
        template_id = current["id"]
        if not template_id:
            return
        removed = remove_override(template_id)
        _populate(keep=template_id)
        status.setText("Reverted to the shipped template." if removed
                       else "Already the shipped template.")

    def _on_take_new():
        """Drop the customisation of the blocks that moved, keeping the rest.

        Deliberately not "revert everything": a user with three customised
        blocks who accepts the new version of one should not silently lose the
        other two.
        """
        template_id = current["id"]
        if not template_id:
            return
        resolved = resolve(template_id)
        if not resolved.stale_blocks:
            return
        blocks = {name: block.text for name, block in resolved.spec.blocks.items()}
        for name in resolved.stale_blocks:
            blocks[name] = resolved.builtin.blocks[name].text
        save_override(template_id, blocks)
        _populate(keep=template_id)
        status.setText(
            f"Took the new default for {', '.join(resolved.stale_blocks)}; "
            f"your other customisations are unchanged."
        )

    tree.itemSelectionChanged.connect(_on_selection)
    validate_button.clicked.connect(_on_validate)
    save_button.clicked.connect(_on_save)
    revert_button.clicked.connect(_on_revert)
    take_button.clicked.connect(_on_take_new)
    _populate()

    right = QWidget()
    right_layout = QVBoxLayout()
    right_layout.setContentsMargins(0, 0, 0, 0)
    right_layout.addWidget(QLabel("Contract"))
    right_layout.addWidget(contract)

    # Default beside Yours, so a customisation is read as a diff rather than as
    # a file with no reference point.
    side_by_side = QSplitter(Qt.Horizontal)
    for label, pane in (("Default (read-only)", source), ("Yours (editable)", editor)):
        holder = QWidget()
        holder_layout = QVBoxLayout()
        holder_layout.setContentsMargins(0, 0, 0, 0)
        holder_layout.addWidget(QLabel(label))
        holder_layout.addWidget(pane)
        holder.setLayout(holder_layout)
        side_by_side.addWidget(holder)

    buttons = QHBoxLayout()
    for button in (validate_button, save_button, revert_button, take_button):
        buttons.addWidget(button)
    buttons.addStretch()
    button_bar = QWidget()
    button_bar.setLayout(buttons)

    lower = QWidget()
    lower_layout = QVBoxLayout()
    lower_layout.setContentsMargins(0, 0, 0, 0)
    preview_label = QLabel("Preview — what would run")
    lower_layout.addWidget(preview_label)
    lower_layout.addWidget(preview)
    lower_layout.addWidget(QLabel("Problems"))
    lower_layout.addWidget(problems_list)
    lower.setLayout(lower_layout)

    panes = QSplitter(Qt.Vertical)
    panes.addWidget(side_by_side)
    panes.addWidget(lower)
    panes.setStretchFactor(0, 2)

    right_layout.addWidget(button_bar)
    right_layout.addWidget(status)
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
