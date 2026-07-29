"""Every recorded node must record code — or declare that it is not code.

Phase 0.3's invariant, enforced statically instead of one tab at a time. A node
whose cell parses to no executable statement replays as a silent no-op:
``allow_errors=False`` sees a cell that ran fine, the notebook "passes", and the
step it claims to document is simply missing. The Tier-2 verification finds
these on a *recorded session* — only for the actions that session happened to
take. This finds them in the source, including the tabs nobody exercised.

The escape hatch is ``kind=NOTE``: viewer state with no notebook equivalent (the
canvas background, an overlay's opacity) says so, renders as markdown, and is
counted separately. Anything else must emit a statement.

The check reads the ``ctx.record_node`` call sites with ``ast``: f-string
placeholders are replaced by a name, so ``f"x = {value}"`` is judged as ``x =
_v``. Only calls whose approximation *parses* are judged — an unparseable
approximation means the substitution was too crude to conclude anything.

Pure ast; no Qt, no imports of the tabs.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

TABS = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer" / "tabs"

# Nodes still recording prose, each with the reason it has not been migrated.
# Shrinking this list is the remaining Phase 0.3 work; growing it needs a reason.
KNOWN_PROSE = {
    # Computes a 2-D transcript histogram from the parquet — real analysis, and
    # a real gap. Needs the transcript loader expressed as plain pyarrow first.
    "viewer:transcript_density",
}


def _approximate(node: ast.AST) -> str | None:
    """Best-effort source text for a code argument, or None if not inferable."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                parts.append(value.value)
            elif isinstance(value, ast.FormattedValue):
                parts.append("_v")
            else:
                return None
        return "".join(parts)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left, right = _approximate(node.left), _approximate(node.right)
        return None if left is None or right is None else left + right
    if isinstance(node, ast.IfExp):
        return _approximate(node.body)          # the branch that adds text
    return None


def _record_node_calls():
    """(module, lineno, node_id, code_text, kind) for every record_node call."""
    for path in sorted(TABS.glob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for call in ast.walk(tree):
            if not isinstance(call, ast.Call):
                continue
            func = call.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
            if name != "record_node" or len(call.args) < 2:
                continue
            node_id = _approximate(call.args[0])
            code = _approximate(call.args[1])
            kind = next(
                (kw.value.id for kw in call.keywords
                 if kw.arg == "kind" and isinstance(kw.value, ast.Name)),
                "artifact",
            )
            yield path.name, call.lineno, node_id, code, kind


def _comment_only(code: str) -> bool:
    try:
        return not ast.parse(code).body
    except SyntaxError:
        return False        # the approximation was too crude to judge


def test_the_scan_finds_the_recorders():
    """Guard the guard: a silent scan that finds nothing would pass forever."""
    calls = list(_record_node_calls())
    assert len(calls) > 20
    assert any(code is not None for _m, _l, _n, code, _k in calls)


def test_no_recorder_emits_prose_where_code_is_expected():
    offenders = [
        f"{module}:{lineno} {node_id or '<computed id>'}"
        for module, lineno, node_id, code, kind in _record_node_calls()
        if kind != "NOTE" and code is not None and _comment_only(code)
        and (node_id or "") not in KNOWN_PROSE
    ]
    assert offenders == [], (
        "these record a comment where the notebook needs a statement — they "
        "replay as silent no-ops. Emit code, or record them with kind=NOTE if "
        "they are viewer state with no notebook equivalent: "
        + ", ".join(offenders)
    )


def test_declared_notes_really_are_comment_only():
    """The other direction: NOTE is for viewer state, not a way to hide code.

    A node marked NOTE whose body contains statements would be dropped from the
    notebook's code entirely — it renders as markdown.
    """
    offenders = [
        f"{module}:{lineno} {node_id}"
        for module, lineno, node_id, code, kind in _record_node_calls()
        if kind == "NOTE" and code is not None and not _comment_only(code)
    ]
    assert offenders == []


def test_the_known_prose_list_is_still_accurate():
    """A migrated node left in KNOWN_PROSE would mask its own regression."""
    recorded = {
        node_id for _m, _l, node_id, code, kind in _record_node_calls()
        if kind != "NOTE" and code is not None and _comment_only(code)
    }
    assert KNOWN_PROSE <= recorded, (
        "no longer prose — remove from KNOWN_PROSE: " + ", ".join(KNOWN_PROSE - recorded)
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
