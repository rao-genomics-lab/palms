"""The H&E registration steps record the code that produced the alignment.

Both fits used to record the 3x3 matrix they produced as a literal, under a
comment explaining where it came from. That passes
``test_recorded_code_is_code.py`` -- the cell does assign a variable -- and it is
still the failure this whole system exists to prevent: the recorded code was not
the code that ran, and a reader could not re-derive the number.

These are source guards in the idiom of ``test_persistence_safety.py``: they
parse the tab rather than driving a GUI, because what they assert is structural.
The behaviour of the fits themselves is measured in ``test_nuclei_align.py`` and
``test_coarse_align.py``.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from palms.utils.step_templates import builtin_ids, builtin_spec

TAB = Path(__file__).resolve().parents[1] / "src/palms/tabs/tab_he_registration.py"

#: The steps this tab runs, and the node kind each must carry. Every one of them
#: binds something a later step can consume -- an image, a flip, a transform --
#: so none of them is a TERMINAL, which ``prov_graph`` defines as "code, and no
#: dependents". They were all terminals until 2026-09-03, which is why nothing
#: declared a dependency on any of them.
_STEPS = {
    "he:load": "he.load",
    "he:flip": "he.flip",
    "he:coarse_align": "he.coarse_align",
    "he:nuclei_register": "he.nuclei_align",
    "he:landmark_register": "he.landmark_align",
}


@pytest.fixture(scope="module")
def tab_ast() -> ast.Module:
    return ast.parse(TAB.read_text())


def _step_calls(tree: ast.Module) -> dict[str, ast.Call]:
    """Every ``Step(id="...", ...)`` in the module, keyed by its node id."""
    out = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name) and node.func.id == "Step"):
            continue
        for kw in node.keywords:
            if kw.arg == "id" and isinstance(kw.value, ast.Constant):
                out[kw.value.value] = node
    return out


def _enclosing_function(tree: ast.Module, target: ast.AST) -> ast.AST:
    """The innermost ``def`` containing *target*."""
    best = None
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if any(sub is target for sub in ast.walk(node)):
            if best is None or node.lineno > best.lineno:
                best = node
    assert best is not None
    return best


def _keyword(call: ast.Call, name: str):
    for kw in call.keywords:
        if kw.arg == name:
            return kw.value
    return None


def test_every_he_step_is_a_step_not_a_record_node(tab_ast):
    """The five nodes are run through ``ctx.run_step``, which is what makes the
    recorded source identical to the executed source by construction."""
    found = _step_calls(tab_ast)
    missing = sorted(set(_STEPS) - set(found))
    assert not missing, (
        f"{missing} are not built as Step(...) in {TAB.name}. A hand-written "
        f"record_node beside the code that runs is two expressions of one "
        f"computation, which is what run_step exists to make impossible."
    )


def test_no_he_node_records_a_matrix_literal(tab_ast):
    """No recorded cell may be `he_affine = np.array([[...]])`.

    The specific regression: recording the *answer* rather than the computation.
    Such a cell replays, succeeds and reproduces nothing, and neither
    ``allow_errors=False`` nor the comment-only guard can see it.
    """
    for node in ast.walk(tab_ast):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "record_node"):
            continue
        source = ast.get_source_segment(TAB.read_text(), node) or ""
        assert "affine_3x3_yx).tolist()" not in source, (
            "a record_node call inlines a fitted transform as a literal; record "
            "the code that produced it (see he.coarse_align.tmpl)"
        )


@pytest.mark.parametrize("node_id", sorted(_STEPS))
def test_every_he_step_is_an_artifact(tab_ast, node_id):
    """Reusable state, not a side effect.

    Qt propagates nothing here -- this is about the graph. ``TERMINAL`` is
    documented as "code, and no dependents", so a terminal is a node nothing may
    declare a dependency on, and the fine fit genuinely does depend on the seed
    the coarse search produced.
    """
    kind = _keyword(_step_calls(tab_ast)[node_id], "kind")
    assert isinstance(kind, ast.Name) and kind.id == "ARTIFACT", (
        f"{node_id} is not recorded as an ARTIFACT"
    )


def test_the_fits_depend_on_the_image_and_the_flip(tab_ast):
    """Staleness has to propagate, which is the point of getting the edges right.

    Re-loading a different H&E, or re-ticking a flip, invalidates both fits --
    the GUI has always *acted* on that (``on_flip_changed`` drops the coarse
    affine) while the graph said every node depended on ``preamble`` alone.
    """
    source = TAB.read_text()
    calls = _step_calls(tab_ast)
    for node_id in ("he:coarse_align", "he:nuclei_register"):
        deps = _keyword(calls[node_id], "deps")
        # The nuclei fit assembles its list at runtime, because the seed edge
        # varies -- so take the names from wherever they are written rather than
        # insisting the literal sit at the call site. The enclosing builder is
        # small and does nothing else, which is what makes reading it safe.
        segment = ast.get_source_segment(source, deps) or ""
        if not segment.strip().startswith("["):
            builder = _enclosing_function(tab_ast, deps)
            segment = ast.get_source_segment(source, builder) or ""
        for needed in ("he:load", "he:flip"):
            assert f'"{needed}"' in segment, (
                f"{node_id} does not declare a dependency on {needed}; it reads "
                f"the image he:load binds and works in the frame he:flip declares"
            )


def test_the_nuclei_fit_is_never_seeded_by_itself(tab_ast):
    """A step that depends on its own output is a cycle, and unreplayable.

    The GUI used to seed a re-run from ``he_state['affine_3x3']`` whatever
    produced it, so pressing Fine Align twice fed the fit its own answer. No
    notebook can express that, so the seed is chosen from a *different* step's
    output or inlined as a restored literal.
    """
    # By AST rather than by slicing the file between two function names: the
    # first draft of this test did the latter and broke the moment the function
    # after `_nuclei_seed` was renamed, which is a guard failing for a reason
    # that has nothing to do with what it guards.
    body = ""
    for node in ast.walk(tab_ast):
        if (isinstance(node, ast.FunctionDef) and node.name == "_nuclei_seed"):
            body = ast.unparse(node)
    assert body, "_nuclei_seed not found"
    assert "'he:nuclei_register'" not in body, (
        "the nuclei fit names itself as a possible seed dependency"
    )
    assert "== 'landmarks'" in body, (
        "the seed_fine branch must check which method produced the fine affine; "
        "without it a previous nuclei fit is indistinguishable from a landmark one"
    )


@pytest.mark.parametrize("template_id", sorted(_STEPS.values()))
def test_every_he_template_is_registered(template_id):
    assert template_id in builtin_ids()
    spec = builtin_spec(template_id)
    assert spec.outputs, f"{template_id} declares no outputs"


def test_the_flip_template_needs_no_palms_import():
    """Two booleans are two booleans.

    The allowance is per-template and granting it where it is not needed would
    make the exemption list stop meaning anything -- so the one H&E template that
    computes nothing does not carry it.
    """
    spec = builtin_spec("he.flip")
    assert not spec.palms_reason
    assert all("palms" not in b.text for b in spec.blocks.values())
