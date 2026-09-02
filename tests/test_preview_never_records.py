"""`preview` executes analysis code that is never recorded. Keep it that way.

`StepExecutor.run` is the one path from a `Step` to executed *and* recorded
code, and that is the property the whole provenance claim rests on:
the recorded source is the executed source, by construction. `preview` exists
so display code can reuse a template's text instead of keeping a second
implementation of the same computation — the drift risk that made CopyKAT's
reconstruction cell a documented exception rather than a pattern.

That makes it the one function in the codebase that can run analysis without
leaving a node behind, so it is fenced on four sides:

* recording is not a *flag* on `run` — there is no way to ask `run` not to record;
* `preview` cannot reach the graph, and executes into a copy of the namespace;
* only declared call sites may call it, and they are listed here with reasons;
* those call sites must label what they draw, because an unlabelled preview
  is precisely the hazard — a picture that looks like a result and is not.

In the idiom of `tests/test_persistence_safety.py`'s
`test_no_source_file_still_does_delete_then_write`.
"""
from __future__ import annotations

import ast
import inspect
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from palms.utils.prov_graph import ProvGraph  # noqa: E402
from palms.utils.steps import Step, StepError, StepExecutor  # noqa: E402

SRC = Path(__file__).resolve().parent.parent / "src" / "palms"

#: Every place in the app that may execute a Step without recording it, with
#: the reason it is allowed to. Compared by equality in both directions: a
#: stale entry rots the list the way `KNOWN_PROSE` was allowed to.
ALLOWED_PREVIEW_CALL_SITES = {
    ("tabs/_helpers.py", "create_shared_helpers._preview_step"):
        "the ctx binding itself — the single door",
    ("tabs/tab_transcripts.py", "build_tab._preview_start._work"):
        "the density preview: display only, labelled, never a node",
}

#: Anything that would make a preview change the analysis rather than draw it.
_RECORDING_NAMES = {
    "run_step", "record_node", "record_code", "record_clustering",
    "record_preamble", "record_environment", "ensure_normalized",
    "ensure_annotations", "ensure_spatial_neighbors", "upsert",
}


def _module_ast(relative: str) -> ast.Module:
    path = SRC / relative
    return ast.parse(path.read_text(), str(path))


def _enclosing_chain(tree: ast.Module) -> dict[ast.AST, str]:
    """Map every node to the dotted chain of functions enclosing it."""
    chain: dict[ast.AST, str] = {}

    def walk(node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                name = f"{prefix}.{child.name}" if prefix else child.name
                chain[child] = name
                walk(child, name)
            else:
                chain[child] = prefix
                walk(child, prefix)

    walk(tree, "")
    return chain


def _preview_calls() -> dict[tuple[str, str], int]:
    """Every call to `preview` / `preview_step` under src/palms, by call site."""
    found: dict[tuple[str, str], int] = {}
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(), str(path))
        chain = _enclosing_chain(tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = (func.attr if isinstance(func, ast.Attribute)
                    else func.id if isinstance(func, ast.Name) else None)
            if name in ("preview", "preview_step"):
                key = (str(path.relative_to(SRC)), chain.get(node, ""))
                found[key] = found.get(key, 0) + 1
    return found


def _function_def(relative: str, qualname: str) -> ast.FunctionDef:
    tree = _module_ast(relative)
    chain = _enclosing_chain(tree)
    for node, name in chain.items():
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and name == qualname:
            return node
    raise AssertionError(f"{qualname} not found in {relative}")


# ── behaviour ────────────────────────────────────────────────────────────────

def _executor(**ns):
    return StepExecutor(namespace=dict(ns), graph=ProvGraph())


def test_a_preview_records_nothing():
    """Not "no node with that id" — no nodes at all."""
    ex = _executor(a=1)
    step = Step(id="t", template="b = a + extra\n", outputs=["b"])

    assert ex.preview(step, bindings={"extra": 41}) == {"b": 42}
    assert list(ex.graph.nodes()) == []
    assert "t" not in ex.graph


def test_a_preview_cannot_touch_the_shared_namespace():
    """The copy is the property, and it is why a preview cannot make a later
    recorded step consume values it did not fetch."""
    sentinel = object()
    ex = _executor(adata=sentinel)
    ex.preview(Step(id="t", template="adata = 'clobbered'\nout = 1\n", outputs=["out"]))

    assert ex.ns["adata"] is sentinel
    assert "out" not in ex.names()


def test_a_preview_may_substitute_a_result_but_not_the_data_behind_it():
    """Substituting `transcript_points` is the whole point, and it is normally
    already bound by the last recorded run — an earlier version refused that
    and so stopped previewing after the first Compute Density. Swapping `adata`
    is the real hazard: a picture drawn from data the recorded step never uses."""
    ex = _executor(adata=object(), transcript_points="from the points element")
    step = Step(id="t", template="out = transcript_points\n", outputs=["out"])

    with pytest.raises(StepError, match="shadow"):
        ex.preview(step, bindings={"adata": object()})

    assert ex.preview(step, bindings={"transcript_points": "from the index"}) == {
        "out": "from the index"}
    assert ex.ns["transcript_points"] == "from the points element"


def test_a_preview_that_does_not_bind_its_output_is_an_error():
    ex = _executor()
    with pytest.raises(StepError, match="did not bind"):
        ex.preview(Step(id="t", template="pass\n", outputs=["missing"]))


# ── the guards ───────────────────────────────────────────────────────────────

def test_recording_is_not_a_mode_of_run():
    """`run` records unconditionally: no flag, and no branch around the upsert.

    A `record=False` parameter would put a condition above the single line that
    guarantees "executed == recorded", and every later reader of `run` would
    have to check the call site before believing its docstring.
    """
    node = _function_def("utils/steps.py", "StepExecutor.run")
    args = node.args
    names = {a.arg for a in args.args + args.kwonlyargs}
    assert not names & {"record", "preview", "dry_run", "no_record"}, (
        "recording must not become a mode of run() — add a separate method")

    upserts = [n for n in ast.walk(node)
               if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
               and n.func.attr == "upsert"]
    assert len(upserts) == 1, "run() must record exactly once"

    guarded = [n for n in ast.walk(node) if isinstance(n, (ast.If, ast.Try))
               for inner in ast.walk(n) if inner is upserts[0]]
    assert guarded == [], "run()'s upsert must not sit inside a branch"


def test_preview_cannot_reach_the_graph():
    """Read the code, not the prose — the docstring says "graph" a lot."""
    node = _function_def("utils/steps.py", "StepExecutor.preview")
    reached = {n.attr for n in ast.walk(node) if isinstance(n, ast.Attribute)}
    assert not reached & {"graph", "upsert"}, (
        f"preview must have no route to the provenance graph, but reaches {reached}")
    assert "dict(self.ns)" in inspect.getsource(StepExecutor.preview), (
        "preview must execute into a copy of the namespace, not the namespace")


def test_only_declared_call_sites_run_a_preview():
    found = set(_preview_calls())
    allowed = set(ALLOWED_PREVIEW_CALL_SITES)

    assert found == allowed, (
        "a preview executes analysis code that is never recorded. A new call "
        "site must be display-only and added to ALLOWED_PREVIEW_CALL_SITES with "
        f"its reason.\n  unexpected: {sorted(found - allowed)}"
        f"\n  stale entries: {sorted(allowed - found)}")


def test_the_scan_finds_the_call_sites():
    """Guard the guard: a scan that silently matches nothing passes forever."""
    found = _preview_calls()
    assert len(found) >= 2 and sum(found.values()) >= 2


def test_a_previewed_step_id_cannot_collide_with_a_recorded_node():
    """Belt and braces: even handed to run_step by mistake, a preview Step
    could not overwrite `transcript_density:<gene>`."""
    node = _function_def("tabs/tab_transcripts.py", "build_tab._preview_bindings")
    ids = [kw.value for call in ast.walk(node)
           if isinstance(call, ast.Call) and getattr(call.func, "id", None) == "Step"
           for kw in call.keywords if kw.arg == "id"]
    assert ids, "the preview must build a Step with an explicit id"
    for value in ids:
        first = value.values[0] if isinstance(value, ast.JoinedStr) else value
        assert isinstance(first, ast.Constant) and first.value.startswith("preview:"), (
            "a preview Step's id must start with 'preview:'")


def test_the_preview_callbacks_do_not_record():
    for qualname in ("build_tab._preview_start", "build_tab._on_preview_ready",
                     "build_tab._preview_bindings"):
        node = _function_def("tabs/tab_transcripts.py", qualname)
        called = {n.func.attr if isinstance(n.func, ast.Attribute) else
                  getattr(n.func, "id", None)
                  for n in ast.walk(node) if isinstance(n, ast.Call)}
        offenders = sorted(called & _RECORDING_NAMES)
        assert offenders == [], (
            f"{qualname} draws a preview; it must not {offenders} — that would "
            "make drawing a picture change the analysis")


def test_the_preview_is_labelled():
    """The labelling is part of the guard: an unlabelled preview is the hazard."""
    from palms.tabs import tab_transcripts

    assert "PREVIEW" in tab_transcripts.PREVIEW_LABEL
    assert "PREVIEW" in tab_transcripts.PREVIEW_LAYER_NAME
    assert tab_transcripts.PREVIEW_LAYER_NAME in (SRC / "app.py").read_text(), (
        "the preview layer app.py creates must carry PREVIEW_LAYER_NAME")

    node = _function_def("tabs/tab_transcripts.py", "build_tab._on_preview_ready")
    names = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
    assert "PREVIEW_LABEL" in names, "the status line must say it is a preview"


def test_a_previewed_template_does_not_write_to_adata_or_sdata():
    """The namespace copy is shallow, so `adata` and `sdata` are the same
    objects: a template that mutated them would mutate them for real."""
    from palms.utils.step_templates import builtin_spec

    for template_id in ("transcripts.density",):
        blocks = builtin_spec(template_id).blocks
        text = "\n".join(b.text for b in blocks.values())
        # `$param` is not valid Python until it is substituted; the shape of the
        # assignment targets is what matters here, not the values.
        tree = ast.parse(re.sub(r"\$\{?(\w+)\}?", r"_param_\1", text))
        for node in ast.walk(tree):
            targets = (node.targets if isinstance(node, ast.Assign)
                       else [node.target] if isinstance(node, (ast.AugAssign, ast.AnnAssign))
                       else [])
            for target in targets:
                root = target
                while isinstance(root, (ast.Subscript, ast.Attribute)):
                    root = root.value
                assert getattr(root, "id", None) not in ("adata", "sdata"), (
                    f"{template_id} writes to {root.id}; it may not be previewed")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
