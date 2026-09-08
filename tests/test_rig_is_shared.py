"""One definition of "drive the running app", used by both callers.

``Rig`` was written inside ``scripts/capture_screenshots.py``. The end-to-end
suite needs the same verbs — navigate to a tab, find a widget by its label,
click it, wait for the worker it started — and the obvious move was to copy them
into ``tests/``. Two copies would disagree the first time either learned
something about a widget the other did not, and the one in ``scripts/`` is
exercised only when somebody re-shoots the documentation.

The same argument applies to the startup seam. The capture script used to call
``app._do_full_init`` with a hand-rolled ``_app`` dict that was missing
``plots_dock`` and ``plots_panel`` — a divergence that existed precisely because
nothing checked.

Source guards in the idiom of ``tests/test_plot_consistency.py``: parsed, not
imported, so this costs nothing and needs no Qt.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
CAPTURE = REPO / "scripts" / "capture_screenshots.py"
RIG = REPO / "src" / "palms" / "testing" / "rig.py"


@pytest.fixture(scope="module")
def capture_tree() -> ast.Module:
    return ast.parse(CAPTURE.read_text())


def _top_level_names(tree: ast.Module) -> set[str]:
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
    return names


def test_the_capture_script_does_not_define_its_own_rig(capture_tree):
    defined = _top_level_names(capture_tree)
    for name in ("Rig", "_process_events", "process_events", "BASE_LAYERS", "FRAME_LAYER"):
        assert name not in defined, (
            f"{name} is defined in scripts/capture_screenshots.py again — "
            f"it belongs in src/palms/testing/rig.py, which tests/e2e/ also uses")


def test_the_capture_script_imports_the_shared_rig(capture_tree):
    imported = {
        alias.asname or alias.name
        for node in ast.walk(capture_tree)
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("palms.testing")
        for alias in node.names
    }
    assert "Rig" in imported
    assert {"BASE_LAYERS", "FRAME_LAYER"} <= imported


def test_the_capture_script_uses_the_startup_seam(capture_tree):
    """It builds the app the way the app builds itself.

    ``create_viewer`` is ``run_viewer`` minus ``napari.run()``, so calling it is
    the only way to get the real wiring — menus, the close guard, the layer
    snapshot hook, a complete ``_app`` — without an event loop.

    Matched on calls rather than on the file's text: the comment explaining why
    ``_do_full_init`` is no longer used names it, and a grep cannot tell that
    apart from using it.
    """
    called = {
        node.func.id if isinstance(node.func, ast.Name) else node.func.attr
        for node in ast.walk(capture_tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, (ast.Name, ast.Attribute))
    }
    assert "create_viewer" in called
    assert "_do_full_init" not in called, (
        "call app.create_viewer(); _do_full_init skips the wiring around it")


def test_the_rig_does_not_import_napari():
    """It is handed a viewer; importing napari would make it a heavyweight import.

    ``palms.testing`` ships in the wheel, and the capture script resolves its
    dataset argument *before* paying for the napari import so that a typo fails
    in a second rather than after ten.
    """
    tree = ast.parse(RIG.read_text())
    modules = {
        node.module.split(".")[0]
        for node in ast.walk(tree) if isinstance(node, ast.ImportFrom) and node.module
    } | {
        alias.name.split(".")[0]
        for node in ast.walk(tree) if isinstance(node, ast.Import)
        for alias in node.names
    }
    assert "napari" not in modules
