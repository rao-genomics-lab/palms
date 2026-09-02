"""Naming clusters must reach the Coloring tab without reselecting the clustering.

Reported from the GUI 2026-09-02, long-standing: after naming clusters with the
LLM, the new names did not appear until the user selected a different clustering
and came back, or opened and closed the label editor.

The annotation paths write ``state['cluster_labels']`` and persist it, so the
data was correct everywhere and only the *widget* was stale — the checkbox
captions are built once, in ``ctx.repopulate_cluster_checkboxes``, and
``tab_gene_analysis`` never called it. The three call sites that did exist
(``app.py`` and two in ``tab_cell_coloring``) are exactly the set of workarounds
users found, which is why this read as a cosmetic lag for so long.

The source guard below is the point of this module. The bug was an *omission*,
and the first attempt at the fix missed one of the three writers — label
transfer — because the report only named the LLM. A test that pins the three
known handlers would have passed on that mistake; one that derives the set from
the source does not.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

TAB = (Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs"
       / "tab_gene_analysis.py")
REFRESH = "refresh_cluster_labels_in_coloring_tab"


def _fake_ctx(shown: str | None):
    calls = []
    return SimpleNamespace(
        repopulate_cluster_checkboxes=lambda: calls.append(1),
        clustering_widget=SimpleNamespace(value=shown),
    ), calls


def test_the_shown_clustering_is_refreshed():
    from palms.tabs.tab_gene_analysis import refresh_cluster_labels_in_coloring_tab

    ctx, calls = _fake_ctx("leiden_r1.0")
    refresh_cluster_labels_in_coloring_tab(ctx, "leiden_r1.0")
    assert calls == [1]


def test_another_clustering_is_left_alone():
    """Rebuilding re-checks every box, so refreshing a clustering the user is not
    looking at would silently clear their cluster filter. Switching to it later
    goes through on_clustering_change, which rebuilds anyway."""
    from palms.tabs.tab_gene_analysis import refresh_cluster_labels_in_coloring_tab

    ctx, calls = _fake_ctx("leiden_r1.0")
    refresh_cluster_labels_in_coloring_tab(ctx, "leiden_r0.5")
    assert calls == []


@pytest.mark.parametrize("missing", ["repopulate_cluster_checkboxes",
                                     "clustering_widget"])
def test_a_half_wired_context_does_not_raise(missing):
    """These run from a worker callback, and ViewerContext declares both as None.

    A tab built before the Coloring tab has wired them would otherwise turn a
    successful annotation into a TypeError the user sees instead of their labels.
    """
    from palms.tabs.tab_gene_analysis import refresh_cluster_labels_in_coloring_tab

    ctx = SimpleNamespace(
        repopulate_cluster_checkboxes=lambda: None,
        clustering_widget=SimpleNamespace(value="leiden_r1.0"),
    )
    setattr(ctx, missing, None)
    refresh_cluster_labels_in_coloring_tab(ctx, "leiden_r1.0")


def _label_writers(tree):
    """Every function that assigns into state['cluster_labels'][...]."""
    writers = {}
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for sub in ast.walk(node):
            if not isinstance(sub, ast.Assign):
                continue
            for target in sub.targets:
                if (isinstance(target, ast.Subscript)
                        and isinstance(target.value, ast.Subscript)
                        and isinstance(target.value.slice, ast.Constant)
                        and target.value.slice.value == "cluster_labels"):
                    writers.setdefault(node.name, node)
    return writers


def test_every_handler_that_names_clusters_announces_it():
    """Derived from the source, so a fourth annotation path cannot forget."""
    tree = ast.parse(TAB.read_text())
    writers = _label_writers(tree)
    assert len(writers) >= 3, (
        f"expected at least the three annotation handlers, found {sorted(writers)}")

    missing = []
    for name, node in writers.items():
        called = {
            n.func.id for n in ast.walk(node)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
        }
        if REFRESH not in called:
            missing.append(name)
    assert not missing, (
        f"{missing} write cluster labels without calling {REFRESH}(ctx, key) — "
        f"the names will not reach the Coloring tab until the clustering is "
        f"reselected")


def test_the_refresh_is_not_hidden_behind_a_branch():
    """A conditional call is the same omission with extra steps."""
    tree = ast.parse(TAB.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        for sub in ast.walk(node):
            if (isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name)
                    and sub.func.id == REFRESH):
                pytest.fail(
                    f"line {sub.lineno}: {REFRESH} is inside a conditional; the "
                    f"decision belongs in the helper, where it is tested")
