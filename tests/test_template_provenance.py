"""Where a node's code came from, and whether that survives a round-trip.

``ProvNode.code`` always records what actually ran, so replay stays correct no
matter what these fields say. What they add is the ability for a *reader* — the
verification report, someone opening the notebook a year later — to tell a stock
run from a customised one, which the rendered source alone cannot reveal.

Two properties matter most and are easy to break:

* **Old graphs must still load.** ``prov_graph.json`` sidecars and zarr session
  attrs written before these fields existed are on disk right now.
* **Template metadata must not invalidate anything.** It describes where the
  same code came from, so changing it must not mark descendants stale — that
  would put a spurious ⚠ on every downstream result.
"""

from __future__ import annotations

import json

import pytest

from xenium_viewer.utils.prov_graph import (
    ARTIFACT,
    TEMPLATE_BUILTIN,
    TEMPLATE_HAND_EDITED,
    TEMPLATE_ORIGINS,
    TEMPLATE_USER,
    ProvGraph,
)
from xenium_viewer.utils.steps import Step, StepExecutor


# ── backward compatibility ───────────────────────────────────────────────────

def test_a_graph_written_before_the_fields_existed_still_loads():
    """The exact shape ``to_list`` produced before this feature."""
    legacy = [
        {"id": "preamble", "code": "import scanpy as sc", "deps": [],
         "kind": "setup", "label": None, "params": {}, "stale": False, "seq": 1},
        {"id": "clustering:k", "code": "x = 1", "deps": ["preamble"],
         "kind": "artifact", "label": "Clustering: k", "params": {"r": 1.0},
         "stale": False, "seq": 2},
    ]
    graph = ProvGraph.from_list(legacy)

    assert len(graph) == 2
    for node in graph.nodes():
        assert node.template_origin == TEMPLATE_BUILTIN, (
            "a graph recorded before templates could be customised is stock by "
            "definition — there was no way for it not to be"
        )
        assert node.template_id is None
        assert node.template_hash is None


def test_a_graph_with_unknown_future_fields_still_loads():
    """``from_list`` enumerates keys rather than splatting, so this must work."""
    forward = [{"id": "a", "code": "x = 1", "kind": "setup",
                "template_origin": TEMPLATE_USER,
                "some_field_from_the_future": {"nested": True}}]
    graph = ProvGraph.from_list(forward)
    assert graph.get("a").template_origin == TEMPLATE_USER


def test_the_fields_round_trip_through_json():
    graph = ProvGraph()
    graph.upsert("a", "x = 1", kind="setup",
                 template_id="clustering.leiden",
                 template_origin=TEMPLATE_USER, template_hash="deadbeef")

    reloaded = ProvGraph.from_list(json.loads(json.dumps(graph.to_list())))
    node = reloaded.get("a")
    assert node.template_id == "clustering.leiden"
    assert node.template_origin == TEMPLATE_USER
    assert node.template_hash == "deadbeef"


def test_a_fresh_node_defaults_to_builtin():
    graph = ProvGraph()
    node = graph.upsert("a", "x = 1", kind="setup")
    assert node.template_origin == TEMPLATE_BUILTIN
    assert node.template_origin in TEMPLATE_ORIGINS


# ── staleness ────────────────────────────────────────────────────────────────

def test_changing_only_template_metadata_does_not_mark_descendants_stale():
    """It describes where the same code came from — the artifact is unaffected."""
    graph = ProvGraph()
    graph.upsert("a", "x = 1", kind="setup")
    graph.upsert("b", "y = x + 1", deps=["a"])
    assert not graph.get("b").stale

    graph.upsert("a", "x = 1", kind="setup",
                 template_id="setup.a", template_origin=TEMPLATE_USER,
                 template_hash="abc123")

    assert graph.get("a").template_origin == TEMPLATE_USER
    assert not graph.get("b").stale, (
        "template metadata entered the `unchanged` comparison; a change of "
        "origin must not invalidate downstream results"
    )


def test_changing_the_code_still_marks_descendants_stale():
    """Guard the guard: the staleness machinery itself must still work."""
    graph = ProvGraph()
    graph.upsert("a", "x = 1", kind="setup")
    graph.upsert("b", "y = x + 1", deps=["a"])

    graph.upsert("a", "x = 2", kind="setup", template_origin=TEMPLATE_USER)
    assert graph.get("b").stale


# ── the Step → node path ─────────────────────────────────────────────────────

def test_a_step_stamps_its_template_metadata_onto_the_node():
    ex = StepExecutor(namespace={})
    ex.run(Step(
        id="a", template="x = $n", params={"n": 1}, kind="setup",
        template_id="demo.step", template_origin=TEMPLATE_USER,
        template_hash="cafe",
    ))
    node = ex.graph.get("a")
    assert (node.template_id, node.template_origin, node.template_hash) == (
        "demo.step", TEMPLATE_USER, "cafe")


def test_a_step_with_no_template_metadata_records_as_builtin():
    ex = StepExecutor(namespace={})
    ex.run(Step(id="a", template="x = 1", kind="setup"))
    assert ex.graph.get("a").template_origin == TEMPLATE_BUILTIN


# ── the hand-edited case ─────────────────────────────────────────────────────

class _FakeCell:
    """Stands in for a ``NotebookCell`` — the parts ``reconcile_edits`` reads."""

    def __init__(self, node_id, code, edited=True):
        self.node_id = node_id
        self._code = code
        self.edited_by_user = edited

    def get_code(self):
        return self._code


def test_reconcile_edits_marks_the_node_hand_edited():
    """A cell edited in the Notebook tab reaches export without being executed.

    That makes ``node.code`` no longer the source that produced the artifact —
    the one guarantee every other node carries. It must say so.
    """
    pytest.importorskip("qtpy")
    from xenium_viewer.tabs.tab_notebook import reconcile_edits

    graph = ProvGraph()
    graph.upsert("a", "x = 1", kind="setup", template_id="demo.step",
                 template_hash="cafe")

    cell = _FakeCell("a", "x = 999  # typed by hand, never run")
    assert reconcile_edits(graph, [cell]) == ["a"]

    node = graph.get("a")
    assert node.code == "x = 999  # typed by hand, never run"
    assert node.template_origin == TEMPLATE_HAND_EDITED
    assert node.template_hash is None, "hand-typed code derives from no template"
    assert node.template_id == "demo.step", (
        "which template it started from is still worth knowing"
    )
    assert not cell.edited_by_user, "the edit was folded in; don't re-apply it"


def test_reconcile_edits_leaves_untouched_cells_alone():
    pytest.importorskip("qtpy")
    from xenium_viewer.tabs.tab_notebook import reconcile_edits

    graph = ProvGraph()
    graph.upsert("a", "x = 1", kind="setup")

    # Not edited, and a free-form cell with no node — neither may be folded in.
    assert reconcile_edits(graph, [
        _FakeCell("a", "x = 2", edited=False),
        _FakeCell(None, "print('scratch')", edited=True),
    ]) == []

    node = graph.get("a")
    assert node.code == "x = 1"
    assert node.template_origin == TEMPLATE_BUILTIN


def test_reconcile_edits_tolerates_no_graph():
    """Called on every export, including before anything has been recorded."""
    pytest.importorskip("qtpy")
    from xenium_viewer.tabs.tab_notebook import reconcile_edits

    assert reconcile_edits(None, [_FakeCell("a", "x = 1")]) == []
