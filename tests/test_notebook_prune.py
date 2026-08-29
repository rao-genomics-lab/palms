"""Dropping stale nodes out of the provenance graph.

``plan_prune`` and ``describe_prune`` are module-level and viewer-free — the
idiom ``test_template_provenance.py`` uses for ``reconcile_edits`` — so the
correctness question can be answered without a napari viewer or a Qt event loop.

The property that matters is the last one: applying a plan through the real
``ProvGraph.remove`` must never raise. ``remove`` refuses any node another still
depends on, and the stale set is *not* dependency-closed (``upsert`` clears
``stale`` on the node it re-records while flagging that node's descendants), so
"remove everything stale" is wrong in two different ways at once. Also the first
coverage ``ProvGraph.remove`` has had.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

pytest.importorskip("qtpy")

from palms.tabs import tab_notebook  # noqa: E402
from palms.utils.prov_graph import SETUP, TERMINAL, ProvGraph  # noqa: E402


def _graph():
    g = ProvGraph()
    g.upsert("preamble", "data_path = ...", kind=SETUP)
    g.upsert("normalize", "normalize(adata)", deps=["preamble"], kind=SETUP)
    g.upsert("clustering:a", "leiden(...)", deps=["normalize"])
    g.upsert("rank_genes:a", "rank(...)", deps=["normalize", "clustering:a"])
    g.upsert("plot:rank_panel:a", "sc.pl.rank_genes_groups(...)",
             deps=["rank_genes:a"], kind=TERMINAL)
    return g


# ── Nothing to do ────────────────────────────────────────────────────────────

def test_a_missing_or_empty_graph_plans_nothing():
    assert tab_notebook.plan_prune(None).is_empty
    assert tab_notebook.plan_prune(ProvGraph()).is_empty


def test_a_graph_with_no_stale_nodes_plans_nothing():
    plan = tab_notebook.plan_prune(_graph())
    assert plan.is_empty and plan.blocked == ()


# ── Ordering ─────────────────────────────────────────────────────────────────

def test_stale_nodes_are_removed_leaves_first():
    """ProvGraph.remove refuses a node with dependents, so order is not cosmetic."""
    g = _graph()
    g.upsert("normalize", "normalize(adata, target_sum=100)",
             deps=["preamble"], kind=SETUP)
    plan = tab_notebook.plan_prune(g)
    order = list(plan.remove)
    assert set(order) == {"clustering:a", "rank_genes:a", "plot:rank_panel:a"}
    assert order.index("plot:rank_panel:a") < order.index("rank_genes:a")
    assert order.index("rank_genes:a") < order.index("clustering:a")


def test_the_re_recorded_node_itself_is_not_dropped():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    assert "normalize" not in tab_notebook.plan_prune(g).remove
    assert "preamble" not in tab_notebook.plan_prune(g).remove


# ── A fresh node holding a stale one ─────────────────────────────────────────

def _graph_with_a_fresh_dependent():
    """The shape "remove everything stale" gets wrong.

    normalize is re-recorded, staling clustering:a and its descendants; a new
    step is then recorded on top of the stale clustering, so it is fresh and
    depends on a stale node.
    """
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    g.upsert("nhood:a", "sq.gr.nhood_enrichment(...)", deps=["clustering:a"])
    return g


def test_a_stale_node_a_fresh_step_depends_on_is_kept():
    plan = tab_notebook.plan_prune(_graph_with_a_fresh_dependent())
    assert "clustering:a" not in plan.remove
    assert ("clustering:a", ("nhood:a",)) in plan.blocked


def test_the_rest_of_the_stale_set_is_still_dropped():
    plan = tab_notebook.plan_prune(_graph_with_a_fresh_dependent())
    assert set(plan.remove) == {"rank_genes:a", "plot:rank_panel:a"}


def test_a_stale_node_held_only_transitively_is_kept():
    """The fresh node's hold reaches through its whole ancestry, not one hop."""
    g = _graph()
    g.upsert("preamble", 'data_path = Path(r"/elsewhere")', kind=SETUP)
    g.upsert("nhood:a", "sq.gr.nhood_enrichment(...)", deps=["clustering:a"])
    plan = tab_notebook.plan_prune(g)
    blocked = dict(plan.blocked)
    # nhood:a is fresh and reaches normalize through clustering:a.
    assert "normalize" in blocked and "clustering:a" in blocked
    assert "normalize" not in plan.remove


# ── The property: a plan is always applicable ────────────────────────────────

@pytest.mark.parametrize("build", [
    _graph_with_a_fresh_dependent,
    lambda: _graph(),
])
def test_applying_a_plan_never_raises(build):
    g = build()
    plan = tab_notebook.plan_prune(g)
    for node_id in plan.remove:
        g.remove(node_id)          # raises ValueError if the order is wrong
    assert all(nid not in g for nid in plan.remove)
    # What is left is still a coherent DAG: every dep still resolves.
    for node in g.nodes():
        for dep in node.deps:
            assert dep in g, f"{node.id} was left depending on a removed {dep}"
    g.topo_sort()


def test_a_pruned_graph_can_still_be_serialised_and_reloaded():
    g = _graph_with_a_fresh_dependent()
    for node_id in tab_notebook.plan_prune(g).remove:
        g.remove(node_id)
    restored = ProvGraph.from_list(g.to_list())
    assert {n.id for n in restored.nodes()} == {n.id for n in g.nodes()}


# ── ProvGraph.remove itself, previously untested ─────────────────────────────

def test_remove_refuses_a_node_something_still_depends_on():
    g = _graph()
    with pytest.raises(ValueError, match="still required by"):
        g.remove("clustering:a")
    assert "clustering:a" in g


def test_remove_of_an_unknown_id_is_a_no_op():
    g = _graph()
    before = len(g)
    g.remove("nothing:here")
    assert len(g) == before


# ── The dialog body ──────────────────────────────────────────────────────────

def test_the_dialog_names_the_orphaned_results_and_where_to_clear_them():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    plan = tab_notebook.plan_prune(g)
    text = tab_notebook.describe_prune(
        plan, {"rank_genes:a": ("uns:table/rank_genes_a",)})
    assert "uns:table/rank_genes_a" in text
    assert "Select Stale Results" in text
    assert "backup" in text.lower()


def test_the_dialog_names_what_it_kept_and_why():
    plan = tab_notebook.plan_prune(_graph_with_a_fresh_dependent())
    text = tab_notebook.describe_prune(plan, {})
    assert "clustering:a" in text and "nhood:a" in text


def test_the_dialog_says_when_it_could_not_check_for_stored_results():
    """An empty `orphans` means "none" or "did not look" — they must not read alike.

    The prune reuses the Dataset tab's scan rather than walking the store itself:
    build_inventory sizes every file under a store that is routinely tens of GB,
    and doing that in a button callback would freeze the GUI to decorate a dialog.
    """
    plan = tab_notebook.plan_prune(_graph_with_a_fresh_dependent())
    unknown = tab_notebook.describe_prune(plan, {}, inventory_known=False)
    assert "was not checked" in unknown and "Scan Dataset" in unknown
    assert "was not checked" not in tab_notebook.describe_prune(plan, {})


def test_the_prune_dialog_does_not_walk_the_dataset_itself():
    """A source guard: the freeze is invisible in a test with a tiny fixture.

    Parsed rather than grepped, so the function's own docstring explaining why
    it must not call build_inventory does not itself trip the check — the idiom
    test_tab_dataset.py uses for the _remove_tree callers.
    """
    import ast
    tree = ast.parse(Path(tab_notebook.__file__).read_text())
    target = next(
        n for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef) and n.name == "_stale_artifacts")
    called = {
        node.func.attr if isinstance(node.func, ast.Attribute) else
        getattr(node.func, "id", "")
        for node in ast.walk(target) if isinstance(node, ast.Call)
    }
    assert "build_inventory" not in called
    assert "_dataset_sections" in ast.unparse(target)


# ── Source guard: the prune must write both records ──────────────────────────

def test_the_prune_updates_the_session_attr_as_well_as_the_sidecar():
    """Both monotonic guards refuse a graph smaller than the stored one.

    app._load_prov_graph_items and session._build_session_attrs each keep the
    larger copy, so a prune that writes only the sidecar is silently undone at
    the next launch. Nothing in a single-process test can catch that, hence a
    source guard.
    """
    source = Path(tab_notebook.__file__).read_text()
    body = source.split("def _persist_pruned_graph")[1].split("\n    def ")[0]
    assert "save_prov_graph" in body
    assert "safe_group_update" in body
    assert '"prov_graph"' in body


def test_the_prune_backs_the_graph_up_before_removing_anything():
    source = Path(tab_notebook.__file__).read_text()
    handler = source.split("def _on_prune_stale")[1].split("\n    def ")[0]
    assert handler.index("_backup_graph") < handler.index("graph.remove")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
