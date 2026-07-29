"""Unit tests for the provenance graph (pure, no Qt/napari/nbformat).

Run standalone:   python tests/test_prov_graph.py
Or with pytest:   pytest tests/test_prov_graph.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xenium_viewer.utils.prov_graph import (  # noqa: E402
    ProvGraph, CycleError, graph_to_cells, graph_to_script,
    graph_to_mermaid, graph_to_dot,
    SETUP, ARTIFACT, TERMINAL, NOTE,
)


def _pipeline() -> ProvGraph:
    """preamble → normalize → clustering → (spatial_neighbors) → nhood."""
    g = ProvGraph()
    g.upsert("preamble", "import scanpy as sc", kind=SETUP)
    g.upsert("normalize", "sc.pp.normalize_total(adata)", deps=["preamble"], kind=SETUP)
    g.upsert("clustering:leiden_r1.0",
             "sc.tl.leiden(adata, resolution=1.0, random_state=0)",
             deps=["normalize"], label="Leiden clustering")
    g.upsert("spatial_neighbors", "sq.gr.spatial_neighbors(adata)", deps=["preamble"])
    g.upsert("nhood:leiden_r1.0", "sq.gr.nhood_enrichment(adata, 'leiden_r1.0')",
             deps=["clustering:leiden_r1.0", "spatial_neighbors"])
    return g


def test_topo_order_respects_deps():
    g = _pipeline()
    order = g.topo_sort()
    assert order.index("preamble") < order.index("normalize")
    assert order.index("normalize") < order.index("clustering:leiden_r1.0")
    assert order.index("clustering:leiden_r1.0") < order.index("nhood:leiden_r1.0")
    assert order.index("spatial_neighbors") < order.index("nhood:leiden_r1.0")


def test_three_session_rerun_revises_in_place():
    """Session 3 re-clusters the SAME column with new params."""
    g = _pipeline()
    before = g.topo_sort()
    nhood = g.get("nhood:leiden_r1.0")
    assert nhood.stale is False

    # Re-run clustering under the same id, different code (new n_neighbors).
    g.upsert("clustering:leiden_r1.0",
             "sc.pp.neighbors(adata, n_neighbors=30)\n"
             "sc.tl.leiden(adata, resolution=1.0, random_state=0)",
             deps=["normalize"], label="Leiden clustering")

    # Not dropped, not appended at the end — same position, same node count.
    after = g.topo_sort()
    assert len(g) == len(before)
    assert after == before, "re-running a step must not reorder the notebook"
    # Downstream nhood is now flagged stale.
    assert g.get("nhood:leiden_r1.0").stale is True
    # The clustering node itself is fresh and carries the new code.
    assert g.get("clustering:leiden_r1.0").stale is False
    assert "n_neighbors=30" in g.get("clustering:leiden_r1.0").code


def test_rerun_identical_is_noop():
    g = _pipeline()
    g.upsert("nhood:leiden_r1.0", "sq.gr.nhood_enrichment(adata, 'leiden_r1.0')",
             deps=["clustering:leiden_r1.0", "spatial_neighbors"])
    # Re-record clustering identically → descendants must NOT go stale.
    g.upsert("clustering:leiden_r1.0",
             "sc.tl.leiden(adata, resolution=1.0, random_state=0)",
             deps=["normalize"], label="Leiden clustering")
    assert g.get("nhood:leiden_r1.0").stale is False


def test_new_key_is_independent_branch():
    g = _pipeline()
    g.upsert("clustering:leiden_r0.5",
             "sc.tl.leiden(adata, resolution=0.5, random_state=0)",
             deps=["normalize"], label="Leiden clustering (r0.5)")
    # nhood on the original clustering is untouched.
    assert g.get("nhood:leiden_r1.0").stale is False
    # Pruning to just the r0.5 branch drops the r1.0 nhood branch.
    keep = g.ancestors_closure(["clustering:leiden_r0.5"])
    assert "clustering:leiden_r0.5" in keep
    assert "nhood:leiden_r1.0" not in keep
    assert "preamble" in keep and "normalize" in keep


def test_ordering_invariance():
    """Same steps recorded in different wall-clock order → identical notebook."""
    g1 = ProvGraph()
    g1.upsert("preamble", "P", kind=SETUP)
    g1.upsert("normalize", "N", deps=["preamble"], kind=SETUP)
    g1.upsert("clustering:c", "C", deps=["normalize"])
    g1.upsert("spatial_neighbors", "S", deps=["preamble"])
    g1.upsert("nhood:c", "H", deps=["clustering:c", "spatial_neighbors"])

    # Different insertion order (neighbors + clustering swapped, as across sessions).
    g2 = ProvGraph()
    g2.upsert("preamble", "P", kind=SETUP)
    g2.upsert("normalize", "N", deps=["preamble"], kind=SETUP)
    g2.upsert("spatial_neighbors", "S", deps=["preamble"])
    g2.upsert("clustering:c", "C", deps=["normalize"])
    g2.upsert("nhood:c", "H", deps=["clustering:c", "spatial_neighbors"])

    assert graph_to_script(g1) == graph_to_script(g2)


def test_missing_dependency_errors_early():
    g = ProvGraph()
    g.upsert("preamble", "P", kind=SETUP)
    try:
        g.upsert("nhood:x", "H", deps=["clustering:x"])  # dep never recorded
    except KeyError as e:
        assert "unknown node" in str(e)
    else:
        raise AssertionError("expected KeyError for missing dependency")


def test_cycle_detected():
    g = ProvGraph()
    g.upsert("a", "A")
    g.upsert("b", "B", deps=["a"])
    try:
        g.upsert("a", "A2", deps=["b"])  # a→b→a cycle
    except CycleError:
        pass
    else:
        raise AssertionError("expected CycleError")


def test_terminals_excluded_from_export():
    g = _pipeline()
    g.upsert("plot:nhood", "sq.pl.nhood_enrichment(adata, 'leiden_r1.0')",
             deps=["nhood:leiden_r1.0"], kind=TERMINAL, label="Nhood heatmap")
    with_term = graph_to_cells(g, include_terminals=True)
    without = graph_to_cells(g, include_terminals=False)
    ids_with = {c.node_id for c in with_term}
    ids_without = {c.node_id for c in without}
    assert "plot:nhood" in ids_with
    assert "plot:nhood" not in ids_without


def test_serialization_roundtrip():
    g = _pipeline()
    g.upsert("clustering:leiden_r1.0", "CHANGED", deps=["normalize"])  # make nhood stale
    g2 = ProvGraph.from_list(g.to_list())
    assert g2.topo_sort() == g.topo_sort()
    assert g2.get("nhood:leiden_r1.0").stale == g.get("nhood:leiden_r1.0").stale
    assert graph_to_script(g2) == graph_to_script(g)


def test_markdown_header_precedes_code():
    g = ProvGraph()
    g.upsert("preamble", "P", kind=SETUP, label="Setup")
    cells = graph_to_cells(g)
    assert cells[0].cell_type == "markdown" and "Setup" in cells[0].source
    assert cells[1].cell_type == "code" and cells[1].source == "P"


def test_mermaid_contains_every_node_and_edge():
    g = _pipeline()
    m = graph_to_mermaid(g)
    assert m.startswith("flowchart TD")
    # Mermaid uses opaque tokens (n0, n1, ...) for node ids, so labels/ids show
    # up in the node-declaration text. Every node's display text must appear.
    for nid in g.topo_sort():
        node = g.get(nid)
        assert (node.label or nid) in m
    # One directed edge per dep; count the arrows.
    n_deps = sum(len(node.deps) for node in g.nodes())
    assert m.count(" --> ") == n_deps
    # Fresh graph → no stale marker.
    assert "⚠" not in m


def test_mermaid_marks_stale_nodes():
    g = _pipeline()
    # Re-record an upstream node → downstream nhood goes stale.
    g.upsert("clustering:leiden_r1.0", "CHANGED", deps=["normalize"])
    assert g.get("nhood:leiden_r1.0").stale is True
    m = graph_to_mermaid(g)
    assert "⚠" in m
    assert "classDef stale" in m


def test_dot_is_valid_digraph_with_all_edges():
    g = _pipeline()
    d = graph_to_dot(g)
    assert d.startswith("digraph provenance {")
    assert d.rstrip().endswith("}")
    for nid in g.topo_sort():
        node = g.get(nid)
        assert (node.label or nid) in d
    n_deps = sum(len(node.deps) for node in g.nodes())
    assert d.count(" -> ") == n_deps
    # Stale annotation only appears once a node is stale.
    assert "(stale)" not in d
    g.upsert("clustering:leiden_r1.0", "CHANGED", deps=["normalize"])
    assert "(stale)" in graph_to_dot(g)


def _run_all():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\n{len(fns)} tests passed.")


if __name__ == "__main__":
    _run_all()


# ── notes: viewer state, declared rather than disguised as code ──────────────

def _with_note() -> ProvGraph:
    g = _pipeline()
    g.upsert("viewer:background",
             "\n# Viewer background set to white (display only)",
             deps=["preamble"], kind=NOTE, label="Viewer background")
    return g


def test_a_note_exports_as_markdown_not_as_an_empty_code_cell():
    """The defect this kind exists for.

    As a TERMINAL its cell was a comment: it parses to nothing, executes
    successfully and documents nothing, which is indistinguishable from an
    analysis step that forgot to record its code.
    """
    cells = graph_to_cells(_with_note())
    note_cells = [c for c in cells if c.node_id == "viewer:background"]

    assert [c.cell_type for c in note_cells] == ["markdown", "markdown"]
    assert "Viewer background" in note_cells[0].source          # the label header
    assert "background set to white" in note_cells[1].source
    assert not note_cells[1].source.lstrip().startswith("#")    # prose, not code


def test_a_note_is_marked_as_viewer_state():
    """A reader must be able to tell it apart from analysis narrative."""
    from xenium_viewer.utils.prov_graph import NOTES_MARKER

    cells = graph_to_cells(_with_note())
    body = [c for c in cells if c.node_id == "viewer:background"][-1].source
    assert NOTES_MARKER in body


def test_a_note_keeps_its_comment_in_the_flat_script():
    """A .py can carry a comment, and dropping it would lose the record."""
    script = graph_to_script(_with_note())
    assert "# Viewer background set to white (display only)" in script


def test_notes_are_dropped_with_the_terminals():
    without = graph_to_cells(_with_note(), include_terminals=False)
    assert "viewer:background" not in {c.node_id for c in without}


def test_a_note_never_sorts_before_the_analysis_it_annotates():
    order = _with_note().topo_sort()
    assert order[-1] == "viewer:background"


def test_the_diagrams_render_a_note():
    g = _with_note()
    assert "classDef note" in graph_to_mermaid(g)
    assert '"viewer:background"' in graph_to_dot(g)
