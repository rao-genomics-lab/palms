"""Unit tests for the provenance graph (pure, no Qt/napari/nbformat).

Run standalone:   python tests/test_prov_graph.py
Or with pytest:   pytest tests/test_prov_graph.py
"""
from __future__ import annotations

import sys

import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from palms.utils.prov_graph import (  # noqa: E402
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


def test_rerunning_a_step_clears_its_own_stale_flag_even_when_the_code_is_identical():
    """Re-running is what the ⚠ badge tells you to do; it has to work.

    Found on a real dataset (crop_6): `clustering:cnv_leiden_res0.2` publishes
    the CNV labels onto the table with a snippet that carries none of the run's
    parameters, so re-running inferCNV with a different reference changed
    `cnv:infercnv` (staling the child) and then re-recorded the child with
    byte-identical code. The early return skipped `stale = False`, so that node
    was stuck stale permanently with no way to clear it.
    """
    g = ProvGraph()
    g.upsert("cnv:infercnv", "CODE_A", kind=SETUP)
    g.upsert("clustering:cnv", "PUBLISH", deps=["cnv:infercnv"])

    g.upsert("cnv:infercnv", "CODE_B", kind=SETUP)          # different settings
    assert g.get("clustering:cnv").stale is True

    g.upsert("clustering:cnv", "PUBLISH", deps=["cnv:infercnv"])   # same code
    assert g.get("clustering:cnv").stale is False


def test_an_identical_rerun_still_does_not_stale_descendants():
    """Only the node's own flag changes — its output did not."""
    g = ProvGraph()
    g.upsert("a", "A", kind=SETUP)
    g.upsert("b", "B", deps=["a"])
    g.upsert("c", "C", deps=["b"])
    g.upsert("b", "B", deps=["a"])          # identical re-run of the middle node
    assert g.get("c").stale is False


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
    from palms.utils.prov_graph import NOTES_MARKER

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


# ── barriers ──────────────────────────────────────────────────────────────────
# A barrier rebinds a name that nodes outside its lineage read. The readers do
# not depend on it, so a deps edge cannot say "run them first"; the barrier
# says it for them by sorting after every node that is not its descendant.

def _two_lineages() -> ProvGraph:
    """Everything recorded before the filter, then the filter, then a second
    lineage under it — plus one node recorded *after* the filter but rooted
    above it, which is the case a recording-order tie-break would get wrong."""
    g = ProvGraph()
    g.upsert("environment", "# versions", kind=SETUP)
    g.upsert("preamble", "adata = sdata['table'].copy()", kind=SETUP)
    g.upsert("normalize", "adata_norm = adata.copy()", deps=["preamble"], kind=SETUP)
    g.upsert("clustering:leiden_r1.0", "sc.tl.leiden(adata_norm)", deps=["normalize"])
    g.upsert("plot:umap:leiden_r1.0", "sc.pl.umap(adata_norm)",
             deps=["clustering:leiden_r1.0"], kind=TERMINAL)
    g.upsert("viewer:background", "# canvas", kind=NOTE)
    g.upsert("qc_filter", "adata = adata[keep].copy()", deps=["preamble"],
             kind=SETUP, barrier=True)
    g.upsert("normalize:qc", "adata_norm = adata.copy()", deps=["qc_filter"], kind=SETUP)
    g.upsert("clustering:leiden_r0.9", "sc.tl.leiden(adata_norm)", deps=["normalize:qc"])
    g.upsert("clustering:leiden_r0.8", "sc.tl.leiden(adata_norm)", deps=["normalize"])
    return g


def test_a_barrier_sorts_after_every_node_that_is_not_its_descendant():
    order = _two_lineages().topo_sort()
    at = order.index("qc_filter")
    for nid in ("environment", "preamble", "normalize", "clustering:leiden_r1.0",
                "plot:umap:leiden_r1.0", "viewer:background",
                "clustering:leiden_r0.8"):
        assert order.index(nid) < at, f"{nid} reads the old binding; it must precede the filter"
    for nid in ("normalize:qc", "clustering:leiden_r0.9"):
        assert order.index(nid) > at


def test_a_node_recorded_after_the_barrier_but_rooted_above_it_still_precedes_it():
    """The reason the rule is structural rather than a recording-order tie-break."""
    order = _two_lineages().topo_sort()
    assert order.index("clustering:leiden_r0.8") < order.index("qc_filter")


def test_a_barrier_flags_nothing_when_inserted():
    g = _two_lineages()
    assert not any(n.stale for n in g.nodes())


def test_the_barrier_keeps_the_ordering_invariance():
    """Same nodes recorded in a different order → the same notebook."""
    a = _two_lineages()
    b = ProvGraph()
    for nid in ("preamble", "qc_filter", "normalize", "environment", "normalize:qc",
                "clustering:leiden_r0.9", "clustering:leiden_r1.0",
                "clustering:leiden_r0.8", "viewer:background", "plot:umap:leiden_r1.0"):
        n = a.get(nid)
        b.upsert(nid, n.code, deps=n.deps, kind=n.kind, barrier=n.barrier)
    assert a.topo_sort() == b.topo_sort()


def test_two_unordered_barriers_are_refused_with_a_reason():
    g = _two_lineages()
    g.upsert("segmentation_swap", "adata = tables['custom_table']",
             deps=["preamble"], kind=SETUP, barrier=True)
    with pytest.raises(CycleError, match="unordered"):
        g.topo_sort()


def test_dep_ordered_barriers_chain():
    g = _two_lineages()
    g.upsert("qc_filter:2", "adata = adata[keep2].copy()", deps=["qc_filter"],
             kind=SETUP, barrier=True)
    g.upsert("clustering:leiden_r0.7", "sc.tl.leiden(adata)", deps=["qc_filter:2"])
    order = g.topo_sort()
    assert order.index("clustering:leiden_r0.9") < order.index("qc_filter:2") \
        < order.index("clustering:leiden_r0.7")


def test_barrier_survives_a_round_trip_and_defaults_off_for_old_graphs():
    g = _two_lineages()
    back = ProvGraph.from_list(g.to_list())
    assert back.get("qc_filter").barrier is True
    assert back.get("normalize").barrier is False
    assert back.topo_sort() == g.topo_sort()

    old = [{"id": "preamble", "code": "x"}, {"id": "n", "code": "y", "deps": ["preamble"]}]
    assert all(n.barrier is False for n in ProvGraph.from_list(old).nodes())


def test_changing_only_the_barrier_flag_stales_nothing():
    """Where a node sorts is not what it computes."""
    g = _two_lineages()
    g.upsert("qc_filter", "adata = adata[keep].copy()", deps=["preamble"],
             kind=SETUP, barrier=False)
    assert g.get("clustering:leiden_r0.9").stale is False
    assert g.get("qc_filter").barrier is False
