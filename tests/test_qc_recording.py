"""What applying and reverting a QC filter does to the provenance graph.

A filter starts a **second lineage**, it does not revise the first. Results
recorded before it were computed on every cell and still were; the graph keeps
saying so (``deps=["preamble"]``, fresh), and work recorded after it roots at
``qc_filter``. Two mechanisms make that honest, and both are checked here.

**Order.** ``qc_filter`` rebinds ``adata``, which every earlier cell-rooted
step read, and ``topo_sort`` breaks ties by ``(kind, id)`` -- SETUP first, so
the filter would otherwise sort ahead of the unfiltered clusterings and the
notebook would run them on the wrong cells. The node is a *barrier*
(``Step(barrier=True)``): the sort places it after every node that is not its
descendant, from the flag alone, so no history has to be rewritten to get the
notebook right.

**Ids.** ``normalize``, ``spatial_neighbors`` and ``roi_deg`` carry no key, so
a filtered run would upsert the node the unfiltered results depend on and flag
them stale for nothing. ``cell_scoped_id`` gives the filtered lineage its own
(``normalize:qc``); a source guard makes sure every dependent asks for it.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")
pytest.importorskip("scanpy")
pytest.importorskip("qtpy")

SRC = Path(__file__).resolve().parent.parent / "src" / "palms"
TABS = SRC / "tabs"


@pytest.fixture
def ctx(tmp_path, qapp):
    """A ViewerContext with the shared helpers bound and a real table."""
    from palms.loader import label_to_obs_for
    from palms.tabs._helpers import create_shared_helpers
    from palms.utils.viewer_context import ViewerContext

    rng = np.random.default_rng(0)
    n, g = 60, 12
    X = sparse.csr_matrix(rng.poisson(0.6, (n, g)).astype("float32"))
    adata = anndata.AnnData(
        X=X,
        obs=pd.DataFrame({"cell_id": np.arange(1, n + 1)},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(g)]),
    )
    adata.obsm["spatial"] = rng.uniform(0, 100, (n, 2))

    context = ViewerContext(
        data_path=tmp_path,
        adata=adata,
        full_adata=adata,
        label_to_obs=label_to_obs_for(adata),
        full_label_to_obs=label_to_obs_for(adata),
        gene_names=list(adata.var_names),
        state={"record_code": True, "code_journal": [],
               "prov_graph_restored": True},
    )
    create_shared_helpers(context)
    context.record_preamble()
    return context


def _apply(ctx, min_counts=5, min_cells=20):
    from palms.tabs._helpers import qc_filter_preview
    ctx.apply_qc_filter(qc_filter_preview(min_counts, min_cells))


def test_applying_records_a_setup_node_rooted_at_the_preamble(ctx):
    from palms.utils.prov_graph import SETUP

    _apply(ctx)
    node = ctx.state["prov_graph"].get("qc_filter")
    assert node is not None
    assert node.kind == SETUP
    assert node.deps == ["preamble"]
    assert "sc.pp.filter_cells" in node.code
    assert node.stale is False


def test_applying_narrows_the_bound_table(ctx):
    full = ctx.full_adata
    _apply(ctx)
    assert ctx.adata is not full
    assert ctx.adata.n_obs < full.n_obs
    assert ctx.full_adata is full, "the unfiltered pair must not move"
    assert len(ctx.label_to_obs) == len(ctx.full_label_to_obs)


def test_normalize_under_the_filter_is_its_own_node(ctx):
    _apply(ctx)
    ctx.ensure_normalized()

    graph = ctx.state["prov_graph"]
    assert "normalize" not in graph
    assert graph.get("normalize:qc").deps == ["qc_filter"]
    order = graph.topo_sort()
    assert order.index("preamble") < order.index("qc_filter") < order.index("normalize:qc")


def test_the_filter_node_is_a_barrier(ctx):
    _apply(ctx)
    assert ctx.state["prov_graph"].get("qc_filter").barrier is True


def test_the_filter_runs_on_the_full_table_not_the_current_one(ctx):
    """Re-running looser must bring cells back; a filter of a filter cannot."""
    _apply(ctx, min_counts=5, min_cells=20)
    narrow = ctx.adata.n_obs
    _apply(ctx, min_counts=1, min_cells=1)
    assert ctx.adata.n_obs == ctx.full_adata.n_obs > narrow


def _two_lineages(ctx):
    """r1.0 on every cell, then a filter, then r0.9 under it -- the crop_6 session."""
    ctx.ensure_normalized()
    ctx.record_node("clustering:leiden_r1.0", "adata.obs['leiden_r1.0'] = 1",
                    deps=[ctx.cell_scoped_id("normalize")])
    _apply(ctx)
    ctx.ensure_normalized()
    ctx.record_node("clustering:leiden_r0.9", "adata.obs['leiden_r0.9'] = 1",
                    deps=[ctx.cell_scoped_id("normalize")])
    return ctx.state["prov_graph"]


def test_applying_leaves_results_recorded_before_it_fresh(ctx):
    """They were computed on every cell, and they still were.

    The first design re-pointed them at the filter and forced them stale; the
    user's reading -- two lineages, both true -- is the one the graph keeps.
    """
    graph = _two_lineages(ctx)

    r1 = graph.get("clustering:leiden_r1.0")
    assert r1.deps == ["normalize"]
    assert r1.stale is False
    assert graph.get("normalize").deps == ["preamble"]
    assert graph.get("normalize").stale is False
    assert graph.get("normalize:qc").deps == ["qc_filter"]
    assert graph.get("clustering:leiden_r0.9").deps == ["normalize:qc"]


def test_the_notebook_runs_unfiltered_work_before_the_filter(ctx):
    """The barrier at work: r1.0 sorts before ``qc_filter`` although SETUP
    sorts first and nothing names r1.0 as the filter's dependency."""
    from palms.utils.prov_graph import graph_to_cells

    graph = _two_lineages(ctx)
    ids = [c.node_id for c in graph_to_cells(graph) if c.cell_type == "code"]
    expect = ["normalize", "clustering:leiden_r1.0", "qc_filter",
              "normalize:qc", "clustering:leiden_r0.9"]
    assert [i for i in ids if i in expect] == expect


def test_changing_a_cutoff_flags_only_the_filtered_lineage(ctx):
    graph = _two_lineages(ctx)

    _apply(ctx, min_counts=9, min_cells=20)

    assert graph.get("normalize:qc").stale is True
    assert graph.get("clustering:leiden_r0.9").stale is True
    assert graph.get("normalize").stale is False
    assert graph.get("clustering:leiden_r1.0").stale is False


def test_reverting_removes_an_unused_node(ctx):
    """There is no code for "un-filter": a notebook without one never filtered."""
    _apply(ctx)
    graph = ctx.state["prov_graph"]
    assert "qc_filter" in graph

    ctx.clear_qc_filter()

    assert "qc_filter" not in graph
    assert ctx.adata is ctx.full_adata
    assert ctx.state["qc_filter"] is None


def test_reverting_keeps_a_node_something_was_computed_under(ctx):
    """r0.9 exists and depends on the filter, so its step must stay -- and
    nothing goes stale, because nothing about either lineage changed."""
    graph = _two_lineages(ctx)

    ctx.clear_qc_filter()

    assert "qc_filter" in graph
    assert graph.get("clustering:leiden_r0.9").deps == ["normalize:qc"]
    assert not [n.id for n in graph.nodes() if n.stale]
    assert ctx.adata is ctx.full_adata
    assert ctx.state["qc_filter"] is None
    assert ctx.cell_root() == "preamble"


def test_work_after_a_revert_roots_above_the_barrier(ctx):
    """No restore node: post-revert work depends on ``preamble`` and the sort
    puts it before the filter, where ``adata`` is still the full table."""
    graph = _two_lineages(ctx)
    ctx.clear_qc_filter()
    ctx.ensure_normalized()
    ctx.record_node("clustering:leiden_r0.5", "adata.obs['leiden_r0.5'] = 1",
                    deps=[ctx.cell_scoped_id("normalize")])

    assert graph.get("clustering:leiden_r0.5").deps == ["normalize"]
    order = graph.topo_sort()
    assert order.index("clustering:leiden_r0.5") < order.index("qc_filter")
    assert not [n.id for n in graph.nodes() if n.stale]


def test_revert_then_reapply_flags_nothing(ctx):
    graph = _two_lineages(ctx)
    ctx.clear_qc_filter()

    _apply(ctx)

    assert not [n.id for n in graph.nodes() if n.stale]
    assert graph.get("qc_filter").barrier is True
    assert ctx.cell_root() == "qc_filter"


def test_reverting_without_a_filter_is_a_no_op(ctx):
    graph_before = len(ctx.state["prov_graph"])
    ctx.clear_qc_filter()
    assert len(ctx.state["prov_graph"]) == graph_before


def test_ensure_is_idempotent_and_records_nothing_new(ctx):
    """What a launch calls: a restored session must not go stale for nothing."""
    _apply(ctx)
    ctx.ensure_normalized()
    graph = ctx.state["prov_graph"]
    before = {n.id: (n.code, list(n.deps), n.stale) for n in graph.nodes()}
    filtered = ctx.adata

    ctx.ensure_qc_filter()

    assert ctx.adata is filtered, "an already-applied filter must not re-run"
    assert {n.id: (n.code, list(n.deps), n.stale)
            for n in graph.nodes()} == before


def test_ensure_applies_a_filter_restored_from_a_session(ctx):
    ctx.state["qc_filter"] = {"min_counts": 5, "min_cells": 20}
    ctx.ensure_qc_filter()
    assert "qc_filter" in ctx.state["prov_graph"]
    assert ctx.adata.n_obs < ctx.full_adata.n_obs


def test_cell_root_answers_preamble_until_a_filter_exists(ctx):
    assert ctx.cell_root() == "preamble"
    _apply(ctx)
    assert ctx.cell_root() == "qc_filter"
    ctx.clear_qc_filter()
    assert ctx.cell_root() == "preamble"


# ── The scoped ids, checked at the source ────────────────────────────────────

_SCOPED = ("normalize", "spatial_neighbors", "roi_deg")


def _literal_scoped_uses() -> list[str]:
    """``Step(...)`` / ``record_node(...)`` sites naming a scoped id literally.

    An id or dep written as ``"normalize"`` is the unfiltered node whatever
    cell set is bound, so a step recorded under a filter would depend on --
    or, worse, revise -- the node the unfiltered results hang from. The dep
    has to come from ``ctx.cell_scoped_id``.
    """
    hits: list[str] = []
    for path in sorted(TABS.glob("*.py")):
        module = ast.parse(path.read_text())
        for node in ast.walk(module):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
            if name not in ("Step", "record_node"):
                continue
            checked = [kw.value for kw in node.keywords if kw.arg in ("deps", "id")]
            if node.args:
                checked.append(node.args[0])
            for expr in checked:
                for const in _constants_outside_scoped_calls(expr):
                    if const in _SCOPED:
                        hits.append(f"{path.name}:{node.lineno} {const!r}")
    return hits


def _constants_outside_scoped_calls(expr):
    """String constants in *expr*, not descending into ``cell_scoped_id(...)``."""
    if (isinstance(expr, ast.Call)
            and (getattr(expr.func, "attr", None) == "cell_scoped_id"
                 or getattr(expr.func, "id", None) == "_cell_scoped_id")):
        return
    if isinstance(expr, ast.Constant) and isinstance(expr.value, str):
        yield expr.value
    for child in ast.iter_child_nodes(expr):
        yield from _constants_outside_scoped_calls(child)


def test_no_step_names_a_scoped_id_literally():
    hits = _literal_scoped_uses()
    assert not hits, (
        "these sites write a cell-scoped id as a literal; use "
        f"ctx.cell_scoped_id(...) so the filtered lineage gets its own node: {hits}"
    )


def test_cell_scoped_id_follows_the_filter(ctx):
    assert ctx.cell_scoped_id("normalize") == "normalize"
    _apply(ctx)
    assert ctx.cell_scoped_id("normalize") == "normalize:qc"
    assert ctx.cell_scoped_id("roi_deg") == "roi_deg:qc"
    ctx.clear_qc_filter()
    assert ctx.cell_scoped_id("normalize") == "normalize"


def test_spatial_neighbors_is_scoped_with_its_dep(ctx):
    _apply(ctx)
    ctx.ensure_spatial_neighbors(6)
    graph = ctx.state["prov_graph"]
    assert graph.get("spatial_neighbors:qc").deps == ["normalize:qc"]
    assert "spatial_neighbors" not in graph


def test_no_template_reads_the_stores_table_directly():
    """``sdata`` keeps the *full* table after a rebind; ``adata`` is the subset.

    A template reaching into ``sdata["table"]`` would therefore quietly opt out
    of the filter — and of a custom segmentation before it.
    """
    builtin = SRC / "utils" / "step_templates" / "builtin"
    offenders = [
        p.name for p in sorted(builtin.glob("*.tmpl"))
        if 'sdata["table"]' in p.read_text() or "sdata['table']" in p.read_text()
    ]
    assert not offenders, (
        f"template(s) {offenders} read the store's table directly, bypassing "
        f"whatever narrowed `adata`"
    )


def test_a_failed_apply_leaves_the_previous_state_intact(ctx, monkeypatch):
    """The intermediate state is the dangerous one, so it must not survive.

    Mid-apply ``adata`` is the full table while ``label_to_obs`` still describes
    the previous filter — a mismatch that paints cells with other cells' values
    instead of raising.

    Injected at ``StepExecutor.run`` rather than at ``ctx.run_step``: the helper
    calls its own closure, so patching the public attribute would prove nothing.
    """
    from palms.utils.steps import StepExecutor

    _apply(ctx)
    filtered, filtered_map = ctx.adata, ctx.label_to_obs

    def _boom(self, step, progress=None):
        raise RuntimeError("template blew up")

    monkeypatch.setattr(StepExecutor, "run", _boom)

    with pytest.raises(RuntimeError):
        _apply(ctx, min_counts=9)

    assert ctx.adata is filtered
    assert ctx.label_to_obs is filtered_map
    assert ctx.executor.ns["adata"] is filtered
