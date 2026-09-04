"""What applying and reverting a QC filter does to the provenance graph.

Two things here are easy to get wrong and invisible once wrong.

**Order.** ``_KIND_ORDER`` puts SETUP first and ``topo_sort`` breaks remaining
ties by ``(kind, id)`` — and ``"normalize" < "qc_filter"``. So without a real
edge between them the exported notebook would normalise the *unfiltered* table
and then filter, which runs cleanly and answers a different question.

**Staleness.** ``upsert`` of a new id has no descendants, so inserting the
filter node flags nothing. Everything already recorded would keep its clean
badge while being about a different set of cells. The retroactive re-rooting is
what closes that, and it is the half a reviewer should look at hardest.
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


def test_normalize_depends_on_the_filter_and_sorts_after_it(ctx):
    _apply(ctx)
    ctx.ensure_normalized()

    graph = ctx.state["prov_graph"]
    assert graph.get("normalize").deps == ["qc_filter"]
    order = graph.topo_sort()
    assert order.index("preamble") < order.index("qc_filter") < order.index("normalize")


def test_the_filter_runs_on_the_full_table_not_the_current_one(ctx):
    """Re-running looser must bring cells back; a filter of a filter cannot."""
    _apply(ctx, min_counts=5, min_cells=20)
    narrow = ctx.adata.n_obs
    _apply(ctx, min_counts=1, min_cells=1)
    assert ctx.adata.n_obs == ctx.full_adata.n_obs > narrow


def test_applying_flags_results_recorded_before_it(ctx):
    """The load-bearing half: a new node has no descendants to flag."""
    ctx.ensure_normalized()
    ctx.record_node("clustering:leiden_r1.0", "adata.obs['leiden_r1.0'] = 1",
                    deps=["normalize"])
    graph = ctx.state["prov_graph"]
    assert graph.get("clustering:leiden_r1.0").stale is False

    _apply(ctx)

    assert graph.get("normalize").deps == ["qc_filter"]
    assert graph.get("normalize").stale is True
    assert graph.get("clustering:leiden_r1.0").stale is True


def test_changing_a_cutoff_flags_everything_downstream(ctx):
    _apply(ctx, min_counts=5, min_cells=20)
    ctx.ensure_normalized()
    graph = ctx.state["prov_graph"]
    assert graph.get("normalize").stale is False

    _apply(ctx, min_counts=9, min_cells=20)
    assert graph.get("normalize").stale is True


def test_reverting_removes_the_node_and_re_roots_what_named_it(ctx):
    """There is no code for "un-filter": a notebook without one never filtered."""
    _apply(ctx)
    ctx.ensure_normalized()
    graph = ctx.state["prov_graph"]
    assert "qc_filter" in graph

    ctx.clear_qc_filter()

    assert "qc_filter" not in graph
    assert graph.get("normalize").deps == ["preamble"]
    assert ctx.adata is ctx.full_adata
    assert ctx.state["qc_filter"] is None
    for node in graph.nodes():
        assert "qc_filter" not in node.deps, "a dangling dep would break the graph"


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


# ── The id bridge, checked in both directions ────────────────────────────────

def _literal_cell_root_ids() -> set[str]:
    """Node ids at a ``cell_root()`` dep site, where the id is written literally.

    Not every site is resolvable — ``tab_roi`` builds its deps in a variable and
    ``tab_crop_dataset`` gets its id back from ``crop_export_note`` — so this
    covers the literal ones and :func:`test_every_declared_rule_matches_a_real_id`
    covers the rest from the other end. An f-string id is reduced to the part
    before the first placeholder, which is exactly what the prefixes match on.
    """
    found: set[str] = set()
    for path in sorted(TABS.glob("*.py")):
        module = ast.parse(path.read_text())
        for node in ast.walk(module):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
            if name not in ("Step", "record_node"):
                continue
            deps = next((kw.value for kw in node.keywords if kw.arg == "deps"), None)
            if deps is None or "cell_root()" not in ast.unparse(deps):
                continue
            ident = next((kw.value for kw in node.keywords if kw.arg == "id"), None)
            if ident is None and node.args:
                ident = node.args[0]
            if isinstance(ident, ast.Constant) and isinstance(ident.value, str):
                found.add(ident.value)
            elif isinstance(ident, ast.JoinedStr) and ident.values:
                head = ident.values[0]
                if isinstance(head, ast.Constant) and isinstance(head.value, str):
                    found.add(head.value)
    return found


def test_every_literal_cell_rooted_id_is_recognised():
    """A miss here costs a stale badge that should have been raised."""
    from palms.tabs._helpers import is_cell_rooted

    ids = _literal_cell_root_ids()
    assert ids, "expected to find the call sites that pass ctx.cell_root()"
    unrecognised = sorted(i for i in ids if not is_cell_rooted(i))
    assert not unrecognised, (
        f"node id(s) {unrecognised} declare ctx.cell_root() but is_cell_rooted() "
        f"does not claim them, so applying a filter would leave their results "
        f"looking fresh"
    )


def test_every_declared_rule_matches_a_real_id():
    """The other direction: a rule matching nothing is a rule that has rotted."""
    from palms.tabs._helpers import _CELL_ROOTED_IDS, _CELL_ROOTED_PREFIXES

    source = "\n".join(p.read_text() for p in sorted(SRC.rglob("*.py")))
    for prefix in _CELL_ROOTED_PREFIXES:
        assert f'"{prefix}' in source or f"'{prefix}" in source, (
            f"prefix {prefix!r} matches no recorded node id anywhere in the app"
        )
    for node_id in _CELL_ROOTED_IDS:
        assert f'"{node_id}"' in source or f"'{node_id}'" in source, (
            f"id {node_id!r} matches no recorded node id anywhere in the app"
        )


def test_reverting_moves_a_node_the_id_rules_do_not_claim(ctx):
    """Revert cannot depend on the id list, so it does not.

    A miss on apply costs a stale badge; a miss on revert leaves a dangling dep
    and ``ProvGraph.remove`` refuses. So revert re-points everything naming the
    node it is about to remove, whatever the id rules say about it.
    """
    from palms.tabs._helpers import is_cell_rooted

    _apply(ctx)
    ctx.record_node("something:unclaimed", "pass", deps=[ctx.cell_root()])
    assert not is_cell_rooted("something:unclaimed")

    ctx.clear_qc_filter()

    graph = ctx.state["prov_graph"]
    assert "qc_filter" not in graph
    assert graph.get("something:unclaimed").deps == ["preamble"]


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
