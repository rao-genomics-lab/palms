"""End-to-end: does the *exported notebook* reproduce the viewer's results?

``test_clustering_step.py`` asserts that the recorded source equals the executed
source, and replays it with ``exec`` in the same process. That is the claim
stated; this is the claim *measured*. Here the real ``Step`` templates are run
through ``StepExecutor``, the resulting provenance graph is written out as a
real ``.ipynb`` by ``notebook_export.write_graph_notebook``, and that notebook is
executed in a **clean kernel** — a separate process, no viewer, no session
state, no imports carried over. The two sets of results are then compared.

Tolerances are deliberately absent where equality is achievable: adjusted Rand
index must be exactly 1.0, and the ranked gene names must match exactly. A
tolerance would let a real regression pass as "close enough".

**One documented substitution.** The ``preamble`` node the viewer records loads
the raw Xenium output with ``spatialdata_io.xenium(data_path)``, and CI has no
Xenium dataset. This test substitutes *only that node* for one reading an h5ad
written to ``tmp_path``. Every other cell is byte-identical to what the viewer
recorded — asserted below by
``test_the_notebook_cells_are_the_recorded_sources_verbatim``. This is the same
preamble exception CLAUDE.md already documents; ``scripts/verify_notebook.py``
exercises the real preamble against a real dataset.

Run standalone:   python tests/test_notebook_replay.py
"""
from __future__ import annotations

import ast
import sys
import types
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
sq = pytest.importorskip("squidpy")
pytest.importorskip("anndata")
pytest.importorskip("nbformat")
pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")
adjusted_rand_score = pytest.importorskip("sklearn.metrics").adjusted_rand_score

from xenium_viewer.utils import notebook_export  # noqa: E402
from xenium_viewer.utils.environment import environment_code  # noqa: E402
from xenium_viewer.utils.prov_graph import SETUP  # noqa: E402
from xenium_viewer.utils.steps import Step, StepExecutor  # noqa: E402

# The real templates, imported as constants so the tab modules' Qt/napari
# imports are never triggered.
from xenium_viewer.tabs.tab_clustering import (  # noqa: E402
    FLAVOR_DEFAULTS, _leiden_template,
)
from xenium_viewer.tabs._helpers import (  # noqa: E402
    _NORMALIZE_TEMPLATE, _SPATIAL_NEIGHBORS_TEMPLATE,
)
from xenium_viewer.tabs.tab_gene_analysis import _RANK_GENES_TEMPLATE  # noqa: E402
from xenium_viewer.utils.gene_analysis import rank_genes_key  # noqa: E402
from xenium_viewer.tabs.tab_nhood import _NHOOD_TEMPLATE  # noqa: E402

CLUSTER_KEYS = ("leiden_igraph_r0.5", "leiden_igraph_r1.0")
PRIMARY_KEY = CLUSTER_KEYS[1]
N_PERMS = 100
TOP_N = 10

# Stands in for the viewer's ``preamble`` node — the one substitution (see the
# module docstring). Same shape as the real one: the imports every later cell
# relies on, then a single expression binding ``adata``.
_PREAMBLE_TEMPLATE = """
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from pathlib import Path

data_path = Path($data_path)
adata = sc.read_h5ad(data_path)"""

# Appended to the exported notebook, never recorded in the graph: the replayed
# kernel has to hand its results back across a process boundary somehow.
_DUMP_TEMPLATE = """
out = Path({out!r})
adata.obs[{keys!r}].to_csv(out / "replay_obs.csv")
rank_df.to_csv(out / "replay_rank.csv", index=False)
np.save(out / "replay_nhood.npy", nhood_zscore)"""


def _leiden_step(resolution: float, flavor: str = "igraph") -> Step:
    """Exactly what ``on_run_leiden`` builds, key naming included."""
    key = f"leiden_{flavor}_r{resolution}"
    n_iterations, directed = FLAVOR_DEFAULTS[flavor]
    return Step(
        id=f"clustering:{key}",
        template=_leiden_template(use_hvg=False, do_scale=False),
        params={
            "key": key, "resolution": resolution, "n_neighbors": 15, "n_pcs": 10,
            "n_top_genes": 30, "flavor": flavor, "n_iterations": n_iterations,
            "directed": directed, "random_state": 0,
        },
        deps=["normalize"],
        label=f"Clustering: {key}",
        outputs=["leiden_labels"],
    )


def _run_the_analysis(h5ad_path: Path) -> StepExecutor:
    """Run the real steps in-process, as the viewer's tab callbacks do.

    The namespace starts *empty*: every name a later step needs must come from
    the preamble's imports, which is the same guarantee the notebook needs.
    """
    ex = StepExecutor(namespace={})
    # Recorded, not executed in-process — exactly as the viewer records it. Its
    # point here is that a clean kernel can run it: it is the first cell of
    # every exported notebook, so if it raises, nothing downstream replays.
    ex.graph.upsert("environment", environment_code(), kind=SETUP,
                    label="Environment & seeds")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(Step(
            id="preamble", template=_PREAMBLE_TEMPLATE,
            params={"data_path": str(h5ad_path)}, kind=SETUP,
            label="Setup & data loading",
        ))
        ex.run(Step(
            id="normalize", template=_NORMALIZE_TEMPLATE, deps=["preamble"],
            kind=SETUP, label="Normalize, log-transform, PCA",
            outputs=["adata_norm"],
        ))
        for key in CLUSTER_KEYS:
            ex.run(_leiden_step(float(key.rsplit("r", 1)[1])))
        ex.run(Step(
            id=f"rank_genes:{PRIMARY_KEY}", template=_RANK_GENES_TEMPLATE,
            params={"groupby": PRIMARY_KEY, "method": "wilcoxon", "n_genes": 25,
                    "rank_key": rank_genes_key(PRIMARY_KEY)},
            deps=["normalize", f"clustering:{PRIMARY_KEY}"],
            label=f"Rank genes: {PRIMARY_KEY}", outputs=["rank_df"],
        ))
        ex.run(Step(
            id="spatial_neighbors", template=_SPATIAL_NEIGHBORS_TEMPLATE,
            params={"n_neighs": 6}, deps=["normalize"], label="Spatial neighbors",
        ))
        ex.run(Step(
            id=f"nhood:{PRIMARY_KEY}", template=_NHOOD_TEMPLATE,
            params={
                "cluster_key": PRIMARY_KEY,
                "uns_key": f"{PRIMARY_KEY}_nhood_enrichment",
                "n_perms": N_PERMS, "seed": 42,
            },
            deps=[f"clustering:{PRIMARY_KEY}", "spatial_neighbors"],
            label=f"Neighborhood enrichment: {PRIMARY_KEY}",
            outputs=["nhood_zscore"],
        ))
    return ex


class _Replay:
    """In-process results, the exported notebook, and the replayed results."""

    def __init__(self, executor, nb_path, out_dir):
        self.executor = executor
        self.graph = executor.graph
        self.nb_path = nb_path
        self.gui_obs = executor.get("adata").obs
        self.gui_rank = executor.get("rank_df")
        self.gui_nhood = executor.get("nhood_zscore")
        self.replay_obs = pd.read_csv(out_dir / "replay_obs.csv", index_col=0)
        self.replay_rank = pd.read_csv(out_dir / "replay_rank.csv")
        self.replay_nhood = np.load(out_dir / "replay_nhood.npy")


@pytest.fixture(scope="module")
def replay(tmp_path_factory, replay_adata):
    """Run the analysis, export the notebook, execute it, collect both sides.

    Module-scoped: starting a kernel and importing scanpy/squidpy in it costs
    far more than every assertion made against the result.
    """
    tmp_path = tmp_path_factory.mktemp("replay")
    h5ad_path = tmp_path / "input.h5ad"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        replay_adata().write_h5ad(h5ad_path)

    executor = _run_the_analysis(h5ad_path)

    nb_path = tmp_path / "analysis_notebook.ipynb"
    cells = notebook_export.graph_to_cells(executor.graph)
    cells.append(("code", _DUMP_TEMPLATE.format(
        out=str(tmp_path), keys=list(CLUSTER_KEYS),
    )))
    notebook_export.write_notebook(cells, nb_path)

    notebook_export.execute_notebook(nb_path, cwd=tmp_path, timeout=900)
    return _Replay(executor, nb_path, tmp_path)


# ── the claim ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("key", CLUSTER_KEYS)
def test_the_replayed_notebook_reproduces_the_clusterings(replay, key):
    """ARI exactly 1.0 — the partitions must be identical, not merely similar."""
    gui = replay.gui_obs[key].astype(str).to_numpy()
    replayed = replay.replay_obs[key].astype(str).to_numpy()

    assert len(set(gui)) >= 2, "the fixture must produce something to cluster"
    assert adjusted_rand_score(gui, replayed) == 1.0
    # Stronger than ARI, which is invariant to relabelling: the notebook must
    # also hand back the same cluster *names*, since downstream cells key on them.
    assert (gui == replayed).all()


def test_the_replayed_notebook_reproduces_the_ranked_genes(replay):
    """Identical top-N gene names per group, in the same order."""
    def top_n(frame):
        return {
            str(group): list(sub.head(TOP_N)["names"])
            for group, sub in frame.groupby("group", observed=True)
        }

    gui, replayed = top_n(replay.gui_rank), top_n(replay.replay_rank)
    assert set(gui) == set(replayed)
    assert gui == replayed


def test_the_replayed_notebook_reproduces_the_nhood_enrichment(replay):
    """Z-scores come from a seeded permutation test; they must land together."""
    assert replay.gui_nhood.shape == replay.replay_nhood.shape
    assert np.allclose(replay.gui_nhood, replay.replay_nhood, equal_nan=True)


# ── the properties that make the claim above mean something ──────────────────

def test_the_notebook_cells_are_the_recorded_sources_verbatim(replay):
    """No rewriting between the graph and the .ipynb.

    Without this, the replay could be passing because the exporter fixed
    something up on its way out — which the viewer's own graph would not carry.
    """
    code_cells = [
        source for kind, source in notebook_export.read_notebook(replay.nb_path)
        if kind == "code"
    ]
    injected = code_cells.pop()  # the results dump appended by this test
    assert "replay_obs.csv" in injected

    recorded = [
        replay.graph.get(nid).code.strip("\n")
        for nid in replay.graph.topo_sort()
    ]
    assert code_cells == recorded


def test_every_exported_cell_is_executable_python(replay):
    """A comment-only cell replays as a silent no-op, not as a failure.

    ``allow_errors=False`` cannot catch a node that records prose instead of
    code — it just does nothing and the notebook still "passes". The steps
    exercised here are all real code; ``scripts/verify_notebook.py`` reports the
    comment-only nodes a *real* session still contains (Phase 0.3).
    """
    import ast

    for nid in replay.graph.topo_sort():
        code = replay.graph.get(nid).code
        tree = ast.parse(code)
        assert tree.body, f"node {nid!r} exports no executable statement"


def test_the_notebook_opens_by_declaring_its_environment(replay):
    """First cell, and it must run — a version pin that raises is worse than none.

    Recorded versions belong at the top of the notebook because that is where
    someone reading a disagreeing result will look first; ``print_header()``
    puts the replay's own versions right beside them in the output.
    """
    order = replay.graph.topo_sort()
    assert order[0] == "environment"

    code_cells = [
        source for kind, source in notebook_export.read_notebook(replay.nb_path)
        if kind == "code"
    ]
    assert "sc.logging.print_header()" in code_cells[0]
    assert "np.random.seed(" in code_cells[0]


def test_the_preamble_is_the_only_substituted_node(replay):
    """Guards the one documented divergence from staying alone.

    If a second node ever has to be swapped out to make CI pass, that is a real
    reproducibility gap and this test should be the thing that says so.
    """
    substituted = {"preamble"}
    for nid in replay.graph.topo_sort():
        if nid in substituted:
            continue
        code = replay.graph.get(nid).code
        assert "read_h5ad" not in code, f"node {nid!r} was rewritten for the test"


# ── scripts/verify_notebook.py: the pure half of the Tier-2 report ───────────

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
verify_notebook = pytest.importorskip("verify_notebook")


def _graph_with(nodes):
    from xenium_viewer.utils.prov_graph import ProvGraph
    graph = ProvGraph()
    for node_id, code in nodes:
        graph.upsert(node_id, code, kind=SETUP)
    return graph


def test_comment_only_nodes_are_the_ones_that_replay_as_no_ops():
    """The Phase-0.3 punch list, derived by parsing rather than by eyeballing.

    ``allow_errors=False`` cannot see these: a cell of comments executes fine
    and does nothing, so the report has to name them explicitly.
    """
    graph = _graph_with([
        ("real", "adata_norm = adata.copy()"),
        ("he:flip", "\n# H&E image flip: flip_vertical=True, flip_horizontal=False"),
        ("viewer:zoom", "# zoomed to (100, 200)"),
        ("blank", "\n\n"),
    ])
    assert verify_notebook.comment_only_nodes(graph) == ["blank", "he:flip", "viewer:zoom"]


def test_declared_viewer_state_is_counted_apart_from_the_punch_list():
    """A note is not a gap, and a gap must not hide among the notes.

    The canvas background has no notebook equivalent; a clustering that
    recorded prose does. Both used to land in the same list, so the list was
    mostly display state and the defects in it were invisible.
    """
    from xenium_viewer.utils.prov_graph import NOTE

    graph = _graph_with([
        ("real", "adata_norm = adata.copy()"),
        ("clustering:k", "# Leiden clustering (shown in the viewer)"),
    ])
    graph.upsert("viewer:background", "# background set to white",
                 kind=NOTE)

    assert verify_notebook.comment_only_nodes(graph) == ["clustering:k"]
    assert verify_notebook.note_nodes(graph) == ["viewer:background"]


def test_a_failing_replay_can_be_traced_back_to_the_node_that_broke(tmp_path):
    """nbclient reports a cell index; the report has to name the recorded step.

    Markdown label cells shift the numbering, so the mapping is built from the
    derived cells rather than assumed to be the topo order.
    """
    graph = _graph_with([("preamble", "x = 1")])
    graph.upsert("clustering:k", "y = 2", deps=["preamble"], label="Clustering: k")
    _nb_path, node_ids = verify_notebook.build_notebook(graph, tmp_path)

    cells = notebook_export.read_notebook(_nb_path)
    assert len(node_ids) == len(cells)
    # the markdown header and the code cell both map back to the same node
    assert node_ids[cells.index(("markdown", "## Clustering: k"))] == "clustering:k"
    assert node_ids[-1] is None, "the injected dump cell belongs to no node"


def test_a_failing_cell_is_reported_against_its_node_and_not_counted_as_run(tmp_path):
    """Regression: nbclient calls ``on_cell_executed`` for the failing cell too.

    Taking that hook as "succeeded" named no node at all and counted the broken
    cell among the ones that ran — which is how the first real run of this
    script reported a `FileNotFoundError` with ``failed_node: None``.
    """
    graph = _graph_with([("preamble", "x = 1")])
    graph.upsert("clustering:missing", 'open("no-such-file-xv-test")', deps=["preamble"])
    nb_path, node_ids = verify_notebook.build_notebook(graph, tmp_path)

    timings, cursor = [], {}
    with pytest.raises(Exception):
        verify_notebook.execute(nb_path, tmp_path, 120, timings, node_ids, cursor)

    assert cursor["node"] == "clustering:missing"
    assert [c for c in timings if c.get("failed")] == [
        c for c in timings if c["node"] == "clustering:missing"
    ]
    assert [c["node"] for c in timings if not c.get("failed")] == ["preamble"]


def test_the_report_flags_a_clustering_the_notebook_never_produced():
    """Silence is the failure mode that matters — a viewer result with no cell
    behind it looks like a pass unless it is named."""
    viewer_obs = pd.DataFrame({"clustering_leiden_r1.0": ["0", "1", "0", "1"]})
    replay_obs = pd.DataFrame({"something_else": [1, 2, 3, 4]})
    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)
    assert entry["status"] == "not_in_replay"


def test_the_report_scores_a_matching_clustering_as_ari_one():
    labels = ["a", "a", "b", "b", "c"]
    viewer_obs = pd.DataFrame({"clustering_leiden_r1.0": labels})
    replay_obs = pd.DataFrame({"leiden_r1.0": labels})
    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)
    assert entry["status"] == "ok"
    assert entry["ari"] == 1.0
    assert entry["identical_labels"]
    assert entry["n_clusters_viewer"] == 3


def test_a_relabelled_clustering_scores_ari_one_but_not_identical():
    """ARI is invariant to relabelling; downstream cells key on the names.

    So the report carries both, and a run that is ARI-1.0 with different names
    is still visibly not the same result.
    """
    viewer_obs = pd.DataFrame({"clustering_k": ["0", "0", "1", "1"]})
    replay_obs = pd.DataFrame({"k": ["1", "1", "0", "0"]})
    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)
    assert entry["ari"] == 1.0
    assert not entry["identical_labels"]


def test_clusterings_are_compared_on_cell_identity_not_row_order():
    viewer_obs = pd.DataFrame(
        {"cell_id": ["c", "a", "b"], "clustering_k": ["2", "0", "1"]},
    )
    replay_obs = pd.DataFrame({"cell_id": ["a", "b", "c"], "k": ["0", "1", "2"]})
    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)
    assert entry["identical_labels"]
    assert entry["n_cells_shared"] == 3


def test_top_n_gene_agreement_is_order_sensitive():
    viewer_names = {"0": ["A", "B", "C"], "1": ["D", "E", "F"]}
    same = pd.DataFrame({
        "group": ["0", "0", "0", "1", "1", "1"],
        "names": ["A", "B", "C", "D", "E", "F"],
    })
    swapped = same.copy()
    swapped.loc[0, "names"], swapped.loc[1, "names"] = "B", "A"

    assert verify_notebook.compare_rank_genes(viewer_names, same, 3)["status"] == "ok"
    reordered = verify_notebook.compare_rank_genes(viewer_names, swapped, 3)
    assert reordered["status"] == "diverged"
    # same genes, different order — the report says so rather than just failing
    assert reordered["groups"]["0"]["n_shared"] == 3


# ── the two ways the comparison itself was wrong ─────────────────────────────

def test_cells_without_a_label_are_excluded_rather_than_crashing():
    """A CNV clustering is computed on a subset and reindexed onto the table.

    The mask tested the *rendered* label against "nan"; under pandas 3 a null in
    a categorical renders as `<NA>`, so it passed straight through and sklearn
    died with "Input contains NaN" — after a ten-minute replay, at the last step.
    """
    viewer_obs = pd.DataFrame({
        "clustering_cnv": pd.Categorical(["0", "1", None, "0", None]),
    })
    replay_obs = pd.DataFrame({"cnv": pd.Categorical(["0", "1", None, "0", None])})

    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)

    assert entry["status"] == "ok"
    assert entry["ari"] == 1.0
    assert entry["n_cells_compared"] == 3        # the labelled ones
    assert entry["n_cells_shared"] == 5
    assert entry["n_cells_unlabelled_viewer"] == 2


def test_a_clustering_labelled_on_neither_side_is_reported_not_scored():
    viewer_obs = pd.DataFrame({"clustering_k": [None, None]})
    replay_obs = pd.DataFrame({"k": [None, None]})

    (entry,) = verify_notebook.compare_clusterings(viewer_obs, replay_obs)
    assert entry["status"] == "empty"


def test_ranked_genes_are_compared_against_the_clustering_they_belong_to():
    """The viewer stores one `uns['rank_genes_groups']`; scanpy overwrites it.

    A session that ranked two clusterings stores one and the notebook binds the
    other last, so comparing them reported every group as diverged when nothing
    had: right genes, group numbering from a different clustering. The replayed
    frame now carries a `groupby` column and the comparison selects on it.
    """
    viewer_names = {"0": ["A", "B"], "1": ["C", "D"]}
    replay = pd.DataFrame({
        "group": ["0", "0", "1", "1", "0", "0", "1", "1"],
        "names": ["A", "B", "C", "D", "Z", "Y", "X", "W"],
        "groupby": ["leiden_r1.0"] * 4 + ["leiden_r0.5"] * 4,
    })

    matched = verify_notebook.compare_rank_genes(viewer_names, replay, 2,
                                                 "leiden_r1.0")
    assert matched["status"] == "ok"

    other = verify_notebook.compare_rank_genes(viewer_names, replay, 2,
                                               "leiden_r0.5")
    assert other["status"] == "diverged"


def test_a_clustering_the_notebook_never_ranked_is_named_as_such():
    """Not "diverged" — nothing diverged; the comparison had nothing to compare."""
    replay = pd.DataFrame({
        "group": ["0"], "names": ["A"], "groupby": ["leiden_r0.5"],
    })
    result = verify_notebook.compare_rank_genes({"0": ["A"]}, replay, 1,
                                                "leiden_r1.0")
    assert result["status"] == "different_groupby"
    assert result["replay_groupbys"] == ["leiden_r0.5"]


def test_every_ranking_survives_into_the_notebook():
    """`rank_df` is rebound by the next ranking; the keyed `uns` slots are not.

    Without this the notebook ends holding markers for whichever clustering was
    ranked last, and the others are simply gone from the replayed state. This
    used to be a notebook-local `rank_results` dict, which kept the frames but
    put them somewhere no `sc.pl.rank_genes_groups*` call could reach — they all
    take `key=`. `key_added=` is the scanpy API for exactly this.
    """
    assert "key_added=$rank_key" in _RANK_GENES_TEMPLATE
    assert "key=$rank_key" in _RANK_GENES_TEMPLATE
    assert "globals()" not in _RANK_GENES_TEMPLATE

    # The harvest cell in verify_notebook enumerates those slots and tags each
    # frame with scanpy's own recorded groupby, so the two sides of
    # `compare_rank_genes` can still be matched up by clustering name.
    assert '"params"]["groupby"]' in verify_notebook._DUMP_CELL
    assert "rank_results" not in verify_notebook._DUMP_CELL
    _harvest_the_rankings(
        {rank_genes_key("a"): {"params": {"groupby": "a"}},
         rank_genes_key("b"): {"params": {"groupby": "b"}}},
        expected=[rank_genes_key("a"), rank_genes_key("b")])


def test_the_harvest_skips_uns_entries_that_only_look_like_rankings():
    """`rank_genes_groupby` shares the prefix and is a plain string.

    A restored session also carries a stale unkeyed `rank_genes_groups`
    alongside the keyed ones. Selecting on the prefix alone put both into
    `sc.get.rank_genes_groups_df(key=...)`, which is a crash in the *injected*
    cell — i.e. a verification run that fails while the analysis is fine.
    """
    _harvest_the_rankings(
        {"rank_genes_groupby": "leiden_r1.0",
         "rank_genes_groups": {"names": ["A"]},          # no params
         rank_genes_key("a"): {"params": {"groupby": "a"}}},
        expected=[rank_genes_key("a")])


def _harvest_the_rankings(uns: dict, expected: list[str]):
    """Run the `_rank_keys = ...` statement of the injected cell over *uns*."""
    source = verify_notebook._DUMP_CELL.format(out="/tmp")
    stmt = next(s for s in ast.parse(source).body
                if isinstance(s, ast.Assign)
                and getattr(s.targets[0], "id", None) == "_rank_keys")
    namespace = {"adata_norm": types.SimpleNamespace(uns=uns)}
    exec(compile(ast.Module([stmt], []), "<harvest>", "exec"), namespace)  # noqa: S102
    assert sorted(namespace["_rank_keys"]) == sorted(expected)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
