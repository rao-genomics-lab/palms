"""The normalize and rank-genes steps: recorded source == executed source.

These pin the second known divergence between the GUI and the notebook: the
viewer normalised with ``target_sum=1e4`` (``gene_analysis.get_normalized_adata``)
while the recorded cell used ``sc.pp.normalize_total(adata)``, i.e. scanpy's
median default — different X, so different PCA, neighbours, clusters and DEG.

Run standalone:   python tests/test_normalize_rank_genes_steps.py
Or with pytest:   pytest tests/test_normalize_rank_genes_steps.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")

from xenium_viewer.utils.steps import Step, StepExecutor, check_step  # noqa: E402
from xenium_viewer.tabs._helpers import _NORMALIZE_TEMPLATE  # noqa: E402
from xenium_viewer.tabs.tab_gene_analysis import _RANK_GENES_TEMPLATE  # noqa: E402
from xenium_viewer.utils.gene_analysis import rank_genes_key  # noqa: E402


def _adata(n_obs: int = 120, n_vars: int = 40):
    rng = np.random.default_rng(1)
    counts = rng.poisson(4, size=(n_obs, n_vars)).astype("float32")
    counts[: n_obs // 2, : n_vars // 4] += 15
    a = anndata.AnnData(counts)
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_vars)]
    a.obs["group"] = pd.Categorical(
        ["a"] * (n_obs // 2) + ["b"] * (n_obs - n_obs // 2)
    )
    return a


def _normalize_step():
    return Step(id="normalize", template=_NORMALIZE_TEMPLATE,
                label="Normalize, log-transform, PCA", outputs=["adata_norm"])


def _rank_step(groupby="group", method="wilcoxon", n_genes=10):
    return Step(
        id=f"rank_genes:{groupby}",
        template=_RANK_GENES_TEMPLATE,
        params={"groupby": groupby, "method": method, "n_genes": n_genes,
                # Derived from the clustering by the tab, exactly as here.
                "rank_key": rank_genes_key(groupby)},
        outputs=["rank_df", "adata_norm"],
    )


def _run_both(adata):
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        ex.run(_rank_step())
    return ex


# ── the divergence this fixes ────────────────────────────────────────────────

def test_normalize_records_the_target_sum_it_uses():
    """The old recorded cell omitted target_sum and silently median-scaled."""
    assert "target_sum=1e4" in _NORMALIZE_TEMPLATE


def test_normalize_does_not_mutate_adata():
    """It binds a copy. The old SETUP node mutated adata in place, so any step
    that copied adata could normalise twice."""
    adata = _adata()
    before = adata.X.copy()
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
    assert np.array_equal(adata.X, before)
    assert not np.array_equal(ex.ns["adata_norm"].X, before)


def test_normalize_matches_get_normalized_adata():
    """The step must reproduce what the viewer has always actually computed."""
    from xenium_viewer.utils.gene_analysis import get_normalized_adata

    adata = _adata()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy = get_normalized_adata(adata)
        ex = StepExecutor(namespace={"sc": sc, "adata": adata})
        ex.run(_normalize_step())
    np.testing.assert_allclose(
        np.asarray(ex.ns["adata_norm"].X), np.asarray(legacy.X), rtol=1e-6
    )


def test_rank_genes_runs_on_the_normalized_copy_not_raw_counts():
    """The old recorded cell ranked on `adata`, which was raw unless a SETUP
    normalize node happened to be in the graph."""
    assert "sc.tl.rank_genes_groups(" in _RANK_GENES_TEMPLATE
    assert "adata_norm, groupby=" in _RANK_GENES_TEMPLATE
    assert "rank_genes_groups(adata," not in _RANK_GENES_TEMPLATE.replace(" ", "")


def test_a_second_ranking_does_not_overwrite_the_first():
    """What the notebook-local ``rank_results`` dict used to work around.

    scanpy overwrites ``uns['rank_genes_groups']`` in place, so a session that
    ranked two clusterings ended holding markers for whichever ran last. The
    dict kept both but no ``sc.pl.rank_genes_groups*`` call could reach them —
    they all take ``key=``. ``key_added=`` gives both properties at once.
    """
    adata = _adata()
    adata.obs["other"] = pd.Categorical(["x", "y"] * (adata.n_obs // 2))

    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        ex.run(_rank_step(groupby="group"))
        first = ex.ns["rank_df"]
        ex.run(_rank_step(groupby="other"))

    uns = ex.ns["adata_norm"].uns
    assert rank_genes_key("group") in uns and rank_genes_key("other") in uns
    # Nothing landed in scanpy's default slot, so nothing could be clobbered.
    assert "rank_genes_groups" not in uns

    recovered = sc.get.rank_genes_groups_df(ex.ns["adata_norm"], group=None,
                                            key=rank_genes_key("group"))
    pd.testing.assert_frame_equal(first, recovered)


# ── the guarantee ────────────────────────────────────────────────────────────

def test_recorded_source_is_executed_source():
    adata = _adata()
    ex = _run_both(adata)
    assert ex.graph.get("normalize").code == _normalize_step().render()
    assert ex.graph.get("rank_genes:group").code == _rank_step().render()


def test_recorded_cells_replay_to_the_same_results():
    adata = _adata()
    ex = _run_both(adata)
    gui_df = ex.ns["rank_df"]

    replay_ns = {"sc": sc, "adata": _adata()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for node_id in ex.graph.topo_sort():
            exec(compile(ex.graph.get(node_id).code, "<replay>", "exec"), replay_ns)  # noqa: S102

    pd.testing.assert_frame_equal(gui_df, replay_ns["rank_df"])


def test_steps_are_self_contained():
    assert check_step(_normalize_step(), available={"sc", "adata"}) == set()
    assert check_step(_rank_step(), available={"sc", "adata", "adata_norm"}) == set()


def test_rank_genes_params_appear_literally():
    source = _rank_step(groupby="leiden_r1.0", method="t-test", n_genes=25).render()
    assert "groupby='leiden_r1.0'" in source
    assert "method='t-test'" in source
    assert "n_genes=25" in source
    assert "key_added='rank_genes_leiden_r1.0'" in source


def test_topo_order_puts_normalize_before_rank_genes():
    """rank_genes declares normalize as a dependency, so adata_norm exists."""
    adata = _adata()
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        rank = _rank_step()
        rank.deps = ["normalize"]
        ex.run(rank)
    order = ex.graph.topo_sort()
    assert order.index("normalize") < order.index("rank_genes:group")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
