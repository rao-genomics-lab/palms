"""The Leiden step: the source the viewer runs is the source the notebook records.

These tests execute the *real* template from ``tabs.tab_clustering`` against a
small synthetic AnnData, then re-run the recorded source in a fresh namespace
and check the labels match. That is the property the paper claims, reduced to
something CI can assert.

Importing the tab module pulls in Qt/napari, so the template constants are
imported directly rather than through ``build_tab``.

Run standalone:   python tests/test_clustering_step.py
Or with pytest:   pytest tests/test_clustering_step.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
sc = pytest.importorskip("scanpy")

from xenium_viewer.utils.steps import Step, StepExecutor, check_step  # noqa: E402

# Imported from the tab module's constants without triggering its Qt imports.
from xenium_viewer.tabs.tab_clustering import _leiden_template  # noqa: E402
from xenium_viewer.tabs._helpers import _NORMALIZE_TEMPLATE  # noqa: E402


def _adata(n_obs: int = 200, n_vars: int = 60):
    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(n_obs, n_vars)).astype("float32")
    # two crude populations so clustering has something to find
    counts[: n_obs // 2, : n_vars // 3] += 12
    a = anndata.AnnData(counts)
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_vars)]
    return a


def _normalize_step():
    return Step(id="normalize", template=_NORMALIZE_TEMPLATE, kind="setup",
                label="Normalize, log-transform, PCA", outputs=["adata_norm"])


def _step(key="leiden_r1.0", use_hvg=False, do_scale=False, n_pcs=10):
    return Step(
        id=f"clustering:{key}",
        template=_leiden_template(use_hvg, do_scale),
        params={
            "key": key, "resolution": 1.0, "n_neighbors": 15, "n_pcs": n_pcs,
            "n_top_genes": 30, "random_state": 0,
        },
        deps=["normalize"],
        label=f"Clustering: {key}",
    )


def _run(step, adata):
    """Run the real normalize -> clustering pair, as the viewer does."""
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        ex.run(step)
    return ex


@pytest.mark.parametrize("use_hvg,do_scale", [
    (False, False), (True, False), (False, True), (True, True),
])
def test_template_executes_and_labels_the_cells(use_hvg, do_scale):
    adata = _adata()
    step = _step(use_hvg=use_hvg, do_scale=do_scale)
    ex = _run(step, adata)

    assert "leiden_r1.0" in adata.obs
    assert adata.obs["leiden_r1.0"].nunique() >= 2
    # recorded source == executed source
    assert ex.graph.get("clustering:leiden_r1.0").code == step.render()


def test_recorded_source_replays_to_the_same_labels():
    """Re-running the recorded cell in a clean namespace reproduces the labels.

    This is the reproducibility claim in miniature: no viewer objects, no
    session state — just the source the notebook would contain.
    """
    adata = _adata()
    step = _step()
    ex = _run(step, adata)
    gui_labels = adata.obs["leiden_r1.0"].astype(str).to_numpy()

    # Replay every recorded cell in dependency order, as the notebook does.
    replay_ns = {"sc": sc, "adata": _adata()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for node_id in ex.graph.topo_sort():
            exec(compile(ex.graph.get(node_id).code, "<replay>", "exec"), replay_ns)  # noqa: S102
    replay_labels = replay_ns["adata"].obs["leiden_r1.0"].astype(str).to_numpy()

    assert (gui_labels == replay_labels).all()


def test_step_is_self_contained_given_its_declared_inputs():
    """The template must not reach for anything the notebook won't have bound."""
    available = {"sc", "adata", "adata_norm"}
    assert check_step(_step(), available=available) == set()
    assert check_step(_step(use_hvg=True, do_scale=True), available=available) == set()


def test_clustering_starts_from_the_shared_normalized_copy():
    """Regression, twice over.

    Originally the step copied ``adata`` and relied on a SETUP ``normalize``
    node having mutated it in place — invisible in the DAG, and wrong whenever
    that node was absent. The over-correction then inlined normalisation into
    the step, which normalised twice in the exported notebook and hid the
    ``normalize -> clustering`` edge a reader expects.

    Now ``normalize`` binds a *copy* (``adata_norm``), so the step can consume
    it: one normalisation, one real edge.
    """
    for use_hvg in (False, True):
        source = _leiden_template(use_hvg=use_hvg, do_scale=False)
        assert "adata_leiden = adata_norm.copy()" in source
        # it must not re-normalise what `normalize` already normalised
        assert "normalize_total" not in source
        assert "log1p" not in source


def test_clustering_declares_normalize_as_a_dependency():
    """The edge must be in the graph, not just implied by the source."""
    assert _step().deps == ["normalize"]


def test_pca_is_recomputed_only_when_the_gene_set_or_scaling_changed():
    """Otherwise ``adata_norm``'s X_pca is exactly what we'd recompute."""
    assert "sc.pp.pca(adata_leiden)" not in _leiden_template(False, False)
    for use_hvg, do_scale in [(True, False), (False, True), (True, True)]:
        assert "sc.pp.pca(adata_leiden)" in _leiden_template(use_hvg, do_scale)


def test_normalize_sorts_before_clustering():
    adata = _adata()
    ex = _run(_step(), adata)
    order = ex.graph.topo_sort()
    assert order.index("normalize") < order.index("clustering:leiden_r1.0")


def test_leiden_call_pins_flavor_and_iterations():
    """scanpy warns that the default flavor will change to igraph.

    Leaving it implicit would silently change clusterings on a scanpy upgrade,
    so the recorded source pins flavor, n_iterations and random_state.
    """
    source = _step().render()
    assert "flavor='igraph'" in source
    assert "n_iterations=2" in source
    assert "random_state=0" in source


def test_params_appear_literally_in_the_recorded_source():
    source = _step(key="leiden_r0.5").render()
    assert "resolution=1.0" in source
    assert "key_added='leiden_r0.5'" in source
    assert "n_neighbors=15" in source
    assert "n_pcs=10" in source


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
