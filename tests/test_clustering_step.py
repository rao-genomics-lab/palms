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

from palms.utils.steps import Step, StepExecutor, check_step  # noqa: E402

# Imported from the tab module's constants without triggering its Qt imports.
from palms.tabs.tab_clustering import (  # noqa: E402
    FLAVOR_DEFAULTS, LEIDEN_FLAVORS, _leiden_template,
)
from palms.tabs._helpers import _NORMALIZE_TEMPLATE  # noqa: E402


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


def _step(key=None, use_hvg=False, do_scale=False, n_pcs=10, flavor="igraph",
          resolution=1.0):
    """Build the step exactly as ``on_run_leiden`` does, key naming included."""
    key = key if key is not None else f"leiden_{flavor}_r{resolution}"
    n_iterations, directed = FLAVOR_DEFAULTS[flavor]
    return Step(
        id=f"clustering:{key}",
        template=_leiden_template(use_hvg, do_scale),
        params={
            "key": key, "resolution": resolution, "n_neighbors": 15, "n_pcs": n_pcs,
            "n_top_genes": 30, "flavor": flavor, "n_iterations": n_iterations,
            "directed": directed, "random_state": 0,
        },
        deps=["normalize"],
        label=f"Clustering: {key}",
        outputs=["leiden_labels"],
    )


def _run(step, adata):
    """Run the real normalize -> clustering pair, as the viewer does."""
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        ex.run(step)
    return ex


def test_the_step_hands_back_the_labels_it_declared():
    """The tab reads labels from the returned outputs, not from ctx.adata.obs.

    Reading them back off ``adata`` worked only while the executor namespace and
    ``ctx.adata`` were the same object — an invariant maintained by hand in
    ``_run_step``, and invisible to anyone editing the template. Declaring the
    output moves the check into ``StepExecutor``, which raises if the template
    stops binding the name.
    """
    adata = _adata()
    step = _step()
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        outputs = ex.run(step)

    labels = outputs["leiden_labels"]
    assert list(labels) == list(adata.obs[step.params["key"]])
    assert labels.nunique() >= 2


def test_a_template_that_stops_binding_the_labels_fails_loudly():
    """The point of the declared output: silent breakage becomes a StepError."""
    from palms.utils.steps import StepError

    adata = _adata()
    step = _step()
    # A user edit that drops the last line — the analysis still "works", but the
    # tab would previously have read a stale or missing obs column.
    step.template = step.template.rsplit("\n", 1)[0]
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.run(_normalize_step())
        with pytest.raises(StepError, match="leiden_labels"):
            ex.run(step)


@pytest.mark.parametrize("flavor", LEIDEN_FLAVORS)
@pytest.mark.parametrize("use_hvg,do_scale", [
    (False, False), (True, False), (False, True), (True, True),
])
def test_template_executes_and_labels_the_cells(use_hvg, do_scale, flavor):
    if flavor == "leidenalg":
        pytest.importorskip("leidenalg")
    adata = _adata()
    step = _step(use_hvg=use_hvg, do_scale=do_scale, flavor=flavor)
    ex = _run(step, adata)
    key = step.params["key"]

    assert key in adata.obs
    assert adata.obs[key].nunique() >= 2
    # recorded source == executed source
    assert ex.graph.get(f"clustering:{key}").code == step.render()


@pytest.mark.parametrize("flavor", LEIDEN_FLAVORS)
def test_recorded_source_replays_to_the_same_labels(flavor):
    """Re-running the recorded cell in a clean namespace reproduces the labels.

    This is the reproducibility claim in miniature: no viewer objects, no
    session state — just the source the notebook would contain. Both backends
    are seeded from ``random_state`` (leidenalg via ``seed=``, igraph via
    ``set_igraph_random_state``), so the claim has to hold for either.
    """
    if flavor == "leidenalg":
        pytest.importorskip("leidenalg")
    adata = _adata()
    step = _step(flavor=flavor)
    key = step.params["key"]
    ex = _run(step, adata)
    gui_labels = adata.obs[key].astype(str).to_numpy()

    # Replay every recorded cell in dependency order, as the notebook does.
    replay_ns = {"sc": sc, "adata": _adata()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for node_id in ex.graph.topo_sort():
            exec(compile(ex.graph.get(node_id).code, "<replay>", "exec"), replay_ns)  # noqa: S102
    replay_labels = replay_ns["adata"].obs[key].astype(str).to_numpy()

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
    step = _step()
    ex = _run(step, adata)
    order = ex.graph.topo_sort()
    assert order.index("normalize") < order.index(step.id)


@pytest.mark.parametrize("flavor", LEIDEN_FLAVORS)
def test_leiden_call_pins_flavor_and_its_backend_specific_arguments(flavor):
    """scanpy's default flavor is scheduled to change, and the two backends
    disagree on ``n_iterations`` (2 vs -1) and ``directed`` (False vs True).

    Leaving any of them implicit would silently change clusterings on a scanpy
    upgrade, so the recorded source pins all four literally.
    """
    n_iterations, directed = FLAVOR_DEFAULTS[flavor]
    source = _step(flavor=flavor).render()
    assert f"flavor={flavor!r}" in source
    assert f"n_iterations={n_iterations}" in source
    assert f"directed={directed}" in source
    assert "random_state=0" in source


def test_the_flavor_is_part_of_the_key():
    """Both backends at one resolution must coexist, not overwrite each other.

    They produce genuinely different partitions, so sharing a key would both
    lose a result and revise the DAG node in place, flagging its descendants
    stale for a clustering the user still wanted.
    """
    igraph_step = _step(flavor="igraph")
    leidenalg_step = _step(flavor="leidenalg")
    assert igraph_step.params["key"] == "leiden_igraph_r1.0"
    assert leidenalg_step.params["key"] == "leiden_leidenalg_r1.0"
    assert igraph_step.id != leidenalg_step.id


def test_params_appear_literally_in_the_recorded_source():
    source = _step(resolution=0.5).render()
    assert "resolution=0.5" in source
    assert "key_added='leiden_igraph_r0.5'" in source
    assert "n_neighbors=15" in source
    assert "n_pcs=10" in source


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
