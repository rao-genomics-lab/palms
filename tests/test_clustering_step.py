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


def _adata(n_obs: int = 200, n_vars: int = 60):
    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(n_obs, n_vars)).astype("float32")
    # two crude populations so clustering has something to find
    counts[: n_obs // 2, : n_vars // 3] += 12
    a = anndata.AnnData(counts)
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_vars)]
    return a


def _step(key="leiden_r1.0", use_hvg=False, do_scale=False, n_pcs=10):
    return Step(
        id=f"clustering:{key}",
        template=_leiden_template(use_hvg, do_scale),
        params={
            "key": key, "resolution": 1.0, "n_neighbors": 15, "n_pcs": n_pcs,
            "n_top_genes": 30, "random_state": 0,
        },
        label=f"Clustering: {key}",
    )


def _run(step, adata):
    ex = StepExecutor(namespace={"sc": sc, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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

    recorded = ex.graph.get("clustering:leiden_r1.0").code
    replay_ns = {"sc": sc, "adata": _adata()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exec(compile(recorded, "<replay>", "exec"), replay_ns)  # noqa: S102
    replay_labels = replay_ns["adata"].obs["leiden_r1.0"].astype(str).to_numpy()

    assert (gui_labels == replay_labels).all()


def test_step_is_self_contained_given_sc_and_adata():
    """The template must not reach for anything the notebook won't have bound."""
    assert check_step(_step(), available={"sc", "adata"}) == set()
    assert check_step(_step(use_hvg=True, do_scale=True),
                      available={"sc", "adata"}) == set()


def test_template_does_not_depend_on_a_shared_normalize_node():
    """Regression: the old HVG branch copied an externally-normalised adata.

    ``normalize`` is a SETUP node, so it sorts ahead of every artifact node and
    mutates ``adata`` in place. The step must normalise its own copy, so that a
    notebook containing both cells cannot double-normalise this one.
    """
    source = _leiden_template(use_hvg=True, do_scale=False)
    assert "adata_leiden = adata.copy()" in source
    assert "sc.pp.normalize_total(adata_leiden, target_sum=1e4)" in source
    # and it must not silently rely on adata already being log-normalised
    assert "sc.pp.log1p(adata_leiden)" in source


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
