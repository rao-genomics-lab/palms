"""The QC filter template: what it produces, and what it must not touch.

The dangerous half of ``qc.filter`` is not the arithmetic, it is the object it
is handed. In the viewer ``adata`` **is** the table inside the open SpatialData
store, while in the exported notebook the preamble binds
``adata = sdata["table"].copy()``. So a template that filtered in place would be
harmless on replay and destructive in the GUI — a divergence no replay test can
see, because the notebook never reaches the failing case.

Hence the two properties pinned hardest below: the step *rebinds* rather than
mutating, and the object it was handed comes out unchanged.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")
sc = pytest.importorskip("scanpy")

TEMPLATE_ID = "qc.filter"


def _adata(n=80, g=16, seed=0):
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(rng.poisson(0.5, (n, g)).astype("float32"))
    obs = pd.DataFrame({"cell_id": np.arange(1, n + 1)},
                       index=[f"c{i}" for i in range(n)])
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(g)])
    return anndata.AnnData(X=X, obs=obs, var=var)


def _run(adata, blocks, params):
    from palms.utils.step_templates import builtin_assemble
    from palms.utils.steps import Step, StepExecutor

    executor = StepExecutor(namespace={"sc": sc, "adata": adata})
    step = Step(id="qc_filter", template=builtin_assemble(TEMPLATE_ID, blocks),
                params=params, outputs=["adata"])
    return executor.run(step)["adata"], executor


def _assemblies():
    from palms.utils.step_templates import builtin_spec
    return list(builtin_spec(TEMPLATE_ID).assemblies)


def _params_for(blocks):
    params = {}
    if "cells" in blocks:
        params["min_counts"] = 5
    if "genes" in blocks:
        params["min_cells"] = 20
    return params


@pytest.mark.parametrize("blocks", _assemblies(),
                         ids=lambda b: "+".join(b))
def test_every_assembly_filters_and_rebinds(blocks):
    adata = _adata()
    before = adata.shape
    result, executor = _run(adata, blocks, _params_for(blocks))

    assert result is not adata, "the step must rebind, not mutate"
    assert executor.ns["adata"] is result, (
        "the namespace must hold the new object — _run_step follows it onto "
        "ctx.adata from there"
    )
    assert adata.shape == before, (
        "the object the step was handed is the live sdata['table'] in the "
        "viewer; filtering it in place would shrink the store's own cells"
    )
    assert result.n_obs <= before[0] and result.n_vars <= before[1]


@pytest.mark.parametrize("blocks", _assemblies(),
                         ids=lambda b: "+".join(b))
def test_the_result_is_materialised_not_a_view(blocks):
    """A view would be converted mid-analysis by the next step's obs write.

    ``clustering.leiden`` ends with ``adata.obs[key] = ...``; against a view
    anndata quietly makes a copy at that point, so the object the viewer holds
    and the one the step wrote to stop being the same.
    """
    result, _ = _run(_adata(), blocks, _params_for(blocks))
    assert result.is_view is False


def test_the_two_filters_are_sequential_not_independent():
    """Genes are counted over the cells that survived, as scanpy's own order.

    Computing both masks against the unfiltered matrix would be a different
    analysis, and one no scanpy tutorial describes.
    """
    adata = _adata()
    by_hand = adata.copy()
    sc.pp.filter_cells(by_hand, min_counts=5)
    sc.pp.filter_genes(by_hand, min_cells=20)

    result, _ = _run(adata, ["cells", "genes", "bind"],
                     {"min_counts": 5, "min_cells": 20})

    assert list(result.obs_names) == list(by_hand.obs_names)
    assert list(result.var_names) == list(by_hand.var_names)


def test_a_permissive_cutoff_keeps_everything():
    adata = _adata()
    result, _ = _run(adata, ["cells", "genes", "bind"],
                     {"min_counts": 1, "min_cells": 1})
    assert result.shape == adata.shape


def test_the_filter_does_not_write_scanpy_qc_columns():
    """``inplace=False`` is what keeps ``obs['n_counts']`` off the table.

    The in-place form annotates as a side effect, onto a table that already
    carries Xenium's own count columns and gets persisted to the store.
    """
    result, _ = _run(_adata(), ["cells", "genes", "bind"],
                     {"min_counts": 5, "min_cells": 20})
    assert "n_counts" not in result.obs.columns
    assert "n_cells" not in result.var.columns


def test_the_shipped_text_uses_inplace_false_and_rebinds():
    """Source guard: the shape above is the whole safety argument."""
    from palms.tabs.tab_qc import _qc_filter_template

    text = _qc_filter_template(filter_cells=True, filter_genes=True)
    assert "sc.pp.filter_cells(adata, min_counts=$min_counts, inplace=False)" in text
    assert "sc.pp.filter_genes(adata, min_cells=$min_cells, inplace=False)" in text
    assert "adata = adata.copy()" in text
    # Comments stripped: they explain why inplace=True is wrong here, and would
    # otherwise trip the check that no *call* uses it.
    code = "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))
    assert "inplace=True" not in code


def test_there_is_no_assembly_that_filters_nothing():
    """"Neither cutoff" is not a step; it is ``clear_qc_filter``."""
    from palms.utils.step_templates import builtin_spec

    for assembly in builtin_spec(TEMPLATE_ID).assemblies:
        assert "cells" in assembly or "genes" in assembly


def test_the_preview_helper_selects_declared_assemblies():
    from palms.tabs._helpers import qc_filter_preview
    from palms.utils.step_templates import builtin_spec

    declared = builtin_spec(TEMPLATE_ID).assemblies
    for min_counts, min_cells in ((10, 3), (10, None), (None, 3)):
        preview = qc_filter_preview(min_counts, min_cells)
        assert tuple(preview.blocks) in declared
        assert set(preview.params) == {
            name for name, value in (("min_counts", min_counts),
                                     ("min_cells", min_cells))
            if value is not None
        }
