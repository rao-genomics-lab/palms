"""The tab's keep/drop readout is a second implementation of the filter.

It has to be, because the readout answers on every spin-box tick and the step
copies a whole AnnData. So the numbers come out of two sorted vectors and a
``searchsorted``, while the filter comes out of scanpy — and a second
implementation needs a gate pinning it to the first, in the shape of
``test_the_feather_preview_bins_identically_to_the_recorded_step``.

The gene half is the one worth watching: it counts detections over *the cells
that survive the cell cutoff*, without slicing the matrix, so it is easy for it
to drift into counting over all cells instead.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")
sc = pytest.importorskip("scanpy")


def _adata(sparse_x=True, n=90, g=25, seed=3):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(0.7, (n, g)).astype("float32")
    X = sparse.csr_matrix(counts) if sparse_x else counts
    return anndata.AnnData(
        X=X,
        obs=pd.DataFrame({"cell_id": np.arange(1, n + 1)},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(g)]),
    )


def _readout_numbers(adata, min_counts, min_cells):
    """What the tab would show, using the tab's own helpers."""
    from palms.tabs.tab_qc import cells_per_gene, counts_per_cell

    counts = np.sort(counts_per_cell(adata))
    n_cells = (len(counts) - int(np.searchsorted(counts, min_counts))
               if min_counts is not None else len(counts))
    mask = (counts_per_cell(adata) >= min_counts if min_counts is not None
            else np.ones(adata.n_obs, dtype=bool))
    per_gene = np.sort(cells_per_gene(adata, mask))
    n_genes = (len(per_gene) - int(np.searchsorted(per_gene, min_cells))
               if min_cells is not None else adata.n_vars)
    return n_cells, n_genes


def _step_numbers(adata, min_counts, min_cells):
    """What the recorded step actually produces."""
    from palms.tabs._helpers import qc_filter_preview
    from palms.utils.step_templates import builtin_assemble
    from palms.utils.steps import Step, StepExecutor

    blocks, params, _ = qc_filter_preview(min_counts, min_cells)
    executor = StepExecutor(namespace={"sc": sc, "adata": adata})
    result = executor.run(Step(
        id="qc_filter", template=builtin_assemble("qc.filter", blocks),
        params=params, outputs=["adata"]))["adata"]
    return result.n_obs, result.n_vars


CUTOFFS = list(itertools.product([None, 1, 8, 14, 25], [None, 1, 20, 45, 200]))


@pytest.mark.parametrize("sparse_x", [True, False], ids=["csr", "dense"])
@pytest.mark.parametrize("min_counts,min_cells", CUTOFFS)
def test_the_readout_matches_the_step_it_predicts(sparse_x, min_counts, min_cells):
    if min_counts is None and min_cells is None:
        pytest.skip("no assembly filters nothing; that case is clear_qc_filter")
    adata = _adata(sparse_x=sparse_x)
    assert _readout_numbers(adata, min_counts, min_cells) == \
        _step_numbers(adata, min_counts, min_cells)


def test_genes_are_counted_over_the_surviving_cells_only():
    """The specific drift this pins: counting detections over all cells.

    Constructed so the two answers differ — a gene detected only in cells the
    count cutoff removes must not survive the gene cutoff.
    """
    from palms.tabs.tab_qc import cells_per_gene

    X = np.zeros((6, 2), dtype="float32")
    X[:2, 0] = 50          # a rich pair, detecting Gene0
    X[2:, 1] = 1           # four near-empty cells, detecting Gene1
    adata = anndata.AnnData(sparse.csr_matrix(X),
                            var=pd.DataFrame(index=["Gene0", "Gene1"]))

    over_all = cells_per_gene(adata)
    rich_only = cells_per_gene(adata, np.asarray(X.sum(1)).ravel() >= 10)

    assert list(over_all) == [2, 4]
    assert list(rich_only) == [2, 0], (
        "Gene1 is detected only in cells the count cutoff drops, so after that "
        "cutoff it is detected in none"
    )


def test_counts_per_cell_matches_scanpys_own():
    from palms.tabs.tab_qc import counts_per_cell

    adata = _adata()
    _, per_cell = sc.pp.filter_cells(adata, min_counts=1, inplace=False)
    np.testing.assert_allclose(np.sort(counts_per_cell(adata)), np.sort(per_cell))


def test_explicit_zeros_in_a_sparse_matrix_do_not_count_as_detections():
    """``csr.indices`` lists stored entries, which need not be non-zero."""
    from palms.tabs.tab_qc import cells_per_gene

    X = sparse.csr_matrix(np.array([[1.0, 0.0], [0.0, 0.0]], dtype="float32"))
    X[1, 1] = 0.0            # store an explicit zero
    X = X.tocsr()
    adata = anndata.AnnData(X, var=pd.DataFrame(index=["Gene0", "Gene1"]))
    assert list(cells_per_gene(adata)) == [1, 0]
