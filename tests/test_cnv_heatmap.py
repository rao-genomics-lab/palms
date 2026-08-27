"""Saving the CNV heatmap must work however many windows the run produced.

``chromosome_heatmap(dendrogram=True)`` ends in ``sc.tl.dendrogram``, which
represents the cells with ``pd.DataFrame(_choose_representation(...))``. For a
matrix with more columns than ``settings.N_PCS`` that representation is a PCA —
dense, fine. For a *narrow* one it is ``.X`` itself, and
``pd.DataFrame(csr_matrix)`` does not densify: it builds a one-column object
frame of 1×n row matrices, so the ``.groupby().mean()`` that follows dies with
``TypeError: agg function failed [how->mean,dtype->object]``.

So the heatmap worked on a big panel and crashed on a small one — which is
exactly when a user is most likely to be looking, since few windows usually
means few genes mapped to the genome.

Run headless: ``MPLBACKEND=Agg pytest tests/test_cnv_heatmap.py``
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
pytest.importorskip("infercnvpy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from palms.utils.cnv_analysis import make_cnv_heatmap  # noqa: E402

GROUPBY = "cnv_leiden_res0.2"


def _adata_with_cnv(n_windows: int, n_obs: int = 60, sparse_cnv: bool = True):
    """An AnnData shaped like the CNV pipeline's output, with *n_windows*."""
    rng = np.random.default_rng(0)
    adata = anndata.AnnData(rng.random((n_obs, 4), dtype="float32"))
    adata.obs_names = [f"cell{i}" for i in range(n_obs)]
    adata.obs[GROUPBY] = pd.Categorical(["0", "1"] * (n_obs // 2))

    values = rng.normal(0, 0.1, size=(n_obs, n_windows))
    values[values < 0.05] = 0.0                      # genuinely sparse
    adata.obsm["X_cnv"] = sparse.csr_matrix(values) if sparse_cnv else values

    # infercnvpy reads chromosome boundaries from uns[use_rep]["chr_pos"].
    per_chrom = max(1, n_windows // 2)
    adata.uns["cnv"] = {"chr_pos": {"chr1": 0, "chr2": per_chrom}}
    return adata


@pytest.mark.parametrize("n_windows", [5, 40])
def test_a_narrow_sparse_cnv_matrix_still_plots(n_windows):
    """The regression: 5 windows (8 genes mapped) crashed on save."""
    adata = _adata_with_cnv(n_windows)
    figure = make_cnv_heatmap(adata, GROUPBY)
    assert figure is not None
    assert figure.get_axes()


def test_the_sparse_matrix_is_put_back_after_plotting():
    """Densifying is a plotting detail; the caller's object must not change.

    ``adata_cnv`` is the live session object — it is persisted, re-plotted at
    other resolutions and handed to the CopyKAT path, so silently swapping a
    dense array into it would leak memory the run never asked for.
    """
    adata = _adata_with_cnv(5)
    original = adata.obsm["X_cnv"]

    make_cnv_heatmap(adata, GROUPBY)

    assert adata.obsm["X_cnv"] is original
    assert sparse.issparse(adata.obsm["X_cnv"])


def test_a_dense_cnv_matrix_is_left_alone():
    adata = _adata_with_cnv(5, sparse_cnv=False)
    original = adata.obsm["X_cnv"]

    make_cnv_heatmap(adata, GROUPBY)

    assert adata.obsm["X_cnv"] is original


def test_the_failure_this_guards_against_is_real():
    """Pin the upstream behaviour, so the workaround is removed for a reason.

    If a future scanpy/pandas densifies here, this fails and the context
    manager in ``cnv_analysis`` can go.
    """
    matrix = sparse.csr_matrix(np.arange(20, dtype="float64").reshape(4, 5))
    frame = pd.DataFrame(matrix)

    assert frame.shape == (4, 1), "pd.DataFrame(csr) no longer builds an object column"
    assert frame.dtypes.iloc[0] == object
    with pytest.raises(TypeError):
        frame.groupby(pd.Categorical(["a", "a", "b", "b"]), observed=True).mean()


# ── why there were only 5 windows in the first place ─────────────────────────

def test_a_panel_that_barely_maps_is_rejected_before_it_produces_a_result():
    """The run that crashed above had 8 of 5006 genes mapped, and said nothing.

    InSituCNV's default gene-position reference is human; the dataset was mouse,
    so only symbols spelled the same in both nomenclatures matched. inferCNV
    reads copy number from runs of neighbouring genes, so 8 scattered genes is
    not a weak signal — it is noise that looks like a result.
    """
    from palms.utils.cnv_analysis import GeneMappingError, check_gene_mapping

    with pytest.raises(GeneMappingError) as excinfo:
        check_gene_mapping(8, 5006, ["A1cf", "A2m", "Aatf", "Abca1"])

    message = str(excinfo.value)
    assert "8 of 5006" in message
    assert "mouse" in message                 # the casing heuristic fired
    assert "reference is human" in message


def test_a_human_panel_with_the_same_mapping_rate_gets_no_species_hint():
    """The species guess is a heuristic; the rejection is not."""
    from palms.utils.cnv_analysis import GeneMappingError, check_gene_mapping

    with pytest.raises(GeneMappingError) as excinfo:
        check_gene_mapping(8, 5006, ["A1CF", "A2M", "AATF", "ABCA1"])

    assert "mouse" not in str(excinfo.value)


@pytest.mark.parametrize("n_mapped,n_total", [
    (4800, 5006),      # a normal human panel
    (300, 5006),       # poor but usable — 6%, above the floor
    (0, 0),            # degenerate input must not raise here
])
def test_a_usable_mapping_is_left_alone(n_mapped, n_total):
    from palms.utils.cnv_analysis import check_gene_mapping
    check_gene_mapping(n_mapped, n_total, ["A1CF"])      # must not raise


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
