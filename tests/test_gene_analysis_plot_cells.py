"""The recorded plot cells must run against the state the notebook has.

Found by replaying a real session: `plot:dotplot:leiden_igraph_r1.0` recorded
`sc.pl.rank_genes_groups_dotplot(adata, ...)` and died at cell 29 of 39 with
`KeyError: 'rank_genes_groups'`. The ranked genes are written to `adata_norm` by
the rank-genes step; these terminals were written back when the viewer
normalised `adata` in place and were never updated.

`tests/test_recorded_code_is_code.py` guards the *shape* of the mistake
statically. This one executes the cell — the only way to know it works.

Run headless: ``MPLBACKEND=Agg pytest tests/test_gene_analysis_plot_cells.py``
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from xenium_viewer.tabs.tab_gene_analysis import dotplot_code  # noqa: E402

GROUPBY = "leiden_igraph_r1.0"


@pytest.fixture(scope="module")
def ranked():
    """An `adata_norm` in the state the notebook has at the dotplot cell."""
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(80, 30)).astype("float32")
    counts[:40, :10] += 8          # something for the ranking to find
    adata_norm = sc.AnnData(counts)
    adata_norm.var_names = [f"GENE{i}" for i in range(30)]
    adata_norm.obs[GROUPBY] = pd.Categorical(["0"] * 40 + ["1"] * 40)
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.tl.rank_genes_groups(adata_norm, groupby=GROUPBY, method="wilcoxon")
    return adata_norm


@pytest.mark.parametrize("dendrogram", [True, False])
def test_the_recorded_dotplot_cell_runs_and_writes_its_figure(
        ranked, tmp_path, monkeypatch, dendrogram):
    monkeypatch.chdir(tmp_path)          # the cell saves to a relative path
    code = dotplot_code(GROUPBY, n_genes=5, dendrogram=dendrogram, fmt="png")

    exec(compile(code, "<plot:dotplot>", "exec"),  # noqa: S102
         {"sc": sc, "np": np, "pd": pd, "adata_norm": ranked})

    assert (tmp_path / "dotplot.png").exists()


def test_the_cell_never_names_raw_adata(ranked):
    """The regression itself: with only `adata` bound, it must not be runnable.

    Executing against a namespace that has *only* the raw object is how the
    notebook failed — the name resolved, the `uns` key did not.
    """
    code = dotplot_code(GROUPBY, n_genes=5, dendrogram=False, fmt="png")
    assert "adata_norm" in code
    assert "(adata," not in code and "(adata " not in code

    with pytest.raises(NameError):
        exec(compile(code, "<plot:dotplot>", "exec"),  # noqa: S102
             {"sc": sc, "adata": ranked})


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
