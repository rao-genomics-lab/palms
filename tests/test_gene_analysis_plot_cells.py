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
from xenium_viewer.utils.gene_analysis import rank_genes_key  # noqa: E402

GROUPBY = "leiden_igraph_r1.0"


@pytest.fixture(scope="module")
def ranked():
    """An `adata_norm` in the state the notebook has at the dotplot cell.

    Ranked with ``key_added``, as the rank-genes step does — the dotplot cell
    passes the matching ``key=``, and ranking here the old unkeyed way would
    test the cell against a state no notebook ever has.
    """
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(80, 30)).astype("float32")
    counts[:40, :10] += 8          # something for the ranking to find
    adata_norm = sc.AnnData(counts)
    adata_norm.var_names = [f"GENE{i}" for i in range(30)]
    adata_norm.obs[GROUPBY] = pd.Categorical(["0"] * 40 + ["1"] * 40)
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.tl.rank_genes_groups(adata_norm, groupby=GROUPBY, method="wilcoxon",
                            key_added=rank_genes_key(GROUPBY))
    return adata_norm


@pytest.mark.parametrize("dendrogram", [True, False])
def test_the_recorded_dotplot_cell_runs_and_writes_its_figure(
        ranked, tmp_path, monkeypatch, dendrogram):
    monkeypatch.chdir(tmp_path)          # the cell saves to relative paths
    # Both formats, as the viewer now writes them: the cell names the files
    # that were actually produced rather than a bare "dotplot.svg".
    code = dotplot_code(GROUPBY, n_genes=5, dendrogram=dendrogram,
                        paths=["plots/dotplot.png", "plots/dotplot.pdf"])

    exec(compile(code, "<plot:dotplot>", "exec"),  # noqa: S102
         {"sc": sc, "np": np, "pd": pd, "Path": Path, "adata_norm": ranked})

    assert (tmp_path / "plots" / "dotplot.png").exists()
    assert (tmp_path / "plots" / "dotplot.pdf").exists()


def test_the_cell_never_names_raw_adata(ranked):
    """The regression itself: with only `adata` bound, it must not be runnable.

    Executing against a namespace that has *only* the raw object is how the
    notebook failed — the name resolved, the `uns` key did not.
    """
    code = dotplot_code(GROUPBY, n_genes=5, dendrogram=False,
                        paths=["plots/dotplot.png"])
    assert "adata_norm" in code
    assert "(adata," not in code and "(adata " not in code

    with pytest.raises(NameError):
        exec(compile(code, "<plot:dotplot>", "exec"),  # noqa: S102
             {"sc": sc, "Path": Path, "adata": ranked})


def test_the_cell_reaches_the_keyed_ranking(ranked):
    """The step writes `uns['rank_genes_<clustering>']`, not scanpy's default.

    Same shape of failure as the `adata`/`adata_norm` one above, one layer in:
    every *name* resolves, and the cell dies on the `uns` lookup. Asserted both
    on the string and by running it against an `adata_norm` that has only the
    unkeyed slot — which is what the notebook would have if the cell and the
    step disagreed about the key.
    """
    code = dotplot_code(GROUPBY, n_genes=5, dendrogram=False,
                        paths=["plots/dotplot.png"])
    assert f'key="{rank_genes_key(GROUPBY)}"' in code

    legacy = ranked.copy()
    legacy.uns["rank_genes_groups"] = legacy.uns.pop(rank_genes_key(GROUPBY))
    with pytest.raises(KeyError):
        exec(compile(code, "<plot:dotplot>", "exec"),  # noqa: S102
             {"sc": sc, "np": np, "pd": pd, "Path": Path, "adata_norm": legacy})


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
