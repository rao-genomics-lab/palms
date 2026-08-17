"""Rankings are keyed per clustering, and older caches still open.

``sc.tl.rank_genes_groups`` overwrites ``uns['rank_genes_groups']`` in place, so
a session that ranked two clusterings kept only the last — and the notebook's
workaround (a local ``rank_results`` dict) put the frames somewhere no
``sc.pl.rank_genes_groups*`` call could reach, since they all take ``key=``.
The step template passes ``key_added=`` instead, which fixes both at once.

That makes the ``uns`` slot name derived rather than fixed, and three things
downstream have to agree about it: what the viewer persists, what counts as
user data worth warning about before a cache rebuild, and what Tools → Dataset
will let you delete. Each of those had a *literal* ``"rank_genes_groups"`` in
it, so each would have silently stopped recognising a ranking.

Every cache written before the keying holds the unkeyed slot and nothing else,
so the read path falls back — the assertions below pin that in both directions.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_rank_genes_keying.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")

from xenium_viewer.utils.gene_analysis import (  # noqa: E402
    LEGACY_RANK_KEY, rank_genes_key, resolve_rank_key,
)

GROUPBY = "leiden_r1.0"


def _ranked(key: str | None):
    """A small ranked AnnData, with the result under *key* (None = unranked)."""
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(40, 12)).astype("float32")
    counts[:20, :4] += 8
    adata = anndata.AnnData(counts)
    adata.var_names = [f"GENE{i}" for i in range(12)]
    adata.obs[GROUPBY] = pd.Categorical(["0"] * 20 + ["1"] * 20)
    if key is not None:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, groupby=GROUPBY, method="wilcoxon",
                                key_added=key)
        adata.uns["rank_genes_groupby"] = GROUPBY
    return adata


# ── the key itself ───────────────────────────────────────────────────────────

def test_the_key_names_the_clustering_it_belongs_to():
    assert rank_genes_key(GROUPBY) == "rank_genes_leiden_r1.0"


def test_a_keyed_result_resolves_to_its_own_slot():
    adata = _ranked(rank_genes_key(GROUPBY))
    assert resolve_rank_key(adata, GROUPBY) == rank_genes_key(GROUPBY)


def test_a_cache_written_before_the_keying_falls_back():
    """The whole point of the fallback: it holds the unkeyed slot and no other."""
    adata = _ranked(LEGACY_RANK_KEY)
    assert rank_genes_key(GROUPBY) not in adata.uns
    assert resolve_rank_key(adata, GROUPBY) == LEGACY_RANK_KEY


def test_an_unranked_object_resolves_to_a_key_it_does_not_have():
    """Callers test membership themselves; resolve never invents a result."""
    adata = _ranked(None)
    assert resolve_rank_key(adata, GROUPBY) == LEGACY_RANK_KEY
    assert LEGACY_RANK_KEY not in adata.uns


def test_resolve_tolerates_no_groupby():
    assert resolve_rank_key(_ranked(None), None) == LEGACY_RANK_KEY


# ── restoring a stored ranking ───────────────────────────────────────────────

@pytest.mark.parametrize("key", [rank_genes_key(GROUPBY), LEGACY_RANK_KEY])
def test_a_stored_ranking_restores_under_either_key(key):
    from xenium_viewer.utils.adata_persistence import load_rank_genes_from_adata

    adata = _ranked(key)
    df, adata_norm, groupby = load_rank_genes_from_adata(adata, None)

    assert groupby == GROUPBY
    assert adata_norm is None                      # no sdata, so no sidecar
    expected = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    pd.testing.assert_frame_equal(df, expected)


def test_an_object_with_no_ranking_restores_nothing():
    from xenium_viewer.utils.adata_persistence import load_rank_genes_from_adata

    assert load_rank_genes_from_adata(_ranked(None), None) == (None, None, None)


# ── the three places that used to hardcode the slot name ─────────────────────

def test_a_keyed_ranking_still_counts_as_user_data():
    """`loader` warns before rebuilding a cache that holds user data.

    A fixed name list would have stopped recognising a ranking the moment the
    keying landed — and "no user data" is what lets a cache be rebuilt with no
    dialog at all.
    """
    from xenium_viewer import loader

    assert loader._is_user_uns(rank_genes_key(GROUPBY))
    assert loader._is_user_uns(LEGACY_RANK_KEY)
    assert loader._is_user_uns("rank_genes_groupby")
    assert not loader._is_user_uns("spatial_neighbors")


def test_a_keyed_ranking_is_still_deletable_in_the_dataset_tab():
    """Unrecognised defaults to *not* deletable, so a missed name is not
    dangerous — it is a result the user is simply unable to remove."""
    from xenium_viewer.utils import store_inventory

    assert store_inventory._is_deletable_uns(rank_genes_key(GROUPBY))
    assert store_inventory._is_deletable_uns(LEGACY_RANK_KEY)
    assert store_inventory._is_deletable_uns("rank_genes_groupby")
    assert not store_inventory._is_deletable_uns("spatial")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
