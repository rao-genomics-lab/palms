"""Results computed while filtered must reach the store, and the store must
keep every cell.

Both halves were live defects the moment ``adata`` could be rebound. Every
``save_*_to_adata`` writes onto ``ctx.adata`` and then calls ``_persist_table``,
which writes the element out of ``sdata`` — the *full* table. While the two are
the same object that works; once a filter makes ``ctx.adata`` a subset, the
column is written to an object nothing persists and is lost at exit, silently.

The opposite fix would have been worse: persisting the subset overwrites the
store's cell set, and ``_persist_custom_table`` would write it straight over
``custom_table`` — the only copy of a cell set the raw output does not contain.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")


def _tables(n=20, g=4, keep=None):
    """A full table and a filtered view of it, as the viewer would hold them."""
    full = anndata.AnnData(
        np.ones((n, g), dtype="float32"),
        obs=pd.DataFrame({"cell_id": np.arange(1, n + 1)},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(g)]),
    )
    if keep is None:
        keep = np.arange(n) % 2 == 0
    return full, full[keep].copy()


class _Ctx:
    """The three attributes the persistence helpers read."""

    def __init__(self, full, sub):
        self.full_adata = full
        self.adata = sub
        self.no_cache = True
        self.sdata = None
        self.segmentation_source = "xenium"


def test_full_table_is_the_unfiltered_one():
    from palms.utils.adata_persistence import full_table

    full, sub = _tables()
    assert full_table(_Ctx(full, sub)) is full


def test_full_table_falls_back_when_nothing_was_filtered():
    """Every dataset opened before this feature has no ``full_adata``."""
    from palms.utils.adata_persistence import full_table

    ctx = _Ctx(None, _tables()[0])
    assert full_table(ctx) is ctx.adata


def test_a_column_computed_while_filtered_lands_on_every_cell():
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, sub = _tables()
    sub.obs["clustering_leiden"] = pd.Categorical(
        [str(i % 3) for i in range(sub.n_obs)])

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    assert full.n_obs == 20, "the store's cell set must not shrink"
    assert "clustering_leiden" in full.obs.columns
    written = full.obs["clustering_leiden"]
    assert written.notna().sum() == sub.n_obs
    assert written.isna().sum() == full.n_obs - sub.n_obs, (
        "a cell the filter dropped has no cluster, and NaN says so"
    )
    for name in sub.obs_names:
        assert written[name] == sub.obs["clustering_leiden"][name]


def test_the_bare_template_written_column_reaches_the_store_too():
    """Leiden writes ``obs[key]`` itself so the notebook reproduces it.

    ``store_inventory._clustering_twin_of`` pairs that bare column with
    ``clustering_<key>`` so a deletion cascades — which only works if the merge
    catches columns a *template* wrote, not just the ones a save function did.
    """
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, sub = _tables()
    sub.obs["leiden_r1.0"] = pd.Categorical(["a"] * sub.n_obs)
    sub.obs["clustering_leiden_r1.0"] = pd.Categorical(["a"] * sub.n_obs)

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    assert "leiden_r1.0" in full.obs.columns
    assert "clustering_leiden_r1.0" in full.obs.columns


def _filtered_ctx(tiny_sdata, keep):
    """A ctx whose ``adata`` is a filtered view of a real, on-disk store's table."""
    full = tiny_sdata["table"]
    ctx = _Ctx(full, full[keep].copy())
    ctx.sdata = tiny_sdata
    ctx.no_cache = False
    return ctx


def test_a_filtered_write_reaches_a_real_store(tmp_path, tiny_sdata):
    """The xv-cbz gate: what ``_persist_table`` does after a filter, end to end.

    The first merge reindexed every column of the filtered ``obs`` -- the
    structural ones included -- so ``instance_id`` was null for every dropped
    cell, ``TableModel`` refused it, and ``_persist_table`` reported the failure
    and carried on. Nothing computed under a filter was ever written; the graph
    kept recording that it had been.
    """
    from spatialdata import read_zarr

    from palms.utils.adata_persistence import _persist_table

    ctx = _filtered_ctx(tiny_sdata, np.array([True, True, False, True, False, True]))
    ctx.adata.obs["clustering_leiden"] = pd.Categorical(["0", "1", "0", "1"])

    _persist_table(ctx)

    obs = read_zarr(tiny_sdata.path)["table"].obs
    assert "clustering_leiden" in obs.columns, "the result never reached the store"
    assert obs["clustering_leiden"].isna().sum() == 2
    assert obs["instance_id"].isna().sum() == 0, "the instance key was nulled"
    assert obs["region"].isna().sum() == 0


def test_structural_columns_keep_their_values_and_dtypes_for_dropped_cells():
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, sub = _tables()
    full.obs["region"] = pd.Categorical(["r"] * full.n_obs)
    full.obs["transcript_counts"] = np.arange(full.n_obs, dtype="int64")
    sub = full[np.arange(full.n_obs) % 2 == 0].copy()
    before = full.obs.copy()

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    pd.testing.assert_frame_equal(full.obs, before)
    assert full.obs["transcript_counts"].dtype == "int64"
    assert full.obs["cell_id"].isna().sum() == 0


def test_a_pre_filter_column_survives_for_dropped_cells():
    """r0.8 was computed on every cell; the filtered view carries only the kept
    rows of it, and reindexing that back used to erase the rest."""
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, _ = _tables()
    full.obs["clustering_r0.8"] = pd.Categorical([str(i % 4) for i in range(full.n_obs)])
    keep = np.arange(full.n_obs) % 2 == 0
    sub = full[keep].copy()
    sub.obs["clustering_r0.9"] = pd.Categorical(["a"] * sub.n_obs)

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    assert full.obs["clustering_r0.8"].isna().sum() == 0
    assert list(full.obs["clustering_r0.8"]) == [str(i % 4) for i in range(full.n_obs)]
    assert full.obs["clustering_r0.9"].isna().sum() == (~keep).sum()


def test_a_column_recomputed_under_the_filter_is_null_for_dropped_cells():
    """Re-run r1.0 under the filter and the dropped cells have no r1.0 any more.

    Not "update the kept rows and leave the rest": that leaves one column
    holding two runs, and nothing downstream could tell.
    """
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, _ = _tables()
    full.obs["clustering_r1.0"] = pd.Categorical(["old"] * full.n_obs)
    keep = np.arange(full.n_obs) % 2 == 0
    sub = full[keep].copy()
    sub.obs["clustering_r1.0"] = pd.Categorical(["new"] * sub.n_obs)

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    col = full.obs["clustering_r1.0"]
    assert (col[keep] == "new").all()
    assert col[~keep].isna().all(), "the old run must not survive on dropped cells"


def test_uns_results_copy_across_whole():
    """They are cluster-level, not per-cell, so there is nothing to reindex."""
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, sub = _tables()
    sub.uns["nhood_enrichment"] = {"zscore": np.eye(3), "clusters": ["1", "2", "3"]}

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    assert "nhood_enrichment" in full.uns
    assert full.uns["nhood_enrichment"]["clusters"] == ["1", "2", "3"]


def test_genes_and_embeddings_are_left_alone():
    """The store keeps the whole panel, even when the analysis dropped genes."""
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, _ = _tables()
    sub = full[:, :2].copy()
    full.obsm["X_umap"] = np.zeros((full.n_obs, 2))

    _sync_filtered_obs_into_full(_Ctx(full, sub))

    assert full.n_vars == 4
    assert full.obsm["X_umap"].shape == (full.n_obs, 2)


def test_syncing_is_a_no_op_when_nothing_is_filtered():
    from palms.utils.adata_persistence import _sync_filtered_obs_into_full

    full, _ = _tables()
    before = list(full.obs.columns)
    _sync_filtered_obs_into_full(_Ctx(full, full))
    assert list(full.obs.columns) == before


def test_the_custom_table_written_is_the_full_one():
    """The destructive case: ``_persist_custom_table`` copies ``ctx.adata``.

    Under a custom segmentation plus a filter, writing the subset would put a
    strict subset of the custom cells into ``custom_table`` — and the raw output
    has no copy of them to rebuild from.
    """
    import ast

    source = (Path(__file__).resolve().parent.parent / "src" / "palms"
              / "utils" / "adata_persistence.py").read_text()
    tree = ast.parse(source)
    func = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_persist_custom_table")
    body = ast.unparse(func)
    assert "full_table(ctx).copy()" in body
    assert "ctx.adata.copy()" not in body


def test_both_persist_paths_sync_first():
    """A save function writes ctx.adata; the element written is the full table."""
    import ast

    source = (Path(__file__).resolve().parent.parent / "src" / "palms"
              / "utils" / "adata_persistence.py").read_text()
    tree = ast.parse(source)
    for name in ("_persist_table", "_persist_custom_table"):
        func = next(n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == name)
        assert "_sync_filtered_obs_into_full(ctx)" in ast.unparse(func), (
            f"{name} must merge the filtered results back before writing, or "
            f"everything computed under a QC filter is lost at exit"
        )


# ── Session round-trip ───────────────────────────────────────────────────────

def _attrs(**state):
    """``_build_session_attrs`` with the arguments it does not care about here."""
    from palms.utils.session import _build_session_attrs

    base = {"segmentation_source": "xenium"}
    base.update(state)
    return _build_session_attrs(state=base, he_state={}, snapshot={}, prev_attrs={})


def test_the_cutoffs_are_written_to_the_session():
    cutoffs = {"min_counts": 10, "min_cells": 3}
    assert _attrs(qc_filter=cutoffs)["qc_filter"] == cutoffs


def test_a_store_written_before_this_feature_saves_no_filter():
    """The opt-in guarantee: an older session must not acquire one."""
    assert _attrs()["qc_filter"] is None


def test_load_session_reads_the_attr_back_defaulting_to_off():
    """``load_session`` does I/O, so the read side is checked at the source.

    The default is what makes every store written before this feature open
    unfiltered rather than raising.
    """
    import ast

    source = (Path(__file__).resolve().parent.parent / "src" / "palms"
              / "utils" / "session.py").read_text()
    tree = ast.parse(source)
    func = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "load_session")
    # ast.unparse normalises quoting, so match its spelling, not the file's.
    assert "'qc_filter': attrs.get('qc_filter')" in ast.unparse(func)


def test_the_cutoffs_are_listed_but_never_deletable():
    """A setting, not data: deleting it would mean nothing."""
    from palms.utils.store_inventory import _BLOCKED_SESSION_ATTRS

    assert "qc_filter" in _BLOCKED_SESSION_ATTRS
