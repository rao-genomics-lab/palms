"""Point every derived viewer object at a new cell set.

Two features rebind ``ctx.adata``: Tools -> Segmentation swaps the Xenium cells
for a custom segmentation, and Tools -> QC drops the cells and genes a filter
rejects. Both leave the same debris behind, because almost everything the viewer
derives from the table is indexed by **obs row position** rather than by cell id:
``label_to_obs``, ``CellColorManager``'s colour scatter, ``centroids_yx``, the
UMAP window's frozen ``n_cells``, and every gene ComboBox bound to
``var_names`` at build time.

Getting that wrong does not raise. ``CellColorManager`` does

    color_arr[valid_labels] = rgba_obs[obs_indices]

so a stale ``label_to_obs`` paints each cell with some *other* cell's value and
reports nothing. That is the reason this module exists as one shared routine
rather than as a checklist each tab remembers to follow.

Qt-free and I/O-free on purpose, so it can be tested without a viewer.
"""
from __future__ import annotations

import numpy as np

#: Cached analysis results that are about the old cell set. Popped rather than
#: recomputed: the tabs rebuild them on demand, and a stale zscore matrix shown
#: beside a new cell set is worse than an empty pane.
STALE_RESULT_KEYS = (
    "nhood_result", "nhood_fig", "co_result", "co_fig",
    "rank_genes_df", "rank_genes_adata_norm", "rank_genes_groupby",
    "ligrec_result", "annot_dist_distances",
)

#: Memo keys whose value is an ``id()`` of an object derived from the old cells.
STALE_MEMO_KEYS = ("_norm_src_id", "_spatial_neighbors_key", "_qc_applied_key")


def repoint_label_to_obs(label_to_obs: np.ndarray, kept: np.ndarray) -> np.ndarray:
    """Re-point a label->obs map at the rows that survived *kept*.

    ``label_to_obs`` is indexed by the raster's *pixel value*, so the result
    keeps the input's length whatever happens to the table. Shortening it to the
    surviving maximum label would make ``CellColorManager`` build a colormap the
    raster runs off the end of.

    A cell that ``kept`` rejects becomes ``-1``, which every consumer already
    reads as "no cell here" -- so a filtered-out cell renders transparent, which
    is also the clearest possible readout of what a filter did.

    Parameters
    ----------
    label_to_obs
        The map for the *full* table: ``arr[label] = obs row``, ``-1`` for none.
    kept
        Boolean mask over the full table's obs rows, ``True`` for survivors.
    """
    label_to_obs = np.asarray(label_to_obs)
    kept = np.asarray(kept, dtype=bool)
    if kept.size and int(label_to_obs.max(initial=-1)) >= kept.size:
        raise ValueError(
            "label_to_obs points past the end of the mask: it must be the map "
            "for the same table the mask indexes"
        )

    # Old row -> new row, or -1. Built once, then gathered through the labels.
    new_row = np.full(kept.shape, -1, dtype=np.int32)
    new_row[kept] = np.arange(int(kept.sum()), dtype=np.int32)

    out = np.full(label_to_obs.shape, -1, dtype=np.int32)
    valid = label_to_obs >= 0
    out[valid] = new_row[label_to_obs[valid]]
    return out


def kept_mask(full_adata, new_adata) -> np.ndarray:
    """Boolean mask over *full_adata*'s rows, ``True`` where *new_adata* kept it.

    By ``obs_names`` rather than by the filter's own mask, because only one of
    the two QC blocks binds a cell mask -- a gene-only filter binds none -- and
    because a segmentation swap has no mask at all.
    """
    return full_adata.obs_names.isin(new_adata.obs_names).astype(bool)


def rebind_cells(ctx, new_adata, new_label_to_obs):
    """Point every derived object at *new_adata*. Returns nothing.

    ``ctx.clusterings`` is deliberately not touched, because the two callers
    want opposite things and only they know which. A QC filter subsets rows of
    the same table, and the clustering Series are indexed by ``cell_id`` and
    realigned by ``reindex`` at every use, so they stay valid as they are; a
    segmentation swap replaces the cells outright, so ``tab_segmentation``
    reloads them from the new table before calling this. Either way the combo
    refresh at the end picks up whatever the caller left.
    """
    from palms.utils.coloring import CellColorManager

    ctx.adata = new_adata
    ctx.label_to_obs = new_label_to_obs

    # A fresh manager rather than a mutated one: it is what discards the gene,
    # cluster and continuous colour caches, which have no invalidate-all.
    ctx.color_manager = CellColorManager(new_adata, new_label_to_obs)

    if "spatial" in getattr(new_adata, "obsm", {}):
        centroids_px = np.asarray(new_adata.obsm["spatial"], dtype=np.float64) / ctx.pixel_size
        ctx.centroids_yx = centroids_px[:, ::-1]

    ctx.gene_names = list(new_adata.var_names)

    _rebuild_umap_viewer(ctx, new_adata)

    state = ctx.state
    for key in STALE_RESULT_KEYS:
        state.pop(key, None)
    state["label_to_cluster"] = None
    state["active_clustering_name"] = None

    _reset_executor_namespace(ctx)

    refresh_genes = getattr(ctx, "refresh_gene_choices", None)
    if callable(refresh_genes):
        refresh_genes()
    refresh_clusterings = getattr(ctx, "refresh_clustering_choices", None)
    if callable(refresh_clusterings):
        refresh_clusterings()

    for listener in list(state.get("qc_listeners") or ()):
        try:
            listener()
        except Exception:  # noqa: BLE001 - a label must never break a rebind
            pass


def _rebuild_umap_viewer(ctx, new_adata) -> None:
    """Rebuild the UMAP window's index, which is frozen at construction.

    ``UMAPViewer.__init__`` stores ``n_cells`` and reindexes the UMAP frame onto
    the cell ids it was handed; ``_labels_color_arr_to_obs_colors`` then indexes
    against that length. Left alone across a rebind it produces wrong colours
    with no error -- the same silent class as a stale ``label_to_obs``.
    """
    umap_df = getattr(ctx, "umap_df", None)
    if umap_df is None or "cell_id" not in new_adata.obs.columns:
        return
    from palms.utils.umap_widget import UMAPViewer

    old = ctx.umap_viewer
    ctx.umap_viewer = UMAPViewer(umap_df, new_adata.obs["cell_id"].values)
    if old is not None:
        try:
            old.close()
        except Exception:  # noqa: BLE001 - the window may already be gone
            pass


def _reset_executor_namespace(ctx) -> None:
    """Drop every executor binding derived from the old cell set.

    ``adata_norm``, ``adata_leiden``, ``annotations``, ``transcript_points`` and
    the rest are still bound at the *old* cell count -- ``run_step`` re-syncs
    only ``adata`` and ``sdata``. Enumerating the derived names would rot as
    templates are added, so the rule is the honest one: everything not in the
    base namespace is gone, which is also the clean-kernel semantics the
    exported notebook has.

    A later step whose producer has not re-run then fails as a loud ``StepError``
    naming the missing binding, which is correct -- exporting a table computed on
    cells that no longer exist is the failure this prevents.
    """
    from palms.utils.step_templates.namespace import EXECUTOR_BASE_NAMES

    executor = getattr(ctx, "executor", None)
    if executor is not None:
        for name in [n for n in executor.ns if n not in EXECUTOR_BASE_NAMES]:
            executor.ns.pop(name, None)
        executor.ns["adata"] = ctx.adata
        executor.ns["sdata"] = ctx.sdata

    for key in STALE_MEMO_KEYS:
        ctx.state.pop(key, None)
