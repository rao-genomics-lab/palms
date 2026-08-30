"""The recorded cell that reloads a stored clustering rather than recomputing it.

Two callers need the identical text and must not drift: ``_record_clustering``
emits it for a column whose producer never recorded code, and a Crop Dataset
export substitutes it for a ``read_csv`` of a file the export does not carry.

**A reloaded clustering is not a recomputed one.** The cell reproduces the
labels; it does not re-derive them, so an ARI quoted off a notebook that
reloaded its clustering is measuring the round trip, not the analysis. The
comment says so in the notebook, where the reader is.
"""

from __future__ import annotations

from palms.utils.adata_persistence import CLUSTERING_PREFIX


def reload_clustering_code(key: str, cache_path: str, *, reason: str) -> str:
    """Code that reads ``clustering_<key>`` back out of a zarr store's obs.

    *reason* is the first line of the comment — why there is no producer to
    replay — and is the only part that differs between the two call sites.
    """
    return (
        f"\n# {reason}\n"
        f"# This RELOADS the stored labels from the viewer's cache —\n"
        f"# it does not recompute them.\n"
        f"import zarr\n"
        f"from anndata.io import read_elem\n"
        f"_cached_obs = read_elem(zarr.open(r\"{cache_path}\", mode='r')['tables/table/obs'])\n"
        f"adata.obs[\"{key}\"] = pd.Categorical(\n"
        f"    _cached_obs[\"{CLUSTERING_PREFIX}{key}\"].reindex(adata.obs_names).astype(str).values)"
    )
