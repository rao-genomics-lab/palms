"""Regression tests for the cluster-filter path used by Cells > Coloring.

Guards the bug where selecting a CNV clustering with the "Filter by cluster"
checkbox engaged blanked *every* cell. CNV clusterings carry **string**
categories ('0','1','2', or 'tumor'/'normal'/'unknown'), whereas ordinary
Leiden clusterings carry **integer** categories. The filter checkboxes
(``_repopulate_cluster_checkboxes``) used to coerce the raw ids to ``int``,
while ``_cluster_raw_to_id`` (built by factorizing the string categorical in
``_get_cluster_ids_per_obs``) is keyed by the raw strings — so
``_translate_selected_ids_to_int`` matched nothing, returned ``[]``, and the
mask ``~np.isin(label_to_cluster, [])`` removed all cells.

These functions are ``ctx``/Qt-bound closures in ``tabs/_helpers.py``, so the
tests replicate the exact three code paths (checkbox id typing, raw->id map,
translate) and assert they stay mutually consistent. Pure numpy/pandas.

Run with:  pytest tests/test_cluster_filter.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ── replicas of the real code paths in tabs/_helpers.py ─────────────────────

def _checkbox_ids(series: pd.Series):
    """Mirror of the FIXED _repopulate_cluster_checkboxes id typing.

    Sort numerically for display when possible, but keep the raw category type.
    """
    raw_ids = series.dropna().unique().tolist()
    try:
        return sorted(raw_ids, key=lambda x: int(x))
    except (ValueError, TypeError):
        return sorted(raw_ids, key=lambda x: str(x))


def _cluster_codes_and_raw_map(series: pd.Series, cell_ids):
    """Mirror of _get_cluster_ids_per_obs: reindex then fast-path or factorize.

    Returns (cluster_values, raw_to_id) where raw_to_id is None on the fast path.
    """
    aligned = series.reindex(cell_ids)
    try:
        cluster_values = aligned.fillna(-1).values.astype(np.int32)
        raw_to_id = None
    except (ValueError, TypeError):
        codes, uniques = pd.factorize(aligned.values)
        cluster_values = codes.astype(np.int32)
        raw_to_id = {u: int(i) for i, u in enumerate(uniques)}
    return cluster_values, raw_to_id


def _translate(selected_ids, raw_to_id):
    """Mirror of _translate_selected_ids_to_int."""
    if raw_to_id is None:
        return list(selected_ids)
    return [raw_to_id[sid] for sid in selected_ids if sid in raw_to_id]


# ── fixtures ────────────────────────────────────────────────────────────────

def _string_categorical_series(cell_ids):
    """A CNV-style clustering: string categories, some cells unassigned (NaN)."""
    vals = ["0", "1", "2", "0", "2", None]
    return pd.Series(pd.Categorical(vals), index=cell_ids, name="cnv_leiden_res0.2")


def _int_categorical_series(cell_ids):
    """An ordinary Leiden clustering: integer categories."""
    vals = [0, 1, 2, 0, 2, 0]
    return pd.Series(pd.Categorical(vals), index=cell_ids, name="leiden_r1.0")


# ── tests ───────────────────────────────────────────────────────────────────

def test_string_clustering_checkbox_ids_match_raw_to_id_keys():
    cell_ids = [f"cell{i}" for i in range(6)]
    s = _string_categorical_series(cell_ids)
    _, raw_to_id = _cluster_codes_and_raw_map(s, cell_ids)

    # String categorical must take the factorize path with a raw->id map.
    assert raw_to_id is not None
    assert all(isinstance(k, str) for k in raw_to_id)

    ids = _checkbox_ids(s)
    # The fix keeps raw string ids (not int-coerced) so they line up with raw_to_id.
    assert all(isinstance(i, str) for i in ids)
    assert set(ids) == set(raw_to_id.keys())


def test_string_clustering_filter_all_selected_keeps_all_labeled_cells():
    """The reported bug: filter engaged -> translate returns [] -> all blanked."""
    cell_ids = [f"cell{i}" for i in range(6)]
    s = _string_categorical_series(cell_ids)
    cluster_values, raw_to_id = _cluster_codes_and_raw_map(s, cell_ids)

    ids = _checkbox_ids(s)               # all checkboxes checked
    int_ids = _translate(ids, raw_to_id)

    # Must resolve to the full set of codes (regression: this used to be []).
    assert int_ids != []
    assert set(int_ids) == set(int(c) for c in cluster_values if c >= 0)

    # np.isin keeps every assigned cell; only the NaN cell (code -1) drops out.
    keep = np.isin(cluster_values, int_ids)
    assert keep.sum() == int((cluster_values >= 0).sum()) == 5


def test_string_clustering_single_cluster_filter_selects_only_that_cluster():
    cell_ids = [f"cell{i}" for i in range(6)]
    s = _string_categorical_series(cell_ids)
    cluster_values, raw_to_id = _cluster_codes_and_raw_map(s, cell_ids)

    int_ids = _translate(["1"], raw_to_id)      # only cluster "1"
    keep = np.isin(cluster_values, int_ids)
    # Exactly the cells whose raw label is "1".
    expected = np.array([str(v) == "1" for v in s.reindex(cell_ids).values])
    assert np.array_equal(keep, expected)
    assert keep.sum() == 1


def test_integer_clustering_uses_fast_path_and_passthrough():
    cell_ids = [f"cell{i}" for i in range(6)]
    s = _int_categorical_series(cell_ids)
    cluster_values, raw_to_id = _cluster_codes_and_raw_map(s, cell_ids)

    # Integer categorical stays on the fast path — no raw->id remap needed.
    assert raw_to_id is None

    ids = _checkbox_ids(s)
    assert ids == [0, 1, 2]
    int_ids = _translate(ids, raw_to_id)        # passthrough
    keep = np.isin(cluster_values, int_ids)
    assert keep.all()                            # every cell assigned + kept
