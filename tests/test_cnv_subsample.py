"""Unit tests for cnv_analysis.subsample_indices — the CopyKAT budget split.

Pure/deterministic (numpy + pandas only). Locks in the invariant fixed in the
recent CopyKAT commits: a large reference cluster must not starve the analyzed
cells out of the subsample.

Run with:  pytest tests/test_cnv_subsample.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from xenium_viewer.utils.cnv_analysis import (
    subsample_indices, _MIN_REFERENCE_CELLS,
)


def _series(labels):
    """Build a (reference_series, obs_index) pair from a list of cluster labels."""
    idx = pd.Index([f"cell{i}" for i in range(len(labels))])
    return pd.Series(labels, index=idx), idx


def test_all_kept_when_under_budget():
    s, idx = _series(["ref", "tumor", "ref", "tumor"])
    keep = subsample_indices(s, ["ref"], ["tumor"], idx, max_cells=100)
    assert keep.dtype == bool
    assert keep.all()


def test_analyzed_cells_are_prioritized_over_large_reference():
    # 9000 reference cells vs 1000 analyzed; budget 2000. The old behaviour let
    # the huge reference consume the whole subsample, starving analyzed cells.
    labels = ["ref"] * 9000 + ["tumor"] * 1000
    s, idx = _series(labels)
    keep = subsample_indices(s, ["ref"], ["tumor"], idx, max_cells=2000, seed=0)

    is_tumor = np.array([lab == "tumor" for lab in labels])
    assert keep[is_tumor].all(), "every analyzed cell must survive the subsample"
    assert keep.sum() == 2000, "budget fully used"
    assert keep[~is_tumor].sum() == 1000, "leftover budget tops up the reference"
    # A modest reference baseline is still reserved.
    assert keep[~is_tumor].sum() >= _MIN_REFERENCE_CELLS


def test_out_of_scope_labels_are_never_selected():
    labels = ["ref"] * 600 + ["tumor"] * 300 + ["stroma"] * 600
    s, idx = _series(labels)
    keep = subsample_indices(s, ["ref"], ["tumor"], idx, max_cells=800, seed=0)

    is_stroma = np.array([lab == "stroma" for lab in labels])
    is_tumor = np.array([lab == "tumor" for lab in labels])
    assert not keep[is_stroma].any(), "cells outside the analysis scope are excluded"
    assert keep[is_tumor].all()
    assert keep.sum() == 800


def test_seeded_rng_is_deterministic():
    labels = ["ref"] * 600 + ["tumor"] * 400
    s, idx = _series(labels)
    keep1 = subsample_indices(s, ["ref"], ["tumor"], idx, max_cells=700, seed=42)
    keep2 = subsample_indices(s, ["ref"], ["tumor"], idx, max_cells=700, seed=42)
    assert np.array_equal(keep1, keep2)
    assert keep1.sum() == 700


def test_no_analyze_ids_uses_all_non_reference_as_scope():
    labels = ["ref"] * 600 + ["other"] * 600
    s, idx = _series(labels)
    keep = subsample_indices(s, ["ref"], None, idx, max_cells=800, seed=0)
    assert keep.sum() == 800  # every cell is in scope; budget respected
