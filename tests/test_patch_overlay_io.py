"""Unit tests for patch_overlay_io size/stride inference (pure regex + numpy).

Run with:  pytest tests/test_patch_overlay_io.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from xenium_viewer.utils.patch_overlay_io import (
    _candidate_sizes, infer_patch_size_from_path, estimate_stride,
)


def test_candidate_sizes_px_tokens_win():
    assert _candidate_sizes("patches_256px") == [256]
    assert _candidate_sizes("run_512px_v3") == [512]


def test_candidate_sizes_bare_ints_only_in_range():
    assert _candidate_sizes("tiles_128") == [128]
    assert _candidate_sizes("scale_16") == []      # below the 32-4096 window
    assert _candidate_sizes("no_numbers_here") == []


def test_infer_from_leaf_name():
    assert infer_patch_size_from_path(Path("/data/phikon_256px")) == 256


def test_infer_falls_back_to_parent_name():
    assert infer_patch_size_from_path(Path("/data/256px/coords.npy")) == 256


def test_infer_prefers_power_of_two_on_conflict():
    assert infer_patch_size_from_path(Path("/data/256_300")) == 256


def test_infer_none_when_sizes_disagree():
    assert infer_patch_size_from_path(Path("/data/300_500")) is None


def test_estimate_stride_regular_grid():
    coords = np.array([[0, 0], [64, 0], [128, 0], [192, 10]], dtype=float)
    assert estimate_stride(coords) == 64


def test_estimate_stride_uses_mode_of_diffs():
    coords = np.array([[0, 0], [64, 0], [128, 0], [500, 0]], dtype=float)
    assert estimate_stride(coords) == 64  # diffs 64, 64, 372 -> mode 64


def test_estimate_stride_too_few_points():
    assert estimate_stride(np.array([[0, 0]], dtype=float)) is None
