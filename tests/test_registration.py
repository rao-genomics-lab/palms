"""Unit tests for registration landmark math and landmark JSON I/O.

`compute_landmark_affine` fits a similarity transform; `save_landmarks`/
`load_landmarks` are a JSON round-trip. Both are pure (numpy/skimage/json).

Run with:  pytest tests/test_registration.py
"""
from __future__ import annotations

import numpy as np

from xenium_viewer.utils.registration import (
    compute_landmark_affine, save_landmarks, load_landmarks,
)


def _similarity_yx(pts_yx, scale, theta_deg, translation_yx):
    """Apply a known similarity (scaled rotation + translation) in (y, x)."""
    t = np.deg2rad(theta_deg)
    R = scale * np.array([[np.cos(t), -np.sin(t)],
                          [np.sin(t),  np.cos(t)]])
    return (R @ pts_yx.T).T + np.asarray(translation_yx, dtype=float)


def test_compute_landmark_affine_recovers_known_similarity():
    he = np.array([[0.0, 0.0], [0.0, 10.0], [10.0, 0.0], [5.0, 7.0]])
    xenium = _similarity_yx(he, scale=1.5, theta_deg=30.0, translation_yx=(3.0, -2.0))

    affine, residuals = compute_landmark_affine(xenium, he)

    assert affine.shape == (3, 3)
    assert residuals.shape == (len(he),)
    assert residuals.max() < 1e-6, "an exact similarity must fit with ~zero residual"

    # The returned affine should map H&E (y, x) points onto the Xenium points.
    he_homo = np.hstack([he, np.ones((len(he), 1))])
    recovered = (affine @ he_homo.T).T[:, :2]
    assert np.allclose(recovered, xenium, atol=1e-6)


def test_landmarks_json_roundtrip_full(tmp_path):
    xenium = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    he = np.array([[0.5, 0.5], [1.5, 2.5], [7.0, 8.0]])
    affine = np.eye(3)
    path = tmp_path / "landmarks.json"

    save_landmarks(path, xenium, he, affine=affine, he_filename="slide_he.ome.tif")
    data = load_landmarks(path)

    assert np.allclose(data["xenium_landmarks_yx"], xenium)
    assert np.allclose(data["he_landmarks_yx"], he)
    assert np.allclose(data["affine_3x3_yx"], affine)
    assert data["he_filename"] == "slide_he.ome.tif"


def test_landmarks_json_roundtrip_optional_fields_absent(tmp_path):
    xenium = np.array([[1.0, 2.0], [3.0, 4.0]])
    he = np.array([[0.0, 0.0], [1.0, 1.0]])
    path = tmp_path / "landmarks_min.json"

    save_landmarks(path, xenium, he)  # no affine, no filename
    data = load_landmarks(path)

    assert np.allclose(data["xenium_landmarks_yx"], xenium)
    assert np.allclose(data["he_landmarks_yx"], he)
    assert "affine_3x3_yx" not in data
    assert "he_filename" not in data
