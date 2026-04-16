"""
Affine linking helpers for external overlays.

Provides utilities to mirror one napari layer's affine onto another so that
re-registering a source image (e.g. H&E) automatically updates any external
image / patch overlay that was linked to it.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np


def _affine_matrix(layer) -> np.ndarray:
    """Return the 3x3 affine matrix of ``layer`` (or identity if none)."""
    try:
        m = np.asarray(layer.affine.affine_matrix, dtype=np.float64)
    except Exception:
        return np.eye(3, dtype=np.float64)
    # napari uses (ndim+1, ndim+1); take the 2-D (y, x) portion.
    if m.shape == (3, 3):
        return m
    if m.shape[0] >= 3 and m.shape[1] >= 3:
        return m[-3:, -3:]
    return np.eye(3, dtype=np.float64)


def _is_identity(layer, atol: float = 1e-6) -> bool:
    return np.allclose(_affine_matrix(layer), np.eye(3), atol=atol)


def list_transformable_layers(viewer, exclude=()) -> list:
    """Return layers whose affine is not (approximately) the identity.

    Parameters
    ----------
    viewer : napari.Viewer
    exclude : iterable of napari.layers.Layer
        Layers to skip (e.g. the target itself).
    """
    excluded_ids = {id(l) for l in exclude}
    out = []
    for layer in viewer.layers:
        if id(layer) in excluded_ids:
            continue
        if not _is_identity(layer):
            out.append(layer)
    return out


def link_affine(target, source, viewer=None) -> Callable[[], None]:
    """Copy ``source.affine`` onto ``target`` and mirror future changes.

    Returns a disconnect callable that tears down both the affine-change
    subscription and (if ``viewer`` is provided) the layer-removed guard.
    """
    # Apply current affine once.
    try:
        target.affine = source.affine
    except Exception:
        target.affine = _affine_matrix(source)

    disconnected = {"value": False}

    def _on_affine_change(event=None):
        if disconnected["value"]:
            return
        try:
            target.affine = source.affine
        except Exception:
            try:
                target.affine = _affine_matrix(source)
            except Exception:
                pass

    def _on_source_removed(event=None):
        if event is not None and getattr(event, "value", None) is not source:
            return
        disconnect()

    try:
        source.events.affine.connect(_on_affine_change)
    except Exception:
        pass

    if viewer is not None:
        try:
            viewer.layers.events.removed.connect(_on_source_removed)
        except Exception:
            pass

    def disconnect():
        if disconnected["value"]:
            return
        disconnected["value"] = True
        try:
            source.events.affine.disconnect(_on_affine_change)
        except Exception:
            pass
        if viewer is not None:
            try:
                viewer.layers.events.removed.disconnect(_on_source_removed)
            except Exception:
                pass

    return disconnect


def find_layer_by_name(viewer, name: Optional[str]):
    """Return the first layer whose ``.name`` matches ``name`` (or None)."""
    if not name:
        return None
    for layer in viewer.layers:
        if layer.name == name:
            return layer
    return None
