"""
Multichannel → RGB composite builder for the External Images tab.

Builds a dask pyramid of (Y, X, 3) uint8 RGB arrays from a (C, Y, X) source
pyramid + per-channel state (visible, color, contrast limits). The composite
is a lazy dask expression — napari evaluates only the tiles it needs.
"""

from __future__ import annotations

import numpy as np
import dask.array as da


# ── Default IF channel palette ───────────────────────────────────────────────
# Order: blue (DAPI), green, red, magenta, cyan, yellow, white, orange,
# then cycle through a secondary set.
_IF_COLORS = [
    [0.0, 0.0, 1.0],   # blue
    [0.0, 1.0, 0.0],   # green
    [1.0, 0.0, 0.0],   # red
    [1.0, 0.0, 1.0],   # magenta
    [0.0, 1.0, 1.0],   # cyan
    [1.0, 1.0, 0.0],   # yellow
    [1.0, 1.0, 1.0],   # white
    [1.0, 0.5, 0.0],   # orange
    [0.5, 0.0, 1.0],   # violet
    [0.0, 0.8, 0.5],   # teal
    [1.0, 0.4, 0.7],   # pink
    [0.6, 0.8, 0.2],   # lime
]


def default_channel_colors(channel_names: list[str]) -> list[list[float]]:
    """Auto-assign RGB colors based on channel names.

    "DAPI" (case-insensitive) → blue; remaining channels cycle through the
    IF palette starting at green.
    """
    n = len(channel_names)
    colors = [None] * n
    # Assign DAPI first
    non_dapi = list(range(n))
    for i, name in enumerate(channel_names):
        if "dapi" in name.lower():
            colors[i] = list(_IF_COLORS[0])  # blue
            non_dapi.remove(i)
            break
    # Cycle remaining through palette (skip blue at index 0)
    palette_idx = 1
    for i in non_dapi:
        colors[i] = list(_IF_COLORS[palette_idx % len(_IF_COLORS)])
        palette_idx += 1
    return colors


def auto_contrast(pyramid_cyx: list, channel_axis: int = 0,
                  low_pct: float = 1.0, high_pct: float = 99.9) -> list[tuple[float, float]]:
    """Compute per-channel contrast limits from the lowest-resolution pyramid level.

    Returns a list of ``(min, max)`` tuples, one per channel.
    """
    # Use the smallest (last) level — should be small enough to compute
    arr = pyramid_cyx[-1]
    if hasattr(arr, "compute"):
        arr = arr.compute()
    arr = np.asarray(arr, dtype=np.float32)
    n_channels = arr.shape[channel_axis]
    limits = []
    for c in range(n_channels):
        ch = np.take(arr, c, axis=channel_axis).ravel()
        lo = float(np.percentile(ch, low_pct))
        hi = float(np.percentile(ch, high_pct))
        if hi <= lo:
            hi = lo + 1.0
        limits.append((lo, hi))
    return limits


def build_composite_pyramid(
    pyramid_cyx: list,
    channel_states: list[dict],
    channel_axis: int = 0,
) -> list:
    """Build a lazy (Y, X, 3) uint8 RGB composite for each pyramid level.

    Parameters
    ----------
    pyramid_cyx : list of dask.array
        Raw pyramid levels, each with shape ``(..., C, ..., Y, X)`` where
        ``channel_axis`` locates C.
    channel_states : list of dict
        Per-channel state: ``{"visible": bool, "color": [r, g, b], "clim": [min, max]}``.
    channel_axis : int
        Axis index of the channel dimension.

    Returns
    -------
    list of dask.array
        One ``(Y, X, 3)`` uint8 dask array per pyramid level.
    """
    composite_pyramid = []
    for level in pyramid_cyx:
        level_da = da.asarray(level)
        h = level_da.shape[channel_axis + 1] if channel_axis == 0 else level_da.shape[0]
        w = level_da.shape[-1]
        # Spatial shape (Y, X) — assumes channel_axis is 0 and remaining are (Y, X)
        spatial_shape = tuple(
            s for i, s in enumerate(level_da.shape) if i != channel_axis
        )

        layers = []
        for i, cs in enumerate(channel_states):
            if not cs.get("visible", True):
                continue
            # Extract single channel: (Y, X) dask
            ch = da.take(level_da, i, axis=channel_axis)
            lo, hi = cs["clim"]
            lo, hi = float(lo), float(hi)
            color = np.array(cs["color"][:3], dtype=np.float32).reshape(1, 1, 3)
            # Normalize and colorize
            norm = da.clip(
                (ch.astype(np.float32) - lo) / (hi - lo + 1e-8), 0.0, 1.0
            )
            rgb = norm[:, :, None] * color  # (Y, X, 3) float32
            layers.append(rgb)

        if layers:
            composite = da.clip(sum(layers) * 255.0, 0, 255).astype(np.uint8)
        else:
            composite = da.zeros(spatial_shape + (3,), dtype=np.uint8,
                                 chunks=(-1, -1, 3))
        composite_pyramid.append(composite)
    return composite_pyramid
