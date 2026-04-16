"""
Loaders for phikon patch-cluster overlays and subclone prediction CSVs.

Both data sources place a regular patch grid over an image and colour each
patch by a cluster id. This module provides lightweight loaders plus a
patch-size inference helper that inspects folder / file names.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class PatchOverlayData:
    coords_xy: np.ndarray            # (N, 2) float64, top-left (x, y)
    patch_size: int                  # pixels
    cluster_columns: dict            # {column_name: np.ndarray int64}
    active_cluster_column: str
    source_path: str
    source_kind: str                 # "phikon" | "subclone"
    confidence: Optional[np.ndarray] = None
    extra: dict = field(default_factory=dict)


# ─── Patch-size inference ────────────────────────────────────────────────────

_PX_TOKEN = re.compile(r"(\d{2,4})\s*px", re.IGNORECASE)
_INT_TOKEN = re.compile(r"(?<!\d)(\d{2,4})(?!\d)")


def _candidate_sizes(name: str) -> list[int]:
    hits_px = [int(m.group(1)) for m in _PX_TOKEN.finditer(name)]
    if hits_px:
        return hits_px
    return [int(m.group(1)) for m in _INT_TOKEN.finditer(name)
            if 32 <= int(m.group(1)) <= 4096]


def infer_patch_size_from_path(path: Path) -> Optional[int]:
    """Guess patch size from a file/folder name.

    Preference order: tokens of the form ``<N>px`` in the leaf, then the
    parent; bare integers in [32, 4096] as a fallback. Returns None if the
    inferred sizes disagree.
    """
    p = Path(path)
    candidates = _candidate_sizes(p.name) or _candidate_sizes(p.parent.name)
    if not candidates:
        return None
    uniq = set(candidates)
    if len(uniq) == 1:
        return candidates[0]
    # Prefer powers of two if multiple are present (common patch convention)
    pow2 = [s for s in uniq if (s & (s - 1)) == 0]
    if len(pow2) == 1:
        return pow2[0]
    return None


def estimate_stride(coords_xy: np.ndarray) -> Optional[int]:
    """Estimate grid stride from neighbouring x-coordinates (mode of diffs)."""
    if coords_xy.shape[0] < 2:
        return None
    xs = np.sort(coords_xy[:, 0].astype(np.int64))
    diffs = np.diff(xs)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return None
    vals, counts = np.unique(diffs, return_counts=True)
    return int(vals[np.argmax(counts)])


# ─── Loaders ─────────────────────────────────────────────────────────────────

def load_phikon_folder(folder: Path) -> PatchOverlayData:
    """Load phikon patches + cluster labels from a results folder.

    Expected layout::

        <folder>/patches/coordinates.npy      # (N, 2) int64, (x, y) top-left
        <folder>/clustering/cluster_labels.npy  # (N,) int64
    """
    folder = Path(folder)
    coord_path = folder / "patches" / "coordinates.npy"
    label_path = folder / "clustering" / "cluster_labels.npy"
    if not coord_path.exists():
        raise FileNotFoundError(f"Missing {coord_path}")
    if not label_path.exists():
        raise FileNotFoundError(f"Missing {label_path}")

    coords = np.asarray(np.load(coord_path), dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f"coordinates.npy must be (N, 2); got {coords.shape}")
    labels = np.asarray(np.load(label_path)).astype(np.int64).ravel()
    if labels.shape[0] != coords.shape[0]:
        raise ValueError(
            f"Label count {labels.shape[0]} != coord count {coords.shape[0]}"
        )

    size = infer_patch_size_from_path(folder)
    return PatchOverlayData(
        coords_xy=coords,
        patch_size=size or 0,
        cluster_columns={"phikon_cluster": labels},
        active_cluster_column="phikon_cluster",
        source_path=str(folder),
        source_kind="phikon",
    )


_SUBCLONE_COLS = ("predicted_genomic_cluster", "morphology_cluster")


def load_subclone_csv(
    csv_path: Path,
    cluster_column: str = "predicted_genomic_cluster",
) -> PatchOverlayData:
    """Load a subclone-prediction CSV with ``x_coord, y_coord`` patch grid."""
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)
    required = {"x_coord", "y_coord"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path.name} missing columns {missing}")

    coords = df[["x_coord", "y_coord"]].to_numpy(dtype=np.float64)

    cluster_columns = {}
    for col in _SUBCLONE_COLS:
        if col not in df.columns:
            continue
        raw = df[col]
        if pd.api.types.is_numeric_dtype(raw):
            vals = raw.fillna(-1).to_numpy(dtype=np.int64)
        else:
            codes, _ = pd.factorize(raw, sort=True)
            vals = codes.astype(np.int64)
        cluster_columns[col] = vals

    if not cluster_columns:
        raise ValueError(
            f"{csv_path.name} has none of {_SUBCLONE_COLS}"
        )
    if cluster_column not in cluster_columns:
        cluster_column = next(iter(cluster_columns))

    # Drop rows whose active cluster is NaN (-1) for cleanliness.
    active = cluster_columns[cluster_column]
    keep = active >= 0
    if not keep.all():
        coords = coords[keep]
        cluster_columns = {k: v[keep] for k, v in cluster_columns.items()}
        df = df.loc[keep].reset_index(drop=True)

    confidence = None
    if "prediction_confidence" in df.columns:
        confidence = df["prediction_confidence"].to_numpy(dtype=np.float32)

    size = infer_patch_size_from_path(csv_path.parent) \
        or infer_patch_size_from_path(csv_path)

    return PatchOverlayData(
        coords_xy=coords,
        patch_size=size or 0,
        cluster_columns=cluster_columns,
        active_cluster_column=cluster_column,
        source_path=str(csv_path),
        source_kind="subclone",
        confidence=confidence,
    )
