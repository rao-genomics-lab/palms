"""
Per-gene transcript loader with feather cache.

Loads transcript x/y coordinates for a single gene either from:
1. A pre-built feather file in transcript_cache/ (fast: <100ms)
2. The full transcripts.parquet (slow: ~5s, fallback if cache missing)

Usage:
    loader = TranscriptLoader()
    df = loader.load_gene("MMSET")  # returns DataFrame with x, y columns
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Default paths (relative to the data root)
_SCRIPTS_DIR = Path(__file__).parent.parent
_DATA_DIR = _SCRIPTS_DIR.parent
_CACHE_DIR = _DATA_DIR / "transcript_cache"
_PARQUET_PATH = _DATA_DIR / "transcripts.parquet"

# Column names in the parquet / feather files
X_COL = "x_location"
Y_COL = "y_location"
QV_COL = "qv"
GENE_COL = "feature_name"  # Xenium 3.x column name
CELL_COL = "cell_id"       # absent from caches built before it was kept

# Default pixel size (µm/pixel) — overridden by experiment.xenium at runtime
DEFAULT_PIXEL_SIZE = 0.2125


class TranscriptLoader:
    """
    Loads per-gene transcript locations.

    Parameters
    ----------
    cache_dir : Path
        Directory containing per-gene .feather files.
        Created by 00_preprocess_transcripts.py.
    parquet_path : Path
        Full transcripts.parquet path (fallback).
    min_qv : int
        Minimum quality value (applied when reading from parquet fallback).
        Feather files are pre-filtered.
    """

    def __init__(
        self,
        cache_dir: Path = _CACHE_DIR,
        parquet_path: Path = _PARQUET_PATH,
        min_qv: int = 20,
        pixel_size: float = DEFAULT_PIXEL_SIZE,
    ):
        self.cache_dir = cache_dir
        self.parquet_path = parquet_path
        self.min_qv = min_qv
        self.pixel_size = pixel_size
        self._cached_genes: Optional[set] = None

    @property
    def cached_genes(self) -> set[str]:
        """Set of gene names available in the feather cache."""
        if self._cached_genes is None:
            if self.cache_dir.exists():
                self._cached_genes = {p.stem for p in self.cache_dir.glob("*.feather")}
            else:
                self._cached_genes = set()
        return self._cached_genes

    def load_gene(self, gene_name: str) -> pd.DataFrame:
        """
        Load transcript coordinates for a single gene.

        Returns a DataFrame with at minimum columns:
            x_location, y_location

        Parameters
        ----------
        gene_name : str

        Returns
        -------
        pd.DataFrame
        """
        if gene_name in self.cached_genes:
            return self._load_from_feather(gene_name)
        else:
            print(f"Warning: '{gene_name}' not in feather cache. Scanning parquet (slow)...")
            return self._load_from_parquet(gene_name)

    def _load_from_feather(self, gene_name: str) -> pd.DataFrame:
        feather_path = self.cache_dir / f"{gene_name}.feather"
        df = pd.read_feather(feather_path)
        return df

    def _load_from_parquet(self, gene_name: str) -> pd.DataFrame:
        import pyarrow.parquet as pq
        import pyarrow.compute as pc

        pf = pq.ParquetFile(self.parquet_path)
        results = []
        for batch in pf.iter_batches(batch_size=1_000_000):
            df = batch.to_pandas()
            mask = df[GENE_COL] == gene_name
            if QV_COL in df.columns:
                mask &= df[QV_COL] >= self.min_qv
            if "is_gene" in df.columns:
                mask &= df["is_gene"].astype(bool)
            sub = df[mask]
            if len(sub):
                results.append(sub[[X_COL, Y_COL] + ([QV_COL] if QV_COL in sub.columns else [])])
        if results:
            return pd.concat(results, ignore_index=True)
        return pd.DataFrame(columns=[X_COL, Y_COL])

    def get_points_array(self, gene_name: str) -> np.ndarray:
        """
        Return transcript coordinates as numpy array suitable for napari Points layer.

        Transcript coordinates in the parquet/feather files are in microns.
        Images and labels are in pixel space. We convert microns → pixels
        by dividing by PIXEL_SIZE (0.2125 µm/px).

        Returns
        -------
        np.ndarray, shape (N, 2), columns = [y, x]  (napari uses row, col ordering)
        """
        df = self.load_gene(gene_name)
        if df.empty:
            return np.empty((0, 2), dtype=np.float32)
        # Convert microns → pixels and swap to (row, col) = (y, x) for napari
        return np.column_stack([
            df[Y_COL].values.astype(np.float32) / self.pixel_size,
            df[X_COL].values.astype(np.float32) / self.pixel_size,
        ])

    def get_multi_gene_points(
        self,
        gene_names: list[str],
        palette: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Load transcripts for multiple genes and assign per-point colors.

        Parameters
        ----------
        gene_names : list[str]
            Up to 10 gene names to overlay.
        palette : np.ndarray, optional
            Shape (K, 4) RGBA palette. Defaults to TRANSCRIPT_PALETTE from coloring.py.

        Returns
        -------
        points : np.ndarray, shape (N_total, 2)
        colors : np.ndarray, shape (N_total, 4)
        """
        if palette is None:
            from palms.utils.coloring import TRANSCRIPT_PALETTE
            palette = TRANSCRIPT_PALETTE

        all_points = []
        all_colors = []
        for i, gene in enumerate(gene_names):
            pts = self.get_points_array(gene)
            if len(pts) == 0:
                continue
            all_points.append(pts)
            color = palette[i % len(palette)]
            all_colors.append(np.broadcast_to(color, (len(pts), 4)).copy())

        if all_points:
            return np.concatenate(all_points), np.concatenate(all_colors)
        return np.empty((0, 2), dtype=np.float32), np.empty((0, 4), dtype=np.float32)


# ── Density preview ──────────────────────────────────────────────────────────
# Display only. The recorded density step reads ``sdata.points['transcripts']``
# through the ``transcripts.gene`` template, and must go on doing so — a
# notebook replaying from raw Xenium output has no feather index. What follows
# produces the *same rows* by the fast route, for a picture the user is only
# looking at while they hunt for a gene. See tab_transcripts.py.


def cached_columns(loader: "TranscriptLoader", gene: str) -> list[str]:
    """Column names in *gene*'s feather file, without decoding any data.

    The schema is the thing to check, not whether the directory exists: caches
    built before cell ids were kept are still on disk and still load.
    """
    import pyarrow as pa

    path = loader.cache_dir / f"{gene}.feather"
    with pa.OSFile(str(path), "rb") as handle:
        return list(pa.ipc.open_file(handle).schema.names)


def points_for_preview(loader: "TranscriptLoader", sdata, gene: str,
                       min_qv: int, need_cell_id: bool):
    """The rows ``transcripts.gene`` would produce, from the feather index.

    Returns ``(frame, "")`` on success or ``(None, reason)`` on refusal, where
    *reason* is shown to the user in place of a picture. **It never falls back
    to the parquet scan**: that path takes ~22 s and drops ``cell_id``
    entirely, so it is neither a preview nor equivalent.

    The micron-to-pixel step is deliberately not arithmetic. The frame is
    parsed as a ``PointsModel`` carrying the transcripts element's *own
    declared* transformation and handed to ``sd.transform`` — which is what
    ``transcripts.gene`` does, so the scale comes off the element and
    ``pixel_size`` never enters. ``get_points_array``'s ``/ pixel_size`` is not
    reused for the same reason, and because it also swaps to napari's
    ``(row, col)`` order, which is wrong for a points frame.

    *need_cell_id* comes from the rendered source the caller is about to
    execute, not from a guess about which blocks are selected: only the cluster
    filter reads ``cell_id``, and refusing every preview on a cache that lacks
    a column the code never touches would be a bad trade.
    """
    import spatialdata as sd
    from spatialdata.transformations import get_transformation

    from palms.preprocess import MIN_QV as CACHE_MIN_QV, cache_is_valid

    if not cache_is_valid(loader.cache_dir, loader.parquet_path):
        # The loader itself globs *.feather and never checks this, which is
        # fine for the overlay but not here: a cache built from an older
        # parquet would preview rows the recorded step would not reproduce,
        # and that is the one way a preview could be a lie about the result.
        return None, ("no preview — the transcript index was not built from "
                      "this dataset's transcripts.parquet")
    if min_qv < CACHE_MIN_QV:
        return None, (f"no preview below Min QV {CACHE_MIN_QV} — the transcript "
                      f"index is filtered at {CACHE_MIN_QV}")
    if gene not in loader.cached_genes:
        return None, f"no preview — {gene} is not in the transcript index"
    element = getattr(sdata, "points", {}).get("transcripts")
    if element is None:
        return None, "no preview — no transcripts element to take the frame from"

    wanted = {X_COL: "x", Y_COL: "y"}
    if need_cell_id:
        wanted[CELL_COL] = "cell_id"
    try:
        available = cached_columns(loader, gene)
    except Exception as exc:                      # unreadable file, mid-rebuild
        return None, f"no preview — the transcript index could not be read ({exc})"
    missing = [c for c in wanted if c not in available]
    if missing:
        return None, (f"no preview — this transcript index has no {missing}; "
                      f"rebuild it with palms-preprocess")

    df = loader._load_from_feather(gene)
    if min_qv > CACHE_MIN_QV:
        if QV_COL not in df.columns:
            return None, "no preview — this transcript index has no qv column"
        df = df[df[QV_COL] >= min_qv]
    frame = df[list(wanted)].rename(columns=wanted).reset_index(drop=True)

    points = sd.transform(
        sd.models.PointsModel.parse(
            frame, transformations=get_transformation(element, get_all=True)),
        to_coordinate_system="global",
    ).compute()
    return points, ""
