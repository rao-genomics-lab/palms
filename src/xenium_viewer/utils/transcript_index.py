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
            from xenium_viewer.utils.coloring import TRANSCRIPT_PALETTE
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
