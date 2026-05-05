"""
One-time transcript preprocessing: split transcripts.parquet into per-gene feather files.

Run this script ONCE before launching the viewer. It produces
scripts/transcript_cache/{gene_name}.feather — one file per gene (~480 files).

This enables <100 ms per-gene transcript loading in the viewer versus
the 4–5 s needed to scan the full 1.3 GB parquet each time.

Usage:
    conda activate xenium_viewer
    python scripts/00_preprocess_transcripts.py

Estimated runtime: 30–60 min (I/O bound on the 1.3 GB parquet).
Output size: ~1 GB total in transcript_cache/
"""

import sys
import json
from pathlib import Path
from collections import defaultdict

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

DATA_PATH = Path(__file__).parent.parent
PARQUET_PATH = DATA_PATH / "transcripts.parquet"
CACHE_DIR = DATA_PATH / "transcript_cache"

# Columns to keep in the feather files
# Xenium 3.x uses "feature_name" for the gene name column
GENE_COL = "feature_name"
KEEP_COLS = ["x_location", "y_location", "qv", GENE_COL, "cell_id"]

# Minimum quality value threshold
MIN_QV = 20
CHUNK_SIZE = 1_000_000  # rows per batch
SENTINEL_FILENAME = ".complete"


def _write_sentinel(cache_dir: Path, parquet_path: Path):
    stat = parquet_path.stat()
    sentinel = cache_dir / SENTINEL_FILENAME
    sentinel.write_text(json.dumps({
        "mtime": stat.st_mtime,
        "size": stat.st_size,
    }))


def _cache_is_valid(cache_dir: Path, parquet_path: Path) -> bool:
    sentinel = cache_dir / SENTINEL_FILENAME
    if not sentinel.exists():
        return False
    try:
        info = json.loads(sentinel.read_text())
        stat = parquet_path.stat()
        return info["mtime"] == stat.st_mtime and info["size"] == stat.st_size
    except Exception:
        return False


def preprocess(parquet_path: Path = PARQUET_PATH,
               cache_dir: Path = CACHE_DIR,
               min_qv: int = MIN_QV,
               chunk_size: int = CHUNK_SIZE):
    cache_dir.mkdir(parents=True, exist_ok=True)

    if _cache_is_valid(cache_dir, parquet_path):
        print(f"Transcript cache is up to date ({cache_dir}). Nothing to do.")
        return

    # Stale or missing sentinel — wipe existing feather files before reprocessing
    for f in cache_dir.glob("*.feather"):
        f.unlink()
    if (cache_dir / SENTINEL_FILENAME).exists():
        (cache_dir / SENTINEL_FILENAME).unlink()
    print("Rebuilding transcript cache...")

    print(f"Reading {parquet_path} in chunks of {chunk_size:,} rows ...")
    print(f"Filtering: is_gene=True, qv >= {min_qv}")
    print(f"Writing feather files to {cache_dir}")
    print()

    # Buffer: gene_name -> list of DataFrames
    buffers: dict[str, list[pd.DataFrame]] = defaultdict(list)
    buffer_rows = 0
    FLUSH_THRESHOLD = 10_000_000  # flush when buffer has this many rows

    parquet_file = pq.ParquetFile(parquet_path)
    total_rows = parquet_file.metadata.num_rows
    processed = 0
    kept = 0

    for batch in parquet_file.iter_batches(batch_size=chunk_size):
        df = batch.to_pandas()
        processed += len(df)

        # Filter: only gene transcripts with sufficient quality
        if "is_gene" in df.columns:
            df = df[df["is_gene"].astype(bool)]
        if "qv" in df.columns:
            df = df[df["qv"] >= min_qv]

        # Keep only needed columns (subset to available)
        available_cols = [c for c in KEEP_COLS if c in df.columns]
        df = df[available_cols].copy()

        kept += len(df)

        # Accumulate per gene
        for gene, group in df.groupby(GENE_COL, sort=False):
            buffers[gene].append(group)
        buffer_rows += len(df)

        pct = 100 * processed / total_rows
        print(f"  {processed:>12,} / {total_rows:,}  ({pct:.1f}%)  kept={kept:,}  "
              f"genes={len(buffers)}  buffer_rows={buffer_rows:,}",
              end="\r")

        # Periodically flush large buffers to avoid memory blow-up
        if buffer_rows >= FLUSH_THRESHOLD:
            _flush_buffers(buffers, cache_dir)
            buffers.clear()
            buffer_rows = 0

    # Final flush
    _flush_buffers(buffers, cache_dir)
    buffers.clear()

    _write_sentinel(cache_dir, parquet_path)

    print(f"\n\nDone. Processed {processed:,} rows, kept {kept:,} gene transcripts.")
    feather_files = list(cache_dir.glob("*.feather"))
    print(f"Wrote {len(feather_files)} feather files to {cache_dir}")


def _flush_buffers(buffers: dict, cache_dir: Path):
    """Merge accumulated per-gene DataFrames and write/append to feather files."""
    for gene, dfs in buffers.items():
        if not dfs:
            continue
        merged = pd.concat(dfs, ignore_index=True)
        feather_path = cache_dir / f"{gene}.feather"
        if feather_path.exists():
            # Append to existing file by reading and re-writing
            existing_df = pd.read_feather(feather_path)
            merged = pd.concat([existing_df, merged], ignore_index=True)
        merged.to_feather(feather_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", nargs="?", default=None,
                        help="Path to Xenium output directory (default: parent of scripts/)")
    parser.add_argument("--min-qv", type=int, default=MIN_QV,
                        help=f"Minimum quality value (default: {MIN_QV})")
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE,
                        help=f"Rows per batch (default: {CHUNK_SIZE:,})")
    args = parser.parse_args()

    data_dir = Path(args.data_dir) if args.data_dir else DATA_PATH
    parquet_path = data_dir / "transcripts.parquet"
    cache_dir = data_dir / "transcript_cache"

    preprocess(parquet_path=parquet_path, cache_dir=cache_dir,
               min_qv=args.min_qv, chunk_size=args.chunk_size)
