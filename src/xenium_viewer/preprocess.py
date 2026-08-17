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
import argparse
import resource
from contextlib import contextmanager
from pathlib import Path
from collections import defaultdict

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Columns to keep in the feather files
# Xenium 3.x uses "feature_name" for the gene name column
GENE_COL = "feature_name"
KEEP_COLS = ["x_location", "y_location", "qv", GENE_COL, "cell_id"]

# Minimum quality value threshold
MIN_QV = 20
CHUNK_SIZE = 1_000_000  # rows per batch
FLUSH_THRESHOLD = 10_000_000  # flush buffered rows to disk once they reach this
SENTINEL_FILENAME = ".complete"

# Feather files are written to "<gene>.feather.partial" and renamed only once every
# writer has been closed. An Arrow IPC file has no footer until it is closed, so an
# interrupted run would otherwise leave unreadable "<gene>.feather" files behind —
# and TranscriptLoader._cached_genes globs "*.feather" without consulting the
# .complete sentinel, so a half-built cache is visible to the viewer.
PARTIAL_SUFFIX = ".partial"

# Genes that cannot hold an open writer (see _writer_budget) are written one
# "<gene>.part<K>.feather" shard per flush instead, then merged at the end.
SHARD_GLOB = "*.part*.feather"

# How many file descriptors to leave for everything else in the process. This runs
# inside the GUI, which already holds the zarr store, dask, Qt and the parquet reader
# open, so the writer pool must not claim the whole allowance.
_FD_RESERVE = 256
# Enough for any realistic panel (this one has 5101 genes) without asking for the
# entire hard limit.
_FD_TARGET = 16384


@contextmanager
def _raised_fd_limit(target: int = _FD_TARGET):
    """Temporarily raise this process's RLIMIT_NOFILE soft limit; yield what we got.

    A desktop-launched viewer inherits systemd's default soft limit of **1024** while
    the hard limit is ~1e6, so a process can raise its own soft limit unprivileged.
    Getting this wrong is what broke the previous version of this module: the soft
    limit was checked in an interactive shell (where conda's init had already raised
    it to ~1e6) rather than in the viewer process that actually runs the code, so a
    pool of one writer per gene died with EMFILE partway through.

    Best-effort by design — this runs inside the GUI and must never take the viewer
    down. Restoring a lower soft limit afterwards is safe: it does not invalidate
    file descriptors that are already open.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    new_soft = soft
    if soft < target:
        try:
            new_soft = min(target, hard) if hard != resource.RLIM_INFINITY else target
            resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))
        except (ValueError, OSError):
            new_soft = soft
    try:
        yield new_soft
    finally:
        if new_soft != soft:
            try:
                resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))
            except (ValueError, OSError):
                pass


def _writer_budget(fd_limit: int) -> int:
    """How many feather writers may be open at once, given the fd limit we have."""
    return max(1, fd_limit - _FD_RESERVE)


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


def preprocess(parquet_path: Path,
               cache_dir: Path,
               min_qv: int = MIN_QV,
               chunk_size: int = CHUNK_SIZE):
    cache_dir.mkdir(parents=True, exist_ok=True)

    if _cache_is_valid(cache_dir, parquet_path):
        print(f"Transcript cache is up to date ({cache_dir}). Nothing to do.")
        return

    # Stale or missing sentinel — wipe existing feather files before reprocessing
    # (including any *.feather.partial left behind by an interrupted run)
    for f in cache_dir.glob("*.feather"):
        f.unlink()
    for f in cache_dir.glob(f"*.feather{PARTIAL_SUFFIX}"):
        f.unlink()
    for f in cache_dir.glob(SHARD_GLOB):
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

    parquet_file = pq.ParquetFile(parquet_path)
    total_rows = parquet_file.metadata.num_rows
    processed = 0
    kept = 0

    # Read only the columns we actually use. Xenium's transcripts.parquet carries ~13
    # columns; without this every one of them is decoded and converted to pandas on
    # every batch just to be dropped a few lines later. Intersected with the file's
    # real schema so a dataset missing a column still works (the guards below already
    # handle absent 'is_gene'/'qv').
    present = set(parquet_file.schema_arrow.names)
    read_cols = [c for c in [*KEEP_COLS, "is_gene"] if c in present]

    # One open Arrow IPC writer per gene, created on that gene's first flush. This is
    # what replaces the old read-modify-write flush: feather has no append, so the
    # previous code read back and rewrote each gene's whole accumulated file on every
    # flush cycle — with ~13 cycles that is ~13x the necessary I/O (and the same
    # multiplier on the pd.concat cost). Here each row is written exactly once.
    #
    # The pool is capped at `budget` (see _writer_budget): a gene that arrives once the
    # pool is full is written as one shard per flush instead, costing a single fd at a
    # time, and its shards are merged in the consolidation step below. Both modes end
    # up as exactly one "<gene>.feather", so readers see no difference.
    # gene -> (writer, sink). The sink is kept and closed explicitly: pyarrow's
    # RecordBatchFileWriter.close() finalizes the IPC stream but does **not** release
    # the underlying file descriptor when pyarrow opened it from a path — the fd only
    # goes away when the object is garbage collected. Measured: 200 writers, all
    # closed, still held 200 fds until gc.collect(). Relying on refcounting here would
    # reintroduce EMFILE non-deterministically, which is precisely the failure mode
    # this pool is meant to make impossible.
    writers: dict[str, tuple] = {}
    sharded: dict[str, int] = {}   # gene -> number of shards written so far
    # All batches for a gene must share one schema — feature_name comes back from
    # to_pandas() as a *dictionary* type whose categories vary batch to batch, which
    # the writer rejects on the first mismatch. Decode dictionary columns to their
    # value type; every other column keeps the parquet's own type, so the cached
    # dtypes match what the previous implementation produced.
    source_schema = parquet_file.schema_arrow

    def _stable_type(col: str) -> pa.DataType:
        t = source_schema.field(col).type
        return t.value_type if pa.types.is_dictionary(t) else t

    schema = pa.schema([(c, _stable_type(c)) for c in KEEP_COLS if c in present])

    with _raised_fd_limit() as fd_limit:
        budget = _writer_budget(fd_limit)
        try:
            for batch in parquet_file.iter_batches(batch_size=chunk_size,
                                                   columns=read_cols):
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

                n_seen = len(writers) + len(sharded) or len(buffers)
                pct = 100 * processed / total_rows
                print(f"  {processed:>12,} / {total_rows:,}  ({pct:.1f}%)  kept={kept:,}  "
                      f"genes={n_seen}  buffer_rows={buffer_rows:,}",
                      end="\r")

                # Periodically flush large buffers to avoid memory blow-up
                if buffer_rows >= FLUSH_THRESHOLD:
                    _flush_buffers(buffers, cache_dir, writers, sharded, schema, budget)
                    buffers.clear()
                    buffer_rows = 0

            # Final flush
            _flush_buffers(buffers, cache_dir, writers, sharded, schema, budget)
            buffers.clear()
        finally:
            # Close every writer *and its sink* even on failure, so we never leak
            # handles (see the note on the writers dict above).
            for w, sink in writers.values():
                for closeable in (w, sink):
                    try:
                        closeable.close()
                    except Exception:
                        pass

        # Only now that every file has its IPC footer, publish under the real names.
        n_genes = 0
        for gene in writers:
            partial = cache_dir / f"{gene}.feather{PARTIAL_SUFFIX}"
            if partial.exists():
                partial.replace(cache_dir / f"{gene}.feather")
                n_genes += 1

        # Merge the shards of any gene that overflowed the writer pool.
        if sharded:
            print(f"\n  merging shards for {len(sharded)} gene(s) "
                  f"(writer pool capped at {budget})...")
            for gene, n_shards in sharded.items():
                paths = [cache_dir / f"{gene}.part{k}.feather" for k in range(n_shards)]
                paths = [p for p in paths if p.exists()]
                if not paths:
                    continue
                merged = pd.concat([pd.read_feather(p) for p in paths],
                                   ignore_index=True)
                merged.to_feather(cache_dir / f"{gene}.feather")
                for p in paths:
                    p.unlink()
                n_genes += 1

    _write_sentinel(cache_dir, parquet_path)

    print(f"\n\nDone. Processed {processed:,} rows, kept {kept:,} gene transcripts.")
    print(f"Wrote {n_genes} feather files to {cache_dir}")


def _flush_buffers(buffers: dict, cache_dir: Path, writers: dict, sharded: dict,
                   schema: "pa.Schema", budget: int):
    """Write each gene's accumulated rows out, without re-reading what is on disk.

    A gene normally keeps one Arrow IPC writer open for the whole run, so its rows are
    written exactly once. Once ``budget`` writers are open, further genes fall back to
    one shard per flush (opened and closed here, so one fd at a time); the caller
    merges those shards at the end. Whichever mode a gene starts in, it stays in.
    """
    for gene, dfs in buffers.items():
        if not dfs:
            continue
        table = pa.Table.from_pandas(pd.concat(dfs, ignore_index=True),
                                     schema=schema, preserve_index=False)

        entry = writers.get(gene)
        if entry is None and gene not in sharded:
            if len(writers) < budget:
                sink = pa.OSFile(
                    str(cache_dir / f"{gene}.feather{PARTIAL_SUFFIX}"), "wb")
                entry = (pa.ipc.new_file(sink, schema), sink)
                writers[gene] = entry
            else:
                sharded[gene] = 0

        if entry is not None:
            entry[0].write_table(table)
        else:
            k = sharded[gene]
            sink = pa.OSFile(str(cache_dir / f"{gene}.part{k}.feather"), "wb")
            try:
                w = pa.ipc.new_file(sink, schema)
                w.write_table(table)
                w.close()
            finally:
                sink.close()
            sharded[gene] = k + 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir",
                        help="Path to Xenium output directory (containing transcripts.parquet)")
    parser.add_argument("--min-qv", type=int, default=MIN_QV,
                        help=f"Minimum quality value (default: {MIN_QV})")
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE,
                        help=f"Rows per batch (default: {CHUNK_SIZE:,})")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    parquet_path = data_dir / "transcripts.parquet"
    cache_dir = data_dir / "transcript_cache"

    if not parquet_path.exists():
        sys.exit(f"Error: {parquet_path} not found. Is this a Xenium output directory?")

    preprocess(parquet_path=parquet_path, cache_dir=cache_dir,
               min_qv=args.min_qv, chunk_size=args.chunk_size)


if __name__ == "__main__":
    main()
