# Sorting the cached transcripts by gene — measured, and rejected

Status: **rejected by measurement (2026-09-01)**, beads `xv-qiq.7`. Kept as the record of
why, because the idea is an obvious one that will occur to the next reader, and the
metadata number that motivated it is genuinely misleading.

## The proposal

After PR #67 the density heatmap reads `sdata.points['transcripts']` instead of the
per-gene feather index, which cost *preprocessed* datasets their fast path (~100 ms → 6.5 s;
un-preprocessed datasets got 3.4× faster). The cache stores transcripts in acquisition
order across 41 parquet parts, so **every** row group's `feature_name` statistics span
`A2M..ZNF683` and can never exclude a gene. Sorting each part by `feature_name` should let
the reader skip most row groups on statistics alone.

The staging said: settle the pushdown question first, and if the filter cannot be stated
cleanly, reconsider rather than push through. Both halves of that came back negative.

## 1. dask does not push the mask down — two independent reasons

The template's mask is `points[(feature_name == g) & (qv >= q) & is_gene]`. Neither the
optimized expression nor the lowered one carries a `filters=` operand:

```
Filter
  Projection
    ReadParquetFSSpec  filters=None
```

**(a) The reader.** Predicate pushdown is `_filter_passthrough`, and in dask 2026.7.1 it is
set on exactly one class:

| class | `_filter_passthrough` |
|---|---|
| `ReadParquet` (base) | `False` |
| `ReadParquetPyarrowFS` | **`True`** |
| `ReadParquetFSSpec` | `False` |

`spatialdata._io.io_points` calls `read_parquet(path)` with no `filesystem=`, and
`read_parquet`'s default is `"fsspec"` — so the collection the template is handed can never
push a filter down, whatever the store looks like. There is no dask config knob for it;
`filesystem` is a plain function default.

**(b) The predicate shape.** Even on the pyarrow reader, a *bare* boolean column kills the
whole conjunction — `_DNF.extract_pq_filters` returns `None` for it, and one `None` term
collapses the `And`:

| predicate | pushed-down `filters` |
|---|---|
| `feature_name == g` | `[('feature_name', '==', g)]` |
| `... & (qv >= 20)` | `[('qv', '>=', 20), ('feature_name', '==', g)]` |
| `... & is_gene` | **`None`** |
| `... & (is_gene == True)` | `[('qv','>=',20), ('is_gene','==',True), ('feature_name','==',g)]` |

So the template's `& _transcripts['is_gene']` would have to become
`& (_transcripts['is_gene'] == True)` — a one-token change, and worth nothing on its own
(see §4).

## 2. Sorting does nothing, because pyarrow does not prune on this column at all

Measured on the real 41-part cache (128.7 M rows, 2.2 GB), against a per-part sorted copy of
it. "Row groups surviving" is the metadata calculation the original proposal rested on:

| store | row groups surviving statistics | wall clock |
|---|---|---|
| as built (`NOSIP`) | 123/123 | 1.16 s |
| sorted (`NOSIP`) | **47/123** | **1.16 s** |
| as built (`BSG`) | 123/123 | 1.18 s |
| sorted (`BSG`) | 41/123 | 1.14 s |

A 2.6–3× structural reduction bought **zero** time. The reason is not that the saving is
small — it is that no skipping happens. A gene name excluded by *every* row group's min/max
costs full price:

```
unsorted  ZZZZZZZ  1.21s        0 rows   (excluded by every row group's statistics)
unsorted  CXCR3    1.15s    3,111 rows   (feature_name alone; 1,964 once qv/is_gene apply)
sorted    ZZZZZZZ  1.22s        0 rows
```

**`feature_name` is stored dictionary-encoded** (`dictionary<values=string, indices=int16>`,
`category` in pandas), and pyarrow's scanner does not apply row-group statistics to it. The
"26/32 → 1/32" figure that motivated this task was a count of row groups whose min/max
*could* exclude the gene — not evidence that the reader acts on them. It never did.

Pruning engages only if the column is a plain string. The 2×2, on 6 parts:

| layout | `ZZZZZZZ` (no match) | `CXCR3` (rare) | `BSG` (abundant) |
|---|---|---|---|
| dict / unsorted | 0.206 s | 0.230 s | 0.225 s |
| dict / **sorted** | 0.209 s | 0.208 s | 0.223 s |
| string / unsorted | **0.003 s** | 0.273 s | 0.277 s |
| string / **sorted** | **0.003 s** | **0.155 s** | **0.146 s** |

Three things to read off it. Sorting a dictionary column is a no-op. Retyping alone is a
*regression* on real genes (0.230 → 0.273 s) — the dictionary is also the compression.
Only the two together win, and only by ~1.5×, in exchange for changing the on-disk dtype
of `feature_name` for every consumer. That is a far larger change than the sort it was
supposed to justify, for less than the free win in §3.

## 3. Where the actual win is

All of it is in the reader, and it needs both fixes together. Measured through
`spatialdata`, on the **unmodified** store, with `io_points.read_parquet` monkeypatched to
pass `filesystem="arrow"`:

| | bare `is_gene` | `is_gene == True` |
|---|---|---|
| stock spatialdata (fsspec) | 5.90 s | 4.21 s |
| pyarrow filesystem | 5.44 s | **1.32 s** |

4.1× against the stock path, with no change to the cache. Neither half works alone.

## 4. Why the template cannot do this itself

The template would have to build its own `dd.read_parquet(..., filesystem="arrow")`, which
needs the store's path. `SpatialData.path` is `None` until a store is read or written, and
the recorded preamble builds `sdata = xenium(data_path)` from the **raw** output — so in a
replayed notebook there is no path to reach for. Reading `data_path/"transcripts.parquet"`
instead trades that for two worse problems: it breaks cache-only datasets (a crop export
has no raw source — `loader.has_raw_xenium_source()`), and it means hand-rolling the column
renames and the micron→pixel transform that `spatialdata_io` and the element's own
transformation already declare, which is exactly the smell rule (e) in CLAUDE.md names.

So the fix belongs upstream: one `filesystem="arrow"` in `spatialdata._io.io_points`, after
which the template's one-token predicate change is worth 4×. Until then the feather index
stays the fast path for preprocessed datasets, and **`palms-preprocess` cannot be retired** —
the possible outcome the original task floated is off the table on this route.

## Reproducing

`pyarrow.ParquetFile(...).metadata.row_group(i).column(j).statistics` for the surviving-row-
group counts; a gene name outside every min/max range (`ZZZZZZZ`) is the one-line test for
whether a reader prunes at all. Per-part sorting of the full cache took 105 s at 1.4 GB peak
RSS — bounded per part, as `sort_values` over ~128.7 M rows would not have been.
