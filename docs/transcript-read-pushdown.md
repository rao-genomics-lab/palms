# Sorting the cached transcripts by gene — measured, and rejected

> **Partly superseded — read §6 and "Correction" before acting on §2, §3 or §5.**
> Sorting is a no-op *for pyarrow's own pruning*, which is what §2 measured and is
> still true. It is **5.8×** when the pruning is done in Python and the row groups
> are made finer. The "floor is decoding `feature_name`" claim in §5 is wrong. And
> the reader fix in §3/§5 **shipped and was then reverted** — §6 says why.

Status: **the sort is rejected by measurement (2026-09-01)**, beads `xv-qiq.7`. The reader
fix it was standing in for shipped the same day (PR #68) and was **reverted on 2026-09-02**
(`xv-wuk`), because it left an object stock spatialdata could not introspect and that broke
every session write on a cached launch — §6. So the transcript read is back on stock
spatialdata, and §§1–4 are again a description of an *unsolved* problem rather than a
solved one. Kept as the record of why the sort was dropped, because it is an obvious idea
that will occur to the next reader, and the metadata number that motivated it is genuinely
misleading.

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

*(Superseded in part — see Correction. Everything measured here holds; what it does
not cover is pruning done by the caller rather than by pyarrow.)*

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

All of it is in the reader, and it needs both fixes together. *(This is the measurement
that motivated the patch reverted in §6. It stands as a measurement — the pyarrow
filesystem really is 4.1× here — but PALMS no longer takes it.)* Measured through
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

So the fix belongs upstream: one `filesystem="arrow"` in `spatialdata._io.io_points`. That
is filed as beads `xv-5n7`; §5 is what PALMS does until it lands.

## 5. What shipped, and what was then withdrawn

Both halves shipped together in PR #68, because neither works alone. **The first half was
reverted on 2026-09-02** — see §6. What remains is `transcripts.gene.tmpl` saying
`(_transcripts['is_gene'] == True)`, which is correct either way and free.

The withdrawn half was `utils/arrow_points_read.py`: a context manager applied around the
`read_zarr` in `loader._open_cache`, rebinding `spatialdata._io.io_points.read_parquet` to
pass `filesystem="arrow"`. It was scoped rather than process-wide because the pyarrow
reader returns identical data (same rows, dtypes, partition count, transformations) in a
different partition **order**, and keeps `__null_dask_index__` as the index name. Nothing
in PALMS depends on points row order — `crop_export` does its own filtered read and then
`reset_index(drop=True)`, the density step histograms, the overlay uses the feather index —
but a caller's own spatialdata objects should not inherit ours. Remote stores
(`simplecache::`, `s3://`, `https://`) were left on fsspec, which cannot parse a chained
fsspec URL.

The speed it bought was real, and is recorded here so the trade is on the record. Measured
through `loader._open_cache`, the same six-gene sequence, warm cache:

| gene | stock | patched |
|---|---|---|
| CXCR3 | 2.15 s | 0.97 s |
| BSG | 2.14 s | 0.94 s |
| CXCR3 *(again)* | 2.26 s | 0.92 s |
| NOSIP | 2.34 s | 0.95 s |
| CTLA4 | 2.31 s | 0.96 s |
| BSG *(again)* | 2.26 s | 0.91 s |
| **total** | **13.45 s** | **5.65 s** |

2.4× warm, and larger cold (5.62 s → 1.49 s on the first gene) because pushdown also
avoids reading the bytes. Row counts are identical throughout. Note the shape change: on
the stock reader a 1,645-row gene costs the same as a 2.5 M-row one, because the cost is
scanning `feature_name`; patched, the cost tracks what actually matches.

This does **not** restore the ~100 ms the per-gene feather index gave. *(The sentence
that stood here — "and nothing on this route will — the floor is decoding
`feature_name`" — is wrong; see Correction. The floor on this route is decoding the
**payload** columns of every scanned row group, and it falls to 0.164 s once the
caller prunes row groups.)* `palms-preprocess` is not retired, and the feather index
stays the fast path for the point overlay.

## 6. Why it was reverted (2026-09-02)

**The patch did not just change the speed of the read. It changed the object.**
`ReadParquetPyarrowFS` builds its dask graph out of `Task` objects, and spatialdata's
`_search_for_backing_files_recursively` tests `if "piece" in v.args[0]` — so every walk of
a points element's graph raised `TypeError: argument of type 'Task' is not iterable`
(spatialdata 0.8.0 / dask 2026.7.1). Isolated by reading one real cache both ways in a
single process:

| call | stock reader | pyarrow reader |
|---|---|---|
| `repr(sdata)` | OK | **TypeError** |
| `sdata.is_self_contained()` | OK | **TypeError** |
| `_backed_elements_contained_in_path(...)` | OK | **TypeError** |
| `points.head()`, `len(points)` | OK | OK |

Data access was never affected. Two consequences were, and the second is what settled it:

1. `loader.load_sdata` printed the summary *inside* the `try` whose `except` means "the
   cache could not be opened", so a healthy store was condemned and a rebuild offered on
   **every launch after the one that built the cache**, with `cache_repair.verify` printing
   `✓ Cache is healthy` two lines later. That is a defect on its own terms and is fixed
   separately — the summary now lives outside that block, because *an unrenderable summary
   is not an unopenable store*.
2. `zarr_safe._assert_not_dask_backed` calls `_backed_elements_contained_in_path(path=live,
   object=sdata)` on the live object, and catches only `ImportError`. So **every
   `safe_write_element` on a cached launch raised** — ROIs, annotations, table persistence,
   H&E and ARMS registration, i.e. all session persistence. Confirmed directly against the
   pancreas cache. (The staging write itself was fine: it builds a fresh `SpatialData`
   holding only the element being written, so upstream never walks the points graph. The
   break was the guard, not the write.)

A fix was available and would have kept the speed: the true backing paths *are* recoverable
from the pyarrow graph, at `task.args[0].kwargs["fragment_wrapper"].fragment.path`
(verified), so a second patch of `_search_for_backing_files_recursively` could have made
the walk return the same answer as the stock reader. It was rejected as too much
private-API surface for the gain — two monkeypatches into `spatialdata`'s private modules,
one of them necessarily process-wide, to speed up a read that has a supported route
(§4/`xv-5n7`).

**The generalisable lesson**, and the reason this section exists: the patch shipped with
tripwires that checked whether *upstream's source had moved*. That is not the same question
as whether *the patched result still works with upstream's own APIs*, and no tripwire of the
first kind can catch a defect of the second. A monkeypatch that changes a returned object
needs a round-trip test — take the patched object and call the library's public
introspection (`repr`, `is_self_contained`, a write) — or it ships with an unmeasured blast
radius. Here that radius was every session write.

## Correction (2026-09-01, later the same day)

**"Sorting does nothing" is right about pyarrow and wrong as a conclusion.** §2 measured
whether *pyarrow prunes row groups by itself* on a dictionary column. It does not, and that
part stands. What it never tested is doing the pruning **in Python** — read each file's
footer statistics, decide which row groups can hold the gene, and hand pyarrow an explicit
`row_groups=` list via `ParquetFileFormat.make_fragment(...)`. The statistics are present
and correct; pyarrow simply declines to act on them, so a caller can.

Re-measured on the same 41-part cache, `NGFR` (40,226 matching rows, a median-abundance
gene), warm, best of 3:

| store | row groups | method | wall clock |
|---|---|---|---|
| as built | 123 | dask + arrow FS (the reverted patch) | 0.944 s |
| as built | 123 | Python-pruned | 1.202 s — nothing to prune |
| per-part sorted, **row groups as built** | 123 | Python-pruned | 0.724 s |
| per-part sorted, **100 k-row row groups** | 1309 | Python-pruned | **0.164 s** |
| per-part sorted, 100 k-row row groups | 1309 | stock pyarrow, no pruning | 0.830 s |

**All three changes are required, and any two of them buy almost nothing.** Sorting at the
cache's existing granularity is 1.3×, which is why §2's conclusion looked right — three row
groups of ~1.05 M rows per part is far too coarse for a 541-value column, so even a
perfectly sorted part cannot exclude much. Fine row groups without pruning are 1.1×.
Pruning without sorting is a regression. Together they are **5.8×**, and they put a gene
query at 0.164 s — within ~1.6× of the per-gene feather index's ~100 ms, from a plain
parquet layout a replayed notebook can read.

**The "floor is decoding `feature_name`" claim above is also wrong**, and wrong in the
direction that matters. Splitting the cost by projection on the store as built:

```
NGFR, columns=[feature_name]     0.349 s
NGFR, columns=[x, y, cell_id]    1.190 s
```

pyarrow decodes every *requested* column of every *scanned* row group and filters
afterwards, so the payload dominates and its cost tracks rows **scanned**, not rows
matched — 40,226 matches cost the same as 6.4 M. That is precisely why touching fewer row
groups is the lever, and why 0.164 s is reachable when §2 concluded ~1 s was a floor.

Two smaller results from the same re-measurement, both worth not re-deriving:

- **dask is *faster* than stock pyarrow here** (0.944 s against 1.190 s on the store as
  built), the opposite of what a synthetic reproduction of this problem showed. dask
  parallelises across the 41 parts; a single `dataset.to_table()` does not. Do not assume
  "go straight to pyarrow" is a win on this store.
- **An absent gene becomes nearly free with pruning on any layout** — 0.016–0.049 s against
  ~0.9 s — because a name outside every row group's min/max is excluded even when the data
  is unsorted. That is the cheapest available improvement and needs no re-layout.

**What is not settled** is whether PALMS should act on this. The pruning reader is bespoke
Python, and CLAUDE.md rule (e) says a template must prefer a library API to a hand-rolled
equivalent — a recorded cell that hand-picks row groups from footer statistics is exactly
the kind of code that rule exists to keep out of the notebook. The layout that gets the
same win through a *stock* call is hive partitioning by `feature_name`, which is what
`palms-preprocess` already does in feather form, but it changes what
`sdata.points['transcripts']` is on disk. Tracked in beads; no change has been made.

## Reproducing

`pyarrow.ParquetFile(...).metadata.row_group(i).column(j).statistics` for the surviving-row-
group counts; a gene name outside every min/max range (`ZZZZZZZ`) is the one-line test for
whether a reader prunes at all. Per-part sorting of the full cache took 105 s at 1.4 GB peak
RSS — bounded per part, as `sort_values` over ~128.7 M rows would not have been.
