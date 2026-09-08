# Test data

## `crop_7/` — the end-to-end fixture

6.7 MB in 378 files. Real 10x-derived data, checked in so `tests/e2e/` runs in
CI with no download and no network. `tests/e2e/conftest.py` stages a copy into
`tmp_path` for every test and never launches the viewer against this directory:
a launch writes `palms.log`, `analysis.py`, `analysis_notebook.ipynb` and
`viewer_cache/prov_graph.json` beside the dataset, and repairs the cache.

### What it is

A rectangular crop of a real Xenium 3.1.0 run (`COMBAT_PROSTATE_CRC_04DEC24`,
Region_2, panel `hMulti_100g`), cut with the viewer's own Tools → Crop Dataset,
then used for a full analysis session. **151 cells × 5099 genes.**

It is deliberately *not* pristine. The session it carries is what makes it worth
having, and reloading it is most of what the tests assert:

- a **31-node provenance graph** with two lineages — a `qc_filter` barrier node
  and the `:qc`-scoped `normalize:qc` / `spatial_neighbors:qc` beside the
  unfiltered `normalize` / `spatial_neighbors`;
- ~14 clustering columns, named cluster labels, `rank_genes`, an inferCNV
  result, `X_umap`;
- H&E and ARMS affines, registration landmark shapes, patch overlays;
- the four `migrated_*` markers, so the one-time migrations stay exercised.

It is a **crop export, so it is cache-only**: `has_raw_xenium_source()` is False
and `.xv_manifest.json` declares `cache_only: true`. The zarr *is* the data.
That makes it a fixture for the never-stale / never-rebuild branch, and it means
these tests do **not** cover the raw `spatialdata_io.xenium` load path — the one
a real first launch takes. Nothing here can; that needs raw 10x output.

### How it was made

```bash
python scripts/prepare_e2e_fixture.py /path/to/demo_data/crop_7
```

Three things happen to the copy, and the script refuses to finish if the last
two did not take:

1. **Stripped** — `analysis.py`, `analysis_notebook.ipynb`, `plots/`,
   `palms.log` (all derived; the tests assert the launch recreates them),
   `transcript_cache/` (15 MB of per-gene feathers — without it the Transcripts
   tab takes its parquet fallback, which is a real path and costs nothing on a
   683 KB file), and the store's `.xv_trash` / `.xv_staging` / `.xv_journal` /
   `.xv.lock`.

2. **Repathed** to the placeholder `/palms-e2e-fixture/crop_7`, via the shipped
   `rename_dataset.repair`. Two reasons, and the second is not cosmetic: the
   recorded paths would otherwise publish the author's home directory, *and*
   `app.py` re-emits the preamble for the current `data_path` on every launch,
   so a graph recorded elsewhere flags all 30 descendants stale for nothing.
   Substitution is scoped to the **parent** directory, because the
   `crop_export:crop_7` node records `output_dir=Path('<parent>')` and no
   substitution scoped to the dataset itself can reach it.

3. **Scrubbed** — `arms_he_filename` replaced with `demo_he_section.svs`. A real
   scan filename is a slide identifier and it reaches the napari layer list.
   Same leak `scripts/prepare_demo_dataset.sh` exists to close.

### Using a different dataset

Point `PALMS_E2E_DATASET` at any Xenium output directory and the E2E suite runs
against that instead. Assertions that name this fixture's own numbers are
skipped; the lifecycle ones are not.
