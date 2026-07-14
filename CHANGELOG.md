# Changelog

## [Unreleased] — 2026-07-14 (d)

### Fixed
- **Datasets with no clusterings crashed at startup** — a dataset with an empty
  clustering list (e.g. a Crop Dataset export, which has no raw Xenium
  `analysis/clustering/` folder) crashed while building the control panel:
  `ValueError: None is not a valid choice. must be in ()`, raised by every
  "Clustering" magicgui `ComboBox` when handed `value=None` against empty
  choices. Added a shared `combo_value_kwargs()` helper
  (`tabs/_helpers.py`) that omits the `value` kwarg when the choice list can't
  supply the requested index (magicgui then safely defaults to the first
  choice, or `None` when empty), and applied it to every clustering selector
  (Cell Coloring, Rank Genes, Markers, Lig-Rec, Nhood Enrich, Co-occur, Annot
  Nhood, Annot Dist) plus the Gene Correlation gene selectors. Such datasets
  now open normally, with the clustering dropdowns simply empty until a
  clustering is computed.

## [Unreleased] — 2026-07-14 (c)

### Fixed
- **Crop Dataset: exported datasets failed to open with `FileNotFoundError:
  ... projection.csv`** — `load_umap()`/`load_clusterings()` unconditionally
  read `analysis/umap/.../projection.csv` and `analysis/clustering/`, which a
  Crop Dataset export never has (no raw Xenium `analysis/` folder). Both now
  tolerate a missing `analysis/` folder (return `None`/`{}` instead of
  raising). `app.py`'s dataset loader also now reconstructs the UMAP tab's
  coordinates from `adata.obsm['X_umap']` when the raw CSV is absent but the
  table already carries embedded UMAP coordinates from the source dataset —
  so a cropped dataset's UMAP view still works if the source had one computed
  at crop time. Verified against a real cropped dataset via the full
  `app._load_dataset()` reload path.

## [Unreleased] — 2026-07-14 (b)

### Fixed
- **Crop Dataset: exported datasets failed to open with `--no-cache`** —
  `load_sdata()` unconditionally rebuilt from raw Xenium files whenever
  `use_cache=False`, but a Crop Dataset export has no raw files (by design —
  it's a lightweight, self-contained zarr package), so this crashed with
  `FileNotFoundError: ... cells.zarr.zip`. Now `load_sdata()` falls back to
  the zarr cache whenever raw Xenium files aren't present, regardless of the
  `use_cache`/`--no-cache` flag, since that's the only way such a directory
  can ever be loaded. Fixes already-exported crop datasets too, not just new
  ones, since the check is based on file presence rather than anything
  written at export time.

## [Unreleased] — 2026-07-14 (a)

### Fixed
- **Crop Dataset: orphan nuclei in `nucleus_labels`** — nucleus IDs are their own
  independent numbering, unrelated to `cell_id`/`cell_labels`, so an ID-based overlap
  check (added, then removed, in the previous fix) couldn't reliably mask them and
  the fallback left `nucleus_labels` completely unmasked. Now masked correctly via
  spatial overlap with the already-masked `cell_labels` crop: a nucleus is kept only
  if it occupies at least one pixel within a kept cell's footprint.

## [Unreleased] — 2026-07-13

### Added
- **Crop Dataset (Tools tab)** — draw one or more polygons in a new "Crop Regions"
  napari layer and export each as its own standalone, independently-openable
  xenium-viewer data directory (cropped morphology image, cell/nucleus labels,
  transcripts, and AnnData table). Images/labels are cropped to each polygon's pixel
  bounding box; cells and transcripts are filtered to the exact drawn polygon via
  true point-in-polygon tests. Output/name for each region is chosen via sequential
  folder-picker + name-prompt dialogs, run in a background thread with a progress
  dialog. New module `src/xenium_viewer/utils/crop_export.py`; new tab
  `src/xenium_viewer/tabs/tab_crop_dataset.py`.

## [Unreleased] — 2026-07-01

### Fixed
- **Permission-denied dialog gaps** — two write paths were missing coverage from the
  read-only zarr dialog: (1) `delete_element_from_disk` in the Patches tab only printed
  to console on failure; (2) `record_code` in `_helpers.py` had no exception handling at
  all when writing `code.py`. Both now call `_maybe_show_permission_dialog` on
  `PermissionError`/`OSError`.

## [Unreleased] — 2026-06-30

### Added
- **tab10 palette in Patches tab** — `tab10` (matplotlib's 10-colour categorical
  palette) is now available in the Patches tab palette dropdown alongside tab20,
  glasbey_dark, Set1, Set3, and ARMS.

## [Unreleased] — 2026-06-26 (d)

### Fixed
- **Leiden clustering UI freeze** — `sc.tl.leiden` held the Python GIL during graph
  construction and partition, causing the progress bar and status-bar spinner to freeze.
  The Leiden step now runs in a subprocess via `ProcessPoolExecutor` (spawn context),
  giving it its own GIL so the main process's Qt event loop stays responsive throughout.
  New module: `src/xenium_viewer/utils/leiden_worker.py`.

## [Unreleased] — 2026-06-26 (c)

### Added
- **Progress bar for long-running analyses** — an indeterminate `QProgressBar` now
  appears inside the control panel directly below the run button while an analysis is
  in progress (Leiden clustering, rank genes, L-R, neighbourhood enrichment,
  co-occurrence, ROI DEG, annotation nhood enrichment, Novae domains). The bar
  disappears automatically when the analysis finishes or errors out. The existing
  napari status-bar spinner/tqdm text is retained unchanged.

## [Unreleased] — 2026-06-26 (b)

### Added
- **Wiki screenshots** — all 52 documentation screenshot placeholders are now filled with
  actual PNGs captured from a running viewer instance. A new script
  `scripts/capture_screenshots.py` automates future recapture by programmatically
  navigating each tab and grabbing the control-panel and full-window views.
  One placeholder (`tutorial-clustering-step5.png`, the matplotlib dotplot window)
  remains a comment for manual capture.

## [Unreleased] — 2026-06-26

### Fixed
- **Zarr cache rebuild no longer silently destroys user data** — when `experiment.xenium`
  is newer than `sdata_cached.zarr` (e.g. after Xenium Explorer opens the dataset, or a
  backup/rsync resets timestamps), a Qt dialog now appears listing all user-generated data
  found in the cache (ROIs, H&E/ARMS landmarks and images, clusterings, analysis results,
  etc.) and offers three choices: *Rebuild and restore my data* (backs up the old cache,
  rebuilds from raw files, merges user elements back, then deletes the backup),
  *Rebuild without restoring* (previous behaviour), or *Keep existing cache* (skip
  rebuild and load from the stale cache as-is — safe when only metadata changed).
  When the cache is unreadable (corrupt), it is now preserved as
  `sdata_cached_corrupt_<timestamp>.zarr` instead of being silently deleted, so data can
  be recovered manually.

## [Unreleased] — 2026-06-25

### Fixed
- **ROI DEG region grouping** — `compute_roi_deg` was reading `uns['rank_genes_groups']`
  (the scanpy default key) instead of the key it had just written (`uns['wilcoxon']` /
  `uns['t-test']`). Because the subset AnnData is a copy of `ctx.adata`, it could
  inherit stale cell-type DEG results from a previous analysis, causing "Run ROI DEG"
  to display cluster-level rather than region-level results. Fixed by passing `key=method`
  to `sc.get.rank_genes_groups_df`.
- **Permission errors on read-only datasets** — write failures (clustering import, H&E
  registration, session save, etc.) previously silenced to the console only. A
  `QMessageBox` warning dialog now fires on the first `PermissionError` / read-only
  `OSError` of a session, telling the user to copy the dataset to a writable location
  or launch with `--no-cache`.

### Added
- **ARMS Overlay: save/load landmarks** — "Save Landmarks..." and "Load Landmarks..."
  buttons added to the ARMS Overlay tab (below "Compute Registration"), mirroring the
  H&E Registration tab. Landmarks and the computed affine are saved to a portable JSON
  file via the existing `save_landmarks` / `load_landmarks` API. The save button
  enables when ≥1 landmark pair is placed and disables on "Clear All".

### Changed
- **H&E Registration: removed Save Affine button** — the "Save Affine..." button has
  been removed. There was no corresponding "Load Affine..." button, making it a dead end.
  The affine is already auto-persisted to the sdata zarr cache on every registration.

## [Unreleased] — 2026-05-05

### Changed
- **Installable Python package** — the repo is now packaged with `pyproject.toml`
  (hatchling build backend) and shipped under `src/xenium_viewer/`. The three
  numbered entry-point scripts have been renamed:
  `00_preprocess_transcripts.py` → `xenium_viewer/preprocess.py`,
  `01_load_sdata.py` → `xenium_viewer/loader.py`,
  `02_xenium_viewer.py` → `xenium_viewer/app.py`. Cross-imports were rewritten
  from `from utils.X` / `from tabs.X` to `from xenium_viewer.utils.X` /
  `from xenium_viewer.tabs.X`.
- **Console-script entry points** — `xenium-viewer`, `xenium-preprocess`,
  `xenium-fetch-references`, and `xenium-build-custom-segmentation` are
  installed on `PATH`.
- **Conda environment file** — `environment.yml` reproduces a working env in
  one step (`conda env create -f environment.yml`); pypi-only deps and the
  package itself (`-e .`) are listed in the embedded `pip:` block. Optional
  extras (`celltypist`, `r`, `gpu`, `references`, `full`) are exposed via
  `[project.optional-dependencies]`.
- **Removed importlib hacks** — `app.py` no longer loads its sibling modules
  via `importlib.util.spec_from_file_location`; they are normal package
  imports. The `sys.path.insert` at module top is gone.

## [Unreleased] — 2026-04-18

### Changed
- **External Images tab: single composite layer** — multichannel images (e.g. 25-channel PhenoCycler)
  now display as a single RGB composite napari layer instead of one layer per channel. Per-channel
  visibility checkboxes, color buttons, and contrast range sliders are in the tab widget.
  Composite is built lazily via dask so large images render efficiently.
- **External Images tab: landmark registration** — external images can now be independently aligned
  to the Xenium image via landmark-based registration (same `compute_landmark_affine` as H&E tab).
  Includes flip V/H checkboxes. The "Apply transform from" dropdown remains as an alternative.
  Landmarks are persisted to sdata.shapes; affine to sdata transformations.

## [Unreleased] — 2026-04-17

### Fixed
- **Overlay affine persistence** — affine transforms for patch overlays and external images are now
  saved to the SpatialData object via `set_transformation()` / `write_transformations()` (same
  pattern as H&E registration). On restore, the affine is read back from sdata directly, so
  overlays are correctly positioned even before the source layer (e.g. H&E) finishes loading.
  A deferred-linking listener also re-establishes live affine mirroring once the source layer
  appears.
- **QCheckBox crash on close** — `_snapshot_layers` now reads `entry["hidden_cluster_ids"]` (a
  plain set maintained by the tab) instead of iterating Qt checkbox widgets that may already be
  destroyed during shutdown.

### Changed
- **ARMS palette for subclone predictions** — subclone prediction overlays now default to the
  ARMS palette (RColorBrewer Set1+Set2+Dark2, 24 colours) to match the R-based ARMS visualisation
  pipeline. Cluster-to-colour mapping is 0-normalised so 1-based genomic cluster IDs (1,2,3) map
  to Set1 red, blue, green — matching the R package exactly.

## [Unreleased] — 2026-04-16

### Added
- **External Images tab** — load arbitrary multichannel OME-TIFF/TIFF/SVS files (e.g. PhenoCycler)
  as lazy-loaded napari layers. Channel axis is detected from `tif.series[0].axes` / OME-XML so
  IF images are never misinterpreted as RGB. Each channel gets its own napari sub-layer with
  native contrast / colormap / visibility controls; the tab adds group opacity, show/hide-all,
  and an "Apply transform from…" selector that live-mirrors another layer's affine (updates when
  the source is re-registered). Pixels and affine are written to `sdata.images["ext_<slug>"]`;
  only per-entry UI state (opacity, affine source) lives in zarr attrs.
- **Patch Overlays tab** — load phikon patch-cluster outputs (folder with `patches/coordinates.npy`
  + `clustering/cluster_labels.npy`) and subclone-prediction CSVs (`x_coord`, `y_coord`,
  `predicted_genomic_cluster`, `morphology_cluster`, `prediction_confidence`). Patches render as
  polygon rectangles in a napari Shapes layer (scales to 10k+ patches). Controls for cluster
  column, palette (tab20 / glasbey_dark / Set1 / Set3 / ARMS), outline-only, edge width, opacity,
  confidence threshold, and per-cluster visibility. Patch size is inferred from folder / file name
  with a stride-based sanity check and confirmation dialog. Geometry, cluster columns, and
  confidence are stored in `sdata.shapes["patch_<slug>"]`; UI settings persist in zarr attrs.
  Subclone prediction overlays default to the ARMS palette (RColorBrewer Set1+Set2+Dark2) to
  match the R-based ARMS visualisation pipeline.
- **Affine mirroring helper** — `utils/affine_linking.py` lists layers with non-identity affine
  and wires up `events.affine` subscriptions so external overlays stay aligned when the source
  image is re-registered.

## [Unreleased] — 2026-04-13

### Added
- **Cluster labels persisted to SpatialData** — user-assigned cluster names (from the label editor,
  reference atlas annotation, and LLM annotation) are now saved as `adata.obs["cluster_labels_<key>"]`
  columns (e.g. `cluster_labels_leiden_r1.0`) inside `sdata.tables["table"]` immediately on
  assignment. Labels are loaded back on the next launch and merged with the session-attrs fallback
  (sdata wins on conflict). The obs columns are readable in any standalone Python session.

## [Unreleased] — 2026-04-12

### Changed
- **ROI DEG and ARMS tile DEG persistence migrated to SpatialData** — results are now saved as
  sidecar parquets inside the zarr cache (`roi_deg_cache.parquet`, `arms_tile_deg_cache.parquet`)
  immediately on computation rather than at session close. Restores automatically on relaunch.
  One-time migration copies old `viewer_session/*.parquet` files to the new location.
- **ARMS tile DEG code recording fixed** — generated code snippet in `code.py` is now fully
  executable: loads tile polygons from `sdata` via `load_arms_tiles_from_sdata`, reconstructs the
  ARMS registration affine (fine × flip), applies it to the tiles, filters by selected tile
  clusters, and optionally applies a Xenium cell cluster mask before calling `compute_arms_tile_deg`.

## [Unreleased] — 2026-04-10

### Added
- **Custom segmentation SpatialData persistence** — custom segmentations are now cached inside the
  SpatialData zarr store (`sdata.labels["custom_cell_labels"]` + `sdata.tables["custom_table"]`)
  after first load, so subsequent opens do not require re-selecting the h5ad file.
  - On "Load Custom Segmentation...": if a cached version is detected in sdata, a dialog offers to
    load from cache (fast) or select a new h5ad file.
  - Auto-restore on launch: if custom segmentation was active when the dataset was closed,
    `_restore_session` reloads it from sdata automatically.
  - All per-action saves (clustering, rank genes, nhood enrichment, etc.) now route to
    `sdata.tables["custom_table"]` when custom segmentation is active, so analysis results are
    preserved across sessions.
  - `session.py` now persists `segmentation_source` ("xenium" | "custom") in the viewer session attrs.
  - **"Update SpatialData on disk"** button in the Segmentation tab force-syncs current in-memory
    state to sdata (custom table + labels in custom mode; xenium table in xenium mode).

## [Unreleased] — 2026-03-26 (2)

### Added
- **Custom cell segmentation** — two-stage pipeline to replace the native Xenium segmentation
  with a custom one from a Seurat v5 RDS file.
  - `scripts/extract_seurat_segmentation.R` — Stage 1: extracts polygon vertices, count matrix,
    and cell metadata from a Seurat RDS via slot access (no Seurat installation required). Outputs
    `segmentation_polygons.csv`, `counts.mtx`, `genes.txt`, `barcodes.txt`, `cell_metadata.csv`.
  - `scripts/build_custom_segmentation.py` — Stage 2 (run in `xenium_viewer` conda env): reads
    Stage 1 outputs, validates coordinate system against `cells.parquet`, rasterizes 300K+ polygons
    using `rasterio.features.rasterize`, writes a multi-scale zarr label store
    (`custom_labels.zarr`) and an AnnData file (`custom_segmentation.h5ad`) with
    `obs['cell_id']` = integer label and `obsm['spatial']` = xy-µm centroids, plus a metadata JSON.
  - **Segmentation tab** — new viewer tab with "Load Custom Segmentation..." and "Revert to Xenium
    Segmentation" buttons. Load flow: select `custom_segmentation.h5ad` → zarr found automatically
    via JSON metadata → cell labels layer replaced, `ctx.adata`, `ctx.color_manager`,
    `ctx.centroids_yx`, `ctx.clusterings` all rebuilt from the new data, clustering ComboBoxes
    refreshed across all tabs, stale analysis caches cleared. Revert restores the native
    spatialdata labels and reloads the original clustering files.
- `segmentation_source: str` field on `ViewerContext` (values: `"xenium"` | `"custom"`).

## [Unreleased] — 2026-04-01

### Added
- **UMAP save** (`tabs/tab_umap.py`) — "Save UMAP Plot..." button saves a scanpy-native `sc.pl.umap`
  figure in PNG or SVG. Uses the current cell-coloring clustering; cluster labels (if set) are
  applied as category names so they appear in the legend and on-data annotation.
- **Marker Genes tab** (`tabs/tab_marker_genes.py`) — new tab accepting a JSON marker-genes dict
  (`{"Cell type": ["Gene1", "Gene2"]}`). Generates dotplot, heatmap, matrixplot, tracksplot, and
  correlation matrix using scanpy's native plotting functions with the chosen clustering. Format
  (PNG/SVG) is selectable. The marker dict is persisted to the zarr session and restored on startup.
- **Volcano plot (revised)** — `make_volcano_plot()` redesigned with EnhancedVolcano aesthetics:
  only genes passing **both** thresholds are coloured (up-regulated red `#DC0000`,
  down-regulated blue `#4DBBD5`); everything else is grey. x-axis is symmetric around zero with
  auto-computed 99th-percentile limits. Points outside the display range are shown as directional
  triangle markers (`>`, `<`, `^`) pinned to the axis edges. Default labelled genes increased from
  10 to 20. Uses seaborn `despine` for theme_classic look and `adjustText` for label repulsion
  (requires `pip install adjusttext`; falls back to plain text if not installed).
- **ROI DEG: volcano plots** — "Save Volcano Plot(s)..." button in the ROI tab generates pairwise
  volcano plots after running ROI differential expression. With 2 ROIs, 1 plot is saved; with N ROIs,
  all C(N, 2) pairwise plots are saved as `roi_volcano_Region_X_vs_Region_Y.png` (300 dpi).
  Thresholds: adjusted p-value < 0.01 and |log2FC| > 1. Significantly up-regulated genes are red,
  down-regulated are blue, with dashed threshold lines. Progress shown in status bar.
  - `compute_roi_deg()` now returns `(df, adata_norm)` tuple so the normalized subset AnnData is
    available for pairwise comparisons (matching the existing `compute_arms_tile_deg()` signature).


- **ARMS tiles: "Outline only" checkbox** — when checked, tiles are rendered as colored outlines
  (edge = cluster color, fill = transparent) instead of filled polygons. The checkbox is enabled
  once tiles are loaded and toggles the existing layer in place without reloading.
- **ARMS tiles: "Tile edge width" slider** — adjusts outline thickness live (range 1–100, default
  20). Enabled once tiles are loaded.

### Changed
- **ARMS tiles: ColorBrewer Set1+Set2 color palette** — replaced the hardcoded 8-color custom
  palette with a concatenation of ColorBrewer Set1 (9 colors) + Set2 (8 colors) = 17 distinct
  colors. Cluster filter checkboxes and the legend now show `C{id}` labels without color names.

## [Unreleased] — 2026-03-27 (3)

### Fixed
- **Annot. Nhood and Annot. Distance clustering dropdowns not refreshed on segmentation swap** —
  both tabs created their `clustering_widget` ComboBoxes without registering them on `ctx`, so
  `refresh_clustering_choices` skipped them. Fix: register as `ctx.annot_nhood_clustering_widget`
  and `ctx.annot_dist_clustering_widget`, add corresponding fields to `ViewerContext`, and include
  them in the refresh loop in `_helpers.py`.

## [Unreleased] — 2026-03-27 (2)

### Fixed
- **Custom segmentation: cluster coloring all-transparent** — obs columns (`cell_type`, etc.)
  added to `ctx.clusterings` from `new_adata.obs` were indexed by obs row numbers (`'0'`, `'1'`,
  …) rather than `cell_id` integer label values. `get_cluster_colors` calls
  `cluster_series.reindex(cell_ids)` using label values as the lookup key, so the index mismatch
  produced all-NaN alignments and every cell got RGBA `(0,0,0,0)`. Fix: wrap obs columns in a new
  `pd.Series` indexed by `adata.obs['cell_id'].values` before storing in `ctx.clusterings`.
  Gene expression coloring was unaffected because it uses `label_to_obs` directly.
- **Custom segmentation: `IndexError` in UMAP hover** — `UMAPViewer._valid` has shape
  `(n_original_cells,)` but `cluster_ids_per_obs` from a custom segmentation has a different
  length. Guard added: skip UMAP hover cluster IDs when sizes don't match.

## [Unreleased] — 2026-03-27

### Added
- **Annotation layer** — new napari Shapes layer "Annotations" for drawing named tissue regions
  (bone, adipocyte, vessel, etc.) with per-shape type labels displayed as text overlays.
  Annotations persist across sessions via `sdata.shapes['annotations']`.
- **Annotations tab** — manage annotation type labels (assign, colour-pick, delete), import from
  GeoJSON (compatible with QuPath exports), export to GeoJSON.
- **Annot. Nhood tab** — neighbourhood enrichment analysis that includes annotation regions as
  virtual cell types. Samples a configurable-density grid of virtual cells inside annotation
  polygons, builds an augmented AnnData, and runs the squidpy neighbourhood enrichment pipeline.
  Annotation types appear as rows/columns in the Z-score heatmap.
- **Annot. Distance tab** — for each real cell computes minimum distance (µm) to the boundary of
  a selected annotation type (using shapely's vectorised distance API). Visualises the distribution
  per cell-type cluster as violin, box, or strip plots. Cells can optionally be coloured on the
  canvas by their distance value using a Points layer.

## [Unreleased] — 2026-03-26

### Added
- **View menu: Show Minimap toggle** — checkable "Show Minimap" action appended to napari's native
  View menu. Enabled and checked when the minimap overlay is available; disabled (grayed out) when
  there is no morphology data. State is reset on dataset reload.

## [Unreleased] — 2026-03-25

### Changed
- **Phase 3 SpatialData storage refactoring** — H&E and ARMS landmark pairs, and ARMS tile
  geometries migrated from custom zarr arrays into native SpatialData shapes elements, making
  the sdata object self-sufficient for Python analysis.
  - **H&E landmarks** → `sdata.shapes['he_xenium_landmarks']` and `sdata.shapes['he_he_landmarks']`
    as GeoDataFrames of shapely Points (xy coords in native pixel spaces). Saved after each
    registration event; cleared on landmark clear.
  - **ARMS landmarks** → `sdata.shapes['arms_xenium_landmarks']` and `sdata.shapes['arms_he_landmarks']`
    — same pattern as H&E landmarks.
  - **ARMS tile polygons** → `sdata.shapes['arms_tiles']` as a GeoDataFrame of shapely Polygons
    with `tile_name` and `cluster_id` columns. Saved after GeoJSON+CSV load; restored at startup
    without requiring the original files.
  - UI-only state (flip toggles, coarse_affine, file paths) remains in zarr attrs as before.
  - `utils/session.py` no longer writes landmark zarr arrays; zarr load paths kept as migration
    fallback for existing datasets.
  - 4 new functions in `utils/adata_persistence.py`: `save_landmarks_to_sdata`,
    `load_landmarks_from_sdata`, `save_arms_tiles_to_sdata`, `load_arms_tiles_from_sdata`.
  - **Automatic migration**: on first launch after upgrade, `migrate_landmarks_to_sdata()`
    migrates existing landmarks into sdata.shapes. Sources checked in order: (1) zarr session
    arrays (snapshot-captured at close time), (2) `landmarks.json` in the dataset folder
    (saved via Save Landmarks button). Also re-parses GeoJSON/CSV tile files into
    sdata.shapes['arms_tiles'] if the files still exist. Marked with
    `migrated_landmarks_to_sdata` flag so it only runs once.
  - All shape GeoDataFrames now constructed via `ShapesModel.parse()` to ensure spatialdata
    compatibility (attaches coordinate transform metadata required by the zarr writer).
- **Phase 2 SpatialData storage refactoring** — rank genes results, normalized AnnData, and
  ROI polygons migrated out of custom `viewer_session/` files into native SpatialData locations.
  - **Rank genes** → `adata.uns['rank_genes_groups']` (scanpy native); DataFrame reconstructed
    via `sc.get.rank_genes_groups_df()` on restore. Removes `rank_genes.parquet`.
  - **Normalized AnnData** → `sdata.tables['adata_norm']` (zarr-backed); replaces the
    `rank_genes_adata_norm.h5ad` file. Required for dotplot and volcano plots after restore.
  - **ROI polygons** → `sdata.shapes['rois']` as a GeoDataFrame of shapely Polygons
    (xy coords); replaces per-polygon zarr arrays in `viewer_session/rois/`. Old zarr arrays
    serve as fallback for datasets not yet saved with the new code.
  - Rank genes results now enable Show/Export buttons at startup (loaded before `restore_fn`,
    same pattern as Phase 1 analyses).
  - Old `rank_genes.parquet` + `rank_genes_adata_norm.h5ad` files are auto-migrated and
    deleted on first launch with the new code.
  - `utils/session.py` no longer handles rank genes or ROI polygon storage.
- **Phase 1 SpatialData storage refactoring** — clusterings, nhood enrichment, co-occurrence,
  L-R interaction results, and UMAP coordinates are now stored in native AnnData locations
  (`adata.obs`, `adata.obsm["X_umap"]`, `adata.uns`) instead of custom zarr/parquet files in
  `viewer_session/`. Persisted via `sdata.delete_element_from_disk` + `sdata.write_element`.
  - Analysis results are now saved immediately after computation (crash-safe), rather than
    deferred to viewer exit.
  - Old `viewer_session/` data is automatically migrated on first launch with the new code.
  - New module: `utils/adata_persistence.py` centralizes all adata read/write/migration logic.
  - `utils/session.py` no longer handles clusterings, nhood, co-occurrence, or L-R data.

### Fixed
- **Index alignment for adata persistence** — `adata.obs.index` uses integer strings (`'0'`, `'1'`, ...)
  while clustering series and UMAP are indexed by cell barcode (`'aaaagflk-1'`, ...). Save/load now
  maps between the two via the `cell_id` column.
- **Show/Export buttons disabled at startup for persisted analyses** — nhood enrichment,
  co-occurrence, and L-R results were loaded from `adata.uns` *after* tab session restore ran,
  so the buttons stayed disabled despite results existing. Fixed by loading `adata.uns` analysis
  results before calling `restore_fn()`.

## [Unreleased] — 2026-03-23

### Added
- **Level slider in Novae tab** — exposes the `level` parameter (hierarchical tree level,
  default 7, range 1–15) for `assign_domains()`, giving finer control over domain granularity.
- **Console variable injection** — key variables (`adata`, `sdata`, `viewer`, `ctx`,
  `clusterings`, `color_manager`, `gene_names`, `data_path`) are now pushed into napari's
  built-in IPython console on dataset load, with a help message listing what's available.
  Variables are refreshed on dataset reload. Enables a GUI-to-code handoff workflow alongside
  the existing code recording feature.

### Fixed
- **Console button hidden by layer widgets** — capped layer controls and layer list dock
  widgets to 200px max height so the console toggle button stays visible.
- **Console not resizable when opened** — the Xenium Controls dock (QTabWidget with 13 tabs)
  had a ~970px minimum height from stacked widget content, leaving no room for the console.
  Fixed by wrapping tab content in `QScrollArea` inside `make_tab()`, dropping the dock's
  minimum height to ~117px. Also defers the `resizeDocks` call via `QTimer.singleShot(0)` to
  ensure it runs after Qt finishes the visibility layout pass.

## [Unreleased] — 2026-03-20

### Added
- **Spatial Domains tab (Novae)** — new tab that runs [Novae](https://mics-lab.github.io/novae/)
  zero-shot spatial domain inference. Select species (human/mouse), optionally specify N domains
  (0 = auto-detect), and click "Run Novae Domains". On completion the cell labels layer is
  automatically recolored by inferred domains, the `novae_domains` key is added to the clustering
  dropdowns, and results are persisted to the session cache so they are restored on relaunch.
  Requires `pip install novae`. The full pipeline is recorded to `code.py`.

## [Unreleased] — 2026-03-19

### Added
- **Session-scoped code file** — each viewer session now writes to a timestamped file
  (`code_YYYYMMDD_HHMMSS.py`) instead of overwriting `code.py`. Opening a second dataset
  creates a fresh timestamped file; prior session files are preserved.
- **"Continue from existing code file" in Preferences** — lets users resume a previous
  session's code journal. Opens a file dialog, warns about potential duplication, then
  re-detects preamble/normalize/clustering tags from the file content so subsequent actions
  do not emit duplicate boilerplate.
- **Cluster label mapping recorded as executable dict** — editing cluster labels via the
  label editor now emits a `cluster_labels_<key> = {...}` dict to the code journal instead
  of a generic comment, making the output directly usable in downstream plotting calls.

### Changed
- **Auto-save plots on Show** — clicking "Show" in the Neighborhood Enrichment,
  Co-occurrence, Ligand-Receptor, Gene Correlation, and Gene Analysis (dotplot) tabs now
  automatically saves the figure to `<data_dir>/plots/<stem>.<format>` and reports the
  path in the status bar. The separate "Save Plot" button has been removed from all tabs.
- **Default plot format is now SVG** — plots are saved as vector graphics out of the box;
  PNG remains available via Preferences → Plot Format.

### Added
- **Gene Correlation tab** — scatter plot of per-cell expression for any two selected
  genes, annotated with Pearson r and Spearman ρ (both with p-values). Optionally
  restricted to the current cluster filter selection. Plot can be saved in the
  user's preferred format (PNG/SVG).
- **Normalisation options for Gene Correlation** — new "Normalisation" ComboBox with three
  choices: *Raw counts* (default off), *Fraction of total* (library-size correction), and
  *Log1p(CPM)* (default; matches normalisation used in DEG / dotplot tabs). Axis labels and
  status-bar output reflect the selected mode. `code.py` emits the correct normalisation
  snippet for each choice.
- **Filter small clusters by cell count** — new "Min cluster size" slider (100–10,000)
  and "Filter Small Clusters" button in the Cell Colouring tab. Clicking the button
  auto-unchecks clusters with fewer than N cells and enables the cluster filter,
  propagating to all downstream analyses (cell colouring, ligand-receptor,
  neighborhood enrichment, co-occurrence).

### Fixed
- **Rank genes results lost on crash** — rank genes results are now auto-saved immediately
  after computation via `save_rank_genes_incremental` in `utils/session.py`. This writes
  `rank_genes.parquet`, `rank_genes_adata_norm.h5ad`, and the `rank_genes_groupby` key
  to the zarr session. On session restore, all downstream buttons (dotplot, rank plot,
  volcano, edit/reset labels, export, LLM annotate) are re-enabled and the results
  preview is repopulated. Previously, only export and volcano were re-enabled and
  dotplot/rank-plot required re-running rank genes. Skipped when `--no-cache` is active.
- **Custom clusterings lost on crash** — Leiden runs and CSV imports now write an
  incremental parquet + update `custom_clustering_names` in the zarr session immediately
  after creation (`save_clustering_incremental` in `utils/session.py`). Previously,
  clusterings were only persisted on clean exit; a force-kill or crash meant they were lost.
  Skipped automatically when `--no-cache` is active or the zarr store does not yet exist.
- **Duplicate File/Preferences menus on dataset reload** — `create_file_menu()` now appends
  Xenium actions to napari's native `viewer.window.file_menu` (called once at startup) instead
  of creating a parallel custom "File" QMenu on every dataset load. `create_preferences_menu()`
  now finds-or-creates a single "Preferences" menu and clears stale actions before repopulating,
  so N dataset loads no longer produce N Preferences menus. The now-redundant
  `_attach_bare_file_menu()` helper and its cleanup code in `_do_full_init()` have been removed.
- **Dotplot crash with `'NoneType' object has no attribute 'get_axes'`** — in a
  `@thread_worker` context `plt.gcf()` returns `None`; the `_relabel_axes` fallback
  now checks `.fig` and `.figure` attributes explicitly and skips relabelling rather
  than crashing when no figure object is found.
- **Dotplot/Rank Genes Plot crash with duplicate cluster labels** — when multiple clusters
  share the same user-assigned label (e.g. several clusters all named "Unknown"), the old
  code called `rename_categories()` which raised `ValueError: Categorical categories must
  be unique`. Fixed by removing the rename + re-run approach entirely: plots now use the
  original integer cluster IDs from `uns['rank_genes_groups']` and apply label substitution
  post-hoc via axis tick replacement (`_relabel_axes`).
- **Minimap persists on dataset reload** — when opening a new dataset via File > Open Dataset,
  the old `MinimapWidget` is now explicitly hidden and scheduled for deletion before the new
  dataset is loaded, preventing stale thumbnails from overlapping the new minimap.
- **Cluster label dialog sorts numerically** — the "Edit Cluster Labels" dialog in the
  Clustering and Gene Analysis tabs now sorts cluster IDs numerically (0, 1, 2, … 19) instead
  of lexicographically (0, 1, 10, 11, … 19, 2, 3); falls back to string sort for non-integer IDs.

### Changed
- **"Edit Cluster Labels..." moved to Cell Colouring tab** — the button is now in the Cell
  Colouring tab (below "Apply Cell Coloring"), where it is more naturally discovered alongside
  the clustering selector and cluster filter controls. Removed from the Clustering tab.

### Added
- **LLM-based cluster annotation** — new "LLM Annotation" section in Gene Analysis tab.
  After running Rank Genes, click "Annotate with LLM" to send top 10 marker genes per
  cluster to a locally installed AI CLI (Claude, Gemini, or Codex). The LLM returns cell
  type annotations as JSON, which are stored as cluster labels and shown in dotplots,
  rank gene plots, and hover text. Runs in a background thread to keep UI responsive.

- **Cluster ID shown on hover** — when hovering over cells with custom labels (from
  CellTypist, label transfer, or LLM annotation), the status bar now shows both the
  numeric cluster ID and the label name, e.g. "Cell 12345 — leiden: 3 (T cells)".

- **Reference dataset auto-discovery** — the label transfer section now discovers
  reference datasets from `.metadata.json` sidecar files in `reference_datasets/`
  instead of using a hardcoded dictionary. Each h5ad gets a companion JSON file with
  structured metadata (paper, authors, journal, year, platform, tissue, annotation
  columns, default column, annotation workflow). Selecting a dataset in the ComboBox
  shows paper/platform/tissue info in the status bar. New datasets can be added by
  placing an h5ad + metadata.json pair in `reference_datasets/`.

- **`scripts/fetch_reference_datasets.py`** — standalone download & conversion script
  for 4 additional public scRNA-seq reference datasets: Prostate Cell Atlas
  (Clatworthy/Tuong 2021, direct h5ad), HuPSA (Cheng et al. 2024, Seurat .rds via
  rpy2), GSE176031 (Song/Huang 2021, GEO TXT matrices), GSE181294 (Mei et al. 2023,
  GEO MTX/CSV). Run with `--all` or `--dataset <name>`. Supports `--force` re-download,
  `--keep-raw` to preserve staging files. Each dataset generates both an h5ad file and a
  `.metadata.json` sidecar.

- **Label transfer via sc.tl.ingest()** — new "Label Transfer" section in Gene Analysis
  tab. Select a reference scRNA-seq dataset (two preconfigured prostate cancer datasets
  included, or browse for any h5ad file), load it, pick an annotation column, and click
  "Run Label Transfer" to project cell type labels onto Xenium clusters. Pipeline: find
  common genes → subset both → preprocess reference (normalize, log1p, HVG, PCA, neighbors,
  UMAP) → preprocess Xenium (normalize, log1p) → `sc.tl.ingest()` → majority vote per
  cluster. Reference datasets are cached in memory after first load. Preconfigured datasets:
  Zhao et al. (222K cells, ct.main/ct.sub/ct.sub.epi) and eBioMedicine (21K cells,
  ID/ID_coarse/refinedID). New utility functions: `load_reference_h5ad()`,
  `get_annotation_columns()`, `run_label_transfer()` in `utils/gene_analysis.py`.

- **CellTypist automated cell type annotation** — new "CellTypist Annotation" section
  in the Gene Analysis tab. Select a CellTypist model from a dropdown, click "Annotate
  with CellTypist" to run per-cell predictions, then majority-vote assigns cell type
  labels to each cluster. Labels are stored in `state["cluster_labels"]` (same system
  as manual labels), so dotplots, rank genes plots, hover info, and exports all pick
  up the new names automatically. Includes "Download Models" button for one-time model
  download, and a **confidence threshold slider** (default 0.5) that filters out
  low-confidence per-cell predictions before the majority vote — prevents immune-focused
  models from labeling non-immune clusters as "T cells". Status bar shows how many cells
  passed the threshold. Gracefully degrades when celltypist is not installed (widgets
  shown but disabled). New `run_celltypist_annotation()` utility in
  `utils/gene_analysis.py`.

- **Reset Labels button** in Gene Analysis tab — quickly reverts cluster labels back to
  original numeric IDs after CellTypist annotation or manual editing. Removes the entry
  from `state["cluster_labels"]` for the selected clustering.

- **Complete code recording coverage** — all analysis actions now record parameters to
  `code.py` via `ctx.record_code()`:
  - **Transcripts**: QV threshold in overlay recording; new density heatmap recording
    (gene, bin_size, normalise, cluster_filter)
  - **Ligand-receptor**: interaction database selection and CellPhoneDB flag
  - **ROI DEG**: cluster filter status and selected clusters
  - **UMAP**: display action recorded as comment
  - **Cell coloring**: background color toggle
  - **H&E registration**: flip state changes
  - **ARMS overlay**: flip state changes; Xenium cluster filter in tile DEG
  - **Co-occurrence**: cluster subset in run recording; filter_targets in plot recording
  - **Cluster labels**: label editing recorded in `_helpers.py`

### Fixed
- **Progress feedback not showing** — two bugs prevented tqdm progress from appearing in
  the status bar during nhood/L-R/co-occurrence analyses:
  - QTimer objects were garbage-collected immediately after `worker.start()` because no
    Python reference was retained; fixed by storing timers in `state['_spinner_timer']`
    and `state['_progress_timer']`
  - `qt_tqdm_context` was not patching `tqdm.std.tqdm` (only `tqdm.tqdm` and
    `tqdm.auto.tqdm`); squidpy's `parallelize()` imports from `tqdm.std` when
    `ipywidgets` is absent, so the monkey-patch had no effect; fixed by also patching
    `tqdm.std.tqdm`
  - `attach_tqdm_progress` now returns `(post_fn, timer)` instead of just `post_fn`

### Added
- **Progress feedback for long-running analyses** — status bar now updates during analyses
  instead of showing a static "Running…" message:
  - **Leiden clustering**: animated spinner (`| / - \`) with stage messages
    (`"Preparing data…"` → `"Computing neighbors…"` → `"Running Leiden algorithm…"`)
  - **Rank genes**: animated spinner while scanpy computes
  - **Neighborhood enrichment**: live permutation counter (`"Enrichment permutations: 450/1000 (45%)"`)
  - **Ligand-receptor**: live permutation counter (`"L-R permutations: 230/1000 (23%)"`)
  - **Co-occurrence**: live interval counter (`"Co-occurrence: 23/50 iterations"`)
- New helpers in `tabs/_helpers.py`: `attach_spinner`, `attach_tqdm_progress`,
  `qt_tqdm_context`, `ProgressMailbox`

### Changed
- **Cursor hover debounce** — `_on_cursor_move` in `02_xenium_viewer.py` now fires
  `get_value()` at most every 80 ms via a `QTimer` debounce, eliminating main-thread
  stutter during rapid pan/zoom.
- **Cluster coloring off main thread** — cluster filter (numpy masking) and
  `build_direct_label_colormap()` (100K dict allocations) now run inside the
  `@thread_worker`; the main-thread callback only sets `colormap` and calls
  `refresh()`, eliminating the UI stall on Apply.
- **Faster colormap dict** — `build_direct_label_colormap()` uses bulk `.tolist()`
  (C-level) before building the `color_dict`, eliminating 100K+ Python object
  allocations per colormap rebuild.
- **Revert vectorized cluster loop** — `get_cluster_colors()` reverted to a simple
  Python `for` loop; result is cached so the vectorization complexity wasn't
  visible to users.

### Fixed
- **"Filter by selected clusters" in Transcript Density** — the checkbox was
  silently ignored when `label_to_cluster` was not yet populated (e.g. before
  clicking "Apply Cell Coloring" or after switching to gene coloring). The fix
  pre-computes filter data on the main thread in `on_compute_density()`: if a
  clustering is selected it builds `label_to_cluster` on-the-fly when needed;
  if no clustering is available it shows "No clustering applied — filter
  skipped" in the status label and returns early instead of silently falling
  through.

### Added
- **Transcript cache completeness check** (`00_preprocess_transcripts.py`) — writes a
  `.complete` sentinel file to `transcript_cache/` on successful completion, storing the
  parquet mtime and size. Re-running the script detects an up-to-date cache and exits
  immediately. If the parquet has changed, stale feather files are wiped before
  reprocessing, preventing silent data duplication.

- **Interactive minimap** (`utils/minimap_widget.py`) — floating overlay in the
  top-right corner of the napari canvas showing the DAPI thumbnail with a white
  viewport rectangle. Updates live as the camera pans/zooms; clicking navigates
  the camera. Repositions itself when the window is resized. Instantiated
  automatically after dataset load in `_do_full_init()`.

- **Transcript density bins** — new "Transcript Density" section in the
  Transcripts tab. Aggregates per-gene transcript counts into spatial bins
  (user-controlled µm bin size via slider) and displays a `transcript_density`
  napari Image layer (hot colormap). Supports optional cluster-based filtering
  using centroids from the Cell Coloring tab. Computation runs in a background
  thread worker; contrast limits are set to the 99th percentile of non-zero bins.
  New `transcript_bins_layer` field added to `ViewerContext`.

- **Transcript density normalisation** — "Normalise by cells per bin" checkbox
  in the Transcript Density section. When ticked, each bin value is divided by
  the number of cells in that bin (transcripts-per-cell), correcting for local
  cell-density variation. Bins with zero cells remain zero. When combined with
  "Filter by selected clusters", the cell count per bin uses only selected-cluster
  cell centroids.


- **Empty-viewer startup** — cancelling the startup folder dialog (or omitting the
  CLI path argument) now opens a bare napari window instead of calling `sys.exit(0)`.
  A minimal File menu (Open Dataset... `Ctrl+O` + Preprocess Dataset...) is attached
  directly to the bare viewer via the new `_attach_bare_file_menu()` helper. When the
  user picks a dataset from File → Open Dataset, `_do_full_init()` builds the full UI
  (layers, dock widget, session restore) and replaces the bare menu. Works identically
  for the first load and subsequent dataset switches.

- **`_do_full_init()` helper** — extracted from `main()`: loads a dataset, creates
  `UMAPViewer`, populates napari layers, builds `ViewerContext`, rebuilds the control
  panel dock widget, and restores the session. Called both at initial startup (when
  `data_path` is given) and from `_on_open_dataset()` for every subsequent switch.
  `main()` is now ~100 lines shorter and `_on_open_dataset()` no longer duplicates the
  load/build sequence.

### Changed
- **Zarr skip logic in `PreprocessWorker`** — `PreprocessWorker.run()` now performs a
  staleness pre-check before calling `load_sdata()`: if `sdata_cached.zarr` already
  exists and its mtime ≥ `experiment.xenium` mtime, the zarr step is skipped and the
  progress dialog shows "Zarr cache already up to date, skipping: …" instead of the
  misleading "Creating zarr cache…" message. Uses the same mtime comparison as
  `load_sdata()`.

- **`_on_open_dataset()` refactored** — now uses `nonlocal ctx` and delegates the
  load/build sequence to `_do_full_init()`. Cleanup path (session save, UMAP close,
  generation counter increment, ref nulling, gc) is guarded by `if ctx is not None`
  so it also handles the first load from an empty viewer correctly.

- **`_on_preprocess_dataset()` hardened** — all `ctx` accesses are guarded by
  `if ctx is not None`, so preprocessing works from an empty viewer. The confirmation
  dialog's "session will be cleared" line is omitted when no dataset is loaded.

- **`_parse_args()`** — dialog cancel now returns `(None, no_cache)` instead of
  calling `sys.exit(0)`. `experiment.xenium` validation only runs when a path was
  actually selected.

### Added (continued from 2026-03-16)
- **File > Preprocess Dataset...** — new menu item that runs both preprocessing
  steps (zarr cache creation via `01_load_sdata.py` and per-gene transcript
  feather splitting via `00_preprocess_transcripts.py`) in a background thread
  so the UI stays responsive. Auto-detects whether the selected folder is a
  single Xenium dataset or a parent folder containing several datasets, and
  processes all found datasets sequentially. A modal progress dialog with a
  live progress bar is shown during preprocessing. The current viewer data is
  cleared before starting to avoid a RAM peak of old dataset + I/O overhead.
  Implemented via `PreprocessWorker(QThread)`, `_find_xenium_datasets()`, and
  `_make_progress_dialog()` in `02_xenium_viewer.py`; `create_file_menu()` in
  `tabs/_helpers.py` updated to accept the new `on_preprocess_dataset` callback.

### Changed
- **Memory: free old dataset before loading new one** — in `_on_open_dataset()`,
  all heavy `ctx` fields (`sdata`, `adata`, `clusterings`, `color_manager`,
  `transcript_loader`, layer references, etc.) are now explicitly set to `None`
  and `gc.collect()` is called after clearing napari layers (step 8b) and before
  `_load_dataset()` is called (step 9). This prevents peak RSS from reaching
  old-dataset + new-dataset simultaneously during a dataset switch.

## [Unreleased] — 2026-03-14

### Added
- **File > Open Dataset** — new menu entry (and `Ctrl+O` shortcut) lets the user
  switch to a different Xenium output directory without restarting the application.
  The current session is saved automatically before the swap; the new dataset's
  session is restored if a zarr cache exists. All napari layers are replaced and
  the full control-panel dock widget is rebuilt in-place.

  Implementation details:
  - `utils/viewer_context.py` — added `dataset_generation: int = 0` field.
  - `tabs/_helpers.py` — new `create_file_menu(ctx, on_open_dataset)` function;
    `_build_control_panel()` now accepts an optional `on_open_dataset` callback.
  - `02_xenium_viewer.py` — extracted `_load_dataset()`, `_populate_viewer()`,
    `_snapshot_layers()`, `_make_initial_state()`, `_make_initial_he_state()`, and
    `_make_initial_arms_state()` as standalone helpers; introduced `_app` container
    in `main()`; implemented `_on_open_dataset()` closure with re-entry guard,
    session save/restore, and full viewer reload sequence.
  - Stale thread-worker guard added to six tab modules (`tab_transcripts.py`,
    `tab_cell_coloring.py`, `tab_clustering.py`, `tab_roi.py`,
    `tab_he_registration.py`, `tab_arms.py`): each worker captures
    `ctx.dataset_generation` at dispatch time and silently drops its result if
    the generation has changed when the worker returns.

## [Unreleased] — 2026-03-10

### Fixed
- **Co-occurrence plot colors now match labels layer** — replaced seaborn `tab20`
  fallback in `make_co_occurrence_plot()` with `CLUSTER_PALETTE` from
  `utils/coloring.py`, ensuring line colors match the cell labels layer colors.

### Changed
- **Custom nhood enrichment heatmap now matches squidpy native style** — rewrote
  `make_nhood_enrichment_plot()` in `spatial_analysis.py` to use `imshow` +
  `make_axes_locatable` instead of `sns.heatmap`. Uses `viridis` colormap with
  symmetric normalization for z-score, colored cluster category bars on left + top
  edges (using `CLUSTER_PALETTE`), cluster labels on the left bar, right-side
  vertical colorbar with `%0.2f` ticks, no gridlines. Added `annotate` parameter
  (default False) with luminance-based white/black text. Figsize scales with cluster
  count matching squidpy's `(2*n//3, 2*n//3)` formula.

- **Decomposed `02_xenium_viewer.py` into per-tab modules** — the 4295-line monolith
  has been split into 14 new files. `_build_control_panel()` shrinks from ~3850 lines
  to ~80 lines (thin orchestrator). Each of the 11 tabs is now a standalone module in
  `scripts/tabs/` with a single `build_tab(ctx)` entry point. Cross-tab state is managed
  through a `ViewerContext` dataclass (`scripts/utils/viewer_context.py`) instead of 18
  closure-captured parameters. Shared helpers (status proxy, code recording, clustering
  refresh, label editor, plot save, etc.) live in `scripts/tabs/_helpers.py`. No
  behavioral changes — all features, session persistence, and code recording work
  identically.

  New files:
  - `scripts/utils/viewer_context.py` — ViewerContext dataclass
  - `scripts/tabs/__init__.py` — re-exports all tab builders
  - `scripts/tabs/_helpers.py` — make_tab, StatusProxy, create_shared_helpers, create_preferences_menu
  - `scripts/tabs/tab_clustering.py` — Tab 0: Leiden, import/export
  - `scripts/tabs/tab_cell_coloring.py` — Tab 1: gene/cluster coloring, filter checkboxes
  - `scripts/tabs/tab_transcripts.py` — Tab 2: multi-gene transcript overlay
  - `scripts/tabs/tab_umap.py` — Tab 3: UMAP viewer
  - `scripts/tabs/tab_roi.py` — Tab 4: ROI analysis + ROI DEG
  - `scripts/tabs/tab_he_registration.py` — Tab 5: H&E registration
  - `scripts/tabs/tab_gene_analysis.py` — Tab 6: rank genes, dotplot, volcanos
  - `scripts/tabs/tab_ligrec.py` — Tab 7: ligand-receptor
  - `scripts/tabs/tab_nhood.py` — Tab 8: neighborhood enrichment
  - `scripts/tabs/tab_co_occurrence.py` — Tab 9: co-occurrence
  - `scripts/tabs/tab_arms.py` — Tab 10: ARMS overlay

### Added
- **ARMS tab code recording** — All ARMS operations (H&E load, landmark registration,
  GeoJSON/CSV load, tile DEG, DEG export, volcano plots) now emit `_record_code()`
  entries so they appear in the reproducible `code.py` journal.

## [2026-03-09]

### Fixed
- **ARMS metadata not persisting across sessions** — `save_session()` overwrote
  real-time-saved ARMS attrs (filename, affine, shape, paths) with empty snapshot
  values, causing ARMS registration to be lost on reload. Fixed by: (1) preserving
  existing ARMS attrs before the overwrite and using them as fallback, (2) saving
  all essential ARMS metadata in `_save_arms_affine_to_sdata()` (not just
  affine/flips), (3) calling `_save_arms_affine_to_sdata()` after ARMS restore to
  repair any previously corrupted attrs.
- **NameError on close** — `_on_viewer_closing` referenced `_arms_state` which is
  defined later in `main()`; wrapped in `try/except NameError` so closing the
  viewer before the ARMS tab initializes no longer crashes.

### Added
- **ARMS Tile DEG analysis** — "Run ARMS Tile DEG" button in ARMS Overlay tab
  performs differential expression between cells grouped by ARMS tile cluster ID.
  Tile polygons are transformed from GeoJSON space to Xenium pixel space via the
  registration affine before point-in-polygon tests. Clusters with <10 cells are
  excluded. Results display in a text widget and can be exported to CSV. Session
  persistence saves/restores DEG results across viewer restarts.
- **ARMS pairwise volcano plots** — "Generate ARMS Volcano Plots..." button
  (enabled after ARMS Tile DEG completes) runs pairwise DEG for every pair of
  ARMS tile clusters and saves volcano PNGs to a user-selected directory.
  Progress updates shown in status bar. `compute_arms_tile_deg()` now returns
  the normalized subset adata as a 3rd element for downstream pairwise analysis.
- **Pairwise volcano plots** — "Generate All Volcano Plots..." button in Gene
  Analysis tab runs DEG for every pairwise cluster comparison and saves volcano
  PNGs (3-color scatter with threshold lines and top gene labels) to a
  user-selected directory. Progress updates shown in status bar.
- **ARMS Overlay tab (Tab 10)** — load a larger ARMS H&E image, align it to
  Xenium coordinates via manual landmark registration (same workflow as the H&E
  Registration tab), then load GeoJSON tile boundaries + CSV cluster IDs to
  display 288 tile polygons colored by genomic cluster (1–8). The affine
  transform is applied to both the H&E image and the polygon shapes layer, so
  re-registering landmarks automatically re-aligns everything.
- **Auto-save reproducible code on exit** — code recording is now on by default
  and automatically saves to `{data_dir}/code.py` when the viewer closes (if
  recording is enabled and journal is non-empty).
- **Persist custom clusterings across sessions** — Leiden and imported clusterings
  are saved as parquet files in `viewer_session/clusterings/` inside the zarr
  cache, and restored on next launch.
- **ARMS overlay session persistence** — ARMS H&E image, registration affine,
  flip state, landmarks, and GeoJSON/CSV file paths are saved to the zarr cache
  and restored on next launch. The H&E image is stored as
  `sdata.images["arms_he_image"]` with a multiscale pyramid; the affine and
  landmarks are saved in `viewer_session/arms/`. GeoJSON/CSV tiles are re-parsed
  from the original file paths on restore (with a warning if files have moved).

### Changed
- `record_code` default changed from `False` to `True` (always-on recording).

## [Unreleased] — 2026-03-06

### Added
- **Import clustering from CSV/TSV** — "Import Clustering..." button in the
  Clustering tab. Reads a file with `cell_id` + `group` columns (auto-detects
  tab vs comma separator). Supports string-valued group names (e.g. imported
  cell type annotations). Imported clusterings appear in all downstream dropdowns.
- **Export clustering to CSV/TSV** — "Export Clustering..." button exports the
  currently selected clustering with cell_id and group columns. Custom cluster
  labels are applied to the exported group names when available.
- **Per-clustering label storage** — cluster labels are now stored per-clustering
  (nested dict `{clustering_name: {cluster_id: label}}`), so labels don't collide
  across different clusterings. Session save/restore updated with backward
  compatibility for the old flat format.
- **Multi-column label editor** — the "Edit Cluster Labels..." dialog now uses a
  multi-column grid layout (up to 3 columns, ~10 rows per column) with a scroll
  area, preventing the dialog from being too tall for clusterings with many groups.
  Handles both integer and string cluster IDs.
- **Label propagation to all plots** — cluster labels now appear in:
  - Neighborhood enrichment heatmap (axis tick labels)
  - Co-occurrence plots (subplot titles + legend entries)
  - Ligand-receptor dotplot (column axis labels)
  - Rank genes panel plot (group titles)
  - Mouse hover status bar (shows label instead of raw cluster ID)
  - Cluster filter checkboxes (show labels when available)
  - Dotplot (already worked, now uses per-clustering lookup)

### Changed
- **Cluster label editor shared between tabs** — Gene Analysis and Clustering
  tabs now share the same `_build_label_editor_dialog()` helper, eliminating
  code duplication. Both editors use the per-clustering label storage.
- **String cluster ID support** — `add_clustering_to_obs()` in gene_analysis.py
  and `_get_cluster_ids_per_obs()` in the viewer now handle string-valued
  cluster IDs (from imported clusterings) via factorization fallback.

## [Previous] — 2026-03-05

### Added
- **Cluster label editor in Clustering tab** — "Edit Cluster Labels..." button
  in the Clustering tab opens a dialog to rename clusters (manual cell type
  annotation) using the Cell Coloring tab's selected clustering.

### Changed
- **Co-occurrence plot colors match napari** — line colors in co-occurrence plots
  now use the same palette as napari cell coloring (CLUSTER_PALETTE) instead of
  seaborn tab20. Falls back to tab20 for clusters without a stored color.
- **Clustering sync across tabs** — selecting a clustering in the Cell Coloring
  tab now auto-sets the same clustering in Gene Analysis, Ligand-Receptor,
  Nhood Enrichment, and Co-occurrence tabs (one-directional sync).

## [Previous] — 2026-03-04

### Added
- **Leiden clustering tab** — new "Clustering" tab (Tab 0) with configurable
  `n_neighbors` (5–50), `n_pcs` (10–50), and `resolution` (0.1–5.0) parameters.
  Runs `sc.pp.neighbors` + `sc.tl.leiden` on a worker thread, stores results as
  `leiden_r{resolution}` in the clusterings dict, and refreshes all downstream
  ComboBoxes (Cell Coloring, Gene Analysis, L-R, Nhood, Co-occurrence).
  Reproducible code recording supported.

## [Previous] — 2026-03-03

### Added
- **Interaction database filtering for L-R analysis** — 4 checkboxes (OmniPath,
  LigRecExtra, PathwayExtra, KinaseExtra) to select which OmniPath interaction
  datasets are queried, plus a "CellPhoneDB only" toggle to restrict to canonical
  ligand-receptor pairs. Selections are passed via `interactions_params` to
  `sq.gr.ligrec()`. Unchecking PathwayExtra + KinaseExtra removes intracellular
  proteins (e.g. TP53) that leak through as false "ligands".

- **Plot format preference (PNG / SVG)** — new Preferences → Plot format menu in the
  napari menu bar with exclusive PNG/SVG radio actions. All 4 plot save handlers
  (dotplot, L-R, nhood enrichment, co-occurrence) now respect the chosen format.
  Default is PNG. `matplotlib.savefig()` infers format from file extension.

- **Native squidpy heatmap** for Neighborhood Enrichment — `on_show_nhood_plot()`
  now uses `sq.pl.nhood_enrichment()` when no cluster filter is active, falling
  back to the custom `make_nhood_enrichment_plot()` when clusters are filtered or
  on session restore. `_adata_norm` and `_cluster_key` are now stored on the result
  dict (matching co-occurrence).

- **"Filter targets" checkbox** in Co-occurrence tab — when checked with a cluster
  filter active, restricts target lines (not just query subplots) to the selected
  clusters. Three-path plot logic: (1) filter targets ON → custom plot with both
  query and target restricted, (2) filter targets OFF → squidpy native plot with
  cluster-filtered subplots, (3) session restore fallback → custom plot.
  `make_co_occurrence_plot()` gains a `target_clusters` parameter.

- **Co-occurrence tab** (Tab 9) — squidpy-based spatial co-occurrence analysis
  (`sq.gr.co_occurrence`). Computes how cluster types co-occur spatially across
  increasing distance radii. Configurable: clustering, distance bins (10-100).
  Displays summary with cluster count and distance range. Line-plot visualization
  showing co-occurrence probability vs distance for selected clusters — one subplot
  per query cluster with colored lines for each target cluster and baseline at y=1.
  Uses Cell Coloring tab's cluster filter to select which clusters get subplots.
  Export long-form CSV, save plot as PNG. Session persistence via zarr. New functions
  in `spatial_analysis.py`: `run_co_occurrence()`, `make_co_occurrence_plot()`.

- **Neighborhood Enrichment tab** (Tab 8) — squidpy-based neighborhood enrichment
  analysis (`sq.gr.nhood_enrichment`). Computes which cluster types are spatially
  enriched or depleted in each other's neighborhoods via permutation testing.
  Configurable: clustering, permutation count (100-1000), neighbor count (3-20).
  Displays summary with top enriched/depleted cluster pairs. Heatmap visualization
  (z-score or count mode) via seaborn with diverging colormap, supports filtering
  by the Cell Coloring tab's cluster selection. Export z-score matrix as CSV, save
  plot as PNG. Session persistence via zarr. New functions in `spatial_analysis.py`:
  `run_nhood_enrichment()`, `make_nhood_enrichment_plot()`.

### Changed
- Renamed `_get_lr_cluster_filter()` to `_get_cluster_filter()` — shared helper for
  cluster filtering in both L-R plot and nhood enrichment plot.

## [Unreleased] — 2026-03-02

### Added
- **Session persistence** — viewer state (ROIs, H&E registration, analysis results,
  cluster labels) is automatically saved to `sdata_cached.zarr/` when the viewer closes.
  The H&E image is stored as a spatialdata multiscale image element (`images/he_image`)
  with its affine transformation, so the next launch restores the H&E overlay with
  registration already applied — no need to re-load or re-register. ROI polygons,
  cluster labels, and analysis results (rank genes, ROI DEG, L-R) are persisted in
  `viewer_session/` and restored on startup. Affine transformations are saved in
  real-time (on each registration/flip change). Skipped when `--no-cache` is used.
  New module `scripts/utils/session.py` with `save_session()` and `load_session()`.

### Previously added
- **Gene Analysis tab** (Tab 6) — rank marker genes per cluster using Wilcoxon, t-test, or
  logreg methods via scanpy's `rank_genes_groups`. Features: configurable top-N genes,
  dotplot visualization with optional dendrogram, editable cluster labels (dialog to rename
  cluster IDs for publication-ready plots), rank genes panel plot, results preview in
  monospace text area, and full CSV export. New module `scripts/utils/gene_analysis.py` with
  `get_normalized_adata()` (cached log-normalization + PCA), `add_clustering_to_obs()`,
  `run_rank_genes()`, `make_rank_genes_dotplot()`, `make_rank_genes_plot()`.
- **ROI differential expression** — extended ROI Analysis tab (Tab 4) with genome-wide
  DEG between drawn ROI polygons. Draw >= 2 ROIs, select method (Wilcoxon/t-test), optional
  cluster filtering, then compute DE across all genes. Uses `compute_roi_deg()` which
  normalizes the subset, runs `rank_genes_groups` with `reference='rest'`, and returns
  full results table. Export to CSV supported.
- **Ligand-Receptor tab** (Tab 7) — squidpy-based spatial ligand-receptor interaction
  analysis. Computes spatial neighbor graph, runs permutation-based L-R testing
  (`sq.gr.ligrec`), displays summary of significant interactions, and generates
  interaction dotplot. Configurable: clustering, permutation count (100-1000), neighbor
  count (3-20), p-value threshold. Export means and p-values as CSV, save plot as PNG.
  Warns when the 480-gene Xenium panel lacks common L-R pairs. New module
  `scripts/utils/spatial_analysis.py` with `compute_spatial_neighbors()`, `run_ligrec()`,
  `make_ligrec_plot()`.

## [Unreleased] — 2026-03-01

### Added
- **H&E image registration** — new "H&E Registration" tab (5th tab). Load an H&E
  OME-TIFF/SVS image as a multiscale overlay, place paired landmark points on both
  Xenium and H&E images, then compute a similarity transform (rotation + uniform scale +
  translation) to align the H&E onto Xenium coordinates. Features: opacity slider,
  per-landmark residual display, save/load landmarks as JSON, save affine matrix.
  Flip checkboxes (vertical/horizontal) let you mirror the H&E before placing
  landmarks, useful when the H&E was scanned in a different orientation.
  New module `scripts/utils/registration.py` with `load_he_pyramid()`,
  `compute_landmark_affine()`, `save_landmarks()`, `load_landmarks()`.
- **Coarse tissue-outline alignment** — "Coarse Align" button in H&E Registration tab
  automatically snaps H&E roughly into position before manual landmark placement.
  Extracts tissue masks from morphology_focus (max-projection across channels) and H&E
  (HSV saturation channel), then computes a similarity transform using OpenCV image
  moments with multi-rotation hypothesis testing (0/90/180/270 deg, best IoU wins).
  The coarse affine is a visual stepping-stone; computing landmark registration replaces
  it entirely. New functions in `registration.py`: `extract_tissue_mask()`,
  `extract_tissue_mask_fluorescence()`, `extract_tissue_mask_he()`,
  `compute_coarse_affine()`.
- **ROI polygon analysis** — new "ROI Analysis" tab (4th tab). Draw polygons on the
  napari Shapes layer ("ROIs"), click "Calculate Expression" to get per-region stats
  (cell count, mean, median, std, min, max) for the selected gene. "Export CSV" saves
  per-cell data (region_id, cell_id, x/y centroid in microns, expression). Uses shapely
  `contains_xy` for fast point-in-polygon queries on precomputed cell centroids.
  Includes pairwise Welch's t-tests between regions with Benjamini-Hochberg
  correction when >1 comparison.

### Changed
- **Multi-cluster filter** — "Filter by cluster" now supports selecting multiple clusters
  simultaneously via checkboxes (3-column grid with Select All / Deselect All). Works
  in both Gene Expression and Cluster coloring modes. Replaces the single-cluster
  ComboBox dropdown.
- **White background toggle** — new checkbox at the top of the Cell Coloring tab sets the
  napari canvas background to white, useful for H&E overlay visibility or light-colored
  cluster inspection.
- **Cluster filter works in both modes** — "Filter by cluster" checkbox and cluster ID
  selector are now available in Cluster coloring mode (not just Gene Expression). When
  enabled, only cells belonging to the selected clusters are shown; all others are
  transparent.
- **UMAP window no longer auto-opens** — `color_by_gene()` and `color_by_cluster()`
  store colors/metadata but only update the viewer if already open. Click "Show UMAP
  Window" to open it manually; colors and title are applied on first open.
- **Min/max colormap scaling** — gene expression coloring now normalizes non-zero cells
  using vmin/vmax of expressing cells (instead of 0-to-max). This spreads the full
  viridis range across expressing cells so low-expressors are visually distinguishable.
  Zero-expression cells remain transparent.

### Previously added
- **Multi-gene transcript overlay** — add up to 10 genes simultaneously with distinct
  colors (yellow, cyan, magenta, orange, green, sky blue, red, violet, pink, brown).
  New gene list widget with Add/Remove/Clear All buttons and color legend.
  `TranscriptLoader.get_multi_gene_points()` merges per-gene points with palette colors.
- **Cluster filter for gene expression** — checkbox "Filter by cluster" in Cell Coloring
  section. When enabled, only cells in the selected cluster are colored; rest are
  transparent. `CellColorManager.get_gene_colors_filtered()` method.
- **Generic dataset loading** — CLI argument `python 02_xenium_viewer.py /path/to/data`
  or Qt file dialog if no argument given. Reads `pixel_size` from `experiment.xenium`.
  All scripts accept `--data-dir` / positional argument. Transcript cache stored alongside
  data (`data_dir/transcript_cache/`) instead of inside `scripts/`.
- **Zarr cache** — `01_load_sdata.py` writes `sdata_cached.zarr` on first load;
  subsequent launches use the cache (~60-70% faster). Staleness detection via
  `experiment.xenium` mtime. `--no-cache` flag to force rebuild.
- **Deferred UMAP scatter** — UMAP widget shows placeholder text until first
  "Apply Cell Coloring" click, avoiding 318K-point scatter at startup.
- **Timing instrumentation** — wall-clock times printed for each startup phase.

### Changed
- `load_sdata()` now skips `cells_boundaries` and `nucleus_boundaries` (were hidden
  anyway; 318K polygon shapes freeze napari). Saves ~30-40% of xenium() load time.
- `load_umap()` and `load_clusterings()` accept a `path` parameter instead of using
  module-level globals.
- Transcript cache default location changed from `scripts/transcript_cache/` to
  `data_dir/transcript_cache/`.

## [0.1.0] — 2026-02-28

### Added
- `scripts/00_preprocess_transcripts.py` — one-time script to split `transcripts.parquet`
  (1.3 GB, 128M rows) into per-gene feather files in `scripts/transcript_cache/`.
  Enables <100 ms per-gene transcript loading in the viewer (vs 4–5 s from parquet).
- `scripts/01_load_sdata.py` — SpatialData loader using `spatialdata_io.xenium()` with
  5-level software image pyramid for the morphology_focus TIFFs (no internal pyramid).
  Also loads UMAP coordinates and all cluster assignments from `analysis/`.
- `scripts/utils/coloring.py` — `CellColorManager`: builds napari `DirectLabelColormap`
  from gene expression (viridis/magma/plasma/RdBu/YlOrRd) or cluster assignments.
  Transparent alpha for zero-expression cells. LRU cache on repeated gene queries.
- `scripts/utils/transcript_index.py` — `TranscriptLoader`: loads per-gene transcript
  x/y coordinates from feather cache with parquet fallback.
- `scripts/utils/umap_widget.py` — `UMAPWindow`: linked matplotlib UMAP scatter
  (318K cells, rasterized=True for performance). Syncs colors with napari Labels layer.
  Handles 91-cell UMAP/adata mismatch via `.reindex()`.
- `scripts/02_xenium_viewer.py` — main entry point. Opens napari with:
  - 4-channel morphology_focus image (DAPI, ATP1A1/CD45/E-Cad, 18S, AlphaSMA/Vim)
  - Cell and nucleus labels layers (raster masks for fast coloring via DirectLabelColormap)
  - Transcript points layer (populated on demand from feather cache)
  - Magicgui control panel: gene/cluster coloring, colormap, transcript toggle, QV slider
  - Linked matplotlib UMAP window
- `.gitignore` — excludes parquet, ome.tif, zarr.zip, h5, and transcript_cache/
