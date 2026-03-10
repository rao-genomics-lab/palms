# Changelog

## [Unreleased] — 2026-03-10

### Changed
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
