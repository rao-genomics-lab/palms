# Changelog

## [Unreleased] — 2026-03-01

### Added
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
