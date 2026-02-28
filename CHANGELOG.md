# Changelog

## [Unreleased] — 2026-02-28

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
