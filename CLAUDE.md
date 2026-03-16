# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Xenium Viewer** is a napari-based spatial transcriptomics viewer for Xenium 3.x output data — a Linux equivalent to the commercial Xenium Explorer. It visualizes high-resolution spatial gene expression data with cell-level resolution.

## Running the Viewer

```bash
# Activate environment
conda activate xenium_viewer

# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Without this, transcript loading falls back to scanning transcripts.parquet (~5s/gene instead of <100ms)
python 00_preprocess_transcripts.py

# Launch viewer (file dialog opens if no path given)
python 02_xenium_viewer.py [/path/to/xenium/output/]

# Launch without SpatialData zarr cache
python 02_xenium_viewer.py /path/to/xenium/output/ --no-cache
```

There is no test suite or CI/CD. All testing is manual/exploratory.

## Architecture

### Entry Points & Load Sequence

1. **`02_xenium_viewer.py`** — Main entry point (~600 lines). Validates data dir, orchestrates loading, builds napari viewer with layers, instantiates all managers, creates `ViewerContext`, builds all tab widgets, then restores session.
2. **`01_load_sdata.py`** — Loads SpatialData from Xenium output. Uses a zarr cache (`sdata_cached.zarr/`) for 60–70% faster subsequent launches; staleness detected via `experiment.xenium` mtime.
3. **`00_preprocess_transcripts.py`** — One-time step that splits the transcript parquet into ~480 per-gene feather files for fast per-gene loading (~100ms vs 4–5s).

### Central State: ViewerContext

**`utils/viewer_context.py`** — `ViewerContext` dataclass is the single shared state object passed to all tab modules. It holds:
- Core data: `viewer`, `adata`, `sdata`, `clusterings`, pixel size
- Layer references: `cell_labels_layer`, `transcript_layer`, `roi_layer`
- Manager objects: `CellColorManager`, `TranscriptLoader`, `UMAPViewer`
- Mutable state dicts: `state` (general), `he_state` (H&E registration), `arms_state` (ARMS overlay)
- Shared callables: `record_code()`, `set_status()`, `refresh_clustering_choices()`

### Tab Modules (`tabs/`)

Each tab follows a consistent pattern:

```python
def build_tab(ctx: ViewerContext) -> tuple[QWidget, dict]:
    # Build UI, register callbacks that access ctx
    # Return (tab_widget, exports_dict)
    # exports_dict may contain "restore_session" callable
```

The 11 tabs cover: Clustering, Cell Coloring, Transcripts, UMAP, ROI Analysis, H&E Registration, Gene Analysis, Ligand-Receptor, Neighborhood Enrichment, Co-occurrence, and ARMS Overlay.

`tabs/_helpers.py` contains shared utilities (e.g., `StatusProxy`, `make_tab()`).

### Key Utilities (`utils/`)

| Module | Purpose |
|---|---|
| `coloring.py` | `CellColorManager` with `DirectLabelColormap` for O(nonzero) raster colorization |
| `gene_analysis.py` | Rank genes, normalization, Leiden clustering |
| `spatial_analysis.py` | Squidpy-based spatial analysis (neighborhood enrichment, co-occurrence, L-R) |
| `registration.py` | Landmark-based similarity affine registration for H&E/ARMS |
| `transcript_index.py` | Per-gene feather loader |
| `session.py` | Zarr-based session persistence (ROIs, H&E/ARMS registration, clusterings, DEG results) |
| `umap_widget.py` | Separate linked napari window for UMAP scatter |

### Session Persistence

Stored in `sdata_cached.zarr/viewer_session/` as zarr arrays and parquet files. Auto-saved on relevant actions and restored at startup. Supports zarr v2 and v3.

### Code Recording

All user actions generate reproducible Python code saved to `data_dir/code.py`. Use `ctx.record_code(snippet)` in tab callbacks. Preamble (imports, data loading) is auto-inserted.

## Key Dependencies

- **napari** + **PyQt5/qtpy** + **magicgui** — UI framework
- **scanpy** / **anndata** — single-cell analysis
- **squidpy** — spatial transcriptomics analysis
- **spatialdata** / **spatialdata_io** — spatial data container and Xenium loader
- **zarr** — caching and session persistence
- **dask** — lazy array loading
- **tifffile**, **opencv**, **scikit-image** — image processing

## Known Compatibility Patches

- **ICE/X11 disconnect** — handled at startup (lines 33–51 of `02_xenium_viewer.py`)
- **pandas 3.0 PyArrow strings** — `_convert_arrow_strings()` in `01_load_sdata.py`
- **NumPy 2.0** — `np.NAN` fallback for omnipath compatibility

## Version History

See `CHANGELOG.md`. The codebase was refactored from a 4295-line monolith into 11 modular tabs in March 2026.
