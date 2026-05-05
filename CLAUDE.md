# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Xenium Viewer** is a napari-based spatial transcriptomics viewer for Xenium 3.x output data — a Linux equivalent to the commercial Xenium Explorer. It visualizes high-resolution spatial gene expression data with cell-level resolution.

## Running the Viewer

```bash
# First-time setup
conda env create -f environment.yml
conda activate xenium_viewer

# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Without this, transcript loading falls back to scanning transcripts.parquet (~5s/gene instead of <100ms)
xenium-preprocess /path/to/xenium/output/

# Launch viewer (file dialog opens if no path given)
xenium-viewer [/path/to/xenium/output/]

# Launch without SpatialData zarr cache
xenium-viewer /path/to/xenium/output/ --no-cache
```

The package is installed as `xenium-viewer` (PyPI name) / `xenium_viewer` (import name) via `pip install -e .` (handled automatically by `environment.yml`). Console scripts: `xenium-viewer`, `xenium-preprocess`, `xenium-fetch-references`, `xenium-build-custom-segmentation`. You can also run `python -m xenium_viewer ...`.

There is no test suite or CI/CD. All testing is manual/exploratory.

## Architecture

### Entry Points & Load Sequence

1. **`src/xenium_viewer/app.py`** — Main entry point (~1300 lines). Validates data dir, orchestrates loading, builds napari viewer with layers, instantiates all managers, creates `ViewerContext`, builds all tab widgets, then restores session. The `main()` function is the `xenium-viewer` console-script entry point.
2. **`src/xenium_viewer/loader.py`** — Loads SpatialData from Xenium output. Uses a zarr cache (`sdata_cached.zarr/`) for 60–70% faster subsequent launches; staleness detected via `experiment.xenium` mtime. Public API: `load_sdata`, `load_umap`, `load_clusterings`, `get_label_to_obs_mapping`.
3. **`src/xenium_viewer/preprocess.py`** — One-time step that splits the transcript parquet into ~480 per-gene feather files for fast per-gene loading (~100ms vs 4–5s). The `main()` function is the `xenium-preprocess` console-script entry point.

### Central State: ViewerContext

**`src/xenium_viewer/utils/viewer_context.py`** — `ViewerContext` dataclass is the single shared state object passed to all tab modules. It holds:
- Core data: `viewer`, `adata`, `sdata`, `clusterings`, pixel size
- Layer references: `cell_labels_layer`, `transcript_layer`, `roi_layer`
- Manager objects: `CellColorManager`, `TranscriptLoader`, `UMAPViewer`
- Mutable state dicts: `state` (general), `he_state` (H&E registration), `arms_state` (ARMS overlay)
- Shared callables: `record_code()`, `set_status()`, `refresh_clustering_choices()`

### Tab Modules (`src/xenium_viewer/tabs/`)

Each tab follows a consistent pattern:

```python
def build_tab(ctx: ViewerContext) -> tuple[QWidget, dict]:
    # Build UI, register callbacks that access ctx
    # Return (tab_widget, exports_dict)
    # exports_dict may contain "restore_session" callable
```

The 11 tabs cover: Clustering, Cell Coloring, Transcripts, UMAP, ROI Analysis, H&E Registration, Gene Analysis, Ligand-Receptor, Neighborhood Enrichment, Co-occurrence, and ARMS Overlay.

`src/xenium_viewer/tabs/_helpers.py` contains shared utilities (e.g., `StatusProxy`, `make_tab()`).

### Key Utilities (`src/xenium_viewer/utils/`)

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

- **ICE/X11 disconnect** — handled at startup of `src/xenium_viewer/app.py`
- **pandas 3.0 PyArrow strings** — `_convert_arrow_strings()` in `src/xenium_viewer/loader.py`
- **NumPy 2.0** — `np.NAN` fallback for omnipath compatibility

## Version History

See `CHANGELOG.md`. The codebase was refactored from a 4295-line monolith into 11 modular tabs in March 2026.
