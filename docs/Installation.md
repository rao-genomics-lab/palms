# Installation

## Prerequisites

- **Operating system**: Linux
- **Memory**: 8 GB RAM minimum; 32 GB or more is recommended for large datasets (500k+ cells)
- **Package manager**: [mamba](https://mamba.readthedocs.io/) or [conda](https://docs.conda.io/)

## Standard Install

```bash
git clone https://github.com/sraorao/xenium_viewer.git
cd xenium_viewer
conda env create -f environment.yml
conda activate xenium_viewer
```

This installs the `xenium-viewer` package in editable mode along with all required dependencies.

## Optional Extras

| Extra | Install command | What it adds |
|-------|----------------|--------------|
| CellTypist | `pip install celltypist` | Automated cell type annotation using pre-trained models |
| GPU support | `pip install rapids-singlecell` | GPU-accelerated UMAP and clustering via RAPIDS |
| R integration | Install R + Seurat (see docs) | Custom segmentation pipeline using Seurat-based workflows |
| Reference datasets | `xenium-fetch-references` | Download pre-built label transfer reference panels |
| CNV inference | `pip install -e ".[cnv]"` | Copy-number variation inference using InSituCNV/infercnvpy — see note below |

Install optional extras after activating the `xenium_viewer` conda environment.

**CNV inference note**: `insitucnv`'s published package metadata pins `anndata<0.12`/`pandas<3`, which conflicts with this app's own `anndata`/`pandas` requirements and will make `pip install -e ".[cnv]"` fail to resolve. Install its real dependencies directly instead, then add `insitucnv` itself without letting pip re-resolve its stale pin:

```bash
pip install infercnvpy scvelo
pip install --no-deps insitucnv
```

The pin is stale — `insitucnv` only uses stable AnnData APIs and runs fine against newer `anndata`/`pandas` in practice (verified against `anndata` 0.13 / `pandas` 3.0).

## Reinstalling

If the conda environment becomes broken or corrupted, remove it and start fresh:

```bash
conda env remove -n xenium_viewer
git clone https://github.com/sraorao/xenium_viewer.git
cd xenium_viewer
conda env create -f environment.yml
conda activate xenium_viewer
```

Your analyses stored in the dataset directory (`sdata_cached.zarr/viewer_session/`) are not affected by reinstalling the environment.

## Launching the Viewer

```bash
# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Speeds up per-gene transcript loading from ~5 s to <100 ms
xenium-preprocess /path/to/xenium/output/

# Launch viewer (a file dialog opens if no path is given)
xenium-viewer [/path/to/xenium/output/]

# Launch without building or reading the SpatialData zarr cache
xenium-viewer /path/to/xenium/output/ --no-cache
```

## Console Scripts

The following commands are available after activating the environment:

| Command | Description |
|---------|-------------|
| `xenium-viewer` | Launch the viewer |
| `xenium-preprocess` | Preprocess transcripts for fast per-gene loading |
| `xenium-fetch-references` | Download pre-built label transfer reference datasets |
| `xenium-build-custom-segmentation` | Run the custom cell segmentation pipeline |

All commands accept `--help` for usage details.

## Troubleshooting

**ICE/X11 disconnect on startup**
Handled automatically at startup. No action required.

**Permission errors on read-only filesystems**
If your Xenium output directory is on a read-only filesystem, copy the dataset to a writable location before launching. Alternatively, use `--no-cache` to skip zarr cache creation, though session persistence and faster loading will be unavailable.

**Corrupt zarr cache**
If the viewer detects a corrupt cache at startup, it will prompt you. The old cache is automatically preserved as `sdata_cached_corrupt_<timestamp>.zarr` in the same directory for manual inspection or recovery before a fresh cache is built.

**Stale zarr cache**
If `experiment.xenium` has been modified since the cache was last built (for example, by running Xenium Explorer), a dialog will appear at startup offering to rebuild the cache. Your existing analyses (clusterings, ROIs, registration) are preserved during the rebuild.
