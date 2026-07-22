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

This installs the `xenium-viewer` package in editable mode along with all required dependencies. It also installs the CNV inference stack (`infercnvpy` and `insitucnv`), so the CNV tab's **inferCNV** backend works immediately — no extra step required.

## Optional Extras

These are the extras declared in `pyproject.toml`. Install them after activating the `xenium_viewer` environment:

| Extra | Install command | What it adds |
|-------|----------------|--------------|
| CellTypist | `pip install -e ".[celltypist]"` | Automated cell type annotation using pre-trained models |
| R integration | `pip install -e ".[r]"` | rpy2 + anndata2ri, for the R-based reference fetcher (needs a system R) |
| GPU support | `pip install -e ".[gpu]"` | torch, torch-geometric, pytorch-lightning, xgboost |
| Reference datasets | `pip install -e ".[references]"` | readfcs + rasterio, for reading reference panel formats |
| CNV inference | `pip install -e ".[cnv]"` | infercnvpy + insitucnv — **already included** in `environment.yml` |
| Everything | `pip install -e ".[full]"` | All of the above |

Each extra is independent — install only what you need. A tab whose optional dependency is missing still appears in the UI; it reports a clear "not installed" error naming the command to run, rather than failing at startup.

The `cnv` extra is the exception to "optional": `environment.yml` already installs it, so you only need it for a plain `pip` install. It pulls `insitucnv` from the [`insituCNV-copykat`](https://github.com/sraorao/insituCNV-copykat) fork, which is public and resolves cleanly against this app's `anndata`/`pandas` versions.

To download pre-built reference panels, use the `xenium-fetch-references` command (see [Console Scripts](#console-scripts) below) — it ships with the package and needs no extra install.

## CopyKAT Backend (Optional Second Environment)

The CNV tab's inferCNV backend runs in the main environment. Its **CopyKAT** backend does not: CopyKAT needs rpy2 with R 4.3, a stack that only builds on python 3.11, which cannot coexist with the main environment's python 3.12. It therefore lives in a separate environment:

```bash
conda env create -f environment-copykat.yml   # creates 'xenium_viewer_copykat'
```

The viewer locates that environment by name and launches CopyKAT there as a detached background job, so the analysis (typically a couple of hours) survives closing the GUI. The `copykat` R package is GitHub-only and installs automatically on first run.

To use a differently-named environment, set one of:

| Variable | Meaning |
|----------|---------|
| `XENIUM_COPYKAT_ENV` | Name of the conda environment to use |
| `XENIUM_COPYKAT_PYTHON` | Full path to that environment's `python` |

Skip this section entirely if you only need inferCNV.

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
