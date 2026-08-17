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
| Test suite | `pip install -e ".[test]"` | pytest + nbformat, nbclient and ipykernel, needed to run the notebook-replay test |
| Development | `pip install -e ".[dev]"` | The `test` extra plus ruff, matching what CI runs |

`full` deliberately does *not* include `test` or `dev` — install those explicitly if you intend to run the suite. A plain `pip` install needs Python 3.10 or newer; the conda environment pins 3.12.

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

# Launch ignoring any customised analysis templates, running the shipped ones
xenium-viewer /path/to/xenium/output/ --no-user-templates
```

`python -m xenium_viewer` is equivalent to `xenium-viewer` and accepts the same arguments.

### Environment variables

| Variable | Meaning |
|----------|---------|
| `XENIUM_COPYKAT_ENV` | Name of the conda environment holding the CopyKAT stack |
| `XENIUM_COPYKAT_PYTHON` | Full path to that environment's `python` |
| `XENIUM_VIEWER_TEMPLATE_PATH` | Colon-separated search path replacing the default location for customised analysis templates. Set it to an empty value to disable overrides entirely. See [Templates](Tab-Templates). |

## Console Scripts

The following commands are available after activating the environment:

| Command | Description |
|---------|-------------|
| `xenium-viewer` | Launch the viewer |
| `xenium-preprocess` | Preprocess transcripts for fast per-gene loading |
| `xenium-fetch-references` | Download pre-built label transfer reference datasets |
| `xenium-build-custom-segmentation` | Run the custom cell segmentation pipeline |

All commands accept `--help` for usage details.

The repository also ships `scripts/verify_notebook.py`, which exports a dataset's recorded analysis as a notebook, replays it against the raw Xenium output, and compares the replayed results against what the viewer saved. Use `--dry-run` for a summary in seconds without replaying.

## Troubleshooting

**ICE/X11 disconnect on startup**
Handled automatically at startup. No action required.

**Permission errors on read-only filesystems**
If your Xenium output directory is on a read-only filesystem, copy the dataset to a writable location before launching. Alternatively, use `--no-cache` to skip zarr cache creation, though session persistence and faster loading will be unavailable.

**Corrupt zarr cache**
The viewer tries to repair the cache itself before asking you anything — replaying any interrupted write, clearing debris, rebuilding the metadata index, and restoring elements from its own backup copies. Most damage is resolved there and you never see a dialog.

If the repair does not succeed, a dialog offers three choices:

| Choice | What happens |
|--------|--------------|
| Rebuild and restore my data | The cache is set aside as `sdata_cached_backup_<timestamp>.zarr`, a fresh one is built from the raw files, and your clusterings, ROIs, registrations and session state are copied across. The backup is removed once the rebuild succeeds. |
| Rebuild without restoring | The cache is set aside as `sdata_cached_prev_<timestamp>.zarr` and **kept**, and a fresh empty cache is built. Use this when you suspect the saved analyses themselves are the problem. |
| Quit | Nothing is changed, so you can back the directory up first. |

The viewer never deletes a cache — every path here is a rename. With no display available it raises an error rather than choosing for you. To inspect or repair a cache deliberately rather than at startup, use the [Cache](Tab-Cache) tab; for the whole procedure, see [Recovering a Cache](Tutorial-Recovering-a-Cache).

**Stale zarr cache**
Freshness is decided by a content hash of `experiment.xenium` recorded when the cache was built, so copying a dataset, re-downloading it or restoring it from backup does not make the cache look stale — only a genuine change to the file does. If a change is detected, a dialog at startup offers to rebuild; your existing analyses are preserved.

Caches built before content hashing existed fall back to comparing modification times. That is treated as an uncertain signal: it prompts rather than rebuilding, and choosing to keep the cache records a hash so it stops asking.

**Where the viewer writes**
Besides the zarr cache, a dataset directory accumulates `viewer_cache/` (normalised expression, CNV profiles, cached DEG tables, the analysis provenance graph), `plots/` (saved figures), `transcript_cache/` (per-gene transcript files), `analysis.py`, `analysis_notebook.ipynb`, and a rotating `xenium_viewer.log`. The [Dataset](Tab-Dataset) tab lists all of it with sizes and can delete the regenerable parts.
