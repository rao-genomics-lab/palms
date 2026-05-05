# xenium_viewer

A napari-based spatial transcriptomics viewer for 10x Genomics Xenium 3.x output —
the Linux equivalent of the commercial Xenium Explorer. Visualises high-resolution
spatial gene expression at cell-level resolution, with cell coloring, transcript
density, ROI / DEG / spatial-stats analysis, H&E and external-image registration,
ARMS / patch overlays, and more.

## Install

The viewer pulls in a heavy scientific stack (napari, scanpy, squidpy, spatialdata,
zarr, dask, Qt). The recommended way to install is via conda:

```bash
git clone https://github.com/sraorao/xenium_viewer.git
cd xenium_viewer

conda env create -f environment.yml
conda activate xenium_viewer
```

The env file installs the core stack via conda-forge and `xenium-viewer` itself
in editable mode, so source edits in this checkout are picked up immediately.

### Optional extras

Add on top of the core install:

```bash
pip install -e ".[celltypist]"   # CellTypist label transfer
pip install -e ".[r]"            # rpy2-based reference fetcher (needs system R)
pip install -e ".[gpu]"          # torch / torch-geometric / xgboost
pip install -e ".[references]"   # rasterio + readfcs
pip install -e ".[full]"         # all of the above
```

## Usage

```bash
# Launch the viewer (file dialog opens if no path given)
xenium-viewer [/path/to/xenium/output]

# One-time per-gene transcript preprocessing (~30–60 min, ~1 GB output)
# Produces transcript_cache/ next to the dataset, used for fast gene loading.
xenium-preprocess /path/to/xenium/output

# Skip the SpatialData zarr cache (force reload from raw output)
xenium-viewer /path/to/xenium/output --no-cache
```

Console scripts shipped:

| Command                            | Purpose                                                      |
| ---------------------------------- | ------------------------------------------------------------ |
| `xenium-viewer`                    | Launch the GUI                                               |
| `xenium-preprocess`                | Build the per-gene transcript feather cache (run once)       |
| `xenium-fetch-references`          | Download public scRNA-seq references for label transfer      |
| `xenium-build-custom-segmentation` | Build custom segmentation assets from Seurat extract output  |

You can also invoke the package directly: `python -m xenium_viewer ...`.

## Repo layout

```
src/xenium_viewer/      # the installable package
├── app.py              # main GUI entry point
├── loader.py           # SpatialData loader with zarr cache
├── preprocess.py       # transcript feather cache builder
├── tabs/               # 22 napari control-panel tabs
├── utils/              # shared utilities (coloring, persistence, registration, ...)
└── scripts/            # standalone CLI utilities
scripts/extract_seurat_segmentation.R   # stage-1 R script (run via Rscript)
reference_datasets/                     # fetched references + metadata sidecars
```

## Documentation

- `CLAUDE.md` — architecture overview and developer notes
- `CHANGELOG.md` — release history

## License

MIT
