# xenium_viewer

![CI](https://github.com/sraorao/xenium_viewer/actions/workflows/ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)

A napari-based spatial transcriptomics viewer for 10x Genomics Xenium 3.x output —
the Linux equivalent of the commercial Xenium Explorer. Visualises high-resolution
spatial gene expression at cell-level resolution with:

- **Cell visualisation** — colour by gene expression, cluster, or metadata; load
  per-gene transcript point clouds and density heatmaps; linked UMAP window
- **Clustering & DEG** — Leiden clustering, rank-genes, dotplots, import/export of
  cluster assignments; ROI-based differential expression
- **Spatial analysis** — neighbourhood enrichment, co-occurrence, ligand-receptor
  (via squidpy/omnipath), spatial domain inference (Novae)
- **Image registration** — landmark-based affine registration for H&E and ARMS
  fluorescence images; external OME-TIFF/SVS loader
- **Annotation tools** — draw/label/export annotation shapes (GeoJSON); annotate
  cells by neighbourhood composition or distance to a reference population
- **Session persistence** — all analyses auto-saved to a zarr cache and restored on
  relaunch; reproducible Python code journal exported to `code.py`

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

Several features depend on heavier or more niche packages that aren't installed
by default. Add them on top of the core install, after `conda activate
xenium_viewer`:

```bash
pip install -e ".[celltypist]"   # CellTypist label transfer (Rank Genes tab)
pip install -e ".[r]"            # rpy2-based reference fetcher (needs system R)
pip install -e ".[gpu]"          # torch / torch-geometric / xgboost
pip install -e ".[references]"   # rasterio + readfcs
pip install -e ".[cnv]"          # InSituCNV/infercnvpy CNV inference (CNV tab)
pip install -e ".[full]"         # all of the above
```

Each extra is independent — install only the ones you need. A tab whose
optional dependency isn't installed still appears in the UI; it just reports
a clear "not installed" error (with the `pip install` command to run) the
first time you try to use it, instead of failing at startup.

**CNV extra note**: the `cnv` extra installs `insitucnv` from the
[`insituCNV-copykat`](https://github.com/sraorao/insituCNV-copykat) fork,
which resolves against this app's `anndata`/`pandas` versions — no
`--no-deps` workaround is needed. (Earlier versions of the fork carried
stale `anndata<0.12`/`pandas<3` upper pins that broke resolution; those
were removed.) The `cnv` extra covers the **inferCNV** backend only —
CopyKAT needs a separate environment, see below.

### Reinstalling

To wipe the environment and start fresh (e.g. after a dependency conflict or a
major update), remove the old environment, re-clone the repo, and reinstall:

```bash
conda env remove -n xenium_viewer

git clone https://github.com/sraorao/xenium_viewer.git
cd xenium_viewer
conda env create -f environment.yml
conda activate xenium_viewer
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
src/xenium_viewer/          # the installable package
├── app.py                  # main GUI entry point (~1300 lines)
├── loader.py               # SpatialData loader with zarr cache
├── preprocess.py           # transcript feather cache builder
├── tabs/                   # 21 control-panel tab modules in 5 groups
│   │                       #   Cells: Clustering, Coloring, Transcripts, UMAP
│   │                       #   Genes: Rank Genes, Markers, Correlation, CNV
│   │                       #   Spatial: ROI DEG, Lig-Rec, Nhood Enrich,
│   │                       #            Co-occur, Domains, Annot Nhood, Annot Dist
│   │                       #   Images: H&E, ARMS, Ext Images, Patches
│   │                       #   Tools: Annotations, Segmentation, Notebook
│   └── _helpers.py         # shared tab utilities (StatusProxy, make_tab, …)
└── utils/
    ├── viewer_context.py   # ViewerContext dataclass — shared state for all tabs
    ├── coloring.py         # CellColorManager with DirectLabelColormap
    ├── gene_analysis.py    # rank genes, normalization, Leiden clustering
    ├── spatial_analysis.py # squidpy spatial analysis wrappers
    ├── registration.py     # landmark-based affine registration
    ├── transcript_index.py # per-gene feather loader
    ├── session.py          # zarr-based session persistence
    ├── adata_persistence.py# AnnData / SpatialData result persistence
    ├── umap_widget.py      # linked UMAP scatter window
    └── …                   # annotation utils, notebook engine, patch I/O, …
scripts/
├── extract_seurat_segmentation.R   # stage-1 R script for custom segmentation
├── capture_screenshots.py          # automated wiki screenshot capture
└── push_to_wiki.sh                 # sync docs/ to GitHub Wiki
reference_datasets/                 # fetched scRNA-seq references + metadata
```

## Documentation

Full documentation with screenshots is available on the [GitHub Wiki](https://github.com/sraorao/xenium_viewer/wiki).

- `CLAUDE.md` — architecture overview and developer notes
- `CHANGELOG.md` — release history

## License

MIT
