# PALMS

**P**rovenance-**A**ware **L**inking of **M**ultimodal **S**patial-omics

![CI](https://github.com/rao-genomics-lab/palms/actions/workflows/ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)

A napari-based viewer that brings spatial transcriptomics, histology and genomic
overlays into one coordinate space — and records every action you take as
replayable code. It reads 10x Genomics Xenium 3.x output, and is an open
alternative to the commercial Xenium Explorer, which does not ship a Linux build.
Runs on Linux, macOS and WSL2. Visualises high-resolution spatial gene
expression at cell-level resolution with:

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
  relaunch; reproducible analysis exported as `analysis.py` and
  `analysis_notebook.ipynb`, accumulated across sessions

## Install

```bash
pip install palms
palms /path/to/xenium/output/
```

That is the whole install for using the viewer: PyQt6 carries Qt inside its wheel, so
there is no system Qt or GL package to add. Python 3.10 or newer. Extras are
`pip install "palms[cnv]"` (inferCNV) or `"palms[full]"`.

The one thing a wheel install cannot reach is the **CopyKAT** CNV backend, which needs a
second conda environment for rpy2 + R 4.3. The **inferCNV** backend runs in-process, so
CNV inference itself is available either way.

### From source

To develop PALMS, run the test suite, or use the CopyKAT backend. The viewer pulls in a
heavy scientific stack (napari, scanpy, squidpy, spatialdata, zarr, dask, Qt), and conda
resolves it more reliably than pip:

```bash
git clone https://github.com/rao-genomics-lab/palms.git
cd palms

./scripts/install.sh
conda activate palms
```

The env file installs the core stack via conda-forge and `palms` itself
in editable mode, so source edits in this checkout are picked up immediately.
It also installs the CNV inference stack (`infercnvpy` and `insitucnv-copykat`), so the
CNV tab's **inferCNV** backend works out of the box — no extra step needed.

`install.sh` is just `conda env create -f environment.yml` plus one OS-dependent
step: on **Linux and WSL** it also applies `environment-linux.yml`, which adds
`libglx-devel`. That package fixes a Qt6 startup abort on remote X displays but
is linux-only, so it cannot live in `environment.yml` — with it there, the file
does not solve on macOS at all. If you prefer to run conda yourself:

```bash
conda env create -f environment.yml                              # all platforms
conda env update -n palms -f environment-linux.yml       # Linux/WSL only
```

Skipping the second line on Linux is not silent: the viewer checks for it at
startup and tells you what to run.

### Optional extras

Several features depend on heavier or more niche packages that aren't installed
by default. Add them on top of the core install, after `conda activate
palms`:

```bash
pip install -e ".[celltypist]"   # CellTypist label transfer (Rank Genes tab)
pip install -e ".[r]"            # rpy2-based reference fetcher (needs system R)
pip install -e ".[gpu]"          # torch / torch-geometric / xgboost
pip install -e ".[references]"   # rasterio + readfcs
pip install -e ".[novae]"        # Novae spatial-domain inference (Domains tab)
pip install -e ".[cnv]"          # InSituCNV/infercnvpy — already in environment.yml
pip install -e ".[full]"         # all of the above
```

Each extra is independent — install only the ones you need. A tab whose
optional dependency isn't installed still appears in the UI; it just reports
a clear "not installed" error (with the `pip install` command to run) the
first time you try to use it, instead of failing at startup.

The `cnv` extra is the exception: it is already covered by `environment.yml`, so
you only need it if you installed with plain `pip` instead of conda. It pulls
`insitucnv-copykat` — the [`insituCNV-copykat`](https://github.com/sraorao/insituCNV-copykat)
fork, published under that name because `insitucnv` on PyPI is upstream's. It imports as
`insitucnv` either way, so don't install both. The fork exists because upstream's release
pins `anndata<0.12` and `pandas<3`, which cannot resolve against this app.

### CopyKAT backend (optional second environment)

The CNV tab's inferCNV backend runs in the main environment. Its **CopyKAT**
backend does not: CopyKAT needs rpy2 with R 4.3, a stack that only builds on
python 3.11 and so cannot coexist with the main env's python 3.12. It therefore
lives in a second environment:

```bash
conda env create -f environment-copykat.yml   # creates 'palms_copykat'
```

**Linux only.** This env is not solvable on Apple Silicon: `r-dlm` (a CopyKAT
dependency) comes from the Anaconda `r` channel, which publishes no `osx-arm64`
builds. inferCNV runs in the main env and is unaffected on every platform.

The viewer finds that environment by name and launches CopyKAT there as a
detached background job, which survives the GUI closing. The `copykat` R package
itself is GitHub-only and installs automatically on first run. To point the
viewer at a differently-named environment, set `PALMS_COPYKAT_ENV` (env name) or
`PALMS_COPYKAT_PYTHON` (full path to its python).

Skip this entirely if you only need inferCNV.

### Reinstalling

To wipe the environment and start fresh (e.g. after a dependency conflict or a
major update), remove the old environment, re-clone the repo, and reinstall:

```bash
conda env remove -n palms

git clone https://github.com/rao-genomics-lab/palms.git
cd palms
./scripts/install.sh
conda activate palms
```

## Usage

```bash
# Launch the viewer (file dialog opens if no path given)
palms [/path/to/xenium/output]

# One-time per-gene transcript preprocessing (~30–60 min, ~1 GB output)
# Produces transcript_cache/ next to the dataset, used for fast gene loading.
palms-preprocess /path/to/xenium/output

# Build the SpatialData zarr cache without starting the GUI (tens of minutes).
# The viewer does this on first launch anyway; this is how to get it out of the
# way over ssh or overnight. `--check` reports on an existing cache instead.
palms-build-cache /path/to/xenium/output

# Skip the SpatialData zarr cache (force reload from raw output)
palms /path/to/xenium/output --no-cache
```

Console scripts shipped:

| Command                            | Purpose                                                      |
| ---------------------------------- | ------------------------------------------------------------ |
| `palms`                    | Launch the GUI                                               |
| `palms-preprocess`                | Build the per-gene transcript feather cache (run once)       |
| `palms-build-cache`               | Build/inspect the SpatialData zarr cache without the GUI     |
| `palms-fetch-references`          | Download public scRNA-seq references for label transfer      |
| `palms-build-custom-segmentation` | Build custom segmentation assets from Seurat extract output  |

You can also invoke the package directly: `python -m palms ...`.

## Repo layout

```
src/palms/          # the installable package
├── app.py                  # main GUI entry point (~1800 lines)
├── loader.py               # SpatialData loader with zarr cache
├── preprocess.py           # transcript feather cache builder
├── tabs/                   # 26 control-panel tab modules in 5 groups
│   │                       #   Cells: Clustering, Coloring, Transcripts, UMAP
│   │                       #   Genes: Rank Genes, Markers, Correlation, CNV
│   │                       #   Spatial: ROI DEG, Lig-Rec, Nhood Enrich,
│   │                       #            Co-occur, Domains, Annot Nhood, Annot Dist
│   │                       #   Images: H&E, ARMS, Ext Images, Patches
│   │                       #   Tools: Annotations, Segmentation, Crop Dataset,
│   │                       #          Notebook, Dataset, Cache, Templates
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

Full documentation with screenshots is available on the [GitHub Wiki](https://github.com/rao-genomics-lab/palms/wiki).

- `CLAUDE.md` — architecture overview and developer notes
- `CHANGELOG.md` — release history

## Citing PALMS

If PALMS contributes to work you publish, please cite it. `CITATION.cff` in the repo root
holds the machine-readable metadata; GitHub's "Cite this repository" button reads it
directly. A Zenodo DOI for v1.0.1 is being minted and will be added to that file.

## License

MIT
