# PALMS

PALMS is an open-source, napari-based viewer for 10x Genomics Xenium spatial transcriptomics data, designed as an open alternative to the commercial Xenium Explorer, which ships no Linux build. It runs on Linux, macOS and WSL2. It provides interactive, high-resolution visualisation of cell-level gene expression data directly from Xenium 3.x output directories, with no data export or conversion required. The project is released under the MIT licence and is currently at version 0.1.0 (alpha).

## Key Features

- Visualise cells coloured by gene expression level or clustering assignment
- Leiden clustering with adjustable resolution, computed directly on the dataset
- Transcript overlay with per-gene point layers and density heatmaps
- Linked UMAP viewer in a separate floating window
- Region-of-interest (ROI) analysis with differential gene expression
- Ligand-receptor interaction analysis (squidpy, over the OmniPath, LigRecExtra and CellPhoneDB databases)
- Neighbourhood enrichment analysis
- Co-occurrence analysis
- Spatial domain inference (Novae integration)
- Copy-number variation (CNV) inference with two backends — inferCNV and CopyKAT
- H&E image registration using landmark-based affine alignment
- ARMS fluorescence image overlay with landmark registration
- Custom segmentation pipeline support
- Annotation tools for labelling cells and regions
- Crop Dataset tool for exporting drawn regions as standalone, independently-openable datasets
- Reproducible analysis: every user action is recorded as a step in a provenance graph, exported as `analysis.py` and `analysis_notebook.ipynb` and carried across sessions
- Customisable analysis templates: read, edit and validate the exact code each analysis button runs, per block, without touching the installation
- Cache health checking and repair, with recovery of user-generated data from a backup store
- On-disk inventory of everything a dataset holds, with selective deletion of viewer-created files

## Requirements

- Linux, macOS, or WSL2 (native Windows is not supported)
- conda or mamba

## Quick Install

```bash
./scripts/install.sh
conda activate palms
palms /path/to/xenium/output/
```

`install.sh` creates the environment from `environment.yml` and, on Linux and WSL, applies the `environment-linux.yml` overlay that a Qt6/GLX fix needs.

For full installation instructions, optional extras, and troubleshooting, see the [Installation](Installation) page.

## Navigation

### Getting Started

| Page | Description |
|------|-------------|
| [Installation](Installation) | Full install guide, launch commands, optional extras, troubleshooting |
| [Interface Overview](Interface-Overview) | Canvas, control panel, layers, UMAP window, session persistence |

### Tab Reference

Listed in the order the tabs appear in the control panel. The name in brackets is the
label on the tab itself, where it is abbreviated to fit.

| Tab | Group | Description |
|-----|-------|-------------|
| [Clustering](Tab-Clustering) | Cells | Leiden clustering and cluster assignment |
| [Cell Coloring](Tab-Cell-Coloring) (Coloring) | Cells | Colour cells by gene expression or metadata |
| [Transcripts](Tab-Transcripts) | Cells | Per-gene transcript point layers and heatmaps |
| [UMAP](Tab-UMAP) | Cells | UMAP dimensionality reduction and linked viewer |
| [Rank Genes](Tab-Rank-Genes) | Genes | Rank marker genes per cluster |
| [Markers](Tab-Markers) | Genes | Marker gene panels and scoring |
| [Gene Correlation](Tab-Gene-Correlation) (Correlation) | Genes | Pairwise gene correlation |
| [CNV](Tab-CNV) | Genes | Copy-number variation inference (inferCNV and CopyKAT) |
| [ROI Analysis](Tab-ROI-Analysis) (ROI DEG) | Spatial | ROI-based differential gene expression |
| [Ligand-Receptor](Tab-Ligand-Receptor) (Lig-Rec) | Spatial | Ligand-receptor interaction analysis |
| [Neighborhood Enrichment](Tab-Neighborhood-Enrichment) (Nhood Enrich) | Spatial | Spatial neighbourhood enrichment |
| [Co-occurrence](Tab-Co-occurrence) (Co-occur) | Spatial | Spatial co-occurrence scoring |
| [Spatial Domains](Tab-Domains) (Domains) | Spatial | Spatial domain inference with Novae |
| [Annot Nhood](Tab-Annot-Nhood) | Spatial | Neighbourhood enrichment around annotation regions |
| [Annot Distance](Tab-Annot-Distance) (Annot Dist) | Spatial | Cell-to-annotation distance analysis |
| [H&E Registration](Tab-HE-Registration) (H&E) | Images | H&E image landmark registration |
| [ARMS Overlay](Tab-ARMS-Overlay) (ARMS) | Images | ARMS fluorescence image overlay and registration |
| [External Images](Tab-External-Images) (Ext Images) | Images | Load and align external image files |
| [Patches](Tab-Patches) | Images | Tile/patch-based overlays |
| [Annotations](Tab-Annotations) | Tools | Draw and label annotation shapes |
| [Segmentation](Tab-Segmentation) | Tools | Custom cell segmentation pipeline |
| [Crop Dataset](Tab-Crop-Dataset) | Tools | Export drawn regions as standalone datasets |
| [Notebook](Tab-Notebook) | Tools | Recorded analysis steps, notebook and script export |
| [Dataset](Tab-Dataset) | Tools | On-disk inventory and selective deletion |
| [Cache](Tab-Cache) | Tools | Cache health check, repair and recovery |
| [Templates](Tab-Templates) | Tools | View and edit the code each analysis button runs |

### Tutorials

| Tutorial | Description |
|----------|-------------|
| [Tutorial-Getting-Started](Tutorial-Getting-Started) | Load a dataset and explore the interface |
| [Tutorial-Clustering](Tutorial-Clustering) | Run Leiden clustering and inspect marker genes |
| [Tutorial-HE-Registration](Tutorial-HE-Registration) | Align an H&E image to the Xenium coordinate space |
| [Tutorial-ARMS-Overlay](Tutorial-ARMS-Overlay) | Overlay and register an ARMS fluorescence image |
| [Tutorial-ROI-Analysis](Tutorial-ROI-Analysis) | Draw ROIs and run differential gene expression |
| [Tutorial-Annotations](Tutorial-Annotations) | Annotate cells and export results |
| [Tutorial-Recovering-a-Cache](Tutorial-Recovering-a-Cache) | Diagnose and repair a damaged zarr cache without losing your work |

### The code underneath

| Page | Description |
|------|-------------|
| [Analysis Templates](Analysis-Templates) | What every analysis actually computes, its contract, and its full default source |
| [API Reference](API-Reference) | Calling the loading, analysis and persistence functions directly from a notebook |
