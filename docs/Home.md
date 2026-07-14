# Xenium Viewer

Xenium Viewer is an open-source, napari-based viewer for 10x Genomics Xenium spatial transcriptomics data, designed as a Linux alternative to the commercial Xenium Explorer. It provides interactive, high-resolution visualisation of cell-level gene expression data directly from Xenium 3.x output directories, with no data export or conversion required. The project is released under the MIT licence and is currently at version 0.1.0 (alpha).

## Key Features

- Visualise cells coloured by gene expression level or clustering assignment
- Leiden clustering with adjustable resolution, computed directly on the dataset
- Transcript overlay with per-gene point layers and density heatmaps
- Linked UMAP viewer in a separate floating window
- Region-of-interest (ROI) analysis with differential gene expression
- Ligand-receptor interaction analysis (via CellChat/squidpy)
- Neighbourhood enrichment analysis
- Co-occurrence analysis
- Spatial domain inference (Novae integration)
- H&E image registration using landmark-based affine alignment
- ARMS fluorescence image overlay with landmark registration
- Custom segmentation pipeline support
- Annotation tools for labelling cells and regions
- Crop Dataset tool for exporting drawn regions as standalone, independently-openable datasets
- Interactive Notebook tab for running code against the loaded dataset
- Reproducible code journal: every user action generates a Python snippet saved to `code.py`

## Requirements

- Linux
- conda or mamba

## Quick Install

```bash
conda env create -f environment.yml
conda activate xenium_viewer
xenium-viewer /path/to/xenium/output/
```

For full installation instructions, optional extras, and troubleshooting, see the [Installation](Installation) page.

## Navigation

### Getting Started

| Page | Description |
|------|-------------|
| [Installation](Installation) | Full install guide, launch commands, optional extras, troubleshooting |
| [Interface Overview](Interface-Overview) | Canvas, control panel, layers, UMAP window, session persistence |

### Tab Reference

| Tab | Description |
|-----|-------------|
| [Tab-Clustering](Tab-Clustering) | Leiden clustering and cluster assignment |
| [Tab-Cell-Coloring](Tab-Cell-Coloring) | Colour cells by gene expression or metadata |
| [Tab-Transcripts](Tab-Transcripts) | Per-gene transcript point layers and heatmaps |
| [Tab-UMAP](Tab-UMAP) | UMAP dimensionality reduction and linked viewer |
| [Tab-Rank-Genes](Tab-Rank-Genes) | Rank marker genes per cluster |
| [Tab-Markers](Tab-Markers) | Marker gene panels and scoring |
| [Tab-Gene-Correlation](Tab-Gene-Correlation) | Pairwise and module gene correlation |
| [Tab-ROI-Analysis](Tab-ROI-Analysis) | ROI-based differential gene expression |
| [Tab-Ligand-Receptor](Tab-Ligand-Receptor) | Ligand-receptor interaction analysis |
| [Tab-Neighborhood-Enrichment](Tab-Neighborhood-Enrichment) | Spatial neighbourhood enrichment |
| [Tab-Co-occurrence](Tab-Co-occurrence) | Spatial co-occurrence scoring |
| [Tab-Domains](Tab-Domains) | Spatial domain inference with Novae |
| [Tab-Annot-Nhood](Tab-Annot-Nhood) | Neighbourhood-based cell annotation |
| [Tab-Annot-Distance](Tab-Annot-Distance) | Distance-based cell annotation |
| [Tab-HE-Registration](Tab-HE-Registration) | H&E image landmark registration |
| [Tab-ARMS-Overlay](Tab-ARMS-Overlay) | ARMS fluorescence image overlay and registration |
| [Tab-External-Images](Tab-External-Images) | Load and align external image files |
| [Tab-Patches](Tab-Patches) | Tile/patch-based image export |
| [Tab-Annotations](Tab-Annotations) | Draw and label annotation shapes |
| [Tab-Segmentation](Tab-Segmentation) | Custom cell segmentation pipeline |
| [Tab-Crop-Dataset](Tab-Crop-Dataset) | Export drawn regions as standalone datasets |
| [Tab-Notebook](Tab-Notebook) | Interactive notebook and code journal |

### Tutorials

| Tutorial | Description |
|----------|-------------|
| [Tutorial-Getting-Started](Tutorial-Getting-Started) | Load a dataset and explore the interface |
| [Tutorial-Clustering](Tutorial-Clustering) | Run Leiden clustering and inspect marker genes |
| [Tutorial-HE-Registration](Tutorial-HE-Registration) | Align an H&E image to the Xenium coordinate space |
| [Tutorial-ARMS-Overlay](Tutorial-ARMS-Overlay) | Overlay and register an ARMS fluorescence image |
| [Tutorial-ROI-Analysis](Tutorial-ROI-Analysis) | Draw ROIs and run differential gene expression |
| [Tutorial-Annotations](Tutorial-Annotations) | Annotate cells and export results |
