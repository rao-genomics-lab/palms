# Interface Overview

![Interface Overview](screenshots/interface-overview.png)

The viewer consists of a **napari canvas** occupying the left portion of the window and a **control panel** docked on the right. The canvas provides interactive pan, zoom, and layer toggle controls using napari's standard mouse and keyboard bindings. The control panel is organised into five top-level tab groups, each with further sub-tabs along the bottom edge.

## Control Panel Tab Groups

Each top-level tab contains sub-tabs positioned along the bottom of the panel.

### Cells

| Sub-tab | Purpose |
|---------|---------|
| Clustering | Run Leiden clustering, adjust resolution, view cluster assignments |
| Coloring | Colour cells by gene expression, metadata column, or cluster |
| Transcripts | Load per-gene transcript point layers and density heatmaps |
| UMAP | Compute and display UMAP; open the linked UMAP window |

### Genes

| Sub-tab | Purpose |
|---------|---------|
| Rank Genes | Rank marker genes per cluster using Wilcoxon or t-test |
| Markers | Score and visualise curated marker gene panels |
| Correlation | Compute pairwise and module-level gene correlations |
| CNV | Infer copy-number variation using the InSituCNV method |

### Spatial

| Sub-tab | Purpose |
|---------|---------|
| ROI DEG | Draw regions of interest and run differential gene expression |
| Lig-Rec | Ligand-receptor interaction analysis |
| Nhood Enrich | Spatial neighbourhood enrichment |
| Co-occur | Spatial co-occurrence scoring between cell types |
| Domains | Spatial domain inference using Novae |
| Annot Nhood | Annotate cells based on neighbourhood composition |
| Annot Dist | Annotate cells based on distance to a reference population |

### Images

| Sub-tab | Purpose |
|---------|---------|
| H&E | Load and register an H&E brightfield image |
| ARMS | Load and register an ARMS fluorescence image |
| Ext Images | Load and align additional external image files |
| Patches | Export tiled image patches for downstream use |

### Tools

| Sub-tab | Purpose |
|---------|---------|
| Annotations | Draw, label, and export annotation shapes |
| Segmentation | Run the custom cell segmentation pipeline |
| Crop Dataset | Export drawn regions as standalone, independently-openable datasets |
| Notebook | Interactive code notebook and code journal viewer |

## Layer List

The napari layer panel (top-left of the canvas) contains the following layers, created at startup:

| Layer | Type | Description |
|-------|------|-------------|
| `morphology_focus` | Image | 4-channel morphology image from the Xenium output |
| `cell_labels` | Labels | Cell segmentation masks (used for cell colouring) |
| `nucleus_labels` | Labels | Nucleus segmentation masks |
| `transcripts` | Points | Transcript coordinates; populated when a gene is loaded |
| `ROIs` | Shapes | Region-of-interest shapes drawn in the ROI Analysis tab |
| `annotations` | Shapes | Annotation shapes drawn in the Annotations tab |
| `Crop Regions` | Shapes | Polygons drawn in the Crop Dataset tab, marking regions to export |

Additional layers are added dynamically as you load genes, images, or run analyses.

## Status Bar

A status bar at the bottom of the control panel displays progress messages for long-running operations such as clustering, spatial analysis, and image registration. It shows the current operation name and, where available, a progress indicator.

## UMAP Window

A separate floating napari window can be opened from the **UMAP** sub-tab. This window displays the UMAP scatter plot with cells coloured to match the main canvas. Clicking a cell in the UMAP window highlights the corresponding cell in the main canvas, and vice versa. The UMAP window can be repositioned independently of the main viewer.

## Reproducible Code

Every user action that modifies the dataset or triggers an analysis is recorded as a step in a provenance graph, where each step carries its own code and its dependencies on earlier steps. Two files are written into the dataset directory: `analysis.py`, a flat script derived from the graph, and `analysis_notebook.ipynb`, the same steps as notebook cells.

Because the notebook is derived by sorting the graph rather than by logging actions in order, it always respects dependencies no matter what order you worked in. Both files open with a preamble containing the imports and data loading, so either can be run independently against the raw Xenium output.

You can view the accumulated code, inspect the graph, and export the notebook from the **Notebook** tab in the Tools group. The graph is saved with the session and restored when you reopen the dataset, so an analysis spanning several sessions accumulates into one notebook.

## Session Persistence

Your analyses are automatically saved to `sdata_cached.zarr/viewer_session/` within the dataset directory. The following state is persisted and restored when you reopen the same dataset:

- Clustering assignments and resolution settings
- ROI shapes and DEG results
- H&E and ARMS registration landmarks and affine matrices
- Custom colour assignments

No explicit save action is required; the session is written whenever a relevant action is performed.
