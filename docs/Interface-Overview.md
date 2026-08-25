# Interface Overview

![Interface Overview](screenshots/interface-overview.png)

The viewer consists of a **napari canvas** occupying the left portion of the window and a **control panel** docked on the right. The canvas provides interactive pan, zoom, and layer toggle controls using napari's standard mouse and keyboard bindings. The control panel is organised into five top-level tab groups, each with further sub-tabs along the bottom edge.

## Control Panel Tab Groups

Each top-level tab contains sub-tabs positioned along the bottom of the panel. The **Sub-tab** column below is the label on the tab itself, abbreviated where the tab bar is too narrow for the full name; the **Reference page** column links to the page documenting it. This table is the authority for which name refers to which tab.

### Cells

| Sub-tab | Purpose | Reference page |
|---------|---------|----------------|
| Clustering | Run Leiden clustering, adjust resolution, view cluster assignments | [Clustering](Tab-Clustering) |
| Coloring | Colour cells by gene expression, metadata column, or cluster | [Cell Coloring](Tab-Cell-Coloring) |
| Transcripts | Load per-gene transcript point layers and density heatmaps | [Transcripts](Tab-Transcripts) |
| UMAP | Display the UMAP embedding; open the linked UMAP window | [UMAP](Tab-UMAP) |

### Genes

| Sub-tab | Purpose | Reference page |
|---------|---------|----------------|
| Rank Genes | Rank marker genes per cluster using Wilcoxon or t-test | [Rank Genes](Tab-Rank-Genes) |
| Markers | Score and visualise curated marker gene panels | [Markers](Tab-Markers) |
| Correlation | Compute pairwise gene correlations | [Gene Correlation](Tab-Gene-Correlation) |
| CNV | Infer copy-number variation with the inferCNV and CopyKAT backends | [CNV](Tab-CNV) |

### Spatial

| Sub-tab | Purpose | Reference page |
|---------|---------|----------------|
| ROI DEG | Draw regions of interest and run differential gene expression | [ROI Analysis](Tab-ROI-Analysis) |
| Lig-Rec | Ligand-receptor interaction analysis | [Ligand-Receptor](Tab-Ligand-Receptor) |
| Nhood Enrich | Spatial neighbourhood enrichment | [Neighborhood Enrichment](Tab-Neighborhood-Enrichment) |
| Co-occur | Spatial co-occurrence scoring between cell types | [Co-occurrence](Tab-Co-occurrence) |
| Domains | Spatial domain inference using Novae | [Spatial Domains](Tab-Domains) |
| Annot Nhood | Neighbourhood enrichment around annotation regions | [Annot Nhood](Tab-Annot-Nhood) |
| Annot Dist | Distance from each cell to an annotation region | [Annot Distance](Tab-Annot-Distance) |

### Images

| Sub-tab | Purpose | Reference page |
|---------|---------|----------------|
| H&E | Load and register an H&E brightfield image | [H&E Registration](Tab-HE-Registration) |
| ARMS | Load and register an ARMS fluorescence image | [ARMS Overlay](Tab-ARMS-Overlay) |
| Ext Images | Load and align additional external image files | [External Images](Tab-External-Images) |
| Patches | Overlay tile-level clustering and subclone predictions | [Patches](Tab-Patches) |

### Tools

| Sub-tab | Purpose | Reference page |
|---------|---------|----------------|
| Annotations | Draw, label, and export annotation shapes | [Annotations](Tab-Annotations) |
| Segmentation | Run the custom cell segmentation pipeline | [Segmentation](Tab-Segmentation) |
| Crop Dataset | Export drawn regions as standalone, independently-openable datasets | [Crop Dataset](Tab-Crop-Dataset) |
| Notebook | Review recorded analysis steps and export the notebook | [Notebook](Tab-Notebook) |
| Dataset | Inventory everything on disk and delete viewer-created files | [Dataset](Tab-Dataset) |
| Cache | Check, repair and recover the zarr cache | [Cache](Tab-Cache) |
| Templates | View and edit the code each analysis button runs | [Templates](Tab-Templates) |

## Menu Bar

The viewer adds three menus to napari's own menu bar.

### File

| Item | Description |
|------|-------------|
| Open Dataset... (Ctrl+O) | Close the current dataset and open another. The session is saved first. |
| Preprocess Dataset... | Run the one-time per-gene transcript preprocessing on a dataset, the same work `xenium-preprocess` does from the command line. |

### View

| Item | Description |
|------|-------------|
| Show Xenium Controls (Ctrl+Shift+X) | Show or hide the control panel dock. On by default. |
| Show Minimap | Show or hide the overview minimap. Off, and disabled until a minimap exists. |

### Scale bar

The canvas carries a scale bar reading in **micrometres**, switching to millimetres
as you zoom out — a Xenium section is a few mm across, so both units get used. It is
on from the moment a dataset loads and is napari's own scale bar, so it also responds
to napari's scale-bar settings.

The conversion comes from `pixel_size` in `experiment.xenium`; nothing is assumed.
Every layer is given that scale and labelled in micrometres, which means napari's
world coordinates — the ones the scale bar, the minimap and the camera all read —
are micrometres rather than image pixels. Layer *data* is still in pixels, and so is
every registration affine stored in the zarr; the two are converted at the napari
boundary (`utils/units.py`). Nothing you have registered, drawn or exported changes
position or meaning as a result.

### Preferences

| Item | Description |
|------|-------------|
| Plot format | **PNG** or **SVG** for every automatically saved figure. Defaults to SVG. |
| Plot font size | Base font size for generated figures, from 1 to 20. Defaults to 10. |
| CPU cores | The core budget for parallel analyses — currently CopyKAT only. Choices scale to the machine, with "half" and "all" labelled; defaults to half. |
| Record reproducible code | Whether user actions are recorded into the provenance graph. On by default. |
| Save recorded code... | Write the derived script to a chosen location. |
| Continue from existing code file... | Load a previously written script so a new session appends to it. |

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

Because the notebook is derived by sorting the graph rather than by logging actions in order, it always respects dependencies no matter what order you worked in. Both files open with a preamble containing the imports and data loading, so either can be run independently against the raw Xenium output. Immediately after it comes an `environment` step recording the versions of every relevant package along with the random seeds, so a result can be read against the software that produced it. That step has no dependents, so a version change never marks your analysis stale.

Steps that have no code equivalent — the canvas background colour, which overlays were visible, a crop export — are recorded as **notes** and rendered as markdown rather than as code. This distinguishes them from a step that was simply never recorded: a comment-only code cell runs successfully and does nothing, which made the two indistinguishable.

If any step ran a customised analysis template (see the [Templates](Tab-Templates) tab), the exported notebook opens with a banner naming those steps, and cells you edited by hand in the Notebook tab are marked as such rather than presented as recorded provenance.

You can view the accumulated code, inspect the graph, and export the notebook from the **Notebook** tab in the Tools group. The graph is written to `viewer_cache/prov_graph.json` after every recorded step, and also saved into the session at exit, so an analysis spanning several sessions accumulates into one notebook — and an interrupted session does not lose the steps it had already recorded.

## Session Persistence

Your analyses are automatically saved to `sdata_cached.zarr/viewer_session/` within the dataset directory. The following state is persisted and restored when you reopen the same dataset:

- Clustering assignments and resolution settings
- ROI shapes and DEG results
- H&E and ARMS registration landmarks and affine matrices
- Custom colour assignments and cluster names
- CNV results and the parameters each run used
- External image and patch overlay setup
- The provenance graph behind the Notebook tab

No explicit save action is required; the session is written whenever a relevant action is performed.

Larger derived results are written beside the store rather than inside it, in `<dataset>/viewer_cache/` — normalised expression, CNV profiles, cached DEG tables and the provenance graph. Keeping them out of the zarr store means rebuilding the cache does not discard hours of computation. Generated figures go to `<dataset>/plots/`, and the session log to `<dataset>/xenium_viewer.log`. All of these are listed, with their sizes, in the [Dataset](Tab-Dataset) tab.
