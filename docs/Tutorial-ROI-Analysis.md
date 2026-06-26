# ROI Expression and Differential Expression Analysis

**Prerequisites:** Viewer loaded with a dataset; a gene of interest selected

**Time required:** ~15 minutes

---

## Overview

The ROI workflow lets you draw free-form polygons on the spatial canvas and compare gene expression or cell composition between the enclosed regions. You can calculate per-region statistics for a single gene or run a full differential expression analysis between two or more regions.

---

## Steps

### 1. Colour cells by a gene of interest

Before drawing ROIs, set the active gene so the viewer knows which gene to use for expression statistics.

1. Open the **Cells** group and click the **Coloring** tab.
2. Set **Colour by** to **Gene Expression**.
3. Choose your gene from the gene dropdown.
4. Click **Apply Cell Coloring**.

The cell layer now reflects per-cell expression of that gene, making it easier to draw ROIs around high- or low-expression regions.

<!-- SCREENSHOT: docs/screenshots/tutorial-roi-analysis-step1.png -->

### 2. Activate the ROI layer

In the **napari layer list** on the left, click on the **ROI polygons** Shapes layer to select it. The layer highlights to show it is active.

### 3. Draw ROI polygons

Select the **polygon tool** in the napari toolbar (or press `P`). Click to place each vertex and press **Enter** to close the polygon. Draw at least 2 polygons to enable pairwise comparison.

- For quick rectangular regions, switch to the **rectangle tool** (or press `R`) and click-drag.
- You can draw as many ROIs as needed; the viewer labels them Region 1, Region 2, etc. in order of creation.
- To delete a polygon, select it with the selection tool and press Delete.

<!-- SCREENSHOT: docs/screenshots/tutorial-roi-analysis-step3.png -->

**Tip:** Zoom into the region of interest first, then draw the polygon. For large regions, zoom out so the entire boundary is visible before closing the polygon.

### 4. Calculate expression statistics

1. Open the **Spatial** group and click the **ROI DEG** tab.
2. Click **Calculate Expression**.

The results text shows per-region statistics for the active gene:

| Statistic | Description |
|---|---|
| Mean | Mean expression among cells in the region |
| Median | Median expression |
| Std | Standard deviation |
| N cells | Number of cells inside the polygon |

Pairwise t-test results between all region pairs are shown below the per-region table.

<!-- SCREENSHOT: docs/screenshots/tutorial-roi-analysis-step4.png -->

### 5. Export per-cell expression data (optional)

Click **Export CSV** to save a table of all cells inside all ROIs with columns for cell ID, coordinates, cluster assignment, and expression value for the active gene. This is useful for downstream analysis in R or Python.

### 6. Run differential expression between regions

Click **Run ROI DEG**. The viewer runs the Wilcoxon rank-sum test comparing all cells inside each ROI against cells in all other ROIs. The top 50 differentially expressed genes per region appear in the results text area.

**Optional:** If you want to compare only cells of a specific type within the ROIs:

1. Enable **Filter by cluster**.
2. Select the clustering key and cluster label from the dropdowns that appear.

Only cells matching that cluster label are included in the DEG analysis.

<!-- SCREENSHOT: docs/screenshots/tutorial-roi-analysis-step6.png -->

### 7. Export DEG results

Click **Export DEG CSV...** and choose a save location. The file contains the full ranked gene list for each pairwise region comparison, including log-fold changes and adjusted p-values.

### 8. Save volcano plots

Click **Save Volcano Plot(s)...** and choose a directory. One PNG volcano plot is generated per region pair (e.g. `Region1_vs_Region2.png`). Genes above the significance threshold are labelled with their names.

<!-- SCREENSHOT: docs/screenshots/tutorial-roi-analysis-step8.png -->

---

## Notes

- ROI polygons are saved automatically when you close the viewer and restored on the next launch. You do not need to redraw them each session.
- To clear all ROIs, select all shapes in the ROI polygons layer and press Delete, or use the Edit menu in napari.
- If you have registered an H&E image, you can draw ROIs using the H&E layer as a visual reference while keeping the ROI polygons layer active.
- ROI polygons use the Xenium pixel coordinate system. If you export them and reimport them in another tool, apply the pixel size conversion (`pixel_size` in `experiment.xenium`, typically 0.2125 µm/pixel).

---

## Next steps

- [Tutorial-Annotations](Tutorial-Annotations) — draw named tissue annotations for spatial analysis
- [Tab-ROI-Analysis](Tab-ROI-Analysis) — full reference for the ROI DEG tab
- [Tab-Neighborhood-Enrichment](Tab-Neighborhood-Enrichment) — neighbourhood enrichment analysis between cell types
