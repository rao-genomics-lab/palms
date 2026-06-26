# UMAP

The UMAP tab opens a linked scatter plot window showing cells projected into two-dimensional UMAP space and lets you export a publication-ready UMAP figure coloured by the active cluster assignment. The scatter plot updates automatically whenever you apply a new colouring in the [Coloring](Tab-Cell-Coloring) tab.

![Umap](screenshots/tab-umap.png)

## Controls

| Control | Description |
|---|---|
| Show UMAP Window | Button. Opens a separate floating napari window displaying the UMAP scatter plot, coloured by the currently active clustering. |
| UMAP pt size | Slider (1–50, default 15). Controls point size in the UMAP window. The slider value is divided by 100 to produce a fractional size argument. |
| Save format | Dropdown: **PNG** or **SVG**. Output format for the saved UMAP figure. |
| Save UMAP Plot... | Button. Generates a scanpy UMAP plot using the current clustering and cluster labels, prompts for a save path, and writes the file. |

## Workflow

1. Click **Show UMAP Window** to open the scatter plot. If the window is already open, clicking again brings it to the foreground.
2. Adjust **UMAP pt size** to find a point size that is legible at your screen resolution.
3. Change cluster colouring at any time in the [Coloring](Tab-Cell-Coloring) tab and click **Apply Cell Coloring**; the UMAP window reflects the update immediately.
4. To save a figure, choose a **Save format** and click **Save UMAP Plot...**. Select a destination in the file dialog.

## Notes

- UMAP coordinates are loaded from `analysis/umap/gene_expression_2_components/projection.csv` in the Xenium output directory and stored in `adata.obsm['X_umap']`. If this file is absent, the UMAP tab will not display data.
- A clustering must be active before saving for the figure to be colour-coded. If no clustering is selected, the plot is saved without cluster colours.
- PNG output is saved at 150 dpi. SVG output is saved in vector format and is recommended for publications or figures requiring further editing.
- Cluster label renames applied in [Gene Analysis](Tab-Gene-Analysis) are reflected in the saved figure.
