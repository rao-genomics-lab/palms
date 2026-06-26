# Transcripts

The Transcripts tab lets you overlay raw transcript locations as colour-coded point clouds (up to 10 genes at once) and compute 2D spatial density heatmaps for any single gene, with optional quality filtering and cell-level subsetting. The overlay and density layers are independent and can be displayed simultaneously.

<!-- SCREENSHOT: docs/screenshots/tab-transcripts.png -->

## Controls

### Transcript Overlay

| Control | Description |
|---|---|
| Transcript gene | Dropdown. Gene to add to the display list. |
| Add Gene | Button. Appends the selected gene to the list (maximum 10 genes; duplicates are ignored). |
| Remove Selected | Button. Removes the highlighted gene from the list. |
| Clear All | Button. Removes all genes from the list. |
| Gene list | List widget. Shows all currently selected genes. Each gene is assigned a fixed colour by list position: Yellow, Cyan, Magenta, Orange, Green, Sky Blue, Red, Violet, Pink, Brown. |
| Show transcripts | Checkbox. Toggles the transcript layer on and off without reloading data. |
| Min QV | Slider (0–40, default 20). Minimum Q-score for transcript inclusion. Transcripts with a quality score below this threshold are excluded. |
| Apply Transcripts | Button. Loads transcript coordinates for all genes in the list and renders them as a napari points layer. |

### Density

| Control | Description |
|---|---|
| Density gene | Dropdown. Gene to use for the density heatmap. |
| Bin size (µm) | Slider (10–500, default 50). Spatial bin size in micrometres. |
| Filter by selected clusters | Checkbox. Restricts the density calculation to transcripts located within cells belonging to the currently visible clusters (requires an active cluster filter in the [Coloring](Tab-Cell-Coloring) tab). |
| Normalise by cells per bin | Checkbox. Divides each bin's transcript count by the number of cells in that bin, producing a per-cell density estimate. |
| Compute Density | Button. Computes a 2D transcript density histogram and adds it as a new napari image layer. |

## Workflow

1. Select a gene from the **Transcript gene** dropdown and click **Add Gene**. Repeat for each additional gene (up to 10).
2. Adjust **Min QV** if you want to tighten or relax quality filtering.
3. Click **Apply Transcripts** to load and display the point overlay. Use **Show transcripts** to toggle visibility without reloading.
4. For a density heatmap, choose a gene in **Density gene**, set **Bin size (µm)**, and enable optional filters.
5. Click **Compute Density**. The heatmap appears as a separate layer and can be adjusted independently in the napari layer list.

## Notes

- Without preprocessing, each gene load scans the full `transcripts.parquet` file, which takes roughly 4–5 seconds per gene. Running `xenium-preprocess` once on the dataset creates per-gene feather files and reduces load time to under 100 ms per gene.
- Colour assignment in the point overlay is determined by list position, not gene name. If you remove and re-add a gene, it may receive a different colour.
- The density heatmap and transcript point overlay are independent napari layers. Adjusting one does not affect the other.
