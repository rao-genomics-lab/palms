# Transcripts

The Transcripts tab lets you overlay raw transcript locations as colour-coded point clouds (up to 10 genes at once) and compute 2D spatial density heatmaps for any single gene, with optional quality filtering and cell-level subsetting. The overlay and density layers are independent and can be displayed simultaneously.

![Transcripts](screenshots/tab-transcripts.png)

## Controls

### Transcript Overlay

| Control | Description |
|---|---|
| Transcript gene | Dropdown. Gene to add to the display list. |
| Add Gene | Button. Appends the selected gene to the list (maximum 10 genes; duplicates are ignored). |
| Remove Selected | Button. Removes the highlighted gene from the list. |
| Clear All | Button. Removes all genes from the list. |
| Gene list | List widget. Shows all currently selected genes. Each gene is assigned a fixed colour by list position: Yellow, Cyan, Magenta, Orange, Green, Sky Blue, Red, Violet, Pink, Brown. |
| Show transcripts | Checkbox. Whether the next **Apply Transcripts** loads the overlay or clears it. It does not act on its own — see the notes. |
| Min QV | Slider (0–40, default 20). A real filter for the **density** (applied as the transcripts are read); provenance-only for the point **overlay**, whose feather files were filtered when `palms-preprocess` built them. See the notes. |
| Apply Transcripts | Button. Loads transcript coordinates for all genes in the list and renders them as a napari points layer. Nothing in this section takes effect until you press it. |
| Colour legend | Label under the gene list showing the colour assigned to each loaded gene, e.g. `Yellow: EPCAM | Cyan: PTPRC`. |

### Density

| Control | Description |
|---|---|
| Gene (density section) | Dropdown. Gene to use for the density heatmap. |
| Bin size (µm) | Slider (10–500, default 50). Spatial bin size in micrometres. |
| Filter by selected clusters | Checkbox. Restricts the density calculation to transcripts located within cells belonging to the currently visible clusters (requires an active cluster filter in the [Coloring](Tab-Cell-Coloring) tab). |
| Normalise by cells per bin | Checkbox. Divides each bin's transcript count by the number of cells in that bin, producing a per-cell density estimate. |
| Preview the density as I change settings | Checkbox (default on). Redraws a density picture in about a tenth of a second whenever you change the gene, bin size, Min QV or either filter, binning the per-gene index `palms-preprocess` built. It is a preview: it draws into its own layer, says so in the status line, and is **never recorded**. Untick it and nothing is drawn until you press Compute Density. |
| Compute Density | Button. Computes a 2D transcript density histogram into the transcript density image layer. This is the button that runs and records the analysis; its result replaces the preview. |

## Workflow

1. Select a gene from the **Transcript gene** dropdown and click **Add Gene**. Repeat for each additional gene (up to 10).
2. Make sure **Show transcripts** is ticked, then click **Apply Transcripts** to load and display the point overlay.
3. To hide the overlay again, either untick the layer in the napari layer list, or untick **Show transcripts** and press **Apply Transcripts** a second time.
4. For a density heatmap, choose a gene in the density section's **Gene** dropdown, set **Bin size (µm)**, and enable optional filters.
5. With the preview on, the heatmap follows your settings as you change them. It is labelled `transcript_density (PREVIEW - not recorded)` in the layer list.
6. Click **Compute Density** to run and record the analysis. The recorded result replaces the preview, in the `transcript_density` layer, and is what lands in `analysis.py` and the exported notebook.

## Notes

- **Nothing in the overlay section happens until you click Apply Transcripts.** Ticking or unticking **Show transcripts** on its own has no effect — the checkbox is only read when **Apply Transcripts** runs, and it decides whether that run loads the points or clears them. To hide the overlay without a reload, untick the layer in the napari layer list instead.
- **Min QV means two different things, depending on the layer.** For the point **overlay** it is provenance-only: the per-gene feather files were already filtered when `palms-preprocess` built them, so the slider records which threshold produced the points but does not change them. For the **density** it is a real filter, applied as the transcripts are read. It cannot go *below* the threshold `palms-preprocess` used (20 by default): the preview refuses and says so, because returning the qv≥20 rows for a qv≥10 question would be the wrong picture. To go below it, re-run `palms-preprocess` with a lower `--min-qv`.
- Without preprocessing, each gene load scans the full `transcripts.parquet` file, which takes roughly 4–5 seconds per gene. Running `palms-preprocess` once on the dataset creates per-gene feather files and reduces load time to under 100 ms per gene.
- Colour assignment in the point overlay is determined by list position, not gene name. If you remove and re-add a gene, it may receive a different colour. The colour legend under the gene list shows the current assignment.
- **Filter by selected clusters is always per cell.** It keeps transcripts whose assigned cell is in the selected clusters, reading that assignment from the transcripts element itself. It used to fall back to filtering whole *bins* — every transcript in a bin containing a selected cell — when the feather cache had been built before cell assignment was stored, and warned `[WARNING: rebuild transcript cache for cell-level filtering]`. Neither the fallback nor the warning exists any more.
- **The density heatmap is recorded**, so it appears in `analysis.py` and the exported notebook and replays from the raw Xenium output. It reads `sdata.points['transcripts']` rather than the per-gene feather index, which is a viewer artifact a notebook would not have — so the *first* density for a gene takes a few seconds even on a preprocessed dataset. Moving the bin-size slider afterwards is instant: fetching a gene and binning it are separate steps ([transcripts.gene](Analysis-Templates#transcriptsgene), [transcripts.density](Analysis-Templates#transcriptsdensity)), and only the histogram re-runs. The point **overlay** is unaffected — it is display, not analysis, and still uses the fast index.
- **The preview is not the result.** It bins the per-gene index rather than the points element, which is why it is fast, and it records nothing at all — browsing twenty genes adds nothing to the notebook. It is safe to trust visually because the two routes were measured to produce identical histograms, and a test keeps them identical. It refuses rather than showing an approximation when it cannot match what Compute Density would draw: below the index's quality floor, for a gene the index does not hold, for an index built before cell ids were kept when the cluster filter is on, or when the clustering has not been applied yet. Each refusal says which in the status line.
- The density heatmap and transcript point overlay are independent napari layers. Adjusting one does not affect the other.
