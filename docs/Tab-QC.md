# QC

The QC tab shows the standard Xenium quality-control panel and applies cell and gene filtering cutoffs to the whole session. It sits beside [Segmentation](Tab-Segmentation) because the two answer the same question — *which cells is this analysis about?* — and both change it for every other tab at once.

Nothing is filtered until you press **Apply filter**. A dataset that has never been through this tab is untouched: no node in the provenance graph, no result marked stale, and the exported notebook is exactly what it would have been.

## Controls

Each control's tooltip names the template parameter its value lands in, so a caption here can be traced to the corresponding argument in the exported notebook.

| Control | Description |
|---|---|
| Filter cells | Checkbox (default ticked). Whether the cell cutoff is part of the filter. Unticking it disables **Min transcripts per cell** and removes that step from the recorded code. |
| Min transcripts per cell | Spinbox (1–100000, default 10). Cells with fewer total transcripts than this are dropped. 10 is the conventional Xenium starting point. Sets `min_counts`. |
| Filter genes | Checkbox (default ticked). Whether the gene cutoff is part of the filter. |
| Min cells per gene | Spinbox (1–1000000, default 3). Genes detected in fewer cells than this are dropped, counted over the cells that survive the cell cutoff. Sets `min_cells`. |
| Keep line | Read-only. Updates as you move the spin boxes: how many cells and genes the current cutoffs would keep, and what percentage of cells that is. Exact, not an estimate. |
| QC metrics & plots | Button. Draws the four-panel QC figure and prints the negative-control rates. Does not filter anything. |
| Apply filter | Button. Runs the filter, records it, and re-points the viewer at the surviving cells. Disabled when neither checkbox is ticked. |
| Revert to all cells | Button. Restores every cell and removes the filter step from the notebook. Enabled only while a filter is in force. |
| Status area | Read-only. The control-probe percentages after **QC metrics**, and the cell and gene counts after **Apply**. |

## Workflow

1. Click **QC metrics & plots**. The figure lands in the [Plots](Plots-Window) dock: transcripts per cell, unique genes per cell, segmented cell area, and nucleus ratio. The status area reports the two negative-control rates.
2. Read a cutoff off the first panel — the left tail is the population you are considering dropping.
3. Set **Min transcripts per cell** and watch the keep line. Adjust until the loss looks acceptable.
4. Optionally set **Min cells per gene**. On a targeted panel this usually removes very little.
5. Click **Apply filter**.
6. Carry on in [Clustering](Tab-Clustering), [Rank Genes](Tab-Rank-Genes) or anywhere else. Every expression-based analysis now runs on the filtered cells, and the [Clustering](Tab-Clustering) tab shows a one-line reminder of what is in force.

## Notes

- **Cells the filter drops are no longer coloured in the image.** They keep their outline in the label layer but have no row to take a colour from, so they render transparent. That is the clearest confirmation that the filter did what you asked; it is not cells disappearing from the data.
- **The store always keeps every cell.** The filter is two integers saved with the session, not a rewritten dataset. Results you compute while filtered are written back onto the full table, with a blank for the cells the filter dropped — so nothing is lost, and a cache rebuild cannot silently un-filter anything.
- **Reopening the dataset restores the filter.** It is re-applied before any other tab restores, so the session comes back exactly as you left it: same cells, same clusterings, ready to continue.
- **Applying a filter marks earlier results stale.** They were computed on cells that are no longer bound, so the ⚠ badges in the [Notebook](Tab-Notebook) tab are accurate. Re-run whichever ones you still want; [Tools → Dataset](Tab-Dataset) can clear the rest.
- **Reverting removes the step rather than adding one.** There is no code for "un-filter" — a notebook without a filter simply never filtered — so the node is deleted and anything that depended on it is re-pointed at the preamble and marked stale.
- **Changing the segmentation clears the filter.** The cutoffs were chosen against a cell set that no longer exists. Re-apply them after the swap if you still want them.
- **The metrics need the 10x output's own columns.** The control rates come from `control_probe_counts` and `control_codeword_counts`, and the area panels from `cell_area` and `nucleus_area`. A table built from a custom segmentation may carry neither; the tab says so and draws the two-panel form instead of failing.
- **The filter runs from a template**, so the exact code — including which of the two cutoffs is applied — can be read and changed in the [Templates](Tab-Templates) tab, and is recorded verbatim into the exported notebook. A user override of that template will mark the whole notebook stale on the next launch, because restoring the session re-runs the step and the code has changed.
