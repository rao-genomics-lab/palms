# Markers

The Markers tab lets you visualise expression of a user-defined set of genes across cell clusters by providing a JSON dictionary that maps cluster names to gene lists, then generating dotplots, heatmaps, matrix plots, tracks plots, or correlation matrices with a single button press.

![Markers](screenshots/tab-markers.png)

## Controls

| Control | Description |
|---|---|
| Clustering | Clustering that defines the cell groups shown in all plots |
| Marker genes editor | Multi-line text box — paste a JSON dictionary of the form `{"Cluster A": ["Gene1", "Gene2"], "Cluster B": ["Gene3"]}`. Genes absent from the dataset are dropped, and the status bar names them |
| Save format | Output file format: `PNG` or `SVG` |
| Dotplot | Generates and saves a scanpy dotplot (mean expression + fraction detected) |
| Heatmap | Generates and saves a mean-expression heatmap |
| Matrix plot | Generates and saves a matrix plot |
| Tracks plot | Generates and saves a stacked-violin (tracks) plot |
| Correlation matrix | Generates and saves a gene-gene Pearson correlation matrix |

Each button opens a save-file dialog before rendering.

## Workflow

1. Select the desired **Clustering** from the dropdown.
2. Paste a valid JSON dictionary into the **Marker genes editor**.
3. Choose a **Save format**.
4. Click whichever plot button you need — each is independent and can be run in any order.
5. Confirm the save path in the dialog that appears.

## Notes

- The JSON is validated on each button press; a status message reports any genes that were not found in the dataset.
- The JSON content is stored in session state and restored automatically on the next load.
- Plot types are independent — you can generate any subset in any order without re-running a computation step.
- A valid clustering must be selected first; without one, the run stops with "Select a valid clustering first".
- PNG output is written at a set resolution; SVG is vector and has no dpi setting.
- All five plots are templated steps, so the exact code behind them can be read and changed in the [Templates](Tab-Templates) tab, and each is recorded into the analysis provenance.
