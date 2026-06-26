# Co-occur

The Co-occur tab analyses how frequently pairs of cell clusters appear together across a range of spatial distances, producing a co-occurrence score as a function of distance that reveals at what scales clusters are associated or mutually exclusive.

![Co Occurrence](screenshots/tab-co-occurrence.png)

## Controls

| Control | Description |
|---|---|
| Clustering | Clustering to analyse |
| Distance bins | Slider (10–100, default 50) — number of distance bins to compute |
| Run Co-occurrence | Computes the co-occurrence array across all distance bins |
| Cluster checkboxes | Three-column grid — select which source clusters to include in the plot |
| Select All | Checks all cluster checkboxes |
| Deselect All | Unchecks all cluster checkboxes |
| Filter targets | When checked, uses the active cluster selection from Cell Coloring to filter which target clusters appear in the plot |
| Show Co-occurrence Plot | Displays subplots for each selected source cluster showing co-occurrence vs. distance for all target clusters; enabled after running |
| Export CSV... | Saves a table with columns: source cluster, target cluster, distance bin, co-occurrence score; enabled after running |

## Workflow

1. Select a **Clustering** from the dropdown.
2. Set the number of **Distance bins** (more bins give finer resolution but take longer).
3. Click **Run Co-occurrence** and wait for completion.
4. Use the cluster checkboxes to choose which source clusters to visualise, then click **Show Co-occurrence Plot**.
5. Optionally enable **Filter targets** to restrict target clusters to those currently selected in Cell Coloring.
6. Click **Export CSV...** to save the full co-occurrence table.

## Notes

- A co-occurrence score greater than 1 indicates the cluster pair appears together more than expected by chance at that distance; a score less than 1 indicates avoidance.
- Results are persisted to `adata.uns['co_occurrence']` and restored automatically on the next load.
