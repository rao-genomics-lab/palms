# Clustering

The Clustering tab lets you run Leiden clustering on the cell expression matrix with adjustable graph and resolution parameters, import pre-computed cluster assignments from an external file, and export the current result for downstream use. Clusters produced here are immediately available to all other tabs that accept a clustering input.

![Clustering](screenshots/tab-clustering.png)

## Controls

Each control's tooltip names the template parameter its value lands in, so a caption here can be traced to the corresponding argument in the exported notebook (for example **Neighbours** sets `n_neighbors`).

| Control | Description |
|---|---|
| Neighbours | Slider (5–50, default 15). Number of nearest neighbours used to build the kNN graph. Higher values produce coarser, more stable clusters. |
| Principal components | Slider (10–50, default 40). Number of principal components used when constructing the kNN graph. |
| Resolution | Float spinbox (0.1–5.0, step 0.1, default 1.0). Leiden resolution parameter. Higher values produce more, smaller clusters. |
| Clustering backend | Dropdown (`igraph` / `leidenalg`, default `igraph`). Which implementation of the Leiden algorithm to use. `igraph` is orders of magnitude faster. `leidenalg` is scanpy's historical backend and optimises a different objective, so it produces a different partition — choose it to reproduce an existing scanpy pipeline. The backend also fixes `directed`, which has no control of its own — see the note below. |
| Iterations | Spinbox (-1–100). How many Leiden iterations to run; `-1` iterates until the partition stops improving. Resets to the selected backend's default when you change **Clustering backend** (`2` for igraph, `-1` for leidenalg). |
| Use HVGs only | Checkbox (default unchecked). When checked, restricts the analysis to highly variable genes and enables the **Highly variable genes** slider. |
| Highly variable genes | Slider (500–4000, default 2000). Number of highly variable genes selected when "Use HVGs only" is active. |
| Scale (max_value=10) | Checkbox (default unchecked). Scales the expression matrix to unit variance, capping values at 10. |
| Run Leiden Clustering | Button. Executes Leiden clustering. The result is stored under the key `leiden_{flavor}_r{resolution}` (e.g. `leiden_igraph_r1.0`). |
| Import Clustering... | Button. Opens a file dialog to load cluster assignments from a CSV or TSV file with columns `cell_id` and `group`. If those column names are absent, the first two columns are used in that order. |
| Export Clustering... | Button. Saves the active clustering to a CSV or TSV file. |
| Status text area | Read-only. Displays the clustering key, number of clusters found, and the parameters used to produce the result. |

## Workflow

1. Adjust **Neighbours**, **Principal components**, and **Resolution** to suit your dataset size and the granularity you need.
2. Leave **Clustering backend** at `igraph` unless you are reproducing a pipeline that used `leidenalg`. Because the flavour is part of the result key, you can run both at one resolution and compare them side by side.
3. Optionally enable **Use HVGs only** and set **Highly variable genes** to speed up computation on large datasets.
4. Optionally enable **Scale (max_value=10)** if you want variance-scaled input.
5. Click **Run Leiden Clustering**. Check the status area for the result key and cluster count.
6. Switch to [Coloring](Tab-Cell-Coloring) or [Rank Genes](Tab-Rank-Genes) to use the new clustering.
7. To use an externally computed partition, click **Import Clustering...** and select a CSV/TSV file with `cell_id` and `group` columns. The filename becomes the clustering key.
8. To save results for sharing or downstream scripts, click **Export Clustering...**.

## Notes

- **`directed` follows the backend and is not a control.** The two are not symmetric: `igraph` **raises** `ValueError` on a directed graph, while `leidenalg` merely defaults to one. So there is one illegal combination and it is the one a checkbox would invite. Both values the viewer writes are documentary rather than functional — under `igraph` the argument is validated and then ignored, since the graph is built undirected regardless; under `leidenalg`, `True` is exactly what omitting it would give you. They appear in the recorded cell so it states its assumptions, and so a future scanpy changing its defaults cannot move your clustering without the code showing it.

- The result key is deterministic: re-running with the same flavour and resolution produces the same key and overwrites the previous result.
- Datasets clustered before the flavour picker existed hold keys of the older form `leiden_r{resolution}`. Those still load and appear in every dropdown, but a new run at the same resolution now writes `leiden_igraph_r{resolution}` alongside them rather than replacing them — delete the older key if you don't want both.
- Imported clusterings appear alongside computed ones in every clustering dropdown across the viewer.
- To rename clusters for display and export, use the **Edit Cluster Labels...** button, available in both the [Rank Genes](Tab-Rank-Genes) and [Coloring](Tab-Cell-Coloring) tabs.
- A progress bar under **Run Leiden Clustering** tracks the run, which is dominated by the neighbour-graph construction on a large section.
- The clustering runs from a template, so the exact code — including the fixed random seed — can be read and changed in the [Templates](Tab-Templates) tab, and is recorded verbatim into the exported notebook.
- Cluster assignments produced before the move to scanpy's own Leiden implementation will not match a new run at the same parameters, because the two optimise slightly different objectives. Existing saved clusterings are left exactly as they were.
